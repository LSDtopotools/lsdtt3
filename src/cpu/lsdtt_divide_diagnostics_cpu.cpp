// lsdtt_divide_diagnostics_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_divide_diagnostics_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "lsdtt_statistics.hpp"
#include "../utils/lsdtt_thread_utils.hpp"   // resolve_thread_count

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

namespace
{

void validate_divide_inputs(
    const lsdtt_flowinfo& flow,
    const lsdtt_raster_cpu& dem,
    const std::vector<float>* class_raster)
{
    if (!flow.has_topological_order())
    {
        throw std::runtime_error(
            "compute_divide_step_stats: flow has no topological order");
    }

    const std::size_t n = flow.full_grid_size();
    const lsdtt_raster_meta& fmeta = flow.meta();
    const lsdtt_raster_meta& dmeta = dem.meta();

    if (fmeta.n_rows != dmeta.n_rows || fmeta.n_cols != dmeta.n_cols)
    {
        throw std::runtime_error(
            "compute_divide_step_stats: flow/dem shape mismatch");
    }

    if (dem.get_data().size() != n || dem.get_mask().size() != n)
    {
        throw std::runtime_error(
            "compute_divide_step_stats: dem arrays do not match flow grid size");
    }

    if (class_raster != nullptr && class_raster->size() != n)
    {
        throw std::runtime_error(
            "compute_divide_step_stats: class_raster size mismatch");
    }
}

bool is_valid_for_divide(
    std::int64_t flat,
    const lsdtt_flowinfo& flow,
    const std::vector<std::uint8_t>& dem_mask,
    const std::vector<float>& dem_data)
{
    if (!flow.is_valid_cell(flat)) return false;
    const std::size_t idx = static_cast<std::size_t>(flat);
    if (dem_mask[idx] == 0) return false;
    return std::isfinite(static_cast<double>(dem_data[idx]));
}

// -----------------------------------------------------------------------------
// divide_pair_tester
//
// Decides whether an adjacent pair (a, b) sits across a divide, i.e. whether
// their flow paths deliver to different points in the network.
//
// The test is the path-set test of the plan (§3.3.1a): store, for every valid
// cell, its forward path
//
//     P_K[i] = ( i, r(i), r^2(i), ..., r^K(i) )
//
// as a fixed-width (K+1) row of one flat array, and call (a, b) a divide pair
// iff P_K[a] and P_K[b] share no node. Paths that reach an outlet (self-loop or
// invalid receiver) are absorbed there and the remaining entries repeat the
// terminal node, which leaves the intersection test unchanged.
//
// This replaces the older terminal-node test anc_K[a] != anc_K[b]. That test is
// strictly stronger than "the paths have not merged": once two paths merge they
// share a node *sequence* but not a step *index*, so two neighbours entering a
// junction at different distances walk to different K-th nodes forever. At
// K = 1 the two tests agree exactly (anc_1[i] = receiver(i), so equality means
// both pixels deliver to the same cell), so K == 1 keeps the cheap label array.
// -----------------------------------------------------------------------------
class divide_pair_tester
{
public:
    divide_pair_tester(const lsdtt_flowinfo& flow, int merge_search_steps, int n_threads)
        : k_(merge_search_steps)
    {
        if (k_ <= 0)
        {
            throw std::runtime_error(
                "divide_pair_tester: merge_search_steps must be > 0");
        }

        if (k_ == 1)
        {
            anc_ = compute_downstream_ancestor_labels(flow, 1);
            return;
        }

        const std::size_t n = flow.full_grid_size();
        const std::size_t width = static_cast<std::size_t>(k_) + 1;
        const auto& mask = flow.mask_array().get_vector();

        paths_.assign(n * width, -1);

        const int threads = resolve_thread_count(n_threads);
        (void)threads;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
        for (std::int64_t idx = 0; idx < static_cast<std::int64_t>(n); ++idx)
        {
            if (mask[static_cast<std::size_t>(idx)] == 0) continue;

            std::int64_t* row = &paths_[static_cast<std::size_t>(idx) * width];
            std::int64_t cur = idx;
            row[0] = cur;

            for (std::size_t step = 1; step < width; ++step)
            {
                const std::int64_t receiver = flow.receiver_of_unchecked(cur);
                if (receiver >= 0 && receiver != cur && flow.is_valid_cell(receiver))
                {
                    cur = receiver;
                }
                row[step] = cur;
            }
        }
    }

    // True if the flow paths of a and b do not meet within k steps.
    [[nodiscard]] bool is_divide_pair(std::int64_t a, std::int64_t b) const
    {
        if (k_ == 1)
        {
            return anc_[static_cast<std::size_t>(a)] != anc_[static_cast<std::size_t>(b)];
        }

        const std::size_t width = static_cast<std::size_t>(k_) + 1;
        const std::int64_t* pa = &paths_[static_cast<std::size_t>(a) * width];
        const std::int64_t* pb = &paths_[static_cast<std::size_t>(b) * width];

        for (std::size_t i = 0; i < width; ++i)
        {
            const std::int64_t va = pa[i];
            for (std::size_t j = 0; j < width; ++j)
            {
                if (va == pb[j]) return false;
            }
        }

        return true;
    }

private:
    int                       k_ = 1;
    std::vector<std::int64_t> anc_;    // K == 1 fast path
    std::vector<std::int64_t> paths_;  // (K+1) * N, cell-major
};

// -----------------------------------------------------------------------------
// for_each_divide_pair_in_row
//
// The one neighbour walk shared by compute_divide_step_stats and
// make_divide_step_raster. For every valid cell in `row` it visits the forward
// half-neighbourhood (E, S, SE, SW), so each unordered adjacent pair is seen
// exactly once, skips pairs that are directly connected by a receiver link, and
// applies the divide test. Surviving pairs are handed to
//
//     fn(a, b, col, nr, nc, is_diag)
//
// where (nr, nc) is the neighbour of the (row, col) cell holding a.
//
// A row only ever touches rows `row` and `row + 1`, so a caller writing into
// per-pixel storage can parallelise safely with a two-colour (even/odd row)
// sweep.
// -----------------------------------------------------------------------------
template <typename PairFn>
void for_each_divide_pair_in_row(
    int                              row,
    const lsdtt_flowinfo&            flow,
    const std::vector<std::uint8_t>& dem_mask,
    const std::vector<float>&        dem_data,
    const divide_pair_tester&        tester,
    bool                             use_d8_diagonals,
    int                              n_rows,
    int                              n_cols,
    PairFn&&                         fn)
{
    for (int col = 0; col < n_cols; ++col)
    {
        const std::int64_t a = flow.flat_index_unchecked(row, col);
        if (!is_valid_for_divide(a, flow, dem_mask, dem_data)) continue;

        for (int k = 0; k < (use_d8_diagonals ? 4 : 2); ++k)
        {
            int nr = row;
            int nc = col;
            if (k == 0) nc = col + 1;
            if (k == 1) nr = row + 1;
            if (k == 2) { nr = row + 1; nc = col + 1; }
            if (k == 3) { nr = row + 1; nc = col - 1; }

            if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols) continue;

            const std::int64_t b = flow.flat_index_unchecked(nr, nc);
            if (!is_valid_for_divide(b, flow, dem_mask, dem_data)) continue;

            if (flow.receiver_of_unchecked(a) == b ||
                flow.receiver_of_unchecked(b) == a)
            {
                continue;
            }

            if (!tester.is_divide_pair(a, b)) continue;

            const bool is_diag = (row != nr) && (col != nc);
            fn(a, b, col, nr, nc, is_diag);
        }
    }
}

} // namespace

std::vector<std::int64_t> compute_downstream_ancestor_labels(
    const lsdtt_flowinfo& flow,
    int merge_search_steps)
{
    if (merge_search_steps <= 0)
    {
        throw std::runtime_error(
            "compute_downstream_ancestor_labels: merge_search_steps must be > 0");
    }
    if (!flow.has_topological_order())
    {
        throw std::runtime_error(
            "compute_downstream_ancestor_labels: flow has no topological order");
    }

    const std::size_t n = flow.full_grid_size();
    const auto& mask = flow.mask_array().get_vector();

    std::vector<std::int64_t> previous(n, -1);
    std::vector<std::int64_t> current(n, -1);

    for (std::size_t idx = 0; idx < n; ++idx)
    {
        if (mask[idx] == 0) continue;

        const std::int64_t flat = static_cast<std::int64_t>(idx);
        const std::int64_t receiver = flow.receiver_of_unchecked(flat);

        if (receiver >= 0 && flow.is_valid_cell(receiver))
        {
            previous[idx] = receiver;
        }
        else
        {
            previous[idx] = flat;
        }
    }

    for (int step = 2; step <= merge_search_steps; ++step)
    {
        current.assign(n, -1);

        for (std::size_t idx = 0; idx < n; ++idx)
        {
            if (mask[idx] == 0) continue;

            const std::int64_t anc = previous[idx];
            if (anc < 0 || !flow.is_valid_cell(anc))
            {
                current[idx] = static_cast<std::int64_t>(idx);
                continue;
            }

            const std::int64_t receiver = flow.receiver_of_unchecked(anc);
            if (receiver >= 0 && flow.is_valid_cell(receiver))
            {
                current[idx] = receiver;
            }
            else
            {
                current[idx] = anc;
            }
        }

        previous.swap(current);
    }

    return previous;
}

lsdtt_divide_step_stats compute_divide_step_stats(
    const lsdtt_flowinfo&                  flow,
    const lsdtt_raster_cpu&                dem,
    const lsdtt_divide_diagnostic_options& opts,
    const std::vector<float>*              class_raster,
    float                                  class_high,
    float                                  class_low)
{
    if (opts.divide_merge_search_steps <= 0)
    {
        throw std::runtime_error(
            "compute_divide_step_stats: divide_merge_search_steps must be > 0");
    }

    validate_divide_inputs(flow, dem, class_raster);

    const divide_pair_tester tester(
        flow, opts.divide_merge_search_steps, opts.n_threads);

    lsdtt_divide_step_stats out;
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dem_data = dem.get_data();
    const std::vector<std::uint8_t>& dem_mask = dem.get_mask();

    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    const double dx = meta.cell_size;
    const double diag = dx * std::sqrt(2.0);

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(flow.full_grid_size());
         ++flat)
    {
        if (is_valid_for_divide(flat, flow, dem_mask, dem_data))
        {
            ++out.n_valid_cells;
        }
    }

    // The pair loop is parallel over rows. Each row accumulates into its own
    // buffers and the rows are merged back in row order, so the result --
    // including which pair wins the step_max tie -- is identical to the serial
    // sweep and independent of the thread count.
    struct row_accumulator
    {
        std::vector<double> steps;
        std::vector<double> gradients;
        std::vector<double> class_steps;
        std::int64_t        n_class_higher = 0;

        double       step_max_m      = 0.0;
        std::int64_t step_max_flat_a = -1;
        std::int64_t step_max_flat_b = -1;
        int          step_max_row    = -1;
        int          step_max_col    = -1;
    };

    std::vector<row_accumulator> row_acc(static_cast<std::size_t>(n_rows));

    const int threads = resolve_thread_count(opts.n_threads);
    (void)threads;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        row_accumulator& acc = row_acc[static_cast<std::size_t>(row)];

        for_each_divide_pair_in_row(
            row, flow, dem_mask, dem_data, tester, opts.use_d8_diagonals,
            n_rows, n_cols,
            [&](std::int64_t a, std::int64_t b, int col, int, int, bool is_diag)
            {
                const double za = static_cast<double>(dem_data[static_cast<std::size_t>(a)]);
                const double zb = static_cast<double>(dem_data[static_cast<std::size_t>(b)]);
                const double dz = std::fabs(za - zb);
                const double sep = is_diag ? diag : dx;

                acc.steps.push_back(dz);
                acc.gradients.push_back((sep > 0.0) ? dz / sep : 0.0);

                if (dz > acc.step_max_m)
                {
                    acc.step_max_m      = dz;
                    acc.step_max_flat_a = a;
                    acc.step_max_flat_b = b;
                    acc.step_max_row    = row;
                    acc.step_max_col    = col;
                }

                if (class_raster != nullptr)
                {
                    const float ca = (*class_raster)[static_cast<std::size_t>(a)];
                    const float cb = (*class_raster)[static_cast<std::size_t>(b)];

                    if (std::isfinite(static_cast<double>(ca)) &&
                        std::isfinite(static_cast<double>(cb)))
                    {
                        const bool a_high_b_low = (ca == class_high && cb == class_low);
                        const bool b_high_a_low = (cb == class_high && ca == class_low);

                        if (a_high_b_low || b_high_a_low)
                        {
                            acc.class_steps.push_back(dz);

                            const double z_high = a_high_b_low ? za : zb;
                            const double z_low  = a_high_b_low ? zb : za;
                            if (z_high > z_low) ++acc.n_class_higher;
                        }
                    }
                }
            });
    }

    std::vector<double> steps;
    std::vector<double> gradients;
    std::vector<double> class_steps;
    std::int64_t n_class_higher = 0;

    for (const row_accumulator& acc : row_acc)
    {
        steps.insert(steps.end(), acc.steps.begin(), acc.steps.end());
        gradients.insert(gradients.end(), acc.gradients.begin(), acc.gradients.end());
        class_steps.insert(class_steps.end(), acc.class_steps.begin(), acc.class_steps.end());
        n_class_higher += acc.n_class_higher;
        out.n_class_pairs += static_cast<std::int64_t>(acc.class_steps.size());

        if (acc.step_max_m > out.step_max_m)
        {
            out.step_max_m      = acc.step_max_m;
            out.step_max_flat_a = acc.step_max_flat_a;
            out.step_max_flat_b = acc.step_max_flat_b;
            out.step_max_row    = acc.step_max_row;
            out.step_max_col    = acc.step_max_col;
        }
    }

    out.n_divide_pairs = static_cast<std::int64_t>(steps.size());
    if (steps.empty()) return out;

    double sum_steps = 0.0;
    for (const double v : steps) sum_steps += v;

    out.step_mean_m = sum_steps / static_cast<double>(steps.size());
    out.step_median_m = compute_percentile(steps, 50.0);
    out.step_p90_m = compute_percentile(steps, 90.0);
    out.gradient_median = compute_percentile(gradients, 50.0);
    out.gradient_p90 = compute_percentile(gradients, 90.0);

    if (!class_steps.empty())
    {
        out.class_step_median_m = compute_percentile(class_steps, 50.0);
        out.class_higher_fraction =
            static_cast<double>(n_class_higher) /
            static_cast<double>(class_steps.size());
    }

    return out;
}

lsdtt_raster_cpu make_divide_step_raster(
    const lsdtt_flowinfo&                  flow,
    const lsdtt_raster_cpu&                dem,
    const lsdtt_divide_diagnostic_options& opts)
{
    if (opts.divide_merge_search_steps <= 0)
    {
        throw std::runtime_error(
            "make_divide_step_raster: divide_merge_search_steps must be > 0");
    }

    validate_divide_inputs(flow, dem, nullptr);

    const divide_pair_tester tester(
        flow, opts.divide_merge_search_steps, opts.n_threads);

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dem_data = dem.get_data();
    const std::vector<std::uint8_t>& dem_mask = dem.get_mask();

    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float> out_data(flow.full_grid_size(), 0.0f);
    std::vector<std::uint8_t> out_mask = dem_mask;

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(flow.full_grid_size());
         ++flat)
    {
        if (!is_valid_for_divide(flat, flow, dem_mask, dem_data))
        {
            // The DEM mask can accept a cell the flow graph rejects. Clear the
            // mask bit as well as writing nodata, so a consumer trusting the
            // mask never reads the no-data value as an elevation.
            out_data[static_cast<std::size_t>(flat)] = nodata;
            out_mask[static_cast<std::size_t>(flat)] = 0;
        }
    }

    const int threads = resolve_thread_count(opts.n_threads);
    (void)threads;

    // Two-colour sweep: a row writes only into rows `row` and `row + 1`, so the
    // even rows are mutually disjoint and so are the odd rows. The per-pixel
    // reduction is a max, so the merge order does not matter.
    for (int colour = 0; colour < 2; ++colour)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
        for (int row = colour; row < n_rows; row += 2)
        {
            for_each_divide_pair_in_row(
                row, flow, dem_mask, dem_data, tester, opts.use_d8_diagonals,
                n_rows, n_cols,
                [&](std::int64_t a, std::int64_t b, int, int, int, bool)
                {
                    const double za =
                        static_cast<double>(dem_data[static_cast<std::size_t>(a)]);
                    const double zb =
                        static_cast<double>(dem_data[static_cast<std::size_t>(b)]);
                    const float dz = static_cast<float>(std::fabs(za - zb));

                    float& out_a = out_data[static_cast<std::size_t>(a)];
                    float& out_b = out_data[static_cast<std::size_t>(b)];
                    if (dz > out_a) out_a = dz;
                    if (dz > out_b) out_b = dz;
                });
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

} // namespace lsdtt
