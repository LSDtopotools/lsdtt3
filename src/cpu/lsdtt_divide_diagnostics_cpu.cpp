// lsdtt_divide_diagnostics_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_divide_diagnostics_cpu.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "lsdtt_statistics.hpp"

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

    const std::vector<std::int64_t> anc =
        compute_downstream_ancestor_labels(flow, opts.divide_merge_search_steps);

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

    std::vector<double> steps;
    std::vector<double> gradients;
    std::vector<double> class_steps;
    std::int64_t n_class_higher = 0;

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::int64_t a = flow.flat_index_unchecked(row, col);
            if (!is_valid_for_divide(a, flow, dem_mask, dem_data)) continue;

            for (int k = 0; k < (opts.use_d8_diagonals ? 4 : 2); ++k)
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

                if (anc[static_cast<std::size_t>(a)] == anc[static_cast<std::size_t>(b)])
                {
                    continue;
                }

                const double za = static_cast<double>(dem_data[static_cast<std::size_t>(a)]);
                const double zb = static_cast<double>(dem_data[static_cast<std::size_t>(b)]);
                const double dz = std::fabs(za - zb);
                const bool is_diag = (row != nr) && (col != nc);
                const double sep = is_diag ? diag : dx;

                steps.push_back(dz);
                gradients.push_back((sep > 0.0) ? dz / sep : 0.0);

                if (dz > out.step_max_m)
                {
                    out.step_max_m = dz;
                    out.step_max_flat_a = a;
                    out.step_max_flat_b = b;
                    out.step_max_row = row;
                    out.step_max_col = col;
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
                            ++out.n_class_pairs;
                            class_steps.push_back(dz);

                            const double z_high = a_high_b_low ? za : zb;
                            const double z_low  = a_high_b_low ? zb : za;
                            if (z_high > z_low) ++n_class_higher;
                        }
                    }
                }
            }
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

} // namespace lsdtt
