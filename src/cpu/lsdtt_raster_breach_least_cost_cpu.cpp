// lsdtt_raster_breach_least_cost_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Least-cost depression breaching (Lindsay 2020).
//
//   Lindsay, J.B. (2020) Pit-centric depression removal methods.
//   Geomorphometry 2020 Conference Proceedings.
//
// A port of WhiteboxTools next-gen BreachDepressionsLeastCost, as shipped in
// whitebox_workflows 2.0.6 (jblindsay/whitebox_next_gen,
// crates/wbtools_oss/src/tools/hydrology/mod.rs: detect_strict_pits_with_raise,
// CostNode, breach_depressions_least_cost_core; MIT OR Apache-2.0). The
// algorithm is described in lsdtt_raster_breach_least_cost_cpu.hpp.
//
// Reproducing Whitebox bit for bit (edge_policy = whitebox, explicit
// flat_increment, fill_remainder = false) depends on every detail below:
//   * a double working copy, rounded once to float32 at the end;
//   * Whitebox's neighbour order (NE, E, SE, S, SW, W, NW, N), which decides
//     first-encounter ties; it is NOT lsdtt3's grid::d8_clockwise order;
//   * a port of Rust's std::collections::BinaryHeap, because equal costs are
//     common on flats and the heap's tie order picks the path (std::
//     priority_queue matches it only on libstdc++);
//   * no fused multiply-add: this file is compiled with -ffp-contract=off
//     (CMakeLists.txt), as Rust never contracts.
// The tests in tests/cpp/test_breach_least_cost.cpp pin all of this.
//
// Departures from Whitebox, none of which change a Whitebox-mode result:
//   * validity comes from the raster mask, not from comparing with nodata;
//   * max_dist outside 1..32767 is rejected (Whitebox clamps below at 1, and
//     its i16 path length wraps above 32767);
//   * path length is int32;
//   * edge_policy = outlet treats an off-grid neighbour as Whitebox treats a
//     nodata neighbour.
// Departures that do change the result, and so are off in every comparison
// with Whitebox:
//   * the automatic flat increment is 2 x ULP_f32(max |z|) (auto_flat_increment
//     below), not Whitebox's digit rule, whose ~1e-10 m steps vanish in
//     float32;
//   * the remainder fill is lsdtt3's fill_pits (Priority-Flood with min_slope,
//     or a Barnes et al. 2014b flat resolution), not Whitebox's fill_deps
//     (pit-centric fill with Garbrecht-Martz flats).
// -----------------------------------------------------------------------------

#include "lsdtt_raster_breach_least_cost_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"  // fill_pits (remainder fill)
#include "lsdtt_fill_flat_resolution_reporting_cpu.hpp"  // barnes2014b remainder fill

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace lsdtt {

namespace {

// Whitebox's neighbour order: NE, E, SE, S, SW, W, NW, N (mod.rs DX / DY).
// k_back[k] is the direction from neighbour k back to the centre cell.
constexpr int k_wb_dc[8]   = { 1, 1, 1, 0, -1, -1, -1,  0 };
constexpr int k_wb_dr[8]   = {-1, 0, 1, 1,  1,  0, -1, -1 };
constexpr int k_wb_back[8] = { 4, 5, 6, 7,  0,  1,  2,  3 };

constexpr int k_max_path_steps = 32767;

// -----------------------------------------------------------------------------
// lc_heap -- a port of Rust's std::collections::BinaryHeap for Whitebox's
// CostNode, whose ordering is by cost only, REVERSED (a min-heap), with equal
// costs comparing Equal.
//
//   push: append, then sift_up, stopping as soon as the new element is <= its
//         parent.
//   pop:  swap_remove the root (the last element takes its place), then
//         sift_down_to_bottom (always descend to a leaf, taking the right
//         child when left <= right) and sift_up from there.
//
// In the max-heap's terms "a <= b" for CostNode means b.cost <= a.cost.
// -----------------------------------------------------------------------------
struct lc_cost_node
{
    double      cost;
    std::size_t i;
};

class lc_heap
{
public:
    bool empty() const { return d_.empty(); }

    void push(const lc_cost_node& n)
    {
        d_.push_back(n);
        sift_up(0, d_.size() - 1);
    }

    lc_cost_node pop()
    {
        lc_cost_node item = d_.back();
        d_.pop_back();
        if (!d_.empty())
        {
            std::swap(item, d_[0]);
            sift_down_to_bottom(0);
        }
        return item;
    }

private:
    std::vector<lc_cost_node> d_;

    static bool le(const lc_cost_node& a, const lc_cost_node& b)
    {
        return b.cost <= a.cost;
    }

    void sift_up(std::size_t start, std::size_t pos)
    {
        const lc_cost_node e = d_[pos];
        while (pos > start)
        {
            const std::size_t parent = (pos - 1) / 2;
            if (le(e, d_[parent])) break;
            d_[pos] = d_[parent];
            pos = parent;
        }
        d_[pos] = e;
    }

    void sift_down_to_bottom(std::size_t pos)
    {
        const std::size_t end   = d_.size();
        const std::size_t start = pos;
        const lc_cost_node e    = d_[pos];

        std::size_t child = 2 * pos + 1;
        while (end >= 2 && child <= end - 2)
        {
            if (le(d_[child], d_[child + 1])) child += 1;
            d_[pos] = d_[child];
            pos     = child;
            child   = 2 * pos + 1;
        }
        if (child == end - 1)
        {
            d_[pos] = d_[child];
            pos     = child;
        }
        d_[pos] = e;
        sift_up(start, pos);
    }
};

// The automatic flat increment used when flat_increment <= 0:
//
//   2 x ULP_f32(max(|z|max, 1))
//
// |z|max is the largest |z| over valid, finite cells (0 if there are none).
// ULP_f32(a) is the gap from the float32 value a to the next float32 value
// ABOVE it, nextafterf(a, +inf) - a, so at an exact power of two it is the
// (larger) gap of the binade above. The floor of 1 keeps the increment at
// 2^-22 m (~2.4e-7 m) on rasters whose elevations are all within 1 m of zero,
// instead of descending into subnormal steps; an all-nodata raster gets the
// same value. A flat raster is not a special case.
//
// Why 2 ULP: every breach target pit_z - L * small lies at or below |z|max in
// magnitude (except below a negative minimum), where the float32 spacing is at
// most one ULP, and rounding moves each value by at most half of it. Adjacent
// path values differ by `small` before rounding, so by at least
// small - ULP = ULP > 0 after it: carved channels stay strictly descending once
// stored as float32.
double auto_flat_increment(const std::vector<float>&   z,
                           const std::vector<uint8_t>& valid)
{
    float zabs = 0.0f;
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        if (!valid[i] || !std::isfinite(z[i])) continue;
        zabs = std::max(zabs, std::fabs(z[i]));
    }
    const float a  = std::max(zabs, 1.0f);
    const float up = std::nextafter(a, std::numeric_limits<float>::infinity());
    // At FLT_MAX there is no float above; use the gap below instead.
    const float ulp = std::isfinite(up) ? (up - a) : (a - std::nextafter(a, 0.0f));
    return 2.0 * static_cast<double>(ulp);
}

// Interior cells (all 8 neighbours in the raster and valid) with no strictly
// lower neighbour: what the remainder fill is for.
std::int64_t count_cells_without_descent(const std::vector<float>&   z,
                                         const std::vector<uint8_t>& valid,
                                         std::int64_t nrows, std::int64_t ncols)
{
    std::int64_t n = 0;
    for (std::int64_t r = 1; r < nrows - 1; ++r)
    {
        for (std::int64_t c = 1; c < ncols - 1; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r * ncols + c);
            if (!valid[i]) continue;
            bool interior = true;
            bool lower    = false;
            for (int k = 0; k < 8; ++k)
            {
                const std::size_t j = static_cast<std::size_t>(
                    (r + k_wb_dr[k]) * ncols + (c + k_wb_dc[k]));
                if (!valid[j]) { interior = false; break; }
                if (z[j] < z[i]) lower = true;
            }
            if (interior && !lower) ++n;
        }
    }
    return n;
}

} // namespace

// -----------------------------------------------------------------------------
// parse_breach_lc_edge_policy / breach_lc_edge_policy_name
// -----------------------------------------------------------------------------
breach_lc_edge_policy parse_breach_lc_edge_policy(const std::string& text,
                                                  const std::string& key)
{
    if (text == "outlet")   return breach_lc_edge_policy::outlet;
    if (text == "whitebox") return breach_lc_edge_policy::whitebox;
    throw std::invalid_argument(
        key + " must be one of outlet, whitebox; got '" + text + "'");
}

const char* breach_lc_edge_policy_name(breach_lc_edge_policy policy)
{
    return policy == breach_lc_edge_policy::whitebox ? "whitebox" : "outlet";
}

// -----------------------------------------------------------------------------
// validate_breach_least_cost_options
// -----------------------------------------------------------------------------
void validate_breach_least_cost_options(const breach_least_cost_options& options)
{
    if (options.max_dist < 1 || options.max_dist > k_max_path_steps)
    {
        throw std::invalid_argument(
            "lsdtt::breach_depressions_least_cost: max_dist must be between 1 and " +
            std::to_string(k_max_path_steps) +
            " (maximum breach path length in cells); got " +
            std::to_string(options.max_dist));
    }
    if (std::isnan(options.max_cost))
    {
        throw std::invalid_argument(
            "lsdtt::breach_depressions_least_cost: max_cost is NaN");
    }
    if (!std::isfinite(options.flat_increment))
    {
        throw std::invalid_argument(
            "lsdtt::breach_depressions_least_cost: flat_increment must be finite "
            "(<= 0 selects the automatic value)");
    }
    if (!std::isfinite(options.fill_min_slope) || options.fill_min_slope < 0.0f)
    {
        throw std::invalid_argument(
            "lsdtt::breach_depressions_least_cost: fill_min_slope must be finite "
            "and >= 0; got " + std::to_string(options.fill_min_slope));
    }
    if (options.fill_remainder &&
        options.fill_flat_resolution == flat_resolution_choice::barnes2014b_gradient &&
        !(options.fill_min_slope > 0.0f))
    {
        throw std::invalid_argument(
            "lsdtt::breach_depressions_least_cost: fill_flat_resolution "
            "barnes2014b_gradient needs fill_min_slope > 0 (it sets the gradient "
            "imposed across flats); got " + std::to_string(options.fill_min_slope));
    }
}

// -----------------------------------------------------------------------------
// breach_depressions_least_cost (options)
// -----------------------------------------------------------------------------
lsdtt_raster_cpu breach_depressions_least_cost(
    const lsdtt_raster_cpu&          raster,
    const breach_least_cost_options& options,
    breach_least_cost_diagnostics*   diagnostics)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::breach_depressions_least_cost: invalid raster metadata"
        );
    }
    validate_breach_least_cost_options(options);

    const std::int64_t nrows = meta.n_rows;
    const std::int64_t ncols = meta.n_cols;
    const std::size_t  n     = static_cast<std::size_t>(nrows * ncols);
    const bool edge_is_outlet = options.edge_policy == breach_lc_edge_policy::outlet;

    // Fixed cells and their bank zone (see the options). Empty = none.
    const std::vector<uint8_t>* fixed = options.fixed_cells;
    if (fixed && fixed->size() != n)
    {
        throw std::invalid_argument(
            "lsdtt::breach_depressions_least_cost: fixed_cells must hold one value per cell");
    }
    std::vector<uint8_t> bank;
    if (fixed)
    {
        bank.assign(n, uint8_t{0});
        for (std::int64_t r = 0; r < nrows; ++r)
            for (std::int64_t c = 0; c < ncols; ++c)
            {
                if (!(*fixed)[static_cast<std::size_t>(r * ncols + c)]) continue;
                for (int k = 0; k < 8; ++k)
                {
                    const std::int64_t rn = r + k_wb_dr[k], cn = c + k_wb_dc[k];
                    if (rn < 0 || cn < 0 || rn >= nrows || cn >= ncols) continue;
                    const std::size_t j = static_cast<std::size_t>(rn * ncols + cn);
                    if (!(*fixed)[j]) bank[j] = uint8_t{1};
                }
            }
    }
    auto is_fixed = [&](std::size_t i) { return fixed && (*fixed)[i] != 0; };
    auto in_bank  = [&](std::size_t i) { return fixed && bank[i] != 0; };

    breach_least_cost_diagnostics diag;

    // ---- inputs: validity from the mask, double working copy ----------------
    const auto& src_data = raster.get_data();
    const auto& src_mask = raster.get_mask();

    std::vector<uint8_t> valid(n, uint8_t{0});
    std::vector<double>  z(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (src_mask[i] != 0 && !meta.is_nodata(src_data[i]))
        {
            valid[i] = uint8_t{1};
            z[i]     = static_cast<double>(src_data[i]);
        }
    }

    const double small = (options.flat_increment > 0.0)
                       ? options.flat_increment
                       : auto_flat_increment(src_data, valid);
    diag.flat_increment = small;

    auto in_grid = [&](std::int64_t r, std::int64_t c)
    {
        return r >= 0 && c >= 0 && r < nrows && c < ncols;
    };
    auto at = [&](std::int64_t r, std::int64_t c)
    {
        return static_cast<std::size_t>(r * ncols + c);
    };

    // ---- detect_strict_pits_with_raise: in place, row-major -----------------
    // A pit is raised (or, on a flat floor, lowered) to its lowest neighbour
    // minus `small` as soon as it is found, so the cells scanned after it see
    // the new value.
    std::vector<std::pair<std::size_t, double>> pits;   // (cell, original z)
    for (std::int64_t r = 0; r < nrows; ++r)
    {
        for (std::int64_t c = 0; c < ncols; ++c)
        {
            const std::size_t i = at(r, c);
            if (!valid[i]) continue;
            if (is_fixed(i)) continue;   // fixed cells are never pits
            const double zi = z[i];

            bool   is_pit = true;
            double min_zn = std::numeric_limits<double>::infinity();
            for (int k = 0; k < 8; ++k)
            {
                const std::int64_t rn = r + k_wb_dr[k];
                const std::int64_t cn = c + k_wb_dc[k];
                if (!in_grid(rn, cn) || !valid[at(rn, cn)]) { is_pit = false; break; }
                const double zn = z[at(rn, cn)];
                if (zn < min_zn) min_zn = zn;
                if (zn < zi) { is_pit = false; break; }
            }
            if (is_pit)
            {
                z[i] = min_zn - small;
                pits.emplace_back(i, zi);
            }
        }
    }
    diag.pits_detected = static_cast<std::int64_t>(pits.size());

    // Stable sort, descending original z; processed from the back, so lowest
    // first and, among equal z, the last in scan order first.
    std::stable_sort(pits.begin(), pits.end(),
                     [](const std::pair<std::size_t, double>& a,
                        const std::pair<std::size_t, double>& b)
                     { return a.second > b.second; });

    // ---- breach_depressions_least_cost_core ---------------------------------
    const double diag_dist    = std::sqrt(2.0);
    const double step_dist[8] = { diag_dist, 1.0, diag_dist, 1.0,
                                  diag_dist, 1.0, diag_dist, 1.0 };
    const std::int64_t max_dist = options.max_dist;
    const double max_cost       = options.max_cost;

    // Search state for the whole raster, reset through `touched` after every
    // pit, so a search costs what it explores, not the size of its window.
    std::vector<uint8_t>      encountered(n, uint8_t{0});
    std::vector<int8_t>       backlink(n, int8_t{-1});
    std::vector<std::int32_t> path_length(n, 0);
    std::vector<std::size_t>  touched;

    auto parent_of = [&](std::size_t i)
    {
        const int b = backlink[i];
        const std::int64_t r = static_cast<std::int64_t>(i) / ncols;
        const std::int64_t c = static_cast<std::int64_t>(i) % ncols;
        return at(r + k_wb_dr[b], c + k_wb_dc[b]);
    };

    while (!pits.empty())
    {
        const std::size_t pit = pits.back().first;
        pits.pop_back();

        const double pit_z = z[pit];
        const std::int64_t r0 = static_cast<std::int64_t>(pit) / ncols;
        const std::int64_t c0 = static_cast<std::int64_t>(pit) % ncols;

        // Skip a pit that an earlier breach has given a lower neighbour.
        bool still_pit = true;
        for (int k = 0; k < 8; ++k)
        {
            const std::int64_t rn = r0 + k_wb_dr[k];
            const std::int64_t cn = c0 + k_wb_dc[k];
            if (!in_grid(rn, cn)) continue;
            const std::size_t j = at(rn, cn);
            if (valid[j] && z[j] < pit_z) { still_pit = false; break; }
        }
        if (!still_pit)
        {
            ++diag.pits_skipped_not_pit;
            continue;
        }
        ++diag.pits_searched;

        encountered[pit] = 1;
        touched.push_back(pit);

        lc_heap heap;
        heap.push(lc_cost_node{0.0, pit});

        bool solved      = false;
        bool hit_maxcost = false;
        while (!heap.empty())
        {
            const lc_cost_node node = heap.pop();
            if (node.cost > max_cost) { hit_maxcost = true; break; }

            const std::int64_t length = path_length[node.i];
            if (length > max_dist) continue;

            const double cost1 = z[node.i] - pit_z + (static_cast<double>(length) * small);

            const std::int64_t r = static_cast<std::int64_t>(node.i) / ncols;
            const std::int64_t c = static_cast<std::int64_t>(node.i) % ncols;
            // From the bank zone only a step onto a fixed cell is allowed.
            const bool from_bank = in_bank(node.i);

            // The outlet's backlink always points at node.i, so a solved path
            // is lowered from node.i back to the pit.
            std::size_t carve_from = 0;
            for (int k = 0; k < 8; ++k)
            {
                const std::int64_t rn = r + k_wb_dr[k];
                const std::int64_t cn = c + k_wb_dc[k];
                if (!in_grid(rn, cn))
                {
                    if (!edge_is_outlet) continue;
                    // Off the grid: an outlet, exactly as a nodata neighbour.
                    carve_from = node.i;
                    solved     = true;
                    break;
                }
                const std::size_t ni = at(rn, cn);
                if (encountered[ni] == 1) continue;
                if (fixed)
                {
                    const double zt = pit_z - (static_cast<double>(length + 1) * small);
                    if (is_fixed(ni))
                    {
                        // A fixed cell low enough is an outlet; otherwise it
                        // is impassable. It is never lowered (the carve stops
                        // at node.i).
                        if (valid[ni] && z[ni] > zt) continue;
                    }
                    else if (from_bank && valid[ni])
                    {
                        continue;
                    }
                }

                touched.push_back(ni);
                const std::int64_t next_steps = length + 1;
                path_length[ni] = static_cast<std::int32_t>(next_steps);
                backlink[ni]    = static_cast<int8_t>(k_wb_back[k]);

                const double zn   = z[ni];
                const double zout = pit_z - (static_cast<double>(next_steps) * small);
                if (valid[ni] && zn > zout)
                {
                    const double cost2 = zn - zout;
                    const double next_cost = options.minimize_dist
                        ? node.cost + ((cost1 + cost2) / 2.0) * step_dist[k]
                        : node.cost + cost2;
                    encountered[ni] = 1;
                    if (next_steps <= max_dist) heap.push(lc_cost_node{next_cost, ni});
                }
                else
                {
                    carve_from = node.i;
                    solved     = true;
                    break;
                }
            }

            if (solved)
            {
                // Lower each path cell to pit_z - L * small; the pit (L = 0,
                // no backlink) is left at its raised value.
                std::size_t ti = carve_from;
                for (;;)
                {
                    const double desired =
                        pit_z - static_cast<double>(path_length[ti]) * small;
                    if (z[ti] > desired) z[ti] = desired;
                    if (backlink[ti] < 0) break;
                    ti = parent_of(ti);
                }
                break;
            }
        }

        if (solved)           ++diag.pits_solved;
        else if (hit_maxcost) ++diag.pits_unsolved_max_cost;
        else                  ++diag.pits_unsolved_exhausted;

        for (const std::size_t i : touched)
        {
            encountered[i] = 0;
            backlink[i]    = -1;
            path_length[i] = 0;
        }
        touched.clear();
    }

    // ---- round once to float32 ----------------------------------------------
    const float nodata = static_cast<float>(meta.no_data_value);
    std::vector<float> out(n, nodata);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (valid[i]) out[i] = static_cast<float>(z[i]);
    }
    std::vector<double>().swap(z);

    diag.cells_without_descent = count_cells_without_descent(out, valid, nrows, ncols);

    std::cout << "lsdtt::breach_depressions_least_cost: " << diag.pits_detected
              << " pits found, " << diag.pits_skipped_not_pit
              << " skipped (no longer pits), " << diag.pits_solved << " breached, "
              << (diag.pits_unsolved_exhausted + diag.pits_unsolved_max_cost)
              << " unsolved. flat_increment = " << small << ", edge_policy = "
              << breach_lc_edge_policy_name(options.edge_policy) << "\n";

    lsdtt_raster_cpu breached(std::move(out), std::move(valid), meta);

    if (diag.cells_without_descent > 0)
    {
        std::cout << "lsdtt::breach_depressions_least_cost: "
                  << diag.cells_without_descent << " interior cells without descent";
        if (options.fill_remainder)
        {
            std::cout << " - filling remainder with Priority-Flood, min_slope = "
                      << options.fill_min_slope << ", flat resolution = "
                      << flat_resolution_choice_name(options.fill_flat_resolution)
                      << ".\n";
            if (options.fill_flat_resolution == flat_resolution_choice::inline_gradient)
            {
                // Priority-Flood with min_slope (Barnes et al. 2014); at
                // min_slope 0 it dispatches to the bit-identical flat fill.
                breached = fill_pits(breached, options.fill_min_slope, 0,
                                     priority_flood_queue_backend::binary_heap);
            }
            else
            {
                // Flat fill, then Barnes et al. (2014b) flat resolution, with
                // the same report the CLIs print for their own fill.
                breached = fill_pits_reporting_flat_resolution(
                    breached, options.fill_min_slope, options.fill_flat_resolution, 0);
            }
            diag.remainder_filled = true;
        }
        else
        {
            std::cout << " - left unfilled (fill_remainder = false).\n";
        }
    }

    if (diagnostics != nullptr) *diagnostics = diag;
    return breached;
}

// -----------------------------------------------------------------------------
// breach_depressions_least_cost (positional)
// -----------------------------------------------------------------------------
lsdtt_raster_cpu breach_depressions_least_cost(
    const lsdtt_raster_cpu& raster,
    int   max_dist,
    float max_cost,
    bool  min_dist,
    float flat_increment,
    bool  fill_remainder,
    int   /*n_threads*/)
{
    breach_least_cost_options options;
    options.max_dist       = max_dist;
    options.max_cost       = (max_cost >= (std::numeric_limits<float>::max)())
                           ? std::numeric_limits<double>::infinity()
                           : static_cast<double>(max_cost);
    options.minimize_dist  = min_dist;
    options.flat_increment = static_cast<double>(flat_increment);
    options.fill_remainder = fill_remainder;
    return breach_depressions_least_cost(raster, options, nullptr);
}

} // namespace lsdtt
