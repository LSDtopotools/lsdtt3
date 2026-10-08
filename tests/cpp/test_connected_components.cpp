// test_connected_components.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for label_connected_components
// (src/utils/lsdtt_connected_components.hpp, consolidation C08).
//
// Four private flood fills were replaced. Each is kept below verbatim (the
// raster-level ones reduced to the vectors they read), and each is compared
// with the shared labeller called the way the migrated site now calls it, on
// what the site actually returns:
//
//   binary_mask   label_binary_components_by_size: output labels (as float) and
//                 mask after the size filter and renumbering, D4 and D8
//   slope_break   label_components: the int label image. Under decision D3 it
//                 now counts only valid cells; it must equal the removed fill
//                 run on the grid with invalid cells zeroed
//   terrace       label_centreline_components_d8_: labels AND the per-component
//                 member lists in their breadth-first order
//   synthetic     keep_largest_land_component: the resulting mask
//
// Inputs: random grids (1x1 up to 40x57) at several fill densities, with
// masked, nodata and NaN cells for the raster-level predicates, plus an empty
// grid, an all-foreground grid and a single-cell grid.
//
// Usage: lsdtt-test-connected-components
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "../../src/utils/lsdtt_connected_components.hpp"
#include "../../src/utils/lsdtt_grid_neighbours.hpp"

using namespace lsdtt;

namespace
{

// ============================================================================
// Removed copies (references)
// ============================================================================

// src/cpu/lsdtt_raster_binary_mask_cpu.cpp:1109-1236 (as of 2026-10-07),
// label_binary_components_by_size, reduced to its data/mask/nodata inputs and
// its out_data/out_mask outputs (the raster wrapping and argument checks are
// unchanged in the migrated function).
void old_binary_mask(const std::vector<float>& in_data, const std::vector<std::uint8_t>& in_mask,
                     float nodata_value, int nrows, int ncols, int min_component_pixels, bool d4,
                     std::vector<float>& out_data, std::vector<std::uint8_t>& out_mask)
{
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const float nodata = 0.0f;
    out_data.assign(npix, nodata);
    out_mask.assign(npix, std::uint8_t{0});
    std::vector<std::uint8_t> visited(npix, std::uint8_t{0});

    auto is_valid = [&](std::size_t idx) -> bool
    {
        return in_mask[idx] != 0 &&
               !(in_data[idx] == nodata_value) &&
               std::isfinite(in_data[idx]);
    };

    auto is_foreground = [&](std::size_t idx) -> bool
    {
        return is_valid(idx) && in_data[idx] > 0.5f;
    };

    int next_label = 1;

    for (int row0 = 0; row0 < nrows; ++row0)
    {
        for (int col0 = 0; col0 < ncols; ++col0)
        {
            const std::size_t seed =
                static_cast<std::size_t>(row0) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(col0);

            if (visited[seed] != 0 || !is_foreground(seed))
            {
                continue;
            }

            std::vector<std::size_t> stack;
            std::vector<std::size_t> component;

            stack.push_back(seed);
            visited[seed] = std::uint8_t{1};

            while (!stack.empty())
            {
                const std::size_t idx = stack.back();
                stack.pop_back();
                component.push_back(idx);

                const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
                const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));

                for (int dr = -1; dr <= 1; ++dr)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;
                        if (d4 && dr != 0 && dc != 0) continue;

                        const int rr = row + dr;
                        const int cc = col + dc;

                        if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols)
                        {
                            continue;
                        }

                        const std::size_t nidx =
                            static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) +
                            static_cast<std::size_t>(cc);

                        if (visited[nidx] == 0 && is_foreground(nidx))
                        {
                            visited[nidx] = std::uint8_t{1};
                            stack.push_back(nidx);
                        }
                    }
                }
            }

            if (static_cast<int>(component.size()) >= min_component_pixels)
            {
                const float label = static_cast<float>(next_label);
                ++next_label;

                for (std::size_t idx : component)
                {
                    out_data[idx] = label;
                    out_mask[idx] = std::uint8_t{1};
                }
            }
        }
    }
}

// src/cpu/lsdtt_slope_break_finder_cpu.cpp:413-444 (as of 2026-10-07)
std::vector<int> old_slope_break_label_components(const std::vector<float>& data,
                                                  int nrows, int ncols)
{
    std::vector<int> labels(nrows * ncols, 0);
    int next_label = 1;
    for (int r = 0; r < nrows; ++r)
    for (int c = 0; c < ncols; ++c)
    {
        int idx = r * ncols + c;
        if (data[idx] <= 0.5f || labels[idx] != 0) continue;
        std::queue<std::pair<int,int>> q;
        q.push({r, c});
        labels[idx] = next_label;
        while (!q.empty())
        {
            auto [cr, cc] = q.front(); q.pop();
            for (int dr = -1; dr <= 1; ++dr)
            for (int dc = -1; dc <= 1; ++dc)
            {
                if (!dr && !dc) continue;
                int nr = cr + dr, nc = cc + dc;
                if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                int ni = nr * ncols + nc;
                if (data[ni] <= 0.5f || labels[ni] != 0) continue;
                labels[ni] = next_label;
                q.push({nr, nc});
            }
        }
        ++next_label;
    }
    return labels;
}

// src/cpu/lsdtt_terrace_edge_tracer_cpu.cpp:50-117 (as of 2026-10-07), with
// raster_cell_is_valid(meta, data, mask, i) written out.
std::vector<int> old_terrace_labels(const std::vector<float>& data, const std::vector<std::uint8_t>& mask,
                                    float nodata_value, int nrows, int ncols,
                                    std::vector<std::vector<std::size_t>>& comp_pixels)
{
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    auto valid = [&](std::size_t i) {
        return mask[i] != 0 && !(data[i] == nodata_value) && std::isfinite(data[i]);
    };

    std::vector<int> labels(npix, 0);
    comp_pixels.clear();
    comp_pixels.push_back({}); // label 0 unused

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    static constexpr const auto& DR8 = grid::d8_rowmajor_row;
    static constexpr const auto& DC8 = grid::d8_rowmajor_col;

    int next_label = 1;

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t seed = idx_of(row, col);

            if (!valid(seed)) continue;
            if (data[seed] <= 0.5f) continue;
            if (labels[seed] != 0) continue;

            labels[seed] = next_label;
            comp_pixels.push_back({});
            std::queue<std::size_t> q;
            q.push(seed);

            while (!q.empty())
            {
                const std::size_t cur = q.front();
                q.pop();
                comp_pixels[next_label].push_back(cur);

                const int cr = static_cast<int>(cur / static_cast<std::size_t>(ncols));
                const int cc = static_cast<int>(cur % static_cast<std::size_t>(ncols));

                for (int d = 0; d < 8; ++d)
                {
                    const int nr = cr + DR8[d];
                    const int nc = cc + DC8[d];
                    if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;

                    const std::size_t nidx = idx_of(nr, nc);
                    if (!valid(nidx)) continue;
                    if (data[nidx] <= 0.5f) continue;
                    if (labels[nidx] != 0) continue;

                    labels[nidx] = next_label;
                    q.push(nidx);
                }
            }

            ++next_label;
        }
    }

    return labels;
}

// src/cpu/lsdtt_synthetic_domain_cpu.cpp:29-85 (as of 2026-10-07)
void old_keep_largest_land_component(std::vector<std::uint8_t>& mask, int n_rows, int n_cols)
{
    const std::size_t n = static_cast<std::size_t>(n_rows) * n_cols;

    // zero the one-pixel border so the coastline is internal
    for (int c = 0; c < n_cols; ++c)
    {
        mask[static_cast<std::size_t>(0) * n_cols + c] = 0;
        mask[static_cast<std::size_t>(n_rows - 1) * n_cols + c] = 0;
    }
    for (int r = 0; r < n_rows; ++r)
    {
        mask[static_cast<std::size_t>(r) * n_cols + 0] = 0;
        mask[static_cast<std::size_t>(r) * n_cols + (n_cols - 1)] = 0;
    }

    std::vector<int> label(n, 0);
    int best_label = 0;
    std::size_t best_size = 0;
    int next_label = 0;

    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            const std::size_t s = static_cast<std::size_t>(r) * n_cols + c;
            if (mask[s] == 0 || label[s] != 0) continue;

            ++next_label;
            std::size_t comp_size = 0;
            std::queue<std::pair<int,int>> q;
            q.push({r, c});
            label[s] = next_label;
            while (!q.empty())
            {
                const auto [cr, cc] = q.front(); q.pop();
                ++comp_size;
                static constexpr const auto& dr = grid::d4_nswe_row;
                static constexpr const auto& dc = grid::d4_nswe_col;
                for (int k = 0; k < 4; ++k)
                {
                    const int nr = cr + dr[k], nc = cc + dc[k];
                    if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols) continue;
                    const std::size_t ns = static_cast<std::size_t>(nr) * n_cols + nc;
                    if (mask[ns] == 0 || label[ns] != 0) continue;
                    label[ns] = next_label;
                    q.push({nr, nc});
                }
            }
            if (comp_size > best_size) { best_size = comp_size; best_label = next_label; }
        }
    }

    for (std::size_t s = 0; s < n; ++s)
        mask[s] = (label[s] == best_label) ? 1 : 0;
}

// ============================================================================
// The migrated sites, written exactly as they now appear in the sources
// ============================================================================

void new_binary_mask(const std::vector<float>& in_data, const std::vector<std::uint8_t>& in_mask,
                     float nodata_value, int nrows, int ncols, int min_component_pixels, bool d4,
                     std::vector<float>& out_data, std::vector<std::uint8_t>& out_mask)
{
    const std::size_t npix = static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    out_data.assign(npix, 0.0f);
    out_mask.assign(npix, std::uint8_t{0});
    auto is_foreground = [&](std::size_t idx) -> bool
    {
        return in_mask[idx] != 0 && !(in_data[idx] == nodata_value) &&
               std::isfinite(in_data[idx]) && in_data[idx] > 0.5f;
    };
    std::vector<std::int32_t> label;
    std::vector<std::int64_t> size;
    const std::int32_t n_comp =
        label_connected_components(nrows, ncols, d4 ? 4 : 8, is_foreground, label, size);
    std::vector<float> kept_label(static_cast<std::size_t>(n_comp) + 1, 0.0f);
    int next_label = 1;
    for (std::int32_t l = 1; l <= n_comp; ++l)
        if (size[static_cast<std::size_t>(l)] >= min_component_pixels)
            kept_label[static_cast<std::size_t>(l)] = static_cast<float>(next_label++);
    for (std::size_t i = 0; i < npix; ++i)
    {
        const float k = kept_label[static_cast<std::size_t>(label[i])];
        if (label[i] != 0 && k != 0.0f) { out_data[i] = k; out_mask[i] = std::uint8_t{1}; }
    }
}

// D3 (Simon, 2026-10-05; applied 2026-10-08): the slope-break labelling now
// takes only VALID cells (mask set, not nodata, finite) above 0.5, as
// raster_cell_is_foreground does. Written out here as the migrated site does it.
std::vector<int> new_slope_break(const std::vector<float>& data, const std::vector<std::uint8_t>& mask,
                                 float nodata_value, int nrows, int ncols)
{
    std::vector<std::int32_t> label;
    std::vector<std::int64_t> size;
    label_connected_components(nrows, ncols, 8,
        [&](std::size_t i) {
            return mask[i] != 0 && !(data[i] == nodata_value) && std::isfinite(data[i]) &&
                   data[i] > 0.5f; },
        label, size);
    return std::vector<int>(label.begin(), label.end());
}

std::vector<int> new_terrace(const std::vector<float>& data, const std::vector<std::uint8_t>& mask,
                             float nodata_value, int nrows, int ncols,
                             std::vector<std::vector<std::size_t>>& comp_pixels)
{
    std::vector<std::int32_t> label;
    std::vector<std::int64_t> size;
    label_connected_components(nrows, ncols, 8,
        [&](std::size_t i) {
            return mask[i] != 0 && !(data[i] == nodata_value) && std::isfinite(data[i]) &&
                   data[i] > 0.5f; },
        label, size, &comp_pixels);
    return std::vector<int>(label.begin(), label.end());
}

void new_keep_largest(std::vector<std::uint8_t>& mask, int n_rows, int n_cols)
{
    for (int c = 0; c < n_cols; ++c)
    {
        mask[static_cast<std::size_t>(0) * n_cols + c] = 0;
        mask[static_cast<std::size_t>(n_rows - 1) * n_cols + c] = 0;
    }
    for (int r = 0; r < n_rows; ++r)
    {
        mask[static_cast<std::size_t>(r) * n_cols + 0] = 0;
        mask[static_cast<std::size_t>(r) * n_cols + (n_cols - 1)] = 0;
    }
    std::vector<std::int32_t> label;
    std::vector<std::int64_t> size;
    const std::int32_t n_comp = label_connected_components(n_rows, n_cols, 4,
        [&](std::size_t s) { return mask[s] != 0; }, label, size);
    std::int32_t best_label = 0;
    std::int64_t best_size = 0;
    for (std::int32_t l = 1; l <= n_comp; ++l)
        if (size[static_cast<std::size_t>(l)] > best_size)
        { best_size = size[static_cast<std::size_t>(l)]; best_label = l; }
    const std::size_t n = static_cast<std::size_t>(n_rows) * n_cols;
    for (std::size_t s = 0; s < n; ++s)
        mask[s] = (label[s] == best_label) ? 1 : 0;
}

// ============================================================================
// Harness
// ============================================================================

int g_failures = 0;

void check(const std::string& what, bool ok, long long n_cases, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << " (" << n_cases << " cases)";
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

struct grid_case
{
    int r, c;
    std::vector<float> data;
    std::vector<std::uint8_t> mask;
};

std::vector<grid_case> make_cases(unsigned seed, float nodata)
{
    std::mt19937 rng(seed);
    std::vector<grid_case> out;
    const std::pair<int, int> dims[] = {{1, 1}, {1, 9}, {7, 1}, {5, 5}, {12, 17}, {40, 57}, {23, 31}};
    for (auto [r, c] : dims)
        for (double density : {0.0, 0.15, 0.45, 0.6, 0.85, 1.0})
            for (int rep = 0; rep < 3; ++rep)
            {
                grid_case g{r, c, {}, {}};
                const std::size_t n = static_cast<std::size_t>(r) * c;
                std::uniform_real_distribution<double> u(0.0, 1.0);
                for (std::size_t i = 0; i < n; ++i)
                {
                    float v = (u(rng) < density) ? 1.0f : 0.0f;
                    std::uint8_t m = 1;
                    const double z = u(rng);
                    if (z < 0.04) m = 0;
                    else if (z < 0.07) v = nodata;
                    else if (z < 0.09) v = std::numeric_limits<float>::quiet_NaN();
                    else if (z < 0.10) v = 0.5f;            // exactly the threshold: background
                    g.data.push_back(v);
                    g.mask.push_back(m);
                }
                out.push_back(std::move(g));
            }
    return out;
}

} // namespace

int main()
{
    std::cout << "== C08 connected-component labelling vs the removed flood fills ==\n";

    for (float nodata : {-9999.0f, 1.0f})
    {
        const auto cases = make_cases(nodata < 0 ? 5u : 6u, nodata);
        const std::string tag = " [nodata=" + std::to_string(nodata) + "]";

        for (bool d4 : {false, true})
            for (int min_px : {1, 2, 5})
            {
                long long bad = 0;
                for (const auto& g : cases)
                {
                    std::vector<float> od, nd;
                    std::vector<std::uint8_t> om, nm;
                    old_binary_mask(g.data, g.mask, nodata, g.r, g.c, min_px, d4, od, om);
                    new_binary_mask(g.data, g.mask, nodata, g.r, g.c, min_px, d4, nd, nm);
                    if (od != nd || om != nm) ++bad;
                }
                check(std::string("binary_mask labels and mask, ") + (d4 ? "D4" : "D8") +
                          " min_pixels=" + std::to_string(min_px) + tag,
                      bad == 0, static_cast<long long>(cases.size()));
            }

        long long bad_sb = 0, bad_te = 0, bad_sd = 0, n_sb_changed_by_d3 = 0;
        for (const auto& g : cases)
        {
            // D3: identical to the removed fill run on the grid with every
            // invalid cell (masked, nodata, NaN) set to 0, i.e. the change is
            // exactly the masking and nothing else.
            std::vector<float> sanitised = g.data;
            for (std::size_t i = 0; i < sanitised.size(); ++i)
                if (!(g.mask[i] != 0 && !(g.data[i] == nodata) && std::isfinite(g.data[i])))
                    sanitised[i] = 0.0f;
            if (old_slope_break_label_components(sanitised, g.r, g.c) !=
                new_slope_break(g.data, g.mask, nodata, g.r, g.c))
                ++bad_sb;
            if (old_slope_break_label_components(g.data, g.r, g.c) !=
                new_slope_break(g.data, g.mask, nodata, g.r, g.c))
                ++n_sb_changed_by_d3;

            std::vector<std::vector<std::size_t>> cp_old, cp_new;
            const auto lo = old_terrace_labels(g.data, g.mask, nodata, g.r, g.c, cp_old);
            const auto ln = new_terrace(g.data, g.mask, nodata, g.r, g.c, cp_new);
            if (lo != ln || cp_old != cp_new) ++bad_te;

            std::vector<std::uint8_t> m_old(g.data.size()), m_new;
            for (std::size_t i = 0; i < g.data.size(); ++i) m_old[i] = (g.data[i] > 0.5f) ? 1 : 0;
            m_new = m_old;
            old_keep_largest_land_component(m_old, g.r, g.c);
            new_keep_largest(m_new, g.r, g.c);
            if (m_old != m_new) ++bad_sd;
        }
        check("slope_break label image == old fill on the masked grid (D3); D3 changed " +
                  std::to_string(n_sb_changed_by_d3) + " of them" + tag,
              bad_sb == 0, static_cast<long long>(cases.size()));
        check("terrace labels and member lists in order" + tag, bad_te == 0, static_cast<long long>(cases.size()));
        check("synthetic_domain keep-largest mask (incl. no-land and ties)" + tag, bad_sd == 0,
              static_cast<long long>(cases.size()));
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
