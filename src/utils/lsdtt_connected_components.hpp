// lsdtt_connected_components.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Connected-component labelling on a row-major grid.
//
// label_connected_components(n_rows, n_cols, connectivity, is_foreground, ...)
//   Labels every 4- or 8-connected component of the cells for which
//   is_foreground(flat_index) is true.
//
//   - Labels are 1, 2, 3, ... in the order of each component's first cell in a
//     row-major scan; 0 is background. They do not depend on traversal order.
//   - Every component is labelled; a caller that drops small components does
//     its own filtering and renumbering from size[].
//   - size[l] is the cell count of label l (size[0] = 0).
//   - If members is given, members[l] lists the cells of label l in
//     breadth-first order from the seed, neighbours visited in row-major order
//     (grid::d8_rowmajor for 8-connectivity, grid::d4_nwes for 4); members[0] is
//     empty.
//
// Returns the number of components.
//
// Replaced private flood fills in lsdtt_raster_binary_mask_cpu
// (label_binary_components_by_size), lsdtt_slope_break_finder_cpu,
// lsdtt_terrace_edge_tracer_cpu and lsdtt_synthetic_domain_cpu
// (tests/cpp/test_connected_components.cpp keeps them as the reference). The
// band-by-band labeller in lsdtt_reference_field_ladder_cpu, which erases and
// reuses labels of undersized components, stays separate.
//
// Header-only: integer traversal, no floating arithmetic.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <queue>
#include <vector>

#include "lsdtt_grid_neighbours.hpp"

namespace lsdtt
{

template <class IsForeground>
std::int32_t label_connected_components(
    int                                     n_rows,
    int                                     n_cols,
    int                                     connectivity,   // 4 or 8
    IsForeground                            is_foreground,
    std::vector<std::int32_t>&              label,
    std::vector<std::int64_t>&              size,
    std::vector<std::vector<std::size_t>>*  members = nullptr)
{
    const std::size_t n =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);
    label.assign(n, 0);
    size.assign(1, 0);
    if (members)
    {
        members->clear();
        members->emplace_back();
    }

    const bool d8 = (connectivity == 8);
    const int n_nbr = d8 ? 8 : 4;
    const int* dr = d8 ? grid::d8_rowmajor_row : grid::d4_nwes_row;
    const int* dc = d8 ? grid::d8_rowmajor_col : grid::d4_nwes_col;

    std::int32_t next_label = 0;
    std::queue<std::size_t> q;

    for (std::size_t seed = 0; seed < n; ++seed)
    {
        if (label[seed] != 0 || !is_foreground(seed)) continue;

        ++next_label;
        std::int64_t count = 0;
        if (members) members->emplace_back();

        label[seed] = next_label;
        q.push(seed);

        while (!q.empty())
        {
            const std::size_t cur = q.front();
            q.pop();
            ++count;
            if (members) members->back().push_back(cur);

            const int r = static_cast<int>(cur / static_cast<std::size_t>(n_cols));
            const int c = static_cast<int>(cur % static_cast<std::size_t>(n_cols));

            for (int k = 0; k < n_nbr; ++k)
            {
                const int nr = r + dr[k];
                const int nc = c + dc[k];
                if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols) continue;
                const std::size_t ni =
                    static_cast<std::size_t>(nr) * static_cast<std::size_t>(n_cols) +
                    static_cast<std::size_t>(nc);
                if (label[ni] != 0 || !is_foreground(ni)) continue;
                label[ni] = next_label;
                q.push(ni);
            }
        }

        size.push_back(count);
    }

    return next_label;
}

} // namespace lsdtt
