// lsdtt_grid_geodesic.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Dijkstra on an 8-connected row-major grid.
//
// grid_dijkstra_d8(n_rows, n_cols, dist, sources, passable, edge_cost,
//                  on_relax, stop_at)
//
//   dist      per-cell cost; the caller sets it to +infinity everywhere and to
//             the starting cost at each source before the call. The cost type
//             (float or double) is the caller's, and all accumulation is done in
//             that type: nd = d + edge_cost(from, to, k).
//   sources   flat indices pushed onto the queue at their dist value.
//   passable(to)              may the search enter cell `to`?
//   edge_cost(from, to, k)    cost of the step; k indexes grid::d8_rowmajor
//                             (diagonal when both offsets are non-zero).
//   on_relax(from, to)        called whenever dist[to] is lowered (record a
//                             predecessor, carry a source id, ...).
//   stop_at(cell)             true ends the search when that cell is popped
//                             (after the stale check, before relaxing).
//
// The queue orders entries by (cost, flat index), so the order cells are
// settled in, and therefore every tie, depends only on the costs and indices,
// not on the order sources were given. Neighbours are relaxed in
// grid::d8_rowmajor order.
//
// Replaced three loops of this shape in lsdtt_raster_binary_mask_cpu (geodesic
// detour distance, nearest-seed propagation) and lsdtt_valley_centreline_cpu
// (centreline path), tests/cpp/test_grid_geodesic.cpp keeping them as the
// reference. Searches with other tie or goal rules (least_cost_path, proximity,
// the skeleton A*) stay separate.
//
// A template, so the arithmetic is the caller's own expression inlined at the
// call site exactly as before.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <queue>
#include <utility>
#include <vector>

#include "lsdtt_grid_neighbours.hpp"

namespace lsdtt
{

template <class Cost, class Passable, class EdgeCost, class OnRelax, class StopAt>
void grid_dijkstra_d8(
    int                              n_rows,
    int                              n_cols,
    std::vector<Cost>&               dist,
    const std::vector<std::int64_t>& sources,
    Passable                         passable,
    EdgeCost                         edge_cost,
    OnRelax                          on_relax,
    StopAt                           stop_at)
{
    using qnode = std::pair<Cost, std::int64_t>;
    std::priority_queue<qnode, std::vector<qnode>, std::greater<qnode>> pq;
    for (const std::int64_t s : sources)
        pq.emplace(dist[static_cast<std::size_t>(s)], s);

    const int* dr = grid::d8_rowmajor_row;
    const int* dc = grid::d8_rowmajor_col;

    while (!pq.empty())
    {
        const auto [d, flat] = pq.top();
        pq.pop();
        const std::size_t fi = static_cast<std::size_t>(flat);
        if (d > dist[fi]) continue;
        if (stop_at(flat)) break;
        const int r = static_cast<int>(flat / n_cols);
        const int c = static_cast<int>(flat % n_cols);
        for (int k = 0; k < 8; ++k)
        {
            const int rr = r + dr[k], cc = c + dc[k];
            if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                continue;
            const std::size_t nf = static_cast<std::size_t>(rr) * n_cols + cc;
            if (!passable(nf)) continue;
            const Cost nd = d + edge_cost(fi, nf, k);
            if (nd < dist[nf])
            {
                dist[nf] = nd;
                on_relax(fi, nf);
                pq.emplace(nd, static_cast<std::int64_t>(nf));
            }
        }
    }
}

} // namespace lsdtt
