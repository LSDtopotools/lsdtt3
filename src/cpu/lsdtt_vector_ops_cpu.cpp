// lsdtt_vector_ops_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_vector_ops_cpu.hpp"

#include <cmath>
#include <map>
#include <variant>

namespace lsdtt
{

lsdtt_snap_key make_snap_key(double x, double y, double tolerance_metres)
{
    const double tol = (tolerance_metres > 0.0) ? tolerance_metres : 1e-6;
    lsdtt_snap_key key;
    key.qx = static_cast<long long>(std::llround(x / tol));
    key.qy = static_cast<long long>(std::llround(y / tol));
    return key;
}

lsdtt_point_snapper::lsdtt_point_snapper(double tolerance_metres)
    : tol_(tolerance_metres > 0.0 ? tolerance_metres : 0.0),
      cell_(tolerance_metres > 0.0 ? tolerance_metres : 1.0)
{
}

long long lsdtt_point_snapper::cell_of(double v) const
{
    return static_cast<long long>(std::floor(v / cell_));
}

int lsdtt_point_snapper::find(const lsdtt_xy& p) const
{
    const long long ci = cell_of(p.x), cj = cell_of(p.y);
    int best = -1;
    double best_d = 0.0;
    for (long long i = ci - 1; i <= ci + 1; ++i)
    {
        for (long long j = cj - 1; j <= cj + 1; ++j)
        {
            const auto it = cells_.find({i, j});
            if (it == cells_.end()) continue;
            for (const int id : it->second)
            {
                const lsdtt_xy& q = points_[static_cast<std::size_t>(id)];
                const double d = std::hypot(q.x - p.x, q.y - p.y);
                if (d <= tol_ && (best < 0 || d < best_d || (d == best_d && id < best)))
                {
                    best = id;
                    best_d = d;
                }
            }
        }
    }
    return best;
}

int lsdtt_point_snapper::id_for(const lsdtt_xy& p)
{
    const int found = find(p);
    if (found >= 0) return found;
    const int id = static_cast<int>(points_.size());
    points_.push_back(p);
    cells_[{cell_of(p.x), cell_of(p.y)}].push_back(id);
    return id;
}

lsdtt_endpoint_graph build_line_endpoint_graph(
    const lsdtt_vector_layer& lines,
    double                    snap_tolerance_metres)
{
    lsdtt_endpoint_graph graph;
    lsdtt_point_snapper snapper(snap_tolerance_metres > 0.0 ? snap_tolerance_metres : 1e-6);

    auto node_for = [&](const lsdtt_xy& p) -> int
    {
        const int id = snapper.id_for(p);
        graph.n_nodes = static_cast<int>(snapper.size());
        return id;
    };

    for (std::size_t li = 0; li < lines.records.size(); ++li)
    {
        const auto* line =
            std::get_if<lsdtt_line_geom>(&lines.records[li].geometry);
        if (!line || line->points.size() < 2)
        {
            continue;
        }

        const lsdtt_xy& front = line->points.front();
        const lsdtt_xy& back  = line->points.back();

        lsdtt_line_endpoint ep_front;
        ep_front.line_index = li;
        ep_front.at_start   = true;
        ep_front.xy         = front;
        ep_front.node_id    = node_for(front);
        graph.endpoints.push_back(ep_front);

        lsdtt_line_endpoint ep_back;
        ep_back.line_index = li;
        ep_back.at_start   = false;
        ep_back.xy         = back;
        ep_back.node_id    = node_for(back);
        graph.endpoints.push_back(ep_back);
    }

    graph.node_degree.assign(static_cast<std::size_t>(graph.n_nodes), 0);
    for (const auto& ep : graph.endpoints)
    {
        ++graph.node_degree[static_cast<std::size_t>(ep.node_id)];
    }
    for (auto& ep : graph.endpoints)
    {
        ep.degree = graph.node_degree[static_cast<std::size_t>(ep.node_id)];
    }

    return graph;
}

bool crs_wkt_equivalent(const std::string& a, const std::string& b)
{
    if (a == b) return true;       // identical (incl. both empty)
    if (a.empty() || b.empty()) return false;
    return false;
}

} // namespace lsdtt
