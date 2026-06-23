// lsdtt_vector_ops_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

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

lsdtt_endpoint_graph build_line_endpoint_graph(
    const lsdtt_vector_layer& lines,
    double                    snap_tolerance_metres)
{
    lsdtt_endpoint_graph graph;
    std::map<lsdtt_snap_key, int> node_lookup;

    auto node_for = [&](const lsdtt_xy& p) -> int
    {
        const lsdtt_snap_key key =
            make_snap_key(p.x, p.y, snap_tolerance_metres);
        auto it = node_lookup.find(key);
        if (it != node_lookup.end())
        {
            return it->second;
        }
        const int id = graph.n_nodes++;
        node_lookup.emplace(key, id);
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
