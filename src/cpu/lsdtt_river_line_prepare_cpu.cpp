// lsdtt_river_line_prepare.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_river_line_prepare.hpp for the pipeline. GDAL-free.
// -----------------------------------------------------------------------------

#include "lsdtt_river_line_prepare.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace lsdtt
{

namespace
{

// Smallest distance treated as "the same point" when the tolerance is zero.
constexpr double k_coincident = 1e-9;

double dist(const lsdtt_xy& a, const lsdtt_xy& b)
{
    return std::hypot(a.x - b.x, a.y - b.y);
}

struct union_find
{
    std::vector<int> parent;
    explicit union_find(int n) : parent(static_cast<size_t>(n))
    {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int a)
    {
        while (parent[static_cast<size_t>(a)] != a)
        {
            parent[static_cast<size_t>(a)] = parent[static_cast<size_t>(parent[static_cast<size_t>(a)])];
            a = parent[static_cast<size_t>(a)];
        }
        return a;
    }
    void unite(int a, int b)
    {
        a = find(a);
        b = find(b);
        if (a != b) parent[static_cast<size_t>(std::max(a, b))] = std::min(a, b);
    }
};

// Uniform grid hash of points.
struct point_grid
{
    double cell;
    std::unordered_map<std::int64_t, std::vector<int>> cells;

    explicit point_grid(double cell_size) : cell(cell_size) {}

    static std::int64_t key(std::int64_t i, std::int64_t j)
    {
        return (i << 32) ^ (j & 0xffffffffLL);
    }
    std::int64_t ix(double v) const { return static_cast<std::int64_t>(std::floor(v / cell)); }

    void insert(const lsdtt_xy& p, int id) { cells[key(ix(p.x), ix(p.y))].push_back(id); }

    // Every id in the cells overlapping the box [x0,x1] x [y0,y1].
    template <typename F>
    void visit(double x0, double y0, double x1, double y1, F f) const
    {
        for (std::int64_t i = ix(x0); i <= ix(x1); ++i)
        {
            for (std::int64_t j = ix(y0); j <= ix(y1); ++j)
            {
                auto it = cells.find(key(i, j));
                if (it == cells.end()) continue;
                for (int id : it->second) f(id);
            }
        }
    }
};

// A part in the working pipeline.
struct work_part
{
    std::vector<lsdtt_xy> v;
    size_t feature = 0;
    bool cut_start = false;
    bool cut_end = false;
    std::string start_id;
    std::string end_id;
    int from_node = -1;
    int to_node = -1;
    int n_splits = 0;
};

std::vector<lsdtt_xy> drop_repeated_vertices(const std::vector<lsdtt_xy>& in)
{
    std::vector<lsdtt_xy> out;
    out.reserve(in.size());
    for (const lsdtt_xy& p : in)
    {
        if (!std::isfinite(p.x) || !std::isfinite(p.y)) continue;
        if (!out.empty() && dist(out.back(), p) <= k_coincident) continue;
        out.push_back(p);
    }
    return out;
}

// Undirected bridges of a multigraph (iterative Tarjan). Self-loops are never
// bridges; parallel edges are never bridges.
std::vector<bool> find_bridges(int n_nodes, const std::vector<std::pair<int, int>>& edges)
{
    std::vector<std::vector<std::pair<int, int>>> adj(static_cast<size_t>(n_nodes));
    for (size_t e = 0; e < edges.size(); ++e)
    {
        const int a = edges[e].first;
        const int b = edges[e].second;
        if (a == b) continue;
        adj[static_cast<size_t>(a)].push_back({b, static_cast<int>(e)});
        adj[static_cast<size_t>(b)].push_back({a, static_cast<int>(e)});
    }

    std::vector<bool> is_bridge(edges.size(), false);
    std::vector<int> disc(static_cast<size_t>(n_nodes), -1);
    std::vector<int> low(static_cast<size_t>(n_nodes), -1);
    int timer = 0;

    struct frame { int node; int parent_edge; size_t next; };
    std::vector<frame> stack;

    for (int s = 0; s < n_nodes; ++s)
    {
        if (disc[static_cast<size_t>(s)] != -1) continue;
        disc[static_cast<size_t>(s)] = low[static_cast<size_t>(s)] = timer++;
        stack.push_back({s, -1, 0});
        while (!stack.empty())
        {
            frame& f = stack.back();
            const int v = f.node;
            const auto& nbrs = adj[static_cast<size_t>(v)];
            if (f.next < nbrs.size())
            {
                const int w = nbrs[f.next].first;
                const int e = nbrs[f.next].second;
                ++f.next;
                if (e == f.parent_edge) continue;
                if (disc[static_cast<size_t>(w)] == -1)
                {
                    disc[static_cast<size_t>(w)] = low[static_cast<size_t>(w)] = timer++;
                    stack.push_back({w, e, 0});
                }
                else
                {
                    low[static_cast<size_t>(v)] = std::min(low[static_cast<size_t>(v)],
                                                           disc[static_cast<size_t>(w)]);
                }
            }
            else
            {
                const int pe = f.parent_edge;
                stack.pop_back();
                if (!stack.empty())
                {
                    const int u = stack.back().node;
                    low[static_cast<size_t>(u)] = std::min(low[static_cast<size_t>(u)],
                                                           low[static_cast<size_t>(v)]);
                    if (low[static_cast<size_t>(v)] > disc[static_cast<size_t>(u)])
                    {
                        is_bridge[static_cast<size_t>(pe)] = true;
                    }
                }
            }
        }
    }
    return is_bridge;
}

} // namespace

double polyline_length(const std::vector<lsdtt_xy>& line)
{
    double s = 0.0;
    for (size_t i = 1; i < line.size(); ++i) s += dist(line[i - 1], line[i]);
    return s;
}

std::vector<std::vector<lsdtt_xy>> clip_polyline_to_box(const std::vector<lsdtt_xy>& line,
                                                        double x_min, double y_min,
                                                        double x_max, double y_max,
                                                        std::vector<bool>& cut_start,
                                                        std::vector<bool>& cut_end)
{
    std::vector<std::vector<lsdtt_xy>> pieces;
    cut_start.clear();
    cut_end.clear();
    std::vector<lsdtt_xy> current;
    bool current_cut_start = false;

    auto inside = [&](const lsdtt_xy& p)
    {
        return p.x >= x_min && p.x <= x_max && p.y >= y_min && p.y <= y_max;
    };
    auto close_piece = [&](bool was_cut)
    {
        if (current.size() >= 2 && polyline_length(current) > k_coincident)
        {
            pieces.push_back(current);
            cut_start.push_back(current_cut_start);
            cut_end.push_back(was_cut);
        }
        current.clear();
    };

    for (size_t i = 0; i + 1 < line.size(); ++i)
    {
        const lsdtt_xy a = line[i];
        const lsdtt_xy b = line[i + 1];
        // Liang-Barsky
        double t0 = 0.0;
        double t1 = 1.0;
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double p[4] = {-dx, dx, -dy, dy};
        const double q[4] = {a.x - x_min, x_max - a.x, a.y - y_min, y_max - a.y};
        bool reject = false;
        for (int k = 0; k < 4 && !reject; ++k)
        {
            if (p[k] == 0.0)
            {
                if (q[k] < 0.0) reject = true;
            }
            else
            {
                const double r = q[k] / p[k];
                if (p[k] < 0.0) { if (r > t1) reject = true; else if (r > t0) t0 = r; }
                else            { if (r < t0) reject = true; else if (r < t1) t1 = r; }
            }
        }
        if (reject)
        {
            if (!current.empty()) close_piece(true);
            continue;
        }
        const lsdtt_xy ca{a.x + t0 * dx, a.y + t0 * dy};
        const lsdtt_xy cb{a.x + t1 * dx, a.y + t1 * dy};
        if (current.empty())
        {
            current.push_back(ca);
            current_cut_start = (t0 > 0.0) || !inside(a);
        }
        else if (dist(current.back(), ca) > k_coincident)
        {
            current.push_back(ca);
        }
        if (dist(current.back(), cb) > k_coincident) current.push_back(cb);
        if (t1 < 1.0) close_piece(true);
    }
    if (!current.empty()) close_piece(false);
    return pieces;
}

river_line_network prepare_river_lines(const std::vector<river_line_feature>& features,
                                       const river_line_prepare_options& opts)
{
    river_line_network net;
    const double tol = std::max(opts.node_snap_tolerance_metres, 0.0);
    const double tol_eff = std::max(tol, k_coincident);

    // ---- 1-3. explode, clean, orient, clip ---------------------------------
    std::vector<work_part> parts;
    for (size_t fi = 0; fi < features.size(); ++fi)
    {
        const river_line_feature& f = features[fi];
        const bool multipart = f.parts.size() > 1;
        for (const auto& raw : f.parts)
        {
            ++net.n_input_parts;
            work_part wp;
            wp.v = drop_repeated_vertices(raw);
            wp.feature = fi;
            if (!multipart)
            {
                wp.start_id = f.start_node_id;
                wp.end_id = f.end_node_id;
            }
            if (wp.v.size() < 2 || polyline_length(wp.v) <= k_coincident)
            {
                ++net.n_parts_dropped_degenerate;
                continue;
            }
            if (f.direction < 0)
            {
                std::reverse(wp.v.begin(), wp.v.end());
                std::swap(wp.start_id, wp.end_id);
                ++net.n_parts_reversed;
            }
            if (!opts.clip_to_box)
            {
                parts.push_back(std::move(wp));
                continue;
            }
            std::vector<bool> cs, ce;
            auto pieces = clip_polyline_to_box(wp.v, opts.clip_x_min, opts.clip_y_min,
                                               opts.clip_x_max, opts.clip_y_max, cs, ce);
            if (pieces.empty())
            {
                ++net.n_parts_clipped_away;
                continue;
            }
            for (size_t k = 0; k < pieces.size(); ++k)
            {
                work_part piece;
                piece.v = std::move(pieces[k]);
                piece.feature = fi;
                piece.cut_start = cs[k];
                piece.cut_end = ce[k];
                piece.start_id = (!cs[k] && k == 0) ? wp.start_id : std::string();
                piece.end_id = (!ce[k] && k + 1 == pieces.size()) ? wp.end_id : std::string();
                parts.push_back(std::move(piece));
            }
        }
    }

    // ---- 4a. merge endpoints ------------------------------------------------
    const int n_end = static_cast<int>(parts.size() * 2);
    union_find uf(n_end);
    auto endpoint = [&](int e) -> lsdtt_xy&
    {
        work_part& p = parts[static_cast<size_t>(e / 2)];
        return (e % 2 == 0) ? p.v.front() : p.v.back();
    };
    {
        std::map<std::string, int> by_id;
        for (int e = 0; e < n_end; ++e)
        {
            const work_part& p = parts[static_cast<size_t>(e / 2)];
            const std::string& id = (e % 2 == 0) ? p.start_id : p.end_id;
            if (id.empty()) continue;
            auto it = by_id.find(id);
            if (it == by_id.end()) by_id[id] = e;
            else uf.unite(it->second, e);
        }
        point_grid g(std::max(tol_eff, 1e-3));
        for (int e = 0; e < n_end; ++e) g.insert(endpoint(e), e);
        for (int e = 0; e < n_end; ++e)
        {
            const lsdtt_xy p = endpoint(e);
            g.visit(p.x - tol_eff, p.y - tol_eff, p.x + tol_eff, p.y + tol_eff,
                    [&](int o)
                    {
                        if (o != e && dist(endpoint(o), p) <= tol_eff) uf.unite(e, o);
                    });
        }
    }

    // Nodes: one per endpoint cluster, at the members' mean.
    std::vector<int> node_of_root(static_cast<size_t>(n_end), -1);
    std::vector<lsdtt_xy> node_xy;
    std::vector<int> node_count;
    std::vector<bool> node_clip;
    for (int e = 0; e < n_end; ++e)
    {
        const int r = uf.find(e);
        int& nid = node_of_root[static_cast<size_t>(r)];
        if (nid < 0)
        {
            nid = static_cast<int>(node_xy.size());
            node_xy.push_back({0.0, 0.0});
            node_count.push_back(0);
            node_clip.push_back(false);
        }
        const lsdtt_xy p = endpoint(e);
        node_xy[static_cast<size_t>(nid)].x += p.x;
        node_xy[static_cast<size_t>(nid)].y += p.y;
        ++node_count[static_cast<size_t>(nid)];
        const work_part& wp = parts[static_cast<size_t>(e / 2)];
        if ((e % 2 == 0) ? wp.cut_start : wp.cut_end) node_clip[static_cast<size_t>(nid)] = true;
    }
    for (size_t n = 0; n < node_xy.size(); ++n)
    {
        node_xy[n].x /= node_count[n];
        node_xy[n].y /= node_count[n];
    }
    net.n_endpoints_merged = n_end - static_cast<int>(node_xy.size());
    for (int e = 0; e < n_end; ++e)
    {
        const int nid = node_of_root[static_cast<size_t>(uf.find(e))];
        endpoint(e) = node_xy[static_cast<size_t>(nid)];
        work_part& wp = parts[static_cast<size_t>(e / 2)];
        if (e % 2 == 0) wp.from_node = nid; else wp.to_node = nid;
    }

    // ---- 4b. split at T-junctions ------------------------------------------
    std::vector<work_part> split_parts;
    {
        point_grid g(std::max(tol_eff, 25.0));
        for (size_t n = 0; n < node_xy.size(); ++n) g.insert(node_xy[n], static_cast<int>(n));

        for (work_part& wp : parts)
        {
            // node -> (segment, t, distance), best segment per node
            std::map<int, std::tuple<size_t, double, double>> hits;
            for (size_t i = 0; i + 1 < wp.v.size(); ++i)
            {
                const lsdtt_xy a = wp.v[i];
                const lsdtt_xy b = wp.v[i + 1];
                const double dx = b.x - a.x;
                const double dy = b.y - a.y;
                const double len2 = dx * dx + dy * dy;
                if (len2 <= 0.0) continue;
                g.visit(std::min(a.x, b.x) - tol_eff, std::min(a.y, b.y) - tol_eff,
                        std::max(a.x, b.x) + tol_eff, std::max(a.y, b.y) + tol_eff,
                        [&](int n)
                        {
                            if (n == wp.from_node || n == wp.to_node) return;
                            const lsdtt_xy p = node_xy[static_cast<size_t>(n)];
                            double t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2;
                            t = std::min(1.0, std::max(0.0, t));
                            const lsdtt_xy q{a.x + t * dx, a.y + t * dy};
                            const double d = dist(p, q);
                            if (d > tol_eff) return;
                            // Not at the line's own ends: those are merged already.
                            if (dist(p, wp.v.front()) <= tol_eff || dist(p, wp.v.back()) <= tol_eff) return;
                            auto it = hits.find(n);
                            if (it == hits.end() || d < std::get<2>(it->second))
                            {
                                hits[n] = std::make_tuple(i, t, d);
                            }
                        });
            }
            if (hits.empty())
            {
                split_parts.push_back(std::move(wp));
                continue;
            }
            // A cut is (vertex or segment index, t, node). A node within the
            // tolerance of an interior vertex REPLACES that vertex (t = -1, sorts
            // first at that index), so no spike is left beside it; otherwise it
            // is inserted inside segment i at parameter t.
            std::vector<std::tuple<size_t, double, int>> cuts;
            for (const auto& h : hits)
            {
                const size_t i = std::get<0>(h.second);
                const double t = std::get<1>(h.second);
                const lsdtt_xy p = node_xy[static_cast<size_t>(h.first)];
                if (i + 2 < wp.v.size() && dist(p, wp.v[i + 1]) <= tol_eff)
                {
                    cuts.push_back({i + 1, -1.0, h.first});
                }
                else if (i > 0 && dist(p, wp.v[i]) <= tol_eff)
                {
                    cuts.push_back({i, -1.0, h.first});
                }
                else
                {
                    cuts.push_back({i, t, h.first});
                }
            }
            std::sort(cuts.begin(), cuts.end());

            work_part piece;
            piece.feature = wp.feature;
            piece.cut_start = wp.cut_start;
            piece.from_node = wp.from_node;
            size_t ci = 0;
            int n_made = 0;
            auto emit = [&](int to_node, bool cut_end)
            {
                piece.to_node = to_node;
                piece.cut_end = cut_end;
                if (piece.v.size() >= 2 && polyline_length(piece.v) > k_coincident)
                {
                    split_parts.push_back(piece);
                    ++n_made;
                }
                piece.v.clear();
                piece.cut_start = false;
            };
            auto push = [&](const lsdtt_xy& p)
            {
                if (piece.v.empty() || dist(piece.v.back(), p) > k_coincident) piece.v.push_back(p);
            };
            auto cut_at = [&](int n)
            {
                const lsdtt_xy p = node_xy[static_cast<size_t>(n)];
                push(p);
                emit(n, false);
                piece.from_node = n;
                piece.v.push_back(p);
            };
            for (size_t i = 0; i < wp.v.size(); ++i)
            {
                bool replaced = false;
                while (ci < cuts.size() && std::get<0>(cuts[ci]) == i && std::get<1>(cuts[ci]) < 0.0)
                {
                    cut_at(std::get<2>(cuts[ci]));
                    replaced = true;
                    ++ci;
                }
                if (!replaced) push(wp.v[i]);
                while (ci < cuts.size() && std::get<0>(cuts[ci]) == i)
                {
                    cut_at(std::get<2>(cuts[ci]));
                    ++ci;
                }
            }
            emit(wp.to_node, wp.cut_end);
            const int n_split = n_made - 1;
            for (size_t k = split_parts.size() - static_cast<size_t>(n_made); k < split_parts.size(); ++k)
            {
                split_parts[k].n_splits = n_split;
            }
            net.n_t_junction_splits += n_split;
        }
    }

    // ---- 5. drop isolated short fragments ----------------------------------
    const int n_nodes_all = static_cast<int>(node_xy.size());
    std::vector<bool> keep(split_parts.size(), true);
    if (opts.min_line_length_metres > 0.0)
    {
        union_find cu(n_nodes_all);
        for (const work_part& p : split_parts) cu.unite(p.from_node, p.to_node);
        std::map<int, int> lines_in;
        for (const work_part& p : split_parts) ++lines_in[cu.find(p.from_node)];
        for (size_t k = 0; k < split_parts.size(); ++k)
        {
            const work_part& p = split_parts[k];
            if (lines_in[cu.find(p.from_node)] == 1 && polyline_length(p.v) < opts.min_line_length_metres)
            {
                keep[k] = false;
                ++net.n_fragments_dropped;
            }
        }
    }

    // Renumber the nodes that are still used.
    std::vector<int> new_id(static_cast<size_t>(n_nodes_all), -1);
    for (size_t k = 0; k < split_parts.size(); ++k)
    {
        if (!keep[k]) continue;
        for (int n : {split_parts[k].from_node, split_parts[k].to_node})
        {
            if (new_id[static_cast<size_t>(n)] < 0)
            {
                new_id[static_cast<size_t>(n)] = static_cast<int>(net.nodes.size());
                river_node rn;
                rn.node_id = static_cast<int>(net.nodes.size());
                rn.xy = node_xy[static_cast<size_t>(n)];
                rn.is_clip_point = node_clip[static_cast<size_t>(n)];
                net.nodes.push_back(rn);
            }
        }
    }
    for (size_t k = 0; k < split_parts.size(); ++k)
    {
        if (!keep[k]) continue;
        const work_part& p = split_parts[k];
        const river_line_feature& f = features[p.feature];
        river_line rl;
        rl.line_id = static_cast<int>(net.lines.size());
        rl.vertices = p.v;
        rl.from_node = new_id[static_cast<size_t>(p.from_node)];
        rl.to_node = new_id[static_cast<size_t>(p.to_node)];
        rl.source_feature_id = f.source_feature_id;
        rl.name = f.name;
        rl.feature_class = f.feature_class;
        rl.direction_source = f.direction_source;
        rl.culverted = f.culverted;
        rl.strahler_order = f.strahler_order;
        rl.n_splits = p.n_splits;
        rl.length_m = polyline_length(p.v);
        net.lines.push_back(std::move(rl));
    }

    // ---- 6. direction consistency from the graph alone ----------------------
    // In each component the candidate outlets are every degree-1 node and
    // every data sink (inflow, no outflow: e.g. a confluence whose outflow was
    // cut off by an extract). For a
    // candidate o, a line u -> v agrees with "toward o" when v is nearer o
    // than u (network distance, Dijkstra on line length), and disagrees when
    // it is farther; equal distances (the two halves of a symmetric braid)
    // are not judged. The outlet is the candidate with the least LENGTH of
    // line disagreeing (then the fewest lines, then a node the data already
    // make a sink, then the lowest id), so the data's own directions vote for
    // it and a long main stem outweighs a short spur.
    std::vector<std::pair<int, int>> chosen_outlets;   // (node, n_against)
    {
        const int nn = static_cast<int>(net.nodes.size());
        std::vector<std::vector<std::pair<int, int>>> adj(static_cast<size_t>(nn));  // (other, line)
        std::vector<int> din(static_cast<size_t>(nn), 0), dout(static_cast<size_t>(nn), 0);
        for (size_t k = 0; k < net.lines.size(); ++k)
        {
            const river_line& l = net.lines[k];
            adj[static_cast<size_t>(l.from_node)].push_back({l.to_node, static_cast<int>(k)});
            adj[static_cast<size_t>(l.to_node)].push_back({l.from_node, static_cast<int>(k)});
            ++dout[static_cast<size_t>(l.from_node)];
            ++din[static_cast<size_t>(l.to_node)];
        }
        union_find cu(nn);
        for (const river_line& l : net.lines) cu.unite(l.from_node, l.to_node);
        std::map<int, std::vector<int>> comp_nodes, comp_lines;
        for (int n = 0; n < nn; ++n) comp_nodes[cu.find(n)].push_back(n);
        for (size_t k = 0; k < net.lines.size(); ++k)
            comp_lines[cu.find(net.lines[k].from_node)].push_back(static_cast<int>(k));

        const double inf = std::numeric_limits<double>::infinity();
        std::vector<double> dist(static_cast<size_t>(nn), inf);
        auto dijkstra = [&](int src, const std::vector<int>& nodes)
        {
            for (int n : nodes) dist[static_cast<size_t>(n)] = inf;
            using item = std::pair<double, int>;
            std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
            dist[static_cast<size_t>(src)] = 0.0;
            pq.push({0.0, src});
            while (!pq.empty())
            {
                const auto [d, u] = pq.top();
                pq.pop();
                if (d > dist[static_cast<size_t>(u)]) continue;
                for (const auto& [v, k] : adj[static_cast<size_t>(u)])
                {
                    const double nd = d + std::max(net.lines[static_cast<size_t>(k)].length_m, 1e-9);
                    if (nd < dist[static_cast<size_t>(v)])
                    {
                        dist[static_cast<size_t>(v)] = nd;
                        pq.push({nd, v});
                    }
                }
            }
        };
        // +1 agrees, -1 disagrees, 0 not judged, for line k toward the current dist.
        auto judge = [&](int k) -> int
        {
            const river_line& l = net.lines[static_cast<size_t>(k)];
            const double du = dist[static_cast<size_t>(l.from_node)];
            const double dv = dist[static_cast<size_t>(l.to_node)];
            if (std::fabs(du - dv) <= 1e-9 * std::max(1.0, std::max(du, dv))) return 0;
            return dv < du ? 1 : -1;
        };

        for (const auto& [root, nodes] : comp_nodes)
        {
            const std::vector<int>& lines = comp_lines[root];
            if (lines.empty()) continue;
            std::vector<int> cands;
            for (int n : nodes)
            {
                const int i = din[static_cast<size_t>(n)], o = dout[static_cast<size_t>(n)];
                if (i + o == 1 || (o == 0 && i >= 1)) cands.push_back(n);
            }
            if (cands.empty()) continue;   // a ring: no outlet to vote for
            if (static_cast<long long>(cands.size()) * static_cast<long long>(lines.size()) >
                opts.max_direction_check_work)
            {
                net.notes.push_back("Direction check skipped for a component of " +
                                    std::to_string(lines.size()) + " lines (too large)");
                continue;
            }
            int best = -1, best_against = 0;
            double best_len = 0.0;
            bool best_is_outlet = false;
            for (int o : cands)
            {
                dijkstra(o, nodes);
                int against = 0;
                double len = 0.0;
                for (int k : lines)
                {
                    if (judge(k) < 0)
                    {
                        ++against;
                        len += net.lines[static_cast<size_t>(k)].length_m;
                    }
                }
                const bool is_outlet = din[static_cast<size_t>(o)] >= 1 && dout[static_cast<size_t>(o)] == 0;
                const double eps = 1e-6;
                bool better = best < 0 || len < best_len - eps;
                if (!better && std::fabs(len - best_len) <= eps)
                {
                    better = against < best_against ||
                             (against == best_against && is_outlet && !best_is_outlet);
                }
                if (better)
                {
                    best = o; best_against = against; best_len = len; best_is_outlet = is_outlet;
                }
            }
            dijkstra(best, nodes);
            for (int k : lines)
            {
                river_line& l = net.lines[static_cast<size_t>(k)];
                const int j = judge(k);
                l.direction_consistent = j > 0 ? 1 : (j < 0 ? 0 : -1);
                if (j < 0)
                {
                    ++net.n_lines_against_network;
                    if (opts.repair_direction)
                    {
                        std::reverse(l.vertices.begin(), l.vertices.end());
                        std::swap(l.from_node, l.to_node);
                        l.direction_source = "network";
                        ++net.n_lines_flipped;
                    }
                }
            }
            chosen_outlets.push_back({best, best_against});
        }
    }

    // ---- 7. topology report -------------------------------------------------
    const int n_nodes = static_cast<int>(net.nodes.size());
    std::vector<std::pair<int, int>> edges;
    edges.reserve(net.lines.size());
    for (const river_line& l : net.lines) edges.push_back({l.from_node, l.to_node});

    const std::vector<bool> bridge = find_bridges(n_nodes, edges);
    for (size_t e = 0; e < net.lines.size(); ++e)
    {
        net.lines[e].on_loop = !bridge[e];
        ++net.nodes[static_cast<size_t>(edges[e].first)].degree_out;
        ++net.nodes[static_cast<size_t>(edges[e].second)].degree_in;
    }

    union_find cu(n_nodes);
    for (const auto& e : edges) cu.unite(e.first, e.second);
    std::map<int, int> comp_of_root;
    for (river_node& n : net.nodes)
    {
        const int r = cu.find(n.node_id);
        auto it = comp_of_root.find(r);
        if (it == comp_of_root.end())
        {
            const int c = static_cast<int>(net.components.size());
            comp_of_root[r] = c;
            river_component_summary s;
            s.component_id = c;
            net.components.push_back(s);
            n.component_id = c;
        }
        else
        {
            n.component_id = it->second;
        }
    }

    std::vector<bool> voted_outlet(net.nodes.size(), false);
    for (const auto& co : chosen_outlets) voted_outlet[static_cast<size_t>(co.first)] = true;
    for (river_node& n : net.nodes)
    {
        const int din = n.degree_in;
        const int dout = n.degree_out;
        if (din + dout == 1 && n.is_clip_point) n.node_type = "clip_edge";
        // A sink chosen as its network's outlet by the direction vote is that
        // network's outlet, even with several inflows (a confluence whose
        // outflow lies outside the data).
        else if (voted_outlet[static_cast<size_t>(n.node_id)] && dout == 0) n.node_type = "outlet";
        else if (din == 0 && dout == 1)          n.node_type = "tip";
        else if (din == 1 && dout == 0)          n.node_type = "outlet";
        else if (din >= 1 && dout >= 2)          n.node_type = "bifurcation";
        else if (din >= 2 && dout == 1)          n.node_type = "confluence";
        else if (din == 1 && dout == 1)          n.node_type = "link";
        else                                     n.node_type = "conflict";

        river_component_summary& s = net.components[static_cast<size_t>(n.component_id)];
        if (n.node_type == "tip") ++s.n_tips;
        else if (n.node_type == "outlet") ++s.n_outlets;
        else if (n.node_type == "clip_edge") ++s.n_clip_edges;
        else if (n.node_type == "bifurcation") ++s.n_bifurcations;
        else if (n.node_type == "conflict") ++s.n_conflicts;
    }

    std::vector<int> comp_nodes(net.components.size(), 0);
    for (const river_node& n : net.nodes) ++comp_nodes[static_cast<size_t>(n.component_id)];
    for (river_line& l : net.lines)
    {
        l.component_id = net.nodes[static_cast<size_t>(l.from_node)].component_id;
        river_component_summary& s = net.components[static_cast<size_t>(l.component_id)];
        ++s.n_lines;
        s.length_m += l.length_m;
    }
    for (river_component_summary& s : net.components)
    {
        s.n_loops = s.n_lines - comp_nodes[static_cast<size_t>(s.component_id)] + 1;
    }

    for (const auto& [o, against] : chosen_outlets)
    {
        river_component_summary& c =
            net.components[static_cast<size_t>(net.nodes[static_cast<size_t>(o)].component_id)];
        c.outlet_node = o;
        c.n_against_network = against;
    }

    // ---- notes ---------------------------------------------------------------
    {
        int loops = 0, tips = 0, outlets = 0, clips = 0, bifs = 0, conflicts = 0;
        for (const auto& s : net.components)
        {
            loops += s.n_loops; tips += s.n_tips; outlets += s.n_outlets;
            clips += s.n_clip_edges; bifs += s.n_bifurcations; conflicts += s.n_conflicts;
        }
        std::ostringstream os;
        os << "Parts read: " << net.n_input_parts << " (degenerate dropped "
           << net.n_parts_dropped_degenerate << ", reversed " << net.n_parts_reversed
           << ", outside the box " << net.n_parts_clipped_away << ")";
        net.notes.push_back(os.str());
        os.str("");
        os << "Noding: " << net.n_endpoints_merged << " endpoints merged, "
           << net.n_t_junction_splits << " T-junction splits, "
           << net.n_fragments_dropped << " isolated fragments dropped";
        net.notes.push_back(os.str());
        os.str("");
        os << "Direction (graph only): " << net.n_lines_against_network
           << " line(s) pointed against their network's outlet"
           << (opts.repair_direction ? ", flipped" : ", not flipped (repair off)");
        net.notes.push_back(os.str());
        os.str("");
        os << "Network: " << net.lines.size() << " lines, " << net.nodes.size() << " nodes, "
           << net.components.size() << " components; tips " << tips << ", outlets " << outlets
           << ", clip-edge ends " << clips << ", bifurcations " << bifs << ", loops " << loops
           << ", direction conflicts " << conflicts;
        net.notes.push_back(os.str());
    }
    return net;
}

} // namespace lsdtt
