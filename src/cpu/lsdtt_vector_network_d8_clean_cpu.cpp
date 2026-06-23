// lsdtt_vector_network_d8_clean_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_vector_network_d8_clean_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <set>
#include <unordered_set>

namespace lsdtt
{

namespace
{

// World coordinate -> grid column / row (floor). Returns false if outside grid.
inline bool xy_to_rc(const lsdtt_raster_meta& m, double x, double y,
                     int& row, int& col)
{
    col = static_cast<int>(std::floor((x - m.x_min) / m.cell_size));
    row = static_cast<int>(std::floor((m.y_max() - y) / m.cell_size));
    return row >= 0 && row < m.n_rows && col >= 0 && col < m.n_cols;
}

inline lsdtt_xy rc_to_centre(const lsdtt_raster_meta& m, int row, int col)
{
    lsdtt_xy p;
    p.x = m.x_min + (static_cast<double>(col) + 0.5) * m.cell_size;
    p.y = m.y_max() - (static_cast<double>(row) + 0.5) * m.cell_size;
    return p;
}

// Densely sample an edge polyline at `step` metres, returning sample points.
std::vector<lsdtt_xy> sample_polyline(const std::vector<lsdtt_xy>& xy,
                                      double step)
{
    std::vector<lsdtt_xy> out;
    if (xy.empty()) return out;
    out.push_back(xy.front());
    for (std::size_t k = 1; k < xy.size(); ++k)
    {
        const lsdtt_xy& a = xy[k - 1];
        const lsdtt_xy& b = xy[k];
        const double seg = std::hypot(b.x - a.x, b.y - a.y);
        const int n = std::max(1, static_cast<int>(seg / step));
        for (int s = 1; s <= n; ++s)
        {
            const double t = static_cast<double>(s) / n;
            out.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)});
        }
    }
    return out;
}

} // namespace

// ---------------------------------------------------------------------------
// clean_vector_network_with_d8
// ---------------------------------------------------------------------------
vector_network_d8_clean_result clean_vector_network_with_d8(
    const lsdtt_vector_to_junction_network& network,
    const lsdtt_flow_base&                  flow,
    const vector_network_d8_clean_options&  options)
{
    vector_network_d8_clean_result result;

    const lsdtt_raster_meta& meta = flow.meta();
    const int    n_rows = meta.n_rows;
    const int    n_cols = meta.n_cols;
    const double cs     = meta.cell_size;

    const double step = (options.edge_sample_step_metres > 0.0)
                            ? options.edge_sample_step_metres
                            : 0.5 * cs;
    const double search_max_m = (options.channel_search_max_metres > 0.0)
                                    ? options.channel_search_max_metres
                                    : 250.0;
    const int search_max_cells =
        std::max(1, static_cast<int>(std::ceil(search_max_m / cs)));

    const auto& junctions = network.junctions();
    const auto& edges     = network.edges();
    const std::size_t n_edges = edges.size();

    // -----------------------------------------------------------------------
    // D8 channel mask from the contributing-pixel threshold.
    // -----------------------------------------------------------------------
    std::vector<std::uint8_t> channel_mask;
    {
        lsdtt_raster_result ncp = flow.get_n_contributing_pixels_raster();
        channel_mask.assign(static_cast<std::size_t>(n_rows) * n_cols, 0);
        const std::size_t ncells = channel_mask.size();
        for (std::size_t i = 0; i < ncells && i < ncp.data.size(); ++i)
        {
            const bool valid = (i < ncp.mask.size()) ? (ncp.mask[i] != 0) : true;
            if (valid &&
                ncp.data[i] >= static_cast<float>(options.channel_threshold_pixels))
            {
                channel_mask[i] = 1;
                ++result.channel_mask_cells;
            }
        }
    }

    // Bounded nearest-channel distance (metres) from a grid cell. Returns a
    // large sentinel if no channel cell lies within search_max_cells.
    const double kFar = static_cast<double>(search_max_cells + 1) * cs;
    auto channel_dist = [&](int row, int col) -> double
    {
        double best = std::numeric_limits<double>::infinity();
        for (int r = 0; r <= search_max_cells; ++r)
        {
            if (static_cast<double>(r) * cs > best) break;
            // Iterate the Chebyshev ring of radius r.
            const int r0 = row - r, r1 = row + r;
            const int c0 = col - r, c1 = col + r;
            auto test = [&](int rr, int cc)
            {
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) return;
                if (!channel_mask[static_cast<std::size_t>(rr) * n_cols + cc])
                    return;
                const double d = std::hypot(static_cast<double>(rr - row),
                                            static_cast<double>(cc - col)) * cs;
                if (d < best) best = d;
            };
            if (r == 0) { test(row, col); continue; }
            for (int cc = c0; cc <= c1; ++cc) { test(r0, cc); test(r1, cc); }
            for (int rr = r0 + 1; rr <= r1 - 1; ++rr) { test(rr, c0); test(rr, c1); }
        }
        return std::isfinite(best) ? best : kFar;
    };

    // Mean channel distance of an edge (sampled).
    auto edge_channel_dist = [&](const lsdtt_vector_network_edge& e) -> double
    {
        const std::vector<lsdtt_xy> samples = sample_polyline(e.xy, step);
        if (samples.empty()) return kFar;
        double sum = 0.0; int n = 0;
        for (const lsdtt_xy& p : samples)
        {
            int row, col;
            if (!xy_to_rc(meta, p.x, p.y, row, col)) { sum += kFar; ++n; continue; }
            sum += channel_dist(row, col); ++n;
        }
        return n > 0 ? sum / n : kFar;
    };

    // -----------------------------------------------------------------------
    // Undirected adjacency over all compact edges.
    // -----------------------------------------------------------------------
    const std::size_t n_j = junctions.size();
    std::vector<std::vector<std::pair<std::int32_t, std::int32_t>>> adj(n_j);
    for (std::size_t e = 0; e < n_edges; ++e)
    {
        const std::int32_t a = edges[e].junction_a;
        const std::int32_t b = edges[e].junction_b;
        if (a < 0 || b < 0) continue;
        adj[static_cast<std::size_t>(a)].push_back({b, static_cast<std::int32_t>(e)});
        adj[static_cast<std::size_t>(b)].push_back({a, static_cast<std::int32_t>(e)});
    }

    std::vector<std::uint8_t> removed(n_edges, 0);   // 1 = removed (loop arm)

    // -----------------------------------------------------------------------
    // Layer 2b — loop-arm resolution.
    // -----------------------------------------------------------------------
    if (options.resolve_loops)
    {
        for (std::size_t e = 0; e < n_edges; ++e)
        {
            if (!edges[e].is_loop_edge) continue;
            if (removed[e]) continue;
            ++result.n_loops_detected;

            const std::int32_t a = edges[e].junction_a;
            const std::int32_t b = edges[e].junction_b;
            if (a < 0 || b < 0) continue;

            // BFS from a to b over kept, non-loop edges to find the tree path.
            std::vector<std::int32_t> prev_j(n_j, -1);
            std::vector<std::int32_t> prev_e(n_j, -1);
            std::vector<std::uint8_t> seen(n_j, 0);
            std::deque<std::int32_t> q;
            q.push_back(a); seen[static_cast<std::size_t>(a)] = 1;
            bool found = false;
            while (!q.empty() && !found)
            {
                const std::int32_t u = q.front(); q.pop_front();
                for (const auto& [v, ei] : adj[static_cast<std::size_t>(u)])
                {
                    if (static_cast<std::size_t>(ei) == e) continue;       // skip the loop edge itself
                    if (removed[static_cast<std::size_t>(ei)]) continue;
                    if (edges[static_cast<std::size_t>(ei)].is_loop_edge) continue;
                    if (seen[static_cast<std::size_t>(v)]) continue;
                    seen[static_cast<std::size_t>(v)] = 1;
                    prev_j[static_cast<std::size_t>(v)] = u;
                    prev_e[static_cast<std::size_t>(v)] = ei;
                    if (v == b) { found = true; break; }
                    q.push_back(v);
                }
            }

            if (!found)
            {
                // No alternative path: this labelled edge is actually a bridge,
                // not a loop. Keep it as a normal channel.
                continue;
            }

            // Reconstruct the junction ring J[0..k] (a..b) and edge ring.
            std::vector<std::int32_t> ring_j;     // junctions a .. b
            std::vector<std::int32_t> ring_e;     // path edges, then loop edge e
            for (std::int32_t v = b; v != a; v = prev_j[static_cast<std::size_t>(v)])
            {
                ring_j.push_back(v);
                ring_e.push_back(prev_e[static_cast<std::size_t>(v)]);
            }
            ring_j.push_back(a);
            std::reverse(ring_j.begin(), ring_j.end());  // J[0]=a .. J[k]=b
            std::reverse(ring_e.begin(), ring_e.end());  // edges J[i]-J[i+1]
            ring_e.push_back(static_cast<std::int32_t>(e)); // closing edge b-a

            const int k = static_cast<int>(ring_j.size());  // junction count
            std::unordered_set<std::int32_t> cycle_edge_set(ring_e.begin(),
                                                            ring_e.end());

            // Hubs: ring junctions with an incident edge outside the cycle.
            std::vector<int> hub_pos;
            for (int i = 0; i < k; ++i)
            {
                const std::int32_t j = ring_j[static_cast<std::size_t>(i)];
                bool has_external = false;
                for (const auto& [v, ei] : adj[static_cast<std::size_t>(j)])
                {
                    (void)v;
                    if (removed[static_cast<std::size_t>(ei)]) continue;
                    if (!cycle_edge_set.count(ei)) { has_external = true; break; }
                }
                if (has_external) hub_pos.push_back(i);
            }

            std::vector<std::int32_t> arm_to_remove;

            if (hub_pos.size() == 2)
            {
                // Two arcs of the edge ring between the two hub positions.
                // Edge ring index i connects ring_j[i] and ring_j[(i+1)%k].
                const int p1 = hub_pos[0];
                const int p2 = hub_pos[1];
                std::vector<std::int32_t> arc_fwd, arc_bwd;
                for (int i = p1; i != p2; i = (i + 1) % k)
                    arc_fwd.push_back(ring_e[static_cast<std::size_t>(i)]);
                for (int i = p2; i != p1; i = (i + 1) % k)
                    arc_bwd.push_back(ring_e[static_cast<std::size_t>(i)]);

                auto arc_score = [&](const std::vector<std::int32_t>& arc)
                {
                    double s = 0.0; int n = 0;
                    for (std::int32_t ei : arc)
                    {
                        s += edge_channel_dist(edges[static_cast<std::size_t>(ei)]);
                        ++n;
                    }
                    return n > 0 ? s / n : kFar;
                };
                const double sf = arc_score(arc_fwd);
                const double sb = arc_score(arc_bwd);
                arm_to_remove = (sf >= sb) ? arc_fwd : arc_bwd;  // farther arm
            }
            else
            {
                // Isolated ring or complex hub set: remove the single farthest
                // edge to break the cycle.
                double worst = -1.0; std::int32_t worst_e = -1;
                for (std::int32_t ei : ring_e)
                {
                    const double d =
                        edge_channel_dist(edges[static_cast<std::size_t>(ei)]);
                    if (d > worst) { worst = d; worst_e = ei; }
                }
                if (worst_e >= 0) arm_to_remove.push_back(worst_e);
            }

            int removed_here = 0;
            for (std::int32_t ei : arm_to_remove)
            {
                if (!removed[static_cast<std::size_t>(ei)])
                {
                    removed[static_cast<std::size_t>(ei)] = 1;
                    ++removed_here;
                }
            }
            if (removed_here > 0)
            {
                ++result.n_loops_resolved;
                result.n_loop_edges_removed += removed_here;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Segment-coverage raster of kept edges (cell -> edge index, else -1).
    // Built after loop removal so dead-end walks ignore removed arms.
    // -----------------------------------------------------------------------
    std::vector<std::int32_t> edge_cell(
        static_cast<std::size_t>(n_rows) * n_cols, -1);
    for (std::size_t e = 0; e < n_edges; ++e)
    {
        if (removed[e]) continue;
        for (const lsdtt_xy& p : sample_polyline(edges[e].xy, step))
        {
            int row, col;
            if (!xy_to_rc(meta, p.x, p.y, row, col)) continue;
            edge_cell[static_cast<std::size_t>(row) * n_cols + col] =
                static_cast<std::int32_t>(e);
        }
    }

    // -----------------------------------------------------------------------
    // Layer 2a — dead-end resolution.
    // -----------------------------------------------------------------------
    std::vector<vector_clean_segment> connectors;
    if (options.resolve_dead_ends)
    {
        const int max_steps =
            std::max(1, static_cast<int>(std::ceil(
                options.dead_end_search_max_metres / cs)) * 2);

        for (std::size_t j = 0; j < n_j; ++j)
        {
            const lsdtt_vector_network_junction& J = junctions[j];
            // Interior downstream tip: drains to the outlet sentinel but is not
            // a genuine DEM-edge outlet.
            if (!J.is_outlet) continue;
            if (J.is_dem_edge_junction) continue;

            // Incident (kept) edges of this junction — used to ignore self-hits.
            std::unordered_set<std::int32_t> own_edges;
            for (const auto& [v, ei] : adj[j])
            {
                (void)v;
                if (!removed[static_cast<std::size_t>(ei)]) own_edges.insert(ei);
            }
            if (own_edges.empty()) continue;   // fully removed -> not a real tip
            ++result.n_dead_ends_detected;

            int row, col;
            if (!xy_to_rc(meta, J.x, J.y, row, col)) continue;

            std::int64_t cur = flow.flat_index_unchecked(row, col);
            std::vector<lsdtt_xy> path;
            path.push_back({J.x, J.y});

            double walked = 0.0;
            bool joined = false;
            for (int s = 0; s < max_steps; ++s)
            {
                if (flow.is_outlet_unchecked(cur)) break;
                const std::int64_t nxt = flow.receiver_of_unchecked(cur);
                if (nxt == cur || nxt < 0) break;

                const int pr = flow.row_from_flat_unchecked(cur);
                const int pc = flow.col_from_flat_unchecked(cur);
                const int nr = flow.row_from_flat_unchecked(nxt);
                const int nc = flow.col_from_flat_unchecked(nxt);
                walked += std::hypot(static_cast<double>(nr - pr),
                                     static_cast<double>(nc - pc)) * cs;
                if (walked > options.dead_end_search_max_metres) break;

                path.push_back(rc_to_centre(meta, nr, nc));

                const std::int32_t hit =
                    edge_cell[static_cast<std::size_t>(nr) * n_cols + nc];
                if (hit >= 0 && !own_edges.count(hit))
                {
                    joined = true;
                    break;
                }
                cur = nxt;
            }

            if (joined && path.size() >= 2)
            {
                vector_clean_segment seg;
                seg.xy = std::move(path);
                seg.role = "dead_end_connector";
                connectors.push_back(std::move(seg));
                ++result.n_dead_ends_resolved;
            }
            else
            {
                ++result.n_dead_ends_unresolved;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Assemble output: every edge (tagged), then dead-end connectors.
    // -----------------------------------------------------------------------
    result.segments.reserve(n_edges + connectors.size());
    for (std::size_t e = 0; e < n_edges; ++e)
    {
        vector_clean_segment seg;
        seg.xy      = edges[e].xy;
        seg.edge_id = edges[e].edge_id;
        seg.role    = removed[e] ? "removed_loop_arm" : "channel";
        if (edges[e].is_loop_edge || removed[e])
            seg.mean_channel_dist_m = edge_channel_dist(edges[e]);
        result.segments.push_back(std::move(seg));
    }
    for (auto& c : connectors) result.segments.push_back(std::move(c));

    return result;
}

// ---------------------------------------------------------------------------
// write_cleaned_network
// ---------------------------------------------------------------------------
void write_cleaned_network(
    const vector_network_d8_clean_result& result,
    const lsdtt_raster_meta&              meta,
    const std::string&                    out_path,
    vector_output_format                  format,
    bool                                  spatial_index,
    const std::string&                    layer_name)
{
    lsdtt_vector_schema schema;
    schema.field_names = {"edge_id", "role", "mean_channel_dist_m"};
    schema.field_types = {lsdtt_vector_field_type::int64_type,
                          lsdtt_vector_field_type::string_type,
                          lsdtt_vector_field_type::double_type};

    lsdtt_data_loader writer;
    writer.begin_vector_write(out_path, layer_name, schema,
                              meta.crs_wkt, format, spatial_index);

    for (const vector_clean_segment& s : result.segments)
    {
        if (s.xy.size() < 2) continue;
        std::vector<lsdtt_vector_field_value> vals;
        vals.emplace_back(static_cast<std::int64_t>(s.edge_id));
        vals.emplace_back(s.role);
        vals.emplace_back(s.mean_channel_dist_m);
        writer.write_line(s.xy, vals);
    }

    writer.finish_vector_write();
}

} // namespace lsdtt
