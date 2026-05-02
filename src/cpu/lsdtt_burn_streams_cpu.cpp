// lsdtt_burn_streams_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_burn_streams_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "lsdtt_least_cost_path.hpp"
#include "lsdtt_proximity_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// find_vector_network_source_tips
// ---------------------------------------------------------------------------
vector_network_endpoints find_vector_network_source_tips(
    const std::string&      vector_path,
    const lsdtt_raster_cpu& dem,
    const lsdtt_flow_base*  flow)
{
    // Open the vector layer
    lsdtt_data_loader vec_loader;
    try
    {
        vec_loader.open_vector(vector_path);
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error(
            std::string("find_vector_network_source_tips: cannot open vector: ")
            + e.what());
    }

    const auto layer_names = vec_loader.get_layer_names();
    if (layer_names.empty())
    {
        throw std::runtime_error(
            "find_vector_network_source_tips: no layers in " + vector_path);
    }

    const lsdtt_vector_layer vlayer =
        vec_loader.read_vector_layer(layer_names[0]);
    vec_loader.close_vector();

    if (vlayer.empty())
    {
        throw std::runtime_error(
            "find_vector_network_source_tips: empty vector layer in "
            + vector_path);
    }

    const lsdtt_raster_meta& meta = dem.get_meta();

    // Endpoints are compared with a tolerance of half a cell size.
    // Coordinates are quantized to this resolution so that nearby endpoints
    // from the same junction match as the same key.
    const double snap_tol = meta.cell_size * 0.5;

    auto quantize = [snap_tol](double v) -> std::int64_t
    {
        return static_cast<std::int64_t>(std::round(v / snap_tol));
    };

    auto make_key = [&](double x, double y) -> std::string
    {
        return std::to_string(quantize(x)) + "_" + std::to_string(quantize(y));
    };

    // endpoint_map: key → number of features that have this as an endpoint
    std::map<std::string, int>                    endpoint_count;
    std::map<std::string, lsdtt_xy>               endpoint_xy;
    // For each endpoint store the full inward vertex sequence — all vertices
    // of its segment ordered from the tip inward (tip excluded).
    // When sampling the "adjacent" elevation we walk this sequence and use the
    // first vertex whose DEM cell is valid.  This handles the case where the
    // immediately adjacent vertex (or several vertices) lie outside the DEM or
    // on a nodata pixel, which happens whenever a segment crosses the raster
    // boundary near the tip.
    std::map<std::string, std::vector<lsdtt_xy>> endpoint_inward;

    for (std::size_t fi = 0; fi < vlayer.size(); ++fi)
    {
        const auto& record = vlayer.records[fi];
        const lsdtt_line_geom* line =
            std::get_if<lsdtt_line_geom>(&record.geometry);
        if (!line || line->points.size() < 2) continue;

        const std::size_t n = line->points.size();

        // front endpoint → inward sequence is points[1 .. n-1]
        {
            const lsdtt_xy&   ep  = line->points[0];
            const std::string key = make_key(ep.x, ep.y);
            endpoint_count[key]++;
            endpoint_xy.emplace(key, ep);
            if (endpoint_inward.find(key) == endpoint_inward.end())
            {
                std::vector<lsdtt_xy> seq(line->points.begin() + 1,
                                          line->points.end());
                endpoint_inward.emplace(key, std::move(seq));
            }
        }

        // back endpoint → inward sequence is points[n-2 .. 0] (reversed)
        {
            const lsdtt_xy&   ep  = line->points[n - 1];
            const std::string key = make_key(ep.x, ep.y);
            endpoint_count[key]++;
            endpoint_xy.emplace(key, ep);
            if (endpoint_inward.find(key) == endpoint_inward.end())
            {
                std::vector<lsdtt_xy> seq(line->points.rbegin() + 1,
                                          line->points.rend());
                endpoint_inward.emplace(key, std::move(seq));
            }
        }
    }

    // Collect degree-1 endpoints (appear as endpoint of exactly 1 feature).
    struct tip_info
    {
        lsdtt_xy      xy;
        double        elevation     = std::numeric_limits<double>::quiet_NaN();
        double        adjacent_elev = std::numeric_limits<double>::quiet_NaN();
        std::uint64_t upslope_count = 0;
    };
    std::vector<tip_info> tips;

    const lsdtt_array2d_f&  dem_data = dem.get_data_array();
    const lsdtt_array2d_u8& dem_mask = dem.get_mask_array();

    // Sample DEM elevation at (x, y); returns NaN when outside or on nodata.
    auto sample_elev = [&](double x, double y) -> double
    {
        const int col = static_cast<int>(
            std::floor((x - meta.x_min) / meta.cell_size));
        const int row = static_cast<int>(
            std::floor((meta.y_max() - y) / meta.cell_size));
        if (row < 0 || row >= meta.n_rows ||
            col < 0 || col >= meta.n_cols ||
            dem_mask(row, col) == 0)
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return static_cast<double>(dem_data(row, col));
    };

    for (const auto& [key, count] : endpoint_count)
    {
        if (count != 1) continue;

        const lsdtt_xy& pt = endpoint_xy.at(key);

        tip_info ti;
        ti.xy        = pt;
        ti.elevation = sample_elev(pt.x, pt.y);

        // Walk all inward vertices and record the minimum valid DEM elevation.
        //
        // Using the minimum (rather than the first valid vertex) makes the
        // gradient comparison robust to two common failure modes:
        //   1. The immediately adjacent vertex lies on a slightly higher DEM
        //      cell due to cell-alignment noise in flat terrain, causing a
        //      spuriously positive gradient at a source tip.
        //   2. The segment crosses the DEM boundary with several nodata
        //      vertices before reaching valid terrain.
        //
        // For a source tip, the segment runs downhill to the DEM exit, so the
        // minimum inward elevation (~exit vicinity) is well below the tip even
        // in nearly flat terrain.  For an outlet tip, all inward vertices are
        // upstream (higher), so the minimum inward elevation is still above the
        // tip, preserving a positive gradient.
        const auto& inward = endpoint_inward.at(key);
        for (const lsdtt_xy& v : inward)
        {
            const double e = sample_elev(v.x, v.y);
            if (!std::isnan(e))
            {
                if (std::isnan(ti.adjacent_elev) || e < ti.adjacent_elev)
                    ti.adjacent_elev = e;
            }
        }

        if (flow && !std::isnan(ti.elevation))
        {
            try
            {
                const int col = static_cast<int>(
                    std::floor((pt.x - meta.x_min) / meta.cell_size));
                const int row = static_cast<int>(
                    std::floor((meta.y_max() - pt.y) / meta.cell_size));
                const std::int64_t flat =
                    static_cast<std::int64_t>(row) * meta.n_cols + col;
                if (flow->is_flat_in_range(flat))
                    ti.upslope_count = flow->s_vector_size(flat);
            }
            catch (...) {}
        }

        tips.push_back(ti);
    }

    if (tips.empty())
    {
        std::cerr << "Warning: find_vector_network_source_tips: "
                  << "no degree-1 endpoints found in " << vector_path << "\n";
        return {};
    }

    std::cout << "find_vector_network_source_tips: found "
              << tips.size() << " degree-1 tip(s)\n";

    // ---------------------------------------------------------------------------
    // Outlet identification — segment gradient criterion only
    //
    // A tip is at the upstream (source) end of its segment when the first
    // valid inward vertex is lower than the tip (water flows away from the
    // tip down the segment).  A tip is at the downstream (outlet) end when
    // the first valid inward vertex is higher (water flows toward the tip).
    //
    // If no tip satisfies the downstream-end condition the real outlet is not
    // among the degree-1 tips — it is either at a junction (degree > 1, not
    // in this list) or entirely outside the DEM.  In that case has_outlet is
    // left false and every degree-1 tip is returned as a source.
    //
    // We deliberately do NOT fall back to minimum elevation or maximum upslope
    // count when the gradient criterion finds nothing.  Both of those fallbacks
    // force an outlet selection even when none exists, which removes a genuine
    // source tip from the source list.
    //
    // Among multiple downstream-end tips (unusual but possible in braided or
    // multi-outlet networks) the one with the greatest upslope_count is chosen;
    // steepest gradient breaks ties when upslope counts are equal.
    // ---------------------------------------------------------------------------
    int outlet_idx = -1;
    {
        std::uint64_t best_upslope   = 0;
        double        best_elev_diff = 0.0;

        for (int i = 0; i < static_cast<int>(tips.size()); ++i)
        {
            const tip_info& t = tips[i];
            if (std::isnan(t.elevation) || std::isnan(t.adjacent_elev)) continue;

            const double diff = t.adjacent_elev - t.elevation; // > 0 → tip is downstream
            if (diff <= 0.0) continue;

            const bool better_upslope = t.upslope_count > best_upslope;
            const bool same_upslope   = t.upslope_count == best_upslope;
            const bool better_grad    = diff > best_elev_diff;

            if (outlet_idx < 0 || better_upslope || (same_upslope && better_grad))
            {
                best_upslope   = t.upslope_count;
                best_elev_diff = diff;
                outlet_idx     = i;
            }
        }
        if (outlet_idx >= 0)
            std::cout << "  Outlet selected by segment gradient"
                      << (tips[outlet_idx].upslope_count > 0 ? " + upslope count" : "")
                      << "\n";
        else
            std::cout << "  No outlet tip detected among degree-1 endpoints "
                         "(outlet is at a junction or outside the DEM); "
                         "all tips treated as sources.\n";
    }

    // Print diagnostics for every tip to help with debugging.
    for (int i = 0; i < static_cast<int>(tips.size()); ++i)
    {
        const tip_info& t = tips[i];
        const bool is_outlet = (i == outlet_idx);
        std::cout << "  tip[" << i << "] "
                  << (is_outlet ? "[OUTLET] " : "[source] ")
                  << "(" << t.xy.x << ", " << t.xy.y << ")";
        if (!std::isnan(t.elevation))
            std::cout << "  elev=" << t.elevation;
        if (!std::isnan(t.adjacent_elev))
            std::cout << "  adj_elev=" << t.adjacent_elev;
        if (t.upslope_count > 0)
            std::cout << "  upslope=" << t.upslope_count;
        std::cout << "\n";
    }

    vector_network_endpoints result;
    for (int i = 0; i < static_cast<int>(tips.size()); ++i)
    {
        if (i == outlet_idx)
        {
            result.outlet_tip = tips[i].xy;
            result.has_outlet = true;
        }
        else
        {
            result.source_tips.push_back(tips[i].xy);
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// find_first_layer_junction_pairs
// ---------------------------------------------------------------------------
std::vector<vector_dem_junction_pair>
find_first_layer_junction_pairs(
    const std::string&                 vector_path,
    const lsdtt_raster_cpu&            dem,
    const vector_network_endpoints&    endpoints,
    const lsdtt_junction_network_base& dem_network)
{
    // ---- 1. Read vector layer -----------------------------------------------
    lsdtt_data_loader vec_loader;
    try { vec_loader.open_vector(vector_path); }
    catch (const std::exception& e)
    {
        throw std::runtime_error(
            std::string("find_first_layer_junction_pairs: cannot open vector: ")
            + e.what());
    }
    const auto layer_names = vec_loader.get_layer_names();
    if (layer_names.empty())
        throw std::runtime_error(
            "find_first_layer_junction_pairs: no layers in " + vector_path);

    const lsdtt_vector_layer vlayer =
        vec_loader.read_vector_layer(layer_names[0]);
    vec_loader.close_vector();

    const lsdtt_raster_meta& meta     = dem.get_meta();
    const lsdtt_array2d_f&   dem_data = dem.get_data_array();
    const lsdtt_array2d_u8&  dem_mask = dem.get_mask_array();

    const double snap_tol = meta.cell_size * 0.5;
    auto quantize = [snap_tol](double v) -> std::int64_t
    {
        return static_cast<std::int64_t>(std::round(v / snap_tol));
    };
    auto make_key = [&](double x, double y) -> std::string
    {
        return std::to_string(quantize(x)) + "_" + std::to_string(quantize(y));
    };

    // Helper: sample DEM elevation at (x, y); NaN when outside or nodata.
    auto sample_elev = [&](double x, double y) -> double
    {
        const int col = static_cast<int>(
            std::floor((x - meta.x_min) / meta.cell_size));
        const int row = static_cast<int>(
            std::floor((meta.y_max() - y) / meta.cell_size));
        if (row < 0 || row >= meta.n_rows ||
            col < 0 || col >= meta.n_cols ||
            dem_mask(row, col) == 0)
            return std::numeric_limits<double>::quiet_NaN();
        return static_cast<double>(dem_data(row, col));
    };

    // ---- 2. Build endpoint degree map and per-segment endpoint list ----------
    std::map<std::string, int>     ep_count;
    std::map<std::string, lsdtt_xy> ep_xy;

    // For each segment store {front_key, back_key}
    std::vector<std::pair<std::string, std::string>> seg_endpoints;
    seg_endpoints.reserve(vlayer.size());

    for (std::size_t fi = 0; fi < vlayer.size(); ++fi)
    {
        const auto& rec  = vlayer.records[fi];
        const auto* line = std::get_if<lsdtt_line_geom>(&rec.geometry);
        if (!line || line->points.size() < 2)
        {
            seg_endpoints.emplace_back("", "");
            continue;
        }
        const lsdtt_xy& front = line->points.front();
        const lsdtt_xy& back  = line->points.back();
        const std::string kf  = make_key(front.x, front.y);
        const std::string kb  = make_key(back.x,  back.y);
        ep_count[kf]++;
        ep_count[kb]++;
        ep_xy.emplace(kf, front);
        ep_xy.emplace(kb, back);
        seg_endpoints.emplace_back(kf, kb);
    }

    // ---- 3. Sample DEM elevation at every node ------------------------------
    std::map<std::string, double> ep_elev;
    for (const auto& [key, xy] : ep_xy)
        ep_elev[key] = sample_elev(xy.x, xy.y);

    // ---- 4. Build set of source-tip keys ------------------------------------
    // Degree-1 nodes that are NOT the outlet.
    std::string outlet_key;
    if (endpoints.has_outlet)
        outlet_key = make_key(endpoints.outlet_tip.x, endpoints.outlet_tip.y);

    std::set<std::string> source_tip_keys;
    for (const auto& [key, cnt] : ep_count)
    {
        if (cnt == 1 && key != outlet_key)
            source_tip_keys.insert(key);
    }

    // ---- 5. Build adjacency: junction key → list of segment indices ---------
    std::map<std::string, std::vector<std::size_t>> junc_segs;
    for (std::size_t fi = 0; fi < seg_endpoints.size(); ++fi)
    {
        const auto& [kf, kb] = seg_endpoints[fi];
        if (kf.empty()) continue;
        if (ep_count.at(kf) > 1) junc_segs[kf].push_back(fi);
        if (ep_count.at(kb) > 1) junc_segs[kb].push_back(fi);
    }

    // ---- 6. Find first-layer junctions --------------------------------------
    // A junction is first-layer when every upstream-connected segment (whose
    // far endpoint is higher than the junction) ends at a source tip.
    std::vector<vector_dem_junction_pair> result;

    for (const auto& [junc_key, segs] : junc_segs)
    {
        const double junc_elev = ep_elev.count(junc_key)
                                 ? ep_elev.at(junc_key)
                                 : std::numeric_limits<double>::quiet_NaN();
        if (std::isnan(junc_elev)) continue; // junction outside DEM

        // Classify connected segments as upstream or downstream.
        std::vector<std::string> upstream_other_keys;
        int n_downstream = 0;

        for (std::size_t fi : segs)
        {
            const auto& [kf, kb] = seg_endpoints[fi];
            const std::string other_key = (kf == junc_key) ? kb : kf;
            const double other_elev = ep_elev.count(other_key)
                                      ? ep_elev.at(other_key)
                                      : std::numeric_limits<double>::quiet_NaN();
            if (std::isnan(other_elev)) continue;

            if (other_elev > junc_elev)
                upstream_other_keys.push_back(other_key);
            else
                ++n_downstream;
        }

        // Need at least two upstream branches converging here.
        if (upstream_other_keys.size() < 2) continue;

        // First-layer: every upstream branch must end at a source tip.
        bool all_sources = true;
        for (const std::string& k : upstream_other_keys)
        {
            if (source_tip_keys.find(k) == source_tip_keys.end())
            {
                all_sources = false;
                break;
            }
        }
        if (!all_sources) continue;

        vector_dem_junction_pair pair;
        pair.pair_id   = static_cast<int>(result.size());
        pair.vector_xy = ep_xy.at(junc_key);
        pair.n_sources = static_cast<int>(upstream_other_keys.size());

        // -- Find nearest DEM junction ----------------------------------------
        double best_dist_sq = std::numeric_limits<double>::max();
        const std::int32_t n_dem_junc = dem_network.n_junctions();
        for (std::int32_t j = 0; j < n_dem_junc; ++j)
        {
            if (!dem_network.is_valid_junction(j)) continue;
            const lsdtt_xy dem_jxy = cell_center_from_row_col(
                meta,
                dem_network.row_of_junction(j),
                dem_network.col_of_junction(j));
            const double dx = dem_jxy.x - pair.vector_xy.x;
            const double dy = dem_jxy.y - pair.vector_xy.y;
            const double d2 = dx * dx + dy * dy;
            if (d2 < best_dist_sq)
            {
                best_dist_sq = d2;
                pair.dem_xy  = dem_jxy;
            }
        }
        pair.distance_m = std::sqrt(best_dist_sq);
        result.push_back(pair);
    }

    std::cout << "find_first_layer_junction_pairs: found "
              << result.size() << " first-layer junction pair(s)\n";
    return result;
}

// ---------------------------------------------------------------------------
// repair_junction_gaps  — internal helpers
// ---------------------------------------------------------------------------
namespace {

// Burn vector segments that have any vertex inside the bounding box
// (x_lo,x_hi,y_lo,y_hi) in world coordinates.
// Returns the number of DEM pixels actually lowered.
int burn_vector_in_window(
    lsdtt_raster_cpu&           dem,
    const std::string&          vector_path,
    double x_lo, double x_hi,
    double y_lo, double y_hi)
{
    lsdtt_data_loader loader;
    try { loader.open_vector(vector_path); }
    catch (...) { return 0; }

    const auto names = loader.get_layer_names();
    if (names.empty()) { loader.close_vector(); return 0; }

    const lsdtt_vector_layer vlayer = loader.read_vector_layer(names[0]);
    loader.close_vector();

    const lsdtt_raster_meta& meta = dem.get_meta();
    std::vector<float>&         elev = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();
    const int    nr = meta.n_rows;
    const int    nc = meta.n_cols;
    const double cs = meta.cell_size;
    const double perp_thr = static_cast<float>(std::sqrt(2.0) * 0.5 * cs);

    int total_burned = 0;

    for (std::size_t fi = 0; fi < vlayer.size(); ++fi)
    {
        const auto* line = std::get_if<lsdtt_line_geom>(&vlayer.records[fi].geometry);
        if (!line || line->points.size() < 2) continue;

        // Check if any vertex falls within the window
        bool in_window = false;
        for (const auto& pt : line->points)
        {
            if (pt.x >= x_lo && pt.x <= x_hi && pt.y >= y_lo && pt.y <= y_hi)
            { in_window = true; break; }
        }
        if (!in_window) continue;

        // Compute swath distances for this feature on the full DEM
        swath_distances swath;
        try { swath = compute_swath_distances(dem, line->points, 0); }
        catch (...) { continue; }

        const auto& along_d = swath.along_axis.get_data();
        const auto& along_m = swath.along_axis.get_mask();
        const auto& perp_d  = swath.perp_distance.get_data();

        // Collect pixels under the line that are also inside the window
        struct pnode { double along; int flat; };
        std::vector<pnode> nodes;

        for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const int idx = r * nc + c;
            if (!mask[idx] || !along_m[idx]) continue;
            if (perp_d[idx] > perp_thr) continue;

            const double x = meta.x_min + (c + 0.5) * cs;
            const double y = meta.y_min + (nr - r - 0.5) * cs;
            if (x < x_lo || x > x_hi || y < y_lo || y > y_hi) continue;

            nodes.push_back({static_cast<double>(along_d[idx]), idx});
        }
        if (nodes.empty()) continue;

        std::sort(nodes.begin(), nodes.end(),
                  [](const pnode& a, const pnode& b){ return a.along < b.along; });

        // Determine flow direction: lower end is outlet
        const double e_front = elev[nodes.front().flat];
        const double e_back  = elev[nodes.back().flat];

        // Build index sequence from upstream (high) to downstream (low)
        std::vector<int> seq;
        seq.reserve(nodes.size());
        if (e_front >= e_back)
            for (const auto& n : nodes) seq.push_back(n.flat);
        else
            for (int i = static_cast<int>(nodes.size())-1; i >= 0; --i)
                seq.push_back(nodes[i].flat);

        // Count pixels that will be lowered
        double rmin = elev[seq[0]];
        for (int f : seq)
        {
            const double e = elev[f];
            if (e < rmin) { rmin = e; }
            else if (e > rmin)
            {
                elev[f] = static_cast<float>(rmin);
                ++total_burned;
            }
        }
    }
    return total_burned;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// repair_junction_gaps
// ---------------------------------------------------------------------------
void repair_junction_gaps(
    lsdtt_raster_cpu&                      dem,
    std::vector<vector_dem_junction_pair>& pairs,
    const std::string&                     vector_path,
    double                                 min_gap_m,
    int                                    max_search_px,
    double                                 uphill_mult,
    double                                 max_climb_m)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int    nr = meta.n_rows;
    const int    nc = meta.n_cols;
    const double cs = meta.cell_size;

    const int goal_radius = std::max(3, max_search_px / 4);

    auto to_rc = [&](const lsdtt_xy& xy, int& r, int& c) -> bool
    {
        c = static_cast<int>((xy.x - meta.x_min) / cs);
        r = static_cast<int>((meta.y_max() - xy.y) / cs);
        return r >= 0 && r < nr && c >= 0 && c < nc;
    };

    int n_bridge = 0, n_vburn = 0, n_unresolved = 0;

    for (auto& pair : pairs)
    {
        if (pair.distance_m < min_gap_m)
        {
            pair.fix_type = GapFixType::ok;
            continue;
        }

        int rv, cv, rd, cd;
        if (!to_rc(pair.vector_xy, rv, cv) || !to_rc(pair.dem_xy, rd, cd))
        {
            std::cerr << "  pair " << pair.pair_id
                      << ": junction outside raster — skipping\n";
            pair.fix_type = GapFixType::unresolved;
            ++n_unresolved;
            continue;
        }

        // Search window: bounding box of V_j and D_j expanded by max_search_px
        lcp_params lcp;
        lcp.uphill_mult = uphill_mult;
        lcp.goal_radius = goal_radius;
        lcp.r_min = std::max(0,    std::min(rv, rd) - max_search_px);
        lcp.r_max = std::min(nr-1, std::max(rv, rd) + max_search_px);
        lcp.c_min = std::max(0,    std::min(cv, cd) - max_search_px);
        lcp.c_max = std::min(nc-1, std::max(cv, cd) + max_search_px);

        // ── Least-cost path: D_j → V_j neighbourhood ─────────────────────────
        lcp_result lcp_res;
        try
        {
            lcp_res = raster_least_cost_path(dem, rd, cd, rv, cv, lcp);
        }
        catch (const std::exception& e)
        {
            std::cerr << "  pair " << pair.pair_id
                      << ": LCP failed (" << e.what() << ")\n";
        }

        if (lcp_res.found && lcp_res.max_climb <= max_climb_m)
        {
            burn_monotonic_path(dem, lcp_res.path);
            pair.fix_type = GapFixType::curvature_bridge;
            pair.n_burned = static_cast<int>(lcp_res.path.size());
            ++n_bridge;
            std::cout << "  pair " << pair.pair_id
                      << ": bridge (" << pair.n_burned << " px"
                      << ", climb=" << lcp_res.max_climb << " m)\n";
            continue;
        }

        if (lcp_res.found)
            std::cout << "  pair " << pair.pair_id
                      << ": climb " << lcp_res.max_climb
                      << " m > " << max_climb_m << " m — vector burn fallback\n";
        else
            std::cout << "  pair " << pair.pair_id
                      << ": no terrain path — vector burn fallback\n";

        // ── Fallback: burn under vector segments within the window ────────────
        const double x_lo = meta.x_min + lcp.c_min * cs;
        const double x_hi = meta.x_min + (lcp.c_max + 1) * cs;
        const double y_lo = meta.y_min + (nr - lcp.r_max - 1) * cs;
        const double y_hi = meta.y_min + (nr - lcp.r_min) * cs;

        const int nb = burn_vector_in_window(dem, vector_path,
                                             x_lo, x_hi, y_lo, y_hi);
        if (nb > 0)
        {
            pair.fix_type = GapFixType::vector_burn;
            pair.n_burned = nb;
            ++n_vburn;
            std::cout << "  pair " << pair.pair_id
                      << ": vector burn (" << nb << " px)\n";
        }
        else
        {
            pair.fix_type = GapFixType::unresolved;
            ++n_unresolved;
            std::cout << "  pair " << pair.pair_id << ": unresolved (tunnel?)\n";
        }
    }

    std::cout << "repair_junction_gaps: bridge=" << n_bridge
              << "  vector_burn=" << n_vburn
              << "  unresolved=" << n_unresolved << "\n";
}

// ---------------------------------------------------------------------------
// write_junction_comparison_points
// ---------------------------------------------------------------------------
void write_junction_comparison_points(
    const std::vector<vector_dem_junction_pair>& pairs,
    const lsdtt_raster_meta&                      meta,
    const std::string&                            filename,
    point_output_format                           fmt,
    bool                                          spatial_index)
{
    if (pairs.empty()) return;

    lsdtt_vector_schema schema;
    schema.field_names = {"pair_id", "is_dem_junction", "distance_m",
                          "n_sources", "fix_type"};
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::string_type
    };

    const vector_output_format vfmt =
        (fmt == point_output_format::csv)
            ? vector_output_format::flatgeobuf
            : convert_point_to_vector_output_format(fmt);

    lsdtt_data_loader writer;
    writer.begin_vector_write(
        filename, "junction_comparison",
        schema, meta.crs_wkt,
        vfmt, spatial_index);

    for (const auto& p : pairs)
    {
        const std::string ft = gap_fix_type_str(p.fix_type);

        // Vector network junction
        writer.write_point(
            p.vector_xy.x, p.vector_xy.y,
            {std::int64_t{p.pair_id},
             std::int64_t{0},
             p.distance_m,
             std::int64_t{p.n_sources},
             ft});

        // DEM junction
        writer.write_point(
            p.dem_xy.x, p.dem_xy.y,
            {std::int64_t{p.pair_id},
             std::int64_t{1},
             p.distance_m,
             std::int64_t{p.n_sources},
             ft});
    }

    writer.finish_vector_write();
}

// ---------------------------------------------------------------------------
// resolve_stream_burn_mode
// ---------------------------------------------------------------------------
stream_burn_mode resolve_stream_burn_mode(const std::string& s)
{
    if (s == "burn_in_place")    return stream_burn_mode::burn_in_place;
    if (s == "snap_to_network")  return stream_burn_mode::snap_to_network;
    std::cerr << "Warning: unknown stream_burn_mode '" << s
              << "'; defaulting to burn_in_place.\n";
    return stream_burn_mode::burn_in_place;
}

// ---------------------------------------------------------------------------
// burn_streams_in_place — internal helpers
// ---------------------------------------------------------------------------

namespace {

// Apply a running-minimum scan from the upstream end (last element, highest
// flow_distance) toward the outlet (first element, flow_distance = 0).
// After this call every node satisfies:
//   burned_elevation[i] <= burned_elevation[i+1]   (non-increasing downstream)
//   burned_elevation[i] <= dem_elevation[i]         (never raised above DEM)
void apply_downhill_constraint(std::vector<burned_stream_node>& nodes)
{
    if (nodes.empty()) return;

    // Nodes are sorted ascending by flow_distance (index 0 = outlet).
    // Scan from upstream (n-1) toward outlet (0), tracking the running minimum.
    double running_min = nodes.back().dem_elevation;
    for (std::size_t i = nodes.size(); i-- > 0; )
    {
        running_min = std::min(nodes[i].dem_elevation, running_min);
        nodes[i].burned_elevation = running_min;
    }
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// burn_streams_in_place
// ---------------------------------------------------------------------------
std::vector<burned_stream_profile> burn_streams_in_place(
    const lsdtt_raster_cpu& dem,
    const std::string&      vector_path,
    int                     n_threads
)
{
    // -----------------------------------------------------------------------
    // Load vector features
    // -----------------------------------------------------------------------
    lsdtt_data_loader vec_loader;
    try
    {
        vec_loader.open_vector(vector_path);
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error(
            std::string("burn_streams_in_place: cannot open vector file: ")
            + e.what());
    }

    const auto layer_names = vec_loader.get_layer_names();
    if (layer_names.empty())
        throw std::runtime_error(
            "burn_streams_in_place: vector file has no layers: " + vector_path);

    const lsdtt_vector_layer vlayer =
        vec_loader.read_vector_layer(layer_names[0]);
    vec_loader.close_vector();

    if (vlayer.empty())
        throw std::runtime_error(
            "burn_streams_in_place: vector layer is empty: " + vector_path);

    // -----------------------------------------------------------------------
    // DEM metadata
    // -----------------------------------------------------------------------
    const lsdtt_raster_meta& meta      = dem.get_meta();
    const int                n_rows    = meta.n_rows;
    const int                n_cols    = meta.n_cols;
    const double             cell_size = meta.cell_size;

    const std::vector<float>&   dem_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask = dem.get_mask();

    // A pixel is considered "on the channel" when the perpendicular distance
    // to the polyline is within half the pixel diagonal.
    const float perp_threshold =
        static_cast<float>(std::sqrt(2.0) * 0.5 * cell_size);

    // -----------------------------------------------------------------------
    // Process each polyline feature independently
    // -----------------------------------------------------------------------
    std::vector<burned_stream_profile> profiles;
    profiles.reserve(vlayer.size());

    for (std::size_t fi = 0; fi < vlayer.size(); ++fi)
    {
        const auto& record = vlayer.records[fi];

        // Extract polyline vertices
        const lsdtt_line_geom* line_ptr =
            std::get_if<lsdtt_line_geom>(&record.geometry);
        if (!line_ptr)
        {
            std::cerr << "Warning: feature " << fi
                      << " is not a polyline; skipping.\n";
            continue;
        }
        if (line_ptr->points.size() < 2)
        {
            std::cerr << "Warning: feature " << fi
                      << " has fewer than 2 vertices; skipping.\n";
            continue;
        }

        // -------------------------------------------------------------------
        // Swath distances: along-axis and perpendicular for every DEM pixel
        // -------------------------------------------------------------------
        swath_distances swath;
        try
        {
            swath = compute_swath_distances(dem, line_ptr->points, n_threads);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: feature " << fi
                      << ": swath computation failed (" << e.what()
                      << "); skipping.\n";
            continue;
        }

        const std::vector<float>&   along_data = swath.along_axis.get_data();
        const std::vector<uint8_t>& along_mask = swath.along_axis.get_mask();
        const std::vector<float>&   perp_data  = swath.perp_distance.get_data();

        // -------------------------------------------------------------------
        // Collect pixels beneath the polyline
        // -------------------------------------------------------------------
        struct raw_node
        {
            double along_axis;
            double dem_elev;
            double x, y;
        };
        std::vector<raw_node> raw_nodes;

        for (int row = 0; row < n_rows; ++row)
        {
            for (int col = 0; col < n_cols; ++col)
            {
                const auto idx = static_cast<std::size_t>(row * n_cols + col);
                if (!dem_mask[idx])   continue;  // nodata
                if (!along_mask[idx]) continue;  // swath nodata (shouldn't happen)
                if (perp_data[idx] > perp_threshold) continue;  // not on channel

                const double x =
                    meta.x_min + (col + 0.5) * cell_size;
                const double y =
                    meta.y_min + (n_rows - row - 0.5) * cell_size;

                raw_nodes.push_back({
                    static_cast<double>(along_data[idx]),
                    static_cast<double>(dem_data[idx]),
                    x, y
                });
            }
        }

        if (raw_nodes.empty())
        {
            std::cerr << "Warning: feature " << fi
                      << ": no valid DEM pixels found beneath the polyline; "
                         "skipping.\n";
            continue;
        }

        // -------------------------------------------------------------------
        // Sort by along-axis distance (polyline start = 0)
        // -------------------------------------------------------------------
        std::sort(raw_nodes.begin(), raw_nodes.end(),
                  [](const raw_node& a, const raw_node& b)
                  { return a.along_axis < b.along_axis; });

        // -------------------------------------------------------------------
        // Determine which end is the outlet (lower mean elevation)
        // Compare mean elevation of the first and last 10% of nodes.
        // -------------------------------------------------------------------
        const std::size_t n_check =
            std::max(std::size_t{1}, raw_nodes.size() / 10);

        double sum_start = 0.0, sum_end = 0.0;
        for (std::size_t k = 0; k < n_check; ++k)
            sum_start += raw_nodes[k].dem_elev;
        for (std::size_t k = raw_nodes.size() - n_check;
             k < raw_nodes.size(); ++k)
            sum_end += raw_nodes[k].dem_elev;

        // outlet_at_end = true  → polyline end (high along_axis) is downstream
        // outlet_at_end = false → polyline start (along_axis ≈ 0) is downstream
        const bool outlet_at_end = (sum_end < sum_start);
        const double max_along   = raw_nodes.back().along_axis;

        // -------------------------------------------------------------------
        // Build profile nodes with flow_distance (0 = outlet, increasing ↑)
        // -------------------------------------------------------------------
        burned_stream_profile profile;
        profile.feature_id = static_cast<std::int32_t>(fi);
        profile.nodes.reserve(raw_nodes.size());

        for (const auto& rn : raw_nodes)
        {
            burned_stream_node node;
            node.flow_distance =
                outlet_at_end ? (max_along - rn.along_axis) : rn.along_axis;
            node.dem_elevation  = rn.dem_elev;
            node.burned_elevation = 0.0;  // filled below
            node.x = rn.x;
            node.y = rn.y;
            profile.nodes.push_back(node);
        }

        // Sort ascending by flow_distance so index 0 = outlet
        std::sort(profile.nodes.begin(), profile.nodes.end(),
                  [](const burned_stream_node& a, const burned_stream_node& b)
                  { return a.flow_distance < b.flow_distance; });

        // -------------------------------------------------------------------
        // Apply downhill constraint: running minimum from upstream to outlet
        // -------------------------------------------------------------------
        apply_downhill_constraint(profile.nodes);

        std::cout << "  Feature " << fi << ": " << profile.nodes.size()
                  << " channel pixels, outlet elevation = "
                  << profile.nodes.front().dem_elevation
                  << " m, headwater elevation = "
                  << profile.nodes.back().dem_elevation << " m\n";

        profiles.push_back(std::move(profile));
    }

    return profiles;
}

// ---------------------------------------------------------------------------
// write_burned_stream_profiles_csv
// ---------------------------------------------------------------------------
void write_burned_stream_profiles_csv(
    const std::vector<burned_stream_profile>& profiles,
    const std::string&                        filename
)
{
    std::ofstream ofs(filename);
    if (!ofs)
        throw std::runtime_error(
            "write_burned_stream_profiles_csv: cannot open output file: "
            + filename);

    ofs << "feature_id,flow_distance,x,y,dem_elevation,burned_elevation\n";
    ofs.precision(6);
    ofs << std::fixed;

    for (const auto& profile : profiles)
    {
        for (const auto& node : profile.nodes)
        {
            ofs << profile.feature_id << ','
                << node.flow_distance    << ','
                << node.x               << ','
                << node.y               << ','
                << node.dem_elevation   << ','
                << node.burned_elevation << '\n';
        }
    }
}

// ---------------------------------------------------------------------------
// write_burned_stream_profiles_points
// ---------------------------------------------------------------------------
void write_burned_stream_profiles_points(
    const std::vector<burned_stream_profile>& profiles,
    const lsdtt_raster_meta&                  meta,
    const std::string&                        filename,
    point_output_format                       fmt,
    bool                                      spatial_index
)
{
    if (fmt == point_output_format::csv)
        throw std::invalid_argument(
            "write_burned_stream_profiles_points: use "
            "write_burned_stream_profiles_csv for CSV output.");

    lsdtt_vector_schema schema;
    schema.field_names = {
        "feature_id", "flow_distance", "dem_elevation", "burned_elevation"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };

    const vector_output_format vfmt =
        convert_point_to_vector_output_format(fmt);

    lsdtt_data_loader writer;
    writer.begin_vector_write(filename, "burned_stream_profiles",
                              schema, meta.crs_wkt, vfmt, spatial_index);

    for (const auto& profile : profiles)
    {
        for (const auto& node : profile.nodes)
        {
            writer.write_point(
                node.x, node.y,
                {
                    static_cast<std::int64_t>(profile.feature_id),
                    node.flow_distance,
                    node.dem_elevation,
                    node.burned_elevation
                }
            );
        }
    }

    writer.finish_vector_write();
}

// ===========================================================================
// build_vector_networks helpers
// ===========================================================================
namespace {

inline bool bbox_inside(double x, double y,
                        double x0, double x1, double y0, double y1) noexcept
{
    return x >= x0 && x <= x1 && y >= y0 && y <= y1;
}

// Segment [a→b]: a is inside the bbox, b is outside.
// Returns the intersection point on the bbox edge.
lsdtt_xy interpolate_to_bbox_edge(lsdtt_xy a, lsdtt_xy b,
                                   double x0, double x1, double y0, double y1)
{
    double t = 1.0;
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    if (dx >  1e-12) t = std::min(t, (x1 - a.x) / dx);
    if (dx < -1e-12) t = std::min(t, (x0 - a.x) / dx);
    if (dy >  1e-12) t = std::min(t, (y1 - a.y) / dy);
    if (dy < -1e-12) t = std::min(t, (y0 - a.y) / dy);
    t = std::max(0.0, std::min(1.0, t));
    return {a.x + t * dx, a.y + t * dy};
}

struct BboxClipResult { lsdtt_xy pt; bool valid = false; bool clipped = false; };

// Clip the FRONT (index-0) endpoint of a polyline to the DEM bbox.
BboxClipResult clip_poly_front(const std::vector<lsdtt_xy>& pts,
                                double x0, double x1, double y0, double y1)
{
    if (pts.empty()) return {};
    if (bbox_inside(pts[0].x, pts[0].y, x0, x1, y0, y1))
        return {pts[0], true, false};
    for (std::size_t i = 1; i < pts.size(); ++i)
        if (bbox_inside(pts[i].x, pts[i].y, x0, x1, y0, y1))
            return {interpolate_to_bbox_edge(pts[i], pts[i-1], x0, x1, y0, y1), true, true};
    return {};
}

// Clip the BACK (last-index) endpoint of a polyline to the DEM bbox.
BboxClipResult clip_poly_back(const std::vector<lsdtt_xy>& pts,
                               double x0, double x1, double y0, double y1)
{
    if (pts.empty()) return {};
    const int N = static_cast<int>(pts.size()) - 1;
    if (bbox_inside(pts[N].x, pts[N].y, x0, x1, y0, y1))
        return {pts[N], true, false};
    for (int i = N - 1; i >= 0; --i)
        if (bbox_inside(pts[i].x, pts[i].y, x0, x1, y0, y1))
            return {interpolate_to_bbox_edge(pts[i], pts[i+1], x0, x1, y0, y1), true, true};
    return {};
}

} // anonymous namespace

// ===========================================================================
// build_vector_networks
// ===========================================================================
std::vector<VectorNetwork>
build_vector_networks(const std::string& vector_path, const lsdtt_raster_cpu& dem)
{
    // ---- Load vector -------------------------------------------------------
    lsdtt_data_loader loader;
    try { loader.open_vector(vector_path); }
    catch (const std::exception& e)
    {
        throw std::runtime_error(
            std::string("build_vector_networks: cannot open vector: ") + e.what());
    }
    const auto layer_names = loader.get_layer_names();
    if (layer_names.empty())
        throw std::runtime_error("build_vector_networks: no layers in " + vector_path);
    const lsdtt_vector_layer vlayer = loader.read_vector_layer(layer_names[0]);
    loader.close_vector();

    // ---- DEM extent and helpers --------------------------------------------
    const lsdtt_raster_meta& meta = dem.get_meta();
    const double x0 = meta.x_min;
    const double x1 = meta.x_min + meta.n_cols * meta.cell_size;
    const double y0 = meta.y_min;
    const double y1 = meta.y_min + meta.n_rows * meta.cell_size;
    const double cs = meta.cell_size;

    const lsdtt_array2d_f&  dem_data = dem.get_data_array();
    const lsdtt_array2d_u8& dem_mask = dem.get_mask_array();
    auto sample_elev = [&](double x, double y) -> double {
        const int col = static_cast<int>(std::floor((x - meta.x_min) / cs));
        const int row = static_cast<int>(std::floor((meta.y_max() - y) / cs));
        if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols
            || dem_mask(row, col) == 0)
            return std::numeric_limits<double>::quiet_NaN();
        return static_cast<double>(dem_data(row, col));
    };

    // Quantised key: snap endpoints to half-cell grid for robust matching
    const double snap = cs * 0.5;
    auto make_key = [snap](double x, double y) -> std::string {
        return std::to_string(static_cast<std::int64_t>(std::round(x / snap)))
             + "_"
             + std::to_string(static_cast<std::int64_t>(std::round(y / snap)));
    };

    // ---- Clip each polyline to DEM bbox ------------------------------------
    struct PolyInfo {
        lsdtt_xy front, back;
        bool front_clipped = false;
        bool back_clipped  = false;
        bool valid         = false;
    };
    std::vector<PolyInfo> polys(vlayer.size());

    for (std::size_t fi = 0; fi < vlayer.size(); ++fi)
    {
        const auto* line = std::get_if<lsdtt_line_geom>(&vlayer.records[fi].geometry);
        if (!line || line->points.size() < 2) continue;
        const auto fe = clip_poly_front(line->points, x0, x1, y0, y1);
        const auto be = clip_poly_back (line->points, x0, x1, y0, y1);
        if (!fe.valid || !be.valid) continue;
        polys[fi] = {fe.pt, be.pt, fe.clipped, be.clipped, true};
    }

    // ---- Build endpoint adjacency map --------------------------------------
    // ep_adj[key] = {(poly_index, is_back_end), ...}
    std::unordered_map<std::string, std::vector<std::pair<int,bool>>> ep_adj;
    std::unordered_map<std::string, lsdtt_xy>                         ep_xy;
    std::unordered_map<std::string, bool>                              ep_clipped;

    for (int i = 0; i < static_cast<int>(polys.size()); ++i)
    {
        if (!polys[i].valid) continue;
        const std::string kf = make_key(polys[i].front.x, polys[i].front.y);
        const std::string kb = make_key(polys[i].back.x,  polys[i].back.y);
        ep_adj[kf].emplace_back(i, false);
        ep_adj[kb].emplace_back(i, true);
        ep_xy[kf]      = polys[i].front;
        ep_xy[kb]      = polys[i].back;
        ep_clipped[kf] = polys[i].front_clipped;
        ep_clipped[kb] = polys[i].back_clipped;
    }

    // ---- Connected components (BFS over endpoint keys) ---------------------
    std::unordered_map<std::string, int>      key_comp;
    std::vector<std::vector<std::string>>     comp_keys;

    for (const auto& [seed, _] : ep_adj)
    {
        if (key_comp.count(seed)) continue;
        const int cid = static_cast<int>(comp_keys.size());
        comp_keys.emplace_back();
        std::queue<std::string> bfsq;
        bfsq.push(seed);
        key_comp[seed] = cid;
        while (!bfsq.empty())
        {
            const std::string cur = bfsq.front(); bfsq.pop();
            comp_keys[cid].push_back(cur);
            for (const auto& [pi, is_back] : ep_adj.at(cur))
            {
                const auto& p = polys[pi];
                const std::string other = is_back
                    ? make_key(p.front.x, p.front.y)
                    : make_key(p.back.x,  p.back.y);
                if (!key_comp.count(other))
                {
                    key_comp[other] = cid;
                    bfsq.push(other);
                }
            }
        }
    }

    // ---- Helper: follow chain of degree-2 nodes to next graph node ---------
    // Returns the key of the first node reached with degree != 2.
    auto follow_chain = [&](const std::string& start,
                             const std::string& from) -> std::string
    {
        std::string cur = start, prev = from;
        for (int step = 0; step < 100000; ++step)
        {
            const auto it = ep_adj.find(cur);
            if (it == ep_adj.end() || it->second.size() != 2) break;
            std::string next;
            for (const auto& [pi, is_back] : it->second)
            {
                const auto& p = polys[pi];
                const std::string other = is_back
                    ? make_key(p.front.x, p.front.y)
                    : make_key(p.back.x,  p.back.y);
                if (other != prev) { next = other; break; }
            }
            if (next.empty() || next == cur) break;
            prev = cur;
            cur  = next;
        }
        return cur;
    };

    // ---- Build one VectorNetwork per connected component -------------------
    std::vector<VectorNetwork> result;

    for (int cid = 0; cid < static_cast<int>(comp_keys.size()); ++cid)
    {
        const auto& all_keys = comp_keys[cid];

        // Graph nodes are keys with degree != 2
        std::vector<std::string> graph_keys;
        for (const auto& k : all_keys)
            if (ep_adj.at(k).size() != 2) graph_keys.push_back(k);

        if (graph_keys.empty()) continue;

        // Assign node IDs
        std::unordered_map<std::string, int> key_nid;
        VectorNetwork net;
        for (const auto& k : graph_keys)
        {
            const int nid = static_cast<int>(net.nodes.size());
            key_nid[k]   = nid;
            VectorNetworkNode node;
            node.node_id   = nid;
            node.xy        = ep_xy.at(k);
            node.is_clipped = ep_clipped.count(k) && ep_clipped.at(k);
            net.nodes.push_back(node);
        }

        // Build undirected graph edges (one per unique pair via chain following)
        struct UEdge { int a, b; };
        std::vector<UEdge> uedges;
        std::set<std::pair<int,int>> seen;

        for (const auto& k : graph_keys)
        {
            const int nid_a = key_nid.at(k);
            for (const auto& [pi, is_back] : ep_adj.at(k))
            {
                const auto& p = polys[pi];
                const std::string imm = is_back
                    ? make_key(p.front.x, p.front.y)
                    : make_key(p.back.x,  p.back.y);
                const std::string far = follow_chain(imm, k);
                if (!key_nid.count(far)) continue;
                const int nid_b = key_nid.at(far);
                if (nid_a == nid_b) continue;
                const auto ek = std::make_pair(std::min(nid_a, nid_b),
                                               std::max(nid_a, nid_b));
                if (!seen.insert(ek).second) continue;
                uedges.push_back({nid_a, nid_b});
            }
        }

        // Assign flow direction to each edge
        // Clipped nodes are at the downstream (outlet) end of their edges.
        // Otherwise compare DEM elevation: higher endpoint is upstream.
        struct DEdge { int up, dn; };
        std::vector<DEdge> dedges;
        dedges.reserve(uedges.size());

        for (const auto& ue : uedges)
        {
            const auto& na = net.nodes[ue.a];
            const auto& nb = net.nodes[ue.b];
            int up = -1, dn = -1;

            if (na.is_clipped && !nb.is_clipped)        { up = ue.b; dn = ue.a; }
            else if (nb.is_clipped && !na.is_clipped)   { up = ue.a; dn = ue.b; }
            else
            {
                const double ea = sample_elev(na.xy.x, na.xy.y);
                const double eb = sample_elev(nb.xy.x, nb.xy.y);
                if (!std::isnan(ea) && !std::isnan(eb))
                {
                    if (ea >= eb) { up = ue.a; dn = ue.b; }
                    else          { up = ue.b; dn = ue.a; }
                }
                else if (!std::isnan(ea)) { up = ue.a; dn = ue.b; }
                else                      { up = ue.b; dn = ue.a; }
            }

            if (up >= 0) dedges.push_back({up, dn});
        }

        // Populate directed adjacency on nodes
        for (const auto& de : dedges)
        {
            net.nodes[de.up].downstream_node_id = de.dn;
            net.nodes[de.dn].upstream_node_ids.push_back(de.up);
        }

        // Classify source / outlet
        for (auto& node : net.nodes)
        {
            if (node.upstream_node_ids.empty() && node.downstream_node_id >= 0)
                node.is_source = true;
            if (node.downstream_node_id == -1 && !node.upstream_node_ids.empty())
                node.is_outlet = true;
        }

        // Identify THE outlet (prefer clipped, then lowest elevation)
        {
            std::vector<int> outlets;
            for (const auto& n : net.nodes)
                if (n.is_outlet) outlets.push_back(n.node_id);

            if (outlets.empty())
            {
                std::cerr << "build_vector_networks: component " << cid
                          << " has no outlet; skipping\n";
                continue;
            }

            int best = outlets[0];
            for (int oi : outlets)
            {
                if (net.nodes[oi].is_clipped)  { best = oi; break; }
                const double ea = sample_elev(net.nodes[best].xy.x, net.nodes[best].xy.y);
                const double eb = sample_elev(net.nodes[oi].xy.x,   net.nodes[oi].xy.y);
                if (!std::isnan(eb) && (std::isnan(ea) || eb < ea)) best = oi;
            }
            net.outlet_node_id = best;
        }

        // Kahn's topological sort — collect junctions in upstream-first order
        std::unordered_map<int,int> in_deg;
        for (const auto& n : net.nodes)
            in_deg[n.node_id] = static_cast<int>(n.upstream_node_ids.size());

        std::queue<int> kq;
        for (const auto& n : net.nodes)
            if (in_deg[n.node_id] == 0) kq.push(n.node_id);

        while (!kq.empty())
        {
            const int cur = kq.front(); kq.pop();
            const auto& cn = net.nodes[cur];
            if (!cn.is_source && !cn.is_outlet)
                net.junction_topo_order.push_back(cur);
            if (cn.downstream_node_id >= 0)
                if (--in_deg.at(cn.downstream_node_id) == 0)
                    kq.push(cn.downstream_node_id);
        }

        result.push_back(std::move(net));
    }

    std::cout << "build_vector_networks: " << result.size() << " network(s)\n";
    for (int i = 0; i < static_cast<int>(result.size()); ++i)
    {
        int ns = 0, nj = 0;
        for (const auto& n : result[i].nodes)
            if (n.is_source) ++ns;
            else if (!n.is_outlet) ++nj;
        std::cout << "  network " << i << ": " << ns << " source(s), "
                  << nj << " junction(s)\n";
    }
    return result;
}

// ===========================================================================
// repair_stream_junctions
// ===========================================================================
std::vector<JunctionRepairRecord>
repair_stream_junctions(
    lsdtt_raster_cpu&                             dem,
    const std::string&                            vector_path,
    const std::vector<std::int64_t>&              source_flats,
    std::unique_ptr<lsdtt_flow_base>&             flow,
    std::unique_ptr<lsdtt_junction_network_base>& dem_network,
    const lsdtt_flowinfo_options&                 fi_opts,
    const lsdtt_junction_network_options&         jn_opts,
    const junction_repair_params&                 params)
{
    const std::vector<VectorNetwork> networks = build_vector_networks(vector_path, dem);
    const lsdtt_raster_meta& meta = dem.get_meta();

    // ---- DEM headwater helpers ---------------------------------------------
    struct DemHW { std::int32_t junc_id; lsdtt_xy xy; };

    auto build_headwaters = [&]() -> std::vector<DemHW>
    {
        std::vector<DemHW> hw;
        const std::int32_t n = dem_network->n_junctions();
        for (std::int32_t j = 0; j < n; ++j)
        {
            if (!dem_network->is_valid_junction(j)) continue;
            if (dem_network->donor_junction_size(j) != 0) continue;
            const lsdtt_xy jxy = cell_center_from_row_col(
                meta, dem_network->row_of_junction(j), dem_network->col_of_junction(j));
            hw.push_back({j, jxy});
        }
        return hw;
    };

    auto nearest_headwater = [](const lsdtt_xy& pt,
                                 const std::vector<DemHW>& hw) -> std::int32_t
    {
        double best = std::numeric_limits<double>::max();
        std::int32_t best_j = -1;
        for (const auto& h : hw)
        {
            const double d2 = (h.xy.x - pt.x) * (h.xy.x - pt.x)
                            + (h.xy.y - pt.y) * (h.xy.y - pt.y);
            if (d2 < best) { best = d2; best_j = h.junc_id; }
        }
        return best_j;
    };

    // LCA of two junctions: walk ja to outlet collecting ancestors,
    // then walk jb until hitting one.
    auto lca_two = [&](std::int32_t ja, std::int32_t jb) -> std::int32_t
    {
        std::unordered_set<std::int32_t> seen;
        for (std::int32_t cur = ja; cur >= 0; )
        {
            seen.insert(cur);
            const std::int32_t recv = dem_network->receiver_of_junction(cur);
            if (recv < 0 || recv == cur) break;
            cur = recv;
        }
        for (std::int32_t cur = jb; cur >= 0; )
        {
            if (seen.count(cur)) return cur;
            const std::int32_t recv = dem_network->receiver_of_junction(cur);
            if (recv < 0 || recv == cur) break;
            cur = recv;
        }
        return -1;
    };

    // Recursive DEM-junction finder via LCA.
    // cache maps vector node_id → DEM junction_id for the CURRENT extraction.
    std::function<std::int32_t(const VectorNetwork&, int,
                               const std::vector<DemHW>&,
                               std::unordered_map<int,std::int32_t>&)>
    find_dem_junction;

    find_dem_junction = [&](const VectorNetwork&                  net,
                             int                                   node_id,
                             const std::vector<DemHW>&             hw,
                             std::unordered_map<int,std::int32_t>& cache)
        -> std::int32_t
    {
        auto it = cache.find(node_id);
        if (it != cache.end()) return it->second;

        const auto& node = net.nodes[node_id];
        std::int32_t result = -1;

        if (node.is_source)
        {
            result = nearest_headwater(node.xy, hw);
        }
        else
        {
            std::vector<std::int32_t> ups;
            for (int uid : node.upstream_node_ids)
            {
                const std::int32_t uj = find_dem_junction(net, uid, hw, cache);
                if (uj >= 0) ups.push_back(uj);
            }
            if (!ups.empty())
            {
                result = ups[0];
                for (std::size_t i = 1; i < ups.size(); ++i)
                {
                    result = lca_two(result, ups[i]);
                    if (result < 0) break;
                }
            }
        }

        cache[node_id] = result;
        return result;
    };

    // ---- Main repair loop --------------------------------------------------
    std::vector<JunctionRepairRecord> records;
    int repair_index = 0;
    auto dem_hw = build_headwaters();

    for (int ni = 0; ni < static_cast<int>(networks.size()); ++ni)
    {
        const auto& net = networks[ni];
        std::unordered_map<int,std::int32_t> cache;

        for (int junc_nid : net.junction_topo_order)
        {
            const auto& junc_node = net.nodes[junc_nid];

            // Find DEM counterpart
            const std::int32_t dem_junc =
                find_dem_junction(net, junc_nid, dem_hw, cache);

            if (dem_junc < 0)
            {
                std::cout << "  [net " << ni << " junc " << junc_nid
                          << "]: no DEM counterpart found; skipping\n";
                continue;
            }

            const lsdtt_xy dem_jxy = cell_center_from_row_col(
                meta,
                dem_network->row_of_junction(dem_junc),
                dem_network->col_of_junction(dem_junc));

            const double dx   = dem_jxy.x - junc_node.xy.x;
            const double dy   = dem_jxy.y - junc_node.xy.y;
            const double dist = std::sqrt(dx * dx + dy * dy);

            std::cout << "  [net " << ni << " junc " << junc_nid
                      << "]: gap=" << dist << " m";

            if (dist < params.min_gap_m)
            {
                std::cout << " (within tolerance)\n";
                continue;
            }

            std::cout << " — repairing\n";

            // Single-pair repair via existing repair_junction_gaps
            vector_dem_junction_pair pair;
            pair.pair_id    = repair_index;
            pair.vector_xy  = junc_node.xy;
            pair.dem_xy     = dem_jxy;
            pair.distance_m = dist;
            pair.n_sources  = static_cast<int>(junc_node.upstream_node_ids.size());
            std::vector<vector_dem_junction_pair> pv = {pair};

            repair_junction_gaps(dem, pv, vector_path,
                                 params.min_gap_m,
                                 params.max_search_px,
                                 params.uphill_mult,
                                 params.max_climb_m);

            JunctionRepairRecord rec;
            rec.repair_index  = repair_index++;
            rec.network_index = ni;
            rec.vector_xy     = junc_node.xy;
            rec.dem_xy        = dem_jxy;
            rec.distance_m    = dist;
            rec.fix_type      = pv[0].fix_type;
            rec.n_burned      = pv[0].n_burned;
            records.push_back(rec);

            std::cout << "    → " << gap_fix_type_str(pv[0].fix_type)
                      << " (" << pv[0].n_burned << " px)"
                      << "  re-extracting network...\n";

            // Re-extract flow routing and junction network
            flow       = make_flow_base(dem.get_data_array(), dem.get_mask_array(),
                                        meta, fi_opts);
            dem_network = make_junction_network_base(*flow, source_flats, jn_opts);

            // Refresh headwaters and invalidate cache
            dem_hw = build_headwaters();
            cache.clear();
        }
    }

    std::cout << "repair_stream_junctions: " << repair_index
              << " repair(s) applied across "
              << networks.size() << " network(s)\n";
    return records;
}

// ===========================================================================
// write_junction_repair_log
// ===========================================================================
void write_junction_repair_log(
    const std::vector<JunctionRepairRecord>& records,
    const std::string&                       filename)
{
    std::ofstream f(filename);
    if (!f)
        throw std::runtime_error(
            "write_junction_repair_log: cannot open " + filename);

    f << "repair_index,network_index,vector_x,vector_y,"
         "dem_x,dem_y,distance_m,fix_type,n_burned\n";

    for (const auto& r : records)
    {
        f << r.repair_index  << ","
          << r.network_index << ","
          << r.vector_xy.x   << ","
          << r.vector_xy.y   << ","
          << r.dem_xy.x      << ","
          << r.dem_xy.y      << ","
          << r.distance_m    << ","
          << gap_fix_type_str(r.fix_type) << ","
          << r.n_burned      << "\n";
    }
}

} // namespace lsdtt
