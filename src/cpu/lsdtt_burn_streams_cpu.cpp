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
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
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

} // namespace lsdtt
