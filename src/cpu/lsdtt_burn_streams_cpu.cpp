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
#include <iomanip>
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
#include "lsdtt_raster_feature_finder_cpu.hpp"

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





struct dem_candidate_junction_view
{
    int      junction_id = -1;
    lsdtt_xy xy;
    double   elevation = 0.0;
    int      n_donors = 0;
};

double resolved_candidate_match_radius(
    const candidate_junction_search_params& params,
    const lsdtt_raster_meta&                meta)
{
    if (params.match_radius_m > 0.0)
    {
        return params.match_radius_m;
    }

    return 2.0 * meta.cell_size;
}

double planar_distance(
    const lsdtt_xy& a,
    const lsdtt_xy& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

std::vector<dem_candidate_junction_view>
build_dem_candidate_junction_views(
    const lsdtt_junction_network_base& dem_network,
    const lsdtt_raster_cpu&            dem)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dem_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask = dem.get_mask();

    auto flat_of = [&](int row, int col) -> std::int64_t
    {
        return static_cast<std::int64_t>(row) * meta.n_cols + col;
    };

    std::vector<dem_candidate_junction_view> out;
    out.reserve(static_cast<std::size_t>(dem_network.n_junctions()));

    for (std::int32_t j = 0; j < dem_network.n_junctions(); ++j)
    {
        if (!dem_network.is_valid_junction(j)) continue;

        const int row = dem_network.row_of_junction(j);
        const int col = dem_network.col_of_junction(j);

        const std::int64_t flat = flat_of(row, col);
        if (flat < 0 || flat >= static_cast<std::int64_t>(dem_data.size()))
            continue;

        if (!dem_mask[static_cast<std::size_t>(flat)])
            continue;

        if (meta.is_nodata(dem_data[static_cast<std::size_t>(flat)]))
            continue;

        if (!std::isfinite(dem_data[static_cast<std::size_t>(flat)]))
            continue;

        dem_candidate_junction_view v;
        v.junction_id = j;
        v.xy = cell_center_from_row_col(meta, row, col);
        v.elevation = static_cast<double>(
            dem_data[static_cast<std::size_t>(flat)]
        );
        v.n_donors = static_cast<int>(dem_network.donor_junction_size(j));

        out.push_back(v);
    }

    return out;
}

std::vector<dem_candidate_junction_view>
filter_dem_headwaters(
    const std::vector<dem_candidate_junction_view>& dem_views)
{
    std::vector<dem_candidate_junction_view> out;
    out.reserve(dem_views.size());

    for (const auto& v : dem_views)
    {
        if (v.n_donors == 0)
        {
            out.push_back(v);
        }
    }

    return out;
}

int find_nearest_dem_junction(
    const lsdtt_xy&                                 xy,
    const std::vector<dem_candidate_junction_view>& dem_views,
    double&                                         best_distance_m)
{
    int best_id = -1;
    best_distance_m = std::numeric_limits<double>::infinity();

    for (const auto& dj : dem_views)
    {
        const double d = planar_distance(xy, dj.xy);

        if (d < best_distance_m)
        {
            best_distance_m = d;
            best_id = dj.junction_id;
        }
    }

    if (best_id < 0)
    {
        best_distance_m = -1.0;
    }

    return best_id;
}

int lca_two_dem_junctions(
    const lsdtt_junction_network_base& dem_network,
    int                                ja,
    int                                jb)
{
    if (ja < 0 || jb < 0) return -1;
    if (!dem_network.is_valid_junction(ja)) return -1;
    if (!dem_network.is_valid_junction(jb)) return -1;

    std::unordered_set<int> ancestors;
    ancestors.reserve(static_cast<std::size_t>(dem_network.n_junctions()));

    int cur = ja;
    for (int safety = 0; safety <= dem_network.n_junctions(); ++safety)
    {
        if (cur < 0) break;
        if (!dem_network.is_valid_junction(cur)) break;

        ancestors.insert(cur);

        const int rec = dem_network.receiver_of_junction(cur);
        if (rec < 0 || rec == cur) break;
        cur = rec;
    }

    cur = jb;
    for (int safety = 0; safety <= dem_network.n_junctions(); ++safety)
    {
        if (cur < 0) break;
        if (!dem_network.is_valid_junction(cur)) break;

        if (ancestors.find(cur) != ancestors.end())
        {
            return cur;
        }

        const int rec = dem_network.receiver_of_junction(cur);
        if (rec < 0 || rec == cur) break;
        cur = rec;
    }

    return -1;
}

int lca_many_dem_junctions(
    const lsdtt_junction_network_base& dem_network,
    const std::vector<int>&            junction_ids)
{
    if (junction_ids.empty()) return -1;

    int lca = junction_ids.front();
    for (std::size_t i = 1; i < junction_ids.size(); ++i)
    {
        lca = lca_two_dem_junctions(dem_network, lca, junction_ids[i]);
        if (lca < 0) return -1;
    }

    return lca;
}








    
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
// find_candidate_junctions_phase
//
// Finds the most upstream mismatched vector junction in each basin by
// comparing vector-network topology against the DEM-extracted junction
// network.
//
// Approach
// --------
// This is a TOPOLOGICAL matching algorithm, not a nearest-junction search.
//
// The DEM network is seeded from source junctions derived from the vector
// network, so after correct conditioning the DEM network should reproduce the
// same branching structure as the vector network. The right DEM counterpart
// of a vector junction is therefore not simply the closest DEM junction in
// space. Instead, it is the junction where the DEM branches corresponding to
// that vector junction's immediate donor branches first meet.
//
// The algorithm works basin-by-basin, where a basin is identified by the
// vector outlet_junction key.
//
// Per-basin workflow
// ------------------
// 1. Initialise source mappings
//    Every vector source junction in the basin is mapped to the nearest DEM
//    headwater junction (DEM junction with zero donors). This is justified
//    because the DEM network was extracted using vector-derived source seeds.
//
// 2. Walk vector junctions in upstream -> downstream topological order
//    By the time a junction is visited, all of its immediate donor vector
//    junctions should already have DEM counterparts assigned.
//
// 3. For each non-source vector junction:
//    - gather the DEM counterparts of its immediate donor vector junctions
//    - compute the DEM lowest common ancestor (LCA) of those donor
//      DEM junctions
//
//    That LCA is the topologically correct DEM counterpart of the vector
//    junction, because it is the first place downstream where those DEM
//    branches unite.
//
// 4. Test for mismatch
//    A vector junction is flagged as a candidate mismatch when any of the
//    following occur:
//
//      donor_mapping_missing
//        One or more donor vector junctions does not yet have a DEM
//        counterpart.
//
//      no_common_dem_ancestor
//        The DEM counterparts of the donor branches do not have a valid
//        common downstream ancestor.
//
//      dem_branches_join_too_early
//        The DEM LCA is itself one of the donor DEM junctions, indicating
//        that branches merged too early, or multiple vector donors collapsed
//        onto the same DEM branch.
//
//      donor_count_mismatch
//        The topologically matched DEM junction exists, but its number of
//        donor junctions differs from that of the vector junction.
//
//      dem_junction_too_far
//        The DEM donor branches do join, but the location of their join
//        (the LCA junction) is farther from the vector junction than the
//        allowed tolerance.
//
// 5. Candidate selection
//    The search stops at the FIRST mismatched candidate-type junction in each
//    basin. "Candidate-type" means vector junctions with at least
//    minimum_vector_donor_junctions donors, which in the first implementation
//    usually means true confluences only.
//
// 6. Accepted matches propagate downstream
//    If a vector junction passes all checks, its DEM LCA is stored as the
//    accepted DEM counterpart for use when processing downstream vector
//    junctions in the same basin.
//
// Why only one candidate per basin per phase?
// -------------------------------------------
// Repairing a junction changes DEM flow routing downstream, which can alter
// DEM junction positions and even the identity of downstream LCAs. Because of
// this, only the most upstream mismatch in a basin should be repaired before
// the DEM flow network is recalculated. Different basins, however, are
// independent and can each contribute one candidate in the same phase.
//
// Output
// ------
// Returns:
//   - complete_basins   : basins with no mismatched candidate junction
//   - incomplete_basins : basins with a candidate requiring repair
//   - candidates        : at most one candidate record per incomplete basin
//
// Note on distance tolerance
// --------------------------
// match_radius_m is used here as a LOCATION TOLERANCE between the vector
// junction and its topologically correct DEM counterpart (the LCA), not as a
// search radius for arbitrary nearby DEM junctions.
// ---------------------------------------------------------------------------
CandidateJunctionPhaseResult
find_candidate_junctions_phase(
    const lsdtt_vector_to_junction_network& vector_network,
    const lsdtt_junction_network_base&      dem_network,
    const lsdtt_raster_cpu&                 dem,
    int                                     phase,
    const candidate_junction_search_params& params)
{
    CandidateJunctionPhaseResult result;

    const lsdtt_raster_meta& meta = dem.get_meta();
    const double match_radius_m =
        resolved_candidate_match_radius(params, meta);

    const std::vector<dem_candidate_junction_view> dem_views =
        build_dem_candidate_junction_views(dem_network, dem);

    const std::vector<dem_candidate_junction_view> dem_headwaters =
        filter_dem_headwaters(dem_views);

    std::unordered_map<int, const dem_candidate_junction_view*> dem_by_id;
    dem_by_id.reserve(dem_views.size());

    for (const auto& dv : dem_views)
    {
        dem_by_id[dv.junction_id] = &dv;
    }

    // ---------------------------------------------------------------------
    // One vector basin at a time.
    // ---------------------------------------------------------------------
    for (const std::int32_t basin_outlet : vector_network.outlets())
    {
        bool basin_complete = true;

        // For this basin and this phase:
        // vector junction id -> corresponding DEM junction id
        std::unordered_map<int, int> vector_to_dem;
        vector_to_dem.reserve(static_cast<std::size_t>(vector_network.n_junctions()));

        // -----------------------------------------------------------------
        // Step 0: map vector source junctions to DEM headwaters.
        // -----------------------------------------------------------------
        for (const std::int32_t jid : vector_network.topological_order())
        {
            const auto& vj =
                vector_network.junctions()[static_cast<std::size_t>(jid)];

            if (vj.outlet_junction != basin_outlet)
                continue;

            if (!vj.donor_junctions.empty())
                continue;

            double best_distance_m = -1.0;
            const int dem_source_jid =
                find_nearest_dem_junction(
                    lsdtt_xy{vj.x, vj.y},
                    dem_headwaters,
                    best_distance_m
                );

            if (dem_source_jid >= 0)
            {
                vector_to_dem[jid] = dem_source_jid;
            }
        }

        // -----------------------------------------------------------------
        // Step 1: walk vector junctions upstream -> downstream.
        // For each vector junction, the topologically correct DEM counterpart
        // is the LCA of the DEM counterparts of its immediate donor junctions.
        // -----------------------------------------------------------------
        for (const std::int32_t jid : vector_network.topological_order())
        {
            const auto& vj =
                vector_network.junctions()[static_cast<std::size_t>(jid)];

            if (vj.outlet_junction != basin_outlet)
                continue;

            // Source junctions were handled above.
            if (vj.donor_junctions.empty())
                continue;

            const bool is_candidate_type =
                static_cast<int>(vj.donor_junctions.size()) >=
                params.minimum_vector_donor_junctions;

            // Gather DEM counterparts of the immediate donor vector junctions.
            std::vector<int> donor_dem_junctions;
            donor_dem_junctions.reserve(vj.donor_junctions.size());

            bool missing_donor_mapping = false;

            for (const int donor_vj : vj.donor_junctions)
            {
                const auto it = vector_to_dem.find(donor_vj);
                if (it == vector_to_dem.end() || it->second < 0)
                {
                    missing_donor_mapping = true;
                    break;
                }

                donor_dem_junctions.push_back(it->second);
            }

            if (missing_donor_mapping)
            {
                if (is_candidate_type)
                {
                    CandidateJunctionRecord rec;
                    rec.phase = phase;
                    rec.basin_outlet_junction = basin_outlet;

                    rec.vector_junction_id = vj.junction_id;
                    rec.vector_xy = lsdtt_xy{vj.x, vj.y};
                    rec.vector_elevation = vj.elevation;
                    rec.vector_n_donors =
                        static_cast<int>(vj.donor_junctions.size());
                    rec.vector_order = vj.junction_strahler_order;

                    rec.dem_junction_id = -1;
                    rec.dem_n_donors = -1;
                    rec.match_distance_m = -1.0;
                    rec.mismatch_reason = "donor_mapping_missing";

                    result.candidates.push_back(std::move(rec));
                    basin_complete = false;
                    break;
                }

                continue;
            }

            const int dem_lca_jid =
                lca_many_dem_junctions(dem_network, donor_dem_junctions);

            if (dem_lca_jid < 0)
            {
                if (is_candidate_type)
                {
                    CandidateJunctionRecord rec;
                    rec.phase = phase;
                    rec.basin_outlet_junction = basin_outlet;

                    rec.vector_junction_id = vj.junction_id;
                    rec.vector_xy = lsdtt_xy{vj.x, vj.y};
                    rec.vector_elevation = vj.elevation;
                    rec.vector_n_donors =
                        static_cast<int>(vj.donor_junctions.size());
                    rec.vector_order = vj.junction_strahler_order;

                    rec.dem_junction_id = -1;
                    rec.dem_n_donors = -1;
                    rec.match_distance_m = -1.0;
                    rec.mismatch_reason = "no_common_dem_ancestor";

                    result.candidates.push_back(std::move(rec));
                    basin_complete = false;
                    break;
                }

                continue;
            }

            const dem_candidate_junction_view& dem_lca =
                *dem_by_id.at(dem_lca_jid);

            // Degenerate case:
            // if the LCA is one of the donor DEM junctions, then those branches
            // have already merged too early, or multiple vector donors have
            // collapsed onto the same DEM branch.
            bool lca_is_a_donor = false;
            for (const int ddj : donor_dem_junctions)
            {
                if (ddj == dem_lca_jid)
                {
                    lca_is_a_donor = true;
                    break;
                }
            }

            if (lca_is_a_donor)
            {
                if (is_candidate_type)
                {
                    CandidateJunctionRecord rec;
                    rec.phase = phase;
                    rec.basin_outlet_junction = basin_outlet;

                    rec.vector_junction_id = vj.junction_id;
                    rec.vector_xy = lsdtt_xy{vj.x, vj.y};
                    rec.vector_elevation = vj.elevation;
                    rec.vector_n_donors =
                        static_cast<int>(vj.donor_junctions.size());
                    rec.vector_order = vj.junction_strahler_order;

                    rec.dem_junction_id = dem_lca.junction_id;
                    rec.dem_xy = dem_lca.xy;
                    rec.dem_elevation = dem_lca.elevation;
                    rec.dem_n_donors = dem_lca.n_donors;
                    rec.match_distance_m =
                        planar_distance(rec.vector_xy, rec.dem_xy);
                    rec.mismatch_reason = "dem_branches_join_too_early";

                    result.candidates.push_back(std::move(rec));
                    basin_complete = false;
                    break;
                }

                continue;
            }

            const double junction_distance_m =
                planar_distance(lsdtt_xy{vj.x, vj.y}, dem_lca.xy);

            if (is_candidate_type &&
                params.require_matching_donor_count &&
                dem_lca.n_donors != static_cast<int>(vj.donor_junctions.size()))
            {
                CandidateJunctionRecord rec;
                rec.phase = phase;
                rec.basin_outlet_junction = basin_outlet;

                rec.vector_junction_id = vj.junction_id;
                rec.vector_xy = lsdtt_xy{vj.x, vj.y};
                rec.vector_elevation = vj.elevation;
                rec.vector_n_donors =
                    static_cast<int>(vj.donor_junctions.size());
                rec.vector_order = vj.junction_strahler_order;

                rec.dem_junction_id = dem_lca.junction_id;
                rec.dem_xy = dem_lca.xy;
                rec.dem_elevation = dem_lca.elevation;
                rec.dem_n_donors = dem_lca.n_donors;
                rec.match_distance_m = junction_distance_m;
                rec.mismatch_reason = "donor_count_mismatch";

                result.candidates.push_back(std::move(rec));
                basin_complete = false;
                break;
            }

            if (is_candidate_type && junction_distance_m > match_radius_m)
            {
                CandidateJunctionRecord rec;
                rec.phase = phase;
                rec.basin_outlet_junction = basin_outlet;

                rec.vector_junction_id = vj.junction_id;
                rec.vector_xy = lsdtt_xy{vj.x, vj.y};
                rec.vector_elevation = vj.elevation;
                rec.vector_n_donors =
                    static_cast<int>(vj.donor_junctions.size());
                rec.vector_order = vj.junction_strahler_order;

                rec.dem_junction_id = dem_lca.junction_id;
                rec.dem_xy = dem_lca.xy;
                rec.dem_elevation = dem_lca.elevation;
                rec.dem_n_donors = dem_lca.n_donors;
                rec.match_distance_m = junction_distance_m;
                rec.mismatch_reason = "dem_junction_too_far";

                result.candidates.push_back(std::move(rec));
                basin_complete = false;
                break;
            }

            // If we got here, this junction is accepted as correctly matched.
            vector_to_dem[jid] = dem_lca_jid;
        }

        if (basin_complete)
        {
            result.complete_basins.push_back(basin_outlet);
        }
        else
        {
            result.incomplete_basins.push_back(basin_outlet);
        }
    }

    std::sort(result.complete_basins.begin(), result.complete_basins.end());
    std::sort(result.incomplete_basins.begin(), result.incomplete_basins.end());

    return result;
}





// ---------------------------------------------------------------------------
// write_candidate_junction_points
//
// Writes candidate junction diagnostics as point features.
//
// For each candidate record:
//   - one point is written at the vector junction location
//   - if a DEM counterpart exists, one point is also written at the DEM
//     junction location
//
// The attribute is_dem_junction distinguishes the two:
//   0 -> vector junction
//   1 -> DEM junction
//
// This makes it easy to visualise candidate pairs in GIS.
// ---------------------------------------------------------------------------
void write_candidate_junction_points(
    const std::vector<CandidateJunctionRecord>& records,
    const lsdtt_raster_meta&                    meta,
    const std::string&                          filename,
    point_output_format                         fmt,
    bool                                        spatial_index)
{
    if (records.empty()) return;

    if (fmt == point_output_format::csv)
    {
        throw std::invalid_argument(
            "write_candidate_junction_points: CSV output is not supported in this implementation; use flatgeobuf or gpkg."
        );
    }

    lsdtt_vector_schema schema;
    schema.field_names = {
        "phase",
        "basin_outlet",
        "vector_jid",
        "vector_ndonors",
        "vector_order",
        "dem_jid",
        "dem_ndonors",
        "distance_m",
        "is_dem_junction",
        "mismatch_reason"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::string_type
    };

    const vector_output_format vfmt =
        convert_point_to_vector_output_format(fmt);

    lsdtt_data_loader writer;
    writer.begin_vector_write(
        filename,
        "candidate_junctions",
        schema,
        meta.crs_wkt,
        vfmt,
        spatial_index
    );

    for (const auto& r : records)
    {
        // Vector-side point
        writer.write_point(
            r.vector_xy.x,
            r.vector_xy.y,
            {
                static_cast<std::int64_t>(r.phase),
                static_cast<std::int64_t>(r.basin_outlet_junction),
                static_cast<std::int64_t>(r.vector_junction_id),
                static_cast<std::int64_t>(r.vector_n_donors),
                static_cast<std::int64_t>(r.vector_order),
                static_cast<std::int64_t>(r.dem_junction_id),
                static_cast<std::int64_t>(r.dem_n_donors),
                r.match_distance_m,
                std::int64_t{0},
                r.mismatch_reason
            }
        );

        // DEM-side point, if a matched DEM junction exists
        if (r.dem_junction_id >= 0)
        {
            writer.write_point(
                r.dem_xy.x,
                r.dem_xy.y,
                {
                    static_cast<std::int64_t>(r.phase),
                    static_cast<std::int64_t>(r.basin_outlet_junction),
                    static_cast<std::int64_t>(r.vector_junction_id),
                    static_cast<std::int64_t>(r.vector_n_donors),
                    static_cast<std::int64_t>(r.vector_order),
                    static_cast<std::int64_t>(r.dem_junction_id),
                    static_cast<std::int64_t>(r.dem_n_donors),
                    r.match_distance_m,
                    std::int64_t{1},
                    r.mismatch_reason
                }
            );
        }
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
// collect_upslope_channel_pixels
//
// Returns all channel pixels in the upslope contributing area of
// dem_junction_flat that lie within max_distance_m of target_xy.
//
// Why this exists
// ---------------
// In the parallel-swale failure mode the DEM junction is displaced
// DOWNSTREAM of the true divergence point.  The correct repair start — the
// pixel where the channel should branch toward the vector junction — lies
// somewhere upstream of the DEM junction.  Traversing the full upslope
// contributing area and keeping only channel pixels gives the candidate set
// from which repair_junction_channel_guided selects the best divergence point.
//
// The distance filter bounds the search to pixels plausibly within reach of
// the repair target.  Without it, distant headwater channels that share the
// same contributing area would inflate the candidate set and distort scores.
//
// Implementation: flow.for_each_upslope_node traverses all pixels that drain
// through dem_junction_flat; dem_network.is_channel filters to the extracted
// channel skeleton only (not every contributing-area pixel).
// ===========================================================================
static std::vector<std::int64_t> collect_upslope_channel_pixels(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& dem_network,
    std::int64_t                       dem_junction_flat,
    const lsdtt_xy&                    target_xy,
    double                             max_distance_m,
    const lsdtt_raster_meta&           meta)
{
    std::vector<std::int64_t> result;
    const double max_d2 = max_distance_m * max_distance_m;

    flow.for_each_upslope_node(dem_junction_flat, [&](std::int64_t flat)
    {
        if (!dem_network.is_channel(flat)) return;

        const int    row = static_cast<int>(flat / static_cast<std::int64_t>(meta.n_cols));
        const int    col = static_cast<int>(flat % static_cast<std::int64_t>(meta.n_cols));
        const lsdtt_xy xy = cell_center_from_row_col(meta, row, col);

        const double dx = xy.x - target_xy.x;
        const double dy = xy.y - target_xy.y;
        if (dx * dx + dy * dy <= max_d2)
            result.push_back(flat);
    });

    return result;
}

// ===========================================================================
// repair_junction_channel_guided
//
// Three-phase repair for a single mismatched vector–DEM junction pair.
//
// Core problem
// ------------
// In the parallel-swale case the DEM junction is downstream of the true
// divergence point.  Simply routing from the DEM junction to the vector
// junction would carve the wrong reach of channel.  The algorithm must
// instead: (a) find the best divergence point in the upstream channel
// network, (b) find a channelised landing zone near the vector junction,
// and (c) carve a terrain-guided path between them.
//
// Phase 1 — divergence-point selection
//   Collects all channel pixels upstream of dem_junction_flat within
//   2× gap_distance_m of the vector junction (via collect_upslope_channel_pixels).
//   Each candidate pixel is scored:
//     score = w_proximity × 1/(1 + dist/gap) + w_channel × channel_score[px]
//   Proximity dominates (default w_proximity = 0.7) so the divergence point
//   lands close to the repair target; channel score (w_channel = 0.3) gives
//   secondary preference for concave terrain.  When no upstream channel pixels
//   are found the DEM junction itself becomes the divergence point.
//
// Phase 2 — landing zone
//   Scans a square window of radius match_radius_m around the vector junction
//   and collects pixels with channel_score ≥ landing_score_threshold.  The
//   repair does not need to hit the exact vector junction pixel — any
//   channelised pixel in the landing window qualifies as a target.  When none
//   qualify (flat terrain, no curvature signal), the pixel containing the
//   vector junction coordinate is used as the sole landing target.  This
//   correctly handles case B (low-relief parallel swales) where the terrain
//   gives no signal and the vector position is the only ground truth.
//
// Phase 3 — channel-score-guided Dijkstra + monotonic burn
//   Runs multi-target Dijkstra (goal_flats = landing zone pixels) from the
//   divergence point.  Step cost:
//     step_distance × (1 + alpha × (1 − score[px])) + max(0, Δelev) × uphill
//   High channel-score pixels are cheaper, steering the path through existing
//   swale geometry.  The cheapest landing pixel is reached first (Dijkstra
//   guarantee).  The path is then burned monotonically so D8 routes through it.
//
// Returns curvature_bridge on success, unresolved if the LCP fails.
// n_burned is set to the path length in pixels on success.
// ===========================================================================
static GapFixType repair_junction_channel_guided(
    lsdtt_raster_cpu&                  dem,
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& dem_network,
    std::int64_t                       dem_junction_flat,
    const lsdtt_xy&                    vector_junction_xy,
    double                             gap_distance_m,
    double                             match_radius_m,
    const std::vector<float>&          channel_score_data,
    double                             channel_score_alpha,
    double                             w_proximity,
    double                             w_channel,
    double                             landing_score_threshold,
    double                             uphill_mult,
    int&                               n_burned)
{
    const lsdtt_raster_meta&  meta     = dem.get_meta();
    const double              cs       = meta.cell_size;
    const std::int64_t        ncols    = static_cast<std::int64_t>(meta.n_cols);

    auto to_row = [&](std::int64_t f) -> int
        { return static_cast<int>(f / ncols); };
    auto to_col = [&](std::int64_t f) -> int
        { return static_cast<int>(f % ncols); };

    // ── Phase 1: find divergence point ───────────────────────────────────────
    const std::vector<std::int64_t> upstream =
        collect_upslope_channel_pixels(
            flow, dem_network,
            dem_junction_flat,
            vector_junction_xy,
            2.0 * gap_distance_m,
            meta);

    std::int64_t divergence_flat = dem_junction_flat;

    if (!upstream.empty())
    {
        double best_score = -1.0;

        for (std::int64_t flat : upstream)
        {
            const lsdtt_xy xy  = cell_center_from_row_col(
                meta, to_row(flat), to_col(flat));
            const double   dx   = xy.x - vector_junction_xy.x;
            const double   dy   = xy.y - vector_junction_xy.y;
            const double   dist = std::sqrt(dx*dx + dy*dy);

            const double prox = 1.0 / (1.0 + dist / std::max(gap_distance_m, cs));
            const double chan = static_cast<double>(channel_score_data[flat]);
            const double s   = w_proximity * prox + w_channel * chan;

            if (s > best_score) { best_score = s; divergence_flat = flat; }
        }

        const lsdtt_xy div_xy = cell_center_from_row_col(
            meta, to_row(divergence_flat), to_col(divergence_flat));
        std::cout << "    divergence point: ("
                  << div_xy.x << ", " << div_xy.y << ")  score=" << best_score
                  << "  upstream candidates=" << upstream.size() << "\n";
    }
    else
    {
        std::cout << "    no upstream channel pixels within 2× gap; "
                     "using DEM junction as divergence point\n";
    }

    // ── Phase 2: landing zone near vector junction ────────────────────────────
    const int vj_row = static_cast<int>(
        std::floor((meta.y_max() - vector_junction_xy.y) / cs));
    const int vj_col = static_cast<int>(
        std::floor((vector_junction_xy.x - meta.x_min) / cs));

    const int landing_px = static_cast<int>(std::ceil(match_radius_m / cs));

    const std::vector<uint8_t>& dem_mask = dem.get_mask();
    std::vector<std::int64_t>   landing_flats;

    for (int dr = -landing_px; dr <= landing_px; ++dr)
    {
        for (int dc = -landing_px; dc <= landing_px; ++dc)
        {
            const int r = vj_row + dr;
            const int c = vj_col + dc;
            if (r < 0 || r >= meta.n_rows ||
                c < 0 || c >= meta.n_cols) continue;

            const double dist_m =
                std::sqrt(static_cast<double>(dr*dr + dc*dc)) * cs;
            if (dist_m > match_radius_m) continue;

            const std::int64_t flat = static_cast<std::int64_t>(r) * ncols + c;
            if (!dem_mask[flat]) continue;

            const double score =
                static_cast<double>(channel_score_data[flat]);
            if (score >= landing_score_threshold)
                landing_flats.push_back(flat);
        }
    }

    // Fallback: pixel containing the vector junction
    if (landing_flats.empty())
    {
        if (vj_row >= 0 && vj_row < meta.n_rows &&
            vj_col >= 0 && vj_col < meta.n_cols)
        {
            const std::int64_t vj_flat =
                static_cast<std::int64_t>(vj_row) * ncols + vj_col;
            if (dem_mask[vj_flat])
                landing_flats.push_back(vj_flat);
        }
    }

    if (landing_flats.empty())
    {
        std::cerr << "    repair_junction_channel_guided: "
                     "no valid landing pixels near vector junction\n";
        return GapFixType::unresolved;
    }

    std::cout << "    landing zone: " << landing_flats.size()
              << " candidate pixel(s) within " << match_radius_m << " m\n";

    // ── Phase 3: channel-score-guided Dijkstra ────────────────────────────────
    const int div_row = to_row(divergence_flat);
    const int div_col = to_col(divergence_flat);
    const int buf     = landing_px + 10;

    lcp_params lcp;
    lcp.uphill_mult         = uphill_mult;
    lcp.channel_score       = &channel_score_data;
    lcp.channel_score_alpha = channel_score_alpha;
    lcp.goal_flats          = landing_flats;
    lcp.r_min = std::max(0,             std::min(div_row, vj_row) - buf);
    lcp.r_max = std::min(meta.n_rows-1, std::max(div_row, vj_row) + buf);
    lcp.c_min = std::max(0,             std::min(div_col, vj_col) - buf);
    lcp.c_max = std::min(meta.n_cols-1, std::max(div_col, vj_col) + buf);

    lcp_result lcp_res;
    try
    {
        // goal_row / goal_col are ignored because goal_flats is non-empty;
        // pass the vector-junction pixel as a safe (in-bounds) value.
        lcp_res = raster_least_cost_path(
            dem, div_row, div_col,
            vj_row, vj_col,
            lcp);
    }
    catch (const std::exception& e)
    {
        std::cerr << "    repair_junction_channel_guided: LCP failed: "
                  << e.what() << "\n";
        return GapFixType::unresolved;
    }

    if (!lcp_res.found)
    {
        std::cout << "    repair_junction_channel_guided: "
                     "no path found in search window\n";
        return GapFixType::unresolved;
    }

    burn_monotonic_path(dem, lcp_res.path);
    n_burned = static_cast<int>(lcp_res.path.size());

    std::cout << "    channel-guided path burned: " << n_burned
              << " px, climb=" << lcp_res.max_climb << " m\n";

    return GapFixType::curvature_bridge;
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

    // ── Channel candidate scores (computed once on the unmodified DEM) ──────
    //
    // Channel score = multiscale plan-curvature persistence (Hessian + depth).
    // Range [0, 1]: 1 = strong valley / channel, 0 = ridge or planar terrain.
    //
    // Two roles in the repair algorithm:
    //   (a) Phase 1 divergence scoring — channel score contributes secondary
    //       preference for concave terrain when ranking upstream channel pixels
    //       as divergence-point candidates.
    //   (b) Phase 3 Dijkstra weighting — high-score pixels are cheaper to
    //       traverse, steering the repair path through existing swale geometry
    //       rather than cutting across hillslopes.
    //
    // Curvature (not flow accumulation) is used because curvature is a terrain
    // geometry property that remains meaningful even at road crossings and
    // embankments where flow accumulation is artificially disrupted by the
    // obstruction.
    //
    // Computed once before the repair loop.  Small monotonic burns lower only
    // a few pixels and do not significantly alter the multi-scale curvature
    // structure; recomputing after every burn would be expensive and would not
    // materially change repair decisions for subsequent junctions.
    bool               have_channel_scores = false;
    std::vector<float> channel_score_data;

    {
        const double cs = meta.cell_size;
        channel_candidate_request cc_req;
        cc_req.channel_lambda_min_metres =
            (params.channel_lambda_min_m > 0)
            ? params.channel_lambda_min_m : cs * 3.0;
        cc_req.channel_lambda_max_metres =
            (params.channel_lambda_max_m > 0)
            ? params.channel_lambda_max_m : cs * 10.0;

        if (cc_req.channel_lambda_max_metres <= cc_req.channel_lambda_min_metres)
            cc_req.channel_lambda_max_metres = cc_req.channel_lambda_min_metres * 3.0;

        cc_req.channel_n_scales          = params.channel_n_scales;
        cc_req.weight_persistent         = 0.7;
        cc_req.weight_depth              = 0.3;
        cc_req.weight_slope              = 0.0;
        cc_req.use_hessian_filter        = true;
        cc_req.positive_curvature_is_good = true;

        std::cout << "repair_stream_junctions: computing channel candidates "
                  << "(lambda " << cc_req.channel_lambda_min_metres
                  << "–"        << cc_req.channel_lambda_max_metres
                  << " m, " << cc_req.channel_n_scales << " scales)...\n";
        try
        {
            channel_candidate_result cc = find_channel_candidates(dem, cc_req);
            channel_score_data  = cc.score.get_data();
            have_channel_scores = true;
            std::cout << "  channel candidates computed.\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "  Warning: channel candidate computation failed: "
                      << e.what() << "\n"
                      << "  All repairs will use terrain-only Dijkstra.\n";
        }
    }

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

    // ── Nearest-junction distance helper ─────────────────────────────────────
    // Scans the current dem_network for the junction closest to target.
    // Always uses the most recent dem_network (captured by reference).
    auto nearest_junction_dist_to = [&](const lsdtt_xy& target) -> double
    {
        double best_d2 = std::numeric_limits<double>::max();
        const std::int32_t nj = dem_network->n_junctions();
        for (std::int32_t j = 0; j < nj; ++j)
        {
            if (!dem_network->is_valid_junction(j)) continue;
            const lsdtt_xy jxy = cell_center_from_row_col(
                meta,
                dem_network->row_of_junction(j),
                dem_network->col_of_junction(j));
            const double d2 = (jxy.x - target.x) * (jxy.x - target.x)
                            + (jxy.y - target.y) * (jxy.y - target.y);
            if (d2 < best_d2) best_d2 = d2;
        }
        return (best_d2 < std::numeric_limits<double>::max())
               ? std::sqrt(best_d2)
               : std::numeric_limits<double>::max();
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

            GapFixType fix_type = GapFixType::unresolved;
            int        n_burned = 0;
            bool       repair_ok = false;

            // ── Attempt 1: three-phase channel-guided repair ──────────────────
            if (have_channel_scores)
            {
                const std::int64_t dem_junc_flat =
                    dem_network->flat_of_junction(dem_junc);

                fix_type = repair_junction_channel_guided(
                    dem, *flow, *dem_network,
                    dem_junc_flat,
                    junc_node.xy,
                    dist,
                    params.min_gap_m,
                    channel_score_data,
                    params.channel_score_alpha,
                    params.w_proximity,
                    params.w_channel,
                    params.landing_score_threshold,
                    params.uphill_mult,
                    n_burned);

                if (fix_type == GapFixType::curvature_bridge)
                {
                    // Re-extract and check whether the junction moved
                    flow        = make_flow_base(dem.get_data_array(),
                                                 dem.get_mask_array(), meta, fi_opts);
                    dem_network = make_junction_network_base(*flow, source_flats, jn_opts);
                    dem_hw      = build_headwaters();
                    cache.clear();

                    const double new_gap = nearest_junction_dist_to(junc_node.xy);
                    if (new_gap <= params.min_gap_m)
                    {
                        repair_ok = true;
                        std::cout << "    → channel-guided: success"
                                  << "  new gap=" << new_gap << " m\n";
                    }
                    else
                    {
                        std::cout << "    → channel-guided: burned but gap still "
                                  << new_gap << " m; trying terrain Dijkstra fallback\n";
                    }
                }
                else
                {
                    std::cout << "    → channel-guided: " << gap_fix_type_str(fix_type)
                              << "; trying terrain Dijkstra fallback\n";
                }
            }

            // ── Attempt 2: terrain-only Dijkstra fallback ─────────────────────
            if (!repair_ok)
            {
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

                fix_type  = pv[0].fix_type;
                n_burned += pv[0].n_burned;

                std::cout << "    → terrain fallback: " << gap_fix_type_str(fix_type)
                          << " (" << pv[0].n_burned << " px)"
                          << "  re-extracting network...\n";

                flow        = make_flow_base(dem.get_data_array(),
                                             dem.get_mask_array(), meta, fi_opts);
                dem_network = make_junction_network_base(*flow, source_flats, jn_opts);
                dem_hw      = build_headwaters();
                cache.clear();
            }

            JunctionRepairRecord rec;
            rec.repair_index  = repair_index++;
            rec.network_index = ni;
            rec.vector_xy     = junc_node.xy;
            rec.dem_xy        = dem_jxy;
            rec.distance_m    = dist;
            rec.fix_type      = fix_type;
            rec.n_burned      = n_burned;
            records.push_back(rec);
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
