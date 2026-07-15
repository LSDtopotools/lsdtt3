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
#include <sstream>
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
#include "../core/lsdtt_vector_utils.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "../utils/lsdtt_linear_algebra.hpp"
#include "lsdtt_least_cost_path.hpp"
#include "lsdtt_proximity_cpu.hpp"
#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_feature_finder_cpu.hpp"
#include "lsdtt_raster_skeleton_graph_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"

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

    // endpoint_map: key -> number of features that have this as an endpoint
    std::map<std::string, int>                    endpoint_count;
    std::map<std::string, lsdtt_xy>               endpoint_xy;
    // For each endpoint store the full inward vertex sequence - all vertices
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

        // front endpoint -> inward sequence is points[1 .. n-1]
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

        // back endpoint -> inward sequence is points[n-2 .. 0] (reversed)
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
    // Outlet identification - segment gradient criterion only
    //
    // A tip is at the upstream (source) end of its segment when the first
    // valid inward vertex is lower than the tip (water flows away from the
    // tip down the segment).  A tip is at the downstream (outlet) end when
    // the first valid inward vertex is higher (water flows toward the tip).
    //
    // If no tip satisfies the downstream-end condition the real outlet is not
    // among the degree-1 tips - it is either at a junction (degree > 1, not
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

            const double diff = t.adjacent_elev - t.elevation; // > 0 -> tip is downstream
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

    // ---- 5. Build adjacency: junction key -> list of segment indices ---------
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
// repair_junction_gaps  - internal helpers
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
                      << ": junction outside raster - skipping\n";
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

        //  Least-cost path: D_j -> V_j neighbourhood 
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
            burn_monotonic_path(dem, lcp_res.path, 0.0);
            pair.fix_type = GapFixType::channel_guided_bridge;
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
                      << " m > " << max_climb_m << " m - vector burn fallback\n";
        else
            std::cout << "  pair " << pair.pair_id
                      << ": no terrain path - vector burn fallback\n";

        //  Fallback: burn under vector segments within the window 
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

    auto precise_string = [](double value) -> std::string
    {
        std::ostringstream ss;
        ss.precision(17);
        ss << value;
        return ss.str();
    };

    if (fmt == point_output_format::csv)
    {
        if (meta.crs_wkt.empty())
            throw std::runtime_error(
                "write_candidate_junction_points: CSV output requires a valid CRS");

        const lsdtt_coordinate_transform transform =
            lsdtt_coordinate_transform::from_epsg(meta.crs_wkt, 4326);

        std::vector<std::string> headers = {
            "latitude", "longitude", "phase", "basin_outlet",
            "vector_jid", "vector_ndonors", "vector_order",
            "dem_jid", "dem_ndonors", "distance_m",
            "is_dem_junction", "mismatch_reason"
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(records.size() * 2);

        for (const auto& r : records)
        {
            const lsdtt_lonlat llv = xy_to_lonlat(r.vector_xy.x, r.vector_xy.y, transform);
            rows.push_back({
                precise_string(llv.latitude), precise_string(llv.longitude),
                std::to_string(r.phase), std::to_string(r.basin_outlet_junction),
                std::to_string(r.vector_junction_id), std::to_string(r.vector_n_donors),
                std::to_string(r.vector_order), std::to_string(r.dem_junction_id),
                std::to_string(r.dem_n_donors), precise_string(r.match_distance_m),
                "0", r.mismatch_reason
            });

            if (r.dem_junction_id >= 0)
            {
                const lsdtt_lonlat lld = xy_to_lonlat(r.dem_xy.x, r.dem_xy.y, transform);
                rows.push_back({
                    precise_string(lld.latitude), precise_string(lld.longitude),
                    std::to_string(r.phase), std::to_string(r.basin_outlet_junction),
                    std::to_string(r.vector_junction_id), std::to_string(r.vector_n_donors),
                    std::to_string(r.vector_order), std::to_string(r.dem_junction_id),
                    std::to_string(r.dem_n_donors), precise_string(r.match_distance_m),
                    "1", r.mismatch_reason
                });
            }
        }

        lsdtt_data_loader loader;
        loader.write_csv(filename, headers, rows);
        return;
    }

    lsdtt_vector_schema schema;
    schema.field_names = {
        "phase", "basin_outlet", "vector_jid", "vector_ndonors", "vector_order",
        "dem_jid", "dem_ndonors", "distance_m", "is_dem_junction", "mismatch_reason"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type, lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type, lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type, lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type, lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type, lsdtt_vector_field_type::string_type
    };

    const vector_output_format vfmt = convert_point_to_vector_output_format(fmt);

    lsdtt_data_loader writer;
    writer.begin_vector_write(filename, "candidate_junctions", schema,
                              meta.crs_wkt, vfmt, spatial_index);

    for (const auto& r : records)
    {
        writer.write_point(r.vector_xy.x, r.vector_xy.y, {
            static_cast<std::int64_t>(r.phase),
            static_cast<std::int64_t>(r.basin_outlet_junction),
            static_cast<std::int64_t>(r.vector_junction_id),
            static_cast<std::int64_t>(r.vector_n_donors),
            static_cast<std::int64_t>(r.vector_order),
            static_cast<std::int64_t>(r.dem_junction_id),
            static_cast<std::int64_t>(r.dem_n_donors),
            r.match_distance_m, std::int64_t{0}, r.mismatch_reason
        });

        if (r.dem_junction_id >= 0)
        {
            writer.write_point(r.dem_xy.x, r.dem_xy.y, {
                static_cast<std::int64_t>(r.phase),
                static_cast<std::int64_t>(r.basin_outlet_junction),
                static_cast<std::int64_t>(r.vector_junction_id),
                static_cast<std::int64_t>(r.vector_n_donors),
                static_cast<std::int64_t>(r.vector_order),
                static_cast<std::int64_t>(r.dem_junction_id),
                static_cast<std::int64_t>(r.dem_n_donors),
                r.match_distance_m, std::int64_t{1}, r.mismatch_reason
            });
        }
    }

    writer.finish_vector_write();
}


// ---------------------------------------------------------------------------
// compute_reach_divergence
// ---------------------------------------------------------------------------
std::vector<ReachDivergenceRecord>
compute_reach_divergence(
    const lsdtt_raster_cpu&                 dem,
    const lsdtt_flow_base&                  flow,
    const lsdtt_junction_network_base&      dem_network,
    const lsdtt_vector_to_junction_network& vector_network,
    const reach_divergence_params&          params,
    const std::vector<float>*               plan_curvature,
    const std::vector<float>*               channel_score)
{
    std::vector<ReachDivergenceRecord> out;

    const lsdtt_raster_meta&  meta  = dem.get_meta();
    const std::vector<float>& zdata = dem.get_data();
    const int    n_cols = meta.n_cols;
    const double match_radius =
        (params.match_radius_m > 0.0) ? params.match_radius_m
                                      : 5.0 * meta.cell_size;
    const double sep_thresh = params.separation_threshold_m;

    const std::vector<dem_candidate_junction_view> dem_views =
        build_dem_candidate_junction_views(dem_network, dem);

    auto flat_to_xy = [&](std::int64_t flat) -> lsdtt_xy
    {
        const int row = static_cast<int>(flat / n_cols);
        const int col = static_cast<int>(flat % n_cols);
        return cell_center_from_row_col(meta, row, col);
    };
    auto z_of = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(zdata[static_cast<std::size_t>(flat)]);
    };
    auto sample = [&](const std::vector<float>* arr, std::int64_t flat) -> double
    {
        if (arr == nullptr) return 0.0;
        if (flat < 0 || flat >= static_cast<std::int64_t>(arr->size())) return 0.0;
        const float v = (*arr)[static_cast<std::size_t>(flat)];
        return std::isfinite(v) ? static_cast<double>(v) : 0.0;
    };
    auto polyline_len = [](const std::vector<lsdtt_xy>& xy) -> double
    {
        double L = 0.0;
        for (std::size_t i = 1; i < xy.size(); ++i)
            L += std::hypot(xy[i].x - xy[i - 1].x, xy[i].y - xy[i - 1].y);
        return L;
    };
    // Minimum distance from p to the polyline xy.
    auto sep_to_polyline = [](const lsdtt_xy& p,
                              const std::vector<lsdtt_xy>& xy) -> double
    {
        if (xy.empty()) return 0.0;
        if (xy.size() == 1) return std::hypot(p.x - xy[0].x, p.y - xy[0].y);
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t i = 1; i < xy.size(); ++i)
        {
            const double d = point_to_segment_distance(
                p.x, p.y, xy[i - 1].x, xy[i - 1].y, xy[i].x, xy[i].y);
            if (d < best) best = d;
        }
        return best;
    };

    const auto& vjns = vector_network.junctions();
    const auto& veds = vector_network.edges();
    const std::int32_t n_vjn = vector_network.n_junctions();

    for (const auto& E : veds)
    {
        if (E.upstream_junction < 0 || E.receiver_junction < 0)   continue;
        if (E.upstream_junction >= n_vjn || E.receiver_junction >= n_vjn) continue;
        if (E.xy.size() < 2) continue;

        const auto& vj_up = vjns[static_cast<std::size_t>(E.upstream_junction)];
        const auto& vj_dn = vjns[static_cast<std::size_t>(E.receiver_junction)];

        ReachDivergenceRecord r;
        r.vector_edge_id  = E.edge_id;
        r.vector_jid_up   = E.upstream_junction;
        r.vector_jid_dn   = E.receiver_junction;
        r.vector_reach_xy = E.xy;
        r.vector_length_m =
            (E.length_metres > 0.0) ? E.length_metres : polyline_len(E.xy);

        const double vdrop = vj_up.elevation - vj_dn.elevation;
        if (r.vector_length_m > 0.0)
            r.vector_gradient = vdrop / r.vector_length_m;

        // Match both vector endpoints to DEM junctions.
        const int dj_up = find_nearest_dem_junction(
            lsdtt_xy{vj_up.x, vj_up.y}, dem_views, r.match_dist_up_m);
        const int dj_dn = find_nearest_dem_junction(
            lsdtt_xy{vj_dn.x, vj_dn.y}, dem_views, r.match_dist_dn_m);
        r.dem_jid_up = dj_up;
        r.dem_jid_dn = dj_dn;

        const bool up_matched = (dj_up >= 0 && r.match_dist_up_m <= match_radius);
        const bool dn_matched = (dj_dn >= 0 && r.match_dist_dn_m <= match_radius);
        r.endpoint_unmatched = !(up_matched && dn_matched);

        // We still need a DEM junction to start the walk. If even the nearest
        // upstream DEM junction does not exist (no DEM network here at all),
        // there is nothing to draw or repair.
        if (dj_up < 0)
        {
            r.status = ReachStatus::unmatched_endpoint;
            out.push_back(std::move(r));
            continue;
        }

        // Walk the DEM channel downstream from the (nearest) upstream DEM
        // junction. We do this even when an endpoint is beyond match_radius:
        // the walk follows where the DEM channel actually goes (to the real
        // downstream junction, an outlet, or off the DEM edge), and that
        // divergent/off-edge path is exactly what a topology repair needs.
        const std::int64_t up_flat = dem_network.flat_of_junction(dj_up);
        const std::int64_t dn_flat =
            (dj_dn >= 0) ? dem_network.flat_of_junction(dj_dn) : -1;

        std::vector<std::int64_t> walk;
        walk.push_back(up_flat);
        std::int64_t cur = up_flat;
        bool reached = (cur == dn_flat);
        std::int64_t steps = 0;
        while (!reached && steps < params.max_walk_pixels)
        {
            if (!flow.is_flat_in_range(cur)) break;
            if (flow.is_outlet(cur))         break;
            const std::int64_t nxt = flow.receiver_of(cur);
            if (nxt < 0 || nxt == cur)       break;
            cur = nxt;
            walk.push_back(cur);
            ++steps;
            if (cur == dn_flat) { reached = true; break; }
        }

        // Metrics over the walked DEM polyline.
        double dem_len = 0.0, sum_sep = 0.0, max_sep = 0.0, div_len = 0.0;
        double sum_curv = 0.0, sum_score = 0.0;
        r.dem_reach_xy.reserve(walk.size());
        lsdtt_xy prev_xy{0.0, 0.0};
        for (std::size_t i = 0; i < walk.size(); ++i)
        {
            const lsdtt_xy p = flat_to_xy(walk[i]);
            r.dem_reach_xy.push_back(p);

            double seg = 0.0;
            if (i > 0)
            {
                seg = std::hypot(p.x - prev_xy.x, p.y - prev_xy.y);
                dem_len += seg;
            }

            const double sep = sep_to_polyline(p, E.xy);
            sum_sep += sep;
            if (sep > max_sep) max_sep = sep;
            if (sep > sep_thresh)
            {
                div_len += seg;
                if (!r.has_divergence_point)
                {
                    r.has_divergence_point  = true;
                    r.divergence_point_xy   = p;
                    r.divergence_point_flat = walk[i];
                }
            }

            sum_curv  += sample(plan_curvature, walk[i]);
            sum_score += sample(channel_score, walk[i]);
            prev_xy = p;
        }

        const double n = static_cast<double>(std::max<std::size_t>(walk.size(), 1));
        r.dem_path_length_m  = dem_len;
        r.detour_ratio       = (r.vector_length_m > 0.0)
                                   ? dem_len / r.vector_length_m : 0.0;
        r.max_separation_m   = max_sep;
        r.mean_separation_m  = sum_sep / n;
        r.divergence_length_m = div_len;
        r.mean_plan_curvature = sum_curv / n;
        r.mean_channel_score  = sum_score / n;

        const double ddrop = z_of(walk.front()) - z_of(walk.back());
        if (dem_len > 0.0) r.dem_gradient = ddrop / dem_len;
        if (std::abs(r.vector_gradient) > 1.0e-12)
            r.gradient_ratio = r.dem_gradient / r.vector_gradient;

        // Dimensionless severity for ranking/symbolising.
        const double sep_term = (sep_thresh > 0.0)
            ? r.max_separation_m / sep_thresh : 0.0;
        r.severity = sep_term +
            params.severity_detour_weight * std::max(0.0, r.detour_ratio - 1.0);

        // Classify - geometric only. (Geometry is retained for every reach
        // above, including the unmatched ones, so they register in QGIS.)
        if (r.endpoint_unmatched)
        {
            // A vector endpoint has no DEM junction within match_radius: a
            // junction gap, the topological (LCA) repair's territory rather
            // than the divergence repair's.
            r.status = ReachStatus::unmatched_endpoint;
        }
        else if (!reached)
        {
            r.status = ReachStatus::topology_break;
        }
        else if (r.divergence_length_m >= params.diverged_min_divergence_m ||
                 r.detour_ratio        >= params.diverged_detour_ratio)
        {
            r.status = ReachStatus::diverged;
        }
        else
        {
            r.status = ReachStatus::ok;
        }

        out.push_back(std::move(r));
    }

    return out;
}


// ---------------------------------------------------------------------------
// write_reach_divergence
// ---------------------------------------------------------------------------
void write_reach_divergence(
    const std::vector<ReachDivergenceRecord>& records,
    const lsdtt_raster_meta&                  meta,
    const std::string&                        prefix,
    point_output_format                       fmt,
    bool                                      spatial_index)
{
    if (records.empty()) return;

    lsdtt_vector_schema schema;
    schema.field_names = {
        "edge_id", "v_jid_up", "v_jid_dn", "dem_jid_up", "dem_jid_dn",
        "match_up_m", "match_dn_m", "vec_len_m", "dem_len_m", "detour_ratio",
        "vec_grad", "dem_grad", "grad_ratio", "max_sep_m", "mean_sep_m",
        "diverg_len_m", "severity", "mean_plan_curv", "mean_chan_score",
        "status", "endpt_unm"};
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,  lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,  lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,  lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type, lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type, lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type, lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type, lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type, lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type, lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type, lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::int64_type};

    auto fields = [](const ReachDivergenceRecord& r)
        -> std::vector<lsdtt_vector_field_value>
    {
        return {
            static_cast<std::int64_t>(r.vector_edge_id),
            static_cast<std::int64_t>(r.vector_jid_up),
            static_cast<std::int64_t>(r.vector_jid_dn),
            static_cast<std::int64_t>(r.dem_jid_up),
            static_cast<std::int64_t>(r.dem_jid_dn),
            r.match_dist_up_m, r.match_dist_dn_m,
            r.vector_length_m, r.dem_path_length_m, r.detour_ratio,
            r.vector_gradient, r.dem_gradient, r.gradient_ratio,
            r.max_separation_m, r.mean_separation_m, r.divergence_length_m,
            r.severity, r.mean_plan_curvature, r.mean_channel_score,
            std::string(reach_status_str(r.status)),
            static_cast<std::int64_t>(r.endpoint_unmatched ? 1 : 0)};
    };

    const vector_output_format vfmt =
        (fmt == point_output_format::csv)
            ? vector_output_format::flatgeobuf
            : convert_point_to_vector_output_format(fmt);

    // A bare path without an extension makes GDAL create a directory datasource;
    // append the format extension so a single file is written.
    const std::string ext = output_extension(vfmt);

    // DEM-walked reach polylines.
    {
        lsdtt_data_loader writer;
        writer.begin_vector_write(prefix + "_dem_reach" + ext, "dem_reach",
                                  schema, meta.crs_wkt, vfmt, spatial_index);
        for (const auto& r : records)
            if (r.dem_reach_xy.size() >= 2)
                writer.write_line(r.dem_reach_xy, fields(r));
        writer.finish_vector_write();
    }

    // Matched vector edge polylines (same attributes for side-by-side QC).
    {
        lsdtt_data_loader writer;
        writer.begin_vector_write(prefix + "_vector_reach" + ext, "vector_reach",
                                  schema, meta.crs_wkt, vfmt, spatial_index);
        for (const auto& r : records)
            if (r.vector_reach_xy.size() >= 2)
                writer.write_line(r.vector_reach_xy, fields(r));
        writer.finish_vector_write();
    }

    // CSV table (one row per reach).
    {
        std::vector<std::string> headers = schema.field_names;
        std::vector<std::vector<std::string>> rows;
        rows.reserve(records.size());
        auto num = [](double v) {
            std::ostringstream os; os << std::setprecision(8) << v; return os.str();
        };
        for (const auto& r : records)
        {
            rows.push_back({
                std::to_string(r.vector_edge_id), std::to_string(r.vector_jid_up),
                std::to_string(r.vector_jid_dn), std::to_string(r.dem_jid_up),
                std::to_string(r.dem_jid_dn), num(r.match_dist_up_m),
                num(r.match_dist_dn_m), num(r.vector_length_m),
                num(r.dem_path_length_m), num(r.detour_ratio),
                num(r.vector_gradient), num(r.dem_gradient), num(r.gradient_ratio),
                num(r.max_separation_m), num(r.mean_separation_m),
                num(r.divergence_length_m), num(r.severity),
                num(r.mean_plan_curvature),
                num(r.mean_channel_score), std::string(reach_status_str(r.status)),
                std::to_string(r.endpoint_unmatched ? 1 : 0)});
        }
        lsdtt_data_loader loader;
        loader.write_csv(prefix + ".csv", headers, rows);
    }
}


// ---------------------------------------------------------------------------
// resolve_stream_burn_mode
// ---------------------------------------------------------------------------
stream_burn_mode resolve_stream_burn_mode(const std::string& s)
{
    if (s == "burn_in_place")        return stream_burn_mode::burn_in_place;
    if (s == "snap_to_network")      return stream_burn_mode::snap_to_network;
    if (s == "snap_to_accumulation") return stream_burn_mode::snap_to_accumulation;
    std::cerr << "Warning: unknown stream_burn_mode '" << s
              << "'; defaulting to burn_in_place.\n";
    return stream_burn_mode::burn_in_place;
}

// ---------------------------------------------------------------------------
// burn_streams_in_place - internal helpers
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

        // outlet_at_end = true  -> polyline end (high along_axis) is downstream
        // outlet_at_end = false -> polyline start (along_axis ~= 0) is downstream
        const bool outlet_at_end = (sum_end < sum_start);
        const double max_along   = raw_nodes.back().along_axis;

        // -------------------------------------------------------------------
        // Build profile nodes with flow_distance (0 = outlet, increasing up)
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

// Segment [a->b]: a is inside the bbox, b is outside.
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

        // Kahn's topological sort - collect junctions in upstream-first order
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
// DOWNSTREAM of the true divergence point.  The correct repair start - the
// pixel where the channel should branch toward the vector junction - lies
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
// Anonymous namespace - internal repair helpers
// ===========================================================================
namespace
{

struct repair_attempt_result
{
    GapFixType fix_type        = GapFixType::unresolved;
    int        n_burned        = 0;
    double     max_climb_m     = -1.0;
    double     total_cost      = -1.0;
    int        n_landing_pixels = 0;
    lsdtt_xy   divergence_xy;   // breach start (on lower channel)
    lsdtt_xy   target_xy;       // path endpoint (on higher channel)
    std::vector<std::int64_t> burned_path;  // flat indices start->end
};

double nearest_junction_distance_to_xy(
    const lsdtt_junction_network_base& dem_network,
    const lsdtt_raster_meta&           meta,
    const lsdtt_xy&                    target_xy,
    int*                               best_junction = nullptr)
{
    double best_d2 = std::numeric_limits<double>::max();
    int best_j = -1;
    const std::int32_t nj = dem_network.n_junctions();
    for (std::int32_t j = 0; j < nj; ++j)
    {
        if (!dem_network.is_valid_junction(j)) continue;
        const lsdtt_xy jxy = cell_center_from_row_col(
            meta, dem_network.row_of_junction(j), dem_network.col_of_junction(j));
        const double dx = jxy.x - target_xy.x;
        const double dy = jxy.y - target_xy.y;
        const double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; best_j = j; }
    }
    if (best_junction) *best_junction = best_j;
    return (best_j >= 0) ? std::sqrt(best_d2)
                         : std::numeric_limits<double>::infinity();
}

bool xy_to_valid_row_col(
    const lsdtt_raster_cpu& dem,
    const lsdtt_xy&         xy,
    int&                    row,
    int&                    col)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    col = static_cast<int>(std::floor((xy.x - meta.x_min) / meta.cell_size));
    row = static_cast<int>(std::floor((meta.y_max() - xy.y) / meta.cell_size));
    if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols)
        return false;
    const std::int64_t flat = static_cast<std::int64_t>(row) * meta.n_cols + col;
    return dem.get_mask()[static_cast<std::size_t>(flat)] != 0;
}

} // anonymous namespace

// ===========================================================================
// repair_junction_parallel_channels
//
// Both repair endpoints are on DEM-extracted channels.  The vector junction
// is used only to identify WHICH two DEM branches need connecting.
//
// The DEM junction (candidate.dem_junction_id) is the confluence where the
// two mis-matched branches currently meet.  Its donor junctions in the DEM
// network are the two upstream branches:
//
//   Branch B (correct)  - the donor whose junction is closest to the vector
//                         junction location.  D8 is already routing this
//                         branch roughly correctly.
//   Branch A (wrong)    - the donor whose junction is farthest from the
//                         vector junction.  D8 is routing this branch down
//                         the wrong valley.
//
// Algorithm:
//   1. Classify the donors of the DEM junction into Branch A (far) and
//      Branch B (close) using distance from the vector junction.
//   2. Collect DEM channel pixels for each branch within 2x gap.
//   3. Find the closest pair of pixels (one from each branch) - that pair
//      sits at the natural col between the two parallel valleys.
//   4. Determine which pixel is lower.  The lower one is the breach start;
//      the higher branch pixels are the goal set.
//   5. Additionally: from the lower breach pixel, search upstream for the
//      farthest channel pixel that still faces the higher branch (via
//      channelised bank score), and use that as the actual breach start.
//      This ensures the burn reaches as far upstream as needed.
//   6. Run LCP from breach start to higher branch (goal_flats); burn.
// ===========================================================================
static repair_attempt_result repair_junction_parallel_channels(
    lsdtt_raster_cpu&                           dem,
    const lsdtt_flow_base&                      flow,
    const lsdtt_junction_network_base&          dem_network,
    std::int64_t                                dem_junction_flat,
    const lsdtt_xy&                             vector_junction_xy,
    int                                         /* vector_junction_id */,
    const lsdtt_vector_to_junction_network&     /* vector_network */,
    double                                      gap_distance_m,
    double                                      match_radius_m,
    const std::vector<float>&                   channel_score_data,
    const junction_repair_params&               params,
    const lsdtt_xy*                             breach_anchor_override = nullptr)
{
    repair_attempt_result out;
    out.fix_type = GapFixType::unresolved;

    // Branch classification (which donor is the wrong vs correct valley) always
    // uses the true vector junction. The BREACH search can be re-anchored: for a
    // diverged-but-rejoining reach we anchor it at the divergence point (where
    // the channel leaves the correct valley) so flow is pulled back THERE, not
    // at the downstream confluence. Null override = anchor at the vector junction.
    const lsdtt_xy breach_anchor =
        breach_anchor_override ? *breach_anchor_override : vector_junction_xy;

    const lsdtt_raster_meta& meta  = dem.get_meta();
    const double  cs    = meta.cell_size;
    const int     nr    = meta.n_rows;
    const int     nc    = meta.n_cols;
    const std::int64_t ncols64 = static_cast<std::int64_t>(nc);

    auto to_row  = [&](std::int64_t f) { return static_cast<int>(f / ncols64); };
    auto to_col  = [&](std::int64_t f) { return static_cast<int>(f % ncols64); };
    auto to_flat = [&](int r, int c)   { return static_cast<std::int64_t>(r) * ncols64 + c; };

    static constexpr int DR[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    static constexpr int DC[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    const auto& elev_data = dem.get_data();
    const auto& mask_data = dem.get_mask();

    //  Step 1: classify donors into Branch A (far) and Branch B (close) 
    const std::int32_t dem_junc32 =
        static_cast<std::int32_t>(dem_junction_flat / ncols64 * ncols64 == dem_junction_flat
            ? dem_junction_flat / ncols64
            : dem_junction_flat / ncols64); // unused; get via junction API

    // Find which junction ID corresponds to dem_junction_flat
    std::int32_t dem_jid = -1;
    {
        const std::int32_t nj = dem_network.n_junctions();
        for (std::int32_t j = 0; j < nj; ++j)
        {
            if (dem_network.flat_of_junction(j) == dem_junction_flat)
            { dem_jid = j; break; }
        }
    }
    if (dem_jid < 0)
    {
        std::cerr << "    parallel: DEM junction not found for flat\n";
        return out;
    }

    const std::int32_t n_donors = dem_network.donor_junction_size(dem_jid);
    if (n_donors < 2)
    {
        std::cerr << "    parallel: DEM junction has fewer than 2 donors ("
                  << n_donors << "); cannot identify two branches\n";
        return out;
    }

    std::int32_t branch_a_jid = -1, branch_b_jid = -1;
    double best_close_d2 = std::numeric_limits<double>::max();
    double best_far_d2   = 0.0;

    for (std::int32_t donor : dem_network.donor_junctions(dem_jid))
    {
        const lsdtt_xy dxy = cell_center_from_row_col(
            meta,
            dem_network.row_of_junction(donor),
            dem_network.col_of_junction(donor));
        const double dx = dxy.x - vector_junction_xy.x;
        const double dy = dxy.y - vector_junction_xy.y;
        const double d2 = dx*dx + dy*dy;
        if (d2 < best_close_d2) { best_close_d2 = d2; branch_b_jid = donor; }
        if (d2 > best_far_d2)   { best_far_d2   = d2; branch_a_jid = donor; }
    }

    if (branch_a_jid == branch_b_jid)
    {
        // All donors at same distance - pick the two with highest basin area
        // (fallback: just pick first two)
        auto donors = dem_network.donor_junctions(dem_jid);
        branch_a_jid = donors[0];
        branch_b_jid = donors[1];
    }

    const lsdtt_xy branch_a_xy = cell_center_from_row_col(
        meta,
        dem_network.row_of_junction(branch_a_jid),
        dem_network.col_of_junction(branch_a_jid));
    const lsdtt_xy branch_b_xy = cell_center_from_row_col(
        meta,
        dem_network.row_of_junction(branch_b_jid),
        dem_network.col_of_junction(branch_b_jid));

    std::cout << "    parallel: Branch_A(wrong) jid=" << branch_a_jid
              << " at (" << branch_a_xy.x << "," << branch_a_xy.y << ")"
              << "  Branch_B(correct) jid=" << branch_b_jid
              << " at (" << branch_b_xy.x << "," << branch_b_xy.y << ")\n";

    //  Step 2: collect and classify channel pixels for each branch 
    //
    // Branch A goes down the WRONG valley - its pixels are far from the
    // vector junction, so collecting from branch_a_flat upslope centred on
    // the vector junction gives almost nothing.
    //
    // Correct approach: collect ALL upslope channel pixels from the DEM
    // junction (which captures both branches), then classify each pixel to
    // Branch A or B by checking which donor's upslope set it belongs to.
    const std::int64_t branch_a_flat = dem_network.flat_of_junction(branch_a_jid);
    const std::int64_t branch_b_flat = dem_network.flat_of_junction(branch_b_jid);

    // Collect the channel pixels of each DEM segment.  A segment runs from its
    // donor junction (here a source/channel head, so its upslope is ~1 pixel)
    // DOWNSTREAM along the channel link to the matched DEM junction.  Walk the
    // D8 receiver chain from the donor junction to the DEM junction, collecting
    // every pixel of the link.
    //
    // (Collecting the donor's *upslope* set was the bug: a source junction has
    // essentially no upslope, so both segments collapsed to a single pixel and
    // the nearest-pixel search degenerated.)
    auto collect_link = [&](std::int64_t donor_flat)
    {
        std::vector<std::int64_t> seg;
        std::int64_t cur = donor_flat;
        seg.push_back(cur);
        const int max_steps = nr * nc;  // safety bound; walk terminates earlier
        int steps = 0;
        while (cur != dem_junction_flat && !flow.is_outlet(cur) && steps++ < max_steps)
        {
            const std::int64_t nxt = flow.receiver_of(cur);
            if (nxt == cur || nxt < 0) break;
            cur = nxt;
            seg.push_back(cur);
        }
        return seg;
    };

    std::vector<std::int64_t> channel_a = collect_link(branch_a_flat);
    std::vector<std::int64_t> channel_b = collect_link(branch_b_flat);

    if (channel_a.empty() || channel_b.empty())
    {
        std::cerr << "    parallel: empty channel set (A="
                  << channel_a.size() << " B=" << channel_b.size() << ")\n";
        return out;
    }

    std::cout << "    parallel: Branch_A pixels=" << channel_a.size()
              << "  Branch_B pixels=" << channel_b.size() << "\n";

    //  Step 3: find Branch A pixel nearest to vector junction 
    // This anchors the search for the breach start near the correct location.
    std::int64_t a_near_vj = branch_a_flat;
    {
        double best_d2 = std::numeric_limits<double>::max();
        for (std::int64_t f : channel_a)
        {
            const lsdtt_xy xy = cell_center_from_row_col(meta, to_row(f), to_col(f));
            const double dx = xy.x - breach_anchor.x;
            const double dy = xy.y - breach_anchor.y;
            const double d2 = dx*dx + dy*dy;
            if (d2 < best_d2) { best_d2 = d2; a_near_vj = f; }
        }
    }

    //  Step 4: find breach start via channel-score skeleton 
    //
    // Strategy:
    //   a) Build a binary mask of the channel-score raster restricted to a
    //      search window around the vector junction.  Use a low threshold so
    //      the small col/depression between the two valleys is included.
    //   b) Remove small blobs (noise).
    //   c) Skeletonize and prune short spurs.
    //   d) Extract degree-1 skeleton endpoints.
    //   e) From those endpoints, keep only ones that:
    //        - are on the Branch A (lower) side of the area (closer to channel_a
    //          than to channel_b)
    //        - point toward the higher channel (Branch B)
    //   f) The endpoint closest to a Branch B pixel is the breach start.
    //      If no suitable endpoint is found, fall back to a_near_vj.

    // Search window around the breach anchor (vector junction, or the
    // divergence point for a re-anchored diverged reach).
    const int vj_r = static_cast<int>((meta.y_max() - breach_anchor.y) / cs);
    const int vj_c = static_cast<int>((breach_anchor.x - meta.x_min) / cs);
    const int win_px = static_cast<int>(std::ceil(gap_distance_m / cs)) + 20;

    const int w_r0 = std::max(0,    vj_r - win_px);
    const int w_r1 = std::min(nr-1, vj_r + win_px);
    const int w_c0 = std::max(0,    vj_c - win_px);
    const int w_c1 = std::min(nc-1, vj_c + win_px);

    // Threshold: calibrate against the channel-score values found on the
    // confirmed DEM channel network.  The previous window-local 20th-percentile
    // rule discarded too many channelised cells (the skeleton came out too
    // sparse).  Instead, sample the channel score at every DEM-network channel
    // pixel and take the 5th percentile (lowest 5%) of those values: this
    // retains ~95% of genuine channel cells, giving a denser binary mask.
    std::vector<float> network_scores;
    network_scores.reserve(static_cast<std::size_t>(dem_network.n_channel_cells()));
    {
        const std::int64_t n_flat = static_cast<std::int64_t>(nr) * ncols64;
        for (std::int64_t f = 0; f < n_flat; ++f)
        {
            const std::size_t idx = static_cast<std::size_t>(f);
            if (mask_data[idx] && dem_network.is_channel(f))
                network_scores.push_back(channel_score_data[idx]);
        }
    }

    float score_thresh = 0.0f;
    if (!network_scores.empty())
    {
        std::sort(network_scores.begin(), network_scores.end());
        const std::size_t p05 = network_scores.size() / 20;  // lowest 5%
        score_thresh = network_scores[p05];
    }

    std::cout << "    parallel: channel score threshold=" << score_thresh
              << " (5th pct of " << network_scores.size()
              << " DEM-network channel pixels)  window="
              << (w_r1-w_r0+1) << "x" << (w_c1-w_c0+1) << " px\n";

    // Build windowed binary mask of channel score
    lsdtt_raster_cpu score_bin(
        lsdtt_array2d<float>(meta.n_rows, meta.n_cols, 0.0f),
        lsdtt_array2d<uint8_t>(dem.get_mask_array()),
        meta
    );
    {
        auto& sd = score_bin.get_data();
        for (int r = w_r0; r <= w_r1; ++r)
        {
            for (int c = w_c0; c <= w_c1; ++c)
            {
                const std::size_t idx = static_cast<std::size_t>(to_flat(r, c));
                if (mask_data[idx] &&
                    channel_score_data[idx] > 0.0f &&
                    channel_score_data[idx] >= score_thresh)
                    sd[idx] = 1.0f;
            }
        }
    }

    // Minimum component size: at least 3 pixels to suppress single-pixel noise
    const int min_blob = 3;
    lsdtt_raster_cpu filtered = filter_binary_components_by_size(score_bin, min_blob);

    // Keep only components that touch a DEM channel pixel.  Hillslope blobs
    // (high channel score locally but with no connection to the mapped channel
    // network) are discarded, which greatly cleans up the mask.
    {
        lsdtt_raster_cpu labelled = label_binary_components_by_size(filtered, 1);
        const auto& lab = labelled.get_data();

        // Component labels that contain at least one DEM channel pixel
        std::unordered_set<std::uint32_t> keep_labels;
        const std::int64_t n_flat = static_cast<std::int64_t>(nr) * ncols64;
        for (std::int64_t f = 0; f < n_flat; ++f)
        {
            const std::size_t idx = static_cast<std::size_t>(f);
            const std::uint32_t L = static_cast<std::uint32_t>(lab[idx]);
            if (L != 0 && dem_network.is_channel(f))
                keep_labels.insert(L);
        }

        auto& fd = filtered.get_data();
        for (std::int64_t f = 0; f < n_flat; ++f)
        {
            const std::size_t idx = static_cast<std::size_t>(f);
            if (fd[idx] > 0.5f &&
                keep_labels.count(static_cast<std::uint32_t>(lab[idx])) == 0)
                fd[idx] = 0.0f;
        }

        std::cout << "    parallel: kept " << keep_labels.size()
                  << " channel-connected component(s)\n";
    }

    // Report the binary mask before skeletonisation
    std::size_t n_mask_fg = 0;
    for (float v : filtered.get_data()) if (v > 0.5f) ++n_mask_fg;
    std::cout << "    parallel: binary mask foreground pixels (pre-skeleton)="
              << n_mask_fg << "\n";

    // Fill interior holes: background pixels enclosed by foreground.  Each such
    // hole skeletonises into a small ring/loop that suppresses endpoint
    // detection on the branch beyond it.  binary_fill_holes fills exactly the
    // background regions not connected to the domain exterior.
    filtered = binary_fill_holes(filtered);

    // Skeletonize.  thin_distance_guided is more robust than plain Zhang-Suen
    // on blocky / staircase masks: it protects the medial axis and produces far
    // fewer loop and diamond artefacts, so branch tips become clean endpoints.
    lsdtt_raster_cpu skeleton = thin_distance_guided(filtered, 1.0);

    // Prune spurs shorter than 2 pixels
    skeleton = prune_skeleton_spurs(skeleton, 2);

    // Break terminal "lollipop" loops.  Where a channel-score blob skeletonises
    // into a small closed ring at the tip, the true spur tip is a degree-2 loop
    // pixel and is never reported as an endpoint.  Pruning these rings exposes
    // the real endpoint pointing toward the col.
    skeleton = prune_terminal_loops(skeleton, 30);

    // Extract endpoints (degree-1 tips)
    const skeleton_node_result snr = classify_skeleton_nodes(skeleton);

    std::cout << "    parallel: skeleton endpoints=" << snr.endpoints.size() << "\n";

    // Write diagnostic rasters if prefix is set
    if (!params.diagnostic_prefix.empty())
    {
        try
        {
            lsdtt_data_loader dwriter;
            dwriter.write_full(score_bin.get_data(), meta,
                               params.diagnostic_prefix + "_score_binary.tif");
            dwriter.write_full(filtered.get_data(), meta,
                               params.diagnostic_prefix + "_binary_filtered.tif");
            dwriter.write_full(skeleton.get_data(), meta,
                               params.diagnostic_prefix + "_skeleton.tif");
            std::cout << "    parallel: wrote score_binary, binary_filtered and "
                         "skeleton diagnostics\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "    parallel: could not write diagnostics: " << e.what() << "\n";
        }

        // Write ALL classified skeleton endpoints (not just the lower-channel
        // spurs) so the skeleton's endpoint detection can be verified directly.
        const int  ep_ph = std::max(0, params.diagnostic_phase);
        const std::string ep_ph_str =
            (ep_ph < 10 ? std::string("0") : std::string()) + std::to_string(ep_ph);
        const std::string ep_file =
            params.diagnostic_prefix + "_skeleton_endpoints_phase_" + ep_ph_str + ".fgb";
        try
        {
            lsdtt_vector_schema schema;
            schema.field_names = { "row", "col", "degree" };
            schema.field_types = {
                lsdtt_vector_field_type::int64_type,
                lsdtt_vector_field_type::int64_type,
                lsdtt_vector_field_type::int64_type
            };
            lsdtt_data_loader ew;
            ew.begin_vector_write(ep_file, "skeleton_endpoints", schema,
                                  meta.crs_wkt,
                                  vector_output_format::flatgeobuf, false);
            for (const skeleton_node& ep : snr.endpoints)
            {
                const lsdtt_xy xy = cell_center_from_row_col(meta, ep.row, ep.col);
                std::vector<lsdtt_vector_field_value> attrs = {
                    static_cast<std::int64_t>(ep.row),
                    static_cast<std::int64_t>(ep.col),
                    static_cast<std::int64_t>(ep.degree)
                };
                ew.write_point(xy.x, xy.y, attrs);
            }
            ew.finish_vector_write();
            std::cout << "    parallel: wrote " << snr.endpoints.size()
                      << " skeleton endpoints to " << ep_file << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "    parallel: could not write skeleton endpoints: "
                      << e.what() << "\n";
        }
    }

    //  Step c: decide which branch is the LOWER channel and which is UPPER 
    //
    // Compare the elevation of each branch at its pixel nearest the vector
    // junction.  The repair carries flow from the upper (higher) channel down
    // into the lower channel.
    std::int64_t b_near_vj = branch_b_flat;
    {
        double best_d2 = std::numeric_limits<double>::max();
        for (std::int64_t f : channel_b)
        {
            const lsdtt_xy xy = cell_center_from_row_col(meta, to_row(f), to_col(f));
            const double dx = xy.x - vector_junction_xy.x;
            const double dy = xy.y - vector_junction_xy.y;
            const double d2 = dx*dx + dy*dy;
            if (d2 < best_d2) { best_d2 = d2; b_near_vj = f; }
        }
    }

    const float elev_a_near = elev_data[static_cast<std::size_t>(a_near_vj)];
    const float elev_b_near = elev_data[static_cast<std::size_t>(b_near_vj)];
    const bool  a_is_lower  = elev_a_near <= elev_b_near;

    const std::vector<std::int64_t>& lower_channel = a_is_lower ? channel_a : channel_b;
    const std::vector<std::int64_t>& upper_channel = a_is_lower ? channel_b : channel_a;

    std::cout << "    parallel: lower channel = Branch_" << (a_is_lower ? "A" : "B")
              << " (elev " << (a_is_lower ? elev_a_near : elev_b_near) << " m);"
              << " upper channel = Branch_" << (a_is_lower ? "B" : "A")
              << " (elev " << (a_is_lower ? elev_b_near : elev_a_near) << " m)\n";

    const std::unordered_set<std::int64_t> lower_set(
        lower_channel.begin(), lower_channel.end());
    const std::unordered_set<std::int64_t> upper_set(
        upper_channel.begin(), upper_channel.end());

    //  Step d: select the breach endpoint by D8 connectivity 
    //
    // Endpoint selection must NOT be keyed on how close the skeleton tip itself
    // is to the vector junction - the correct tip can sit on a long spur whose
    // end lies well outside the search radius.  Instead we use where each
    // endpoint's flow *joins the lower channel*:
    //
    //   1. From each skeleton endpoint, follow D8 downstream until it reaches
    //      one of the two channels.
    //   2. Discard endpoints whose flow reaches the UPPER channel; keep only
    //      those that drain to the LOWER channel.
    //   3. Record where each kept endpoint joins the lower channel.
    //   4. Keep only endpoints whose join point lies within the search radius
    //      of the (vector) target junction.
    //   5/6. Among those, choose the endpoint closest to a pixel on the upper
    //      channel; burn from that nearest upper pixel to the endpoint, then
    //      follow D8 from the endpoint down to the lower channel.
    const double search_radius_m =
        std::max(match_radius_m > 0.0 ? match_radius_m : 0.0, 10.0 * cs);
    const double radius2 = search_radius_m * search_radius_m;

    // Trace D8 from a start flat to the first lower/upper channel pixel.
    // Returns 0 if it reached the lower channel, 1 if the upper channel,
    // -1 if neither.  On a hit, join_flat holds the channel pixel reached.
    auto trace_to_channel = [&](std::int64_t start, std::int64_t& join_flat) -> int
    {
        std::int64_t cur = start;
        const int max_steps = nr + nc + 8 * win_px;  // generous bound
        int steps = 0;
        while (steps++ < max_steps)
        {
            if (lower_set.count(cur)) { join_flat = cur; return 0; }
            if (upper_set.count(cur)) { join_flat = cur; return 1; }
            if (flow.is_outlet(cur)) break;
            const std::int64_t nxt = flow.receiver_of(cur);
            if (nxt == cur || nxt < 0) break;
            cur = nxt;
        }
        return -1;
    };

    struct breach_candidate
    {
        std::int64_t ep_flat;
        lsdtt_xy     ep_xy;
        std::int64_t join_flat;           // where its D8 path meets the lower channel
        double       join_dist_to_vj_m;
        std::int64_t nearest_upper_flat;
        double       dist_to_upper_m;
    };
    std::vector<breach_candidate> cands;

    for (const skeleton_node& ep : snr.endpoints)
    {
        const std::int64_t ep_flat = to_flat(ep.row, ep.col);

        // (1)(2) D8 from the endpoint; keep only lower-channel-draining tips
        std::int64_t join_flat = -1;
        const int which = trace_to_channel(ep_flat, join_flat);
        if (which != 0) continue;         // no channel, or drains to the upper channel

        // (3)(4) Proximity gate. Normally the JOIN (where the tip's D8 reaches the
        // lower channel) must lie within the search radius of the target junction.
        // BUT for a diverged reach the valley below the obstruction is DRY by
        // definition of being an obstruction, so its re-established channel - and
        // hence the join - is far downstream. There we instead require the TIP itself
        // (where the breach lands) to be near the obstruction (the breach anchor): the
        // burn notches the obstruction onto the dry-valley skeleton and D8 then carries
        // the flow down to the re-established channel (the join, wherever it is).
        // (3)(4) join point must lie within the search radius of the junction
        const lsdtt_xy join_xy =
            cell_center_from_row_col(meta, to_row(join_flat), to_col(join_flat));
        const double jdx = join_xy.x - breach_anchor.x;
        const double jdy = join_xy.y - breach_anchor.y;
        const double join_d2 = jdx*jdx + jdy*jdy;
        if (join_d2 > radius2) continue;

        // nearest upper-channel pixel to the endpoint
        const lsdtt_xy ep_xy = cell_center_from_row_col(meta, ep.row, ep.col);
        double       d2_upper      = std::numeric_limits<double>::max();
        std::int64_t nearest_upper = -1;
        for (std::int64_t f : upper_channel)
        {
            const lsdtt_xy xy = cell_center_from_row_col(meta, to_row(f), to_col(f));
            const double dx = ep_xy.x - xy.x, dy = ep_xy.y - xy.y;
            const double d2 = dx*dx + dy*dy;
            if (d2 < d2_upper) { d2_upper = d2; nearest_upper = f; }
        }

        breach_candidate bc;
        bc.ep_flat            = ep_flat;
        bc.ep_xy              = ep_xy;
        bc.join_flat          = join_flat;
        bc.join_dist_to_vj_m  = std::sqrt(join_d2);
        bc.nearest_upper_flat = nearest_upper;
        bc.dist_to_upper_m    = std::sqrt(d2_upper);
        cands.push_back(bc);
    }

    std::cout << "    parallel: lower-channel endpoints joining within "
              << search_radius_m << " m of vector junction = "
              << cands.size() << "\n";

    // (6) choose the endpoint closest to the upper channel (narrowest breach)
    int    best_idx = -1;
    double best_d   = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < cands.size(); ++i)
        if (cands[i].dist_to_upper_m < best_d)
        { best_d = cands[i].dist_to_upper_m; best_idx = static_cast<int>(i); }

    //  Write all breach candidates to a diagnostic FGB for bug checking 
    if (!params.diagnostic_prefix.empty())
    {
        const int  ph_i = std::max(0, params.diagnostic_phase);
        const std::string ph_str =
            (ph_i < 10 ? std::string("0") : std::string()) + std::to_string(ph_i);
        const std::string spur_file =
            params.diagnostic_prefix + "_spur_endpoints_phase_" + ph_str + ".fgb";
        try
        {
            lsdtt_vector_schema schema;
            schema.field_names = {
                "cand_index", "is_selected", "dist_to_upper_m", "join_dist_to_vj_m"
            };
            schema.field_types = {
                lsdtt_vector_field_type::int64_type,
                lsdtt_vector_field_type::int64_type,
                lsdtt_vector_field_type::double_type,
                lsdtt_vector_field_type::double_type
            };
            lsdtt_data_loader sw;
            sw.begin_vector_write(spur_file, "spur_endpoints", schema,
                                  meta.crs_wkt,
                                  vector_output_format::flatgeobuf, false);
            for (std::size_t i = 0; i < cands.size(); ++i)
            {
                std::vector<lsdtt_vector_field_value> attrs = {
                    static_cast<std::int64_t>(i),
                    std::int64_t{ static_cast<int>(i) == best_idx ? 1 : 0 },
                    cands[i].dist_to_upper_m,
                    cands[i].join_dist_to_vj_m
                };
                sw.write_point(cands[i].ep_xy.x, cands[i].ep_xy.y, attrs);
            }
            sw.finish_vector_write();
            std::cout << "    parallel: wrote " << cands.size()
                      << " breach candidates to " << spur_file << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "    parallel: could not write spur endpoints: "
                      << e.what() << "\n";
        }
    }

    if (best_idx < 0)
    {
        std::cout << "    parallel: no lower-channel endpoint joins near the "
                     "vector junction; leaving unresolved\n";
        return out;
    }

    const breach_candidate& chosen   = cands[static_cast<std::size_t>(best_idx)];
    const std::int64_t      spur_flat = chosen.ep_flat;
    const std::int64_t      burn_init = chosen.nearest_upper_flat;  // on upper channel

    const int spur_r = to_row(spur_flat),  spur_c = to_col(spur_flat);
    const int init_r = to_row(burn_init),  init_c = to_col(burn_init);

    std::cout << "    parallel: chosen endpoint=("
              << chosen.ep_xy.x << "," << chosen.ep_xy.y << ")"
              << " dist_to_upper=" << chosen.dist_to_upper_m << " m;"
              << " join_dist_to_vj=" << chosen.join_dist_to_vj_m << " m;"
              << " burn init (upper)=("
              << cell_center_from_row_col(meta, init_r, init_c).x << ","
              << cell_center_from_row_col(meta, init_r, init_c).y << ")\n";

    //  Burn from the upper-channel initiation point down to the spur 
    const int buf = std::max(params.max_search_px, win_px);

    lcp_params lcp;
    lcp.uphill_mult         = params.uphill_mult;
    lcp.channel_score       = &channel_score_data;
    lcp.channel_score_alpha = params.channel_score_alpha;
    lcp.goal_flats          = std::vector<std::int64_t>{ spur_flat };
    lcp.r_min = std::max(0,    std::min(init_r, spur_r) - buf);
    lcp.r_max = std::min(nr-1, std::max(init_r, spur_r) + buf);
    lcp.c_min = std::max(0,    std::min(init_c, spur_c) - buf);
    lcp.c_max = std::min(nc-1, std::max(init_c, spur_c) + buf);

    lcp_result lcp_res;
    try { lcp_res = raster_least_cost_path(dem, init_r, init_c, spur_r, spur_c, lcp); }
    catch (const std::exception& e)
    { std::cerr << "    parallel LCP failed: " << e.what() << "\n"; return out; }

    if (!lcp_res.found)
    { std::cout << "    parallel: no path from upper channel to spur\n"; return out; }

    out.max_climb_m = lcp_res.max_climb;
    out.total_cost  = lcp_res.total_cost;

    if (lcp_res.max_climb > params.max_climb_m)
    {
        std::cout << "    parallel: climb=" << lcp_res.max_climb
                  << " m > " << params.max_climb_m << " m; rejected\n";
        return out;
    }

    burn_monotonic_path(dem, lcp_res.path, params.repair_min_slope);

    //  Follow D8 from the spur endpoint down to the lower channel 
    std::vector<std::int64_t> d8_tail;
    {
        std::int64_t cur     = spur_flat;
        bool         reached = lower_set.count(cur) > 0;
        int          steps   = 0;
        const int    max_steps = std::max(1000, buf * 8);
        while (!reached && steps++ < max_steps && !flow.is_outlet(cur))
        {
            const std::int64_t nxt = flow.receiver_of(cur);
            if (nxt == cur || nxt < 0) break;
            cur = nxt;
            d8_tail.push_back(cur);
            if (lower_set.count(cur)) { reached = true; break; }
        }
        if (!reached)
            std::cout << "    parallel: warning - D8 from spur did not reach the "
                         "lower channel within " << max_steps << " steps\n";
    }

    // Full repair path: upper channel -> spur (burned) -> (D8) -> lower channel
    out.burned_path = lcp_res.path;
    out.burned_path.insert(out.burned_path.end(), d8_tail.begin(), d8_tail.end());

    out.n_burned         = static_cast<int>(lcp_res.path.size());
    out.n_landing_pixels = static_cast<int>(lower_channel.size());
    out.fix_type         = GapFixType::channel_guided_bridge;

    // divergence = upper-channel initiation; target = lower-channel landing
    out.divergence_xy = cell_center_from_row_col(meta, init_r, init_c);
    if (!out.burned_path.empty())
    {
        const std::int64_t end_flat = out.burned_path.back();
        out.target_xy = cell_center_from_row_col(
            meta, to_row(end_flat), to_col(end_flat));
    }

    std::cout << "    parallel: burned " << out.n_burned
              << " px (climb=" << lcp_res.max_climb << " m), D8 tail="
              << d8_tail.size() << " px"
              << "  start=(" << out.divergence_xy.x << "," << out.divergence_xy.y << ")"
              << "  end=("   << out.target_xy.x     << "," << out.target_xy.y     << ")\n";
    return out;
}


// ===========================================================================
// repair_junction_channel_guided
// ===========================================================================
static repair_attempt_result repair_junction_channel_guided(
    lsdtt_raster_cpu&                  dem,
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& dem_network,
    std::int64_t                       dem_junction_flat,
    const lsdtt_xy&                    vector_junction_xy,
    double                             gap_distance_m,
    double                             match_radius_m,
    const std::vector<float>&          channel_score_data,
    const junction_repair_params&      params)
{
    repair_attempt_result out;
    out.fix_type = GapFixType::unresolved;

    const lsdtt_raster_meta& meta  = dem.get_meta();
    const double cs = meta.cell_size;
    const std::int64_t ncols = static_cast<std::int64_t>(meta.n_cols);

    auto to_row = [&](std::int64_t f) -> int { return static_cast<int>(f / ncols); };
    auto to_col = [&](std::int64_t f) -> int { return static_cast<int>(f % ncols); };

    //  Phase 1: find divergence point 
    //
    // Score upstream channel pixels against the MIDPOINT of the gap, not the
    // vector junction end.  Scoring against the vector junction end selects
    // pixels that are already near the landing zone, producing a trivially
    // short path that does not redirect the channel.  The midpoint target
    // selects a pixel somewhere along the mis-routed reach, from which a
    // meaningful LCP can bridge across to the target valley.
    //
    // Candidates within match_radius_m of the vector junction are excluded:
    // those pixels are landing-zone candidates, not upstream divergence points.
    const std::vector<std::int64_t> upstream =
        collect_upslope_channel_pixels(
            flow, dem_network, dem_junction_flat,
            vector_junction_xy, 2.0 * gap_distance_m, meta);

    // Midpoint between DEM junction and vector junction
    const lsdtt_xy dem_junction_xy =
        cell_center_from_row_col(meta,
            to_row(dem_junction_flat), to_col(dem_junction_flat));

    const lsdtt_xy midpoint{
        (dem_junction_xy.x + vector_junction_xy.x) * 0.5,
        (dem_junction_xy.y + vector_junction_xy.y) * 0.5
    };

    std::int64_t divergence_flat = dem_junction_flat;

    if (!upstream.empty())
    {
        double best_score = -1.0;
        int    n_eligible = 0;

        for (std::int64_t flat : upstream)
        {
            const lsdtt_xy xy = cell_center_from_row_col(meta, to_row(flat), to_col(flat));

            // Exclude pixels already in the landing zone
            const double dv_x = xy.x - vector_junction_xy.x;
            const double dv_y = xy.y - vector_junction_xy.y;
            if (std::sqrt(dv_x*dv_x + dv_y*dv_y) < match_radius_m)
                continue;

            ++n_eligible;

            // Proximity to midpoint - scale denominator by half the gap
            const double dm_x = xy.x - midpoint.x;
            const double dm_y = xy.y - midpoint.y;
            const double dist = std::sqrt(dm_x*dm_x + dm_y*dm_y);
            const double half_gap = std::max(gap_distance_m * 0.5, cs);
            const double prox = 1.0 / (1.0 + dist / half_gap);
            const double chan = static_cast<double>(channel_score_data[flat]);
            const double s = params.w_proximity * prox + params.w_channel * chan;
            if (s > best_score) { best_score = s; divergence_flat = flat; }
        }

        if (n_eligible == 0)
        {
            // All upstream pixels are already at the landing end - fall back
            // to the DEM junction as the divergence point
            std::cout << "    all upstream candidates in landing zone; "
                         "using DEM junction as divergence point\n";
        }
        else
        {
            std::cout << "    divergence point selected  candidates="
                      << upstream.size() << " eligible=" << n_eligible << "\n";
        }
    }
    else
    {
        std::cout << "    no upstream channel pixels within 2x gap; "
                     "using DEM junction as divergence point\n";
    }

    out.divergence_xy = cell_center_from_row_col(
        meta, to_row(divergence_flat), to_col(divergence_flat));

    //  Phase 2: landing zone near vector junction 
    const int vj_row = static_cast<int>(
        std::floor((meta.y_max() - vector_junction_xy.y) / cs));
    const int vj_col = static_cast<int>(
        std::floor((vector_junction_xy.x - meta.x_min) / cs));
    const int landing_px = static_cast<int>(std::ceil(match_radius_m / cs));

    const std::vector<uint8_t>& dem_mask = dem.get_mask();
    std::vector<std::int64_t> landing_flats;

    for (int dr = -landing_px; dr <= landing_px; ++dr)
    {
        for (int dc = -landing_px; dc <= landing_px; ++dc)
        {
            const int r = vj_row + dr;
            const int c = vj_col + dc;
            if (r < 0 || r >= meta.n_rows || c < 0 || c >= meta.n_cols) continue;
            const double dist_m = std::sqrt(static_cast<double>(dr*dr + dc*dc)) * cs;
            if (dist_m > match_radius_m) continue;
            const std::int64_t flat = static_cast<std::int64_t>(r) * ncols + c;
            if (!dem_mask[static_cast<std::size_t>(flat)]) continue;
            if (static_cast<double>(channel_score_data[flat]) >= params.landing_score_threshold)
                landing_flats.push_back(flat);
        }
    }

    if (landing_flats.empty())
    {
        if (vj_row >= 0 && vj_row < meta.n_rows &&
            vj_col >= 0 && vj_col < meta.n_cols)
        {
            const std::int64_t vj_flat = static_cast<std::int64_t>(vj_row) * ncols + vj_col;
            if (dem_mask[static_cast<std::size_t>(vj_flat)])
                landing_flats.push_back(vj_flat);
        }
    }

    out.n_landing_pixels = static_cast<int>(landing_flats.size());

    if (landing_flats.empty())
    {
        std::cerr << "    repair_junction_channel_guided: no valid landing pixels\n";
        return out;
    }

    //  Phase 3: channel-score-guided Dijkstra 
    const int div_row = to_row(divergence_flat);
    const int div_col = to_col(divergence_flat);
    const int buf = std::max(params.max_search_px, landing_px + 10);

    lcp_params lcp;
    lcp.uphill_mult         = params.uphill_mult;
    lcp.channel_score       = &channel_score_data;
    lcp.channel_score_alpha = params.channel_score_alpha;
    lcp.goal_flats          = landing_flats;
    lcp.r_min = std::max(0,             std::min(div_row, vj_row) - buf);
    lcp.r_max = std::min(meta.n_rows-1, std::max(div_row, vj_row) + buf);
    lcp.c_min = std::max(0,             std::min(div_col, vj_col) - buf);
    lcp.c_max = std::min(meta.n_cols-1, std::max(div_col, vj_col) + buf);

    lcp_result lcp_res;
    try
    {
        lcp_res = raster_least_cost_path(dem, div_row, div_col, vj_row, vj_col, lcp);
    }
    catch (const std::exception& e)
    {
        std::cerr << "    channel-guided LCP failed: " << e.what() << "\n";
        return out;
    }

    if (!lcp_res.found)
    {
        std::cout << "    channel-guided: no path found in search window\n";
        return out;
    }

    out.max_climb_m = lcp_res.max_climb;
    out.total_cost  = lcp_res.total_cost;

    if (lcp_res.max_climb > params.max_climb_m)
    {
        std::cout << "    channel-guided path rejected: climb=" << lcp_res.max_climb
                  << " m > " << params.max_climb_m << " m\n";
        return out;
    }

    burn_monotonic_path(dem, lcp_res.path, params.repair_min_slope);
    out.n_burned    = static_cast<int>(lcp_res.path.size());
    out.fix_type    = GapFixType::channel_guided_bridge;
    out.burned_path = lcp_res.path;
    if (!lcp_res.path.empty())
    {
        const std::int64_t end_flat = lcp_res.path.back();
        const lsdtt_raster_meta& m = dem.get_meta();
        out.target_xy = cell_center_from_row_col(
            m, static_cast<int>(end_flat / m.n_cols),
               static_cast<int>(end_flat % m.n_cols));
    }

    std::cout << "    channel-guided: burned " << out.n_burned
              << " px, climb=" << lcp_res.max_climb << " m\n";
    return out;
}

// ===========================================================================
// repair_junction_terrain_fallback
// ===========================================================================
static repair_attempt_result repair_junction_terrain_fallback(
    lsdtt_raster_cpu&              dem,
    const std::string&             vector_path,
    const CandidateJunctionRecord& candidate,
    const junction_repair_params&  params)
{
    repair_attempt_result out;
    out.fix_type = GapFixType::unresolved;

    const lsdtt_raster_meta& meta = dem.get_meta();

    int rv = -1, cv = -1, rd = -1, cd = -1;
    if (!xy_to_valid_row_col(dem, candidate.vector_xy, rv, cv) ||
        !xy_to_valid_row_col(dem, candidate.dem_xy,    rd, cd))
        return out;

    lcp_params lcp;
    lcp.uphill_mult = params.uphill_mult;
    lcp.goal_radius = std::max(3, params.max_search_px / 4);
    lcp.r_min = std::max(0,             std::min(rv, rd) - params.max_search_px);
    lcp.r_max = std::min(meta.n_rows-1, std::max(rv, rd) + params.max_search_px);
    lcp.c_min = std::max(0,             std::min(cv, cd) - params.max_search_px);
    lcp.c_max = std::min(meta.n_cols-1, std::max(cv, cd) + params.max_search_px);

    lcp_result lcp_res;
    try { lcp_res = raster_least_cost_path(dem, rd, cd, rv, cv, lcp); }
    catch (const std::exception& e)
    { std::cerr << "    terrain fallback LCP failed: " << e.what() << "\n"; }

    if (lcp_res.found)
    {
        out.max_climb_m = lcp_res.max_climb;
        out.total_cost  = lcp_res.total_cost;
        if (lcp_res.max_climb <= params.max_climb_m)
        {
            burn_monotonic_path(dem, lcp_res.path, params.repair_min_slope);
            out.fix_type    = GapFixType::terrain_bridge;
            out.n_burned    = static_cast<int>(lcp_res.path.size());
            out.burned_path = lcp_res.path;
            if (!lcp_res.path.empty())
            {
                const std::int64_t ef = lcp_res.path.back();
                const lsdtt_raster_meta& m = dem.get_meta();
                out.target_xy = cell_center_from_row_col(
                    m, static_cast<int>(ef / m.n_cols),
                       static_cast<int>(ef % m.n_cols));
            }
            return out;
        }
    }

    // Vector burn fallback
    const double cs = meta.cell_size;
    const double x_lo = meta.x_min + lcp.c_min * cs;
    const double x_hi = meta.x_min + (lcp.c_max + 1) * cs;
    const double y_lo = meta.y_min + (meta.n_rows - lcp.r_max - 1) * cs;
    const double y_hi = meta.y_min + (meta.n_rows - lcp.r_min) * cs;

    const int nb = burn_vector_in_window(dem, vector_path, x_lo, x_hi, y_lo, y_hi);
    if (nb > 0)
    {
        out.fix_type = GapFixType::vector_burn;
        out.n_burned = nb;
    }
    return out;
}

// ===========================================================================
// repair_candidate_junctions
// ===========================================================================
std::vector<JunctionRepairRecord>
repair_candidate_junctions(
    lsdtt_raster_cpu&                             dem,
    const lsdtt_vector_to_junction_network&       vector_network,
    const std::string&                            vector_path,
    const std::vector<CandidateJunctionRecord>&   candidates,
    const std::vector<std::int64_t>&              source_flats,
    std::unique_ptr<lsdtt_flow_base>&             flow,
    std::unique_ptr<lsdtt_junction_network_base>& dem_network,
    const lsdtt_flowinfo_options&                 fi_opts,
    const lsdtt_junction_network_options&         jn_opts,
    const junction_repair_params&                 params_in,
    int                                           phase,
    int                                           repair_index_offset,
    const std::vector<float>*                     precomputed_channel_scores)
{
    std::vector<JunctionRepairRecord> records;
    records.reserve(candidates.size());
    if (candidates.empty()) return records;

    // Local copy so per-phase diagnostic outputs can be tagged with the phase.
    junction_repair_params params = params_in;
    params.diagnostic_phase = phase;

    const lsdtt_raster_meta& meta = dem.get_meta();

    const bool use_channel_guidance =
        params.channel_score_alpha > 0.0 && params.w_channel > 0.0;

    bool               have_channel_scores = false;
    std::vector<float> channel_score_data;

    if (precomputed_channel_scores && !precomputed_channel_scores->empty())
    {
        // Caller already computed scores - reuse them
        channel_score_data  = *precomputed_channel_scores;
        have_channel_scores = true;
        std::cout << "repair_candidate_junctions: using precomputed channel scores\n";
    }
    else if (use_channel_guidance)
    {
        const double cs = meta.cell_size;
        channel_candidate_request cc_req;
        cc_req.channel_lambda_min_metres =
            (params.channel_lambda_min_m > 0.0) ? params.channel_lambda_min_m : cs * 3.0;
        cc_req.channel_lambda_max_metres =
            (params.channel_lambda_max_m > 0.0) ? params.channel_lambda_max_m : cs * 10.0;
        if (cc_req.channel_lambda_max_metres <= cc_req.channel_lambda_min_metres)
            cc_req.channel_lambda_max_metres = cc_req.channel_lambda_min_metres * 3.0;
        cc_req.channel_n_scales           = params.channel_n_scales;
        cc_req.weight_persistent          = 0.7;
        cc_req.weight_depth               = 0.3;
        cc_req.weight_slope               = 0.0;
        cc_req.use_hessian_filter         = true;
        cc_req.positive_curvature_is_good = true;

        std::cout << "repair_candidate_junctions: computing channel scores "
                  << "(lambda " << cc_req.channel_lambda_min_metres
                  << "-" << cc_req.channel_lambda_max_metres
                  << " m, " << cc_req.channel_n_scales << " scales)\n";
        try
        {
            channel_candidate_result cc = find_channel_candidates(dem, cc_req);
            channel_score_data  = cc.score.get_data();
            have_channel_scores = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "  Warning: channel score computation failed: "
                      << e.what() << "\n"
                      << "  Repairs will use terrain fallback only.\n";
        }
    }

    int repair_index = repair_index_offset;

    for (const CandidateJunctionRecord& c : candidates)
    {
        JunctionRepairRecord rec;
        rec.repair_index          = repair_index++;
        rec.phase                 = phase;
        rec.basin_outlet_junction = c.basin_outlet_junction;
        rec.vector_junction_id    = c.vector_junction_id;
        rec.dem_junction_id       = c.dem_junction_id;
        rec.vector_xy             = c.vector_xy;
        rec.dem_xy                = c.dem_xy;
        rec.distance_m            = c.match_distance_m;
        rec.mismatch_reason       = c.mismatch_reason;
        rec.fix_type              = GapFixType::unresolved;

        if (c.dem_junction_id < 0 || c.match_distance_m < 0.0)
        {
            std::cout << "  candidate vector_jid=" << c.vector_junction_id
                      << ": no DEM counterpart; skipping\n";
            records.push_back(rec);
            continue;
        }

        if (c.match_distance_m < params.min_gap_m)
        {
            rec.fix_type = GapFixType::ok;
            rec.success  = true;
            records.push_back(rec);
            continue;
        }

        // Gap too large to be a physical obstruction: the DEM flow is correct
        // and the vector network is a coarse local representation. Report and
        // skip rather than carve a long, spurious channel.
        if (params.max_repair_gap_m > 0.0 &&
            c.match_distance_m > params.max_repair_gap_m)
        {
            std::cout << "  skipping candidate vector_jid=" << c.vector_junction_id
                      << " basin=" << c.basin_outlet_junction
                      << " gap=" << c.match_distance_m
                      << " m > max_repair_gap_m=" << params.max_repair_gap_m
                      << " m (topology mismatch, not an obstruction; reason="
                      << c.mismatch_reason << ")\n";
            rec.fix_type = GapFixType::unresolved;
            rec.success  = false;
            records.push_back(rec);
            continue;
        }

        std::cout << "  repairing candidate vector_jid=" << c.vector_junction_id
                  << " basin=" << c.basin_outlet_junction
                  << " gap=" << c.match_distance_m
                  << " m reason=" << c.mismatch_reason << "\n";

        repair_attempt_result attempt;

        const std::int64_t dem_junc_flat =
            dem_network->flat_of_junction(
                static_cast<std::int32_t>(c.dem_junction_id));

        // Attempt 1: parallel channel repair (handles parallel-swale case)
        if (have_channel_scores)
        {
            try
            {
                attempt = repair_junction_parallel_channels(
                    dem, *flow, *dem_network,
                    dem_junc_flat, c.vector_xy,
                    c.vector_junction_id, vector_network,
                    c.match_distance_m, params.min_gap_m,
                    channel_score_data, params);
            }
            catch (const std::exception& e)
            {
                std::cerr << "    parallel channel repair exception: "
                          << e.what() << "\n";
                attempt.fix_type = GapFixType::unresolved;
            }
        }

        // Attempt 2: channel-guided fallback (barrier/crossing case)
        if (attempt.fix_type == GapFixType::unresolved && have_channel_scores)
        {
            try
            {
                attempt = repair_junction_channel_guided(
                    dem, *flow, *dem_network,
                    dem_junc_flat, c.vector_xy,
                    c.match_distance_m, params.min_gap_m,
                    channel_score_data, params);
            }
            catch (const std::exception& e)
            {
                std::cerr << "    channel-guided repair exception: "
                          << e.what() << "\n";
                attempt.fix_type = GapFixType::unresolved;
            }
        }

        // Attempt 3: terrain fallback
        if (attempt.fix_type == GapFixType::unresolved)
        {
            attempt = repair_junction_terrain_fallback(dem, vector_path, c, params);
        }

        rec.fix_type          = attempt.fix_type;
        rec.n_burned          = attempt.n_burned;
        rec.max_climb_m       = attempt.max_climb_m;
        rec.total_cost        = attempt.total_cost;
        rec.n_landing_pixels  = attempt.n_landing_pixels;
        rec.divergence_xy     = attempt.divergence_xy;
        rec.target_xy         = attempt.target_xy;
        rec.burned_path       = std::move(attempt.burned_path);

        if (attempt.fix_type != GapFixType::unresolved &&
            attempt.fix_type != GapFixType::ok)
        {
            // A carve can leave a residual pit/flat. make_flow_base is OpenMP
            // and TERMINATES the process on an unresolved sink (uncatchable),
            // so FILL any residual depression the burn left BEFORE reflowing -
            // a capturing burn that daylights is unchanged by the fill; one
            // that accidentally pitted is resolved instead of crashing.
            const float fill_slope = (params.repair_min_slope > 0.0)
                ? static_cast<float>(params.repair_min_slope) : 1.0e-4f;
            lsdtt_raster_cpu filled = fill_pits(dem, fill_slope, fi_opts.n_threads);
            dem.get_data_array() = filled.get_data_array();

            flow = make_flow_base(
                dem.get_data_array(), dem.get_mask_array(), meta, fi_opts);
            dem_network = make_junction_network_base(*flow, source_flats, jn_opts);

            int nearest_j = -1;
            rec.post_distance_m = nearest_junction_distance_to_xy(
                *dem_network, meta, c.vector_xy, &nearest_j);
            rec.success = rec.post_distance_m <= params.min_gap_m;

            std::cout << "    result: " << gap_fix_type_str(rec.fix_type)
                      << " burned=" << rec.n_burned
                      << " post_gap=" << rec.post_distance_m
                      << " m success=" << (rec.success ? "true" : "false") << "\n";
        }
        else
        {
            rec.post_distance_m = c.match_distance_m;
            rec.success = false;
            std::cout << "    repair unresolved\n";
        }

        records.push_back(rec);
    }

    return records;
}

// ===========================================================================
// repair_diverged_reach_capture
//
// Capturing-burn repair for a diverged-but-rejoining reach (road/ditch capture
// where the parallel-channel col-finder fails because the road sits high on the
// valley side, not at a col). The vector edge is lower-resolution (off by
// 10-20 m) and must NOT be burned under; it is used only as a CORRIDOR guide.
// Inside a narrow corridor around the vector edge (narrower than the road's
// separation, so the road is excluded) we trace the local DEM thalweg (lowest
// ground) from the divergence point down to a DEM channel near the downstream
// junction, then burn_capturing_path along it so D8 leaves the road and follows
// the true valley. Every burned pixel is a DEM pixel (a local low), never
// vector geometry. burn_capturing_path daylights at the terminus, so no sink.
// ===========================================================================
static repair_attempt_result repair_diverged_reach_capture(
    lsdtt_raster_cpu&                  dem,
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& dem_network,
    const ReachDivergenceRecord&       reach,
    const lsdtt_xy&                    downstream_vector_xy,
    double                             corridor_half_width_m,
    const std::vector<float>&          channel_score_data,
    const junction_repair_params&      params)
{
    repair_attempt_result out;
    out.fix_type = GapFixType::unresolved;
    if (!reach.has_divergence_point) return out;

    const lsdtt_raster_meta& meta = dem.get_meta();
    const double cs = meta.cell_size;
    const int    nr = meta.n_rows;
    const int    nc = meta.n_cols;
    const std::int64_t ncols64 = nc;
    const auto& mask = dem.get_mask();

    auto to_row  = [&](std::int64_t f){ return static_cast<int>(f / ncols64); };
    auto to_col  = [&](std::int64_t f){ return static_cast<int>(f % ncols64); };
    auto to_flat = [&](int r,int c){ return static_cast<std::int64_t>(r)*ncols64 + c; };
    auto valid   = [&](int r,int c){
        return r>=0 && r<nr && c>=0 && c<nc &&
               mask[static_cast<std::size_t>(to_flat(r,c))]!=0; };
    auto rc_of_xy = [&](const lsdtt_xy& p, int& r, int& c){
        c = static_cast<int>(std::floor((p.x - meta.x_min)/cs));
        r = static_cast<int>(std::floor((meta.y_max() - p.y)/cs)); };

    // Start = divergence-point pixel.
    int sr, sc; rc_of_xy(reach.divergence_point_xy, sr, sc);
    if (!valid(sr,sc)) return out;
    const std::int64_t start_flat = to_flat(sr,sc);

    int gr0, gc0; rc_of_xy(downstream_vector_xy, gr0, gc0);

    // Contributing-area raster (D8 flow accumulation) - also used below to snap the
    // skeleton to the D8 channel and to detect the downstream confluence.
    const std::vector<float> area = flow.get_n_contributing_pixels_raster().data;
    const float div_area = area.empty()?0.0f:area[static_cast<std::size_t>(start_flat)];

    // Road pixels: the D8 channel the diverged reach currently follows from the
    // divergence point. CRITICAL: stop at the REJOIN - once this D8 path merges into
    // a much-larger channel (the trunk / correct downstream reach), area jumps and
    // everything beyond is the legitimate downstream channel, NOT road. Walking all
    // the way to the outlet (the old behaviour) put the trunk into road_set, which
    // blocked the repair from ever following D8 onto the trunk and stopping there.
    std::vector<std::int64_t> road;
    {
        std::int64_t cur = start_flat; int steps=0; const int maxs = nr+nc+1;
        while (steps++ < maxs){
            if (cur!=start_flat &&
                area[static_cast<std::size_t>(cur)] >= 5.0f*std::max(div_area,1.0f))
                break;  // merged into the trunk: the rest is the downstream reach
            road.push_back(cur);
            if (flow.is_outlet(cur)) break;
            const std::int64_t nxt = flow.receiver_of(cur);
            if (nxt<0 || nxt==cur) break;
            cur = nxt;
        }
    }
    std::unordered_set<std::int64_t> road_set(road.begin(), road.end());

    // Corridor window: around the vector edge from the divergence point to the
    // downstream junction, padded by the corridor half-width.
    const int corr_px = std::max(4, static_cast<int>(std::ceil(corridor_half_width_m/cs)));
    int r_min=std::min(sr,gr0), r_max=std::max(sr,gr0);
    int c_min=std::min(sc,gc0), c_max=std::max(sc,gc0);
    for (const auto& p : reach.vector_reach_xy){
        int pr,pc; rc_of_xy(p,pr,pc);
        r_min=std::min(r_min,pr); r_max=std::max(r_max,pr);
        c_min=std::min(c_min,pc); c_max=std::max(c_max,pc);
    }
    const int pad = corr_px + 4;
    const int wr_min=std::max(0,r_min-pad), wr_max=std::min(nr-1,r_max+pad);
    const int wc_min=std::max(0,c_min-pad), wc_max=std::min(nc-1,c_max+pad);

    // --- Channel-score SKELETON of the window (the obvious channel, even where
    // DEM flow went down the road). We extract segments, keep those PROXIMATE to
    // the vector channel (this drops road-side ditches, which run with the road
    // away from the vector), bridge the ~50 m breaks, then route the connected
    // proximate channel from the divergence point to the junction. ---
    lsdtt_raster_cpu skel(
        lsdtt_array2d<float>(nr, nc, 0.0f),
        lsdtt_array2d<uint8_t>(dem.get_mask_array()), meta);
    {
        const std::int64_t nflat = static_cast<std::int64_t>(nr)*ncols64;
        std::vector<float> netsc;
        for (std::int64_t f=0; f<nflat; ++f){
            const std::size_t idx=static_cast<std::size_t>(f);
            if (mask[idx] && dem_network.is_channel(f) && channel_score_data[idx]>0.0f)
                netsc.push_back(channel_score_data[idx]);
        }
        float thr=0.0f;
        if(!netsc.empty()){ std::sort(netsc.begin(),netsc.end());
            thr=netsc[netsc.size()/2]; }   // median: keep the channel complete,
                                           // D8-bridging below handles the gaps
        lsdtt_raster_cpu score_bin(
            lsdtt_array2d<float>(nr, nc, 0.0f),
            lsdtt_array2d<uint8_t>(dem.get_mask_array()), meta);
        {
            auto& sb=score_bin.get_data();
            for(int r=wr_min;r<=wr_max;++r) for(int c=wc_min;c<=wc_max;++c){
                const std::size_t idx=static_cast<std::size_t>(to_flat(r,c));
                if(mask[idx] && channel_score_data[idx]>0.0f && channel_score_data[idx]>=thr)
                    sb[idx]=1.0f;
            }
        }
        lsdtt_raster_cpu filt=filter_binary_components_by_size(score_bin,3);
        filt=binary_fill_holes(filt);
        skel=thin_distance_guided(filt,1.0);
        skel=prune_skeleton_spurs(skel,2);
        skel=prune_terminal_loops(skel,30);
        std::cout << "    capture: skeleton score threshold=" << thr << "\n";
    }

    const std::vector<float>& elevd = dem.get_data();

    // Distance from (x,y) to the vector edge polyline.
    auto sep_to_poly=[&](double x,double y)->double{
        const auto& P=reach.vector_reach_xy;
        double best=std::numeric_limits<double>::max();
        for(std::size_t i=1;i<P.size();++i){
            const double d=point_to_segment_distance(x,y,P[i-1].x,P[i-1].y,P[i].x,P[i].y);
            if(d<best) best=d;
        }
        return best;
    };

    // Keep only segments PROXIMATE to the vector channel (drops the road ditch).
    const double prox_thresh_m = 50.0;
    lsdtt_raster_cpu prox(
        lsdtt_array2d<float>(nr, nc, 0.0f),
        lsdtt_array2d<uint8_t>(dem.get_mask_array()), meta);
    int n_seg=0, n_keep=0;
    {
        const std::vector<skeleton_segment> segs = extract_skeleton_segments(skel);
        n_seg=static_cast<int>(segs.size());
        auto& pd=prox.get_data();
        for(const auto& s: segs){
            if(s.pixels.empty()) continue;
            double sum=0.0; int n=0;
            for(const auto& px: s.pixels){
                const lsdtt_xy xy=cell_center_from_row_col(meta,px.row,px.col);
                sum+=sep_to_poly(xy.x,xy.y); ++n;
            }
            if(n==0 || sum/n > prox_thresh_m) continue;
            ++n_keep;
            for(const auto& px: s.pixels)
                pd[static_cast<std::size_t>(to_flat(px.row,px.col))]=1.0f;
        }
    }
    std::cout << "    capture: skeleton segments=" << n_seg
              << " kept(proximate to vector, <" << prox_thresh_m << " m)="
              << n_keep << "\n";

    // Bridge the ~50 m breaks among the proximate segments (tangent-based).
    lsdtt_raster_cpu conn = prox;
    try {
        std::vector<skeleton_segment> psegs = extract_skeleton_segments(prox);
        std::vector<skeleton_endpoint_tangent> tans =
            estimate_skeleton_endpoint_tangents(psegs, 6.0, 2);
        endpoint_pair_request req;
        req.max_gap_pixels = std::max(60.0, std::ceil(70.0/cs));   // ~70 m
        req.max_facing_angle_degrees = 70.0;        // relax: channels bend at gaps
        req.max_heading_misalignment_degrees = 80.0;
        req.reject_same_component = true;
        std::vector<skeleton_endpoint_pair_candidate> cands =
            find_skeleton_endpoint_pair_candidates(psegs, tans, req);
        std::vector<skeleton_endpoint_pair_candidate> acc =
            select_skeleton_endpoint_pairs(cands, tans.size(),
                                           endpoint_pair_selection_request());
        conn = rasterise_skeleton_endpoint_bridges(prox, acc).connected_skeleton;
        std::cout << "    capture: bridged " << acc.size() << " gap(s)\n";
    } catch (const std::exception& e){
        std::cerr << "    capture: bridging failed: " << e.what() << "\n";
    }

    if (!params.diagnostic_prefix.empty()){
        try {
            lsdtt_data_loader w;
            const std::string base = params.diagnostic_prefix + "_diverged_skeleton_edge"
                + std::to_string(reach.vector_edge_id);
            w.write_full(skel.get_data(), meta, base + ".tif");
            w.write_full(conn.get_data(), meta, base + "_proximate.tif");
        } catch (...) {}
    }

    // Connected proximate skeleton pixel set.
    std::unordered_set<std::int64_t> conn_set;
    {
        const auto& cd=conn.get_data();
        const std::int64_t nflat=static_cast<std::int64_t>(nr)*ncols64;
        for(std::int64_t f=0; f<nflat; ++f)
            if(cd[static_cast<std::size_t>(f)]>0.5f) conn_set.insert(f);
    }
    if (conn_set.empty()){
        std::cout << "    capture: no proximate skeleton in window\n"; return out; }

    // Augment connectivity with D8 flow (Simon's idea): where flow from one
    // proximate skeleton pixel reaches another within a gap budget WITHOUT
    // passing through the road, add the D8 path. This bridges the fragmented
    // skeleton along the real drainage; bridges that would run into the road are
    // rejected (those are the capture points the breach itself must cross).
    {
        const int maxgap = std::max(20, static_cast<int>(std::ceil(60.0/cs)));
        const std::vector<std::int64_t> seeds(conn_set.begin(), conn_set.end());
        std::vector<std::int64_t> add;
        for (std::int64_t f : seeds){
            std::int64_t cur=f; std::vector<std::int64_t> seg; int steps=0;
            bool bridged=false;
            while (steps++ < maxgap){
                if (flow.is_outlet(cur)) break;
                const std::int64_t nxt=flow.receiver_of(cur);
                if (nxt<0 || nxt==cur || road_set.count(nxt)) break;
                cur=nxt;
                if (conn_set.count(cur)){ bridged=true; break; }
                seg.push_back(cur);
            }
            if (bridged) add.insert(add.end(), seg.begin(), seg.end());
        }
        for (std::int64_t p : add) conn_set.insert(p);
        std::cout << "    capture: connected proximate skeleton pixels="
                  << conn_set.size() << " (+ " << add.size() << " D8-bridge px)\n";
    }

    // Source = upstream end of the reach (the dem_jid_up channel head). The
    // repair follows the D8 reach from the source and recarves only the divergent
    // span (below) along the valley skeleton.
    std::int64_t src_flat = start_flat;
    if (!reach.dem_reach_xy.empty()){
        int rr,cc; rc_of_xy(reach.dem_reach_xy.front(), rr, cc);
        if (valid(rr,cc)) src_flat = to_flat(rr,cc);
    }

    // Build the DEM reach (the D8 walk from the source) as flat indices, and find
    // the contiguous DIVERGENT span - the part whose lateral separation from the
    // vector exceeds the corridor (the road-capture detour). Simon's model: FOLLOW
    // the D8 reach (it already tracks the vector and ends ~15 m from the junction);
    // recarve ONLY this divergent span along the true valley. Everything outside the
    // span is kept verbatim from the D8 reach.
    std::vector<std::int64_t> reach_flats; reach_flats.reserve(reach.dem_reach_xy.size());
    for(const auto& p: reach.dem_reach_xy){ int r,c; rc_of_xy(p,r,c);
        reach_flats.push_back(valid(r,c)?to_flat(r,c):-1); }
    const double sep_thresh = corridor_half_width_m;     // outside the corridor
    int i_div0=-1, i_div1=-1;
    for(std::size_t i=0;i<reach.dem_reach_xy.size();++i){
        const auto& p=reach.dem_reach_xy[i];
        if(sep_to_poly(p.x,p.y) > sep_thresh){ if(i_div0<0) i_div0=(int)i; i_div1=(int)i; }
    }
    // Span endpoints on the reach: last good pixel BEFORE the detour and first good
    // pixel AFTER it (these sit in the valley, near the vector).
    auto valid_reach_idx=[&](int i)->std::int64_t{
        while(i>=0 && i<(int)reach_flats.size() && reach_flats[i]<0) --i;
        return (i>=0 && i<(int)reach_flats.size())? reach_flats[i] : -1; };
    auto valid_reach_idx_fwd=[&](int i)->std::int64_t{
        while(i>=0 && i<(int)reach_flats.size() && reach_flats[i]<0) ++i;
        return (i>=0 && i<(int)reach_flats.size())? reach_flats[i] : -1; };
    int entry_idx = (i_div0<0)? -1 : std::max(0, i_div0-1);
    int exit_idx  = (i_div1<0)? -1 : std::min((int)reach_flats.size()-1, i_div1+1);
    const std::int64_t span_entry = (entry_idx<0)? src_flat : valid_reach_idx(entry_idx);
    const std::int64_t span_exit  = (exit_idx<0)? -1        : valid_reach_idx_fwd(exit_idx);
    std::cout << "    capture: reach n="<<reach_flats.size()
              <<" divergent span idx ["<<i_div0<<","<<i_div1<<"] sep_thresh="
              <<sep_thresh<<" m\n";

    // head = conn pixel nearest the span ENTRY (where the reach leaves the valley);
    // goalc = conn pixel nearest the span EXIT (where it rejoins). The BFS then
    // routes ONLY across the detour, along the skeleton (true valley).
    const int er=to_row(span_entry), ec=to_col(span_entry);
    const int xr=(span_exit<0?gr0:to_row(span_exit)), xc=(span_exit<0?gc0:to_col(span_exit));
    std::int64_t head=-1, goalc=-1; double hb=1e18, gb=1e18;
    for(std::int64_t f: conn_set){
        const double dh=std::hypot(double(to_row(f)-er),double(to_col(f)-ec));
        if(dh<hb){hb=dh;head=f;}
        const double dg=std::hypot(double(to_row(f)-xr),double(to_col(f)-xc));
        if(dg<gb){gb=dg;goalc=f;}
    }
    const double max_head_px=std::max(double(corr_px)*3.0,120.0);
    if(head<0 || hb>max_head_px){
        std::cout << "    capture: no proximate skeleton within "
                  << max_head_px << " px of the span entry (nearest="
                  << hb << ")\n"; return out; }

    // BFS along the connected proximate skeleton from head to goalc (robust to
    // the meanders that broke the steepest-descent walk).
    std::vector<std::int64_t> bfspath;
    {
        std::unordered_map<std::int64_t,std::int64_t> par;
        std::vector<std::int64_t> q; q.push_back(head); par[head]=head;
        std::size_t qh=0; bool found=false;
        while(qh<q.size()){
            const std::int64_t cur=q[qh++];
            if(cur==goalc){found=true;break;}
            const int cr=to_row(cur),cc=to_col(cur);
            for(int dr=-1;dr<=1;++dr)for(int dc=-1;dc<=1;++dc){
                if(!dr&&!dc) continue;
                const int r=cr+dr,c=cc+dc; if(!valid(r,c))continue;
                const std::int64_t nf=to_flat(r,c);
                if(!conn_set.count(nf)||par.count(nf)) continue;
                par[nf]=cur; q.push_back(nf);
            }
        }
        if(!found){
            std::cout << "    capture: proximate skeleton not connected "
                         "across the divergent span; leaving unresolved\n";
            return out;
        }
        std::int64_t cur=goalc;
        while(cur!=head){ bfspath.push_back(cur); cur=par[cur]; }
        bfspath.push_back(head);
        std::reverse(bfspath.begin(), bfspath.end());
    }

    auto bres = [&](std::int64_t f0,std::int64_t f1,std::vector<std::int64_t>& p){
        int r0=to_row(f0),c0=to_col(f0),r1=to_row(f1),c1=to_col(f1);
        const int dr=std::abs(r1-r0),dc=std::abs(c1-c0);
        const int srr=r0<r1?1:-1,scc=c0<c1?1:-1; int err=dc-dr,r=r0,c=c0;
        while(true){
            if(valid(r,c)){ const std::int64_t f=to_flat(r,c);
                if(p.empty()||p.back()!=f) p.push_back(f); }
            if(r==r1&&c==c1) break;
            const int e2=2*err;
            if(e2>-dr){err-=dr;c+=scc;}
            if(e2<dc){err+=dc;r+=srr;}
        } };
    // Carve the divergent span along the valley skeleton: start at the span ENTRY
    // (on the D8 reach), follow the BFS skeleton path to goalc, and finish at the
    // span EXIT (back on the D8 reach). This recarves ONLY the road-capture detour,
    // along the true valley.
    std::vector<std::int64_t> carve;
    {
        bres(span_entry, head, carve);            // span entry -> skeleton
        for(std::int64_t f: bfspath)              // along the valley skeleton
            if(carve.empty()||carve.back()!=f) carve.push_back(f);
        if(span_exit>=0) bres(goalc, span_exit, carve);  // skeleton -> span exit
    }

    // Assemble the final route: kept upstream D8 reach [0..entry_idx] + the valley
    // carve + kept downstream D8 reach [exit_idx..end]. The downstream reach already
    // tracks the vector and ends at the trunk confluence ~15 m from the vector
    // junction, so no synthetic stop is needed - we keep the real D8 channel.
    std::vector<std::int64_t> route;
    auto push=[&](std::int64_t f){
        if(f>=0 && (route.empty()||route.back()!=f)) route.push_back(f); };
    if(entry_idx>=0){ for(int i=0;i<=entry_idx;++i) push(reach_flats[i]); }
    else            { push(src_flat); }
    for(std::int64_t f: carve) push(f);
    if(exit_idx>=0){ for(int i=exit_idx;i<(int)reach_flats.size();++i) push(reach_flats[i]); }
    if (route.size()<2){ std::cout<<"    capture: degenerate route\n"; return out; }

    // Burn the WHOLE route as one capturing path. burn_capturing_path runs a
    // running-minimum from the higher (upstream) end, so the result is globally
    // monotone-descending - it cuts the road embankment and the small noise lips
    // and daylights at the downstream channel (which the profile confirms sits
    // well below the upstream side). fill_pits in the reflow mops up any residual.
    try { burn_capturing_path(dem, route, std::max(params.repair_min_slope,1.0e-4)); }
    catch (const std::exception& e){
        std::cerr << "    capture: burn failed: " << e.what() << "\n"; return out; }

    out.burned_path   = route;
    out.n_burned      = static_cast<int>(route.size());
    out.fix_type      = GapFixType::channel_guided_bridge;
    out.divergence_xy = cell_center_from_row_col(meta, sr, sc);
    out.target_xy     = cell_center_from_row_col(meta, to_row(route.back()), to_col(route.back()));
    std::cout << "    capture: burned whole route " << route.size()
              << " px (head " << hb << " px from source)\n";
    return out;
}

// ===========================================================================
// repair_diverged_reaches
//
// Reach-driven repair for the road-capture-then-rejoin case the junction-LCA
// matcher misses. A diverged reach is repaired with a corridor capturing burn
// (repair_diverged_reach_capture) anchored at its divergence point; a
// topology_break is reconnected at its downstream DEM confluence
// (repair_junction_parallel_channels). The anchor confluence is
// re-matched from the downstream vector junction against the CURRENT network so
// it stays correct across reflows.
// ===========================================================================
std::vector<JunctionRepairRecord>
repair_diverged_reaches(
    lsdtt_raster_cpu&                             dem,
    const lsdtt_vector_to_junction_network&       vector_network,
    const std::string&                            vector_path,
    const std::vector<ReachDivergenceRecord>&     reaches,
    const std::vector<std::int64_t>&              source_flats,
    std::unique_ptr<lsdtt_flow_base>&             flow,
    std::unique_ptr<lsdtt_junction_network_base>& dem_network,
    const lsdtt_flowinfo_options&                 fi_opts,
    const lsdtt_junction_network_options&         jn_opts,
    const junction_repair_params&                 params_in,
    int                                           phase,
    int                                           repair_index_offset,
    const std::vector<float>*                     precomputed_channel_scores)
{
    std::vector<JunctionRepairRecord> records;
    if (reaches.empty()) return records;

    junction_repair_params params = params_in;
    params.diagnostic_phase = phase;

    const lsdtt_raster_meta& meta = dem.get_meta();

    // Channel scores are required for the narrows-finding repair. Reuse the
    // caller's precomputed scores (same as repair_candidate_junctions); without
    // them the parallel-channel / channel-guided repairs cannot run.
    std::vector<float> channel_score_data;
    bool have_channel_scores = false;
    if (precomputed_channel_scores && !precomputed_channel_scores->empty())
    {
        channel_score_data  = *precomputed_channel_scores;
        have_channel_scores = true;
        std::cout << "repair_diverged_reaches: using precomputed channel scores\n";
    }
    else
    {
        std::cerr << "repair_diverged_reaches: no channel scores supplied; "
                  << "skipping divergence repairs\n";
        return records;
    }

    const auto& vjns = vector_network.junctions();
    const std::int32_t n_vjn = vector_network.n_junctions();

    // Window scale must not run away on topology_break reaches, whose channel
    // diverts to a different basin so max_separation_m can be kilometres. Cap
    // it: the narrows we want is the local col near the confluence.
    const double max_gap_scale_m = 200.0;

    // Dedup: a topology_break is keyed by its downstream confluence (one repair
    // per confluence); a diverged reach is keyed by its divergence point.
    std::unordered_set<int>          handled_break_vjid;
    std::unordered_set<std::int64_t> handled_div_points;

    // Corridor half-width for the diverged capturing burn: narrower than the
    // road's separation so the road is excluded, but wide enough to span the
    // vector's own 10-20 m offset to the true thalweg.
    const double corridor_half_width_m = 20.0;

    int repair_index = repair_index_offset;

    for (const ReachDivergenceRecord& r : reaches)
    {
        const bool is_div = (r.status == ReachStatus::diverged);
        const bool is_brk = (r.status == ReachStatus::topology_break);
        if (!is_div && !is_brk) continue;
        if (r.vector_jid_dn < 0 || r.vector_jid_dn >= n_vjn) continue;

        if (is_brk)
        {
            if (!handled_break_vjid.insert(r.vector_jid_dn).second) continue;
        }
        else
        {
            if (!r.has_divergence_point) continue;
            if (!handled_div_points.insert(r.divergence_point_flat).second) continue;
        }

        const auto& vj_dn = vjns[static_cast<std::size_t>(r.vector_jid_dn)];
        const lsdtt_xy vxy_dn{vj_dn.x, vj_dn.y};

        JunctionRepairRecord rec;
        rec.repair_index       = repair_index++;
        rec.phase              = phase;
        rec.vector_junction_id = r.vector_jid_dn;
        rec.vector_xy          = is_div ? r.divergence_point_xy : vxy_dn;
        rec.distance_m         = r.max_separation_m;
        rec.mismatch_reason    = std::string("reach_") + reach_status_str(r.status);
        rec.fix_type           = GapFixType::unresolved;

        // Snapshot the DEM so a burn whose reflow fails (leaves an unresolved
        // sink) can be rolled back instead of aborting the run.
        lsdtt_array2d_f dem_backup = dem.get_data_array();
        repair_attempt_result attempt;

        if (is_div)
        {
            // Diverged-but-rejoining: follow the D8 reach, recarve only the divergent
            // span along the valley skeleton (repair_diverged_reach_capture).
            std::cout << "  repairing diverged reach edge=" << r.vector_edge_id
                      << " v_jid_dn=" << r.vector_jid_dn
                      << " max_sep=" << r.max_separation_m
                      << " m capture@(" << r.divergence_point_xy.x << ","
                      << r.divergence_point_xy.y << ")\n";
            try
            {
                attempt = repair_diverged_reach_capture(
                    dem, *flow, *dem_network, r, vxy_dn,
                    corridor_half_width_m, channel_score_data, params);
            }
            catch (const std::exception& e)
            {
                std::cerr << "    capture exception: " << e.what() << "\n";
                attempt.fix_type = GapFixType::unresolved;
            }
        }
        else
        {
            // topology_break: reconnect at the downstream DEM confluence,
            // re-matched against the CURRENT network (reflow may have renumbered
            // junctions since the reach diagnostic ran).
            const std::vector<dem_candidate_junction_view> dem_views =
                build_dem_candidate_junction_views(*dem_network, dem);
            double anchor_dist = 0.0;
            const int anchor_jid =
                find_nearest_dem_junction(vxy_dn, dem_views, anchor_dist);
            rec.dem_junction_id = anchor_jid;

            if (anchor_jid < 0)
            {
                std::cout << "  topology_break reach edge=" << r.vector_edge_id
                          << ": no DEM confluence near downstream vector junction; "
                          << "skipping\n";
                records.push_back(rec);
                continue;
            }
            rec.dem_xy = cell_center_from_row_col(
                meta,
                dem_network->row_of_junction(static_cast<std::int32_t>(anchor_jid)),
                dem_network->col_of_junction(static_cast<std::int32_t>(anchor_jid)));
            const std::int64_t dem_junc_flat =
                dem_network->flat_of_junction(static_cast<std::int32_t>(anchor_jid));
            const double gap_scale = std::clamp(
                r.max_separation_m, params.min_gap_m, max_gap_scale_m);

            std::cout << "  repairing topology_break reach edge=" << r.vector_edge_id
                      << " v_jid_dn=" << r.vector_jid_dn
                      << " anchor_dem_jid=" << anchor_jid
                      << " max_sep=" << r.max_separation_m << " m\n";
            try
            {
                attempt = repair_junction_parallel_channels(
                    dem, *flow, *dem_network,
                    dem_junc_flat, vxy_dn,
                    r.vector_jid_dn, vector_network,
                    gap_scale, params.min_gap_m,
                    channel_score_data, params);
            }
            catch (const std::exception& e)
            {
                std::cerr << "    parallel channel repair exception: "
                          << e.what() << "\n";
                attempt.fix_type = GapFixType::unresolved;
            }
        }

        rec.fix_type         = attempt.fix_type;
        rec.n_burned         = attempt.n_burned;
        rec.max_climb_m      = attempt.max_climb_m;
        rec.total_cost       = attempt.total_cost;
        rec.n_landing_pixels = attempt.n_landing_pixels;
        rec.divergence_xy    = attempt.divergence_xy;
        rec.target_xy        = attempt.target_xy;
        rec.burned_path      = std::move(attempt.burned_path);

        if (attempt.fix_type != GapFixType::unresolved &&
            attempt.fix_type != GapFixType::ok)
        {
            // make_flow_base is OpenMP and TERMINATES the process on an
            // unresolved sink (even single-threaded), so we cannot catch our way
            // out. FILL any residual pit the burn left BEFORE reflowing - a
            // capturing burn that daylights is unchanged by the fill; one that
            // accidentally pitted is resolved instead of crashing.
            lsdtt_flowinfo_options fi_single = fi_opts;
            fi_single.n_threads = 1;
            const float fill_slope = (params.repair_min_slope > 0.0)
                ? static_cast<float>(params.repair_min_slope) : 1.0e-4f;

            bool reflow_ok = true;
            try
            {
                lsdtt_raster_cpu filled = fill_pits(dem, fill_slope, fi_opts.n_threads);
                dem.get_data_array() = filled.get_data_array();
                flow = make_flow_base(
                    dem.get_data_array(), dem.get_mask_array(), meta, fi_single);
                dem_network = make_junction_network_base(*flow, source_flats, jn_opts);
            }
            catch (const std::exception& e)
            {
                // Should not happen now (fill resolves sinks); roll back to be safe.
                reflow_ok = false;
                dem.get_data_array() = dem_backup;
                lsdtt_raster_cpu filled = fill_pits(dem, fill_slope, fi_opts.n_threads);
                dem.get_data_array() = filled.get_data_array();
                flow = make_flow_base(
                    dem.get_data_array(), dem.get_mask_array(), meta, fi_single);
                dem_network = make_junction_network_base(*flow, source_flats, jn_opts);
                rec.fix_type = GapFixType::unresolved;
                rec.success  = false;
                std::cerr << "    reverted divergence burn for edge="
                          << r.vector_edge_id << " (reflow failed: "
                          << e.what() << ")\n";
            }

            if (reflow_ok)
            {
                int nearest_j = -1;
                rec.post_distance_m = nearest_junction_distance_to_xy(
                    *dem_network, meta, vxy_dn, &nearest_j);
                rec.success = true;
                std::cout << "    result: " << gap_fix_type_str(rec.fix_type)
                          << " burned=" << rec.n_burned << "\n";
            }
        }
        else
        {
            rec.post_distance_m = r.max_separation_m;
            rec.success = false;
            std::cout << "    repair unresolved\n";
        }

        records.push_back(rec);
    }

    return records;
}

// ===========================================================================
// repair_stream_junctions  [DEPRECATED]
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

    //  Channel candidate scores (computed once on the unmodified DEM) 
    //
    // Channel score = multiscale plan-curvature persistence (Hessian + depth).
    // Range [0, 1]: 1 = strong valley / channel, 0 = ridge or planar terrain.
    //
    // Two roles in the repair algorithm:
    //   (a) Phase 1 divergence scoring - channel score contributes secondary
    //       preference for concave terrain when ranking upstream channel pixels
    //       as divergence-point candidates.
    //   (b) Phase 3 Dijkstra weighting - high-score pixels are cheaper to
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
                  << "-"        << cc_req.channel_lambda_max_metres
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
    // cache maps vector node_id -> DEM junction_id for the CURRENT extraction.
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

    //  Nearest-junction distance helper 
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

            std::cout << " - repairing\n";

            GapFixType fix_type = GapFixType::unresolved;
            int        n_burned = 0;
            bool       repair_ok = false;

            //  Attempt 1: three-phase channel-guided repair 
            if (have_channel_scores)
            {
                const std::int64_t dem_junc_flat =
                    dem_network->flat_of_junction(dem_junc);

                repair_attempt_result attempt = repair_junction_channel_guided(
                    dem, *flow, *dem_network,
                    dem_junc_flat,
                    junc_node.xy,
                    dist,
                    params.min_gap_m,
                    channel_score_data,
                    params);

                fix_type = attempt.fix_type;
                n_burned = attempt.n_burned;

                if (fix_type == GapFixType::channel_guided_bridge)
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
                        std::cout << "    -> channel-guided: success"
                                  << "  new gap=" << new_gap << " m\n";
                    }
                    else
                    {
                        std::cout << "    -> channel-guided: burned but gap still "
                                  << new_gap << " m; trying terrain Dijkstra fallback\n";
                    }
                }
                else
                {
                    std::cout << "    -> channel-guided: " << gap_fix_type_str(fix_type)
                              << "; trying terrain Dijkstra fallback\n";
                }
            }

            //  Attempt 2: terrain-only Dijkstra fallback 
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

                std::cout << "    -> terrain fallback: " << gap_fix_type_str(fix_type)
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

    f << "repair_index,phase,basin_outlet_junction,vector_junction_id,"
         "dem_junction_id,vector_x,vector_y,dem_x,dem_y,"
         "distance_m,post_distance_m,fix_type,n_burned,success,"
         "max_climb_m,total_cost,n_landing_pixels,"
         "divergence_x,divergence_y,mismatch_reason\n";

    for (const auto& r : records)
    {
        f << r.repair_index           << ","
          << r.phase                  << ","
          << r.basin_outlet_junction  << ","
          << r.vector_junction_id     << ","
          << r.dem_junction_id        << ","
          << r.vector_xy.x            << ","
          << r.vector_xy.y            << ","
          << r.dem_xy.x               << ","
          << r.dem_xy.y               << ","
          << r.distance_m             << ","
          << r.post_distance_m        << ","
          << gap_fix_type_str(r.fix_type) << ","
          << r.n_burned               << ","
          << (r.success ? 1 : 0)      << ","
          << r.max_climb_m            << ","
          << r.total_cost             << ","
          << r.n_landing_pixels       << ","
          << r.divergence_xy.x        << ","
          << r.divergence_xy.y        << ","
          << r.mismatch_reason        << "\n";
    }
}

} // namespace lsdtt
