// lsdtt_vector_line_repair_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_vector_line_repair_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <variant>

#include "lsdtt_vector_ops_cpu.hpp"
#include "../utils/lsdtt_linear_algebra.hpp"
#include "../utils/lsdtt_polygon_geometry.hpp"

namespace lsdtt
{
namespace
{

double dist_xy(const lsdtt_xy& a, const lsdtt_xy& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Outward unit tangent at a polyline tip, estimated over up to fit_length of
// inward arc length. Outward = pointing away from the line body (the direction
// the line would extend if continued past the tip).
void tip_outward_tangent(
    const std::vector<lsdtt_xy>& pts,
    bool                         at_start,
    double                       fit_length,
    double&                      tx,
    double&                      ty)
{
    tx = 0.0;
    ty = 0.0;
    const std::size_t n = pts.size();
    if (n < 2) return;

    const lsdtt_xy tip = at_start ? pts.front() : pts.back();
    lsdtt_xy inward    = at_start ? pts[1] : pts[n - 2];

    double acc = 0.0;
    if (at_start)
    {
        for (std::size_t i = 1; i < n; ++i)
        {
            acc += dist_xy(pts[i - 1], pts[i]);
            inward = pts[i];
            if (acc >= fit_length) break;
        }
    }
    else
    {
        for (std::size_t i = n - 1; i > 0; --i)
        {
            acc += dist_xy(pts[i], pts[i - 1]);
            inward = pts[i - 1];
            if (acc >= fit_length) break;
        }
    }

    const double dx = tip.x - inward.x;
    const double dy = tip.y - inward.y;
    const double len = std::sqrt(dx * dx + dy * dy);
    if (len > 0.0)
    {
        tx = dx / len;
        ty = dy / len;
    }
}

// Flat road segment (world coords) tagged with the road record it came from.
struct road_segment
{
    double      x1, y1, x2, y2;
    std::size_t road_index;
};

std::vector<road_segment> build_road_segments(const lsdtt_vector_layer& roads)
{
    std::vector<road_segment> segs;
    for (std::size_t ri = 0; ri < roads.records.size(); ++ri)
    {
        const auto* line =
            std::get_if<lsdtt_line_geom>(&roads.records[ri].geometry);
        if (!line || line->points.size() < 2) continue;
        for (std::size_t k = 1; k < line->points.size(); ++k)
        {
            const lsdtt_xy& a = line->points[k - 1];
            const lsdtt_xy& b = line->points[k];
            segs.push_back(road_segment{a.x, a.y, b.x, b.y, ri});
        }
    }
    return segs;
}

std::string read_string_field(
    const lsdtt_vector_layer&  layer,
    const lsdtt_vector_record& rec,
    const std::string&         name)
{
    const int idx = layer.schema.index_of(name);
    if (idx < 0 || idx >= static_cast<int>(rec.field_values.size())) return "";
    if (const auto* s =
            std::get_if<std::string>(&rec.field_values[static_cast<std::size_t>(idx)]))
        return *s;
    return "";
}

std::int64_t read_int_field(
    const lsdtt_vector_layer&  layer,
    const lsdtt_vector_record& rec,
    const std::string&         name)
{
    const int idx = layer.schema.index_of(name);
    if (idx < 0 || idx >= static_cast<int>(rec.field_values.size())) return 0;
    if (const auto* v =
            std::get_if<std::int64_t>(&rec.field_values[static_cast<std::size_t>(idx)]))
        return *v;
    return 0;
}

std::string road_name_for(
    const lsdtt_vector_layer& roads,
    std::size_t               road_index,
    const std::string&        primary_field)
{
    if (road_index >= roads.records.size()) return "";
    const lsdtt_vector_record& rec = roads.records[road_index];
    std::string name = read_string_field(roads, rec, primary_field);
    if (name.empty()) name = read_string_field(roads, rec, "roadNumber");
    return name;
}

} // anonymous namespace


// ---------------------------------------------------------------------------
// extract_water_line_tips
// ---------------------------------------------------------------------------
std::vector<water_line_tip> extract_water_line_tips(
    const lsdtt_vector_layer&        water,
    const water_line_repair_options& opt)
{
    const lsdtt_endpoint_graph graph =
        build_line_endpoint_graph(water, opt.snap_tolerance_metres);

    std::vector<water_line_tip> tips;
    for (const lsdtt_line_endpoint& ep : graph.endpoints)
    {
        if (ep.degree != 1) continue;   // only dangling tips

        const auto* line =
            std::get_if<lsdtt_line_geom>(&water.records[ep.line_index].geometry);
        if (!line || line->points.size() < 2) continue;

        water_line_tip tip;
        tip.line_index = ep.line_index;
        tip.at_start   = ep.at_start;
        tip.xy         = ep.xy;
        tip_outward_tangent(line->points, ep.at_start,
                            opt.tangent_fit_length_metres, tip.tx, tip.ty);
        tips.push_back(tip);
    }
    return tips;
}

// ---------------------------------------------------------------------------
// find_water_bridge_candidates
// ---------------------------------------------------------------------------
std::vector<water_bridge_candidate> find_water_bridge_candidates(
    const std::vector<water_line_tip>& tips,
    const lsdtt_vector_layer&          roads,
    const water_line_repair_options&   opt)
{
    std::vector<water_bridge_candidate> out;

    const std::vector<road_segment> rsegs = build_road_segments(roads);

    const double max_facing  = std::max(opt.max_facing_angle_degrees, 1e-9);
    const double max_heading = std::max(opt.max_heading_misalignment_degrees, 1e-9);
    const double road_tol    = std::max(opt.road_intersection_tolerance_metres, 0.0);

    // A pair is bridged if it crosses a road within max_gap_length_metres, OR is
    // a short gap (<= max_nonroad_gap_metres) regardless of roads. The overall
    // distance cap is the larger of the two; the score normaliser uses it too.
    const double gap_cap = std::max({opt.max_gap_length_metres,
                                     opt.max_nonroad_gap_metres, 1e-9});

    for (std::size_t i = 0; i < tips.size(); ++i)
    {
        for (std::size_t j = i + 1; j < tips.size(); ++j)
        {
            const water_line_tip& ta = tips[i];
            const water_line_tip& tb = tips[j];

            // Never bridge a single line's two ends to itself.
            if (ta.line_index == tb.line_index) continue;

            const double gap = dist_xy(ta.xy, tb.xy);
            if (gap <= 0.0 || gap > gap_cap) continue;

            const double dirx = tb.xy.x - ta.xy.x;
            const double diry = tb.xy.y - ta.xy.y;

            const double fa = angle_between_vectors_degrees(ta.tx, ta.ty, dirx, diry);
            const double fb = angle_between_vectors_degrees(tb.tx, tb.ty, -dirx, -diry);
            const double heading =
                angle_between_vectors_degrees(ta.tx, ta.ty, tb.tx, tb.ty);
            if (!(std::isfinite(fa) && std::isfinite(fb) && std::isfinite(heading)))
                continue;
            const double misalign = std::abs(180.0 - heading);

            // Road association: connector within road_tol of a road centreline
            // (strict crossing has distance 0; a staggered severance running
            // ~parallel to the road stays within tolerance).
            bool        near_road = false;
            std::size_t road_idx  = static_cast<std::size_t>(-1);
            double      ix = 0.0, iy = 0.0;
            for (const road_segment& rs : rsegs)
            {
                if (segment_segment_intersection(
                        ta.xy.x, ta.xy.y, tb.xy.x, tb.xy.y,
                        rs.x1, rs.y1, rs.x2, rs.y2, ix, iy))
                {
                    near_road = true; road_idx = rs.road_index; break;
                }
            }
            if (!near_road && road_tol > 0.0)
            {
                for (const road_segment& rs : rsegs)
                {
                    const double d1 = point_to_segment_distance(
                        ta.xy.x, ta.xy.y, rs.x1, rs.y1, rs.x2, rs.y2);
                    const double d2 = point_to_segment_distance(
                        tb.xy.x, tb.xy.y, rs.x1, rs.y1, rs.x2, rs.y2);
                    const double d3 = point_to_segment_distance(
                        rs.x1, rs.y1, ta.xy.x, ta.xy.y, tb.xy.x, tb.xy.y);
                    const double d4 = point_to_segment_distance(
                        rs.x2, rs.y2, ta.xy.x, ta.xy.y, tb.xy.x, tb.xy.y);
                    const double dmin =
                        std::min(std::min(d1, d2), std::min(d3, d4));
                    if (dmin <= road_tol)
                    {
                        near_road = true; road_idx = rs.road_index; break;
                    }
                }
            }

            // Tiered acceptance:
            //  - ROAD tier: near a road, within max_gap_length, AT LEAST ONE tip
            //    facing (the channel kinks at the crossing so the dangling tip's
            //    own tangent is often unreliable, but the partner faces well).
            //    No heading constraint.
            //  - NON-ROAD tier: no road excuse, so require BOTH tips facing,
            //    tangents ~opposite, and a short gap (max_nonroad_gap).
            const bool strict_angles =
                (fa <= opt.max_facing_angle_degrees &&
                 fb <= opt.max_facing_angle_degrees &&
                 misalign <= opt.max_heading_misalignment_degrees);
            const bool road_angle_ok =
                std::min(fa, fb) <= opt.road_max_facing_angle_degrees;

            bool accept  = false;
            bool is_road = false;
            if (near_road && gap <= opt.max_gap_length_metres && road_angle_ok)
            { accept = true; is_road = true; }
            else if (gap <= opt.max_nonroad_gap_metres && strict_angles)
            { accept = true; is_road = false; }
            if (!accept) continue;

            water_bridge_candidate c;
            c.tip_index_a = i;
            c.tip_index_b = j;
            c.gap_length_metres = gap;
            c.facing_angle_a_degrees = fa;
            c.facing_angle_b_degrees = fb;
            c.heading_misalignment_degrees = misalign;
            c.crosses_road = is_road;
            c.crossed_road_index = is_road ? road_idx
                                           : static_cast<std::size_t>(-1);
            c.crossed_road_name  = is_road
                ? road_name_for(roads, road_idx, opt.road_name_field)
                : std::string();
            c.score =
                opt.weight_distance * (gap / gap_cap) +
                opt.weight_facing   * (0.5 * (fa + fb) / max_facing) +
                opt.weight_heading  * (misalign / max_heading);

            out.push_back(std::move(c));
        }
    }

    return out;
}

// ---------------------------------------------------------------------------
// select_water_bridges
// ---------------------------------------------------------------------------
std::vector<water_bridge_candidate> select_water_bridges(
    const std::vector<water_bridge_candidate>& candidates,
    std::size_t                                n_tips,
    const water_line_repair_options&           opt)
{
    std::vector<water_bridge_candidate> sorted = candidates;
    std::sort(sorted.begin(), sorted.end(),
              [](const water_bridge_candidate& a, const water_bridge_candidate& b)
              { return a.score < b.score; });

    // Best (minimum) score seen at each tip, for the mutual-best test.
    std::vector<double> best_score(
        n_tips, std::numeric_limits<double>::infinity());
    for (const auto& c : sorted)
    {
        best_score[c.tip_index_a] = std::min(best_score[c.tip_index_a], c.score);
        best_score[c.tip_index_b] = std::min(best_score[c.tip_index_b], c.score);
    }

    std::vector<char> used(n_tips, 0);
    std::vector<water_bridge_candidate> accepted;

    for (const auto& c : sorted)
    {
        if (c.score > opt.max_score) continue;
        if (used[c.tip_index_a] || used[c.tip_index_b]) continue;
        if (opt.require_mutual_best &&
            !(c.score == best_score[c.tip_index_a] &&
              c.score == best_score[c.tip_index_b]))
            continue;

        accepted.push_back(c);
        used[c.tip_index_a] = 1;
        used[c.tip_index_b] = 1;
    }

    return accepted;
}

// ---------------------------------------------------------------------------
// build_repaired_water_layer  (keep originals + add bridge lines)
// ---------------------------------------------------------------------------
lsdtt_vector_layer build_repaired_water_layer(
    const lsdtt_vector_layer&                  water,
    const std::vector<water_line_tip>&         tips,
    const std::vector<water_bridge_candidate>& bridges,
    const water_line_repair_options&           opt)
{
    if (opt.merge_chains)
    {
        throw std::runtime_error(
            "build_repaired_water_layer: merge_chains mode is not yet "
            "implemented; use keep+bridge (merge_chains = false). Coincident "
            "endpoints already make the network topologically continuous for "
            "lsdtt_vector_to_junction_network.");
    }

    lsdtt_vector_layer out;
    out.schema.field_names = {
        "repair_type", "bridge_length_m", "crossed_road_name",
        "src_fid", "src_featcode"
    };
    out.schema.field_types = {
        lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::int64_type
    };

    // Originals (carry provenance).
    for (const lsdtt_vector_record& rec : water.records)
    {
        const auto* line = std::get_if<lsdtt_line_geom>(&rec.geometry);
        if (!line || line->points.size() < 2) continue;

        lsdtt_vector_record out_rec;
        out_rec.geometry = lsdtt_line_geom{line->points};
        out_rec.field_values = {
            std::string("original"),
            std::monostate{},                              // bridge_length_m NULL
            std::monostate{},                              // crossed_road_name NULL
            read_string_field(water, rec, "Fid"),
            read_int_field(water, rec, "FeatCode")
        };
        out.records.push_back(std::move(out_rec));
    }

    // Bridges.
    for (const water_bridge_candidate& b : bridges)
    {
        const water_line_tip& ta = tips[b.tip_index_a];
        const water_line_tip& tb = tips[b.tip_index_b];

        lsdtt_vector_record out_rec;
        out_rec.geometry = lsdtt_line_geom{ std::vector<lsdtt_xy>{ ta.xy, tb.xy } };
        out_rec.field_values = {
            std::string(b.crosses_road ? "bridge_road" : "bridge_gap"),
            b.gap_length_metres,
            b.crossed_road_name.empty()
                ? lsdtt_vector_field_value{std::monostate{}}
                : lsdtt_vector_field_value{b.crossed_road_name},
            std::string(""),
            std::monostate{}                               // src_featcode NULL
        };
        out.records.push_back(std::move(out_rec));
    }

    return out;
}

// ---------------------------------------------------------------------------
// repair_water_lines
// ---------------------------------------------------------------------------
water_line_repair_result repair_water_lines(
    const lsdtt_vector_layer&        water,
    const lsdtt_vector_layer&        roads,
    const water_line_repair_options& opt)
{
    water_line_repair_result result;

    result.tips = extract_water_line_tips(water, opt);
    result.n_dangling_tips = static_cast<int>(result.tips.size());

    const std::vector<water_bridge_candidate> candidates =
        find_water_bridge_candidates(result.tips, roads, opt);
    result.n_candidate_pairs = static_cast<int>(candidates.size());

    result.accepted_bridges =
        select_water_bridges(candidates, result.tips.size(), opt);
    result.n_bridges_accepted = static_cast<int>(result.accepted_bridges.size());

    result.output_layer = build_repaired_water_layer(
        water, result.tips, result.accepted_bridges, opt);

    return result;
}

} // namespace lsdtt
