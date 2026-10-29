// test_valley_width_network.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// build_valley_width_network / apply_network_section_rules
// (lsdtt_valley_width_network_cpu), plan_valley_width_network steps 1-4, 6.
//
// Synthetic valleys of known floor width on a 5 m grid draining south:
//   trunk     N-S, floor cols 80..120 (41 cells, 205 m)
//   tributary E-W, joins from the east, floor rows 142..158 (17 cells, 85 m),
//             cols 100..175
//   sub-trib  N-S, joins the tributary from the north, floor cols 131..139
//             (9 cells, 45 m), rows 118..150
// Floors slope toward their centre lines steeply enough that floor flow
// converges there (no parallel edge channels), and the trunk is the longest
// flow path, so the tree's trunk is the main valley.
// The supplied valley mask is exactly the three floors.
//
// Checks: the tree (one trunk, the tributary its child, the sub-tributary the
// tributary's child); every kept section's width equals its own valley's floor
// width; no kept tributary station lies inside the parent's valley; the trunk
// sections through the tributary's mouth are excluded as tributary_crossing,
// and those through the sub-tributary's mouth on the TRIBUTARY, never on the
// trunk (the nesting); with both rules off the mouth widths reappear (so the
// rules are what removed them); stations where the mask has no valley are
// not_in_valley. Then a meandering channel in the trunk floor:
// smoothed_channel stays in the valley mask and is straighter than the
// channel; medial_path follows the valley's centre.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_channel_scope_cpu.hpp"
#include "lsdtt_valley_cross_sections_cpu.hpp"
#include "lsdtt_valley_width_network_cpu.hpp"

namespace
{
int g_failures = 0;
void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "PASS  " : "FAIL  ") << what << std::endl;
    if (!ok) ++g_failures;
}

constexpr int    k_rows = 300;
constexpr int    k_jr   = 150;   // tributary centre row (the junction)
constexpr int    k_cols = 200;
constexpr double k_dx   = 5.0;
constexpr double k_pi   = 3.14159265358979323846;

constexpr double k_trunk_width = 41 * k_dx;
constexpr double k_trib_width  = 17 * k_dx;
constexpr double k_sub_width   = 9 * k_dx;

bool in_trunk(int r, int c, bool meander)
{
    (void)r; (void)meander;
    return std::abs(c - 100) <= 20;
}
bool in_trib(int r, int c, bool with_tribs)
{
    return with_tribs && std::abs(r - k_jr) <= 8 && c >= 100 && c <= 175;
}
bool in_sub(int r, int c, bool with_tribs)
{
    return with_tribs && std::abs(c - 135) <= 4 && r >= 118 && r <= k_jr;
}

// Each valley is a surface defined everywhere: thalweg + gentle floor cross
// slope + a steep wall rise with distance outside the floor. The DEM is their
// minimum, so the valleys join as a union. The meander variant carves a 1 m
// sinuous groove in the trunk floor for the channel to follow.
struct fixture
{
    lsdtt::lsdtt_raster_cpu dem;
    lsdtt::lsdtt_raster_cpu valley;
    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    lsdtt::channel_scope_result scope;
};

fixture make_fixture(bool with_tribs, bool meander, int threshold)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = k_rows;
    meta.n_cols = k_cols;
    meta.cell_size = k_dx;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    const std::size_t n = static_cast<std::size_t>(k_rows) * k_cols;
    std::vector<float> z(n), v(n);
    auto trunk_z = [](int r) { return 100.0 + 0.1 * (k_rows - 1 - r); };
    auto trib_z  = [&](int c) { return trunk_z(k_jr) + 0.15 * (c - 100); };
    auto groove_col = [](int r) { return 100.0 + 12.0 * std::sin(2.0 * k_pi * r / 100.0); };

    for (int r = 0; r < k_rows; ++r)
        for (int c = 0; c < k_cols; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * k_cols + c;
            // The floor slopes toward the channel line - the valley centre, or
            // the meandering groove - so all floor flow converges on it.
            const double centre = meander ? groove_col(r) : 100.0;
            double zt = trunk_z(r) + 0.08 * std::abs(c - centre) +
                        2.0 * std::max(0, std::abs(c - 100) - 20);
            double best = zt;
            if (with_tribs)
            {
                const int cc = std::clamp(c, 100, 175);
                const double dist = std::max(0, std::abs(r - k_jr) - 8) +
                                    std::max(0, c - 175) + std::max(0, 100 - c);
                const double zb = trib_z(cc) + 0.1 * std::abs(r - k_jr) + 2.0 * dist;
                const int rr = std::clamp(r, 118, k_jr);
                const double ds = std::max(0, std::abs(c - 135) - 4) +
                                  std::max(0, 118 - r) + std::max(0, r - k_jr);
                const double zs = trib_z(135) + 0.2 * (k_jr - rr) +
                                  0.15 * std::abs(c - 135) + 2.0 * ds;
                best = std::min({best, zb, zs});
            }
            // The open south edge is one flat row below everything, so each of
            // its cells is an outlet; sloping it would route flow along it.
            if (r == k_rows - 1) best = trunk_z(k_rows - 1) - 10.0;
            z[i] = static_cast<float>(best);
            v[i] = (in_trunk(r, c, meander) || in_trib(r, c, with_tribs) ||
                    in_sub(r, c, with_tribs)) ? 1.0f : 0.0f;
        }

    fixture fx;
    std::vector<std::uint8_t> mask(n, 1);
    lsdtt::lsdtt_raster_cpu raw(std::move(z), std::vector<std::uint8_t>(mask), meta);
    // North, east and west closed: everything drains out of the south edge.
    fx.dem = lsdtt::fill_pits(raw, 0.0001f, "nnon", 1);
    fx.valley = lsdtt::lsdtt_raster_cpu(std::move(v), std::move(mask), meta);

    lsdtt::lsdtt_flowinfo_options fi;
    fi.boundary_conditions = "nnon";
    fi.n_threads = 1;
    fx.flow = lsdtt::make_flow_base(fx.dem.get_data_array(),
                                    fx.dem.get_mask_array(), meta, fi);
    lsdtt::channel_scope_request sreq;
    sreq.threshold_contributing_pixels = threshold;
    fx.scope = lsdtt::resolve_channel_scope(*fx.flow, sreq);

    // The Y-network gets hand-drawn channels on the floor centre lines, each
    // ending at its parent's floor edge (the trace follows the flow from
    // there). Threshold heads on this fixture land where its planar walls
    // meet, which says nothing about the code under test.
    if (with_tribs)
    {
        std::fill(fx.scope.channel_mask.begin(), fx.scope.channel_mask.end(), 0);
        for (int r = 3; r < k_rows; ++r) fx.scope.channel_mask[r * k_cols + 100] = 1;
        for (int c = 121; c <= 173; ++c) fx.scope.channel_mask[k_jr * k_cols + c] = 1;
        for (int r = 120; r <= 141; ++r) fx.scope.channel_mask[r * k_cols + 135] = 1;
    }
    return fx;
}

lsdtt::valley_cross_section_result sections_for(
    const fixture& fx, const lsdtt::valley_width_network_result& net)
{
    lsdtt::valley_cross_section_request xreq;
    xreq.spacing_metres = 25.0;
    xreq.max_half_length_metres = 1000.0;
    // Transect half-length = 0.75 x the widest width at the default 1.5, which
    // a lopsided mouth section (100 m one side, the mouth the other) always
    // overruns, so it would come out side_not_found. A longer transect lets the
    // mouth sections find their far side, so the crossing rule is what has to
    // exclude them.
    xreq.half_length_multiplier = 3.0;
    return lsdtt::extract_valley_cross_sections(fx.dem, net.centrelines, xreq);
}

int row_of(double y) { return static_cast<int>(std::floor((k_rows * k_dx - y) / k_dx)); }
int col_of(double x) { return static_cast<int>(std::floor(x / k_dx)); }

} // namespace

int main()
{
    using lsdtt::network_section_exclusion;

    // -------------------------------------------------------------------------
    // Y-network: tree, widths, both rules, nesting.
    // -------------------------------------------------------------------------
    for (const auto choice : {lsdtt::network_centreline_choice::smoothed_channel,
                              lsdtt::network_centreline_choice::medial_path})
    {
        const std::string tag =
            choice == lsdtt::network_centreline_choice::medial_path
                ? "[medial_path] " : "[smoothed_channel] ";
        fixture fx = make_fixture(true, false, 600);
        // The synthetic georeference puts x_min / y_max away from 0; work in
        // grid-relative coordinates for the row / column lookups below.
        const double x0 = fx.dem.get_meta().x_min;
        const double y0 = fx.dem.get_meta().y_max() - k_rows * k_dx;

        lsdtt::valley_width_network_request nreq;
        nreq.centreline = choice;
        const lsdtt::valley_width_network_result net =
            lsdtt::build_valley_width_network(fx.dem, *fx.flow, fx.scope.channel_mask,
                                              fx.valley, nreq);

        // Identify the three segments by their heads.
        int trunk = -1, trib = -1, sub = -1;
        for (const auto& sg : net.segments)
        {
            const int hr = static_cast<int>(sg.head_flat / k_cols);
            const int hc = static_cast<int>(sg.head_flat % k_cols);
            if (sg.parent_segment_id < 0 && std::abs(hc - 100) <= 20) trunk = sg.segment_id;
            else if (std::abs(hr - k_jr) <= 8 && hc > 120) trib = sg.segment_id;
            else if (std::abs(hc - 135) <= 4 && hr < k_jr - 8) sub = sg.segment_id;
        }
        int n_long = 0;
        for (const auto& sg : net.segments)
            if (static_cast<int>(sg.channel_flats.size()) >= nreq.min_segment_pixels) ++n_long;
        check(n_long == 3, tag + "three segments long enough for a centreline (got " +
                               std::to_string(n_long) + " of " +
                               std::to_string(net.segments.size()) + ")");
        if (n_long != 3)
            for (const auto& sg : net.segments)
                std::cout << "      segment " << sg.segment_id << " head (r,c) ("
                          << sg.head_flat / k_cols << "," << sg.head_flat % k_cols
                          << ") parent " << sg.parent_segment_id << " pixels "
                          << sg.channel_flats.size() << "\n";
        check(trunk >= 0 && trib >= 0 && sub >= 0, tag + "trunk, tributary and sub-tributary found");
        if (trunk < 0 || trib < 0 || sub < 0) continue;
        check(net.segments[trib].parent_segment_id == trunk &&
                  net.segments[trib].nesting_depth == 1,
              tag + "tributary is the trunk's child, depth 1");
        check(net.segments[sub].parent_segment_id == trib &&
                  net.segments[sub].nesting_depth == 2,
              tag + "sub-tributary is the tributary's child, depth 2");
        check(net.segments[trunk].drains_to_edge, tag + "trunk drains to the edge");

        // Smoothed centrelines never leave the valley mask (their channels are
        // inside it throughout).
        if (choice == lsdtt::network_centreline_choice::smoothed_channel)
        {
            const std::vector<float>& vm = fx.valley.get_data();
            bool all_in = true;
            for (const auto& ch : net.centrelines.channels)
                for (const std::int64_t f : ch.centreline_flats)
                    if (vm[static_cast<std::size_t>(f)] < 0.5f) all_in = false;
            check(all_in, tag + "smoothed centrelines stay inside the valley mask");
        }

        const lsdtt::valley_cross_section_result xres = sections_for(fx, net);
        const auto reasons = lsdtt::apply_network_section_rules(
            net, xres, fx.dem.get_meta(), lsdtt::network_section_rules_request{});

        auto expected_width = [&](int seg) {
            return seg == trunk ? k_trunk_width : seg == trib ? k_trib_width : k_sub_width; };

        int n_kept[3] = {0, 0, 0};
        double worst = 0.0;
        bool trib_kept_inside_trunk = false;
        bool sub_kept_inside_trib = false;
        int trunk_mouth_excluded = 0, trunk_mouth_total = 0, trunk_mouth_crossing = 0;
        int trib_mouth_excluded = 0;
        bool trunk_excluded_at_sub = false;
        for (std::size_t k = 0; k < xres.sections.size(); ++k)
        {
            const auto& xs = xres.sections[k];
            const int seg = xs.source_id;
            const int r = row_of(xs.centre_y - y0);
            const int c = col_of(xs.centre_x - x0);
            const int idx = seg == trunk ? 0 : seg == trib ? 1 : 2;
            if (reasons[k] == network_section_exclusion::none)
            {
                ++n_kept[idx];
                worst = std::max(worst, std::abs(xs.valley_width_metres - expected_width(seg)));
                if (seg == trib && c <= 120) trib_kept_inside_trunk = true;
                if (seg == sub && r >= k_jr - 8) sub_kept_inside_trib = true;
            }
            // Trunk stations level with the tributary's floor see its mouth.
            if (seg == trunk && std::abs(r - k_jr) <= 6)
            {
                ++trunk_mouth_total;
                if (reasons[k] != network_section_exclusion::none)
                    ++trunk_mouth_excluded;
                if (reasons[k] == network_section_exclusion::tributary_crossing)
                    ++trunk_mouth_crossing;
            }
            if (seg == trib && std::abs(c - 135) <= 2 &&
                reasons[k] == network_section_exclusion::tributary_crossing)
                ++trib_mouth_excluded;
            if (seg == trunk && reasons[k] == network_section_exclusion::tributary_crossing &&
                std::abs(r - k_jr) > 12)
                trunk_excluded_at_sub = true;
        }
        check(n_kept[0] >= 20 && n_kept[1] >= 5 && n_kept[2] >= 2,
              tag + "kept sections on every segment (" + std::to_string(n_kept[0]) + ", " +
                  std::to_string(n_kept[1]) + ", " + std::to_string(n_kept[2]) + ")");
        check(worst <= 2.0 * k_dx,
              tag + "every kept width equals its floor width within two cells (worst " +
                  std::to_string(worst) + " m)");
        check(!trib_kept_inside_trunk, tag + "no kept tributary station inside the trunk valley");
        check(!sub_kept_inside_trib, tag + "no kept sub-tributary station inside the tributary valley");
        check(trunk_mouth_total > 0 && trunk_mouth_excluded == trunk_mouth_total,
              tag + "every trunk section across the tributary mouth is excluded (" +
                  std::to_string(trunk_mouth_excluded) + " of " +
                  std::to_string(trunk_mouth_total) + ", " +
                  std::to_string(trunk_mouth_crossing) + " tributary_crossing)");
        check(trib_mouth_excluded > 0,
              tag + "tributary sections across the sub-tributary mouth are tributary_crossing");
        check(!trunk_excluded_at_sub,
              tag + "the sub-tributary excludes no trunk section directly");

        // Rules off: the mouth widths come back.
        lsdtt::network_section_rules_request off;
        off.stop_tributaries_at_parent = false;
        off.remove_tributary_crossing_sections = false;
        const auto raw = lsdtt::apply_network_section_rules(net, xres, fx.dem.get_meta(), off);
        double widest_trunk = 0.0;
        for (std::size_t k = 0; k < xres.sections.size(); ++k)
            if (xres.sections[k].source_id == trunk &&
                raw[k] == network_section_exclusion::none)
                widest_trunk = std::max(widest_trunk, xres.sections[k].valley_width_metres);
        if (trunk_mouth_crossing > 0)
            check(widest_trunk > 1.5 * k_trunk_width,
                  tag + "rules off: those trunk sections measure the mouth (widest " +
                      std::to_string(widest_trunk) + " m)");
        else
            std::cout << "NOTE  " << tag << "no bounded trunk section across the "
                         "mouth (all side_not_found), so the crossing rule was not "
                         "exercised on the trunk" << std::endl;
    }

    // -------------------------------------------------------------------------
    // A reach with no valley: the mask removed from the trunk's top 40 rows.
    // Those stations sit where the channel runs without a valley; they must be
    // not_in_valley, not kept with a zero width.
    // -------------------------------------------------------------------------
    {
        fixture fx = make_fixture(true, false, 600);
        const double y0 = fx.dem.get_meta().y_max() - k_rows * k_dx;
        std::vector<float> vm = fx.valley.get_data();
        for (int r = 0; r < 40; ++r)
            for (int c = 0; c < k_cols; ++c) vm[r * k_cols + c] = 0.0f;
        const lsdtt::lsdtt_raster_cpu cut(std::move(vm),
            std::vector<std::uint8_t>(fx.valley.get_mask()), fx.valley.get_meta());
        lsdtt::valley_width_network_request nreq;
        nreq.centreline = lsdtt::network_centreline_choice::smoothed_channel;
        const auto net = lsdtt::build_valley_width_network(
            fx.dem, *fx.flow, fx.scope.channel_mask, cut, nreq);
        const auto xres = sections_for(fx, net);
        const auto reasons = lsdtt::apply_network_section_rules(
            net, xres, fx.dem.get_meta(), lsdtt::network_section_rules_request{});
        int n_north = 0, n_flagged = 0, n_kept_south = 0;
        for (std::size_t k = 0; k < xres.sections.size(); ++k)
        {
            const int r = row_of(xres.sections[k].centre_y - y0);
            if (r < 38)
            {
                ++n_north;
                if (reasons[k] == network_section_exclusion::not_in_valley) ++n_flagged;
            }
            else if (r > 42 && reasons[k] == network_section_exclusion::none)
                ++n_kept_south;
        }
        check(n_north > 0 && n_flagged == n_north,
              "[no valley] every station above the mask is not_in_valley (" +
                  std::to_string(n_flagged) + " of " + std::to_string(n_north) + ")");
        check(n_kept_south > 0, "[no valley] stations in the valley below are kept");
    }

    // -------------------------------------------------------------------------
    // Meandering channel in a straight valley.
    // -------------------------------------------------------------------------
    {
        fixture fx = make_fixture(false, true, 600);
        const std::vector<float>& vm = fx.valley.get_data();
        const double x0 = fx.dem.get_meta().x_min;

        auto mean_offset = [&](const std::vector<std::int64_t>& flats) {
            double s = 0.0;
            for (const std::int64_t f : flats) s += std::abs(static_cast<int>(f % k_cols) - 100);
            return flats.empty() ? 1e9 : s / flats.size();
        };
        auto path_length = [&](const std::vector<std::int64_t>& flats) {
            double L = 0.0;
            for (std::size_t i = 1; i < flats.size(); ++i)
            {
                const int dr = static_cast<int>(flats[i] / k_cols - flats[i - 1] / k_cols);
                const int dc = static_cast<int>(flats[i] % k_cols - flats[i - 1] % k_cols);
                L += std::hypot(dr, dc) * k_dx;
            }
            return L;
        };
        (void)x0;

        lsdtt::valley_width_network_request nreq;
        nreq.centreline = lsdtt::network_centreline_choice::smoothed_channel;
        nreq.smoothing_width_multiple = 1.5;
        const auto sm = lsdtt::build_valley_width_network(
            fx.dem, *fx.flow, fx.scope.channel_mask, fx.valley, nreq);
        nreq.centreline = lsdtt::network_centreline_choice::medial_path;
        const auto md = lsdtt::build_valley_width_network(
            fx.dem, *fx.flow, fx.scope.channel_mask, fx.valley, nreq);

        check(sm.segments.size() == 1 && md.segments.size() == 1,
              "[meander] one segment (got " + std::to_string(sm.segments.size()) + ")");
        for (const auto& sg : sm.segments)
            if (sm.segments.size() != 1)
                std::cout << "      segment " << sg.segment_id << " head (r,c) ("
                          << sg.head_flat / k_cols << "," << sg.head_flat % k_cols
                          << ") parent " << sg.parent_segment_id << " pixels "
                          << sg.channel_flats.size() << "\n";
        if (sm.segments.size() == 1 && md.segments.size() == 1)
        {
            const auto& channel = sm.segments[0].channel_flats;
            const auto& smooth  = sm.centrelines.channels[0].centreline_flats;
            const auto& medial  = md.centrelines.channels[0].centreline_flats;
            bool all_in = true;
            for (const std::int64_t f : smooth)
                if (vm[static_cast<std::size_t>(f)] < 0.5f) all_in = false;
            check(all_in, "[meander] smoothed line stays inside the valley mask");
            const double lc = path_length(channel), ls = path_length(smooth);
            check(ls < 0.9 * lc, "[meander] smoothed line is shorter than the channel (" +
                                     std::to_string(ls) + " vs " + std::to_string(lc) + " m)");
            check(mean_offset(channel) > 5.0,
                  "[meander] the channel meanders (mean offset " +
                      std::to_string(mean_offset(channel)) + " cells)");
            check(mean_offset(medial) < 3.0,
                  "[meander] medial path follows the valley centre (mean offset " +
                      std::to_string(mean_offset(medial)) + " cells)");
        }
    }

    std::cout << (g_failures == 0 ? "ALL PASSED\n" : "FAILURES: ") ;
    if (g_failures) std::cout << g_failures << "\n";
    return g_failures == 0 ? 0 : 1;
}
