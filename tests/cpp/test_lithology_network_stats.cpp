// test_lithology_network_stats.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Channel-link bearings and junction angles by lithologic unit
// (lsdtt_lithology_channel_stats_cpu, network part, with the category-keyed
// statistics of lsdtt_channel_bearing_statistics_cpu).
//
// Fixture: the N-S trellis of test_channel_bearing_statistics.cpp (300 x 160 at
// 10 m; four N-S trunk valleys draining south, E-W tributaries; D8, threshold
// 8 pixels), with lithologic units given in process on the grid:
//   layout A: unit 1 west of column 80, unit 2 from column 80. Column 80 is
//             the divide between the second and third valleys, so no channel
//             crosses the contact.
//   layout B: unit 1 north of row 150, unit 2 from row 150: every trunk
//             crosses the contact.
//   layout U: layout A with columns >= 140 unmapped (-1).
//
//   T12 links: on A every link takes the unit of its side with purity 1; on B
//       a crossing trunk link's purity equals the length-weighted fraction
//       computed here independently; a fixed_n_pixels sample at the
//       upstream end moves the tag to the sampled nodes; unmapped links
//       form their own group (never excluded).
//   T13 junctions and summaries: interior <= all, with exact counts; the
//       axial 10-degree rose by unit equals counts made here; the junction
//       median bootstrap (reproducible, fills only its own rows, matches the
//       row's n); appended record columns and the quoting of unit names with
//       commas; pixels_draining_through against a brute-force walk.
//
// Usage: lsdtt-test-lithology-network-stats <work_dir>
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_csv_table.hpp"
#include "../../src/core/lsdtt_flow_base.hpp"
#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/core/lsdtt_junction_network_base.hpp"
#include "../../src/cpu/lsdtt_channel_bearing_statistics_cpu.hpp"
#include "../../src/cpu/lsdtt_channel_network_analysis_cpu.hpp"
#include "../../src/cpu/lsdtt_channel_network_grouping_cpu.hpp"
#include "../../src/cpu/lsdtt_lithology_channel_stats_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0, g_checks = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    ++g_checks;
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

std::string fmt(double v)
{
    std::ostringstream s;
    s << std::setprecision(10) << v;
    return s.str();
}

// ---------------------------------------------------------------------------
// The N-S trellis (as in test_channel_bearing_statistics.cpp)
// ---------------------------------------------------------------------------
constexpr int k_rows = 300, k_cols = 160, k_threshold = 8;
constexpr double k_dx = 10.0, k_slope = 0.02, k_relief = 50.0;
struct valley { int col, west_px, east_px; };
const std::vector<valley> k_valleys = {{6, 6, 34}, {46, 6, 34}, {114, 34, 6}, {154, 34, 6}};

double flank_height(int c)
{
    const double pi = std::acos(-1.0);
    for (const valley& v : k_valleys)
    {
        if (c >= v.col - v.west_px && c <= v.col)
            return k_relief * 0.5 * (1.0 - std::cos(pi * (v.col - c) / v.west_px));
        if (c >= v.col && c <= v.col + v.east_px)
            return k_relief * 0.5 * (1.0 - std::cos(pi * (c - v.col) / v.east_px));
    }
    throw std::runtime_error("flank_height: column outside every valley");
}

struct routed
{
    lsdtt_raster_cpu dem;
    std::unique_ptr<lsdtt_flow_base> flow;
    std::unique_ptr<lsdtt_junction_network_base> network;
};

routed route_trellis()
{
    lsdtt_raster_meta m;
    m.n_rows = k_rows;
    m.n_cols = k_cols;
    m.cell_size = k_dx;
    m.no_data_value = -9999.0;
    m.set_default_synthetic_georeference();
    std::vector<float> z(static_cast<std::size_t>(k_rows) * k_cols);
    for (int r = 0; r < k_rows; ++r)
        for (int c = 0; c < k_cols; ++c)
            z[static_cast<std::size_t>(r) * k_cols + c] =
                static_cast<float>(k_slope * (k_rows - 1 - r) * k_dx + flank_height(c));
    routed rn;
    rn.dem = fill_pits(lsdtt_raster_cpu(std::move(z), std::vector<std::uint8_t>(z.size(), 1), m),
                       0.0001f, "nnon", 1);
    lsdtt_flowinfo_options fi;
    fi.boundary_conditions = "nnon";
    fi.interior_nodata_is_outlet = true;
    fi.require_acyclic_graph = true;
    fi.build_topological_order = true;
    fi.build_s_vector = true;
    fi.n_threads = 1;
    rn.flow = make_flow_base(rn.dem.get_data_array(), rn.dem.get_mask_array(), rn.dem.get_meta(), fi);
    lsdtt_junction_network_options jn;
    jn.build_topological_order = true;
    jn.build_s_vector = true;
    rn.network = make_junction_network_base(*rn.flow, rn.flow->get_sources_index_threshold(k_threshold), jn);
    return rn;
}

// Units: code_at(row, col); names "west"/"east" (or "north"/"south").
template <class F>
lsdtt_lithology_units_on_grid make_units(F code_at, const std::string& name1, const std::string& name2)
{
    lsdtt_lithology_units_on_grid u;
    u.code_of_pixel.resize(static_cast<std::size_t>(k_rows) * k_cols);
    std::size_t n1 = 0, n2 = 0;
    for (int r = 0; r < k_rows; ++r)
        for (int c = 0; c < k_cols; ++c)
        {
            const int code = code_at(r, c);
            u.code_of_pixel[static_cast<std::size_t>(r) * k_cols + c] = code;
            n1 += code == 1;
            n2 += code == 2;
            u.n_unmapped_valid += code == -1;
        }
    u.units = {{1, name1, {name1}, n1}, {2, name2, {name2}, n2}};
    u.n_valid = u.code_of_pixel.size();
    return u;
}

int col_of(std::int64_t f) { return static_cast<int>(f % k_cols); }
int row_of(std::int64_t f) { return static_cast<int>(f / k_cols); }

// ---------------------------------------------------------------------------
// T12 links
// ---------------------------------------------------------------------------
void test_links(const routed& t)
{
    std::cout << "-- T12 link tags --\n";
    lsdtt_link_sampling_options full;   // full_link
    const auto links = calculate_channel_link_bearings(*t.flow, *t.network, {}, full);
    const auto units_a = make_units([](int, int c) { return c < 80 ? 1 : 2; }, "west", "east");
    const auto tags_a = tag_channel_link_bearings_with_lithology(*t.flow, links, units_a);

    std::size_t n_west = 0, n_east = 0, n_bad = 0;
    for (std::size_t i = 0; i < links.size(); ++i)
    {
        const int side = col_of(links[i].sampled_flats.front()) < 80 ? 1 : 2;
        bool one_side = true;
        for (const std::int64_t f : links[i].sampled_flats) one_side = one_side && ((col_of(f) < 80 ? 1 : 2) == side);
        if (!one_side || tags_a[i].code != side || tags_a[i].majority_fraction != 1.0 || tags_a[i].n_units != 1)
            ++n_bad;
        (side == 1 ? n_west : n_east) += 1;
    }
    std::cout << "   layout A: " << links.size() << " links, " << n_west << " west, " << n_east << " east\n";
    check("T12 layout A: every link takes the unit of its side with purity 1 (2396 links, 1198 / 1198)",
          n_bad == 0 && links.size() == 2396 && n_west == 1198 && n_east == 1198,
          std::to_string(n_bad) + " mismatches");
    const auto cats_a = lithology_categories_of_links(tags_a, units_a, 0.8);
    std::size_t c1 = 0, c2 = 0;
    for (auto c : cats_a.category) { c1 += c == 1; c2 += c == 2; }
    check("T12 layout A: categories 1 and 2 hold every link, none excluded",
          c1 == n_west && c2 == n_east && cats_a.n_excluded() == 0 && cats_a.key == "lithology" &&
              cats_a.name_of(1) == "west" && cats_a.name_of(2) == "east");

    // Layout B: the contact crosses every trunk.
    const auto units_b = make_units([](int r, int) { return r < 150 ? 1 : 2; }, "north", "south");
    const auto tags_b = tag_channel_link_bearings_with_lithology(*t.flow, links, units_b);
    std::size_t n_cross = 0, n_cross_ok = 0;
    double example = -1.0;
    for (std::size_t i = 0; i < links.size(); ++i)
    {
        // Independent length-weighted fraction: each node weighted by its D8
        // step to its receiver.
        std::map<int, double> w;
        double total = 0.0;
        for (const std::int64_t f : links[i].sampled_flats)
        {
            double step = 0.0;
            if (!t.flow->is_outlet(f))
            {
                const std::int64_t rc = t.flow->receiver_of(f);
                step = (row_of(rc) != row_of(f) && col_of(rc) != col_of(f)) ? k_dx * std::sqrt(2.0) : k_dx;
            }
            w[row_of(f) < 150 ? 1 : 2] += step;
            total += step;
        }
        if (w.size() < 2) continue;
        ++n_cross;
        const int maj = w[1] >= w[2] ? 1 : 2;
        const double frac = w[maj] / total;
        if (tags_b[i].code == maj && std::fabs(tags_b[i].majority_fraction - frac) < 1e-12 && tags_b[i].n_units == 2)
            ++n_cross_ok;
        if (example < 0.0) example = frac;
    }
    std::cout << "   layout B: " << n_cross << " links cross the contact; first purity " << fmt(example) << "\n";
    check("T12 layout B: crossing links' unit and purity equal the independent length-weighted values",
          n_cross == 12 && n_cross_ok == n_cross, std::to_string(n_cross_ok) + " of " + std::to_string(n_cross));
    const auto cats_b = lithology_categories_of_links(tags_b, units_b, 0.8);
    std::size_t expect_mixed = 0, n_mixed = 0;
    for (const auto& tg : tags_b) expect_mixed += tg.majority_fraction < 0.8;
    for (const auto c : cats_b.category) n_mixed += c == k_lithology_category_mixed;
    check("T12 layout B: the 4 links below purity 0.8 are 'mixed (<80% one unit)', none excluded",
          n_mixed == expect_mixed && n_mixed == 4 && cats_b.n_excluded() == 0 &&
              cats_b.name_of(k_lithology_category_mixed) == "mixed (<80% one unit)" &&
              cats_b.is_descriptive_only(k_lithology_category_mixed) &&
              cats_b.n_descriptive_only() == 4,
          std::to_string(n_mixed));
    check("mixed name follows the threshold: 0.75 -> 'mixed (<75% one unit)'",
          lithology_mixed_category_name(0.75) == "mixed (<75% one unit)");

    // Mixed is summarised but never compared, and no record is lost.
    std::vector<double> bearing_b, length_b;
    for (const auto& l : links) { bearing_b.push_back(l.bearing_deg); length_b.push_back(l.total_length_metres); }
    const std::vector<std::int64_t> no_blocks(links.size(), -1);
    channel_bearing_statistics_options ob;
    ob.do_directional = false;
    ob.min_group_size = 1;
    const auto rows_b = summarise_channel_link_bearings_by_category(bearing_b, length_b, no_blocks, cats_b, nullptr, ob);
    std::size_t n_sum = 0;
    bool mixed_row = false;
    for (const auto& rw : rows_b)
    {
        n_sum += rw.summary.n;
        if (rw.group_value == "mixed (<80% one unit)")
            mixed_row = rw.summary.n == 4 && rw.note.find("descriptive only") != std::string::npos;
    }
    check("mixed has its own summary row (n 4, note 'descriptive only'); the group n sum to every link",
          mixed_row && n_sum == links.size(), std::to_string(n_sum) + " of " + std::to_string(links.size()));
    std::size_t n_north = 0, n_south = 0;
    for (const auto c : cats_b.category) { n_north += c == 1; n_south += c == 2; }
    bool none_mixed = true, rest_ok = false;
    for (const auto mode : {channel_bearing_comparison_mode::pairwise, channel_bearing_comparison_mode::one_vs_rest})
    {
        ob.comparison_mode = mode;
        const auto cmp = compare_channel_link_bearings_by_category(bearing_b, no_blocks, cats_b, nullptr, ob);
        for (const auto& rw : cmp.rows)
        {
            none_mixed = none_mixed && rw.group_a.rfind("mixed", 0) != 0 && rw.group_b.rfind("mixed", 0) != 0;
            if (mode == channel_bearing_comparison_mode::one_vs_rest && rw.group_a == "north")
                rest_ok = rw.n_a == n_north && rw.n_b == n_south;   // rest = south only, no mixed
        }
    }
    check("mixed is in no comparison row and not in 'rest' (north vs rest: rest = the south links only)",
          none_mixed && rest_ok, std::to_string(n_north) + " / " + std::to_string(n_south));

    // The tag follows the SAMPLED nodes. The crossing links are E-W
    // tributaries that cross the contact as they turn into the trunk, at
    // their downstream end; 4 nodes sampled at the UPSTREAM end miss it.
    lsdtt_link_sampling_options fixed;
    fixed.mode = lsdtt_link_sampling_mode::fixed_n_pixels;
    fixed.fixed_n_pixels = 4;
    const auto links_f = calculate_channel_link_bearings(*t.flow, *t.network, {}, fixed,
                                                         lsdtt_link_sample_anchor::upstream_end);
    const auto tags_f = tag_channel_link_bearings_with_lithology(*t.flow, links_f, units_b);
    std::map<std::int32_t, std::size_t> full_index;
    for (std::size_t i = 0; i < links.size(); ++i) full_index[links[i].upstream_junction] = i;
    std::size_t n_changed = 0;
    for (std::size_t i = 0; i < links_f.size(); ++i)
    {
        const auto it = full_index.find(links_f[i].upstream_junction);
        if (it == full_index.end()) continue;
        if (tags_b[it->second].n_units == 2 && tags_f[i].n_units == 1 && tags_f[i].majority_fraction == 1.0)
            ++n_changed;
    }
    std::size_t n_long_cross = 0;
    for (std::size_t i = 0; i < links.size(); ++i)
        n_long_cross += tags_b[i].n_units == 2 && links[i].sampled_flats.size() > 4;
    check("T12 fixed_n_pixels 4 at the upstream end: every crossing link longer than 4 nodes "
          "becomes pure (8 links)", n_changed == n_long_cross && n_changed == 8,
          std::to_string(n_changed) + " of " + std::to_string(n_long_cross));

    // Unmapped is a group, never excluded.
    const auto units_u = make_units([](int, int c) { return c >= 140 ? -1 : (c < 80 ? 1 : 2); }, "west", "east");
    const auto tags_u = tag_channel_link_bearings_with_lithology(*t.flow, links, units_u);
    const auto cats_u = lithology_categories_of_links(tags_u, units_u, 0.8);
    std::size_t n_unmapped = 0, n_unmapped_tags = 0;
    for (std::size_t i = 0; i < links.size(); ++i)
    {
        n_unmapped += cats_u.category[i] == k_lithology_category_unmapped;
        n_unmapped_tags += tags_u[i].code == -1 && tags_u[i].majority_fraction >= 0.8;
    }
    check("T12 unmapped links form the group 'unmapped' (category " +
              std::to_string(k_lithology_category_unmapped) + "), not excluded",
          n_unmapped > 0 && n_unmapped == n_unmapped_tags && cats_u.name_of(k_lithology_category_unmapped) == "unmapped",
          std::to_string(n_unmapped));

    // Appended columns: three, one row per record, and the writer appends them.
    const auto cols = lithology_columns_of_links(tags_a, units_a);
    check("T12 link columns: lith_code, lith_name, lith_majority_fraction",
          cols.names == std::vector<std::string>({"lith_code", "lith_name", "lith_majority_fraction"}) &&
              cols.values.size() == links.size());
}

// ---------------------------------------------------------------------------
// T13 junctions, rose, bootstrap, writers
// ---------------------------------------------------------------------------
void test_junctions_and_summaries(const routed& t, const std::string& dir)
{
    std::cout << "-- T13 junctions and summaries --\n";
    lsdtt_link_sampling_options full;
    const auto angles = calculate_junction_angles(*t.flow, *t.network, {}, full);
    const auto units_b = make_units([](int r, int) { return r < 150 ? 1 : 2; }, "north", "south");
    const auto jt = tag_junction_angles_with_lithology(*t.flow, *t.network, angles, units_b, 0.8);
    const auto all = lithology_categories_of_junctions(jt, units_b, false, "lithology");
    const auto interior = lithology_categories_of_junctions(jt, units_b, true, "lithology_interior");
    std::size_t n_all = 0, n_int = 0, n_agree = 0;
    for (std::size_t i = 0; i < angles.size(); ++i)
    {
        n_all += all.category[i] >= 0;
        n_int += interior.category[i] >= 0;
        // Independent agreement: all three arms' majority = the junction's unit.
        const int jc = units_b.code_of_pixel[static_cast<std::size_t>(t.network->flat_of_junction(angles[i].junction))];
        bool ok = jt[i].valid && jt[i].junction_code == jc;
        for (const auto* a : {&jt[i].donor_1, &jt[i].donor_2, &jt[i].receiver})
            ok = ok && a->code == jc && a->majority_fraction >= 0.8;
        n_agree += ok;
    }
    std::cout << "   layout B: " << angles.size() << " junction records, lithology " << n_all
              << ", lithology_interior " << n_int << "\n";
    std::size_t n_disagree = 0;
    for (const auto c : interior.category) n_disagree += c == k_lithology_category_arms_disagree;
    std::size_t n_int_units = n_int - n_disagree;
    check("T13 every valid junction is in 'lithology' (1192) and in 'lithology_interior' (1180 in "
          "their unit, 12 'arms disagree', descriptive only)",
          n_all == angles.size() && n_all == 1192 && n_int == 1192 && n_int_units == n_agree &&
              n_int_units == 1180 && n_disagree == 12 &&
              interior.is_descriptive_only(k_lithology_category_arms_disagree) &&
              interior.name_of(k_lithology_category_arms_disagree) == "arms disagree",
          std::to_string(n_int_units) + " / " + std::to_string(n_disagree));

    // On layout A no arm crosses the contact, so every junction is interior.
    const auto units_a = make_units([](int, int c) { return c < 80 ? 1 : 2; }, "west", "east");
    const auto jt_a = tag_junction_angles_with_lithology(*t.flow, *t.network, angles, units_a, 0.8);
    std::size_t n_int_a = 0;
    for (const auto& x : jt_a) n_int_a += x.all_arms_agree;
    check("T13 layout A: every junction is interior", n_int_a == angles.size(),
          std::to_string(n_int_a) + " of " + std::to_string(angles.size()));

    // Axial 10-degree rose by unit (layout A) against counts made here.
    const auto links = calculate_channel_link_bearings(*t.flow, *t.network, {}, full);
    const auto tags = tag_channel_link_bearings_with_lithology(*t.flow, links, units_a);
    const auto groups = lithology_categories_of_links(tags, units_a, 0.8);
    std::vector<double> bearing, length;
    for (const auto& l : links) { bearing.push_back(l.bearing_deg); length.push_back(l.total_length_metres); }
    channel_bearing_statistics_options o;
    o.do_directional = false;
    const auto rose = channel_link_bearing_rose_by_category(bearing, length, groups, nullptr, o);
    std::map<std::pair<std::string, int>, std::size_t> want;
    for (std::size_t i = 0; i < links.size(); ++i)
    {
        double ax = std::fmod(bearing[i], 180.0);
        if (ax < 0) ax += 180.0;
        int bin = static_cast<int>(std::floor((ax + 5.0) / 10.0)) % 18;
        ++want[{groups.name_of(groups.category[i]), bin}];
    }
    std::size_t n_bins_ok = 0, n_bins = 0;
    for (const auto& b : rose)
    {
        if (b.stratum_value != "all") continue;
        ++n_bins;
        const int bin = static_cast<int>(std::lround(b.bin.bin_centre_deg / 10.0)) % 18;
        n_bins_ok += b.bin.count == want[{b.group_value, bin}];
    }
    std::cout << "   rose: west bin 0 " << want[{"west", 0}] << ", bin 9 " << want[{"west", 9}]
              << "; east bin 0 " << want[{"east", 0}] << ", bin 9 " << want[{"east", 9}] << "\n";
    check("T13 axial rose by unit: every bin equals the count made here (2 units x 18 bins)",
          n_bins == 36 && n_bins_ok == 36 && want[{"west", 0}] == 596 && want[{"west", 9}] == 602 &&
              want[{"east", 0}] == 596 && want[{"east", 9}] == 602,
          std::to_string(n_bins_ok) + " of " + std::to_string(n_bins));

    // Junction median bootstrap.
    std::vector<double> alpha;
    for (const auto& a : angles) alpha.push_back(a.donor_opening_angle_deg);
    std::vector<std::int64_t> blocks(angles.size(), -1);
    auto rows = summarise_junction_angles_by_category("alpha", alpha, blocks, all, nullptr, o);
    auto rows_order = summarise_junction_angles_by_category(
        "alpha", alpha, blocks, strahler_order_categories([&] {
            std::vector<std::int32_t> v;
            for (const auto& a : angles) v.push_back(a.junction_order);
            return v;
        }(), 0), nullptr, o);
    rows.insert(rows.end(), rows_order.begin(), rows_order.end());
    auto rows2 = rows;
    const auto rows_fresh = rows;
    std::mt19937_64 g1(5489), g2(5489);
    add_junction_angle_median_bootstrap(rows, "alpha", alpha, all, nullptr, 500, 5, g1);
    add_junction_angle_median_bootstrap(rows2, "alpha", alpha, all, nullptr, 500, 5, g2);
    bool same = true, filled_lith = true, untouched = true, brackets = true;
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        same = same && std::memcmp(&rows[i].median_ci_low_deg, &rows2[i].median_ci_low_deg, sizeof(double)) == 0 &&
               std::memcmp(&rows[i].median_ci_high_deg, &rows2[i].median_ci_high_deg, sizeof(double)) == 0;
        if (rows[i].group_key == "lithology")
        {
            filled_lith = filled_lith && rows[i].n_bootstrap == 500;
            brackets = brackets && rows[i].median_ci_low_deg <= rows[i].summary.median_deg &&
                       rows[i].summary.median_deg <= rows[i].median_ci_high_deg;
        }
        else
            untouched = untouched && rows[i].n_bootstrap == 0 && std::isnan(rows[i].median_ci_low_deg);
    }
    check("T13 junction median bootstrap: same seed identical, lithology rows filled and bracket the "
          "median, other keys untouched",
          same && filled_lith && untouched && brackets);
    std::mt19937_64 g3(5489);
    auto wrong = rows_fresh;
    bool threw = false;
    try
    {
        add_junction_angle_median_bootstrap(
            wrong, "alpha", std::vector<double>(alpha.size(), std::nan("")), all, nullptr, 10, 5, g3);
    }
    catch (const std::logic_error&) { threw = true; }
    check("T13 bootstrap refuses rows that do not match the values (all NaN here)", threw);
    std::mt19937_64 g4(5489);
    auto few = rows_fresh;
    add_junction_angle_median_bootstrap(few, "alpha", alpha, all, nullptr, 100, 1000000, g4);
    bool none = true;
    for (const auto& r : few) none = none && r.n_bootstrap == 0 && std::isnan(r.median_ci_high_deg);
    check("T13 rows below min_n keep NaN and n_bootstrap 0", none);

    // Writers: the junction summary gains the interval columns; names with a
    // comma are quoted and read back intact.
    auto named = rows;
    for (auto& r : named) if (r.group_value == "north") r.group_value = "Tonoloway, Wills Creek";
    const std::string f = dir + "/junction_summary.csv";
    write_junction_angle_group_summary_csv(named, f);
    const lsdtt_csv_table tab = load_csv_table(f);
    bool found = false;
    for (const auto& r : tab.rows) found = found || csv_row_string(r, "group_value") == "Tonoloway, Wills Creek";
    check("T13 junction summary: interval columns at the end, a comma in a unit name quoted",
          tab.headers.size() == 22 && tab.headers[14] == "median_ci_low_deg" &&
              tab.headers[16] == "n_bootstrap" && tab.headers[21] == "note" && found);

    // Junction columns appended by the record writer.
    const auto jcols = lithology_columns_of_junctions(jt, units_b);
    const std::string jf = dir + "/junction_angles.csv";
    write_junction_angle_data(angles, *t.flow, *t.network, jf, point_output_format::csv, "junction_angles", false, &jcols);
    const lsdtt_csv_table jtab = load_csv_table(jf);
    const std::vector<std::string> tail(jtab.headers.end() - 9, jtab.headers.end());
    check("T13 _junction_angles gains the nine lithology columns at the end",
          tail == jcols.names && jtab.rows.size() == angles.size());
    const std::string jf0 = dir + "/junction_angles_plain.csv";
    write_junction_angle_data(angles, *t.flow, *t.network, jf0, point_output_format::csv, "junction_angles", false);
    check("T13 without extra columns the writer is unchanged (9 fewer columns)",
          load_csv_table(jf0).headers.size() + 9 == jtab.headers.size());
}

// ---------------------------------------------------------------------------
// pixels_draining_through
// ---------------------------------------------------------------------------
void test_draining_mask(const routed& t)
{
    std::cout << "-- analysed-area mask --\n";
    std::vector<std::int64_t> seeds;
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(t.network->n_junctions()); j += 7)
        seeds.push_back(t.network->flat_of_junction(j));
    const auto mask = pixels_draining_through(*t.flow, seeds);
    const std::set<std::int64_t> seed_set(seeds.begin(), seeds.end());
    std::size_t n_bad = 0, n_in = 0;
    for (std::int64_t p = 0; p < static_cast<std::int64_t>(mask.size()); ++p)
    {
        bool hit = false;
        std::int64_t cur = p;
        for (int guard = 0; guard < k_rows * k_cols; ++guard)
        {
            if (seed_set.count(cur)) { hit = true; break; }
            if (t.flow->is_outlet(cur)) break;
            cur = t.flow->receiver_of(cur);
        }
        n_bad += hit != (mask[static_cast<std::size_t>(p)] == 1);
        n_in += hit;
    }
    check("pixels_draining_through = a brute-force walk to the seeds", n_bad == 0 && n_in > 0,
          std::to_string(n_in) + " pixels, " + std::to_string(n_bad) + " differ");
}

// ---------------------------------------------------------------------------
// Block order rule and block permutations
// ---------------------------------------------------------------------------
double axial_v(const std::vector<double>& a, const std::vector<double>& b, int)
{
    return compare_channel_bearing_samples(a, b, channel_bearing_angular_type::axial).kuiper_V;
}

void test_block_rule_and_permutation(const routed& t)
{
    std::cout << "-- block order rule and permutations --\n";
    // Synthetic order counts: order 1 120, 2 45, 3 31, 4 12, 5 3.
    const std::vector<std::size_t> counts = {120, 45, 31, 12, 3};
    const auto a = choose_channel_block_order(counts, 0, 30);
    check("auto: the largest k with >= 30 basins is 3", a.automatic && a.block_order == 3 && a.note.empty());
    const auto b = choose_channel_block_order(counts, 0, 40);
    check("auto with min 40: k = 2", b.block_order == 2);
    const auto o = choose_channel_block_order(counts, 4, 30);
    check("override 4 is used as given, with a note (12 < 30)",
          !o.automatic && o.block_order == 4 && o.note.find("12 basin(s)") != std::string::npos, o.note);
    const auto n = choose_channel_block_order({20, 6, 1}, 0, 30);
    check("nothing reaches min: no blocks, note names the largest count",
          n.block_order == 0 &&
              n.note == "network too small for block permutation: largest order-k basin count 20 < min_blocks 30",
          n.note);

    // With no blocks every permutation p is nan with that note.
    std::vector<double> v = {0, 1, 2, 90, 91, 92};
    channel_record_categories g;
    g.key = "region";
    g.category = {1, 1, 1, 2, 2, 2};
    channel_bearing_statistics_options op;
    op.do_directional = false;
    op.min_group_size = 1;
    op.n_permutations = 99;
    op.permutation_unavailable_note = n.note;
    const auto cn = compare_channel_link_bearings_by_category(v, std::vector<std::int64_t>(6, -1), g, nullptr, op);
    check("no blocks: permutation_p nan, note carries the reason",
          cn.rows.size() == 1 && std::isnan(cn.rows[0].permutation_p) &&
              cn.rows[0].note.find("network too small") != std::string::npos);

    // Known extreme separation: 40 blocks, 20 all at ~0 deg (group a), 20 at ~90.
    std::vector<double> val;
    std::vector<std::int64_t> blk;
    std::vector<std::size_t> ia, ib;
    for (int bk = 0; bk < 40; ++bk)
        for (int q = 0; q < 3; ++q)
        {
            (bk < 20 ? ia : ib).push_back(val.size());
            val.push_back((bk < 20 ? 0.0 : 90.0) + q + 0.1 * bk);
            blk.push_back(bk);
        }
    std::mt19937_64 r1(7), r2(7);
    const auto p1 = channel_block_permutation_test(val, blk, ia, ib, false, 199, r1, &axial_v, 1);
    const auto p2 = channel_block_permutation_test(val, blk, ia, ib, false, 199, r2, &axial_v, 1);
    check("extreme separation between blocks: p = 1 / (1 + 199)", p1.p == 1.0 / 200.0 && p1.n_blocks == 40,
          fmt(p1.p));
    check("same seed, identical p and statistic", p1.p == p2.p && p1.observed == p2.observed);
    std::mt19937_64 r3(7);
    const auto pw = channel_block_permutation_test(val, blk, ia, ib, true, 199, r3, &axial_v, 1);
    check("within blocks when no block holds both groups: nan, noted",
          std::isnan(pw.p) && pw.note == "no block holds both groups", pw.note);

    // The trellis by Strahler order, blocks = order-1 basins: order-2 links
    // have no block, so order 1 vs 2 has no blocked records in group_b.
    lsdtt_link_sampling_options full;
    const auto links = calculate_channel_link_bearings(*t.flow, *t.network, {}, full);
    std::vector<double> bear, len;
    std::vector<std::int32_t> js, ords;
    for (const auto& l : links)
    {
        bear.push_back(l.bearing_deg);
        len.push_back(l.total_length_metres);
        js.push_back(l.upstream_junction);
        ords.push_back(l.stream_order);
    }
    const auto groups = strahler_order_categories(ords, 0);
    channel_bearing_statistics_options ot;
    ot.do_directional = false;
    ot.n_permutations = 99;
    ot.permute_within_blocks = true;
    ot.block_order = 1;
    const auto blocks1 = block_ids_of_records(js, order_k_basin_of_junctions(*t.network, 1));
    const auto c1 = compare_channel_link_bearings_by_category(bear, blocks1, groups, nullptr, ot);
    check("order 1 vs 2 with order-1 blocks: no blocked order-2 record, p nan with a note",
          c1.rows.size() == 1 && std::isnan(c1.rows[0].permutation_p) &&
              c1.rows[0].note.find("no blocked records in group_b") != std::string::npos,
          c1.rows.empty() ? "" : c1.rows[0].note);
    const auto counts_t = count_order_k_basins_with_records(*t.network, js);
    std::cout << "   trellis order-k basins with links: " << counts_t[0] << ", " << counts_t[1] << "\n";
    ot.block_order = 2;
    const auto blocks2 = block_ids_of_records(js, order_k_basin_of_junctions(*t.network, 2));
    const auto c2a = compare_channel_link_bearings_by_category(bear, blocks2, groups, nullptr, ot);
    const auto c2b = compare_channel_link_bearings_by_category(bear, blocks2, groups, nullptr, ot);
    std::vector<double> b1, b2;
    std::size_t n_unblocked = 0;
    for (std::size_t i = 0; i < links.size(); ++i)
    {
        if (blocks2[i] < 0) { ++n_unblocked; continue; }
        (ords[i] == 1 ? b1 : b2).push_back(bear[i]);
    }
    const double v_blocked = compare_channel_bearing_samples(b1, b2, channel_bearing_angular_type::axial).kuiper_V;
    check("order-2 blocks: p in (0, 1], identical for the same seed; the statistic is Kuiper V on "
          "the blocked records (4 order-1 links without a block left out, as the note says)",
          c2a.rows.size() == 1 && std::isfinite(c2a.rows[0].permutation_p) &&
              c2a.rows[0].permutation_p > 0.0 && c2a.rows[0].permutation_p <= 1.0 &&
              c2a.rows[0].permutation_p == c2b.rows[0].permutation_p &&
              c2a.rows[0].permutation_statistic == v_blocked && n_unblocked == 4 &&
              c2a.rows[0].note.find("4 record(s) without a block left out") != std::string::npos,
          c2a.rows.empty() ? "" : fmt(c2a.rows[0].permutation_p) + ", " + c2a.rows[0].note);
}

} // namespace

int main(int argc, char* argv[])
{
    const std::string dir = argc > 1 ? argv[1] : ".";
    std::filesystem::create_directories(dir);
    gdal_init();
    std::cout << "== channel-network records by lithologic unit ==\n";
    const routed t = route_trellis();
    std::cout << "trellis: " << t.network->n_junctions() << " junctions\n";
    test_links(t);
    test_junctions_and_summaries(t, dir);
    test_draining_mask(t);
    test_block_rule_and_permutation(t);
    std::cout << (g_checks - g_failures) << " / " << g_checks << " checks passed\n"
              << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
