// test_lithology_channel_stats.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Tests for k_sn categorised by lithologic unit (lsdtt_lithology_channel_stats_cpu).
//
//   T-S1 compute_weighted_quantile (lower weighted quantile; exact)
//   T1  contact distance on a two-unit grid (exact)
//   T3  node -> segment on a hand-built fit result (exact)
//   T4  lithology of a set of nodes: majority, purity, ties, unmapped (exact)
//   T5  lithologic reaches on a hand-built network: counts and end types,
//       with and without a contact buffer (exact)
//   T6  reach OLS on z = 3 + 7 chi; min_reach_pixels (exact)
//   T7  summary and histogram bookkeeping (counts add up, rows present,
//       low_sample below 5 qualifying reaches with a NaN interval)
//   T8  bootstrap interval: seed reproducibility, constant values, clusters
//   T9  observation only: the fit result is unchanged
//   T10 known answer, n = 1: a steady-state landscape with K1 = 1e-5 west of a
//       vertical contact and K2 = 4e-5 east of it, U = 2e-4, map area, A0 = 1:
//       k_sn = (U/K)^(1/n) = 20 and 5, ratio 4
//   T11 the same with m = 0.9, n = 2: k_sn = sqrt(20) and sqrt(5), ratio 2
//
// The steady state integrates chi with the area of each cell (left end), the
// chi of lsdtt-chi-analysis is trapezoidal, so the fitted k_sn sits slightly
// off (U/K)^(1/n). Measured 2026-10-17 (120 x 160 at 30 m, seed 11): medians
// 1.14 % and 1.25 % high for both n, ratio within 0.108 %. Tolerances: the
// ratio at twice the measured error (0.25 %); the medians at 2 %, the cap
// (twice the measured 1.25 % would exceed it). Those are the node medians; the
// reach medians (OLS per lithologic reach) and the layout with the contact
// across the drainage (B) have their own measured pins, in main().
//
// Usage:
//   lsdtt-test-lithology-channel-stats <work_dir>
//   lsdtt-test-lithology-channel-stats --write-fixtures <dir> <tiny_dem.asc>
//   lsdtt-test-lithology-channel-stats --check-ksn-by-lithology <summary.csv>
//       <unit_a> <unit_b> <expected_ratio> <relative_tolerance>
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../../src/core/lsdtt_csv_table.hpp"
#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../../src/cpu/lsdtt_chi_network_analysis_cpu.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_lithology_channel_stats_cpu.hpp"
#include "../../src/cpu/lsdtt_lithology_map_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_steady_state_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../../src/utils/lsdtt_statistics.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;
constexpr double k_nan_d = std::numeric_limits<double>::quiet_NaN();

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

std::string fmt(double v)
{
    std::ostringstream s;
    s << std::setprecision(6) << v;
    return s.str();
}

lsdtt_raster_meta utm_grid(int nr, int nc, double dx)
{
    lsdtt_raster_meta m;
    m.n_rows = nr; m.n_cols = nc; m.cell_size = dx;
    m.x_min = 500000.0; m.y_min = 5000000.0; m.no_data_value = -9999.0;
    m.set_crs_from_epsg(32632);
    return m;
}

lsdtt_raster_cpu blank(const lsdtt_raster_meta& m)
{
    const std::size_t n = static_cast<std::size_t>(m.n_rows) * m.n_cols;
    return lsdtt_raster_cpu(std::vector<float>(n, 0.0f), std::vector<std::uint8_t>(n, 1), m);
}

// Units by column: code 1 for col < split, else 2.
lsdtt_lithology_units_on_grid column_units(const lsdtt_raster_meta& m, int split)
{
    lsdtt_lithology_units_on_grid u;
    u.code_of_pixel.resize(static_cast<std::size_t>(m.n_rows) * m.n_cols);
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c)
            u.code_of_pixel[static_cast<std::size_t>(r) * m.n_cols + c] = c < split ? 1 : 2;
    u.units = {{1, "west", {"west"}, 0}, {2, "east", {"east"}, 0}};
    return u;
}

// Units by row: code 1 for row < split, else 2.
lsdtt_lithology_units_on_grid row_units(const lsdtt_raster_meta& m, int split)
{
    lsdtt_lithology_units_on_grid u = column_units(m, 0);
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c)
            u.code_of_pixel[static_cast<std::size_t>(r) * m.n_cols + c] = r < split ? 1 : 2;
    return u;
}

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

const lsdtt_chi_ksn_lithology_summary_row* row_for(const lsdtt_chi_ksn_lithology_summary& s, int code)
{
    for (const auto& r : s.rows)
        if (r.lith_code == code) return &r;
    return nullptr;
}

// ---------------------------------------------------------------------------
// T-S1 weighted quantile (lsdtt_statistics)
// ---------------------------------------------------------------------------
void test_weighted_quantile()
{
    std::cout << "-- T-S1 weighted quantile --\n";
    const std::vector<double> v = {5, 1, 3, 2, 4}, eq(5, 1.0);
    check("equal weights: q = 0, 0.25, 0.5, 1 give 1, 2, 3, 5",
          compute_weighted_quantile(v, eq, 0.0) == 1.0 && compute_weighted_quantile(v, eq, 0.25) == 2.0 &&
              compute_weighted_quantile(v, eq, 0.5) == 3.0 && compute_weighted_quantile(v, eq, 1.0) == 5.0);
    const std::vector<double> v2 = {1, 2, 3, 4}, w2 = {1, 1, 1, 7};
    check("weights {1,1,1,7}: q = 0.25 -> 3 (cumulative 3 >= 2.5), q = 0.5 -> 4, q = 0 -> 1",
          compute_weighted_quantile(v2, w2, 0.25) == 3.0 && compute_weighted_quantile(v2, w2, 0.5) == 4.0 &&
              compute_weighted_quantile(v2, w2, 0.0) == 1.0);
    const std::vector<double> w3 = {0, 2, 0, 2};
    check("zero weights: q = 1 -> 4 (largest with weight), q = 0.5 -> 2",
          compute_weighted_quantile(v2, w3, 1.0) == 4.0 && compute_weighted_quantile(v2, w3, 0.5) == 2.0);
    check("empty input is NaN", std::isnan(compute_weighted_quantile({}, {}, 0.5)));
    auto throws = [](const std::vector<double>& a, const std::vector<double>& w, double q) {
        try { (void)compute_weighted_quantile(a, w, q); } catch (const std::invalid_argument&) { return true; }
        return false;
    };
    check("zero total weight throws", throws({1, 2}, {0, 0}, 0.5));
    check("negative weight, size mismatch and q outside [0, 1] throw",
          throws({1, 2}, {1, -1}, 0.5) && throws({1, 2}, {1}, 0.5) && throws({1, 2}, {1, 1}, 1.5));

    // Equal weights agree with compute_quantile_floor_in_place at q = 0, 0.5, 1.
    std::mt19937_64 g(7);
    std::uniform_int_distribution<int> d(0, 1000);
    bool agree = true;
    for (std::size_t n = 1; n <= 12; ++n)
    {
        std::vector<float> f(n);
        for (auto& x : f) x = static_cast<float>(d(g)) * 0.25f;
        const std::vector<double> dv(f.begin(), f.end()), w(n, 1.0);
        for (double q : {0.0, 0.5, 1.0})
        {
            std::vector<float> c = f;
            agree = agree && compute_weighted_quantile(dv, w, q) == compute_quantile_floor_in_place(c, q);
        }
    }
    check("equal weights = compute_quantile_floor_in_place at q = 0, 0.5, 1 (n = 1..12)", agree);
}

// ---------------------------------------------------------------------------
// T1 contact distance
// ---------------------------------------------------------------------------
void test_contact_distance()
{
    std::cout << "-- T1 contact distance --\n";
    const lsdtt_raster_meta m = utm_grid(20, 20, 10.0);
    lsdtt_raster_cpu g = blank(m);
    std::vector<std::uint8_t> mask = g.get_mask();
    mask[static_cast<std::size_t>(5) * 20 + 2] = 0;   // a nodata pixel in the west
    g = lsdtt_raster_cpu(std::vector<float>(g.get_data()), std::move(mask), m);
    const auto u = column_units(m, 10);
    const std::vector<float> d = compute_lithology_contact_distance(g, u, 1);
    auto at = [&](int r, int c) { return d[static_cast<std::size_t>(r) * 20 + c]; };
    check("contact pixels (cols 9 and 10) are 0", at(3, 9) == 0.0f && at(3, 10) == 0.0f);
    check("col 0 is 9 dx = 90 m, col 19 is 90 m", at(0, 0) == 90.0f && at(19, 19) == 90.0f,
          fmt(at(0, 0)) + ", " + fmt(at(19, 19)));
    check("col 4 is 50 m", at(12, 4) == 50.0f, fmt(at(12, 4)));
    check("grid nodata is NaN", std::isnan(at(5, 2)));

    lsdtt_lithology_units_on_grid one = u;
    for (auto& c : one.code_of_pixel) c = 1;
    const auto none = compute_lithology_contact_distance(blank(m), one, 1);
    check("no contact: all NaN", std::isnan(none[0]) && std::isnan(none[399]));

    lsdtt_lithology_units_on_grid hole = one;
    hole.code_of_pixel[static_cast<std::size_t>(10) * 20 + 10] = -1;
    const auto h = compute_lithology_contact_distance(blank(m), hole, 1);
    check("unmapped counts as a unit (its neighbours are contacts)",
          h[static_cast<std::size_t>(10) * 20 + 11] == 0.0f && h[static_cast<std::size_t>(10) * 20 + 13] == 20.0f);
}

// ---------------------------------------------------------------------------
// T3 node -> segment
// ---------------------------------------------------------------------------
void test_node_segments()
{
    std::cout << "-- T3 node -> segment --\n";
    // Main stem 0 <- 1 <- 2 <- 3 <- 4 <- 5 (0 = outlet), tributary 7 <- 6 -> joins at 3
    // (6 drains to 3). Segments: main [5..3], [3..0]; tributary [7..3].
    lsdtt_chi_network_fit_result r;
    r.flats = {0, 1, 2, 3, 4, 5, 6, 7};
    r.receiver_flats = {k_outlet_receiver, 0, 1, 2, 3, 4, 3, 6};
    r.chi.assign(8, 0.0f);
    r.observed_elevation.assign(8, 0.0f);
    lsdtt_chi_network_fit_segment a, b, c;
    a.segment_id = 0; a.upstream_flat = 3; a.downstream_flat = 0; a.k_sn = 1.0f;
    b.segment_id = 1; b.upstream_flat = 5; b.downstream_flat = 3; b.k_sn = 2.0f;
    c.segment_id = 2; c.upstream_flat = 7; c.downstream_flat = 3; c.k_sn = 3.0f;
    r.segments = {a, b, c};
    const auto s = assign_chi_fit_nodes_to_segments(r);
    const std::vector<std::int32_t> want = {0, 0, 0, 0, 1, 1, 2, 2};
    std::string got;
    for (auto v : s) got += std::to_string(v) + " ";
    check("each node on one segment; confluence on the downstream segment; outlet on the last",
          s == want, got);
}

// ---------------------------------------------------------------------------
// T4 lithology of a set of nodes
// ---------------------------------------------------------------------------
void test_flats_tag()
{
    std::cout << "-- T4 flats tag --\n";
    const std::vector<std::int32_t> code = {1, 1, 2, 2, -1};
    const double dx = 10.0, dd = 10.0 * std::sqrt(2.0);
    const auto a = summarise_lithology_of_flats({0, 1, 2}, {dx, dd, dx}, code);
    check("cardinal + diagonal + cardinal: unit 1, fraction (dx + dx sqrt2) / (2 dx + dx sqrt2)",
          a.code == 1 && a.n_units == 2 && a.majority_fraction == (dx + dd) / (dx + dd + dx),
          fmt(a.majority_fraction));
    const auto t = summarise_lithology_of_flats({0, 1, 2, 3}, {dx, dd, dd, dx}, code);
    check("equal weight in two units: the lower code, fraction 0.5",
          t.code == 1 && t.majority_fraction == 0.5 && t.n_units == 2);
    const auto u = summarise_lithology_of_flats({4, 4}, {dx, dx}, code);
    check("all unmapped: code -1, fraction 1", u.code == -1 && u.majority_fraction == 1.0 && u.n_units == 1);
    const auto z = summarise_lithology_of_flats({0, 2, 3}, {0.0, 0.0, 0.0}, code);
    check("zero weights: nodes counted, unit 2 with 2/3", z.code == 2 && z.majority_fraction == 2.0 / 3.0);
    const auto e = summarise_lithology_of_flats({}, {}, code);
    check("no nodes: code -1, fraction NaN, 0 units", e.code == -1 && std::isnan(e.majority_fraction) && e.n_units == 0);
}

// ---------------------------------------------------------------------------
// T5-T7 reaches on a hand-built network
//
// 40 x 40 grid, dx 10 m, unit 1 west of column 12, unit 2 from column 12.
// Main stem (branch 0) on row 10, columns 0..29, draining west; column 0 is the
// outlet; node k has chi 0.5 k. Tributary (branch 1) on column 20, rows 11..19,
// draining north into the main stem at (10, 20). z = 3 + 7 chi everywhere, so
// every reach has k_sn_ols = 7 exactly (chi and z exact in float).
// ---------------------------------------------------------------------------
struct hand_network
{
    lsdtt_raster_meta meta;
    lsdtt_lithology_units_on_grid units;
    lsdtt_chi_network_fit_collection_result fits;
    std::vector<float> contact;
};

hand_network make_hand_network()
{
    hand_network h;
    h.meta = utm_grid(40, 40, 10.0);
    h.units = column_units(h.meta, 12);
    auto flat = [&](int r, int c) { return static_cast<std::int64_t>(r) * 40 + c; };
    lsdtt_chi_network_fit_result r;
    r.outlet_flat = flat(10, 0);
    for (int k = 0; k < 30; ++k)
    {
        r.flats.push_back(flat(10, k));
        r.receiver_flats.push_back(k == 0 ? k_outlet_receiver : flat(10, k - 1));
        r.chi.push_back(0.5f * k);
        r.branch_id_of_node.push_back(0);
    }
    for (int j = 0; j < 9; ++j)
    {
        r.flats.push_back(flat(11 + j, 20));
        r.receiver_flats.push_back(flat(10 + j, 20));
        r.chi.push_back(10.0f + 0.5f * (j + 1));
        r.branch_id_of_node.push_back(1);
    }
    for (float c : r.chi) r.observed_elevation.push_back(3.0f + 7.0f * c);
    lsdtt_chi_network_fit_segment a, b;
    a.segment_id = 0; a.branch_id = 0; a.upstream_flat = flat(10, 29); a.downstream_flat = flat(10, 0); a.k_sn = 7.0f;
    b.segment_id = 1; b.branch_id = 1; b.upstream_flat = flat(19, 20); b.downstream_flat = flat(10, 20); b.k_sn = 7.0f;
    r.segments = {a, b};
    h.fits.basins.push_back(r);
    h.contact = compute_lithology_contact_distance(blank(h.meta), h.units, 1);
    return h;
}

std::string ends(const lsdtt_chi_lithology_reach_row& r)
{
    return chi_lithology_reach_end_name(r.upstream_end) + "/" + chi_lithology_reach_end_name(r.downstream_end);
}

void test_reaches()
{
    std::cout << "-- T5-T7 reaches --\n";
    const hand_network h = make_hand_network();
    const lsdtt_chi_network_fit_result before = h.fits.basins[0];
    lsdtt_chi_lithology_options o;
    o.min_reach_pixels = 10;
    const auto t = build_chi_lithology_tables(h.fits, h.meta, h.units, h.contact, {}, o);

    check("T5 three reaches", t.reaches.size() == 3, std::to_string(t.reaches.size()));
    if (t.reaches.size() == 3)
    {
        const auto& r0 = t.reaches[0];
        const auto& r1 = t.reaches[1];
        const auto& r2 = t.reaches[2];
        check("T5 reach 0: unit 1, 12 nodes, contact upstream, outlet downstream",
              r0.lith_code == 1 && r0.n_nodes == 12 && ends(r0) == "contact/outlet", ends(r0));
        check("T5 reach 1: unit 2, 18 nodes on branch 0, source upstream, contact downstream",
              r1.lith_code == 2 && r1.n_nodes == 18 && r1.branch_id == 0 && ends(r1) == "source/contact", ends(r1));
        check("T5 reach 2: the tributary, 9 nodes, source upstream, junction downstream",
              r2.lith_code == 2 && r2.n_nodes == 9 && r2.branch_id == 1 && ends(r2) == "source/junction", ends(r2));
        check("T6 z = 3 + 7 chi: k_sn_ols 7, R^2 1, SE 0 (to 1e-12)",
              r0.k_sn_ols == 7.0 && r1.k_sn_ols == 7.0 && std::fabs(r0.r_squared - 1.0) < 1e-12 &&
                  std::fabs(r0.k_sn_ols_se) < 1e-12 && std::fabs(r1.k_sn_ols_se) < 1e-12,
              fmt(r0.k_sn_ols) + " " + fmt(r0.r_squared) + " " + fmt(r0.k_sn_ols_se));
        check("T6 chi span, dz, length and chi-weighted fit k_sn of reach 0",
              r0.chi_span == 5.5 && r0.dz_m == 38.5 && r0.length_m == 110.0 && r0.k_sn_fit_chi_weighted == 7.0,
              fmt(r0.chi_span) + " " + fmt(r0.dz_m) + " " + fmt(r0.length_m));
        check("T6 the 9-node reach does not qualify (min 10), the others do",
              r0.qualifies && r1.qualifies && !r2.qualifies);
    }
    std::size_t n_assigned = 0;
    bool ids_ok = true;
    for (const auto& n : t.nodes.rows)
    {
        if (n.reach_id >= 0) ++n_assigned;
        ids_ok = ids_ok && !n.near_contact;
    }
    check("T5 every node in one reach, none near a contact", n_assigned == 39 && ids_ok);

    // Segments: segment 0 spans the contact.
    check("segments: two rows", t.segments.size() == 2);
    if (t.segments.size() == 2)
    {
        const auto& s0 = t.segments[0];
        const auto& s1 = t.segments[1];
        check("segment 0: majority unit 2 with 180/290 of its length, 2 units, 30 nodes",
              s0.lith.code == 2 && s0.lith.n_units == 2 && s0.lith.majority_fraction == 180.0 / 290.0 &&
                  s0.n_nodes == 30 && s0.length_m == 290.0,
              fmt(s0.lith.majority_fraction) + ", " + fmt(s0.length_m));
        check("segment 1: unit 2 pure, 9 nodes, 90 m",
              s1.lith.code == 2 && s1.lith.majority_fraction == 1.0 && s1.n_nodes == 9 && s1.length_m == 90.0);
    }

    // T7 bookkeeping on the reach summary.
    std::mt19937_64 rng(5489);
    const auto s = summarise_chi_ksn_by_lithology(t, h.units, 0.45, 1.0, o, rng);
    const auto* w = row_for(s, 1);
    const auto* e = row_for(s, 2);
    const auto* un = row_for(s, -1);
    const auto* all = row_for(s, -2);
    check("T7 rows: west, east, unmapped, all", s.rows.size() == 4 && w && e && un && all);
    if (w && e && un && all)
    {
        check("T7 n_nodes add up (12 + 27 + 0 = 39)",
              w->n_nodes + e->n_nodes + un->n_nodes == all->n_nodes && all->n_nodes == 39,
              std::to_string(w->n_nodes) + " + " + std::to_string(e->n_nodes));
        check("T7 channel lengths add up",
              w->channel_length_m + e->channel_length_m + un->channel_length_m == all->channel_length_m);
        check("T7 reaches: 1 / 2 / 0 / 3, qualifying 1 / 1 / 0 / 2",
              w->n_reaches == 1 && e->n_reaches == 2 && un->n_reaches == 0 && all->n_reaches == 3 &&
                  w->n_reaches_qualifying == 1 && e->n_reaches_qualifying == 1 && all->n_reaches_qualifying == 2);
        check("T7 fewer than 5 qualifying reaches: low_sample, NaN interval, no replicates",
              w->low_sample && e->low_sample && all->low_sample && std::isnan(w->reach_k_sn_median_ci_low) &&
                  std::isnan(all->reach_k_sn_median_ci_high) && all->n_bootstrap == 0);
        check("T7 reach median k_sn 7 in both units", w->reach_k_sn_median == 7.0 && e->reach_k_sn_median == 7.0);
        check("T7 unmapped row present with no nodes and NaN statistics",
              un->n_nodes == 0 && std::isnan(un->reach_k_sn_median) && std::isnan(un->node_k_sn_length_weighted_median));
        check("T7 m/n, A0, buffer and min reach copied",
              s.m_over_n == 0.45 && s.A_0 == 1.0 && s.contact_buffer_m == 0.0 && s.min_reach_pixels == 10);
        check("T7 node length-weighted median k_sn 7", all->node_k_sn_length_weighted_median == 7.0);
    }
    std::size_t basin_rows = 0;
    for (const auto& b : s.basin_rows) basin_rows += b.basin_index == 0 ? 1 : 0;
    check("T7 unit x basin: one row per unit with nodes (2)", s.basin_rows.size() == 2 && basin_rows == 2);

    // T5 with a buffer of 2 dx: columns 10..13 are nearer than 20 m to the contact.
    o.contact_buffer_m = 20.0;
    const auto tb = build_chi_lithology_tables(h.fits, h.meta, h.units, h.contact, {}, o);
    std::size_t n_near = 0;
    for (const auto& n : tb.nodes.rows) n_near += n.near_contact ? 1 : 0;
    check("T5 buffer 20 m: 4 nodes excluded", n_near == 4, std::to_string(n_near));
    if (tb.reaches.size() == 3)
    {
        check("T5 buffer: reach 0 is 10 nodes, buffer/outlet",
              tb.reaches[0].n_nodes == 10 && ends(tb.reaches[0]) == "buffer/outlet", ends(tb.reaches[0]));
        check("T5 buffer: reach 1 is 16 nodes, source/buffer",
              tb.reaches[1].n_nodes == 16 && ends(tb.reaches[1]) == "source/buffer", ends(tb.reaches[1]));
        check("T5 buffer: the tributary is unchanged", ends(tb.reaches[2]) == "source/junction");
    }
    else
        check("T5 buffer: three reaches", false, std::to_string(tb.reaches.size()));
    std::mt19937_64 rng2(5489);
    const auto sb = summarise_chi_ksn_by_lithology(tb, h.units, 0.45, 1.0, o, rng2);
    const auto* wb = row_for(sb, 1);
    const auto* eb = row_for(sb, 2);
    const auto* ab = row_for(sb, -2);
    check("T5 buffer: excluded 2 + 2 = 4 in the summary, first-cut n_nodes unchanged",
          wb && eb && ab && wb->n_nodes_excluded_near_contact == 2 && eb->n_nodes_excluded_near_contact == 2 &&
              ab->n_nodes_excluded_near_contact == 4 && ab->n_nodes == 39);

    // T9 observation only.
    const auto& after = h.fits.basins[0];
    bool same = after.flats == before.flats && after.receiver_flats == before.receiver_flats &&
                after.chi == before.chi && after.observed_elevation == before.observed_elevation &&
                after.branch_id_of_node == before.branch_id_of_node && after.segments.size() == before.segments.size();
    for (std::size_t i = 0; same && i < after.segments.size(); ++i)
        same = after.segments[i].k_sn == before.segments[i].k_sn &&
               after.segments[i].upstream_flat == before.segments[i].upstream_flat;
    check("T9 the fit result is unchanged by the tables and the summary", same);
}

// ---------------------------------------------------------------------------
// T8 bootstrap interval of the median
// ---------------------------------------------------------------------------
void test_bootstrap()
{
    std::cout << "-- T8 bootstrap --\n";
    std::mt19937_64 g(1);
    std::normal_distribution<double> nd(10.0, 2.0);
    std::vector<double> v(40);
    for (auto& x : v) x = nd(g);
    std::mt19937_64 a(5489), b(5489), c(5490);
    const auto ca = bootstrap_lithology_median_ci(v, {}, 1000, a);
    const auto cb = bootstrap_lithology_median_ci(v, {}, 1000, b);
    const auto cc = bootstrap_lithology_median_ci(v, {}, 1000, c);
    check("same seed: bitwise-identical interval", same_bits(ca.low, cb.low) && same_bits(ca.high, cb.high));
    check("different seed: a different interval", !(same_bits(ca.low, cc.low) && same_bits(ca.high, cc.high)),
          fmt(ca.low) + "," + fmt(ca.high) + " vs " + fmt(cc.low) + "," + fmt(cc.high));
    std::vector<double> m = v;
    const double med = compute_median_in_place(m);
    check("the interval brackets the point median", ca.low <= med && med <= ca.high && ca.low < ca.high,
          fmt(ca.low) + " <= " + fmt(med) + " <= " + fmt(ca.high));
    check("1000 replicates reported", ca.n_replicates == 1000);
    std::mt19937_64 d(5489);
    const auto cst = bootstrap_lithology_median_ci(std::vector<double>(20, 3.5), {}, 500, d);
    check("constant values: interval [v, v]", cst.low == 3.5 && cst.high == 3.5);
    std::mt19937_64 e(5489);
    const auto cl = bootstrap_lithology_median_ci({1, 1, 1, 10, 10, 10}, {0, 0, 0, 1, 1, 1}, 1000, e);
    check("basin clusters: replicate medians only 1, 5.5 or 10, interval inside [1, 10]",
          cl.low >= 1.0 && cl.high <= 10.0 && (cl.low == 1.0 || cl.low == 5.5) && (cl.high == 10.0 || cl.high == 5.5),
          fmt(cl.low) + "," + fmt(cl.high));
    std::mt19937_64 f(5489);
    const auto off = bootstrap_lithology_median_ci(v, {}, 0, f);
    check("0 replicates: NaN, none reported", std::isnan(off.low) && std::isnan(off.high) && off.n_replicates == 0);
}

// ---------------------------------------------------------------------------
// The steady-state fixture
// ---------------------------------------------------------------------------
struct steady_fixture
{
    lsdtt_raster_cpu surface;      // the noisy surface the flow is routed on
    lsdtt_raster_cpu elevation;    // the steady state
    lsdtt_lithology_units_on_grid units;
};

const int k_nr = 120, k_nc = 160, k_split = 80;
const double k_dx = 30.0;
const float k_U = 2e-4f, k_K1 = 1e-5f, k_K2 = 4e-5f;

lsdtt_flowinfo d8_flow(const lsdtt_raster_cpu& dem)
{
    lsdtt_flowinfo_options fo;
    fo.boundary_conditions = "oooo";
    fo.n_threads = 1;
    fo.build_topological_order = true;
    fo.build_s_vector = true;
    return lsdtt_flowinfo(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fo);
}

lsdtt_raster_cpu tilted_surface()
{
    lsdtt_raster_cpu z = generate_red_noise_surface(k_nr, k_nc, k_dx, 11, 0.0, 3.0, 6, 1);
    std::vector<float> d = z.get_data();
    for (int r = 0; r < k_nr; ++r)
        for (int c = 0; c < k_nc; ++c)
            d[static_cast<std::size_t>(r) * k_nc + c] += static_cast<float>(0.02 * k_dx * r);
    lsdtt_raster_meta m = z.get_meta();
    m.x_min = 500000.0; m.y_min = 5000000.0;
    m.set_crs_from_epsg(32632);
    return fill_pits(lsdtt_raster_cpu(std::move(d), std::vector<std::uint8_t>(z.get_mask()), m),
                     1e-4f, 1);
}

// Layout A: the contact is a column (parallel to the drainage, which runs
// north); layout B: a row (across it), so every trunk crosses it.
enum class contact_layout { strike_parallel, cross_strike };
const int k_split_row = 60;

steady_fixture make_steady(float m_exp, float n_exp, contact_layout layout = contact_layout::strike_parallel)
{
    const bool by_col = layout == contact_layout::strike_parallel;
    steady_fixture f;
    f.surface = tilted_surface();
    const lsdtt_flowinfo flow = d8_flow(f.surface);
    lsdtt_fluvial_steady_state_options o;
    o.U = k_U;
    o.K = k_K1;
    o.m = m_exp;
    o.n = n_exp;
    o.divide_area_by_dx = false;   // map area, as chi with A0 = 1
    o.erodibility_field.resize(static_cast<std::size_t>(k_nr) * k_nc);
    for (int r = 0; r < k_nr; ++r)
        for (int c = 0; c < k_nc; ++c)
            o.erodibility_field[static_cast<std::size_t>(r) * k_nc + c] =
                (by_col ? c < k_split : r < k_split_row) ? k_K1 : k_K2;
    f.elevation = compute_d8_fluvial_steady_state(flow, o).elevation;
    f.units = by_col ? column_units(f.surface.get_meta(), k_split)
                     : row_units(f.surface.get_meta(), k_split_row);
    return f;
}

lsdtt_chi_network_fit_collection_result fit_on(const lsdtt_raster_cpu& surface,
                                               const lsdtt_raster_cpu& elevation, double mn)
{
    const lsdtt_flowinfo flow = d8_flow(surface);
    const lsdtt_raster_result npx = flow.get_n_contributing_pixels_raster();
    std::vector<std::int64_t> outlets;
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(flow.full_grid_size()); ++i)
        if (flow.is_valid_cell(i) && flow.is_outlet(i) && npx.data[static_cast<std::size_t>(i)] >= 50.0f)
            outlets.push_back(i);
    lsdtt_chi_options co;
    co.m_over_n = mn;
    co.A_0 = 1.0;
    lsdtt_chi_collection_result chi = get_chi_from_outlets_zero_base_threshold(flow, outlets, 50, co);
    (void)reject_chi_basins_with_few_channel_pixels(chi, 30);
    lsdtt_chi_network_fit_collection_result fits;
    for (const auto& b : chi.basins) fits.basins.push_back(fit_chi_network(b, elevation, {}));
    return fits;
}

struct reach_pins
{
    double tol_rel = 0.0, tol_ratio = 0.0;     // reach medians and their ratio
    std::size_t n_reaches_west = 0, n_reaches_east = 0;   // pinned exactly
};

std::size_t count_contact_reaches(const lsdtt_chi_lithology_tables& t)
{
    std::size_t n = 0;
    for (const auto& r : t.reaches)
        if (r.upstream_end == lsdtt_chi_lithology_reach_end::contact ||
            r.downstream_end == lsdtt_chi_lithology_reach_end::contact)
            ++n;
    return n;
}

// Reach statistics on a steady state: the per-unit reach medians against
// (U/K)^(1/n), their ratio, and the pinned reach counts. Returns the tables.
lsdtt_chi_lithology_tables check_reach_answer(const std::string& tag, const steady_fixture& f, double mn,
                                              double n_exp, const reach_pins& pins)
{
    const auto fits = fit_on(f.surface, f.elevation, mn);
    const std::vector<float> cd = compute_lithology_contact_distance(f.surface, f.units, 1);
    lsdtt_chi_lithology_options o;   // defaults: buffer 0, 10 nodes, 1000 replicates by reach
    const auto t = build_chi_lithology_tables(fits, f.surface.get_meta(), f.units, cd, {}, o);
    std::mt19937_64 rng(5489);
    const auto s = summarise_chi_ksn_by_lithology(t, f.units, mn, 1.0, o, rng);
    const double want1 = std::pow(k_U / k_K1, 1.0 / n_exp), want2 = std::pow(k_U / k_K2, 1.0 / n_exp);
    const auto* w = row_for(s, 1);
    const auto* e = row_for(s, 2);
    if (!w || !e || w->n_reaches_qualifying == 0 || e->n_reaches_qualifying == 0)
    {
        check(tag + "both units have qualifying reaches", false);
        return t;
    }
    const double err1 = std::fabs(w->reach_k_sn_median / want1 - 1.0);
    const double err2 = std::fabs(e->reach_k_sn_median / want2 - 1.0);
    const double ratio = w->reach_k_sn_median / e->reach_k_sn_median;
    const double err_ratio = std::fabs(ratio / (want1 / want2) - 1.0);
    std::cout << "   measured (reaches): west " << std::setprecision(8) << w->reach_k_sn_median
              << " (want " << want1 << ", rel err " << err1 << ", " << w->n_reaches << " reaches, "
              << w->n_reaches_qualifying << " qualifying, CI [" << w->reach_k_sn_median_ci_low << ", "
              << w->reach_k_sn_median_ci_high << "]), east " << e->reach_k_sn_median << " (want "
              << want2 << ", rel err " << err2 << ", " << e->n_reaches << " reaches, "
              << e->n_reaches_qualifying << " qualifying, CI [" << e->reach_k_sn_median_ci_low << ", "
              << e->reach_k_sn_median_ci_high << "]), ratio " << ratio << " (rel err " << err_ratio
              << "), contact-ended reaches " << count_contact_reaches(t) << std::setprecision(6) << "\n";
    check(tag + "west reach median k_sn = (U/K1)^(1/n)", err1 <= pins.tol_rel,
          "rel err " + fmt(err1) + ", tolerance " + fmt(pins.tol_rel));
    check(tag + "east reach median k_sn = (U/K2)^(1/n)", err2 <= pins.tol_rel,
          "rel err " + fmt(err2) + ", tolerance " + fmt(pins.tol_rel));
    check(tag + "reach median ratio = (K2/K1)^(1/n)", err_ratio <= pins.tol_ratio,
          fmt(ratio) + ", rel err " + fmt(err_ratio) + ", tolerance " + fmt(pins.tol_ratio));
    check(tag + "reach counts pinned (west " + std::to_string(pins.n_reaches_west) + ", east " +
              std::to_string(pins.n_reaches_east) + ")",
          w->n_reaches == pins.n_reaches_west && e->n_reaches == pins.n_reaches_east,
          std::to_string(w->n_reaches) + ", " + std::to_string(e->n_reaches));
    check(tag + "intervals bracket the reach medians",
          w->reach_k_sn_median_ci_low <= w->reach_k_sn_median && w->reach_k_sn_median <= w->reach_k_sn_median_ci_high &&
              e->reach_k_sn_median_ci_low <= e->reach_k_sn_median && e->reach_k_sn_median <= e->reach_k_sn_median_ci_high);
    std::mt19937_64 rng_again(5489);
    const auto s2 = summarise_chi_ksn_by_lithology(t, f.units, mn, 1.0, o, rng_again);
    const auto* w2 = row_for(s2, 1);
    check(tag + "same seed, same interval (bitwise)",
          w2 && same_bits(w2->reach_k_sn_median_ci_low, w->reach_k_sn_median_ci_low) &&
              same_bits(w2->reach_k_sn_median_ci_high, w->reach_k_sn_median_ci_high));
    return t;
}

void test_known_answer(float m_exp, float n_exp, double tol_rel, double tol_ratio, const reach_pins& pins_a,
                       const reach_pins& pins_b)
{
    const std::string tag = "n = " + fmt(n_exp) + ": ";
    std::cout << "-- T10/T11 steady state, m = " << m_exp << ", n = " << n_exp << " --\n";
    const steady_fixture f = make_steady(m_exp, n_exp);
    const double mn = static_cast<double>(m_exp) / n_exp;
    const auto fits = fit_on(f.surface, f.elevation, mn);
    const std::vector<float> cd = compute_lithology_contact_distance(f.surface, f.units, 1);
    const auto before = fits.basins.empty() ? lsdtt_chi_network_fit_result{} : fits.basins[0];
    const auto nodes = build_chi_lithology_node_table(fits, f.surface.get_meta(), f.units, cd);
    check(tag + "fits untouched by the tagging (observation only)",
          !fits.basins.empty() && fits.basins[0].flats == before.flats &&
              fits.basins[0].segments.size() == before.segments.size());
    const auto s = summarise_chi_ksn_by_lithology(nodes, f.units, mn, 1.0);

    const double want1 = std::pow(k_U / k_K1, 1.0 / n_exp), want2 = std::pow(k_U / k_K2, 1.0 / n_exp);
    const auto* w = row_for(s, 1);
    const auto* e = row_for(s, 2);
    if (!w || !e || w->n_nodes == 0 || e->n_nodes == 0)
    {
        check(tag + "both units have channel nodes", false);
        return;
    }
    const double err1 = std::fabs(w->k_sn_median / want1 - 1.0);
    const double err2 = std::fabs(e->k_sn_median / want2 - 1.0);
    const double ratio = w->k_sn_median / e->k_sn_median, want_ratio = want1 / want2;
    const double err_ratio = std::fabs(ratio / want_ratio - 1.0);
    std::cout << "   measured: west median " << w->k_sn_median << " (want " << want1 << ", "
              << w->n_nodes << " nodes, " << w->n_segments << " segments), east "
              << e->k_sn_median << " (want " << want2 << ", " << e->n_nodes << " nodes, "
              << e->n_segments << " segments), ratio " << ratio << "\n";
    check(tag + "west (K1) median k_sn = (U/K1)^(1/n)", err1 <= tol_rel,
          "rel err " + fmt(err1) + ", tolerance " + fmt(tol_rel));
    check(tag + "east (K2) median k_sn = (U/K2)^(1/n)", err2 <= tol_rel,
          "rel err " + fmt(err2) + ", tolerance " + fmt(tol_rel));
    check(tag + "ratio = (K2/K1)^(1/n)", err_ratio <= tol_ratio,
          fmt(ratio) + " vs " + fmt(want_ratio) + ", rel err " + fmt(err_ratio));

    // T7 bookkeeping on the same tables.
    const auto* all = row_for(s, -2);
    const auto* un = row_for(s, -1);
    check(tag + "all = west + east + unmapped (nodes, length)",
          all && un && all->n_nodes == w->n_nodes + e->n_nodes + un->n_nodes &&
              std::fabs(all->channel_length_m -
                        (w->channel_length_m + e->channel_length_m + un->channel_length_m)) < 1e-6 * all->channel_length_m);
    check(tag + "no unmapped nodes", un && un->n_nodes == 0 && s.n_nodes_unmapped == 0);
    check(tag + "m/n and A0 copied", s.m_over_n == mn && s.A_0 == 1.0);

    lsdtt_chi_ksn_histogram_options ho;
    ho.n_bins = 25;
    const auto h = make_chi_ksn_histogram_by_lithology(nodes, f.units, mn, 1.0, ho);
    std::map<int, std::size_t> sum, rows;
    for (const auto& r : h.rows) { sum[r.lith_code] += r.n_nodes; ++rows[r.lith_code]; }
    check(tag + "histogram: 27 rows (25 bins + below + above) per unit, unmapped and all",
          rows.size() == 4 && rows[1] == 27 && rows[2] == 27 && rows[-1] == 27 && rows[-2] == 27);
    check(tag + "histogram counts = summary node counts",
          sum[1] == w->n_nodes && sum[2] == e->n_nodes && sum[-2] == all->n_nodes);
    ho.log_bins = true;
    const auto hl = make_chi_ksn_histogram_by_lithology(nodes, f.units, mn, 1.0, ho);
    std::size_t lsum = 0;
    for (const auto& r : hl.rows) if (r.lith_code == -2) lsum += r.n_nodes;
    check(tag + "log histogram keeps every node", lsum == all->n_nodes);

    std::cout << "   layout A (contact parallel to the drainage):\n";
    const auto ta = check_reach_answer(tag + "layout A: ", f, mn, n_exp, pins_a);
    std::cout << "   layout B (contact across the drainage):\n";
    const steady_fixture fb = make_steady(m_exp, n_exp, contact_layout::cross_strike);
    const auto tb = check_reach_answer(tag + "layout B: ", fb, mn, n_exp, pins_b);
    check(tag + "layout B has more contact-ended reaches than layout A",
          count_contact_reaches(tb) > count_contact_reaches(ta),
          std::to_string(count_contact_reaches(tb)) + " vs " + std::to_string(count_contact_reaches(ta)));
    // The units of the nodes each segment owns (from the node table) against
    // the segment tag: crossing segments (two units) have n_lith_units = 2.
    std::map<std::pair<int, int>, std::set<int>> units_of_segment;
    for (const auto& r : tb.nodes.rows)
        if (r.segment_id >= 0) units_of_segment[{r.basin_index, r.segment_id}].insert(r.lith_code);
    std::size_t n_cross = 0, n_agree = 0;
    for (const auto& sg : tb.segments)
    {
        const std::size_t nu = units_of_segment[{sg.basin_index, sg.segment_id}].size();
        n_cross += nu == 2 ? 1 : 0;
        n_agree += static_cast<std::size_t>(sg.lith.n_units) == nu ? 1 : 0;
    }
    std::cout << "   layout B: " << n_cross << " of " << tb.segments.size() << " segments cross the contact\n";
    check(tag + "layout B: segments crossing the contact exist and every segment's n_lith_units "
                "matches its nodes (crossing = 2)",
          n_cross > 0 && n_agree == tb.segments.size(),
          std::to_string(n_cross) + " crossing, " + std::to_string(n_agree) + " of " +
              std::to_string(tb.segments.size()) + " agree");
}

// ---------------------------------------------------------------------------
// --check-ksn-by-lithology: the ratio of two units' reach_k_sn_median in a
// written _chi_ksn_by_lithology_ CSV (read by header) against an expected
// ratio. Both units must have enough reaches (low_sample = 0). Returns the
// measured relative error, or NaN with `why` set.
// ---------------------------------------------------------------------------
double ksn_by_lithology_ratio_error(const std::string& csv, const std::string& unit_a,
                                    const std::string& unit_b, double expect_ratio, std::string& why)
{
    const lsdtt_csv_table t = load_csv_table(csv);
    double a = k_nan_d, b = k_nan_d;
    for (std::size_t i = 0; i < t.rows.size(); ++i)
    {
        const auto& r = t.rows[i];
        const std::string name = csv_row_string(r, "lith_name");
        if (name != unit_a && name != unit_b) continue;
        if (csv_row_double_required(r, "low_sample", csv, static_cast<int>(i) + 1) != 0.0)
        {
            why = name + " has low_sample = 1";
            return k_nan_d;
        }
        (name == unit_a ? a : b) = csv_row_double_required(r, "reach_k_sn_median", csv, static_cast<int>(i) + 1);
    }
    if (!std::isfinite(a) || !std::isfinite(b) || b == 0.0)
    {
        why = "no finite reach_k_sn_median for " + unit_a + " and " + unit_b;
        return k_nan_d;
    }
    std::cout << "reach_k_sn_median " << unit_a << " " << std::setprecision(10) << a << ", " << unit_b
              << " " << b << ", ratio " << a / b << " (expected " << expect_ratio << "), rel err "
              << std::fabs(a / b / expect_ratio - 1.0) << std::setprecision(6) << "\n";
    return std::fabs(a / b / expect_ratio - 1.0);
}

int check_ksn_by_lithology(int argc, char* argv[])
{
    if (argc != 7)
    {
        std::cerr << "usage: --check-ksn-by-lithology <summary.csv> <unit_a> <unit_b> "
                     "<expected_ratio> <relative_tolerance>\n";
        return 2;
    }
    std::string why;
    const double tol = std::stod(argv[6]);
    try
    {
        const double err = ksn_by_lithology_ratio_error(argv[2], argv[3], argv[4], std::stod(argv[5]), why);
        if (std::isnan(err))
        {
            std::cerr << "FAIL: " << why << "\n";
            return 1;
        }
        std::cout << (err <= tol ? "PASS" : "FAIL") << " ratio within " << tol << "\n";
        return err <= tol ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}

// The check mode on a summary written here from the steady state.
void test_check_mode(const std::string& dir)
{
    std::cout << "-- --check-ksn-by-lithology on a written summary --\n";
    const steady_fixture f = make_steady(0.45f, 1.0f);
    const auto fits = fit_on(f.surface, f.elevation, 0.45);
    const auto cd = compute_lithology_contact_distance(f.surface, f.units, 1);
    lsdtt_chi_lithology_options o;
    o.bootstrap_replicates = 50;
    const auto t = build_chi_lithology_tables(fits, f.surface.get_meta(), f.units, cd, {}, o);
    std::mt19937_64 rng(5489);
    const auto s = summarise_chi_ksn_by_lithology(t, f.units, 0.45, 1.0, o, rng);
    const std::string csv = dir + "/check_mode_summary.csv";
    write_chi_ksn_lithology_summary_csv(csv, s);
    std::string why;
    const double err = ksn_by_lithology_ratio_error(csv, "west", "east", 4.0, why);
    const double want = std::fabs(row_for(s, 1)->reach_k_sn_median / row_for(s, 2)->reach_k_sn_median / 4.0 - 1.0);
    check("check mode reads the ratio back from the CSV (to 1e-8)", std::fabs(err - want) < 1e-8,
          fmt(err) + " vs " + fmt(want));
    const double wrong = ksn_by_lithology_ratio_error(csv, "west", "nowhere", 4.0, why);
    check("check mode: an absent unit is reported", std::isnan(wrong) && !why.empty(), why);
}

// ---------------------------------------------------------------------------
// --write-fixtures: a steady-state DEM and its unit map for the CLI, and a
// two-unit code raster on tiny_dem.asc's grid for the observation-only check.
// ---------------------------------------------------------------------------
void write_polygon_units(const std::string& path, const lsdtt_raster_meta& m_in, int split,
                         double shift_x = 0.0)
{
    lsdtt_raster_meta m = m_in;
    m.x_min += shift_x;
    gdal_init();
    std::filesystem::remove(path);
    GDALDataset* ds = GetGDALDriverManager()->GetDriverByName("GPKG")->Create(
        path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    OGRSpatialReference srs;
    srs.importFromEPSG(32632);
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRLayer* layer = ds->CreateLayer("units", &srs, wkbPolygon, nullptr);
    OGRFieldDefn fd("lith_primary", OFTString);
    (void)layer->CreateField(&fd);
    const double xs = m.x_min + split * m.cell_size, x1 = m.x_min + m.n_cols * m.cell_size;
    const double y0 = m.y_min, y1 = m.y_min + m.n_rows * m.cell_size;
    for (int k = 0; k < 2; ++k)
    {
        OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetField("lith_primary", k == 0 ? "sandstone" : "shale");
        const double a = k == 0 ? m.x_min : xs, b = k == 0 ? xs : x1;
        OGRLinearRing ring;
        ring.addPoint(a, y0); ring.addPoint(b, y0); ring.addPoint(b, y1); ring.addPoint(a, y1);
        ring.addPoint(a, y0);
        OGRPolygon p;
        p.addRing(&ring);
        f->SetGeometry(&p);
        if (layer->CreateFeature(f) != OGRERR_NONE) std::cerr << "could not write a unit polygon\n";
        OGRFeature::DestroyFeature(f);
    }
    GDALClose(ds);
}

int write_fixtures(const std::string& dir, const std::string& tiny_dem)
{
    std::filesystem::create_directories(dir);
    gdal_init();
    const steady_fixture f = make_steady(0.45f, 1.0f);
    lsdtt_data_loader w;
    w.write_full(f.elevation.get_data(), f.elevation.get_meta(), dir + "/lithks_steady_dem.tif",
                 raster_output_format::gtiff);
    write_polygon_units(dir + "/lithks_units.gpkg", f.elevation.get_meta(), k_split);
    // The same units 200 km east, and a code raster there: maps that label no
    // valid DEM pixel (same CRS, no overlap).
    write_polygon_units(dir + "/lithks_units_far.gpkg", f.elevation.get_meta(), k_split, 200000.0);
    {
        lsdtt_raster_meta far = f.elevation.get_meta();
        far.x_min += 200000.0;
        std::vector<float> c(static_cast<std::size_t>(far.n_rows) * far.n_cols, 1.0f);
        w.write_full(c, far, dir + "/lithks_codes_far.tif", raster_output_format::gtiff);
    }
    {
        // A crosswalk with unit names, and one missing the value "shale".
        std::ofstream cw(dir + "/lithks_crosswalk.csv");
        cw << "value,code,name\nsandstone,10,quartz sandstone\nshale,30,grey shale\n";
        std::ofstream gap(dir + "/lithks_crosswalk_gap.csv");
        gap << "value,code\nsandstone,10\n";
    }

    lsdtt_data_loader l;
    l.open_raster(tiny_dem);
    const lsdtt_raster_meta tm = l.get_meta();
    std::vector<float> codes(static_cast<std::size_t>(tm.n_rows) * tm.n_cols);
    for (int r = 0; r < tm.n_rows; ++r)
        for (int c = 0; c < tm.n_cols; ++c)
            codes[static_cast<std::size_t>(r) * tm.n_cols + c] = c < tm.n_cols / 2 ? 1.0f : 2.0f;
    w.write_full(codes, tm, dir + "/lithks_tiny_codes.tif", raster_output_format::gtiff);
    // tiny_dem beside its code raster, so one read_path serves both.
    const std::filesystem::path src(tiny_dem);
    for (const char* ext : {".asc", ".prj"})
    {
        std::filesystem::path from = src;
        from.replace_extension(ext);
        if (std::filesystem::exists(from))
            std::filesystem::copy_file(from, std::filesystem::path(dir) / from.filename(),
                                       std::filesystem::copy_options::overwrite_existing);
    }
    std::cout << "fixtures written to " << dir << "\n";
    return 0;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc > 3 && std::string(argv[1]) == "--write-fixtures") return write_fixtures(argv[2], argv[3]);
    if (argc > 1 && std::string(argv[1]) == "--check-ksn-by-lithology") return check_ksn_by_lithology(argc, argv);
    const std::string dir = argc > 1 ? argv[1] : ".";
    std::filesystem::create_directories(dir);
    gdal_init();
    std::cout << "== k_sn by lithologic unit ==\n";
    test_weighted_quantile();
    test_contact_distance();
    test_node_segments();
    test_flats_tag();
    test_reaches();
    test_bootstrap();
    // Reach pins, measured 2026-10-17 on the same fixture (identical relative
    // errors for n = 1 and n = 2): layout A reach medians +1.759 % (west) and
    // +1.310 % (east), ratio 0.443 %; layout B +1.192 % and +1.581 %, ratio
    // 0.383 %. Medians pinned at the 2 % cap (twice the measured error would
    // exceed it); ratios at about twice the measured error, 0.9 % (A) and
    // 0.8 % (B), under the 1 % cap. Reach counts pinned exactly.
    const reach_pins pa{0.02, 0.009, 46, 42}, pb{0.02, 0.008, 57, 38};
    test_known_answer(0.45f, 1.0f, 0.02, 0.0025, pa, pb);
    test_known_answer(0.9f, 2.0f, 0.02, 0.0025, pa, pb);
    test_check_mode(dir);
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
