// test_channel_bearing_statistics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Channel-link bearing and junction-angle statistics
// (lsdtt_channel_bearing_statistics_cpu, lsdtt_channel_network_grouping_cpu)
// and the CLI fixture for lsdtt-channel-network-analysis.
//
// Modes, one per CTest test:
//
//   --circular-two-sample  The analytic cases of the circular two-sample
//                          functions in lsdtt_statistical_tests (unchanged
//                          numerics). The first four are the cases the
//                          retired lsdtt-bearing-comparison CLI tests pinned,
//                          with the same 8 bearings, now native.
//   --group-summaries      Circular mean, R-bar, axial fold, length weights,
//                          Rayleigh, D8 fraction, rose bins.
//   --comparisons          Holm, the pairwise / one_vs_rest / reference row
//                          counts, the minimum group size, the reference-CSV
//                          reader.
//   --junction-angle-stats Linear statistics of junction angles.
//   --trellis-network      Known-answer synthetic trellis networks (flow
//                          routing, junction network, bearings, angles) and
//                          the sampled_flats on the records.
//   --grouping-keys        Region (with purity), order pooling, donor order
//                          pairs, order-k basins, innermost-wins outlet basins.
//   --isotropic-null       A diamond-square surface: the all-link axial R-bar
//                          measured against the 0.15 bound.
//   --write-fixtures DIR [REFERENCE_CSV]
//                          The CLI fixture (DEMs, region rasters, the
//                          reference CSV copied in).
//
// THE TRELLIS. 10 m cells, 300 rows x 160 columns. z = S (distance north of
// the south edge) + h(x), S = 0.02, h a sum of half-cosine flanks, A = 50 m.
// Four N-S valleys at columns 6, 46, 114, 154. Each valley has one SHORT
// flank (6 px, too short to reach the 8-pixel channel threshold) and one LONG
// flank (34 px, which sheds one tributary per row). The short flank is on the
// west for the two western valleys and on the east for the two eastern ones,
// so tributaries flow west (270) in half the valleys and east (90) in the
// other half: their directional mean cancels and their axial mean does not.
// One channelled flank per valley also makes every trunk confluence a
// two-donor junction (with both flanks channelled, each trunk pixel would take
// a tributary from each side and have three donors, and no junction angle).
// Boundary conditions nnon: only the south edge is open, so every valley is a
// whole basin.
//
// Flank flow is cardinal (E-W) except within ~2 pixels of the valley floor and
// of the crest, where the 0.02 down-valley gradient makes the diagonal step
// steeper, so tributary bearings are a little off 90 / 270. The tolerances
// below allow for that; they are derived from the geometry, not tuned.
//
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_gdal_init.hpp"
#include "lsdtt_flow_base.hpp"
#include "lsdtt_junction_network_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "lsdtt_channel_network_analysis_cpu.hpp"
#include "lsdtt_channel_bearing_statistics_cpu.hpp"
#include "lsdtt_channel_network_grouping_cpu.hpp"
#include "lsdtt_statistical_tests.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS: " : "FAIL: ") << name;
    if (!detail.empty()) std::cout << "  [" << detail << "]";
    std::cout << "\n";
}

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(10) << v;
    return os.str();
}

void check_near(const std::string& name, double value, double expected, double tol)
{
    const bool ok = std::isfinite(value) && std::abs(value - expected) <= tol;
    check(name, ok, "got " + num(value) + ", expected " + num(expected) + " +/- " + num(tol));
}

void check_nan(const std::string& name, double value)
{
    check(name, std::isnan(value), "got " + num(value) + ", expected nan");
}

// Smallest angular difference on a circle of `period`.
double circular_diff(double a, double b, double period)
{
    double d = std::fmod(std::abs(a - b), period);
    return std::min(d, period - d);
}

void check_circular_near(const std::string& name, double value, double expected,
                         double period, double tol)
{
    const bool ok = std::isfinite(value) && circular_diff(value, expected, period) <= tol;
    check(name, ok, "got " + num(value) + ", expected " + num(expected) + " +/- " + num(tol) +
                    " (period " + num(period) + ")");
}

using lsdtt::channel_bearing_angular_type;

// =============================================================================
// --circular-two-sample
// =============================================================================
const std::vector<double> k_sample = {10, 50, 100, 170, 200, 260, 300, 345};

void run_circular_two_sample_checks()
{
    constexpr double tol = 1e-9;
    std::vector<double> plus20, plus180;
    for (const double v : k_sample)
    {
        plus20.push_back(std::fmod(v + 20.0, 360.0));
        plus180.push_back(std::fmod(v + 180.0, 360.0));
    }

    // Old smoke_bearing_comparison_identical.
    const auto w_same = lsdtt::compute_circular_wasserstein_distance(k_sample, k_sample, 360.0, 0.0);
    const auto k_same = lsdtt::compute_kuiper_two_sample_test(k_sample, k_sample, 360.0, 0.0);
    check_near("CT1 identical: W1 = 0", w_same.distance, 0.0, tol);
    check_near("CT1 identical: Kuiper V = 0", k_same.statistic, 0.0, tol);
    check_near("CT1 identical: Kuiper p = 1", k_same.p_value_approx, 1.0, tol);

    // Old smoke_bearing_comparison_rotated.
    const auto w_rot = lsdtt::compute_circular_wasserstein_distance(k_sample, plus20, 360.0, 0.0);
    check_near("CT2 rotated by 20: W1 = 20", w_rot.distance, 20.0, tol);

    // Old smoke_bearing_comparison_rotation_invariant.
    const auto w_ri = lsdtt::compute_rotation_invariant_circular_wasserstein_distance(
        k_sample, plus20, 360.0, 0.0);
    check_near("CT3 rotation-invariant: W1 = 0", w_ri.distance, 0.0, tol);
    check_near("CT3 rotation-invariant: optimal rotation = 340", w_ri.optimal_rotation, 340.0, tol);

    // Old smoke_bearing_comparison_axial.
    const auto w_ax = lsdtt::compute_axial_circular_wasserstein_distance(k_sample, plus180, 180.0, 0.0);
    check_near("CT4 +180 as axial (period 180): W1 = 0", w_ax.distance, 0.0, tol);

    // New closed forms.
    const auto w_wrap = lsdtt::compute_circular_wasserstein_distance({10.0}, {350.0}, 360.0, 0.0);
    check_near("CT5 point mass 10 vs 350: W1 = 20 (wraps, not 340)", w_wrap.distance, 20.0, tol);
    const auto w_ax90 = lsdtt::compute_axial_circular_wasserstein_distance({0.0}, {90.0}, 180.0, 0.0);
    check_near("CT6 axial point mass 0 vs 90: W1 = 90", w_ax90.distance, 90.0, tol);
    const auto w_ax10 = lsdtt::compute_axial_circular_wasserstein_distance({10.0}, {170.0}, 180.0, 0.0);
    check_near("CT7 axial point mass 10 vs 170: W1 = 20", w_ax10.distance, 20.0, tol);
    const auto k_disjoint = lsdtt::compute_kuiper_two_sample_test(
        {10.0, 20.0}, {190.0, 200.0}, 360.0, 0.0);
    check_near("CT8 disjoint arcs: Kuiper V = 1", k_disjoint.statistic, 1.0, tol);

    // The library wrapper calls exactly these functions.
    const auto c = lsdtt::compare_channel_bearing_samples(
        k_sample, plus20, channel_bearing_angular_type::directional);
    check_near("CT9 compare_channel_bearing_samples: W1 = 20", c.circular_wasserstein_deg, 20.0, tol);
    const auto c_ri = lsdtt::compare_channel_bearing_samples(
        k_sample, plus20, channel_bearing_angular_type::directional, true);
    check_near("CT10 wrapper rotation-invariant: optimal rotation = 340",
               c_ri.optimal_rotation_deg, 340.0, tol);
    const auto c_guard = lsdtt::compare_channel_bearing_samples(
        k_sample, plus20, channel_bearing_angular_type::directional, true, 63);
    check("CT11 rotation-invariant guard: 64 pairs > 63 is skipped with nan",
          c_guard.rotation_invariant_skipped && std::isnan(c_guard.circular_wasserstein_deg) &&
          std::isfinite(c_guard.kuiper_V));
}

// =============================================================================
// --group-summaries
// =============================================================================
// Angles on the half-open domains [0, 360) and [0, 180): a value just below
// the period (359.99999999999579, measured on a WV link) rounds to the period at
// the writers' 10 significant digits and must be written 0, never 360 / 180.
void run_angle_format_checks()
{
    using lsdtt::format_channel_bearing_angle;
    check("angle format: 359.99999999999579 -> 0 on [0, 360)",
          format_channel_bearing_angle(359.99999999999579, 360.0) == "0");
    check("angle format: 179.99999999999579 -> 0 on [0, 180)",
          format_channel_bearing_angle(179.99999999999579, 180.0) == "0");
    check("angle format: ordinary values unchanged",
          format_channel_bearing_angle(359.9999, 360.0) == "359.9999" &&
          format_channel_bearing_angle(207.3826203, 360.0) == num(207.3826203) &&
          format_channel_bearing_angle(0.0, 360.0) == "0" &&
          format_channel_bearing_angle(179.99999, 180.0) == "179.99999" &&
          format_channel_bearing_angle(std::nan(""), 360.0) == "nan");

    // The membership writer: the two links of wvbig_order that printed 360 / 180.
    std::vector<lsdtt::lsdtt_channel_link_bearing> recs(3);
    const double b[3] = {359.99999999999579, 359.99999999999721, 270.0};
    for (int i = 0; i < 3; ++i)
    {
        recs[i].upstream_junction = 10 + i;
        recs[i].downstream_junction = 20 + i;
        recs[i].stream_order = 1;
        recs[i].total_length_metres = 100.0;
        recs[i].bearing_deg = b[i];
    }
    lsdtt::channel_record_categories groups;
    groups.key = "strahler_order";
    groups.category = {1, 1, 1};
    groups.names[1] = "1";
    const std::string path = (std::filesystem::temp_directory_path() /
                              ("lsdtt_cbs_angle_format_" +
                               std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) +
                               ".csv")).string();
    lsdtt::write_channel_bearing_group_membership_csv(recs, groups, nullptr, {0, 0, 0}, {-1, -1, -1},
                                                      nullptr, nullptr, path);
    std::ifstream in(path);
    std::string header, l1, l2, l3;
    std::getline(in, header);
    std::getline(in, l1);
    std::getline(in, l2);
    std::getline(in, l3);
    in.close();
    std::filesystem::remove(path);
    check("membership writer: bearings just below 360 are written 0, axial 0",
          l1 == "10,20,1,100,0,0,1,all,0,-1,,nan" && l2 == "11,21,1,100,0,0,1,all,0,-1,,nan", l1 + " | " + l2);
    check("membership writer: other rows unchanged", l3 == "12,22,1,100,270,90,1,all,0,-1,,nan", l3);
}

void run_group_summary_checks()
{
    constexpr double tol = 1e-9;
    const double pi = std::acos(-1.0);

    const auto d1 = lsdtt::compute_channel_bearing_circular_mean({350.0, 10.0});
    check_circular_near("GS1 {350, 10} directional: mean 0 (not 180, not 360)",
                        d1.mean_deg, 0.0, 360.0, tol);
    check("GS1b mean is written in [0, 360)", d1.mean_deg >= 0.0 && d1.mean_deg < 360.0,
          num(d1.mean_deg));
    check_near("GS1c R = cos 10", d1.mean_resultant_length, std::cos(10.0 * pi / 180.0), tol);

    const auto a1 = lsdtt::compute_channel_bearing_axial_mean({170.0, 10.0});
    check_circular_near("GS2 {170, 10} axial: mean 0", a1.mean_deg, 0.0, 180.0, tol);
    check("GS2b axial mean in [0, 180)", a1.mean_deg >= 0.0 && a1.mean_deg < 180.0,
          num(a1.mean_deg));
    check_near("GS2c R = cos 20", a1.mean_resultant_length, std::cos(20.0 * pi / 180.0), tol);

    const auto a2 = lsdtt::compute_channel_bearing_axial_mean({0.0, 90.0});
    check_near("GS3 {0, 90} axial: R = 0", a2.mean_resultant_length, 0.0, 1e-12);
    check_nan("GS3b {0, 90} axial: mean nan", a2.mean_deg);

    const auto w1 = lsdtt::compute_channel_bearing_circular_mean({0.0, 90.0}, {1.0, 3.0});
    check_near("GS4 {0, 90} weights {1, 3}: mean atan2(3, 1)",
               w1.mean_deg, std::atan2(3.0, 1.0) * 180.0 / pi, tol);

    const std::vector<double> ten(10, 37.0);
    const auto s10 = lsdtt::summarise_channel_bearing_group(
        ten, std::vector<double>(10, 1.0), channel_bearing_angular_type::directional);
    check_near("GS5 10 identical: R = 1", s10.mean_resultant_length_count, 1.0, tol);
    check("GS5b 10 identical: Rayleigh p < 1e-4", s10.rayleigh_p_iid_approx < 1e-4,
          num(s10.rayleigh_p_iid_approx));
    check_near("GS5c 10 identical: circular SD = 0", s10.circular_sd_deg_count, 0.0, 1e-6);
    check_near("GS5d 10 identical: Rayleigh Z = n", s10.rayleigh_z, 10.0, tol);

    std::vector<double> d8;
    for (int rep = 0; rep < 25; ++rep)
        for (const double v : {0.0, 45.0, 90.0, 135.0}) d8.push_back(v);
    const std::vector<double> d8w(d8.size(), 1.0);
    const auto sd8a = lsdtt::summarise_channel_bearing_group(d8, d8w, channel_bearing_angular_type::axial);
    const auto sd8d = lsdtt::summarise_channel_bearing_group(d8, d8w, channel_bearing_angular_type::directional);
    check_near("GS6 {0,45,90,135}x25 axial: R = 0 (D8 cancellation)",
               sd8a.mean_resultant_length_count, 0.0, 1e-12);
    check_nan("GS6b axial mean nan", sd8a.circular_mean_deg_count);
    check_nan("GS6c axial circular SD nan where the mean is undefined", sd8a.circular_sd_deg_count);
    // Directional: 0, 45, 90, 135 do NOT cancel as directions (they span a half
    // circle); R = |1 + e^{i45} + e^{i90} + e^{i135}| / 4.
    const double rdir = std::hypot(1.0 + std::cos(pi / 4) + 0.0 + std::cos(3 * pi / 4),
                                   0.0 + std::sin(pi / 4) + 1.0 + std::sin(3 * pi / 4)) / 4.0;
    check_near("GS6d directional R of {0,45,90,135} = |sum|/4", sd8d.mean_resultant_length_count,
               rdir, tol);
    check_near("GS6e fraction_on_d8_multiples = 1", sd8a.fraction_on_d8_multiples, 1.0, tol);
    // The D8 directions that DO cancel directionally: all eight.
    std::vector<double> d8full;
    for (int k = 0; k < 8; ++k) d8full.push_back(45.0 * k);
    const auto d8f = lsdtt::compute_channel_bearing_circular_mean(d8full);
    check_near("GS6f all eight D8 directions: directional R = 0", d8f.mean_resultant_length,
               0.0, 1e-12);
    check_near("GS6g half-degree off D8 is not on a multiple",
               lsdtt::channel_bearing_fraction_on_d8_multiples({0.0, 0.6, 44.6, 90.4}), 0.75, tol);

    // Kish n_eff and total length.
    const auto kish = lsdtt::summarise_channel_bearing_group(
        {10.0, 20.0, 30.0}, {1.0, 1.0, 2.0}, channel_bearing_angular_type::directional);
    check_near("GS7 Kish n_eff = (sum w)^2 / sum w^2 = 16/6", kish.n_eff_length, 16.0 / 6.0, tol);
    check_near("GS7b total length = sum w", kish.total_length_m, 4.0, tol);

    // Axial SD is halved from the doubled circle.
    const auto ax_sd = lsdtt::summarise_channel_bearing_group(
        {170.0, 10.0}, {1.0, 1.0}, channel_bearing_angular_type::axial);
    const double expected_sd =
        std::sqrt(-2.0 * std::log(std::cos(20.0 * pi / 180.0))) * 180.0 / pi / 2.0;
    check_near("GS8 axial circular SD = sqrt(-2 ln R2) / 2", ax_sd.circular_sd_deg_count,
               expected_sd, tol);

    // Rose.
    auto bin_of = [](const std::vector<lsdtt::channel_bearing_rose_bin>& bins) {
        for (const auto& b : bins) if (b.count > 0) return b.bin_centre_deg;
        return -1.0;
    };
    const auto r355 = lsdtt::channel_bearing_rose_histogram(
        {355.0}, {}, channel_bearing_angular_type::directional, 10.0);
    check("RS1 w = 10: 36 directional bins", r355.size() == 36, std::to_string(r355.size()));
    check_near("RS2 355 -> bin centre 0", bin_of(r355), 0.0, tol);
    const auto r5 = lsdtt::channel_bearing_rose_histogram(
        {5.0}, {}, channel_bearing_angular_type::directional, 10.0);
    check_near("RS3 5.0 -> bin centre 10 (left-closed)", bin_of(r5), 10.0, tol);
    const auto r175 = lsdtt::channel_bearing_rose_histogram(
        {175.0}, {}, channel_bearing_angular_type::axial, 10.0);
    check("RS4 w = 10: 18 axial bins", r175.size() == 18, std::to_string(r175.size()));
    check_near("RS5 axial 175 -> bin centre 0 (wraps)", bin_of(r175), 0.0, tol);
    check_near("RS5b bin 0 lower edge = -w/2", r175.front().bin_lower_deg, -5.0, tol);

    const std::vector<double> many = {1, 44, 91, 133, 222, 359.9, 180, 270.2, 15, 97};
    const std::vector<double> lens = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    for (const auto type : {channel_bearing_angular_type::directional,
                            channel_bearing_angular_type::axial})
    {
        const auto bins = lsdtt::channel_bearing_rose_histogram(many, lens, type, 15.0);
        std::size_t c = 0;
        double l = 0.0;
        for (const auto& b : bins) { c += b.count; l += b.length_m; }
        check("RS6 " + lsdtt::channel_bearing_angular_type_name(type) +
              ": bin counts sum to n and lengths to total",
              c == many.size() && std::abs(l - 55.0) < tol);
    }
    check("RS7 width 7 rejected (180/7 not an integer)",
          !lsdtt::channel_bearing_rose_bin_width_is_valid(7.0));
    check("RS7b width 7.5 accepted", lsdtt::channel_bearing_rose_bin_width_is_valid(7.5));

    // Category summaries: excluded records are not summarised; "all" pools.
    lsdtt::channel_record_categories g;
    g.key = "strahler_order";
    g.category = {1, 1, 2, -1, 2};
    g.names = {{1, "1"}, {2, "2"}};
    lsdtt::channel_bearing_statistics_options opt;
    opt.min_group_size = 2;
    const auto rows = lsdtt::summarise_channel_link_bearings_by_category(
        {0.0, 10.0, 90.0, 45.0, 92.0}, {1, 1, 1, 1, 1}, {7, 7, 8, -1, 9}, g, nullptr, opt);
    check("CS1 2 groups x 2 angular types = 4 rows", rows.size() == 4, std::to_string(rows.size()));
    bool n_ok = true;
    for (const auto& r : rows)
    {
        n_ok = n_ok && r.summary.n == 2 && r.stratum_key == "none" && r.stratum_value == "all";
        if (r.group_value == "1") n_ok = n_ok && r.n_blocks == 1;
        if (r.group_value == "2") n_ok = n_ok && r.n_blocks == 2;
    }
    check("CS2 excluded record dropped; n_blocks counts distinct blocks", n_ok);

    // Combined key (order * 100000 + code) needs nothing special.
    lsdtt::channel_record_categories combo;
    combo.key = "order_x_unit";
    combo.category = {1 * 100000 + 7, 2 * 100000 + 7, 2 * 100000 + 7};
    combo.names = {{100007, "1|sandstone"}, {200007, "2|sandstone"}};
    const auto crow = lsdtt::summarise_channel_link_bearings_by_category(
        {0.0, 0.0, 10.0}, {1, 1, 1}, {-1, -1, -1}, combo, nullptr, opt);
    check("CS3 combined key: rows labelled by name, ascending category",
          crow.size() == 4 && crow.front().group_value == "1|sandstone" &&
          crow.back().group_value == "2|sandstone");
}

// =============================================================================
// --comparisons
// =============================================================================
void run_comparison_checks()
{
    constexpr double tol = 1e-12;
    const auto holm = lsdtt::holm_adjusted_p_values({0.01, 0.04, 0.03});
    check_near("HO1 Holm 0.01 -> 0.03", holm[0], 0.03, tol);
    check_near("HO2 Holm 0.04 -> 0.06", holm[1], 0.06, tol);
    check_near("HO3 Holm 0.03 -> 0.06 (monotone)", holm[2], 0.06, tol);
    const auto holm_nan = lsdtt::holm_adjusted_p_values({0.2, std::nan(""), 0.6});
    check("HO4 nan stays nan and is not counted", std::isnan(holm_nan[1]) &&
          std::abs(holm_nan[0] - 0.4) < tol && std::abs(holm_nan[2] - 0.6) < tol);

    // Four groups of 5 records, plus one group of 2 below the minimum.
    std::vector<double> bearings;
    lsdtt::channel_record_categories g;
    g.key = "basin";
    for (int grp = 0; grp < 4; ++grp)
    {
        for (int i = 0; i < 5; ++i)
        {
            bearings.push_back(20.0 * grp + 3.0 * i);
            g.category.push_back(grp);
        }
        g.names[grp] = "b" + std::to_string(grp);
    }
    bearings.push_back(300.0); g.category.push_back(9);
    bearings.push_back(310.0); g.category.push_back(9);
    g.names[9] = "small";
    const std::vector<std::int64_t> blocks(bearings.size(), -1);

    lsdtt::channel_bearing_statistics_options opt;
    opt.min_group_size = 5;
    opt.do_axial = false;

    const auto pw = lsdtt::compare_channel_link_bearings_by_category(
        bearings, blocks, g, nullptr, opt);
    check("CO1 pairwise, K = 4 eligible: 6 rows", pw.rows.size() == 6,
          std::to_string(pw.rows.size()));
    bool no_small = true;
    for (const auto& r : pw.rows)
        no_small = no_small && r.group_a != "small" && r.group_b != "small" &&
                   r.n_comparisons_in_family == 6 && r.comparison_mode == "pairwise";
    check("CO2 the group below min size takes no part; family size 6", no_small);

    opt.comparison_mode = lsdtt::channel_bearing_comparison_mode::one_vs_rest;
    const auto ovr = lsdtt::compare_channel_link_bearings_by_category(
        bearings, blocks, g, nullptr, opt);
    check("CO3 one_vs_rest: 4 rows", ovr.rows.size() == 4, std::to_string(ovr.rows.size()));
    check("CO3b rest includes the small group (n_b = 17)",
          !ovr.rows.empty() && ovr.rows.front().n_b == 17 && ovr.rows.front().group_b == "rest");

    const std::vector<double> reference = {0.0, 5.0, 10.0};
    opt.comparison_mode = lsdtt::channel_bearing_comparison_mode::pairwise;
    const auto withref = lsdtt::compare_channel_link_bearings_by_category(
        bearings, blocks, g, nullptr, opt, &reference, "reference");
    std::size_t n_ref = 0;
    for (const auto& r : withref.rows) if (r.comparison_mode == "reference") ++n_ref;
    check("CO4 a reference adds K = 4 rows", withref.rows.size() == 10 && n_ref == 4,
          std::to_string(withref.rows.size()) + " rows, " + std::to_string(n_ref) + " reference");

    const auto summaries = lsdtt::summarise_channel_link_bearings_by_category(
        bearings, std::vector<double>(bearings.size(), 1.0), blocks, g, nullptr, opt);
    bool small_flag = false;
    for (const auto& r : summaries)
        if (r.group_value == "small") small_flag = !r.meets_min_group_size && r.summary.n == 2;
    check("CO5 group below min size: meets_min_group_size = 0", small_flag);

    // Stratified: comparisons only within a stratum, plus the pooled "all".
    lsdtt::channel_record_categories strata;
    strata.key = "strahler_order";
    for (std::size_t i = 0; i < bearings.size(); ++i)
        strata.category.push_back(static_cast<std::int32_t>(i % 2) + 1);
    strata.names = {{1, "1"}, {2, "2"}};
    opt.min_group_size = 2;
    const auto st = lsdtt::compare_channel_link_bearings_by_category(
        bearings, blocks, g, &strata, opt);
    std::set<std::string> stratum_values;
    for (const auto& r : st.rows) stratum_values.insert(r.stratum_value);
    check("CO6 stratified rows carry strata 1, 2 and all",
          stratum_values == std::set<std::string>{"1", "2", "all"});

    // Reference-CSV reader.
    const std::filesystem::path dir =
        std::filesystem::temp_directory_path() / "lsdtt_test_channel_bearing_statistics";
    std::filesystem::create_directories(dir);
    const std::string csv = (dir / "ref.csv").string();
    {
        std::ofstream os(csv);
        os << "bearing_deg,stream_order\n"
           << "10,1\n" << "20,2\n" << "abc,2\n" << ",2\n" << "40,2\n" << "50,x\n";
    }
    std::string message;
    try
    {
        (void)lsdtt::read_channel_bearing_reference_csv(csv, "x");
    }
    catch (const std::exception& e)
    {
        message = e.what();
    }
    check("RF1 missing column: exact message", message == "Column 'x' not found in " + csv,
          message);
    const auto all = lsdtt::read_channel_bearing_reference_csv(csv, "bearing_deg");
    check("RF2 non-numeric and blank cells skipped: 4 values", all.size() == 4,
          std::to_string(all.size()));
    const auto filtered =
        lsdtt::read_channel_bearing_reference_csv(csv, "bearing_deg", "stream_order", 2);
    check("RF3 order filter keeps only matching rows: {20, 40}",
          filtered == std::vector<double>{20.0, 40.0});
    std::filesystem::remove_all(dir);
}

// =============================================================================
// --junction-angle-stats
// =============================================================================
void run_junction_angle_checks()
{
    constexpr double tol = 1e-12;
    const auto same = lsdtt::compare_junction_angle_samples({60.0, 90.0}, {60.0, 90.0});
    check_near("JA1 {60, 90} vs itself: KS D = 0", same.ks_D, 0.0, tol);
    check_near("JA2 {60, 90} vs itself: linear W1 = 0", same.linear_wasserstein_deg, 0.0, tol);
    const auto apart = lsdtt::compare_junction_angle_samples({30.0}, {90.0});
    check_near("JA3 {30} vs {90}: linear W1 = 60", apart.linear_wasserstein_deg, 60.0, tol);
    // Linear, not circular: 10 and 350 are 340 apart as sectors.
    const auto lin = lsdtt::compare_junction_angle_samples({10.0}, {350.0});
    check_near("JA4 sectors 10 vs 350 are 340 apart (linear)", lin.linear_wasserstein_deg,
               340.0, tol);

    const std::vector<double> alpha = {80, 70, 100};
    const std::vector<double> b1 = {140, 200, 120};
    const std::vector<double> b2 = {140, 90, 140};
    const auto m = lsdtt::make_junction_angle_measurands(alpha, b1, b2);
    const auto swapped = lsdtt::make_junction_angle_measurands(alpha, b2, b1);
    check("JA5 beta min/max invariant to swapping donor labels",
          m.beta_min == swapped.beta_min && m.beta_max == swapped.beta_max &&
          m.beta_min == std::vector<double>({140, 90, 120}) &&
          m.beta_max == std::vector<double>({140, 200, 140}));

    const auto s = lsdtt::summarise_junction_angle_group({10.0, 20.0, 30.0, 40.0});
    check("JA6 linear summary: n 4, mean 25, median 25",
          s.n == 4 && std::abs(s.mean_deg - 25.0) < tol && std::abs(s.median_deg - 25.0) < tol);
}

// =============================================================================
// Trellis fixture
// =============================================================================
constexpr int    k_trellis_rows  = 300;
constexpr int    k_trellis_cols  = 160;
constexpr double k_trellis_dx    = 10.0;
constexpr double k_trellis_slope = 0.02;
constexpr double k_trellis_relief = 50.0;
constexpr int    k_trellis_threshold = 8;

struct trellis_valley
{
    int col;
    int west_px;
    int east_px;
};

const std::vector<trellis_valley> k_valleys = {
    {6, 6, 34}, {46, 6, 34}, {114, 34, 6}, {154, 34, 6}};

// Flank height at cross-valley column c.
double trellis_flank_height(int c)
{
    for (const trellis_valley& v : k_valleys)
    {
        if (c >= v.col - v.west_px && c <= v.col)
        {
            const double d = v.col - c;
            return k_trellis_relief * 0.5 * (1.0 - std::cos(std::acos(-1.0) * d / v.west_px));
        }
        if (c >= v.col && c <= v.col + v.east_px)
        {
            const double d = c - v.col;
            return k_trellis_relief * 0.5 * (1.0 - std::cos(std::acos(-1.0) * d / v.east_px));
        }
    }
    throw std::runtime_error("trellis_flank_height: column outside every valley span");
}

// transposed = false: N-S valleys draining south (boundary nnon).
// transposed = true : the same surface with rows and columns swapped, so the
//                     valleys run E-W and drain west (boundary nnno).
lsdtt::lsdtt_raster_cpu make_trellis_dem(bool transposed)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = transposed ? k_trellis_cols : k_trellis_rows;
    meta.n_cols = transposed ? k_trellis_rows : k_trellis_cols;
    meta.cell_size = k_trellis_dx;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> z(static_cast<std::size_t>(meta.n_rows) * meta.n_cols);
    for (int r = 0; r < meta.n_rows; ++r)
    {
        for (int c = 0; c < meta.n_cols; ++c)
        {
            double down_valley;   // distance up-valley from the outlet edge
            double flank;
            if (!transposed)
            {
                down_valley = (meta.n_rows - 1 - r) * k_trellis_dx;   // rises northward
                flank = trellis_flank_height(c);
            }
            else
            {
                down_valley = c * k_trellis_dx;                         // rises eastward
                flank = trellis_flank_height(r);
            }
            z[static_cast<std::size_t>(r) * meta.n_cols + c] =
                static_cast<float>(k_trellis_slope * down_valley + flank);
        }
    }
    std::vector<std::uint8_t> mask(z.size(), 1);
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta);
}

struct routed_network
{
    lsdtt::lsdtt_raster_cpu dem;
    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    std::unique_ptr<lsdtt::lsdtt_junction_network_base> network;
};

routed_network route(lsdtt::lsdtt_raster_cpu dem, const std::string& boundary, int threshold)
{
    routed_network rn;
    rn.dem = lsdtt::fill_pits(dem, 0.0001f, boundary, 1);

    lsdtt::lsdtt_flowinfo_options fi;
    fi.boundary_conditions       = boundary;
    fi.interior_nodata_is_outlet = true;
    fi.require_acyclic_graph     = true;
    fi.build_topological_order   = true;
    fi.build_s_vector            = true;
    fi.n_threads                 = 1;
    rn.flow = lsdtt::make_flow_base(
        rn.dem.get_data_array(), rn.dem.get_mask_array(), rn.dem.get_meta(), fi);

    const std::vector<std::int64_t> sources = rn.flow->get_sources_index_threshold(threshold);
    lsdtt::lsdtt_junction_network_options jn;
    jn.build_topological_order = true;
    jn.build_s_vector          = true;
    rn.network = lsdtt::make_junction_network_base(*rn.flow, sources, jn);
    return rn;
}

std::vector<double> bearings_of_order(const std::vector<lsdtt::lsdtt_channel_link_bearing>& recs,
                                      int order)
{
    std::vector<double> out;
    for (const auto& r : recs) if (r.stream_order == order) out.push_back(r.bearing_deg);
    return out;
}

double median_of(std::vector<double> v)
{
    if (v.empty()) return std::nan("");
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

void run_trellis_checks()
{
    lsdtt::lsdtt_link_sampling_options sampling;   // full_link, min_nodes 2

    // ---------------------------------------------------------------- N-S trellis
    const routed_network ns = route(make_trellis_dem(false), "nnon", k_trellis_threshold);
    const auto links = lsdtt::calculate_channel_link_bearings(*ns.flow, *ns.network, {}, sampling);
    const auto angles = lsdtt::calculate_junction_angles(*ns.flow, *ns.network, {}, sampling);

    std::cout << "trellis: " << ns.network->n_junctions() << " junctions, " << links.size()
              << " links, " << angles.size() << " two-donor junctions\n";
    check("TR0 the trellis has links and two-donor junctions",
          links.size() > 1000 && angles.size() > 500);

    // sampled_flats are exactly the nodes the bearing was fitted on.
    bool flats_ok = true;
    for (const auto& r : links)
    {
        const auto traced =
            lsdtt::trace_channel_link_between_junctions(*ns.flow, *ns.network, r.upstream_junction);
        const auto sampled = lsdtt::sample_traced_channel_link(
            traced, sampling, lsdtt::lsdtt_link_sample_anchor::upstream_end);
        flats_ok = flats_ok && r.sampled_flats.size() == r.n_nodes_used &&
                   r.sampled_flats == sampled.flats;
    }
    check("SF1 link sampled_flats: count = n_nodes_used and equal to the fitted sample", flats_ok);
    bool arm_ok = true;
    for (const auto& a : angles)
    {
        arm_ok = arm_ok && a.donor_1_sampled_flats.size() == a.donor_1_nodes_used &&
                 a.donor_2_sampled_flats.size() == a.donor_2_nodes_used &&
                 a.receiver_sampled_flats.size() == a.receiver_nodes_used &&
                 !a.donor_1_sampled_flats.empty() &&
                 a.donor_1_sampled_flats.back() == ns.network->flat_of_junction(a.junction) &&
                 a.receiver_sampled_flats.front() == ns.network->flat_of_junction(a.junction);
    }
    check("SF2 junction arms: sampled_flats count = nodes_used; donors end and receiver "
          "starts at the junction", arm_ok);
    // A non-default sampling: the record still holds the nodes it was fitted on.
    lsdtt::lsdtt_link_sampling_options fixed;
    fixed.mode = lsdtt::lsdtt_link_sampling_mode::fixed_n_pixels;
    fixed.fixed_n_pixels = 4;
    const auto links_fixed = lsdtt::calculate_channel_link_bearings(
        *ns.flow, *ns.network, {}, fixed, lsdtt::lsdtt_link_sample_anchor::downstream_end);
    bool fixed_ok = !links_fixed.empty();
    for (const auto& r : links_fixed)
        fixed_ok = fixed_ok && r.sampled_flats.size() == r.n_nodes_used && r.n_nodes_used <= 4 &&
                   r.sampled_flats.back() == r.downstream_flat;
    check("SF3 fixed_n_pixels 4 at downstream_end: sampled_flats hold the last <= 4 nodes",
          fixed_ok);

    const std::vector<double> o1 = bearings_of_order(links, 1);
    const std::vector<double> o2 = bearings_of_order(links, 2);
    std::size_t max_order = 0;
    for (const auto& r : links) max_order = std::max<std::size_t>(max_order, r.stream_order);
    std::cout << "trellis: order-1 links " << o1.size() << ", order-2 links " << o2.size()
              << ", max order " << max_order << "\n";
    check("TR1 orders are 1 and 2 only", max_order == 2 && !o1.empty() && !o2.empty());

    const auto o2_ax = lsdtt::compute_channel_bearing_axial_mean(o2);
    const auto o2_dir = lsdtt::compute_channel_bearing_circular_mean(o2);
    check_circular_near("TR2 order-2 (trunks) axial mean ~ 0", o2_ax.mean_deg, 0.0, 180.0, 2.0);
    check("TR3 order-2 axial R > 0.98", o2_ax.mean_resultant_length > 0.98,
          num(o2_ax.mean_resultant_length));
    check_circular_near("TR4 order-2 directional mean ~ 180 (draining south)",
                        o2_dir.mean_deg, 180.0, 360.0, 2.0);

    const auto o1_ax = lsdtt::compute_channel_bearing_axial_mean(o1);
    const auto o1_dir = lsdtt::compute_channel_bearing_circular_mean(o1);
    check_circular_near("TR5 order-1 axial mean ~ 90 (tributaries E-W)", o1_ax.mean_deg,
                        90.0, 180.0, 5.0);
    std::size_t n_trib = 0;
    std::size_t n_head = 0;
    for (const double b : o1)
    {
        const double ax = std::fmod(b, 180.0);
        if (std::abs(ax - 90.0) < 45.0) ++n_trib; else ++n_head;
    }
    const double expected_r =
        (static_cast<double>(n_trib) - static_cast<double>(n_head)) /
        static_cast<double>(n_trib + n_head);
    std::cout << "trellis: order-1 E-W tributaries " << n_trib << ", N-S head links "
              << n_head << "\n";
    check_near("TR6 order-1 axial R = (n_trib - n_head) / (n_trib + n_head)",
               o1_ax.mean_resultant_length, expected_r, 0.05);
    check("TR7 order-1 directional R < 0.2 (W- and E-flowing tributaries cancel)",
          o1_dir.mean_resultant_length < 0.2, num(o1_dir.mean_resultant_length));

    std::vector<double> alpha;
    for (const auto& a : angles) alpha.push_back(a.donor_opening_angle_deg);
    check_near("TR8 junction alpha median ~ 90 (trunk 180, tributary 90 or 270)",
               median_of(alpha), 90.0, 5.0);

    const auto k12 = lsdtt::compare_channel_bearing_samples(
        o1, o2, channel_bearing_angular_type::axial);
    check("TR9 Kuiper order 1 vs order 2, axial: V > 0.5", k12.kuiper_V > 0.5, num(k12.kuiper_V));

    // ---------------------------------------------------------------- rotated
    const routed_network ew = route(make_trellis_dem(true), "nnno", k_trellis_threshold);
    const auto links_ew = lsdtt::calculate_channel_link_bearings(*ew.flow, *ew.network, {}, sampling);
    const std::vector<double> ew2 = bearings_of_order(links_ew, 2);
    const auto ew2_ax = lsdtt::compute_channel_bearing_axial_mean(ew2);
    const auto ew2_dir = lsdtt::compute_channel_bearing_circular_mean(ew2);
    check_circular_near("RT1 rotated order-2 axial mean ~ 90", ew2_ax.mean_deg, 90.0, 180.0, 2.0);
    check_circular_near("RT2 rotated order-2 directional mean ~ 270 (draining west)",
                        ew2_dir.mean_deg, 270.0, 360.0, 2.0);
    const auto w_ax = lsdtt::compare_channel_bearing_samples(
        o2, ew2, channel_bearing_angular_type::axial);
    const auto w_dir = lsdtt::compare_channel_bearing_samples(
        o2, ew2, channel_bearing_angular_type::directional);
    check_near("RT3 axial W1, trellis order 2 vs rotated order 2 = 90", w_ax.circular_wasserstein_deg,
               90.0, 3.0);
    check_near("RT4 directional W1 (180 vs 270) = 90", w_dir.circular_wasserstein_deg, 90.0, 3.0);
}

// =============================================================================
// --grouping-keys
// =============================================================================
void run_grouping_checks()
{
    lsdtt::lsdtt_link_sampling_options sampling;
    const routed_network ns = route(make_trellis_dem(false), "nnon", k_trellis_threshold);
    const auto links = lsdtt::calculate_channel_link_bearings(*ns.flow, *ns.network, {}, sampling);
    const auto& meta = ns.dem.get_meta();

    // Region: north half (rows < 150) = index 0, south half = index 1.
    lsdtt::channel_network_zone_labels halves;
    halves.key_name = "region";
    halves.values = {"north", "south"};
    halves.value_index.resize(static_cast<std::size_t>(meta.n_rows) * meta.n_cols);
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            halves.value_index[static_cast<std::size_t>(r) * meta.n_cols + c] = r < 150 ? 0 : 1;

    const auto zones = lsdtt::assign_zones_to_channel_links(*ns.flow, links, halves);
    const auto regions = lsdtt::categories_from_zone_assignment("region", zones, halves, 0.0);
    std::size_t n_north = 0, n_south = 0, n_none = 0;
    for (const std::int32_t c : regions.category)
    {
        if (c == 0) ++n_north; else if (c == 1) ++n_south; else ++n_none;
    }
    std::cout << "grouping: north " << n_north << ", south " << n_south << ", none " << n_none
              << " of " << links.size() << " links\n";
    check("GK1 n_north + n_south + n_none = n_links",
          n_north + n_south + n_none == links.size() && n_north > 0 && n_south > 0);

    // Purity of the links crossing the split equals the length fraction of
    // their SAMPLED nodes (each weighted by its D8 step to its receiver).
    bool purity_ok = true;
    std::size_t n_crossing = 0;
    std::size_t n_crossing_below_06 = 0;
    const double dx = ns.flow->cell_size();
    for (std::size_t i = 0; i < links.size(); ++i)
    {
        double north = 0.0, south = 0.0;
        for (const std::int64_t f : links[i].sampled_flats)
        {
            double w = 0.0;
            if (!ns.flow->is_outlet(f))
            {
                const std::int64_t rc = ns.flow->receiver_of(f);
                const bool diag = ns.flow->row_from_flat(rc) != ns.flow->row_from_flat(f) &&
                                  ns.flow->col_from_flat(rc) != ns.flow->col_from_flat(f);
                w = rc == f ? 0.0 : (diag ? dx * std::sqrt(2.0) : dx);
            }
            (ns.flow->row_from_flat(f) < 150 ? north : south) += w;
        }
        if (north > 0.0 && south > 0.0)
        {
            ++n_crossing;
            const double expected = std::max(north, south) / (north + south);
            purity_ok = purity_ok && std::abs(zones.purity[i] - expected) < 1e-12;
            if (expected < 0.6) ++n_crossing_below_06;
        }
        else
        {
            purity_ok = purity_ok && zones.purity[i] == 1.0;
        }
    }
    check("GK2 purity = majority sampled length fraction (" + std::to_string(n_crossing) +
          " links cross the split)", purity_ok && n_crossing > 0);
    const auto regions_06 = lsdtt::categories_from_zone_assignment("region", zones, halves, 0.6);
    check("GK3 min purity 0.6 excludes exactly the crossing links below 0.6",
          regions_06.n_excluded() == n_crossing_below_06,
          std::to_string(regions_06.n_excluded()) + " vs " + std::to_string(n_crossing_below_06));

    // Order-2 axial W1 between halves < 1 degree.
    std::vector<double> north2, south2;
    for (std::size_t i = 0; i < links.size(); ++i)
    {
        if (links[i].stream_order != 2) continue;
        if (regions.category[i] == 0) north2.push_back(links[i].bearing_deg);
        if (regions.category[i] == 1) south2.push_back(links[i].bearing_deg);
    }
    const auto w = lsdtt::compare_channel_bearing_samples(
        north2, south2, channel_bearing_angular_type::axial);
    check("GK4 order-2 axial W1 north vs south < 1 deg", w.circular_wasserstein_deg < 1.0,
          num(w.circular_wasserstein_deg));

    // Order pooling and donor pairs.
    std::vector<std::int32_t> orders;
    std::vector<std::int32_t> upstream;
    for (const auto& r : links)
    {
        orders.push_back(r.stream_order);
        upstream.push_back(r.upstream_junction);
    }
    const auto pooled = lsdtt::strahler_order_categories(orders, 2);
    check("GK5 max_order_group 2: orders >= 2 pool into '2+'",
          pooled.names.size() == 2 && pooled.names.at(2) == "2+" && pooled.names.at(1) == "1");
    const auto pairs = lsdtt::donor_order_pair_categories({2, 1, 3}, {1, 1, 0});
    check("GK6 donor_order_pair: '1_2', '1_1', excluded",
          pairs.category == std::vector<std::int32_t>({1002, 1001, -1}) &&
          pairs.names.at(1002) == "1_2");

    // Order-2 basins. The open south edge adds one per valley: the last row's
    // flank drains ALONG row 299 and enters the outlet pixel beside the trunk,
    // so the outlet junction is order 3 and each valley has two order-2 mouths:
    // the trunk (row 298) and that edge channel (row 299). Re-derived from the
    // geometry; the trunks are the whole valleys.
    const auto basin2 = lsdtt::order_k_basin_of_junctions(*ns.network, 2);
    std::set<std::int32_t> distinct;
    for (const std::int32_t b : basin2) if (b >= 0) distinct.insert(b);
    std::set<std::int32_t> trunk_basins;
    std::size_t edge_mouths = 0;
    for (const std::int32_t b : distinct)
    {
        const int row = ns.network->row_of_junction(b);
        if (row == k_trellis_rows - 1) ++edge_mouths;
    }
    for (std::int32_t j = 0; j < ns.network->n_junctions(); ++j)
    {
        const int col = ns.network->col_of_junction(j);
        const bool on_trunk = std::any_of(k_valleys.begin(), k_valleys.end(),
            [&](const trellis_valley& v) { return v.col == col; });
        if (on_trunk && ns.network->stream_order_of_junction(j) == 2 &&
            ns.network->row_of_junction(j) < k_trellis_rows - 2)
            trunk_basins.insert(basin2[static_cast<std::size_t>(j)]);
    }
    check("GK7 the trunk junctions fall in one order-2 basin per whole valley (4)",
          trunk_basins.size() == k_valleys.size() && trunk_basins.count(-1) == 0,
          std::to_string(trunk_basins.size()));
    check("GK7a the other order-2 basins are the 4 south-edge channels (row 299)",
          distinct.size() == 2 * k_valleys.size() && edge_mouths == k_valleys.size(),
          std::to_string(distinct.size()) + " basins, " + std::to_string(edge_mouths) +
          " with an edge mouth");
    bool order_ok = true;
    for (const std::int32_t b : distinct)
        order_ok = order_ok && ns.network->stream_order_of_junction(b) == 2;
    check("GK7b each order-2 basin mouth is an order-2 junction", order_ok);
    // Order-1 basins: every order-1 junction is its own basin's member and
    // order-2 junctions have none.
    const auto basin1 = lsdtt::order_k_basin_of_junctions(*ns.network, 1);
    bool b1_ok = true;
    for (std::int32_t j = 0; j < ns.network->n_junctions(); ++j)
    {
        const int o = ns.network->stream_order_of_junction(j);
        if (o > 1) b1_ok = b1_ok && basin1[static_cast<std::size_t>(j)] == -1;
        else b1_ok = b1_ok && basin1[static_cast<std::size_t>(j)] >= 0;
    }
    check("GK8 k = 1: order > 1 junctions have no order-1 basin", b1_ok);

    // Innermost wins: an outer outlet and a junction half way up its trunk.
    std::int32_t outer = -1;
    for (std::int32_t i = 0; i < ns.network->n_outlet_junctions(); ++i)
    {
        const std::int32_t o = ns.network->outlet_junction(i);
        if (ns.network->col_of_junction(o) == k_valleys.front().col) outer = o;
    }
    std::int32_t inner = -1;
    for (std::int32_t j = 0; j < ns.network->n_junctions(); ++j)
    {
        if (ns.network->col_of_junction(j) == k_valleys.front().col &&
            ns.network->row_of_junction(j) == k_trellis_rows / 2 &&
            ns.network->stream_order_of_junction(j) == 2)
            inner = j;
    }
    check("GK9 found an outer outlet and an inner trunk junction", outer >= 0 && inner >= 0);
    if (outer >= 0 && inner >= 0)
    {
        const auto ob = lsdtt::outlet_basin_of_junctions(*ns.network, {inner, outer});
        bool inner_ok = ob[static_cast<std::size_t>(inner)] == inner;
        bool outer_ok = true;
        std::size_t n_inner = 0;
        std::size_t n_outer = 0;
        for (std::int32_t j = 0; j < ns.network->n_junctions(); ++j)
        {
            const bool up_of_inner = j == inner || ns.network->is_upslope_of(j, inner);
            const bool up_of_outer = j == outer || ns.network->is_upslope_of(j, outer);
            const std::int32_t b = ob[static_cast<std::size_t>(j)];
            if (up_of_inner) { inner_ok = inner_ok && b == inner; ++n_inner; }
            else if (up_of_outer) { outer_ok = outer_ok && b == outer; ++n_outer; }
            else outer_ok = outer_ok && b == -1;
        }
        check("GK10 innermost wins: junctions upslope of the inner outlet get it (" +
              std::to_string(n_inner) + "), the rest of the valley the outer (" +
              std::to_string(n_outer) + "), elsewhere none",
              inner_ok && outer_ok && n_inner > 0 && n_outer > 0);
        const auto bc = lsdtt::basin_categories(upstream, ob);
        check("GK11 basin categories: two basins, the other valleys excluded",
              bc.names.size() == 2 && bc.n_excluded() > 0);
    }

    // The grouped summary on the trellis: stratified by region.
    const auto groups = lsdtt::strahler_order_categories(orders, 0);
    const auto blocks = lsdtt::block_ids_of_records(upstream, basin2);
    std::vector<double> lengths;
    for (const auto& r : links) lengths.push_back(r.total_length_metres);
    lsdtt::channel_bearing_statistics_options opt;
    const auto rows = lsdtt::summarise_channel_link_bearings_by_category(
        [&] { std::vector<double> b; for (const auto& r : links) b.push_back(r.bearing_deg); return b; }(),
        lengths, blocks, groups, &regions, opt);
    std::size_t blocks_o2 = 0;
    for (const auto& r : rows)
        if (r.stratum_value == "all" && r.group_value == "2") blocks_o2 = r.n_blocks;
    check("GK12 order-2 links span every order-2 block (" + std::to_string(distinct.size()) + ")",
          blocks_o2 == distinct.size(), std::to_string(blocks_o2));
}

// =============================================================================
// --isotropic-null
// =============================================================================
// Seed choice, measured 2026-10-17 (257 x 257 at 30 m, roughness 0.55, sd 100 m,
// threshold 20 px, full_link): all-link axial R-bar for seeds 1-8 was 0.0751,
// 0.0063, 0.0523, 0.0191, 0.0457, 0.0272, 0.0369, 0.0149. Seed 1 sits 1.998x
// below the 0.15 bound (under the 2x margin asked for), so seed 2 is used.
// Directional R-bar on the same surfaces was 0.09-0.36: a D8 network on an
// isotropic surface is not directionally uniform, which is why the axial form is
// the one bounded here.
constexpr int k_null_seed = 2;

void run_isotropic_null_checks(int seed)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 257;
    meta.n_cols = 257;
    meta.cell_size = 30.0;
    meta.set_default_synthetic_georeference();
    const lsdtt::lsdtt_raster_cpu dem = lsdtt::generate_diamond_square_surface(
        meta, seed, 0.55, 1.0, 500.0, 100.0, 1);
    const routed_network rn = route(dem, "oooo", 20);
    lsdtt::lsdtt_link_sampling_options sampling;
    const auto links = lsdtt::calculate_channel_link_bearings(*rn.flow, *rn.network, {}, sampling);
    std::vector<double> b;
    for (const auto& r : links) b.push_back(r.bearing_deg);
    const auto ax = lsdtt::compute_channel_bearing_axial_mean(b);
    const auto dir = lsdtt::compute_channel_bearing_circular_mean(b);
    std::cout << "isotropic null: diamond-square seed " << seed << ", 257 x 257 @ 30 m, "
              << "threshold 20 px: " << links.size() << " links, axial R = "
              << std::setprecision(6) << ax.mean_resultant_length
              << ", directional R = " << dir.mean_resultant_length
              << ", fraction on D8 multiples = "
              << lsdtt::channel_bearing_fraction_on_d8_multiples(b) << "\n";
    check("NU1 isotropic null: all-link axial R < 0.15", ax.mean_resultant_length < 0.15,
          num(ax.mean_resultant_length));
}

// =============================================================================
// --write-fixtures DIR [REFERENCE_CSV]
// =============================================================================
void write_fixtures(const std::string& output_dir, const std::string& reference_csv)
{
    lsdtt::gdal_init();
    std::filesystem::create_directories(output_dir);
    lsdtt::lsdtt_data_loader loader;

    const lsdtt::lsdtt_raster_cpu ns = make_trellis_dem(false);
    loader.write_full(ns.get_data(), ns.get_meta(), output_dir + "/cbs_trellis_dem.tif",
                      lsdtt::raster_output_format::gtiff);
    const lsdtt::lsdtt_raster_cpu ew = make_trellis_dem(true);
    loader.write_full(ew.get_data(), ew.get_meta(), output_dir + "/cbs_trellis_rotated_dem.tif",
                      lsdtt::raster_output_format::gtiff);

    // Region code raster on the trellis grid: 1 = north half, 2 = south half.
    const auto& m = ns.get_meta();
    std::vector<float> codes(static_cast<std::size_t>(m.n_rows) * m.n_cols);
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c)
            codes[static_cast<std::size_t>(r) * m.n_cols + c] = r < m.n_rows / 2 ? 1.0f : 2.0f;
    loader.write_full(codes, m, output_dir + "/cbs_region_halves.tif",
                      lsdtt::raster_output_format::gtiff);

    // A region raster that is off the grid and has no CRS to warp it with.
    lsdtt::lsdtt_raster_meta off = m;
    off.n_rows = 50;
    off.n_cols = 40;
    off.crs_wkt.clear();
    off.epsg_code = 0;
    loader.write_full(std::vector<float>(static_cast<std::size_t>(50) * 40, 1.0f), off,
                      output_dir + "/cbs_region_off_grid.tif",
                      lsdtt::raster_output_format::gtiff);

    // Region / lithology maps 200 km east of the trellis (same CRS, no
    // overlap): a code raster and a one-polygon GeoPackage (field "unit").
    {
        lsdtt::lsdtt_raster_meta far = m;
        far.x_min += 200000.0;
        loader.write_full(std::vector<float>(codes.size(), 1.0f), far,
                          output_dir + "/cbs_region_far.tif", lsdtt::raster_output_format::gtiff);
        const std::string gpkg = output_dir + "/cbs_region_far.gpkg";
        std::filesystem::remove(gpkg);
        GDALDataset* ds = GetGDALDriverManager()->GetDriverByName("GPKG")->Create(
            gpkg.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        OGRSpatialReference srs;
        srs.importFromWkt(m.crs_wkt.c_str());
        srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        OGRLayer* layer = ds->CreateLayer("units", &srs, wkbPolygon, nullptr);
        OGRFieldDefn fd("unit", OFTString);
        (void)layer->CreateField(&fd);
        OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetField("unit", "far");
        const double x0 = far.x_min, x1 = far.x_min + far.n_cols * far.cell_size;
        const double y0 = far.y_min, y1 = far.y_min + far.n_rows * far.cell_size;
        OGRLinearRing ring;
        ring.addPoint(x0, y0); ring.addPoint(x1, y0); ring.addPoint(x1, y1); ring.addPoint(x0, y1);
        ring.addPoint(x0, y0);
        OGRPolygon poly;
        poly.addRing(&ring);
        f->SetGeometry(&poly);
        if (layer->CreateFeature(f) != OGRERR_NONE) std::cerr << "could not write the far polygon\n";
        OGRFeature::DestroyFeature(f);
        GDALClose(ds);
    }

    if (!reference_csv.empty())
    {
        std::filesystem::copy_file(reference_csv,
            std::filesystem::path(output_dir) / std::filesystem::path(reference_csv).filename(),
            std::filesystem::copy_options::overwrite_existing);
    }
    std::cout << "fixture written to " << output_dir << "\n";
}

int finish()
{
    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}

} // namespace

int main(int argc, char* argv[])
{
    try
    {
        const std::string mode = argc >= 2 ? argv[1] : "";
        if (mode == "--write-fixtures" && (argc == 3 || argc == 4))
        {
            write_fixtures(argv[2], argc == 4 ? argv[3] : "");
            return 0;
        }
        if (mode == "--isotropic-null" && argc == 3)
        {
            // Measurement only: the same check on another seed.
            run_isotropic_null_checks(std::stoi(argv[2]));
            return finish();
        }
        if (argc != 2)
        {
            std::cerr << "usage: lsdtt-test-channel-bearing-statistics "
                         "--circular-two-sample | --group-summaries | --comparisons | "
                         "--junction-angle-stats | --trellis-network | --grouping-keys | "
                         "--isotropic-null | --write-fixtures DIR [REFERENCE_CSV]\n";
            return 2;
        }
        if (mode == "--circular-two-sample")       run_circular_two_sample_checks();
        else if (mode == "--group-summaries")      { run_group_summary_checks(); run_angle_format_checks(); }
        else if (mode == "--comparisons")          run_comparison_checks();
        else if (mode == "--junction-angle-stats") run_junction_angle_checks();
        else if (mode == "--trellis-network")      run_trellis_checks();
        else if (mode == "--grouping-keys")        run_grouping_checks();
        else if (mode == "--isotropic-null")       run_isotropic_null_checks(k_null_seed);
        else
        {
            std::cerr << "unknown mode: " << mode << "\n";
            return 2;
        }
        return finish();
    }
    catch (const std::exception& e)
    {
        std::cerr << "EXCEPTION: " << e.what() << "\n";
        return 1;
    }
}
