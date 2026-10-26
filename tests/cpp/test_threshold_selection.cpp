// test_threshold_selection.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Checks lsdtt_threshold_selection, the shared threshold picker
// (fixed | qq_gaussian | cdf_divergence).
//
// The qq_gaussian checks use samples whose answer is known in advance rather
// than random mixtures, which have no exact answer to compare with:
//   * a Gaussian sample built from exact normal quantiles, with everything below
//     the 20th percentile compressed towards it - the threshold must land just
//     below the 20th percentile, where the compressed tail rejoins the line;
//   * the same sample uncompressed - every quantile agrees, so the threshold is
//     the first one (a rule that "always finds something in the middle" fails);
//   * QQ tables written directly, to pin the persistence rule exactly: the
//     candidate plus persistence_points quantiles after it, all inside the
//     arrays.
//
// cdf_divergence is not re-derived here: its numbers are pinned end to end by
// rigorous_flat_surface_finder_cdf_divergence and
// rigorous_mountain_front_cdf_divergence, recorded before the code moved.

#include "lsdtt_histogram.hpp"
#include "lsdtt_statistics.hpp"
#include "lsdtt_threshold_selection.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int n_failed = 0;
int n_checks = 0;

void check(bool condition, const std::string& label)
{
    ++n_checks;
    if (condition)
    {
        std::cout << "PASS: " << label << "\n";
    }
    else
    {
        std::cout << "FAIL: " << label << "\n";
        ++n_failed;
    }
}

// n values at the exact standard-normal quantiles of the plotting positions
// compute_qq uses, scaled to mean mu and standard deviation sigma.
std::vector<double> exact_gaussian(std::size_t n, double mu, double sigma)
{
    std::vector<double> v(n);
    for (std::size_t k = 0; k < n; ++k)
    {
        const double p = (static_cast<double>(k) + 0.5) / static_cast<double>(n);
        v[k] = mu + sigma * lsdtt::normal_quantile(p);
    }
    return v;
}

lsdtt::threshold_selection_request qq_request()
{
    lsdtt::threshold_selection_request r;
    r.method = lsdtt::threshold_selection_choice::qq_gaussian;
    return r;   // paper defaults: 25/75, 1 %, 10000 points, 200 persistence
}

void test_compressed_tail()
{
    std::cout << "\n--- qq_gaussian: a compressed low tail ---\n";
    const std::size_t n = 10000;
    std::vector<double> v = exact_gaussian(n, 10.0, 2.0);
    // Everything below the 20th percentile is pulled to a fifth of its distance
    // from it: a low population that is not on the Gaussian line. The reference
    // line (through Q25 and Q75) is untouched. Worked out by hand, the
    // compressed points rejoin the line - within 1 % of the range - at about
    // the 18.3rd percentile.
    const double q20 = v[2000];
    for (std::size_t k = 0; k < 2000; ++k)
        v[k] = q20 - 0.2 * (q20 - v[k]);

    const lsdtt::threshold_selection_result r =
        lsdtt::select_threshold(v, qq_request());
    std::cout << "  qq_index " << r.qq_index << ", threshold " << r.threshold
              << " (20th percentile value " << q20 << ")\n";
    check(r.found, "QQ1: a threshold is found");
    check(r.qq_index >= 1800 && r.qq_index <= 2000,
          "QQ2: it sits within 2 percentiles below the break (index 1800-2000 "
          "of 10000)");
    check(r.threshold <= q20 && r.threshold > v[1700],
          "QQ3: the threshold value lies just below the 20th percentile");
}

void test_pure_gaussian()
{
    std::cout << "\n--- qq_gaussian: a pure Gaussian ---\n";
    const std::vector<double> v = exact_gaussian(10000, 10.0, 2.0);
    const lsdtt::threshold_selection_result r =
        lsdtt::select_threshold(v, qq_request());
    // Not exactly index 0: the lowest quantile is read by type-7 interpolation
    // between the two lowest sample values, which sit 0.27 sigma apart in the
    // extreme tail, and misses the line by 1.8 % of the range. From the second
    // quantile on every point agrees. What matters is that the threshold is the
    // bottom of the distribution, not somewhere inside it.
    std::cout << "  qq_index " << r.qq_index << "\n";
    check(r.found && r.qq_index <= 10,
          "QQ4: every quantile agrees, so the threshold is at the very bottom "
          "(within the first 0.1 %)");
}

void test_persistence_rule()
{
    std::cout << "\n--- the persistence rule on written QQ tables ---\n";
    // reference 0 everywhere; a 100 at index 0 sets the range to 100, so a
    // value agrees when it is below 1 (1 % of 100). Disagreements are 50.
    const std::size_t P = 200;
    auto table = [](std::size_t n, const std::vector<std::size_t>& disagree)
    {
        std::vector<double> dq(n, 0.0);
        dq[0] = 100.0;
        for (std::size_t i : disagree) dq[i] = 50.0;
        return dq;
    };
    const std::vector<double> ref_1000(1000, 0.0);

    // Disagreements every 201 positions: every agreeing run is exactly 200
    // long - one short of a candidate plus 200 more.
    {
        std::vector<std::size_t> d;
        for (std::size_t i = 201; i < 1000; i += 201) d.push_back(i);
        d.push_back(999);   // the tail run is short too
        const long i = lsdtt::qq_first_persistent_agreement(
            table(1000, d), ref_1000, 0.01, P);
        check(i == -1, "QQ5: runs of 200 agreeing quantiles never qualify");
    }
    // Widen one gap by one: a single run of 201, starting at 403.
    {
        const std::vector<std::size_t> d = {201, 402, 604, 805, 999};
        const long i = lsdtt::qq_first_persistent_agreement(
            table(1000, d), ref_1000, 0.01, P);
        check(i == 403, "QQ6: the first run of 201 qualifies, at its start");
    }
    // At the top of the array: agreement only over the last 201 quantiles
    // qualifies; over the last 200 it does not (the rule never looks past the
    // end - the legacy code read beyond it here).
    {
        std::vector<std::size_t> d;
        for (std::size_t i = 1; i < 1000 - 201; ++i) d.push_back(i);
        check(lsdtt::qq_first_persistent_agreement(
                  table(1000, d), ref_1000, 0.01, P) == 1000 - 201,
              "QQ7: a run of 201 ending at the last quantile qualifies");
        d.push_back(1000 - 201);
        check(lsdtt::qq_first_persistent_agreement(
                  table(1000, d), ref_1000, 0.01, P) == -1,
              "QQ8: a run of 200 ending at the last quantile does not");
    }
    // A constant sample has no range: no threshold, reported as not found.
    {
        const std::vector<double> flat(5000, 3.0);
        const lsdtt::threshold_selection_result r =
            lsdtt::select_threshold(flat, qq_request());
        check(!r.found, "QQ9: a constant sample reports found = false");
    }
}

void test_compute_qq_percentiles()
{
    std::cout << "\n--- compute_qq percentile arguments ---\n";
    std::vector<double> v = exact_gaussian(3000, 0.0, 1.0);
    for (std::size_t k = 0; k < 400; ++k) v[k] *= 3.0;   // a heavy low tail
    const lsdtt::lsdtt_qq_data a = lsdtt::compute_qq(v, 500);
    const lsdtt::lsdtt_qq_data b = lsdtt::compute_qq(v, 500, 25.0, 75.0);
    check(a.normal_reference == b.normal_reference &&
              a.data_quantile == b.data_quantile,
          "QQ10: the default is exactly the 25 / 75 line");
    const lsdtt::lsdtt_qq_data c = lsdtt::compute_qq(v, 500, 10.0, 90.0);
    check(c.normal_reference != a.normal_reference,
          "QQ11: other percentiles move the reference line");
    bool threw = false;
    try { (void)lsdtt::compute_qq(v, 500, 75.0, 25.0); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "QQ12: lower >= upper throws");
}

void test_fixed_cdf_and_names()
{
    std::cout << "\n--- fixed, cdf_divergence guards, names ---\n";
    lsdtt::threshold_selection_request r;
    r.method      = lsdtt::threshold_selection_choice::fixed;
    r.fixed_value = 2.5;
    const lsdtt::threshold_selection_result f = lsdtt::select_threshold({}, r);
    check(f.found && f.threshold == 2.5, "TS1: fixed returns the value");

    r.method = lsdtt::threshold_selection_choice::cdf_divergence;
    bool threw = false;
    try { (void)lsdtt::select_threshold({1.0, 2.0, 3.0}, r); }
    catch (const std::runtime_error&) { threw = true; }
    check(threw, "TS2: cdf_divergence with fewer than 10 values throws");

    bool round_trip = true;
    for (const char* name : {"fixed", "qq_gaussian", "cdf_divergence"})
        round_trip = round_trip &&
            lsdtt::threshold_selection_choice_name(
                lsdtt::parse_threshold_selection_choice(name)) == name;
    check(round_trip, "TS3: choice names round-trip");
    threw = false;
    try { (void)lsdtt::parse_threshold_selection_choice("qq"); }
    catch (const std::invalid_argument& e)
    {
        threw = std::string(e.what()).find("qq_gaussian") != std::string::npos;
    }
    check(threw, "TS4: an unknown name throws, listing the valid ones");
}

} // namespace

int main()
{
    test_compressed_tail();
    test_pure_gaussian();
    test_persistence_rule();
    test_compute_qq_percentiles();
    test_fixed_cdf_and_names();

    std::cout << "\n" << (n_checks - n_failed) << " / " << n_checks
              << " checks passed\n";
    return n_failed == 0 ? 0 : 1;
}
