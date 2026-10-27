// test_surface_method_comparison.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// compare_surface_classes / surface_class_agreement_bits
// (lsdtt_surface_method_comparison_cpu), fluvial-surface methods plan 7b/7d.
// Hand-checkable class maps: a method compared with itself, two methods that
// disagree on known cells, invalid cells, and the bit layout.
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_surface_method_comparison_cpu.hpp"

namespace
{
int g_failures = 0;
void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "PASS  " : "FAIL  ") << what << "\n";
    if (!ok) ++g_failures;
}
} // namespace

int main()
{
    using lsdtt::method_class_map;
    //                 0  1  2  3  4  5  6  7  8  9
    const method_class_map a{"a", {1, 1, 1, 2, 2, 0, 0, 0, 1, 2}};
    const method_class_map b{"b", {1, 1, 0, 2, 1, 0, 0, 2, 1, 2}};
    std::vector<std::uint8_t> valid(10, 1);

    {
        const auto r = lsdtt::compare_surface_classes({a, a}, valid);
        check(r.size() == 1 && r[0].floodplain_iou == 1.0 &&
                  r[0].terrace_iou == 1.0 && r[0].classified_agreement == 1.0,
              "a method against itself: IoU 1, agreement 1");
        std::int64_t off = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (i != j) off += r[0].confusion[i][j];
        check(off == 0, "a method against itself: no off-diagonal cells");
    }
    {
        const auto r = lsdtt::compare_surface_classes({a, b}, valid);
        const auto& c = r[0].confusion;
        std::int64_t sum = 0;
        for (const auto& row : c) for (const auto v : row) sum += v;
        check(sum == 10, "confusion sums to the valid cell count");
        // floodplain: a {0,1,2,8}, b {0,1,4,8} -> inter 3, union 5
        check(std::fabs(r[0].floodplain_iou - 3.0 / 5.0) < 1e-12,
              "floodplain IoU 3/5");
        // terrace: a {3,4,9}, b {3,7,9} -> inter 2, union 4
        check(std::fabs(r[0].terrace_iou - 2.0 / 4.0) < 1e-12, "terrace IoU 2/4");
        // classified by either: cells 0,1,2,3,4,7,8,9 = 8; same class: 0,1,3,8,9 = 5
        check(std::fabs(r[0].classified_agreement - 5.0 / 8.0) < 1e-12,
              "classified agreement 5/8");
        check(c[1][0] == 1 && c[2][1] == 1 && c[0][2] == 1,
              "the three disagreements land in [1][0], [2][1], [0][2]");
    }
    {
        valid[0] = 0;   // cell 0 (floodplain in both) leaves the comparison
        const auto r = lsdtt::compare_surface_classes({a, b}, valid);
        check(std::fabs(r[0].floodplain_iou - 2.0 / 4.0) < 1e-12,
              "an invalid cell takes no part (floodplain IoU 2/4)");
        valid[0] = 1;
    }
    {
        const auto bits = lsdtt::surface_class_agreement_bits({a, b}, valid, 1);
        check(bits[0] == 3 && bits[2] == 1 && bits[4] == 2 && bits[5] == 0,
              "agreement bits: both 3, only a 1, only b 2, neither 0");
    }
    {
        const auto r = lsdtt::compare_surface_classes(
            {method_class_map{"x", {0, 0}}, method_class_map{"y", {0, 0}}},
            std::vector<std::uint8_t>(2, 1));
        check(std::isnan(r[0].floodplain_iou) && std::isnan(r[0].classified_agreement),
              "no class in either method: IoU and agreement are NaN");
    }
    {
        bool threw = false;
        try { (void)lsdtt::compare_surface_classes({a, method_class_map{"c", {1}}}, valid); }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "a class map off the grid is refused");
        threw = false;
        try { (void)lsdtt::compare_surface_classes({method_class_map{"d", std::vector<std::uint8_t>(10, 3)}}, valid); }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "a class above 2 is refused");
    }
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
