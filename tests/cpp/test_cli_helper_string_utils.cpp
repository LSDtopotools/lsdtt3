// test_cli_helper_string_utils.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Phase 2 CLI-consistency helper pins.
//
// Covers:
//   B1          lsdtt::to_lower_ascii            (src/utils/lsdtt_string_utils.hpp)
//   Category A  lsdtt::planar_distance           (src/utils/lsdtt_polygon_geometry.hpp)
//   regression  lsdtt::polyline_length           (same header)
//   regression  lsdtt::cell_centre_x / cell_centre_y
//                                                (src/utils/lsdtt_raster_geometry.hpp)
//
// The two regression pins exist because Stage 2 of the remediation routes CLI
// copies of those onto these library versions; if their behaviour moved, the
// migration would silently change CLI output.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/core/lsdtt_types.hpp"
#include "../../src/utils/lsdtt_polygon_geometry.hpp"
#include "../../src/utils/lsdtt_raster_geometry.hpp"
#include "../../src/utils/lsdtt_string_utils.hpp"

namespace
{

int g_failures = 0;

void check_string(
    const std::string& what,
    const std::string& got,
    const std::string& expected)
{
    if (got != expected)
    {
        std::printf("FAIL %s: expected \"%s\", got \"%s\"\n",
                    what.c_str(), expected.c_str(), got.c_str());
        ++g_failures;
    }
}

void check_close(
    const std::string& what,
    double got,
    double expected,
    double tol = 1e-12)
{
    if (!(std::abs(got - expected) <= tol))
    {
        std::printf("FAIL %s: expected %.17g, got %.17g (tol %.3g)\n",
                    what.c_str(), expected, got, tol);
        ++g_failures;
    }
}

void check_true(const std::string& what, bool ok)
{
    if (!ok)
    {
        std::printf("FAIL %s: expected true\n", what.c_str());
        ++g_failures;
    }
}

// -----------------------------------------------------------------------------
// B1: to_lower_ascii
// -----------------------------------------------------------------------------
void test_to_lower_ascii()
{
    check_string("to_lower_ascii(\"ABC\")",
                 lsdtt::to_lower_ascii("ABC"), "abc");

    check_string("to_lower_ascii mixed case",
                 lsdtt::to_lower_ascii("Same_As_Input"), "same_as_input");

    check_string("to_lower_ascii already lowercase",
                 lsdtt::to_lower_ascii("d_infinity"), "d_infinity");

    check_string("to_lower_ascii empty string",
                 lsdtt::to_lower_ascii(""), "");

    // Digits, punctuation and whitespace must pass through untouched.
    check_string("to_lower_ascii digits/punctuation",
                 lsdtt::to_lower_ascii("D8 / M-FD_3.14 [TIF]"),
                 "d8 / m-fd_3.14 [tif]");

    // A realistic parameter-string use: a file extension.
    check_string("to_lower_ascii extension",
                 lsdtt::to_lower_ascii(".TIF"), ".tif");

    // Bytes >= 0x80. The point of the unsigned char cast: passing a negative
    // char to std::tolower is UB. Assert we do not crash and that the bytes
    // come back unchanged.
    {
        std::string high;
        high.push_back(static_cast<char>(0x80));
        high.push_back(static_cast<char>(0xC3));  // UTF-8 lead byte
        high.push_back(static_cast<char>(0xA9));  // UTF-8 continuation ("e-acute")
        high.push_back(static_cast<char>(0xFF));
        high.push_back('A');

        const std::string got = lsdtt::to_lower_ascii(high);

        check_true("to_lower_ascii high bytes: length preserved",
                   got.size() == high.size());
        for (std::size_t i = 0; i + 1 < high.size(); ++i)
        {
            if (got[i] != high[i])
            {
                std::printf("FAIL to_lower_ascii high byte %zu: expected 0x%02X, "
                            "got 0x%02X\n",
                            i,
                            static_cast<unsigned>(
                                static_cast<unsigned char>(high[i])),
                            static_cast<unsigned>(
                                static_cast<unsigned char>(got[i])));
                ++g_failures;
            }
        }
        check_true("to_lower_ascii high bytes: trailing ASCII still folded",
                   got.back() == 'a');
    }
}

// -----------------------------------------------------------------------------
// Category A: planar_distance
// -----------------------------------------------------------------------------
void test_planar_distance()
{
    check_close("planar_distance 3-4-5",
                lsdtt::planar_distance(lsdtt::lsdtt_xy{0.0, 0.0},
                                       lsdtt::lsdtt_xy{3.0, 4.0}),
                5.0);

    check_close("planar_distance zero",
                lsdtt::planar_distance(lsdtt::lsdtt_xy{7.5, -2.25},
                                       lsdtt::lsdtt_xy{7.5, -2.25}),
                0.0);

    // Symmetric in its arguments (the burn-streams copy computes a - b rather
    // than b - a; the two agree because the result is squared).
    check_close("planar_distance symmetry",
                lsdtt::planar_distance(lsdtt::lsdtt_xy{3.0, 4.0},
                                       lsdtt::lsdtt_xy{0.0, 0.0}),
                5.0);

    check_close("planar_distance negative coordinates",
                lsdtt::planar_distance(lsdtt::lsdtt_xy{-3.0, -4.0},
                                       lsdtt::lsdtt_xy{-6.0, -8.0}),
                5.0);

    // Realistic projected (UTM metre) magnitudes: no precision loss.
    check_close("planar_distance UTM magnitudes",
                lsdtt::planar_distance(lsdtt::lsdtt_xy{500000.0, 5000000.0},
                                       lsdtt::lsdtt_xy{500030.0, 5000040.0}),
                50.0,
                1e-9);

    // Large-magnitude behaviour. The implementation is deliberately plain
    // sqrt(dx*dx + dy*dy), matching the three copies it replaces, so the
    // squares overflow to infinity here. std::hypot would return a finite
    // value; this test pins the CURRENT behaviour so that a later switch to
    // hypot is a deliberate, visible change rather than a silent one. No
    // projected coordinate pair comes within many orders of magnitude of this.
    {
        const double big = 1.0e200;
        const double d   = lsdtt::planar_distance(lsdtt::lsdtt_xy{0.0, 0.0},
                                                  lsdtt::lsdtt_xy{big, big});
        check_true("planar_distance 1e200: plain-sqrt overflow to inf "
                   "(NOT hypot semantics)",
                   std::isinf(d));
    }
}

// -----------------------------------------------------------------------------
// Regression pins for the helpers Stage 2 will route CLIs onto
// -----------------------------------------------------------------------------
void test_polyline_length()
{
    check_close("polyline_length empty",
                lsdtt::polyline_length(std::vector<lsdtt::lsdtt_xy>{}), 0.0);

    check_close("polyline_length single vertex",
                lsdtt::polyline_length(
                    std::vector<lsdtt::lsdtt_xy>{{1.0, 2.0}}), 0.0);

    const std::vector<lsdtt::lsdtt_xy> line{
        {0.0, 0.0}, {3.0, 4.0}, {3.0, 14.0}};
    check_close("polyline_length 3 vertices", lsdtt::polyline_length(line), 15.0);

    // OPEN polyline: no implied closing edge (that is polygon_perimeter's job).
    const std::vector<lsdtt::lsdtt_xy> tri{{0.0, 0.0}, {3.0, 0.0}, {3.0, 4.0}};
    check_close("polyline_length is open (no closing edge)",
                lsdtt::polyline_length(tri), 7.0);
    check_close("polygon_perimeter closes the ring",
                lsdtt::polygon_perimeter(tri), 12.0);
}

void test_cell_centres()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows    = 4;
    meta.n_cols    = 5;
    meta.cell_size = 10.0;
    meta.x_min     = 1000.0;
    meta.y_min     = 2000.0;

    check_close("cell_centre_x col 0", lsdtt::cell_centre_x(meta, 0), 1005.0);
    check_close("cell_centre_x col 4", lsdtt::cell_centre_x(meta, 4), 1045.0);

    // Row 0 is the TOP row: y decreases with increasing row index.
    check_close("cell_centre_y row 0 (top)",
                lsdtt::cell_centre_y(meta, 0), 2035.0);
    check_close("cell_centre_y row 3 (bottom)",
                lsdtt::cell_centre_y(meta, 3), 2005.0);

    // The US-spelling aliases must agree.
    check_close("cell_center_x alias",
                lsdtt::cell_center_x(meta, 2), lsdtt::cell_centre_x(meta, 2));
    check_close("cell_center_y alias",
                lsdtt::cell_center_y(meta, 2), lsdtt::cell_centre_y(meta, 2));
}

} // namespace

int main()
{
    std::printf("test_cli_helper_string_utils: B1 to_lower_ascii + "
                "Category A planar_distance\n");

    test_to_lower_ascii();
    test_planar_distance();
    test_polyline_length();
    test_cell_centres();

    if (g_failures != 0)
    {
        std::printf("test_cli_helper_string_utils: %d FAILURE(S)\n", g_failures);
        return 1;
    }

    std::printf("test_cli_helper_string_utils: all checks passed\n");
    return 0;
}
