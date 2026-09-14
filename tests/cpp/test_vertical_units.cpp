// test_vertical_units.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Vertical units, and the plausibility warning that makes them discoverable.
//
// The defect this pins, found 2026-08-25 on ../data/IndianCreek_1m: the DEM was
// NAVD88 FEET on a 1 m UTM grid. No raster format records a vertical unit, so
// nothing detected it, and every slope, curvature, relief, E* and R* in the run
// was too large by 3.28 WITHOUT LOOKING BROKEN -- the landscape merely appeared
// steeper. It was diagnosed instead as "S_c = 0.8 is too low for this
// catchment", and a whole analysis was written around that wrong conclusion.
//
// So the tests are as much about the WARNING as the conversion:
//
//   V1  parsing accepts the documented spellings
//   V2  A TYPO IS FATAL. Falling back to metres would turn a misspelling into a
//       silent 3.28x error, which is the exact failure mode above
//   V3  conversion scales valid pixels and LEAVES NODATA ALONE
//   V4  metres is bit-for-bit a no-op
//   V5  us survey feet differ from international feet, but immeasurably
//   V6  the plausibility check fires on an impossibly steep surface and stays
//       quiet on a gentle one -- the MEDIAN is what is tested, so a few cliffs
//       in gentle terrain must not trip it
//   V7  the round trip: a feet-DEM trips the check, and converting clears it
//
// Usage: lsdtt-test-vertical-units
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{

int g_checks   = 0;
int g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

void check_close(double got, double expected, double tol, const std::string& what)
{
    ++g_checks;
    if (!(std::fabs(got - expected) <= tol))
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n        got " << got
                  << ", expected " << expected << "\n";
    }
    else std::cout << "  ok  : " << what << "  (" << got << ")\n";
}

// A raster of `n` x `n` cells on a 1 m grid, elevation = slope * x, so the
// surface has a known, uniform gradient.
lsdtt::lsdtt_raster_cpu ramp(int n, double gradient, float nodata_value = -9999.0f)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n; meta.n_cols = n; meta.cell_size = 1.0;
    meta.no_data_value = nodata_value;

    std::vector<float>   d(static_cast<std::size_t>(n) * n);
    std::vector<uint8_t> m(static_cast<std::size_t>(n) * n, 1);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            d[static_cast<std::size_t>(r) * n + c] =
                static_cast<float>(gradient * c);

    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

}  // namespace


int main()
{
    using namespace lsdtt;

    std::cout << "=====================================================\n"
              << " Vertical units\n"
              << "=====================================================\n";

    // -------------------------------------------------------------------------
    std::cout << "\nV1: the documented spellings parse\n";
    check(parse_vertical_unit("metres") == vertical_unit::metres, "metres");
    check(parse_vertical_unit("meters") == vertical_unit::metres, "meters");
    check(parse_vertical_unit("M")      == vertical_unit::metres, "M (case)");
    check(parse_vertical_unit("feet")   == vertical_unit::feet,   "feet");
    check(parse_vertical_unit(" ft ")   == vertical_unit::feet,   "ft (spaces)");
    check(parse_vertical_unit("us_survey_feet") == vertical_unit::us_survey_feet,
          "us_survey_feet");

    // -------------------------------------------------------------------------
    std::cout << "\nV2: a typo is FATAL, never a silent fallback to metres\n";
    for (const char* bad : {"fee", "metre s", "", "yards", "0.3048"})
    {
        bool threw = false;
        try { (void)parse_vertical_unit(bad); }
        catch (const std::exception&) { threw = true; }
        check(threw, std::string("rejects '") + bad + "'");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nV3: conversion scales valid pixels, leaves nodata alone\n";
    {
        lsdtt_raster_cpu r = ramp(8, 1.0);
        // Punch a nodata hole and mark it invalid, as a loader would.
        r.get_data()[10]      = -9999.0f;
        r.get_mask_array().get_vector()[10] = 0;

        const lsdtt_raster_cpu out =
            convert_vertical_units_to_metres(r, vertical_unit::feet);

        check_close(out.get_data()[3], 3.0 * 0.3048, 1e-6, "valid pixel x 0.3048");
        check_close(out.get_data()[10], -9999.0, 1e-9,
                    "NODATA IS NOT SCALED (would move off the sentinel)");
        check(out.get_mask()[10] == 0, "mask preserved");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nV4: metres is a bit-for-bit no-op\n";
    {
        const lsdtt_raster_cpu r = ramp(8, 0.37);
        const lsdtt_raster_cpu out =
            convert_vertical_units_to_metres(r, vertical_unit::metres);
        bool identical = true;
        for (std::size_t i = 0; i < r.get_data().size(); ++i)
            if (r.get_data()[i] != out.get_data()[i]) identical = false;
        check(identical, "every pixel bit-identical");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nV5: us survey feet differ from international feet\n";
    {
        const lsdtt_raster_cpu r = ramp(8, 1000.0);
        const double intl = convert_vertical_units_to_metres(
                                r, vertical_unit::feet).get_data()[3];
        const double usft = convert_vertical_units_to_metres(
                                r, vertical_unit::us_survey_feet).get_data()[3];
        check(intl != usft, "the two feet are not the same factor");
        // 2 ppm: real, but never measurable on a catchment. Documented so nobody
        // reaches for it to explain a discrepancy.
        check(std::fabs(intl - usft) / intl < 5e-6, "they differ by ~2 ppm");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nV6: the plausibility check tests the MEDIAN\n";
    {
        // Gentle: 15 degrees, comfortably plausible.
        const auto gentle = check_vertical_plausibility(ramp(40, 0.268));
        check(!gentle.implausible, "15 deg surface is not flagged");
        check_close(gentle.median_slope_deg, 15.0, 0.5, "median slope recovered");

        // Impossibly steep: 42 degrees everywhere, as Indian Creek appeared.
        const auto steep = check_vertical_plausibility(ramp(40, 0.900));
        check(steep.implausible, "42 deg surface IS flagged");
        check_close(steep.median_slope_deg_if_feet,
                    std::atan(0.900 * 0.3048) * 180.0 / M_PI, 0.5,
                    "reports what the median would be if feet");

        // A few cliffs in gentle terrain must NOT trip it: that is the whole
        // reason the median is tested rather than a maximum.
        lsdtt_raster_cpu spotty = ramp(40, 0.268);
        for (int r = 5; r < 9; ++r)
            for (int c = 5; c < 9; ++c)
                spotty.get_data()[static_cast<std::size_t>(r) * 40 + c] += 60.0f;
        check(!check_vertical_plausibility(spotty).implausible,
              "isolated cliffs in gentle terrain do not flag");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nV7: round trip -- feet trips the check, converting clears it\n";
    {
        // A 15 degree landscape whose elevations are recorded in feet reads as
        // 0.268 * 3.28084 = 0.879, i.e. 41 degrees. Exactly Indian Creek.
        const lsdtt_raster_cpu as_feet = ramp(40, 0.268 / 0.3048);

        std::ostringstream warned;
        const auto before = warn_if_vertically_implausible(as_feet, warned);
        check(before.implausible, "feet-on-a-metre-grid is flagged");
        check(warned.str().find("vertical_units: feet") != std::string::npos,
              "the warning names the parameter that fixes it");
        check(warned.str().find("WILL NOT CONVERT THIS FOR YOU") != std::string::npos,
              "the warning is explicit that nothing is auto-converted");

        std::ostringstream quiet;
        const auto after = warn_if_vertically_implausible(
            convert_vertical_units_to_metres(as_feet, vertical_unit::feet), quiet);
        check(!after.implausible, "converted DEM is not flagged");
        check(quiet.str().empty(), "and prints nothing at all");
        check_close(after.median_slope_deg, 15.0, 0.5, "back to 15 degrees");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
