// test_wetness_index.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_wetness_index_cpu -- the shared Topographic
// Wetness Index.
//
// WHY THIS EXISTS ALONGSIDE rigorous_drainage_twi_refactor_identical
// ------------------------------------------------------------------
// The two guard DIFFERENT claims and neither subsumes the other.
//
//   rigorous_drainage_twi_refactor_identical pins the raster that
//   lsdtt-drainage-analysis writes end to end. That catches a broken hookup --
//   wrong flow graph, wrong window radius, wrong nodata -- but it is a 5x5 DEM
//   with nine valid pixels, and it was MEASURED not to discriminate a float
//   computation from a double one: perturbing the core arithmetic to double
//   left its checksum unchanged. A checksum over nine well-conditioned pixels
//   is simply too coarse for that.
//
//   This driver pins the core function directly, including the precision
//   contract, using values chosen BECAUSE float and double disagree on them.
//
// On the precision contract: the point is not that float is more accurate than
// double here -- it is not, and the difference is at most one ulp of the float32
// the raster stores either way, which is hydrologically nothing. The point is
// that TWI is an existing published output and STI is derived from it, so the
// numbers must not drift silently underneath either. Changing to double is a
// decision to be taken deliberately, with these tests updated, not a side effect
// of tidying an expression.
//
// Usage: lsdtt-test-wetness-index
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_wetness_index_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

constexpr double k_cell      = 5.0;      // so pixel_area == 25 m^2
constexpr double k_nodata    = -9999.0;

lsdtt::lsdtt_raster_meta make_meta(int n_rows, int n_cols)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = k_cell;
    meta.no_data_value = k_nodata;
    meta.set_default_synthetic_georeference();
    return meta;
}

// Build a raster from explicit values; a value of exactly k_nodata is masked off.
lsdtt::lsdtt_raster_cpu make_raster(int n_rows, int n_cols,
                                    const std::vector<float>& values)
{
    lsdtt::lsdtt_array2d_f  data(n_rows, n_cols, static_cast<float>(k_nodata));
    lsdtt::lsdtt_array2d_u8 mask(n_rows, n_cols, std::uint8_t(0));

    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
        {
            const float v = values[static_cast<std::size_t>(r) * n_cols + c];
            data(r, c) = v;
            mask(r, c) = (v == static_cast<float>(k_nodata))
                             ? std::uint8_t(0) : std::uint8_t(1);
        }

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask),
                                   make_meta(n_rows, n_cols));
}

bool bit_equal(float a, float b)
{
    std::uint32_t ua = 0;
    std::uint32_t ub = 0;
    std::memcpy(&ua, &a, sizeof ua);
    std::memcpy(&ub, &b, sizeof ub);
    return ua == ub;
}

// The reference float evaluation, written out in full so the test states the
// contract rather than just calling the thing it is testing.
float twi_float_reference(float area, float slope, double pixel_area)
{
    return std::log((area / static_cast<float>(pixel_area)) / slope);
}

// =============================================================================
// 1. The formula itself
// =============================================================================
void test_formula()
{
    // pixel_area = 25. area = 2500 -> a = 100. slope = 1 -> TWI = ln(100).
    const auto area  = make_raster(1, 1, { 2500.0f });
    const auto slope = make_raster(1, 1, { 1.0f });

    const auto twi = lsdtt::compute_twi_from_area_and_slope(
        area, slope, make_meta(1, 1));

    const float got = twi.get_data_array()(0, 0);
    check("formula: ln(a/tanB) at a=100, tanB=1",
          std::fabs(got - static_cast<float>(std::log(100.0))) < 1.0e-5f,
          "got " + std::to_string(got));
    check("formula: cell is marked valid",
          twi.get_mask_array()(0, 0) == 1);

    // Halving the slope must raise TWI by exactly ln(2).
    const auto slope_half = make_raster(1, 1, { 0.5f });
    const auto twi_half   = lsdtt::compute_twi_from_area_and_slope(
        area, slope_half, make_meta(1, 1));
    const float got_half = twi_half.get_data_array()(0, 0);
    check("formula: halving slope adds ln(2)",
          std::fabs((got_half - got) - static_cast<float>(std::log(2.0))) < 1.0e-5f,
          "delta " + std::to_string(got_half - got));

    // Doubling the area must also raise TWI by exactly ln(2).
    const auto area_double = make_raster(1, 1, { 5000.0f });
    const auto twi_double  = lsdtt::compute_twi_from_area_and_slope(
        area_double, slope, make_meta(1, 1));
    const float got_double = twi_double.get_data_array()(0, 0);
    check("formula: doubling area adds ln(2)",
          std::fabs((got_double - got) - static_cast<float>(std::log(2.0))) < 1.0e-5f,
          "delta " + std::to_string(got_double - got));
}

// =============================================================================
// 2. Validity rules
//
// Zero and negative slope, and zero and negative area, must all leave the cell
// nodata AND masked off rather than clamped to a sentinel. A flat cell is
// exactly where a zero slope arises, so this is the common case, not an edge
// case, and writing a finite number there would be a lie a downstream model
// would happily consume.
// =============================================================================
void test_validity()
{
    const int n = 6;
    const std::vector<float> areas  = {
        2500.0f, 2500.0f, 2500.0f, 0.0f, -5.0f, static_cast<float>(k_nodata)
    };
    const std::vector<float> slopes = {
        1.0f, 0.0f, -0.5f, 1.0f, 1.0f, 1.0f
    };

    const auto twi = lsdtt::compute_twi_from_area_and_slope(
        make_raster(1, n, areas), make_raster(1, n, slopes), make_meta(1, n));

    const auto& d = twi.get_data_array();
    const auto& m = twi.get_mask_array();

    check("validity: positive area and slope is valid", m(0, 0) == 1);
    check("validity: zero slope is masked off",         m(0, 1) == 0);
    check("validity: negative slope is masked off",     m(0, 2) == 0);
    check("validity: zero area is masked off",          m(0, 3) == 0);
    check("validity: negative area is masked off",      m(0, 4) == 0);
    check("validity: nodata area is masked off",        m(0, 5) == 0);

    bool all_nodata = true;
    for (int c = 1; c < n; ++c)
        if (d(0, c) != static_cast<float>(k_nodata)) all_nodata = false;
    check("validity: every invalid cell carries the nodata value", all_nodata);
}

// A cell masked off in the SLOPE raster must be masked off in the output even
// though its area is perfectly good, and vice versa. This is the pair of
// early-continue guards, and getting one of them backwards is the sort of
// transcription slip a refactor can introduce.
void test_input_masks()
{
    const auto area_ok   = make_raster(1, 2, { 2500.0f, 2500.0f });
    const auto slope_gap = make_raster(1, 2, { 1.0f, static_cast<float>(k_nodata) });
    const auto twi_a = lsdtt::compute_twi_from_area_and_slope(
        area_ok, slope_gap, make_meta(1, 2));
    check("masks: masked slope masks the output",
          twi_a.get_mask_array()(0, 0) == 1 && twi_a.get_mask_array()(0, 1) == 0);

    const auto area_gap = make_raster(1, 2, { 2500.0f, static_cast<float>(k_nodata) });
    const auto slope_ok = make_raster(1, 2, { 1.0f, 1.0f });
    const auto twi_b = lsdtt::compute_twi_from_area_and_slope(
        area_gap, slope_ok, make_meta(1, 2));
    check("masks: masked area masks the output",
          twi_b.get_mask_array()(0, 0) == 1 && twi_b.get_mask_array()(0, 1) == 0);
}

// =============================================================================
// 3. The precision contract
//
// These three (area, slope) pairs were found by search BECAUSE the float and
// double evaluations of ln((A/pixel_area)/S) round to DIFFERENT float32 values
// on them. The CLI checksum test does not discriminate here; this does.
// =============================================================================
void test_precision_contract()
{
    struct case_t { float area; float slope; };
    const case_t cases[] = {
        { 4900874.5f,   0.23621974885463715f },
        { 2090614.75f,  1.5143061876296997f  },
        { 947283.3125f, 0.949573814868927f   },
    };

    int n_differ = 0;
    int n_match  = 0;

    for (const case_t& k : cases)
    {
        const auto twi = lsdtt::compute_twi_from_area_and_slope(
            make_raster(1, 1, { k.area }),
            make_raster(1, 1, { k.slope }),
            make_meta(1, 1));

        const float got     = twi.get_data_array()(0, 0);
        const float as_f32  = twi_float_reference(k.area, k.slope, k_cell * k_cell);
        const float as_f64  = static_cast<float>(std::log(
            (static_cast<double>(k.area) / (k_cell * k_cell)) /
            static_cast<double>(k.slope)));

        // Guard the guard: if these ever stop disagreeing (a libm change, a
        // different rounding mode) the case has gone blunt and must be replaced,
        // so say so loudly rather than passing vacuously.
        if (!bit_equal(as_f32, as_f64)) ++n_differ;

        if (bit_equal(got, as_f32)) ++n_match;

        check("precision: float path reproduced bit-for-bit "
              "(area=" + std::to_string(k.area) + ")",
              bit_equal(got, as_f32),
              "got " + std::to_string(got) +
              ", float ref " + std::to_string(as_f32) +
              ", double ref " + std::to_string(as_f64));
    }

    check("precision: the chosen cases still discriminate float from double",
          n_differ == 3,
          std::to_string(n_differ) + " of 3 still differ -- if this drops, "
          "replace the cases, do not relax the test");
    check("precision: every case matched the float reference", n_match == 3);
}

// =============================================================================
// 4. Dimension mismatch
// =============================================================================
void test_dimension_mismatch()
{
    bool threw = false;
    try
    {
        (void)lsdtt::compute_twi_from_area_and_slope(
            make_raster(2, 2, { 1.0f, 1.0f, 1.0f, 1.0f }),
            make_raster(1, 2, { 1.0f, 1.0f }),
            make_meta(2, 2));
    }
    catch (const std::runtime_error&) { threw = true; }
    check("mismatch: disagreeing raster sizes throw", threw);
}

// =============================================================================
// 5. Pixel area really is used
//
// A fixed area raster on a coarser grid must give a SMALLER index, because the
// same contributing area spread over a bigger cell is less specific area. If
// pixel_area were dropped from the expression this check is what notices.
// =============================================================================
void test_pixel_area_scaling()
{
    lsdtt::lsdtt_raster_meta coarse = make_meta(1, 1);
    coarse.cell_size = 10.0;   // pixel_area 100 rather than 25, a factor of 4

    const auto area  = make_raster(1, 1, { 2500.0f });
    const auto slope = make_raster(1, 1, { 1.0f });

    const float fine_v = lsdtt::compute_twi_from_area_and_slope(
        area, slope, make_meta(1, 1)).get_data_array()(0, 0);

    // The inputs must carry the coarse meta too, or the size check is all that
    // is being exercised.
    lsdtt::lsdtt_array2d_f  a_d(1, 1, 2500.0f);
    lsdtt::lsdtt_array2d_u8 a_m(1, 1, std::uint8_t(1));
    lsdtt::lsdtt_array2d_f  s_d(1, 1, 1.0f);
    lsdtt::lsdtt_array2d_u8 s_m(1, 1, std::uint8_t(1));
    const lsdtt::lsdtt_raster_cpu area_c(std::move(a_d), std::move(a_m), coarse);
    const lsdtt::lsdtt_raster_cpu slope_c(std::move(s_d), std::move(s_m), coarse);

    const float coarse_v = lsdtt::compute_twi_from_area_and_slope(
        area_c, slope_c, coarse).get_data_array()(0, 0);

    check("pixel area: coarser cells lower the index by ln(4)",
          std::fabs((fine_v - coarse_v) - static_cast<float>(std::log(4.0))) < 1.0e-5f,
          "fine " + std::to_string(fine_v) + ", coarse " + std::to_string(coarse_v));
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-wetness-index\n"
              << "------------------------\n";

    try
    {
        test_formula();
        test_validity();
        test_input_masks();
        test_precision_contract();
        test_dimension_mismatch();
        test_pixel_area_scaling();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed;
        ++g_n_checks;
    }

    std::cout << "------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;

    return g_n_failed == 0 ? 0 : 1;
}
