// test_cosmo_ratio.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// The isotope-ratio output policy (make_isotope_ratio_raster).
//
// A ratio raster looks like the most trivial thing in the codebase -- one
// division per pixel -- which is exactly why it is worth pinning. The failure
// modes are all silent: an inf written into a raster that a basin average then
// consumes, a nodata numerator divided and emitted as a number, a mask that
// disagrees with its own data.
//
// The policy is promoted from a lambda inside lsdtt-cosmo-prediction so that
// the landslide coupling and the basin-accumulated path share it rather than
// each re-deciding what a zero denominator means.
//
//   R1  a plain ratio is the plain quotient
//   R2  zero denominator -> nodata and COUNTED, never inf
//   R3  one input masked -> nodata and counted. The original lambda tested only
//       the DENOMINATOR's mask, so it would divide a nodata numerator and emit
//       the result; harmless for MultiNuclideRasters, wrong in general
//   R4  both inputs masked -> outside the domain, counted as nothing
//   R5  geometry mismatch throws rather than reading out of bounds
//   R6  the output mask is exactly the set of pixels with a finite ratio
//   R7  a non-finite quotient is caught and counted separately -- unreachable
//       given R2/R3, so if it ever fires something upstream produced inf/NaN
//   R8  THE PROMOTION IS FAITHFUL: on MultiNuclideRasters-shaped input (all
//       three nuclides share a mask) the helper reproduces the original lambda
//       bit for bit. This is what makes the refactor safe at unit level rather
//       than only through the CLI oracle.
//
// Usage: lsdtt-test-cosmo-ratio
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_cosmo_raster_cpu.hpp"
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

lsdtt::lsdtt_raster_meta make_meta(int n_rows, int n_cols)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = n_rows;
    m.n_cols = n_cols;
    m.cell_size = 10.0;
    m.no_data_value = -9999.0;
    m.set_default_synthetic_georeference();
    return m;
}

// A raster with per-pixel values and a per-pixel mask, both given as flat rows.
lsdtt::lsdtt_raster_cpu make_raster(const lsdtt::lsdtt_raster_meta& m,
                                    const std::vector<float>& vals,
                                    const std::vector<std::uint8_t>& mask)
{
    lsdtt::lsdtt_array2d_f  d(m.n_rows, m.n_cols, static_cast<float>(m.no_data_value));
    lsdtt::lsdtt_array2d_u8 k(m.n_rows, m.n_cols, 0);
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * m.n_cols + c;
            d(r, c) = vals[i];
            k(r, c) = mask[i];
        }
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(k), m);
}

}  // namespace


int main()
{
    std::cout << "=====================================================\n"
              << " Isotope ratio output policy\n"
              << "=====================================================\n";

    const lsdtt::lsdtt_raster_meta meta = make_meta(2, 3);
    //                 idx:      0     1     2      3      4     5
    const std::vector<float> num{ 6.0f, 0.0f, 5.0f,  9.0f,  1.0f, 2.0f};
    const std::vector<float> den{ 3.0f, 4.0f, 0.0f,  3.0f,  1.0f, 8.0f};

    // -------------------------------------------------------------------------
    std::cout << "\nR1/R2/R6: plain quotient, zero denominator, output mask\n";
    {
        const std::vector<std::uint8_t> all{1, 1, 1, 1, 1, 1};
        const auto rr = lsdtt::make_isotope_ratio_raster(
            make_raster(meta, num, all), make_raster(meta, den, all));

        const auto& d = rr.ratio.get_data_array();
        const auto& m = rr.ratio.get_mask_array();

        check(std::fabs(d(0, 0) - 2.0f) < 1e-6f, "6/3 = 2");
        check(std::fabs(d(0, 1) - 0.0f) < 1e-6f,
              "a ZERO NUMERATOR is a legitimate ratio of 0, not nodata");
        check(m(0, 2) == 0, "zero denominator -> nodata");
        check(rr.n_zero_denominator == 1, "the zero denominator is COUNTED");
        check(std::isfinite(static_cast<double>(d(0, 2))) == false ||
              d(0, 2) == static_cast<float>(meta.no_data_value),
              "the zero-denominator pixel carries the nodata value, not inf");
        check(rr.n_valid == 5, "five valid pixels");

        int n_masked = 0;
        for (int r = 0; r < meta.n_rows; ++r)
            for (int c = 0; c < meta.n_cols; ++c)
                if (m(r, c)) ++n_masked;
        check(n_masked == rr.n_valid,
              "the output MASK agrees with the reported valid count");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nR3/R4: partial and total masking\n";
    {
        // pixel 0: numerator masked out, denominator valid  -> one_input_masked
        // pixel 5: both masked out                          -> silent, outside
        const std::vector<std::uint8_t> nmask{0, 1, 1, 1, 1, 0};
        const std::vector<std::uint8_t> dmask{1, 1, 1, 1, 1, 0};
        const auto rr = lsdtt::make_isotope_ratio_raster(
            make_raster(meta, num, nmask), make_raster(meta, den, dmask));

        check(rr.n_one_input_masked == 1,
              "numerator masked with a valid denominator is counted, not divided");
        check(rr.ratio.get_mask_array()(0, 0) == 0,
              "and it is nodata in the output");
        check(rr.n_valid == 3,
              "both-masked pixels are outside the domain, not failures "
              "(3 valid: idx 1,3,4; idx 2 has a zero denominator)");
        check(rr.n_zero_denominator == 1, "the zero denominator still counts");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nR7: a non-finite quotient is caught, not emitted\n";
    {
        const std::vector<std::uint8_t> all{1, 1, 1, 1, 1, 1};
        std::vector<float> bad = num;
        bad[3] = std::numeric_limits<float>::infinity();
        const auto rr = lsdtt::make_isotope_ratio_raster(
            make_raster(meta, bad, all), make_raster(meta, den, all));

        check(rr.n_non_finite == 1, "an inf numerator is counted as non-finite");
        check(rr.ratio.get_mask_array()(1, 0) == 0,
              "and is nodata rather than an inf written into the raster");
        check(rr.summary("Al26/Be10").find("NON-FINITE") != std::string::npos,
              "summary() says so loudly -- it means an inf/NaN upstream");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nR5: geometry mismatch throws\n";
    {
        const lsdtt::lsdtt_raster_meta other = make_meta(3, 2);
        const std::vector<std::uint8_t> all{1, 1, 1, 1, 1, 1};
        bool threw = false;
        try
        {
            (void)lsdtt::make_isotope_ratio_raster(
                make_raster(other, num, all), make_raster(meta, den, all));
        }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "mismatched geometry throws invalid_argument");
    }

    // -------------------------------------------------------------------------
    // R8 -- the promotion is faithful.
    //
    // Reproduces the ORIGINAL lambda from lsdtt-cosmo-prediction verbatim and
    // requires bit-for-bit agreement on the input shape it was written for:
    // MultiNuclideRasters, where all three nuclides share one mask. Anything
    // else and the promotion would have quietly changed a shipped output.
    // -------------------------------------------------------------------------
    std::cout << "\nR8: bit-for-bit agreement with the original lambda\n";
    {
        const std::vector<std::uint8_t> shared{1, 1, 1, 0, 1, 1};
        const lsdtt::lsdtt_raster_cpu N = make_raster(meta, num, shared);
        const lsdtt::lsdtt_raster_cpu D = make_raster(meta, den, shared);

        // --- the original, as it stood on 2026-08-05 ---
        const float ndv = static_cast<float>(meta.no_data_value);
        lsdtt::lsdtt_array2d_f  r_d(meta.n_rows, meta.n_cols, ndv);
        lsdtt::lsdtt_array2d_u8 r_m(meta.n_rows, meta.n_cols, 0);
        {
            const auto& be_d = D.get_data_array();
            const auto& be_m = D.get_mask_array();
            const auto& n_d  = N.get_data_array();
            for (int r = 0; r < meta.n_rows; ++r)
                for (int c = 0; c < meta.n_cols; ++c)
                {
                    if (!be_m(r, c)) continue;
                    if (be_d(r, c) == 0.0f) continue;
                    r_d(r, c) = n_d(r, c) / be_d(r, c);
                    r_m(r, c) = 1;
                }
        }

        const auto rr = lsdtt::make_isotope_ratio_raster(N, D);
        const auto& g_d = rr.ratio.get_data_array();
        const auto& g_m = rr.ratio.get_mask_array();

        bool same = true;
        for (int r = 0; r < meta.n_rows && same; ++r)
            for (int c = 0; c < meta.n_cols && same; ++c)
            {
                if (r_m(r, c) != g_m(r, c)) same = false;
                // Bitwise on the value, not "close": a promotion that shifts a
                // shipped raster by an ulp is still a shipped raster that moved.
                else if (r_m(r, c) &&
                         std::memcmp(&r_d(r, c), &g_d(r, c), sizeof(float)) != 0)
                    same = false;
            }
        check(same, "helper reproduces the original lambda BITWISE on "
                    "shared-mask input");
    }

    std::cout << "\n=====================================================\n"
              << " checks: " << g_checks << ", failures: " << g_failures << "\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
