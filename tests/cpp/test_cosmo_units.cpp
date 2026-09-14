// test_cosmo_units.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// The metre <-> g/cm2 conversions (lsdtt_crn_units.hpp).
//
// These are four-line functions, and pinning them looks like ceremony until you
// remember what they are for. The CRN core works in MASS DEPTH; metres appear
// only at the boundaries; every conversion carries a density. Until 2026-08-05
// `x * rho * 0.1` was written out by hand at eight sites, each independently
// choosing a density -- and the density split of 2026-08-04 exists precisely
// because one of those choices was wrong.
//
//   U1  round trips, per material density
//   U2  the documented constants, computed rather than asserted from memory
//   U3  THE FACTOR THE DENSITY SPLIT WAS ABOUT: a soil scar carried at bedrock
//       density is 1.77x too deep in mass. This is the number in the docs and
//       it should be reproducible from the conversions alone
//   U4  rate conversions compose the way the depth ones do
//   U5  mm/kyr <-> m/yr, which carries no density and must not acquire one
//   U6  g/cm2/yr -> kg/m2/yr is exactly 10, the delivery-weight conversion
//   U7  EXPRESSION ORDER is preserved bit-for-bit against the hand-written
//       form the functions replaced -- the migration claim, tested
//
// Usage: lsdtt-test-cosmo-units
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

#include "../../src/cpu/lsdtt_crn_units.hpp"

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

void check_close(double got, double expected, double rel_tol,
                 const std::string& what)
{
    const double denom = std::fabs(expected) > 0.0 ? std::fabs(expected) : 1.0;
    const double rel = std::fabs(got - expected) / denom;
    ++g_checks;
    if (!(rel <= rel_tol))
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n        got " << got
                  << ", expected " << expected << " (rel " << rel << ")\n";
    }
    else
    {
        std::cout << "  ok  : " << what << "  (" << got << ")\n";
    }
}

// The four densities the codebase names, per the density convention.
constexpr double kRock     = 2650.0;
constexpr double kSoil     = 1500.0;
constexpr double kAlluvium = 2200.0;

}  // namespace


int main()
{
    using namespace lsdtt;

    std::cout << "=====================================================\n"
              << " CRN unit conversions\n"
              << "=====================================================\n";

    // -------------------------------------------------------------------------
    std::cout << "\nU1: depth round trips, per material\n";
    for (double rho : {kRock, kSoil, kAlluvium})
        for (double m : {0.05, 0.5, 3.0, 40.0})
        {
            const double back = g_cm2_to_metres(metres_to_g_cm2(m, rho), rho);
            check_close(back, m, 1.0e-12,
                        "m -> g/cm2 -> m at rho = " + std::to_string(int(rho)));
        }

    // -------------------------------------------------------------------------
    std::cout << "\nU2: the documented constants, computed\n";
    {
        // 1 m of rock at 2650 kg/m3 = 265 g/cm2. This is the number that makes
        // "a metre of rock removes essentially all spallation" true, since the
        // spallation attenuation length is 160 g/cm2.
        check_close(metres_to_g_cm2(1.0, kRock), 265.0, 1.0e-12,
                    "1 m of rock = 265 g/cm2");
        check_close(metres_to_g_cm2(1.0, kSoil), 150.0, 1.0e-12,
                    "1 m of soil at 1500 = 150 g/cm2");
        // 10 m of ice at 917 kg/m3 = 917 g/cm2, quoted in the ice section.
        check_close(metres_to_g_cm2(10.0, 917.0), 917.0, 1.0e-12,
                    "10 m of ice at 917 = 917 g/cm2");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nU3: the factor the density split was about\n";
    {
        // A soil scar carried at the BEDROCK density is 2650/1500 = 1.7667x too
        // deep in mass. The docs quote 1.77; it should fall out of the
        // conversions rather than being a remembered number.
        const double d_m = 1.0;
        const double as_soil = metres_to_g_cm2(d_m, kSoil);
        const double as_rock = metres_to_g_cm2(d_m, kRock);
        check_close(as_rock / as_soil, 1.7667, 1.0e-3,
                    "carrying a soil scar at rock density overstates mass depth "
                    "by 1.77x");

        // And the consequence: spallation at the scar base, exp(-d/160).
        const double p_soil = std::exp(-as_soil / 160.0);
        const double p_rock = std::exp(-as_rock / 160.0);
        check_close(p_soil / p_rock, 2.05, 0.02,
                    "which changes spallation production at the scar base by "
                    "2.05x");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nU4: rate conversions compose like the depth ones\n";
    for (double rho : {kRock, kSoil})
        for (double r : {1.0e-5, 1.0e-4, 1.0e-3})
        {
            const double back = eff_erate_to_m_yr(m_yr_to_eff_erate(r, rho), rho);
            check_close(back, r, 1.0e-12, "m/yr -> g/cm2/yr -> m/yr");
        }
    {
        // A rate conversion is the same arithmetic as a depth conversion, and
        // must stay that way: an erosion rate is a depth per unit time.
        check(m_yr_to_eff_erate(0.001, kRock) == metres_to_g_cm2(0.001, kRock),
              "m_yr_to_eff_erate and metres_to_g_cm2 are the same arithmetic");
        // 0.1 mm/yr of rock = 1e-4 m/yr = 0.0265 g/cm2/yr, the value that
        // appears throughout the cosmogenic fixtures.
        check_close(m_yr_to_eff_erate(1.0e-4, kRock), 0.0265, 1.0e-12,
                    "0.1 mm/yr of rock = 0.0265 g/cm2/yr");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nU5: mm/kyr <-> m/yr carries NO density\n";
    {
        check_close(mm_kyr_to_m_yr(1000.0), 1.0e-3, 1.0e-15,
                    "1000 mm/kyr = 1e-3 m/yr");
        check_close(m_yr_to_mm_kyr(mm_kyr_to_m_yr(25.0)), 25.0, 1.0e-12,
                    "mm/kyr round trip");
        // The landslide series prescribes a rate in mm/kyr and the CRN side
        // wants g/cm2/yr; the density enters ONLY at the second step. Composing
        // them the other way round is meaningless, and this pins the split.
        check_close(m_yr_to_eff_erate(mm_kyr_to_m_yr(25.0), kSoil), 0.00375,
                    1.0e-12,
                    "25 mm/kyr of soil = 0.00375 g/cm2/yr (the Niemi fixture)");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nU6: the delivery-weight conversion is exactly 10\n";
    {
        check(eff_erate_to_kg_m2_yr(1.0) == 10.0,
              "1 g/cm2/yr = 10 kg/m2/yr, exactly");
        check_close(eff_erate_to_kg_m2_yr(0.002), 0.02, 1.0e-15,
                    "and it scales linearly");
    }

    // -------------------------------------------------------------------------
    // U7 -- the migration claim.
    //
    // Every site replaced on 2026-08-05 wrote the multiplication in a specific
    // order. `0.1 * rho * m` and `m * rho * 0.1` are not the same double, and
    // the refactor oracle compares file hashes, so "the same value" is not good
    // enough -- it has to be the same BITS.
    // -------------------------------------------------------------------------
    std::cout << "\nU7: expression order preserved bit-for-bit\n";
    {
        bool all_same = true;
        for (double rho : {kRock, kSoil, kAlluvium, 917.0})
            for (double m : {0.017, 0.5, 1.0/3.0, 2.718281828, 137.035})
            {
                const double hand = m * rho * 0.1;      // as it was written
                const double via  = metres_to_g_cm2(m, rho);
                if (std::memcmp(&hand, &via, sizeof(double)) != 0)
                    all_same = false;

                const double hand_inv = (m * rho * 0.1) / (rho * 0.1);
                const double via_inv  = g_cm2_to_metres(metres_to_g_cm2(m, rho),
                                                        rho);
                if (std::memcmp(&hand_inv, &via_inv, sizeof(double)) != 0)
                    all_same = false;
            }
        check(all_same,
              "metres_to_g_cm2 / g_cm2_to_metres reproduce the hand-written "
              "expressions BITWISE");

        bool ten_same = true;
        for (double e : {0.001, 0.0265, 1.0/7.0})
        {
            const double hand = e * 10.0;
            const double via  = eff_erate_to_kg_m2_yr(e);
            if (std::memcmp(&hand, &via, sizeof(double)) != 0) ten_same = false;
        }
        check(ten_same, "eff_erate_to_kg_m2_yr reproduces `x * 10.0` bitwise");

        bool mm_same = true;
        for (double r : {25.0, 7.5, 1.0/3.0})
        {
            const double hand = r * 1.0e-6;
            const double via  = mm_kyr_to_m_yr(r);
            if (std::memcmp(&hand, &via, sizeof(double)) != 0) mm_same = false;
        }
        check(mm_same, "mm_kyr_to_m_yr reproduces `x * 1.0e-6` bitwise");
    }

    std::cout << "\n=====================================================\n"
              << " checks: " << g_checks << ", failures: " << g_failures << "\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
