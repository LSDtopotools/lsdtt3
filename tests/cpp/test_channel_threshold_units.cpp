// test_channel_threshold_units.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// What `incision_threshold_area` is actually compared against.
//
// It was documented as "Minimum drainage area (m^2)". It is not, on the default
// setting: with `divide_area_by_dx` true the comparison is against SPECIFIC
// catchment area a = A/dx, which is a LENGTH IN METRES. The two conventions are
// dimensionally different, and the default is the one the documentation did not
// describe.
//
// (An earlier note in the plan called the behaviour "a pixel count". That is
// also wrong -- a pixel count would be A/dx^2. It is A/dx.)
//
// This is a units error that cannot show up as a crash or a NaN: it selects a
// different channel network, which looks like a modelling choice rather than a
// mistake. So the semantics are pinned here arithmetically, on a fixture whose
// upstream pixel counts are known exactly.
//
//   U1  the compared quantity, both conventions, from a known area
//   U2  THE GRID-DEPENDENCE TRAP: the same threshold selects a different
//       number of upstream pixels at a different cell size
//   U3  a value chosen as though it were m^2 is wrong by a factor of dx
//   U4  the mask agrees with the quantity, cell by cell
//   U5  the report converts consistently, and round-trips
//   U6  a non-positive threshold selects nothing, in both conventions
//
// Usage: lsdtt-test-channel-threshold-units
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_stream_power_cpu.hpp"

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

}  // namespace


int main()
{
    using namespace lsdtt;

    std::cout << "=====================================================\n"
              << " incision_threshold_area: what it is compared against\n"
              << "=====================================================\n";

    // -------------------------------------------------------------------------
    std::cout << "\nU1: the compared quantity\n";
    {
        const double dx = 10.0;
        const double n_upstream = 250.0;
        const double A = n_upstream * dx * dx;      // 25000 m^2

        check_close(channel_threshold_quantity(A, dx, false), 25000.0, 1e-9,
                    "divide_area_by_dx OFF -> map area in m^2");
        check_close(channel_threshold_quantity(A, dx, true), 2500.0, 1e-9,
                    "divide_area_by_dx ON  -> SPECIFIC area A/dx, in METRES");

        // Not a pixel count, which is the other thing it gets mistaken for.
        check(std::fabs(channel_threshold_quantity(A, dx, true) - n_upstream) > 1.0,
              "and it is NOT a pixel count (that would be A/dx^2 = 250)");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nU2: THE GRID-DEPENDENCE TRAP\n";
    {
        // One threshold, two cell sizes. Under the default convention the
        // number of upstream PIXELS it demands scales with 1/dx, so a .param
        // file moved between grids silently asks for a different network.
        const double T = 1000.0;

        const auto coarse = describe_channel_threshold(T, 10.0, true);
        const auto fine   = describe_channel_threshold(T,  1.0, true);

        check_close(coarse.equivalent_pixels, 100.0, 1e-9,
                    "T=1000 on a 10 m grid demands 100 upstream pixels");
        check_close(fine.equivalent_pixels, 1000.0, 1e-9,
                    "the SAME T=1000 on a 1 m grid demands 1000 pixels");
        check_close(fine.equivalent_pixels / coarse.equivalent_pixels, 10.0, 1e-9,
                    "a 10x finer grid asks for a 10x larger network in pixels");

        // Under the total-area convention the pixel count scales as 1/dx^2, so
        // the trap is worse, not absent -- the conventions differ in HOW they
        // depend on the grid, not in whether they do.
        const auto coarse_A = describe_channel_threshold(T, 10.0, false);
        const auto fine_A   = describe_channel_threshold(T,  1.0, false);
        check_close(fine_A.equivalent_pixels / coarse_A.equivalent_pixels, 100.0, 1e-9,
                    "under map area the same move costs a factor of 100");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nU3: a value picked as though it were m^2\n";
    {
        // Someone wants 10 upstream pixels on a 10 m grid: 10 * 100 = 1000 m^2.
        // Entered as 1000 under the default convention it means 1000 m of
        // specific area = 100 pixels. Ten times too many.
        const double dx = 10.0;
        const auto got = describe_channel_threshold(1000.0, dx, true);
        check_close(got.equivalent_pixels, 100.0, 1e-9,
                    "'1000' means 100 pixels, not the 10 that 1000 m^2 implies");
        check_close(got.equivalent_pixels / 10.0, dx, 1e-9,
                    "the error is exactly a factor of dx");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nU4: the mask agrees with the quantity, cell by cell\n";
    {
        const int nr = 21, nc = 21;
        const double dx = 10.0;

        lsdtt_raster_meta meta;
        meta.n_rows = nr; meta.n_cols = nc; meta.cell_size = dx;
        meta.x_min = 0.0; meta.y_min = 0.0; meta.no_data_value = -9999.0;

        std::vector<float>        z(static_cast<std::size_t>(nr) * nc);
        std::vector<std::uint8_t> m(z.size(), 1u);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[static_cast<std::size_t>(r) * nc + c] =
                    static_cast<float>(500.0 - 1.0 * r + 0.6 * std::fabs(c - 10.0));

        const lsdtt_raster_cpu dem(std::move(z), std::move(m), meta);
        lsdtt_flowinfo_options fi;
        fi.boundary_conditions = "oooo";
        const lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                                  dem.get_meta(), fi);
        const std::vector<double> area = flow.accumulate();

        for (bool by_dx : {false, true})
        {
            const double T = by_dx ? 200.0 : 2000.0;
            const auto mask = channel_mask_from_area_threshold(
                flow, area, T, dx, by_dx);

            bool agrees = true;
            std::int64_t n_chan = 0;
            for (std::size_t i = 0; i < mask.size(); ++i)
            {
                if (!flow.is_valid_cell(static_cast<std::int64_t>(i)))
                {
                    if (mask[i]) agrees = false;
                    continue;
                }
                const double q = channel_threshold_quantity(area[i], dx, by_dx);
                if ((q >= T) != (mask[i] != 0u)) agrees = false;
                if (mask[i]) ++n_chan;
            }
            check(agrees, std::string("mask == (quantity >= T) everywhere, "
                                      "divide_area_by_dx ") +
                          (by_dx ? "on" : "off"));
            check(n_chan > 0, "and the fixture actually has channels");
        }
    }

    // -------------------------------------------------------------------------
    std::cout << "\nU5: the report converts consistently\n";
    {
        const double dx = 5.0;
        for (bool by_dx : {false, true})
        {
            const auto r = describe_channel_threshold(400.0, dx, by_dx);
            check_close(r.map_area_m2, r.specific_area_m * dx, 1e-9,
                        std::string("map = specific * dx (by_dx ") +
                        (by_dx ? "on)" : "off)"));
            check_close(r.equivalent_pixels, r.map_area_m2 / (dx * dx), 1e-9,
                        "pixels = map / dx^2");
        }

        // The two conventions describe the SAME network when the threshold is
        // converted -- which is the check that the conversion is the right one
        // rather than merely self-consistent.
        const auto as_specific = describe_channel_threshold(1000.0, dx, true);
        const auto as_map      = describe_channel_threshold(
                                     as_specific.map_area_m2, dx, false);
        check_close(as_map.equivalent_pixels, as_specific.equivalent_pixels, 1e-9,
                    "converting the threshold between conventions preserves the "
                    "network");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nU6: a non-positive threshold selects nothing\n";
    {
        const int nr = 9, nc = 9;
        lsdtt_raster_meta meta;
        meta.n_rows = nr; meta.n_cols = nc; meta.cell_size = 1.0;
        meta.x_min = 0.0; meta.y_min = 0.0; meta.no_data_value = -9999.0;
        std::vector<float>        z(static_cast<std::size_t>(nr) * nc);
        std::vector<std::uint8_t> m(z.size(), 1u);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[static_cast<std::size_t>(r) * nc + c] =
                    static_cast<float>(100.0 - r);
        const lsdtt_raster_cpu dem(std::move(z), std::move(m), meta);
        lsdtt_flowinfo_options fi;
        fi.boundary_conditions = "oooo";
        const lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                                  dem.get_meta(), fi);
        const auto area = flow.accumulate();

        for (bool by_dx : {false, true})
        {
            const auto mask = channel_mask_from_area_threshold(
                flow, area, 0.0, 1.0, by_dx);
            std::int64_t n = 0;
            for (auto v : mask) n += (v != 0u);
            check(n == 0, std::string("threshold 0 selects no channels (by_dx ") +
                          (by_dx ? "on)" : "off)"));
        }
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
