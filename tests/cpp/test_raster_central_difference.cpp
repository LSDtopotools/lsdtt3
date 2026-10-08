// test_raster_central_difference.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation of the two scalar users of the four-cardinal centred-
// difference gradient after consolidation C18 moved the stencil into
// lsdtt_raster_central_difference_cpu.hpp:
//
//   - check_vertical_plausibility (lsdtt_raster_cpu), every report field;
//   - median_interior_gradient (lsdtt_lem_timestep_cpu).
//
// Both functions as they were before the change are kept below VERBATIM as the
// reference, and the results are compared bit for bit on DEMs with nodata holes,
// nodata edges, flat areas and steep areas.
//
// Usage: lsdtt-test-raster-central-difference
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_lem_timestep_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/utils/lsdtt_random.hpp"
#include "../../src/utils/lsdtt_statistics.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

bool same_bits(double a, double b)
{
    std::uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

// ---- reference: lsdtt_raster_cpu.cpp check_vertical_plausibility, verbatim ----
vertical_plausibility_report reference_check_vertical_plausibility(
    const lsdtt_raster_cpu& dem,
    double                  median_slope_limit_deg)
{
    vertical_plausibility_report rep;

    const int    nr = dem.n_rows();
    const int    nc = dem.n_cols();
    const double cs = dem.cell_size();
    if (nr < 3 || nc < 3 || cs <= 0.0) return rep;

    const std::vector<float>&   d = dem.get_data();
    const std::vector<uint8_t>& m = dem.get_mask();

    std::vector<float> slopes;
    slopes.reserve(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc) / 4);

    float zmin =  std::numeric_limits<float>::max();
    float zmax = -std::numeric_limits<float>::max();

    for (int r = 1; r < nr - 1; ++r)
    {
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            if (!m[i]) continue;

            if (d[i] < zmin) zmin = d[i];
            if (d[i] > zmax) zmax = d[i];

            const std::size_t n = i - nc, s = i + nc, w = i - 1, e = i + 1;
            if (!m[n] || !m[s] || !m[w] || !m[e]) continue;

            const double dzdx = (static_cast<double>(d[e]) - d[w]) / (2.0 * cs);
            const double dzdy = (static_cast<double>(d[s]) - d[n]) / (2.0 * cs);
            slopes.push_back(static_cast<float>(std::hypot(dzdx, dzdy)));
        }
    }

    if (slopes.empty()) return rep;

    const auto pct = [&slopes](double p) -> double
    {
        return compute_quantile_floor_in_place(slopes, p / 100.0);
    };

    const auto deg = [](double g) { return std::atan(g) * 180.0 / M_PI; };

    rep.n_sampled                = slopes.size();
    const double med             = pct(50.0);
    rep.median_slope_deg         = deg(med);
    rep.p90_slope_deg            = deg(pct(90.0));
    rep.p99_slope_deg            = deg(pct(99.0));
    rep.median_slope_deg_if_feet = deg(med * 0.3048);
    rep.relief_m                 = static_cast<double>(zmax - zmin);
    rep.implausible              = rep.median_slope_deg > median_slope_limit_deg;

    return rep;
}

// ---- reference: lsdtt_lem_timestep_cpu.cpp median_interior_gradient, verbatim ----
double reference_median_interior_gradient(const lsdtt_raster_cpu& dem)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nrows = meta.n_rows, ncols = meta.n_cols;
    const double dx = meta.cell_size;
    if (nrows < 3 || ncols < 3 || !(dx > 0.0)) return 0.0;

    const std::vector<float>&   z    = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    std::vector<double> g;
    g.reserve(static_cast<std::size_t>(nrows) * ncols / 2);

    for (int r = 1; r < nrows - 1; ++r)
        for (int c = 1; c < ncols - 1; ++c)
        {
            const std::size_t f = static_cast<std::size_t>(r) * ncols + c;
            const std::size_t fn = f - ncols, fs = f + ncols, fw = f - 1, fe = f + 1;
            if (!mask[f] || !mask[fn] || !mask[fs] || !mask[fw] || !mask[fe]) continue;
            const double gx = (static_cast<double>(z[fe]) - static_cast<double>(z[fw]))
                            / (2.0 * dx);
            const double gy = (static_cast<double>(z[fs]) - static_cast<double>(z[fn]))
                            / (2.0 * dx);
            g.push_back(std::sqrt(gx * gx + gy * gy));
        }

    return compute_median_in_place(g, 0.0);
}

// ---- fixtures ----

lsdtt_raster_cpu make_dem(int nr, int nc, double cell, std::uint64_t seed, double hole_fraction,
                          double relief)
{
    lsdtt_raster_meta meta;
    meta.n_rows = nr;
    meta.n_cols = nc;
    meta.cell_size = cell;
    meta.no_data_value = -9999.0;
    random::splitmix64_stream rng(seed);
    std::vector<float> z(static_cast<std::size_t>(nr) * nc);
    std::vector<std::uint8_t> m(z.size(), 1);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            z[i] = static_cast<float>(relief * (std::sin(r * 0.11) * std::cos(c * 0.07) + rng.uniform01()));
            if (rng.uniform01() < hole_fraction) m[i] = 0;
            // A flat block and a nodata strip on the west edge.
            if (r > nr / 2 && r < nr / 2 + 6 && c > nc / 3 && c < nc / 3 + 8) z[i] = 12.5f;
            if (c < 2 && r % 5 == 0) m[i] = 0;
        }
    return lsdtt_raster_cpu(std::move(z), std::move(m), meta);
}

void compare(const std::string& label, const lsdtt_raster_cpu& dem)
{
    const vertical_plausibility_report a = check_vertical_plausibility(dem, 30.0);
    const vertical_plausibility_report b = reference_check_vertical_plausibility(dem, 30.0);
    check(label + ": check_vertical_plausibility report bit-identical (n_sampled " +
              std::to_string(b.n_sampled) + ")",
          same_bits(a.median_slope_deg, b.median_slope_deg) && same_bits(a.p90_slope_deg, b.p90_slope_deg) &&
              same_bits(a.p99_slope_deg, b.p99_slope_deg) && same_bits(a.relief_m, b.relief_m) &&
              same_bits(a.median_slope_deg_if_feet, b.median_slope_deg_if_feet) &&
              a.n_sampled == b.n_sampled && a.implausible == b.implausible);
    check(label + ": median_interior_gradient bit-identical",
          same_bits(median_interior_gradient(dem), reference_median_interior_gradient(dem)));
}

} // namespace

int main()
{
    std::cout << "== C18 four-cardinal centred-difference gradient ==\n";
    for (std::uint64_t seed = 1; seed <= 6; ++seed)
    {
        compare("200x170 gentle, seed " + std::to_string(seed), make_dem(200, 170, 10.0, seed, 0.02, 3.0));
        compare("200x170 steep, seed " + std::to_string(seed), make_dem(200, 170, 2.0, seed, 0.05, 40.0));
    }
    compare("odd cell size, many holes", make_dem(97, 61, 0.3048, 99, 0.30, 5.0));
    compare("3x3", make_dem(3, 3, 1.0, 7, 0.0, 1.0));
    compare("all holes", make_dem(20, 20, 1.0, 8, 1.0, 1.0));
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
