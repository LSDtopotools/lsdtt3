// test_slope_stability_sinmap.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit test for the SINMAP component (src/cpu/lsdtt_slope_stability_sinmap_cpu).
//
// Every expected value comes from the model's own equations evaluated in the
// test, or from surfaces whose geometry is known in closed form. SINMAP's own
// software is not run.
//
//   * cells: FS_min / FS_max / saturation against the equations written out
//     here; P(FS > 1) against 1e6-draw Monte Carlo of C, t and R/T in each
//     saturation region and at the region boundaries; the SI definition
//     (FS_min where > 1, 0 where FS_max < 1); the deterministic limit; SI and
//     P(stable) never increase with wetter (larger a) cells; flats;
//   * classes at the Table 4 break points; parameter validation;
//   * planar slope raster: D-infinity a is exactly (i + 1) b in row i below
//     the divide, every output equals the single-cell function at that a, and
//     the first saturated row is where a passes (T/R)_min sin(theta);
//   * cone and convergent bowl: D-infinity a against a(r) = r/2 and
//     (R^2 - r^2)/(2r), converging as the grid is refined;
//   * edge contamination: all of a plane draining from the edge, only the
//     edge ring of a cone;
//   * regions, region nodata, exclusion, a missing region id; one thread and
//     four give identical rasters; class areas sum to the region;
//   * SA-plot curves: every point lies on its curve (FS = target, or
//     x a / sin = 1), the vertical part is at the saturated root, and a curve
//     whose target cohesion alone exceeds has no points.
//
// Usage: lsdtt-test-slope-stability-sinmap
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_sinmap_cpu.hpp"
#include "../../src/cpu/lsdtt_wetness_index_cpu.hpp"
#include "../../src/utils/lsdtt_array2d.hpp"
#include "../../src/utils/lsdtt_random.hpp"

using namespace lsdtt;
using lsdtt::random::splitmix64_stream;

namespace
{

int g_failures = 0;
int g_checks = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_checks;
    if (!ok) ++g_failures;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
}

std::string num(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.4g", v);
    return buf;
}

template <class F>
bool throws(F&& f)
{
    try { f(); } catch (const std::invalid_argument&) { return true; }
    return false;
}

// The model written out (manual eqs. 3 and 8).
double fs_equation(double S, double a, double C, double t, double x, double r)
{
    const double th = std::atan(S);
    const double w = std::min(x * a / std::sin(th), 1.0);
    return (C + std::cos(th) * (1.0 - w * r) * t) / std::sin(th);
}

double monte_carlo_p_stable(double S, double a, const sinmap_parameters& p, int n, std::uint64_t seed,
                            double& se)
{
    splitmix64_stream rng(seed);
    const double x1 = 1.0 / p.transmissivity_over_recharge_max_metres;
    const double x2 = 1.0 / p.transmissivity_over_recharge_min_metres;
    int stable = 0;
    for (int i = 0; i < n; ++i)
    {
        const double C = p.dimensionless_cohesion_min
            + (p.dimensionless_cohesion_max - p.dimensionless_cohesion_min) * rng.uniform01();
        const double t = p.friction_coefficient_min
            + (p.friction_coefficient_max - p.friction_coefficient_min) * rng.uniform01();
        const double x = x1 + (x2 - x1) * rng.uniform01();
        if (fs_equation(S, a, C, t, x, p.water_to_soil_density_ratio) > 1.0) ++stable;
    }
    const double pr = static_cast<double>(stable) / n;
    se = std::sqrt(std::max(pr * (1.0 - pr), 1.0 / n) / n);
    return pr;
}

lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.x_min = 500000.0;
    meta.y_min = 4000000.0;
    meta.no_data_value = -9999.0;
    return meta;
}

struct surface
{
    lsdtt_raster_meta meta;
    lsdtt_array2d_f dem;
    lsdtt_array2d_u8 mask;
};

// z(x, y) on cell centres; cells for which keep(x, y) is false are nodata.
template <class Z, class Keep>
surface make_surface(int n, double cs, Z z, Keep keep)
{
    surface s{make_meta(n, n, cs), lsdtt_array2d_f(n, n, -9999.0f), lsdtt_array2d_u8(n, n, 0)};
    for (int row = 0; row < n; ++row)
        for (int col = 0; col < n; ++col)
        {
            const double x = (col + 0.5) * cs;              // from the west edge
            const double y = (n - row - 0.5) * cs;          // from the south edge
            if (!keep(x, y)) continue;
            s.dem(row, col) = static_cast<float>(z(x, y));
            s.mask(row, col) = 1;
        }
    return s;
}

lsdtt_flowgraph_cpu dinf_graph(const surface& s, bool interior_nodata_is_outlet = true)
{
    lsdtt_flowgraph_options o;
    o.method = lsdtt_flowgraph_method::d_infinity;
    o.interior_nodata_is_outlet = interior_nodata_is_outlet;
    o.n_threads = 1;
    return lsdtt_flowgraph_cpu(s.dem, s.mask, s.meta, o);
}

// Mean of (A / b) / a_exact over cells with r_lo < r < r_hi from centre (cx, cy).
template <class Exact>
double mean_area_ratio(const surface& s, const lsdtt_flowgraph_cpu& g, double cx, double cy,
                       double r_lo, double r_hi, Exact exact)
{
    const std::vector<double> A = g.accumulate();
    const int n = s.meta.n_rows;
    const double cs = s.meta.cell_size;
    double sum = 0.0;
    int count = 0;
    for (int row = 0; row < n; ++row)
        for (int col = 0; col < n; ++col)
        {
            if (!s.mask(row, col)) continue;
            const double x = (col + 0.5) * cs, y = (n - row - 0.5) * cs;
            const double r = std::hypot(x - cx, y - cy);
            if (r <= r_lo || r >= r_hi) continue;
            sum += A[static_cast<std::size_t>(row) * n + col] / cs / exact(r);
            ++count;
        }
    return count ? sum / count : 1.0;
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-slope-stability-sinmap\n";
    const sinmap_parameters P;   // SINMAP's defaults
    const double r = P.water_to_soil_density_ratio;
    const double x1 = 1.0 / P.transmissivity_over_recharge_max_metres;
    const double x2 = 1.0 / P.transmissivity_over_recharge_min_metres;

    // ------------------------------------------------------- parameters
    {
        auto bad = [](auto edit) { sinmap_parameters p; edit(p); return throws([&] { validate_sinmap_parameters(p); }); };
        check("parameters: defaults validate", !throws([&] { validate_sinmap_parameters(P); }));
        check("parameters: reversed, negative, zero T/R, r outside [0,1] and NaN all throw",
              bad([](sinmap_parameters& p) { p.dimensionless_cohesion_min = 0.3; })
              && bad([](sinmap_parameters& p) { p.friction_coefficient_min = -0.1; })
              && bad([](sinmap_parameters& p) { p.transmissivity_over_recharge_min_metres = 0.0; })
              && bad([](sinmap_parameters& p) { p.water_to_soil_density_ratio = 1.5; })
              && bad([](sinmap_parameters& p) { p.friction_coefficient_max = std::nan(""); }));
        check("classes: Table 4 break points",
              sinmap_stability_class(1.5001) == 1 && sinmap_stability_class(1.5) == 2
              && sinmap_stability_class(1.25) == 3 && sinmap_stability_class(1.0) == 4
              && sinmap_stability_class(0.5) == 5 && sinmap_stability_class(1e-9) == 5
              && sinmap_stability_class(0.0) == 6
              && sinmap_stability_class(std::numeric_limits<double>::infinity()) == 1);
    }

    // ------------------------------------------------------- single cells
    {
        // FS_min, FS_max and saturation against the equations, over a sweep.
        double worst = 0.0;
        bool sat_ok = true;
        for (double S : {0.05, 0.2, 0.5, 0.84, 1.0, 1.5, 3.0})
            for (double a : {1.0, 30.0, 300.0, 1000.0, 3000.0, 1e4, 1e5})
            {
                const sinmap_cell_result c = sinmap_evaluate_cell(S, a, P);
                const double s = std::sin(std::atan(S));
                const double e_min = fs_equation(S, a, P.dimensionless_cohesion_min, P.friction_coefficient_min, x2, r);
                const double e_max = fs_equation(S, a, P.dimensionless_cohesion_max, P.friction_coefficient_max, x1, r);
                worst = std::max({worst, std::fabs(c.fs_min - e_min) / e_min, std::fabs(c.fs_max - e_max) / e_max});
                const double q1 = x1 * a / s, q2 = x2 * a / s;
                const double sat = (q1 >= 1.0) ? 3.0 : (q2 > 1.0) ? 2.0 : q2;
                sat_ok &= c.saturation == sat;
            }
        check("cell: FS_min and FS_max equal eqs. 3 and 8 at the extreme parameter sets", worst < 1e-14,
              "max relative " + num(worst));
        check("cell: saturation coding (3 always, 2 threshold, else x2 a / sin)", sat_ok);

        // P(FS > 1) against Monte Carlo, in each saturation region and near boundaries.
        struct mc_case { double S, a; const char* where; };
        const double S = 0.84, s = std::sin(std::atan(S));
        const mc_case cases[] = {
            {S, 0.5 * s / x2, "never saturated"},
            {S, 0.999 * s / x2, "just below the x2 saturation line"},
            {S, 1.3 * s / x2, "threshold region"},
            {S, 1.001 * s / x1, "just above the x1 saturation line"},
            {S, 4.0 * s / x1, "always saturated"},
            {0.62, 1200.0, "steeper, threshold"},
            {1.1, 200.0, "steep, never saturated"},
        };
        double worst_z = 0.0;
        std::string detail;
        int n_partial = 0;
        for (const mc_case& mc : cases)
        {
            const sinmap_cell_result c = sinmap_evaluate_cell(mc.S, mc.a, P);
            double se = 0.0;
            const double p_mc = monte_carlo_p_stable(mc.S, mc.a, P, 1000000, 99, se);
            worst_z = std::max(worst_z, std::fabs(c.probability_stable - p_mc) / se);
            if (c.probability_stable > 0.02 && c.probability_stable < 0.98) ++n_partial;
            detail += std::string(mc.where) + " " + num(c.probability_stable) + "/" + num(p_mc) + "; ";
        }
        check("cell: P(FS > 1) matches 1e6-draw Monte Carlo within 5 standard errors in every region",
              worst_z < 5.0 && n_partial >= 4, "worst " + num(worst_z) + " SE; " + detail);

        // SI definition and the 0/1 ends of P.
        bool si_ok = true;
        for (double S2 = 0.02; S2 < 3.0; S2 += 0.037)
            for (double a = 1.0; a < 3e4; a *= 1.7)
            {
                const sinmap_cell_result c = sinmap_evaluate_cell(S2, a, P);
                if (c.fs_min > 1.0) si_ok &= (c.stability_index == c.fs_min) && std::fabs(c.probability_stable - 1.0) < 1e-12;
                else if (c.fs_max < 1.0) si_ok &= (c.stability_index == 0.0) && c.probability_stable < 1e-12;
                else si_ok &= (c.stability_index == c.probability_stable);
                si_ok &= c.stability_class == sinmap_stability_class(c.stability_index);
            }
        check("cell: SI = FS_min above 1, 0 below FS_max < 1, else P(stable); P is 1 / 0 at those ends", si_ok);

        // Wetter cells are never more stable.
        bool mono = true;
        for (double S2 : {0.3, 0.84, 1.4})
        {
            double prev_si = std::numeric_limits<double>::infinity(), prev_p = 2.0;
            for (double a = 1.0; a < 1e5; a *= 1.15)
            {
                const sinmap_cell_result c = sinmap_evaluate_cell(S2, a, P);
                mono &= c.stability_index <= prev_si * (1.0 + 1e-12) + 1e-12;
                mono &= c.probability_stable <= prev_p + 1e-12;
                prev_si = c.stability_index;
                prev_p = c.probability_stable;
            }
        }
        check("cell: SI and P(stable) never increase as a grows (wetter)", mono);

        // Deterministic limit: SI is FS where FS > 1, else 0, and P is 0 or 1.
        sinmap_parameters d;
        d.dimensionless_cohesion_min = d.dimensionless_cohesion_max = 0.0;
        d.friction_coefficient_min = d.friction_coefficient_max = 0.7;
        d.transmissivity_over_recharge_min_metres = d.transmissivity_over_recharge_max_metres = 2500.0;
        bool det = true;
        for (double S2 = 0.1; S2 < 1.5; S2 += 0.013)
            for (double a : {10.0, 500.0, 5000.0})
            {
                const sinmap_cell_result c = sinmap_evaluate_cell(S2, a, d);
                const double fs = fs_equation(S2, a, 0.0, 0.7, 1.0 / 2500.0, 0.5);
                det &= (fs > 1.0) ? std::fabs(c.stability_index - fs) < 1e-12 * fs : c.stability_index == 0.0;
                det &= c.probability_stable == 0.0 || c.probability_stable == 1.0;
            }
        check("cell: deterministic limit (zero-width bounds) gives SI = FS or 0", det);

        const sinmap_cell_result f = sinmap_evaluate_cell(0.0, 100.0, P);
        check("cell: a flat cell is unconditionally stable (class 1, P = 1, saturated)",
              f.flat && f.stability_class == 1 && f.probability_stable == 1.0 && f.saturation == 3.0);
    }

    // ------------------------------------------------------- planar slope
    {
        const int n = 60;
        const double cs = 10.0, S0 = 0.6;
        // Divide along the north edge; z falls to the south.
        const surface s = make_surface(n, cs, [&](double, double y) { return S0 * y; },
                                       [](double, double) { return true; });
        const lsdtt_flowgraph_cpu g = dinf_graph(s);
        lsdtt_raster_cpu slope(g.dinf_slopes());
        lsdtt_raster_cpu area(g.accumulate_raster());
        sinmap_parameters p = P;
        p.transmissivity_over_recharge_min_metres = 200.0;   // saturates within the grid
        p.transmissivity_over_recharge_max_metres = 300.0;
        const sinmap_rasters out = compute_sinmap(slope, area, {}, {{0, p}});

        bool a_exact = true, same = true;
        int first_threshold_row = -1;
        for (int row = 0; row < n - 1; ++row)          // the south row is the outlet
            for (int col = 1; col < n - 1; ++col)
            {
                const std::size_t i = static_cast<std::size_t>(row) * n + col;
                const double a = out.specific_catchment_area.get_data()[i];
                a_exact &= std::fabs(a - (row + 1) * cs) < 1e-4 * (row + 1) * cs;
                const sinmap_cell_result c = sinmap_evaluate_cell(S0, (row + 1) * cs, p);
                same &= std::fabs(out.probability_stable.get_data()[i] - static_cast<float>(c.probability_stable)) < 1e-6
                     && out.stability_class.get_data()[i] == static_cast<float>(c.stability_class)
                     && std::fabs(out.fs_min.get_data()[i] - static_cast<float>(c.fs_min)) < 1e-6f * static_cast<float>(c.fs_min);
                if (col == n / 2 && first_threshold_row < 0 && out.saturation.get_data()[i] >= 2.0f)
                    first_threshold_row = row;
            }
        const double s_theta = std::sin(std::atan(S0));
        const double a_star = p.transmissivity_over_recharge_min_metres * s_theta;   // x2 a / sin = 1
        const int expected_row = static_cast<int>(std::floor(a_star / cs + 1e-12));  // first (i+1) cs > a*
        check("plane: D-infinity a is (i + 1) b in row i below the divide", a_exact);
        check("plane: every output equals the single-cell function at the analytic a", same);
        check("plane: first threshold-saturated row is where a passes (T/R)_min sin(theta)",
              first_threshold_row == expected_row,
              "row " + std::to_string(first_threshold_row) + ", expected " + std::to_string(expected_row));

        const std::vector<std::uint8_t> edge = dinf_edge_contaminated_mask(g, s.meta, s.mask.get_vector());
        bool all = true;
        for (std::size_t i = 0; i < edge.size(); ++i) all &= edge[i] == 1;
        check("edge: every cell of a plane draining from the edge is contaminated", all);
    }

    // ------------------------------------------------ cone and convergent bowl
    {
        // SINMAP's a is A / b with b the cell size, whatever the flow direction
        // (manual Part I section 3.4). A D-infinity flow tube on a plane at angle
        // alpha to the grid is b max(|cos|, |sin|) wide (at 45 degrees each
        // diagonal chain of cells is b / sqrt 2 from the next), so on a plane
        // A / b = L max(|cos alpha|, |sin alpha|), L the flow distance from the
        // divide. Exact for D-infinity on planes; asserted at four angles.
        double worst_plane = 0.0, worst_corrected = 0.0;
        std::string plane_detail, corrected_detail;
        for (double deg : {0.0, 15.0, 30.0, 45.0})
        {
            const double al = deg * M_PI / 180.0, cs = 10.0;
            const int n = 101;
            const surface s = make_surface(n, cs, [&](double x, double y) { return -0.3 * (x * std::cos(al) + y * std::sin(al)); },
                                           [](double, double) { return true; });
            const lsdtt_flowgraph_cpu g = dinf_graph(s);
            const std::vector<double> A = g.accumulate();
            const std::vector<double> width = dinf_contour_width(g);
            const double factor = std::max(std::fabs(std::cos(al)), std::fabs(std::sin(al)));
            double sum = 0.0, sum_corrected = 0.0;
            int cnt = 0;
            for (int row = 20; row < 80; ++row)
                for (int col = 20; col < 80; ++col)
                {
                    const double x = (col + 0.5) * cs, y = (n - row - 0.5) * cs;
                    const double Lx = std::cos(al) > 1e-12 ? x / std::cos(al) : 1e300;
                    const double Ly = std::sin(al) > 1e-12 ? y / std::sin(al) : 1e300;
                    const double L = std::min(Lx, Ly);
                    const std::size_t k = static_cast<std::size_t>(row) * n + col;
                    sum += A[k] / cs / (L * factor);
                    sum_corrected += A[k] / width[k] / L;
                    ++cnt;
                }
            const double ratio = sum / cnt;
            worst_plane = std::max(worst_plane, std::fabs(ratio - 1.0));
            plane_detail += num(deg) + " deg: " + num(ratio) + "; ";
            const double rc = sum_corrected / cnt;
            worst_corrected = std::max(worst_corrected, std::fabs(rc - 1.0));
            corrected_detail += num(deg) + " deg: " + num(rc) + "; ";
        }
        check("oblique plane: A / b = L max(|cos|, |sin|) at 0, 15, 30 and 45 degrees (within 3 %)",
              worst_plane < 0.03, plane_detail);
        check("oblique plane: with the D-infinity contour width, A / w = L at every angle (within 3 %)",
              worst_corrected < 0.03, corrected_detail);

        // Divergent cone: dispersion mixes directions, so A / b tends to
        // (r/2) times the mean of max(|cos|, |sin|) over angles, 2 sqrt(2)/pi.
        const double k_cone = 2.0 * std::sqrt(2.0) / M_PI;
        auto cone_ratio = [](int n)
        {
            const double cs = 1000.0 / n, c = 500.0;
            const surface s = make_surface(n, cs, [&](double x, double y) { return 300.0 - 0.5 * std::hypot(x - c, y - c); },
                                           [](double, double) { return true; });
            const lsdtt_flowgraph_cpu g = dinf_graph(s);
            return mean_area_ratio(s, g, c, c, 0.25 * c, 0.8 * c, [](double rr) { return rr / 2.0; });
        };
        const double c51 = cone_ratio(51), c201 = cone_ratio(201);
        // With the contour width, the mean over directions is 0.9003 times the
        // mean of 1 / max(|cos|, |sin|), (4/pi) ln(1 + sqrt 2) = 1.122, i.e.
        // 1.010 if the uncorrected ratio is direction-independent (as measured).
        auto cone_corrected = [](int n)
        {
            const double cs = 1000.0 / n, c = 500.0;
            const surface s = make_surface(n, cs, [&](double x, double y) { return 300.0 - 0.5 * std::hypot(x - c, y - c); },
                                           [](double, double) { return true; });
            const lsdtt_flowgraph_cpu g = dinf_graph(s);
            const std::vector<double> A = g.accumulate();
            const std::vector<double> w = dinf_contour_width(g);
            double sum = 0.0;
            int cnt = 0;
            for (int row = 0; row < n; ++row)
                for (int col = 0; col < n; ++col)
                {
                    const double x = (col + 0.5) * cs, y = (n - row - 0.5) * cs;
                    const double rr = std::hypot(x - c, y - c);
                    if (rr <= 0.25 * c || rr >= 0.8 * c) continue;
                    const std::size_t k = static_cast<std::size_t>(row) * n + col;
                    sum += A[k] / w[k] / (rr / 2.0);
                    ++cnt;
                }
            return sum / cnt;
        };
        const double cc201 = cone_corrected(201);
        const double k_corrected = k_cone * 4.0 / M_PI * std::log(1.0 + std::sqrt(2.0));
        check("cone: with the contour width the mean (A / w) / (r/2) is within 2 % of 1 (exact)",
              std::fabs(cc201 - 1.0) < 0.02,
              "201 cells: " + num(cc201) + "; " + num(k_corrected) + " predicted if the uncorrected ratio were direction-independent");
        check("cone: mean (A / b) / (r/2) converges to 2 sqrt(2)/pi = " + num(k_cone),
              std::fabs(c201 - k_cone) < std::fabs(c51 - k_cone) && std::fabs(c201 - k_cone) < 0.015,
              "51 cells: " + num(c51) + ", 201 cells: " + num(c201));

        // Convergent bowl: every cell drains to the central sink, so the sink
        // collects the whole domain (exact). The ratio to (R^2 - r^2)/(2r) is
        // direction dependent under convergence and is printed, not asserted.
        {
            const int n = 101;
            const double cs = 10.0, c = 505.0, R = 480.0;
            const surface s = make_surface(n, cs, [&](double x, double y) { return 0.5 * std::hypot(x - c, y - c); },
                                           [&](double x, double y) { return std::hypot(x - c, y - c) < R; });
            const lsdtt_flowgraph_cpu g = dinf_graph(s, false);
            const std::vector<double> A = g.accumulate();
            double total = 0.0, at_sink = 0.0;
            for (std::size_t i = 0; i < A.size(); ++i) if (s.mask.get_vector()[i]) total += cs * cs;
            for (std::int64_t f : g.outlet_indices()) at_sink += A[static_cast<std::size_t>(f)];
            const double bowl = mean_area_ratio(s, g, c, c, 0.25 * R, 0.8 * R,
                                                [R](double rr) { return (R * R - rr * rr) / (2.0 * rr); });
            // Tolerance 1e-6: the D-infinity receiver weights are float32, so each
            // cell's split sums to 1 only to about 1e-7.
            check("bowl: the central sink collects the whole domain (D-infinity conserves area)",
                  std::fabs(at_sink / total - 1.0) < 1e-6 && g.outlet_indices().size() == 1,
                  "sink share " + num(at_sink / total) + ", outlets " + std::to_string(g.outlet_indices().size()) + "; mean (A / b) / ((R^2 - r^2)/(2r)) = " + num(bowl) + " (measured)");
        }

        // Edge contamination on a cone: only cells on the edge ring.
        const int n = 41;
        const surface s = make_surface(n, 10.0, [](double x, double y) { return 100.0 - 0.3 * std::hypot(x - 205.0, y - 205.0); },
                                       [](double, double) { return true; });
        const lsdtt_flowgraph_cpu g = dinf_graph(s);
        const std::vector<std::uint8_t> edge = dinf_edge_contaminated_mask(g, s.meta, s.mask.get_vector());
        bool ring = true;
        for (int row = 0; row < n; ++row)
            for (int col = 0; col < n; ++col)
            {
                const bool on_edge = row == 0 || col == 0 || row == n - 1 || col == n - 1;
                ring &= edge[static_cast<std::size_t>(row) * n + col] == (on_edge ? 1 : 0);
            }
        check("edge: on a cone only the edge ring is contaminated", ring);
    }

    // ------------------------------------------------ relation to TOPMODEL's TWI
    {
        // With the same slope and area, SINMAP's unsaturated wetness at R/T = x is
        //   w = x a / sin(theta),  a = A / b,
        // and lsdtt3's TWI = ln((A / b^2) / tan(theta)) (lsdtt_wetness_index_cpu), so
        //   w = x b exp(TWI) / cos(theta)  exactly: the indices differ only by
        // sin vs tan (the cos factor) and constants. TWI is evaluated in float
        // (its precision contract), hence the 1e-6 tolerance.
        const int n = 40;
        const surface s = make_surface(n, 10.0, [](double x, double y) { return 0.6 * y + 0.1 * x + 4.0 * std::sin(x / 41.0) * std::cos(y / 29.0); },
                                       [](double, double) { return true; });
        const lsdtt_flowgraph_cpu g = dinf_graph(s);
        lsdtt_raster_cpu slope(g.dinf_slopes());
        lsdtt_raster_cpu area(g.accumulate_raster());
        const lsdtt_raster_cpu twi = compute_twi_from_area_and_slope(area, slope, s.meta);
        sinmap_parameters dry = P;   // dry enough that no cell saturates
        dry.transmissivity_over_recharge_min_metres = 1e6;
        dry.transmissivity_over_recharge_max_metres = 2e6;
        const double xw = 1.0 / dry.transmissivity_over_recharge_min_metres;
        double worst = 0.0;
        int compared = 0;
        for (std::size_t i = 0; i < slope.get_data().size(); ++i)
        {
            if (!twi.get_mask()[i] || !(slope.get_data()[i] > 0.0f)) continue;
            const double S = slope.get_data()[i];
            const sinmap_cell_result c = sinmap_evaluate_cell(S, area.get_data()[i] / 10.0, dry);
            if (c.saturation_region != k_sinmap_never_saturated) continue;
            const double from_twi = xw * 10.0 * std::exp(static_cast<double>(twi.get_data()[i])) / std::cos(std::atan(S));
            worst = std::max(worst, std::fabs(c.saturation / from_twi - 1.0));
            ++compared;
        }
        check("TOPMODEL: SINMAP wetness = (R/T) b exp(TWI) / cos(theta) on the same slope and area",
              worst < 1e-6 && compared > 1000, std::to_string(compared) + " cells, max relative " + num(worst));
    }

    // ------------------------------------------------ regions, exclusion, threads
    {
        const int n = 40;
        const surface s = make_surface(n, 10.0, [](double x, double y) { return 0.4 * y + 0.2 * x + 3.0 * std::sin(x / 37.0) * std::cos(y / 23.0); },
                                       [](double, double) { return true; });
        const lsdtt_flowgraph_cpu g = dinf_graph(s);
        lsdtt_raster_cpu slope(g.dinf_slopes());
        lsdtt_raster_cpu area(g.accumulate_raster());
        sinmap_parameters wet = P;
        wet.transmissivity_over_recharge_min_metres = 50.0;
        wet.transmissivity_over_recharge_max_metres = 120.0;
        const std::size_t nn = static_cast<std::size_t>(n) * n;
        std::vector<std::int32_t> region(nn, 1);
        for (std::size_t i = 0; i < nn; ++i) if ((i / n) >= 20) region[i] = 7;
        region[5] = k_sinmap_region_nodata;
        std::vector<std::uint8_t> exclude(nn, 0);
        exclude[n * 10 + 10] = 1;

        const std::map<std::int32_t, sinmap_parameters> table{{1, P}, {7, wet}};
        const sinmap_rasters a1 = compute_sinmap(slope, area, region, table, exclude, 1);
        const sinmap_rasters a4 = compute_sinmap(slope, area, region, table, exclude, 4);

        bool per_region = true;
        for (std::size_t i = 0; i < nn; ++i)
        {
            if (!a1.stability_class.get_mask()[i]) continue;
            const sinmap_cell_result c = sinmap_evaluate_cell(slope.get_data()[i], area.get_data()[i] / 10.0,
                                                              region[i] == 1 ? P : wet);
            per_region &= a1.stability_class.get_data()[i] == static_cast<float>(c.stability_class)
                       && a1.probability_stable.get_data()[i] == static_cast<float>(c.probability_stable);
        }
        check("regions: each cell uses its own region's parameters", per_region);
        check("regions: a region-nodata cell and an excluded cell are nodata; the exclusion is counted",
              !a1.stability_class.get_mask()[5] && !a1.stability_class.get_mask()[n * 10 + 10]
              && a1.stability_class.get_data()[5] == -9999.0f && a1.n_excluded == 1);
        check("threads: 1 and 4 threads give identical rasters",
              a1.stability_index.get_data() == a4.stability_index.get_data()
              && a1.probability_stable.get_data() == a4.probability_stable.get_data()
              && a1.saturation.get_data() == a4.saturation.get_data()
              && a1.stability_class.get_mask() == a4.stability_class.get_mask());

        std::map<std::int32_t, std::int64_t> cells;
        std::map<std::int32_t, double> frac;
        for (const auto& row : sinmap_class_areas(a1, region)) { cells[row.region_id] += row.n_cells; frac[row.region_id] += row.fraction_of_region; }
        std::int64_t valid1 = 0, valid7 = 0;
        for (std::size_t i = 0; i < nn; ++i)
            if (a1.stability_class.get_mask()[i]) (region[i] == 1 ? valid1 : valid7)++;
        check("class areas: counts sum to each region's cells and fractions to 1",
              cells[1] == valid1 && cells[7] == valid7 && std::fabs(frac[1] - 1.0) < 1e-12 && std::fabs(frac[7] - 1.0) < 1e-12);

        std::vector<std::int32_t> bad = region;
        bad[100] = 3;
        check("regions: a valid cell whose region id has no parameter row throws",
              throws([&] { (void)compute_sinmap(slope, area, bad, table); }));
    }

    // ------------------------------------------------ SA-plot curves
    {
        const double slope_max = 2.0, a_min = 10.0, a_max = 1e5;
        const auto pts = sinmap_slope_area_curves(P, slope_max, a_min, a_max, 400);
        double worst_fs = 0.0, worst_sat = 0.0;
        std::map<std::string, int> count;
        for (const auto& pt : pts)
        {
            ++count[pt.curve];
            const double s = std::sin(std::atan(pt.slope));
            if (pt.curve == "always_saturated_above") worst_sat = std::max(worst_sat, std::fabs(x1 * pt.specific_catchment_area_metres / s - 1.0));
            else if (pt.curve == "never_saturated_below") worst_sat = std::max(worst_sat, std::fabs(x2 * pt.specific_catchment_area_metres / s - 1.0));
            else
            {
                const bool is_max = pt.curve == "fs_max_1.0";
                const double target = (pt.curve == "fs_min_1.5") ? 1.5 : 1.0;
                const double fs = is_max
                    ? fs_equation(pt.slope, pt.specific_catchment_area_metres, P.dimensionless_cohesion_max, P.friction_coefficient_max, x1, r)
                    : fs_equation(pt.slope, pt.specific_catchment_area_metres, P.dimensionless_cohesion_min, P.friction_coefficient_min, x2, r);
                worst_fs = std::max(worst_fs, std::fabs(fs - target));
            }
        }
        check("curves: every FS-curve point has FS equal to its target", worst_fs < 1e-9 && count["fs_min_1.0"] > 50,
              "max |FS - target| " + num(worst_fs));
        check("curves: every saturation-line point has x a / sin(theta) = 1", worst_sat < 1e-12 && count["never_saturated_below"] > 50);
        check("curves: all five curves present for SINMAP's defaults",
              count.size() == 5 && count["fs_max_1.0"] > 10 && count["fs_min_1.5"] > 10);

        const double Ssat = sinmap_saturated_slope_for_safety_factor(1.2, 0.1, 0.8, 0.5);
        const double th = std::atan(Ssat);
        check("curves: saturated root satisfies FS(w = 1) = target",
              std::fabs((0.1 + std::cos(th) * 0.5 * 0.8) / std::sin(th) - 1.2) < 1e-12);
        check("curves: no root when cohesion alone exceeds the target",
              std::isinf(sinmap_saturated_slope_for_safety_factor(1.0, 1.2, 0.8, 0.5)));

        sinmap_parameters strong = P;
        strong.dimensionless_cohesion_min = 2.0;
        strong.dimensionless_cohesion_max = 2.5;
        int n_fs = 0;
        for (const auto& pt : sinmap_slope_area_curves(strong, slope_max, a_min, a_max, 200))
            if (pt.curve.rfind("fs_", 0) == 0) ++n_fs;
        check("curves: a target that cohesion alone exceeds has no points", n_fs == 0);
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
