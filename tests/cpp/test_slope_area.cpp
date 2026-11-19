// test_slope_area.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// d8_flow_path_gradient (lsdtt_slope_area_cpu) and the options form of
// channel_slope_at it uses (plan_slope_area_gradient_and_low_concavity.md,
// step 1b).
//
// Cases:
//   1. A plane falling east at S = 0.05 (0.5 m per 10 m cell, exact in float):
//      the D8 path is the row, and the gradient is S at every interior pixel
//      for windows of 4, 6 and 10 cells.
//   2. A plane falling south-east: the D8 path is the diagonal, and the
//      gradient along it (steps of 10 sqrt 2 m) is the plane's diagonal
//      gradient S.
//   3. One-sided windows: at the divide (west edge, no donors) and at the
//      outlet (east edge) the gradient is S; with a window longer than twice
//      the whole path neither side spans half of it, so the pixel is nodata.
//   4. Exclusion: a 1 m bump raised on three cells of the row changes the
//      gradient next to it; with those cells excluded the gradient is S again
//      and the bump cells themselves are nodata.
//   5. Donor ranking (the valley of test_channel_gradient: the floor falls
//      east at 0.01, a tributary at 0.05 joins at (10, 30)): with the
//      contributing-pixel ranking the walk up from (10, 32) follows the floor
//      (gradient 0.01); with a rank raster that favours the tributary it goes
//      up the tributary and the gradient is steeper.
//   6. Refusals: a window under four cell sizes; a raster off the grid.
//   7. n_threads 1 and 4 give bit-identical rasters.
//   8. conditioned_pixel_mask: on a plane with a 3 x 3 depression, a pit
//      fill marks exactly the cells it raised as filled (1); a least-cost
//      breach marks the cells it lowered as carved (2) and those it raised as
//      filled; equal surfaces mark nothing; a raster off the grid is refused.
//   9. F1 (plan section 10): on 30 m cells with a 150 m window a one-sided
//      window (divide, outlet) reaches 60 m of the 75 m half-window; with
//      the span rule of window/2 less one diagonal step it is valid, and the
//      gradient is exact.
//  10. fit_log_slope_area on arrays: S = 2 A^-0.5 gives theta 0.5 (r2 1);
//      constant S gives theta 0 (r2 NaN); S = A^0.1 gives theta -0.1; equal
//      areas or fewer than min_pixels pairs are invalid; S <= 0 is dropped.
//  11. find_low_concavity with given area and gradient rasters on an
//      east-draining plane (every row a channel): the flag rule |theta| <
//      0.15 for theta 0.5, 0.1, -0.1, -0.3; the window is cut at the outlet
//      (too few pixels: flag -1); an excluded pixel has no fit; excluded
//      window pixels are left out of the fit.
//  12. junctions: on the valley with a tributary, the downstream window of a
//      floor pixel upstream of the confluence counts one junction.
//  13. integration: a valley whose floor area grows downstream, with the long
//      profile built so S = k A^-theta; theta from the D8 flow-path gradient
//      (window 4 cells) is within 0.05 of 0.5 and of 0.05 on the middle of
//      the floor.
//  14. refusals: window_pixels < 2, min_pixels < 3 or > window, negative
//      threshold, rasters off the grid.
//  15. segments (build_low_concavity_segments) on the valley with a
//      tributary, flags set by hand: the segments partition the channel
//      pixels; three links (floor above the confluence, tributary, floor
//      below); a flag change splits a link; a flagged run shorter than
//      min_run_pixels is relabelled 0 and merged; Strahler order 1 above
//      and 2 below the confluence; pixels outside `include` are left out;
//      a segment's length counts the step to the next pixel down.
//
// Usage: lsdtt-test-slope-area
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_channel_gradient_cpu.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_breach_least_cost_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_area_cpu.hpp"
#include "../../src/utils/lsdtt_regression.hpp"

namespace
{
int g_n = 0, g_fail = 0;
void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n;
    if (!ok) ++g_fail;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}
std::string s(double v) { return std::to_string(v); }

constexpr double k_cs = 10.0;

lsdtt::lsdtt_raster_meta make_meta(int nr, int nc)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = nr;
    m.n_cols = nc;
    m.x_min = 0.0;
    m.y_min = 0.0;
    m.cell_size = k_cs;
    m.no_data_value = -9999.0f;
    return m;
}

lsdtt::lsdtt_raster_cpu make_raster(int nr, int nc, const std::vector<float>& v)
{
    std::vector<float> data = v;
    std::vector<std::uint8_t> mask(data.size(), 1);
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), make_meta(nr, nc));
}

std::unique_ptr<lsdtt::lsdtt_flow_base> route(const lsdtt::lsdtt_raster_cpu& dem)
{
    lsdtt::lsdtt_flowinfo_options fo;
    fo.boundary_conditions = "oooo";
    fo.build_s_vector = true;
    fo.build_topological_order = true;
    fo.n_threads = 1;
    return lsdtt::make_flow_base(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fo);
}

// D8 drainage area (m^2) from the flow, as the run's area raster would be.
lsdtt::lsdtt_raster_cpu d8_area(const lsdtt::lsdtt_flow_base& flow, const lsdtt::lsdtt_raster_cpu& dem)
{
    const auto& m = dem.get_meta();
    std::vector<float> a(dem.get_data().size(), 0.0f);
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = static_cast<float>(flow.s_vector_size(static_cast<std::int64_t>(i)) * m.cell_size * m.cell_size);
    return make_raster(m.n_rows, m.n_cols, a);
}

float at(const lsdtt::lsdtt_raster_cpu& r, int row, int col)
{
    return r.get_data()[static_cast<std::size_t>(row * r.get_meta().n_cols + col)];
}
bool valid(const lsdtt::lsdtt_raster_cpu& r, int row, int col)
{
    return r.get_mask()[static_cast<std::size_t>(row * r.get_meta().n_cols + col)] != 0;
}
bool near(double a, double b, double tol) { return std::abs(a - b) <= tol; }

} // namespace

int main()
{
    using namespace lsdtt;
    const double S = 0.05;

    // ---- 1. east-falling plane ------------------------------------------------
    {
        const int nr = 9, nc = 30;
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(200.0 - S * k_cs * c);
        const auto dem = make_raster(nr, nc, z);
        const auto flow = route(dem);
        const auto area = d8_area(*flow, dem);
        check("plane: the D8 path is the row", flow->receiver_of(4 * nc + 10) == 4 * nc + 11);
        for (const int cells : {4, 6, 10})
        {
            lsdtt_d8_flow_path_gradient_options o;
            o.window_metres = cells * k_cs;
            o.n_threads = 1;
            const auto g = d8_flow_path_gradient(*flow, dem, area, o);
            double worst = 0.0;
            bool all_valid = true;
            for (int c = 1; c < nc - 1; ++c)
            {
                if (!valid(g, 4, c)) { all_valid = false; continue; }
                worst = std::max(worst, std::abs(at(g, 4, c) - S));
            }
            check("plane, window " + std::to_string(cells) + " cells: gradient S on the whole interior row",
                  all_valid && worst < 1e-6, "max |g - S| " + s(worst));
        }

        // ---- 3. one-sided windows ---------------------------------------------
        lsdtt_d8_flow_path_gradient_options o;
        o.window_metres = 8 * k_cs;
        o.n_threads = 1;
        const auto g = d8_flow_path_gradient(*flow, dem, area, o);
        check("divide (west edge): one-sided window, gradient S",
              valid(g, 4, 0) && near(at(g, 4, 0), S, 1e-6), s(at(g, 4, 0)));
        check("outlet (east edge): one-sided window, gradient S",
              valid(g, 4, nc - 1) && near(at(g, 4, nc - 1), S, 1e-6), s(at(g, 4, nc - 1)));
        o.window_metres = 2.0 * nc * k_cs + 4 * k_cs;   // half of it longer than the whole row
        const auto g_long = d8_flow_path_gradient(*flow, dem, area, o);
        bool none = true;
        for (int c = 0; c < nc; ++c) none = none && !valid(g_long, 4, c);
        check("window longer than twice the path: every pixel of the row is nodata", none);

        // ---- 6. refusals --------------------------------------------------------
        o.window_metres = 3.9 * k_cs;
        bool threw = false;
        try { (void)d8_flow_path_gradient(*flow, dem, area, o); } catch (const std::invalid_argument&) { threw = true; }
        check("a window under four cell sizes is refused", threw);
        o.window_metres = 6 * k_cs;
        threw = false;
        try { (void)d8_flow_path_gradient(*flow, dem, make_raster(3, 3, std::vector<float>(9, 1.0f)), o); }
        catch (const std::invalid_argument&) { threw = true; }
        check("an area raster off the grid is refused", threw);
        std::vector<std::uint8_t> short_mask(5, 0);
        o.exclude = &short_mask;
        threw = false;
        try { (void)d8_flow_path_gradient(*flow, dem, area, o); } catch (const std::invalid_argument&) { threw = true; }
        check("an exclusion mask off the grid is refused", threw);
        o.exclude = nullptr;

        // ---- 7. threads ---------------------------------------------------------
        o.window_metres = 6 * k_cs;
        o.n_threads = 1;
        const auto g1 = d8_flow_path_gradient(*flow, dem, area, o);
        o.n_threads = 4;
        const auto g4 = d8_flow_path_gradient(*flow, dem, area, o);
        const bool same = g1.get_mask() == g4.get_mask() &&
                          std::memcmp(g1.get_data().data(), g4.get_data().data(), g1.get_data().size() * sizeof(float)) == 0;
        check("n_threads 1 and 4: bit-identical", same);
    }

    // ---- 4. exclusion -----------------------------------------------------------
    {
        const int nr = 9, nc = 30;
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(200.0 - S * k_cs * c);
        const auto plain = make_raster(nr, nc, z);
        // the routing (D8 receivers) of the plain plane, so the bump does not
        // re-route: it stands for cells a fill raised along the same path
        const auto flow = route(plain);
        const auto area = d8_area(*flow, plain);
        std::vector<float> bumped = z;
        std::vector<std::uint8_t> ex(z.size(), 0);
        for (int c = 14; c <= 16; ++c)
        {
            bumped[static_cast<std::size_t>(4 * nc + c)] += 1.0f;
            ex[static_cast<std::size_t>(4 * nc + c)] = 1;
        }
        const auto dem = make_raster(nr, nc, bumped);
        lsdtt_d8_flow_path_gradient_options o;
        o.window_metres = 10 * k_cs;
        o.n_threads = 1;
        const auto g_in = d8_flow_path_gradient(*flow, dem, area, o);
        check("bump kept in the fit: the gradient next to it is not S",
              valid(g_in, 4, 12) && std::abs(at(g_in, 4, 12) - S) > 1e-3, s(at(g_in, 4, 12)));
        o.exclude = &ex;
        const auto g_ex = d8_flow_path_gradient(*flow, dem, area, o);
        double worst = 0.0;
        bool ok = true;
        for (int c = 1; c < nc - 1; ++c)
        {
            if (c >= 14 && c <= 16) { ok = ok && !valid(g_ex, 4, c); continue; }
            if (!valid(g_ex, 4, c)) { ok = false; continue; }
            worst = std::max(worst, std::abs(at(g_ex, 4, c) - S));
        }
        check("bump excluded: gradient S elsewhere on the row, the bump cells nodata",
              ok && worst < 1e-6, "max |g - S| " + s(worst));
    }

    // ---- 2. south-east plane ------------------------------------------------------
    {
        const int nr = 30, nc = 30;
        // z falls by S per metre of distance along the diagonal (x + y) / sqrt 2;
        // rows count south, so y grows with r.
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[static_cast<std::size_t>(r * nc + c)] =
                    static_cast<float>(500.0 - S * k_cs * (c + r) / std::sqrt(2.0));
        const auto dem = make_raster(nr, nc, z);
        const auto flow = route(dem);
        const auto area = d8_area(*flow, dem);
        check("diagonal plane: the D8 path is the diagonal", flow->receiver_of(10 * nc + 10) == 11 * nc + 11);
        lsdtt_d8_flow_path_gradient_options o;
        o.window_metres = 6 * k_cs * std::sqrt(2.0);
        o.n_threads = 1;
        const auto g = d8_flow_path_gradient(*flow, dem, area, o);
        check("diagonal plane: gradient along the diagonal is S (sqrt 2 step lengths)",
              valid(g, 15, 15) && near(at(g, 15, 15), S, 1e-5), s(at(g, 15, 15)));
    }

    // ---- 5. donor ranking ------------------------------------------------------------
    {
        const int nr = 21, nc = 60;
        const double Sf = 0.01;
        auto floor_z = [&](int c) { return 100.0 - Sf * k_cs * c; };
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                double h = floor_z(c) + 1.0 * std::abs(r - 10);
                if (r < 10) h = std::min(h, floor_z(30) + 0.5 * (10 - r) + 1.0 * std::abs(c - 30));
                z[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(h);
            }
        const auto dem = make_raster(nr, nc, z);
        const auto flow = route(dem);
        const auto area = d8_area(*flow, dem);
        auto F = [&](int r, int c) { return static_cast<std::int64_t>(r) * nc + c; };
        lsdtt_channel_slope_options so;
        so.donor_rank = &area;
        const auto by_area = channel_slope_at(*flow, dem, F(10, 32), 6 * k_cs, so);
        const auto legacy = channel_slope_at(*flow, dem, F(10, 32), 6 * k_cs, std::uint64_t{1});
        check("rank by the D8 area raster = the contributing-pixel walk (floor gradient 0.01)",
              by_area.slope_valid && legacy.slope_valid && by_area.slope_m_per_m == legacy.slope_m_per_m &&
                  near(by_area.slope_m_per_m, Sf, 1e-6),
              s(by_area.slope_m_per_m) + " vs " + s(legacy.slope_m_per_m));
        // a rank raster that favours the tributary column
        std::vector<float> fav(area.get_data());
        for (int r = 0; r < 10; ++r) fav[static_cast<std::size_t>(F(r, 30))] = 1.0e9f;
        const auto rank_trib = make_raster(nr, nc, fav);
        so.donor_rank = &rank_trib;
        const auto by_trib = channel_slope_at(*flow, dem, F(10, 32), 6 * k_cs, so);
        check("rank raster favouring the tributary: the walk goes up it and the gradient is steeper",
              by_trib.slope_valid && by_trib.slope_m_per_m > Sf + 1e-3, s(by_trib.slope_m_per_m));
    }

    // ---- 8. filled_pixel_mask ----------------------------------------------------
    {
        const int nr = 15, nc = 15;
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(100.0 - 0.5 * c);
        for (int r = 6; r <= 8; ++r)
            for (int c = 6; c <= 8; ++c) z[static_cast<std::size_t>(r * nc + c)] -= 5.0f;   // depression
        const auto input = make_raster(nr, nc, z);
        const auto filled = fill_pits(input, 0.0001f, "oooo", 1);
        const auto m = conditioned_pixel_mask(input, filled);
        bool exact = true;
        int n = 0;
        for (std::size_t i = 0; i < m.size(); ++i)
        {
            const bool raised = filled.get_data()[i] > input.get_data()[i];
            exact = exact && (m[i] == (raised ? k_pixel_filled : k_pixel_untouched));
            n += m[i] == k_pixel_filled ? 1 : 0;
        }
        bool depression_marked = true;
        for (int r = 6; r <= 8; ++r)
            for (int c = 6; c <= 8; ++c) depression_marked = depression_marked && m[static_cast<std::size_t>(r * nc + c)];
        check("fill: the mask is exactly the raised cells (code 1), and covers the depression",
              exact && depression_marked && n >= 9, std::to_string(n) + " marked");
        breach_least_cost_options bo;
        bo.fill_remainder = false;
        const auto breached = breach_depressions_least_cost(input, bo);
        const auto mb = conditioned_pixel_mask(input, breached);
        bool codes_right = true;
        int n_lowered = 0;
        for (std::size_t i = 0; i < mb.size(); ++i)
        {
            const float a = input.get_data()[i], b = breached.get_data()[i];
            const std::uint8_t want = b < a ? k_pixel_carved : (b > a ? k_pixel_filled : k_pixel_untouched);
            codes_right = codes_right && mb[i] == want;
            n_lowered += b < a ? 1 : 0;
        }
        check("least-cost breach: lowered cells carved (2), raised cells filled (1), others 0",
              n_lowered > 0 && codes_right, std::to_string(n_lowered) + " lowered");
        const auto same = conditioned_pixel_mask(input, input);
        bool none = true;
        for (const auto v : same) none = none && v == 0;
        check("equal surfaces: nothing marked", none);
        bool threw = false;
        try { (void)conditioned_pixel_mask(input, make_raster(3, 3, std::vector<float>(9, 1.0f))); }
        catch (const std::invalid_argument&) { threw = true; }
        check("conditioned_pixel_mask: a raster off the grid is refused", threw);
    }

    // ---- 9. F1: one-sided windows on 30 m cells ------------------------------------
    {
        const int nr = 9, nc = 20;
        const double cs30 = 30.0;
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(400.0 - S * cs30 * c);
        auto meta = make_meta(nr, nc);
        meta.cell_size = cs30;
        std::vector<float> data = z;
        std::vector<std::uint8_t> mask(data.size(), 1);
        const lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);
        const auto flow = route(dem);
        const auto area = d8_area(*flow, dem);
        lsdtt_d8_flow_path_gradient_options o;
        o.window_metres = 150.0;
        o.n_threads = 1;
        const auto g = d8_flow_path_gradient(*flow, dem, area, o);
        check("30 m cells, 150 m window: the divide pixel (one-sided, 60 m of 75 m) is valid, gradient S",
              valid(g, 4, 0) && near(at(g, 4, 0), S, 1e-6), valid(g, 4, 0) ? s(at(g, 4, 0)) : "nodata");
        check("30 m cells, 150 m window: the outlet pixel (one-sided) is valid, gradient S",
              valid(g, 4, nc - 1) && near(at(g, 4, nc - 1), S, 1e-6));
        // the old rule (min span = window/2) would have rejected both:
        lsdtt_channel_slope_options so;
        so.donor_rank = &area;
        const auto old_rule = channel_slope_at(*flow, dem, 4 * nc, 75.0, so);
        check("the unrelaxed rule (span >= window/2) rejects the divide pixel on 30 m cells",
              !old_rule.slope_valid);
    }

    // ---- 10. fit_log_slope_area ------------------------------------------------------
    {
        std::vector<double> A, S5, S0, Sm;
        for (int i = 0; i < 10; ++i)
        {
            const double a = std::pow(10.0, 2.0 + 0.2 * i);
            A.push_back(a);
            S5.push_back(2.0 * std::pow(a, -0.5));
            S0.push_back(0.3);
            Sm.push_back(std::pow(a, 0.1));
        }
        const auto f5 = fit_log_slope_area(A, S5, 5);
        check("S = 2 A^-0.5: theta 0.5, r2 1",
              f5.valid && near(f5.theta, 0.5, 1e-9) && near(f5.r2, 1.0, 1e-9) && f5.n_pixels == 10 &&
                  near(f5.log10_area_range, 1.8, 1e-9),
              s(f5.theta));
        const auto f0 = fit_log_slope_area(A, S0, 5);
        check("constant S: theta 0, r2 NaN", f0.valid && f0.theta == 0.0 && std::isnan(f0.r2), s(f0.theta));
        const auto fm = fit_log_slope_area(A, Sm, 5);
        check("S = A^0.1: theta -0.1", fm.valid && near(fm.theta, -0.1, 1e-9), s(fm.theta));
        check("exact power law: theta_se 0; constant S: theta_se 0",
              near(f5.theta_se, 0.0, 1e-9) && f0.theta_se == 0.0, s(f5.theta_se));
        // log10 A = 0, 1, 2, 3 and log10 S = 0, -1, -1, -3: b = -0.9, SS_res = 0.7,
        // Sxx = 5, so theta 0.9 and theta_se = sqrt(0.7 / 2 / 5) = 0.264575...
        const auto fse = fit_log_slope_area(std::vector<double>{1, 10, 100, 1000},
                                            std::vector<double>{1, 0.1, 0.1, 0.001}, 3);
        // the same case straight through the shared regression utility
        const auto reg = compute_linear_regression(std::vector<double>{0, 1, 2, 3}, std::vector<double>{0, -1, -1, -3});
        check("compute_linear_regression: slope -0.9, slope_standard_error sqrt(0.07)",
              near(reg.slope, -0.9, 1e-12) && near(reg.slope_standard_error, std::sqrt(0.07), 1e-12),
              s(reg.slope_standard_error));
        const auto reg2 = compute_linear_regression(std::vector<double>{0, 1}, std::vector<double>{0, 1});
        check("compute_linear_regression: two points, slope_standard_error NaN", std::isnan(reg2.slope_standard_error));
        check("hand case: theta 0.9, theta_se sqrt(0.07)",
              fse.valid && near(fse.theta, 0.9, 1e-12) && near(fse.theta_se, std::sqrt(0.07), 1e-12),
              s(fse.theta) + ", " + s(fse.theta_se));
        const auto feq = fit_log_slope_area(std::vector<double>(6, 500.0), std::vector<double>{1, 2, 3, 4, 5, 6}, 5);
        check("equal areas: invalid", !feq.valid && feq.n_pixels == 6);
        const auto ffew = fit_log_slope_area(std::vector<double>{1, 10, 100, 1000}, std::vector<double>{1, 1, 1, 1}, 5);
        check("4 pairs with min_pixels 5: invalid", !ffew.valid && ffew.n_pixels == 4);
        std::vector<double> S5bad = S5;
        S5bad[3] = 0.0;
        S5bad[4] = -1.0;
        const auto fdrop = fit_log_slope_area(A, S5bad, 5);
        check("S <= 0 dropped: 8 pairs, theta still 0.5", fdrop.valid && fdrop.n_pixels == 8 && near(fdrop.theta, 0.5, 1e-9));
    }

    // ---- 11. flag rule, outlet cut, exclusion -----------------------------------------
    {
        const int nr = 5, nc = 30;
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c) z[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(300.0 - 1.0 * c);
        const auto dem = make_raster(nr, nc, z);
        const auto flow = route(dem);
        auto rasters_for = [&](double theta) {
            std::vector<float> a(z.size()), g(z.size());
            for (int r = 0; r < nr; ++r)
                for (int c = 0; c < nc; ++c)
                {
                    const double av = 1000.0 * std::pow(1.3, c);       // log A grows along the row
                    a[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(av);
                    g[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(0.5 * std::pow(av, -theta));
                }
            return std::make_pair(make_raster(nr, nc, a), make_raster(nr, nc, g));
        };
        lsdtt_low_concavity_options lo;
        lo.threshold_contributing_pixels = 1;
        lo.window_pixels = 10;
        lo.min_pixels = 5;
        lo.max_abs_concavity = 0.15;
        lo.n_threads = 1;
        auto flag_at = [&](const lsdtt_low_concavity_result& res, int r, int c) -> int {
            for (std::size_t k = 0; k < res.flat.size(); ++k)
                if (res.flat[k] == static_cast<std::int64_t>(r) * nc + c) return res.flag[k];
            return -9;
        };
        auto idx_of = [&](const lsdtt_low_concavity_result& res, int r, int c) -> std::size_t {
            for (std::size_t k = 0; k < res.flat.size(); ++k)
                if (res.flat[k] == static_cast<std::int64_t>(r) * nc + c) return k;
            return res.flat.size();
        };
        for (const double th : {0.5, 0.1, -0.1, -0.3})
        {
            const auto rs = rasters_for(th);
            const auto res = find_low_concavity(*flow, rs.first, rs.second, lo);
            const std::size_t k = idx_of(res, 2, 5);
            const int want = std::abs(th) < 0.15 ? 1 : 0;
            check("theta " + s(th) + ": fitted exactly, flag " + std::to_string(want),
                  k < res.flat.size() && near(res.theta[k], th, 1e-5) && res.flag[k] == want,
                  k < res.flat.size() ? s(res.theta[k]) : "missing");
        }
        const auto rs = rasters_for(0.1);
        const auto res = find_low_concavity(*flow, rs.first, rs.second, lo);
        check("every pixel of the plane is a channel pixel (threshold 1)", res.flat.size() == static_cast<std::size_t>(nr * nc));
        check("window cut at the outlet: 3 pixels left, no fit (-1)", flag_at(res, 2, nc - 3) == -1);
        check("window cut at the outlet: 5 pixels left, a fit", flag_at(res, 2, nc - 5) == 1);
        std::vector<std::uint8_t> ex(z.size(), 0);
        ex[static_cast<std::size_t>(2 * nc + 5)] = 2;               // carved
        ex[static_cast<std::size_t>(2 * nc + 8)] = 1;               // filled
        lo.exclude = &ex;
        const auto rx = find_low_concavity(*flow, rs.first, rs.second, lo);
        check("an excluded (carved) pixel gets no fit", flag_at(rx, 2, 5) == -1);
        const std::size_t k4 = idx_of(rx, 2, 4);
        check("excluded window pixels are left out: 8 of 10 pairs, theta still exact",
              k4 < rx.flat.size() && rx.n_pixels[k4] == 8 && near(rx.theta[k4], 0.1, 1e-5),
              k4 < rx.flat.size() ? std::to_string(rx.n_pixels[k4]) : "missing");
        lo.exclude = nullptr;
        lo.n_threads = 4;
        const auto r4 = find_low_concavity(*flow, rs.first, rs.second, lo);
        lo.n_threads = 1;
        bool same = r4.flat == res.flat && r4.flag == res.flag && r4.n_pixels == res.n_pixels;
        for (std::size_t k = 0; same && k < res.theta.size(); ++k)
            same = (std::isnan(res.theta[k]) && std::isnan(r4.theta[k])) || res.theta[k] == r4.theta[k];
        check("find_low_concavity: n_threads 1 and 4 identical", same);

        // ---- 14. refusals ----------------------------------------------------------------
        auto refused = [&](lsdtt_low_concavity_options o, const lsdtt_raster_cpu& a) {
            try { (void)find_low_concavity(*flow, a, rs.second, o); } catch (const std::invalid_argument&) { return true; }
            return false;
        };
        lsdtt_low_concavity_options bad = lo;
        bad.window_pixels = 1;
        check("window_pixels 1 refused", refused(bad, rs.first));
        bad = lo;
        bad.min_pixels = 2;
        check("min_pixels 2 refused", refused(bad, rs.first));
        bad = lo;
        bad.min_pixels = 11;
        check("min_pixels above the window refused", refused(bad, rs.first));
        bad = lo;
        bad.max_abs_concavity = -0.1;
        check("negative max_abs_concavity refused", refused(bad, rs.first));
        check("an area raster off the grid refused", refused(lo, make_raster(3, 3, std::vector<float>(9, 1.0f))));
    }

    // ---- 12. junctions --------------------------------------------------------------------
    {
        const int nr = 21, nc = 60;
        auto floor_z = [&](int c) { return 100.0 - 0.1 * c; };
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                double h = floor_z(c) + 1.0 * std::abs(r - 10);
                if (r < 10) h = std::min(h, floor_z(30) + 0.5 * (10 - r) + 1.0 * std::abs(c - 30));
                z[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(h);
            }
        const auto dem = make_raster(nr, nc, z);
        const auto flow = route(dem);
        const auto area = d8_area(*flow, dem);
        auto F = [&](int r, int c) { return static_cast<std::int64_t>(r) * nc + c; };
        lsdtt_d8_flow_path_gradient_options go;
        go.window_metres = 4 * k_cs;
        go.n_threads = 1;
        const auto g = d8_flow_path_gradient(*flow, dem, area, go);
        lsdtt_low_concavity_options lo;
        lo.threshold_contributing_pixels = flow->s_vector_size(F(9, 30));   // the tributary's mouth is a channel
        lo.n_threads = 1;
        const auto res = find_low_concavity(*flow, area, g, lo);
        std::size_t k = res.flat.size(), k2 = res.flat.size();
        for (std::size_t i = 0; i < res.flat.size(); ++i)
        {
            if (res.flat[i] == F(10, 25)) k = i;
            if (res.flat[i] == F(10, 40)) k2 = i;
        }
        check("junctions: the tributary mouth is a channel pixel, a side cell is not",
              flow->s_vector_size(F(9, 20)) < lo.threshold_contributing_pixels);
        check("junctions: the window from (10, 25) passes the confluence at (10, 30): 1 junction",
              k < res.flat.size() && res.n_junctions[k] == 1,
              k < res.flat.size() ? std::to_string(res.n_junctions[k]) : "missing");
        check("junctions: the window from (10, 40) has none", k2 < res.flat.size() && res.n_junctions[k2] == 0);
    }

    // ---- 13. integration: S = k A^-theta on a valley floor ---------------------------------
    for (const double th : {0.5, 0.05})
    {
        const int nr = 21, nc = 80;
        // floor along row 10; sides rise 5 m per row, so every side cell drains
        // straight to the floor and the floor area grows by 21 cells per column
        std::vector<double> floor(static_cast<std::size_t>(nc), 0.0);
        floor[static_cast<std::size_t>(nc - 1)] = 10.0;
        for (int c = nc - 2; c >= 0; --c)
        {
            const double a = (c + 1) * 21.0 * k_cs * k_cs;            // area at the floor cell c
            floor[static_cast<std::size_t>(c)] = floor[static_cast<std::size_t>(c + 1)] + 0.2 * std::pow(a / 1.0e4, -th) * k_cs;
        }
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(floor[static_cast<std::size_t>(c)] + 5.0 * std::abs(r - 10));
        const auto dem = make_raster(nr, nc, z);
        const auto flow = route(dem);
        const auto area = d8_area(*flow, dem);
        lsdtt_d8_flow_path_gradient_options go;
        go.window_metres = 4 * k_cs;
        go.n_threads = 1;
        const auto g = d8_flow_path_gradient(*flow, dem, area, go);
        lsdtt_low_concavity_options lo;
        lo.threshold_contributing_pixels = 200;                         // the floor from column 9 on
        lo.window_pixels = 10;
        lo.n_threads = 1;
        const auto res = find_low_concavity(*flow, area, g, lo);
        double worst = 0.0;
        int n_mid = 0;
        for (std::size_t k = 0; k < res.flat.size(); ++k)
        {
            const int r = static_cast<int>(res.flat[k] / nc), c = static_cast<int>(res.flat[k] % nc);
            if (r != 10 || c < 20 || c > 60) continue;
            ++n_mid;
            worst = std::max(worst, std::abs(res.theta[k] - th));
        }
        std::cout << "  (valley floor, theta " << th << ": max |error| " << worst << ")\n";
        check("valley floor built with theta " + s(th) + ": fitted theta within 0.05 on columns 20-60",
              n_mid == 41 && worst < 0.05, "max |theta - " + s(th) + "| = " + s(worst) + ", n " + std::to_string(n_mid));
    }

    // ---- 15. segments --------------------------------------------------------------------
    {
        const int nr = 21, nc = 60;
        auto floor_z = [&](int c) { return 100.0 - 0.1 * c; };
        std::vector<float> z(static_cast<std::size_t>(nr * nc));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                double h = floor_z(c) + 1.0 * std::abs(r - 10);
                if (r < 10) h = std::min(h, floor_z(30) + 0.5 * (10 - r) + 1.0 * std::abs(c - 30));
                z[static_cast<std::size_t>(r * nc + c)] = static_cast<float>(h);
            }
        const auto dem = make_raster(nr, nc, z);
        const auto flow = route(dem);
        const auto area = d8_area(*flow, dem);
        auto F = [&](int r, int c) { return static_cast<std::int64_t>(r) * nc + c; };
        lsdtt_d8_flow_path_gradient_options go;
        go.window_metres = 4 * k_cs;
        go.n_threads = 1;
        const auto g = d8_flow_path_gradient(*flow, dem, area, go);
        lsdtt_low_concavity_options lo;
        lo.threshold_contributing_pixels = flow->s_vector_size(F(9, 30));
        lo.n_threads = 1;
        auto res = find_low_concavity(*flow, area, g, lo);
        // flags by hand: floor columns 40-44 flagged (a run of 5), column 50
        // flagged alone, everything else 0
        for (std::size_t k = 0; k < res.flat.size(); ++k)
        {
            const int r = static_cast<int>(res.flat[k] / nc), c = static_cast<int>(res.flat[k] % nc);
            res.flag[k] = (r == 10 && ((c >= 40 && c <= 44) || c == 50)) ? 1 : 0;
        }
        const auto segs = build_low_concavity_segments(*flow, res, 1);
        std::size_t total = 0;
        std::vector<int> seen(res.flat.size(), 0);
        std::int32_t max_link = -1;
        for (const auto& sg : segs)
        {
            total += sg.flats.size();
            max_link = std::max(max_link, sg.link_id);
            for (const auto f : sg.flats)
                for (std::size_t k = 0; k < res.flat.size(); ++k)
                    if (res.flat[k] == f) ++seen[k];
        }
        bool once = true;
        for (const int v : seen) once = once && v == 1;
        check("segments partition the channel pixels (each exactly once)", total == res.flat.size() && once,
              std::to_string(total) + " vs " + std::to_string(res.flat.size()));
        check("three links: floor above the confluence, tributary, floor below", max_link == 2,
              std::to_string(max_link + 1));
        auto seg_of = [&](std::int64_t f) -> const lsdtt_low_concavity_segment* {
            for (const auto& sg : segs)
                for (const auto q : sg.flats)
                    if (q == f) return &sg;
            return nullptr;
        };
        const auto* run = seg_of(F(10, 42));
        check("the flagged run 40-44 is one segment of 5 pixels, flag 1",
              run && run->flag == 1 && run->flats.size() == 5 && run->flats.front() == F(10, 40),
              run ? std::to_string(run->flats.size()) : "missing");
        const auto* lone = seg_of(F(10, 50));
        check("the lone flagged pixel is its own flag-1 segment", lone && lone->flag == 1 && lone->flats.size() == 1);
        check("a one-pixel segment's length is the step to the next pixel down",
              lone && near(lone->length_m, k_cs, 1e-9) && lone->end_flat == F(10, 51));
        const auto* above = seg_of(F(10, 20));
        const auto* trib = seg_of(F(9, 30));
        const auto* below = seg_of(F(10, 35));
        check("Strahler: floor above and tributary 1, floor below the confluence 2",
              above && trib && below && above->stream_order == 1 && trib->stream_order == 1 &&
                  below->stream_order == 2);
        check("the confluence pixel starts the link below it",
              seg_of(F(10, 30)) && seg_of(F(10, 30))->link_id == below->link_id &&
                  seg_of(F(10, 29))->link_id != below->link_id);
        const auto merged = build_low_concavity_segments(*flow, res, 3);
        const lsdtt_low_concavity_segment* lone3 = nullptr;
        for (const auto& sg : merged)
            for (const auto q : sg.flats)
                if (q == F(10, 50)) lone3 = &sg;
        check("min_run_pixels 3: the lone flagged pixel is relabelled 0 and merged",
              lone3 && lone3->flag == 0 && lone3->flats.size() > 1);
        std::vector<std::uint8_t> inc(z.size(), 1);
        for (int c = 0; c < 25; ++c) inc[static_cast<std::size_t>(F(10, c))] = 0;
        const auto clipped = build_low_concavity_segments(*flow, res, 1, &inc);
        bool none_outside = true;
        for (const auto& sg : clipped)
            for (const auto q : sg.flats) none_outside = none_outside && inc[static_cast<std::size_t>(q)];
        check("pixels outside include are left out of every segment", none_outside);
    }

    std::cout << "RESULT: " << (g_fail == 0 ? "PASS" : "FAIL") << " (" << g_n << " checks)\n";
    return g_fail == 0 ? 0 : 1;
}
