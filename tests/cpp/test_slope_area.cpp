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
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_area_cpu.hpp"

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

    std::cout << "RESULT: " << (g_fail == 0 ? "PASS" : "FAIL") << " (" << g_n << " checks)\n";
    return g_fail == 0 ? 0 : 1;
}
