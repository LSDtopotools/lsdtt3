// test_chi_step_changes.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Step-change detection on chi network branches (lsdtt_chi_step_changes_cpu;
// docs/plan_docs/plan_chi_knickpoints.md, Part B).
//
// Path cases (gradient threshold 0.40, minimum drop 5 m, pairs given as
// gradients on a background of 0.01):
//  1. one 0.50 pair at 30 m: one reach of that pair (15 m drop) with a 30,
//     90 and 300 m window;
//  2. pairs 0.45, 0.35, 0.45 at 30 m: a 30 m window gives two one-pair
//     reaches of 13.5 m; a 90 m window (mean 0.417) gives one reach of
//     37.5 m, the middle pair bridged;
//  3. pairs 0.90, 0.30, 0.05 at 30 m, 90 m window (mean 0.42): all three
//     marked, the trim leaves the 0.90 pair (27 m);
//  4. cases 2 and 3 at 10 m pixels (each pair three pixels) with a 90 m
//     window: the same reaches in metres (9 and 3 pairs);
//  5. a 0.45 pair at 10 m (4.5 m drop, under 5 m): no reach;
//  6. pairs 0.80, -0.05, 0.80 at 30 m, 90 m window: one reach of three
//     pairs, the uphill pair bridged and flagged, drop 46.5 m;
//  7. bad options throw.
// Network case (3 x 40 grid, 30 m, fitted by fit_chi_network): a main stem
// along row 1 (outlet at column 0) with a 20 m step on its first pair and a
// 15 m step at column 30, and a tributary along row 0 from column 20 whose
// first pair is a diagonal 30 m step onto the main stem at column 19.
// Expected: three reaches; the outlet one flagged at_basin_outlet, the
// tributary one at_tributary_mouth with length 30 sqrt 2 and its bottom at
// the junction; the column 30 one unflagged, 870 m above the base.
//
// Usage: lsdtt-test-chi-step-changes
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_chi_network_analysis_cpu.hpp"
#include "../../src/cpu/lsdtt_chi_step_changes_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

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
bool near(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol; }

using namespace lsdtt;

struct path
{
    std::vector<double> dz, len;
};

// 10 background pairs, the given gradients (each repeated `rep` times), 10
// background pairs; pixel length px.
path make_path(const std::vector<double>& grads, double px, int rep = 1)
{
    path p;
    auto add = [&](double g) { p.dz.push_back(g * px); p.len.push_back(px); };
    for (int i = 0; i < 10 * rep; ++i) add(0.01);
    for (const double g : grads)
        for (int r = 0; r < rep; ++r) add(g);
    for (int i = 0; i < 10 * rep; ++i) add(0.01);
    return p;
}

std::vector<lsdtt_chi_step_reach_on_path> run(const path& p, double window)
{
    lsdtt_chi_step_options o;
    o.max_window_length_metres = window;
    return detect_chi_step_reaches(p.dz, p.len, o);
}

std::string describe(const std::vector<lsdtt_chi_step_reach_on_path>& r)
{
    std::string d = s(static_cast<double>(r.size())) + " reaches:";
    for (const auto& x : r)
        d += " [" + s(static_cast<double>(x.first_pair)) + "-" + s(static_cast<double>(x.last_pair)) +
             " drop " + s(x.drop_metres) + "]";
    return d;
}
} // namespace

int main()
{
    // ---- 1. one steep pair ------------------------------------------------------
    {
        const path p = make_path({0.50}, 30.0);
        for (const double w : {30.0, 90.0, 300.0})
        {
            const auto r = run(p, w);
            check("1. one 0.50 pair, window " + s(w) + ": one reach of that pair, 15 m",
                  r.size() == 1 && r[0].first_pair == 10 && r[0].last_pair == 10 &&
                      near(r[0].drop_metres, 15.0) && near(r[0].max_window_gradient, 0.5),
                  describe(r));
        }
    }

    // ---- 2. 0.45, 0.35, 0.45 --------------------------------------------------------
    {
        const path p = make_path({0.45, 0.35, 0.45}, 30.0);
        const auto a = run(p, 30.0);
        check("2. 30 m window: two one-pair reaches of 13.5 m",
              a.size() == 2 && a[0].n_pairs == 1 && a[1].n_pairs == 1 && near(a[0].drop_metres, 13.5) &&
                  near(a[1].drop_metres, 13.5) && a[0].first_pair == 10 && a[1].first_pair == 12,
              describe(a));
        const auto b = run(p, 90.0);
        check("2. 90 m window: one reach of 37.5 m, middle pair bridged",
              b.size() == 1 && b[0].n_pairs == 3 && near(b[0].drop_metres, 37.5) &&
                  near(b[0].flow_length_metres, 90.0) && !b[0].contains_uphill_pairs,
              describe(b));
    }

    // ---- 3. 0.90, 0.30, 0.05 ---------------------------------------------------------
    {
        const auto r = run(make_path({0.90, 0.30, 0.05}, 30.0), 90.0);
        check("3. 90 m window: trimmed to the 0.90 pair, 27 m",
              r.size() == 1 && r[0].first_pair == 10 && r[0].last_pair == 10 && near(r[0].drop_metres, 27.0),
              describe(r));
    }

    // ---- 4. the same at 10 m pixels ----------------------------------------------------
    {
        const auto a = run(make_path({0.45, 0.35, 0.45}, 10.0, 3), 90.0);
        check("4. 10 m, 0.45/0.35/0.45, 90 m window: one reach of 9 pairs, 37.5 m",
              a.size() == 1 && a[0].n_pairs == 9 && near(a[0].drop_metres, 37.5, 1e-9) &&
                  near(a[0].flow_length_metres, 90.0),
              describe(a));
        const auto b = run(make_path({0.90, 0.30, 0.05}, 10.0, 3), 90.0);
        check("4. 10 m, 0.90/0.30/0.05, 90 m window: trimmed to the 3 0.90 pairs, 27 m",
              b.size() == 1 && b[0].n_pairs == 3 && b[0].first_pair == 30 && near(b[0].drop_metres, 27.0, 1e-9),
              describe(b));
    }

    // ---- 5. drop under the minimum ---------------------------------------------------
    {
        const auto r = run(make_path({0.45}, 10.0), 90.0);
        check("5. a 4.5 m step (minimum 5 m): no reach", r.empty(), describe(r));
    }

    // ---- 6. uphill pair inside a reach --------------------------------------------------
    {
        const auto r = run(make_path({0.80, -0.05, 0.80}, 30.0), 90.0);
        check("6. 0.80/-0.05/0.80: one reach of 3 pairs, uphill flagged, 46.5 m",
              r.size() == 1 && r[0].n_pairs == 3 && r[0].contains_uphill_pairs && near(r[0].drop_metres, 46.5),
              describe(r));
        const auto w = run(make_path({0.80, -0.05, 0.80}, 30.0), 30.0);
        check("6. same with a 30 m window: two one-pair reaches, no uphill flag",
              w.size() == 2 && !w[0].contains_uphill_pairs && !w[1].contains_uphill_pairs, describe(w));
    }

    // ---- 7. bad options -------------------------------------------------------------------
    {
        auto throws = [](lsdtt_chi_step_options o)
        {
            try { (void)detect_chi_step_reaches({1.0}, {30.0}, o); }
            catch (const std::invalid_argument&) { return true; }
            return false;
        };
        lsdtt_chi_step_options g; g.min_gradient = 0.0;
        lsdtt_chi_step_options w; w.max_window_length_metres = -1.0;
        lsdtt_chi_step_options d; d.min_drop_metres = -1.0;
        check("7. bad options throw", throws(g) && throws(w) && throws(d));
        bool size_throws = false;
        try { (void)detect_chi_step_reaches({1.0, 2.0}, {30.0}, {}); }
        catch (const std::invalid_argument&) { size_throws = true; }
        check("7. mismatched sizes throw", size_throws);
    }

    // ---- 8. network ---------------------------------------------------------------------------
    {
        const int nc = 40, nr = 3;
        const double cs = 30.0;
        lsdtt_raster_meta m;
        m.n_rows = nr;
        m.n_cols = nc;
        m.cell_size = cs;
        m.x_min = 500000.0;
        m.y_min = 5000000.0;
        m.set_crs_from_epsg(32632);
        std::vector<float> z(static_cast<std::size_t>(nr * nc), 0.0f);
        std::vector<std::uint8_t> valid(static_cast<std::size_t>(nr * nc), 0u);
        auto flat = [&](int row, int col) { return static_cast<std::int64_t>(row) * nc + col; };
        auto z_main = [&](int c) { return 2.0 * c + (c >= 1 ? 20.0 : 0.0) + (c >= 30 ? 15.0 : 0.0); };
        auto z_trib = [&](int c) { return z_main(19) + 30.0 + 2.0 + 2.0 * (c - 20); };

        lsdtt_chi_basin_result basin;
        basin.outlet_flat = flat(1, 0);
        basin.outlet_chi = 0.0f;
        const double dchi = 0.02;
        for (int c = 0; c < nc; ++c)
        {
            const std::int64_t f = flat(1, c);
            z[static_cast<std::size_t>(f)] = static_cast<float>(z_main(c));
            valid[static_cast<std::size_t>(f)] = 1u;
            basin.flats.push_back(f);
            basin.chi.push_back(static_cast<float>(dchi * c));
            basin.receiver_flats.push_back(c == 0 ? k_outlet_receiver : flat(1, c - 1));
        }
        for (int c = 20; c < 30; ++c)
        {
            const std::int64_t f = flat(0, c);
            z[static_cast<std::size_t>(f)] = static_cast<float>(z_trib(c));
            valid[static_cast<std::size_t>(f)] = 1u;
            basin.flats.push_back(f);
            basin.chi.push_back(static_cast<float>(dchi * 19 + dchi * std::sqrt(2.0) + dchi * (c - 20)));
            basin.receiver_flats.push_back(c == 20 ? flat(1, 19) : flat(0, c - 1));
        }
        lsdtt_raster_cpu dem(std::move(z), std::move(valid), m);
        lsdtt_chi_network_fit_collection_result fits;
        fits.basins.push_back(fit_chi_network(basin, dem, {}));
        lsdtt_chi_collection_result nets;
        nets.basins.push_back(basin);

        lsdtt_chi_step_options o;  // 0.40, 90 m, 5 m
        const auto rows = build_chi_step_changes(fits, nets, {}, m, o);
        check("8. network: three reaches", rows.size() == 3, s(static_cast<double>(rows.size())));
        const lsdtt_chi_step_row* outlet = nullptr;
        const lsdtt_chi_step_row* mid = nullptr;
        const lsdtt_chi_step_row* mouth = nullptr;
        for (const auto& r : rows)
        {
            if (r.bottom_flat == flat(1, 0)) outlet = &r;
            if (r.bottom_flat == flat(1, 29)) mid = &r;
            if (r.bottom_flat == flat(1, 19)) mouth = &r;
        }
        check("8. outlet step: bottom at the outlet, 22 m, at_basin_outlet",
              outlet && outlet->top_flat == flat(1, 1) && near(outlet->drop_m, 22.0, 1e-4) &&
                  outlet->at_basin_outlet == 1 && outlet->at_tributary_mouth == 0 &&
                  near(outlet->distance_to_branch_base_m, 0.0),
              outlet ? s(outlet->drop_m) : "missing");
        check("8. column 30 step: 17 m, 870 m above the base, no flags",
              mid && mid->top_flat == flat(1, 30) && near(mid->drop_m, 17.0, 1e-4) &&
                  near(mid->distance_to_branch_base_m, 870.0, 1e-6) && mid->at_basin_outlet == 0 &&
                  mid->at_tributary_mouth == 0,
              mid ? s(mid->distance_to_branch_base_m) : "missing");
        check("8. tributary mouth step: diagonal (30 sqrt 2), 32 m, at_tributary_mouth, bottom at the junction",
              mouth && mouth->top_flat == flat(0, 20) && near(mouth->flow_length_m, 30.0 * std::sqrt(2.0), 1e-9) &&
                  near(mouth->drop_m, 32.0, 1e-4) && mouth->at_tributary_mouth == 1 &&
                  mouth->at_basin_outlet == 0 && mouth->branch_id != 0 &&
                  near(mouth->bottom_chi, dchi * 19, 1e-6),
              mouth ? s(mouth->flow_length_m) + " branch " + s(mouth->branch_id) : "missing");
        check("8. main-stem steps are on branch 0", outlet && mid && outlet->branch_id == 0 && mid->branch_id == 0);
        check("8. coordinates are cell centres",
              outlet && near(outlet->bottom_x, m.x_min + 0.5 * cs) && near(outlet->bottom_y, m.y_min + 1.5 * cs));
        check("8. no keys: source_key -1, flow distance NaN",
              outlet && outlet->source_key == -1 && std::isnan(outlet->bottom_flow_distance_m));
    }

    std::cout << "\n" << (g_n - g_fail) << " / " << g_n << " checks passed\n";
    return g_fail == 0 ? 0 : 1;
}
