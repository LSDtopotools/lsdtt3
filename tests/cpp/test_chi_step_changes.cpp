// test_chi_step_changes.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Step-change detection on chi network branches (lsdtt_chi_step_changes_cpu).
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

    // ---- 9. Part C: step-aware fit, one row ------------------------------------------
    // k_sn 50 (1 m per 0.02 chi per 10 m pixel), one 8 m step on pair 150
    // (nodes 149 -> 150, 9 m over 10 m), critical divergence 10 m: stage 1
    // fits one bent segment (the step's misfit is under 10 m); stage 2 gives
    // k_sn 50 either side and one step segment.
    {
        const int n = 301;
        const double dchi = 0.02;
        lsdtt_raster_meta m;
        m.n_rows = 1;
        m.n_cols = n;
        m.cell_size = 10.0;
        m.x_min = 500000.0;
        m.y_min = 5000000.0;
        m.set_crs_from_epsg(32632);
        std::vector<float> z(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) z[static_cast<std::size_t>(i)] = static_cast<float>(i * 1.0 + (i >= 150 ? 8.0 : 0.0));
        lsdtt_raster_cpu dem(std::move(z), std::vector<std::uint8_t>(static_cast<std::size_t>(n), 1u), m);
        lsdtt_chi_basin_result basin;
        basin.outlet_flat = 0;
        basin.outlet_chi = 0.0f;
        for (int i = 0; i < n; ++i)
        {
            basin.flats.push_back(i);
            basin.chi.push_back(static_cast<float>(dchi * i));
            basin.receiver_flats.push_back(i == 0 ? k_outlet_receiver : i - 1);
        }
        lsdtt_chi_network_fit_options o1;
        o1.critical_divergence_m = 10.0f;
        lsdtt_chi_network_fit_collection_result s1;
        s1.basins.push_back(fit_chi_network(basin, dem, o1));
        lsdtt_chi_collection_result nets;
        nets.basins.push_back(basin);
        const auto steps = build_chi_step_changes(s1, nets, {}, m, {});
        check("9. one step reach, pair 149 -> 150, 9 m",
              steps.size() == 1 && steps[0].bottom_flat == 149 && steps[0].top_flat == 150 &&
                  near(steps[0].drop_m, 9.0, 1e-4),
              s(static_cast<double>(steps.size())));
        check("9. stage 1: one bent segment, k_sn 50 + 8 / 6",
              s1.basins[0].segments.size() == 1 && near(s1.basins[0].segments[0].k_sn, 50.0 + 8.0 / 6.0, 1e-3) &&
                  !s1.basins[0].step_aware && s1.basins[0].is_step.empty(),
              s(static_cast<double>(s1.basins[0].segments.size())));

        lsdtt_chi_network_fit_options o2 = o1;
        o2.step_aware = true;
        o2.step_reaches = chi_step_reaches_for_fit(steps);
        const lsdtt_chi_network_fit_result r2 = fit_chi_network(basin, dem, o2);
        int n_step = 0, n_gentle_50 = 0;
        bool step_ok = false;
        for (const auto& g : r2.segments)
        {
            if (g.is_step)
            {
                ++n_step;
                step_ok = g.downstream_flat == 149 && g.upstream_flat == 150 && near(g.k_sn, 9.0 / 0.02, 1e-2);
            }
            else if (near(g.k_sn, 50.0, 1e-3))
                ++n_gentle_50;
        }
        check("9. stage 2: three segments, one step (149 -> 150), two of k_sn 50",
              r2.segments.size() == 3 && n_step == 1 && step_ok && n_gentle_50 == 2,
              s(static_cast<double>(r2.segments.size())) + " segments, " + s(n_step) + " step, " +
                  s(n_gentle_50) + " at 50");
        long n_step_nodes = 0;
        for (const auto v : r2.is_step) n_step_nodes += v;
        check("9. stage 2 is_step on nodes 149 and 150 only",
              r2.step_aware && r2.is_step.size() == static_cast<std::size_t>(n) && n_step_nodes == 2 &&
                  r2.is_step[149] == 1 && r2.is_step[150] == 1);
        double max_res = 0.0;
        for (const auto v : r2.abs_residual_elevation) max_res = std::max(max_res, static_cast<double>(v));
        check("9. stage 2 fits the profile exactly", max_res < 1e-3, s(max_res));

        // tables: is_step column only for the step-aware fit
        lsdtt_chi_network_fit_collection_result s2;
        s2.basins.push_back(r2);
        const auto nt1 = export_chi_network_fit_nodes(s1);
        const auto nt2 = export_chi_network_fit_nodes(s2);
        const auto st2 = export_chi_network_fit_segments(s2);
        check("9. node / segment tables carry is_step only when step-aware",
              !nt1.step_aware && nt1.is_step.empty() && nt2.step_aware && nt2.is_step.size() == nt2.flat.size() &&
                  st2.step_aware && st2.is_step.size() == st2.segment_id.size());

        // stage 2 knickpoints exclude the boundaries at the step
        const auto tk2 = build_chi_tributaries_and_knickpoints(s2, nets, {}, m);
        check("9. stage 2 knickpoints: none at the step's ends", tk2.knickpoints.empty(),
              s(static_cast<double>(tk2.knickpoints.size())));
        // stage 1 at critical divergence 2 m: knickpoints at the step, cross-referenced
        lsdtt_chi_network_fit_options o3;
        o3.critical_divergence_m = 2.0f;
        lsdtt_chi_network_fit_collection_result s3;
        s3.basins.push_back(fit_chi_network(basin, dem, o3));
        auto tk3 = build_chi_tributaries_and_knickpoints(s3, nets, {}, m);
        cross_reference_chi_knickpoints_with_steps(tk3.knickpoints, steps, s3);
        bool all_ref = !tk3.knickpoints.empty(), any_inside = false;
        for (const auto& k : tk3.knickpoints)
        {
            all_ref = all_ref && k.nearest_step_id == 0 && k.step_in_adjacent_segments == 1;
            any_inside = any_inside || near(k.nearest_step_chi_distance, 0.0);
        }
        check("9. stage 1 knickpoints at the step: nearest step 0, step in adjacent segments, one inside it",
              all_ref && any_inside, s(static_cast<double>(tk3.knickpoints.size())) + " rows");
    }

    // ---- 10. Part C: tributary anchored on a step segment -----------------------------
    // 2 x 300 grid, 10 m: main stem along row 1 (k_sn 50) with a 3-pair step
    // (+6 m on each of pairs 100, 101, 102: nodes 99 -> 102); a tributary
    // along row 0 from column 100 joins the main stem at column 100, inside
    // the step. Stage 2 must fit (consistent elevations) and anchor the
    // tributary on the step segment's fitted elevation there.
    {
        const int nc = 300, nr = 2;
        const double dchi = 0.02;
        lsdtt_raster_meta m;
        m.n_rows = nr;
        m.n_cols = nc;
        m.cell_size = 10.0;
        m.x_min = 500000.0;
        m.y_min = 5000000.0;
        m.set_crs_from_epsg(32632);
        auto flat = [&](int row, int col) { return static_cast<std::int64_t>(row) * nc + col; };
        auto z_main = [&](int c) { return 1.0 * c + 6.0 * ((c >= 100) + (c >= 101) + (c >= 102)); };
        std::vector<float> z(static_cast<std::size_t>(nr * nc), 0.0f);
        std::vector<std::uint8_t> valid(static_cast<std::size_t>(nr * nc), 0u);
        lsdtt_chi_basin_result basin;
        basin.outlet_flat = flat(1, 0);
        basin.outlet_chi = 0.0f;
        for (int c = 0; c < nc; ++c)
        {
            const std::int64_t f = flat(1, c);
            z[static_cast<std::size_t>(f)] = static_cast<float>(z_main(c));
            valid[static_cast<std::size_t>(f)] = 1u;
            basin.flats.push_back(f);
            basin.chi.push_back(static_cast<float>(dchi * c));
            basin.receiver_flats.push_back(c == 0 ? k_outlet_receiver : flat(1, c - 1));
        }
        for (int c = 100; c < 150; ++c)
        {
            const std::int64_t f = flat(0, c);
            z[static_cast<std::size_t>(f)] = static_cast<float>(z_main(100) + 1.0 + (c - 100));
            valid[static_cast<std::size_t>(f)] = 1u;
            basin.flats.push_back(f);
            basin.chi.push_back(static_cast<float>(dchi * 100 + dchi * (c - 99)));
            basin.receiver_flats.push_back(c == 100 ? flat(1, 100) : flat(0, c - 1));
        }
        lsdtt_raster_cpu dem(std::move(z), std::move(valid), m);
        lsdtt_chi_collection_result nets;
        nets.basins.push_back(basin);
        lsdtt_chi_network_fit_collection_result s1;
        s1.basins.push_back(fit_chi_network(basin, dem, {}));
        const auto steps = build_chi_step_changes(s1, nets, {}, m, {});
        const bool one_main_step = steps.size() == 1 && steps[0].branch_id == 0 &&
                                   steps[0].bottom_flat == flat(1, 99) && steps[0].top_flat == flat(1, 102);
        check("10. one step reach on the main stem, nodes 99 -> 102", one_main_step,
              steps.empty() ? "none" : s(static_cast<double>(steps[0].bottom_flat)) + " -> " +
                                           s(static_cast<double>(steps[0].top_flat)));

        lsdtt_chi_network_fit_options o2;
        o2.step_aware = true;
        o2.step_reaches = chi_step_reaches_for_fit(steps);
        bool threw = false;
        lsdtt_chi_network_fit_result r2;
        try { r2 = fit_chi_network(basin, dem, o2); }
        catch (const std::exception& e) { threw = true; std::cout << "  (" << e.what() << ")\n"; }
        check("10. stage 2 fits with a tributary anchored inside the step", !threw);
        if (!threw)
        {
            const lsdtt_chi_network_fit_segment* step = nullptr;
            const lsdtt_chi_network_fit_segment* trib_first = nullptr;
            for (const auto& g : r2.segments)
            {
                if (g.is_step) step = &g;
                if (g.branch_id != 0 && g.downstream_flat == flat(1, 100)) trib_first = &g;
            }
            const double z_line = step ? step->downstream_elevation +
                                             (step->upstream_elevation - step->downstream_elevation) / 3.0
                                       : 0.0;
            check("10. the tributary's first segment starts at the junction on the step line",
                  step && trib_first && step->downstream_flat == flat(1, 99) &&
                      step->upstream_flat == flat(1, 102) && near(trib_first->downstream_elevation, z_line, 1e-3),
                  trib_first ? s(trib_first->downstream_elevation) + " vs " + s(z_line) : "missing");
        }
    }

    // ---- 11. Part C: two steps sharing a bottom node -------------------------------------
    // As 8, plus a 20 m main-stem step on pair 19 -> 20: the junction (1, 19)
    // is the bottom of both that step and the tributary-mouth step. Both must
    // become step segments (one per branch).
    {
        const int nc = 40, nr = 3;
        const double cs = 30.0, dchi = 0.02;
        lsdtt_raster_meta m;
        m.n_rows = nr;
        m.n_cols = nc;
        m.cell_size = cs;
        m.x_min = 500000.0;
        m.y_min = 5000000.0;
        m.set_crs_from_epsg(32632);
        auto flat = [&](int row, int col) { return static_cast<std::int64_t>(row) * nc + col; };
        auto z_main = [&](int c) { return 2.0 * c + (c >= 20 ? 20.0 : 0.0); };
        std::vector<float> z(static_cast<std::size_t>(nr * nc), 0.0f);
        std::vector<std::uint8_t> valid(static_cast<std::size_t>(nr * nc), 0u);
        lsdtt_chi_basin_result basin;
        basin.outlet_flat = flat(1, 0);
        basin.outlet_chi = 0.0f;
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
            z[static_cast<std::size_t>(f)] = static_cast<float>(z_main(19) + 32.0 + 2.0 * (c - 20));
            valid[static_cast<std::size_t>(f)] = 1u;
            basin.flats.push_back(f);
            basin.chi.push_back(static_cast<float>(dchi * 19 + dchi * std::sqrt(2.0) + dchi * (c - 20)));
            basin.receiver_flats.push_back(c == 20 ? flat(1, 19) : flat(0, c - 1));
        }
        lsdtt_raster_cpu dem(std::move(z), std::move(valid), m);
        lsdtt_chi_collection_result nets;
        nets.basins.push_back(basin);
        lsdtt_chi_network_fit_collection_result s1;
        s1.basins.push_back(fit_chi_network(basin, dem, {}));
        const auto steps = build_chi_step_changes(s1, nets, {}, m, {});
        int shared = 0;
        for (const auto& st : steps) shared += st.bottom_flat == flat(1, 19);
        check("11. two steps with their bottom at the junction", shared == 2, s(shared));
        lsdtt_chi_network_fit_options o2;
        o2.step_aware = true;
        o2.step_reaches = chi_step_reaches_for_fit(steps);
        const auto r2 = fit_chi_network(basin, dem, o2);
        int n_step = 0, from_junction = 0;
        for (const auto& g : r2.segments)
        {
            n_step += g.is_step;
            if (g.is_step && g.downstream_flat == flat(1, 19)) ++from_junction;
        }
        check("11. stage 2: every step is a step segment, both junction ones",
              n_step == static_cast<int>(steps.size()) && from_junction == 2,
              s(n_step) + " of " + s(static_cast<double>(steps.size())) + ", " + s(from_junction) + " at the junction");
    }

    std::cout << "\n" << (g_n - g_fail) << " / " << g_n << " checks passed\n";
    return g_fail == 0 ? 0 : 1;
}
