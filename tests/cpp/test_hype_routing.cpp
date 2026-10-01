// test_hype_routing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// unit_hype_routing: HYPE river routing (WP4, Stage A).
//
// Usage: lsdtt-test-hype-routing <oracle_dir>   (tests/expected/hype)
//
// Analytic checks (no oracle):
//   * the damping box equals the exact solution of dS/dt = I - S/k over one
//     step, for inflow spread evenly over the step;
//   * the translation queue delays by the integer part and splits by the
//     fractional part, and conserves volume;
//   * a zero-length river passes flow straight on.
// Oracle checks: the `routing` case (cout, cloc, crun for two sub-basins in
// series) and `default_all` (cout, crun), with the HYPE column producing the
// runoff. Tolerance as in unit_hype_column: 5e-4 + 2e-5 x |value|.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_csv_table.hpp"
#include "../../src/core/lsdtt_distributed_hydrology_forcing.hpp"
#include "../../src/cpu/lsdtt_hype_column_cpu.hpp"
#include "../../src/cpu/lsdtt_hype_parameters.hpp"
#include "../../src/cpu/lsdtt_hype_routing_cpu.hpp"

namespace
{

int g_checks = 0;
int g_fails  = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok)
    {
        ++g_fails;
        std::cout << "  FAIL: " << what << "\n";
    }
}

constexpr double k_abs_tol = 5.0e-4;
constexpr double k_rel_tol = 2.0e-5;
constexpr double k_day     = 86400.0;

int day_of_year(const std::string& iso)
{
    const int y = std::stoi(iso.substr(0, 4));
    const int m = std::stoi(iso.substr(5, 2));
    const int d = std::stoi(iso.substr(8, 2));
    static const int cum[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    const bool leap = (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
    return cum[m - 1] + d + ((leap && m > 2) ? 1 : 0);
}

void analytic_checks()
{
    std::cout << "\n--- analytic checks ---\n";

    // Damping only (damp = 1 -> no translation): k = L / (v dt) timesteps.
    {
        const double L = 60000.0, v = 0.5;
        lsdtt::hype_river r = lsdtt::make_hype_river(L, v, 1.0, k_day);
        const double k = L / v / k_day;
        check(r.transtime == 0.0, "damp = 1: no translation");
        double S = 0.0, worst = 0.0;
        const double inflows[] = {2.0, 0.0, 5.0, 5.0, 0.0, 0.0, 1.0};
        for (double qin : inflows)
        {
            const double I = qin * k_day;                    // m3 over the step
            // Exact: dS/dt = I/1 - S/k over a unit step, I spread evenly.
            const double S1 = S * std::exp(-1.0 / k) + I * k * (1.0 - std::exp(-1.0 / k));
            const double Q_exact = (S + I - S1) / k_day;     // m3/s
            const double Q = lsdtt::hype_river_step(r, qin, k_day);
            worst = std::max(worst, std::fabs(Q - Q_exact));
            S = S1;
        }
        std::printf("  damping box vs exact linear reservoir: max |diff| %.3e m3/s\n", worst);
        check(worst < 1e-12, "damping box equals the exact linear-reservoir step");
    }

    // Translation only (damp = 0): transtime = L/(v dt). Pick 2.25 steps.
    {
        const double v = 1.0, L = 2.25 * v * k_day;
        lsdtt::hype_river r = lsdtt::make_hype_river(L, v, 0.0, k_day);
        check(r.ttstep == 2 && std::fabs(r.ttpart - 0.25) < 1e-12, "translation: 2.25 steps split into 2 + 0.25");
        std::vector<double> out;
        const double pulse[] = {4.0, 0, 0, 0, 0, 0};
        double vin = 0.0, vout = 0.0;
        for (double q : pulse)
        {
            out.push_back(lsdtt::hype_river_translate(r, q, k_day));
            vin += q;
            vout += out.back();
        }
        check(out[0] == 0.0 && out[1] == 0.0, "translation: nothing leaves before 2 steps");
        check(std::fabs(out[2] - 3.0) < 1e-12 && std::fabs(out[3] - 1.0) < 1e-12,
              "translation: pulse leaves 75% after 2 steps and 25% after 3");
        check(std::fabs(vin - vout) < 1e-12, "translation conserves volume");
    }

    // Zero-length river passes flow straight on.
    {
        lsdtt::hype_river r = lsdtt::make_hype_river(0.0, 1.0, 0.5, k_day);
        check(lsdtt::hype_river_step(r, 3.5, k_day) == 3.5, "zero-length river passes flow on");
    }
}

void run_case(const std::string& oracle_dir, const std::string& name,
              const std::vector<std::string>& basin_vars)
{
    std::cout << "\n--- case " << name << " ---\n";
    const std::string in = oracle_dir + "/" + name + "/lsdtt_input";
    const lsdtt::hype_setup setup = lsdtt::load_hype_setup_directory(in);
    const std::string fpath = in + "/forcing.csv";
    const auto prec_m = lsdtt::load_distributed_hydrology_depth_series(fpath, "precipitation_metres", "precipitation_mm");
    const auto temp_c = lsdtt::load_distributed_hydrology_temperature_series(fpath, "temperature_c");
    const auto ftab   = lsdtt::load_csv_table(fpath);
    const std::size_t n_steps = prec_m.size();

    const std::size_t nsub = setup.subbasins.size();
    std::map<int, std::size_t> index_of;
    for (std::size_t s = 0; s < nsub; ++s) index_of[setup.subbasins[s].subid] = s;

    // Units grouped by sub-basin.
    struct unit { std::size_t sub; double fraction; lsdtt::hype_unit_constants c; lsdtt::hype_unit_state st; };
    std::vector<unit> units;
    for (const auto& sc : setup.subbasin_classes)
    {
        if (!(sc.fraction > 0.0)) continue;
        lsdtt::hype_class_definition cd;
        for (const auto& c : setup.classes) if (c.class_id == sc.class_id) cd = c;
        lsdtt::hype_landuse_parameters lu;
        for (const auto& l : setup.landuse) if (l.code == cd.landuse) lu = l;
        lsdtt::hype_soil_parameters so;
        for (const auto& x : setup.soils) if (x.code == cd.soil) so = x;
        const std::size_t s = index_of.at(sc.subid);
        unit u{s, sc.fraction,
               lsdtt::make_hype_unit_constants(setup.general, lu, so, cd,
                                               setup.subbasins[s].slope_mean_percent, sc.dh_m), {}};
        u.st = lsdtt::initial_hype_unit_state(u.c);
        units.push_back(u);
    }

    // Rivers.
    lsdtt::hype_river_network net;
    for (const auto& sb : setup.subbasins)
    {
        net.local.push_back(lsdtt::make_hype_river(sb.loc_rivlen_m, setup.general.rivvel, setup.general.damp, k_day));
        net.main.push_back(lsdtt::make_hype_river(sb.rivlen_m, setup.general.rivvel, setup.general.damp, k_day));
        net.downstream.push_back(sb.maindown == 0 ? -1 : static_cast<int>(index_of.at(sb.maindown)));
        net.area_m2.push_back(sb.area_m2);
    }
    net.order = lsdtt::hype_subbasin_downstream_order(setup.subbasins);

    std::map<std::string, std::vector<std::vector<double>>> sim;   // var -> [sub][step]
    for (const auto& v : basin_vars) sim[v].assign(nsub, std::vector<double>(n_steps, 0.0));

    for (std::size_t t = 0; t < n_steps; ++t)
    {
        const lsdtt::hype_unit_forcing forcing{prec_m[t] * 1000.0, temp_c[t],
                                               day_of_year(lsdtt::csv_row_string(ftab.rows[t], "date"))};
        std::vector<double> runoff(nsub, 0.0);
        for (auto& u : units)
        {
            lsdtt::hype_unit_fluxes f;
            lsdtt::hype_unit_step(setup.general, u.c, forcing, u.st, f);
            runoff[u.sub] += f.runoff_mm * u.fraction;          // model_hype.f90:2361
        }
        lsdtt::hype_network_flows flows;
        lsdtt::hype_route_network_step(net, runoff, k_day, flows);
        for (std::size_t s = 0; s < nsub; ++s)
        {
            if (sim.count("cout")) sim["cout"][s][t] = flows.outflow_m3s[s];
            if (sim.count("cloc")) sim["cloc"][s][t] = flows.local_outflow_m3s[s];
            if (sim.count("crun")) sim["crun"][s][t] = runoff[s];   // land fraction 1
        }
    }

    for (const auto& v : basin_vars)
    {
        const std::string path = oracle_dir + "/" + name + "/expected/basin_" + v + ".csv";
        const auto ref = lsdtt::load_csv_table(path);
        check(ref.rows.size() == n_steps, name + "/" + v + ": oracle has one row per step");
        double max_abs = 0.0;
        int n_out = 0;
        std::string first_bad;
        for (std::size_t s = 0; s < nsub; ++s)
        {
            const std::string col = std::to_string(setup.subbasins[s].subid);
            for (std::size_t t = 0; t < n_steps && t < ref.rows.size(); ++t)
            {
                const double h = lsdtt::csv_row_double_required(ref.rows[t], col, path, static_cast<int>(t) + 1);
                const double l = sim[v][s][t];
                const double d = std::fabs(l - h);
                max_abs = std::max(max_abs, d);
                if (d > k_abs_tol + k_rel_tol * std::fabs(h))
                {
                    if (n_out == 0)
                    {
                        char buf[200];
                        std::snprintf(buf, sizeof(buf), "sub-basin %s step %zu: lsdtt %.9g vs HYPE %.9g",
                                      col.c_str(), t + 1, l, h);
                        first_bad = buf;
                    }
                    ++n_out;
                }
            }
        }
        std::printf("  %-5s max |diff| %.3e  out of tolerance %d\n", v.c_str(), max_abs, n_out);
        check(n_out == 0, name + "/" + v + ": matches HYPE within tolerance" +
                          (first_bad.empty() ? "" : " (first: " + first_bad + ")"));
    }
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-hype-routing <oracle_dir>\n";
        return 2;
    }
    try
    {
        analytic_checks();
        run_case(argv[1], "routing", {"cout", "cloc", "crun"});
        run_case(argv[1], "default_all", {"cout", "crun"});
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL: exception: " << e.what() << "\n";
        ++g_fails;
    }
    std::cout << "\n" << (g_checks - g_fails) << " / " << g_checks << " checks passed\n";
    std::cout << (g_fails == 0 ? "RESULT: PASS" : "RESULT: FAIL") << "\n";
    return g_fails == 0 ? 0 : 1;
}
