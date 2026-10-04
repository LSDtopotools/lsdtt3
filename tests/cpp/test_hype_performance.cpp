// test_hype_performance.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// performance_hype_model: timing harness for the HYPE model (WP7 B0,
// plan_potential_shype_port.md section 7.7.1).
//
// Usage: lsdtt-test-hype-performance <output_dir> [--quick] [--skip-large] [--repeats N] [--no-balance] [--no-unit-fluxes]
//   --quick       small workload only (the registered CTest test uses this)
//   --skip-large  small and medium only
//   --no-balance  timed repeats run with balance tracking off (WP7 B1); the
//                 phase split's "balance" share is then not part of step time
//   --repeats N   timed repeats per workload (default 5), after one warm-up
//
// Workloads (synthetic, in memory; tests/cpp/hype_synthetic_setup.hpp):
//   small   100 sub-basins x 5 classes,   1 year
//   medium  2,000 x 8,                    5 years
//   large   10,000 x 10,                  10 years
//
// Per workload, each the median of the timed repeats:
//   setup_s           constructing lsdtt_hype_model_cpu from the set-up
//   step_total_s      every model.step() of the run
//   per_step_ms, per_unit_step_ns
//   peak_rss_mb       VmHWM of /proc/self/status after the model runs (the
//                     peak is reset per workload through /proc/self/clear_refs
//                     where the kernel allows it); -1 where not available
//
// PHASE SPLIT, MEASURED FROM OUTSIDE THE MODEL. The model has no internal
// timers and B0 does not change model code, so the split is reconstructed:
//   land     a copy of the freshly constructed model's units is stepped with
//            hype_unit_step over all units, every step, with the same forcing
//   route    a hype_river_network built the same way as the model's
//            (make_hype_river from the same lengths, the same downstream
//            indices and order) is routed with hype_route_network_step on the
//            runoff that the land copy produces, every step
//   balance  model.storage_m3() called n_steps times on the stepped model
//            (its cost does not depend on the values held)
//   agg      the remainder: step_total - land - route - balance. It holds
//            the per-sub-basin aggregation, the outlet-flow sum and anything
//            not captured by the three replicas (call overhead, cache effects
//            of running the phases back to back inside step()).
// Each share is that phase's median time over the median step_total. The
// replica's own aggregation is also timed directly and printed as a check.
//
// PASS/FAIL: the domain water balance closes (|closure error| <= 1e-12 x
// cumulative precipitation) and every flow and state is finite. Speed is
// reported, never asserted.
//
// Writes hype_performance.csv to output_dir.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../../src/cpu/lsdtt_hype_column_cpu.hpp"
#include "../../src/cpu/lsdtt_hype_model_cpu.hpp"
#include "../../src/cpu/lsdtt_hype_parameters.hpp"
#include "../../src/cpu/lsdtt_hype_routing_cpu.hpp"
#include "hype_synthetic_setup.hpp"

namespace
{

using clock_type = std::chrono::steady_clock;

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

double seconds_since(clock_type::time_point t0)
{
    return std::chrono::duration<double>(clock_type::now() - t0).count();
}

double median(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return (n % 2 == 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Peak resident set [MB] from /proc/self/status (VmHWM); -1 off Linux.
double peak_rss_mb()
{
#ifdef __linux__
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line))
    {
        if (line.rfind("VmHWM:", 0) == 0)
        {
            const double kb = std::strtod(line.c_str() + 6, nullptr);
            return kb / 1024.0;
        }
    }
#endif
    return -1.0;
}

// Reset the VmHWM peak (Linux >= 4.0: write "5" to clear_refs). Returns
// whether it worked.
bool reset_peak_rss()
{
#ifdef __linux__
    std::ofstream f("/proc/self/clear_refs");
    if (!f) return false;
    f << "5";
    f.flush();
    return static_cast<bool>(f);
#else
    return false;
#endif
}

struct workload
{
    std::string name;
    std::size_t n_subbasins;
    int         n_classes;
    int         n_years;
};

struct result
{
    std::string name;
    std::size_t n_subbasins = 0, n_units = 0;
    int n_steps = 0;
    double setup_s = 0, step_total_s = 0, land_s = 0, route_s = 0, balance_s = 0, agg_direct_s = 0;
    double replica_total_s = 0;
    double peak_rss_mb = -1;
    bool   rss_reset = false;
    double max_rel_closure = 0;
};

bool all_finite(const std::vector<double>& v)
{
    for (double x : v)
        if (!std::isfinite(x)) return false;
    return true;
}

// Checks after a model run: balance closure and finite flows and states.
double check_model(const lsdtt::lsdtt_hype_model_cpu& model, const std::string& label)
{
    const auto& b = model.balance();
    const double rel = std::abs(b.closure_error_m3()) / b.precipitation_m3;
    check(b.precipitation_m3 > 0.0 && rel <= 1.0e-12,
          label + ": water balance closes (|closure| / precipitation = " + std::to_string(rel) + ")");
    check(all_finite(model.outflow_m3s()) && all_finite(model.local_outflow_m3s()) &&
              all_finite(model.runoff_mm()),
          label + ": final flows finite");
    bool states_ok = std::isfinite(b.storage_m3) && std::isfinite(b.evaporation_m3) &&
                     std::isfinite(b.outflow_m3);
    for (const auto& u : model.units())
    {
        states_ok = states_ok && std::isfinite(u.state.snow_mm) && std::isfinite(u.state.snowliq_mm) &&
                    std::isfinite(u.state.water_mm[0]) && std::isfinite(u.state.water_mm[1]) &&
                    std::isfinite(u.state.water_mm[2]);
    }
    check(states_ok, label + ": final states and balance terms finite");
    return rel;
}

// Process activity over a run, for the record (not asserted).
struct activity
{
    long snow = 0, aet = 0, macropore = 0, surface_div = 0, perc = 0, soil_runoff = 0, sat_excess = 0;
};

void count_activity(const lsdtt::lsdtt_hype_model_cpu& model, activity& a)
{
    for (std::size_t k = 0; k < model.units().size(); ++k)
    {
        const auto& u = model.units()[k];
        const auto& f = model.unit_fluxes()[k];
        if (u.state.snow_mm > 0.0) ++a.snow;
        if (f.soil_evap_mm[0] + f.soil_evap_mm[1] > 0.0) ++a.aet;
        if (f.macropore_mm > 0.0) ++a.macropore;
        if (f.excess_infiltration_mm > 0.0) ++a.surface_div;
        if (f.percolation_mm[0] + f.percolation_mm[1] > 0.0) ++a.perc;
        if (f.soil_runoff_mm[0] + f.soil_runoff_mm[1] + f.soil_runoff_mm[2] > 0.0) ++a.soil_runoff;
        if (f.saturation_excess_mm > 0.0) ++a.sat_excess;
    }
}

// --no-balance: timed repeats build the model with balance tracking off (WP7
// B1). The untimed warm-up keeps it on, so closure is still checked once.
bool g_track_balance = true;
// --no-unit-fluxes: timed repeats do not store per-unit flux records (WP7 B4b).
bool g_store_unit_fluxes = true;

result run_workload(const workload& w, int repeats)
{
    result r;
    r.name        = w.name;
    r.n_subbasins = w.n_subbasins;
    r.n_steps     = 365 * w.n_years;
    r.rss_reset   = reset_peak_rss();

    const lsdtt::hype_setup setup = hype_synthetic::make_hype_synthetic_setup(w.n_subbasins, w.n_classes);
    const hype_synthetic::hype_synthetic_forcing forcing = hype_synthetic::make_hype_synthetic_forcing(r.n_steps);
    const std::size_t n_steps = static_cast<std::size_t>(r.n_steps);

    std::cout << "\n[" << w.name << "] " << w.n_subbasins << " sub-basins x " << w.n_classes
              << " classes, " << r.n_steps << " daily steps, 1 warm-up + " << repeats << " timed repeats\n";

    // ---- The model itself: set-up and stepping -----------------------------
    std::vector<double> setup_t, step_t;
    for (int rep = 0; rep <= repeats; ++rep)
    {
        const bool warmup = (rep == 0);
        auto t0 = clock_type::now();
        lsdtt::lsdtt_hype_model_cpu model(setup, 86400.0, warmup ? true : g_track_balance);
        if (!warmup) model.set_store_unit_fluxes(g_store_unit_fluxes);
        const double ts = seconds_since(t0);
        r.n_units = model.units().size();

        activity act;
        bool flows_finite = true;
        t0 = clock_type::now();
        if (warmup)
        {
            // Untimed checks on every step of the warm-up run.
            for (std::size_t i = 0; i < n_steps; ++i)
            {
                model.step(forcing.prec_mm[i], forcing.temp_c[i], forcing.doy[i]);
                flows_finite = flows_finite && all_finite(model.outflow_m3s()) &&
                               all_finite(model.local_outflow_m3s()) && all_finite(model.runoff_mm());
                count_activity(model, act);
            }
        }
        else
        {
            for (std::size_t i = 0; i < n_steps; ++i)
                model.step(forcing.prec_mm[i], forcing.temp_c[i], forcing.doy[i]);
        }
        const double tt = seconds_since(t0);

        if (model.balance_tracked())
        {
            const double rel = check_model(model, w.name + (warmup ? " warm-up" : " repeat " + std::to_string(rep)));
            r.max_rel_closure = std::max(r.max_rel_closure, rel);
        }
        if (warmup)
        {
            check(flows_finite, w.name + " warm-up: flows finite on every step");
            const double us = static_cast<double>(r.n_units) * static_cast<double>(n_steps);
            std::cout << "  process activity (fraction of unit-steps): snow on ground "
                      << act.snow / us << ", AET " << act.aet / us << ", macropore "
                      << act.macropore / us << ", excess infiltration " << act.surface_div / us
                      << ", percolation " << act.perc / us << ", soil runoff " << act.soil_runoff / us
                      << ", saturation excess " << act.sat_excess / us << "\n";
        }
        else
        {
            setup_t.push_back(ts);
            step_t.push_back(tt);
        }
    }
    r.peak_rss_mb = peak_rss_mb();   // before the replicas below add their own copies
    r.setup_s      = median(setup_t);
    r.step_total_s = median(step_t);

    // ---- Phase replicas ----------------------------------------------------
    std::vector<double> land_t, agg_t, route_t, bal_t, replica_t;
    for (int rep = 0; rep <= repeats; ++rep)
    {
        lsdtt::lsdtt_hype_model_cpu model(setup);
        std::vector<lsdtt::hype_model_unit> units = model.units();   // initial state
        const std::vector<std::size_t>& offset = model.unit_offset();

        // The model's network, rebuilt the same way (lsdtt_hype_model_cpu ctor).
        lsdtt::hype_river_network net;
        for (const auto& sb : setup.subbasins)
        {
            net.local.push_back(lsdtt::make_hype_river(sb.loc_rivlen_m, setup.general.rivvel,
                                                       setup.general.damp, 86400.0));
            net.main.push_back(lsdtt::make_hype_river(sb.rivlen_m, setup.general.rivvel,
                                                      setup.general.damp, 86400.0));
            net.downstream.push_back(sb.maindown == 0 ? -1 : sb.maindown - 1);   // subid = index + 1
            net.area_m2.push_back(sb.area_m2);
        }
        net.order = lsdtt::hype_subbasin_downstream_order(setup.subbasins);
        lsdtt::hype_network_flows flows;
        std::vector<double> runoff(w.n_subbasins, 0.0);

        double t_land = 0, t_agg = 0, t_route = 0;
        double prec_m3 = 0, evap_m3 = 0;
        lsdtt::hype_unit_forcing f;
        std::vector<lsdtt::hype_unit_fluxes> unit_fluxes(units.size());
        const auto t_all = clock_type::now();
        for (std::size_t i = 0; i < n_steps; ++i)
        {
            f.precipitation_mm = forcing.prec_mm[i];
            f.temperature_c    = forcing.temp_c[i];
            f.day_of_year      = forcing.doy[i];
            f.pet_seasonal_factor = lsdtt::hype_pet_seasonal_factor(setup.general.cevpam, setup.general.cevpph, f.day_of_year);   // as the model does (B3)

            auto t0 = clock_type::now();
            for (std::size_t k = 0; k < units.size(); ++k)
                lsdtt::hype_unit_step(setup.general, units[k].constants, f, units[k].state, unit_fluxes[k]);
            auto t1 = clock_type::now();
            double step_prec = 0.0, step_evap = 0.0;
            for (std::size_t s = 0; s + 1 < offset.size(); ++s)
            {
                double ro = 0.0;
                for (std::size_t k = offset[s]; k < offset[s + 1]; ++k)
                {
                    const auto& u = units[k];
                    ro += unit_fluxes[k].runoff_mm * u.fraction;
                    const double area_m2 = u.fraction * net.area_m2[s];
                    step_prec += unit_fluxes[k].class_precipitation_mm * 1.0e-3 * area_m2;
                    step_evap += unit_fluxes[k].evap_mm * 1.0e-3 * area_m2;
                }
                runoff[s] = ro;
            }
            prec_m3 += step_prec;
            evap_m3 += step_evap;
            auto t2 = clock_type::now();
            lsdtt::hype_route_network_step(net, runoff, 86400.0, flows);
            auto t3 = clock_type::now();
            t_land  += std::chrono::duration<double>(t1 - t0).count();
            t_agg   += std::chrono::duration<double>(t2 - t1).count();
            t_route += std::chrono::duration<double>(t3 - t2).count();
        }
        const double t_replica = seconds_since(t_all);

        // The replica must be the model's run: same runoff and outflow at the end.
        for (std::size_t i = 0; i < n_steps; ++i)
            model.step(forcing.prec_mm[i], forcing.temp_c[i], forcing.doy[i]);
        if (rep == 0)
        {
            check(runoff == model.runoff_mm() && flows.outflow_m3s == model.outflow_m3s() &&
                      prec_m3 == model.balance().precipitation_m3 &&
                      evap_m3 == model.balance().evaporation_m3,
                  w.name + ": phase replica reproduces the model run bitwise");
        }

        // Balance share: storage_m3() once per step.
        volatile double sink = 0.0;
        const auto tb = clock_type::now();
        for (std::size_t i = 0; i < n_steps; ++i) sink = sink + model.storage_m3();
        const double t_bal = seconds_since(tb);
        (void)sink;

        if (rep > 0)
        {
            land_t.push_back(t_land);
            agg_t.push_back(t_agg);
            route_t.push_back(t_route);
            bal_t.push_back(t_bal);
            replica_t.push_back(t_replica);
        }
    }
    r.land_s          = median(land_t);
    r.agg_direct_s    = median(agg_t);
    r.route_s         = median(route_t);
    r.balance_s       = median(bal_t);
    r.replica_total_s = median(replica_t);
    return r;
}

} // namespace

int main(int argc, char** argv)
{
    std::string out_dir = ".";
    bool quick = false, skip_large = false;
    int repeats = 5;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--quick") quick = true;
        else if (a == "--skip-large") skip_large = true;
        else if (a == "--repeats" && i + 1 < argc) repeats = std::max(1, std::atoi(argv[++i]));
        else if (a == "--no-balance") g_track_balance = false;
        else if (a == "--no-unit-fluxes") g_store_unit_fluxes = false;
        else out_dir = a;
    }

    const char* omp_env = std::getenv("OMP_NUM_THREADS");
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
    const std::string openmp_note = "OpenMP compiled in";
#else
    const int threads = 1;
    const std::string openmp_note = "OpenMP not compiled in";
#endif

    std::cout << "performance_hype_model (plan_potential_shype_port.md 7.7.1)\n"
              << "  OMP_NUM_THREADS = " << (omp_env ? omp_env : "(unset)") << "; " << openmp_note
              << ", omp_get_max_threads() = " << threads
              << ". The model is serial at B0, so this is recorded, not used.\n"
              << "  Phase split is reconstructed outside the model: land = hype_unit_step over a\n"
              << "  copy of the units; route = hype_route_network_step on a network built the same\n"
              << "  way; balance = model.storage_m3() n_steps times; agg = step_total - the three.\n";

    std::vector<workload> loads = {{"small", 100, 5, 1}};
    if (!quick)
    {
        loads.push_back({"medium", 2000, 8, 5});
        if (!skip_large) loads.push_back({"large", 10000, 10, 10});
    }

    std::vector<result> results;
    for (const auto& w : loads) results.push_back(run_workload(w, repeats));

    std::error_code mkdir_error;
    std::filesystem::create_directories(out_dir, mkdir_error);   // the CSV needs its directory
    const std::string csv_path = out_dir + "/hype_performance.csv";
    std::ofstream csv(csv_path);
    csv << "workload,n_subbasins,n_units,n_steps,threads,setup_s,step_total_s,per_step_ms,"
           "per_unit_step_ns,land_pct,agg_pct,route_pct,balance_pct,peak_rss_mb\n";

    std::cout << "\nSummary (medians):\n";
    for (const auto& r : results)
    {
        const double unit_steps = static_cast<double>(r.n_units) * r.n_steps;
        const double land_pct = 100.0 * r.land_s / r.step_total_s;
        const double route_pct = 100.0 * r.route_s / r.step_total_s;
        const double bal_pct = 100.0 * r.balance_s / r.step_total_s;
        const double agg_pct = 100.0 - land_pct - route_pct - bal_pct;
        const double per_step_ms = 1.0e3 * r.step_total_s / r.n_steps;
        const double per_unit_ns = 1.0e9 * r.step_total_s / unit_steps;
        csv << std::setprecision(6) << r.name << ',' << r.n_subbasins << ',' << r.n_units << ','
            << r.n_steps << ',' << threads << ',' << r.setup_s << ',' << r.step_total_s << ','
            << per_step_ms << ',' << per_unit_ns << ',' << land_pct << ',' << agg_pct << ','
            << route_pct << ',' << bal_pct << ',' << r.peak_rss_mb << '\n';

        std::cout << std::fixed << std::setprecision(3)
                  << "  " << r.name << ": " << r.n_subbasins << " sub-basins, " << r.n_units
                  << " units, " << r.n_steps << " steps\n"
                  << "    set-up " << 1.0e3 * r.setup_s << " ms; stepping " << r.step_total_s << " s; "
                  << per_step_ms << " ms/step; " << std::setprecision(1) << per_unit_ns
                  << " ns/unit-step\n"
                  << "    split land/agg/route/balance = " << land_pct << " / " << agg_pct << " / "
                  << route_pct << " / " << bal_pct << " %\n"
                  << "    check: replica aggregation timed directly = "
                  << 100.0 * r.agg_direct_s / r.step_total_s << " %; replica land+agg+route = "
                  << 100.0 * r.replica_total_s / r.step_total_s << " % of model stepping\n"
                  << "    peak RSS " << r.peak_rss_mb << " MB"
                  << (r.rss_reset ? " (peak reset for this workload)" : " (cumulative peak; reset unavailable)")
                  << "; max |closure|/precipitation " << std::scientific << std::setprecision(2)
                  << r.max_rel_closure << std::defaultfloat << "\n";
    }
    std::cout << "\nWrote " << csv_path << "\n"
              << g_checks << " checks, " << g_fails << " failures\n";
    return g_fails == 0 ? 0 : 1;
}
