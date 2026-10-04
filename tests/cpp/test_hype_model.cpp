// test_hype_model.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// unit_hype_model: the assembled HYPE model (WP5) against HYPE.
//
// Usage: lsdtt-test-hype-model <oracle_dir>   (tests/expected/hype)
//
//   * default_all and routing through lsdtt_hype_model_cpu: every per-class
//     variable and every per-sub-basin variable HYPE wrote, every step.
//     Tolerance as in unit_hype_column: 5e-4 + 2e-5 x |value|.
//   * Domain water balance (snow, soil, river boxes, water in transit in the
//     translation queues) closes on every step.
//   * Observed PET: when supplied it replaces the calculated value in every
//     unit (HYPE PET option 0); actual ET still needs temperature above ttmp
//     and never exceeds it.
//   * A non-daily timestep is rejected.
//   * Balance tracking off (WP7 B1): flows and unit states bitwise identical
//     to tracking on; storage terms NaN, cumulative fluxes still kept.
//   * A column failure inside the parallel land loop (WP7 B2) comes out of
//     step() as an exception, at 1 thread and at the maximum thread count.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_csv_table.hpp"
#include "../../src/core/lsdtt_distributed_hydrology_forcing.hpp"
#include "../../src/cpu/lsdtt_hype_model_cpu.hpp"
#include "../../src/cpu/lsdtt_hype_parameters.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

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

int day_of_year(const std::string& iso)
{
    const int y = std::stoi(iso.substr(0, 4));
    const int m = std::stoi(iso.substr(5, 2));
    const int d = std::stoi(iso.substr(8, 2));
    static const int cum[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    const bool leap = (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
    return cum[m - 1] + d + ((leap && m > 2) ? 1 : 0);
}

struct forcing_series
{
    std::vector<double> prec_mm, temp_c;
    std::vector<int> doy;
};

forcing_series load_forcing(const std::string& in)
{
    const std::string fpath = in + "/forcing.csv";
    forcing_series f;
    for (double p : lsdtt::load_distributed_hydrology_depth_series(fpath, "precipitation_metres", "precipitation_mm"))
        f.prec_mm.push_back(p * 1000.0);
    f.temp_c = lsdtt::load_distributed_hydrology_temperature_series(fpath, "temperature_c");
    for (const auto& row : lsdtt::load_csv_table(fpath).rows)
        f.doy.push_back(day_of_year(lsdtt::csv_row_string(row, "date")));
    return f;
}

double class_value(const std::string& v, const lsdtt::hype_model_unit& u, const lsdtt::hype_unit_fluxes& f)
{
    if (v == "snow") return u.state.snow_mm;
    if (v == "melt") return f.out_melt_mm;
    if (v == "evap") return f.evap_mm;
    if (v == "infi") return f.out_infi_mm;
    if (v == "ros1") return f.saturation_excess_mm;
    if (v == "ros2") return f.excess_infiltration_mm;
    if (v == "crun") return f.runoff_mm;
    if (v == "sml1") return f.out_sml_mm[0];
    if (v == "sml2") return f.out_sml_mm[1];
    if (v == "sml3") return f.out_sml_mm[2];
    throw std::runtime_error("unknown class variable " + v);
}

void compare(const std::string& label, const std::string& path, const std::string& col,
             const std::vector<double>& sim)
{
    const auto ref = lsdtt::load_csv_table(path);
    double max_abs = 0.0;
    int n_out = 0;
    for (std::size_t t = 0; t < sim.size() && t < ref.rows.size(); ++t)
    {
        const double h = lsdtt::csv_row_double_required(ref.rows[t], col, path, static_cast<int>(t) + 1);
        const double d = std::fabs(sim[t] - h);
        max_abs = std::max(max_abs, d);
        if (d > k_abs_tol + k_rel_tol * std::fabs(h)) ++n_out;
    }
    check(ref.rows.size() == sim.size() && n_out == 0,
          label + " column " + col + ": matches HYPE (" + std::to_string(n_out) + " out of tolerance)");
    std::printf("  %-28s max |diff| %.3e\n", (label + " " + col).c_str(), max_abs);
}

void run_case(const std::string& oracle, const std::string& name,
              const std::vector<std::string>& class_vars, const std::vector<std::string>& basin_vars)
{
    std::cout << "\n--- case " << name << " ---\n";
    const std::string in = oracle + "/" + name + "/lsdtt_input";
    lsdtt::lsdtt_hype_model_cpu model(lsdtt::load_hype_setup_directory(in));
    const forcing_series fs = load_forcing(in);
    const std::size_t n = fs.prec_mm.size();

    std::map<std::string, std::vector<std::vector<double>>> cls, bas;   // var -> [unit|sub][t]
    for (const auto& v : class_vars) cls[v].assign(model.units().size(), std::vector<double>(n));
    for (const auto& v : basin_vars) bas[v].assign(model.n_subbasins(), std::vector<double>(n));

    double worst_balance = 0.0, scale = 0.0;
    for (std::size_t t = 0; t < n; ++t)
    {
        model.step(fs.prec_mm[t], fs.temp_c[t], fs.doy[t]);
        for (const auto& v : class_vars)
            for (std::size_t k = 0; k < model.units().size(); ++k) cls[v][k][t] = class_value(v, model.units()[k], model.unit_fluxes()[k]);
        for (std::size_t s = 0; s < model.n_subbasins(); ++s)
        {
            if (bas.count("cout")) bas["cout"][s][t] = model.outflow_m3s()[s];
            if (bas.count("cloc")) bas["cloc"][s][t] = model.local_outflow_m3s()[s];
            if (bas.count("crun")) bas["crun"][s][t] = model.runoff_mm()[s];
        }
        worst_balance = std::max(worst_balance, std::fabs(model.balance().closure_error_m3()));
        scale = std::max(scale, model.balance().precipitation_m3);
    }
    std::printf("  domain balance: largest |closure error| %.3e m3 (cumulative precipitation %.3e m3)\n",
                worst_balance, scale);
    check(worst_balance <= 1.0e-12 * std::max(scale, 1.0),
          name + ": domain water balance closes to 1e-12 of cumulative precipitation");

    for (const auto& v : class_vars)
        for (std::size_t k = 0; k < model.units().size(); ++k)
        {
            const auto& u = model.units()[k];
            compare(name + "/" + v, oracle + "/" + name + "/expected/class_" + v + ".csv",
                    std::to_string(model.subids()[u.subbasin]) + "_" + std::to_string(u.class_id), cls[v][k]);
        }
    for (const auto& v : basin_vars)
        for (std::size_t s = 0; s < model.n_subbasins(); ++s)
            compare(name + "/" + v, oracle + "/" + name + "/expected/basin_" + v + ".csv",
                    std::to_string(model.subids()[s]), bas[v][s]);
}

void observed_pet_checks(const std::string& oracle)
{
    std::cout << "\n--- observed PET ---\n";
    const std::string in = oracle + "/default_all/lsdtt_input";
    const auto setup = lsdtt::load_hype_setup_directory(in);
    lsdtt::lsdtt_hype_model_cpu model(setup);
    const forcing_series fs = load_forcing(in);
    bool epot_replaced = true, aet_bounded = true, aet_off_when_cold = true;
    int warm_steps = 0, cold_steps = 0;
    for (std::size_t t = 0; t < fs.prec_mm.size(); ++t)
    {
        const double obs = 1.0 + 0.5 * static_cast<double>(t % 7);   // 1.0 .. 4.0 mm
        model.step(fs.prec_mm[t], fs.temp_c[t], fs.doy[t], obs);
        for (std::size_t k = 0; k < model.units().size(); ++k)
        {
            const auto& u = model.units()[k];
            const auto& f = model.unit_fluxes()[k];
            if (f.epot_mm != obs) epot_replaced = false;
            if (f.evap_mm > obs + 1e-12) aet_bounded = false;
            if (!(f.class_temperature_c > u.constants.landuse.ttmp))
            {
                ++cold_steps;
                if (f.evap_mm != 0.0) aet_off_when_cold = false;
            }
            else ++warm_steps;
        }
    }
    check(epot_replaced, "observed PET replaces the calculated PET in every unit on every step");
    check(aet_bounded, "actual ET never exceeds the observed PET");
    check(cold_steps > 0 && warm_steps > 0, "the case has unit-steps both above and below ttmp");
    check(aet_off_when_cold, "no actual ET when the class temperature is not above ttmp");
    check(std::fabs(model.balance().closure_error_m3()) <= 1e-12 * model.balance().precipitation_m3,
          "domain balance closes with observed PET");

    // NaN means "not supplied": identical to a model stepped without it.
    lsdtt::lsdtt_hype_model_cpu a(setup), b(setup);
    bool same = true;
    for (std::size_t t = 0; t < fs.prec_mm.size(); ++t)
    {
        a.step(fs.prec_mm[t], fs.temp_c[t], fs.doy[t]);
        b.step(fs.prec_mm[t], fs.temp_c[t], fs.doy[t], std::numeric_limits<double>::quiet_NaN());
        for (std::size_t s = 0; s < a.n_subbasins(); ++s)
            if (a.outflow_m3s()[s] != b.outflow_m3s()[s]) same = false;
    }
    check(same, "NaN observed PET gives exactly the calculated-PET run");
}

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

bool same_bits(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

bool same_state(const lsdtt::hype_unit_state& a, const lsdtt::hype_unit_state& b)
{
    return same_bits(a.snow_mm, b.snow_mm) && same_bits(a.snowliq_mm, b.snowliq_mm) &&
           same_bits(a.water_mm[0], b.water_mm[0]) && same_bits(a.water_mm[1], b.water_mm[1]) &&
           same_bits(a.water_mm[2], b.water_mm[2]);
}

// WP7 B1: switching the balance off skips only the storage pass.
void balance_off_checks(const std::string& oracle)
{
    std::cout << "\n--- balance tracking off ---\n";
    const std::string in = oracle + "/default_all/lsdtt_input";
    const auto setup = lsdtt::load_hype_setup_directory(in);
    const forcing_series fs = load_forcing(in);
    lsdtt::lsdtt_hype_model_cpu on(setup);
    lsdtt::lsdtt_hype_model_cpu off(setup, 86400.0, false);
    check(on.balance_tracked() && !off.balance_tracked(), "balance_tracked() reports the constructor flag");

    bool flows_same = true, states_same = true, fluxes_same = true;
    for (std::size_t t = 0; t < fs.prec_mm.size(); ++t)
    {
        on.step(fs.prec_mm[t], fs.temp_c[t], fs.doy[t]);
        off.step(fs.prec_mm[t], fs.temp_c[t], fs.doy[t]);
        if (!same_bits(on.runoff_mm(), off.runoff_mm()) ||
            !same_bits(on.local_outflow_m3s(), off.local_outflow_m3s()) ||
            !same_bits(on.outflow_m3s(), off.outflow_m3s())) flows_same = false;
        for (std::size_t k = 0; k < on.units().size(); ++k)
        {
            const auto& a = on.units()[k];
            const auto& b = off.units()[k];
            if (!same_state(a.state, b.state)) states_same = false;
            const auto& fa = on.unit_fluxes()[k];
            const auto& fb = off.unit_fluxes()[k];
            if (!same_bits(fa.runoff_mm, fb.runoff_mm) ||
                !same_bits(fa.evap_mm, fb.evap_mm) ||
                !same_bits(fa.class_precipitation_mm, fb.class_precipitation_mm)) fluxes_same = false;
        }
    }
    check(flows_same, "tracking off: sub-basin runoff, local outflow and outflow bitwise identical, every step");
    check(states_same, "tracking off: every unit state bitwise identical, every step");
    check(fluxes_same, "tracking off: unit runoff, evaporation and precipitation bitwise identical, every step");

    const auto& bo = on.balance();
    const auto& bf = off.balance();
    check(bo.storage_tracked && !bf.storage_tracked, "storage_tracked follows the flag");
    check(same_bits(bo.precipitation_m3, bf.precipitation_m3) && same_bits(bo.evaporation_m3, bf.evaporation_m3) &&
          same_bits(bo.outflow_m3, bf.outflow_m3),
          "tracking off: cumulative precipitation, evaporation and outflow bitwise identical");
    check(std::isnan(bf.storage_start_m3) && std::isnan(bf.storage_m3) && std::isnan(bf.closure_error_m3()),
          "tracking off: storage terms and closure error are NaN");
    check(std::isfinite(bo.closure_error_m3()), "tracking on: closure error is finite");
    check(same_bits(on.storage_m3(), off.storage_m3()), "storage_m3() on request is the same with tracking off");
}


// WP7 B4b: not storing per-unit fluxes must not change any result.
void unit_flux_storage_checks(const std::string& oracle)
{
    std::cout << "\n--- per-unit flux storage off ---\n";
    const std::string in = oracle + "/default_all/lsdtt_input";
    const auto setup = lsdtt::load_hype_setup_directory(in);
    const forcing_series fs = load_forcing(in);
    lsdtt::lsdtt_hype_model_cpu on(setup), off(setup);
    off.set_store_unit_fluxes(false);
    bool flows_same = true, states_same = true;
    for (std::size_t t = 0; t < fs.prec_mm.size(); ++t)
    {
        on.step(fs.prec_mm[t], fs.temp_c[t], fs.doy[t]);
        off.step(fs.prec_mm[t], fs.temp_c[t], fs.doy[t]);
        if (!same_bits(on.runoff_mm(), off.runoff_mm()) ||
            !same_bits(on.outflow_m3s(), off.outflow_m3s()) ||
            !same_bits(on.local_outflow_m3s(), off.local_outflow_m3s())) flows_same = false;
        for (std::size_t k = 0; k < on.units().size(); ++k)
            if (!same_state(on.units()[k].state, off.units()[k].state)) states_same = false;
    }
    check(off.unit_fluxes().empty() && !off.store_unit_fluxes(), "flux storage off: unit_fluxes() is empty");
    check(flows_same, "flux storage off: sub-basin flows bitwise identical, every step");
    check(states_same, "flux storage off: unit states bitwise identical, every step");
    check(same_bits(on.balance().closure_error_m3(), off.balance().closure_error_m3()),
          "flux storage off: balance bitwise identical");
}

// WP7 B2: an exception thrown by a unit inside the parallel land loop must
// leave step() as an exception (escaping an OpenMP region would terminate),
// and the reported failure must be the lowest-index failing unit at any
// thread count. Validation leaves no input that makes the column fail, so the
// model's test hook makes every unit from index k_first upward throw; with
// many failing units spread over the threads, the report must still name
// unit k_first.
void parallel_throw_checks(const std::string& oracle)
{
    std::cout << "\n--- exception from the parallel land loop ---\n";
    const std::string in = oracle + "/default_all/lsdtt_input";
    const auto setup = lsdtt::load_hype_setup_directory(in);
    const forcing_series fs = load_forcing(in);

    auto run = [&](int n_threads, std::size_t first, std::string& message) -> bool
    {
#ifdef _OPENMP
        omp_set_num_threads(n_threads);
#else
        (void)n_threads;
#endif
        lsdtt::lsdtt_hype_model_cpu model(setup);
        model.step(fs.prec_mm[0], fs.temp_c[0], fs.doy[0]);           // one clean step
        model.inject_failure_from_unit_for_testing(first);
        try
        {
            model.step(fs.prec_mm[1], fs.temp_c[1], fs.doy[1]);
        }
        catch (const std::runtime_error& e)
        {
            message = e.what();
            return true;
        }
        return false;
    };

#ifdef _OPENMP
    const int n_max = std::max(2, omp_get_max_threads());
#else
    const int n_max = 1;
#endif
    const std::size_t first = 1;   // default_all has 6 units; units 1..5 all throw
    std::string msg_one, msg_many;
    const bool threw_one  = run(1, first, msg_one);
    const bool threw_many = run(n_max, first, msg_many);
#ifdef _OPENMP
    omp_set_num_threads(n_max);
#endif
    std::cout << "  " << n_max << " threads: " << msg_many << "\n";
    const std::string want = "injected test failure in unit " + std::to_string(first);
    check(threw_one, "a unit failure propagates out of step() as std::runtime_error at 1 thread");
    check(threw_many, "a unit failure propagates out of step() as std::runtime_error at " +
                          std::to_string(n_max) + " threads");
    check(msg_many.find(want) != std::string::npos && msg_one == msg_many,
          "the lowest-index failing unit is reported, the same at 1 and " + std::to_string(n_max) + " threads");
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-hype-model <oracle_dir>\n";
        return 2;
    }
    const std::string oracle = argv[1];
    try
    {
        run_case(oracle, "default_all",
                 {"snow", "melt", "evap", "infi", "ros1", "ros2", "crun", "sml1", "sml2", "sml3"},
                 {"cout", "crun"});
        run_case(oracle, "routing", {}, {"cout", "cloc", "crun"});
        observed_pet_checks(oracle);
        balance_off_checks(oracle);
        parallel_throw_checks(oracle);
        unit_flux_storage_checks(oracle);

        bool threw = false;
        try { lsdtt::lsdtt_hype_model_cpu m(lsdtt::load_hype_setup_directory(oracle + "/routing/lsdtt_input"), 3600.0); }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "a non-daily timestep is rejected");
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
