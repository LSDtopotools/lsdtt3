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
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
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

double class_value(const std::string& v, const lsdtt::hype_model_unit& u)
{
    const auto& f = u.fluxes;
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
            for (std::size_t k = 0; k < model.units().size(); ++k) cls[v][k][t] = class_value(v, model.units()[k]);
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
        for (const auto& u : model.units())
        {
            if (u.fluxes.epot_mm != obs) epot_replaced = false;
            if (u.fluxes.evap_mm > obs + 1e-12) aet_bounded = false;
            if (!(u.fluxes.class_temperature_c > u.constants.landuse.ttmp))
            {
                ++cold_steps;
                if (u.fluxes.evap_mm != 0.0) aet_off_when_cold = false;
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
