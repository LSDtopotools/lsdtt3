// =============================================================================
//  DO NOT CHANGE THE NUMBERS THIS CODE PRODUCES — ask Simon Mudd first.
//
//  The cosmogenic calculations in this file (the CAIRN erosion-rate method and
//  the step-change concentration solutions) were checked by hand over several
//  years. To stop anyone changing those numbers by accident, the test suite
//  keeps a stored list of known-correct output values in
//  tests/expected/cosmo/numerics_sentinel.csv. The tests in
//  tests/cpp/test_cosmo_physics.cpp re-run these calculations and compare the
//  results against that stored list to 12 decimal places.
//
//  So if you edit this file and any output value changes — even from a refactor
//  or an algebra rewrite you believe is harmless — those tests will fail. When
//  that happens, do NOT edit the stored list to make the tests pass again. A
//  changed value has to be reviewed and approved by Simon Mudd before it is
//  committed.
// =============================================================================
#include "lsdtt_crn_column.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace lsdtt
{

// 
void lsdtt_crn_column::read_erosion_time_series(const std::string& csv_path)
{
    std::ifstream ifs(csv_path);
    if (!ifs.is_open())
        throw std::runtime_error("lsdtt_crn_column: cannot open \"" + csv_path + "\"");

    std::vector<double> t_v, r_v, rem_v;
    std::string line;
    std::getline(ifs, line);   // skip header

    while (std::getline(ifs, line))
    {
        if (line.empty()) continue;
        // strip control chars
        line.erase(std::remove_if(line.begin(), line.end(),
                                  [](char c){ return c == '\r' || c == '\n'; }),
                   line.end());
        std::istringstream ss(line);
        std::string t_s, r_s, rem_s;
        if (!std::getline(ss, t_s, ',') ||
            !std::getline(ss, r_s, ',') ||
            !std::getline(ss, rem_s, ','))
            throw std::runtime_error("lsdtt_crn_column: erosion history file must have "
                                     "three comma-separated columns: time,rate,removal");
        t_v.push_back(std::stod(t_s));
        r_v.push_back(std::stod(r_s));
        rem_v.push_back(std::stod(rem_s));
    }

    set_erosion_time_series(t_v, r_v, rem_v);
}

// 
void lsdtt_crn_column::set_erosion_time_series(const std::vector<double>& time_yr,
                                                const std::vector<double>& rate_m_yr,
                                                const std::vector<double>& removal_m)
{
    time_vec_    = time_yr;
    rate_vec_    = rate_m_yr;
    removal_vec_ = removal_m;
    check_erosion_time_series();
}

// 
void lsdtt_crn_column::check_erosion_time_series() const
{
    int n = static_cast<int>(time_vec_.size());
    if (n == 0)
        throw std::runtime_error("lsdtt_crn_column: erosion time series is empty.");
    if (time_vec_[n-1] != -9999.0)
        throw std::runtime_error("lsdtt_crn_column: last row must have time = -9999 "
                                 "(long-term background erosion rate).");
    if (rate_vec_[n-1] <= 0.0)
        throw std::runtime_error("lsdtt_crn_column: background erosion rate (last row) "
                                 "must be > 0.");
    for (int i = 0; i < n-1; ++i)
    {
        if (rate_vec_[i] != 0.0 && removal_vec_[i] != 0.0)
            throw std::runtime_error("lsdtt_crn_column: row " + std::to_string(i) +
                                     " has both nonzero rate and nonzero removal.");
        if (removal_vec_[i] != 0.0 && i > 0 && time_vec_[i] != time_vec_[i-1])
            throw std::runtime_error("lsdtt_crn_column: removal event at row " +
                                     std::to_string(i) + " must have the same time "
                                     "as the preceding row.");
    }
    if (rate_vec_[0] == 0.0 && time_vec_[0] != 0.0)
        throw std::runtime_error("lsdtt_crn_column: first removal event must be at time=0.");

    // TIME ORDERING (added 2026-07-28). run_transient walks the record with
    // dt = time_vec_[ts] - time_vec_[ts-1], and dt = time_vec_[0] for the most
    // recent step, so index 0 must be the YOUNGEST row and yr BP must increase
    // with index. Without this check a reversed record fed negative dt straight
    // into update_conc and returned a confidently wrong number -- measured
    // 72357 forward against 178467 reversed -- with no error, no warning and no
    // NaN. Silence was the whole problem, so this is a hard failure.
    //
    // Equal times are legal ONLY for a removal event, which is instantaneous
    // and is required above to share the time of the row before it.
    if (time_vec_[0] < 0.0)
        throw std::runtime_error("lsdtt_crn_column: first row time must be >= 0 "
                                 "(yr BP of the most recent step); got " +
                                 std::to_string(time_vec_[0]) + ".");
    for (int i = 1; i < n-1; ++i)
    {
        if (time_vec_[i] < time_vec_[i-1])
            throw std::runtime_error(
                "lsdtt_crn_column: erosion history must be ordered with index 0 "
                "the MOST RECENT row and yr BP increasing with index, but row " +
                std::to_string(i) + " (time " + std::to_string(time_vec_[i]) +
                ") is younger than row " + std::to_string(i-1) + " (time " +
                std::to_string(time_vec_[i-1]) + "). A reversed record would "
                "produce negative timesteps and a silently wrong result. "
                "Reverse the rows, keeping time = -9999 last.");
        if (time_vec_[i] == time_vec_[i-1] && removal_vec_[i] == 0.0)
            throw std::runtime_error(
                "lsdtt_crn_column: rows " + std::to_string(i-1) + " and " +
                std::to_string(i) + " share time " +
                std::to_string(time_vec_[i]) + " but row " + std::to_string(i) +
                " is not a removal event, so it would contribute a zero-length "
                "timestep.");
    }
}

// 
double lsdtt_crn_column::calculate_depth_change() const
{
    int n = static_cast<int>(time_vec_.size());
    double d = 0.0;

    // First step
    if (rate_vec_[0] != 0.0)
        d += rate_vec_[0] * time_vec_[0];
    else
        d += removal_vec_[0];

    // Remaining steps (excluding background, the last row)
    for (int i = 1; i < n-1; ++i)
    {
        if (rate_vec_[i] != 0.0)
            d += rate_vec_[i] * (time_vec_[i] - time_vec_[i-1]);
        else
            d += removal_vec_[i];
    }
    return d;
}

// 
lsdtt_crn_parameters lsdtt_crn_column::make_params(double production_scaling,
                                                    const std::string& muon_scheme) const
{
    lsdtt_crn_parameters p;
    p.reset_scaling(muon_scheme);
    // Apply Stone / site production scaling by multiplying P0 directly.
    // This handles both factors < 1 and > 1 (high-latitude enhancement).
    // scale_F_values is reserved for shielding < 1 in the simple modes.
    if (production_scaling != 1.0)
        p.apply_production_scaling(production_scaling);
    return p;
}

// 
std::vector<CrnColumnSample>
lsdtt_crn_column::run_transient(const std::vector<double>& sampling_depths_m,
                                 double production_scaling,
                                 const std::string& muon_scheme,
                                 double initial_10Be,
                                 double initial_26Al,
                                 double initial_14C) const
{
    int n_steps   = static_cast<int>(time_vec_.size());

    if (n_steps < 2)
        throw std::runtime_error("lsdtt_crn_column::run_transient: erosion time series "
                                 "not loaded or has fewer than 2 rows.");

    return run_transient_impl(sampling_depths_m,
                              make_params(production_scaling, muon_scheme),
                              initial_10Be, initial_26Al, initial_14C);
}

// Raster-convention overload: production scaling and shielding go through
// scale_F_values as one factor (see make_params_shielded).
std::vector<CrnColumnSample>
lsdtt_crn_column::run_transient(const std::vector<double>& sampling_depths_m,
                                 double production_scaling,
                                 double total_shielding,
                                 const std::string& muon_scheme,
                                 double initial_10Be,
                                 double initial_26Al,
                                 double initial_14C) const
{
    return run_transient_impl(sampling_depths_m,
                              make_params_shielded(production_scaling,
                                                   total_shielding, muon_scheme),
                              initial_10Be, initial_26Al, initial_14C);
}

// Shared body. The ONLY thing the two overloads differ in is the parameter
// object, so the evolution itself lives in exactly one place.
std::vector<CrnColumnSample>
lsdtt_crn_column::run_transient_impl(const std::vector<double>& sampling_depths_m,
                                      const lsdtt_crn_parameters& params,
                                      double initial_10Be,
                                      double initial_26Al,
                                      double initial_14C) const
{
    int n_steps   = static_cast<int>(time_vec_.size());
    int n_samples = static_cast<int>(sampling_depths_m.size());

    if (n_steps < 2)
        throw std::runtime_error("lsdtt_crn_column::run_transient: erosion time series "
                                 "not loaded or has fewer than 2 rows.");

    // Background (long-term) erosion rate
    double ss_rate_m_yr  = rate_vec_[n_steps-1];
    double eff_ss_rate   = to_eff(ss_rate_m_yr);   // g/cm^2/yr

    // Depth change: how much deeper each particle was at the start of the record
    double d_change = calculate_depth_change();

    bool use_initial_10Be = (initial_10Be >= 0.0);
    bool use_initial_26Al = (initial_26Al >= 0.0);
    bool use_initial_14C  = (initial_14C  >= 0.0);

    std::vector<CrnColumnSample> results;
    results.reserve(static_cast<std::size_t>(n_samples));

    for (int s = 0; s < n_samples; ++s)
    {
        double start_depth_m  = sampling_depths_m[s] + d_change;
        double start_d_eff    = to_eff(start_depth_m);   // g/cm^2

        lsdtt_crn_particle part(start_d_eff);

        //  Initial concentration 
        if (use_initial_10Be || use_initial_26Al || use_initial_14C)
        {
            // User-specified initial concentration (e.g. pre-existing inventory)
            part.set_conc_10Be(use_initial_10Be ? initial_10Be : 0.0);
            part.set_conc_26Al(use_initial_26Al ? initial_26Al : 0.0);
            part.set_conc_14C (use_initial_14C  ? initial_14C  : 0.0);
        }
        else
        {
            // Initialise at steady state with background rate at the starting depth
            part.init_ss_all(eff_ss_rate, params);
        }

        //  Evolve through erosion history (oldest -> youngest) 
        // Loop from second-to-last row down to row 1, then handle row 0 last.
        for (int ts = n_steps-2; ts >= 1; --ts)
        {
            if (rate_vec_[ts] == 0.0)
            {
                // Instantaneous removal
                part.remove_mass(removal_vec_[ts], rock_density_);
            }
            else
            {
                double dt          = time_vec_[ts] - time_vec_[ts-1];
                double eff_erate   = to_eff(rate_vec_[ts]);
                part.update_conc(dt, eff_erate, params);
                part.remove_mass(rate_vec_[ts] * dt, rock_density_);
            }
        }

        // First row (most recent step before present)
        {
            double dt        = time_vec_[0];
            double eff_erate = to_eff(rate_vec_[0]);
            part.update_conc(dt, eff_erate, params);
            part.remove_mass(rate_vec_[0] * dt, rock_density_);
        }

        CrnColumnSample out;
        out.depth_m         = sampling_depths_m[s];
        out.effective_depth = part.get_effective_depth();
        out.conc_10Be       = part.get_conc_10Be();
        out.conc_26Al       = part.get_conc_26Al();
        out.conc_14C        = part.get_conc_14C();
        results.push_back(out);
    }

    return results;
}

// 
CrnColumnSample lsdtt_crn_column::steady_state(double d_eff,
                                                double eff_erate,
                                                double production_scaling,
                                                const std::string& muon_scheme) const
{
    lsdtt_crn_parameters params = make_params(production_scaling, muon_scheme);
    lsdtt_crn_particle part(d_eff);
    part.init_ss_all(eff_erate, params);

    CrnColumnSample out;
    out.depth_m         = d_eff / (rock_density_ * 0.1);
    out.effective_depth = d_eff;
    out.conc_10Be       = part.get_conc_10Be();
    out.conc_26Al       = part.get_conc_26Al();
    out.conc_14C        = part.get_conc_14C();
    return out;
}

// 
CrnColumnSample lsdtt_crn_column::step_change(double d_eff,
                                               double eff_erate_old,
                                               double eff_erate_new,
                                               double t_step,
                                               double production_scaling,
                                               const std::string& muon_scheme) const
{
    lsdtt_crn_parameters params = make_params(production_scaling, muon_scheme);

    // Starting depth at beginning of phase 2 (t_step yr ago)
    double d_start = d_eff + eff_erate_new * t_step;
    lsdtt_crn_particle part(d_start);

    // Initial condition: SS at old erosion rate, starting depth
    part.init_ss_all(eff_erate_old, params);

    // Evolve for t_step at new rate (depth goes from d_start to d_eff)
    part.update_conc(t_step, eff_erate_new, params);

    CrnColumnSample out;
    out.depth_m         = d_eff / (rock_density_ * 0.1);
    out.effective_depth = d_eff;
    out.conc_10Be       = part.get_conc_10Be();
    out.conc_26Al       = part.get_conc_26Al();
    out.conc_14C        = part.get_conc_14C();
    return out;
}

//
// ADDITIVE overloads taking shielding separately (2026-07-27). See the header
// for why: shielding must go through the Vermeesch redistribution, not a flat
// P0 multiply, to match the raster pipeline. The two functions above are
// unchanged and remain the single-factor form.
lsdtt_crn_parameters
lsdtt_crn_column::make_params_shielded(double production_scaling,
                                        double total_shielding,
                                        const std::string& muon_scheme) const
{
    // THE RASTER'S RULE, EXACTLY (changed 2026-07-28).
    // calculate_crn_concentration_raster* forms S_tot = production_scaling x
    // shielding and puts the PRODUCT through scale_F_values -- the Vermeesch
    // virtual-attenuation redistribution -- with no guard, including the
    // S_tot > 1 negative-effective-depth branch that reproduces CRONUS v3 at
    // altitude. That is the convention the numerics sentinel pins and every
    // validated comparison exercises, so it is the reference.
    //
    // This previously applied production scaling through a flat P0 multiply
    // (make_params) and only shielding through scale_F_values, which made the
    // column disagree with the raster by -14.87 % at S_tot = 0.5 and +5.39 % at
    // 1.4 -- and, because of a `< 1.0` guard, made it unable to express
    // altitude enhancement at all. Both are gone: one factor, one mechanism.
    //
    // make_params and the single-factor overloads are deliberately NOT changed;
    // they remain the flat-P0 form for callers that want it, and run_transient
    // still uses them unless the shielded overload is called.
    lsdtt_crn_parameters p;
    p.reset_scaling(muon_scheme);
    const double s_tot = production_scaling * total_shielding;
    p.scale_F_values(s_tot, "Be10");
    p.scale_F_values(s_tot, "Al26");
    p.scale_F_values(s_tot, "C14");
    return p;
}

CrnColumnSample lsdtt_crn_column::steady_state(double d_eff,
                                                double eff_erate,
                                                double production_scaling,
                                                double total_shielding,
                                                const std::string& muon_scheme) const
{
    lsdtt_crn_parameters params =
        make_params_shielded(production_scaling, total_shielding, muon_scheme);
    lsdtt_crn_particle part(d_eff);
    part.init_ss_all(eff_erate, params);

    CrnColumnSample out;
    out.depth_m         = d_eff / (rock_density_ * 0.1);
    out.effective_depth = d_eff;
    out.conc_10Be       = part.get_conc_10Be();
    out.conc_26Al       = part.get_conc_26Al();
    out.conc_14C        = part.get_conc_14C();
    return out;
}

CrnColumnSample lsdtt_crn_column::step_change(double d_eff,
                                               double eff_erate_old,
                                               double eff_erate_new,
                                               double t_step,
                                               double production_scaling,
                                               double total_shielding,
                                               const std::string& muon_scheme) const
{
    lsdtt_crn_parameters params =
        make_params_shielded(production_scaling, total_shielding, muon_scheme);

    // Identical evolution to the single-factor overload: only the parameter
    // object differs. Kept as a separate body rather than delegating, because
    // the existing overload takes the combined factor and must not change.
    double d_start = d_eff + eff_erate_new * t_step;
    lsdtt_crn_particle part(d_start);
    part.init_ss_all(eff_erate_old, params);
    part.update_conc(t_step, eff_erate_new, params);

    CrnColumnSample out;
    out.depth_m         = d_eff / (rock_density_ * 0.1);
    out.effective_depth = d_eff;
    out.conc_10Be       = part.get_conc_10Be();
    out.conc_26Al       = part.get_conc_26Al();
    out.conc_14C        = part.get_conc_14C();
    return out;
}

} // namespace lsdtt
