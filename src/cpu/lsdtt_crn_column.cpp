#include "lsdtt_crn_column.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace lsdtt
{

// ─────────────────────────────────────────────────────────────────────────────
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

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_column::set_erosion_time_series(const std::vector<double>& time_yr,
                                                const std::vector<double>& rate_m_yr,
                                                const std::vector<double>& removal_m)
{
    time_vec_    = time_yr;
    rate_vec_    = rate_m_yr;
    removal_vec_ = removal_m;
    check_erosion_time_series();
}

// ─────────────────────────────────────────────────────────────────────────────
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
}

// ─────────────────────────────────────────────────────────────────────────────
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

// ─────────────────────────────────────────────────────────────────────────────
lsdtt_crn_parameters lsdtt_crn_column::make_params(double production_scaling,
                                                    const std::string& muon_scheme) const
{
    lsdtt_crn_parameters p;
    p.reset_scaling(muon_scheme);
    p.scale_F_values(production_scaling, "Be10");
    // All three nuclides use the same Γ and F ratios after scaling;
    // scale_F_values writes only the chosen nuclide's F[].  We need a unified
    // approach: scale all three with the same shielding factor.
    // Reset and re-scale individually, then merge.
    // Actually in lsdtt3, scale_F_values applies the SAME effective depth to
    // the chosen nuclide.  For the column, we apply the same production_scaling
    // to all three.  We need to call it three times.
    // Re-do cleanly:
    p.reset_scaling(muon_scheme);
    // scale_F_values modifies F for the named nuclide in-place; call for each.
    // But reset_scaling wipes all, so we have to call each without resetting.
    // Simplest: scale_F_values should be idempotent; call each separately
    // after reset.  Here we must NOT reset between calls.
    // We need a separate approach: just set production_scaling into the
    // CRNParameters by abusing scale_F_values three times.
    // (scale_F_values(s, nuclide) multiplies F[nuclide][] by effective depth
    //  correction that maps scalar shielding s to an equivalent F[].)
    // For a column simulation production_scaling == Stone factor × topo shielding,
    // a dimensionless number ≤ 1.  scale_F_values expects total_shielding in [0,1].
    // Cap at 1 to satisfy the guard.
    double s = std::min(1.0, std::max(0.0, production_scaling));
    if (s < 1.0)
    {
        p.scale_F_values(s, "Be10");
        p.scale_F_values(s, "Al26");
        p.scale_F_values(s, "C14");
    }
    return p;
}

// ─────────────────────────────────────────────────────────────────────────────
std::vector<CrnColumnSample>
lsdtt_crn_column::run_transient(const std::vector<double>& sampling_depths_m,
                                 double production_scaling,
                                 const std::string& muon_scheme,
                                 double initial_10Be,
                                 double initial_26Al,
                                 double initial_14C) const
{
    int n_steps   = static_cast<int>(time_vec_.size());
    int n_samples = static_cast<int>(sampling_depths_m.size());

    if (n_steps < 2)
        throw std::runtime_error("lsdtt_crn_column::run_transient: erosion time series "
                                 "not loaded or has fewer than 2 rows.");

    lsdtt_crn_parameters params = make_params(production_scaling, muon_scheme);

    // Background (long-term) erosion rate
    double ss_rate_m_yr  = rate_vec_[n_steps-1];
    double eff_ss_rate   = to_eff(ss_rate_m_yr);   // g/cm²/yr

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
        double start_d_eff    = to_eff(start_depth_m);   // g/cm²

        lsdtt_crn_particle part(start_d_eff);

        // ── Initial concentration ──────────────────────────────────────────
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

        // ── Evolve through erosion history (oldest → youngest) ─────────────
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

// ─────────────────────────────────────────────────────────────────────────────
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

// ─────────────────────────────────────────────────────────────────────────────
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

} // namespace lsdtt
