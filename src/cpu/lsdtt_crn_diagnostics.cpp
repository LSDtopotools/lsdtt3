// lsdtt_crn_diagnostics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

// =============================================================================
//  DO NOT CHANGE THE NUMBERS THIS CODE PRODUCES — ask Simon Mudd first.
//  See lsdtt_crn_diagnostics.hpp for what this file is and how it is pinned.
// =============================================================================
#include "lsdtt_crn_diagnostics.hpp"
#include "lsdtt_crn_particle.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lsdtt
{

const char* to_string(CrnSolveStatus s)
{
    switch (s)
    {
        case CrnSolveStatus::Converged:    return "converged";
        case CrnSolveStatus::NoSolution:   return "no solution";
        case CrnSolveStatus::NotUnique:    return "not unique";
        case CrnSolveStatus::NotConverged: return "not converged";
    }
    return "unknown";
}

const char* crn_solve_status_token(CrnSolveStatus s)
{
    switch (s)
    {
        case CrnSolveStatus::Converged:    return "Converged";
        case CrnSolveStatus::NoSolution:   return "NoSolution";
        case CrnSolveStatus::NotUnique:    return "NotUnique";
        case CrnSolveStatus::NotConverged: return "NotConverged";
    }
    return "Unknown";
}

int crn_solve_status_code(CrnSolveStatus s)
{
    switch (s)
    {
        case CrnSolveStatus::Converged:    return 0;
        case CrnSolveStatus::NoSolution:   return 1;
        case CrnSolveStatus::NotUnique:    return 2;
        case CrnSolveStatus::NotConverged: return 3;
    }
    return 4;
}

CrnUncertaintySummary summarise_crn_uncertainty(
    const std::vector<double>& values,
    const double center,
    const bool   add_production_uncertainty,
    const double production_rel_sigma,
    const bool   add_decay_uncertainty,
    const double decay_rel_sigma)
{
    CrnUncertaintySummary u;
    if (values.empty()) return u;
    double vmin = values[0], vmax = values[0];
    for (double v : values)
    {
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    u.lo = vmin;
    u.hi = vmax;
    u.unc_int = 0.5 * std::fabs(vmax - vmin);
    const double up = add_production_uncertainty
                          ? std::fabs(center) * production_rel_sigma : 0.0;
    const double ud = add_decay_uncertainty
                          ? std::fabs(center) * decay_rel_sigma : 0.0;
    u.unc_ext = std::sqrt(u.unc_int * u.unc_int + up * up + ud * ud);
    u.n_valid = static_cast<int>(values.size());
    return u;
}

namespace
{

// Per-nuclide accessors, so the three-way branch lives in one place instead of
// being repeated at every use.
struct NuclideView
{
    double P0;
    double lambda;
    double F[4];
    double Gamma[4];
};

NuclideView view_of(const std::string& nuclide, const lsdtt_crn_parameters& p)
{
    NuclideView v{};
    for (int j = 0; j < 4; ++j) v.Gamma[j] = p.get_Gamma(j);
    if (nuclide == "Al26")
    {
        v.P0 = p.get_P0_26Al(); v.lambda = p.get_lambda_26Al();
        for (int j = 0; j < 4; ++j) v.F[j] = p.get_F_26Al(j);
    }
    else if (nuclide == "C14")
    {
        v.P0 = p.get_P0_14C(); v.lambda = p.get_lambda_14C();
        for (int j = 0; j < 4; ++j) v.F[j] = p.get_F_14C(j);
    }
    else if (nuclide == "Be10")
    {
        v.P0 = p.get_P0_10Be(); v.lambda = p.get_lambda_10Be();
        for (int j = 0; j < 4; ++j) v.F[j] = p.get_F_10Be(j);
    }
    else
    {
        throw std::runtime_error("lsdtt_crn_diagnostics: unknown nuclide \"" +
                                 nuclide + "\". Use Be10, Al26 or C14.");
    }
    return v;
}


// Steady-state concentration from a prebuilt view.  The bisections below call
// this tens of thousands of times per sample, so the view is hoisted out of the
// loop rather than rebuilt (three std::string comparisons) on each evaluation.
double ss_conc_view(const NuclideView& v, double eff_erate_g_cm2_yr)
{
    double sum = 0.0;
    for (int j = 0; j < 4; ++j)
        sum += v.F[j] * v.Gamma[j] / (eff_erate_g_cm2_yr + v.Gamma[j] * v.lambda);
    return v.P0 * sum;
}

double ss_ratio_view(const NuclideView& num, const NuclideView& den,
                     double eff_erate_g_cm2_yr)
{
    const double d = ss_conc_view(den, eff_erate_g_cm2_yr);
    if (d <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return ss_conc_view(num, eff_erate_g_cm2_yr) / d;
}

lsdtt_crn_parameters scaled_params(double total_shielding,
                                   const std::string& muon_scheme)
{
    lsdtt_crn_parameters p;
    p.reset_scaling(muon_scheme);
    // Same convention as the raster pipeline and the column: one factor,
    // through the Vermeesch redistribution, no guard.
    p.scale_F_values(total_shielding, "Be10");
    p.scale_F_values(total_shielding, "Al26");
    p.scale_F_values(total_shielding, "C14");
    return p;
}

// Bounds on the erosion-rate search, in g/cm2/yr. 1e-7 is ~4e-5 mm/kyr and 10.0
// is ~38 m/kyr; nothing geomorphic lies outside, and the concentration is
// monotone across the whole range, so these only ever fail to bracket when the
// measurement itself is outside the model's reach.
constexpr double k_erate_lo = 1.0e-7;
constexpr double k_erate_hi = 10.0;
constexpr int    k_bisect_iters = 80;    // 10/2^80 is far below the ULP of any
                                         // erate in [1e-7,10]; 200 spent 120 iterations
                                         // sitting on the same double.

// Invert a monotone-decreasing function f(eps) = target by bisection.
// Returns false when the target is not bracketed, which is a real answer.
template <typename F>
bool bisect_decreasing_in(F&& f, double target, double lo, double hi, double& out)
{
    const double f_lo = f(lo), f_hi = f(hi);
    // f decreases with eps, so a bracketed target satisfies f(hi) <= target <= f(lo)
    if (!(target <= f_lo && target >= f_hi)) return false;
    for (int i = 0; i < k_bisect_iters; ++i)
    {
        const double mid = 0.5 * (lo + hi);
        if (f(mid) > target) lo = mid; else hi = mid;
    }
    out = 0.5 * (lo + hi);
    return true;
}

// Whole-range convenience wrapper.
template <typename F>
bool bisect_decreasing(F&& f, double target, double& out)
{
    return bisect_decreasing_in(std::forward<F>(f), target,
                                k_erate_lo, k_erate_hi, out);
}

} // namespace

//
double crn_steady_state_conc(const std::string& nuclide,
                             double eff_erate_g_cm2_yr,
                             const lsdtt_crn_parameters& params)
{
    return ss_conc_view(view_of(nuclide, params), eff_erate_g_cm2_yr);
}

//
double crn_steady_state_ratio(const std::string& nuclide_num,
                              const std::string& nuclide_den,
                              double eff_erate_g_cm2_yr,
                              const lsdtt_crn_parameters& params)
{
    const double den = crn_steady_state_conc(nuclide_den, eff_erate_g_cm2_yr, params);
    if (den <= 0.0) return 0.0;
    return crn_steady_state_conc(nuclide_num, eff_erate_g_cm2_yr, params) / den;
}

//
double crn_apparent_erate_from_conc(const std::string& nuclide,
                                    double conc_atoms_g,
                                    double total_scaling,
                                    const std::string& muon_scheme,
                                    CrnSolveStatus& status)
{
    status = CrnSolveStatus::NoSolution;
    if (!(conc_atoms_g > 0.0)) return 0.0;

    const lsdtt_crn_parameters params = scaled_params(total_scaling, muon_scheme);
    const NuclideView v = view_of(nuclide, params);   // throws on a bad nuclide

    double eps = 0.0;
    if (!bisect_decreasing([&](double e) { return ss_conc_view(v, e); },
                           conc_atoms_g, eps))
        return 0.0;

    status = CrnSolveStatus::Converged;
    return eps;
}

//
CrnBurialResult crn_burial_diagnostics_pair(const std::string& nuclide_fast,
                                            const std::string& nuclide_slow,
                                            double conc_fast_atoms_g,
                                            double conc_slow_atoms_g,
                                            double total_shielding,
                                            const std::string& muon_scheme)
{
    CrnBurialResult out{};
    out.status = CrnSolveStatus::NoSolution;
    out.nuclide_fast = nuclide_fast;
    out.nuclide_slow = nuclide_slow;
    out.ratio_measured = (conc_slow_atoms_g > 0.0)
        ? conc_fast_atoms_g / conc_slow_atoms_g : 0.0;

    if (!(conc_slow_atoms_g > 0.0) || !(conc_fast_atoms_g > 0.0))
        return out;

    const lsdtt_crn_parameters params = scaled_params(total_shielding, muon_scheme);
    const NuclideView v_fast = view_of(nuclide_fast, params);
    const NuclideView v_slow = view_of(nuclide_slow, params);
    const double lam_fast = v_fast.lambda;
    const double lam_slow = v_slow.lambda;
    const double dlam  = lam_fast - lam_slow;      // should be > 0 for burial solve
    if (!(dlam > 0.0))
    {
        out.status = CrnSolveStatus::NoSolution;
        return out;
    }

    // The single-nuclide answer, for comparison: what a slow-nuclide-only study
    // reports.
    double eps_apparent = 0.0;
    const bool app_ok = bisect_decreasing(
        [&](double e) { return ss_conc_view(v_slow, e); },
        conc_slow_atoms_g, eps_apparent);
    out.apparent_erate_no_burial = app_ok ? eps_apparent : 0.0;

    // Is burial required at all? Compare the measured ratio with the
    // steady-state ratio at the apparent rate.
    if (app_ok)
    {
        out.ratio_ss_at_solution =
            ss_ratio_view(v_fast, v_slow, eps_apparent);
        if (out.ratio_measured >= out.ratio_ss_at_solution)
        {
            // On or above the steady-state line: a simple steady-erosion history
            // explains the pair. Not an error -- the common case for an actively
            // eroding surface.
            out.status = CrnSolveStatus::Converged;
            out.erate_preburial_g_cm2_yr = eps_apparent;
            out.burial_age_yr = 0.0;
            out.burial_required = false;
            return out;
        }
    }
    out.burial_required = true;

    // Burial required. Eliminate t analytically:
    //   R_meas = R_ss(eps) * exp(-dlam * t)  =>  t(eps) = ln(R_ss/R_meas)/dlam
    //   N10_meas = N10_ss(eps) * exp(-lam10 * t(eps))
    // The right-hand side decreases monotonically in eps (both factors do), so a
    // single bracketed bisection resolves it.
    auto predicted_N10 = [&](double e) {
        const double r_ss = ss_ratio_view(v_fast, v_slow, e);
        if (!(r_ss > out.ratio_measured)) return 0.0;   // no burial time exists
        const double t = std::log(r_ss / out.ratio_measured) / dlam;
        return ss_conc_view(v_slow, e) * std::exp(-lam_slow * t);
    };

    // The search must start where a burial time EXISTS. R_ss rises with erosion
    // rate (less time to decay), so below some e_min the steady-state ratio is
    // already at or under the measured one and no positive burial time can
    // reproduce it -- predicted_N10 is clipped to 0 there, which destroys the
    // monotonicity bisection needs. Locate that lower edge first, then bisect
    // on [e_min, hi] where the function really is monotone decreasing.
    //
    // Found by the B1 round trip failing at a slow rate with a short burial,
    // where e_min sits well inside the search range; without this the solver
    // reported NoSolution for a sample it could in fact resolve exactly.
    double e_min = k_erate_lo;
    {
        double lo = k_erate_lo, hi = k_erate_hi;
        if (ss_ratio_view(v_fast, v_slow, hi)
              <= out.ratio_measured)
        {
            // Even the fastest erosion cannot produce a ratio above the measured
            // one: the pair is off the model surface entirely.
            out.status = CrnSolveStatus::NoSolution;
            return out;
        }
        for (int i = 0; i < k_bisect_iters; ++i)
        {
            const double mid = 0.5 * (lo + hi);
            if (ss_ratio_view(v_fast, v_slow, mid)
                    <= out.ratio_measured) lo = mid; else hi = mid;
        }
        e_min = hi;   // smallest rate admitting a positive burial time
    }

    double eps = 0.0;
    if (!bisect_decreasing_in(predicted_N10, conc_slow_atoms_g, e_min,
                              k_erate_hi, eps))
    {
        // The pair lies off the model surface entirely: no (eps, t) reproduces
        // both. Report it rather than returning a nearest fit.
        out.status = CrnSolveStatus::NoSolution;
        return out;
    }

    const double r_ss = ss_ratio_view(v_fast, v_slow, eps);
    out.erate_preburial_g_cm2_yr = eps;
    out.ratio_ss_at_solution     = r_ss;
    out.burial_age_yr            = std::log(r_ss / out.ratio_measured) / dlam;
    out.status                   = CrnSolveStatus::Converged;
    return out;
}

CrnBurialResult crn_burial_diagnostics(double conc_26Al_atoms_g,
                                       double conc_10Be_atoms_g,
                                       double total_shielding,
                                       const std::string& muon_scheme)
{
    return crn_burial_diagnostics_pair("Al26", "Be10",
                                       conc_26Al_atoms_g,
                                       conc_10Be_atoms_g,
                                       total_shielding,
                                       muon_scheme);
}

//
CrnStepPairResult crn_invert_step_change_pair(
    const std::string& nuclide_fast,
    const std::string& nuclide_slow,
    double conc_fast_atoms_g,
    double conc_slow_atoms_g,
    double eff_erate_old_g_cm2_yr,
    double total_shielding,
    const std::string& muon_scheme)
{
    CrnStepPairResult out{};
    out.status = CrnSolveStatus::NoSolution;
    if (!(conc_fast_atoms_g > 0.0) || !(conc_slow_atoms_g > 0.0))
        return out;

    const lsdtt_crn_parameters params = scaled_params(total_shielding, muon_scheme);

    // Forward model: the SAME step-change evolution the raster kernel uses --
    // start at the depth the new rate will have exhumed, initialise at steady
    // state under the old rate, evolve. Reusing lsdtt_crn_particle rather than
    // restating the algebra is deliberate: a diagnostic that inverted a
    // *different* forward model than the predictor would be worse than useless.
    auto forward = [&](const std::string& nuc, double e_new, double t) {
        const double d_start = e_new * t;
        lsdtt_crn_particle part(d_start);
        part.set_conc_10Be(0.0); part.set_conc_26Al(0.0); part.set_conc_14C(0.0);
        if (nuc == "Al26")      part.update_26Al_ss(eff_erate_old_g_cm2_yr, params, d_start, d_start);
        else if (nuc == "C14")  part.update_14C_ss (eff_erate_old_g_cm2_yr, params, d_start, d_start);
        else                    part.update_10Be_ss(eff_erate_old_g_cm2_yr, params, d_start, d_start);
        part.update_conc(t, e_new, params);
        return nuc == "Al26" ? part.get_conc_26Al()
             : nuc == "C14"  ? part.get_conc_14C()
                             : part.get_conc_10Be();
    };

    // Inner solve: at a fixed t, find the rate reproducing the SLOW nuclide.
    // Monotone decreasing in e_new at fixed t, so bisection is safe.
    auto erate_for_time = [&](double t, double& e_out) {
        return bisect_decreasing(
            [&](double e) { return forward(nuclide_slow, e, t); },
            conc_slow_atoms_g, e_out);
    };

    // Outer residual: having matched the slow nuclide, how far off is the fast?
    auto residual = [&](double t, bool& ok) {
        double e = 0.0;
        ok = erate_for_time(t, e);
        if (!ok) return 0.0;
        return forward(nuclide_fast, e, t) - conc_fast_atoms_g;
    };

    // Scan the time axis on a log grid to BRACKET sign changes before bisecting.
    // A 2-D Newton iteration would happily converge on one root and never reveal
    // that another existed; counting brackets is what makes non-uniqueness
    // detectable instead of invisible.
    const double t_lo = 10.0, t_hi = 1.0e7;
    const int n_scan = 240;
    std::vector<std::pair<double, double>> brackets;
    double t_prev = t_lo;
    bool ok_prev = false;
    double r_prev = residual(t_prev, ok_prev);
    for (int i = 1; i <= n_scan; ++i)
    {
        const double frac = static_cast<double>(i) / n_scan;
        const double t = t_lo * std::pow(t_hi / t_lo, frac);
        bool ok = false;
        const double r = residual(t, ok);
        if (ok && ok_prev && ((r_prev <= 0.0 && r >= 0.0) || (r_prev >= 0.0 && r <= 0.0)))
            brackets.emplace_back(t_prev, t);
        t_prev = t; r_prev = r; ok_prev = ok;
    }

    out.n_brackets = static_cast<int>(brackets.size());
    if (brackets.empty())  { out.status = CrnSolveStatus::NoSolution; return out; }
    if (brackets.size() > 1)
    {
        // Genuinely ambiguous: the pair does not resolve (rate, time) here.
        // Reporting this is the correct answer -- returning the first root would
        // present a coin flip as a measurement.
        out.status = CrnSolveStatus::NotUnique;
        return out;
    }

    double lo = brackets[0].first, hi = brackets[0].second;
    bool ok_lo = false;
    const double r_lo = residual(lo, ok_lo);
    for (int i = 0; i < 100; ++i)
    {
        const double mid = std::sqrt(lo * hi);      // geometric: the axis is log
        bool ok = false;
        const double r = residual(mid, ok);
        if (!ok) break;
        if ((r <= 0.0) == (r_lo <= 0.0)) lo = mid; else hi = mid;
    }
    const double t_sol = std::sqrt(lo * hi);
    double e_sol = 0.0;
    if (!erate_for_time(t_sol, e_sol))
    { out.status = CrnSolveStatus::NotConverged; return out; }

    out.status               = CrnSolveStatus::Converged;
    out.time_since_change_yr = t_sol;
    out.erate_new_g_cm2_yr   = e_sol;
    out.conc_fast_predicted  = forward(nuclide_fast, e_sol, t_sol);
    out.conc_slow_predicted  = forward(nuclide_slow, e_sol, t_sol);
    return out;
}

} // namespace lsdtt
