// =============================================================================
//  DO NOT CHANGE THE NUMBERS THIS CODE PRODUCES — ask Simon Mudd first.
//
//  These are cosmogenic diagnostics built on the locked numeric core. They add
//  no new physics: every concentration comes from the existing steady-state and
//  step-change solutions in lsdtt_crn_particle / lsdtt_crn_parameters. What is
//  new is the INVERSION — solving those same forward models backwards — so a
//  change here changes published numbers just as surely as a change to the core.
//  The tests in tests/cpp/test_cosmo_physics.cpp pin them by round trip: a
//  forward model at known inputs, inverted, must return those inputs.
// =============================================================================
//
// Two-nuclide diagnostics.
//
// A single nuclide gives one number and therefore one unknown: an apparent
// erosion rate, on the assumption of steady erosion and no burial. A PAIR of
// nuclides with different half-lives gives two numbers, and can therefore
// resolve two unknowns. This file provides the two cases that matter:
//
//   1. BURIAL / INHERITANCE (Al-26 + Be-10). If a sample was exposed, then
//      buried deeply enough to stop production, both inventories decay, but
//      Al-26 decays roughly twice as fast. The Al26/Be10 ratio therefore falls
//      below the steady-state erosion line by an amount that dates the burial.
//      Solves for (pre-burial erosion rate, burial duration).
//
//   2. STEP CHANGE (any pair, but C-14 + Be-10 is the useful one). After a
//      change in erosion rate the two nuclides re-equilibrate at different
//      speeds, so the pair resolves both the new rate and the time since the
//      change. Solves for (new erosion rate, time since change).
//
// BOTH RETURN A STATUS, AND NON-CONVERGENCE IS NOT A FAILURE OF THE CALLER'S
// DATA. A pair can be genuinely uninformative: outside the resolvable window a
// solution may not exist or may not be unique, and that is a property of the
// physics, not a bug. Reporting "not identifiable" is the correct answer and
// must never be papered over with a nearest-fit.

#ifndef LSDTT_CRN_DIAGNOSTICS_HPP
#define LSDTT_CRN_DIAGNOSTICS_HPP

#include <limits>
#include <string>
#include <vector>
#include "lsdtt_crn_parameters.hpp"

namespace lsdtt
{

enum class CrnSolveStatus
{
    Converged,        // a unique solution was found
    NoSolution,       // the measurements lie outside the forward model's range
    NotUnique,        // more than one solution brackets; the pair cannot resolve it
    NotConverged      // bracketed but the solver did not tighten (should not happen)
};

// Human-readable form, e.g. "no solution". For prose and log lines.
const char* to_string(CrnSolveStatus s);

// Machine-readable single-word token for CSV columns, e.g. "NoSolution". A tidy
// table wants a token a loader can compare without normalising case or
// whitespace, which to_string() deliberately is not.
const char* crn_solve_status_token(CrnSolveStatus s);

// Small integer status for callers that want a numeric column or an exit code.
// Converged is 0 so that "nonzero means look at it" holds.
int crn_solve_status_code(CrnSolveStatus s);

//  Steady-state concentration of one nuclide at a surface sample
// C(eps) = P0 * sum_j F_j * Gamma_j / (eps + Gamma_j * lambda), the same
// analytic form the inversion in lsdtt_cosmo_raster_cpu uses. Exposed here so
// diagnostics and tests can reference one implementation rather than restating
// the algebra.
//
// nuclide: "Be10", "Al26" or "C14". params must already carry the site scaling
// (see lsdtt_crn_parameters::scale_F_values).
double crn_steady_state_conc(const std::string& nuclide,
                             double eff_erate_g_cm2_yr,
                             const lsdtt_crn_parameters& params);

//  Steady-state nuclide ratio at a surface sample
// R_ss(eps) = C_num(eps) / C_den(eps). Monotone in eps over the physical range
// because the shorter-lived nuclide equilibrates faster, which is exactly the
// property the burial solve relies on.
double crn_steady_state_ratio(const std::string& nuclide_num,
                              const std::string& nuclide_den,
                              double eff_erate_g_cm2_yr,
                              const lsdtt_crn_parameters& params);

//  Apparent (steady-state, no-burial) effective erosion rate from a SURFACE
//  concentration
//
// This is the number a single-nuclide study reports. It uses the same monotone
// bisection over crn_steady_state_conc that
// crn_burial_diagnostics_pair::apparent_erate_no_burial already uses, so it adds
// no new physics -- it only makes the existing primitive callable.
//
// SURFACE ONLY. crn_steady_state_conc is the d = 0 form, so a shielded sample
// must not be passed through here: it would be inverted as if it were exposed.
//
// total_scaling: production_scaling x shielding, ONE factor, exactly as
// lsdtt_crn_parameters::scale_F_values expects and as
// lsdtt_crn_column::make_params_shielded forms it. Passing the same product the
// forward model used makes the inversion exactly consistent with it.
//
// Returns 0 and sets status != Converged when the concentration is outside the
// model's reachable range (erosion rates 1e-7 to 10 g/cm2/yr). That is a real
// answer about the measurement, not a solver failure, and must not be papered
// over with a nearest fit.
double crn_apparent_erate_from_conc(const std::string& nuclide,
                                    double conc_atoms_g,
                                    double total_scaling,
                                    const std::string& muon_scheme,
                                    CrnSolveStatus& status);

//  Burial / inheritance diagnostics
struct CrnBurialResult
{
    CrnSolveStatus status;
    std::string nuclide_fast;
    std::string nuclide_slow;
    double erate_preburial_g_cm2_yr; // steady erosion rate before burial
    double burial_age_yr;            // duration of shielded burial; 0 if none
    double ratio_measured;           // N_fast / N_slow as supplied
    double ratio_ss_at_solution;     // the steady-state ratio at the solved rate
    bool   burial_required;          // measured ratio sits below the SS line
    double apparent_erate_no_burial; // rate from the slow nuclide alone, ignoring
                                     // burial -- the number a single-nuclide
                                     // study reports
};

// Solves for (pre-burial erosion rate, burial duration) from a two-nuclide
// pair, assuming: steady erosion to a surface sample, then instantaneous burial
// deep enough that production stops, then no further exposure.
//
// The burial time is eliminated analytically -- from the ratio equation
// t = ln(R_ss(eps)/R_meas) / (lambda_fast - lambda_slow) -- which reduces the
// nominally two-dimensional problem to a ONE-dimensional monotone root find in
// eps. That is why this is robust where a 2-D Newton solve would not be.
//
// nuclide_fast should be the shorter-lived member (larger lambda), so
// lambda_fast > lambda_slow and burial lowers N_fast/N_slow below the
// steady-state line. If lambdas are equal, the pair cannot date burial.
CrnBurialResult crn_burial_diagnostics_pair(const std::string& nuclide_fast,
                                            const std::string& nuclide_slow,
                                            double conc_fast_atoms_g,
                                            double conc_slow_atoms_g,
                                            double total_shielding,
                                            const std::string& muon_scheme);

// Legacy Al-26 / Be-10 convenience wrapper.
// Equivalent to crn_burial_diagnostics_pair("Al26", "Be10", ...)
// and retained for compatibility with existing callers/tests.
//
// Solves for (pre-burial erosion rate, burial duration) from an Al-26 / Be-10
// pair.
//
// If the measured ratio sits ON or ABOVE the steady-state line, no burial is
// required: burial_required is false, burial_age_yr is 0, and the pre-burial
// rate equals the apparent rate. That is a legitimate result, not an error.
//
// total_shielding: production scaling x shielding, as elsewhere (see
// scale_F_values); applied to both nuclides.
CrnBurialResult crn_burial_diagnostics(double conc_26Al_atoms_g,
                                       double conc_10Be_atoms_g,
                                       double total_shielding,
                                       const std::string& muon_scheme);

//  Step-change pair inversion
struct CrnStepPairResult
{
    CrnSolveStatus status;
    double erate_new_g_cm2_yr;   // erosion rate since the change
    double time_since_change_yr; // when the change happened
    double conc_fast_predicted;  // forward-modelled at the solution, for checking
    double conc_slow_predicted;
    int    n_brackets;           // how many sign changes were found; > 1 => NotUnique
};

// Solves for (new erosion rate, time since change) from a pair of nuclides
// measured at the same surface sample, given an assumed prior steady rate.
//
// Reduces the 2-D problem to NESTED 1-D solves rather than a 2-D Newton
// iteration: for a trial time, the slow nuclide is inverted for the rate that
// reproduces its concentration (monotone, bracketed, robust), then the fast
// nuclide's residual is a function of time alone. The time axis is scanned on a
// log grid to find sign changes before bisecting, so a non-unique solution is
// DETECTED and reported rather than silently returning whichever root the
// iteration happened to fall into.
//
// nuclide_fast should be the shorter-lived member (C14, or Al26 against Be10);
// it is the one carrying the timing information.
CrnStepPairResult crn_invert_step_change_pair(
    const std::string& nuclide_fast,
    const std::string& nuclide_slow,
    double conc_fast_atoms_g,
    double conc_slow_atoms_g,
    double eff_erate_old_g_cm2_yr,
    double total_shielding,
    const std::string& muon_scheme);

//  Uncertainty summary over a set of inverted values
// Summarises the spread of a set of solutions obtained by re-inverting at the
// measurement's uncertainty bounds, and combines it with the systematic terms.
//
//   unc_int : INTERNAL (measurement) uncertainty, the half-range of `values`
//   unc_ext : EXTERNAL uncertainty, unc_int combined in quadrature with the
//             production-rate and decay-constant terms, each taken as a
//             relative sigma on |center|
//
// The internal/external split is the usual cosmogenic reporting convention:
// quote internal when comparing samples measured the same way, external when
// comparing against other methods or other laboratories.
struct CrnUncertaintySummary
{
    double lo      = std::numeric_limits<double>::quiet_NaN();
    double hi      = std::numeric_limits<double>::quiet_NaN();
    double unc_int = std::numeric_limits<double>::quiet_NaN();
    double unc_ext = std::numeric_limits<double>::quiet_NaN();
    int    n_valid = 0;
};

// An empty `values` returns the default-constructed summary (all NaN, n_valid
// 0) rather than throwing: a sample whose bounds could not be inverted still
// has a central estimate worth reporting.
CrnUncertaintySummary summarise_crn_uncertainty(
    const std::vector<double>& values,
    double center,
    bool   add_production_uncertainty,
    double production_rel_sigma,
    bool   add_decay_uncertainty,
    double decay_rel_sigma);

} // namespace lsdtt

#endif // LSDTT_CRN_DIAGNOSTICS_HPP
