// lsdtt_landslide_size_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Landslide SIZE DISTRIBUTION: sampling areas and converting them to volumes.
//
// This is the distribution half of the triggered landslide module. It has no
// landscape in it at all -- no DEM, no flow routing, no placement. Feed it a
// distribution and a cell size, get pixel counts and volumes. Placement lives
// in a separate file so that the statistics can be tested on their own.
//
// Design: the size distribution is IMPOSED, not emergent. The user prescribes a
// trigger time and an event volume; this file answers "how many landslides, of
// what sizes". The alternative -- deriving sizes from slope stability, as
// HyLands does (Campforts et al. 2020, GMD 13, 3863-3886) -- is the more
// fundamental approach and is the intended eventual comparison, but
// it makes the landslide flux an outcome of rock-strength calibration, which is
// exactly the variable we do not want in the loop when asking what landslides
// do to cosmogenic nuclide concentrations.
//
// ---------------------------------------------------------------------------
// AREAS ARE DISCRETE. Landslides excavate whole pixels, so the continuous
// probability density is turned into a probability mass function over integer
// pixel counts ONCE, at construction, and sampled by inverse-CDF lookup.
//
// Do NOT sample the continuous density and round. Rounding distorts the small
// end, which is precisely where the rollover lives and where n_min = 1 puts most
// of the probability mass. The pmf construction integrates the density across
// each pixel-count bin instead, which is correct at any cell size.
// ---------------------------------------------------------------------------
//
// Two distribution families:
//
//   1. INVERSE GAMMA (default) -- Malamud et al. (2004), ESPL 29, 687-711,
//      Eq. 3. Three parameters. Power-law tail with exponent -(rho+1) AND an
//      exponential rollover at small areas. The rollover is included
//      deliberately: Malamud argues it is real rather than a mapping artefact
//      (their inventories are complete well below the modal area of ~400 m2),
//      and HyLands independently produces one from mechanics.
//
//   2. TRUNCATED PARETO -- Yanites et al. (2009), JGR-ES 114, F01007, Eq. 1.
//      Pure power law between hard bounds, no rollover. Kept because it is what
//      Niemi et al. (2005) and Yanites use, so it is needed to reproduce their
//      published cosmogenic results.
//
// EXPONENT CONVENTION -- a live trap. Niemi's b = 1.18 is a CUMULATIVE
// exponent; Yanites' and Valagussa's are PDF exponents, with
// alpha_pdf = b_cumulative + 1. Everything in this file is a PDF exponent, and
// the parameter is named pareto_alpha_pdf to say so.
//
// AREA-TO-VOLUME: V = alpha * A^gamma, with A in m^2 and V in m^3, so alpha
// carries units L^(3-2*gamma) with L in metres. Defaults are from Larsen,
// Montgomery & Korup (2010), Nature Geoscience 3, 247-251, SUPPLEMENTARY
// TABLE S1 -- which separates SCAR from DEPOSIT geometry. This module excavates
// scars, so the scar rows are the material-specific presets. The commonly
// quoted gamma_soil = 1.145 / gamma_bedrock = 1.35 are the COMBINED scar +
// deposit fits and are the wrong rows for this purpose; they are provided for
// comparison with literature that quotes them.
//
// The unit convention is cross-checked two ways in the tests, because getting
// alpha's units wrong is silent: Simonett (1967) in km units converts to
// log10(alpha) = -0.828 against Table S1's -0.82, and Hovius eps = 0.05 gives
// -1.301 against Table S1's footnoted -1.30.
// -----------------------------------------------------------------------------

#ifndef LSDTT_LANDSLIDE_SIZE_CPU_HPP
#define LSDTT_LANDSLIDE_SIZE_CPU_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Distribution family
// -----------------------------------------------------------------------------
enum class LandslideSizeModel
{
    InverseGamma,     // Malamud et al. (2004) Eq. 3 -- default
    TruncatedPareto   // Yanites et al. (2009) Eq. 1
};

std::string to_string(LandslideSizeModel model);

// Parse a parameter-file spelling ("inverse_gamma", "truncated_pareto").
// Returns false on an unknown name rather than guessing.
bool parse_landslide_size_model(const std::string& name, LandslideSizeModel& out);


// -----------------------------------------------------------------------------
// Area-to-volume presets, V = alpha * A^gamma  (A in m^2, V in m^3)
//
// From Larsen et al. (2010) Supplementary Table S1 unless noted. alpha is
// 10^(log10 alpha) as tabulated there.
// -----------------------------------------------------------------------------
struct AreaVolumeScaling
{
    double alpha = 0.1459;   // larsen_global
    double gamma = 1.332;

    // DENSITY OF THE MATERIAL THE SCAR EXCAVATES, kg/m3.
    //
    // Larsen's V = alpha A^gamma is fitted to MEASURED scar volumes, so
    // d = V/A is a true geometric thickness of whatever failed -- and for a
    // shallow landslide that is soil, not rock. Converting it to a mass depth
    // with 2650 overstates the shielding by roughly a factor of two, and since
    // CRN production falls off exponentially with mass depth, a factor of two
    // there is a several-fold error in the concentration of the material
    // delivered.
    //
    // So the density travels WITH the preset: choosing larsen_soil_scar is
    // already a statement that the failures are in soil.
    double scar_density_kg_m3 = 2650.0;

    // False for presets fitted across BOTH materials (larsen_global) or with no
    // material attached (hovius_*). Those have no defensible density, and the
    // 2650 above is a placeholder rather than a value. Callers should say so
    // rather than pick silently -- see area_volume_density_note().
    bool density_is_material_specific = false;

    // Table S1 rows. Prefer the SCAR presets: this module excavates scars.
    static AreaVolumeScaling larsen_global();        // g 1.332, a 0.1459, n=4231
    static AreaVolumeScaling larsen_soil_scar();     // g 1.262, a 0.2244, n=1617
    static AreaVolumeScaling larsen_bedrock_scar();  // g 1.410, a 0.2344, n=168
    static AreaVolumeScaling larsen_soil_all();      // g 1.145, a 0.3631  (scar+deposit)
    static AreaVolumeScaling larsen_bedrock_all();   // g 1.350, a 0.1862  (scar+deposit)

    // d = 0.05*sqrt(A) as used by Hovius (1997), Niemi (2005), Yanites (2009)
    // and Malamud (2004). Present so those results can be reproduced. Larsen
    // shows gamma = 1.5 overestimates total volume by >= 2x versus gamma = 1.4,
    // growing with maximum landslide area -- and the slab form is itself twice
    // the parabolic form at equal eps, so a careless choice can be 4x out.
    static AreaVolumeScaling hovius_slab();          // g 1.5, a 0.05
    static AreaVolumeScaling hovius_parabolic();     // g 1.5, a 0.025

    static bool from_name(const std::string& name, AreaVolumeScaling& out);

    double volume_m3(double area_m2) const;

    // Mass depth (g/cm2) of a scar of this geometric depth, using the preset's
    // own material density. This is the ONLY correct way to turn a scar
    // thickness into the quantity the CRN physics works in.
    double scar_effective_depth_g_cm2(double depth_m) const;

    // Geometric thickness of ROCK carrying the same mass as depth_m of this
    // preset's material. Needed because lsdtt_crn_column carries a single
    // density for a whole erosion history: a soil removal is handed to it as
    // the rock-equivalent thickness, which removes the right mass without the
    // column needing to know two densities.
    double rock_equivalent_depth_m(double depth_m, double rock_density_kg_m3) const;
};

// Empty when the preset carries a defensible material density; otherwise a
// sentence naming the problem, suitable for printing as a warning. A caller
// that has been given an explicit density by the user should not print it.
[[nodiscard]] std::string area_volume_density_note(const AreaVolumeScaling& av);


// -----------------------------------------------------------------------------
// Distribution parameters
// -----------------------------------------------------------------------------
struct LandslideSizeParams
{
    LandslideSizeModel model = LandslideSizeModel::InverseGamma;

    // --- Inverse gamma (Malamud et al. 2004 Eq. 3), areas in m^2 ---
    // Maximum-likelihood fit to three substantially complete inventories:
    // Umbria snowmelt, Northridge earthquake, Guatemala rainfall. r^2 = 0.965.
    // NOTE s is NEGATIVE, and the density is normalised over s <= A < infinity.
    double rho     = 1.40;      // tail exponent is -(rho+1) = -2.40
    double a_param = 1280.0;    // 1.28e-3 km^2; sets the modal position
    double s_param = -132.0;    // -1.32e-4 km^2; sets the small-area rollover

    // --- Truncated Pareto (Yanites et al. 2009 Eq. 1) ---
    // PDF exponent, NOT cumulative. Valagussa et al. (2019) measure 1.75-3.11
    // across lithologies for coseismic inventories; Niemi's b = 1.18 cumulative
    // is alpha_pdf = 2.18.
    double pareto_alpha_pdf = 2.3;

    // --- Bounds (both families) ---
    int    n_min_pixels = 1;    // smallest landslide, in pixels
    double area_max_m2  = 0.0;  // must be > 0; capped by the channel threshold

    AreaVolumeScaling av{};

    // Validate. Returns false and fills `error` rather than throwing, so a CLI
    // can report cleanly.
    bool validate(double cell_size_m, std::string& error) const;
};


// -----------------------------------------------------------------------------
// The discrete sampler
//
// Builds the pmf over pixel counts [n_min, n_max] at construction and samples
// by inverse CDF. Construction is O(n_max); sampling is O(log n_max).
//
// n_max is floor(area_max_m2 / cell_area), so the caller controls it through
// area_max_m2 -- which the placement layer pins to the channel initiation
// threshold, since a landslide cannot be larger than the hillslope it comes
// from.
// -----------------------------------------------------------------------------
class LandslideSizeSampler
{
  public:
    // Throws std::invalid_argument if params.validate() fails.
    LandslideSizeSampler(const LandslideSizeParams& params, double cell_size_m);

    // --- Sampling ---
    int    sample_pixels(double u01) const;          // u01 in [0,1)
    double sample_area_m2(double u01) const;

    // --- Distribution properties, computed from the TRUNCATED pmf ---
    // mean_volume_m3 is the number that converts a volume budget into a
    // landslide count, and it depends strongly on n_max: the volume density of
    // the Malamud tail goes as A^(gamma - rho - 1) ~ A^-1.068, so each
    // logarithmic decade of area contributes roughly equally to total volume.
    // Changing the channel threshold changes the landslide count for a fixed
    // budget by tens of per cent. That is physical, not a defect -- Niemi:
    // "large, but infrequent landslides dominate the overall sediment flux".
    double mean_area_m2()   const noexcept { return mean_area_m2_; }
    double mean_volume_m3() const noexcept { return mean_volume_m3_; }

    // Volume of a landslide of n pixels, through this sampler's own
    // area-volume scaling. Callers must use this rather than reconstructing
    // alpha*A^gamma themselves: V is non-linear in A, so a mean-area shortcut
    // is wrong, and duplicating the scaling invites the two copies to drift.
    double volume_of_pixels(int n_pixels) const
    {
        return params_.av.volume_m3(static_cast<double>(n_pixels) * cell_area_m2_);
    }

    const AreaVolumeScaling& area_volume_scaling() const noexcept { return params_.av; }

    int    n_min_pixels() const noexcept { return n_min_; }
    int    n_max_pixels() const noexcept { return n_max_; }
    double cell_area_m2() const noexcept { return cell_area_m2_; }

    // pmf / cdf access, for tests and for reporting the realised distribution.
    double pmf(int n_pixels) const;
    const std::vector<double>& cdf() const noexcept { return cdf_; }

    // Fraction of the UNTRUNCATED distribution lying above area_max, by number
    // and by volume. Expect the volume fraction to be large (~59% at a 1000-px
    // threshold on a 30 m DEM). Reported so the
    // truncation is visible, never silently absorbed.
    double truncated_fraction_number() const noexcept { return trunc_frac_n_; }
    double truncated_fraction_volume() const noexcept { return trunc_frac_v_; }

    // True when the cell is too coarse to resolve the modal area, so the
    // rollover cannot be represented and the sampled distribution is
    // effectively a truncated power law whatever the caller selected. Callers
    // should warn. Always false for TruncatedPareto, which has no rollover.
    bool rollover_unresolved() const noexcept { return rollover_unresolved_; }

  private:
    LandslideSizeParams params_;
    double cell_area_m2_   = 0.0;
    int    n_min_          = 1;
    int    n_max_          = 1;
    std::vector<double> cdf_;          // cumulative, size n_max_-n_min_+1
    double mean_area_m2_   = 0.0;
    double mean_volume_m3_ = 0.0;
    double trunc_frac_n_   = 0.0;
    double trunc_frac_v_   = 0.0;
    bool   rollover_unresolved_ = false;
};


// -----------------------------------------------------------------------------
// Free functions: the continuous densities, exposed for testing against the
// published equations and for the truncation diagnostics.
// -----------------------------------------------------------------------------

// Malamud et al. (2004) Eq. 3. Returns 0 for A <= s.
double inverse_gamma_pdf(double area_m2, double rho, double a_param, double s_param);

// P(A > area) for the inverse gamma, via the regularised incomplete gamma.
double inverse_gamma_survival(double area_m2, double rho, double a_param, double s_param);

// Malamud et al. (2004) Eq. 14: mean of the UNTRUNCATED distribution,
// a/(rho-1) + s. Requires rho > 1; returns infinity otherwise. With the
// published parameters this is 3070 m^2, independent of event size.
double inverse_gamma_mean_area_m2(double rho, double a_param, double s_param);

// Yanites et al. (2009) Eq. 1, normalised over [a_min, a_max].
double truncated_pareto_pdf(double area_m2, double alpha_pdf,
                            double a_min_m2, double a_max_m2);

// Closed-form inverse transform for the truncated Pareto. Handles alpha -> 1.
double truncated_pareto_sample(double u01, double alpha_pdf,
                               double a_min_m2, double a_max_m2);

// Closed-form mean of the truncated Pareto. Handles alpha -> 1 and alpha -> 2.
double truncated_pareto_mean_area_m2(double alpha_pdf,
                                     double a_min_m2, double a_max_m2);

}  // namespace lsdtt

#endif  // LSDTT_LANDSLIDE_SIZE_CPU_HPP
