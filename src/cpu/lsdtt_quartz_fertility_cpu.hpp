// lsdtt_quartz_fertility_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Quartz fertility, and the ice mask that shares its machinery.
//
// WHY THIS EXISTS. A basin-averaged cosmogenic sample is not a mixture of ROCK,
// it is a mixture of the QUARTZ separated from that rock. Two pixels eroding at
// identical rates contribute unequally if their rock carries different quartz,
// and in lithologically complex terrain some source rocks carry almost none --
// limestone, basalt, many volcanics and mafics. Route concentrations pixel by
// pixel without accounting for that and the answer is simply WRONG: a
// quartz-free unit eroding fast contributes heavily to the sediment flux and
// nothing at all to the sample.
//
// This is a correctness issue in mixed-lithology catchments, not a refinement.
// It is also the reason the module is lithology-first: chemical weathering also
// enriches quartz in regolith relative to fresh bedrock, and the depth-dependent
// option below can express that, but it is a second-order effect and not the
// reason the module exists.
//
// WHAT IT PRODUCES. One field: a mixing WEIGHT in kg of quartz per square metre
// per year,
//
//     w(x) = mass_flux(x) * f_q(x)
//
// which is handed to the existing accumulation kernel in place of the erosion
// rate. Nothing in the CRN physics changes; only the weights do. With
// f_q == 1 everywhere the weight is proportional to the mass flux and the
// accumulation is unchanged, which is why the whole feature can default to off
// and be a no-op.
//
// THE ICE MASK. An ice-covered pixel differs from a bare one in exactly three
// ways, and lsdtt3 already has a field for each:
//
//   1. cosmic rays are attenuated by the ice column  -> SNOW SHIELDING, which is
//      already an effective depth in g/cm2 (ice at 917 kg/m3: 10 m = 917 g/cm2,
//      i.e. essentially total);
//   2. the erosion rate differs, possibly zero       -> the EROSION RATE raster;
//   3. little or no quartz reaches the sample        -> QUARTZ FERTILITY.
//
// So an ice mask is a helper that writes three existing rasters, NOT a module of
// its own. A module whose entire content is "set three fields" earns nothing and
// creates a place for the three to drift apart.
//
// LIMITATION, STATED RATHER THAN HIDDEN: this models PRESENT-DAY ice only.
// Concentrations depend on the shielding HISTORY, and real glaciated landscapes
// advanced and retreated repeatedly. lsdtt3 can express a per-pixel erosion
// history but not a per-pixel time-varying shielding, so a glacial-interglacial
// cycle cannot be represented at all. Concentrations under formerly glaciated
// ground will be overestimated.
// -----------------------------------------------------------------------------

#ifndef LSDTT_QUARTZ_FERTILITY_CPU_HPP
#define LSDTT_QUARTZ_FERTILITY_CPU_HPP

#include <cstdint>
#include <string>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// QuartzFertilityConfig
//
// Fertility is a mass fraction, 0 to 1. Zero is legitimate and meaningful: the
// pixel contributes sediment but no sample.
// -----------------------------------------------------------------------------
struct QuartzFertilityConfig
{
    // Uniform fertility, used where no lithology raster is supplied. 1.0 makes
    // the whole feature a no-op, which is the default.
    double uniform_fertility = 1.0;

    // Depth-dependent fertility, enabled when regolith_depth_m > 0.
    //
    // Chemical weathering removes soluble minerals and leaves quartz behind, so
    // regolith is enriched relative to the fresh bedrock a deep landslide scar
    // exposes. If that holds, landslide debris is UNDER-represented in the
    // quartz separate relative to its mass, and the dilution a landslide causes
    // is weaker than a mass-weighted calculation predicts. That is a testable
    // claim, and this is the switch that tests it.
    //
    // Negative values mean "same as the surface fertility", i.e. no contrast.
    double regolith_fertility = -1.0;
    double bedrock_fertility  = -1.0;
    double regolith_depth_m   = 0.0;

    [[nodiscard]] bool depth_dependent() const
    {
        return regolith_depth_m > 0.0
            && regolith_fertility >= 0.0 && bedrock_fertility >= 0.0;
    }
};

// -----------------------------------------------------------------------------
// IceMaskConfig
// -----------------------------------------------------------------------------
struct IceMaskConfig
{
    // Uniform ice thickness in metres, used where no thickness raster is given.
    double ice_thickness_m = 50.0;

    // Density of ice. 917 kg/m3 for solid ice; lower for firn.
    double ice_density_kg_m3 = 917.0;

    // Effective erosion rate beneath the ice, g/cm2/yr. Zero is the cold-based
    // idealisation and is numerically safe: the steady-state kernel forms its
    // denominator as (eps + Gamma*lambda), never divides by eps, and returns
    // the saturation concentration P/lambda.
    double erate_under_ice_g_cm2_yr = 0.0;

    // Quartz fertility beneath the ice. Zero means the ice contributes no
    // material to the sample at all, which is the usual assumption when the
    // erosion rate under it is also zero.
    double quartz_fertility_under_ice = 0.0;
};

// -----------------------------------------------------------------------------
// QuartzWeightSummary
//
// Reported, not silently absorbed. A catchment whose upslope area is entirely
// quartz-free yields no sample -- a meaningful and, in volcanic terrain, common
// answer that must not be returned as nodata without explanation.
// -----------------------------------------------------------------------------
struct QuartzWeightSummary
{
    std::int64_t n_valid            = 0;
    std::int64_t n_zero_fertility   = 0;   // contributes sediment, no quartz
    std::int64_t n_zero_weight      = 0;   // contributes nothing at all
    double total_quartz_kg_yr       = 0.0;
    double mean_fertility           = 0.0; // over valid pixels, unweighted

    [[nodiscard]] bool yields_no_quartz() const { return total_quartz_kg_yr <= 0.0; }
    [[nodiscard]] std::string summary() const;
};

// -----------------------------------------------------------------------------
// make_quartz_fertility_raster
//
// The SURFACE fertility field: what the material eroding off each pixel carries.
// lithology_raster may be null, in which case the uniform value is used
// everywhere. Values outside [0, 1] are clamped, and the number clamped is
// reported through n_clamped so a mis-scaled lithology map (percentages rather
// than fractions, say) is caught rather than silently truncated.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu make_quartz_fertility_raster(
    const lsdtt_raster_meta&     meta,
    const QuartzFertilityConfig& config,
    const lsdtt_raster_cpu*      lithology_raster,
    std::int64_t*                n_clamped = nullptr);

// -----------------------------------------------------------------------------
// delivered_quartz_fertility
//
// Fertility of the material a landslide of geometric depth depth_m delivers:
// the depth average over the excavated slab of a two-layer profile, regolith
// above regolith_depth_m and bedrock below.
//
//   f(d) = [ f_reg * min(d, D) + f_bed * max(0, d - D) ] / d
//
// Returns surface_fertility unchanged when the depth-dependent option is off,
// so the two-layer model is strictly opt-in.
// -----------------------------------------------------------------------------
[[nodiscard]] double delivered_quartz_fertility(
    const QuartzFertilityConfig& config,
    double                       surface_fertility,
    double                       depth_m);

// -----------------------------------------------------------------------------
// make_quartz_weight_raster
//
// weight = mass_flux * fertility, kg of quartz m^-2 yr^-1. This is what the
// accumulation kernel should mix on.
//
// When excavation_depth_m and the depth-dependent option are both supplied, a
// failed pixel's fertility is the DELIVERED fertility of its slab rather than
// its surface value -- that is the only place the two-layer model bites.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu make_quartz_weight_raster(
    const lsdtt_raster_cpu&      mass_flux_kg_m2_yr,
    const lsdtt_raster_cpu&      fertility,
    const QuartzFertilityConfig& config,
    const lsdtt_raster_cpu*      excavation_depth_m,
    QuartzWeightSummary*         summary = nullptr);

// -----------------------------------------------------------------------------
// apply_ice_mask
//
// Writes the three fields ice actually changes. ice_mask is treated as a
// boolean: any valid, non-zero, non-nodata value is ice.
//
// ice_thickness_raster may be null, in which case the uniform thickness is used.
// Any of the three output rasters may be empty (n_rows == 0) and is then skipped,
// so a caller that does not carry, say, a fertility field need not invent one.
//
// Returns the number of pixels marked as ice.
// -----------------------------------------------------------------------------
std::int64_t apply_ice_mask(
    const lsdtt_raster_cpu&  ice_mask,
    const lsdtt_raster_cpu*  ice_thickness_raster,
    const IceMaskConfig&     config,
    lsdtt_raster_cpu&        snow_shield_g_cm2,
    lsdtt_raster_cpu&        erosion_rate_g_cm2_yr,
    lsdtt_raster_cpu&        quartz_fertility);

// -----------------------------------------------------------------------------
// ice_effective_depth_g_cm2
//
// Mass depth of an ice column. Separate and named because the conversion is the
// same 0.1 factor used everywhere else and getting it wrong is silent.
// -----------------------------------------------------------------------------
[[nodiscard]] double ice_effective_depth_g_cm2(double thickness_m,
                                               double ice_density_kg_m3);

} // namespace lsdtt

#endif // LSDTT_QUARTZ_FERTILITY_CPU_HPP
