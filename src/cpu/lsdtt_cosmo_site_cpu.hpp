// lsdtt_cosmo_site_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Everything between "I have a DEM and some parameters" and "I have a
// production raster and three shielding rasters".
//
// WHY THIS EXISTS. That stretch -- loading crn_parameters.nc, resolving the
// atmosphere and scaling scheme, building the pressure raster and the CRONUS
// emulator inputs, computing or loading production, and resolving shielding
// from either a raster or a scalar -- was written out separately in
// lsdtt-cosmo-prediction and lsdtt-cosmo-tool, roughly 150 lines each. Giving
// lsdtt-cosmo-landslides real production rates would have made it a third copy.
//
// CAPABILITIES, NOT A UNION. The two CLIs are NOT the same shape and unifying
// them naively would have been a silent feature change:
//
//   lsdtt-cosmo-tool has NO shielding-raster parameters and NO
//   production-from-disk path. It declares only the scalars and always computes.
//
// Handing it those parameters as a side effect of a refactor would have added
// an untested code path to an inversion tool while passing every numeric check,
// because the prefixes default to NULL. So a caller states which capabilities
// it has, and anything outside that set is not merely unused -- it is not
// available. Changing a CLI's capability set is a feature decision with its own
// review, never a consequence of sharing code.
//
// VALIDATION IS SEPARATE FROM CONSTRUCTION, and callable earlier. A bad
// scaling_scheme should fail before a DEM is read and a step change generated,
// not after. validate_cosmo_site_config() does no I/O and is safe to call
// immediately after parsing.
//
// PRODUCTION IS RE-COMPUTABLE, not just readable. lsdtt-cosmo-tool's
// self-consistent-age mode rebuilds production over ONE BASIN's pixels at each
// iterate, via the eval_only overload; its own comment records that evaluating
// the whole DEM every iterate "was the dominant cost of this mode". A component
// that owned only a finished raster could not express that, and the tool would
// have had to bypass it -- putting the duplication straight back.
// -----------------------------------------------------------------------------

#ifndef LSDTT_COSMO_SITE_CPU_HPP
#define LSDTT_COSMO_SITE_CPU_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_cosmo_raster_cpu.hpp"
#include "lsdtt_crn_parameters.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// What a given CLI is allowed to do. See the header note: this is a capability
// statement, not a convenience filter.
// -----------------------------------------------------------------------------
enum class CosmoSiteCaps : unsigned
{
    None               = 0u,
    Production         = 1u << 0,  // scaling scheme, atmosphere, crn_nc, sample age
    ProductionFromDisk = 1u << 1,  // load_production_raster + prefix
    ShieldingScalars   = 1u << 2,  // topo_shield_value, snow/self eff thickness
    ShieldingRasters   = 1u << 3,  // the three *_raster_prefix parameters

    // Build the pressure raster even when no muon scheme needs it.
    // lsdtt-cosmo-tool reports a basin-mean effective pressure in its output
    // CSV and the CRONUS comparison workflow requires it, so for that tool the
    // pressure raster is a RESULT, not an intermediate. Where it is only an
    // intermediate, a failure to build it is fatal; where it is a reported
    // value, a failure degrades to nan with a warning. Both behaviours were
    // already in the tree and neither is wrong.
    PressureAlways     = 1u << 4
};

constexpr CosmoSiteCaps operator|(CosmoSiteCaps a, CosmoSiteCaps b)
{
    return static_cast<CosmoSiteCaps>(static_cast<unsigned>(a) |
                                      static_cast<unsigned>(b));
}
constexpr bool has_cap(CosmoSiteCaps set, CosmoSiteCaps c)
{
    return (static_cast<unsigned>(set) & static_cast<unsigned>(c)) != 0u;
}


// -----------------------------------------------------------------------------
// The parameters, exactly as the CLIs already declare them. Field names match
// the parameter names so a reader can map one to the other without a table.
// -----------------------------------------------------------------------------
struct CosmoSiteConfig
{
    std::string nuclide        = "Be10";
    std::string muon_scheme    = "BraucherBorchers";
    std::string scaling_scheme = "stone";
    std::string atmosphere     = "ERA40";
    std::string crn_nc_file;                 // empty -> the bundled nc
    double      sample_age_yr  = 0.0;

    // ProductionFromDisk
    bool        load_production_raster   = false;
    std::string production_raster_prefix = "NULL";

    // ShieldingRasters (each "NULL" -> use the matching scalar)
    std::string topographic_shielding_raster_prefix = "NULL";
    std::string snow_shielding_raster_prefix        = "NULL";
    std::string self_shielding_raster_prefix        = "NULL";

    // ShieldingScalars
    double topo_shield_value         = 1.0;
    double snow_shield_eff_thickness = 0.0;
    double self_shield_eff_thickness = 0.0;

    // Directory that raster prefixes are resolved against.
    std::string read_path;
};


// -----------------------------------------------------------------------------
// Parameter validation. NO I/O, so it is safe -- and intended -- to call
// immediately after parsing, before a DEM is read or any expensive scenario is
// generated.
//
// Throws std::runtime_error whose what() is the complete diagnostic WITHOUT a
// "FATAL: " prefix; the caller adds that, matching every other CLI message.
//
// The lsdn-grids check is NOT here: it needs crn_parameters.nc loaded, so it
// lives in build_cosmo_site().
// -----------------------------------------------------------------------------
void validate_cosmo_site_config(const CosmoSiteConfig& cfg, CosmoSiteCaps caps);


// -----------------------------------------------------------------------------
// The built site.
//
// Owns its rasters. `emulator_inputs()` is built on demand rather than stored,
// because CronusEmulatorInputs holds a RAW POINTER to the pressure raster and a
// stored copy would dangle the moment the site were moved.
// -----------------------------------------------------------------------------
class CosmoSite
{
public:
    const lsdtt_raster_cpu& production()            const { return production_; }
    const lsdtt_raster_cpu& topographic_shielding() const { return topo_shield_; }
    const lsdtt_raster_cpu& snow_shielding()        const { return snow_shield_; }
    const lsdtt_raster_cpu& self_shielding()        const { return self_shield_; }

    const CrnNcData&  crn_nc()            const { return crn_nc_; }
    const std::string& crn_nc_path()      const { return crn_nc_path_; }
    AtmSource         atmosphere_source() const { return atm_source_; }
    ScalingScheme     scaling_scheme()    const { return scaling_scheme_; }
    bool              uses_lsdn()         const
    { return scaling_scheme_ == ScalingScheme::LSDn; }

    // "NCEP", "ERA40" or "standard" -- what was ACTUALLY used after the
    // documented absent-grid fallback, not what was asked for.
    const std::string& atmosphere_report() const { return atm_report_; }

    // Empty unless an emulator muon scheme is in play or PressureAlways was
    // requested. has_pressure() is false if it was wanted but could not be
    // built -- which is a warning, not a failure, when nothing depends on it.
    const lsdtt_raster_cpu& pressure()     const { return pressure_; }
    bool                    has_pressure() const { return has_pressure_; }

    CronusEmulatorInputs emulator_inputs() const
    {
        CronusEmulatorInputs in;
        if (has_pressure_) in.pressure_raster = &pressure_;
        in.lsdn_slhl_reference = lsdn_slhl_reference_;
        return in;
    }

    // Rebuild production at a different age, optionally over a SUBSET of pixels
    // (flat row*n_cols + col indices). The subset form is what makes
    // lsdtt-cosmo-tool's self-consistent-age iteration affordable; see the
    // header note. `dem` must be the same raster the site was built from.
    [[nodiscard]] lsdtt_raster_cpu recompute_production(
        const lsdtt_raster_cpu& dem,
        double                  age_yr,
        const std::vector<std::int64_t>* eval_only = nullptr) const;

    [[nodiscard]] std::string summary() const;

private:
    friend CosmoSite build_cosmo_site(const lsdtt_raster_cpu&,
                                      const CosmoSiteConfig&, CosmoSiteCaps);

    lsdtt_raster_cpu production_, topo_shield_, snow_shield_, self_shield_;
    lsdtt_raster_cpu pressure_;
    CrnNcData        crn_nc_;
    std::string      crn_nc_path_;
    std::string      atm_report_;
    AtmSource        atm_source_     = AtmSource::ERA40;
    ScalingScheme    scaling_scheme_ = ScalingScheme::Stone;
    double           lsdn_slhl_reference_ = 0.0;
    bool             needs_pressure_ = false;
    bool             has_pressure_   = false;
    bool             production_from_disk_ = false;
    std::string      nuclide_;
};


// Build the site. Loads crn_parameters.nc, resolves the atmosphere, builds the
// pressure raster and emulator inputs where the muon scheme needs them,
// computes or loads production, and resolves the three shielding fields.
//
// Throws std::runtime_error (message without a "FATAL: " prefix) on any fatal
// condition: missing LSDn grids, an unreadable raster, a production failure.
[[nodiscard]] CosmoSite build_cosmo_site(const lsdtt_raster_cpu& dem,
                                         const CosmoSiteConfig&  cfg,
                                         CosmoSiteCaps           caps);


// -----------------------------------------------------------------------------
// load_raster_or_uniform
//
// "A raster from disk if a prefix was given, otherwise a uniform field." Four
// call sites in lsdtt-cosmo-prediction use this shape -- the three shielding
// fields and the background erosion rate -- and it was a lambda.
//
// MaskedFill::Value reproduces that lambda: the uniform branch filled the whole
// array including pixels outside the mask. See lsdtt_raster_cpu.hpp for why
// that choice is a parameter rather than a convention.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu load_raster_or_uniform(
    const lsdtt_raster_cpu& reference,
    const std::string&      prefix,      // "NULL" -> uniform
    const std::string&      read_path,
    double                  uniform_value);

}  // namespace lsdtt

#endif  // LSDTT_COSMO_SITE_CPU_HPP
