// lsdtt_landcover_scenario_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Land-cover scenarios: one class raster, three derived fields
//
// Phase 9. A scenario RECLASSIFIES pixels -- grazing becomes woodland inside an
// intervention footprint -- and three things follow from that one change:
//
//   n     hillslope and floodplain roughness, from the class table
//   LAI   and hence interception capacity and free throughfall
//   T0    transmissivity, by RE-KEYING the measured soil lookup
//
// Changing one without the others is not a scenario, it is a bug: planting
// woodland raises roughness, raises LAI, and raises conductivity through
// reduced compaction (Marshall et al. 2009, 2014). They move together because
// they are consequences of the same physical change.
//
// WHY TRANSMISSIVITY IS RE-KEYED RATHER THAN MULTIPLIED
// -----------------------------------------------------
// This is the substantive difference from SD-TOPMODEL, and the reason it is
// possible here at all.
//
// Zhu et al. (2025) represent woodland planting as `Ks x 1.5`: a calibration
// choice, defensible but arbitrary, and the same multiplier wherever it is
// applied. lsdtt3 has the measured alternative. `lsdtt_sti_cpu`'s lookup is
// keyed on (soil series, land use) and its `ln_kd` field IS ln(T0), so
// reclassifying a pixel re-keys that lookup and returns the MEASURED
// transmissivity for THAT SOIL under the new land use.
//
// The consequences are real, not cosmetic. A multiplier gives every soil the
// same proportional response to planting; re-keying gives a peaty gley and a
// freely drained brown earth the responses their profiles actually show, which
// may differ in size and can differ in sign. And no calibration parameter is
// introduced, because nothing was fitted.
//
// Where the measured pair is absent the STI fallback chain applies and the
// substitution is COUNTED, so a scenario resting on fallbacks is visible rather
// than silently plausible.
//
// WHAT A SCENARIO MUST NOT TOUCH
// ------------------------------
// Channel roughness. That is grain-size controlled (see
// lsdtt_hydro_roughness_cpu.hpp); planting a field does not change the bed
// material. The roughness this module produces is the LAND COVER field, which
// compute_roughness_field() applies to hillslope and floodplain cells only.
//
// References
// ----------
//   Marshall M.R. et al. (2009, 2014) -- grazing to woodland raises Ks
//   Zhu et al. (2025) Hydrol. Process. 39:e70122 -- the multiplier approach
//   Boorman D.B., Hollis J.M. & Lilly A. (1995) -- HOST, the usual UK source
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_sti_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// One land-cover class and the properties it implies.
// ---------------------------------------------------------------------------
struct landcover_class
{
    std::string name;
    int         code = 0;
    double      lai  = 0.0;
    double      manning_n = 0.0;

    // The land-use key used against the soil lookup: "CULT", "SEMI", "OTHER"
    // in the HOST/SSKIB vocabulary. This is what makes the transmissivity
    // measured rather than assumed.
    std::string soil_land_use;
};

// ---------------------------------------------------------------------------
// The class table, parsed from a parameter string:
//
//   "grazing:code=1,lai=1.5,n=0.30,soil=CULT; woodland:code=2,lai=5.0,n=0.60,soil=SEMI"
//
// Throws std::invalid_argument on a malformed entry, a duplicate code or name,
// a negative LAI or roughness, or a missing required field.
// ---------------------------------------------------------------------------
struct landcover_table
{
    std::vector<landcover_class> classes;

    [[nodiscard]] static landcover_table parse(const std::string& spec);

    [[nodiscard]] const landcover_class* by_code(int code) const noexcept;
    [[nodiscard]] const landcover_class* by_name(const std::string& name) const noexcept;
    [[nodiscard]] std::size_t size() const noexcept { return classes.size(); }
};

// ---------------------------------------------------------------------------
// A scenario: a set of reclassification rules, applied where an intervention
// mask is set.
//
//   "grazing->woodland; arable->woodland"
//
// A rule with source "*" matches any class. Rules are applied in order and a
// pixel is reclassified at most once, so an earlier rule wins -- which makes
// "grazing->woodland; *->hedgerow" mean what it reads like.
//
// Throws std::invalid_argument on a malformed rule or an unknown class name.
// ---------------------------------------------------------------------------
struct landcover_scenario_rule
{
    std::string from;   // "*" matches any
    std::string to;
};

struct landcover_scenario
{
    std::vector<landcover_scenario_rule> rules;

    [[nodiscard]] static landcover_scenario parse(const std::string& spec);
    [[nodiscard]] bool empty() const noexcept { return rules.empty(); }
};

// ---------------------------------------------------------------------------
// Applying a scenario.
// ---------------------------------------------------------------------------
struct scenario_application_result
{
    lsdtt_raster_cpu landcover;      // the reclassified class raster

    std::int64_t n_considered   = 0; // valid cells inside the intervention mask
    std::int64_t n_reclassified = 0;
    std::int64_t n_unmatched    = 0; // in the mask but matching no rule

    // Reclassification counts, keyed "from->to".
    std::unordered_map<std::string, std::int64_t> transitions;
};

// `intervention_mask` may be null, in which case the scenario applies
// everywhere. Cells where it is zero or nodata keep their class.
//
// Throws std::invalid_argument for a grid mismatch, or a class code in the
// raster that the table does not define.
[[nodiscard]]
scenario_application_result apply_landcover_scenario(
    const lsdtt_raster_cpu&   landcover,
    const lsdtt_raster_cpu*   intervention_mask,
    const landcover_table&    table,
    const landcover_scenario& scenario
);

// ---------------------------------------------------------------------------
// Deriving the three fields from a (possibly reclassified) class raster.
// ---------------------------------------------------------------------------
struct landcover_field_inputs
{
    // Required.
    const lsdtt_raster_cpu* landcover = nullptr;

    // Required only for the transmissivity field. Soil series codes, with the
    // string lookup that came out of rasterise_polygon_string_field().
    const lsdtt_raster_cpu* soil_series_codes = nullptr;
    const std::vector<std::pair<std::uint32_t, std::string>>* soil_series_lookup = nullptr;
    const sti_lnkd_table*   lnkd_table = nullptr;
};

struct landcover_field_result
{
    lsdtt_raster_cpu lai;
    lsdtt_raster_cpu roughness;       // land-cover n: hillslope and floodplain only
    lsdtt_raster_cpu ln_t0;           // ln(T0) in m^2/h; empty when no soil inputs

    std::int64_t n_cells               = 0;
    std::int64_t n_t0_direct           = 0;  // exact (series, land use) match
    std::int64_t n_t0_fallback         = 0;  // a related land use supplied it
    std::int64_t n_t0_missing          = 0;  // nothing did; left nodata
};

// Derives LAI and roughness always, and ln(T0) when the soil inputs are
// supplied. The ln(T0) values are converted from the STI table's `ln_kd` with
// sti_ln_kd_to_ln_t0_m2_per_hour(), so they are ready to use as
// topmodel_parameters::ln_t0 without a further factor.
//
// Throws std::invalid_argument for a missing required input, a grid mismatch,
// or a class code the table does not define.
[[nodiscard]]
landcover_field_result derive_landcover_fields(
    const landcover_field_inputs& inputs,
    const landcover_table&        table
);

} // namespace lsdtt
