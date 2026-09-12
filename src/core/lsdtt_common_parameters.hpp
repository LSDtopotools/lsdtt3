// lsdtt_common_parameters.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_common_parameters
//
// WHY THIS FILE EXISTS
//
// Until now there was no shared declaration layer between the CLIs. Every
// `cli/lsdtt-*.cpp` hand-typed its own parameter block into its own
// `*_default_map`, so a concept used by twenty tools was spelled twenty times
// and drifted. The channel-initiation threshold acquired five different names;
// the same fill flag acquired two defaults; the same output-format parameter
// acquired three help texts. None of that is a coding error - it is the
// inevitable consequence of there being no single place to put the name.
//
// This file is that place. Each function below declares a COHERENT GROUP of
// parameters into the CALLER'S default maps and the CALLER'S help registry, so
// that after adoption a shared parameter name, type, default and help text are
// written down EXACTLY ONCE in the codebase. Renaming a shared parameter then
// becomes a one-line edit here plus one deprecated-alias entry per CLI, rather
// than a twenty-CLI sweep. That is what makes the naming standard in
// `docs/cli_parameter_conventions.md` self-enforcing rather than aspirational.
//
// This is Phase 0a of `docs/plan_docs/cli_consistency_remediation_plan.md`.
// Phase 0 changes no CLI and no behaviour: these functions are purely additive
// and nothing calls them yet. Adoption by the 38 CLIs is Phase 2+ work.
//
// USAGE
//
//   std::map<std::string,bool>        bool_default_map;
//   std::map<std::string,int>         int_default_map;
//   std::map<std::string,float>       float_default_map;
//   std::map<std::string,std::string> string_default_map;
//   lsdtt::lsdtt_help_registry        help_map;
//
//   lsdtt::declare_routing_core_parameters(
//       int_default_map, string_default_map, help_map);
//   lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);
//   lsdtt::declare_fill_min_slope_parameter(float_default_map, help_map);
//   lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);
//
//   // A CLI that genuinely wants a different default overwrites AFTER the call:
//   int_default_map["threshold_contributing_pixels"] = 200;
//   help_map["threshold_contributing_pixels"][1]     = "200";
//
// OVERWRITE POLICY
//
// The helper ALWAYS writes its key into the default map and the help registry,
// unconditionally. It does not check whether the caller has already declared
// the key, and it does not preserve a caller value set beforehand. A CLI that
// needs a different default therefore calls the helper first and overwrites
// afterwards, as shown above, and must also update element [1] of the help
// entry so that the documented default and the actual default cannot diverge.
// The rule was chosen because the alternative (first-writer-wins) makes the
// effective default depend on declaration order, which is exactly the kind of
// invisible divergence this file exists to remove.
//
// ORDERING
//
// `lsdtt_help_registry` records first-insertion order and the terminal help is
// printed in that order, so call these functions at the point in the CLI's
// parameter block where the group belongs in the pipeline, not all together at
// the top.
//
// PREFIX RULE
//
// Every function takes an optional `prefix` so a CLI that runs a group twice
// (two flow-routing stages, a coarse and a fine synthetic domain) can declare
// both. The rule is:
//
//   - empty prefix (the default) -> keys are used verbatim:
//       `raster_is_filled`
//   - non-empty prefix          -> `prefix` + `_` + key:
//       prefix "coarse" -> `coarse_raster_is_filled`
//   - a prefix that already ends in `_` does not gain a second underscore:
//       prefix "coarse_" -> `coarse_raster_is_filled`
//
// The help-entry GROUP label (field 4) also gains the prefix, as
// `"<group> (<prefix>)"`, so the two stages do not collapse into one section in
// the help output. Cross-references inside the help TEXT name the unprefixed
// key, since prefixing them would require string surgery on every sentence.
//
// PROVENANCE OF THE NAMES AND DEFAULTS
//
// Canonical names and defaults come from §3 of the remediation plan where that
// table pins them. Where it does not, the name and default were taken from the
// existing CLI that already spells it in the most conformant way:
//
//   raster_is_filled, fill_min_slope,
//   boundary_conditions, n_threads                cli/lsdtt-channel-extraction.cpp
//   fill_before_routing                           cli/lsdtt-hillslope-channel-coupling.cpp
//   threshold_contributing_pixels                 plan D1/D3 (int, 1000)
//   write_output_format, write_vector_format      cli/lsdtt-fluvial-surfaces.cpp
//   calculate_spatial_index                       cli/lsdtt-basin-extraction.cpp
//   basin_outlet_* (all five)                     cli/lsdtt-basin-extraction.cpp
//   polyfit_lambda_metres, polyfit_metrics,
//   compute_backend                               cli/lsdtt-surface-metrics.cpp
//   n_rows, n_cols, cell_size_metres              cli/lsdtt-raster-synthetic-tests.cpp
//   random_seed                                   cli/lsdtt-synthetic-networks.cpp
//   synthetic_crs_wkt, noise_type,
//   noise_target_stddev_metres                           cli/lsdtt-lem.cpp
//   noise_relief_metres, diamond_square_roughness,
//   diamond_square_initial_scale                  cli/lsdtt-raster-synthetic-tests.cpp
//
// Two names declared here are NOT yet spelled this way by any CLI, and are
// deliberate renames that adoption will carry:
//
//   cell_size_metres  (CLIs say `cell_size` / `grid_resolution`; §3 requires
//       units spelled out)
//   noise_type / noise_relief_metres  (plan §9.4 names; the closest existing keys are
//       lem's `initial_surface_type` and raster-synthetic-tests'
//       `diamond_square_noise_relief_metres`)
//
// THE TWO FILL FLAGS ARE NOT THE SAME PARAMETER
//
// `raster_is_filled` and `fill_before_routing` are one word apart and ask
// different questions (plan section 2.2 / decision D2). They have separate
// declaration functions, and the help text of each names the other:
//
//   raster_is_filled     - "has the input ALREADY been filled?" Filling is
//                          MANDATORY for the operation; this flag only lets the
//                          CLI skip work it would otherwise do. Setting it
//                          wrongly gives you BAD routing, not different routing.
//   fill_before_routing  - "SHOULD I fill at all?" Filling is OPTIONAL and this
//                          is a genuine analytical choice; false means routing
//                          is built on the unfilled surface and flow paths
//                          terminate at sinks.
//
// A CLI declares whichever of the two matches its operation. It should declare
// both only if it genuinely offers both questions.
// -----------------------------------------------------------------------------

#pragma once

#include <map>
#include <string>

#include "lsdtt_param_parser.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Applies the prefix rule described in the file header to a bare parameter key.
// Exposed so a CLI can build the key it needs to read back out of the parsed
// maps after declaring a prefixed group.
// -----------------------------------------------------------------------------
[[nodiscard]] std::string prefixed_parameter_key(
    const std::string& prefix,
    const std::string& key);

// -----------------------------------------------------------------------------
// Flow-routing core.
//
//   boundary_conditions           string  "oooo"
//   n_threads                     int     0
// -----------------------------------------------------------------------------
void declare_routing_core_parameters(
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

// Mandatory-fill assertion: filling is required; true skips work already done.
void declare_raster_is_filled_parameter(
    std::map<std::string, bool>& bool_default_map,
    lsdtt_help_registry&         help_map,
    const std::string&           prefix = "");

// Optional-fill choice: false deliberately routes into sinks.
void declare_fill_before_routing_parameter(
    std::map<std::string, bool>& bool_default_map,
    lsdtt_help_registry&         help_map,
    const std::string&           prefix = "");

// Fill slope is independent because not every caller exposing a fill concept
// performs the fill itself or spells the stage slope without a prefix.
void declare_fill_min_slope_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    const std::string&            prefix = "");

// Channel-initiation threshold (plan D1/D3).
void declare_channel_initiation_parameters(
    std::map<std::string, int>& int_default_map,
    lsdtt_help_registry&        help_map,
    const std::string&          prefix = "");

// -----------------------------------------------------------------------------
// Output formats. Declare only the capabilities a CLI actually exposes.
// -----------------------------------------------------------------------------
void declare_raster_output_format_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

void declare_vector_output_format_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

void declare_spatial_index_parameter(
    std::map<std::string, bool>& bool_default_map,
    lsdtt_help_registry&         help_map,
    const std::string&           prefix = "");

// Genuine per-output overrides inherit the canonical selector by default.
inline constexpr const char* inherit_write_vector_format =
    "same_as_write_vector_format";
inline constexpr const char* inherit_write_output_format =
    "same_as_write_output_format";

[[nodiscard]] std::string resolve_output_format_override(
    const std::string& override_value,
    const std::string& canonical_value);

// -----------------------------------------------------------------------------
// Basin outlet input and outlet snapping.
//
//   basin_outlet_fname                          string  "NULL"
//   basin_outlet_snapping_method                string  "downstream"
//   basin_outlet_contributing_pixels_threshold  int     1000
//   basin_outlet_snap_min_search_radius_pixels  int     0
//   basin_outlet_snap_max_search_radius_pixels  int     25
// -----------------------------------------------------------------------------
void declare_basin_outlet_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

// -----------------------------------------------------------------------------
// Polynomial surface fitting.
//
//   polyfit_lambda_metres  float   70.0    (HALF-WAVELENGTH, not a radius)
//   polyfit_metrics        string  ""
//   compute_backend        string  "auto"
// -----------------------------------------------------------------------------
void declare_surface_fitting_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

// -----------------------------------------------------------------------------
// Synthetic model domain (plan §9.4). Each concept is independently selectable
// because the synthetic CLIs do not share one CRS model or every domain knob.
//
//   n_rows            int     128
//   n_cols            int     128
//   cell_size_metres  float   1.0
//   random_seed       int     33
//   synthetic_crs_wkt string  ""
// -----------------------------------------------------------------------------
void declare_synthetic_dimensions_parameters(
    std::map<std::string, int>&         int_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

void declare_synthetic_cell_size_parameter(
    std::map<std::string, float>&       float_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

void declare_synthetic_random_seed_parameter(
    std::map<std::string, int>&         int_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

void declare_synthetic_crs_wkt_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

// Convenience composition for callers that genuinely expose every concept.
void declare_synthetic_domain_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

// -----------------------------------------------------------------------------
// Synthetic noise surfaces (plan §9.4).
//
//   noise_type                    string  "flat"
//   noise_target_stddev_metres           float   1.0
//   noise_relief_metres                  float   0.0
//   diamond_square_roughness      float   0.55
//   diamond_square_initial_scale  float   1.0
// -----------------------------------------------------------------------------
void declare_synthetic_noise_type_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

void declare_synthetic_noise_target_stddev_parameter(
    std::map<std::string, float>&       float_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

void declare_synthetic_noise_relief_parameter(
    std::map<std::string, float>&       float_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

void declare_diamond_square_parameters(
    std::map<std::string, float>&       float_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

// Convenience composition for callers that genuinely expose every concept.
void declare_synthetic_noise_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

} // namespace lsdtt
