// lsdtt_common_parameters.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
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
// This is Phase 0a of the CLI consistency clean-up.
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
// Canonical names and defaults follow the CLI parameter conventions. Where
// those do not pin them, the name and default were taken from the existing CLI that already spells it in the most conformant way:
//
//   raster_is_filled, fill_min_slope,
//   boundary_conditions, n_threads                cli/lsdtt-channel-extraction.cpp
//   fill_before_routing                           cli/lsdtt-hillslope-channel-coupling.cpp
//   threshold_contributing_pixels                 (int, 1000)
//   write_output_format, write_vector_format      cli/lsdtt-fluvial-surfaces.cpp
//   calculate_spatial_index                       cli/lsdtt-basin-extraction.cpp
//   basin_outlet_* (all five)                     cli/lsdtt-basin-extraction.cpp
//   polyfit_lambda_metres, polyfit_metrics,
//   compute_backend                               cli/lsdtt-surface-metrics.cpp
//   n_rows, n_cols, cell_size_metres              cli/lsdtt-synthetic-surfaces.cpp
//   random_seed                                   cli/lsdtt-synthetic-networks.cpp
//   synthetic_crs_wkt, noise_type,
//   noise_target_stddev_metres                           cli/lsdtt-lem.cpp
//   noise_relief_metres, diamond_square_roughness,
//   diamond_square_initial_scale                  cli/lsdtt-synthetic-surfaces.cpp
//
// Two names declared here are NOT yet spelled this way by any CLI, and are
// deliberate renames that adoption will carry:
//
//   cell_size_metres  (CLIs say `cell_size` / `grid_resolution`; §3 requires
//       units spelled out)
//   noise_type / noise_relief_metres  (the closest existing keys are
//       lem's `initial_surface_type` and raster-synthetic-tests'
//       `diamond_square_noise_relief_metres`)
//
// THE TWO FILL FLAGS ARE NOT THE SAME PARAMETER
//
// `raster_is_filled` and `fill_before_routing` are one word apart and ask
// different questions. They have separate
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
#include "../utils/lsdtt_threshold_selection.hpp"

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

// How flats are made drainable after filling (fill_flat_resolution_choice,
// string, default inline_gradient = the existing Priority-Flood behaviour).
// Parse the value with lsdtt::parse_flat_resolution_choice.
// min_slope_key names the CLI's fill slope parameter in the help text (a few
// CLIs call it fill_min_slope_for_routing).
void declare_fill_flat_resolution_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "",
    const std::string&                  min_slope_key = "fill_min_slope");

// Interior nodata: nodata pixels not connected to the raster edge through
// nodata (isolated voids, moulins, lakes masked out). Edge-connected nodata,
// such as sea removed by min_elevation, is always an open boundary.
//
//   interior_nodata_method      string  "sink"       sink | fill
//   nodata_fill_method_choice   string  "laplacian"  laplacian | iterative_weighted_average
//   nodata_fill_max_iter        int     10000
//
// Apply with lsdtt::apply_interior_nodata_method (lsdtt_raster_nodata_filling_cpu).
void declare_interior_nodata_parameters(
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix = "");

// Channel-initiation threshold.
void declare_channel_initiation_parameters(
    std::map<std::string, int>& int_default_map,
    lsdtt_help_registry&        help_map,
    const std::string&          prefix = "");

// -----------------------------------------------------------------------------
// Flat-ground admission and floodplain bounds, shared by
// lsdtt-flat-surface-finder, lsdtt-fluvial-surfaces, lsdtt-mountain-front and
// lsdtt-hillslope-channel-coupling. One name per aim (2026-10-26); each tool
// declares only the concepts it uses. Defaults differ legitimately between the
// tools, so each is an argument rather than an overwrite after the call: the
// help text then always states the default the tool actually uses.
//
//   relief_window_radius_metres       float   0 = default_cells x cell size
//   max_relief_metres                 float   default_metres
//   relief_threshold_algorithm_choice string  fixed | cdf_divergence
//   max_admission_gradient            float   m/m
//   max_floodplain_hand_metres        float   m
//   max_floodplain_gradient           float   m/m
//
// Gradients are m/m, never angles.
// -----------------------------------------------------------------------------
void declare_relief_window_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    int                           default_cells,
    const std::string&            parent = "");

void declare_max_relief_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    float                         default_metres,
    const std::string&            parent = "");

void declare_relief_threshold_algorithm_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  default_choice);

void declare_max_admission_gradient_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    float                         default_gradient,
    const std::string&            parent = "");

void declare_max_floodplain_hand_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    float                         default_metres,
    const std::string&            parent = "");

// Morphological clean-up of a valley mask, shared by lsdtt-fluvial-surfaces
// (on the mask it writes) and lsdtt-valley-centreline (on the mask it reads):
// valley_mask_closing_radius_metres (dilation then erosion: closes narrow
// non-valley intrusions and slivers), valley_mask_opening_radius_metres
// (erosion then dilation: removes narrow valley tendrils) and
// fill_valley_mask_holes. Applied in that order.
void declare_valley_mask_morphology_parameters(
    std::map<std::string, float>& float_default_map,
    std::map<std::string, bool>&  bool_default_map,
    lsdtt_help_registry&          help_map,
    const std::string&            group,
    const std::string&            parent,
    bool                          fill_holes_default);

// Valley centrelines and the relative elevation measured from them, shared by
// lsdtt-valley-centreline (run_valley_width_network, write_relative_elevation_
// raster) and lsdtt-fluvial-surfaces (reference_field_algorithm_choice
// relative_elevation): network_centreline_algorithm_choice and its two
// settings, centreline_baseline_algorithm_choice and
// relative_elevation_max_distance_metres.
void declare_relative_elevation_parameters(
    std::map<std::string, std::string>& string_default_map,
    std::map<std::string, float>&       float_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  group,
    const std::string&                  parent);

void declare_max_floodplain_gradient_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    float                         default_gradient,
    const std::string&            parent = "");

// The ladder's selection rule - the largest level count whose band centres are
// at least kappa * E_v apart - as used by lsdtt-fluvial-surfaces on HAND and by
// lsdtt-valley-centreline on cross-section break heights:
//   ladder_break_algorithm_choice  string  antimode   gap | antimode | jenks
//   ladder_max_level_count         int     12
//   vertical_accuracy_metres       float   0          E_v; 0 disables the floor
//   ladder_resolvability_kappa     float   3
void declare_ladder_selection_parameters(
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map);

// How the floodplain bounds (max_floodplain_hand_metres,
// max_floodplain_gradient) are chosen: fixed | qq_gaussian | cdf_divergence.
// Shared by lsdtt-fluvial-surfaces and lsdtt-hillslope-channel-coupling.
void declare_floodplain_threshold_algorithm_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  default_choice);

// The Clubb et al. (2017) QQ rule's controls, for one quantity:
//   <prefix>_qq_lower_percentile    float  25
//   <prefix>_qq_upper_percentile    float  75
//   <prefix>_qq_threshold_fraction  float  0.01  (of the data range)
// e.g. prefix "floodplain" or "relief".
void declare_qq_threshold_parameters(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    const std::string&            prefix,
    const std::string&            parent);

// A threshold_selection_request from a *_threshold_algorithm_choice value, the
// fixed value to use under `fixed`, and the <prefix>_qq_* parameters. Throws
// std::invalid_argument on an unknown choice, naming choice_key and the valid
// values.
[[nodiscard]] threshold_selection_request threshold_request_from_parameters(
    const std::string&                  choice,
    const std::string&                  choice_key,
    double                              fixed_value,
    const std::map<std::string, float>& float_parameters,
    const std::string&                  prefix);

// Window radius in whole pixels from relief_window_radius_metres: 0 gives
// default_cells; otherwise the nearest whole number of cells, at least 1.
// Throws std::invalid_argument on a negative or non-finite radius.
[[nodiscard]] int relief_window_radius_pixels(double radius_metres,
                                              double cell_size,
                                              int    default_cells);

// Conversions for compatibility aliases (lsdtt_deprecated_compatibility_alias),
// applied by the CLI after parsing. Each throws std::invalid_argument naming
// old_key when the stored text is not a valid old value.
//
// gradient_from_legacy_degrees: tan of an angle in [0, 90); 0 stays exactly 0.
// relief_threshold_choice_from_legacy_bool: true -> cdf_divergence,
//   false -> fixed, accepting the parser's boolean spellings.
[[nodiscard]] double gradient_from_legacy_degrees(const std::string& text,
                                                  const std::string& old_key);
[[nodiscard]] std::string relief_threshold_choice_from_legacy_bool(
    const std::string& text, const std::string& old_key);

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
// Default of basin_outlet_contributing_pixels_threshold. Every CLI that takes
// the key, and FlowInfo.snap_points in the Python bindings, use this one value
// so the same (absent) parameter snaps to the same pixel everywhere.
// -----------------------------------------------------------------------------
inline constexpr int default_basin_outlet_contributing_pixels_threshold = 1000;

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
// Synthetic model domain. Each concept is independently selectable
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
// Synthetic noise surfaces.
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

// -----------------------------------------------------------------------------
// Fetch tools (lsdtt-fetch-raster, lsdtt-fetch-lithology)
//
// declare_fetch_bbox_parameters: west, east, south, north (WGS84 degrees, as
//   strings so the parser keeps double precision).
// declare_fetch_bounds_file_parameters: bounds_fname, bounds_path, bounds_layer.
// declare_http_request_parameters: http_timeout_seconds,
//   http_connect_timeout_seconds, http_max_retry.
// declare_linear_feature_fetch_parameters: see its declaration below.
// -----------------------------------------------------------------------------
void declare_fetch_bbox_parameters(std::map<std::string, std::string>& string_default_map,
                                   lsdtt_help_registry&                help_map);
void declare_fetch_bounds_file_parameters(std::map<std::string, std::string>& string_default_map,
                                          lsdtt_help_registry&                help_map);
void declare_http_request_parameters(std::map<std::string, int>& int_default_map,
                                     lsdtt_help_registry&        help_map);

// declare_linear_feature_fetch_parameters: fetch_from_linear_feature (the
//   switch), fetch_feature_fname, fetch_feature_path, fetch_feature_layer,
//   linear_feature_fetch_width, linear_feature_fetch_end_buffer,
//   rotate_after_fetch. Fetch along a line, two points or a polygon into an
//   oblique frame (lsdtt_oblique_frame.hpp).
void declare_linear_feature_fetch_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map);

// -----------------------------------------------------------------------------
// Oblique frames in lsdtt-raster-preprocessing
//
// declare_linear_feature_rotation_parameters: rotate_to_linear_feature,
//   rotation_feature_fname, rotation_feature_path, rotation_feature_layer,
//   linear_feature_rotation_width, linear_feature_rotation_end_buffer,
//   rotate_to_data_footprint. Rotate an existing DEM into the oblique frame
//   of a feature file or of its own valid cells (lsdtt_oblique_frame.hpp).
// declare_rotate_back_parameters: run_rotate_back, rotation_frame_fname,
//   rotate_back_file_list, rotate_back_resampling_method, rotate_back_epsg.
//   Warp rasters in a frame back onto the frame file's target lattice.
// -----------------------------------------------------------------------------
void declare_linear_feature_rotation_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map);
void declare_rotate_back_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map);

// -----------------------------------------------------------------------------
// Orientations (every CLI that writes a bearing, aspect or orientation)
//
// declare_true_north_bearing_parameters: bearings_relative_to_true_north.
//   Off: bearings are clockwise from grid north (grid up). On: each output is
//   converted at its own location with the grid convergence
//   (lsdtt_true_north.hpp). Hillshade illumination is never converted.
// -----------------------------------------------------------------------------
void declare_true_north_bearing_parameters(
    std::map<std::string, bool>& bool_default_map,
    lsdtt_help_registry&         help_map);

// -----------------------------------------------------------------------------
// Lithology input (lsdtt-lithology-erodibility and analyses by lithologic unit)
//
// declare_lithology_map_parameters: lithology_fname, lithology_source,
//   lithology_layer_name, lithology_field_name, lithology_crosswalk_fname. The
//   help states the codes-only rules of read_lithology_units_on_grid; a CLI
//   with extra rules (lsdtt-lithology-erodibility's table-name matching)
//   overwrites the guidance element [3] after the call.
// declare_write_lithology_code_raster_parameter: write_lithology_code_raster
//   (default true; analysis CLIs overwrite the default to false).
// -----------------------------------------------------------------------------
void declare_lithology_map_parameters(std::map<std::string, std::string>& string_default_map,
                                      lsdtt_help_registry&                help_map,
                                      const std::string&                  prefix = "");
void declare_write_lithology_code_raster_parameter(std::map<std::string, bool>& bool_default_map,
                                                   lsdtt_help_registry&         help_map,
                                                   const std::string&           prefix = "");

// -----------------------------------------------------------------------------
// Channel region input (grouping channel statistics by user-drawn regions)
//
// declare_channel_region_parameters: channel_region_fname,
//   channel_region_source, channel_region_layer_name, channel_region_field_name.
//   Regions are read like a lithology map (polygon vector or integer code
//   raster). Used by lsdtt-channel-network-analysis to group link bearings and
//   junction angles by region.
// -----------------------------------------------------------------------------
void declare_channel_region_parameters(std::map<std::string, std::string>& string_default_map,
                                       lsdtt_help_registry&                help_map,
                                       const std::string&                  prefix = "");

// -----------------------------------------------------------------------------
// Lithology summaries (per-unit statistics of channel metrics)
//
// Split into three so that no CLI advertises a key it does not use.
//
// declare_lithology_ksn_summary_parameters: lithology_contact_buffer_metres
//   (float 0), lithology_min_reach_pixels (int 10), lithology_bootstrap_cluster
//   (string "reach"). Per-unit k_sn summaries (lsdtt-chi-analysis).
// declare_lithology_network_summary_parameters: lithology_min_majority_fraction
//   (float 0.8). Per-unit link / junction statistics
//   (lsdtt-channel-network-analysis).
// declare_lithology_bootstrap_parameters: lithology_bootstrap_replicates
//   (int 1000), lithology_bootstrap_random_seed (int 5489). Shared by both.
// -----------------------------------------------------------------------------
void declare_lithology_ksn_summary_parameters(std::map<std::string, float>&       float_default_map,
                                              std::map<std::string, int>&         int_default_map,
                                              std::map<std::string, std::string>& string_default_map,
                                              lsdtt_help_registry&                help_map,
                                              const std::string&                  prefix = "");
void declare_lithology_network_summary_parameters(std::map<std::string, float>& float_default_map,
                                                  lsdtt_help_registry&          help_map,
                                                  const std::string&            prefix = "");
void declare_lithology_bootstrap_parameters(std::map<std::string, int>& int_default_map,
                                            lsdtt_help_registry&        help_map,
                                            const std::string&          prefix = "");

// -----------------------------------------------------------------------------
// Erosion-rate samples and erodibility calibration (stream-power K from
// basin-averaged erosion rates; lsdtt_erodibility_calibration_cpu)
//
// declare_erosion_rate_input_parameters (group "Erosion-rate samples"):
//   erosion_rate_fname, erosion_rate_layer_name,
//   erosion_rate_sample_id_column, erosion_rate_latitude_column,
//   erosion_rate_longitude_column, erosion_rate_column,
//   erosion_rate_uncertainty_column, erosion_rate_units (required, no
//   default), erosion_rate_reported_area_column, erosion_rate_snapping_method,
//   erosion_rate_snap_max_distance_metres, erosion_rate_area_mismatch_tolerance,
//   erosion_rate_area_mismatch_action, erosion_rate_nested_basins.
// declare_erodibility_calibration_parameters (group "Erodibility calibration"):
//   erosion_rate_ksn_statistic, erosion_rate_min_reach_pixels,
//   erosion_rate_min_reaches_per_basin, erosion_rate_fit_n,
//   erosion_rate_min_groups_for_n_fit, erosion_rate_min_ksn_range_ratio,
//   erosion_rate_mono_lithology_min_area_fraction,
//   erosion_rate_min_unit_reaches_per_basin, erosion_rate_unit_k_methods,
//   erosion_rate_property_table_method,
//   erosion_rate_systematic_relative_uncertainty,
//   erosion_rate_bootstrap_replicates, erosion_rate_bootstrap_random_seed,
//   write_erodibility_calibration_basin_polygons,
//   write_erodibility_calibration_property_table.
// Every key is effective only with erosion_rate_fname. The list of n values
// is not declared here.
// -----------------------------------------------------------------------------
// declare_stream_power_n_values_parameter: stream_power_n_values (string "1"),
// the slope exponent(s) n shared by relative erodibility and the erodibility
// calibration in lsdtt-chi-analysis (relative_erodibility_n_values is its
// deprecated alias there).
void declare_stream_power_n_values_parameter(std::map<std::string, std::string>& string_default_map,
                                             lsdtt_help_registry&                help_map,
                                             const std::string&                  prefix = "");
void declare_erosion_rate_input_parameters(std::map<std::string, float>&       float_default_map,
                                           std::map<std::string, std::string>& string_default_map,
                                           lsdtt_help_registry&                help_map,
                                           const std::string&                  prefix = "");
void declare_erodibility_calibration_parameters(std::map<std::string, bool>&        bool_default_map,
                                                std::map<std::string, int>&         int_default_map,
                                                std::map<std::string, float>&       float_default_map,
                                                std::map<std::string, std::string>& string_default_map,
                                                lsdtt_help_registry&                help_map,
                                                const std::string&                  prefix = "");

} // namespace lsdtt
