// lsdtt-fluvial-surfaces.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-fluvial-surfaces
//
// Stratifies flat ground into fluvial levels -- floodplain, terrace treads, fan
// lobes -- by banding a REFERENCE FIELD and labelling connected components
// within each band.
//
// WHERE THIS SITS. lsdtt-flat-surface-finder answers "which ground is flat?"
// from morphometry alone: no flow routing, no channel, no height reference.
// This tool answers the next question, "which fluvial level is each flat
// surface on?", which is the first one that needs a channel network. It takes
// the flat mask as input (flat_mask_fname) and does everything that needs a
// height reference. lsdtt-valley-centreline then takes the valley mask this
// tool writes and traces the centreline through it.
//
// THE REFERENCE FIELD. Two choices, one ladder:
//
//   hand            height above the in-scope channel network. Terraces along a
//                   reach sit at a roughly constant height above the channel,
//                   so HAND makes a terrace flight a step function.
//   fitted_surface  the signed residual to a surface fitted to the flat ground.
//                   A fan is NOT a step function in HAND - it climbs
//                   continuously from toe to apex, so a HAND ladder cuts it
//                   into arcs that correspond to nothing. Its inset lobes sit a
//                   roughly constant depth below the projection of the older
//                   surface, which is what the residual measures.
//
// BAND FIRST, THEN LABEL. Connected components are labelled once per band, not
// once over the whole mask with height consulted afterwards. Two pixels can
// only share a label if they share a band, so two treads at different heights
// cannot merge - and the riser between them never has to be found, and never
// has to be continuous. A road or gully cutting a riser is simply irrelevant.
//
// THERE IS NO VALLEY MODE AND NO TERRACE MODE. Finding a valley floor and
// segmenting a terrace flight are the same computation at different band
// counts: the valley floor is the lowest band of a two-band cut. Set
// ladder_level_count to 2 for the valley cut, to n for an n-level ladder, or
// leave it at 0 and let the resolvability floor decide how many levels the DEM
// can actually carry.
//
// Main outputs:
//   <stem>_band_index.tif         Int16 band per pixel, -1 = not banded
//   <stem>_banded_labels.tif      UInt32 patch label, a global running counter
//   <stem>_valley_mask.tif        band 0: the input to lsdtt-valley-centreline
//   <stem>_surface_class.tif      floodplain / terrace class (hand field only)
//   <stem>_surface_patches.csv    one row per patch - the join key downstream
//   <stem>_ladder_windows.csv     the resolvability report: read this first
//   <stem>_reference_field.tif    the field that was actually banded
//
// Usage:
//   lsdtt-fluvial-surfaces /path/to/ params.driver
//   lsdtt-fluvial-surfaces -h
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_true_north.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/cpu/lsdtt_raster_vectorise_cpu.hpp"
#include "../src/cpu/lsdtt_surface_method_comparison_cpu.hpp"
#include "../src/cpu/lsdtt_surface_training_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_channel_scope_cpu.hpp"
#include "../src/cpu/lsdtt_flat_surface_channel_classify_cpu.hpp"
#include "../src/core/lsdtt_stat_csv.hpp"
#include "../src/cpu/lsdtt_floodplain_thresholds_cpu.hpp"
#include "../src/cpu/lsdtt_proximity_cpu.hpp"
#include "../src/cpu/lsdtt_raster_binary_mask_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_reference_field_ladder_cpu.hpp"
#include "../src/cpu/lsdtt_valley_cross_sections_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-fluvial-surfaces";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-fluvial-surfaces                    ||\n"
        << "|| I sort flat ground into fluvial levels for you:      ||\n"
        << "|| floodplains, terrace treads and fan lobes, by        ||\n"
        << "|| banding a height reference and labelling within      ||\n"
        << "|| each band.                                           ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();
    lsdtt::lsdtt_cli_args cli_args =
        lsdtt::lsdtt_param_parser::parse_args_deferred(
            argc, argv, k_tool_name, k_status, k_citation);

    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry         help_map;

    // --- General ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of threads. 0 = all available.",
        "OpenMP thread count for every parallel stage.",
        "General","io",""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Raster output format: cog, gtiff or same_as_input.",
        "cog is a Cloud Optimised GeoTIFF and is the project default.",
        "General","io",""};

    string_default_map["write_vector_format"] = "flatgeobuf";
    help_map["write_vector_format"] = {"string","flatgeobuf",
        "Vector output format: flatgeobuf, gpkg or csv.",
        "Governs every vector layer this tool writes, including the channel "
        "network. FlatGeobuf is preferred: it is compact, carries its own "
        "georeferencing, and does not have GeoPackage's problems writing "
        "across drives (a common WSL2 situation).",
        "General","io",""};

    // --- Method: a preset for the stage parameters ---
    // Every stage parameter defaults to "preset" and takes the preset's value
    // unless the parameter file names it; the run log says which is which.
    string_default_map["fluvial_surface_algorithm_choice"] = "ladder";
    help_map["fluvial_surface_algorithm_choice"] = {"string","ladder",
        "Method preset: ladder, clubb2017 or cross_section_levels.",
        "Fills every stage parameter left at `preset`. ladder (default): the "
        "admitted ground (flat_mask_fname, or the built-in relief admission) "
        "is banded by a data-driven ladder on HAND, components labelled within "
        "each band. clubb2017: Clubb et al. (2017, Earth Surf. Dynam. "
        "5:369-385) - floodplain where HAND above a Strahler-order channel and "
        "gradient are both below bounds read off quantile-quantile plots, "
        "floodplain patches must contain a channel pixel, terraces above the "
        "floodplain bound measured from the nearest channel (swath baseline). "
        "cross_section_levels: the ladder, but with the levels read from "
        "supplied_levels_fname (normally <stem>_cross_section_levels.csv from "
        "lsdtt-valley-centreline) instead of found. "
        "Any stage named in the parameter file overrides the preset, so e.g. "
        "clubb2017 with level_algorithm_choice: ladder bands the Clubb "
        "admission with the ladder. Stage values (ladder / clubb2017 / "
        "cross_section_levels): admission_algorithm_choice flat_mask / "
        "floodplain_thresholds / flat_mask; "
        "floodplain_threshold_algorithm_choice fixed / qq_gaussian / fixed; "
        "channel_definition_algorithm_choice contributing_pixels / "
        "stream_order / contributing_pixels; level_algorithm_choice ladder / "
        "fixed_thresholds / supplied_levels; "
        "terrace_reference_field_algorithm_choice same / swath_baseline / "
        "same; segmentation_algorithm_choice within_bands / whole_mask / "
        "within_bands. Several presets in one run: list them comma-separated "
        "with no spaces (ladder,clubb2017); each writes under "
        "<write_fname>_<preset>, then the comparison outputs "
        "(write_method_record_csv).",
        "Method","switch",""};

    string_default_map["level_algorithm_choice"] = "preset";
    help_map["level_algorithm_choice"] = {"string","preset",
        "How surfaces are put on levels: preset, ladder, fixed_thresholds or "
        "supplied_levels.",
        "ladder: data-driven breaks in the reference field, as many levels as "
        "the DEM can resolve (the ladder_* parameters). fixed_thresholds: two "
        "levels from two bounds - floodplain below max_floodplain_hand_metres "
        "(and touching the channel), terrace above it by "
        "min_terrace_height_above_floodplain_metres - with the bounds chosen "
        "by floodplain_threshold_algorithm_choice. fixed_thresholds needs "
        "reference_field_algorithm_choice hand. supplied_levels: the ladder's "
        "machinery with the levels read from supplied_levels_fname - band "
        "breaks at the midpoints between consecutive levels, per channel - "
        "instead of found; band centres and the resolvability floor are still "
        "measured and reported (status supplied in _ladder_windows.csv).",
        "Method","param","fluvial_surface_algorithm_choice"};

    string_default_map["supplied_levels_fname"] = "";
    help_map["supplied_levels_fname"] = {"string","",
        "Levels file for level_algorithm_choice supplied_levels, in read_path.",
        "CSV with at least source_id and level_hand_m columns (heights above "
        "the channel, m): <stem>_cross_section_levels.csv from "
        "lsdtt-valley-centreline, or a file written by hand from survey. A "
        "source_id of -1 applies to every channel. Each channel's source_id "
        "must match this run's channel scope - use the same source points and "
        "threshold_contributing_pixels as the run that wrote it - and a "
        "channel with no levels stops the run, naming it.",
        "Method","io","level_algorithm_choice"};

    string_default_map["segmentation_algorithm_choice"] = "preset";
    help_map["segmentation_algorithm_choice"] = {"string","preset",
        "How admitted ground becomes patches: preset, within_bands or "
        "whole_mask.",
        "within_bands (ladder and cross_section_levels presets): connected "
        "components inside each level, so two treads at different heights can "
        "never share a patch. whole_mask (clubb2017 preset): components of the "
        "admitted ground minus the risers in riser_mask_fname, each taking its "
        "majority level. whole_mask is the riser segmentation, kept so it can "
        "be compared: where a riser is broken (a road, a gully) the two treads "
        "either side join into one patch on one level. Under fixed_thresholds "
        "there is one terrace level, so whole_mask only removes the riser "
        "pixels.",
        "Method","param","fluvial_surface_algorithm_choice"};

    string_default_map["riser_mask_fname"] = "";
    help_map["riser_mask_fname"] = {"string","",
        "Riser raster for segmentation_algorithm_choice whole_mask, in "
        "read_path.",
        "Cells > 0.5 are risers and split the admitted ground. Normally "
        "<stem>_riser_centreline.tif from lsdtt-feature-finder (feature_type: "
        "find_risers), on the DEM grid; gap closing there "
        "(riser_gap_closing_algorithm_choice: astar) decides whether broken "
        "risers still separate treads. Empty = no risers: whole_mask then "
        "labels the admitted ground as it stands.",
        "Method","io","segmentation_algorithm_choice"};

    string_default_map["channel_definition_algorithm_choice"] = "preset";
    help_map["channel_definition_algorithm_choice"] = {"string","preset",
        "What counts as a channel: preset, contributing_pixels or stream_order.",
        "contributing_pixels: the whole network down from heads at "
        "threshold_contributing_pixels. stream_order: only its pixels of "
        "Strahler order >= min_channel_stream_order; HAND, the floodplain "
        "contact test and the swath baseline then all use that channel. The "
        "threshold still decides which first-order channels exist.",
        "Method","param","fluvial_surface_algorithm_choice"};

    int_default_map["min_channel_stream_order"] = 3;
    help_map["min_channel_stream_order"] = {"int","3",
        "Lowest Strahler order counted as channel (stream_order only).",
        "Clubb et al. (2017) used third order at all eight of their sites. "
        "LSDTopoTools2 called this threshold_SO.",
        "Method","param","channel_definition_algorithm_choice"};

    string_default_map["terrace_reference_field_algorithm_choice"] = "preset";
    help_map["terrace_reference_field_algorithm_choice"] = {"string","preset",
        "What terrace height is measured from: preset, same or swath_baseline.",
        "same: the main reference field (HAND). swath_baseline: elevation "
        "above the NEAREST channel pixel within swath_half_width_metres, which "
        "does not follow drainage, so a terrace above a confluence is "
        "referenced to the nearest channel rather than the one its flow path "
        "reaches (Clubb et al. 2017). fixed_thresholds only.",
        "Method","param","fluvial_surface_algorithm_choice"};

    float_default_map["swath_half_width_metres"] = 500.0f;
    help_map["swath_half_width_metres"] = {"float","500",
        "Half-width [m] of the swath baseline: terraces further than this from "
        "the channel are not measured.",
        "A HALF-width, as half_width_metres in lsdtt-swath-profile. "
        "LSDTopoTools2's swath_width (default 1000) was a full width.",
        "Method","param","terrace_reference_field_algorithm_choice"};

    // --- Admission: which ground may take part ---
    string_default_map["flat_mask_fname"] = "";
    help_map["flat_mask_fname"] = {"string","",
        "Flat mask raster from lsdtt-flat-surface-finder, in read_path.",
        "Cells > 0.5 are admitted. This is normally "
        "<stem>_combined_flat_mask.tif from lsdtt-flat-surface-finder, which "
        "is the intended two-step workflow: detect flat ground there, sort it "
        "into levels here. Leave empty to use the built-in relief or gradient "
        "admission (admission_algorithm_choice) instead, which is coarser and is meant for a quick look rather than "
        "for a result you intend to publish.",
        "Admission","io",""};

    string_default_map["admission_algorithm_choice"] = "preset";
    help_map["admission_algorithm_choice"] = {"string","preset",
        "Which ground is admitted: preset, flat_mask, relief, gradient or "
        "floodplain_thresholds.",
        "floodplain_thresholds (the clubb2017 preset): cells whose gradient, "
        "read over floodplain_gradient_window_radius_metres, is below "
        "max_floodplain_gradient as chosen by "
        "floodplain_threshold_algorithm_choice, among cells with a defined "
        "HAND. flat_mask (the ladder preset) reads flat_mask_fname, normally "
        "<stem>_combined_flat_mask.tif from lsdtt-flat-surface-finder. relief "
        "and gradient are a coarse built-in admission for a quick look when "
        "there is no flat mask: relief admits cells whose local relief over "
        "relief_window_radius_metres is below max_relief_metres; gradient "
        "admits cells whose gradient is below max_admission_gradient. "
        "flat_mask with flat_mask_fname empty falls back to relief, as "
        "earlier versions did. Replaces fallback_admission_algorithm_choice "
        "(whose `slope` is now `gradient`).",
        "Admission","param",""};

    // Built-in admission, shared names (lsdtt_common_parameters): the same
    // window and bound as lsdtt-flat-surface-finder and lsdtt-mountain-front.
    // Five cells is the radius this tool has always used (an 11 x 11 window).
    lsdtt::declare_relief_window_parameter(float_default_map, help_map, 5,
                                           "admission_algorithm_choice");
    lsdtt::declare_max_relief_parameter(float_default_map, help_map, 1.0f,
                                        "admission_algorithm_choice");
    lsdtt::declare_max_admission_gradient_parameter(
        float_default_map, help_map, 0.05f, "admission_algorithm_choice");
    // Hidden storage for compatibility aliases, resolved after the DEM loads.
    string_default_map["legacy_relief_window_radius_pixels"] = "";
    string_default_map["legacy_fallback_admission_algorithm_choice"] = "";

    bool_default_map["fill_admitted_holes"] = false;
    help_map["fill_admitted_holes"] = {"bool","false",
        "Fill enclosed holes in the admitted ground (admission flat_mask, "
        "relief or gradient).",
        "Roughness - DEM noise, crops, buildings - leaves enclosed gaps in a "
        "relief- or gradient-based admission, which become holes in every "
        "surface built on it. A gap is filled only if it is entirely "
        "surrounded by admitted ground. Not applied to the "
        "floodplain_thresholds admission, whose floodplain has "
        "fill_floodplain_holes.",
        "Admission","param",""};

    // Training polygons (fluvial-surface methods plan 7i).
    string_default_map["training_surfaces_fname"] = "";
    help_map["training_surfaces_fname"] = {"string","",
        "Polygons (gpkg, fgb, shp) of floodplain and terrace drawn by hand, used "
        "to choose the HAND thresholds.",
        "Class values: floodplain or fp; terrace, terraceN, tN or 'terrace N' "
        "(case and spaces ignored). Terrace numbers are ignored: the levels come "
        "from the HAND inside the polygons. Each polygon is summarised by its "
        "median and its 10-90 percentile core; per reach window (the ladder's "
        "windows) the floodplain / terrace bound and the breaks between "
        "terrace levels are chosen by best separation (Youden's J) of those "
        "cores. clubb2017: the bound replaces the data-chosen HAND bound per "
        "window (and the second pass is skipped). ladder: the trained breaks "
        "are the lowest breaks and the ladder's own search runs above the top "
        "one. A window with too little training takes the fit over all "
        "training cells. Writes <stem>_training_levels.tif, "
        "<stem>_training_fit.csv and <stem>_training_polygons.csv. Needs "
        "reference_field_choice hand.",
        "Training","param",""};
    string_default_map["training_surfaces_layer_name"] = "";
    help_map["training_surfaces_layer_name"] = {"string","",
        "Layer of training_surfaces_fname; empty = the first layer.",
        "", "Training","param",""};
    string_default_map["training_surfaces_class_field_name"] = "class";
    help_map["training_surfaces_class_field_name"] = {"string","class",
        "Field of training_surfaces_fname holding each polygon's class.",
        "", "Training","param",""};
    float_default_map["training_spread_factor"] = 3.0f;
    help_map["training_spread_factor"] = {"float","3.0",
        "A training polygon whose HAND p90 - p10 exceeds this many required "
        "gaps (kappa x E_v) is not a single level.",
        "Such a polygon (a sloping surface, or one drawn across a riser) still "
        "counts as terrace for the floodplain / terrace bound but takes no part "
        "in the terrace level breaks. Flagged in <stem>_training_polygons.csv. "
        "With E_v 0 nothing is flagged.",
        "Training","param",""};
    int_default_map["min_training_cells_per_class"] = 20;
    help_map["min_training_cells_per_class"] = {"int","20",
        "Core cells a training population needs in a window to take part.",
        "Below it the window falls back to the fit over all training cells, "
        "and below it there the threshold is the data-chosen one.",
        "Training","param",""};

    int_default_map["min_band_component_pixels"] = 50;
    help_map["min_band_component_pixels"] = {"int","50",
        "Minimum patch size in pixels, applied PER BAND.",
        "Components are formed inside each band, so this is a floor on the "
        "size of one surface on one level - not the whole-mask filter that "
        "min_component_pixels is in lsdtt-flat-surface-finder, which is why "
        "the name differs. Replaces min_component_pixels here.",
        "Admission","param",""};

    // --- Reference field ---
    string_default_map["reference_field_algorithm_choice"] = "hand";
    help_map["reference_field_algorithm_choice"] = {"string","hand",
        "Field to band: hand or fitted_surface.",
        "hand = height above the in-scope channel network; the right field for "
        "a terrace flight, where treads sit at a constant height above the "
        "channel. fitted_surface = the signed residual to a surface fitted to "
        "the flat ground; the right field for a fan, which climbs continuously "
        "in HAND and would be cut into meaningless arcs by a HAND ladder.",
        "Reference field","switch",""};

    string_default_map["hand_fname"] = "";
    help_map["hand_fname"] = {"string","",
        "Precomputed height-above-channel raster, in read_path.",
        "THE COMPOSABILITY LEVER. HAND is the expensive stage of this tool and "
        "the level count is what you iterate on, so compute HAND once, write "
        "it with write_hand_raster, then re-run with hand_fname set and try "
        "different ladder settings for the cost of the ladder alone. Flow "
        "routing is still built when a channel scope or network output is "
        "requested, but the HAND computation is skipped.",
        "Reference field","io",""};

    string_default_map["fitted_surface_algorithm_choice"] = "cone";
    help_map["fitted_surface_algorithm_choice"] = {"string","cone",
        "Surface to fit: cone, plane_ransac or poly2.",
        "cone is the fan option and the default: z falling linearly with "
        "distance from an apex. A fan is radial, and a cone is not a low-order "
        "polynomial in x and y - on the synthetic fan fixture the older "
        "surface's own residual spread is about 34 m under a RANSAC plane and "
        "10 m under a poly2, against a 4 m lobe inset, so the misfit swamps "
        "the signal; under a cone it is about 0.3 m. Use plane_ransac for a "
        "piedmont or a terrace tread and poly2 for a broadly curved surface "
        "with no single centre.",
        "Reference field","param","reference_field_algorithm_choice"};

    string_default_map["fitted_surface_mask_fname"] = "";
    help_map["fitted_surface_mask_fname"] = {"string","",
        "Optional raster restricting where the surface is fitted, in read_path.",
        "Cells > 0.5 take part in the fit. Use it to fit on the older fan "
        "surface alone. Optional: a cone fitted to everything, inset lobes "
        "included, still separates the populations - the lobes drag the whole "
        "fit down rather than blurring the boundary, and the ladder does not "
        "care about the absolute offset. Empty = fit over the admitted mask.",
        "Reference field","io","reference_field_algorithm_choice"};

    float_default_map["fitted_surface_ransac_tolerance_metres"] = 5.0f;
    help_map["fitted_surface_ransac_tolerance_metres"] = {"float","5.0",
        "RANSAC inlier tolerance for plane_ransac [m].",
        "Set this to roughly the roughness of the surface being fitted, NOT to "
        "the depth of the features being isolated: a tolerance as large as a "
        "lobe inset lets the lobes join the fit and flattens the signal.",
        "Reference field","param","reference_field_algorithm_choice"};

    int_default_map["fitted_surface_ransac_iterations"] = 1000;
    help_map["fitted_surface_ransac_iterations"] = {"int","1000",
        "Number of RANSAC hypotheses for plane_ransac.",
        "Increase when the inlier fraction is expected to be below about 30%.",
        "Reference field","param","reference_field_algorithm_choice"};

    int_default_map["fitted_surface_seed"] = 1;
    help_map["fitted_surface_seed"] = {"int","1",
        "RANSAC random seed, for reproducibility.",
        "Same seed and same input give the same plane.",
        "Reference field","param","reference_field_algorithm_choice"};

    float_default_map["fitted_surface_apex_x"] = 0.0f;
    help_map["fitted_surface_apex_x"] = {"float","0.0",
        "Fan apex easting in the DEM's projected CRS [m]; cone only.",
        "Leave both apex coordinates at 0 to locate the apex automatically, as "
        "the centroid of the highest 0.1% of fit-mask cells. The apex is "
        "supplied rather than solved for because it is a strongly non-linear "
        "parameter and a user can point at a fan head more reliably than a "
        "solver can find it.",
        "Reference field","param","reference_field_algorithm_choice"};

    float_default_map["fitted_surface_apex_y"] = 0.0f;
    help_map["fitted_surface_apex_y"] = {"float","0.0",
        "Fan apex northing in the DEM's projected CRS [m]; cone only.",
        "See fitted_surface_apex_x.",
        "Reference field","param","reference_field_algorithm_choice"};

    // --- The ladder ---
    int_default_map["ladder_level_count"] = 0;
    help_map["ladder_level_count"] = {"int","0",
        "Number of levels: 0 = auto, 2 = the valley cut, n = an n-level ladder.",
        "0 lets the resolvability floor choose the largest level count the DEM "
        "can actually carry. 2 gives the valley/wall separation that "
        "lsdtt-valley-centreline consumes. An explicit count is honoured and "
        "reported honestly: if that partition is not resolvable the window "
        "still records the count you asked for, with status unresolved, rather "
        "than silently substituting a different one.",
        "Ladder","param",""};



    // Shared with lsdtt-valley-centreline, which clusters cross-section break
    // levels with the same rule (lsdtt_common_parameters).
    lsdtt::declare_ladder_selection_parameters(
        int_default_map, float_default_map, string_default_map, help_map);
    help_map["ladder_break_algorithm_choice"][3] +=
        " Here: antimode bins on log10 and therefore cannot band a field with "
        "negative values, which every fitted_surface residual has. This tool "
        "switches to jenks automatically when you select fitted_surface "
        "without naming a break algorithm, and rejects the combination if you "
        "ask for it explicitly.";

    string_default_map["ladder_gap_break_position"] = "midpoint";
    help_map["ladder_gap_break_position"] = {"string","midpoint",
        "Where in the empty span between clusters the break goes: lower, "
        "midpoint or upper.",
        "The span is nearly empty by construction, so any position separates "
        "the same two clusters - but the sparse pixels INSIDE it (the bank, "
        "the riser, the floodplain margin) go to whichever side the break "
        "leaves them on. Use `upper` when the low band IS the valley floor and "
        "you want the whole margin in it; `lower` for the strictest reading of "
        "the active floodplain; `midpoint` splits the difference.",
        "Ladder","param","ladder_break_algorithm_choice"};

    string_default_map["ladder_gap_rank"] = "lowest";
    help_map["ladder_gap_rank"] = {"string","lowest",
        "Order the gap finder offers its candidate breaks in: lowest or width.",
        "gap algorithm only, and it decides which break wins when more than one "
        "is admissible. `lowest` walks the ladder from the bottom, so the "
        "valley cut is taken before any terrace or hillslope cut. `width` "
        "prefers the widest empty span, which is only a good guide where the "
        "whole field range is a few metres: in a mountain catchment the widest "
        "empty span sits in the sparse tail a hundred metres above the river, "
        "and taking it bands the entire terrace staircase as valley floor.",
        "Ladder","param","ladder_break_algorithm_choice"};

    float_default_map["ladder_gap_trough_max_ratio"] = 0.6f;
    help_map["ladder_gap_trough_max_ratio"] = {"float","0.6",
        "Offer a local minimum as a candidate break when its floor is below "
        "this fraction of the smaller flanking peak. 0 = empty spans only.",
        "gap algorithm only. `ladder_gap_empty_fraction` measures emptiness "
        "against the TALLEST bin, so where one level dominates the population "
        "a real dip above it never reads as empty and the finder has no "
        "candidate anywhere near the valley floor. This offers those dips as "
        "candidates; `ladder_max_separation_ratio` still has to accept them "
        "before any becomes a break, which is why this is the looser of the "
        "two numbers.",
        "Ladder","param","ladder_break_algorithm_choice"};

    float_default_map["ladder_gap_empty_fraction"] = 0.02f;
    help_map["ladder_gap_empty_fraction"] = {"float","0.02",
        "A histogram bin below this fraction of the tallest bin counts as empty.",
        "gap algorithm only. Raise it to treat sparsely-populated ground as a "
        "gap; lower it to require a near-total absence of pixels before calling "
        "a space between clusters.",
        "Ladder","param","ladder_break_algorithm_choice"};

    float_default_map["ladder_gap_min_width_metres"] = 0.5f;
    help_map["ladder_gap_min_width_metres"] = {"float","0.5",
        "Narrowest empty span the gap finder will call a separation [m].",
        "gap algorithm only. Below this a run of empty bins is noise, not the "
        "space between two surfaces.",
        "Ladder","param","ladder_break_algorithm_choice"};

    int_default_map["ladder_gap_hist_bins"] = 80;
    help_map["ladder_gap_hist_bins"] = {"int","80",
        "Linear histogram bins the gap finder searches for empty runs.",
        "Too few and a real gap is bridged by a single populated bin; too many "
        "and ordinary sampling noise fragments a cluster into false gaps.",
        "Ladder","param","ladder_break_algorithm_choice"};

    float_default_map["ladder_max_separation_ratio"] = 0.0f;
    help_map["ladder_max_separation_ratio"] = {"float","0.0",
        "Reject a partition whose worst break sits where the density is above "
        "this fraction of the populations either side. 0 disables.",
        "A SECOND CONSTRAINT alongside the resolvability floor, never blended "
        "with it into a score. The floor asks whether the DEM can MEASURE a "
        "separation; this asks whether there IS one. Both are necessary and "
        "neither implies the other: two lidar terraces 30 cm apart pass this "
        "and need a small vertical_accuracy_metres to pass the floor, while a "
        "26 m cut through one broad mode passes any floor and fails this "
        "outright. 0.6 matches the antimode finder's own acceptance limit. "
        "Off by default because turning it on changes which partitions are "
        "admissible.",
        "Ladder","param",""};

    float_default_map["ladder_min_break_metres"] = 0.0f;
    help_map["ladder_min_break_metres"] = {"float","0.0",
        "Reject any break below this height [m]. 0 = unbounded.",
        "A hard bound applied to every candidate before selection, for when the "
        "lowest plausible level boundary is known independently.",
        "Ladder","param",""};

    float_default_map["ladder_max_break_metres"] = 0.0f;
    help_map["ladder_max_break_metres"] = {"float","0.0",
        "Reject any break above this height [m]. 0 = unbounded.",
        "The companion to ladder_min_break_metres, and the direct lever against "
        "a break landing far above the valley in a landscape whose flat uplands "
        "dominate the population.",
        "Ladder","param",""};

    int_default_map["ladder_k_bootstrap"] = 0;
    help_map["ladder_k_bootstrap"] = {"int","0",
        "Resample each window and re-run the break search this many times, to "
        "report how stable the chosen level count is. 0 = off.",
        "`k_selected` on its own says nothing about whether a different draw "
        "of the same landscape would have given the same answer, and a ladder "
        "that returns k = 7 every time is not the same claim as one that "
        "returns it once in fifteen. This costs B times the break search per "
        "window, and the break search is where the cost sits, not the pixel "
        "count - every finder bins before it searches. Measured on a 3.58 M "
        "cell DEM with 41 windows, 15 replicates added 4.9 s under `gap` "
        "(16.4 to 21.3 s). The same DEM's whole ladder is 6.2 s under "
        "`antimode` and 37.7 s under `jenks`, so scale accordingly: affordable "
        "on gap, think first on jenks. Results go to _ladder_windows.csv as "
        "k_bootstrap_mode, k_stability, k_ci_lo, k_ci_hi and n_bootstrap.",
        "Ladder", "param", ""};

    int_default_map["ladder_k_bootstrap_seed"] = 1;
    help_map["ladder_k_bootstrap_seed"] = {"int","1",
        "Seed for the k bootstrap, so a stability figure is reproducible.",
        "Each window is seeded from this plus its own source and index, so a "
        "window's stability does not shift because a neighbouring reach "
        "changed. The generator is specified in the source rather than taken "
        "from the standard library, so the number reproduces across platforms.",
        "Ladder", "param", "ladder_k_bootstrap"};

    bool_default_map["write_patch_terrain_aspect"] = false;
    help_map["write_patch_terrain_aspect"] = {"bool","false",
        "Add the terrain aspect and its resultant length to the patch table.",
        "The existing `aspect_degrees` column is the direction the BANDED "
        "FIELD climbs, fitted as a plane through it - on a HAND field that is "
        "roughly `away from the channel`, not the direction the ground faces. "
        "This adds the topographic aspect as a CIRCULAR mean (the arithmetic "
        "mean of 350 and 10 degrees is 180, the opposite way) together with "
        "`aspect_r_bar`, the mean resultant length in [0, 1]. R_bar is the "
        "number that says whether the mean direction is worth quoting: 1 means "
        "every cell in the patch faces the same way, 0 means they face every "
        "way equally and the aspect is meaningless. Costs one polynomial fit.",
        "Outputs", "output", ""};

    float_default_map["patch_aspect_window_radius_metres"] = 0.0f;
    help_map["patch_aspect_window_radius_metres"] = {"float","0",
        "Window radius [m] of the polynomial fit the terrain aspect is "
        "measured over. 0 = five cells.",
        "Stated in cells by default because a length tuned on lidar is "
        "sub-pixel on 30 m COP30.",
        "Outputs", "output", "write_patch_terrain_aspect"};

    string_default_map["vertical_accuracy_algorithm_choice"] = "fixed";
    help_map["vertical_accuracy_algorithm_choice"] = {"string","fixed",
        "How E_v is arrived at: fixed or slope_dependent.",
        "Photogrammetric and InSAR DEMs degrade on slopes - TanDEM-X relative "
        "LE90 is roughly 2 m on gentle ground and 4 m on steep - so a single "
        "number for a raster spanning a valley floor and a 40-degree wall "
        "states an accuracy that is wrong in both places. Under "
        "slope_dependent each along-channel WINDOW reads one E_v off the "
        "median gradient of its own admitted cells. Per window, not per pixel: "
        "the resolvability floor constrains the separation between band "
        "CENTRES, and a band centre is a property of a partition rather than "
        "of a place, so it has no slope of its own. The value each window used "
        "is written to _ladder_windows.csv.",
        "Ladder", "param", ""};

    float_default_map["vertical_accuracy_gentle_metres"] = 0.0f;
    help_map["vertical_accuracy_gentle_metres"] = {"float","0.0",
        "E_v [m] for windows at or below the break gradient. "
        "slope_dependent only.",
        "For TanDEM-X this is the quoted relative LE90 on gentle ground, "
        "about 2 m.",
        "Ladder", "param", ""};

    float_default_map["vertical_accuracy_steep_metres"] = 0.0f;
    help_map["vertical_accuracy_steep_metres"] = {"float","0.0",
        "E_v [m] for windows above the break gradient. slope_dependent only.",
        "About 4 m for TanDEM-X. It must be the LARGER of the two or the mode "
        "is doing nothing useful.",
        "Ladder", "param", ""};

    float_default_map["vertical_accuracy_break_gradient"] = 0.20f;
    help_map["vertical_accuracy_break_gradient"] = {"float","0.2",
        "Gradient [m/m] separating the gentle and steep E_v.",
        "0.2 is a 20 percent grade, which is where the TanDEM-X accuracy "
        "specification splits. Replaces vertical_accuracy_break_degrees; an "
        "old value is converted with gradient = tan(angle).",
        "Ladder", "param", ""};
    string_default_map["legacy_vertical_accuracy_break_degrees"] = "";

    float_default_map["vertical_accuracy_gradient_window_radius_metres"] = 0.0f;
    help_map["vertical_accuracy_gradient_window_radius_metres"] = {"float","0",
        "Window radius [m] of the polynomial fit the window gradient is "
        "measured over. 0 = five cells. slope_dependent only.",
        "Stated in cells by default because a length tuned on lidar is "
        "sub-pixel on 30 m COP30. This is the same measurement "
        "lsdtt-flat-surface-finder makes for max_admission_gradient, and "
        "it is worth using the same window in both so the gradient a surface "
        "was admitted on is the gradient its E_v was read from.",
        "Ladder", "param", ""};



    float_default_map["ladder_reach_window_metres"] = 0.0f;
    help_map["ladder_reach_window_metres"] = {"float","0.0",
        "Along-channel window length for reach-relative breaks [m].",
        "0 = auto: the larger of 2 km and an eighth of each channel's length. "
        "Negative = one global partition over the whole population. Windowing "
        "lets a break level track downstream changes in incision instead of "
        "being averaged into a single basin-wide ladder. Ignored for the "
        "fitted_surface field, which has no along-channel coordinate.",
        "Ladder","param",""};

    int_default_map["ladder_min_window_pixels"] = 1000;
    help_map["ladder_min_window_pixels"] = {"int","1000",
        "Minimum admitted pixels for a window to find its own breaks.",
        "A sparser window falls back to the global partition rather than "
        "inventing structure from too few pixels.",
        "Ladder","param",""};

    int_default_map["ladder_hist_bins"] = 64;
    help_map["ladder_hist_bins"] = {"int","64",
        "Histogram bin count for the antimode break finder.",
        "Clamped to at least 16.",
        "Ladder","param",""};

    float_default_map["ladder_max_valley_fraction"] = 0.6f;
    help_map["ladder_max_valley_fraction"] = {"float","0.6",
        "Antimode acceptance: dip depth relative to the smaller flanking peak.",
        "A dip is accepted as a break when its smoothed density is below this "
        "fraction of the smaller of its two flanking modal peaks.",
        "Ladder","param",""};

    float_default_map["ladder_min_field_metres"] = 0.05f;
    help_map["ladder_min_field_metres"] = {"float","0.05",
        "Log floor for the antimode histogram [m].",
        "Field values below this are clamped before taking log10.",
        "Ladder","param",""};

    // --- Floodplain and terrace bounds (fixed_thresholds; clubb2017) ---
    lsdtt::declare_floodplain_threshold_algorithm_parameter(
        string_default_map, help_map, "preset");
    help_map["floodplain_threshold_algorithm_choice"][2] =
        "How the floodplain bounds are chosen: preset, fixed, qq_gaussian or "
        "cdf_divergence.";
    lsdtt::declare_qq_threshold_parameters(
        float_default_map, help_map, "floodplain",
        "floodplain_threshold_algorithm_choice");
    lsdtt::declare_max_floodplain_gradient_parameter(
        float_default_map, help_map, 0.1f, "floodplain_threshold_algorithm_choice");
    help_map["max_floodplain_gradient"][3] +=
        " The gradient bound of admission floodplain_thresholds, used as given "
        "under floodplain_threshold_algorithm_choice fixed.";

    string_default_map["floodplain_bound_algorithm_choice"] = "auto";
    help_map["floodplain_bound_algorithm_choice"] = {"string","auto",
        "Which bounds the floodplain admission uses: auto, hand_and_gradient "
        "or hand_only (admission floodplain_thresholds).",
        "hand_and_gradient: HAND below its bound AND gradient below its bound, "
        "as Clubb et al. (2017). hand_only: HAND alone - for low-relief "
        "landscapes, where the floodplain is cut into a surface that is "
        "itself nearly flat, and field-, road- and levee-scale relief puts "
        "much of the true floodplain above a data-chosen gradient bound. "
        "auto (default): hand_only when more than low_relief_min_area_fraction "
        "of the DEM's floodplain gradient is below low_relief_max_gradient, "
        "else hand_and_gradient. The fraction and the branch are in the run "
        "log, and the branch is recorded in _ingested_params.param so a rerun "
        "repeats it. Measured 2026-10-27 on four tiles (fraction below 0.1): "
        "Royalton 0.97 (gradient bound removed the floodplain), Mayoco 0.53, "
        "Salinas 0.27, Freiburg 0.17.",
        "Floodplain","param","admission_algorithm_choice"};

    float_default_map["low_relief_max_gradient"] = 0.1f;
    help_map["low_relief_max_gradient"] = {"float","0.1",
        "Gradient [m/m] below which a cell counts as low relief "
        "(floodplain_bound_algorithm_choice auto).",
        "Read on the floodplain gradient (floodplain_gradient_window_radius_"
        "metres). Gradient depends on resolution and window, so the fraction "
        "is too.",
        "Floodplain","param","floodplain_bound_algorithm_choice"};

    float_default_map["low_relief_min_area_fraction"] = 0.8f;
    help_map["low_relief_min_area_fraction"] = {"float","0.8",
        "Fraction of the DEM below low_relief_max_gradient above which auto "
        "drops the gradient bound.", "",
        "Floodplain","param","floodplain_bound_algorithm_choice"};

    float_default_map["floodplain_gradient_window_radius_metres"] = 0.0f;
    help_map["floodplain_gradient_window_radius_metres"] = {"float","0",
        "Window radius [m] of the polynomial fit the floodplain gradient is "
        "read over. 0 = three cells.",
        "LSDTopoTools2's surface_fitting_radius (30 m by default, three cells "
        "on 10 m data).",
        "Floodplain","param","admission_algorithm_choice"};

    bool_default_map["run_terrace_classification"] = true;
    help_map["run_terrace_classification"] = {"bool","true",
        "Find terraces as well as the floodplain (fixed_thresholds).",
        "false = the floodplain alone, as LSDTopoTools2 allowed: the terrace "
        "field (and its swath baseline, the most memory- and time-hungry "
        "stage on a large DEM) is never built and no terrace class is "
        "assigned. The floodplain mask is identical either way: it is decided "
        "before the terraces and independently of them.",
        "Floodplain","switch","level_algorithm_choice"};

    float_default_map["min_terrace_height_above_floodplain_metres"] = 0.0f;
    help_map["min_terrace_height_above_floodplain_metres"] = {"float","0",
        "How far above the floodplain bound a terrace starts [m] "
        "(fixed_thresholds). Must be >= 0.",
        "Terraces are ground whose terrace field exceeds max_floodplain_hand_"
        "metres plus this, so the two never overlap. LSDTopoTools2's "
        "minimum_terrace_height (whose window could overlap the floodplain). "
        "Different from min_terrace_hand_metres, the ladder's absolute bound.",
        "Floodplain","param","level_algorithm_choice"};

    float_default_map["max_terrace_height_metres"] = 0.0f;
    help_map["max_terrace_height_metres"] = {"float","0",
        "Highest a terrace may stand on the terrace field [m] "
        "(fixed_thresholds). 0 = no ceiling.",
        "With no ceiling the swath half-width and the gradient bound are the "
        "only limits on how high terraces reach.",
        "Floodplain","param","level_algorithm_choice"};

    bool_default_map["fill_floodplain_holes"] = false;
    help_map["fill_floodplain_holes"] = {"bool","false",
        "Fill enclosed holes in the floodplain (fixed_thresholds).",
        "LSDTopoTools2's fill_floodplain (lsdtt-valley-metrics), done with "
        "binary_fill_holes rather than OpenCV: a background region not "
        "connected to the outside of the valid DEM is a hole. Default false, "
        "but TRUE under the clubb2017 preset unless the parameter file names "
        "it; the value used is in the run log's stage table.",
        "Floodplain","param","level_algorithm_choice"};

    string_default_map["floodplain_second_pass_algorithm_choice"] = "none";
    help_map["floodplain_second_pass_algorithm_choice"] = {"string","none",
        "A second HAND bound read off the floodplain itself: none or "
        "qq_upper_departure (fixed_thresholds).",
        "qq_upper_departure: the QQ plot of HAND over the first-pass "
        "floodplain cells, against the same Gaussian line, and the HIGHEST "
        "quantile at which the data still agree with it (the mirror of the "
        "Clubb rule, same floodplain_qq_* settings). Above it HAND departs "
        "upward: wide flats standing above the valley floor - fans at "
        "tributary mouths, upland basins. The floodplain is then rebuilt "
        "with that HAND bound (same gradient bound, channel contact and size "
        "filter) and its holes are always filled; terraces are measured from "
        "the second-pass bound. Not part of Clubb et al. (2017). The first "
        "pass is still written (_floodplain_first_pass_mask.tif, "
        "_floodplain_first_pass_polygons) and the second QQ table is "
        "<stem>_floodplain_hand_qq_second_pass.csv.",
        "Floodplain","param","level_algorithm_choice"};

    bool_default_map["write_floodplain_polygons"] = false;
    help_map["write_floodplain_polygons"] = {"bool","false",
        "Write the floodplain as polygons (fixed_thresholds).",
        "Output file: <stem>_floodplain_polygons.<fgb|gpkg> (write_vector_format; "
        "one polygon per connected floodplain component, 8-connected), and "
        "<stem>_floodplain_first_pass_polygons.<fgb|gpkg> when a second pass "
        "ran.",
        "Outputs","output","level_algorithm_choice"};

    bool_default_map["write_floodplain_threshold_csv"] = true;
    help_map["write_floodplain_threshold_csv"] = {"bool","true",
        "Write the table each data-chosen floodplain bound was read from.",
        "Output files: <stem>_floodplain_hand_qq.csv and "
        "<stem>_floodplain_gradient_qq.csv (qq_gaussian), or "
        "_adaptive_cdf.csv (cdf_divergence). A threshold is only "
        "interpretable against the distribution it was cut from.",
        "Outputs","output","floodplain_threshold_algorithm_choice"};

    bool_default_map["write_terrace_field_raster"] = false;
    help_map["write_terrace_field_raster"] = {"bool","false",
        "Write the terrace reference field (fixed_thresholds).",
        "Output file: <stem>_terrace_field.tif - swath-baseline relief, or "
        "HAND under same.",
        "Outputs","output","terrace_reference_field_algorithm_choice"};

    bool_default_map["write_floodplain_gradient_raster"] = false;
    help_map["write_floodplain_gradient_raster"] = {"bool","false",
        "Write the gradient raster (m/m) the floodplain bound was applied to.",
        "Output file: <stem>_floodplain_gradient.tif.",
        "Outputs","output","admission_algorithm_choice"};

    // --- Classification (hand field) ---
    // Shared name with lsdtt-hillslope-channel-coupling's floodplain rule.
    lsdtt::declare_max_floodplain_hand_parameter(float_default_map, help_map,
                                                 1.0f);
    help_map["max_floodplain_hand_metres"][3] =
        "Here an optional absolute override on the lowest band's upper bound. "
        "Under band-then-label the ladder already decides where one level ends "
        "and the next begins; this is a hard cap for cases where you know the "
        "floodplain height independently. Only meaningful for the hand field. "
        "Replaces floodplain_max_height_metres.";

    float_default_map["min_terrace_hand_metres"] = 0.0f;
    help_map["min_terrace_hand_metres"] = {"float","0.0",
        "Lower HAND bound of the terrace class [m].",
        "As max_floodplain_hand_metres, at the other end. Only meaningful "
        "for the hand field. Replaces terrace_min_height_metres.",
        "Classification","param",""};

    float_default_map["channel_adjacency_radius_metres"] = 0.0f;
    help_map["channel_adjacency_radius_metres"] = {"float","0.0",
        "Distance within which a surface counts as channel-adjacent [m].",
        "0 = touching only.",
        "Classification","param",""};

    float_default_map["valley_corridor_half_width_metres"] = 5000.0f;
    help_map["valley_corridor_half_width_metres"] = {"float","5000.0",
        "Valley-corridor half-width [m]. 0 or negative = no limit.",
        "A terrace must have at least one pixel within this Euclidean "
        "distance of an in-scope channel. This is what keeps the regional "
        "upland out of the terrace class: an interfluve is also flat and also "
        "elevated above the trunk channel, but it is not a fluvial terrace.",
        "Classification","param",""};

    bool_default_map["clip_corridor_per_pixel"] = true;
    help_map["clip_corridor_per_pixel"] = {"bool","true",
        "Apply the corridor limit per pixel as well as per surface.",
        "Stops a surface that merely touches the valley edge from carrying its "
        "class across the whole interfluve.",
        "Classification","param","valley_corridor_half_width_metres"};

    float_default_map["max_band_spread_metres"] = 50.0f;
    help_map["max_band_spread_metres"] = {"float","50.0",
        "QC tripwire on the field spread within one surface [m]. 0 disables.",
        "A surface whose robust field spread still exceeds this cannot take "
        "the terrace class. Under band-then-label a band's spread is bounded "
        "by its own breaks, so this should almost never fire; it is retained "
        "as a tripwire, not as a mechanism.",
        "Classification","param",""};

    // --- Flow routing ---
    lsdtt::declare_fill_before_routing_parameter(bool_default_map, help_map);

    float_default_map["fill_min_slope_for_routing"] = 0.0001f;
    help_map["fill_min_slope_for_routing"] = {"float","0.0001",
        "Minimum imposed slope when filling for routing [m/m].",
        "Used only when fill_before_routing is true. inline_gradient "
        "(default): imposed during the Barnes et al. 2014 Priority-Flood. "
        "barnes2014b_gradient: the minimum gradient across flats. "
        "barnes2014b_float_step: ignored.",
        "Flow routing","param","fill_before_routing"};

    lsdtt::declare_fill_flat_resolution_parameter(
        string_default_map, help_map, "", "fill_min_slope_for_routing");

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Edge boundary codes, north east south west.",
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        "Flow routing","param",""};

    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    // --- Channel scope: source and outlet points ---
    bool_default_map["break_flats_with_channels"] = false;
    help_map["break_flats_with_channels"] = {"bool","false",
        "Remove in-scope channel pixels from the admitted set before banding.",
        "Default FALSE, changed from the old flat-surface-finder behaviour. It "
        "used to be needed to stop one surface spanning a channel; banding by "
        "height does that on its own, since the two banks of an incised channel "
        "are not in the same band as the channel. All it does now is keep the "
        "channel itself out of a surface, which matters for area statistics and "
        "little else.",
        "Channel scope","switch",""};

    bool_default_map["connect_surfaces_to_channels"] = false;
    help_map["connect_surfaces_to_channels"] = {"bool","false",
        "Record which in-scope channel each surface is linked to.",
        "Adds the connection diagnostics used to attribute a surface to a "
        "reach. The patch table already carries source_id from the height "
        "reference's donor, so this is for the cases where the linkage itself "
        "is the question.",
        "Channel scope","switch",""};

    string_default_map["drainage_injection_points_fname"] = "";
    help_map["drainage_injection_points_fname"] = {"string","",
        "Vector of points injecting extra drainage area, in read_path.",
        "For a channel entering the DEM from outside it, where the "
        "contributing-pixel threshold alone would not start a channel at the "
        "boundary. Empty = no injection.",
        "Channel scope","io",""};

    string_default_map["source_points_fname"] = "";
    help_map["source_points_fname"] = {"string","",
        "Channel source points (CSV or vector), in read_path.",
        "Empty = every channel above the threshold is in scope. Supplying "
        "sources restricts the analysis to those channels, which is how you "
        "measure terraces along one trunk rather than every rill in the DEM.",
        "Channel scope","io",""};

    string_default_map["source_points_input_format"] = "auto_detect";
    help_map["source_points_input_format"] = {"string","auto_detect",
        "Source point format: auto_detect, csv or ogr_vector_layer.",
        "auto_detect infers from the file extension.",
        "Channel scope","io","source_points_fname"};

    string_default_map["source_points_layer_name"] = "";
    help_map["source_points_layer_name"] = {"string","",
        "Layer name within a source point vector file.",
        "Empty = the first layer.",
        "Channel scope","io","source_points_fname"};

    string_default_map["source_points_latitude_field_name"] = "latitude";
    help_map["source_points_latitude_field_name"] = {"string","latitude",
        "Latitude column name in a source point CSV.",
        "CSV points are read as EPSG:4326 and transformed to the DEM CRS.",
        "Channel scope","io","source_points_fname"};

    string_default_map["source_points_longitude_field_name"] = "longitude";
    help_map["source_points_longitude_field_name"] = {"string","longitude",
        "Longitude column name in a source point CSV.",
        "See source_points_latitude_field_name.",
        "Channel scope","io","source_points_fname"};

    string_default_map["outlet_points_fname"] = "";
    help_map["outlet_points_fname"] = {"string","",
        "Channel outlet points (CSV or vector), in read_path.",
        "Truncates each traced channel at its matched outlet, so the analysis "
        "covers a reach rather than everything down to the DEM edge.",
        "Channel scope","io",""};

    string_default_map["outlet_points_input_format"] = "auto_detect";
    help_map["outlet_points_input_format"] = {"string","auto_detect",
        "Outlet point format: auto_detect, csv or ogr_vector_layer.",
        "As source_points_input_format.",
        "Channel scope","io","outlet_points_fname"};

    string_default_map["outlet_points_layer_name"] = "";
    help_map["outlet_points_layer_name"] = {"string","",
        "Layer name within an outlet point vector file.",
        "Empty = the first layer.",
        "Channel scope","io","outlet_points_fname"};

    string_default_map["outlet_points_latitude_field_name"] = "latitude";
    help_map["outlet_points_latitude_field_name"] = {"string","latitude",
        "Latitude column name in an outlet point CSV.",
        "See source_points_latitude_field_name.",
        "Channel scope","io","outlet_points_fname"};

    string_default_map["outlet_points_longitude_field_name"] = "longitude";
    help_map["outlet_points_longitude_field_name"] = {"string","longitude",
        "Longitude column name in an outlet point CSV.",
        "See source_points_latitude_field_name.",
        "Channel scope","io","outlet_points_fname"};

    bool_default_map["auto_detect_outlet_points"] = true;
    help_map["auto_detect_outlet_points"] = {"bool","true",
        "Infer source and outlet roles from a single combined point file.",
        "A supplied point that snaps strictly downstream of another is treated "
        "as that one's outlet. Ignored when outlet_points_fname is given, "
        "because then the roles are stated directly.",
        "Channel scope","switch","source_points_fname"};

    bool_default_map["snap_sources_to_channel_threshold"] = true;
    help_map["snap_sources_to_channel_threshold"] = {"bool","true",
        "Snap supplied source points onto the threshold channel network.",
        "A hand-digitised point rarely lands exactly on the routed channel.",
        "Channel scope","switch","source_points_fname"};

    int_default_map["source_snap_max_search_radius_pixels"] = 20;
    help_map["source_snap_max_search_radius_pixels"] = {"int","20",
        "Search radius when snapping a source point [pixels].",
        "A point that finds no threshold channel within this radius is "
        "reported as unsnapped in the point snap diagnostics.",
        "Channel scope","param","snap_sources_to_channel_threshold"};

    float_default_map["outlet_snap_max_distance_metres"] = 500.0f;
    help_map["outlet_snap_max_distance_metres"] = {"float","500.0",
        "Maximum distance an outlet may be moved onto a traced channel [m].",
        "Outlets snap onto the traced flow paths, not onto the channel "
        "network, so a slightly misplaced outlet still truncates the intended "
        "channel.",
        "Channel scope","param","outlet_points_fname"};

    // --- Outputs ---
    bool_default_map["write_reference_field_raster"] = true;
    help_map["write_reference_field_raster"] = {"bool","true",
        "Write the field the ladder actually ran on.",
        "Output file: <stem>_reference_field.tif. Whichever field was "
        "selected, so a run is self-documenting.",
        "Outputs","switch",""};

    bool_default_map["write_hand_raster"] = false;
    help_map["write_hand_raster"] = {"bool","false",
        "Write the height-above-channel raster.",
        "Output file: <stem>_hand.tif. Feed it back as hand_fname to iterate "
        "on the ladder without recomputing HAND.",
        "Outputs","switch",""};

    bool_default_map["write_band_index_raster"] = true;
    help_map["write_band_index_raster"] = {"bool","true",
        "Write the per-pixel band index.",
        "Output file: <stem>_band_index.tif (Int16, nodata -1). The "
        "authoritative statement of which level each pixel is on; the band is "
        "NOT packed into the patch label.",
        "Outputs","switch",""};

    bool_default_map["write_banded_labels_raster"] = true;
    help_map["write_banded_labels_raster"] = {"bool","true",
        "Write the per-pixel patch label.",
        "Output file: <stem>_banded_labels.tif (UInt32). A global running "
        "counter 1..N assigned band by band; 0 is background. Join it to "
        "<stem>_surface_patches.csv on the label.",
        "Outputs","switch",""};

    bool_default_map["write_valley_mask_raster"] = true;
    help_map["write_valley_mask_raster"] = {"bool","true",
        "Write band 0 as a valley mask.",
        "Output file: <stem>_valley_mask.tif. This is the input "
        "lsdtt-valley-centreline expects as valley_mask_fname. With "
        "ladder_level_count 2 it is exactly the valley/wall cut.",
        "Outputs","switch",""};

    // --- Valley mask cleaning ---
    float_default_map["valley_mask_closing_width_fraction"] = 0.0f;
    help_map["valley_mask_closing_width_fraction"] = {"float","0.0",
        "Set the closing radius as a fraction of the measured valley width. "
        "0 = use valley_mask_closing_radius_metres instead.",
        "A radius in metres does not travel: a value that seals the slivers in "
        "a 400 m floodplain erases a 100 m one, so it has to be re-guessed for "
        "every DEM. This measures the width instead - cross sections cast "
        "perpendicular to each channel, the CUMULATIVE length of valley cells "
        "along each (so a bar or an island does not count as floor), median "
        "across sections - and takes the radius as this fraction of it. About "
        "0.2 is a reasonable starting point. Overrides "
        "valley_mask_closing_radius_metres when set.",
        "Outputs","param","clean_valley_mask"};

    float_default_map["valley_width_section_spacing_metres"] = 200.0f;
    help_map["valley_width_section_spacing_metres"] = {"float","200.0",
        "Spacing of the cross sections used to measure the valley width [m].",
        "Only used when valley_mask_closing_width_fraction is set. Closer "
        "spacing gives a better-sampled median at more cost.",
        "Outputs","param","valley_mask_closing_width_fraction"};

    float_default_map["valley_width_max_half_length_metres"] = 3000.0f;
    help_map["valley_width_max_half_length_metres"] = {"float","3000.0",
        "How far each cross section looks either side before giving up [m].",
        "A section reaching this on BOTH sides is discarded rather than "
        "counted, because it is measuring the search limit and not the valley. "
        "Set it comfortably wider than the widest valley expected.",
        "Outputs","param","valley_mask_closing_width_fraction"};

    float_default_map["valley_mask_closing_radius_metres"] = 0.0f;
    help_map["valley_mask_closing_radius_metres"] = {"float","0.0",
        "Seal narrow gaps in the valley mask with a morphological closing [m]. "
        "0 = off.",
        "A dilation followed by an erosion of the same radius. It bridges thin "
        "OPEN-ENDED gaps - the parallel slivers left by scroll bars, levees, "
        "abandoned channels and field boundaries lying a little above the "
        "surrounding floor - which hole filling cannot touch, because a finger "
        "open at one end is not an enclosed hole. Because the erosion undoes "
        "the dilation, the outer extent is left about where it was rather than "
        "grown. Set it to roughly the width of the widest sliver you want "
        "closed; too large and it will bridge across a genuine terrace riser.",
        "Outputs","param","clean_valley_mask"};

    bool_default_map["fill_valley_mask_holes"] = true;
    help_map["fill_valley_mask_holes"] = {"bool","true",
        "Fill enclosed holes in the valley mask before writing it.",
        "A hole in the valley floor - a rise the ladder put in a higher band, a "
        "patch the size filter dropped - is a hole in the DOMAIN the centreline "
        "is traced through, and the least-cost walk goes around it. The "
        "centreline then hugs one wall instead of running down the middle. "
        "Filling enclosed holes is what makes the mask a valley rather than a "
        "valley with islands. Only affects <stem>_valley_mask.tif; the band "
        "index and the patch table keep the holes, because there the rise is a "
        "real observation.",
        "Outputs","switch","write_valley_mask_raster"};

    // --- Forced valley floor: state the threshold instead of finding it ---
    float_default_map["force_valley_hand_metres"] = 0.0f;
    help_map["force_valley_hand_metres"] = {"float","0.0",
        "Define the valley floor directly as height above channel below this "
        "[m], instead of using the ladder's break. 0 = use the ladder.",
        "The ladder finds a break from the data; sometimes you already know the "
        "answer - from field mapping, a flood stage, a published terrace "
        "height - and want it applied rather than rediscovered. Combined with "
        "force_valley_max_distance_metres this defines the valley floor "
        "outright: HAND below the threshold AND within the distance. The band "
        "index, labels and patch table are UNAFFECTED and still come from the "
        "ladder; only <stem>_valley_mask.tif is forced, so the two can be "
        "compared.",
        "Valley floor","param",""};

    float_default_map["force_valley_max_distance_metres"] = 0.0f;
    help_map["force_valley_max_distance_metres"] = {"float","0.0",
        "Euclidean distance to an in-scope channel bounding the forced valley "
        "floor [m]. 0 = no distance limit.",
        "The companion to force_valley_hand_metres, and the reason the pair is "
        "useful: a HAND threshold alone admits every low patch in the DEM, "
        "including ground belonging to other drainages entirely. The distance "
        "limit is what confines it to this channel's valley. Usable on its own "
        "as well, if the ladder is finding the break but reaching too far.",
        "Valley floor","param",""};

    bool_default_map["clean_valley_mask"] = true;
    help_map["clean_valley_mask"] = {"bool","true",
        "Tidy band 0 before writing it as the valley mask.",
        "Band 0 straight off the ladder is speckled and contains flat patches "
        "that have nothing to do with the channel, so it is not a drop-in for "
        "lsdtt-valley-centreline until it has been cleaned the way that tool "
        "used to clean its own. Erode, drop patches that do not reach the "
        "channel, then regrow. Only affects <stem>_valley_mask.tif; the band "
        "index, the labels and the patch table are untouched.",
        "Outputs","switch","write_valley_mask_raster"};

    int_default_map["valley_dilation_pixels"] = 0;
    help_map["valley_dilation_pixels"] = {"int","0",
        "Erosion / regrow radius used when cleaning the valley mask [pixels]. "
        "0 = no erosion.",
        "DEFAULT 0, and erosion is usually the wrong tool here. It exists to "
        "SEVER thin threads of flat ground joining one valley floor to a "
        "neighbour's across a col, so that connected-component filtering can "
        "tell them apart - it is a means of making connectivity discriminating, "
        "not a way to shrink the mask. The cost is that a reach narrower than "
        "2 * this + 1 pixels erodes to nothing and cannot regrow: on the "
        "Royalton 10 m DEM a value of 4 (64 m) destroyed 30% of a valley floor "
        "that connectivity alone had got right. Raise it only when two valleys "
        "are genuinely merging through a col.",
        "Outputs","param","clean_valley_mask"};

    float_default_map["valley_max_detour_ratio"] = 0.0f;
    help_map["valley_max_detour_ratio"] = {"float","0.0",
        "Maximum along-ground path to the channel, as a multiple of the "
        "straight-line distance.",
        "A flat that sits NEAR the channel but is only reachable the long way "
        "round is a different landform, not this channel's valley floor. 0 or "
        "negative disables the test, which is the DEFAULT: on a meandering "
        "floodplain the along-ground path around a bend legitimately exceeds "
        "the straight line several times over, and the test removes real "
        "valley floor. Turn it on where flats from separate drainages are "
        "connected by a long corridor.",
        "Outputs","param","clean_valley_mask"};

    bool_default_map["restrict_valley_to_line_of_sight"] = false;
    help_map["restrict_valley_to_line_of_sight"] = {"bool","false",
        "Require a straight line from each valley pixel to the channel to stay "
        "inside the mask.",
        "A STRONG filter, written for separating one valley floor from a "
        "neighbouring one across a broad col, where connectivity alone is not "
        "enough. Default FALSE here: on a meandering floodplain a great deal of "
        "genuine valley floor has no straight line to the channel that stays "
        "inside the mask. On the Royalton 10 m DEM it removed 82% of band 0 on "
        "its own, with every other filter disabled - which made the valley mask "
        "SMALLER than the floodplain class computed from the same field. Turn "
        "it on where two valleys meet around a col and connectivity alone "
        "merges them.",
        "Outputs","param","clean_valley_mask"};

    bool_default_map["restrict_valley_to_channel_connected"] = true;
    help_map["restrict_valley_to_channel_connected"] = {"bool","true",
        "Drop valley patches with no pixel on or beside the channel.",
        "The valley floor of a channel touches that channel. A patch that "
        "touches none is flat ground somewhere else.",
        "Outputs","param","clean_valley_mask"};

    bool_default_map["write_surface_class_raster"] = true;
    help_map["write_surface_class_raster"] = {"bool","true",
        "Write the floodplain / terrace class raster.",
        "Output file: <stem>_surface_class.tif. Only written for the hand "
        "field: floodplain and terrace are statements about height above a "
        "channel, and mean nothing on a fitted-surface residual.",
        "Outputs","switch",""};

    bool_default_map["write_floodplain_mask_raster"] = false;
    help_map["write_floodplain_mask_raster"] = {"bool","false",
        "Write the floodplain class on its own.",
        "Output file: <stem>_floodplain_mask.tif. hand field only.",
        "Outputs","switch",""};

    bool_default_map["write_method_record_csv"] = true;
    help_map["write_method_record_csv"] = {"bool","true",
        "Write the method record: every stage, every data-chosen bound, and "
        "the measured counts.",
        "Output file: <stem>_method_record.csv (method, item, value, source). "
        "source says whether a stage came from the preset or the parameter "
        "file, or that a value was measured. With several methods in "
        "fluvial_surface_algorithm_choice (comma-separated, no spaces) each "
        "has its own, and <write_fname>_method_summary.csv, "
        "_method_comparison.csv (per pair: floodplain and terrace "
        "intersection-over-union, agreement on classified cells, the 3 x 3 "
        "confusion in cells and m^2) and _method_agreement_floodplain.tif / "
        "_terrace.tif (UInt8, bit k = method k in list order) are written.",
        "Outputs","output",""};

    bool_default_map["write_surface_class_polygons"] = false;
    help_map["write_surface_class_polygons"] = {"bool","false",
        "Write the floodplain and terrace classes as polygons (hand field).",
        "Output file: <stem>_surface_class_polygons.<fgb|gpkg> "
        "(write_vector_format): one polygon per connected area of one class "
        "(8-connected), field class (1 floodplain, 2 terrace), label_id, "
        "pixel_count, area_m2. Outer boundaries only.",
        "Outputs","output",""};

    bool_default_map["write_surface_patch_polygons"] = false;
    help_map["write_surface_patch_polygons"] = {"bool","false",
        "Write every labelled patch as a polygon, with its level.",
        "Output file: <stem>_surface_patch_polygons.<fgb|gpkg> "
        "(write_vector_format). Fields: label_id (joins "
        "_surface_patches.csv and _banded_labels.tif), pixel_count, area_m2, "
        "band (the ladder level, or 0 floodplain / 1 terrace under "
        "fixed_thresholds). Outer boundaries only: a patch enclosing another "
        "level's patch is drawn without the hole.",
        "Outputs","output",""};

    bool_default_map["write_surface_patches_csv"] = true;
    help_map["write_surface_patches_csv"] = {"bool","true",
        "Write one row per labelled patch.",
        "Output file: <stem>_surface_patches.csv with label, band, class, "
        "source_id, n_pixels, area_m2, median_field_m, aspect_degrees and "
        "window_centre_m. The join key for everything downstream.",
        "Outputs","switch",""};

    bool_default_map["write_surface_class_csv"] = true;
    help_map["write_surface_class_csv"] = {"bool","true",
        "Write one row per classification record.",
        "Output file: <stem>_surface_class.csv. hand field only.",
        "Outputs","switch",""};

    bool_default_map["write_reference_field_profile_csv"] = false;
    help_map["write_reference_field_profile_csv"] = {"bool","false",
        "Write a per-pixel table of the reference field along the channel.",
        "Output file: <stem>_reference_field_profile.csv with x, y, field_m, "
        "channel_flow_distance_m, source_id, class, label and band. One row per "
        "labelled pixel: the input for a height-against-distance profile "
        "figure, and the table that proves the raster labels join to the patch "
        "and class tables. Requires the hand field with a computed HAND, since "
        "the along-channel distance comes from the pixel's HAND donor node.",
        "Outputs","switch",""};

    bool_default_map["write_nearest_source_id_raster"] = false;
    help_map["write_nearest_source_id_raster"] = {"bool","false",
        "Write which channel each pixel's height is measured to.",
        "Output file: <stem>_nearest_source_id.tif (Int32, -1 off-network). "
        "The value is the source id of the pixel's HAND DONOR - the channel "
        "node it actually drains to - not the nearest channel by straight-line "
        "distance, so it partitions the DEM the same way the height reference "
        "does. Use it to look at one channel's results in isolation. Needs a "
        "computed HAND: a supplied hand_fname carries no donor nodes.",
        "Outputs","switch",""};

    bool_default_map["write_ladder_break_metric_csv"] = false;
    help_map["write_ladder_break_metric_csv"] = {"bool","false",
        "Write the break DECISION, not just the break.",
        "Output file: <stem>_ladder_break_metric.csv, one row per histogram bin "
        "per window: the raw count, the smoothed density the separation ratio "
        "is measured on, which detected gap the bin falls in, and whether a "
        "candidate or the chosen break lands there. The chosen break on a raw "
        "histogram looks arbitrary, and whether it beat several candidates or "
        "was the only one on offer cannot be recovered from the answer alone.",
        "Outputs","switch",""};

    bool_default_map["write_ladder_window_id_raster"] = false;
    help_map["write_ladder_window_id_raster"] = {"bool","false",
        "Write which along-channel window each pixel's breaks came from.",
        "Output file: <stem>_ladder_window_id.tif (Int16, -1 where the global "
        "partition was used or the pixel is not banded). The value is the "
        "window's ordinal within its channel, counting downstream - a reach "
        "number. Use it to trace a break that looks wrong back to the ground "
        "that produced it.",
        "Outputs","switch",""};

    bool_default_map["write_ladder_histogram_csv"] = false;
    help_map["write_ladder_histogram_csv"] = {"bool","false",
        "Write the reference-field density each window's breaks were chosen "
        "from.",
        "Output file: <stem>_ladder_histogram.csv, one row per histogram bin "
        "per window, with the bin containing each break flagged. THIS is the "
        "diagnostic that says whether a break separates two populations or "
        "merely cuts through one - a deep, wide dip between two modes means two "
        "surfaces; a shallow scoop through one broad mode means one surface and "
        "an arbitrary cut. Neither the resolvability floor nor the chosen break "
        "value can tell you which you have. Plot it with lsdviz3-ladder-density.",
        "Outputs","switch",""};

    int_default_map["ladder_hist_plot_bins"] = 80;
    help_map["ladder_hist_plot_bins"] = {"int","80",
        "Bin count for the per-window density in write_ladder_histogram_csv.",
        "Linear in field units, not log10, so it can be read directly against "
        "the breaks. Distinct from ladder_hist_bins, which is the log10 "
        "histogram the antimode finder searches for dips.",
        "Outputs","param","write_ladder_histogram_csv"};

    bool_default_map["write_ladder_windows_csv"] = true;
    help_map["write_ladder_windows_csv"] = {"bool","true",
        "Write the per-window resolvability report.",
        "Output file: <stem>_ladder_windows.csv with the level count chosen, "
        "the status, the breaks, the band centres and the gap that had to be "
        "cleared. READ THIS FIRST: what may be claimed about a ladder depends "
        "on whether its windows resolved, and that is not visible in any "
        "raster.",
        "Outputs","switch",""};

    bool_default_map["write_channel_scope_raster"] = false;
    help_map["write_channel_scope_raster"] = {"bool","false",
        "Write the in-scope channel network as a raster.",
        "Output file: <stem>_channel_scope.tif (UInt32, source_id + 1).",
        "Outputs","switch",""};

    bool_default_map["write_channel_scope_nodes_csv"] = false;
    help_map["write_channel_scope_nodes_csv"] = {"bool","false",
        "Write the in-scope channel nodes as a table.",
        "Output file: <stem>_channel_scope_nodes.csv.",
        "Outputs","switch",""};

    bool_default_map["write_point_snap_csv"] = false;
    help_map["write_point_snap_csv"] = {"bool","false",
        "Write where each supplied point snapped to.",
        "Outputs: <stem>_point_snaps.csv plus raw and snapped FlatGeobuf "
        "point layers, so both can be overlaid in one QC figure.",
        "Outputs","switch",""};

    bool_default_map["write_channel_network"] = false;
    help_map["write_channel_network"] = {"bool","false",
        "Write the channel network as vector layers.",
        "Master switch for the three layers below. The network this tool "
        "measured its surfaces against is written here so you do not have to "
        "chain a second CLI to see it.",
        "Outputs","switch",""};

    bool_default_map["write_channel_network_lines"] = true;
    help_map["write_channel_network_lines"] = {"bool","true",
        "Write channel links as lines.",
        "Lines are the compact default representation of the network used to fit "
        "the surfaces; per-node points remain opt-in. Output file: "
        "<stem>_channel_network_lines.fgb.",
        "Outputs","switch","write_channel_network"};

    bool_default_map["write_channel_network_points"] = false;
    help_map["write_channel_network_points"] = {"bool","false",
        "Write channel nodes as points.",
        "Output file: <stem>_channel_network_points.fgb.",
        "Outputs","switch","write_channel_network"};

    bool_default_map["write_junction_network_points"] = false;
    help_map["write_junction_network_points"] = {"bool","false",
        "Write junctions as points.",
        "Output file: <stem>_junction_network_points.fgb.",
        "Outputs","switch","write_channel_network"};

    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);
    lsdtt::declare_true_north_bearing_parameters(bool_default_map, help_map);

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    // Records the run: closing message, run summary file, failure box.
    lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp(cli_args.path, cli_args.fname);

    // Renames. Since 2026-10-26 the admission and floodplain names are shared
    // with lsdtt-flat-surface-finder, lsdtt-mountain-front and
    // lsdtt-hillslope-channel-coupling (one name per aim), bounds take min_ /
    // max_ prefixes, and gradients are m/m.
    pp.set_deprecated_aliases({
        {"source_snap_max_radius_pixels", "source_snap_max_search_radius_pixels"},
        {"fill_dem_for_flow_routing", "fill_before_routing"},
        {"channel_threshold_contributing_pixels",
         "threshold_contributing_pixels"},
        {"fallback_relief_threshold_metres", "max_relief_metres"},
        {"fallback_slope_threshold", "max_admission_gradient"},
        {"min_component_pixels", "min_band_component_pixels"},
        {"vertical_accuracy_mode", "vertical_accuracy_algorithm_choice"},
        {"vertical_accuracy_slope_window_radius_metres",
         "vertical_accuracy_gradient_window_radius_metres"},
        {"floodplain_max_height_metres", "max_floodplain_hand_metres"},
        {"terrace_min_height_metres", "min_terrace_hand_metres"} });
    pp.set_deprecated_numeric_aliases({
        {"patch_aspect_window_metres",
         {"patch_aspect_window_radius_metres", 0.5, 0.0, true}},
        {"vertical_accuracy_slope_window_metres",
         {"vertical_accuracy_gradient_window_radius_metres", 0.5, 0.0, true}} });
    // Old values that cannot be carried by a plain or affine alias: converted
    // after parsing (a window in cells needs the cell size; degrees need a
    // tangent; the old admission selector only applied without a flat mask).
    pp.set_deprecated_compatibility_aliases({
        {"fallback_relief_window_radius_pixels",
         {"legacy_relief_window_radius_pixels", "relief_window_radius_metres",
          "The radius in cells is converted to metres with the DEM cell size; "
          "the radius used is in the run log."}},
        {"fallback_relief_window_pixels",
         {"legacy_relief_window_radius_pixels", "relief_window_radius_metres",
          "The radius in cells is converted to metres with the DEM cell size; "
          "the radius used is in the run log."}},
        {"fallback_admission_algorithm_choice",
         {"legacy_fallback_admission_algorithm_choice",
          "admission_algorithm_choice",
          "It applies only when flat_mask_fname is empty, as it always did; "
          "its `slope` is now `gradient`."}},
        {"vertical_accuracy_break_degrees",
         {"legacy_vertical_accuracy_break_degrees",
          "vertical_accuracy_break_gradient",
          "Converted with gradient = tan(angle)."}} });

    try
    {
        pp.parse_all_parameters(float_default_map, int_default_map,
                                bool_default_map, string_default_map);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto ip = pp.get_int_parameters();
    auto sf = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter resolution and validation
    // -------------------------------------------------------------------------
    lsdtt::raster_output_format out_format = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt = sf["write_output_format"];
        if (fmt == "gtiff" || fmt == "geotiff")
            out_format = lsdtt::raster_output_format::gtiff;
        else if (fmt == "same_as_input")
            out_format = lsdtt::raster_output_format::same_as_input;
    }

    // -------------------------------------------------------------------------
    // Method preset -> stage choices
    //
    // Each stage parameter defaults to "preset" and is filled from the preset
    // unless the parameter file named it. The table printed below says which,
    // so a run can be read back without the parameter file.
    // -------------------------------------------------------------------------
    // Methods to run (plan 7b): one preset, or several, comma-separated with
    // no spaces (ladder,clubb2017). Each runs the whole pipeline below under
    // its own stem, <write_fname>_<method>; stage parameters named in the file
    // apply to every one. One method keeps the plain <write_fname> stem.
    std::vector<std::string> method_list;
    {
        const std::string& ml = sf["fluvial_surface_algorithm_choice"];
        std::size_t start = 0;
        while (start <= ml.size())
        {
            const std::size_t comma = ml.find(',', start);
            const std::string m = ml.substr(start, comma == std::string::npos
                                                       ? std::string::npos
                                                       : comma - start);
            if (m.empty())
            {
                std::cerr << "FATAL: fluvial_surface_algorithm_choice '" << ml
                          << "' has an empty entry.\n";
                return 1;
            }
            if (std::find(method_list.begin(), method_list.end(), m) !=
                method_list.end())
            {
                std::cerr << "FATAL: fluvial_surface_algorithm_choice lists '"
                          << m << "' twice.\n";
                return 1;
            }
            method_list.push_back(m);
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
        if (method_list.size() > 8)
        {
            std::cerr << "FATAL: at most 8 methods in one run (the agreement "
                         "rasters hold one bit per method).\n";
            return 1;
        }
        if (std::find(method_list.begin(), method_list.end(),
                      "cross_section_levels") != method_list.end() &&
            sf["supplied_levels_fname"].empty() &&
            (sf["level_algorithm_choice"] == "preset" ||
             sf["level_algorithm_choice"] == "supplied_levels"))
        {
            std::cerr << "FATAL: cross_section_levels needs supplied_levels_fname "
                         "(the levels CSV from lsdtt-valley-centreline).\n";
            return 1;
        }
    }

    // What each method's run hands back for the comparison and the summary.
    struct method_run_output
    {
        std::vector<std::uint8_t> cls;     // surface_class per cell (hand only)
        std::vector<std::uint8_t> valid;
        lsdtt::lsdtt_raster_meta  meta;
        std::vector<std::pair<std::string, std::string>> summary;
    };

    auto run_method = [&](const std::string& method, const std::string& stem,
                          const std::string& fname_stem,
                          method_run_output& run_out) -> int
    {
    if (method != "ladder" && method != "clubb2017" &&
        method != "cross_section_levels")
    {
        std::cerr << "FATAL: unknown fluvial_surface_algorithm_choice '"
                  << method << "'. Valid: ladder, clubb2017, "
                     "cross_section_levels.\n";
        return 1;
    }
    const bool clubb_preset = (method == "clubb2017");
    const bool xs_preset    = (method == "cross_section_levels");
    struct stage_value { std::string key, value, source; };
    std::vector<stage_value> stages;
    auto resolve_stage = [&](const std::string& key,
                             const std::string& ladder_value,
                             const std::string& clubb_value,
                             const std::string& xs_value) -> std::string
    {
        const std::string& v = sf[key];
        if (v != "preset")
        {
            stages.push_back({key, v, "parameter file"});
            return v;
        }
        const std::string& pv = clubb_preset ? clubb_value
                              : xs_preset    ? xs_value
                                             : ladder_value;
        stages.push_back({key, pv, "preset " + method});
        return pv;
    };

    // Admission. A legacy fallback_admission_algorithm_choice only ever chose
    // the built-in admission when there was no flat mask, so it is honoured
    // only then, and only under the ladder preset it predates.
    std::string admission_choice;
    {
        const bool have_mask = !sf["flat_mask_fname"].empty();
        const bool named = sf["admission_algorithm_choice"] != "preset";
        const std::string& legacy =
            sf["legacy_fallback_admission_algorithm_choice"];
        if (!named && !clubb_preset && !xs_preset && !legacy.empty() &&
            !have_mask)
        {
            admission_choice = (legacy == "slope") ? "gradient" : legacy;
            stages.push_back({"admission_algorithm_choice", admission_choice,
                              "fallback_admission_algorithm_choice"});
            pp.record_resolved_parameter("admission_algorithm_choice",
                                         admission_choice);
        }
        else
        {
            admission_choice = resolve_stage("admission_algorithm_choice",
                                             "flat_mask", "floodplain_thresholds",
                                             "flat_mask");
        }
        if (admission_choice != "flat_mask" && admission_choice != "relief" &&
            admission_choice != "gradient" &&
            admission_choice != "floodplain_thresholds")
        {
            std::cerr << "FATAL: unknown admission_algorithm_choice '"
                      << admission_choice << "'. Valid: preset, flat_mask, "
                         "relief, gradient, floodplain_thresholds.\n";
            return 1;
        }
        if (have_mask && admission_choice != "flat_mask")
        {
            // Named in the file: a contradiction, refuse it. From a preset
            // (clubb2017 in a list run beside ladder): that preset's admission
            // wins for this method and the mask is for the others.
            if (named)
            {
                std::cerr << "FATAL: admission_algorithm_choice is "
                          << admission_choice << " but flat_mask_fname is set. "
                             "Use flat_mask to read the mask, or clear "
                             "flat_mask_fname.\n";
                return 1;
            }
            std::cout << "NOTE: the " << method << " preset admits by "
                      << admission_choice << "; flat_mask_fname is not used "
                         "by this method.\n";
        }
        if (admission_choice == "flat_mask" && !have_mask)
            admission_choice = "relief";   // no mask: the built-in fallback
    }

    const std::string level_choice =
        resolve_stage("level_algorithm_choice", "ladder", "fixed_thresholds",
                      "supplied_levels");
    const std::string floodplain_threshold_choice =
        resolve_stage("floodplain_threshold_algorithm_choice",
                      "fixed", "qq_gaussian", "fixed");
    const std::string channel_definition =
        resolve_stage("channel_definition_algorithm_choice",
                      "contributing_pixels", "stream_order",
                      "contributing_pixels");
    const std::string terrace_field_choice =
        resolve_stage("terrace_reference_field_algorithm_choice",
                      "same", "swath_baseline", "same");
    const std::string segmentation_choice =
        resolve_stage("segmentation_algorithm_choice",
                      "within_bands", "whole_mask", "within_bands");

    auto check_values = [&](const std::string& key, const std::string& v,
                            std::initializer_list<const char*> valid) -> bool
    {
        for (const char* ok : valid) if (v == ok) return true;
        std::cerr << "FATAL: unknown " << key << " '" << v << "'. Valid: preset";
        for (const char* ok : valid) std::cerr << ", " << ok;
        std::cerr << ".\n";
        return false;
    };
    if (!check_values("level_algorithm_choice", level_choice,
                      {"ladder", "fixed_thresholds", "supplied_levels"}) ||
        !check_values("floodplain_threshold_algorithm_choice",
                      floodplain_threshold_choice,
                      {"fixed", "qq_gaussian", "cdf_divergence"}) ||
        !check_values("channel_definition_algorithm_choice", channel_definition,
                      {"contributing_pixels", "stream_order"}) ||
        !check_values("terrace_reference_field_algorithm_choice",
                      terrace_field_choice, {"same", "swath_baseline"}) ||
        !check_values("segmentation_algorithm_choice", segmentation_choice,
                      {"within_bands", "whole_mask"}))
        return 1;
    const bool segment_whole_mask = (segmentation_choice == "whole_mask");

    // fill_floodplain_holes: false by default, true under clubb2017 unless
    // the file names it (a bool cannot carry the "preset" sentinel).
    bool fill_floodplain_holes = bf["fill_floodplain_holes"];
    if (pp.parameter_supplied("fill_floodplain_holes"))
    {
        stages.push_back({"fill_floodplain_holes",
                          fill_floodplain_holes ? "true" : "false",
                          "parameter file"});
    }
    else if (clubb_preset)
    {
        fill_floodplain_holes = true;
        stages.push_back({"fill_floodplain_holes", "true", "preset " + method});
        pp.record_resolved_parameter("fill_floodplain_holes", "true");
    }

    const bool level_fixed     = (level_choice == "fixed_thresholds");
    const bool level_supplied  = (level_choice == "supplied_levels");

    // Supplied levels -> per-channel band breaks, read before any work so a
    // missing or malformed file fails fast.
    std::map<std::int32_t, std::vector<double>> supplied_breaks;
    if (level_supplied)
    {
        if (sf["supplied_levels_fname"].empty())
        {
            std::cerr << "FATAL: level_algorithm_choice supplied_levels needs "
                         "supplied_levels_fname.\n";
            return 1;
        }
        try
        {
            const std::string lp =
                pp.add_read_path_if_needed(sf["supplied_levels_fname"]);
            const std::map<std::int32_t, std::vector<double>> levels =
                lsdtt::read_cross_section_levels_csv(lp);
            supplied_breaks = lsdtt::breaks_from_levels(levels);
            std::cout << "Supplied levels: " << lp << "\n";
            for (const auto& [sid, lv] : levels)
            {
                std::cout << "  channel " << sid << ":";
                for (double v : lv) std::cout << " " << v;
                std::cout << " m\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: supplied levels: " << e.what() << "\n";
            return 1;
        }
    }
    const bool admit_by_bounds = (admission_choice == "floodplain_thresholds");
    // The bounds are chosen when either stage reads them.
    const bool use_floodplain_bounds = level_fixed || admit_by_bounds;

    std::cout << "--- Method: " << method << " ---\n";
    for (const stage_value& sv : stages)
        std::cout << "  " << std::left << std::setw(42) << sv.key << std::right
                  << sv.value << "   (" << sv.source << ")\n";
    std::cout << "\n";

    // Parameters the resolved stages never read: say so rather than let a
    // setting silently do nothing.
    {
        std::vector<std::pair<std::string, bool>> unused = {
            {"min_channel_stream_order", channel_definition != "stream_order"},
            {"swath_half_width_metres",
             !(level_fixed && terrace_field_choice == "swath_baseline")},
            {"min_terrace_height_above_floodplain_metres", !level_fixed},
            {"max_terrace_height_metres", !level_fixed},
            {"fill_floodplain_holes", !level_fixed},
            {"supplied_levels_fname", !level_supplied},
            {"training_surfaces_layer_name", sf["training_surfaces_fname"].empty()},
            {"training_surfaces_class_field_name", sf["training_surfaces_fname"].empty()},
            {"training_spread_factor", sf["training_surfaces_fname"].empty()},
            {"min_training_cells_per_class", sf["training_surfaces_fname"].empty()},
            {"run_terrace_classification", !level_fixed},
            {"floodplain_second_pass_algorithm_choice", !level_fixed},
            {"write_floodplain_polygons", !level_fixed},
            {"riser_mask_fname", !segment_whole_mask},
            {"max_floodplain_gradient", !admit_by_bounds},
            {"floodplain_gradient_window_radius_metres", !admit_by_bounds},
            {"floodplain_bound_algorithm_choice", !admit_by_bounds},
            {"low_relief_max_gradient", !admit_by_bounds},
            {"low_relief_min_area_fraction", !admit_by_bounds},
            {"floodplain_qq_lower_percentile",
             !(use_floodplain_bounds && floodplain_threshold_choice == "qq_gaussian")},
            {"floodplain_qq_upper_percentile",
             !(use_floodplain_bounds && floodplain_threshold_choice == "qq_gaussian")},
            {"floodplain_qq_threshold_fraction",
             !(use_floodplain_bounds && floodplain_threshold_choice == "qq_gaussian")},
        };
        for (const auto& [key, is_unused] : unused)
            if (is_unused && pp.parameter_supplied(key))
                std::cerr << "WARNING: " << key << " is set but the resolved "
                             "method stages do not use it.\n";
    }

    // vertical_accuracy_break_gradient, or a legacy angle converted to one.
    // Held as float, the parameter's type, so the value recorded for
    // _ingested_params.param reads back exactly.
    float ev_break_gradient = ff["vertical_accuracy_break_gradient"];
    if (!sf["legacy_vertical_accuracy_break_degrees"].empty() &&
        !pp.parameter_supplied("vertical_accuracy_break_gradient"))
    {
        try
        {
            ev_break_gradient = static_cast<float>(
                lsdtt::gradient_from_legacy_degrees(
                    sf["legacy_vertical_accuracy_break_degrees"],
                    "vertical_accuracy_break_degrees"));
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
        std::ostringstream os;
        os.precision(9);
        os << ev_break_gradient;
        pp.record_resolved_parameter("vertical_accuracy_break_gradient",
                                     os.str());
        std::cout << "vertical_accuracy_break_degrees "
                  << sf["legacy_vertical_accuracy_break_degrees"]
                  << " -> vertical_accuracy_break_gradient " << os.str()
                  << "\n";
    }

    lsdtt::reference_field_choice field_choice =
        lsdtt::reference_field_choice::hand;
    {
        const std::string& rf = sf["reference_field_algorithm_choice"];
        if (rf == "hand")
            field_choice = lsdtt::reference_field_choice::hand;
        else if (rf == "fitted_surface")
            field_choice = lsdtt::reference_field_choice::fitted_surface;
        else
        {
            std::cerr << "FATAL: unknown reference_field_algorithm_choice '"
                      << rf << "'. Valid: hand, fitted_surface.\n";
            return 1;
        }
    }
    const bool use_hand =
        field_choice == lsdtt::reference_field_choice::hand;
    if (use_floodplain_bounds && !use_hand)
    {
        std::cerr << "FATAL: " << (level_fixed ? "level_algorithm_choice "
                                                 "fixed_thresholds"
                                               : "admission_algorithm_choice "
                                                 "floodplain_thresholds")
                  << " needs reference_field_algorithm_choice hand: "
                     "floodplain bounds are heights above a channel.\n";
        return 1;
    }

    lsdtt::fitted_surface_choice fit_choice =
        lsdtt::fitted_surface_choice::cone;
    {
        const std::string& fs = sf["fitted_surface_algorithm_choice"];
        if (fs == "cone")              fit_choice = lsdtt::fitted_surface_choice::cone;
        else if (fs == "plane_ransac") fit_choice = lsdtt::fitted_surface_choice::plane_ransac;
        else if (fs == "poly2")        fit_choice = lsdtt::fitted_surface_choice::poly2;
        else
        {
            std::cerr << "FATAL: unknown fitted_surface_algorithm_choice '"
                      << fs << "'. Valid: cone, plane_ransac, poly2.\n";
            return 1;
        }
    }

    // The antimode finder bins on log10, so it cannot band a signed field, and
    // every fitted-surface residual is signed. Switch silently only when the
    // user did not name an algorithm; reject an explicit request, because
    // quietly running something other than what was asked for is worse than
    // saying no.
    lsdtt::ladder_break_choice break_choice =
        lsdtt::ladder_break_choice::antimode;
    {
        const bool named = pp.parameter_supplied("ladder_break_algorithm_choice");
        const std::string& ba = sf["ladder_break_algorithm_choice"];
        if (ba == "antimode")   break_choice = lsdtt::ladder_break_choice::antimode;
        else if (ba == "jenks") break_choice = lsdtt::ladder_break_choice::jenks;
        else if (ba == "gap")   break_choice = lsdtt::ladder_break_choice::gap;
        else
        {
            std::cerr << "FATAL: unknown ladder_break_algorithm_choice '" << ba
                      << "'. Valid: antimode, jenks, gap.\n";
            return 1;
        }

        if (!use_hand && break_choice == lsdtt::ladder_break_choice::antimode)
        {
            if (named)
            {
                std::cerr
                    << "FATAL: ladder_break_algorithm_choice antimode cannot "
                       "be used with reference_field_algorithm_choice "
                       "fitted_surface.\n"
                       "  The antimode finder bins on log10, so every negative "
                       "residual is floored into the lowest bin and the whole "
                       "below-surface population - the very thing a fitted "
                       "surface is used to find - collapses into one spike.\n"
                       "  Use ladder_break_algorithm_choice: jenks.\n";
                return 1;
            }
            std::cout << "NOTE: reference_field_algorithm_choice is "
                         "fitted_surface, so the break finder is jenks.\n"
                         "  antimode bins on log10 and cannot band a field "
                         "with negative values.\n";
            break_choice = lsdtt::ladder_break_choice::jenks;
        }
    }

    const std::string& ev_mode_name = sf["vertical_accuracy_algorithm_choice"];
    const bool ev_slope_dependent = (ev_mode_name == "slope_dependent");
    if (ev_mode_name != "fixed" && !ev_slope_dependent)
    {
        std::cerr << "FATAL: unknown vertical_accuracy_algorithm_choice '" << ev_mode_name
                  << "'. Valid: fixed, slope_dependent.\n";
        return 1;
    }
    if (ev_slope_dependent)
    {
        const double eg = static_cast<double>(ff["vertical_accuracy_gentle_metres"]);
        const double es = static_cast<double>(ff["vertical_accuracy_steep_metres"]);
        if (!(eg > 0.0) || !(es > 0.0))
        {
            std::cerr << "FATAL: vertical_accuracy_algorithm_choice slope_dependent needs "
                         "both vertical_accuracy_gentle_metres and "
                         "vertical_accuracy_steep_metres > 0.\n";
            return 1;
        }
        if (es < eg)
        {
            std::cerr << "FATAL: vertical_accuracy_steep_metres (" << es
                      << ") is below vertical_accuracy_gentle_metres (" << eg
                      << "). A DEM does not get MORE accurate on a slope; if "
                         "that is really what you mean, use fixed.\n";
            return 1;
        }
    }

    if (!level_fixed && !ev_slope_dependent && ff["vertical_accuracy_metres"] <= 0.0f)
    {
        std::cerr << "WARNING: vertical_accuracy_metres is 0, so the "
                     "resolvability floor is disabled.\n"
                     "  Every window will report `resolved` and the ladder "
                     "carries no instrument guarantee.\n"
                     "  Set it to your DEM's vertical accuracy to have levels "
                     "checked against what the DEM can measure.\n"
                     "  (A DEM whose vertical units are feet on a metre grid "
                     "has an effective accuracy 3.28x larger than it looks.)\n";
    }

    // How flats are made drainable (Barnes et al. 2014b options). Checked
    // here, before any raster work. Only the routing fill uses it: with
    // fill_before_routing: false the DEM is routed as it is.
    lsdtt::flat_resolution_choice flat_choice =
        lsdtt::flat_resolution_choice::inline_gradient;
    try
    {
        flat_choice = lsdtt::parse_flat_resolution_choice(
            sf["fill_flat_resolution_choice"]);
        if (bf["fill_before_routing"])
        {
            lsdtt::check_fill_flat_resolution_min_slope(
                flat_choice, ff["fill_min_slope_for_routing"],
                "fill_min_slope_for_routing");
        }
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    if (!bf["fill_before_routing"] &&
        flat_choice != lsdtt::flat_resolution_choice::inline_gradient)
    {
        std::cout << "NOTE: fill_flat_resolution_choice ("
                  << lsdtt::flat_resolution_choice_name(flat_choice)
                  << ") is ignored because fill_before_routing is false.\n";
    }

    // -------------------------------------------------------------------------
    // Load DEM
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_raster_cpu  dem;
    {
        const std::string in_path = pp.get_full_read_path();
        std::cout << "Loading DEM: " << in_path << "\n";
        try { loader.open_raster(in_path); }
        catch (const std::exception& e)
        { std::cerr << "FATAL: loading DEM: " << e.what() << "\n"; return 1; }
        dem = lsdtt::lsdtt_raster_cpu(loader);
        try { dem.get_meta().check_projected_metres(); }
        catch (const std::exception& e)
        { std::cerr << "FATAL: projection: " << e.what() << "\n"; return 1; }
        std::cout << "  rows=" << dem.get_meta().n_rows
                  << "  cols=" << dem.get_meta().n_cols
                  << "  cell=" << dem.get_meta().cell_size << " m\n\n";
    }

    // Patch aspects: grid-relative unless bearings_relative_to_true_north.
    std::unique_ptr<lsdtt::lsdtt_grid_convergence> true_north;
    try
    {
        true_north = lsdtt::start_true_north_bearings(
            dem.get_meta(), bf["bearings_relative_to_true_north"], std::cout);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: bearings_relative_to_true_north: " << e.what() << "\n";
        return 1;
    }
    // Interior nodata: kept as sinks or filled (interior_nodata_method).
    try
    {
        std::cout << lsdtt::apply_interior_nodata_method(
                         dem, sf.at("interior_nodata_method"),
                         sf.at("nodata_fill_method_choice"),
                         ip.at("nodata_fill_max_iter"), ip.at("n_threads"))
                  << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: interior nodata: " << e.what() << "\n";
        return 1;
    }
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n_pixels =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    auto load_companion_raster = [&](const std::string& fname,
                                     const char* what,
                                     lsdtt::lsdtt_raster_cpu& out) -> bool
    {
        const std::string path = pp.add_read_path_if_needed(fname);
        std::cout << "Loading " << what << ": " << path << "\n";
        lsdtt::lsdtt_data_loader companion;
        companion.open_raster(path);
        out = lsdtt::lsdtt_raster_cpu(companion);
        const lsdtt::lsdtt_raster_meta& m = out.get_meta();
        if (m.n_rows != meta.n_rows || m.n_cols != meta.n_cols)
        {
            std::cerr << "FATAL: " << what << " is " << m.n_rows << "x"
                      << m.n_cols << " but the DEM is " << meta.n_rows << "x"
                      << meta.n_cols << ".\n";
            return false;
        }
        return true;
    };

    // Risers for segmentation_algorithm_choice whole_mask (1/0 on the DEM grid).
    lsdtt::lsdtt_raster_cpu riser_mask;
    const lsdtt::lsdtt_raster_cpu* riser_barrier = nullptr;
    if (segment_whole_mask && !sf["riser_mask_fname"].empty())
    {
        lsdtt::lsdtt_raster_cpu in;
        try
        {
            if (!load_companion_raster(sf["riser_mask_fname"], "riser mask", in))
                return 1;
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: loading riser mask: " << e.what() << "\n";
            return 1;
        }
        std::vector<float> rd(n_pixels, 0.0f);
        std::int64_t n_riser = 0;
        for (std::size_t i = 0; i < n_pixels; ++i)
            if (in.get_mask()[i] && in.get_data()[i] > 0.5f) { rd[i] = 1.0f; ++n_riser; }
        riser_mask = lsdtt::lsdtt_raster_cpu(std::move(rd),
                                             std::vector<uint8_t>(n_pixels, 1), meta);
        riser_barrier = &riser_mask;
        std::cout << "  " << n_riser << " riser cell(s)\n";
    }
    else if (segment_whole_mask)
    {
        std::cout << "segmentation whole_mask with no riser_mask_fname: the "
                     "admitted ground is labelled as it stands.\n";
    }

    // -------------------------------------------------------------------------
    // Admission: which ground may take part
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu admitted;
    {
        std::cout << "--- Admission ---\n";
        if (admit_by_bounds)
        {
            // floodplain_thresholds: every valid cell for now. The gradient
            // bound is read off cells with a defined HAND, so it is applied
            // after the reference field (below).
            std::cout << "Admission by floodplain bounds: the gradient bound "
                         "is applied once HAND is known.\n";
            admitted = lsdtt::lsdtt_raster_cpu(
                std::vector<float>(n_pixels, 1.0f),
                std::vector<uint8_t>(dem.get_mask()), meta);
        }
        else if (admission_choice == "flat_mask")
        {
            lsdtt::lsdtt_raster_cpu in_mask;
            try
            {
                if (!load_companion_raster(sf["flat_mask_fname"],
                                           "flat mask", in_mask))
                    return 1;
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: loading flat mask: " << e.what() << "\n";
                return 1;
            }
            const std::vector<float>&   mv = in_mask.get_data();
            const std::vector<uint8_t>& mm = in_mask.get_mask();
            std::vector<float>   ad(n_pixels, 0.0f);
            std::vector<uint8_t> am(n_pixels, 1);
            for (std::size_t i = 0; i < n_pixels; ++i)
                if (mm[i] && mv[i] > 0.5f) ad[i] = 1.0f;
            admitted = lsdtt::lsdtt_raster_cpu(
                std::move(ad), std::move(am), meta);
        }
        else
        {
            // No flat mask: a coarse built-in admission so the tool can be run
            // on its own for a first look. This is deliberately not as good as
            // lsdtt-flat-surface-finder's texture metrics.
            const std::string& how = admission_choice;
            std::cout << "No flat_mask_fname; using the built-in " << how
                      << " admission.\n"
                      << "  For a publishable result run "
                         "lsdtt-flat-surface-finder first and pass its "
                         "_combined_flat_mask.tif as flat_mask_fname.\n";
            std::vector<float>   ad(n_pixels, 0.0f);
            std::vector<uint8_t> am(n_pixels, 1);
            const std::vector<uint8_t>& dm = dem.get_mask();

            if (how == "relief")
            {
                // Radius in cells: relief_window_radius_metres, or a legacy
                // fallback_relief_window_radius_pixels, or the default 5.
                int radius_px = 0;
                const std::string& legacy_px =
                    sf["legacy_relief_window_radius_pixels"];
                try
                {
                    if (!legacy_px.empty() &&
                        !pp.parameter_supplied("relief_window_radius_metres"))
                    {
                        std::size_t used = 0;
                        radius_px = std::stoi(legacy_px, &used);
                        if (used != legacy_px.size() || radius_px < 1)
                            throw std::invalid_argument(
                                "fallback_relief_window_radius_pixels must be "
                                "a whole number >= 1; got '" + legacy_px + "'.");
                        std::ostringstream os;
                        os.precision(9);
                        os << static_cast<float>(radius_px * meta.cell_size);
                        pp.record_resolved_parameter(
                            "relief_window_radius_metres", os.str());
                    }
                    else
                    {
                        radius_px = lsdtt::relief_window_radius_pixels(
                            ff["relief_window_radius_metres"],
                            meta.cell_size, 5);
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: " << e.what() << "\n";
                    return 1;
                }
                const double radius =
                    static_cast<double>(radius_px) * meta.cell_size;
                std::cout << "  Relief window radius: " << radius_px
                          << " cell(s) = " << radius << " m\n";
                const lsdtt::lsdtt_raster_cpu relief =
                    lsdtt::moving_window_local_relief(dem, radius);
                const std::vector<float>& rv = relief.get_data();
                const float thr = ff["max_relief_metres"];
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (dm[i] && std::isfinite(rv[i]) && rv[i] < thr)
                        ad[i] = 1.0f;
            }
            else if (how == "gradient")
            {
                const lsdtt::poly_fit_coeffs pc =
                    lsdtt::direct_polyfit(dem, 2.0 * meta.cell_size,
                                          ip["n_threads"]);
                lsdtt::polyfit_metric_request mreq;
                mreq.metrics = {lsdtt::polyfit_metric::slope};
                const lsdtt::polyfit_metric_result mres =
                    lsdtt::compute_polyfit_metrics(pc, mreq);
                const std::vector<float>& sv =
                    mres.at(lsdtt::polyfit_metric::slope).get_data();
                const float thr = ff["max_admission_gradient"];
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (dm[i] && std::isfinite(sv[i]) && sv[i] < thr)
                        ad[i] = 1.0f;
            }
            else
            {
                std::cerr << "FATAL: unknown admission_algorithm_choice '"
                          << how << "'. Valid: flat_mask, relief, gradient.\n";
                return 1;
            }
            admitted = lsdtt::lsdtt_raster_cpu(
                std::move(ad), std::move(am), meta);
        }

        // fill_admitted_holes: roughness (noise, crops, buildings) leaves
        // enclosed gaps in a relief / gradient / flat-mask admission that
        // become holes in every surface built on it (Mayoco's flat mask:
        // 988 holes). The floodplain_thresholds admission is not filled here:
        // its floodplain has its own fill (fill_floodplain_holes).
        if (!admit_by_bounds && bf["fill_admitted_holes"])
        {
            std::int64_t before = 0, after = 0;
            for (const float v : admitted.get_data()) if (v > 0.5f) ++before;
            admitted = lsdtt::binary_fill_holes(admitted);
            for (const float v : admitted.get_data()) if (v > 0.5f) ++after;
            std::cout << "  fill_admitted_holes: " << (after - before)
                      << " enclosed cell(s) admitted\n";
        }
        std::int64_t n_admitted = 0;
        for (const float v : admitted.get_data()) if (v > 0.5f) ++n_admitted;
        std::cout << "  " << n_admitted << " admitted cell(s)\n\n";
        if (n_admitted == 0)
        {
            std::cerr << "FATAL: no cells admitted; nothing to band.\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Flow routing and channel scope
    //
    // Needed for the hand field, for any channel-linked output, and for the
    // channel network. A fitted_surface run with no channel outputs needs none
    // of it and skips the whole block.
    // -------------------------------------------------------------------------
    const bool need_flow =
        use_hand ||
        bf["write_channel_scope_raster"] || bf["write_channel_scope_nodes_csv"] ||
        bf["write_point_snap_csv"] || bf["write_channel_network"];

    lsdtt::lsdtt_raster_cpu routing_dem = dem;
    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    lsdtt::channel_scope_result scope;
    bool scope_resolved = false;

    if (need_flow)
    {
        std::cout << "--- D8 flow routing ---\n";
        std::cout << "Fill before routing  : "
                  << (bf["fill_before_routing"] ? "yes" : "no") << "\n";
        try
        {
            if (bf["fill_before_routing"])
            {
                std::cout << "Filling DEM for flow routing...\n";
                if (flat_choice ==
                    lsdtt::flat_resolution_choice::inline_gradient)
                {
                    routing_dem = lsdtt::fill_pits(
                        dem, ff["fill_min_slope_for_routing"], ip["n_threads"]);
                }
                else
                {
                    routing_dem = lsdtt::fill_pits_reporting_flat_resolution(
                        dem, ff["fill_min_slope_for_routing"], flat_choice,
                        ip["n_threads"]);
                }
            }
            else
            {
                std::cout << "Routing on the unfilled input DEM.\n";
            }
            lsdtt::lsdtt_flowinfo_options fo;
            fo.boundary_conditions       = sf["boundary_conditions"];
            fo.interior_nodata_is_outlet = true;
            fo.require_acyclic_graph     = true;
            fo.build_topological_order   = true;
            fo.build_s_vector            = true;
            fo.n_threads                 = ip["n_threads"];
            flow = lsdtt::make_flow_base(
                routing_dem.get_data_array(), routing_dem.get_mask_array(),
                routing_dem.get_meta(), fo);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
            return 1;
        }

        std::cout << "\n--- Channel scope resolution ---\n";
        try
        {
            // The five-key convention lives in lsdtt_channel_scope_cpu; this CLI
            // spells its paths as <role>_points_fname and resolves them against
            // read_path, which the accessor below does.
            auto role_flats = [&](const std::string& role)
            {
                return lsdtt::load_role_point_flats_by_key(
                    role, "_fname",
                    [&](const std::string& key) {
                        const std::string v = sf[key];
                        return (key.size() >= 6 &&
                                key.compare(key.size() - 6, 6, "_fname") == 0 &&
                                !v.empty())
                            ? pp.add_read_path_if_needed(v)
                            : v;
                    },
                    routing_dem.get_meta());
            };

            lsdtt::channel_scope_request sreq;
            sreq.threshold_contributing_pixels =
                ip["threshold_contributing_pixels"];
            sreq.source_flats = role_flats("source");
            sreq.outlet_flats = role_flats("outlet");
            sreq.snap_sources = bf["snap_sources_to_channel_threshold"];
            sreq.source_snap_max_radius_pixels =
                ip["source_snap_max_search_radius_pixels"];
            sreq.outlet_snap_max_distance_metres =
                static_cast<double>(ff["outlet_snap_max_distance_metres"]);
            sreq.auto_detect_roles =
                bf["auto_detect_outlet_points"] &&
                sf["outlet_points_fname"].empty() &&
                !sreq.source_flats.empty();
            if (channel_definition == "stream_order")
            {
                sreq.channel_definition =
                    lsdtt::channel_definition_choice::stream_order;
                sreq.min_stream_order = ip["min_channel_stream_order"];
            }

            scope = lsdtt::resolve_channel_scope(*flow, sreq);
            scope_resolved = true;
            std::cout << "  " << scope.resolved_sources.size()
                      << " channel source(s) in scope\n";
            if (channel_definition == "stream_order")
            {
                std::int64_t n_chan = 0;
                std::int32_t max_order = 0;
                for (std::size_t i = 0; i < scope.channel_mask.size(); ++i)
                {
                    if (scope.channel_mask[i]) ++n_chan;
                    max_order = std::max(max_order, scope.stream_order[i]);
                }
                std::cout << "  channel = Strahler order >= "
                          << ip["min_channel_stream_order"] << ": " << n_chan
                          << " channel pixel(s); highest order on the network "
                          << max_order << "\n";
                if (n_chan == 0)
                {
                    std::cerr << "FATAL: no channel pixel of order >= "
                              << ip["min_channel_stream_order"]
                              << " (the network reaches order " << max_order
                              << "). Lower min_channel_stream_order or "
                                 "threshold_contributing_pixels.\n";
                    return 1;
                }
            }

            lsdtt::channel_scope_diagnostic_request diag;
            diag.write_scope_raster    = bf["write_channel_scope_raster"];
            diag.write_scope_nodes_csv = bf["write_channel_scope_nodes_csv"];
            diag.write_point_snap_csv  = bf["write_point_snap_csv"];
            diag.raster_format         = out_format;
            lsdtt::write_channel_scope_diagnostics(
                scope, *flow, dem, loader, stem, diag);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: channel scope resolution failed: "
                      << e.what() << "\n";
            return 1;
        }
        std::cout << "\n";
    }

    // -------------------------------------------------------------------------
    // The reference field
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu reference_field;
    lsdtt::lsdtt_raster_cpu flow_distance;
    std::vector<std::int64_t> donor_flat;
    bool have_donors = false;

    if (use_hand)
    {
        std::cout << "--- Reference field: height above channel ---\n";
        if (!sf["hand_fname"].empty())
        {
            // The supplied field is used exactly as given. Donors are not
            // recoverable from a raster, so a supplied HAND has no reach
            // coordinate and the ladder runs globally: an intentional and
            // stated cost of the composability lever.
            try
            {
                if (!load_companion_raster(sf["hand_fname"], "HAND",
                                           reference_field))
                    return 1;
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: loading HAND: " << e.what() << "\n";
                return 1;
            }
            std::cout << "  Using the supplied HAND raster. Reach-relative "
                         "windowing is unavailable (a raster carries no donor "
                         "nodes), so the ladder is global.\n";
        }
        else
        {
            if (!scope_resolved)
            {
                std::cerr << "FATAL: no channel scope, so HAND cannot be "
                             "computed.\n";
                return 1;
            }
            try
            {
                std::cout << "Computing height above the in-scope channel "
                             "network...\n";
                std::vector<uint8_t> cmask(scope.channel_mask.begin(),
                                           scope.channel_mask.end());
                lsdtt::height_above_channel_donor_result hd =
                    flow->get_height_above_channel_mask_donor_raster(
                        dem.get_data(), dem.get_mask(), cmask);
                reference_field = lsdtt::lsdtt_raster_cpu(std::move(hd.height));
                donor_flat      = std::move(hd.donor_flat);
                flow_distance   = lsdtt::lsdtt_raster_cpu(
                    flow->get_distance_from_outlet_raster());
                have_donors     = true;
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: HAND computation failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (bf["write_hand_raster"])
        {
            const std::string p = stem + "_hand.tif";
            std::cout << "Writing HAND: " << p << "\n";
            loader.write_full(reference_field.get_data(),
                              reference_field.get_meta(), p, out_format);
        }
    }
    else
    {
        std::cout << "--- Reference field: residual to a fitted surface ---\n";
        lsdtt::lsdtt_raster_cpu fit_mask = admitted;
        if (!sf["fitted_surface_mask_fname"].empty())
        {
            try
            {
                if (!load_companion_raster(sf["fitted_surface_mask_fname"],
                                           "fitted surface mask", fit_mask))
                    return 1;
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: loading fitted surface mask: "
                          << e.what() << "\n";
                return 1;
            }
        }

        lsdtt::fitted_surface_request freq;
        freq.algorithm               = fit_choice;
        freq.ransac_tolerance_metres =
            static_cast<double>(ff["fitted_surface_ransac_tolerance_metres"]);
        freq.ransac_iterations       = ip["fitted_surface_ransac_iterations"];
        freq.seed                    = ip["fitted_surface_seed"];
        freq.n_threads               = ip["n_threads"];
        freq.apex_x = static_cast<double>(ff["fitted_surface_apex_x"]);
        freq.apex_y = static_cast<double>(ff["fitted_surface_apex_y"]);
        freq.apex_is_set = (freq.apex_x != 0.0 || freq.apex_y != 0.0);

        lsdtt::fitted_surface_field fs;
        try
        {
            fs = lsdtt::build_fitted_surface_reference_field(
                dem, fit_mask, freq);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: surface fitting failed: " << e.what() << "\n";
            return 1;
        }
        if (!fs.valid)
        {
            std::cerr << "FATAL: the surface fit failed (too few usable cells, "
                         "or a degenerate geometry).\n";
            return 1;
        }

        std::cout << "  " << sf["fitted_surface_algorithm_choice"]
                  << " fit: rmse " << fs.rmse << " m over "
                  << fs.n_pixels_used << " cell(s)\n";
        if (fit_choice == lsdtt::fitted_surface_choice::cone)
            std::cout << "  apex (" << fs.cone_coeffs.apex_x << ", "
                      << fs.cone_coeffs.apex_y << "), radial gradient "
                      << fs.cone_coeffs.gradient << " m/m\n";
        reference_field = std::move(fs.field);
    }

    if (bf["write_reference_field_raster"])
    {
        const std::string p = stem + "_reference_field.tif";
        std::cout << "Writing reference field: " << p << "\n";
        loader.write_full(reference_field.get_data(),
                          reference_field.get_meta(), p, out_format);
    }
    std::cout << "\n";

    // -------------------------------------------------------------------------
    // Floodplain bounds (admission floodplain_thresholds, level
    // fixed_thresholds)
    //
    // The bounds are read off the cells that could be floodplain: admitted,
    // with a defined HAND (a pixel that reaches no in-scope channel has no
    // height above one), and for the gradient a defined gradient.
    // -------------------------------------------------------------------------
    // One row per data-chosen floodplain bound, written as
    // <stem>_floodplain_thresholds.csv: the numbers a plot of the QQ tables
    // needs, so nothing downstream has to recompute them.
    std::vector<std::vector<std::string>> floodplain_threshold_rows;
    auto num9 = [](double v)
    {
        std::ostringstream os;
        os.precision(9);
        os << v;
        return os.str();
    };
    float floodplain_hand_bound     = ff["max_floodplain_hand_metres"];
    float floodplain_gradient_bound = ff["max_floodplain_gradient"];
    lsdtt::lsdtt_raster_cpu floodplain_gradient;
    // floodplain_bound_algorithm_choice: false = HAND alone (hand_only).
    bool use_floodplain_gradient_bound = true;
    double low_relief_fraction = -1.0;   // auto only; -1 = not measured
    if (use_floodplain_bounds)
    {
        std::cout << "--- Floodplain bounds (" << floodplain_threshold_choice
                  << ") ---\n";
        if (admit_by_bounds)
        {
            const double cell = dem.get_meta().cell_size;
            double radius = static_cast<double>(
                ff["floodplain_gradient_window_radius_metres"]);
            if (!(radius > 0.0)) radius = 3.0 * cell;
            std::cout << "Floodplain gradient: polyfit window radius " << radius
                      << " m\n";
            try
            {
                // direct_polyfit takes the full physical window diameter.
                const lsdtt::poly_fit_coeffs fit =
                    lsdtt::direct_polyfit(dem, 2.0 * radius, ip["n_threads"]);
                lsdtt::polyfit_metric_request mreq;
                mreq.metrics   = {lsdtt::polyfit_metric::slope};
                mreq.n_threads = ip["n_threads"];
                lsdtt::polyfit_metric_result mres =
                    lsdtt::compute_polyfit_metrics(fit, mreq);
                floodplain_gradient = mres.at(lsdtt::polyfit_metric::slope);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: floodplain gradient failed: " << e.what()
                          << "\n";
                return 1;
            }
            // floodplain_bound_algorithm_choice: is the gradient bound used?
            {
                const std::string& bc = sf["floodplain_bound_algorithm_choice"];
                if (bc != "auto" && bc != "hand_and_gradient" && bc != "hand_only")
                {
                    std::cerr << "FATAL: unknown floodplain_bound_algorithm_choice '"
                              << bc << "'. Valid: auto, hand_and_gradient, "
                                 "hand_only.\n";
                    return 1;
                }
                if (bc == "auto")
                {
                    const float g_low = ff["low_relief_max_gradient"];
                    std::int64_t n_def = 0, n_low = 0;
                    const std::vector<float>&   gd = floodplain_gradient.get_data();
                    const std::vector<uint8_t>& gm = floodplain_gradient.get_mask();
                    for (std::size_t i = 0; i < n_pixels; ++i)
                    {
                        if (!gm[i] || !std::isfinite(gd[i]) ||
                            floodplain_gradient.get_meta().is_nodata(gd[i]))
                            continue;
                        ++n_def;
                        if (gd[i] < g_low) ++n_low;
                    }
                    const double frac = n_def > 0 ? double(n_low) / double(n_def) : 0.0;
                    low_relief_fraction = frac;
                    use_floodplain_gradient_bound =
                        !(frac > static_cast<double>(ff["low_relief_min_area_fraction"]));
                    std::cout << "  low relief test: " << frac
                              << " of the DEM's floodplain gradient is below "
                              << g_low << " (switch above "
                              << ff["low_relief_min_area_fraction"] << ") -> "
                              << (use_floodplain_gradient_bound ? "hand_and_gradient"
                                                                : "hand_only")
                              << "\n";
                    pp.record_resolved_parameter(
                        "floodplain_bound_algorithm_choice",
                        use_floodplain_gradient_bound ? "hand_and_gradient"
                                                      : "hand_only");
                }
                else
                {
                    use_floodplain_gradient_bound = (bc == "hand_and_gradient");
                }
            }
            if (bf["write_floodplain_gradient_raster"])
            {
                const std::string p = stem + "_floodplain_gradient.tif";
                std::cout << "Writing floodplain gradient: " << p << "\n";
                loader.write_full(floodplain_gradient.get_data(),
                                  floodplain_gradient.get_meta(), p, out_format);
            }
        }

        const std::vector<float>&   rf  = reference_field.get_data();
        const std::vector<uint8_t>& rfm = reference_field.get_mask();
        const std::vector<float>&   ad  = admitted.get_data();
        auto hand_defined = [&](std::size_t i)
        {
            return ad[i] > 0.5f && rfm[i] && std::isfinite(rf[i]);
        };
        auto grad_defined = [&](std::size_t i)
        {
            return floodplain_gradient.get_mask()[i] &&
                   std::isfinite(floodplain_gradient.get_data()[i]);
        };

        if (floodplain_threshold_choice != "fixed")
        {
            // stem suffix for the table a bound was read from
            const std::string table_suffix =
                floodplain_threshold_choice == "qq_gaussian" ? "_qq.csv"
                                                             : "_adaptive_cdf.csv";
            auto choose = [&](const std::vector<double>& values,
                              double fixed_value, const char* what,
                              const std::string& out_name,
                              float& bound) -> bool
            {
                lsdtt::threshold_selection_result sel;
                try
                {
                    sel = lsdtt::select_threshold(values,
                        lsdtt::threshold_request_from_parameters(
                            floodplain_threshold_choice,
                            "floodplain_threshold_algorithm_choice",
                            fixed_value, ff, "floodplain"));
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: floodplain " << what << " bound: "
                              << e.what() << "\n";
                    return false;
                }
                if (!sel.found)
                {
                    std::cerr << "FATAL: " << floodplain_threshold_choice
                              << " found no floodplain " << what << " bound in "
                              << values.size() << " value(s): the data never "
                                 "settle onto the reference line. Look at the "
                                 "table, or set the bound and use fixed.\n";
                    return false;
                }
                bound = static_cast<float>(sel.threshold);
                floodplain_threshold_rows.push_back({
                    "1", out_name, floodplain_threshold_choice,
                    num9(sel.threshold),
                    floodplain_threshold_choice == "qq_gaussian"
                        ? std::to_string(sel.qq_index) : "",
                    std::to_string(values.size()),
                    fname_stem + "_floodplain_" + out_name + table_suffix});
                if (bf["write_floodplain_threshold_csv"])
                {
                    const std::string p = stem + "_floodplain_" + out_name +
                                          table_suffix;
                    if (floodplain_threshold_choice == "qq_gaussian")
                        lsdtt::write_qq_csv(loader, p, sel.qq);
                    else
                        lsdtt::write_smoothed_cdf_csv(loader, p, sel.cdf,
                                                      sel.fits);
                    std::cout << "  Writing " << p << "\n";
                }
                if (floodplain_threshold_choice == "qq_gaussian")
                    std::cout << "  " << what << ": quantile " << sel.qq_index
                              << " of " << sel.qq.data_quantile.size()
                              << " (data range " << sel.qq_range << ")\n";
                else
                    std::cout << "  " << what << ": best fit "
                              << sel.best_family << ", percentile "
                              << sel.threshold_percentile << "\n";
                return true;
            };

            if (level_fixed)
            {
                std::vector<double> hv;
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (hand_defined(i)) hv.push_back(rf[i]);
                if (!choose(hv, floodplain_hand_bound, "HAND", "hand",
                            floodplain_hand_bound))
                    return 1;
                std::ostringstream os; os.precision(9); os << floodplain_hand_bound;
                pp.record_resolved_parameter("max_floodplain_hand_metres", os.str());
            }
            if (admit_by_bounds && use_floodplain_gradient_bound)
            {
                std::vector<double> gv;
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (hand_defined(i) && grad_defined(i))
                        gv.push_back(floodplain_gradient.get_data()[i]);
                if (!choose(gv, floodplain_gradient_bound, "gradient",
                            "gradient", floodplain_gradient_bound))
                    return 1;
                std::ostringstream os; os.precision(9); os << floodplain_gradient_bound;
                pp.record_resolved_parameter("max_floodplain_gradient", os.str());
            }
        }

        if (level_fixed)
            std::cout << "  floodplain HAND bound     : " << floodplain_hand_bound
                      << " m\n";
        if (admit_by_bounds)
        {
            if (use_floodplain_gradient_bound)
                std::cout << "  floodplain gradient bound : "
                          << floodplain_gradient_bound << " (m/m)\n";
            else
                std::cout << "  floodplain gradient bound : none (hand_only)\n";
            std::vector<float>& ad_w = admitted.get_data();
            std::int64_t n_in = 0;
            for (std::size_t i = 0; i < n_pixels; ++i)
            {
                const bool keep = hand_defined(i) &&
                    (!use_floodplain_gradient_bound ||
                     (grad_defined(i) &&
                      floodplain_gradient.get_data()[i] < floodplain_gradient_bound));
                ad_w[i] = keep ? 1.0f : 0.0f;
                if (keep) ++n_in;
            }
            std::cout << (use_floodplain_gradient_bound
                              ? "  admitted below the gradient bound: "
                              : "  admitted (HAND defined): ")
                      << n_in << " cell(s)\n";
        }
        std::cout << "\n";
    }

    // Negative HAND -> 0 for the antimode ladder (Simon, 2026-10-27). The
    // antimode finder bins log10 HAND and refuses a negative value, and every
    // real DEM tried had some: cells routed to a channel cell higher than
    // themselves (fill and flat-resolution artefacts). Below the channel is
    // height 0 above it; the clamp count is in the log and the method record.
    std::int64_t n_hand_clamped = 0;
    if (use_hand && !level_fixed &&
        break_choice == lsdtt::ladder_break_choice::antimode)
    {
        std::vector<float>&         rfd = reference_field.get_data();
        const std::vector<uint8_t>& rfm = reference_field.get_mask();
        for (std::size_t i = 0; i < rfd.size(); ++i)
            if (rfm[i] && std::isfinite(rfd[i]) &&
                !reference_field.get_meta().is_nodata(rfd[i]) && rfd[i] < 0.0f)
            {
                rfd[i] = 0.0f;
                ++n_hand_clamped;
            }
        if (n_hand_clamped > 0)
            std::cout << "  " << n_hand_clamped << " negative HAND cell(s) set "
                         "to 0 for the antimode ladder\n";
    }

    // ---- training polygons (plan 7i) ----------------------------------------
    const bool training_on = !sf["training_surfaces_fname"].empty();
    lsdtt::training_set training;
    std::vector<lsdtt::training_polygon>  training_polygon_stats;
    std::vector<std::vector<std::string>> training_fit_rows;
    const double training_spread_factor =
        static_cast<double>(ff["training_spread_factor"]);
    const int training_min_cells = ip["min_training_cells_per_class"];
    if (training_on)
    {
        std::cout << "--- Training polygons ---\n";
        if (!use_hand)
        {
            std::cerr << "FATAL: training_surfaces_fname sets HAND bounds and "
                         "breaks; reference_field_choice "
                      << sf["reference_field_choice"] << " has no HAND.\n";
            return 1;
        }
        const std::string tp =
            pp.add_read_path_if_needed(sf["training_surfaces_fname"]);
        try
        {
            training = lsdtt::read_training_surfaces(
                reference_field, tp, sf["training_surfaces_layer_name"],
                sf["training_surfaces_class_field_name"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: training polygons: " << e.what() << "\n";
            return 1;
        }
        std::int64_t n_poly[2] = {0, 0}, n_cell[2] = {0, 0};
        for (const lsdtt::training_polygon& tpoly : training.polygons)
        {
            ++n_poly[tpoly.cls];
            n_cell[tpoly.cls] += tpoly.n_cells;
        }
        std::cout << "  " << tp << ": " << n_poly[0] << " floodplain polygon(s) ("
                  << n_cell[0] << " cells), " << n_poly[1]
                  << " terrace polygon(s) (" << n_cell[1] << " cells)\n";
        if (level_supplied)
            std::cout << "  NOTE: supplied levels are used as given; the "
                         "training polygons do not set the ladder breaks.\n";

        std::vector<float> tl(n_pixels, -1.0f);
        for (std::size_t i = 0; i < n_pixels; ++i)
        {
            const std::int32_t f = training.polygon_of_pixel[i];
            if (f >= 0)
                tl[i] = static_cast<float>(
                    training.polygons[static_cast<std::size_t>(f)].cls);
        }
        lsdtt::lsdtt_raster_meta tmeta = meta;
        tmeta.no_data_value = -1.0;
        const std::string p = stem + "_training_levels.tif";
        std::cout << "  writing training levels (0 floodplain, 1 terrace): "
                  << p << "\n\n";
        loader.write_full(tl, tmeta, p, out_format);
    }
    auto add_training_fit_rows = [&](const std::string& scope,
                                     std::int32_t source, double centre,
                                     const std::string& used,
                                     const lsdtt::training_window_fit& fit)
    {
        auto row = [&](const std::string& name, const lsdtt::training_split& sp)
        {
            training_fit_rows.push_back({
                scope, std::to_string(source), num9(centre), used, name,
                num9(sp.threshold), num9(sp.youden_j),
                std::to_string(sp.n_lower), std::to_string(sp.n_upper),
                std::to_string(sp.lower_below), std::to_string(sp.upper_above),
                fit.note});
        };
        if (fit.has_bound) row("floodplain_bound", fit.bound);
        for (std::size_t b = 0; b < fit.break_fits.size(); ++b)
            row("break_" + std::to_string(b + 1), fit.break_fits[b]);
        if (!fit.has_bound && fit.break_fits.empty())
            training_fit_rows.push_back({scope, std::to_string(source),
                num9(centre), used, "none", "", "", "", "", "", "", fit.note});
    };

    // -------------------------------------------------------------------------
    // The ladder, band-then-label, and classification
    //
    // The hand field goes through classify_flat_surfaces_by_channel, which runs
    // the ladder and then assigns floodplain / terrace classes from height
    // above the channel. The fitted_surface field runs the ladder directly:
    // floodplain and terrace are statements about height above a channel and
    // mean nothing on a residual, so no class is invented for one.
    // -------------------------------------------------------------------------
    lsdtt::ladder_result                          ladder;
    lsdtt::banded_label_result                    banded;
    std::vector<lsdtt::flat_surface_label_class>  classes;
    lsdtt::lsdtt_raster_cpu                       class_raster;
    lsdtt::lsdtt_raster_cpu                       floodplain_mask;
    lsdtt::lsdtt_raster_cpu                       floodplain_first_pass_mask;

    std::cout << (level_fixed ? "--- Floodplain and terraces (fixed thresholds) ---\n"
                              : "--- Ladder and band-then-label ---\n");
    if (use_hand)
    {
        lsdtt::flat_surface_classify_request creq;
        creq.floodplain_max_height_metres =
            static_cast<double>(ff["max_floodplain_hand_metres"]);
        creq.terrace_min_height_metres =
            static_cast<double>(ff["min_terrace_hand_metres"]);
        creq.channel_adjacency_radius_metres =
            static_cast<double>(ff["channel_adjacency_radius_metres"]);
        creq.terrace_max_distance_to_channel_metres =
            static_cast<double>(ff["valley_corridor_half_width_metres"]);
        creq.terrace_corridor_per_pixel  = bf["clip_corridor_per_pixel"];
        creq.terrace_max_hand_spread_metres =
            static_cast<double>(ff["max_band_spread_metres"]);
        creq.n_threads                   = ip["n_threads"];

        creq.label_within_bands          = true;
        creq.segment_whole_mask          = segment_whole_mask;
        creq.whole_mask_barrier          = riser_barrier;
        creq.ladder_break_algorithm      = break_choice;
        creq.ladder_level_count          = ip["ladder_level_count"];
        creq.ladder_max_level_count      = ip["ladder_max_level_count"];
        creq.vertical_accuracy_metres =
            static_cast<double>(ff["vertical_accuracy_metres"]);
        creq.k_bootstrap      = ip["ladder_k_bootstrap"];
        creq.k_bootstrap_seed = ip["ladder_k_bootstrap_seed"];
        creq.vertical_accuracy_mode = ev_slope_dependent
            ? lsdtt::vertical_accuracy_choice::slope_dependent
            : lsdtt::vertical_accuracy_choice::fixed;
        creq.vertical_accuracy_gentle_metres =
            static_cast<double>(ff["vertical_accuracy_gentle_metres"]);
        creq.vertical_accuracy_steep_metres =
            static_cast<double>(ff["vertical_accuracy_steep_metres"]);
        creq.vertical_accuracy_break_gradient =
            static_cast<double>(ev_break_gradient);
        creq.resolvability_kappa =
            static_cast<double>(ff["ladder_resolvability_kappa"]);
        creq.band_reach_window_metres =
            static_cast<double>(ff["ladder_reach_window_metres"]);
        creq.band_min_window_pixels      = ip["ladder_min_window_pixels"];
        creq.band_hist_bins              = ip["ladder_hist_bins"];
        creq.band_min_hand_metres =
            static_cast<double>(ff["ladder_min_field_metres"]);
        creq.band_max_valley_fraction =
            static_cast<double>(ff["ladder_max_valley_fraction"]);
        {
            const std::string& gp = sf["ladder_gap_break_position"];
            if (gp == "lower")
                creq.ladder_gap_break_at = lsdtt::gap_break_position::lower;
            else if (gp == "upper")
                creq.ladder_gap_break_at = lsdtt::gap_break_position::upper;
            else if (gp == "midpoint")
                creq.ladder_gap_break_at = lsdtt::gap_break_position::midpoint;
            else
            {
                std::cerr << "FATAL: unknown ladder_gap_break_position '" << gp
                          << "'. Valid: lower, midpoint, upper.\n";
                return 1;
            }
        }
        {
            const std::string& gr = sf["ladder_gap_rank"];
            if (gr == "lowest")     creq.ladder_gap_rank = lsdtt::gap_rank_choice::lowest;
            else if (gr == "width") creq.ladder_gap_rank = lsdtt::gap_rank_choice::width;
            else
            {
                std::cerr << "FATAL: unknown ladder_gap_rank '" << gr
                          << "'. Valid: lowest, width.\n";
                return 1;
            }
        }
        creq.ladder_gap_trough_max_ratio =
            static_cast<double>(ff["ladder_gap_trough_max_ratio"]);
        creq.ladder_gap_empty_fraction =
            static_cast<double>(ff["ladder_gap_empty_fraction"]);
        creq.ladder_gap_min_width_metres =
            static_cast<double>(ff["ladder_gap_min_width_metres"]);
        creq.ladder_gap_hist_bins    = ip["ladder_gap_hist_bins"];
        creq.ladder_max_separation_ratio =
            static_cast<double>(ff["ladder_max_separation_ratio"]);
        creq.ladder_min_break_metres =
            static_cast<double>(ff["ladder_min_break_metres"]);
        creq.ladder_max_break_metres =
            static_cast<double>(ff["ladder_max_break_metres"]);
        creq.band_min_component_pixels   = ip["min_band_component_pixels"];
        creq.ladder_supplied_breaks      = supplied_breaks;
        creq.ladder_training             = training_on ? &training : nullptr;
        creq.ladder_training_spread_factor       = training_spread_factor;
        creq.ladder_training_min_cells_per_class = training_min_cells;
        creq.record_window_histograms    = bf["write_ladder_histogram_csv"] ||
                                           bf["write_ladder_break_metric_csv"];
        creq.window_hist_bins            = ip["ladder_hist_plot_bins"];
        creq.true_north                  = true_north.get();

        std::vector<std::int32_t> src_ids =
            scope_resolved ? scope.channel_source_id
                           : std::vector<std::int32_t>(n_pixels, -1);
        // break_flats_with_channels: drop channel cells from the admitted set.
        if (bf["break_flats_with_channels"] && scope_resolved)
        {
            std::vector<float>& ad = admitted.get_data();
            std::int64_t n_dropped = 0;
            for (std::size_t i = 0; i < n_pixels; ++i)
                if (scope.channel_mask[i] && ad[i] > 0.5f)
                { ad[i] = 0.0f; ++n_dropped; }
            std::cout << "  break_flats_with_channels: " << n_dropped
                      << " channel cell(s) removed from the admitted set\n";
        }

        std::vector<std::uint8_t> cmask =
            scope_resolved ? std::vector<std::uint8_t>(
                                 scope.channel_mask.begin(),
                                 scope.channel_mask.end())
                           : std::vector<std::uint8_t>(n_pixels, 0);

        // Slope-dependent E_v needs a gradient. Computed here rather than in
        // the engine so the fit window is a stated parameter and the same
        // measurement lsdtt-flat-surface-finder makes for admission.
        lsdtt::lsdtt_raster_cpu ev_slope;
        bool have_ev_slope = false;
        if (ev_slope_dependent)
        {
            const double cell = dem.get_meta().cell_size;
            // direct_polyfit takes the full physical window diameter.
            double lambda = 2.0 * static_cast<double>(
                ff["vertical_accuracy_gradient_window_radius_metres"]);
            if (!(lambda > 0.0)) lambda = 5.0 * cell;
            std::cout << "Window gradient for E_v: polyfit window " << lambda
                      << " m, break at "
                      << ev_break_gradient << " (gradient)\n";
            try
            {
                const lsdtt::poly_fit_coeffs fit =
                    lsdtt::direct_polyfit(dem, lambda, ip["n_threads"]);
                lsdtt::polyfit_metric_request mreq;
                mreq.metrics   = {lsdtt::polyfit_metric::slope};
                mreq.n_threads = ip["n_threads"];
                lsdtt::polyfit_metric_result mres =
                    lsdtt::compute_polyfit_metrics(fit, mreq);
                ev_slope = mres.at(lsdtt::polyfit_metric::slope);
                have_ev_slope = true;
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: gradient for slope-dependent E_v failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        lsdtt::lsdtt_raster_cpu patch_aspect;
        bool have_patch_aspect = false;
        if (bf["write_patch_terrain_aspect"])
        {
            const double cell = dem.get_meta().cell_size;
            // direct_polyfit takes the full physical window diameter.
            double lambda = 2.0 *
                static_cast<double>(ff["patch_aspect_window_radius_metres"]);
            if (!(lambda > 0.0)) lambda = 5.0 * cell;
            std::cout << "Terrain aspect for the patch table: polyfit window "
                      << lambda << " m\n";
            try
            {
                const lsdtt::poly_fit_coeffs fit =
                    lsdtt::direct_polyfit(dem, lambda, ip["n_threads"]);
                lsdtt::polyfit_metric_request mreq;
                mreq.metrics   = {lsdtt::polyfit_metric::aspect_degrees};
                mreq.n_threads = ip["n_threads"];
                lsdtt::polyfit_metric_result mres =
                    lsdtt::compute_polyfit_metrics(fit, mreq);
                patch_aspect = mres.at(lsdtt::polyfit_metric::aspect_degrees);
                have_patch_aspect = true;
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: terrain aspect failed: " << e.what()
                          << "\n";
                return 1;
            }
        }

        if (level_fixed)
        {
            // ---- fixed_thresholds: two levels from two bounds ----------------
            std::cout << "  (level_algorithm_choice fixed_thresholds: no ladder)\n";
            if (bf["break_flats_with_channels"])
                std::cout << "  NOTE: the channel pixels are seeded into the "
                             "floodplain, as Clubb et al. do; "
                             "break_flats_with_channels only affected admission.\n";

            // run_terrace_classification false: floodplain only. No terrace
            // field is built (not even a copy of HAND), and the engine never
            // reads one.
            const bool run_terraces = bf["run_terrace_classification"];
            if (!run_terraces)
                std::cout << "  run_terrace_classification false: floodplain "
                             "only, no terrace field.\n";
            lsdtt::lsdtt_raster_cpu terrace_field;
            if (run_terraces) terrace_field = reference_field;
            std::vector<std::int64_t> swath_nearest;
            if (run_terraces && terrace_field_choice == "swath_baseline")
            {
                std::vector<std::int64_t> baseline;
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (cmask[i]) baseline.push_back(static_cast<std::int64_t>(i));
                try
                {
                    lsdtt::baseline_relief_result br =
                        lsdtt::relief_above_nearest_baseline(
                            dem, baseline,
                            static_cast<double>(ff["swath_half_width_metres"]),
                            ip["n_threads"]);
                    terrace_field = std::move(br.relief);
                    swath_nearest = std::move(br.nearest_baseline_flat);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: swath baseline failed: " << e.what()
                              << "\n";
                    return 1;
                }
                std::cout << "  terrace field: relief above the nearest of "
                          << baseline.size() << " channel pixel(s) within "
                          << ff["swath_half_width_metres"] << " m\n";
            }
            if (run_terraces && bf["write_terrace_field_raster"])
            {
                const std::string p = stem + "_terrace_field.tif";
                std::cout << "Writing terrace field: " << p << "\n";
                loader.write_full(terrace_field.get_data(),
                                  terrace_field.get_meta(), p, out_format);
            }

            // Training: a HAND bound per reach window (the ladder's windows,
            // pooled k-1..k+1), else the fit over all training cells, else
            // the data-chosen bound.
            std::vector<float> trained_bound;
            if (training_on)
            {
                lsdtt::training_set tset = training;
                const double ev = ev_slope_dependent
                    ? static_cast<double>(ff["vertical_accuracy_gentle_metres"])
                    : static_cast<double>(ff["vertical_accuracy_metres"]);
                const double gap =
                    static_cast<double>(ff["ladder_resolvability_kappa"]) * ev;
                lsdtt::compute_training_polygon_stats(
                    tset, reference_field, training_spread_factor * gap);
                training_polygon_stats = tset.polygons;

                const std::vector<float>&   rf  = reference_field.get_data();
                const std::vector<uint8_t>& rfm = reference_field.get_mask();
                auto cell_ok = [&](std::size_t i)
                { return rfm[i] && std::isfinite(rf[i]); };
                std::vector<std::pair<std::int32_t, double>> all_cells;
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (tset.polygon_of_pixel[i] >= 0 && cell_ok(i))
                        all_cells.emplace_back(tset.polygon_of_pixel[i],
                                               static_cast<double>(rf[i]));
                const lsdtt::training_window_fit gfit = lsdtt::fit_training_window(
                    all_cells, tset, gap, training_min_cells);
                add_training_fit_rows("global", -1, 0.0,
                                      gfit.has_bound ? "fallback" : "no", gfit);
                const float gbound = gfit.has_bound
                    ? static_cast<float>(gfit.bound.threshold)
                    : std::numeric_limits<float>::quiet_NaN();

                trained_bound.assign(n_pixels, gbound);
                std::int64_t n_window = 0, n_global = 0;
                const double reach_w =
                    static_cast<double>(ff["ladder_reach_window_metres"]);
                if (have_donors && scope_resolved && reach_w >= 0.0)
                {
                    const lsdtt::donor_reach_coordinate rc =
                        lsdtt::donor_reach_coordinates(donor_flat, flow_distance,
                                                       src_ids);
                    std::vector<std::uint8_t> defined(n_pixels, 0);
                    for (std::size_t i = 0; i < n_pixels; ++i)
                        defined[i] = cell_ok(i) ? 1 : 0;
                    const lsdtt::reach_window_set rw = lsdtt::assign_reach_windows(
                        rc.source_of_pixel, rc.distance, defined, reach_w);
                    std::vector<std::vector<std::pair<std::int32_t, double>>>
                        wcells(rw.windows.size());
                    for (std::size_t i = 0; i < n_pixels; ++i)
                    {
                        const std::int32_t w = rw.window_of_pixel[i];
                        if (w >= 0 && tset.polygon_of_pixel[i] >= 0)
                            wcells[static_cast<std::size_t>(w)].emplace_back(
                                tset.polygon_of_pixel[i], static_cast<double>(rf[i]));
                    }
                    std::vector<float> wbound(rw.windows.size(), gbound);
                    for (std::size_t w = 0; w < rw.windows.size(); ++w)
                    {
                        std::vector<std::pair<std::int32_t, double>> pooled;
                        for (const std::int32_t k : rw.windows[w].neighbours)
                            pooled.insert(pooled.end(),
                                          wcells[static_cast<std::size_t>(k)].begin(),
                                          wcells[static_cast<std::size_t>(k)].end());
                        if (pooled.empty()) continue;
                        const lsdtt::training_window_fit fit =
                            lsdtt::fit_training_window(pooled, tset, gap,
                                                       training_min_cells);
                        if (fit.has_bound)
                            wbound[w] = static_cast<float>(fit.bound.threshold);
                        add_training_fit_rows(
                            "window", rw.windows[w].source, rw.windows[w].centre_m,
                            fit.has_bound ? "window"
                                          : (gfit.has_bound ? "global" : "no"),
                            fit);
                        fit.has_bound ? ++n_window : ++n_global;
                    }
                    for (std::size_t i = 0; i < n_pixels; ++i)
                    {
                        const std::int32_t w = rw.window_of_pixel[i];
                        if (w >= 0) trained_bound[i] = wbound[static_cast<std::size_t>(w)];
                    }
                }
                std::cout << "  training: floodplain / terrace HAND bound ";
                if (gfit.has_bound)
                    std::cout << gfit.bound.threshold << " m over all training "
                                 "(J " << gfit.bound.youden_j << ")";
                else
                    std::cout << "not fitted over all training ("
                              << gfit.note << ")";
                std::cout << "; " << n_window << " window(s) with their own "
                          << "bound, " << n_global << " with training but "
                             "too little of it\n";
                if (!gfit.has_bound && n_window == 0)
                {
                    std::cout << "  NOTE: the training sets no bound; the "
                                 "data-chosen bound is used.\n";
                    trained_bound.clear();
                }
            }

            lsdtt::floodplain_threshold_request freq;
            freq.max_floodplain_hand_metres = floodplain_hand_bound;
            if (!trained_bound.empty()) freq.hand_bound_per_cell = &trained_bound;
            freq.max_floodplain_gradient    = 0.0;   // applied in admission
            freq.seed_channel_pixels        = true;
            freq.require_channel_contact    = true;
            freq.min_component_pixels       = ip["min_band_component_pixels"];
            freq.fill_holes                 = fill_floodplain_holes;
            freq.connectivity               = lsdtt::pixel_connectivity::d8;
            freq.find_terraces              = run_terraces;
            freq.min_terrace_height_above_floodplain_metres = static_cast<double>(
                ff["min_terrace_height_above_floodplain_metres"]);
            freq.max_terrace_height_metres =
                static_cast<double>(ff["max_terrace_height_metres"]);

            lsdtt::floodplain_threshold_result fres;
            try
            {
                fres = lsdtt::classify_floodplain_terraces_by_thresholds(
                    admitted, reference_field, reference_field, cmask,
                    run_terraces ? &terrace_field : nullptr, freq);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: floodplain / terrace thresholds: "
                          << e.what() << "\n";
                return 1;
            }
            // ---- optional second pass on the first-pass floodplain's HAND ----
            const std::string& second_pass =
                sf["floodplain_second_pass_algorithm_choice"];
            if (second_pass != "none" && second_pass != "qq_upper_departure")
            {
                std::cerr << "FATAL: unknown floodplain_second_pass_algorithm_"
                             "choice '" << second_pass << "'. Valid: none, "
                             "qq_upper_departure.\n";
                return 1;
            }
            if (second_pass == "qq_upper_departure" && !trained_bound.empty())
                std::cout << "  NOTE: the training sets the HAND bound; the "
                             "second pass is skipped.\n";
            else if (second_pass == "qq_upper_departure")
            {
                std::cout << "  first pass: " << fres.n_floodplain_pixels
                          << " floodplain cell(s) at HAND <= "
                          << floodplain_hand_bound << " m\n";
                const std::vector<float>& fpd = fres.floodplain_mask.get_data();
                const std::vector<float>& rf  = reference_field.get_data();
                const std::vector<uint8_t>& rfm = reference_field.get_mask();
                std::vector<double> fp_hand;
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (fpd[i] > 0.5f && rfm[i] && std::isfinite(rf[i]))
                        fp_hand.push_back(static_cast<double>(rf[i]));
                lsdtt::threshold_selection_result sel2;
                try
                {
                    sel2 = lsdtt::select_qq_upper_departure(fp_hand,
                        lsdtt::threshold_request_from_parameters(
                            "qq_gaussian", "floodplain_second_pass_algorithm_choice",
                            0.0, ff, "floodplain"));
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: floodplain second pass: " << e.what()
                              << "\n";
                    return 1;
                }
                if (bf["write_floodplain_threshold_csv"] && !sel2.qq.data_quantile.empty())
                {
                    const std::string p = stem + "_floodplain_hand_qq_second_pass.csv";
                    lsdtt::write_qq_csv(loader, p, sel2.qq);
                    std::cout << "  Writing " << p << "\n";
                }
                if (!sel2.found || !(sel2.threshold < floodplain_hand_bound))
                {
                    std::cout << "  second pass: "
                              << (sel2.found ? "the upper departure is not below "
                                               "the first-pass bound"
                                             : "no upper departure found")
                              << "; the first-pass floodplain is kept.\n";
                }
                else
                {
                    std::cout << "  second pass (qq_upper_departure): HAND bound "
                              << sel2.threshold << " m (quantile " << sel2.qq_index
                              << " of " << sel2.qq.data_quantile.size()
                              << " over " << fp_hand.size()
                              << " first-pass floodplain cell(s))\n";
                    floodplain_threshold_rows.push_back({
                        "2", "hand", "qq_upper_departure",
                        num9(sel2.threshold), std::to_string(sel2.qq_index),
                        std::to_string(fp_hand.size()),
                        fname_stem +
                            "_floodplain_hand_qq_second_pass.csv"});
                    floodplain_first_pass_mask = fres.floodplain_mask;
                    floodplain_hand_bound = static_cast<float>(sel2.threshold);
                    freq.max_floodplain_hand_metres = floodplain_hand_bound;
                    freq.fill_holes = true;   // the second pass is always filled
                    try
                    {
                        fres = lsdtt::classify_floodplain_terraces_by_thresholds(
                            admitted, reference_field, reference_field, cmask,
                            run_terraces ? &terrace_field : nullptr, freq);
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << "FATAL: floodplain second pass: "
                                  << e.what() << "\n";
                        return 1;
                    }
                }
            }

            std::cout << "  floodplain: " << fres.n_floodplain_candidates
                      << " candidate cell(s) (channel included); kept "
                      << fres.n_floodplain_components_kept
                      << " component(s), dropped "
                      << fres.n_floodplain_dropped_small << " small and "
                      << fres.n_floodplain_dropped_no_channel
                      << " not touching the channel; "
                      << fres.n_floodplain_pixels << " cell(s)\n"
                      << "  terraces: " << fres.n_terrace_components
                      << " component(s), " << fres.n_terrace_pixels
                      << " cell(s), above " << floodplain_hand_bound +
                             ff["min_terrace_height_above_floodplain_metres"]
                      << " m on the terrace field\n";

            // Each patch's own field: HAND on the floodplain, the terrace field
            // on a terrace - the quantity its level was decided on.
            const std::vector<float>& bi = fres.band_index.get_data();
            std::vector<float>   pf(n_pixels, 0.0f);
            std::vector<uint8_t> pfm(n_pixels, 0);
            std::vector<float>   on(n_pixels, 0.0f);
            std::int64_t n_riser_removed = 0;
            for (std::size_t i = 0; i < n_pixels; ++i)
            {
                if (bi[i] < -0.5f) continue;
                // whole_mask: one terrace level, so the risers only split it.
                if (riser_barrier != nullptr && riser_mask.get_data()[i] > 0.5f)
                {
                    ++n_riser_removed;
                    continue;
                }
                on[i] = 1.0f;
                const lsdtt::lsdtt_raster_cpu& f =
                    (bi[i] < 0.5f) ? reference_field : terrace_field;
                if (f.get_mask()[i] && std::isfinite(f.get_data()[i]))
                {
                    pf[i]  = f.get_data()[i];
                    pfm[i] = 1;
                }
            }
            if (riser_barrier != nullptr)
                std::cout << "  whole_mask: " << n_riser_removed
                          << " floodplain / terrace cell(s) removed as riser\n";
            const lsdtt::lsdtt_raster_cpu patch_field(
                std::move(pf), std::move(pfm), meta);
            const lsdtt::lsdtt_raster_cpu on_level(
                std::move(on), std::vector<uint8_t>(dem.get_mask()), meta);
            try
            {
                // Size was already filtered on the whole floodplain / terrace
                // masks, so labelling keeps every component, 8-connected as
                // the engine counted them.
                banded = lsdtt::label_components_within_bands(
                    fres.band_index, on_level, patch_field, 1,
                    lsdtt::pixel_connectivity::d8,
                    have_patch_aspect ? &patch_aspect : nullptr,
                    true_north.get());
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: labelling failed: " << e.what() << "\n";
                return 1;
            }

            // Class records, one per patch, in the classifier's schema.
            std::vector<float> chan_bin(n_pixels, 0.0f);
            for (std::size_t i = 0; i < n_pixels; ++i)
                if (cmask[i]) chan_bin[i] = 1.0f;
            const lsdtt::lsdtt_raster_cpu chan_dist =
                lsdtt::euclidean_distance_transform(
                    lsdtt::lsdtt_raster_cpu(std::move(chan_bin),
                                            std::vector<uint8_t>(dem.get_mask()),
                                            meta),
                    ip["n_threads"]);
            const std::vector<float>& lab = banded.labels.get_data();
            const std::size_t n_patch = banded.patches.size();
            std::vector<std::vector<double>> heights(n_patch), cdist(n_patch);
            std::vector<std::map<std::int32_t, std::int64_t>> src_votes(n_patch);
            std::vector<double> min_dist(n_patch,
                                         std::numeric_limits<double>::infinity());
            for (std::size_t i = 0; i < n_pixels; ++i)
            {
                if (!(lab[i] > 0.5f)) continue;
                const std::size_t k = static_cast<std::size_t>(lab[i]) - 1;
                if (k >= n_patch) continue;
                if (patch_field.get_mask()[i])
                    heights[k].push_back(patch_field.get_data()[i]);
                if (chan_dist.get_mask()[i])
                    min_dist[k] = std::min(min_dist[k],
                        static_cast<double>(chan_dist.get_data()[i]));
                // The channel a pixel belongs to: the swath node it was
                // measured from on a terrace, its HAND donor otherwise.
                std::int64_t node = -1;
                if (bi[i] > 0.5f && !swath_nearest.empty())
                    node = swath_nearest[i];
                else if (have_donors)
                    node = donor_flat[i];
                if (node >= 0 && static_cast<std::size_t>(node) < n_pixels)
                {
                    const std::int32_t sid =
                        scope.channel_source_id[static_cast<std::size_t>(node)];
                    if (sid >= 0) ++src_votes[k][sid];
                    if (have_donors && donor_flat[i] >= 0)
                        cdist[k].push_back(flow_distance.get_data()[
                            static_cast<std::size_t>(donor_flat[i])]);
                }
            }
            classes.clear();
            classes.reserve(n_patch);
            const double cell = meta.cell_size;
            for (std::size_t k = 0; k < n_patch; ++k)
            {
                const lsdtt::patch_record& pr = banded.patches[k];
                lsdtt::flat_surface_label_class c;
                c.label     = pr.label;
                c.cls       = pr.band == 0 ? lsdtt::surface_class::floodplain
                                           : lsdtt::surface_class::terrace;
                c.n_pixels  = pr.n_pixels;
                c.hand_band = pr.band;
                std::int64_t best = 0;
                for (const auto& [sid, nv] : src_votes[k])
                    if (nv > best) { best = nv; c.source_id = sid; }
                std::vector<double>& h = heights[k];
                if (!h.empty())
                {
                    std::sort(h.begin(), h.end());
                    c.has_height    = true;
                    c.min_height    = h.front();
                    c.max_height    = h.back();
                    c.median_height = (h.size() % 2)
                        ? h[h.size() / 2]
                        : 0.5 * (h[h.size() / 2 - 1] + h[h.size() / 2]);
                    double sum = 0.0;
                    for (double v : h) sum += v;
                    c.mean_height = sum / static_cast<double>(h.size());
                    c.hand_spread_metres = c.max_height - c.min_height;
                }
                if (std::isfinite(min_dist[k]))
                {
                    c.min_distance_to_channel_metres = min_dist[k];
                    c.adjacent_to_channel = min_dist[k] <= cell * 1.5;
                }
                std::vector<double>& cd = cdist[k];
                if (!cd.empty())
                {
                    std::sort(cd.begin(), cd.end());
                    c.median_channel_distance_metres = cd[cd.size() / 2];
                }
                c.ladder_window_centre_metres = -1.0;
                classes.push_back(c);
            }

            std::vector<float> cls(n_pixels, 0.0f);
            for (std::size_t i = 0; i < n_pixels; ++i)
            {
                if (bi[i] > -0.5f && bi[i] < 0.5f)
                    cls[i] = static_cast<float>(lsdtt::surface_class::floodplain);
                else if (bi[i] > 0.5f)
                    cls[i] = static_cast<float>(lsdtt::surface_class::terrace);
            }
            class_raster = lsdtt::lsdtt_raster_cpu(
                std::move(cls), std::vector<uint8_t>(dem.get_mask()), meta);
            floodplain_mask   = fres.floodplain_mask;
            ladder.band_index = fres.band_index;
            ladder.max_band_index = 1;

            std::int64_t n_fp = 0, n_tr = 0;
            for (const auto& c : classes)
                (c.cls == lsdtt::surface_class::floodplain ? n_fp : n_tr)++;
            std::cout << "  " << n_fp << " floodplain, " << n_tr
                      << " terrace surface(s)\n";
        }
        else
        {
        lsdtt::flat_surface_classify_result cres;
        try
        {
            cres = lsdtt::classify_flat_surfaces_by_channel(
                admitted, reference_field, cmask, src_ids, creq,
                have_donors ? &donor_flat  : nullptr,
                have_donors ? &flow_distance : nullptr,
                have_ev_slope ? &ev_slope : nullptr,
                have_patch_aspect ? &patch_aspect : nullptr);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: classification failed: " << e.what() << "\n";
            return 1;
        }

        ladder.windows       = cres.ladder_windows;
        ladder.max_band_index = 0;
        if (training_on && !level_supplied)
        {
            training_polygon_stats = cres.ladder_training_polygons;
            const lsdtt::training_window_fit& g = cres.ladder_training_global_fit;
            add_training_fit_rows("global", -1, 0.0,
                                  g.breaks.empty() ? "no" : "fallback", g);
            for (const lsdtt::ladder_window& w : cres.ladder_windows)
                if (w.source_id >= 0)
                    add_training_fit_rows(
                        "window", w.source_id, w.window_centre_m,
                        w.training_source.empty() ? "no" : w.training_source,
                        w.training_fit);
            std::cout << "  training: " << cres.ladder_n_windows_trained
                      << " of " << cres.ladder_windows.size()
                      << " window(s) banded with trained breaks";
            if (!g.breaks.empty())
            {
                std::cout << "; over all training:";
                for (const double b : g.breaks) std::cout << " " << b;
                std::cout << " m";
            }
            std::cout << "\n";
        }
        banded.labels        = cres.banded_labels_raster;
        banded.patches       = cres.patches;
        classes              = cres.labels;
        class_raster         = cres.surface_class_raster;
        floodplain_mask      = cres.floodplain_mask;
        ladder.band_index    = cres.hand_band_raster;
        ladder.window_index  = cres.ladder_window_index;

        if (segment_whole_mask)
            std::cout << "  whole_mask: " << cres.n_whole_mask_barrier_pixels
                      << " admitted cell(s) removed as riser, "
                      << cres.n_whole_mask_reassigned_pixels
                      << " cell(s) moved to their patch's majority band\n";
        std::cout << "  " << cres.n_floodplain_labels << " floodplain, "
                  << cres.n_terrace_labels << " terrace surface(s)\n";
        }
    }
    else
    {
        lsdtt::ladder_request lreq;
        lreq.reference_field  = field_choice;
        lreq.break_algorithm  = break_choice;
        lreq.level_count      = ip["ladder_level_count"];
        lreq.max_level_count  = ip["ladder_max_level_count"];
        lreq.vertical_accuracy_metres =
            static_cast<double>(ff["vertical_accuracy_metres"]);
        lreq.k_bootstrap      = ip["ladder_k_bootstrap"];
        lreq.k_bootstrap_seed = ip["ladder_k_bootstrap_seed"];
        lreq.vertical_accuracy_mode = ev_slope_dependent
            ? lsdtt::vertical_accuracy_choice::slope_dependent
            : lsdtt::vertical_accuracy_choice::fixed;
        lreq.vertical_accuracy_gentle_metres =
            static_cast<double>(ff["vertical_accuracy_gentle_metres"]);
        lreq.vertical_accuracy_steep_metres =
            static_cast<double>(ff["vertical_accuracy_steep_metres"]);
        lreq.vertical_accuracy_break_gradient =
            static_cast<double>(ev_break_gradient);
        lreq.resolvability_kappa =
            static_cast<double>(ff["ladder_resolvability_kappa"]);
        lreq.reach_window_metres = -1.0;   // no reach on a fitted surface (C6)
        lreq.supplied_breaks     = supplied_breaks;
        lreq.min_window_pixels   = ip["ladder_min_window_pixels"];
        lreq.hist_bins           = ip["ladder_hist_bins"];
        lreq.min_field_metres =
            static_cast<double>(ff["ladder_min_field_metres"]);
        lreq.max_valley_fraction =
            static_cast<double>(ff["ladder_max_valley_fraction"]);
        {
            const std::string& gr = sf["ladder_gap_rank"];
            if (gr == "lowest")     lreq.gap_rank = lsdtt::gap_rank_choice::lowest;
            else if (gr == "width") lreq.gap_rank = lsdtt::gap_rank_choice::width;
            else
            {
                std::cerr << "FATAL: unknown ladder_gap_rank '" << gr
                          << "'. Valid: lowest, width.\n";
                return 1;
            }
        }
        lreq.gap_trough_max_ratio =
            static_cast<double>(ff["ladder_gap_trough_max_ratio"]);
        lreq.gap_empty_fraction =
            static_cast<double>(ff["ladder_gap_empty_fraction"]);
        lreq.gap_min_width_metres =
            static_cast<double>(ff["ladder_gap_min_width_metres"]);
        lreq.gap_hist_bins        = ip["ladder_gap_hist_bins"];
        lreq.max_separation_ratio =
            static_cast<double>(ff["ladder_max_separation_ratio"]);
        lreq.min_break_metres =
            static_cast<double>(ff["ladder_min_break_metres"]);
        lreq.max_break_metres =
            static_cast<double>(ff["ladder_max_break_metres"]);
        lreq.n_threads = ip["n_threads"];
        lreq.record_window_histograms = bf["write_ladder_histogram_csv"] ||
                                        bf["write_ladder_break_metric_csv"];
        lreq.window_hist_bins         = ip["ladder_hist_plot_bins"];

        try
        {
            ladder = lsdtt::build_reference_field_ladder(
                reference_field, admitted, lreq);
            lsdtt::lsdtt_raster_cpu label_admitted = admitted;
            if (segment_whole_mask)
            {
                lsdtt::whole_mask_segmentation seg = lsdtt::segment_whole_mask(
                    ladder.band_index, admitted, riser_barrier,
                    lsdtt::pixel_connectivity::d4);
                std::cout << "  whole_mask: " << seg.n_barrier_pixels
                          << " admitted cell(s) removed as riser, "
                          << seg.n_reassigned_pixels
                          << " cell(s) moved to their patch's majority band\n";
                ladder.band_index = std::move(seg.band_index);
                label_admitted    = std::move(seg.admitted);
            }
            banded = lsdtt::label_components_within_bands(
                ladder.band_index, label_admitted, reference_field,
                std::max(1, ip["min_band_component_pixels"]),
                lsdtt::pixel_connectivity::d4, nullptr, true_north.get());
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: ladder failed: " << e.what() << "\n";
            return 1;
        }
        std::cout << "  A fitted-surface residual carries no floodplain / "
                     "terrace class: those are statements about height above "
                     "a channel.\n"
                     "  The band index and the patch table are the result.\n";
    }

    // The resolvability report, on stdout as well as in the CSV: which branch
    // the ladder took decides what may be claimed about it.
    if (level_fixed)
    {
        std::cout << "  " << banded.patches.size() << " patch(es)\n";
    }
    else
    {
        std::int64_t n_res = 0, n_unres = 0, n_deg = 0, n_sup = 0, n_sup_low = 0;
        for (const lsdtt::ladder_window& w : ladder.windows)
        {
            if (w.status == lsdtt::ladder_status::resolved ||
                w.status == lsdtt::ladder_status::forced) ++n_res;
            else if (w.status == lsdtt::ladder_status::unresolved) ++n_unres;
            else if (w.status == lsdtt::ladder_status::supplied)
            {
                ++n_sup;
                // Two measured bands at least (an empty band has no centre).
                if (w.min_band_centre_gap_m > 0.0 &&
                    w.min_band_centre_gap_m < w.required_gap_m) ++n_sup_low;
            }
            else ++n_deg;
        }
        std::cout << "  " << banded.patches.size() << " patch(es)\n"
                  << "  ladder: " << n_res << " window(s) resolved, "
                  << n_unres << " unresolved, " << n_deg << " degenerate";
        if (n_sup > 0)
            std::cout << ", " << n_sup << " with supplied levels ("
                      << n_sup_low << " below the resolvability floor)";
        if (!ladder.windows.empty())
        {
            // Under slope_dependent the floor is per window, so reporting one
            // number would be a lie about the other windows. Print the range
            // when they differ and the single value when they do not.
            double gmin = ladder.windows[0].required_gap_m;
            double gmax = gmin;
            for (const auto& w : ladder.windows)
            {
                gmin = std::min(gmin, w.required_gap_m);
                gmax = std::max(gmax, w.required_gap_m);
            }
            if (gmax > gmin)
                std::cout << "  (required band-centre gap = " << gmin
                          << " to " << gmax << " m across windows)";
            else
                std::cout << "  (required band-centre gap = " << gmin << " m)";
        }
        std::cout << "\n";
        for (const std::string& note : ladder.notes)
            std::cout << "  NOTE: " << note << "\n";
        if (n_unres > 0)
            std::cout << "  NOTE: an unresolved window means the DEM cannot "
                         "separate the requested levels there; do not describe "
                         "them as resolved levels.\n";
    }
    std::cout << "\n";

    // -------------------------------------------------------------------------
    // Outputs
    // -------------------------------------------------------------------------
    std::cout << "--- Writing outputs ---\n";

    if (bf["write_band_index_raster"] && !ladder.band_index.get_data().empty())
    {
        const std::string p = stem + "_band_index.tif";
        std::cout << "Writing band index (Int16): " << p << "\n";
        loader.write_full(ladder.band_index.get_data(),
                          ladder.band_index.get_meta(), p, out_format,
                          lsdtt::lsdtt_pixel_type::int16);
    }

    if (bf["write_banded_labels_raster"] && !banded.labels.get_data().empty())
    {
        const std::string p = stem + "_banded_labels.tif";
        std::cout << "Writing banded labels (UInt32): " << p << "\n";
        loader.write_full(banded.labels.get_data(), banded.labels.get_meta(),
                          p, out_format, lsdtt::lsdtt_pixel_type::uint32);
    }

    if (bf["write_valley_mask_raster"] && !ladder.band_index.get_data().empty())
    {
        // Band 0 is the valley floor by construction: the ladder is ordered, so
        // the lowest band is the lowest surface. With ladder_level_count 2 this
        // is exactly the valley/wall cut lsdtt-valley-centreline wants.
        const std::vector<float>& bi = ladder.band_index.get_data();
        std::vector<float>   vm(n_pixels, 0.0f);
        std::vector<uint8_t> vmm(n_pixels, 1);
        std::int64_t n_valley = 0;

        const double force_h = static_cast<double>(ff["force_valley_hand_metres"]);
        const double force_d =
            static_cast<double>(ff["force_valley_max_distance_metres"]);

        if (force_h > 0.0 || force_d > 0.0)
        {
            // Stated, not discovered. The ladder still ran, and the band index,
            // labels and patch table are untouched - only this raster is
            // forced, so the two answers can be laid side by side.
            const std::vector<float>&   rf  = reference_field.get_data();
            const std::vector<uint8_t>& rfm = reference_field.get_mask();
            const std::vector<float>&   am  = admitted.get_data();

            std::vector<float> dist;
            if (force_d > 0.0 && scope_resolved)
            {
                std::vector<float>   cdv(n_pixels, 0.0f);
                std::vector<uint8_t> cmv(n_pixels, 1u);
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (scope.channel_mask[i]) cdv[i] = 1.0f;
                const lsdtt::lsdtt_raster_cpu chan(
                    std::move(cdv), std::move(cmv), meta);
                dist = lsdtt::euclidean_distance_transform(
                    chan, ip["n_threads"]).get_data();
            }
            else if (force_d > 0.0)
            {
                std::cerr << "WARNING: force_valley_max_distance_metres needs a "
                             "channel scope; the distance limit is ignored.\n";
            }

            for (std::size_t i = 0; i < n_pixels; ++i)
            {
                if (!(am[i] > 0.5f)) continue;
                if (!rfm[i] || !std::isfinite(rf[i])) continue;
                if (force_h > 0.0 && static_cast<double>(rf[i]) > force_h) continue;
                if (!dist.empty() && static_cast<double>(dist[i]) > force_d) continue;
                vm[i] = 1.0f; ++n_valley;
            }
            std::cout << "  valley floor FORCED: ";
            if (force_h > 0.0) std::cout << "field <= " << force_h << " m";
            if (force_h > 0.0 && force_d > 0.0) std::cout << " and ";
            if (force_d > 0.0) std::cout << "within " << force_d
                                         << " m of a channel";
            std::cout << " -> " << n_valley << " cell(s)\n";
        }
        else
        for (std::size_t i = 0; i < n_pixels; ++i)
        {
            if (bi[i] >= -0.5f && bi[i] < 0.5f) { vm[i] = 1.0f; ++n_valley; }
        }
        lsdtt::lsdtt_raster_cpu valley_raster(
            std::move(vm), std::move(vmm), meta);

        if (bf["clean_valley_mask"] && scope_resolved)
        {
            // Generic constrained morphology seeded by the channel network -
            // erode, drop what the channel cannot reach, regrow. Nothing about
            // it is valley-specific, so it lives with the rest of the binary
            // mask operations rather than here.
            lsdtt::seed_reachability_options sropt;
            sropt.dilation_pixels = ip["valley_dilation_pixels"];
            sropt.require_seed_adjacent =
                bf["restrict_valley_to_channel_connected"];
            sropt.max_detour_ratio =
                static_cast<double>(ff["valley_max_detour_ratio"]);
            sropt.require_seed_line_of_sight =
                bf["restrict_valley_to_line_of_sight"];
            sropt.n_threads = ip["n_threads"];
            try
            {
                std::vector<float>   sd(n_pixels, 0.0f);
                std::vector<uint8_t> smk(n_pixels, 1u);
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (scope.channel_mask[i]) sd[i] = 1.0f;
                const lsdtt::lsdtt_raster_cpu seed_raster(
                    std::move(sd), std::move(smk), meta);

                valley_raster =
                    lsdtt::filter_binary_components_by_seed_reachability(
                        valley_raster, seed_raster, sropt);
                std::int64_t n_after = 0;
                for (const float v : valley_raster.get_data())
                    if (v > 0.5f) ++n_after;
                std::cout << "  valley mask cleaned: " << n_valley << " -> "
                          << n_after << " cell(s)\n";
                n_valley = n_after;
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: valley mask cleaning failed: "
                          << e.what() << "\n";
            }
        }
        else if (bf["clean_valley_mask"])
        {
            std::cout << "  valley mask not cleaned: cleaning needs a channel "
                         "network, and no channel scope was resolved.\n";
        }

        // Closing BEFORE hole filling, and both after the connectivity filter:
        // closing turns an open-ended sliver into either solid ground or an
        // enclosed hole, and the fill then takes whatever is left.
        {
            double cr =
                static_cast<double>(ff["valley_mask_closing_radius_metres"]);

            // A fraction of the MEASURED width, when asked for. Measured on the
            // mask as it stands now - after connectivity, before closing - so
            // the radius is derived from the valley the run actually found.
            const double wfrac =
                static_cast<double>(ff["valley_mask_closing_width_fraction"]);
            if (wfrac > 0.0 && scope_resolved)
            {
                // Measure on a HOLE-FILLED copy, not on the mask as it stands.
                // The slivers this closing exists to seal are exactly what
                // makes the cumulative width undercount: a section crossing
                // floor / sliver / floor counts only the floor, so the width
                // comes out narrower than the valley and the radius derived
                // from it comes out too small to seal the slivers that caused
                // the undercount. Filling first breaks that circularity.
                lsdtt::lsdtt_raster_cpu width_mask = valley_raster;
                try { width_mask = lsdtt::binary_fill_holes(width_mask); }
                catch (const std::exception&) { /* measure on what we have */ }

                std::vector<double> widths;
                for (const std::vector<std::int64_t>& nodes :
                         scope.channel_node_lists)
                {
                    const lsdtt::mask_width_result mw =
                        lsdtt::measure_mask_width_across_polyline(
                            width_mask, nodes,
                            static_cast<double>(
                                ff["valley_width_section_spacing_metres"]),
                            static_cast<double>(
                                ff["valley_width_max_half_length_metres"]));
                    if (mw.n_sections > 0)
                    {
                        std::cout << "  valley width, channel "
                                  << widths.size() << ": floor "
                                  << mw.median_width_metres
                                  << " m, envelope "
                                  << mw.median_envelope_metres << " m (broken "
                                  << (mw.median_width_metres > 0.0
                                      ? mw.median_envelope_metres /
                                        mw.median_width_metres : 0.0)
                                  << "x, " << mw.n_sections << " sections)\n";
                        // The closing is trying to reconstitute the ENVELOPE,
                        // so that is the width its radius is a fraction of.
                        widths.push_back(mw.median_envelope_metres);
                    }
                }
                if (!widths.empty())
                {
                    std::sort(widths.begin(), widths.end());
                    const double w = widths[widths.size() / 2];
                    cr = wfrac * w;
                    std::cout << "  closing radius = " << wfrac << " x " << w
                              << " m = " << cr << " m\n";
                }
                else
                {
                    std::cerr << "WARNING: no usable cross section, so the "
                                 "valley width could not be measured; falling "
                                 "back to valley_mask_closing_radius_metres.\n";
                }
            }
            if (cr > 0.0)
            {
                try
                {
                    valley_raster =
                        lsdtt::binary_morphological_closing(valley_raster, cr);
                    std::int64_t n_closed = 0;
                    for (const float v : valley_raster.get_data())
                        if (v > 0.5f) ++n_closed;
                    std::cout << "  valley mask closed (radius " << cr
                              << " m): " << n_valley << " -> " << n_closed
                              << " cell(s)\n";
                    n_valley = n_closed;
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: closing the valley mask failed: "
                              << e.what() << "\n";
                }
            }
        }

        if (bf["fill_valley_mask_holes"])
        {
            try
            {
                valley_raster = lsdtt::binary_fill_holes(valley_raster);
                std::int64_t n_filled = 0;
                for (const float v : valley_raster.get_data())
                    if (v > 0.5f) ++n_filled;
                std::cout << "  valley mask holes filled: " << n_valley
                          << " -> " << n_filled << " cell(s)\n";
                n_valley = n_filled;
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: filling valley mask holes failed: "
                          << e.what() << "\n";
            }
        }

        const std::string p = stem + "_valley_mask.tif";
        std::cout << "Writing valley mask (band 0, " << n_valley
                  << " cell(s)): " << p << "\n";
        loader.write_full(valley_raster.get_data(), meta, p, out_format,
                          lsdtt::lsdtt_pixel_type::uint8);
    }

    if (use_hand && bf["write_surface_class_polygons"] &&
        !class_raster.get_data().empty())
    {
        lsdtt::vector_output_format pfmt;
        try
        {
            pfmt = lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
        }
        catch (const std::exception&)
        {
            std::cerr << "FATAL: write_surface_class_polygons needs "
                         "write_vector_format flatgeobuf or gpkg, not '"
                      << sf["write_vector_format"] << "'.\n";
            return 1;
        }
        // One label per connected area of each class; the class travels as
        // the polygon attribute.
        const std::vector<float>& cv = class_raster.get_data();
        std::vector<float>        labels(n_pixels, 0.0f);
        std::vector<std::int32_t> class_of(1, -1);
        std::int64_t next = 0;
        for (const int c : {1, 2})
        {
            std::vector<float> m(n_pixels, 0.0f);
            for (std::size_t i = 0; i < n_pixels; ++i)
                if (static_cast<int>(cv[i] + 0.5f) == c) m[i] = 1.0f;
            const lsdtt::lsdtt_raster_cpu lab = lsdtt::label_binary_components_by_size(
                lsdtt::lsdtt_raster_cpu(std::move(m),
                                        std::vector<uint8_t>(n_pixels, 1), meta),
                1, lsdtt::pixel_connectivity::d8);
            std::int64_t max_here = 0;
            for (std::size_t i = 0; i < n_pixels; ++i)
            {
                const float l = lab.get_data()[i];
                if (l > 0.5f)
                {
                    const std::int64_t li = static_cast<std::int64_t>(l + 0.5f);
                    labels[i] = static_cast<float>(next + li);
                    max_here = std::max(max_here, li);
                }
            }
            for (std::int64_t k = 0; k < max_here; ++k) class_of.push_back(c);
            next += max_here;
        }
        lsdtt::lsdtt_vectorise_options vopt;
        vopt.layer_name           = "surface_classes";
        vopt.format               = pfmt;
        vopt.background           = 0.0f;
        vopt.spatial_index        = true;
        vopt.label_attribute_name = "class";
        const std::string p = stem + "_surface_class_polygons" +
                              lsdtt::output_extension(pfmt);
        std::cout << "Writing surface class polygons: " << p << "\n";
        try
        {
            lsdtt::lsdtt_data_loader vec_loader;
            lsdtt::trace_labelled_polygons(
                lsdtt::lsdtt_raster_cpu(std::move(labels),
                                        std::vector<uint8_t>(n_pixels, 1), meta),
                vec_loader, p, vopt, class_of);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: surface class polygons: " << e.what() << "\n";
            return 1;
        }
    }

    if (use_hand && bf["write_surface_class_raster"] &&
        !class_raster.get_data().empty())
    {
        const std::string p = stem + "_surface_class.tif";
        std::cout << "Writing surface class raster: " << p << "\n";
        loader.write_full(class_raster.get_data(), class_raster.get_meta(),
                          p, out_format, lsdtt::lsdtt_pixel_type::uint8);
    }

    if (use_hand && bf["write_floodplain_mask_raster"] &&
        !floodplain_mask.get_data().empty())
    {
        const std::string p = stem + "_floodplain_mask.tif";
        std::cout << "Writing floodplain mask: " << p << "\n";
        loader.write_full(floodplain_mask.get_data(),
                          floodplain_mask.get_meta(), p, out_format, lsdtt::lsdtt_pixel_type::uint8);
    }

    if (use_hand && bf["write_floodplain_mask_raster"] &&
        !floodplain_first_pass_mask.get_data().empty())
    {
        const std::string p = stem + "_floodplain_first_pass_mask.tif";
        std::cout << "Writing first-pass floodplain mask: " << p << "\n";
        loader.write_full(floodplain_first_pass_mask.get_data(),
                          floodplain_first_pass_mask.get_meta(), p, out_format,
                          lsdtt::lsdtt_pixel_type::uint8);
    }

    if (training_on)
    {
        const std::string pf = stem + "_training_fit.csv";
        std::cout << "Writing training fit: " << pf << "\n";
        loader.write_csv(pf, {"scope", "source_id", "window_centre_m", "used",
                              "threshold", "threshold_m", "youden_j",
                              "n_lower", "n_upper", "n_lower_below",
                              "n_upper_above", "note"},
                         training_fit_rows);
        std::vector<std::vector<std::string>> prow;
        for (std::size_t k = 0; k < training_polygon_stats.size(); ++k)
        {
            const lsdtt::training_polygon& tp = training_polygon_stats[k];
            const bool has = tp.n_with_field > 0;
            prow.push_back({std::to_string(k), tp.label,
                            tp.cls == 0 ? "floodplain" : "terrace",
                            std::to_string(tp.number), std::to_string(tp.n_cells),
                            std::to_string(tp.n_with_field),
                            has ? num9(tp.median) : "", has ? num9(tp.p10) : "",
                            has ? num9(tp.p90) : "",
                            has ? num9(tp.p90 - tp.p10) : "",
                            tp.flagged ? "1" : "0"});
        }
        const std::string pp2 = stem + "_training_polygons.csv";
        std::cout << "Writing training polygons: " << pp2 << "\n";
        loader.write_csv(pp2, {"polygon_id", "label", "class", "terrace_number",
                               "n_cells", "n_with_field", "median_m", "p10_m",
                               "p90_m", "spread_m", "flagged"},
                         prow);
    }

    if (bf["write_floodplain_threshold_csv"] && !floodplain_threshold_rows.empty())
    {
        const std::string p = stem + "_floodplain_thresholds.csv";
        std::cout << "Writing floodplain thresholds: " << p << "\n";
        loader.write_csv(p, {"pass", "quantity", "algorithm", "threshold",
                             "qq_index", "n_values", "table_fname"},
                         floodplain_threshold_rows);
    }

    if (level_fixed && bf["write_floodplain_polygons"] &&
        !floodplain_mask.get_data().empty())
    {
        lsdtt::vector_output_format pfmt;
        try
        {
            pfmt = lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
        }
        catch (const std::exception&)
        {
            std::cerr << "FATAL: write_floodplain_polygons needs "
                         "write_vector_format flatgeobuf or gpkg, not '"
                      << sf["write_vector_format"] << "'.\n";
            return 1;
        }
        auto write_polygons = [&](const lsdtt::lsdtt_raster_cpu& mask,
                                  const std::string& suffix,
                                  const std::string& layer)
        {
            const lsdtt::lsdtt_raster_cpu labels =
                lsdtt::label_binary_components_by_size(
                    mask, 1, lsdtt::pixel_connectivity::d8);
            lsdtt::lsdtt_vectorise_options vopt;
            vopt.layer_name = layer;
            vopt.format     = pfmt;
            vopt.background = 0.0f;
            vopt.spatial_index = true;   // as lsdtt-marsh-topography
            const std::string p = stem + suffix + lsdtt::output_extension(pfmt);
            std::cout << "Writing " << layer << ": " << p << "\n";
            lsdtt::lsdtt_data_loader vec_loader;
            lsdtt::trace_labelled_polygons(labels, vec_loader, p, vopt);
        };
        try
        {
            write_polygons(floodplain_mask, "_floodplain_polygons",
                           "floodplain");
            if (!floodplain_first_pass_mask.get_data().empty())
                write_polygons(floodplain_first_pass_mask,
                               "_floodplain_first_pass_polygons",
                               "floodplain_first_pass");
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: floodplain polygons: " << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_surface_patch_polygons"] && !banded.labels.get_data().empty())
    {
        lsdtt::vector_output_format pfmt;
        try
        {
            pfmt = lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
        }
        catch (const std::exception&)
        {
            std::cerr << "FATAL: write_surface_patch_polygons needs "
                         "write_vector_format flatgeobuf or gpkg, not '"
                      << sf["write_vector_format"] << "'.\n";
            return 1;
        }
        std::int64_t max_label = 0;
        for (const auto& pr : banded.patches) max_label = std::max(max_label, pr.label);
        std::vector<std::int32_t> band_of(static_cast<std::size_t>(max_label) + 1, -1);
        for (const auto& pr : banded.patches)
            band_of[static_cast<std::size_t>(pr.label)] = pr.band;
        lsdtt::lsdtt_vectorise_options vopt;
        vopt.layer_name           = "surface_patches";
        vopt.format               = pfmt;
        vopt.background           = 0.0f;
        vopt.spatial_index        = true;
        vopt.label_attribute_name = "band";
        const std::string p = stem + "_surface_patch_polygons" +
                              lsdtt::output_extension(pfmt);
        std::cout << "Writing surface patch polygons: " << p << "\n";
        try
        {
            lsdtt::lsdtt_data_loader vec_loader;
            lsdtt::trace_labelled_polygons(banded.labels, vec_loader, p, vopt, band_of);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: surface patch polygons: " << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_surface_patches_csv"])
    {
        const std::string p = stem + "_surface_patches.csv";
        std::cout << "Writing surface patches: " << p << "\n";
        lsdtt::write_surface_patches_csv(loader, p, banded.patches, classes);
    }

    if (use_hand && bf["write_surface_class_csv"])
    {
        const std::string p = stem + "_surface_class.csv";
        std::cout << "Writing surface class table: " << p << "\n";
        lsdtt::write_surface_class_csv(loader, p, classes);
    }

    if (bf["write_ladder_windows_csv"] && level_fixed)
    {
        std::cout << "No ladder windows: level_algorithm_choice is "
                     "fixed_thresholds.\n";
    }
    else if (bf["write_ladder_windows_csv"])
    {
        const std::string p = stem + "_ladder_windows.csv";
        std::cout << "Writing ladder windows: " << p << "\n";
        lsdtt::write_ladder_windows_csv(loader, p, ladder.windows);
    }

    if (bf["write_nearest_source_id_raster"])
    {
        if (!have_donors || !scope_resolved)
        {
            std::cerr << "WARNING: write_nearest_source_id_raster needs a "
                         "computed HAND (a supplied hand_fname carries no donor "
                         "nodes), so no raster was written.\n";
        }
        else
        {
            std::vector<float>   sid(n_pixels, -1.0f);
            std::vector<uint8_t> sidm(n_pixels, 1u);
            for (std::size_t i = 0; i < n_pixels; ++i)
            {
                const std::int64_t d = donor_flat[i];
                if (d < 0) continue;
                const std::size_t di = static_cast<std::size_t>(d);
                if (di >= n_pixels) continue;
                sid[i] = static_cast<float>(scope.channel_source_id[di]);
            }
            const std::string p = stem + "_nearest_source_id.tif";
            std::cout << "Writing nearest source id (Int32): " << p << "\n";
            loader.write_full(sid, meta, p, out_format, lsdtt::lsdtt_pixel_type::int32);
        }
    }

    if (bf["write_ladder_window_id_raster"] &&
        !ladder.window_index.get_data().empty())
    {
        const std::string p = stem + "_ladder_window_id.tif";
        std::cout << "Writing ladder window id (Int16): " << p << "\n";
        loader.write_full(ladder.window_index.get_data(),
                          ladder.window_index.get_meta(), p, out_format,
                          lsdtt::lsdtt_pixel_type::int16);
    }

    if (bf["write_ladder_break_metric_csv"])
    {
        const std::string p = stem + "_ladder_break_metric.csv";
        std::cout << "Writing break decision metric: " << p << "\n";
        lsdtt::write_ladder_break_metric_csv(loader, p, ladder.windows);
    }

    if (bf["write_ladder_histogram_csv"])
    {
        const std::string p = stem + "_ladder_histogram.csv";
        std::cout << "Writing ladder densities: " << p << "\n";
        lsdtt::write_ladder_histogram_csv(loader, p, ladder.windows);
    }

    if (bf["write_reference_field_profile_csv"])
    {
        // Per-pixel field against along-channel distance. The distance is read
        // at the pixel's HAND DONOR node, not at the pixel, which is why this
        // needs a computed HAND: a supplied HAND raster carries no donors, and
        // a fitted surface has no channel to measure along at all.
        if (!have_donors || !scope_resolved)
        {
            std::cerr << "WARNING: write_reference_field_profile_csv needs a "
                         "computed HAND (it reads the along-channel distance at "
                         "each pixel's donor node), so no profile was written.\n";
        }
        else
        {
            const std::string p = stem + "_reference_field_profile.csv";
            std::cout << "Writing reference field profile: " << p << "\n";

            const std::vector<float>&   fv  = reference_field.get_data();
            const std::vector<uint8_t>& fm  = reference_field.get_mask();
            const std::vector<float>&   fdd = flow_distance.get_data();
            const std::vector<uint8_t>& fdm = flow_distance.get_mask();
            const std::vector<float>&   lv  = banded.labels.get_data();
            const std::vector<float>&   bv  = ladder.band_index.get_data();
            const std::vector<float>*   cv  =
                class_raster.get_data().empty() ? nullptr
                                                : &class_raster.get_data();

            std::vector<std::vector<std::string>> rows;
            for (std::size_t i = 0; i < n_pixels; ++i)
            {
                if (!(lv[i] > 0.0f)) continue;          // not a labelled patch
                if (!fm[i] || !std::isfinite(fv[i])) continue;
                const std::int64_t donor = donor_flat[i];
                if (donor < 0) continue;
                const std::size_t d = static_cast<std::size_t>(donor);
                if (d >= n_pixels || !fdm[d] || !std::isfinite(fdd[d])) continue;

                const int r = static_cast<int>(i / static_cast<std::size_t>(meta.n_cols));
                const int c = static_cast<int>(i % static_cast<std::size_t>(meta.n_cols));
                const lsdtt::lsdtt_xy xy =
                    lsdtt::cell_center_from_row_col(meta, r, c);
                rows.push_back({
                    std::to_string(xy.x),
                    std::to_string(xy.y),
                    std::to_string(fv[i]),
                    std::to_string(fdd[d]),
                    std::to_string(scope.channel_source_id[d]),
                    cv ? std::to_string(static_cast<int>((*cv)[i]))
                       : std::string("-1"),
                    std::to_string(static_cast<std::int64_t>(std::llround(lv[i]))),
                    std::to_string(static_cast<int>(bv[i])) });
            }
            loader.write_csv(p,
                {"x","y","field_m","channel_flow_distance_m",
                 "source_id","class","label","band"},
                rows);
        }
    }

    if (bf["write_channel_network"] && scope_resolved && flow)
    {
        std::cout << "\n--- Writing channel network ---\n";
        try
        {
            lsdtt::lsdtt_junction_network_options jo;
            jo.build_topological_order = true;
            jo.build_s_vector          = true;
            std::unique_ptr<lsdtt::lsdtt_junction_network_base> jnet =
                lsdtt::make_junction_network_base(
                    *flow, scope.resolved_sources, jo);

            lsdtt::channel_network_write_request nreq;
            nreq.write_lines           = bf["write_channel_network_lines"];
            nreq.write_channel_points  = bf["write_channel_network_points"];
            nreq.write_junction_points = bf["write_junction_network_points"];
            nreq.format = lsdtt::resolve_vector_output_format(
                sf["write_vector_format"]);
            nreq.channel_source_id = scope.channel_source_id;

            lsdtt::write_channel_network_layers(*jnet, stem, nreq);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: channel network failed: "
                      << e.what() << "\n";
        }
    }

    // ---- method record (plan 7a) and what the comparison needs -------------
    {
        std::vector<std::vector<std::string>> rec;
        auto add = [&](const std::string& item, const std::string& value,
                       const std::string& source)
        { rec.push_back({method, item, value, source}); };
        for (const stage_value& sv : stages) add(sv.key, sv.value, sv.source);
        for (const auto& r : floodplain_threshold_rows)
            add("pass" + r[0] + "_" + r[1] + "_bound", r[3], r[2]);
        if (low_relief_fraction >= 0.0)
        {
            std::ostringstream os; os.precision(6); os << low_relief_fraction;
            add("low_relief_fraction", os.str(), "measured");
            add("floodplain_bound_algorithm_choice",
                use_floodplain_gradient_bound ? "hand_and_gradient" : "hand_only",
                "auto");
        }
        if (n_hand_clamped > 0)
            add("negative_hand_cells_set_to_0", std::to_string(n_hand_clamped),
                "antimode ladder");
        std::int64_t n_res = 0, n_unres = 0, n_deg = 0, n_sup = 0;
        for (const lsdtt::ladder_window& w : ladder.windows)
        {
            if (w.status == lsdtt::ladder_status::resolved ||
                w.status == lsdtt::ladder_status::forced) ++n_res;
            else if (w.status == lsdtt::ladder_status::unresolved) ++n_unres;
            else if (w.status == lsdtt::ladder_status::supplied) ++n_sup;
            else ++n_deg;
        }
        if (!ladder.windows.empty())
        {
            add("ladder_windows_resolved", std::to_string(n_res), "measured");
            add("ladder_windows_unresolved", std::to_string(n_unres), "measured");
            add("ladder_windows_degenerate", std::to_string(n_deg), "measured");
            if (n_sup > 0)
                add("ladder_windows_supplied", std::to_string(n_sup), "measured");
        }
        int max_band = -1;
        for (const auto& pr : banded.patches) max_band = std::max(max_band, pr.band);
        add("patches", std::to_string(banded.patches.size()), "measured");
        add("levels", std::to_string(max_band + 1), "measured");

        std::int64_t n_fp = 0, n_tr = 0;
        run_out.meta = meta;
        run_out.valid.assign(n_pixels, 0);
        for (std::size_t i = 0; i < n_pixels; ++i)
            run_out.valid[i] = dem.get_mask()[i] ? 1 : 0;
        if (use_hand && !class_raster.get_data().empty())
        {
            run_out.cls.assign(n_pixels, 0);
            const std::vector<float>& cv = class_raster.get_data();
            for (std::size_t i = 0; i < n_pixels; ++i)
            {
                const int c = static_cast<int>(cv[i] + 0.5f);
                if (c == 1 || c == 2) run_out.cls[i] = static_cast<std::uint8_t>(c);
                if (c == 1) ++n_fp;
                if (c == 2) ++n_tr;
            }
        }
        const double cell_area = meta.cell_size * meta.cell_size;
        add("floodplain_cells", std::to_string(n_fp), "measured");
        add("terrace_cells", std::to_string(n_tr), "measured");

        if (bf["write_method_record_csv"])
        {
            const std::string p = stem + "_method_record.csv";
            std::cout << "Writing method record: " << p << "\n";
            loader.write_csv(p, {"method", "item", "value", "source"}, rec);
        }
        auto num = [](double v) { std::ostringstream os; os.precision(9); os << v; return os.str(); };
        run_out.summary = {
            {"floodplain_cells", std::to_string(n_fp)},
            {"floodplain_area_m2", num(n_fp * cell_area)},
            {"terrace_cells", std::to_string(n_tr)},
            {"terrace_area_m2", num(n_tr * cell_area)},
            {"patches", std::to_string(banded.patches.size())},
            {"levels", std::to_string(max_band + 1)},
            {"ladder_windows_resolved", std::to_string(n_res)},
            {"ladder_windows_unresolved", std::to_string(n_unres)},
            {"ladder_windows_degenerate", std::to_string(n_deg)}};
        for (const auto& r : floodplain_threshold_rows)
            run_out.summary.push_back({"pass" + r[0] + "_" + r[1] + "_bound", r[3]});
    }
    return 0;
    };   // run_method

    std::vector<method_run_output> runs(method_list.size());
    const bool several = method_list.size() > 1;
    for (std::size_t k = 0; k < method_list.size(); ++k)
    {
        const std::string& m = method_list[k];
        const std::string run_stem  = several
            ? pp.get_full_write_path() + "_" + m : pp.get_full_write_path();
        const std::string run_fname = several
            ? pp.get_write_fname() + "_" + m : pp.get_write_fname();
        if (several)
            std::cout << "\n================ method " << (k + 1) << " of "
                      << method_list.size() << ": " << m
                      << " ================\n";
        const int rc = run_method(m, run_stem, run_fname, runs[k]);
        if (rc != 0) return rc;
    }

    // ---- comparison of the methods (plan 7b) -----------------------------------
    if (several)
    {
        const std::string base = pp.get_full_write_path();
        lsdtt::lsdtt_data_loader out_loader;
        {
            std::vector<std::string> head = {"method"};
            for (const auto& kv : runs[0].summary) head.push_back(kv.first);
            std::vector<std::vector<std::string>> rows;
            for (std::size_t k = 0; k < runs.size(); ++k)
            {
                std::vector<std::string> row = {method_list[k]};
                for (const std::string& h : std::vector<std::string>(head.begin() + 1, head.end()))
                {
                    std::string v;
                    for (const auto& kv : runs[k].summary) if (kv.first == h) v = kv.second;
                    row.push_back(v);
                }
                rows.push_back(row);
            }
            // Bounds only some methods have go in extra columns.
            for (std::size_t k = 0; k < runs.size(); ++k)
                for (const auto& kv : runs[k].summary)
                    if (std::find(head.begin(), head.end(), kv.first) == head.end())
                    {
                        head.push_back(kv.first);
                        for (auto& r : rows) r.push_back("");
                    }
            for (std::size_t k = 0; k < runs.size(); ++k)
                for (const auto& kv : runs[k].summary)
                {
                    const auto it = std::find(head.begin(), head.end(), kv.first);
                    rows[k][static_cast<std::size_t>(it - head.begin())] = kv.second;
                }
            const std::string p = base + "_method_summary.csv";
            std::cout << "\nWriting method summary: " << p << "\n";
            out_loader.write_csv(p, head, rows);
        }

        std::vector<lsdtt::method_class_map> maps;
        for (std::size_t k = 0; k < runs.size(); ++k)
            if (!runs[k].cls.empty()) maps.push_back({method_list[k], runs[k].cls});
        if (maps.size() < 2)
        {
            std::cout << "Fewer than two methods produced floodplain / terrace "
                         "classes (hand field only); no comparison.\n";
        }
        else
        {
            const auto pairs = lsdtt::compare_surface_classes(maps, runs[0].valid);
            const double ca = runs[0].meta.cell_size * runs[0].meta.cell_size;
            auto num = [](double v) { std::ostringstream os; os.precision(6); os << v; return os.str(); };
            std::vector<std::vector<std::string>> rows;
            const char* cname[3] = {"none", "floodplain", "terrace"};
            for (const auto& pc : pairs)
            {
                std::vector<std::string> row = {maps[pc.a].name, maps[pc.b].name,
                    num(pc.floodplain_iou), num(pc.terrace_iou),
                    num(pc.classified_agreement)};
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        row.push_back(std::to_string(pc.confusion[i][j]));
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        row.push_back(num(pc.confusion[i][j] * ca));
                rows.push_back(row);
            }
            std::vector<std::string> head = {"method_a", "method_b",
                "floodplain_iou", "terrace_iou", "classified_agreement"};
            for (const char* suffix : {"_cells", "_m2"})
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        head.push_back(std::string("a_") + cname[i] + "_b_" +
                                       cname[j] + suffix);
            const std::string p = base + "_method_comparison.csv";
            std::cout << "Writing method comparison: " << p << "\n";
            out_loader.write_csv(p, head, rows);

            std::cout << "Agreement rasters: bit k = method k in this order:";
            for (std::size_t k = 0; k < maps.size(); ++k)
                std::cout << " " << k << "=" << maps[k].name;
            std::cout << "\n";
            for (const auto& [cls, name] : std::vector<std::pair<std::uint8_t, std::string>>{
                     {1, "floodplain"}, {2, "terrace"}})
            {
                const std::vector<std::uint8_t> bits =
                    lsdtt::surface_class_agreement_bits(maps, runs[0].valid, cls);
                std::vector<float> f(bits.begin(), bits.end());
                const std::string rp = base + "_method_agreement_" + name + ".tif";
                std::cout << "Writing " << rp << "\n";
                out_loader.write_full(f, runs[0].meta, rp, out_format,
                                      lsdtt::lsdtt_pixel_type::uint8);
            }
        }
    }

    std::cout << "\n";
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-fluvial-surfaces complete.\n";
    return run_record.succeed();
}
