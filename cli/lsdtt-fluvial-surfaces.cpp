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

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_true_north.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_channel_scope_cpu.hpp"
#include "../src/cpu/lsdtt_flat_surface_channel_classify_cpu.hpp"
#include "../src/cpu/lsdtt_proximity_cpu.hpp"
#include "../src/cpu/lsdtt_raster_binary_mask_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_reference_field_ladder_cpu.hpp"
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

    // --- Admission: which ground may take part ---
    string_default_map["flat_mask_fname"] = "";
    help_map["flat_mask_fname"] = {"string","",
        "Flat mask raster from lsdtt-flat-surface-finder, in read_path.",
        "Cells > 0.5 are admitted. This is normally "
        "<stem>_combined_flat_mask.tif from lsdtt-flat-surface-finder, which "
        "is the intended two-step workflow: detect flat ground there, sort it "
        "into levels here. Leave empty to use the built-in fallback admission "
        "instead, which is coarser and is meant for a quick look rather than "
        "for a result you intend to publish.",
        "Admission","io",""};

    string_default_map["fallback_admission_algorithm_choice"] = "relief";
    help_map["fallback_admission_algorithm_choice"] = {"string","relief",
        "Built-in admission when flat_mask_fname is empty: relief or slope.",
        "relief admits cells whose local relief over "
        "fallback_relief_window_radius_pixels is below "
        "fallback_relief_threshold_metres. slope admits cells whose gradient "
        "magnitude is below fallback_slope_threshold.",
        "Admission","param","flat_mask_fname"};

    int_default_map["fallback_relief_window_radius_pixels"] = 5;
    help_map["fallback_relief_window_radius_pixels"] = {"int","5",
        "Moving-window RADIUS in pixels for the fallback relief admission.",
        "This is a radius, not a window side: 5 means an 11x11 window.  Note "
        "that lsdtt-flat-surface-finder and lsdtt-mountain-front spell their "
        "equivalent relief_window_side_pixels, which is a full SIDE LENGTH and "
        "so means something different for the same number.  Only used when "
        "flat_mask_fname is empty and the fallback is relief.",
        "Admission","param","flat_mask_fname"};

    float_default_map["fallback_relief_threshold_metres"] = 1.0f;
    help_map["fallback_relief_threshold_metres"] = {"float","1.0",
        "Local relief below which a cell is admitted as flat [m].",
        "Only used when flat_mask_fname is empty and the fallback is relief.",
        "Admission","param","flat_mask_fname"};

    float_default_map["fallback_slope_threshold"] = 0.05f;
    help_map["fallback_slope_threshold"] = {"float","0.05",
        "Gradient magnitude below which a cell is admitted as flat [m/m].",
        "Only used when flat_mask_fname is empty and the fallback is slope.",
        "Admission","param","flat_mask_fname"};

    int_default_map["min_component_pixels"] = 50;
    help_map["min_component_pixels"] = {"int","50",
        "Minimum patch size in pixels, applied PER BAND.",
        "Note the change of meaning from lsdtt-flat-surface-finder, where the "
        "same name filters components of the whole mask. Here components are "
        "formed inside each band, so this is a floor on the size of one "
        "surface on one level.",
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

    int_default_map["ladder_max_level_count"] = 12;
    help_map["ladder_max_level_count"] = {"int","12",
        "Ceiling on the level count considered under auto selection.",
        "A ceiling on what one reach can plausibly carry, not a storage limit. "
        "A reach resolving more than a dozen distinct fluvial levels is more "
        "likely over-fitting noise than recording twelve terraces. Raise it "
        "freely for a genuine long Quaternary staircase; there is no numeric "
        "cliff behind it.",
        "Ladder","param",""};

    string_default_map["ladder_break_algorithm_choice"] = "antimode";
    help_map["ladder_break_algorithm_choice"] = {"string","antimode",
        "Break finder: gap, antimode or jenks.",
        "antimode finds dips in a smoothed log10 density and is shared with "
        "the valley engine. jenks is Fisher-Jenks on the linear field; it "
        "takes the level count directly and is the source of the reported "
        "goodness-of-variance fit. NOTE: antimode bins on log10 and therefore "
        "cannot band a field with negative values, which every "
        "fitted_surface residual has. This tool switches to jenks "
        "automatically when you select fitted_surface without naming a break "
        "algorithm, and rejects the combination if you ask for it explicitly.",
        "Ladder","param",""};

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

    string_default_map["vertical_accuracy_mode"] = "fixed";
    help_map["vertical_accuracy_mode"] = {"string","fixed",
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

    float_default_map["vertical_accuracy_break_degrees"] = 11.31f;
    help_map["vertical_accuracy_break_degrees"] = {"float","11.31",
        "Gradient [degrees] separating the gentle and steep E_v.",
        "11.31 degrees is a 20 percent grade, which is where the TanDEM-X "
        "accuracy specification splits.",
        "Ladder", "param", ""};

    float_default_map["vertical_accuracy_slope_window_radius_metres"] = 0.0f;
    help_map["vertical_accuracy_slope_window_radius_metres"] = {"float","0",
        "Window radius [m] of the polynomial fit the window gradient is "
        "measured over. 0 = five cells. slope_dependent only.",
        "Stated in cells by default because a length tuned on lidar is "
        "sub-pixel on 30 m COP30. This is the same measurement "
        "lsdtt-flat-surface-finder makes for admission_slope_cap_degrees, and "
        "it is worth using the same window in both so the gradient a surface "
        "was admitted on is the gradient its E_v was read from.",
        "Ladder", "param", ""};

    float_default_map["vertical_accuracy_metres"] = 0.0f;
    help_map["vertical_accuracy_metres"] = {"float","0.0",
        "Vertical accuracy of the DEM, E_v [m]. 0 disables the floor.",
        "Two levels are only claimed as distinct when their band centres are "
        "at least ladder_resolvability_kappa * E_v apart. This is metadata you "
        "must supply: it cannot be derived from the DEM. 0 disables the check "
        "and every window then reports resolved, carrying no instrument "
        "guarantee at all. Note that a DEM whose vertical units are feet on a "
        "metre grid has an effective E_v 3.28 times larger than you think.",
        "Ladder","param",""};

    float_default_map["ladder_resolvability_kappa"] = 3.0f;
    help_map["ladder_resolvability_kappa"] = {"float","3.0",
        "Multiplier on E_v in the resolvability floor.",
        "Required band-centre gap = kappa * vertical_accuracy_metres.",
        "Ladder","param",""};

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

    // --- Classification (hand field) ---
    float_default_map["floodplain_max_height_metres"] = 1.0f;
    help_map["floodplain_max_height_metres"] = {"float","1.0",
        "Upper HAND bound of the floodplain class [m].",
        "An optional absolute override on the lowest band's upper bound. Under "
        "band-then-label the ladder already decides where one level ends and "
        "the next begins; this is a hard cap for cases where you know the "
        "floodplain height independently. Only meaningful for the hand field.",
        "Classification","param",""};

    float_default_map["terrace_min_height_metres"] = 0.0f;
    help_map["terrace_min_height_metres"] = {"float","0.0",
        "Lower HAND bound of the terrace class [m].",
        "As floodplain_max_height_metres, at the other end. Only meaningful "
        "for the hand field.",
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

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp(cli_args.path, cli_args.fname);

    // A pure rename: the value carries across unchanged. The new spelling says
    // what the number is - a radius, so 5 is an 11x11 window - because the same
    // idea is a full window SIDE in lsdtt-flat-surface-finder and
    // lsdtt-mountain-front.
    pp.set_deprecated_aliases({
        {"fallback_relief_window_pixels", "fallback_relief_window_radius_pixels"},
        {"source_snap_max_radius_pixels", "source_snap_max_search_radius_pixels"},
        {"fill_dem_for_flow_routing", "fill_before_routing"},
        {"channel_threshold_contributing_pixels",
         "threshold_contributing_pixels"} });
    pp.set_deprecated_numeric_aliases({
        {"patch_aspect_window_metres",
         {"patch_aspect_window_radius_metres", 0.5, 0.0, true}},
        {"vertical_accuracy_slope_window_metres",
         {"vertical_accuracy_slope_window_radius_metres", 0.5, 0.0, true}} });

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

    const std::string& ev_mode_name = sf["vertical_accuracy_mode"];
    const bool ev_slope_dependent = (ev_mode_name == "slope_dependent");
    if (ev_mode_name != "fixed" && !ev_slope_dependent)
    {
        std::cerr << "FATAL: unknown vertical_accuracy_mode '" << ev_mode_name
                  << "'. Valid: fixed, slope_dependent.\n";
        return 1;
    }
    if (ev_slope_dependent)
    {
        const double eg = static_cast<double>(ff["vertical_accuracy_gentle_metres"]);
        const double es = static_cast<double>(ff["vertical_accuracy_steep_metres"]);
        if (!(eg > 0.0) || !(es > 0.0))
        {
            std::cerr << "FATAL: vertical_accuracy_mode slope_dependent needs "
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

    if (!ev_slope_dependent && ff["vertical_accuracy_metres"] <= 0.0f)
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
    const std::string stem = pp.get_full_write_path();

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

    // -------------------------------------------------------------------------
    // Admission: which ground may take part
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu admitted;
    {
        std::cout << "--- Admission ---\n";
        if (!sf["flat_mask_fname"].empty())
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
            const std::string& how = sf["fallback_admission_algorithm_choice"];
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
                const double radius =
                    static_cast<double>(ip["fallback_relief_window_radius_pixels"]) *
                    meta.cell_size;
                const lsdtt::lsdtt_raster_cpu relief =
                    lsdtt::moving_window_local_relief(dem, radius);
                const std::vector<float>& rv = relief.get_data();
                const float thr = ff["fallback_relief_threshold_metres"];
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (dm[i] && std::isfinite(rv[i]) && rv[i] < thr)
                        ad[i] = 1.0f;
            }
            else if (how == "slope")
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
                const float thr = ff["fallback_slope_threshold"];
                for (std::size_t i = 0; i < n_pixels; ++i)
                    if (dm[i] && std::isfinite(sv[i]) && sv[i] < thr)
                        ad[i] = 1.0f;
            }
            else
            {
                std::cerr << "FATAL: unknown "
                             "fallback_admission_algorithm_choice '" << how
                          << "'. Valid: relief, slope.\n";
                return 1;
            }
            admitted = lsdtt::lsdtt_raster_cpu(
                std::move(ad), std::move(am), meta);
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

            scope = lsdtt::resolve_channel_scope(*flow, sreq);
            scope_resolved = true;
            std::cout << "  " << scope.resolved_sources.size()
                      << " channel source(s) in scope\n";

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

    std::cout << "--- Ladder and band-then-label ---\n";
    if (use_hand)
    {
        lsdtt::flat_surface_classify_request creq;
        creq.floodplain_max_height_metres =
            static_cast<double>(ff["floodplain_max_height_metres"]);
        creq.terrace_min_height_metres =
            static_cast<double>(ff["terrace_min_height_metres"]);
        creq.channel_adjacency_radius_metres =
            static_cast<double>(ff["channel_adjacency_radius_metres"]);
        creq.terrace_max_distance_to_channel_metres =
            static_cast<double>(ff["valley_corridor_half_width_metres"]);
        creq.terrace_corridor_per_pixel  = bf["clip_corridor_per_pixel"];
        creq.terrace_max_hand_spread_metres =
            static_cast<double>(ff["max_band_spread_metres"]);
        creq.n_threads                   = ip["n_threads"];

        creq.label_within_bands          = true;
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
        creq.vertical_accuracy_break_degrees =
            static_cast<double>(ff["vertical_accuracy_break_degrees"]);
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
        creq.band_min_component_pixels   = ip["min_component_pixels"];
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
                ff["vertical_accuracy_slope_window_radius_metres"]);
            if (!(lambda > 0.0)) lambda = 5.0 * cell;
            std::cout << "Window gradient for E_v: polyfit window " << lambda
                      << " m, break at "
                      << ff["vertical_accuracy_break_degrees"] << " degrees\n";
            try
            {
                const lsdtt::poly_fit_coeffs fit =
                    lsdtt::direct_polyfit(dem, lambda, ip["n_threads"]);
                lsdtt::polyfit_metric_request mreq;
                mreq.metrics   = {lsdtt::polyfit_metric::slope_degrees};
                mreq.n_threads = ip["n_threads"];
                lsdtt::polyfit_metric_result mres =
                    lsdtt::compute_polyfit_metrics(fit, mreq);
                ev_slope = mres.at(lsdtt::polyfit_metric::slope_degrees);
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
        banded.labels        = cres.banded_labels_raster;
        banded.patches       = cres.patches;
        classes              = cres.labels;
        class_raster         = cres.surface_class_raster;
        floodplain_mask      = cres.floodplain_mask;
        ladder.band_index    = cres.hand_band_raster;
        ladder.window_index  = cres.ladder_window_index;

        std::cout << "  " << cres.n_floodplain_labels << " floodplain, "
                  << cres.n_terrace_labels << " terrace surface(s)\n";
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
        lreq.vertical_accuracy_break_degrees =
            static_cast<double>(ff["vertical_accuracy_break_degrees"]);
        lreq.resolvability_kappa =
            static_cast<double>(ff["ladder_resolvability_kappa"]);
        lreq.reach_window_metres = -1.0;   // no reach on a fitted surface (C6)
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
            banded = lsdtt::label_components_within_bands(
                ladder.band_index, admitted, reference_field,
                std::max(1, ip["min_component_pixels"]),
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
    {
        std::int64_t n_res = 0, n_unres = 0, n_deg = 0;
        for (const lsdtt::ladder_window& w : ladder.windows)
        {
            if (w.status == lsdtt::ladder_status::resolved ||
                w.status == lsdtt::ladder_status::forced) ++n_res;
            else if (w.status == lsdtt::ladder_status::unresolved) ++n_unres;
            else ++n_deg;
        }
        std::cout << "  " << banded.patches.size() << " patch(es)\n"
                  << "  ladder: " << n_res << " window(s) resolved, "
                  << n_unres << " unresolved, " << n_deg << " degenerate";
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

    if (bf["write_ladder_windows_csv"])
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

    std::cout << "\n";
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-fluvial-surfaces complete.\n";
    return 0;
}
