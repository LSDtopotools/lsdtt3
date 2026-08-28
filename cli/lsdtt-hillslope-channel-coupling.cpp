// lsdtt-hillslope-channel-coupling.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-hillslope-channel-coupling
//
// Couples hillslope form to channel form. For every ridgeline pixel it measures
// hilltop curvature, traces the flow path down to the channel that pixel drains
// to, and reports the hillslope length, relief and the non-dimensional pair
// E* / R*.
//
//   Hurst, M. D., Mudd, S. M., Walcott, R., Attal, M. and Yoo, K. (2012).
//     Using hilltop curvature to derive the spatial distribution of erosion
//     rates. JGR Earth Surface 117:F02017.
//   Grieve, S. W. D., Mudd, S. M., Hurst, M. D. and Milodowski, D. T. (2016).
//     A nondimensional framework for exploring the relief structure of
//     landscapes. Earth Surface Dynamics 4:309-325.
//
// -----------------------------------------------------------------------------
// WHAT A HILLTOP PIXEL IS HERE
//
// Drainage divides are traced on the pixel-corner "divide grid", so a divide
// runs BETWEEN DEM pixels. The hilltop pixels are the pixels a ridgeline
// SEPARATES, one on each side of every edge of the polyline, and because the
// ridgeline is ORDERED each pixel carries its position along the divide and
// which side of it the pixel sits on. That side is the useful part: it
// partitions a divide's hilltops into the two OPPOSING hillslopes, so E*/R* can
// be compared ACROSS a divide rather than only pooled.
//
// -----------------------------------------------------------------------------
// THREE THINGS THE OUTPUT WILL NOT LET YOU FORGET
//
// 1. EVERY hilltop pixel gets a row, including ones that were never traced.
//    The tracer enumerates sources strictly inside the outer ring, so hilltop
//    pixels on the DEM border produce no trace at all; they appear with
//    trace_status = not_traced_border rather than vanishing. On a test fixture
//    that was 35% of the filtered set, which is not a rounding error.
//
// 2. Lh is measured from the source pixel CENTRE to where the path ENTERS the
//    channel cell, not to the channel cell's centre. Short by up to ~0.7 of a
//    cell.
//
// 3. Quantities that are undefined are EMPTY, never a magic number. E*/R* are
//    undefined when Lh == 0 or the curvature cell was nodata; every
//    trace-derived field is undefined when the trace did not reach a channel.
//
// Usage:
//   lsdtt-hillslope-channel-coupling <path> <parameter_file>
//   lsdtt-hillslope-channel-coupling <parameter_file>   (path assumed ./)
//   lsdtt-hillslope-channel-coupling -h                 (help)
//   lsdtt-hillslope-channel-coupling -v                 (version)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <map>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_csv_table.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_basin_cpu.hpp"
#include "../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../src/cpu/lsdtt_chi_network_analysis_cpu.hpp"
#include "../src/cpu/lsdtt_basin_outline_cpu.hpp"
#include "../src/cpu/lsdtt_divide_diagnostics_cpu.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_hilltop_flow_routing_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_raster_binary_mask_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_ridgeline_cpu.hpp"
#include "../src/cpu/lsdtt_streamlines_cpu.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-hillslope-channel-coupling";

namespace
{

std::string lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

void print_banner()
{
    std::cout
        << "\n=========================================================\n"
        << "|| Welcome to the lsdtt3 hillslope-channel coupling tool ||\n"
        << "|| I link hilltop curvature to the channels below it.    ||\n"
        << "|| CLI version: " << k_version << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "=========================================================\n\n";
}

lsdtt::vector_output_format parse_vector_format(const std::string& raw,
                                                const std::string& which)
{
    const std::string v = lower(raw);
    if (v == "flatgeobuf" || v == "fgb") return lsdtt::vector_output_format::flatgeobuf;
    if (v == "gpkg" || v == "geopackage") return lsdtt::vector_output_format::gpkg;
    throw std::invalid_argument(
        which + " must be flatgeobuf, fgb or gpkg (got '" + raw + "')");
}

std::string raster_ext_of(lsdtt::raster_output_format fmt)
{
    return (fmt == lsdtt::raster_output_format::png) ? ".png" : ".tif";
}

lsdtt::raster_output_format parse_raster_format(const std::string& raw)
{
    const std::string v = lower(raw);
    if (v == "cog")   return lsdtt::raster_output_format::cog;
    if (v == "gtiff" || v == "tif" || v == "tiff")
        return lsdtt::raster_output_format::gtiff;
    if (v == "png")   return lsdtt::raster_output_format::png;
    throw std::invalid_argument(
        "write_output_format must be cog, gtiff or png (got '" + raw + "')");
}

// Build a raster of one float per cell, valid everywhere the DEM is valid.
lsdtt::lsdtt_raster_cpu make_scalar_raster(const lsdtt::lsdtt_raster_meta& meta,
                                           std::vector<float> data,
                                           const lsdtt::lsdtt_array2d_u8& valid)
{
    std::vector<std::uint8_t> mask(valid.get_vector());
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

double easting_of(const lsdtt::lsdtt_raster_meta& m, int col)
{
    return m.x_min + (col + 0.5) * m.cell_size;
}
double northing_of(const lsdtt::lsdtt_raster_meta& m, int row)
{
    return m.y_min + (m.n_rows - row - 1 + 0.5) * m.cell_size;
}

// A field value that is EMPTY when the quantity is undefined. Writing a NaN or
// a magic number here is how a sentinel ends up being read as data.
lsdtt::lsdtt_vector_field_value opt_double(double v)
{
    if (std::isnan(v)) return lsdtt::lsdtt_vector_field_value{std::monostate{}};
    return lsdtt::lsdtt_vector_field_value{v};
}

} // namespace

int main(int argc, char* argv[])
{
    print_banner();
    const lsdtt::lsdtt_cli_args cli_args =
        lsdtt::lsdtt_param_parser::parse_args_deferred(
            argc, argv, k_tool_name, k_version, k_citation);

    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry         help_map;

    // ---- IO ------------------------------------------------------------------
    string_default_map["boundary_conditions"] = "nnnn";
    help_map["boundary_conditions"] = {"string", "nnnn",
        "Four-character N,E,S,W boundary code for the flow router.",
        "'n' no flux, 'p' periodic, 'o'/'b' open or base level. MUST match how "
        "the DEM was made: a grid with no draining edge makes the router throw. "
        "Note 'o' and 'b' are inconsistent between the router and the hillslope "
        "solvers, so divides within one hillslope length of such an edge are not "
        "trustworthy.",
        "Input", "param", ""};

    string_default_map["vertical_units"] = "metres";
    help_map["vertical_units"] = {"string","metres",
        "Vertical unit of the input elevations: metres, feet or us_survey_feet.",
        "NO raster format records the vertical unit, so this CANNOT be detected "
        "and is NEVER applied automatically. Much US state lidar is NAVD88 FEET "
        "on a metre grid; left uncorrected, hilltop curvature, relief, E* and R* "
        "are all wrong by 3.28 and NOTHING downstream looks broken - the "
        "landscape merely appears far steeper than it is. Establish the vertical "
        "datum from the data provider before setting this.",
        "Input", "param", "metres|feet|us_survey_feet"};

    float_default_map["minimum_elevation"] = -9999.0f;
    help_map["minimum_elevation"] = {"float","-9999",
        "Pixels at or below this elevation become nodata. -9999 disables.",
        "MANY DEMs CARRY NODATA AS ZERO regardless of what their header "
        "declares, and a block of zeros beside real terrain is a cliff hundreds "
        "of metres high: it corrupts the polyfit, drags every flow path into it, "
        "and inflates hillslope relief and gradient. Check the elevation "
        "histogram before trusting a run. Set to 0 to drop sea-level and "
        "nodata-as-zero pixels.",
        "Input", "param", ""};

    bool_default_map["fill_before_routing"] = true;
    help_map["fill_before_routing"] = {"bool", "true",
        "Pit-fill the DEM before building the flow network.",
        "The router and the channel network need a pit-free surface. The TRACER "
        "still runs on the RAW DEM so paths follow true topography -- both "
        "surfaces are held at once and this switch only affects the router.",
        "Input", "switch", ""};

    // ---- channels ------------------------------------------------------------
    int_default_map["source_pixels_threshold"] = 1000;
    help_map["source_pixels_threshold"] = {"int", "1000",
        "Contributing pixels defining a channel head.",
        "Sets where hillslopes end, so it sets the hillslope length. This is the "
        "single most consequential parameter in the tool.",
        "Channel network", "param", ""};

    string_default_map["channel_network_file"] = "";
    help_map["channel_network_file"] = {"string", "",
        "Optional channel network to ingest instead of extracting one.",
        "A point file written by lsdtt-channel-extraction "
        "(write_channel_network_points), carrying row, col, flat and "
        "stream_order. It MUST have been extracted on this DEM's grid: 'flat' is "
        "row-major on those dimensions and a mismatch is rejected.",
        "Channel network", "param", ""};

    // ---- basin selection -----------------------------------------------------
    string_default_map["basin_selection_choice"] = "sibling";
    help_map["basin_selection_choice"] = {"string", "sibling",
        "How basins are chosen; ridgelines come from their outlines.",
        "sibling: every junction with two or more donors becomes a parent and "
        "its donor sub-basins are the selection. This makes an interfluve "
        "STRUCTURAL -- a divide's two flanking basins share a parent by "
        "construction -- and is the recommended choice for hilltop work. "
        "all_junctions: one basin per junction; network-wide but deeply nested, "
        "so divides are traced many times and must be deduplicated.",
        "Basin selection", "param", ""};

    int_default_map["minimum_basin_pixels"] = 100;
    help_map["minimum_basin_pixels"] = {"int", "100",
        "Skip basins smaller than this.",
        "Tiny basins produce degenerate outlines and no usable hillslopes.",
        "Basin selection", "param", ""};

    // ---- ridge extraction ----------------------------------------------------
    string_default_map["ridge_extraction_algorithm_choice"] = "ridgeline";
    help_map["ridge_extraction_algorithm_choice"] = {"string", "ridgeline",
        "How hilltop pixels are found.",
        "ridgeline: hilltop pixels are the DEM pixels that an ORDERED ridgeline "
        "separates, so each one carries its position along the divide and which "
        "side it is on. divide_pairs: a flow-path merge test over the whole DEM "
        "that needs no basins, but yields no ordering and no divide identity. "
        "Use it as a fallback when basin outlines are slow or fragile.",
        "Ridge extraction", "param", ""};

    string_default_map["diagonal_step_policy"] = "skip";
    help_map["diagonal_step_policy"] = {"string", "skip",
        "What to do with a non-4-connected step in a ridgeline.",
        "A diagonal corner-to-corner step lies along no lattice edge, so it "
        "separates no well-defined pixel pair. skip / bridge / fail. Measured on "
        "real basin outlines: 1808 of 1808 edges were cardinal, so this is a "
        "safety net rather than a routine path.",
        "Ridge extraction", "param", ""};

    int_default_map["divide_merge_search_steps"] = 5;
    help_map["divide_merge_search_steps"] = {"int", "5",
        "Merge tolerance K for the divide_pairs fallback route.",
        "Two adjacent pixels are a divide pair if their flow paths do not meet "
        "within K steps. Ignored unless ridge_extraction_algorithm_choice is "
        "divide_pairs.",
        "Ridge extraction", "param", ""};

    int_default_map["stream_network_padding_pixels"] = 2;
    help_map["stream_network_padding_pixels"] = {"int", "2",
        "Drop hilltop pixels within this many pixels of a channel.",
        "NOT optional in practice. A basin outline includes the segments running "
        "beside the channel down to the outlet, so an unfiltered ridge set "
        "contains pixels sitting ON the channel with no hillslope at all -- "
        "measured at 39 percent on a test fixture. 0 disables; 1 drops only "
        "pixels that are themselves channel pixels.",
        "Ridge extraction", "param", ""};

    float_default_map["threshold_hilltop_gradient"] = 0.4f;
    help_map["threshold_hilltop_gradient"] = {"float", "0.4",
        "Drop hilltop pixels steeper than this. <= 0 disables.",
        "Steep 'ridges' are usually landslide scarps or DEM noise.",
        "Ridge extraction", "param", ""};

    bool_default_map["remove_positive_curvature_hilltops"] = true;
    help_map["remove_positive_curvature_hilltops"] = {"bool", "true",
        "Drop hilltop pixels whose curvature is positive.",
        "A hilltop is convex, and under lsdtt3's sign convention a convex-up "
        "ridge has NEGATIVE curvature. A positive value is therefore not a "
        "hilltop. Note this is the opposite sign to some other packages.",
        "Ridge extraction", "param", ""};

    // ---- surface metrics -----------------------------------------------------
    float_default_map["surface_fitting_window_radius"] = 0.0f;
    help_map["surface_fitting_window_radius"] = {"float", "0",
        "Polynomial surface-fitting window radius in metres. 0 = 3 * cell size.",
        "Cht is SCALE DEPENDENT, so this is a scientific parameter, not a "
        "smoothing detail. It is written into the output summary because an E* "
        "quoted without its window is not reproducible.",
        "Surface metrics", "param", ""};

    // ---- floodplain ----------------------------------------------------------
    string_default_map["floodplain_mask_file"] = "";
    help_map["floodplain_mask_file"] = {"string", "",
        "Optional binary floodplain raster to ingest.",
        "It does TWO things, and both matter: floodplain pixels are added to the "
        "SINK set, so a trace entering the valley floor terminates there, AND "
        "removed from the ridge set, so a floodplain-interior 'hilltop' cannot "
        "produce a degenerate zero-length row.",
        "Floodplain", "param", ""};

    bool_default_map["find_floodplain"] = false;
    help_map["find_floodplain"] = {"bool", "false",
        "Detect a floodplain with a simple gradient + HAND rule.",
        "floodplain = gradient < threshold AND HAND < threshold, hole-filled. "
        "This is a CONVENIENCE, not a classifier: it has no reach-relative "
        "banding and no terrace separation, so it will swallow terraces and "
        "plateaus. Use lsdtt-flat-surface-finder for real work. Default off so "
        "the tool's behaviour never silently depends on two arbitrary numbers. "
        "HAND is height above nearest drainage (Nobre et al. 2011, J. Hydrology "
        "404:13-29).",
        "Floodplain", "switch", ""};

    float_default_map["floodplain_gradient_threshold"] = 0.1f;
    help_map["floodplain_gradient_threshold"] = {"float", "0.1",
        "Gradient below which a pixel may be floodplain.",
        "Used only when find_floodplain is true.",
        "Floodplain", "param", "find_floodplain"};

    float_default_map["floodplain_hand_threshold"] = 5.0f;
    help_map["floodplain_hand_threshold"] = {"float", "5.0",
        "HAND (m) below which a pixel may be floodplain.",
        "Used only when find_floodplain is true.",
        "Floodplain", "param", "find_floodplain"};

    // ---- tracing -------------------------------------------------------------
    string_default_map["trace_surface_choice"] = "raw";
    help_map["trace_surface_choice"] = {"string", "raw",
        "Which surface the flow paths are traced on: raw or fitted.",
        "raw (default): the unmodified DEM, so paths follow true topography. "
        "That is the right choice in principle, but on high-resolution lidar "
        "micro-topography makes the Tarboton edge-routing branch cycle, and a "
        "trace that cycles yields no hillslope length at all -- measured at 48 "
        "percent of traces on a 1 m dataset. "
        "fitted: the polynomial fitted surface at surface_fitting_window_radius, "
        "which is the SAME surface the curvature is measured on, so aspect and "
        "Cht become mutually consistent. It CHANGES WHAT Lh MEANS -- shorter, "
        "smoother paths -- so it is a declared choice and is reported in the "
        "run header, never applied silently.",
        "Tracing", "param", ""};

    int_default_map["hilltop_thinning"] = 1;
    help_map["hilltop_thinning"] = {"int", "1",
        "Analyse every Nth hilltop pixel. 1 = all of them.",
        "This THINS THE ANALYSIS and therefore changes every distribution in the "
        "output table. It is not a rendering option -- for that use "
        "trace_thinning. Leave at 1 unless you know why you are subsampling.",
        "Tracing", "param", ""};

    int_default_map["trace_thinning"] = 20;
    help_map["trace_thinning"] = {"int", "20",
        "Write every Nth trace geometry. Affects the FILE only.",
        "Full-resolution trace geometry on a large DEM is millions of "
        "LineStrings and tens of gigabytes, and no viewer will draw it. The "
        "hilltop table always comes from every analysed pixel regardless of "
        "this setting.",
        "Tracing", "param", "write_hillslope_traces"};

    int_default_map["trace_max_steps"] = 100000;
    help_map["trace_max_steps"] = {"int", "100000",
        "Hard cap on cell steps per trace.",
        "Stops a runaway trace on flat or noisy topography.",
        "Tracing", "param", ""};

    float_default_map["flat_terminus_gradient"] = 0.02f;
    help_map["flat_terminus_gradient"] = {"float", "0.02",
        "Gradient below which a trace that revisits a cell is reported as "
        "stalled_on_flat (5) rather than loop_detected (2).",
        "These are different failures and want different responses. A stall "
        "means the path descended CORRECTLY and then found no channel to "
        "terminate on - a statement about the channel network's coverage, not "
        "about the routing; extend the network, or read Lh as covering only the "
        "part of the DEM the network reaches. A loop means the trace revisited a "
        "cell on ground it could still have descended, which is a genuine "
        "routing failure and makes Lh unreliable. Measured as the steepest "
        "descent to any of the eight neighbours at the terminating cell. The "
        "default is below any hillslope and above the residual gradient a fill "
        "leaves on a flat; it is calibrated on one 1 m landscape, where 74% of "
        "reported loops turned out to be stalls.",
        "Tracing", "param", ""};

    float_default_map["plan_curvature_threshold"] = 0.001f;
    help_map["plan_curvature_threshold"] = {"float", "0.001",
        "abs(plan curvature) above this counts a path cell as curved.",
        "Drives divergent_count and planar_count, which together sum to the "
        "number of visited cells. Note 'divergent' counts CURVED cells of "
        "either sign, because the test is on the magnitude.",
        "Tracing", "param", ""};

    // ---- chi ------------------------------------------------------------------
    bool_default_map["run_chi_analysis"] = false;
    help_map["run_chi_analysis"] = {"bool", "false",
        "Fit the chi network and attach channel steepness to every hilltop row.",
        "This is the COUPLING the tool is named for: hillslope form (Cht, Lh, "
        "E*, R*) alongside the k_sn of the channel that hillslope actually "
        "drains to. k_sn is fitted per SEGMENT, so each channel node is matched "
        "to the segment whose chi interval contains it. Off by default because "
        "the fit costs real time and is not needed for the hillslope metrics "
        "alone.",
        "Chi", "switch", ""};

    float_default_map["m_over_n"] = 0.45f;
    help_map["m_over_n"] = {"float", "0.45",
        "Concavity index m/n used for the chi integration.",
        "Chi is only comparable between basins computed with the SAME m/n, so "
        "this value is written into the summary.",
        "Chi", "param", "run_chi_analysis"};

    float_default_map["chi_A_0"] = 1.0f;
    help_map["chi_A_0"] = {"float", "1.0",
        "Reference drainage area A_0 for the chi integration.",
        "Changes the units of chi and therefore of k_sn.",
        "Chi", "param", "run_chi_analysis"};

    float_default_map["chi_critical_divergence_m"] = 5.0f;
    help_map["chi_critical_divergence_m"] = {"float", "5.0",
        "Segment-fitting tolerance: max divergence (m) before a new segment.",
        "Larger values give fewer, longer segments and therefore smoother k_sn.",
        "Chi", "param", "run_chi_analysis"};

    // ---- non-dimensional -----------------------------------------------------
    float_default_map["critical_hillslope_slope"] = 0.8f;
    help_map["critical_hillslope_slope"] = {"float", "0.8",
        "Critical gradient S_c used to non-dimensionalise E* and R*.",
        "E* = 2|Cht|Lh/S_c and R* = R/(Lh S_c), so both are only as meaningful "
        "as this number and the tool cannot check it for you. 0.8 reproduces the "
        "value hard-coded in the legacy tool, which is the value to use for any "
        "comparison against legacy output.",
        "Non-dimensional", "param", ""};

    // ---- output --------------------------------------------------------------
    string_default_map["hilltop_points_output_format"] = "flatgeobuf";
    help_map["hilltop_points_output_format"] = {"string", "flatgeobuf",
        "Format for the hilltop point output.",
        "flatgeobuf (or fgb), gpkg, or csv. FlatGeobuf is the default and is "
        "recommended: GeoPackage is unreliable on mounted filesystems such as "
        "WSL2's /mnt/c.",
        "Output", "io", ""};

    string_default_map["hillslope_traces_output_format"] = "flatgeobuf";
    help_map["hillslope_traces_output_format"] = {"string", "flatgeobuf",
        "Format for the traced flow paths.",
        "flatgeobuf (or fgb) or gpkg.",
        "Output", "io", "write_hillslope_traces"};

    string_default_map["ridgelines_output_format"] = "flatgeobuf";
    help_map["ridgelines_output_format"] = {"string", "flatgeobuf",
        "Format for the ridgeline polylines.",
        "flatgeobuf (or fgb) or gpkg.",
        "Output", "io", "write_ridgelines"};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string", "cog",
        "Raster output format: cog, gtiff or png.",
        "Shared by every raster this tool writes.",
        "Output", "io", ""};

    bool_default_map["write_hilltop_points"] = true;
    help_map["write_hilltop_points"] = {"bool", "true",
        "Write the hilltop table (one row per hilltop pixel).",
        "The main product.",
        "Output", "switch", ""};

    bool_default_map["write_hillslope_traces"] = false;
    help_map["write_hillslope_traces"] = {"bool", "false",
        "Write the traced flow paths as lines, thinned by trace_thinning.",
        "Off by default because the geometry is large and is a diagnostic "
        "rather than a result.",
        "Output", "switch", ""};

    bool_default_map["write_ridgelines"] = false;
    help_map["write_ridgelines"] = {"bool", "false",
        "Write the ordered ridgelines as lines.",
        "Cheap, and it is what makes ridgeline_id and "
        "distance_along_ridgeline_m on the points usable -- without the lines "
        "those columns reference something the user cannot see.",
        "Output", "switch", ""};

    bool_default_map["write_ridge_mask_raster"] = false;
    help_map["write_ridge_mask_raster"] = {"bool", "false",
        "Write a raster of the hilltop pixels, coded by side.",
        "0 not a hilltop, 1 left of travel, 2 right of travel.",
        "Output", "switch", ""};

    bool_default_map["write_hilltop_curvature_raster"] = false;
    help_map["write_hilltop_curvature_raster"] = {"bool", "false",
        "Write Cht at the hilltop pixels.", "", "Output", "switch", ""};

    bool_default_map["write_hillslope_length_raster"] = false;
    help_map["write_hillslope_length_raster"] = {"bool", "false",
        "Write Lh at the hilltop pixels.", "", "Output", "switch", ""};

    bool_default_map["write_relief_raster"] = false;
    help_map["write_relief_raster"] = {"bool", "false",
        "Write hillslope relief at the hilltop pixels.", "", "Output", "switch", ""};

    bool_default_map["write_hand_raster"] = false;
    help_map["write_hand_raster"] = {"bool", "false",
        "Write the HAND raster used by the floodplain detector.",
        "Only produced when find_floodplain is true.",
        "Output", "switch", "find_floodplain"};

    bool_default_map["write_channel_mask_raster"] = false;
    help_map["write_channel_mask_raster"] = {"bool", "false",
        "Write the channel mask actually used as sinks.",
        "Worth having when the network was INGESTED: it is the only way to see "
        "what was really used, as opposed to what the file was supposed to "
        "contain.",
        "Output", "switch", ""};

    bool_default_map["write_floodplain_mask_raster"] = false;
    help_map["write_floodplain_mask_raster"] = {"bool", "false",
        "Write the floodplain mask actually used.",
        "Whether ingested or detected.",
        "Output", "switch", ""};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int", "0",
        "Thread count. 0 = hardware default.",
        "Measured scaling on this workload saturates near 2x on four physical "
        "cores; hyperthreads add nothing.",
        "Performance", "param", ""};

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

    try
    {
        lsdtt::lsdtt_param_parser pp(cli_args.path, cli_args.fname);
        pp.parse_all_parameters(float_default_map, int_default_map,
                                bool_default_map, string_default_map);
        const auto bf = pp.get_bool_parameters();
        const auto ff = pp.get_float_parameters();
        const auto ip = pp.get_int_parameters();
        const auto sf = pp.get_string_parameters();

        // ---- validate up front, before any expensive work -------------------
        const std::string ridge_choice = lower(sf.at("ridge_extraction_algorithm_choice"));
        if (ridge_choice != "ridgeline" && ridge_choice != "divide_pairs")
            throw std::invalid_argument(
                "ridge_extraction_algorithm_choice must be ridgeline or divide_pairs");

        const std::string basin_choice = lower(sf.at("basin_selection_choice"));
        if (basin_choice != "sibling" && basin_choice != "all_junctions")
            throw std::invalid_argument(
                "basin_selection_choice must be sibling or all_junctions");

        const std::string diag_choice = lower(sf.at("diagonal_step_policy"));
        if (diag_choice != "skip" && diag_choice != "bridge" && diag_choice != "fail")
            throw std::invalid_argument(
                "diagonal_step_policy must be skip, bridge or fail");

        if (ff.at("critical_hillslope_slope") <= 0.0f)
            throw std::invalid_argument(
                "critical_hillslope_slope must be > 0: E* and R* are undefined otherwise");
        if (ip.at("hilltop_thinning") < 1)
            throw std::invalid_argument("hilltop_thinning must be >= 1");
        if (ip.at("trace_thinning") < 1)
            throw std::invalid_argument("trace_thinning must be >= 1");
        if (ip.at("source_pixels_threshold") <= 0)
            throw std::invalid_argument("source_pixels_threshold must be > 0");

        const lsdtt::raster_output_format raster_format =
            parse_raster_format(sf.at("write_output_format"));
        const std::string raster_ext = raster_ext_of(raster_format);

        const std::string points_fmt_raw = lower(sf.at("hilltop_points_output_format"));
        const bool points_as_csv = (points_fmt_raw == "csv");
        lsdtt::vector_output_format points_format =
            lsdtt::vector_output_format::flatgeobuf;
        if (!points_as_csv)
            points_format = parse_vector_format(points_fmt_raw,
                                                "hilltop_points_output_format");

        const lsdtt::vector_output_format traces_format =
            parse_vector_format(sf.at("hillslope_traces_output_format"),
                                "hillslope_traces_output_format");
        const lsdtt::vector_output_format ridge_format =
            parse_vector_format(sf.at("ridgelines_output_format"),
                                "ridgelines_output_format");

        const int n_threads = ip.at("n_threads");

        // ---- load ------------------------------------------------------------
        lsdtt::lsdtt_data_loader loader;
        // read_ext is split off read_fname by the parser; default it so a bare
        // stem works, matching the other CLIs.
        const std::string input_path =
            pp.get_read_path() + pp.get_read_fname() +
            (pp.get_read_ext().empty() ? std::string(".tif") : pp.get_read_ext());
        loader.open_raster(input_path);
        loader.print_raster_info();
        loader.get_meta().check_projected_metres();

        lsdtt::lsdtt_raster_cpu dem_raw(loader);

        // Vertical units BEFORE anything else: every hillslope metric this tool
        // reports scales linearly with elevation, so a foot/metre mix-up is a
        // uniform 3.28x error that leaves no visible symptom. Opt-in only.
        const lsdtt::vertical_unit vunit =
            lsdtt::parse_vertical_unit(sf.at("vertical_units"));
        if (vunit != lsdtt::vertical_unit::metres)
        {
            std::cout << "Converting elevations from " << sf.at("vertical_units")
                      << " to metres.\n";
            dem_raw = lsdtt::convert_vertical_units_to_metres(dem_raw, vunit);
        }
        // No plausibility warning here: it lives in lsdtt-raster-preprocessing,
        // which is where a DEM gets inspected and corrected, and where the 3DEP
        // check can be answered by rerunning one cheap tool rather than this one.

        if (ff.at("minimum_elevation") > -9998.0f)
        {
            std::cout << "Masking pixels at or below "
                      << ff.at("minimum_elevation") << " m to nodata\n";
            dem_raw = lsdtt::mask_to_nodata_using_threshold(
                dem_raw, ff.at("minimum_elevation"), true);
        }
        const lsdtt::lsdtt_raster_meta meta = dem_raw.get_meta();
        const int nr = meta.n_rows, nc = meta.n_cols;
        const double cell = meta.cell_size;
        const std::size_t ncell = static_cast<std::size_t>(nr) * nc;

        // ---- flow routing -----------------------------------------------------
        // The router needs a pit-free surface; the TRACER uses the raw one.
        std::cout << "\n--- Flow routing ---\n";
        const lsdtt::lsdtt_raster_cpu dem_routed =
            bf.at("fill_before_routing") ? lsdtt::fill_pits(dem_raw) : dem_raw;

        lsdtt::lsdtt_flowinfo_options fopts;
        fopts.boundary_conditions     = sf.at("boundary_conditions");
        fopts.build_topological_order = true;
        fopts.build_s_vector          = true;

        std::unique_ptr<lsdtt::lsdtt_flowinfo> flow_owner;
        try
        {
            flow_owner = std::make_unique<lsdtt::lsdtt_flowinfo>(
                dem_routed.get_data_array(), dem_routed.get_mask_array(), meta, fopts);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: flow routing failed with boundary_conditions '"
                      << sf.at("boundary_conditions") << "': " << e.what()
                      << "\n       These must match how the DEM was made; a grid "
                         "with no draining edge cannot be routed.\n";
            return 1;
        }
        const lsdtt::lsdtt_flowinfo& flow = *flow_owner;

        // ---- channel network ---------------------------------------------------
        std::cout << "\n--- Channel network ---\n";
        const int chan_threshold = ip.at("source_pixels_threshold");

        std::vector<float>        chan_data(ncell, 0.0f);
        std::int64_t n_channel = 0;
        std::vector<std::int64_t> sources;

        const std::string chan_file = sf.at("channel_network_file");
        if (!chan_file.empty())
        {
            // INGEST. A channel network extracted elsewhere -- typically by
            // lsdtt-channel-extraction with a Wiener source detector, which
            // finds channel heads from curvature rather than from a drainage
            // area threshold. The point of supporting this is that the channel
            // head placement is the single thing that sets hillslope length,
            // so it should be possible to use a better detector than a
            // threshold.
            std::cout << "Ingesting channel network from " << chan_file << "\n";

            std::vector<std::pair<int,int>> pts;
            const std::string ext = chan_file.size() >= 4
                ? lower(chan_file.substr(chan_file.size() - 4)) : std::string();

            if (ext == ".csv")
            {
                const lsdtt::lsdtt_csv_table t = lsdtt::load_csv_table(chan_file);
                if (!lsdtt::csv_has_column(t, "row") || !lsdtt::csv_has_column(t, "col"))
                    throw std::invalid_argument(
                        "channel_network_file " + chan_file + " has no 'row'/'col' "
                        "columns. Write it with lsdtt-channel-extraction's "
                        "write_channel_network_points.");
                int line = 1;
                for (const auto& r : t.rows)
                {
                    ++line;
                    pts.emplace_back(
                        static_cast<int>(lsdtt::csv_row_double_required(r, "row", chan_file, line)),
                        static_cast<int>(lsdtt::csv_row_double_required(r, "col", chan_file, line)));
                }
            }
            else
            {
                lsdtt::lsdtt_data_loader vl;
                vl.open_vector(chan_file);
                const auto names = vl.get_layer_names();
                if (names.empty())
                    throw std::invalid_argument(
                        "channel_network_file " + chan_file + " has no layers");
                const lsdtt::lsdtt_vector_layer layer = vl.read_vector_layer(names.front());
                const int i_row = layer.schema.index_of("row");
                const int i_col = layer.schema.index_of("col");
                if (i_row < 0 || i_col < 0)
                    throw std::invalid_argument(
                        "channel_network_file " + chan_file + " has no 'row'/'col' "
                        "fields. Write it with lsdtt-channel-extraction's "
                        "write_channel_network_points.");
                for (const auto& rec : layer.records)
                    pts.emplace_back(
                        static_cast<int>(std::get<std::int64_t>(rec.field_values[i_row])),
                        static_cast<int>(std::get<std::int64_t>(rec.field_values[i_col])));
            }

            if (pts.empty())
                throw std::runtime_error(
                    "channel_network_file " + chan_file + " contained no points");

            // The row/col indices are row-major on the grid the network was
            // extracted from. Reject an out-of-range index rather than wrapping
            // it into a wrong pixel -- a network from a different DEM must not
            // be silently accepted.
            std::int64_t n_bad = 0, n_invalid = 0;
            for (const auto& rc : pts)
            {
                if (rc.first < 0 || rc.first >= nr || rc.second < 0 || rc.second >= nc)
                { ++n_bad; continue; }
                const std::size_t k =
                    static_cast<std::size_t>(rc.first) * nc + rc.second;
                if (!flow.is_valid_cell(static_cast<std::int64_t>(k)))
                { ++n_invalid; continue; }
                if (chan_data[k] == 0.0f) { chan_data[k] = 1.0f; ++n_channel; }
            }
            if (n_bad > 0)
                throw std::runtime_error(
                    "channel_network_file has " + std::to_string(n_bad) +
                    " points outside this DEM (" + std::to_string(nr) + "x" +
                    std::to_string(nc) + "). It was extracted on a different "
                    "grid; row/col are row-major on that grid and cannot be "
                    "reused here.");

            if (n_invalid > 0)
                std::cout << "Dropped " << n_invalid << " ingested channel points "
                             "that fall on nodata in this DEM. If that count is "
                             "large the network was extracted BEFORE the DEM was "
                             "masked, and should be re-extracted.\n";
            if (n_channel == 0)
            {
                // Two different causes, and blaming the wrong one sends the user
                // to re-extract a network that was never the problem. If NO
                // point had usable geometry the file is the issue; if points
                // were read but all landed on nodata, the DEM is.
                if (pts.empty())
                    throw std::runtime_error(
                        "channel_network_file: no usable POINT geometry was read "
                        "from '" + chan_file + "'. Multipart points are supported "
                        "and split automatically, so this is either an empty "
                        "layer or a geometry type that is not points - check that "
                        "the file is the *_channel_network_points output and not "
                        "the *_channel_network_lines one.");

                throw std::runtime_error(
                    "channel_network_file: all " + std::to_string(pts.size()) +
                    " ingested channel points fell on nodata in this DEM. The "
                    "network was extracted on a different (or unmasked) DEM - "
                    "re-extract it on the same one.");
            }

            // Sources are the channel cells with no channel donor. Derived
            // rather than ingested, because the points file carries the whole
            // network and not its heads.
            for (int r = 0; r < nr; ++r)
                for (int c = 0; c < nc; ++c)
                {
                    const std::size_t k = static_cast<std::size_t>(r) * nc + c;
                    if (chan_data[k] < 0.5f) continue;
                    const std::int64_t flat = static_cast<std::int64_t>(k);
                    bool has_channel_donor = false;
                    for (int dr = -1; dr <= 1 && !has_channel_donor; ++dr)
                        for (int dc = -1; dc <= 1 && !has_channel_donor; ++dc)
                        {
                            if (dr == 0 && dc == 0) continue;
                            const int rr = r + dr, cc = c + dc;
                            if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
                            const std::int64_t nf =
                                static_cast<std::int64_t>(rr) * nc + cc;
                            if (chan_data[static_cast<std::size_t>(nf)] < 0.5f) continue;
                            if (!flow.is_valid_cell(nf)) continue;
                            if (flow.receiver_of(nf) == flat) has_channel_donor = true;
                        }
                    if (!has_channel_donor) sources.push_back(flat);
                }

            std::cout << "Ingested " << pts.size() << " channel points -> "
                      << n_channel << " channel pixels, " << sources.size()
                      << " derived heads\n";
        }
        else
        {
            const lsdtt::lsdtt_raster_result cp =
                flow.get_n_contributing_pixels_raster();
            for (std::size_t i = 0; i < ncell; ++i)
                if (cp.mask[i] && cp.data[i] >= static_cast<float>(chan_threshold))
                { chan_data[i] = 1.0f; ++n_channel; }

            if (n_channel == 0)
                throw std::runtime_error(
                    "no channel pixels at source_pixels_threshold "
                    + std::to_string(chan_threshold) + " - lower it");

            sources = flow.get_sources_index_threshold(chan_threshold);
        }

        if (sources.empty())
            throw std::runtime_error("no channel sources");

        lsdtt::lsdtt_junction_network_options jopts;
        lsdtt::lsdtt_junction_network_cpu jn(flow, sources, jopts);
        std::cout << "Channel pixels: " << n_channel
                  << "   sources: " << sources.size()
                  << "   junctions: " << jn.n_junctions() << "\n";

        // ---- surface metrics ---------------------------------------------------
        std::cout << "\n--- Surface metrics ---\n";
        const double window = ff.at("surface_fitting_window_radius") > 0.0f
            ? static_cast<double>(ff.at("surface_fitting_window_radius"))
            : 3.0 * cell;
        std::cout << "Polynomial window radius: " << window << " m\n";

        // -------------------------------------------------------------------
        // Cht IS THE LAPLACIAN. This is not a choice, and it used to be one.
        //
        // E* = 2|Cht|Lh/S_c is defined with the Laplacian in Roering et al.
        // (2007) and Hurst et al. (2012). lsdtt3's `mean_curvature` is the
        // DIFFERENTIAL-GEOMETRY mean curvature,
        //
        //     H = ((1+q^2)r - 2pqs + (1+p^2)t) / (2 (1+p^2+q^2)^(3/2))
        //
        // which is a different quantity: it reduces to grad^2(z)/2 only where
        // p = q = 0, and at finite gradient the metric terms make the ratio
        // slope-dependent. Feeding it to E* therefore biases the result by
        // about a factor of two AND biases gentle and steep hilltops
        // differently.
        //
        // Measured on a synthetic landscape with KNOWN forcing (U = 4e-4 m/yr,
        // D = 0.05 m^2/yr, so the divide curvature must be U/D = 0.0080 1/m):
        // mean_curvature gave 0.0035 (0.43x) and the Laplacian gave 0.0071
        // (0.89x, the remainder being a hillslope solve short of steady state).
        // Nothing about the wrong answer looked wrong -- R* stayed below 1, the
        // cloud was tight, no warning was printed. It simply sat in the wrong
        // place on the curve.
        //
        // Offering the choice invited exactly that, so it is gone.
        // -------------------------------------------------------------------
        constexpr lsdtt::polyfit_metric k_hilltop_curvature =
            lsdtt::polyfit_metric::laplacian;

        // The parameter parser does not report unrecognised keys, so a
        // parameter file written against the old CLI would have its setting
        // dropped without a word -- and silently getting a DIFFERENT answer
        // than the file asks for is the failure mode this change exists to
        // remove. Say so.
        if (pp.parameter_supplied("curvature_metric_choice"))
            std::cout << "\nNOTE: curvature_metric_choice is no longer a "
                         "parameter and has been ignored.\n"
                         "      Hilltop curvature is always the LAPLACIAN, "
                         "which is what E* is defined with.\n"
                         "      mean_curvature is the differential-geometry "
                         "mean curvature, a different\n      quantity: it "
                         "biases E* by roughly a factor of two, and by a "
                         "slope-dependent\n      amount, without anything in "
                         "the output looking wrong.\n\n";

        const std::string trace_surface = lower(sf.at("trace_surface_choice"));
        if (trace_surface != "raw" && trace_surface != "fitted")
            throw std::invalid_argument(
                "trace_surface_choice must be raw or fitted");

        std::unique_ptr<lsdtt::lsdtt_raster_cpu> curv_r, grad_r, plan_r, fitted_r;
        {
            const auto coeffs = lsdtt::direct_polyfit(dem_raw, window, n_threads);
            lsdtt::polyfit_metric_request req;
            req.metrics = std::vector<lsdtt::polyfit_metric>{
                k_hilltop_curvature,
                lsdtt::polyfit_metric::slope,
                lsdtt::polyfit_metric::plan_curvature,
                lsdtt::polyfit_metric::fitted_surface};
            req.n_threads = n_threads;
            const auto m = lsdtt::compute_polyfit_metrics(coeffs, req);

            curv_r = std::make_unique<lsdtt::lsdtt_raster_cpu>(
                m.rasters.at(k_hilltop_curvature));
            grad_r = std::make_unique<lsdtt::lsdtt_raster_cpu>(
                m.rasters.at(lsdtt::polyfit_metric::slope));
            plan_r = std::make_unique<lsdtt::lsdtt_raster_cpu>(
                m.rasters.at(lsdtt::polyfit_metric::plan_curvature));
            fitted_r = std::make_unique<lsdtt::lsdtt_raster_cpu>(
                m.rasters.at(lsdtt::polyfit_metric::fitted_surface));
        }

        const lsdtt::lsdtt_raster_cpu& trace_dem =
            (trace_surface == "fitted") ? *fitted_r : dem_raw;
        std::cout << "Tracing surface: " << trace_surface
                  << (trace_surface == "fitted"
                        ? "  (polyfit fitted surface -- Lh is smoother and shorter "
                          "than on the raw DEM)"
                        : "  (raw DEM)")
                  << "\n";

        // ---- floodplain ---------------------------------------------------------
        std::vector<float> floodplain(ncell, 0.0f);
        bool have_floodplain = false;

        if (!sf.at("floodplain_mask_file").empty())
        {
            lsdtt::lsdtt_data_loader fp_loader;
            fp_loader.open_raster(sf.at("floodplain_mask_file"));
            const lsdtt::lsdtt_raster_cpu fp(fp_loader);
            if (fp.get_meta().n_rows != nr || fp.get_meta().n_cols != nc)
                throw std::invalid_argument(
                    "floodplain_mask_file does not match the DEM grid");
            for (int r = 0; r < nr; ++r)
                for (int c = 0; c < nc; ++c)
                    if (fp.get_mask_array().at(r, c) &&
                        fp.get_data_array().at(r, c) > 0.5f)
                        floodplain[static_cast<std::size_t>(r) * nc + c] = 1.0f;
            have_floodplain = true;
            std::cout << "\n--- Floodplain ---\nIngested " << sf.at("floodplain_mask_file") << "\n";
        }
        else if (bf.at("find_floodplain"))
        {
            std::cout << "\n--- Floodplain (simple detector) ---\n"
                      << "This is a convenience rule, not a classifier: it will "
                         "include terraces and plateaus.\n";
            std::vector<std::uint8_t> chan_u8(ncell, 0);
            for (std::size_t i = 0; i < ncell; ++i)
                chan_u8[i] = chan_data[i] > 0.5f ? 1 : 0;

            const lsdtt::lsdtt_raster_result hand =
                flow.get_height_above_channel_mask_raster(
                    dem_routed.get_data_array().get_vector(),
                    dem_routed.get_mask_array().get_vector(),
                    chan_u8);

            const float g_thr = ff.at("floodplain_gradient_threshold");
            const float h_thr = ff.at("floodplain_hand_threshold");
            std::vector<float> raw_fp(ncell, 0.0f);
            for (int r = 0; r < nr; ++r)
                for (int c = 0; c < nc; ++c)
                {
                    const std::size_t k = static_cast<std::size_t>(r) * nc + c;
                    if (!hand.mask[k]) continue;
                    if (!grad_r->get_mask_array().at(r, c)) continue;
                    if (grad_r->get_data_array().at(r, c) < g_thr &&
                        hand.data[k] < h_thr)
                        raw_fp[k] = 1.0f;
                }

            lsdtt::lsdtt_raster_cpu fp_raster(
                std::move(raw_fp),
                std::vector<std::uint8_t>(dem_raw.get_mask_array().get_vector()),
                meta);
            const lsdtt::lsdtt_raster_cpu filled = lsdtt::binary_fill_holes(fp_raster);
            floodplain = filled.get_data_array().get_vector();
            have_floodplain = true;

            if (bf.at("write_hand_raster"))
                loader.write_full(hand.data, meta,
                                  pp.make_output_path("_HAND", raster_ext),
                                  raster_format);
        }

        std::int64_t n_floodplain = 0;
        if (have_floodplain)
        {
            for (std::size_t i = 0; i < ncell; ++i)
                if (floodplain[i] > 0.5f) ++n_floodplain;
            std::cout << "Floodplain pixels: " << n_floodplain
                      << " (added to sinks, removed from ridges)\n";
        }

        // Sinks = channels plus floodplain.
        std::vector<float> sink_data = chan_data;
        for (std::size_t i = 0; i < ncell; ++i)
            if (floodplain[i] > 0.5f) sink_data[i] = 1.0f;

        const lsdtt::lsdtt_raster_cpu channel_raster = make_scalar_raster(
            meta, chan_data, dem_raw.get_mask_array());
        const lsdtt::lsdtt_raster_cpu sink_raster = make_scalar_raster(
            meta, sink_data, dem_raw.get_mask_array());

        // ---- ridge extraction ---------------------------------------------------
        std::cout << "\n--- Ridge extraction ---\n";
        std::vector<lsdtt::lsdtt_ridgeline_cpu> ridgelines;
        std::vector<std::int64_t>               basin_flats;
        lsdtt::lsdtt_ridge_extraction_result    ridge;

        if (ridge_choice == "ridgeline")
        {
            std::vector<std::int64_t> outlets;
            if (basin_choice == "sibling")
            {
                // Sibling sets: a divide's two flanking basins share a parent by
                // construction, so an interfluve is structural rather than
                // accidental.
                for (std::int32_t j = 0; j < jn.n_junctions(); ++j)
                {
                    std::vector<std::int32_t> donors;
                    jn.for_each_donor_junction(
                        j, [&](std::int32_t d) { donors.push_back(d); });
                    if (donors.size() < 2) continue;
                    for (const std::int32_t d : donors)
                    {
                        try { outlets.push_back(jn.basin_outlet_flat_of_junction(d)); }
                        catch (const std::exception&) {}
                    }
                }
            }
            else
            {
                for (std::int32_t j = 0; j < jn.n_junctions(); ++j)
                    try { outlets.push_back(jn.basin_outlet_flat_of_junction(j)); }
                    catch (const std::exception&) {}
            }

            std::sort(outlets.begin(), outlets.end());
            outlets.erase(std::unique(outlets.begin(), outlets.end()), outlets.end());

            int n_basins = 0;
            for (const std::int64_t outlet : outlets)
            {
                if (outlet < 0) continue;
                try
                {
                    lsdtt::lsdtt_basin_cpu basin(outlet, flow);
                    if (basin.get_number_of_basin_pixels() < ip.at("minimum_basin_pixels"))
                        continue;
                    lsdtt::lsdtt_basin_outline_cpu outline(basin);
                    if (outline.outline_divide_rc().size() < 4) continue;
                    for (auto& rl : outline.split_outline_into_ridgelines_by_equal_length())
                    {
                        if (rl.size() < 2) continue;
                        ridgelines.push_back(std::move(rl));
                        basin_flats.push_back(outlet);
                    }
                    ++n_basins;
                }
                catch (const std::exception&) { continue; }
            }

            if (ridgelines.empty())
                throw std::runtime_error(
                    "no ridgelines extracted - try basin_selection_choice: "
                    "all_junctions, or lower minimum_basin_pixels");

            lsdtt::lsdtt_ridge_extraction_options ropts;
            ropts.diagonal_policy =
                diag_choice == "bridge" ? lsdtt::diagonal_step_policy::bridge
              : diag_choice == "fail"   ? lsdtt::diagonal_step_policy::fail
                                        : lsdtt::diagonal_step_policy::skip;

            ridge = lsdtt::extract_hilltop_pixels(ridgelines, meta, ropts, basin_flats);

            std::cout << "Basins: " << n_basins
                      << "   ridgelines: " << ridgelines.size()
                      << "   edges: " << ridge.n_edges_total
                      << " (cardinal " << ridge.n_edges_cardinal
                      << ", diagonal " << ridge.n_edges_diagonal << ")\n";
        }
        else
        {
            lsdtt::lsdtt_divide_diagnostic_options dopts;
            dopts.divide_merge_search_steps = ip.at("divide_merge_search_steps");
            dopts.n_threads = n_threads;
            const auto pairs = lsdtt::compute_divide_pairs(flow, dem_routed, dopts);

            std::vector<std::int64_t> seen;
            seen.reserve(pairs.size() * 2);
            for (const auto& pr : pairs) { seen.push_back(pr.flat_a); seen.push_back(pr.flat_b); }
            std::sort(seen.begin(), seen.end());
            seen.erase(std::unique(seen.begin(), seen.end()), seen.end());

            for (const std::int64_t f : seen)
            {
                lsdtt::lsdtt_hilltop_pixel px;
                px.flat = f;
                px.row  = static_cast<int>(f / nc);
                px.col  = static_cast<int>(f % nc);
                ridge.pixels.push_back(px);
            }
            std::cout << "divide_pairs route: " << pairs.size() << " pairs -> "
                      << ridge.pixels.size() << " hilltop pixels "
                      << "(no ordering or side; those columns will be empty)\n";
        }

        std::cout << "Hilltop pixels before filtering: " << ridge.pixels.size()
                  << "   (duplicates dropped " << ridge.n_duplicates_dropped
                  << ", off-grid " << ridge.n_off_grid_dropped << ")\n";

        // ---- filters -------------------------------------------------------------
        lsdtt::lsdtt_hilltop_filter_options filt;
        filt.stream_network_padding_pixels      = ip.at("stream_network_padding_pixels");
        filt.threshold_hilltop_gradient         = ff.at("threshold_hilltop_gradient");
        filt.remove_positive_curvature_hilltops = bf.at("remove_positive_curvature_hilltops");

        // Floodplain pixels are removed from the ridge set by folding them into
        // the mask the padding filter tests against.
        lsdtt::lsdtt_raster_cpu filter_mask = have_floodplain ? sink_raster : channel_raster;

        lsdtt::lsdtt_hilltop_filter_report frep;
        std::vector<lsdtt::lsdtt_hilltop_pixel> filtered =
            lsdtt::filter_hilltop_pixels(ridge.pixels, meta, filt,
                                         &filter_mask, grad_r.get(), curv_r.get(), frep);

        std::cout << "Filters: on a channel " << frep.n_on_channel_before
                  << ", dropped near channel " << frep.n_dropped_near_channel
                  << ", dropped steep " << frep.n_dropped_steep
                  << ", dropped convex-up " << frep.n_dropped_convex_up
                  << "  ->  " << frep.n_out << " hilltop pixels\n";

        if (ip.at("hilltop_thinning") > 1)
        {
            std::vector<lsdtt::lsdtt_hilltop_pixel> thinned;
            const int step = ip.at("hilltop_thinning");
            for (std::size_t i = 0; i < filtered.size(); i += static_cast<std::size_t>(step))
                thinned.push_back(filtered[i]);
            std::cout << "hilltop_thinning=" << step << " -> " << thinned.size()
                      << " analysed (THIS CHANGES EVERY DISTRIBUTION BELOW)\n";
            filtered = std::move(thinned);
        }

        if (filtered.empty())
            throw std::runtime_error(
                "every hilltop pixel was filtered out - relax "
                "stream_network_padding_pixels, threshold_hilltop_gradient or "
                "remove_positive_curvature_hilltops");

        // ---- metrics --------------------------------------------------------------
        std::cout << "\n--- Hillslope metrics ---\n";
        lsdtt::lsdtt_hillslope_metrics_options mopts;
        mopts.critical_slope_Sc        = ff.at("critical_hillslope_slope");
        mopts.plan_curvature_threshold = ff.at("plan_curvature_threshold");
        mopts.max_steps                = ip.at("trace_max_steps");
        mopts.flat_terminus_gradient   = ff.at("flat_terminus_gradient");
        mopts.n_threads                = n_threads;

        auto metrics = lsdtt::compute_hillslope_metrics(
            filtered, trace_dem, sink_raster,
            curv_r.get(), grad_r.get(), plan_r.get(), mopts);

        std::cout << "Rows: " << metrics.rows.size()
                  << "   reached channel " << metrics.n_reached_channel
                  << ", hit edge " << metrics.n_hit_edge
                  << ", looped " << metrics.n_loop_detected
                  << ", stalled on flat " << metrics.n_stalled_on_flat
                  << ", no aspect " << metrics.n_no_aspect
                  << ", never traced (DEM border) " << metrics.n_not_traced_border
                  << "\n";
        if (metrics.n_not_traced_border > 0)
            std::cout << "Note: " << metrics.n_not_traced_border
                      << " hilltop pixels sit on the outer ring and are never "
                         "offered to the tracer. They are in the table with "
                         "trace_status = 4 and empty metrics.\n";
        // These two say different things and get different advice. A stall is a
        // statement about where the CHANNEL NETWORK reaches; a loop is a
        // statement about the routing.
        if (metrics.n_stalled_on_flat > metrics.n_reached_channel / 4 &&
            metrics.n_reached_channel > 0)
            std::cout << "NOTE: " << metrics.n_stalled_on_flat
                      << " traces (trace_status = 5) descended correctly and "
                         "then ran out of\n      descent with no channel to "
                         "terminate on. That is the channel network's\n"
                         "      coverage, not a routing failure: those hilltops "
                         "drain to valleys the\n      extraction did not put a "
                         "channel in. Extend the network, or read the Lh\n"
                         "      distribution as covering only the part of the "
                         "DEM the network reaches.\n";
        if (metrics.n_loop_detected >
            metrics.n_reached_channel / 2 && metrics.n_reached_channel > 0)
            std::cout << "WARNING: a large fraction of traces cycled ON "
                         "DESCENDING GROUND (trace_status = 2).\n         That "
                         "is a routing failure, and Lh distributions from this "
                         "run are not reliable.\n";

        // ---- chi join (Phase 5) -------------------------------------------------
        // The channel a hilltop drains to is already known from the trace; this
        // attaches that channel's steepness. Done IN PROCESS: a post-hoc CSV
        // join on floating-point coordinates drops rows silently, and a dropped
        // row is indistinguishable from a hilltop that genuinely has no channel.
        if (bf.at("run_chi_analysis"))
        {
            std::cout << "\n--- Chi ---\n";
            lsdtt::lsdtt_chi_options chi_opts;
            chi_opts.m_over_n = ff.at("m_over_n");
            chi_opts.A_0      = ff.at("chi_A_0");

            // Seed from the basin outlets of the junction network, so chi is
            // computed on the same channels the traces terminate on.
            std::vector<std::int64_t> outlet_flats;
            for (std::int32_t j = 0; j < jn.n_junctions(); ++j)
                try { outlet_flats.push_back(jn.basin_outlet_flat_of_junction(j)); }
                catch (const std::exception&) {}
            std::sort(outlet_flats.begin(), outlet_flats.end());
            outlet_flats.erase(std::unique(outlet_flats.begin(), outlet_flats.end()),
                               outlet_flats.end());

            // The chi extractor REJECTS an outlet below its contributing-pixel
            // threshold, and that throw would abort a whole run for one small
            // basin. With an INGESTED network the junction outlets have nothing
            // to do with source_pixels_threshold, so filter here rather than
            // letting the extractor decide.
            const lsdtt::lsdtt_raster_result cp_chi =
                flow.get_n_contributing_pixels_raster();
            const std::int64_t n_before = static_cast<std::int64_t>(outlet_flats.size());
            std::vector<std::int64_t> kept;
            kept.reserve(outlet_flats.size());
            for (const std::int64_t f : outlet_flats)
            {
                const std::size_t k = static_cast<std::size_t>(f);
                if (k < cp_chi.data.size() && cp_chi.mask[k] &&
                    cp_chi.data[k] >= static_cast<float>(chan_threshold))
                    kept.push_back(f);
            }
            outlet_flats.swap(kept);
            if (n_before != static_cast<std::int64_t>(outlet_flats.size()))
                std::cout << "Dropped " << (n_before - static_cast<std::int64_t>(outlet_flats.size()))
                          << " of " << n_before << " basin outlets below the chi "
                             "contributing-pixel threshold (" << chan_threshold << ")\n";

            if (outlet_flats.empty())
            {
                std::cout << "No basin outlet meets the chi threshold; skipping the "
                             "chi join. Lower source_pixels_threshold to enable it.\n";
                goto chi_done;
            }

            {
            const lsdtt::lsdtt_chi_collection_result chi_networks =
                lsdtt::get_chi_from_outlets_zero_base_threshold(
                    flow, outlet_flats, chan_threshold, chi_opts);

            lsdtt::lsdtt_chi_network_fit_options fit_opts;
            fit_opts.critical_divergence_m = ff.at("chi_critical_divergence_m");

            const auto fits = lsdtt::fit_chi_networks(chi_networks, dem_routed, fit_opts);
            const auto chi_nodes    = lsdtt::export_chi_network_fit_nodes(fits);
            const auto chi_segments = lsdtt::export_chi_network_fit_segments(fits);

            std::vector<float> flow_dist;
            std::vector<float> drainage_area;
            {
                const lsdtt::lsdtt_raster_result fd =
                    flow.get_distance_from_outlet_raster();
                flow_dist = fd.data;
                const lsdtt::lsdtt_raster_result da = flow.get_drainage_area_raster();
                drainage_area = da.data;
            }

            const auto chan_props = lsdtt::build_channel_property_map(
                chi_nodes, chi_segments, &flow_dist, &drainage_area);

            const std::int64_t n_joined =
                lsdtt::attach_channel_properties(metrics, chan_props, nc);

            std::cout << "m/n = " << chi_opts.m_over_n
                      << "   fitted channel nodes " << chi_nodes.flat.size()
                      << ", segments " << chi_segments.segment_id.size() << "\n"
                      << "Joined k_sn onto " << n_joined << " of "
                      << metrics.n_reached_channel
                      << " hilltop rows that reached a channel";
            if (metrics.n_reached_channel > 0)
                std::cout << "  ("
                          << (100.0 * static_cast<double>(n_joined)
                              / static_cast<double>(metrics.n_reached_channel))
                          << " %)";
            std::cout << "\n";
            if (n_joined < metrics.n_reached_channel)
                std::cout << "Note: rows whose channel pixel was not part of a "
                             "fitted network keep an empty channel_k_sn. They are "
                             "not zero and are not dropped.\n";
            }
            chi_done: ;
        }

        // ---- outputs ---------------------------------------------------------------
        std::cout << "\n--- Output ---\n";

        if (bf.at("write_hilltop_points"))
        {
            lsdtt::lsdtt_vector_schema schema;
            const auto add = [&](const std::string& n, lsdtt::lsdtt_vector_field_type t)
            { schema.field_names.push_back(n); schema.field_types.push_back(t); };

            using FT = lsdtt::lsdtt_vector_field_type;
            add("row", FT::int64_type);
            add("col", FT::int64_type);
            add("Cht", FT::double_type);
            add("hilltop_gradient", FT::double_type);
            add("mean_path_gradient", FT::double_type);
            add("relief_m", FT::double_type);
            add("hillslope_length_m", FT::double_type);
            add("euclidean_distance_m", FT::double_type);
            add("channel_row", FT::int64_type);
            add("channel_col", FT::int64_type);
            add("E_star", FT::double_type);
            add("R_star", FT::double_type);
            add("ridgeline_id", FT::int64_type);
            add("ridgeline_node_index", FT::int64_type);
            add("distance_along_ridgeline_m", FT::double_type);
            add("ridge_side", FT::int64_type);
            add("basin_outlet_flat", FT::int64_type);
            add("divergent_count", FT::int64_type);
            add("planar_count", FT::int64_type);
            add("trace_status", FT::int64_type);
            add("channel_k_sn", FT::double_type);
            add("channel_chi", FT::double_type);
            add("channel_elevation", FT::double_type);
            add("channel_flow_distance_m", FT::double_type);
            add("channel_drainage_area_m2", FT::double_type);
            add("channel_basin_index", FT::int64_type);
            add("channel_segment_id", FT::int64_type);

            const std::string ext = points_as_csv
                ? ".csv" : lsdtt::output_extension(points_format);
            const std::string path = pp.make_output_path("_hilltop_points", ext);

            // CSV goes through write_csv, NOT the vector session: handing a
            // .csv path to the OGR vector writer makes a DIRECTORY of that name
            // with a .fgb inside it.
            std::vector<std::string> headers;
            std::vector<std::vector<std::string>> csv_rows;
            const auto num = [](double v) {
                if (std::isnan(v)) return std::string();
                std::ostringstream os; os << std::setprecision(10) << v;
                return os.str();
            };

            if (points_as_csv)
            {
                headers = {"latitude_proxy_northing", "easting", "row", "col",
                           "Cht", "hilltop_gradient", "mean_path_gradient",
                           "relief_m", "hillslope_length_m",
                           "euclidean_distance_m", "channel_row", "channel_col",
                           "E_star", "R_star", "ridgeline_id",
                           "ridgeline_node_index", "distance_along_ridgeline_m",
                           "ridge_side", "basin_outlet_flat", "divergent_count",
                           "planar_count", "trace_status", "channel_k_sn",
                           "channel_chi", "channel_elevation",
                           "channel_flow_distance_m", "channel_drainage_area_m2",
                           "channel_basin_index", "channel_segment_id"};
                csv_rows.reserve(metrics.rows.size());
                for (const auto& r : metrics.rows)
                    csv_rows.push_back({
                        num(northing_of(meta, r.ridge.row)),
                        num(easting_of(meta, r.ridge.col)),
                        std::to_string(r.ridge.row), std::to_string(r.ridge.col),
                        num(r.hilltop_curvature), num(r.hilltop_gradient),
                        num(r.mean_path_gradient), num(r.relief_m),
                        num(r.hillslope_length_m), num(r.euclidean_distance_m),
                        std::to_string(r.channel_row), std::to_string(r.channel_col),
                        num(r.E_star), num(r.R_star),
                        std::to_string(r.ridge.ridgeline_id),
                        std::to_string(r.ridge.ridgeline_node_index),
                        num(r.ridge.distance_along_ridgeline_m),
                        std::to_string(static_cast<int>(r.ridge.ridge_side)),
                        std::to_string(r.ridge.basin_outlet_flat),
                        std::to_string(r.divergent_count),
                        std::to_string(r.planar_count),
                        std::to_string(static_cast<int>(r.status)),
                        num(r.channel_k_sn), num(r.channel_chi),
                        num(r.channel_elevation), num(r.channel_flow_distance_m),
                        num(r.channel_drainage_area_m2),
                        std::to_string(r.channel_basin_index),
                        std::to_string(r.channel_segment_id)});
                loader.write_csv(path, headers, csv_rows);
                std::cout << "Wrote " << path << "  (" << metrics.rows.size()
                          << " rows)\n";
            }
            else
            {
            loader.begin_vector_write(path, "hilltop_points", schema,
                                      meta.crs_wkt, points_format);
            for (const auto& r : metrics.rows)
            {
                std::vector<lsdtt::lsdtt_vector_field_value> v;
                v.emplace_back(static_cast<std::int64_t>(r.ridge.row));
                v.emplace_back(static_cast<std::int64_t>(r.ridge.col));
                v.push_back(opt_double(r.hilltop_curvature));
                v.push_back(opt_double(r.hilltop_gradient));
                v.push_back(opt_double(r.mean_path_gradient));
                v.push_back(opt_double(r.relief_m));
                v.push_back(opt_double(r.hillslope_length_m));
                v.push_back(opt_double(r.euclidean_distance_m));
                v.emplace_back(static_cast<std::int64_t>(r.channel_row));
                v.emplace_back(static_cast<std::int64_t>(r.channel_col));
                v.push_back(opt_double(r.E_star));
                v.push_back(opt_double(r.R_star));
                v.emplace_back(static_cast<std::int64_t>(r.ridge.ridgeline_id));
                v.emplace_back(static_cast<std::int64_t>(r.ridge.ridgeline_node_index));
                v.emplace_back(r.ridge.distance_along_ridgeline_m);
                v.emplace_back(static_cast<std::int64_t>(r.ridge.ridge_side));
                v.emplace_back(r.ridge.basin_outlet_flat);
                v.emplace_back(static_cast<std::int64_t>(r.divergent_count));
                v.emplace_back(static_cast<std::int64_t>(r.planar_count));
                v.emplace_back(static_cast<std::int64_t>(r.status));
                v.push_back(opt_double(r.channel_k_sn));
                v.push_back(opt_double(r.channel_chi));
                v.push_back(opt_double(r.channel_elevation));
                v.push_back(opt_double(r.channel_flow_distance_m));
                v.push_back(opt_double(r.channel_drainage_area_m2));
                v.emplace_back(static_cast<std::int64_t>(r.channel_basin_index));
                v.emplace_back(static_cast<std::int64_t>(r.channel_segment_id));
                loader.write_point(easting_of(meta, r.ridge.col),
                                   northing_of(meta, r.ridge.row), v);
            }
            loader.finish_vector_write();
            std::cout << "Wrote " << path << "  (" << metrics.rows.size()
                      << " rows)\n";
            }
        }

        if (bf.at("write_ridgelines") && !ridgelines.empty())
        {
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {"ridgeline_id", "n_nodes", "total_length_m",
                                  "basin_outlet_flat"};
            schema.field_types = {lsdtt::lsdtt_vector_field_type::int64_type,
                                  lsdtt::lsdtt_vector_field_type::int64_type,
                                  lsdtt::lsdtt_vector_field_type::double_type,
                                  lsdtt::lsdtt_vector_field_type::int64_type};

            const std::string path = pp.make_output_path(
                "_ridgelines", lsdtt::output_extension(ridge_format));
            loader.begin_vector_write(path, "ridgelines", schema,
                                      meta.crs_wkt, ridge_format);
            for (std::size_t i = 0; i < ridgelines.size(); ++i)
            {
                const auto& rl = ridgelines[i];
                if (rl.size() < 2) continue;
                std::vector<lsdtt::lsdtt_xy> pts(rl.positions().begin(),
                                                 rl.positions().end());
                std::vector<lsdtt::lsdtt_vector_field_value> v;
                v.emplace_back(static_cast<std::int64_t>(i));
                v.emplace_back(static_cast<std::int64_t>(rl.size()));
                v.emplace_back(rl.total_length());
                v.emplace_back(i < basin_flats.size() ? basin_flats[i]
                                                      : static_cast<std::int64_t>(-1));
                loader.write_line(pts, v);
            }
            loader.finish_vector_write();
            std::cout << "Wrote " << path << "  (" << ridgelines.size() << " ridgelines)\n";
        }

        if (bf.at("write_hillslope_traces"))
        {
            // Re-trace only the thinned subset WITH geometry: tracing everything
            // with save_nodes on would multiply peak memory for output nobody
            // can draw.
            const int thin = ip.at("trace_thinning");
            std::vector<float> src(ncell, 0.0f);
            std::int64_t n_written_src = 0;
            for (std::size_t i = 0; i < filtered.size(); i += static_cast<std::size_t>(thin))
            { src[static_cast<std::size_t>(filtered[i].flat)] = 1.0f; ++n_written_src; }

            const lsdtt::lsdtt_raster_cpu src_raster =
                make_scalar_raster(meta, src, dem_raw.get_mask_array());

            lsdtt::lsdtt_streamline_options topts;
            topts.save_nodes = true;
            topts.max_steps  = ip.at("trace_max_steps");
            topts.n_threads  = n_threads;

            const auto traces =
                lsdtt::trace_streamlines(trace_dem, src_raster, sink_raster, topts);

            const std::string path = pp.make_output_path(
                "_hillslope_traces", lsdtt::output_extension(traces_format));
            lsdtt::write_streamlines(traces, meta, loader, path,
                                     "hillslope_traces", traces_format, true);
            // Report what was WRITTEN, not what was requested: border sources
            // and failed traces produce no geometry, so the two differ.
            std::cout << "Wrote " << path << "  (" << n_written_src
                      << " sources offered, thinned by " << thin << "; "
                      << traces.streamlines.size() << " traced, "
                      << traces.n_reached_sink << " reached a channel)\n";
        }

        // ---- rasters ---------------------------------------------------------------
        const auto write_hilltop_raster =
            [&](bool enabled, const std::string& suffix,
                const std::function<double(const lsdtt::lsdtt_hilltop_metrics_row&)>& get)
        {
            if (!enabled) return;
            std::vector<float> v(ncell, static_cast<float>(meta.no_data_value));
            for (const auto& r : metrics.rows)
            {
                const double d = get(r);
                if (!std::isnan(d))
                    v[static_cast<std::size_t>(r.ridge.flat)] = static_cast<float>(d);
            }
            const std::string path = pp.make_output_path(suffix, raster_ext);
            loader.write_full(v, meta, path, raster_format);
            std::cout << "Wrote " << path << "\n";
        };

        write_hilltop_raster(bf.at("write_ridge_mask_raster"), "_ridge_mask",
            [](const lsdtt::lsdtt_hilltop_metrics_row& r)
            { return static_cast<double>(r.ridge.ridge_side + 1); });
        write_hilltop_raster(bf.at("write_hilltop_curvature_raster"), "_Cht",
            [](const lsdtt::lsdtt_hilltop_metrics_row& r) { return r.hilltop_curvature; });
        write_hilltop_raster(bf.at("write_hillslope_length_raster"), "_Lh",
            [](const lsdtt::lsdtt_hilltop_metrics_row& r) { return r.hillslope_length_m; });
        write_hilltop_raster(bf.at("write_relief_raster"), "_hillslope_relief",
            [](const lsdtt::lsdtt_hilltop_metrics_row& r) { return r.relief_m; });

        if (bf.at("write_channel_mask_raster"))
        {
            const std::string path = pp.make_output_path("_channel_mask", raster_ext);
            loader.write_full(chan_data, meta, path, raster_format);
            std::cout << "Wrote " << path << "\n";
        }

        if (bf.at("write_floodplain_mask_raster") && have_floodplain)
        {
            const std::string path = pp.make_output_path("_floodplain_mask", raster_ext);
            loader.write_full(floodplain, meta, path, raster_format);
            std::cout << "Wrote " << path << "\n";
        }

        pp.print_help(help_map, pp.get_write_fname() + "_README",
                      k_version, k_citation);
        pp.print_parameters();

        std::cout << "\n" << k_tool_name << ": done.\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }
}
