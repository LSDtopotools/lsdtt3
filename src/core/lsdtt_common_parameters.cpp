// lsdtt_common_parameters.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Shared parameter-group declarations. See the header for why this file
// exists, the prefix rule, the overwrite policy and the provenance of every
// name and default declared here.

#include "lsdtt_common_parameters.hpp"

#include <algorithm>
#include <cmath>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace lsdtt
{

namespace
{

// Group label for the help entry. Prefixed groups get their own section so a
// CLI declaring the same group twice does not collapse into one heading.
std::string prefixed_group_label(const std::string& prefix,
                                 const std::string& group)
{
    if (prefix.empty())
    {
        return group;
    }
    std::string bare = prefix;
    while (!bare.empty() && bare.back() == '_')
    {
        bare.pop_back();
    }
    if (bare.empty())
    {
        return group;
    }
    return group + " (" + bare + ")";
}

// Shortest text that reads back to the same float, for the default column of
// a help entry ("0.05", not "0.050000").
std::string round_trip_text(float value)
{
    for (int precision = 1; precision <= 9; ++precision)
    {
        std::ostringstream os;
        os.imbue(std::locale::classic());
        os.precision(precision);
        os << value;
        std::istringstream is(os.str());
        is.imbue(std::locale::classic());
        float back = 0.0f;
        is >> back;
        if (back == value) return os.str();
    }
    std::ostringstream os;
    os.imbue(std::locale::classic());
    os.precision(9);
    os << value;
    return os.str();
}

} // namespace

// -----------------------------------------------------------------------------
std::string prefixed_parameter_key(const std::string& prefix,
                                   const std::string& key)
{
    if (prefix.empty())
    {
        return key;
    }
    if (prefix.back() == '_')
    {
        return prefix + key;
    }
    return prefix + "_" + key;
}

// -----------------------------------------------------------------------------
// Flow routing core
// -----------------------------------------------------------------------------
void declare_routing_core_parameters(
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Flow routing");
    const std::string k_boundary_conditions =
        prefixed_parameter_key(prefix, "boundary_conditions");
    string_default_map[k_boundary_conditions] = "oooo";
    help_map[k_boundary_conditions] = {"string", "oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        group, "param", ""};

    const std::string k_n_threads = prefixed_parameter_key(prefix, "n_threads");
    int_default_map[k_n_threads] = 0;
    help_map[k_n_threads] = {"int", "0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. A positive value is honoured "
        "verbatim and OVERRIDES OMP_NUM_THREADS, so setting 1 here runs serial "
        "on any machine.",
        group, "io", ""};
}

// -----------------------------------------------------------------------------
// Mandatory-fill assertion
// -----------------------------------------------------------------------------
void declare_raster_is_filled_parameter(
    std::map<std::string, bool>& bool_default_map,
    lsdtt_help_registry&         help_map,
    const std::string&           prefix)
{
    const std::string group = prefixed_group_label(prefix, "Pit filling");
    const std::string key = prefixed_parameter_key(prefix, "raster_is_filled");
    bool_default_map[key] = false;
    help_map[key] = {"bool", "false",
        "Assert that the input raster has ALREADY been pit-filled.",
        "ASKS: has the input already been filled? Filling is REQUIRED for this "
        "analysis, so this switch does not choose whether to fill -- it only "
        "lets the tool skip filling a surface that is filled already. Setting "
        "it true for a DEM that still contains pits gives BAD routing, not "
        "different routing. Do not confuse this with fill_before_routing, "
        "which asks whether to fill at all.",
        group, "switch", ""};
}

// -----------------------------------------------------------------------------
// Optional-fill analytical choice
// -----------------------------------------------------------------------------
void declare_fill_before_routing_parameter(
    std::map<std::string, bool>& bool_default_map,
    lsdtt_help_registry&         help_map,
    const std::string&           prefix)
{
    const std::string group = prefixed_group_label(prefix, "Pit filling");
    const std::string key = prefixed_parameter_key(prefix, "fill_before_routing");
    bool_default_map[key] = true;
    help_map[key] = {"bool", "true",
        "Choose whether to pit-fill the DEM before building the flow network.",
        "ASKS: should the DEM be filled at all? Filling is OPTIONAL here and "
        "this is a genuine analytical choice. True (the default) routes over a "
        "pit-free surface. False builds the flow network on the UNFILLED "
        "surface, so flow paths terminate at sinks instead of being carried "
        "through them. Do not confuse this with raster_is_filled, which is an "
        "assertion that the input is already filled.",
        group, "switch", ""};
}

// -----------------------------------------------------------------------------
void declare_interior_nodata_parameters(
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Interior nodata");
    const std::string method_key = prefixed_parameter_key(prefix, "interior_nodata_method");
    string_default_map[method_key] = "sink";
    help_map[method_key] = {"string", "sink",
        "How nodata pixels inside the DEM are treated: sink or fill.",
        "Interior nodata are nodata pixels not connected to the raster edge "
        "through other nodata (isolated voids, moulins, masked lakes). sink (the "
        "default, unchanged behaviour): they are outlets, so flow that reaches "
        "them ends there, as at the DEM edge; use this for real internal "
        "drainage (e.g. moulins on an ice sheet). fill: they are interpolated "
        "with nodata_fill_method_choice before any filling or routing, so flow "
        "passes across them; use this for data voids. Nodata connected to the "
        "edge (e.g. sea removed by min_elevation) is never filled and always "
        "acts as an open boundary.",
        group, "param", ""};

    const std::string choice_key = prefixed_parameter_key(prefix, "nodata_fill_method_choice");
    string_default_map[choice_key] = "laplacian";
    help_map[choice_key] = {"string", "laplacian",
        "Method for filling interior nodata pixels.",
        "laplacian: smooth Dirichlet interpolation (recommended). "
        "iterative_weighted_average: fast, well-suited to thin voids. Used when "
        + method_key + " is fill.",
        group, "param", method_key + "=fill"};

    const std::string iter_key = prefixed_parameter_key(prefix, "nodata_fill_max_iter");
    int_default_map[iter_key] = 10000;
    help_map[iter_key] = {"int", "10000",
        "Maximum iterations for iterative nodata fill methods.",
        "Used when " + method_key + " is fill.",
        group, "param", method_key + "=fill"};
}

// -----------------------------------------------------------------------------
void declare_fill_min_slope_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    const std::string&            prefix)
{
    const std::string group = prefixed_group_label(prefix, "Pit filling");
    const std::string key = prefixed_parameter_key(prefix, "fill_min_slope");
    float_default_map[key] = 0.0001f;
    help_map[key] = {"float", "0.0001",
        "Minimum slope (m/m) imposed by pit filling or flat resolution.",
        "inline_gradient (default): imposed during the Priority-Flood fill "
        "across filled depressions and flats, so used only when "
        "raster_is_filled is false. barnes2014b_gradient: the minimum gradient "
        "across flats, also on a raster_is_filled input. "
        "barnes2014b_float_step: ignored.",
        group, "param", ""};
}

// -----------------------------------------------------------------------------
void declare_fill_flat_resolution_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix,
    const std::string&                  min_slope_key)
{
    const std::string group = prefixed_group_label(prefix, "Pit filling");
    const std::string key =
        prefixed_parameter_key(prefix, "fill_flat_resolution_choice");
    string_default_map[key] = "inline_gradient";
    help_map[key] = {"string", "inline_gradient",
        "How flats are made drainable: inline_gradient, barnes2014b_gradient or "
        "barnes2014b_float_step.",
        "inline_gradient (default, unchanged): the Priority-Flood fill imposes "
        + min_slope_key + " during the flood. The barnes2014b options make flow "
        "converge across flats (Barnes et al. 2014, "
        "doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses "
        + min_slope_key + " (m/m, must be > 0) as the minimum gradient; flats are "
        "raised about 3-4x more, capped below their rims, and capped flats "
        "are reported. barnes2014b_float_step is equivalent to RichDEM "
        "ResolveFlats on float32 in its first pass: it ignores " + min_slope_key +
        " and its slopes depend on elevation (of order ULP(z) / cell size). With "
        "either barnes2014b option the result is always routable: cells the first "
        "pass leaves without a lower neighbour are refilled and resolved again (up "
        "to 4 rounds), and any still left are raised by a final priority-flood pass "
        "by one float32 step each; the report gives the rounds and the cells "
        "changed after the first pass.",
        group, "param", ""};
}

// -----------------------------------------------------------------------------
void declare_channel_initiation_parameters(
    std::map<std::string, int>& int_default_map,
    lsdtt_help_registry&        help_map,
    const std::string&          prefix)
{
    const std::string group = prefixed_group_label(prefix, "Channel network");
    const std::string key =
        prefixed_parameter_key(prefix, "threshold_contributing_pixels");
    int_default_map[key] = 1000;
    help_map[key] = {"int", "1000",
        "Contributing-pixel threshold defining a channel head.",
        "The total number of upslope contributing pixels from the flow-routing "
        "object, not the number of direct donors to a pixel. Reduce this value "
        "to extract a denser network. This is a PIXEL COUNT and therefore "
        "depends on the grid resolution: the equivalent drainage area is "
        "threshold_contributing_pixels x cell area. Where a genuine area is "
        "wanted instead, the convention is <what>_area_threshold_m2.",
        group, "param", ""};
}

// -----------------------------------------------------------------------------
// Flat-ground admission and floodplain bounds
// -----------------------------------------------------------------------------
void declare_relief_window_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    int                           default_cells,
    const std::string&            parent)
{
    const std::string cells = std::to_string(default_cells);
    float_default_map["relief_window_radius_metres"] = 0.0f;
    help_map["relief_window_radius_metres"] = {"float", "0",
        "Radius (m) of the square window local relief (max - min) is measured "
        "over. 0 = " + cells + " cells.",
        "Converted to the nearest whole number of cells, at least one; the run "
        "log states the radius actually used. A radius r cells is a (2r+1) x "
        "(2r+1) window, so the default is a " +
        std::to_string(2 * default_cells + 1) + " x " +
        std::to_string(2 * default_cells + 1) + " window. Replaces "
        "relief_window_side_pixels (a full side) and "
        "fallback_relief_window_radius_pixels (a radius in cells), which "
        "spelled the same idea two ways.",
        "Admission", "param", parent};
}

void declare_max_relief_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    float                         default_metres,
    const std::string&            parent)
{
    float_default_map["max_relief_metres"] = default_metres;
    help_map["max_relief_metres"] = {"float", round_trip_text(default_metres),
        "Local relief (m) below which a pixel is admitted as flat.",
        "Read a value off the relief CDF the tool writes. Where the tool has "
        "relief_threshold_algorithm_choice, this is the value used by `fixed`.",
        "Admission", "param", parent};
}

void declare_relief_threshold_algorithm_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  default_choice)
{
    string_default_map["relief_threshold_algorithm_choice"] = default_choice;
    help_map["relief_threshold_algorithm_choice"] = {"string", default_choice,
        "How the relief bound is chosen: fixed, cdf_divergence or qq_gaussian.",
        "fixed uses max_relief_metres. cdf_divergence computes the relief CDF, "
        "fits the best of a family of parametric distributions (normal, "
        "exponential, Weibull, Pareto) and takes the relief at which the "
        "empirical CDF most exceeds the fit: flat pixels inflate the empirical "
        "distribution at low relief. qq_gaussian is the Clubb et al. (2017) "
        "rule: the lowest relief at which the quantiles rejoin a Gaussian drawn "
        "through two percentiles of the QQ plot (relief_qq_*) and stay on it. "
        "The table and the chosen value are written out. Replaces the boolean "
        "adaptive_relief_threshold.",
        "Admission", "param", ""};
}

void declare_max_admission_gradient_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    float                         default_gradient,
    const std::string&            parent)
{
    float_default_map["max_admission_gradient"] = default_gradient;
    help_map["max_admission_gradient"] = {"float",
        round_trip_text(default_gradient),
        "Gradient (m/m) below which a pixel is admitted as flat.",
        "A gradient, not an angle: 0.087 is 5 degrees, 0.176 is 10, 0.268 is "
        "15 (gradient = tan(angle)). Measured from a fitted polynomial surface "
        "rather than a 3x3 difference.",
        "Admission", "param", parent};
}

void declare_max_floodplain_hand_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    float                         default_metres,
    const std::string&            parent)
{
    float_default_map["max_floodplain_hand_metres"] = default_metres;
    help_map["max_floodplain_hand_metres"] = {"float",
        round_trip_text(default_metres),
        "HAND (m) above which ground is not floodplain.",
        "HAND is height above nearest drainage (Nobre et al. 2011, J. "
        "Hydrology 404:13-29), measured along the flow path to the channel.",
        "Floodplain", "param", parent};
}

void declare_max_floodplain_gradient_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    float                         default_gradient,
    const std::string&            parent)
{
    float_default_map["max_floodplain_gradient"] = default_gradient;
    help_map["max_floodplain_gradient"] = {"float",
        round_trip_text(default_gradient),
        "Gradient (m/m) above which ground is not floodplain.",
        "A gradient, not an angle (gradient = tan(angle)).",
        "Floodplain", "param", parent};
}

void declare_floodplain_threshold_algorithm_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  default_choice)
{
    string_default_map["floodplain_threshold_algorithm_choice"] = default_choice;
    help_map["floodplain_threshold_algorithm_choice"] = {"string", default_choice,
        "How the floodplain bounds are chosen: fixed, qq_gaussian or "
        "cdf_divergence.",
        "fixed uses max_floodplain_hand_metres and max_floodplain_gradient as "
        "given. qq_gaussian is Clubb et al. (2017): for each quantity, the "
        "lowest value at which the data rejoin a Gaussian drawn through two "
        "percentiles of the quantile-quantile plot (floodplain_qq_*), and stay "
        "on it. cdf_divergence fits a family of distributions and takes the "
        "value at which the empirical CDF most exceeds the best fit. The chosen "
        "values are printed and written to the ingested-parameters file, so a "
        "run can be repeated with fixed.",
        "Floodplain", "param", ""};
}

void declare_qq_threshold_parameters(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    const std::string&            prefix,
    const std::string&            parent)
{
    const std::string lo = prefix + "_qq_lower_percentile";
    const std::string hi = prefix + "_qq_upper_percentile";
    const std::string fr = prefix + "_qq_threshold_fraction";
    float_default_map[lo] = 25.0f;
    help_map[lo] = {"float", "25",
        "Lower percentile of the QQ cloud the Gaussian reference line is drawn "
        "through (qq_gaussian).",
        "25 and 75 (the IQR) are the values Clubb et al. (2017) used.",
        "Threshold selection", "param", parent};
    float_default_map[hi] = 75.0f;
    help_map[hi] = {"float", "75",
        "Upper percentile of the QQ cloud the Gaussian reference line is drawn "
        "through (qq_gaussian).",
        "Must be above the lower percentile.",
        "Threshold selection", "param", parent};
    float_default_map[fr] = 0.01f;
    help_map[fr] = {"float", "0.01",
        "Fraction of the data range within which a quantile counts as on the "
        "Gaussian line (qq_gaussian).",
        "0.01 (1%) is the paper's value. The threshold is the lowest quantile "
        "that is on the line and stays on it for the next 200 quantiles of "
        "10000. Larger values admit more ground.",
        "Threshold selection", "param", parent};
}

threshold_selection_request threshold_request_from_parameters(
    const std::string&                  choice,
    const std::string&                  choice_key,
    double                              fixed_value,
    const std::map<std::string, float>& float_parameters,
    const std::string&                  prefix)
{
    threshold_selection_request r;
    if (choice == "fixed")               r.method = threshold_selection_choice::fixed;
    else if (choice == "qq_gaussian")    r.method = threshold_selection_choice::qq_gaussian;
    else if (choice == "cdf_divergence") r.method = threshold_selection_choice::cdf_divergence;
    else
        throw std::invalid_argument(
            "unknown " + choice_key + " '" + choice +
            "'. Valid: fixed, qq_gaussian, cdf_divergence.");
    r.fixed_value = fixed_value;
    auto get = [&](const std::string& key, double fallback)
    {
        const auto it = float_parameters.find(key);
        return it == float_parameters.end() ? fallback
                                            : static_cast<double>(it->second);
    };
    r.qq_lower_percentile   = get(prefix + "_qq_lower_percentile", 25.0);
    r.qq_upper_percentile   = get(prefix + "_qq_upper_percentile", 75.0);
    r.qq_threshold_fraction = get(prefix + "_qq_threshold_fraction", 0.01);
    return r;
}

int relief_window_radius_pixels(double radius_metres,
                                double cell_size,
                                int    default_cells)
{
    if (!std::isfinite(radius_metres) || radius_metres < 0.0)
        throw std::invalid_argument(
            "relief_window_radius_metres must be a finite value >= 0");
    if (radius_metres == 0.0) return default_cells;
    if (!(cell_size > 0.0))
        throw std::invalid_argument("relief window: cell size must be > 0");
    const long r = std::lround(radius_metres / cell_size);
    return static_cast<int>(std::max(1L, r));
}

double gradient_from_legacy_degrees(const std::string& text,
                                    const std::string& old_key)
{
    std::size_t used = 0;
    double deg = 0.0;
    try { deg = std::stod(text, &used); }
    catch (const std::exception&) { used = 0; }
    if (used == 0 || used != text.size() || !std::isfinite(deg) ||
        deg < 0.0 || deg >= 90.0)
        throw std::invalid_argument(
            "'" + old_key + "' must be an angle in [0, 90) degrees to convert "
            "to a gradient; got '" + text + "'.");
    if (deg == 0.0) return 0.0;
    return std::tan(deg * 3.14159265358979323846 / 180.0);
}

std::string relief_threshold_choice_from_legacy_bool(const std::string& text,
                                                     const std::string& old_key)
{
    for (const char* t : {"true", "True", "TRUE", "t", "T", "1"})
        if (text == t) return "cdf_divergence";
    for (const char* f : {"false", "False", "FALSE", "f", "F", "0"})
        if (text == f) return "fixed";
    throw std::invalid_argument(
        "'" + old_key + "' must be true or false; got '" + text + "'.");
}

// -----------------------------------------------------------------------------
// Output formats
// -----------------------------------------------------------------------------
void declare_raster_output_format_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Output");
    const std::string k_write_output_format =
        prefixed_parameter_key(prefix, "write_output_format");
    string_default_map[k_write_output_format] = "cog";
    help_map[k_write_output_format] = {"string", "cog",
        "Raster output format: cog, gtiff, png or same_as_input.",
        "cog is a Cloud Optimised GeoTIFF and is the project default. It falls "
        "back to gtiff for small rasters, where the tiling overhead is not "
        "worth it. There is ONE raster format parameter per CLI: it governs "
        "every raster the tool writes. Vector outputs are governed separately "
        "by write_vector_format.",
        group, "io", ""};
}

void declare_vector_output_format_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Output");
    const std::string k_write_vector_format =
        prefixed_parameter_key(prefix, "write_vector_format");
    string_default_map[k_write_vector_format] = "flatgeobuf";
    help_map[k_write_vector_format] = {"string", "flatgeobuf",
        "Vector output format: flatgeobuf or gpkg.",
        "Governs every vector layer this tool writes. FlatGeobuf is preferred: "
        "it is compact, carries its own georeferencing, and does not have "
        "GeoPackage's problems writing across drives (a common WSL2 "
        "situation). CSV is available only through selectors for point/table "
        "outputs. Per-output overrides exist only where a container "
        "genuinely differs. Raster outputs are governed separately by "
        "write_output_format.",
        group, "io", ""};
}

void declare_spatial_index_parameter(
    std::map<std::string, bool>& bool_default_map,
    lsdtt_help_registry&         help_map,
    const std::string&           prefix)
{
    const std::string group = prefixed_group_label(prefix, "Output");
    const std::string k_calculate_spatial_index =
        prefixed_parameter_key(prefix, "calculate_spatial_index");
    bool_default_map[k_calculate_spatial_index] = true;
    help_map[k_calculate_spatial_index] = {"bool", "true",
        "Whether to calculate a spatial index for vector data writes.",
        "A spatial index helps a GIS speed up spatial queries, so it is often "
        "helpful to have one. But calculating it requires the entire vector "
        "layer to sit in memory, preventing feature-by-feature writes and "
        "adding a lot of overhead. For large datasets it may therefore be "
        "preferable to write without an index; one can always be built later "
        "with e.g. gdal.",
        group, "io", ""};
}

std::string resolve_output_format_override(
    const std::string& override_value,
    const std::string& canonical_value)
{
    if (override_value == inherit_write_vector_format ||
        override_value == inherit_write_output_format)
        return canonical_value;
    return override_value;
}

// -----------------------------------------------------------------------------
// Basin outlet input and snapping
// -----------------------------------------------------------------------------
void declare_basin_outlet_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    (void)bool_default_map;
    (void)float_default_map;

    const std::string group = prefixed_group_label(prefix, "Basin outlets");

    const std::string k_basin_outlet_fname =
        prefixed_parameter_key(prefix, "basin_outlet_fname");
    string_default_map[k_basin_outlet_fname] = "NULL";
    help_map[k_basin_outlet_fname] = {"string", "NULL",
        "File of requested basin outlet points (CSV, FlatGeobuf or GeoPackage).",
        "Format is auto-detected from the file extension: *.csv is read as a "
        "text file with latitude and longitude columns (WGS84 / EPSG:4326 "
        "decimal degrees); *.fgb or *.gpkg is read as a point vector layer in "
        "its own CRS, which is then transformed to the raster CRS. For CSV "
        "inputs the file must contain columns named exactly 'latitude' and "
        "'longitude'. For vector inputs all features must be point geometry.",
        group, "param", ""};

    const std::string k_basin_outlet_snapping_method =
        prefixed_parameter_key(prefix, "basin_outlet_snapping_method");
    string_default_map[k_basin_outlet_snapping_method] = "downstream";
    help_map[k_basin_outlet_snapping_method] = {"string", "downstream",
        "How outlet points are snapped onto the channel network.",
        "Allowed values: none, expanding_square, downstream, largest_area. "
        "none: use the mapped raster pixel directly. "
        "expanding_square: search outward ring-by-ring until a pixel reaching "
        "basin_outlet_contributing_pixels_threshold is found. "
        "downstream: follow the D8 receiver chain until a qualifying pixel is "
        "found. "
        "largest_area: take the pixel with the most contributing pixels in the "
        "square of half-width basin_outlet_snap_max_search_radius_pixels, if it "
        "reaches basin_outlet_contributing_pixels_threshold (a large radius can "
        "jump onto a bigger neighbouring river).",
        group, "param", ""};

    const std::string k_basin_outlet_contributing_pixels_threshold =
        prefixed_parameter_key(prefix,
                               "basin_outlet_contributing_pixels_threshold");
    int_default_map[k_basin_outlet_contributing_pixels_threshold] =
        default_basin_outlet_contributing_pixels_threshold;
    help_map[k_basin_outlet_contributing_pixels_threshold] = {"int",
        std::to_string(default_basin_outlet_contributing_pixels_threshold),
        "Contributing-pixel threshold a pixel must reach (>=) to accept an outlet.",
        "This is the total number of upslope contributing pixels from the "
        "flow-routing object, not the number of direct donors to a pixel. Used "
        "only when basin_outlet_snapping_method is not none. It is a separate "
        "quantity from threshold_contributing_pixels, which sets where the "
        "channel network begins.",
        group, "param", ""};

    const std::string k_basin_outlet_snap_min_search_radius_pixels =
        prefixed_parameter_key(prefix, "basin_outlet_snap_min_search_radius_pixels");
    int_default_map[k_basin_outlet_snap_min_search_radius_pixels] = 0;
    help_map[k_basin_outlet_snap_min_search_radius_pixels] = {"int", "0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "0 means the mapped input pixel is tested first.",
        group, "param", ""};

    const std::string k_basin_outlet_snap_max_search_radius_pixels =
        prefixed_parameter_key(prefix, "basin_outlet_snap_max_search_radius_pixels");
    int_default_map[k_basin_outlet_snap_max_search_radius_pixels] = 25;
    help_map[k_basin_outlet_snap_max_search_radius_pixels] = {"int", "25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square or "
        "largest_area. For expanding_square the search stops at the first "
        "radius containing one or more qualifying pixels; for largest_area it "
        "is the half-width of the square searched.",
        group, "param", ""};
}

// -----------------------------------------------------------------------------
// Polynomial surface fitting
// -----------------------------------------------------------------------------
void declare_surface_fitting_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    (void)bool_default_map;
    (void)int_default_map;

    const std::string group =
        prefixed_group_label(prefix, "Polynomial surface fit");

    const std::string k_polyfit_lambda_metres =
        prefixed_parameter_key(prefix, "polyfit_lambda_metres");
    float_default_map[k_polyfit_lambda_metres] = 70.0f;
    help_map[k_polyfit_lambda_metres] = {"float", "70.0",
        "Characteristic lengthscale (metres) of the quadratic polynomial fit.",
        "This is a HALF-WAVELENGTH, not a radius: it is the scale of the "
        "topographic feature the fit is tuned to resolve, and it sets the size "
        "of the moving window. The default (70 m) is tuned for Copernicus 30 m "
        "data; reduce it for lidar (e.g. 10-30 m). Parameters that really are "
        "radii are spelled <stage>_<quantity>_radius_metres.",
        group, "param", ""};

    const std::string k_polyfit_metrics =
        prefixed_parameter_key(prefix, "polyfit_metrics");
    string_default_map[k_polyfit_metrics] = "";
    help_map[k_polyfit_metrics] = {"string", "",
        "Comma-separated list of derived polyfit metrics to write.",
        "Examples: slope,aspect_degrees,plan_curvature,profile_curvature. "
        "Use all to write every supported metric. "
        "Available metrics: fitted_surface, dzdx, dzdy, d2zdx2, d2zdxdy, "
        "d2zdy2, slope, slope_degrees, aspect_radians, aspect_degrees, "
        "laplacian, hessian_determinant, mean_curvature, gaussian_curvature, "
        "maximum_curvature, minimum_curvature, curvedness, shape_index, "
        "profile_curvature, plan_curvature, tangential_curvature, "
        "difference_curvature, horizontal_excess_curvature, "
        "vertical_excess_curvature, accumulation_curvature, ring_curvature, "
        "rotor.",
        group, "param", ""};

    const std::string k_compute_backend =
        prefixed_parameter_key(prefix, "compute_backend");
    string_default_map[k_compute_backend] = "auto";
    help_map[k_compute_backend] = {"string", "auto",
        "Preferred compute backend for kernels that have a GPU port.",
        "auto: use the GPU when a device and GPU kernels are both available, "
        "otherwise CPU (default). cpu: force CPU. gpu: prefer the GPU, falling "
        "back to CPU per-function when no GPU kernel exists or no device is "
        "present. Currently only direct_polyfit has a GPU port; all other "
        "metrics run on the CPU regardless of this setting. Has no effect in a "
        "build without GPU support (-DLSDTT_GPU=OFF).",
        group, "param", ""};
}

// -----------------------------------------------------------------------------
// Synthetic model domain
// -----------------------------------------------------------------------------
void declare_synthetic_dimensions_parameters(
    std::map<std::string, int>&         int_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Synthetic domain");

    const std::string k_n_rows = prefixed_parameter_key(prefix, "n_rows");
    int_default_map[k_n_rows] = 128;
    help_map[k_n_rows] = {"int", "128",
        "Number of rows in the synthetic model domain.",
        "Used only when the domain is generated rather than read from a DEM.",
        group, "param", ""};

    const std::string k_n_cols = prefixed_parameter_key(prefix, "n_cols");
    int_default_map[k_n_cols] = 128;
    help_map[k_n_cols] = {"int", "128",
        "Number of columns in the synthetic model domain.",
        "Used only when the domain is generated rather than read from a DEM.",
        group, "param", ""};
}

void declare_synthetic_cell_size_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    const std::string&            prefix)
{
    const std::string group = prefixed_group_label(prefix, "Synthetic domain");
    const std::string k_cell_size_metres =
        prefixed_parameter_key(prefix, "cell_size_metres");
    float_default_map[k_cell_size_metres] = 1.0f;
    help_map[k_cell_size_metres] = {"float", "1.0",
        "Cell size of the synthetic model domain in metres.",
        "Sets the map scale of the generated surface. Anything expressed as a "
        "pixel count (for example threshold_contributing_pixels) corresponds to "
        "a different ground area as this changes.",
        group, "param", ""};
}

void declare_synthetic_random_seed_parameter(
    std::map<std::string, int>& int_default_map,
    lsdtt_help_registry&        help_map,
    const std::string&          prefix)
{
    const std::string group = prefixed_group_label(prefix, "Synthetic domain");
    const std::string k_random_seed =
        prefixed_parameter_key(prefix, "random_seed");
    int_default_map[k_random_seed] = 33;
    help_map[k_random_seed] = {"int", "33",
        "Seed for every stochastic generator used to build the domain.",
        "A fixed seed makes a run reproducible. Change it to draw a different "
        "realisation of the same statistical surface.",
        group, "param", ""};
}

void declare_synthetic_crs_wkt_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Synthetic domain");
    const std::string k_synthetic_crs_wkt =
        prefixed_parameter_key(prefix, "synthetic_crs_wkt");
    string_default_map[k_synthetic_crs_wkt] = "";
    help_map[k_synthetic_crs_wkt] = {"string", "",
        "Coordinate Reference System WKT for the synthetic surface.",
        "A generated surface has no natural location, but downstream tools "
        "expect a projected CRS in metres. If left empty (the default) the "
        "project's default synthetic georeference is assigned and the "
        "lower-left corner is placed at a fixed easting and northing.",
        group, "param", ""};
}

void declare_synthetic_domain_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    (void)bool_default_map;
    declare_synthetic_dimensions_parameters(int_default_map, help_map, prefix);
    declare_synthetic_cell_size_parameter(float_default_map, help_map, prefix);
    declare_synthetic_random_seed_parameter(int_default_map, help_map, prefix);
    declare_synthetic_crs_wkt_parameter(string_default_map, help_map, prefix);
}

// -----------------------------------------------------------------------------
// Synthetic noise surfaces
// -----------------------------------------------------------------------------
void declare_synthetic_noise_type_parameter(
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Synthetic noise");

    const std::string k_noise_type =
        prefixed_parameter_key(prefix, "noise_type");
    string_default_map[k_noise_type] = "flat";
    help_map[k_noise_type] = {"string", "flat",
        "Type of noise surface to generate: flat, white_noise, red_noise, "
        "diamond_square.",
        "flat = all zeros, i.e. no noise; white_noise = uniform random, equal "
        "power at every wavelength; red_noise = low-frequency dominated; "
        "diamond_square = fractal terrain whose roughness is set by "
        "diamond_square_roughness.",
        group, "param", ""};
}

void declare_synthetic_noise_target_stddev_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    const std::string&            prefix)
{
    const std::string group = prefixed_group_label(prefix, "Synthetic noise");
    const std::string k_noise_target_stddev_metres =
        prefixed_parameter_key(prefix, "noise_target_stddev_metres");
    float_default_map[k_noise_target_stddev_metres] = 1.0f;
    help_map[k_noise_target_stddev_metres] = {"float", "1.0",
        "Target standard deviation of the generated noise surface (metres).",
        "The raw surface is rescaled to this standard deviation. Used when "
        "noise_type is red_noise or diamond_square.",
        group, "param", ""};
}

void declare_synthetic_noise_relief_parameter(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    const std::string&            prefix)
{
    const std::string group = prefixed_group_label(prefix, "Synthetic noise");
    const std::string k_noise_relief_metres =
        prefixed_parameter_key(prefix, "noise_relief_metres");
    float_default_map[k_noise_relief_metres] = 0.0f;
    help_map[k_noise_relief_metres] = {"float", "0.0",
        "Amplitude (metres) of white noise added on top of the surface. "
        "0 disables it.",
        "Roughens whatever closed depressions survive, so that when they are "
        "filled the result still has a gradient to route on rather than a "
        "flat. Applied everywhere, not only inside depressions. Size it "
        "against the mean fill depth rather than against the total relief.",
        group, "param", ""};
}

void declare_diamond_square_parameters(
    std::map<std::string, float>& float_default_map,
    lsdtt_help_registry&          help_map,
    const std::string&            prefix)
{
    const std::string group = prefixed_group_label(prefix, "Synthetic noise");
    const std::string k_diamond_square_roughness =
        prefixed_parameter_key(prefix, "diamond_square_roughness");
    float_default_map[k_diamond_square_roughness] = 0.55f;
    help_map[k_diamond_square_roughness] = {"float", "0.55",
        "Roughness decay per level for the diamond-square surface (0-1).",
        "0.5 gives Brownian noise; lower values are smoother, higher values "
        "rougher. Used only when noise_type is diamond_square.",
        group, "param", ""};

    const std::string k_diamond_square_initial_scale =
        prefixed_parameter_key(prefix, "diamond_square_initial_scale");
    float_default_map[k_diamond_square_initial_scale] = 1.0f;
    help_map[k_diamond_square_initial_scale] = {"float", "1.0",
        "Initial perturbation amplitude for the diamond-square surface.",
        "Dimensionless: the surface is normalised to noise_target_stddev_metres afterwards, so only the ratio of this to diamond_square_roughness affects the result -- its absolute value is washed out. That is why it carries no unit suffix."
        "subsequent level by diamond_square_roughness. Used only when "
        "noise_type is diamond_square.",
        group, "param", ""};
}

void declare_synthetic_noise_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    (void)bool_default_map;
    (void)int_default_map;
    declare_synthetic_noise_type_parameter(string_default_map, help_map, prefix);
    declare_synthetic_noise_target_stddev_parameter(float_default_map, help_map, prefix);
    declare_synthetic_noise_relief_parameter(float_default_map, help_map, prefix);
    declare_diamond_square_parameters(float_default_map, help_map, prefix);
}

// -----------------------------------------------------------------------------
// Fetch tools
// -----------------------------------------------------------------------------
void declare_fetch_bbox_parameters(std::map<std::string, std::string>& string_default_map,
                                   lsdtt_help_registry&                help_map)
{
    string_default_map["west"] = "";
    help_map["west"] = {"string","",
        "Western edge of the bounding box, WGS84 decimal degrees.",
        "Required. Must be less than east. Range -180 to 180.",
        "Bounding box", "param", ""};

    string_default_map["east"] = "";
    help_map["east"] = {"string","",
        "Eastern edge of the bounding box, WGS84 decimal degrees.",
        "Required. Must be greater than west. Range -180 to 180.",
        "Bounding box", "param", ""};

    string_default_map["south"] = "";
    help_map["south"] = {"string","",
        "Southern edge of the bounding box, WGS84 decimal degrees.",
        "Required. Must be less than north. Range -90 to 90.",
        "Bounding box", "param", ""};

    string_default_map["north"] = "";
    help_map["north"] = {"string","",
        "Northern edge of the bounding box, WGS84 decimal degrees.",
        "Required. Must be greater than south. Range -90 to 90.",
        "Bounding box", "param", ""};
}

void declare_fetch_bounds_file_parameters(std::map<std::string, std::string>& string_default_map,
                                          lsdtt_help_registry&                help_map)
{
    string_default_map["bounds_fname"] = "";
    help_map["bounds_fname"] = {"string","",
        "Vector file (or raster) defining the area, instead of west/south/east/north.",
        "FlatGeobuf (.fgb), GeoPackage (.gpkg) or any other OGR-readable "
        "format. The layer must contain either one or more polygons (the "
        "envelope of all of them is used), or exactly two points marking "
        "opposite corners. For two points either diagonal works -- "
        "lower-left/upper-right and lower-right/upper-left give the same box, "
        "so the orientation is worked out automatically. If the layer is in a "
        "projected CRS such as UTM it is reprojected to WGS84 for you, so a "
        "box drawn in a GIS can be used directly. A raster with a "
        "georeference (e.g. an existing DEM GeoTIFF) works too: its footprint "
        "is reprojected to WGS84. Set this or the four edge parameters, not "
        "both.",
        "Bounding box", "param", ""};

    string_default_map["bounds_path"] = "";
    help_map["bounds_path"] = {"string","",
        "Optional directory for bounds_fname.",
        "If empty, read_path is used. Ignored when bounds_fname is an "
        "absolute path.",
        "Bounding box", "param", ""};

    string_default_map["bounds_layer"] = "";
    help_map["bounds_layer"] = {"string","",
        "Layer name within bounds_fname.",
        "If empty, the first layer is used (with a warning when the file has "
        "more than one).",
        "Bounding box", "param", ""};
}

void declare_linear_feature_fetch_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map)
{
    const std::string group = "Linear feature frame";
    const std::string parent = "fetch_from_linear_feature";

    bool_default_map["fetch_from_linear_feature"] = false;
    help_map["fetch_from_linear_feature"] = {"bool","false",
        "Fetch along a feature (a line, two points or a polygon) in fetch_feature_fname.",
        "The area is a rectangle aligned with the feature: for a line or two "
        "points, the geodesic between the end points (only the end points of a "
        "curved or multipart line are used) lengthened by the end buffer and "
        "linear_feature_fetch_width wide; for polygons, their minimum-area "
        "enclosing rectangle, whose longest side is the axis. With "
        "rotate_after_fetch the DEM is written in an oblique Mercator frame "
        "with the feature along +x. Cannot be combined with west/south/east/"
        "north, bounds_fname, bounds_mode vector or outlets, or size_pixels.",
        group, "switch", ""};

    string_default_map["fetch_feature_fname"] = "";
    help_map["fetch_feature_fname"] = {"string","",
        "Vector file holding the line, two points or polygon(s).",
        "Any OGR-readable format (GeoPackage, FlatGeobuf, shapefile, ...) with a "
        "CRS; it is reprojected to WGS84. Refused before any download: no CRS, "
        "a mix of geometry kinds, one point or more than two, a closed line "
        "(use a polygon), or end points less than three cells apart.",
        group, "io", parent};

    string_default_map["fetch_feature_path"] = "";
    help_map["fetch_feature_path"] = {"string","",
        "Optional directory for fetch_feature_fname.",
        "If empty, read_path is used. Ignored when fetch_feature_fname is an "
        "absolute path.",
        group, "io", parent};

    string_default_map["fetch_feature_layer"] = "";
    help_map["fetch_feature_layer"] = {"string","",
        "Layer name within fetch_feature_fname.",
        "If empty, the first layer is used (with a warning when the file has "
        "more than one).",
        group, "io", parent};

    float_default_map["linear_feature_fetch_width"] = 0.0f;
    help_map["linear_feature_fetch_width"] = {"float","0.0",
        "Full width (m) of the area across a line or two points, centred on it.",
        "Required (> 0) for a line or two points. Ignored, with a warning, for "
        "polygons, whose own short side is the width.",
        group, "param", parent};

    float_default_map["linear_feature_fetch_end_buffer"] = -1.0f;
    help_map["linear_feature_fetch_end_buffer"] = {"float","-1.0",
        "Buffer (m) added beyond each end of the axis.",
        "-1 means 10 grid cells. For polygons it is added on all four sides.",
        group, "param", parent};

    bool_default_map["rotate_after_fetch"] = true;
    help_map["rotate_after_fetch"] = {"bool","true",
        "Write the DEM in the feature's oblique frame (feature along +x).",
        "true: the DEM is warped straight from the source tiles into the frame "
        "(one resampling) and written as <write_fname>_oblique_DEM.tif with "
        "<write_fname>_oblique_DEM_rotation_frame.txt, which records the frame "
        "and the lattice to rotate products back onto; target_epsg (or the "
        "analysis CRS of the frame centre) is that rotate-back target. false: "
        "an ordinary fetch in target_epsg (or UTM) covering the frame's "
        "footprint, with no frame file. Read only when "
        "fetch_from_linear_feature is true.",
        group, "param", parent};
}

void declare_linear_feature_rotation_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map)
{
    const std::string group = "Oblique frame rotation";
    const std::string parent = "rotate_to_linear_feature";

    bool_default_map["rotate_to_linear_feature"] = false;
    help_map["rotate_to_linear_feature"] = {"bool","false",
        "Rotate the DEM into the oblique frame of a feature (a line, two points or a polygon).",
        "Replaces the reprojection step: the DEM is warped once, straight from "
        "its own CRS, into an oblique Mercator frame with the feature along +x "
        "(on the DEM's own datum), and every later step runs on the rotated DEM. "
        "The frame is a rectangle aligned with the feature: for a line or two "
        "points, the geodesic between the end points lengthened by the end "
        "buffer and linear_feature_rotation_width wide; for polygons, their "
        "minimum-area enclosing rectangle. Writes <write_fname>_oblique.tif and "
        "<write_fname>_oblique_rotation_frame.txt, which records the frame and "
        "the lattice to rotate products back onto (the DEM's own grid when it "
        "is projected in metres, otherwise the analysis CRS of the frame "
        "centre). Cannot be combined with rotate_to_data_footprint, "
        "force_reproject or run_rotate_back.",
        group, "switch", ""};

    string_default_map["rotation_feature_fname"] = "";
    help_map["rotation_feature_fname"] = {"string","",
        "Vector file holding the line, two points or polygon(s).",
        "Any OGR-readable format with a CRS; it is reprojected to the DEM's own "
        "geographic CRS. Refused before any warp: no CRS, a mix of geometry "
        "kinds, one point or more than two, a closed line (use a polygon), or "
        "end points less than three cells apart.",
        group, "io", parent};

    string_default_map["rotation_feature_path"] = "";
    help_map["rotation_feature_path"] = {"string","",
        "Optional directory for rotation_feature_fname.",
        "If empty, read_path is used. Ignored when rotation_feature_fname is an "
        "absolute path.",
        group, "io", parent};

    string_default_map["rotation_feature_layer"] = "";
    help_map["rotation_feature_layer"] = {"string","",
        "Layer name within rotation_feature_fname.",
        "If empty, the first layer is used (with a warning when the file has "
        "more than one).",
        group, "io", parent};

    float_default_map["linear_feature_rotation_width"] = 0.0f;
    help_map["linear_feature_rotation_width"] = {"float","0.0",
        "Full width (m) of the frame across a line or two points, centred on it.",
        "Required (> 0) for a line or two points. Ignored, with a warning, for "
        "polygons, whose own short side is the width.",
        group, "param", parent};

    float_default_map["linear_feature_rotation_end_buffer"] = -1.0f;
    help_map["linear_feature_rotation_end_buffer"] = {"float","-1.0",
        "Buffer (m) added beyond each end of the frame axis.",
        "-1 means 10 grid cells. For polygons and rotate_to_data_footprint it "
        "is added on all four sides.",
        group, "param", ""};

    bool_default_map["rotate_to_data_footprint"] = false;
    help_map["rotate_to_data_footprint"] = {"bool","false",
        "Rotate the DEM into the frame of its own valid cells.",
        "The frame is the minimum-area rectangle enclosing the DEM's valid "
        "cells (the polygon rule of rotate_to_linear_feature applied to their "
        "convex hull), so an elongated DEM is turned to lie along +x with as "
        "little nodata as possible. Outputs and frame file as "
        "rotate_to_linear_feature. Cannot be combined with "
        "rotate_to_linear_feature, force_reproject or run_rotate_back.",
        group, "switch", ""};
}

void declare_rotate_back_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map)
{
    const std::string group = "Rotate back";
    const std::string parent = "run_rotate_back";

    bool_default_map["run_rotate_back"] = false;
    help_map["run_rotate_back"] = {"bool","false",
        "Warp rasters from an oblique frame back onto the frame file's target lattice.",
        "A separate mode: no other preprocessing step runs, and switching one on "
        "in the parameter file is refused. Each raster in rotate_back_file_list "
        "is written to write_path as <its name>_rotated_back.tif. A raster at "
        "the frame's cell size goes to the target cell size, any other keeps its "
        "own; the output lattice is anchored at the target origin, so after "
        "rotate_to_linear_feature or rotate_to_data_footprint on a projected DEM "
        "the products land on the source DEM's own pixels. A raster not in the "
        "frame CRS is refused. Vector outputs need no rotating back: they carry "
        "their CRS.",
        group, "switch", ""};

    string_default_map["rotation_frame_fname"] = "";
    help_map["rotation_frame_fname"] = {"string","",
        "The rotation frame file (<stem>_rotation_frame.txt).",
        "Relative paths are taken from read_path. If empty, the frame is rebuilt "
        "from each raster's own CRS and the target is the analysis CRS of the "
        "frame centre on the lattice anchored at (0, 0), with a warning: a "
        "source DEM origin off that lattice cannot be recovered without the "
        "frame file.",
        group, "io", parent};

    string_default_map["rotate_back_file_list"] = "";
    help_map["rotate_back_file_list"] = {"string","",
        "Text file listing the rasters to rotate back, one per line.",
        "Each line is a raster path, optionally followed by a comma and a "
        "resampling method (path,near); use near for categorical rasters "
        "(basins, junctions, masks, lithology). Blank lines and lines starting "
        "with # are skipped. The list and relative raster paths are taken from "
        "read_path.",
        group, "io", parent};

    string_default_map["rotate_back_resampling_method"] = "bilinear";
    help_map["rotate_back_resampling_method"] = {"string","bilinear",
        "Resampling method for list lines that do not give one.",
        "near, bilinear (default), cubic, cubicspline, lanczos, average, mode, "
        "max, min, med, q1, q3, sum.",
        group, "param", parent};

    int_default_map["rotate_back_epsg"] = 0;
    help_map["rotate_back_epsg"] = {"int","0",
        "EPSG code overriding the frame file's target CRS (0: the frame file's).",
        "A CRS other than the frame file's target is gridded on the lattice "
        "anchored at (0, 0), with a warning.",
        group, "param", parent};
}

void declare_true_north_bearing_parameters(
    std::map<std::string, bool>& bool_default_map,
    lsdtt_help_registry&         help_map)
{
    bool_default_map["bearings_relative_to_true_north"] = false;
    help_map["bearings_relative_to_true_north"] = {"bool","false",
        "Write bearings, aspects and orientations clockwise from true north.",
        "Off (default): they are clockwise from grid north (grid up), so in an "
        "oblique frame with the feature along +x, 90 degrees means along the "
        "feature. On: each value is converted with the grid convergence at its "
        "own location (each cell of a raster, each row of a table), which also "
        "removes the up to about 3 degree difference between UTM grid north and "
        "true north. Angles between two directions (junction angles, "
        "deflections, contrasts) are unchanged, and hillshade illumination "
        "always stays relative to grid up. Needs a raster with a CRS.",
        "Orientations", "param", ""};
}

void declare_http_request_parameters(std::map<std::string, int>& int_default_map,
                                     lsdtt_help_registry&        help_map)
{
    int_default_map["http_timeout_seconds"] = 60;
    help_map["http_timeout_seconds"] = {"int","60",
        "Maximum seconds for an individual HTTP request.",
        "Bounds how long an unresponsive server can stall the fetch.",
        "Network", "param", ""};

    int_default_map["http_connect_timeout_seconds"] = 15;
    help_map["http_connect_timeout_seconds"] = {"int","15",
        "Maximum seconds to establish a connection.",
        "Keep this short so an offline machine fails promptly.",
        "Network", "param", ""};

    int_default_map["http_max_retry"] = 3;
    help_map["http_max_retry"] = {"int","3",
        "Number of retries on a transient HTTP failure.",
        "Cloud storage occasionally returns 5xx responses; retrying absorbs "
        "those without failing the run.",
        "Network", "param", ""};
}

void declare_lithology_map_parameters(std::map<std::string, std::string>& string_default_map,
                                      lsdtt_help_registry&                help_map,
                                      const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Lithology input");
    const std::string k_fname = prefixed_parameter_key(prefix, "lithology_fname");

    string_default_map[k_fname] = "";
    help_map[k_fname] = {"string","",
        "Lithology map (in read_path): a polygon vector file or a categorical code raster.",
        "Vectors: any OGR polygon layer (GeoPackage, FlatGeobuf, GeoJSON, Shapefile), e.g. "
        "lsdtt-fetch-lithology output; it is reprojected to the DEM's CRS. Rasters: integer "
        "lithology codes, warped onto the DEM grid with nearest neighbour if not on it already.",
        group, "param", ""};

    const std::string k_source = prefixed_parameter_key(prefix, "lithology_source");
    string_default_map[k_source] = "auto";
    help_map[k_source] = {"string","auto",
        "auto, vector_map or code_raster.",
        "auto decides from the file extension (.tif, .tiff, .asc, .bil, .img, .vrt, .flt, .nc "
        "are rasters).",
        group, "param", ""};

    const std::string k_layer = prefixed_parameter_key(prefix, "lithology_layer_name");
    string_default_map[k_layer] = "";
    help_map[k_layer] = {"string","",
        "Vector layer to read. Empty (default) = the first layer.", "",
        group, "param", "lithology_fname"};

    const std::string k_field = prefixed_parameter_key(prefix, "lithology_field_name");
    string_default_map[k_field] = "lith_primary";
    help_map[k_field] = {"string","lith_primary",
        "Attribute of the vector layer that holds the lithology.",
        "Each distinct value is a lithologic unit unless lithology_crosswalk_fname groups them. "
        "lith_primary is the first lithology of each unit in lsdtt-fetch-lithology output; "
        "name (the map unit) or lith_class group more finely or more coarsely. Units are "
        "numbered within a run; compare runs by unit name.",
        group, "param", "lithology_fname"};

    const std::string k_cw = prefixed_parameter_key(prefix, "lithology_crosswalk_fname");
    string_default_map[k_cw] = "";
    help_map[k_cw] = {"string","",
        "Optional CSV (in read_path) with columns value,code[,name] mapping map values to units.",
        "When given, every map value present must be in it. Use it to group many map values "
        "into a few units; the optional name column names each code (else the values joined "
        "with '; ').",
        group, "param", "lithology_fname"};
}

void declare_write_lithology_code_raster_parameter(std::map<std::string, bool>& bool_default_map,
                                                   lsdtt_help_registry&         help_map,
                                                   const std::string&           prefix)
{
    const std::string group = prefixed_group_label(prefix, "Output");
    const std::string k = prefixed_parameter_key(prefix, "write_lithology_code_raster");
    bool_default_map[k] = true;
    help_map[k] = {"bool","true",
        "Write the lithologic unit code of every pixel.",
        "Output: <write_fname>_lithology_code.tif.", group, "switch", ""};
}

void declare_channel_region_parameters(std::map<std::string, std::string>& string_default_map,
                                       lsdtt_help_registry&                help_map,
                                       const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Channel region input");
    const std::string k_fname = prefixed_parameter_key(prefix, "channel_region_fname");

    string_default_map[k_fname] = "";
    help_map[k_fname] = {"string","",
        "Region map (in read_path): a polygon vector file or an integer code raster.",
        "User-drawn regions such as physiographic provinces, used to group channel statistics "
        "by region. Vectors: any OGR polygon layer, reprojected to the DEM's CRS. Rasters: "
        "integer region codes, warped onto the DEM grid with nearest neighbour if not on it "
        "already. Empty (default) = no region grouping.",
        group, "param", ""};

    const std::string k_source = prefixed_parameter_key(prefix, "channel_region_source");
    string_default_map[k_source] = "auto";
    help_map[k_source] = {"string","auto",
        "auto, vector_map or code_raster.",
        "auto decides from the file extension (.tif, .tiff, .asc, .bil, .img, .vrt, .flt, .nc "
        "are rasters).",
        group, "param", "channel_region_fname"};

    const std::string k_layer = prefixed_parameter_key(prefix, "channel_region_layer_name");
    string_default_map[k_layer] = "";
    help_map[k_layer] = {"string","",
        "Vector layer to read. Empty (default) = the first layer.", "",
        group, "param", "channel_region_fname"};

    const std::string k_field = prefixed_parameter_key(prefix, "channel_region_field_name");
    string_default_map[k_field] = "";
    help_map[k_field] = {"string","",
        "Attribute of the vector layer whose values name the regions.",
        "Required for a vector map; each distinct value is one region. Ignored for a code "
        "raster.",
        group, "param", "channel_region_fname"};
}

void declare_lithology_ksn_summary_parameters(std::map<std::string, float>&       float_default_map,
                                              std::map<std::string, int>&         int_default_map,
                                              std::map<std::string, std::string>& string_default_map,
                                              lsdtt_help_registry&                help_map,
                                              const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Lithology summaries");

    const std::string k_buffer = prefixed_parameter_key(prefix, "lithology_contact_buffer_metres");
    float_default_map[k_buffer] = 0.0f;
    help_map[k_buffer] = {"float","0.0",
        "Exclude channel nodes closer than this to a mapped lithologic contact from the "
        "per-unit k_sn statistics.",
        "Excluded nodes are still written, flagged near_contact = 1. 0 keeps every node. "
        "Mapped contacts carry a positional error of about 0.5 mm at map scale (about 125 m "
        "at 1:250,000).",
        group, "param", "lithology_fname"};

    const std::string k_min_reach = prefixed_parameter_key(prefix, "lithology_min_reach_pixels");
    int_default_map[k_min_reach] = 10;
    help_map[k_min_reach] = {"int","10",
        "Minimum number of channel nodes for a lithologic reach to enter the per-unit "
        "statistics.",
        "A lithologic reach is a run of channel nodes in one unit on one branch. It must also "
        "span a chi interval greater than 0.",
        group, "param", "lithology_fname"};

    const std::string k_cluster = prefixed_parameter_key(prefix, "lithology_bootstrap_cluster");
    string_default_map[k_cluster] = "reach";
    help_map[k_cluster] = {"string","reach",
        "Resampling unit of the bootstrap interval of the per-unit median k_sn: reach or basin.",
        "reach resamples lithologic reaches within a unit; basin resamples whole basins. "
        "Reaches on one trunk channel are correlated, so the reach interval is a lower bound "
        "on the uncertainty.",
        group, "param", "lithology_fname"};
}

void declare_lithology_network_summary_parameters(std::map<std::string, float>& float_default_map,
                                                  lsdtt_help_registry&          help_map,
                                                  const std::string&            prefix)
{
    const std::string group = prefixed_group_label(prefix, "Lithology summaries");
    const std::string k = prefixed_parameter_key(prefix, "lithology_min_majority_fraction");
    float_default_map[k] = 0.8f;
    help_map[k] = {"float","0.8",
        "Minimum fraction of a channel link or junction arm in its majority lithologic unit "
        "for it to count as that unit.",
        "The fraction is length-weighted over the sampled nodes. No record is dropped: links "
        "below it form the group 'mixed (<x% one unit)' and junctions whose arms fall short "
        "the group 'arms disagree'; both are summarised but left out of every comparison "
        "(never compared, never part of 'rest'). Every record keeps its majority unit and "
        "fraction in the per-record columns.",
        group, "param", "lithology_fname"};
}

void declare_lithology_bootstrap_parameters(std::map<std::string, int>& int_default_map,
                                            lsdtt_help_registry&        help_map,
                                            const std::string&          prefix)
{
    const std::string group = prefixed_group_label(prefix, "Lithology summaries");

    const std::string k_reps = prefixed_parameter_key(prefix, "lithology_bootstrap_replicates");
    int_default_map[k_reps] = 1000;
    help_map[k_reps] = {"int","1000",
        "Number of bootstrap replicates for the per-unit confidence intervals.",
        "0 = no bootstrap; the interval columns are written as NaN.",
        group, "param", "lithology_fname"};

    const std::string k_seed = prefixed_parameter_key(prefix, "lithology_bootstrap_random_seed");
    int_default_map[k_seed] = 5489;
    help_map[k_seed] = {"int","5489",
        "Random seed for the per-unit bootstrap.",
        "lsdtt-chi-analysis: each m/n gets its own std::mt19937_64 from "
        "seed_seq{seed, round(m/n * 10000)}, so an interval does not depend on the other m/n "
        "values in the list. lsdtt-channel-network-analysis: one std::mt19937_64 per run. A "
        "fixed seed reproduces the intervals. Same default as bootstrap_random_seed.",
        group, "param", "lithology_fname"};
}

void declare_stream_power_n_values_parameter(std::map<std::string, std::string>& string_default_map,
                                             lsdtt_help_registry&                help_map,
                                             const std::string&                  prefix)
{
    const std::string k = prefixed_parameter_key(prefix, "stream_power_n_values");
    string_default_map[k] = "1";
    help_map[k] = {"string", "1",
        "Slope exponent(s) n of the stream power law E = K A^m S^n, a comma list (e.g. 1,2).",
        "Used by write_relative_erodibility and the erodibility calibration (erosion_rate_fname): "
        "each n gives its own rows, with m = (m/n) n for the m/n of the chi run. Every value must "
        "be > 0. Replaces relative_erodibility_n_values (still accepted, with a warning).",
        prefixed_group_label(prefix, "Stream power"), "param", ""};
}

void declare_erosion_rate_input_parameters(std::map<std::string, float>&       float_default_map,
                                           std::map<std::string, std::string>& string_default_map,
                                           lsdtt_help_registry&                help_map,
                                           const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Erosion-rate samples");
    const std::string dep = "erosion_rate_fname";
    auto str = [&](const char* key, const char* def, const char* text, const char* guidance,
                   const char* depends) {
        const std::string k = prefixed_parameter_key(prefix, key);
        string_default_map[k] = def;
        help_map[k] = {"string", def, text, guidance, group, "param", depends};
    };
    auto flt = [&](const char* key, float def, const char* def_text, const char* text,
                   const char* guidance) {
        const std::string k = prefixed_parameter_key(prefix, key);
        float_default_map[k] = def;
        help_map[k] = {"float", def_text, text, guidance, group, "param", dep};
    };

    str("erosion_rate_fname", "",
        "Basin-averaged erosion-rate samples: a CSV with latitude/longitude (EPSG:4326) or an "
        "OGR point layer in any CRS. Empty = off.",
        "One row per sample: the basin outlet, the erosion rate and its 1-sigma uncertainty, "
        "optionally the reported drainage area. lsdtt-fetch-erosion-rates writes this format "
        "from OCTOPUS.",
        "");
    str("erosion_rate_layer_name", "", "Layer of a vector erosion-rate file.",
        "Empty takes the only (or first) layer. Not used for CSV.", dep.c_str());
    str("erosion_rate_sample_id_column", "sample_id", "Column holding the sample identifier.",
        "When the column is missing or a value is blank the sample is called row_<n> (n counts "
        "data rows from 1).", dep.c_str());
    str("erosion_rate_latitude_column", "latitude", "CSV column holding the outlet latitude.",
        "WGS84 decimal degrees. Not used for vector input.", dep.c_str());
    str("erosion_rate_longitude_column", "longitude", "CSV column holding the outlet longitude.",
        "WGS84 decimal degrees. Not used for vector input.", dep.c_str());
    str("erosion_rate_column", "erosion_rate", "Column holding the basin-averaged erosion rate.",
        "In the units of erosion_rate_units. A missing column is an error.", dep.c_str());
    str("erosion_rate_uncertainty_column", "erosion_rate_uncertainty",
        "Column holding the 1-sigma uncertainty of the erosion rate.",
        "Same units as the rate. A missing column is an error; a blank or non-positive value "
        "marks that sample missing_uncertainty.", dep.c_str());
    str("erosion_rate_units", "",
        "Units of the erosion rate and its uncertainty (required, no default): mm_per_kyr, "
        "m_per_myr, cm_per_kyr, mm_per_yr or m_per_yr.",
        "mm_per_kyr equals m_per_myr (1e-6 m/yr). Mass fluxes (g cm^-2 yr^-1) are refused: "
        "divide by the rock density first. When the file has an erosion_rate_units column every "
        "row must agree with this parameter.", dep.c_str());
    str("erosion_rate_reported_area_column", "reported_area_km2",
        "Column holding each sample's reported drainage area in km^2.",
        "Used for area-matched snapping and the DEM-vs-reported area check. Empty, or a column "
        "not in the file, skips both.", dep.c_str());
    str("erosion_rate_snapping_method", "reported_area",
        "How a sample is snapped to the channel network: reported_area, expanding_square, "
        "downstream, largest_area or none.",
        "reported_area takes the nearest channel cell whose drainage area is within "
        "erosion_rate_area_mismatch_tolerance of the reported area, else the cell closest in "
        "area. Without a reported area it falls back to expanding_square (nearest channel "
        "cell). largest_area takes the cell with the most contributing pixels within "
        "erosion_rate_snap_max_distance_metres.", dep.c_str());
    flt("erosion_rate_snap_max_distance_metres", 300.0f, "300",
        "Search radius for snapping a sample, in metres.",
        "Converted to pixels with the cell size and rounded up, so a parameter file stays valid "
        "across resolutions.");
    flt("erosion_rate_area_mismatch_tolerance", 1.25f, "1.25",
        "Ratio tolerance (> 1) between the DEM drainage area and the reported area.",
        "A basin with A_dem / A_reported outside [1/tol, tol] is area_mismatch. Also the "
        "matching tolerance of reported_area snapping.");
    str("erosion_rate_area_mismatch_action", "exclude",
        "What an area_mismatch basin does: exclude (written, kept out of every fit) or flag "
        "(kept in the fits and flagged).", "", dep.c_str());
    str("erosion_rate_nested_basins", "keep",
        "Nested sample basins: keep, outermost or innermost.",
        "keep uses all basins and resamples nesting groups in the bootstrap; outermost keeps "
        "only basins nested in no other sample; innermost only basins containing no other "
        "sample.", dep.c_str());
}

void declare_erodibility_calibration_parameters(std::map<std::string, bool>&        bool_default_map,
                                                std::map<std::string, int>&         int_default_map,
                                                std::map<std::string, float>&       float_default_map,
                                                std::map<std::string, std::string>& string_default_map,
                                                lsdtt_help_registry&                help_map,
                                                const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Erodibility calibration");
    const std::string dep = "erosion_rate_fname";
    auto str = [&](const char* key, const char* def, const char* text, const char* guidance) {
        const std::string k = prefixed_parameter_key(prefix, key);
        string_default_map[k] = def;
        help_map[k] = {"string", def, text, guidance, group, "param", dep};
    };
    auto integer = [&](const char* key, int def, const char* text, const char* guidance) {
        const std::string k = prefixed_parameter_key(prefix, key);
        int_default_map[k] = def;
        help_map[k] = {"int", std::to_string(def), text, guidance, group, "param", dep};
    };
    auto flt = [&](const char* key, float def, const char* def_text, const char* text,
                   const char* guidance) {
        const std::string k = prefixed_parameter_key(prefix, key);
        float_default_map[k] = def;
        help_map[k] = {"float", def_text, text, guidance, group, "param", dep};
    };
    auto flag = [&](const char* key, bool def, const char* text, const char* guidance) {
        const std::string k = prefixed_parameter_key(prefix, key);
        bool_default_map[k] = def;
        help_map[k] = {"bool", def ? "true" : "false", text, guidance, group, "param", dep};
    };

    str("erosion_rate_ksn_statistic", "reach_median",
        "Basin k_sn that drives the calibration: reach_median, chi_z_slope or segment_median.",
        "reach_median is the median k_sn of the basin's qualifying reaches (split at "
        "confluences only), with a bootstrap interval. The others are written as diagnostics "
        "either way.");
    integer("erosion_rate_min_reach_pixels", 10,
            "Minimum channel nodes for a basin reach to qualify.",
            "Reaches run between confluences, sources and the outlet.");
    integer("erosion_rate_min_reaches_per_basin", 5,
            "Minimum qualifying reaches for a basin to enter the calibration.",
            "Basins with fewer are written with status too_few_reaches and enter no fit.");
    flag("erosion_rate_fit_n", true, "Fit n and K together (York regression of ln E on ln k_sn).",
         "Attempted only with enough independent basin groups and a wide enough k_sn range; "
         "K at each fixed n is always written.");
    integer("erosion_rate_min_groups_for_n_fit", 8,
            "Minimum independent basin groups for fitting n.",
            "Nested basins form one group. Below this only K at fixed n is reported.");
    flt("erosion_rate_min_ksn_range_ratio", 2.0f, "2.0",
        "Minimum max/min k_sn ratio over the used basins for fitting n.", "");
    flt("erosion_rate_mono_lithology_min_area_fraction", 0.8f, "0.8",
        "Minimum area fraction of a basin's dominant unit for the basin to count as "
        "mono-lithologic.", "Needs lithology_fname.");
    integer("erosion_rate_min_unit_reaches_per_basin", 3,
            "Minimum qualifying lithologic reaches of a unit in a basin for the unit-in-basin "
            "estimate.", "Needs lithology_fname.");
    str("erosion_rate_unit_k_methods", "mono,unit_in_basin",
        "Per-unit K estimators, comma separated: mono, unit_in_basin.", "Needs lithology_fname.");
    str("erosion_rate_property_table_method", "unit_in_basin",
        "Per-unit estimator written to the draft lithology property table.",
        "Needs lithology_fname.");
    flt("erosion_rate_systematic_relative_uncertainty", 0.0f, "0",
        "Relative systematic uncertainty of every erosion rate (e.g. the production-rate "
        "term), added to the K interval in quadrature.", "0 adds nothing.");
    integer("erosion_rate_bootstrap_replicates", 1000,
            "Bootstrap replicates of the calibration intervals.", "0 = no intervals.");
    integer("erosion_rate_bootstrap_random_seed", 5489,
            "Seed of the calibration bootstrap.",
            "One stream per m/n, separate from the k_sn summary and relative-erodibility "
            "streams.");
    flag("write_erodibility_calibration_basin_polygons", true,
         "Write the sample basins as polygons (write_vector_format).",
         "For checking the snapping against the published basins in a GIS.");
    flag("write_erodibility_calibration_property_table", true,
         "Write a draft lithology property table from the per-unit K.",
         "Needs lithology_fname.");
}

} // namespace lsdtt
