// lsdtt-raster-preprocessing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-raster-preprocessing
//
// Raster preprocessing pipeline:
//   1. Projection check and reprojection (if needed)
//   2. Threshold masking
//   2b. Binary mask application (optional)
//   2c. Clip raster to polygon (optional)
//   3. Interior nodata filling
//   4. Clip to maximum valid rectangle
//   5. Depression breaching
//   6. Pit filling
//   7. Depression polygons
//   8. Hillshade
//   8b. Multidirectional hillshade
//   9. Vector clipping
//   10. Sample rasters to points
//
// All steps after (1) operate on the reprojected raster when reprojection
// occurs.  Output rasters receive a projection suffix when reprojected:
//   _utm     - UTM zone of raster centre
//   _polar   - NSIDC polar stereographic (|lat| > 75 deg)
//   _lambert - custom Lambert Azimuthal Equal Area (wide-extent rasters)
//
// Usage:
//   lsdtt-raster-preprocessing <path> <parameter_file>
//   lsdtt-raster-preprocessing <parameter_file>   (path assumed to be ./)
//   lsdtt-raster-preprocessing -h                 (print help and exit)
//   lsdtt-raster-preprocessing -v                 (print version and exit)
//   lsdtt-raster-preprocessing -c                 (print citation and exit)
//
// Parameter file format:
//   key: value
//   Lines beginning with # are ignored.
//   All keys use snake_case. No spaces in key names.
//
// Example parameter file:
//   read_path: /data/dems/
//   read_fname: my_dem
//   write_path: /data/output/
//   write_fname: my_dem_preprocessed
//   enforce_projection: true
//   resampling_method: cubic
//   grid_spacing: -1
//   allow_lambert_projection: false
//   write_hillshade_raster: true
//   write_fill_raster: true
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <sstream>


#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_raster_reproject.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_usgs_elevation.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_breach_least_cost_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_raster_linear_features_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_raster_mosaic_cpu.hpp"
#include "../src/cpu/lsdtt_depression_hierarchy_cpu.hpp"
#include "../src/cpu/lsdtt_internal_drainage_cpu.hpp"
#include "../src/cpu/lsdtt_fill_spill_merge_cpu.hpp"
#include "../src/cpu/lsdtt_vector_rasterise_cpu.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_csv_table.hpp"

#include <chrono>
#include <filesystem>


static const std::string k_version  = "1.0.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-raster-preprocessing";

// this is for sampling raster points
struct raster_attribute_request
{
    std::string raster_path;
    std::string attribute_name;
};




static std::vector<raster_attribute_request> read_raster_attribute_csv(
    const std::string& csv_path,
    const lsdtt::lsdtt_param_parser& pp)
{
    std::ifstream in(csv_path);
    if (!in.is_open())
    {
        throw std::runtime_error(
            "Cannot open raster_attribute_fname: " + csv_path
        );
    }

    std::string header_line;
    if (!std::getline(in, header_line))
    {
        throw std::runtime_error(
            "raster_attribute_fname is empty: " + csv_path
        );
    }

    auto header = lsdtt::csv_split_row_trimmed(header_line);

    int raster_col = -1;
    int attr_col   = -1;

    for (int i = 0; i < static_cast<int>(header.size()); ++i)
    {
        if (header[i] == "raster")    raster_col = i;
        if (header[i] == "attribute") attr_col   = i;
    }

    if (raster_col < 0 || attr_col < 0)
    {
        throw std::runtime_error(
            "raster_attribute_fname must contain columns named "
            "'raster' and 'attribute'"
        );
    }

    std::vector<raster_attribute_request> requests;

    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty()) continue;

        auto tokens = lsdtt::csv_split_row_trimmed(line);
        if (static_cast<int>(tokens.size()) <= std::max(raster_col, attr_col))
        {
            continue;
        }

        raster_attribute_request req;
        req.raster_path    = pp.add_read_path_if_needed(tokens[raster_col]);
        req.attribute_name = tokens[attr_col];

        if (!req.raster_path.empty() && !req.attribute_name.empty())
        {
            requests.push_back(req);
        }
    }

    return requests;
}

static std::vector<lsdtt::lsdtt_vector_field_value>
sample_raster_to_vector_points(
    const std::string& raster_path,
    const std::string& point_path,
    const std::string& point_layer_name)
{
    lsdtt::lsdtt_data_loader raster_loader;
    raster_loader.open_raster(raster_path);

    const lsdtt::lsdtt_raster_meta& meta = raster_loader.get_meta();
    lsdtt::lsdtt_raster_tile tile = raster_loader.read_full();

    std::vector<lsdtt::lsdtt_input_point> mapped_points =
        lsdtt::load_vector_point_layer_in_raster_crs(
            point_path, meta, point_layer_name, true
        );

    std::vector<lsdtt::lsdtt_vector_field_value> sampled;
    sampled.reserve(mapped_points.size());

    for (const auto& pt : mapped_points)
    {
        if (pt.status != lsdtt::lsdtt_point_input_status::success ||
            pt.flat < 0 ||
            pt.flat >= static_cast<std::int64_t>(tile.data.size()))
        {
            sampled.emplace_back(std::monostate{});
            continue;
        }

        if (tile.mask[static_cast<std::size_t>(pt.flat)] == 0)
        {
            sampled.emplace_back(std::monostate{});
            continue;
        }

        sampled.emplace_back(
            static_cast<double>(tile.data[static_cast<std::size_t>(pt.flat)])
        );
    }

    return sampled;
}



// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt3 raster preprocessing              ||\n"
        << "|| Projection, filling, breaching, nodata repair       ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "|| If you need to condition your DEM for flow routing, ||\n"
        << "|| which mosty users will need to do, this CLI is your ||\n"
        << "|| huckleberry.                                        ||\n"
        << "=========================================================\n"
        << "|| -v or --version    : full lsdtt3 release version\n"
        << "|| This CLI's version : " << k_version << "\n"
        << "|| If CLI version number is lower than 1.0.0 it means\n"
        << "|| this CLI is still experimental and not fully tested.\n"
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
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter declarations
    //
    // Each parameter has its default value and help entry declared together
    // so that a developer can read the full contract for each parameter in
    // one place.  Parameters are grouped by the pipeline step that uses them.
    // None are const - some booleans conditionally modify others below.
    //
    // help_map entry format:
    //   { type, default, short_desc, long_desc, group, role, parent }
    // (see lsdtt_help_registry in lsdtt_param_parser.hpp)
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    lsdtt::lsdtt_help_registry                         help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = single-threaded. "
        "Applies to filling, breaching, nodata repair, and wavelet steps.",
        "Run control", "io", ""};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format.",
        "Run control", "io", ""};

    // --- Projection ---
    bool_default_map["enforce_projection"] = true;
    help_map["enforce_projection"] = {"bool","true",
        "Check and reproject the raster to a suitable projected CRS.",
        "If true and the raster is geographic or Web Mercator, it is "
        "reprojected automatically before any analysis. "
        "If false and the raster is not projected, the program exits with "
        "an error.",
        "Projection", "switch", ""};

    string_default_map["resampling_method"] = "cubic";
    help_map["resampling_method"] = {"string","cubic",
        "GDAL resampling algorithm used during reprojection.",
        "Cubic is the default here because this workflow reprojects continuous "
        "elevation surfaces. Use near or mode for categorical rasters. Accepted "
        "values: near, bilinear, cubic (default), cubicspline, "
        "lanczos, average, mode, max, min, med, q1, q3, sum.",
        "Projection", "param", "enforce_projection"};

    float_default_map["grid_spacing"] = -1.0f;
    help_map["grid_spacing"] = {"float","-1",
        "Target grid spacing in metres after reprojection.",
        "If negative (default), spacing is inferred from the input pixel size "
        "and rounded: >=20 m -> nearest 5 m; 5-20 m -> nearest 2 m; "
        "<5 m -> nearest 1 m. Set a positive value to enforce a specific spacing.",
        "Projection", "param", "enforce_projection"};

    bool_default_map["allow_lambert_projection"] = false;
    help_map["allow_lambert_projection"] = {"bool","false",
        "Use Lambert Azimuthal Equal Area for rasters spanning 3+ UTM zones.",
        "If true and the raster is very wide (>=3 UTM zones), a LAEA projection "
        "centred on the raster is used instead of UTM. "
        "Output files receive the suffix _lambert.",
        "Projection", "param", "enforce_projection"};

    int_default_map["warp_mem_mb"] = 512;
    help_map["warp_mem_mb"] = {"int","512",
        "GDALWarp working memory in MB.",
        "Larger values reduce source tile re-reads for large rasters. "
        "Reduce on memory-constrained systems (minimum ~64 MB).",
        "Projection", "param", "enforce_projection"};

    // projection forcing
    bool_default_map["force_reproject"] = false;
    help_map["force_reproject"] = {"bool","false",
        "Force reprojection to the recommended analysis CRS.",
        "If true, the raster is reprojected to the recommended CRS even when "
        "the input raster is already projected. This allows conversion from an "
        "existing projected CRS to the preferred UTM, polar stereographic, or "
        "Lambert Azimuthal Equal Area analysis projection.",
        "Projection", "switch", ""};

    bool_default_map["force_resampling"] = false;
    help_map["force_resampling"] = {"bool","false",
        "Force a user-specified grid spacing during reprojection.",
        "If true, forced_grid_spacing is used as the target grid spacing "
        "instead of the automatically estimated or user-supplied grid_spacing.",
        "Projection", "switch", ""};

    float_default_map["forced_grid_spacing"] = -1.0f;
    help_map["forced_grid_spacing"] = {"float","-1",
        "Target grid spacing in metres when force_resampling is true.",
        "Must be positive when force_resampling is true. Ignored otherwise.",
        "Projection", "param", "force_resampling"};

    // --- Vertical units ---
    string_default_map["vertical_units"] = "metres";
    help_map["vertical_units"] = {"string","metres",
        "Vertical unit of the input elevations: metres, feet or us_survey_feet.",
        "NO raster format records the vertical unit, so this CANNOT be detected "
        "and is NEVER applied automatically. Much US state lidar is delivered as "
        "NAVD88 FEET on a metre UTM grid; left uncorrected every slope, "
        "curvature and relief is too large by 3.28. Establish the vertical datum "
        "from the data provider before setting this - a wrong setting silently "
        "corrupts correct data. The tool warns when slopes look implausible.",
        // No parent: this depends on nothing. The allowed values belong in
        // the description (they are stated there) and in the type, never in
        // the parent slot, which names the parameter this one depends on.
        "Vertical units", "param", ""};

    bool_default_map["check_vertical_units_against_usgs"] = true;
    help_map["check_vertical_units_against_usgs"] = {"bool","true",
        "When the slope-plausibility warning fires, query USGS 3DEP to settle "
        "whether the elevations are metres or feet.",
        "Plausibility can only FLAG a suspicious DEM; an external reference "
        "SETTLES it, because nothing but a unit error puts a DEM 3.28x off a "
        "national elevation dataset. Runs ONLY if the warning fires AND the DEM "
        "sits in US coverage - which is not the limitation it appears, since "
        "feet-on-a-metre-grid is a US delivery convention. Sends 5 pixel "
        "coordinates to epqs.nationalmap.gov. Set false to keep the run offline. "
        "Nothing is ever converted automatically.",
        "Vertical units", "param", ""};

    // --- Threshold masking ---
    float_default_map["min_elevation"] = 0.0f;
    help_map["min_elevation"] = {"float","0.0",
        "Pixels at or below this elevation (metres) are set to nodata.",
        "Default of 0.0 removes ocean and below-sea-level artefacts common "
        "in global DEMs. Set to a large negative number to disable.",
        "Threshold masking", "param", ""};

    float_default_map["maximum_elevation"] = 30000.0f;
    help_map["maximum_elevation"] = {"float","30000",
        "Pixels above this elevation (metres) are set to nodata.",
        "Default covers the full range of terrestrial elevations. "
        "Reduce to remove artefacts in low-relief study areas.",
        "Threshold masking", "param", ""};

    // --- Binary mask application ---
    bool_default_map["run_binary_mask"] = false;
    help_map["run_binary_mask"] = {"bool","false",
        "Apply an external binary mask raster to the DEM.",
        "Loads the raster at binary_mask_fname and calls "
        "raster_mask_by_binary_keep: pixels where the mask > 0.5 are kept; "
        "all others are set to nodata. "
        "The mask must cover the same spatial extent and have the same "
        "cell size and origin as the (possibly reprojected) DEM. "
        "Runs after threshold masking. "
        "Output suffix: _binary_masked.tif",
        "Binary mask", "switch", ""};

    string_default_map["binary_mask_fname"] = "";
    help_map["binary_mask_fname"] = {"string","",
        "Full path to the binary mask raster.",
        "Any GDAL-readable single-band raster is accepted. "
        "Pixels > 0.5 are treated as keep (foreground); "
        "all others set the DEM pixel to nodata. "
        "Required when run_binary_mask = true.",
        "Binary mask", "param", "run_binary_mask"};

    // --- Interior nodata filling ---
    bool_default_map["write_filled_nodata_raster"] = false;
    help_map["write_filled_nodata_raster"] = {"bool","false",
        "Write a DEM with interior nodata pixels filled.",
        "Only fills nodata pixels completely surrounded by valid data. "
        "Does not fill nodata that touches the raster edge.",
        "Interior nodata filling", "switch", ""};

    string_default_map["nodata_fill_method_choice"] = "laplacian";
    help_map["nodata_fill_method_choice"] = {"string","laplacian",
        "Method for filling interior nodata pixels.",
        "laplacian: smooth Dirichlet interpolation (recommended). "
        "iterative_weighted_average: fast, well-suited to thin voids.",
        "Interior nodata filling", "param", "write_filled_nodata_raster"};

    int_default_map["nodata_fill_max_iter"] = 10000;
    help_map["nodata_fill_max_iter"] = {"int","10000",
        "Maximum iterations for iterative nodata fill methods.",
        "Only used when nodata_fill_method_choice is laplacian or "
        "iterative_weighted_average.",
        "Interior nodata filling", "param", "write_filled_nodata_raster"};

    // --- Clip to rectangle ---
    bool_default_map["write_clipped_rectangle_raster"] = false;
    help_map["write_clipped_rectangle_raster"] = {"bool","false",
        "Write a DEM clipped to the largest rectangle of all-valid pixels.",
        "Run write_filled_nodata_raster first to maximise the rectangle size.",
        "Clip to rectangle", "switch", ""};

    // --- Depression breaching ---
    bool_default_map["write_breached_raster"] = false;
    help_map["write_breached_raster"] = {"bool","false",
        "Write a depression-breached DEM.",
        "Use breach_algorithm_choice to select the algorithm.",
        "Depression breaching", "switch", ""};

    // --- Smoothing / conditioning ---
    string_default_map["smooth_algorithm_choice"] = "none";
    help_map["smooth_algorithm_choice"] = {"string","none",
        "Edge-preserving smoothing applied BEFORE breaching and filling: none, "
        "perona_malik or coherence_diffusion.",
        "For work that reads heights off surfaces - terrace treads, valley "
        "floors, anything the fluvial-surface ladder bands - the riser is the "
        "signal, not the noise. A mean or Gaussian filter removes DEM noise "
        "and the tread/riser break in equal measure, which moves every height "
        "break downstream. Perona-Malik (1990) diffuses within a facet and "
        "stops at a gradient above K, so the riser survives. It is a "
        "CONDITIONING step: the smoothed raster is written as its own file and "
        "everything downstream should read that, so which DEM an analysis ran "
        "on is visible rather than inferred.",
        "Smoothing", "param", ""};

    int_default_map["perona_malik_iterations"] = 10;
    help_map["perona_malik_iterations"] = {"int","10",
        "Explicit diffusion timesteps.",
        "More iterations remove more within-facet noise. The riser is "
        "preserved at every count, so this trades runtime against how quiet "
        "the treads become rather than against how much structure survives.",
        "Smoothing", "param", ""};

    float_default_map["perona_malik_lambda"] = 0.20f;
    help_map["perona_malik_lambda"] = {"float","0.20",
        "Diffusion timestep. Must be in (0, 0.25].",
        "0.25 is the stability limit of the explicit 4-neighbour stencil. A "
        "value above it is REJECTED rather than warned about: an unstable run "
        "does not give a worse DEM, it gives a meaningless one.",
        "Smoothing", "param", ""};

    string_default_map["perona_malik_conductance"] = "exponential";
    help_map["perona_malik_conductance"] = {"string","exponential",
        "Conductance function: exponential or quadratic.",
        "exponential is exp(-(s/K)^2), Perona & Malik eq. 3, and falls away "
        "faster above K, so it preserves a sharp riser more strictly. "
        "quadratic is 1/(1+(s/K)^2), eq. 4, which keeps diffusing weakly "
        "across strong edges and is the gentler choice on a noisy DEM.",
        "Smoothing", "param", ""};

    float_default_map["perona_malik_k_quantile"] = 0.90f;
    help_map["perona_malik_k_quantile"] = {"float","0.90",
        "Quantile of the DEM's own neighbour-difference distribution used to "
        "set the edge threshold K. 0 = state K directly instead.",
        "K is in metres of rise between adjacent cells, the same units as the "
        "differences the stencil takes - not a central-difference gradient, "
        "which averages across two cells and so cannot see the one-cell noise "
        "the diffusion does. A K carried over from another landscape means "
        "nothing, so following Black et al. (1998) it is read off the "
        "distribution of the DEM in hand. The value actually used is printed.",
        "Smoothing", "param", ""};

    float_default_map["perona_malik_k_metres_per_cell"] = 0.0f;
    help_map["perona_malik_k_metres_per_cell"] = {"float","0",
        "Edge threshold K stated directly, in metres of rise per cell. Only "
        "read when perona_malik_k_quantile is 0.",
        "Use this to pin K across a set of DEMs that must be conditioned "
        "identically; the quantile would otherwise give each one its own.",
        "Smoothing", "param", ""};

    bool_default_map["write_stencil_difference_raster"] = false;
    help_map["write_stencil_difference_raster"] = {"bool","false",
        "Write the neighbour-difference raster (metres) K was read from.",
        "Per cell, the largest absolute elevation difference to a valid "
        "4-neighbour - the exact quantity the conductance is evaluated on. "
        "Write it the first time a DEM is conditioned: a quantile is only "
        "interpretable against the distribution it was taken from.",
        "Smoothing", "output", ""};

    int_default_map["coherence_iterations"] = 10;
    help_map["coherence_iterations"] = {"int","10",
        "Coherence-enhancing diffusion steps.",
        "More iterations bridge longer breaks along a feature and smooth "
        "further within a facet.",
        "Smoothing", "param", ""};

    float_default_map["coherence_time_step"] = 0.15f;
    help_map["coherence_time_step"] = {"float","0.15",
        "Diffusion timestep tau. Must be in (0, 0.25].",
        "Same explicit-stencil stability limit as perona_malik_lambda, and "
        "rejected the same way.",
        "Smoothing", "param", ""};

    float_default_map["coherence_derivative_radius_metres"] = 0.0f;
    help_map["coherence_derivative_radius_metres"] = {"float","0",
        "Gradient pre-smoothing radius before the structure tensor is formed. "
        "0 = one cell.",
        "Defaults are stated in CELLS, not metres, because a radius that is "
        "sensible on 1 m lidar is sub-pixel on 30 m COP30 and the structure "
        "tensor then has nothing to integrate over.",
        "Smoothing", "param", ""};

    float_default_map["coherence_integration_radius_metres"] = 0.0f;
    help_map["coherence_integration_radius_metres"] = {"float","0",
        "Structure-tensor integration radius. 0 = three cells.",
        "This sets the length scale over which a feature has to look linear "
        "before the filter will smooth along it. Too small and every noisy "
        "streak reads as coherent; too large and a curving riser stops "
        "registering as one feature. See the note on cells vs metres under "
        "coherence_derivative_radius_metres.",
        "Smoothing", "param", ""};

    float_default_map["coherence_contrast"] = 0.2f;
    help_map["coherence_contrast"] = {"float","0.2",
        "Coherence contrast C: how selective the along-feature smoothing is.",
        "The coherence it is compared against is normalised to [0, 1], so C is "
        "scale-invariant and transfers between DEMs. Raise it to smooth only "
        "the most coherent threads; lower it to smooth more of the field.",
        "Smoothing", "param", ""};

    float_default_map["coherence_min_diffusivity"] = 0.001f;
    help_map["coherence_min_diffusivity"] = {"float","0.001",
        "Across-feature diffusivity floor, in (0, 1).",
        "The eigenvalue perpendicular to the feature. Near zero means no "
        "smoothing at all across a riser, which is the point; raise it only to "
        "deliberately soften the across-feature behaviour.",
        "Smoothing", "param", ""};

    bool_default_map["write_carve_diff_raster"] = false;
    help_map["write_carve_diff_raster"] = {"bool","false",
        "Write the difference raster: breached minus original.",
        "Negative values show where breaching lowered the surface.",
        "Depression breaching", "switch", ""};

    bool_default_map["breach_before_fill"] = false;
    help_map["breach_before_fill"] = {"bool","false",
        "Breach depressions before filling pits.",
        "Recommended for lidar DEMs where channelised flow paths "
        "must be carved before residual pits are filled.",
        "Depression breaching", "switch", ""};

    bool_default_map["fill_remainder_after_breach"] = true;
    help_map["fill_remainder_after_breach"] = {"bool","true",
        "Fill pits that could not be fully breached.",
        "Only used when breach_algorithm_choice is least_cost. If cells "
        "without a downslope neighbour remain after the breach, they are "
        "filled with fill_min_slope and fill_flat_resolution_choice, as the "
        "barnes2014 pit fill would, so the breached DEM is routable.",
        "Depression breaching", "param", "breach_before_fill"};

    string_default_map["breach_algorithm_choice"] = "lindsay2016";
    help_map["breach_algorithm_choice"] = {"string","lindsay2016",
        "Algorithm for depression breaching.",
        "lindsay2016: Lindsay (2016) efficient breaching. "
        "least_cost: Lindsay (2020) pit-centric least-cost breaching.",
        "Depression breaching", "param", "write_breached_raster"};

    // The declared default is the outlet-mode one; when the file does not
    // set the key, the whitebox edge policy uses 100 instead (see the
    // options block below and default_breach_lc_max_dist).
    int_default_map["breach_lc_max_dist"] = lsdtt::k_breach_lc_default_max_dist_outlet;
    help_map["breach_lc_max_dist"] = {"int","1000 (outlet) / 100 (whitebox)",
        "Maximum breach path length in cells for least-cost breaching.",
        "Counted in steps from the pit (a diagonal step counts as one), not "
        "as a radius. Pits that need a longer path are left unbreached (and "
        "filled if fill_remainder_after_breach is true). 1 to 32767. If not "
        "set, 1000 with breach_lc_edge_policy outlet (long enough to bridge "
        "large water bodies, so a far-off pit is not carved to the nearest "
        "raster edge instead) and 100 with whitebox (the WhiteboxTools "
        "default; longer paths there can carve across divides). Only used "
        "when breach_algorithm_choice is least_cost.",
        "Depression breaching", "param", "write_breached_raster"};

    float_default_map["breach_lc_max_cost"] = std::numeric_limits<float>::max();
    help_map["breach_lc_max_cost"] = {"float","unlimited",
        "Maximum cost of a least-cost breach path.",
        "The cost of a path is the depth each cell on it is lowered, summed: "
        "metres when breach_lc_minimize_dist is false; with "
        "breach_lc_minimize_dist true each step is weighted by its length in "
        "cells (1 or 1.414), so the units are metres x cells. A pit whose "
        "cheapest path costs more is left unbreached. Only used when "
        "breach_algorithm_choice is least_cost.",
        "Depression breaching", "param", "write_breached_raster"};

    bool_default_map["breach_lc_minimize_dist"] = false;
    help_map["breach_lc_minimize_dist"] = {"bool","false",
        "Weight least-cost breach paths by their length.",
        "false: a path costs the total depth lowered (m), so the search "
        "prefers shallow cuts however long. true: each step is weighted by "
        "its length in cells, favouring shorter paths; breach_lc_max_cost is "
        "then in metres x cells. Only used when breach_algorithm_choice is "
        "least_cost.",
        "Depression breaching", "param", "write_breached_raster"};

    string_default_map["breach_lc_edge_policy"] = "outlet";
    help_map["breach_lc_edge_policy"] = {"string","outlet",
        "Whether a least-cost breach may drain off the raster edge.",
        "outlet: the raster edge is an outlet, like a nodata cell, so a pit "
        "near the edge can breach straight off the grid. whitebox: the edge "
        "is not an outlet (as in WhiteboxTools), so a breach must reach an "
        "edge cell that is already low enough, which can carve deep cuts "
        "across divides near the edge. Only used when breach_algorithm_choice "
        "is least_cost.",
        "Depression breaching", "param", "write_breached_raster"};

    float_default_map["breach_lc_flat_increment"] = 0.0f;
    help_map["breach_lc_flat_increment"] = {"float","0.0",
        "Elevation drop (m) per cell along a least-cost breach channel.",
        "0 or less selects the automatic value: twice the float32 spacing at "
        "the largest absolute elevation in the DEM (e.g. 2.4e-4 m for "
        "elevations of 1024-2047 m), the smallest drop that keeps the channel "
        "strictly downhill. Only used when breach_algorithm_choice is "
        "least_cost.",
        "Depression breaching", "param", "write_breached_raster"};

    // --- Pit filling ---
    bool_default_map["write_fill_raster"] = false;
    help_map["write_fill_raster"] = {"bool","false",
        "Write a pit-filled DEM.",
        "Use fill_algorithm_choice to select the algorithm.",
        "Pit filling", "switch", ""};

    bool_default_map["write_fill_diff_raster"] = false;
    help_map["write_fill_diff_raster"] = {"bool","false",
        "Write the difference raster: filled minus original.",
        "Positive values show where filling raised the surface."
        "If fill_agorithm_choice is barnesetal2021, this represents"
        " the water depth standing in topographic depressions.",
        "Pit filling", "switch", ""};

    string_default_map["fill_algorithm_choice"] = "barnes2014";
    help_map["fill_algorithm_choice"] = {"string","barnes2014",
        "Depression filling algorithm.",
        "barnes2014: Improved Priority-Flood (Barnes et al. 2014). The default. "
        "Imposes fill_min_slope across flats so the filled surface is routable. "
        "vincent1993: grayscale morphological reconstruction (Vincent 1993). "
        "Produces a BIT-IDENTICAL result to barnes2014 and is 2-3x faster, but "
        "fills to exactly the spill level and so requires fill_min_slope: 0 -- "
        "it cannot impose a gradient across flats, and lsdtt3's D8 router will "
        "reject a surface with true flats. barnesetal2021: Fill-Spill-Merge "
        "(https://doi.org/10.5194/esurf-9-105-2021); note this distributes a "
        "finite rainfall volume and does not guarantee all depressions are "
        "filled.",
        "Filling", "param", ""};

    string_default_map["priority_flood_queue"] = "binary";
    help_map["priority_flood_queue"] = {"string", "binary",
        "Spill-queue backend for Barnes 2014 filling.",
        "binary (default) uses the historic heap. radix uses the equivalent "
        "hybrid radix queue. Only used when fill_algorithm_choice is barnes2014.",
        "Pit filling", "param", "write_fill_raster"};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum gradient imposed between adjacent filled cells.",
        "Prevents flat regions that cause ambiguous flow direction. "
        "Increase slightly for very coarse DEMs. "
        "This does not apply to fill-spill-merge. Under "
        "fill_flat_resolution_choice barnes2014b_gradient it is the minimum "
        "gradient across flats; barnes2014b_float_step ignores it.",
        "Pit filling", "param", "write_fill_raster"};

    // inline_gradient (default) | barnes2014b_gradient | barnes2014b_float_step.
    // The barnes2014b options apply to the barnes2014 and vincent1993 fills.
    lsdtt::declare_fill_flat_resolution_parameter(string_default_map, help_map);

    float_default_map["fill_spill_merge_rainfall_rate"] = 0.1f;
    help_map["fill_spill_merge_rainfall_rate"] = {"float", "0.1",
        "Uniform rainfall depth in elevation units per cell for "
        "fill-spill-merge flow routing (Barnes et al., 2021).",
        "Spatial patterns in rainfall to route into and over "
        "depressions are not currently supported.",
        "Pit filling", "param", "write_fill_raster"    
    };
    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map);

    string_default_map["fsm_boundary_conditions"] = "oooo";
    help_map["fsm_boundary_conditions"] = {"string","oooo",
        "Boundary conditions for the depression hierarchy (north east south west).",
        "Four characters, one per edge: "
        "o = open (cells on this edge drain freely to the ocean, default). "
        "b = forced baselevel outlet (same as o for hierarchy purposes). "
        "n = no-flux (edge is NOT an ocean boundary). "
        "p = periodic (edge wraps; NOT an ocean boundary). "
        "Example: oooo = all four edges are open ocean (typical for a DEM tile). "
        "Example: nnnn = closed basin (no edge is ocean; requires interior pits).",
        "Pit filling", "param", "write_fill_raster"};
 
    // --- Pixel replacement ---
    //
    // A CSV-driven per-pixel edit: set named cells to arbitrary values. The
    // raster-driven equivalent is run_binary_mask; this covers the case where
    // you know the coordinates of a handful of bad cells and want to blank or
    // correct them without building a mask raster.
    string_default_map["replace_pixels_fname"] = "";
    help_map["replace_pixels_fname"] = {"string","",
        "Point file listing pixels to overwrite with new values.",
        "CSV, FlatGeobuf or GeoPackage. Positions come from geometry, or from "
        "latitude/longitude, or from easting/northing in the raster CRS. The "
        "replacement value is read from the column named by "
        "replace_pixels_value_field. Setting a cell to the raster's nodata value "
        "marks it nodata rather than merely writing the number. Applied before "
        "breaching, filling and nodata filling.",
        "Pixel replacement", "param", ""};

    bool_default_map["write_replaced_pixels_raster"] = false;
    help_map["write_replaced_pixels_raster"] = {"bool","false",
        "Write the DEM after pixel replacement.",
        "Without this the corrected DEM is still used by every later stage, but "
        "is not written on its own. Set it when the corrected DEM is the thing "
        "you want.",
        "Pixel replacement", "switch", "replace_pixels_fname"};

    string_default_map["replace_pixels_value_field"] = "value";
    help_map["replace_pixels_value_field"] = {"string","value",
        "Column holding the replacement value.",
        "Named field in replace_pixels_fname giving the new value for each "
        "listed pixel.",
        "Pixel replacement", "param", "replace_pixels_fname"};

    // --- Sink candidates ---
    //
    // The same detection lsdtt-internal-drainage runs, surfaced here so a
    // reviewer can see closed depressions while conditioning a DEM, before
    // deciding whether any of them should be treated as real sinks. This CLI
    // only WRITES the candidates; routing into them is lsdtt-internal-drainage's
    // job.
    bool_default_map["write_sink_candidates"] = false;
    help_map["write_sink_candidates"] = {"bool","false",
        "Write a point layer of candidate internally drained sinks.",
        "Finds closed depressions that never drain to the DEM edge and writes "
        "them for review in a GIS. Pass the reviewed file to "
        "lsdtt-internal-drainage to route flow into the accepted sinks. Output "
        "format follows write_vector_format.",
        "Sink candidates", "switch", ""};

    int_default_map["sink_candidate_min_pixels"] = 1;
    help_map["sink_candidate_min_pixels"] = {"int","1",
        "Minimum sink size in pixels.",
        "Closed depressions smaller than this are not reported. Raise it to "
        "suppress DEM noise on rough surfaces.",
        "Sink candidates", "param", "write_sink_candidates"};

    // --- Depression polygons ---
    bool_default_map["write_depression_polygons"] = false;
    help_map["write_depression_polygons"] = {"bool","false",
        "Write a vector polygon layer of topographic depressions, "
        "based on the fill diff raster.",
        "Each polygon carries depression_id, area_m2, and volume_m3 attributes. "
        "Output format is controlled by write_vector_format.",
        "Depression polygons", "switch", ""};

    int_default_map["depression_min_pixels"] = 1;
    help_map["depression_min_pixels"] = {"int","1",
        "Minimum depression size in pixels.",
        "Depressions with fewer pixels than this are removed before "
        "polygon tracing. Useful for suppressing single-cell fill artefacts.",
        "Depression polygons", "param", "write_depression_polygons"};

    float_default_map["depression_min_depth"] = 0.0f;
    help_map["depression_min_depth"] = {"float","0.0",
        "Minimum fill depth in metres for a cell to be considered a depression.",
        "Cells where fill_diff <= this value are treated as background. "
        "A small positive value (e.g. 0.01) suppresses near-zero fill artefacts.",
        "Depression polygons", "param", "write_depression_polygons"};

    string_default_map["legacy_depression_polygons_write_vector_format"] =
        lsdtt::inherit_write_vector_format;


    // --- Hillshade ---
    bool_default_map["write_hillshade_raster"] = false;
    help_map["write_hillshade_raster"] = {"bool","false",
        "Write a hillshade raster.",
        "Uses the Zevenbergen and Thorne (1987) algorithm. "
        "Computed from the filled DEM when write_fill_raster is also true.",
        "Hillshade", "switch", ""};

    float_default_map["hillshade_azimuth"] = 315.0f;
    help_map["hillshade_azimuth"] = {"float","315",
        "Sun azimuth in degrees clockwise from north.",
        "315 (northwest) is the cartographic convention.",
        "Hillshade", "param", "write_hillshade_raster"};

    float_default_map["hillshade_altitude"] = 45.0f;
    help_map["hillshade_altitude"] = {"float","45",
        "Sun altitude above the horizon in degrees.",
        "45 degrees is a typical default.",
        "Hillshade", "param", "write_hillshade_raster"};

    // --- Mosaic ---
    string_default_map["mosaic_file_list"] = "";
    help_map["mosaic_file_list"] = {"string","",
        "Path to a text file listing rasters to mosaic (one path per line).",
        "If set, the listed rasters are merged into a single output raster "
        "before any other preprocessing step. Rasters with a different CRS, "
        "cell size, or pixel-misaligned origin are skipped with a warning. "
        "Blend semantics: last valid pixel wins (later files take priority "
        "over earlier ones in overlap regions, matching gdal_merge.py default). "
        "The mosaic is written to write_path/write_fname_mosaic.tif and used "
        "as the input DEM for all subsequent pipeline steps. "
        "Lines starting with # and blank lines are ignored in the file list.",
        "Mosaic", "switch", ""};

    float_default_map["hillshade_z_factor"] = 1.0f;
    help_map["hillshade_z_factor"] = {"float","1.0",
        "Vertical exaggeration factor applied before computing hillshade.",
        "Values greater than 1 emphasise relief in flat terrain.",
        "Hillshade", "param", "write_hillshade_raster"};

    bool_default_map["write_multidirectional_hillshade_raster"] = false;
    help_map["write_multidirectional_hillshade_raster"] = {"bool","false",
        "Write a multidirectional hillshade raster.",
        "Computes hillshade from multiple azimuths and writes the mean "
        "illumination raster to reduce directional bias.",
        "Hillshade", "switch", ""};

    int_default_map["multidirectional_hillshade_n_directions"] = 8;
    help_map["multidirectional_hillshade_n_directions"] = {"int","8",
        "Number of evenly spaced azimuth directions used for multidirectional hillshade.",
        "Must be >= 2. Typical values are 8 to 16.",
        "Hillshade", "param", "write_multidirectional_hillshade_raster"};

    float_default_map["multidirectional_hillshade_start_azimuth"] = 0.0f;
    help_map["multidirectional_hillshade_start_azimuth"] = {"float","0.0",
        "Starting azimuth in degrees clockwise from north for multidirectional hillshade.",
        "Subsequent azimuths are sampled uniformly around 360 degrees.",
        "Hillshade", "param", "write_multidirectional_hillshade_raster"};

    // --- Vector clipping ---
    bool_default_map["clip_vector_data_to_raster"] = false;
    help_map["clip_vector_data_to_raster"] = {"bool","false",
        "Clip a vector file to the bounding box of the output raster.",
        "The vector CRS need not match the raster CRS - the raster bbox is "
        "reprojected into the vector's native CRS before clipping. "
        "Set vector_data_fname to the path of the input vector file.",
        "Vector clipping", "switch", ""};

    string_default_map["vector_data_fname"] = "";
    help_map["vector_data_fname"] = {"string","",
        "Path to the input vector file to clip.",
        "Any OGR-readable format is accepted (Shapefile, GeoPackage, "
        "FlatGeobuf, GeoJSON, etc.). Required when "
        "clip_vector_data_to_raster is true.",
        "Vector clipping", "param", "clip_vector_data_to_raster"};

    string_default_map["legacy_vector_clip_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // --- Clip raster to polygon ---
    bool_default_map["clip_raster_to_polygon"] = false;
    help_map["clip_raster_to_polygon"] = {"bool","false",
        "Crop the DEM to the minimum bounding box of a polygon.",
        "The DEM is cropped to the smallest axis-aligned rectangle that "
        "contains the polygon, in the raster's own CRS. "
        "Pixels outside the polygon are also set to nodata unless "
        "force_rectangular_raster is true. "
        "If the file contains more than one polygon, a warning is issued and "
        "only the first polygon is used. "
        "Set polygon_clip_fname to the path of the polygon file. "
        "Output suffix: _polygon_clipped.tif",
        "Clip raster to polygon", "switch", ""};

    string_default_map["polygon_clip_fname"] = "";
    help_map["polygon_clip_fname"] = {"string","",
        "Path to the polygon vector file used to clip the DEM.",
        "Any OGR-readable polygon format is accepted (GeoPackage, FlatGeobuf, "
        "Shapefile, GeoJSON, etc.). The first layer and first polygon are used. "
        "Required when clip_raster_to_polygon is true.",
        "Clip raster to polygon", "param", "clip_raster_to_polygon"};

    bool_default_map["force_rectangular_raster"] = false;
    help_map["force_rectangular_raster"] = {"bool","false",
        "When clipping to a polygon, keep all data inside the bounding box.",
        "If false (default), pixels that fall inside the bounding box but "
        "outside the polygon shape are set to nodata. "
        "If true, all pixels within the bounding box are kept regardless of "
        "whether they fall inside the polygon - the output is a plain "
        "rectangular crop. Only used when clip_raster_to_polygon is true.",
        "Clip raster to polygon", "param", "clip_raster_to_polygon"};

    string_default_map["polygon_clip_rule"] = "majority_coverage";
    help_map["polygon_clip_rule"] = {"string","majority_coverage",
        "Pixel inclusion rule when masking by the polygon shape.",
        "Only used when clip_raster_to_polygon is true and "
        "force_rectangular_raster is false. "
        "any_coverage: pixel is kept if the polygon touches any part of it. "
        "majority_coverage: pixel is kept if the polygon centre falls inside "
        "(approximately 50% coverage, default). "
        "full_coverage: pixel is kept only if the polygon covers it fully.",
        "Clip raster to polygon", "param", "clip_raster_to_polygon"};

    // --- Point sampling ---
    bool_default_map["sample_rasters_to_points"] = false;
    help_map["sample_rasters_to_points"] = {"bool","false",
        "Sample a list of rasters onto a point dataset.",
        "Reads a point dataset and a CSV file with columns raster,attribute. "
        "Each raster is sampled at each point and written as a new attribute "
        "in an output point file.",
        "Point sampling", "switch", ""};

    string_default_map["point_data_fname"] = "";
    help_map["point_data_fname"] = {"string","",
        "Path to the input point dataset.",
        "Recommended input is an OGR-readable point vector file such as "
        "FlatGeobuf, GeoPackage, Shapefile, or GeoJSON. "
        "CSV latitude/longitude support can also be added, but vector input is "
        "the cleanest first implementation because existing attributes are preserved.",
        "Point sampling", "param", "sample_rasters_to_points"};

    string_default_map["point_data_layer_name"] = "";
    help_map["point_data_layer_name"] = {"string","",
        "Optional layer name for the input point dataset.",
        "If empty, the first layer is used, or the only layer if the file has one.",
        "Point sampling", "param", "sample_rasters_to_points"};

    string_default_map["raster_attribute_fname"] = "";
    help_map["raster_attribute_fname"] = {"string","",
        "CSV file listing rasters to sample, with columns raster,attribute.",
        "The raster column gives the raster filename or path. "
        "The attribute column gives the name of the output attribute to append "
        "to the point file.",
        "Point sampling", "param", "sample_rasters_to_points"};

    string_default_map["legacy_sampled_points_write_vector_format"] =
        lsdtt::inherit_write_vector_format;




    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

    // -------------------------------------------------------------------------
    // Load parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(cli_args.path, cli_args.fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parse - fills maps from the parameter file, falling back to defaults
    // -------------------------------------------------------------------------
    pp.set_deprecated_aliases({
        {"apply_binary_mask", "run_binary_mask"},
        {"raster_attribute_csv", "raster_attribute_fname"},
        {"replace_pixels_from_csv", "replace_pixels_fname"},
        {"minimum_elevation", "min_elevation"} });
    pp.set_deprecated_compatibility_aliases({
        {"depression_polygons_output_format",
         {"legacy_depression_polygons_write_vector_format", "write_vector_format"}},
        {"vector_clip_output_format",
         {"legacy_vector_clip_write_vector_format", "write_vector_format"}},
        {"sampled_points_output_format",
         {"legacy_sampled_points_write_vector_format", "write_vector_format"}} });
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (bf["force_resampling"] && ff["forced_grid_spacing"] <= 0.0f)
    {
        std::cerr << "FATAL: force_resampling is true but forced_grid_spacing "
                  << "is not positive.\n";
        return 1;
    }

    if (bf["force_resampling"])
    {
        std::cout << "Force resampling enabled: using forced_grid_spacing = "
                  << ff["forced_grid_spacing"] << " m\n";
    }

    if (bf["force_reproject"])
    {
        std::cout << "Force reprojection enabled: projected inputs will also "
                  << "be warped to the recommended analysis CRS.\n";
    }



    // -------------------------------------------------------------------------
    // Derived / cross-parameter logic
    //
    // These adjustments are the reason parameters are not const.
    // If breach_before_fill is set, we always need to run breaching even if
    // the user did not explicitly request the breached raster output.
    // -------------------------------------------------------------------------
    if (bf["breach_before_fill"])
    {
        // Ensure the breaching step runs even if no breached-raster output
        // was requested - breach_before_fill implies we want a breached dem
        // to feed into the fill step.
        bf["write_breached_raster"] = bf["write_breached_raster"] || false;
    }

    if (bf["write_fsm_water_depth_raster"])
    {
        bf["write_fsm_raster"] = true;
    }

    // Depression polygons require fill_diff to have been computed in memory.
    if (bf["write_depression_polygons"])
    {
        bf["write_fill_diff_raster"] = true;
    }
 

    int n_threads = ifs["n_threads"];

    // Flat resolution for pit filling. Checked here, before any raster work,
    // so a bad value or combination fails without loading anything.
    // inline_gradient keeps the fill exactly as it was; the barnes2014b
    // options fill to the spill level and then resolve the flats (Barnes et
    // al. 2014b), so they only make sense for the barnes2014 and vincent1993
    // fills.
    lsdtt::flat_resolution_choice flat_choice =
        lsdtt::flat_resolution_choice::inline_gradient;
    try
    {
        flat_choice = lsdtt::parse_flat_resolution_choice(
            sf["fill_flat_resolution_choice"]);
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    if ((bf["write_fill_raster"] || bf["write_fill_diff_raster"]) &&
        flat_choice != lsdtt::flat_resolution_choice::inline_gradient)
    {
        if (sf["fill_algorithm_choice"] == "barnesetal2021")
        {
            std::cerr << "FATAL: fill_flat_resolution_choice = "
                      << lsdtt::flat_resolution_choice_name(flat_choice)
                      << " applies to the barnes2014 and vincent1993 fills, "
                         "not barnesetal2021\n";
            return 1;
        }

        try
        {
            lsdtt::check_fill_flat_resolution_min_slope(
                flat_choice, ff["fill_min_slope"]);
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
    }

    // Least-cost breaching options. Built and validated here, before any
    // raster work, so a bad value fails with the library's message and a
    // non-zero exit rather than being caught as a breach failure later.
    lsdtt::breach_least_cost_options breach_lc_options;
    const bool run_breach = bf["breach_before_fill"] || bf["write_breached_raster"] ||
                            bf["write_carve_diff_raster"];
    if (run_breach && sf["breach_algorithm_choice"] == "least_cost")
    {
        try
        {
            breach_lc_options.max_cost       =
                (ff["breach_lc_max_cost"] >= std::numeric_limits<float>::max())
                    ? std::numeric_limits<double>::infinity()
                    : static_cast<double>(ff["breach_lc_max_cost"]);
            breach_lc_options.minimize_dist  = bf["breach_lc_minimize_dist"];
            breach_lc_options.flat_increment = static_cast<double>(ff["breach_lc_flat_increment"]);
            breach_lc_options.fill_remainder = bf["fill_remainder_after_breach"];
            breach_lc_options.fill_min_slope = ff["fill_min_slope"];
            breach_lc_options.fill_flat_resolution = flat_choice;
            breach_lc_options.edge_policy    = lsdtt::parse_breach_lc_edge_policy(
                sf["breach_lc_edge_policy"], "breach_lc_edge_policy");
            // max_dist defaults by edge policy (1000 outlet, to bridge large
            // water bodies; 100 whitebox, Whitebox's own) unless the file
            // sets it. See default_breach_lc_max_dist.
            const bool max_dist_given = pp.parameter_supplied("breach_lc_max_dist");
            breach_lc_options.max_dist = max_dist_given
                ? ifs["breach_lc_max_dist"]
                : lsdtt::default_breach_lc_max_dist(breach_lc_options.edge_policy);
            lsdtt::validate_breach_least_cost_options(breach_lc_options);
            std::cout << "Least-cost breaching: breach_lc_max_dist = "
                      << breach_lc_options.max_dist
                      << (max_dist_given ? " (from the parameter file)"
                                         : " (default for this edge policy)")
                      << ", breach_lc_edge_policy = "
                      << lsdtt::breach_lc_edge_policy_name(breach_lc_options.edge_policy)
                      << "\n";
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
    }

    // Resolve output format
    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt_str = sf["write_output_format"];
        if (fmt_str == "gtiff")
        {
            out_fmt = lsdtt::raster_output_format::gtiff;
        }
        else if (fmt_str == "same_as_input")
        {
            out_fmt = lsdtt::raster_output_format::same_as_input;
        }
        else if (fmt_str != "cog")
        {
            std::cerr << "Warning: unknown write_output_format '"
                      << fmt_str << "' - defaulting to cog.\n";
        }
        std::cout << "Output format: " << fmt_str << "\n";
    }

    // -------------------------------------------------------------------------
    // Step 0 - Mosaic (optional)
    //
    // If mosaic_file_list is set, read the listed rasters, mosaic them into
    // a single raster, write it, and use it as the DEM for all subsequent
    // pipeline steps.  When not set, the normal read_fname is used.
    //
    // read_fname is ignored when mosaic_file_list is set - warn the user if
    // both are present so there is no silent data-loss confusion.
    // write_fname must be non-empty (either explicitly set or derived from
    // read_fname) so the mosaic output path is well-defined.
    // -------------------------------------------------------------------------
    std::string active_path;

    {
        const std::string& list_path = sf["mosaic_file_list"];
        if (!list_path.empty())
        {
            // Warn if read_fname was also set - it will be ignored.
            if (!pp.get_read_fname().empty())
            {
                std::cerr << "Warning: both mosaic_file_list and read_fname are "
                          << "set. read_fname ('" << pp.get_read_fname()
                          << "') is ignored when mosaic_file_list is used.\n";
            }

            // write_fname must be non-empty for the output path to be usable.
            if (pp.get_write_fname().empty())
            {
                std::cerr << "FATAL: mosaic_file_list is set but write_fname is "
                          << "empty. Set write_fname in the parameter file to "
                          << "define the output filename stem.\n";
                return 1;
            }

            std::cout << "\n--- Mosaic ---\n";
            std::cout << "Reading file list: " << list_path << "\n";

            std::vector<std::string> mosaic_paths;
            {
                std::ifstream list_file(list_path);
                if (!list_file.is_open())
                {
                    std::cerr << "FATAL: cannot open mosaic_file_list '"
                              << list_path << "'.\n";
                    return 1;
                }
                std::string line;
                while (std::getline(list_file, line))
                {
                    const auto first = line.find_first_not_of(" \t\r\n");
                    if (first == std::string::npos) continue;
                    const auto last = line.find_last_not_of(" \t\r\n");
                    line = line.substr(first, last - first + 1);
                    if (line.empty() || line[0] == '#') continue;
                    mosaic_paths.push_back(line);
                }
            }

            if (mosaic_paths.empty())
            {
                std::cerr << "FATAL: mosaic_file_list contains no file paths.\n";
                return 1;
            }

            std::cout << "  " << mosaic_paths.size() << " path(s) found.\n";

            lsdtt::lsdtt_raster_cpu mosaic_raster;
            try
            {
                mosaic_raster = lsdtt::raster_mosaic(mosaic_paths);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: mosaic failed: " << e.what() << "\n";
                return 1;
            }

            const std::string mosaic_out =
                pp.get_write_path() + pp.get_write_fname() + "_mosaic.tif";
            std::cout << "Writing mosaic: " << mosaic_out << "\n";

            try
            {
                lsdtt::lsdtt_data_loader wl;
                wl.write_full(
                    mosaic_raster.get_data(),
                    mosaic_raster.get_meta(),
                    mosaic_out,
                    out_fmt
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing mosaic failed: " << e.what() << "\n";
                return 1;
            }

            std::cout << "Mosaic written.\n";
            active_path = mosaic_out;
        }
        else
        {
            active_path = pp.get_full_read_path();
        }
    }

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    std::cout << "\nOpening DEM: " << active_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(active_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    // -------------------------------------------------------------------------
    // Step 1 - Projection check and reprojection
    //
    // If the raster is geographic or Web Mercator, reproject it to the
    // appropriate projected CRS.  All subsequent steps use the reprojected
    // raster.  The reprojected file is written to the output directory with
    // a projection suffix (_utm, _polar, _lambert).
    // -------------------------------------------------------------------------
    {
        std::cout << "\n--- Projection ---\n";

        lsdtt::lsdtt_reproject_options proj_opts;
        proj_opts.enforce_projection       = bf["enforce_projection"];
        proj_opts.force_reproject          = bf["force_reproject"];
        proj_opts.resampling_method        = sf["resampling_method"];
        proj_opts.grid_spacing             = static_cast<double>(ff["grid_spacing"]);
        proj_opts.force_resampling         = bf["force_resampling"];
        proj_opts.forced_grid_spacing      = static_cast<double>(ff["forced_grid_spacing"]);
        proj_opts.allow_lambert_projection = bf["allow_lambert_projection"];
        proj_opts.n_threads                = n_threads;
        proj_opts.warp_mem_mb              = ifs["warp_mem_mb"];



        
        lsdtt::lsdtt_reproject_result proj_result;
        try
        {
            proj_result = lsdtt::lsdtt_reproject_if_needed(
                active_path,
                pp.get_write_path(),
                pp.get_write_fname(),
                loader.get_meta(),
                proj_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "\nFATAL: " << e.what() << "\n";
            return 1;
        }

        if (proj_result.reprojected)
        {
            active_path = proj_result.out_path;
            loader.close_raster();
            try
            {
                loader.open_raster(active_path);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error reopening reprojected DEM: "
                          << e.what() << "\n";
                return 1;
            }
            std::cout << "Downstream pipeline will use: "
                      << active_path << "\n";
            loader.print_raster_info();
        }

        // Final CRS check - must be projected in metres at this point
        try
        {
            loader.get_meta().check_projected_metres();
        }
        catch (const std::exception& e)
        {
            std::cerr << "\nFATAL: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Load into CPU raster object
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);

    const std::string write_stem = pp.get_write_fname();

    // -------------------------------------------------------------------------
    // Step 1b - Vertical units
    //
    // Before ANY masking or analysis, so that min_elevation and every
    // derived quantity below are in metres. Opt-in only; see the header of
    // lsdtt_raster_cpu.hpp for why this is never inferred.
    // -------------------------------------------------------------------------
    const lsdtt::vertical_unit vunit =
        lsdtt::parse_vertical_unit(sf["vertical_units"]);
    if (vunit != lsdtt::vertical_unit::metres)
    {
        std::cout << "\nConverting elevations from " << sf["vertical_units"]
                  << " to metres.\n";
        dem = lsdtt::convert_vertical_units_to_metres(dem, vunit);
    }
    else
    {
        // Only worth warning when the user has NOT already told us the unit.
        const auto plaus = lsdtt::warn_if_vertically_implausible(dem, std::cout);

        // The warning can only say "this looks impossible". 3DEP can say which
        // unit it actually is, so it runs only when there is a question to
        // answer -- never as a routine network call on a healthy DEM.
        if (plaus.implausible && bf["check_vertical_units_against_usgs"])
        {
            std::cout << "Checking the vertical unit against USGS 3DEP...\n";
            const lsdtt::usgs_unit_check chk =
                lsdtt::check_vertical_units_against_usgs(dem);
            std::cout << "  " << chk.message << "\n\n";

            if (chk.verdict == lsdtt::vertical_unit_verdict::feet)
            {
                // Deliberately NOT applied here. The user asked for this DEM as
                // it stands; converting behind their back would make the output
                // disagree with the input for a reason nothing in the file
                // records. They are told exactly what to set.
                std::cout << "  Rerun with 'vertical_units: feet' to write a "
                             "corrected raster.\n\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 2 - Threshold masking
    // -------------------------------------------------------------------------
    float min_elev = ff["min_elevation"];
    float max_elev = ff["maximum_elevation"];

    // Count valid pixels before masking so we can detect any change.
    const auto count_valid = [](const lsdtt::lsdtt_raster_cpu& r) -> std::size_t {
        std::size_t n = 0;
        for (uint8_t b : r.get_mask()) n += b;
        return n;
    };
    const std::size_t valid_before_masking = count_valid(dem);

    if (min_elev >= 0.0f)
    {
        std::cout << "\nMasking pixels at or below " << min_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, min_elev, true);
    }
    if (max_elev < 30000.0f)
    {
        std::cout << "Masking pixels above " << max_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, max_elev, false);
    }

    if (count_valid(dem) < valid_before_masking)
    {
        const std::string o = pp.make_output_path( "_remove_seas.tif");
        std::cout << "Threshold masking changed "
                  << (valid_before_masking - count_valid(dem))
                  << " pixel(s) - writing: " << o << "\n";
        loader.write_full(dem.get_data(), dem.get_meta(), o, out_fmt);
    }

    // -------------------------------------------------------------------------
    // Step 2b - Binary mask application
    // -------------------------------------------------------------------------
    if (bf["run_binary_mask"])
    {
        std::cout << "\n--- Binary mask ---\n";

        const std::string& mask_path = sf["binary_mask_fname"];
        if (mask_path.empty())
        {
            std::cerr << "Warning: run_binary_mask is true but "
                         "binary_mask_fname is empty - skipping.\n";
        }
        else
        {
            std::cout << "Mask file: " << mask_path << "\n";
            try
            {
                lsdtt::lsdtt_data_loader mask_loader;
                mask_loader.open_raster(mask_path);
                const lsdtt::lsdtt_raster_cpu mask_raster(mask_loader);

                dem = lsdtt::algebra::raster_mask_by_binary_keep(dem, mask_raster);

                const std::size_t valid_after = count_valid(dem);
                std::cout << "  Valid pixels after masking: " << valid_after << "\n";

                const std::string o = pp.make_output_path( "_binary_masked.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(dem.get_data(), dem.get_meta(), o, out_fmt);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: binary mask application failed: "
                          << e.what() << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 2c - Crop to polygon bounding box (and optionally mask by shape)
    // -------------------------------------------------------------------------
    if (bf["clip_raster_to_polygon"])
    {
        std::cout << "\n--- Clip raster to polygon ---\n";

        const std::string& poly_path = sf["polygon_clip_fname"];
        if (poly_path.empty())
        {
            std::cerr << "Warning: clip_raster_to_polygon is true but "
                         "polygon_clip_fname is empty - skipping.\n";
        }
        else
        {
            std::cout << "Polygon file: " << poly_path << "\n";
            try
            {
                // Warn if the file contains multiple polygons.
                const int n_features =
                    lsdtt::lsdtt_data_loader::count_features_in_first_layer(
                        poly_path);
                if (n_features > 1)
                {
                    std::cerr << "Warning: polygon file contains "
                              << n_features
                              << " polygon(s) - only the first will be used.\n";
                }
                else
                {
                    std::cout << "  " << n_features << " polygon(s) found.\n";
                }

                // Get the polygon bbox in the raster's CRS.
                const lsdtt::polygon_bbox bbox =
                    lsdtt::get_first_polygon_bbox_in_raster_crs(
                        poly_path, dem.get_meta().crs_wkt);

                std::cout << "  Polygon bbox (raster CRS):\n"
                          << "    xmin=" << bbox.xmin << "  ymin=" << bbox.ymin
                          << "  xmax=" << bbox.xmax << "  ymax=" << bbox.ymax
                          << "\n";

                // Crop the DEM to the polygon bbox.
                dem = lsdtt::clip_raster_to_bbox(
                          dem, bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax);

                // Optionally mask pixels outside the polygon shape.
                if (!bf["force_rectangular_raster"])
                {
                    lsdtt::polygon_rasterise_rule poly_rule =
                        lsdtt::polygon_rasterise_rule::majority_coverage;
                    {
                        const std::string& rs = sf["polygon_clip_rule"];
                        if (rs == "any_coverage")
                            poly_rule =
                                lsdtt::polygon_rasterise_rule::any_coverage;
                        else if (rs == "full_coverage")
                            poly_rule =
                                lsdtt::polygon_rasterise_rule::full_coverage;
                        else if (rs != "majority_coverage")
                            std::cerr << "Warning: unknown polygon_clip_rule '"
                                      << rs
                                      << "' - using majority_coverage.\n";
                    }

                    const lsdtt::lsdtt_raster_cpu poly_mask =
                        lsdtt::rasterise_polygons(dem, poly_path,
                                                  poly_rule, 1);

                    dem = lsdtt::algebra::raster_mask_by_binary_keep(
                              dem, poly_mask);
                }

                const std::string o = pp.make_output_path( "_polygon_clipped.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(dem.get_data(), dem.get_meta(), o, out_fmt);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: polygon clip failed: "
                          << e.what() << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 2c - Pixel replacement
    //
    // Runs before breaching, filling and nodata filling: this is a correction to
    // the source data, and anything downstream should see the corrected DEM.
    // -------------------------------------------------------------------------
    if (!sf["replace_pixels_fname"].empty())
    {
        std::cout << "\n--- Pixel replacement ---\n";

        const std::filesystem::path given(sf["replace_pixels_fname"]);
        const std::string points_path =
            given.is_absolute() ? sf["replace_pixels_fname"]
                                : pp.get_read_path() + sf["replace_pixels_fname"];

        const std::string value_field = sf["replace_pixels_value_field"];

        const lsdtt::lsdtt_point_attribute_read_result points =
            lsdtt::load_point_layer_with_attributes(
                points_path, dem.get_meta(), std::vector<std::string>{});

        lsdtt::lsdtt_array2d_f&  data = dem.get_data_array();
        lsdtt::lsdtt_array2d_u8& mask = dem.get_mask_array();

        const double nodata = dem.get_meta().no_data_value;

        std::size_t replaced = 0, blanked = 0;
        std::size_t outside = 0, no_value = 0;

        for (const lsdtt::lsdtt_input_point_with_attributes& record :
             points.records)
        {
            if (record.point.status != lsdtt::lsdtt_point_input_status::success)
            {
                ++outside;
                continue;
            }

            const auto it = record.attributes.find(value_field);
            double new_value = 0.0;
            if (it == record.attributes.end() ||
                !it->second.as_double(new_value))
            {
                ++no_value;
                continue;
            }

            const int r = record.point.row;
            const int c = record.point.col;

            // A value equal to the raster's nodata must also clear the validity
            // mask, or the cell would carry the nodata number while still being
            // treated as data by everything downstream.
            if (new_value == nodata)
            {
                data(r, c) = static_cast<float>(nodata);
                mask(r, c) = 0u;
                ++blanked;
            }
            else
            {
                data(r, c) = static_cast<float>(new_value);
                mask(r, c) = 1u;
                ++replaced;
            }
        }

        std::cout << "Replaced " << replaced << " pixel(s); blanked "
                  << blanked << " to nodata.\n";

        if (outside > 0)
            std::cerr << "Warning: " << outside
                      << " listed pixel(s) fall outside the raster and were "
                         "skipped.\n";
        if (no_value > 0)
            std::cerr << "Warning: " << no_value
                      << " listed pixel(s) had no readable '" << value_field
                      << "' value and were skipped.\n";

        if (bf["write_replaced_pixels_raster"])
        {
            const std::string o = pp.make_output_path("_pixels_replaced.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(dem.get_data(), dem.get_meta(), o, out_fmt);
        }
    }

    // -------------------------------------------------------------------------
    // Step 3 - Interior nodata filling
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu nodata_filled_dem = dem;
    bool did_nodata_fill = false;

    if (bf["write_filled_nodata_raster"])
    {
        std::cout << "\n--- Interior nodata filling ---\n";

        lsdtt::nodata_fill_method nfm = lsdtt::nodata_fill_method::laplacian;
        if (sf["nodata_fill_method_choice"] == "iterative_weighted_average")
        {
            nfm = lsdtt::nodata_fill_method::iterative_weighted_average;
        }

        try
        {
            nodata_filled_dem = lsdtt::fill_interior_nodata(
                dem, nfm, ifs["nodata_fill_max_iter"], 1e-6f, n_threads
            );

            std::string o = pp.make_output_path( "_nodata_filled.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(nodata_filled_dem.get_data(),
                              nodata_filled_dem.get_meta(), o, out_fmt);
            did_nodata_fill = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: nodata filling failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 4 - Clip to maximum valid rectangle
    // -------------------------------------------------------------------------
    if (bf["write_clipped_rectangle_raster"])
    {
        std::cout << "\n--- Clip to maximum rectangle ---\n";
        const lsdtt::lsdtt_raster_cpu& src = did_nodata_fill
                                           ? nodata_filled_dem
                                           : dem;
        try
        {
            lsdtt::lsdtt_raster_cpu clipped =
                lsdtt::clip_to_maximum_rectangle(src);
            std::string o = pp.make_output_path( "_clipped_rect.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(clipped.get_data(), clipped.get_meta(),
                              o, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: rectangle clipping failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 4b - Edge-preserving smoothing (conditioning)
    //
    // Before breaching and filling, because it is conditioning: everything
    // downstream should see one surface, and which one should be visible in the
    // file name rather than inferred from a parameter file.
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu working_dem = dem;

    {
        const std::string& smooth_choice = sf["smooth_algorithm_choice"];
        if (smooth_choice != "none" && smooth_choice != "perona_malik" &&
            smooth_choice != "coherence_diffusion")
        {
            std::cerr << "FATAL: unknown smooth_algorithm_choice '"
                      << smooth_choice
                      << "'. Valid: none, perona_malik, coherence_diffusion.\n";
            return 1;
        }

        if (smooth_choice == "perona_malik")
        {
            std::cout << "\n--- Perona-Malik conditioning ---\n";

            lsdtt::perona_malik_options pmo;
            pmo.n_iterations = ifs["perona_malik_iterations"];
            pmo.lambda = static_cast<double>(ff["perona_malik_lambda"]);
            pmo.k_gradient_quantile =
                static_cast<double>(ff["perona_malik_k_quantile"]);
            pmo.k_metres_per_cell =
                static_cast<double>(ff["perona_malik_k_metres_per_cell"]);
            pmo.n_threads = ifs["n_threads"];

            const std::string& cond = sf["perona_malik_conductance"];
            if (cond == "exponential")
                pmo.conductance = lsdtt::pm_conductance_choice::exponential;
            else if (cond == "quadratic")
                pmo.conductance = lsdtt::pm_conductance_choice::quadratic;
            else
            {
                std::cerr << "FATAL: unknown perona_malik_conductance '"
                          << cond << "'. Valid: exponential, quadratic.\n";
                return 1;
            }

            if (bf["write_stencil_difference_raster"])
            {
                try
                {
                    const lsdtt::lsdtt_raster_cpu g =
                        lsdtt::stencil_difference_magnitude(dem);
                    const std::string o =
                        pp.make_output_path("_stencil_difference.tif");
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(g.get_data(), g.get_meta(), o, out_fmt);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: neighbour-difference raster failed: "
                              << e.what() << "\n";
                }
            }

            const auto pm_t0 = std::chrono::steady_clock::now();
            lsdtt::perona_malik_result pmr;
            try
            {
                pmr = lsdtt::perona_malik_diffusion(dem, pmo);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: Perona-Malik failed: " << e.what() << "\n";
                return 1;
            }
            const auto pm_t1 = std::chrono::steady_clock::now();

            std::cout << "  " << pmr.n_iterations << " iteration(s), lambda "
                      << pmo.lambda << ", c = " << cond << ", K = "
                      << pmr.k_metres_per_cell << " m/cell";
            if (pmo.k_gradient_quantile > 0.0)
                std::cout << " (quantile " << pmo.k_gradient_quantile << ")";
            std::cout << "\n";
            std::cout << "  Compute time: "
                      << std::chrono::duration<double>(pm_t1 - pm_t0).count()
                      << " s\n";

            const std::string o = pp.make_output_path("_smoothed.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(pmr.smoothed.get_data(),
                              pmr.smoothed.get_meta(), o, out_fmt);

            working_dem = pmr.smoothed;
        }

        if (smooth_choice == "coherence_diffusion")
        {
            std::cout << "\n--- Coherence-enhancing diffusion ---\n";

            const double cell = dem.get_meta().cell_size;
            double deriv_r =
                static_cast<double>(ff["coherence_derivative_radius_metres"]);
            double integ_r =
                static_cast<double>(ff["coherence_integration_radius_metres"]);
            // Stated in cells, not metres: a radius tuned on 1 m lidar is
            // sub-pixel on 30 m COP30, and the structure tensor then has
            // nothing to integrate over.
            if (!(deriv_r > 0.0)) deriv_r = cell;
            if (!(integ_r > 0.0)) integ_r = 3.0 * cell;

            const int    iters = ifs["coherence_iterations"];
            const double tau   = static_cast<double>(ff["coherence_time_step"]);
            const double contr = static_cast<double>(ff["coherence_contrast"]);
            const double mind  =
                static_cast<double>(ff["coherence_min_diffusivity"]);

            std::cout << "  " << iters << " iteration(s), tau " << tau
                      << ", derivative radius " << deriv_r
                      << " m, integration radius " << integ_r
                      << " m, C = " << contr
                      << ", min diffusivity " << mind << "\n";

            const auto ce_t0 = std::chrono::steady_clock::now();
            lsdtt::lsdtt_raster_cpu diffused;
            try
            {
                diffused = lsdtt::coherence_enhancing_diffusion(
                    dem, iters, tau, deriv_r, integ_r, contr, mind,
                    ifs["n_threads"]);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: coherence-enhancing diffusion failed: "
                          << e.what() << "\n";
                return 1;
            }
            const auto ce_t1 = std::chrono::steady_clock::now();
            std::cout << "  Compute time: "
                      << std::chrono::duration<double>(ce_t1 - ce_t0).count()
                      << " s\n";

            const std::string o = pp.make_output_path("_smoothed.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(diffused.get_data(), diffused.get_meta(),
                              o, out_fmt);

            working_dem = diffused;
        }
    }

    // -------------------------------------------------------------------------
    // Step 5 - Depression breaching
    // -------------------------------------------------------------------------

    if (bf["breach_before_fill"] || bf["write_breached_raster"] ||
        bf["write_carve_diff_raster"])
    {
        std::cout << "\n--- Depression breaching ---\n";
        std::string algo = sf["breach_algorithm_choice"];

        // See the note on the pit-filling timer below: this separates the
        // breach itself from process startup and I/O so a benchmark can quote
        // it against an in-process reference timing.
        const auto breach_compute_t0 = std::chrono::steady_clock::now();

        try
        {
            lsdtt::lsdtt_raster_cpu breached;

            if (algo == "least_cost")
            {
                std::cout << "Algorithm: Lindsay (2020) least-cost "
                             "breaching\n";
                breached = lsdtt::breach_depressions_least_cost(
                    working_dem, breach_lc_options);
            }
            else
            {
                std::cout << "Algorithm: Lindsay (2016) efficient "
                             "breaching\n";
                breached = lsdtt::breach_depressions(working_dem,
                                                     n_threads);
            }

            const auto breach_compute_t1 = std::chrono::steady_clock::now();
            std::cout << "  Compute time: "
                      << std::chrono::duration<double>(
                             breach_compute_t1 - breach_compute_t0).count()
                      << " s\n";

            if (bf["write_breached_raster"])
            {
                std::string o = pp.make_output_path( "_breached.tif");
                std::cout << "Writing: " << o << "\n";
                const auto breach_write_t0 = std::chrono::steady_clock::now();
                loader.write_full(breached.get_data(),
                                  breached.get_meta(), o, out_fmt);
                const auto breach_write_t1 = std::chrono::steady_clock::now();
                std::cout << "  Write time: "
                          << std::chrono::duration<double>(
                                 breach_write_t1 - breach_write_t0).count()
                          << " s\n";
            }

            if (bf["write_carve_diff_raster"])
            {
                std::cout << "Computing carve difference raster "
                             "(breached - original)\n";
                try
                {
                    lsdtt::lsdtt_raster_cpu carve_diff =
                        lsdtt::algebra::raster_subtract(breached, working_dem);
                    std::string o = pp.make_output_path( "_carve_diff.tif");
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(carve_diff.get_data(),
                                      carve_diff.get_meta(), o, out_fmt);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: carve diff raster failed: "
                              << e.what() << "\n";
                }
            }

            if (bf["breach_before_fill"])
            {
                working_dem = breached;
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: breaching failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 - Pit filling
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu filled_dem = working_dem;
    bool did_fill = false;

    lsdtt::lsdtt_raster_cpu fill_diff_raster;
    bool have_fill_diff = false;

    if (bf["write_fill_raster"] || bf["write_fill_diff_raster"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::string algo = sf["fill_algorithm_choice"];
        const std::string& queue_name = sf["priority_flood_queue"];
        const lsdtt::priority_flood_queue_backend queue_backend =
            queue_name == "binary" ? lsdtt::priority_flood_queue_backend::binary_heap :
            queue_name == "radix"  ? lsdtt::priority_flood_queue_backend::radix_hybrid :
            throw std::invalid_argument(
                "priority_flood_queue must be 'binary' or 'radix'");

        // Under a barnes2014b option the fill runs to the spill level by
        // reconstruction (identical to barnes2014 at min_slope 0) and the
        // flats are then resolved, so the spill queue is never used and
        // barnes2014b_float_step has no use for fill_min_slope.
        const bool barnes2014b_flats =
            flat_choice != lsdtt::flat_resolution_choice::inline_gradient;
        if (barnes2014b_flats && algo != "barnesetal2021")
        {
            if (queue_name != "binary")
            {
                std::cout << "NOTE: priority_flood_queue is ignored by "
                             "fill_flat_resolution_choice = "
                          << lsdtt::flat_resolution_choice_name(flat_choice)
                          << ".\n";
            }
            if (flat_choice ==
                    lsdtt::flat_resolution_choice::barnes2014b_float_step &&
                ff["fill_min_slope"] != 0.0f)
            {
                std::cout << "NOTE: fill_min_slope (" << ff["fill_min_slope"]
                          << ") is ignored by fill_flat_resolution_choice = "
                             "barnes2014b_float_step; flats are raised by "
                             "float32 steps.\n";
            }
        }

        // Fill to the spill level, resolve the flats (Barnes et al. 2014b)
        // and print what the resolver did. Used by barnes2014 and vincent1993
        // under a barnes2014b option.
        auto fill_and_resolve_flats = [&]() -> lsdtt::lsdtt_raster_cpu
        {
            return lsdtt::fill_pits_reporting_flat_resolution(
                working_dem, ff["fill_min_slope"], flat_choice, n_threads);
        };

        // Timed so the cost of the fill itself can be separated from process
        // startup, GDAL raster I/O and every other stage this CLI runs. The
        // "Compute time:"/"Write time:" wording matches lsdtt-surface-metrics
        // and is what tests/scripts/run_benchmark_suite.py parses; without it
        // the only figure available to a benchmark is whole-process wall clock,
        // which on a small raster is mostly process startup and is not
        // comparable with an in-process reference timing.
        const auto fill_compute_t0 = std::chrono::steady_clock::now();

        try
        {
            if (algo == "barnes2014")
            {
                std::cout << "Algorithm: Barnes et al. (2014) "
                             "Priority-Flood\n";

                if (barnes2014b_flats)
                {
                    filled_dem = fill_and_resolve_flats();
                }
                else
                {
                    filled_dem = lsdtt::fill_pits_priority_flood(
                        working_dem,
                        ff["fill_min_slope"],
                        n_threads,
                        queue_backend
                    );
                }
            }
            else if (algo == "vincent1993")
            {
                // Reconstruction fills to exactly the spill level and has no
                // mechanism for a gradient across flats. Refusing is better
                // than silently ignoring fill_min_slope and handing back a
                // surface with undrainable flats that the D8 router will then
                // reject. A barnes2014b option supplies that mechanism.
                if (!barnes2014b_flats && ff["fill_min_slope"] > 0.0f)
                {
                    std::cerr << "FATAL: fill_algorithm_choice = vincent1993 "
                                 "cannot impose a gradient across flats, but "
                                 "fill_min_slope is "
                              << ff["fill_min_slope"] << ".\n"
                              << "Set fill_min_slope: 0 to use it, or use "
                                 "barnes2014 (the default), which is what "
                                 "imposes the gradient, or set "
                                 "fill_flat_resolution_choice: "
                                 "barnes2014b_gradient.\n";
                    return 1;
                }

                std::cout << "Algorithm: Vincent (1993) grayscale morphological "
                             "reconstruction\n";

                if (barnes2014b_flats)
                {
                    filled_dem = fill_and_resolve_flats();
                }
                else
                {
                    filled_dem = lsdtt::fill_pits_reconstruction(
                        working_dem,
                        n_threads
                    );
                }
            }
            else if (algo == "barnesetal2021"){
                std::cout << "Algorithm: Barnes et al. (2021), "
                             "https://doi.org/10.5194/esurf-9-105-2021: "
                             "Fill-Spill-Merge\n";

                lsdtt::depression_hierarchy_options dh_opts;
                dh_opts.boundary_conditions = sf["fsm_boundary_conditions"];
    
                std::vector<lsdtt::dep_label_t> fsm_label;
                std::vector<uint8_t>            fsm_flowdir;
    
                std::cout << "I will first build the depression hierarchy "
                             "following Barnes et al. (2020), "
                             "https://doi.org/10.5194/esurf-8-431-2020, using "
                             "boundary conditions: " 
                             << sf["fsm_boundary_conditions"] << "\n";
    
                lsdtt::lsdtt_depression_hierarchy dh = lsdtt::build_depression_hierarchy(
                    working_dem, fsm_label, fsm_flowdir, dh_opts
                );
                
                std::cout << "I am now going to Fill-Spill-Merge this depression "
                            "hierarchy using a uniform rainfall rate of: "
                             << ff["fill_spill_merge_rainfall_rate"] << "\n";
 
                filled_dem = lsdtt::fsm_filled_dem(
                    working_dem, fsm_label, fsm_flowdir, dh,
                    ff["fill_spill_merge_rainfall_rate"]
                );

            }
            else {
                std::cout << "Warning: no correct fill_algorithm choice supplied, "
                          << "choose one of: \n"
                          << "barnes2014, vincent1993, barnesetal2021 \n"
                          << "I will now simply take the non-filled dem as filled dem. \n";

                filled_dem = working_dem;
            }

            const auto fill_compute_t1 = std::chrono::steady_clock::now();
            std::cout << "  Compute time: "
                      << std::chrono::duration<double>(
                             fill_compute_t1 - fill_compute_t0).count()
                      << " s\n";

            if (bf["write_fill_raster"])
            {
                std::string o = pp.make_output_path( "_fill.tif");
                std::cout << "Writing: " << o << "\n";
                const auto fill_write_t0 = std::chrono::steady_clock::now();
                loader.write_full(filled_dem.get_data(),
                                  filled_dem.get_meta(), o, out_fmt);
                const auto fill_write_t1 = std::chrono::steady_clock::now();
                std::cout << "  Write time: "
                          << std::chrono::duration<double>(
                                 fill_write_t1 - fill_write_t0).count()
                          << " s\n";
            }

            if (bf["write_fill_diff_raster"])
            {
                std::cout << "Computing fill difference raster "
                             "(filled - original)\n";
                try
                {
                    fill_diff_raster =
                        lsdtt::algebra::raster_subtract(filled_dem, working_dem);
                    have_fill_diff = true;

                    if (bf["write_fill_diff_raster"])
                    {
                        std::string o = pp.make_output_path( "_fill_diff.tif");
                        std::cout << "Writing: " << o << "\n";
                        loader.write_full(fill_diff_raster.get_data(),
                                          fill_diff_raster.get_meta(), o, out_fmt);
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: fill diff raster failed: "
                              << e.what() << "\n";
                }
            }
            did_fill = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: pit filling failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 7 - Depression polygons
    // -------------------------------------------------------------------------
    if (bf["write_sink_candidates"])
    {
        std::cout << "\n--- Sink candidates ---\n";

        try
        {
            lsdtt::sink_detection_options sink_opts;
            // fill_difference, not depression_hierarchy: on a rough surface the
            // hierarchy's leaves are individually tiny and a pixel threshold can
            // remove every one of them. See docs/flow_routing.md.
            sink_opts.method = lsdtt::sink_detection_method::fill_difference;
            sink_opts.min_pixels_in_patch = ifs["sink_candidate_min_pixels"];

            const std::vector<lsdtt::sink_candidate> sinks =
                lsdtt::find_internal_drainage(working_dem, sink_opts);

            std::cout << "Found " << sinks.size()
                      << " candidate sink(s).\n";

            const lsdtt::vector_output_format vec_fmt =
                lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
            const std::string o =
                pp.make_output_path("_sink_candidates" +
                                    lsdtt::output_extension(vec_fmt));

            lsdtt::write_sink_candidates(
                loader, o, sinks, working_dem.get_meta(),
                lsdtt::resolve_sink_io_format(o, "same_as_extension"));

            std::cout << "Review these in a GIS, delete the ones that are not "
                         "real sinks, then pass the saved file to "
                         "lsdtt-internal-drainage as sink_candidate_fname.\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: sink candidate detection failed: "
                      << e.what() << "\n";
        }
    }

    if (bf["write_depression_polygons"])
    {
        std::cout << "\n--- Depression polygons ---\n";

        if (!have_fill_diff)
        {
            std::cerr << "Warning: fill diff not available - "
                         "skipping depression polygons.\n";
        }
        else
        {
            try
            {
                const lsdtt::vector_output_format vec_fmt =
                    lsdtt::resolve_vector_output_format(
                        lsdtt::resolve_output_format_override(
                            sf["legacy_depression_polygons_write_vector_format"],
                            sf["write_vector_format"]));

                const std::string o =
                    pp.make_output_path( "_depressions" +
                             lsdtt::output_extension(vec_fmt));

                std::cout << "Tracing depression outlines...\n";

                lsdtt::depression_polygons_result dep_result =
                    lsdtt::depressions_to_polygons(
                        fill_diff_raster,
                        ifs["depression_min_pixels"],
                        ff["depression_min_depth"]
                    );

                std::cout << "Found " << dep_result.records.size()
                          << " depression(s).\n";

                lsdtt::write_vector_records(
                    loader, o, "depressions",
                    dep_result.schema, dep_result.records,
                    filled_dem.get_meta().crs_wkt, vec_fmt);

                std::cout << "Writing: " << o << "\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: depression polygons failed: "
                          << e.what() << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 8 - Hillshade
    // -------------------------------------------------------------------------
    if (bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Hillshade ---\n";
        std::cout << "Azimuth: "  << ff["hillshade_azimuth"]  << " deg\n";
        std::cout << "Altitude: " << ff["hillshade_altitude"] << " deg\n";
        std::cout << "Z-factor: " << ff["hillshade_z_factor"] << "\n";

        const lsdtt::lsdtt_raster_cpu& hs_src = did_fill
            ? filled_dem
            : working_dem;

        try
        {
            lsdtt::lsdtt_raster_cpu hillshade = lsdtt::compute_hillshade(
                hs_src,
                ff["hillshade_azimuth"],
                ff["hillshade_altitude"],
                ff["hillshade_z_factor"]
            );

            std::string o = pp.make_output_path( "_hillshade.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(hillshade.get_data(), hillshade.get_meta(),
                              o, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: hillshade failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 8b - Multidirectional hillshade
    // -------------------------------------------------------------------------
    if (bf["write_multidirectional_hillshade_raster"])
    {
        std::cout << "\n--- Multidirectional hillshade ---\n";
        std::cout << "N directions : "
                  << ifs["multidirectional_hillshade_n_directions"] << "\n";
        std::cout << "Start azimuth: "
                  << ff["multidirectional_hillshade_start_azimuth"] << " deg\n";
        std::cout << "Altitude     : " << ff["hillshade_altitude"] << " deg\n";
        std::cout << "Z-factor     : " << ff["hillshade_z_factor"] << "\n";

        const lsdtt::lsdtt_raster_cpu& hs_src = did_fill
            ? filled_dem
            : working_dem;

        try
        {
            lsdtt::lsdtt_raster_cpu hillshade_md =
                lsdtt::compute_multidirectional_hillshade(
                    hs_src,
                    ifs["multidirectional_hillshade_n_directions"],
                    ff["hillshade_altitude"],
                    ff["hillshade_z_factor"],
                    ff["multidirectional_hillshade_start_azimuth"]
                );

            std::string o = pp.make_output_path("_hillshade_multidirectional.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(
                hillshade_md.get_data(),
                hillshade_md.get_meta(),
                o,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: multidirectional hillshade failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 9 - Vector clipping
    // -------------------------------------------------------------------------
    if (bf["clip_vector_data_to_raster"])
    {
        std::cout << "\n--- Vector clipping ---\n";

        const std::string& vec_in = sf["vector_data_fname"];
        if (vec_in.empty())
        {
            std::cerr << "Warning: clip_vector_data_to_raster is true but "
                         "vector_data_fname is empty - skipping.\n";
        }
        else
        {
            try
            {
                const lsdtt::vector_output_format vec_fmt =
                    lsdtt::resolve_vector_output_format(
                    lsdtt::resolve_output_format_override(
                        sf["legacy_vector_clip_write_vector_format"],
                        sf["write_vector_format"])
                    );

                const std::string vec_out =
                    pp.make_output_path( "_clipped_vector" +
                             lsdtt::output_extension(vec_fmt));

                std::cout << "Input vector : " << vec_in  << "\n";
                std::cout << "Output vector: " << vec_out << "\n";

                lsdtt::lsdtt_data_loader::clip_vector_to_raster_bbox(
                    vec_in,
                    vec_out,
                    dem.get_meta(),
                    vec_fmt
                );

                std::cout << "Vector clipping done.\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: vector clipping failed: "
                          << e.what() << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 10 - Sample rasters to points
    // -------------------------------------------------------------------------
    if (bf["sample_rasters_to_points"])
    {
        std::cout << "\n--- Sample rasters to points ---\n";

        const std::string point_path =
            pp.add_read_path_if_needed(sf["point_data_fname"]);
        const std::string list_csv =
            pp.add_read_path_if_needed(sf["raster_attribute_fname"]);
        const std::string layer_name = sf["point_data_layer_name"];

        if (sf["point_data_fname"].empty())
        {
            std::cerr << "Warning: sample_rasters_to_points is true but "
                        "point_data_fname is empty - skipping.\n";
        }
        else if (sf["raster_attribute_fname"].empty())
        {
            std::cerr << "Warning: sample_rasters_to_points is true but "
                        "raster_attribute_fname is empty - skipping.\n";
        }
        else
        {
            try
            {
                // Read the input point layer
                lsdtt::lsdtt_data_loader point_loader;
                point_loader.open_vector(point_path);

                std::string use_layer = layer_name;
                if (use_layer.empty())
                {
                    auto layer_names = point_loader.get_layer_names();
                    if (layer_names.empty())
                    {
                        throw std::runtime_error("Point dataset has no layers.");
                    }
                    use_layer = layer_names[0];
                }

                lsdtt::lsdtt_vector_layer points =
                    point_loader.read_vector_layer(use_layer);

                // You need this extra loader method:
                const std::string point_crs_wkt =
                    point_loader.get_vector_layer_crs_wkt(use_layer);

                // Read raster list
                std::vector<raster_attribute_request> requests =
                    read_raster_attribute_csv(list_csv, pp);

                if (requests.empty())
                {
                    throw std::runtime_error(
                        "No raster requests found in raster_attribute_fname."
                    );
                }

                // Sample each raster
                std::vector<std::vector<lsdtt::lsdtt_vector_field_value>> sampled_columns;
                sampled_columns.reserve(requests.size());

                for (const auto& req : requests)
                {
                    std::cout << "Sampling raster: " << req.raster_path
                            << " -> attribute: " << req.attribute_name << "\n";

                    sampled_columns.emplace_back(
                        sample_raster_to_vector_points(
                            req.raster_path, point_path, use_layer
                        )
                    );

                    if (sampled_columns.back().size() != points.records.size())
                    {
                        throw std::runtime_error(
                            "Point count mismatch while sampling raster: " +
                            req.raster_path
                        );
                    }
                }

                // Build output schema
                lsdtt::lsdtt_vector_schema out_schema = points.schema;
                for (const auto& req : requests)
                {
                    out_schema.field_names.push_back(req.attribute_name);
                    out_schema.field_types.push_back(
                        lsdtt::lsdtt_vector_field_type::double_type
                    );
                }

                const lsdtt::point_output_format pt_fmt =
                    lsdtt::resolve_point_output_format(
                        lsdtt::resolve_output_format_override(
                            sf["legacy_sampled_points_write_vector_format"],
                            sf["write_vector_format"])
                    );

                if (pt_fmt == lsdtt::point_output_format::csv)
                {
                    throw std::runtime_error(
                        "CSV output for vector point input is not yet implemented "
                        "in this first version."
                    );
                }

                const lsdtt::vector_output_format vec_fmt =
                    lsdtt::convert_point_to_vector_output_format(pt_fmt);

                const std::string out_file =
                    pp.make_output_path( "_sampled_points" +
                                lsdtt::output_extension(pt_fmt));

                point_loader.begin_vector_write(
                    out_file,
                    "sampled_points",
                    out_schema,
                    point_crs_wkt,
                    vec_fmt
                );

                for (std::size_t i = 0; i < points.records.size(); ++i)
                {
                    // GUARD, not a bare std::get. These records come from a
                    // user-supplied file, and a geometry this reader cannot
                    // interpret arrives as an empty variant -- an unguarded get
                    // would throw std::bad_variant_access with a message naming
                    // neither the file nor the feature.
                    if (!std::holds_alternative<lsdtt::lsdtt_point_geom>(
                            points.records[i].geometry))
                    {
                        throw std::runtime_error(
                            "point_data_fname: feature " + std::to_string(i) +
                            " carries no usable POINT geometry. Multipart "
                            "points are supported and are split automatically, "
                            "so this is a geometry type the reader does not "
                            "handle (a curve, or a line or polygon in a file "
                            "expected to hold points).");
                    }

                    const auto& geom =
                        std::get<lsdtt::lsdtt_point_geom>(points.records[i].geometry);

                    std::vector<lsdtt::lsdtt_vector_field_value> out_fields =
                        points.records[i].field_values;

                    for (std::size_t j = 0; j < sampled_columns.size(); ++j)
                    {
                        out_fields.emplace_back(sampled_columns[j][i]);
                    }

                    point_loader.write_point(geom.x, geom.y, out_fields);
                }

                point_loader.finish_vector_write();
                std::cout << "Wrote sampled points: " << out_file << "\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: point sampling failed: "
                        << e.what() << "\n";
            }
        }
    }


    // -------------------------------------------------------------------------
    // Parameter echo - always printed so runs are reproducible from the log
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-raster-preprocessing: done.\n";
    return 0;
}
