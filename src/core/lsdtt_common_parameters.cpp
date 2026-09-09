// lsdtt_common_parameters.cpp
//
// Shared parameter-group declarations. See the header for why this file
// exists, the prefix rule, the overwrite policy and the provenance of every
// name and default declared here.

#include "lsdtt_common_parameters.hpp"

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
// Flow routing and pit filling
// -----------------------------------------------------------------------------
void declare_flow_routing_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    const std::string group = prefixed_group_label(prefix, "Flow routing");

    // --- Pit filling -----------------------------------------------------
    // Concept A (plan §2.2): filling is MANDATORY, this only skips the work.
    const std::string k_raster_is_filled =
        prefixed_parameter_key(prefix, "raster_is_filled");
    bool_default_map[k_raster_is_filled] = false;
    help_map[k_raster_is_filled] = {"bool", "false",
        "Assert that the input raster has ALREADY been pit-filled.",
        "ASKS: has the input already been filled? Filling is REQUIRED for this "
        "analysis, so this switch does not choose whether to fill -- it only "
        "lets the tool skip filling a surface that is filled already. Setting "
        "it true for a DEM that still contains pits gives BAD routing, not "
        "different routing. Do not confuse this with fill_before_routing, "
        "which asks whether to fill at all.",
        group, "switch", ""};

    // Concept B (plan §2.2): filling is OPTIONAL and is an analytical choice.
    const std::string k_fill_before_routing =
        prefixed_parameter_key(prefix, "fill_before_routing");
    bool_default_map[k_fill_before_routing] = true;
    help_map[k_fill_before_routing] = {"bool", "true",
        "Choose whether to pit-fill the DEM before building the flow network.",
        "ASKS: should the DEM be filled at all? Filling is OPTIONAL here and "
        "this is a genuine analytical choice. True (the default) routes over a "
        "pit-free surface. False builds the flow network on the UNFILLED "
        "surface, so flow paths terminate at sinks instead of being carried "
        "through them. Do not confuse this with raster_is_filled, which is an "
        "assertion that the input is already filled.",
        group, "switch", ""};

    const std::string k_fill_min_slope =
        prefixed_parameter_key(prefix, "fill_min_slope");
    float_default_map[k_fill_min_slope] = 0.0001f;
    help_map[k_fill_min_slope] = {"float", "0.0001",
        "Minimum slope (m/m) imposed by pit filling.",
        "Passed to fill_pits() to gently enforce downstream descent across "
        "filled depressions and flats. Only used when filling actually "
        "happens.",
        group, "param", ""};

    // --- Routing ---------------------------------------------------------
    const std::string k_boundary_conditions =
        prefixed_parameter_key(prefix, "boundary_conditions");
    string_default_map[k_boundary_conditions] = "oooo";
    help_map[k_boundary_conditions] = {"string", "oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        group, "param", ""};

    const std::string k_threshold_contributing_pixels =
        prefixed_parameter_key(prefix, "threshold_contributing_pixels");
    int_default_map[k_threshold_contributing_pixels] = 1000;
    help_map[k_threshold_contributing_pixels] = {"int", "1000",
        "Contributing-pixel threshold defining a channel head.",
        "The total number of upslope contributing pixels from the flow-routing "
        "object, not the number of direct donors to a pixel. Reduce this value "
        "to extract a denser network. This is a PIXEL COUNT and therefore "
        "depends on the grid resolution: the equivalent drainage area is "
        "threshold_contributing_pixels x cell area. Where a genuine area is "
        "wanted instead, the convention is <what>_area_threshold_m2.",
        group, "param", ""};

    // --- Run control -----------------------------------------------------
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
// Output formats
// -----------------------------------------------------------------------------
void declare_output_format_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    (void)int_default_map;
    (void)float_default_map;

    const std::string group = prefixed_group_label(prefix, "Output");

    const std::string k_write_output_format =
        prefixed_parameter_key(prefix, "write_output_format");
    string_default_map[k_write_output_format] = "cog";
    help_map[k_write_output_format] = {"string", "cog",
        "Raster output format: cog, gtiff, bil, png or same_as_input.",
        "cog is a Cloud Optimised GeoTIFF and is the project default. It falls "
        "back to gtiff for small rasters, where the tiling overhead is not "
        "worth it. There is ONE raster format parameter per CLI: it governs "
        "every raster the tool writes. Vector outputs are governed separately "
        "by write_vector_format.",
        group, "io", ""};

    const std::string k_write_vector_format =
        prefixed_parameter_key(prefix, "write_vector_format");
    string_default_map[k_write_vector_format] = "flatgeobuf";
    help_map[k_write_vector_format] = {"string", "flatgeobuf",
        "Vector output format: flatgeobuf, gpkg or csv.",
        "Governs every vector layer this tool writes. FlatGeobuf is preferred: "
        "it is compact, carries its own georeferencing, and does not have "
        "GeoPackage's problems writing across drives (a common WSL2 "
        "situation). Per-output overrides exist only where a container "
        "genuinely differs. Raster outputs are governed separately by "
        "write_output_format.",
        group, "io", ""};

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
        "Allowed values: none, expanding_square, downstream. "
        "none: use the mapped raster pixel directly. "
        "expanding_square: search outward ring-by-ring until a pixel exceeding "
        "basin_outlet_contributing_pixels_threshold is found. "
        "downstream: follow the D8 receiver chain until a qualifying pixel is "
        "found.",
        group, "param", ""};

    const std::string k_basin_outlet_contributing_pixels_threshold =
        prefixed_parameter_key(prefix,
                               "basin_outlet_contributing_pixels_threshold");
    int_default_map[k_basin_outlet_contributing_pixels_threshold] = 1000;
    help_map[k_basin_outlet_contributing_pixels_threshold] = {"int", "1000",
        "Contributing-pixel threshold a pixel must exceed to accept an outlet.",
        "This is the total number of upslope contributing pixels from the "
        "flow-routing object, not the number of direct donors to a pixel. Used "
        "only when basin_outlet_snapping_method is not none. It is a separate "
        "quantity from threshold_contributing_pixels, which sets where the "
        "channel network begins.",
        group, "param", ""};

    const std::string k_basin_outlet_min_search_radius_pixels =
        prefixed_parameter_key(prefix, "basin_outlet_min_search_radius_pixels");
    int_default_map[k_basin_outlet_min_search_radius_pixels] = 0;
    help_map[k_basin_outlet_min_search_radius_pixels] = {"int", "0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "0 means the mapped input pixel is tested first.",
        group, "param", ""};

    const std::string k_basin_outlet_max_search_radius_pixels =
        prefixed_parameter_key(prefix, "basin_outlet_max_search_radius_pixels");
    int_default_map[k_basin_outlet_max_search_radius_pixels] = 25;
    help_map[k_basin_outlet_max_search_radius_pixels] = {"int", "25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "The search stops at the first radius containing one or more "
        "qualifying pixels.",
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
void declare_synthetic_domain_parameters(
    std::map<std::string, bool>&        bool_default_map,
    std::map<std::string, int>&         int_default_map,
    std::map<std::string, float>&       float_default_map,
    std::map<std::string, std::string>& string_default_map,
    lsdtt_help_registry&                help_map,
    const std::string&                  prefix)
{
    (void)bool_default_map;

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

    const std::string k_cell_size_metres =
        prefixed_parameter_key(prefix, "cell_size_metres");
    float_default_map[k_cell_size_metres] = 1.0f;
    help_map[k_cell_size_metres] = {"float", "1.0",
        "Cell size of the synthetic model domain in metres.",
        "Sets the map scale of the generated surface. Anything expressed as a "
        "pixel count (for example threshold_contributing_pixels) corresponds to "
        "a different ground area as this changes.",
        group, "param", ""};

    const std::string k_random_seed =
        prefixed_parameter_key(prefix, "random_seed");
    int_default_map[k_random_seed] = 33;
    help_map[k_random_seed] = {"int", "33",
        "Seed for every stochastic generator used to build the domain.",
        "A fixed seed makes a run reproducible. Change it to draw a different "
        "realisation of the same statistical surface.",
        group, "param", ""};

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

// -----------------------------------------------------------------------------
// Synthetic noise surfaces
// -----------------------------------------------------------------------------
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

    const std::string k_noise_target_stddev_metres =
        prefixed_parameter_key(prefix, "noise_target_stddev_metres");
    float_default_map[k_noise_target_stddev_metres] = 1.0f;
    help_map[k_noise_target_stddev_metres] = {"float", "1.0",
        "Target standard deviation of the generated noise surface (metres).",
        "The raw surface is rescaled to this standard deviation. Used when "
        "noise_type is red_noise or diamond_square.",
        group, "param", ""};

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

} // namespace lsdtt
