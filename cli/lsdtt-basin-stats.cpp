// lsdtt-basin-stats.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-basin-stats
//
// Reads a DEM and a CSV of basin outlet points, extracts basins, and computes
// a suite of hydromorphic statistics for each basin.
//
// Basin outlet CSV requirements:
//   - parameter: basin_outlet_csv
//   - must contain columns named: latitude, longitude
//   - optional columns: gauge_id, gauge_name
//   - coordinates are interpreted as WGS84 / EPSG:4326 decimal degrees
//
// Workflow:
//   1.  Load raster
//   2.  Check raster is projected in metres
//   3.  Remove seas by masking pixels at or below minimum_elevation to nodata
//   4.  Optionally write a hillshade of the sea-masked DEM
//   5.  Fill pits before routing (default: true)
//   6.  Build D8 flow routing
//   7.  Read outlet CSV (including optional gauge_id / gauge_name columns)
//   8.  Transform each outlet from WGS84 to raster CRS
//   9.  Convert to raster pixel; optionally snap to channel
//  10.  Compute gradient magnitude raster (for slope percentiles)
//  11.  Compute flow-distance-from-outlet raster (for flow distance percentiles)
//  12.  Optionally build junction/channel network (for channel length by order)
//  13.  Per basin:
//        a. Check whether any upslope pixel borders a nodata cell
//        b. Extract DEM elevations, slope magnitudes, flow distances
//        c. Compute descriptive stats and percentiles for each quantity
//        d. Accumulate channel length by Strahler order (1-6)
//  14.  Write basin statistics to CSV
//
// IMPORTANT:
//   - statistics functions assume the supplied flat data are already valid
//   - this CLI is responsible for excluding nodata before calling statistics
//
// Usage:
//   lsdtt-basin-stats <path> <parameter_file>
//   lsdtt-basin-stats <parameter_file>
//   lsdtt-basin-stats -h
//   lsdtt-basin-stats -v
//   lsdtt-basin-stats -c
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>


#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_basin_cpu.hpp"
#include "../src/cpu/lsdtt_basin_outline_cpu.hpp"
#include "../src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.hpp"
#include "../src/cpu/lsdtt_channel_network_analysis_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/utils/lsdtt_polygon_geometry.hpp"
#include "../src/utils/lsdtt_statistics.hpp"
#include "../src/utils/lsdtt_histogram.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/core/lsdtt_csv_table.hpp"


static const std::string k_version  = "0.3.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-basin-stats";

// -----------------------------------------------------------------------------
// Local record types
// -----------------------------------------------------------------------------

struct basin_stats_record
{
    int          basin_id       = -1;
    std::string  gauge_id       = "";
    std::string  gauge_name     = "";
    double       latitude       = 0.0;
    double       longitude      = 0.0;
    double       outlet_x       = 0.0;
    double       outlet_y       = 0.0;
    int          outlet_row     = -1;
    int          outlet_col     = -1;
    std::int64_t outlet_flat    = -1;
    std::size_t  n_pixels       = 0;
    double       drainage_area_km2 = 0.0;
    bool         nodata_on_edge = false;

    lsdtt::lsdtt_descriptive_stats stats;  // elevation statistics

    // Centroid in projected CRS and WGS84
    double centroid_x         = 0.0;
    double centroid_y         = 0.0;
    double centroid_latitude  = 0.0;
    double centroid_longitude = 0.0;

    // Outlet -> centre-of-mass bearing (degrees from North, clockwise) + distance
    double com_bearing_deg    = std::numeric_limits<double>::quiet_NaN();
    double com_distance_m     = std::numeric_limits<double>::quiet_NaN();

    // Slope (gradient magnitude) percentiles
    double slope_25th = 0.0;
    double slope_50th = 0.0;
    double slope_75th = 0.0;

    // Flow distance from outlet percentiles (metres)
    double flow_distance_25th = 0.0;
    double flow_distance_50th = 0.0;
    double flow_distance_75th = 0.0;

    // Channel length (m) by Strahler stream order (keys 1-6)
    // -1.0 indicates channel stats were not calculated
    std::map<int, double> channel_length_by_order;

    // Basin planform shape metrics (see docs/tex_technical/basin_shape_metrics.tex)
    // basin_length_m is L_b under the chosen basin_length_definition.
    bool   shape_valid           = false;
    double basin_length_m        = 0.0;
    double elongation_ratio      = 0.0;
    double form_factor           = 0.0;
    double circularity_ratio     = 0.0;
    double gravelius_compactness = 0.0;
    double solidity              = 0.0;
    double eccentricity          = 0.0;
    double orientation_deg       = 0.0;
    double hypsometric_integral  = 0.0;
};


// -----------------------------------------------------------------------------
// basin_pdf_cdf_path
// -----------------------------------------------------------------------------
static std::string basin_pdf_cdf_path(
    const lsdtt::lsdtt_param_parser& pp,
    int basin_id)
{
    return pp.get_full_write_path() +
           "_basin_" + std::to_string(basin_id) +
           "_pdf_cdf.csv";
}

// -----------------------------------------------------------------------------
// write_basin_pdf_cdf_csv
// -----------------------------------------------------------------------------
static void write_basin_pdf_cdf_csv(
    const std::string& csv_path,
    int basin_id,
    const lsdtt::lsdtt_histogram& hist)
{
    std::ofstream out(csv_path);
    if (!out)
    {
        throw std::runtime_error(
            "Could not open basin PDF/CDF CSV for writing: " + csv_path
        );
    }

    out
        << "basin_id"
        << ",bin_index"
        << ",bin_lower_edge"
        << ",bin_upper_edge"
        << ",bin_centre"
        << ",bin_width"
        << ",count"
        << ",probability"
        << ",pdf"
        << ",cdf"
        << "\n";

    for (std::size_t i = 0; i < hist.n_bins; ++i)
    {
        out
            << basin_id << ","
            << i << ","
            << hist.bin_edges[i] << ","
            << hist.bin_edges[i + 1] << ","
            << hist.bin_centres[i] << ","
            << hist.bin_widths[i] << ","
            << hist.counts[i] << ","
            << hist.probabilities[i] << ","
            << hist.pdf[i] << ","
            << hist.cdf[i]
            << "\n";
    }
}


// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt3 basin statistics                  ||\n"
        << "|| DEM loading, routing, basin hydromorphic stats      ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "|| The original intent of this CLI                     ||\n"
        << "|| was to build hydromorphic data for the CAMELGB      ||\n"
        << "|| dataset to see if morphometric variables could be   ||\n"
        << "|| used for lsdtm hydrologic prediction but this has   ||\n"
        << "|| been extended to get more routine basin stats for   ||\n"
        << "|| standard morphometric analysis                      ||\n"
        << "|| NOTE: For stats on surface metrics (gradient, curv) ||\n"
        << "|| Use lsdtt-basin-surface-metrics                     ||\n"
        << "=========================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// join_boundary_conditions
// -----------------------------------------------------------------------------
static std::string join_boundary_conditions(
    const std::vector<std::string>& bc)
{
    std::string out;
    for (const auto& s : bc)
    {
        out += s;
    }
    return out;
}
// -----------------------------------------------------------------------------
// write_basin_stats_csv
// -----------------------------------------------------------------------------
static void write_basin_stats_csv(
    const std::string& csv_path,
    const std::vector<basin_stats_record>& rows,
    bool has_channel_stats,
    bool has_shape_metrics,
    bool has_com_bearing)
{
    std::ofstream out(csv_path);
    if (!out)
    {
        throw std::runtime_error(
            "Could not open output CSV for writing: " + csv_path
        );
    }

    out
        << "basin_id"
        << ",gauge_id"
        << ",gauge_name"
        << ",latitude"
        << ",longitude"
        << ",outlet_x"
        << ",outlet_y"
        << ",outlet_row"
        << ",outlet_col"
        << ",outlet_flat"
        << ",n_pixels"
        << ",drainage_area_km2"
        << ",nodata_on_edge"
        << ",centroid_x"
        << ",centroid_y"
        << ",centroid_latitude"
        << ",centroid_longitude";

    if (has_com_bearing)
    {
        out << ",com_bearing_deg" << ",com_distance_m";
    }

    out
        << ",elevation_min"
        << ",elevation_max"
        << ",elevation_range"
        << ",elevation_mean"
        << ",elevation_25th"
        << ",elevation_50th"
        << ",elevation_75th"
        << ",elevation_iqr"
        << ",elevation_stddev"
        << ",elevation_mad"
        << ",slope_25th"
        << ",slope_50th"
        << ",slope_75th"
        << ",flow_distance_25th"
        << ",flow_distance_50th"
        << ",flow_distance_75th";

    if (has_shape_metrics)
    {
        out
            << ",basin_length_m"
            << ",elongation_ratio"
            << ",form_factor"
            << ",circularity_ratio"
            << ",gravelius_compactness"
            << ",solidity"
            << ",eccentricity"
            << ",orientation_deg"
            << ",hypsometric_integral";
    }

    if (has_channel_stats)
    {
        for (int o = 1; o <= 6; ++o)
        {
            out << ",channel_length_O" << o;
        }
    }

    out << "\n";

    for (const auto& row : rows)
    {
        out
            << row.basin_id << ","
            << row.gauge_id << ","
            << row.gauge_name << ","
            << row.latitude << ","
            << row.longitude << ","
            << row.outlet_x << ","
            << row.outlet_y << ","
            << row.outlet_row << ","
            << row.outlet_col << ","
            << row.outlet_flat << ","
            << row.n_pixels << ","
            << row.drainage_area_km2 << ","
            << (row.nodata_on_edge ? 1 : 0) << ","
            << row.centroid_x << ","
            << row.centroid_y << ","
            << row.centroid_latitude << ","
            << row.centroid_longitude;

        if (has_com_bearing)
        {
            out << "," << row.com_bearing_deg << "," << row.com_distance_m;
        }

        out
            << "," << row.stats.min << ","
            << row.stats.max << ","
            << row.stats.range << ","
            << row.stats.mean << ","
            << row.stats.q25 << ","
            << row.stats.median << ","
            << row.stats.q75 << ","
            << row.stats.iqr << ","
            << row.stats.stddev_sample << ","
            << row.stats.median_absolute_deviation << ","
            << row.slope_25th << ","
            << row.slope_50th << ","
            << row.slope_75th << ","
            << row.flow_distance_25th << ","
            << row.flow_distance_50th << ","
            << row.flow_distance_75th;

        if (has_shape_metrics)
        {
            out
                << "," << row.basin_length_m
                << "," << row.elongation_ratio
                << "," << row.form_factor
                << "," << row.circularity_ratio
                << "," << row.gravelius_compactness
                << "," << row.solidity
                << "," << row.eccentricity
                << "," << row.orientation_deg
                << "," << row.hypsometric_integral;
        }

        if (has_channel_stats)
        {
            for (int o = 1; o <= 6; ++o)
            {
                auto it = row.channel_length_by_order.find(o);
                out << "," << (it != row.channel_length_by_order.end()
                                   ? it->second : 0.0);
            }
        }

        out << "\n";
    }
}

// -----------------------------------------------------------------------------
// write_basin_polygons
//
// Writes one basin-outline polygon per basin to a vector layer (FlatGeobuf /
// GeoPackage), with the per-basin hydromorphic and (optionally) shape and
// centre-of-mass metrics attached as typed attributes - i.e. the basin metrics
// joined to their associated basin polygons.
// -----------------------------------------------------------------------------
static void write_basin_polygons(
    const std::string&                      out_path,
    const std::vector<basin_stats_record>&  rows,
    const lsdtt::lsdtt_flow_base&           flow,
    const lsdtt::lsdtt_raster_meta&         meta,
    lsdtt::vector_output_format             fmt,
    bool                                    has_channel_stats,
    bool                                    has_shape_metrics,
    bool                                    has_com_bearing,
    bool                                    spatial_index)
{
    using FT = lsdtt::lsdtt_vector_field_type;

    lsdtt::lsdtt_vector_schema sch;
    auto add = [&](const std::string& name, FT t)
    {
        sch.field_names.push_back(name);
        sch.field_types.push_back(t);
    };

    add("basin_id", FT::int64_type);
    add("gauge_id", FT::string_type);
    add("gauge_name", FT::string_type);
    add("latitude", FT::double_type);
    add("longitude", FT::double_type);
    add("drainage_area_km2", FT::double_type);
    add("n_pixels", FT::int64_type);
    add("nodata_on_edge", FT::int64_type);
    add("centroid_latitude", FT::double_type);
    add("centroid_longitude", FT::double_type);
    add("elevation_min", FT::double_type);
    add("elevation_max", FT::double_type);
    add("elevation_mean", FT::double_type);
    add("elevation_range", FT::double_type);
    add("slope_50th", FT::double_type);
    add("flow_distance_50th", FT::double_type);
    if (has_com_bearing)
    {
        add("com_bearing_deg", FT::double_type);
        add("com_distance_m", FT::double_type);
    }
    if (has_shape_metrics)
    {
        add("basin_length_m", FT::double_type);
        add("elongation_ratio", FT::double_type);
        add("form_factor", FT::double_type);
        add("circularity_ratio", FT::double_type);
        add("gravelius_compactness", FT::double_type);
        add("solidity", FT::double_type);
        add("eccentricity", FT::double_type);
        add("orientation_deg", FT::double_type);
        add("hypsometric_integral", FT::double_type);
    }
    if (has_channel_stats)
    {
        for (int o = 1; o <= 6; ++o)
            add("channel_length_O" + std::to_string(o), FT::double_type);
    }

    lsdtt::lsdtt_data_loader pw;
    pw.begin_vector_write(out_path, "basin_polygons", sch,
                          meta.crs_wkt, fmt, spatial_index);

    for (const auto& row : rows)
    {
        try
        {
            lsdtt::lsdtt_basin_cpu basin(row.outlet_flat, flow);
            lsdtt::lsdtt_basin_outline_cpu outline = basin.get_basin_outline();
            const lsdtt::lsdtt_vector_record rec =
                outline.get_outline_as_polygon_record(row.basin_id);

            if (std::holds_alternative<std::monostate>(rec.geometry)) continue;
            const auto& poly = std::get<lsdtt::lsdtt_polygon_geom>(rec.geometry);
            if (poly.exterior.size() < 4) continue;

            std::vector<lsdtt::lsdtt_vector_field_value> v;
            v.emplace_back(static_cast<std::int64_t>(row.basin_id));
            v.emplace_back(row.gauge_id);
            v.emplace_back(row.gauge_name);
            v.emplace_back(row.latitude);
            v.emplace_back(row.longitude);
            v.emplace_back(row.drainage_area_km2);
            v.emplace_back(static_cast<std::int64_t>(row.n_pixels));
            v.emplace_back(static_cast<std::int64_t>(row.nodata_on_edge ? 1 : 0));
            v.emplace_back(row.centroid_latitude);
            v.emplace_back(row.centroid_longitude);
            v.emplace_back(row.stats.min);
            v.emplace_back(row.stats.max);
            v.emplace_back(row.stats.mean);
            v.emplace_back(row.stats.range);
            v.emplace_back(row.slope_50th);
            v.emplace_back(row.flow_distance_50th);
            if (has_com_bearing)
            {
                v.emplace_back(row.com_bearing_deg);
                v.emplace_back(row.com_distance_m);
            }
            if (has_shape_metrics)
            {
                v.emplace_back(row.basin_length_m);
                v.emplace_back(row.elongation_ratio);
                v.emplace_back(row.form_factor);
                v.emplace_back(row.circularity_ratio);
                v.emplace_back(row.gravelius_compactness);
                v.emplace_back(row.solidity);
                v.emplace_back(row.eccentricity);
                v.emplace_back(row.orientation_deg);
                v.emplace_back(row.hypsometric_integral);
            }
            if (has_channel_stats)
            {
                for (int o = 1; o <= 6; ++o)
                {
                    auto it = row.channel_length_by_order.find(o);
                    v.emplace_back(it != row.channel_length_by_order.end()
                                       ? it->second : 0.0);
                }
            }

            pw.write_polygon(poly.exterior, v);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: basin polygon for basin_id=" << row.basin_id
                      << " failed: " << e.what() << "\n";
        }
    }

    pw.finish_vector_write();
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
    // Load parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = single-threaded.",
        "Threading", "io", ""};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format.",
        "Output format", "io", ""};

    // --- Basin outlets ---
    string_default_map["basin_outlet_csv"] = "NULL";
    help_map["basin_outlet_csv"] = {"string","NULL",
        "CSV file of requested basin outlet points.",
        "Required. CSV must contain columns named latitude and longitude. "
        "Optional columns: gauge_id, gauge_name. "
        "Coordinates are interpreted as WGS84 (EPSG:4326) decimal degrees.",
        "Basin outlets", "param", ""};

    // --- Basin outlet snapping ---
    string_default_map["basin_outlet_snapping_method"] = "none";
    help_map["basin_outlet_snapping_method"] = {"string","none",
        "Optional snapping of basin outlet CSV points to pixels that exceed "
        "a contributing-pixel threshold.",
        "Allowed values: none, expanding_square, downstream. "
        "none: use the mapped raster pixel directly. "
        "expanding_square: search outward ring-by-ring until a qualifying pixel "
        "is found. downstream: follow the D8 receiver chain until a qualifying "
        "pixel is found.",
        "Basin outlet snapping", "param", ""
    };

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for optional basin outlet snapping.",
        "This is the total number of upslope contributing pixels from the "
        "flow-routing object, not the number of direct donors to a pixel.",
        "Basin outlet snapping", "param", ""
    };

    int_default_map["basin_outlet_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "0 means the mapped input pixel is tested first.",
        "Basin outlet snapping", "param", ""
    };

    int_default_map["basin_outlet_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "The search stops at the first radius containing one or more qualifying pixels.",
        "Basin outlet snapping", "param", ""
    };

    // --- Sea removal ---
    float_default_map["minimum_elevation"] = 0.0f;
    help_map["minimum_elevation"] = {"float","0.0",
        "Pixels at or below this elevation are set to nodata.",
        "Default of 0.0 removes sea and ocean pixels from many DEM products.",
        "Sea removal", "param", ""};

    // --- Routing preparation ---
    bool_default_map["fill_before_routing"] = true;
    help_map["fill_before_routing"] = {"bool","true",
        "Fill pits before constructing flow routing.",
        "Recommended for basin extraction from DEMs.",
        "Routing preparation", "param", ""};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed during pit filling.",
        "Used when fill_before_routing is true.",
        "Routing preparation", "param", ""};

    // --- Basin quality control ---
    bool_default_map["reject_nodata_edge_basins"] = true;
    help_map["reject_nodata_edge_basins"] = {"bool","true",
        "Reject basins whose upslope area touches any nodata or raster-edge pixel.",
        "When true, basins that drain from the edge of the DEM are skipped. "
        "This prevents incorrect drainage area and statistics for incomplete basins.",
        "Basin quality control", "param", ""};

    // --- Hillshade ---
    bool_default_map["write_hillshade_raster"] = false;
    help_map["write_hillshade_raster"] = {"bool","false",
        "Write a hillshade raster after sea masking.",
        "Useful for quick visual QC of the masked DEM.",
        "Hillshade", "switch", ""};

    float_default_map["hillshade_azimuth"] = 315.0f;
    help_map["hillshade_azimuth"] = {"float","315",
        "Sun azimuth in degrees clockwise from north.",
        "315 is the standard cartographic northwest illumination.",
        "Hillshade", "param", "write_hillshade_raster"};

    float_default_map["hillshade_altitude"] = 45.0f;
    help_map["hillshade_altitude"] = {"float","45",
        "Sun altitude above horizon in degrees.",
        "45 degrees is a common default.",
        "Hillshade", "param", "write_hillshade_raster"};

    float_default_map["hillshade_z_factor"] = 1.0f;
    help_map["hillshade_z_factor"] = {"float","1.0",
        "Vertical exaggeration factor for hillshade.",
        "Values greater than 1 emphasise relief.",
        "Hillshade", "param", "write_hillshade_raster"};

    // --- Channel statistics ---
    bool_default_map["calculate_channel_stats"] = true;
    help_map["calculate_channel_stats"] = {"bool","true",
        "Build a channel/junction network and compute channel length by stream order.",
        "Requires threshold_contributing_pixels to define channel sources.",
        "Channel statistics", "switch", ""};

    int_default_map["threshold_contributing_pixels"] = 500;
    help_map["threshold_contributing_pixels"] = {"int","500",
        "Minimum number of contributing pixels to define a channel source.",
        "Used when calculate_channel_stats is true. Larger values give a sparser network.",
        "Channel statistics", "param", ""};

    bool_default_map["calculate_junction_area_ratios"] = false;
    help_map["calculate_junction_area_ratios"] = {"bool","false",
        "Compute per-junction drainage-area ratios (complement to Tokunaga).",
        "Opt-in. Builds a channel network (using threshold_contributing_pixels) "
        "if one is not already built, then writes the ratio r = A_trib/A_main at "
        "every interior confluence plus a distribution summary. Outputs: "
        "<write_fname>_junction_area_ratios.csv, "
        "<write_fname>_junction_area_ratio_summary.csv (with branch-bootstrap "
        "confidence intervals), and a _junction_area_ratio_decades.csv collapse "
        "table.",
        "Channel statistics", "switch", ""};

    int_default_map["junction_area_ratio_bootstrap"] = 1000;
    help_map["junction_area_ratio_bootstrap"] = {"int","1000",
        "Branch-bootstrap replicates for junction area-ratio confidence intervals.",
        "Resampling is over non-mainstem tributary branches. 0 disables it.",
        "Channel statistics", "param", ""};

    bool_default_map["write_basin_polygons"] = false;
    help_map["write_basin_polygons"] = {"bool","false",
        "Write basin outline polygons with the per-basin metrics as attributes.",
        "Attaches the hydromorphic stats (and shape metrics / centre-of-mass "
        "bearing when those are enabled) to each basin's outline polygon as a "
        "vector layer. Output: <write_fname>_basin_polygons.<ext>. "
        "Format set by basin_polygon_format.",
        "Channel statistics", "switch", ""};

    string_default_map["basin_polygon_format"] = "flatgeobuf";
    help_map["basin_polygon_format"] = {"string","flatgeobuf",
        "Vector format for the basin-polygon output.",
        "flatgeobuf (default) or gpkg.",
        "Channel statistics", "param", ""};

    bool_default_map["calculate_centre_of_mass_bearing"] = false;
    help_map["calculate_centre_of_mass_bearing"] = {"bool","false",
        "Compute the bearing from each basin outlet to its centre of mass.",
        "Opt-in. The basin centre of mass is the planimetric centroid of all "
        "upslope pixels (already reported as centroid_x/centroid_y). This adds "
        "com_bearing_deg (geographic bearing from North, clockwise, outlet -> "
        "centre of mass) and com_distance_m columns. The outlet is the user-"
        "supplied (optionally snapped) outlet pixel; the channel is not followed.",
        "Channel statistics", "switch", ""};

    // --- Basin shape metrics ---
    bool_default_map["calculate_shape_metrics"] = false;
    help_map["calculate_shape_metrics"] = {"bool","false",
        "Trace each basin outline polygon and compute planform shape metrics.",
        "Opt-in; off by default so the original curated hydromorphic stat set is "
        "reproduced unless explicitly requested. When true, adds "
        "basin_length_m, elongation_ratio, form_factor, circularity_ratio, "
        "gravelius_compactness, solidity, eccentricity, orientation_deg, and "
        "hypsometric_integral columns to the output CSV. When false, the outline "
        "is not traced and those columns are omitted, reproducing the original "
        "basin-stats output. See docs/tex_technical/basin_shape_metrics.tex.",
        "Basin shape metrics", "switch", ""};

    string_default_map["basin_length_definition"] = "flow_path";
    help_map["basin_length_definition"] = {"string","flow_path",
        "Definition of basin length L_b used in elongation ratio and form factor.",
        "flow_path: the longest flow-path distance from the outlet (hydrologic "
        "basin length; the default and the usual geomorphic convention). "
        "geometric: the maximum Feret diameter of the basin outline polygon "
        "(longest straight-line caliper). The two are not interchangeable and "
        "Re, Rf both scale with L_b. See docs/tex_technical/basin_shape_metrics.tex.",
        "Basin shape metrics", "param", "calculate_shape_metrics"};

    float_default_map["shape_perimeter_simplify_cells"] = 1.0f;
    help_map["shape_perimeter_simplify_cells"] = {"float","1.0",
        "Douglas-Peucker simplification tolerance for the perimeter, in cells.",
        "The basin outline is traced on the pixel-corner grid and is staircased, "
        "which inflates the perimeter and biases circularity_ratio and "
        "gravelius_compactness (the coastline paradox). The perimeter for those "
        "two metrics is measured on a ring simplified with tolerance "
        "(this value * cell_size) in metres. Set to 0 to disable simplification "
        "and use the raw staircased perimeter. Area-based metrics (elongation, "
        "form factor, solidity) are unaffected.",
        "Basin shape metrics", "param", "calculate_shape_metrics"};

    // --- Basin distribution output ---
    bool_default_map["write_basin_pdf_cdf_csv"] = false;
    help_map["write_basin_pdf_cdf_csv"] = {"bool","false",
        "Write PDF/CDF CSV output for each basin.",
        "If true, writes one CSV per basin containing histogram counts, "
        "probabilities, PDF, and CDF of basin elevations.",
        "Basin distribution output", "switch", ""};

    int_default_map["n_histogram_bins"] = 100;
    help_map["n_histogram_bins"] = {"int","100",
        "Number of bins used for basin elevation histogram/PDF/CDF output.",
        "Used when write_basin_pdf_cdf_csv is true. Default is 100 bins.",
        "Basin distribution output", "param", "write_basin_pdf_cdf_csv"};


    // -------------------------------------------------------------------------
    // Parse parameters
    // -------------------------------------------------------------------------
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

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    auto to_lower_copy = [](std::string s)
    {
        for (char& ch : s)
        {
            ch = static_cast<char>(
                std::tolower(static_cast<unsigned char>(ch))
            );
        }
        return s;
    };

    const std::string basin_outlet_snapping_method =
        to_lower_copy(sf["basin_outlet_snapping_method"]);

    if (basin_outlet_snapping_method != "none" &&
        basin_outlet_snapping_method != "expanding_square" &&
        basin_outlet_snapping_method != "downstream")
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                  << "none, expanding_square, downstream.\n";
        return 1;
    }

    if (ifs["basin_outlet_contributing_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: basin_outlet_contributing_pixels_threshold "
                  << "must be > 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_min_search_radius_pixels"] < 0)
    {
        std::cerr << "FATAL: basin_outlet_min_search_radius_pixels "
                  << "must be >= 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_max_search_radius_pixels"] <
        ifs["basin_outlet_min_search_radius_pixels"])
    {
        std::cerr << "FATAL: basin_outlet_max_search_radius_pixels must be >= "
                  << "basin_outlet_min_search_radius_pixels.\n";
        return 1;
    }

    if (bf["calculate_channel_stats"] &&
        ifs["threshold_contributing_pixels"] <= 0)
    {
        std::cerr << "FATAL: threshold_contributing_pixels must be > 0 "
                  << "when calculate_channel_stats is true.\n";
        return 1;
    }

    const std::string basin_length_definition =
        to_lower_copy(sf["basin_length_definition"]);

    if (basin_length_definition != "flow_path" &&
        basin_length_definition != "geometric")
    {
        std::cerr << "FATAL: basin_length_definition must be one of: "
                  << "flow_path, geometric.\n";
        return 1;
    }

    if (basin_outlet_snapping_method != "none")
    {
        std::cout << "Outlet snapping method                 : "
                  << basin_outlet_snapping_method << "\n";
        std::cout << "Outlet snapping threshold (pixels)     : "
                  << ifs["basin_outlet_contributing_pixels_threshold"] << "\n";

        if (basin_outlet_snapping_method == "expanding_square")
        {
            std::cout << "Outlet snapping min search radius      : "
                      << ifs["basin_outlet_min_search_radius_pixels"] << "\n";
            std::cout << "Outlet snapping max search radius      : "
                      << ifs["basin_outlet_max_search_radius_pixels"] << "\n";
        }
    }

    if (sf["basin_outlet_csv"] == "NULL" || sf["basin_outlet_csv"].empty())
    {
        std::cerr << "FATAL: basin_outlet_csv is required.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Resolve output format
    // -------------------------------------------------------------------------
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
    // Load raster
    // -------------------------------------------------------------------------
    const std::string read_path = pp.get_full_read_path();
    std::cout << "Opening DEM: " << read_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(read_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    // -------------------------------------------------------------------------
    // Require projected metres
    // -------------------------------------------------------------------------
    try
    {
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n"
                  << "Use lsdtt-raster-preprocessing first if reprojection is needed.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Load raster into CPU representation
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);

    // -------------------------------------------------------------------------
    // Remove seas
    // -------------------------------------------------------------------------
    if (ff["minimum_elevation"] >= 0.0f)
    {
        std::cout << "\nMasking pixels at or below "
                  << ff["minimum_elevation"] << " m\n";

        dem = lsdtt::mask_to_nodata_using_threshold(
            dem,
            ff["minimum_elevation"],
            true
        );
    }

    // -------------------------------------------------------------------------
    // Optional hillshade
    // -------------------------------------------------------------------------
    if (bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Hillshade ---\n";

        try
        {
            lsdtt::lsdtt_raster_cpu hillshade = lsdtt::compute_hillshade(
                dem,
                ff["hillshade_azimuth"],
                ff["hillshade_altitude"],
                ff["hillshade_z_factor"]
            );

            const std::string o = pp.make_output_path( "_hillshade.tif");
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
    // Prepare routing DEM
    //
    // Statistics are computed from the sea-masked DEM (dem).
    // Flow routing is built on routing_dem, which may optionally be pit-filled.
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu routing_dem = dem;

    if (bf["fill_before_routing"])
    {
        std::cout << "\n--- Pit filling before routing ---\n";

        try
        {
            // Periodic-aware overload. This CLI routes with the user's
            // boundary_conditions (see get_boundary_conditions below), so
            // filling with open edges and then routing periodically leaves
            // cells with no downhill neighbour at the seam and routing aborts
            // with "could not assign a receiver at row N, col 0".
            routing_dem = lsdtt::fill_pits(
                dem,
                ff["fill_min_slope"],
                join_boundary_conditions(pp.get_boundary_conditions()),
                ifs["n_threads"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Build flow routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- Build flow routing ---\n";

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.n_threads = ifs["n_threads"];

    {
        const std::vector<std::string> bc = pp.get_boundary_conditions();
        if (!bc.empty())
        {
            fi_opts.boundary_conditions = join_boundary_conditions(bc);
        }
    }

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try
    {
        flow = lsdtt::make_flow_base(
            routing_dem.get_data_array(),
            routing_dem.get_mask_array(),
            routing_dem.get_meta(),
            fi_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Gradient raster for slope percentiles
    // -------------------------------------------------------------------------
    std::cout << "\n--- Computing slope (gradient magnitude) raster ---\n";

    lsdtt::lsdtt_raster_cpu gradient_raster;
    try
    {
        gradient_raster = lsdtt::compute_gradient_magnitude_3x3_square(dem);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: gradient computation failed: " << e.what() << "\n";
        return 1;
    }

    const lsdtt::lsdtt_array2d_f&  grad_data = gradient_raster.get_data_array();
    const lsdtt::lsdtt_array2d_u8& grad_mask = gradient_raster.get_mask_array();

    // -------------------------------------------------------------------------
    // Flow-distance-from-outlet raster
    // -------------------------------------------------------------------------
    std::cout << "\n--- Computing flow distance raster ---\n";

    lsdtt::lsdtt_raster_cpu flow_dist_raster;
    try
    {
        flow_dist_raster = lsdtt::lsdtt_raster_cpu(
            flow->get_distance_from_outlet_raster()
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow distance computation failed: "
                  << e.what() << "\n";
        return 1;
    }

    const lsdtt::lsdtt_array2d_f&  fd_data = flow_dist_raster.get_data_array();
    const lsdtt::lsdtt_array2d_u8& fd_mask = flow_dist_raster.get_mask_array();

    // -------------------------------------------------------------------------
    // Optional channel network for stream-order statistics
    // -------------------------------------------------------------------------
    const bool do_shape_metrics = bf["calculate_shape_metrics"];
    const bool do_com_bearing   = bf["calculate_centre_of_mass_bearing"];

    std::unique_ptr<lsdtt::lsdtt_junction_network_base> channel_network;
    const bool do_channel_stats         = bf["calculate_channel_stats"];
    const bool do_junction_area_ratios  = bf["calculate_junction_area_ratios"];

    if (do_channel_stats || do_junction_area_ratios)
    {
        std::cout << "\n--- Building channel network "
                  << "(threshold = " << ifs["threshold_contributing_pixels"]
                  << " px) ---\n";

        try
        {
            std::vector<std::int64_t> source_flats =
                flow->get_sources_index_threshold(
                    ifs["threshold_contributing_pixels"]
                );

            if (source_flats.empty())
            {
                std::cerr << "Warning: source extraction returned no pixels. "
                          << "Channel stats will be skipped. "
                          << "Try lowering threshold_contributing_pixels.\n";
            }
            else
            {
                std::cout << "Source pixels: " << source_flats.size() << "\n";

                lsdtt::lsdtt_junction_network_options jn_opts;
                jn_opts.build_topological_order = true;
                jn_opts.build_s_vector          = false;

                channel_network = lsdtt::make_junction_network_base(
                    *flow, source_flats, jn_opts
                );

                std::cout << "Channel cells : "
                          << channel_network->n_channel_cells() << "\n";
                std::cout << "Junctions     : "
                          << channel_network->n_junctions() << "\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: channel network construction failed: "
                      << e.what() << "\n"
                      << "Channel stats will be skipped.\n";
        }
    }

    const bool has_channel_network = (channel_network != nullptr);

    // -------------------------------------------------------------------------
    // Junction drainage-area ratios (network-wide; complement to Tokunaga)
    // -------------------------------------------------------------------------
    if (do_junction_area_ratios && has_channel_network)
    {
        std::cout << "\n--- Junction area-ratio analysis ---\n";
        try
        {
            const lsdtt::lsdtt_junction_area_ratio_result jar =
                lsdtt::compute_junction_area_ratios(
                    *flow, *channel_network,
                    static_cast<std::int32_t>(ifs["junction_area_ratio_bootstrap"]),
                    1u);

            std::cout << "Confluence junctions: " << jar.n_junctions
                      << "  median ratio: " << jar.median_ratio
                      << "  lambda_tail: " << jar.lambda_tail << "\n";

            const std::string f_data =
                pp.make_output_path("_junction_area_ratios.csv");
            std::cout << "Writing: " << f_data << "\n";
            lsdtt::write_junction_area_ratio_data(jar, f_data);

            const std::string f_sum =
                pp.make_output_path("_junction_area_ratio_summary.csv");
            std::cout << "Writing: " << f_sum << "\n";
            lsdtt::write_junction_area_ratio_summary(jar, f_sum);

            const std::string f_dec =
                pp.make_output_path("_junction_area_ratio_decades.csv");
            std::cout << "Writing: " << f_dec << "\n";
            lsdtt::write_junction_area_ratio_decades(jar, f_dec);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: junction area-ratio analysis failed: "
                      << e.what() << "\n";
        }
    }
    else if (do_junction_area_ratios && !has_channel_network)
    {
        std::cerr << "Warning: calculate_junction_area_ratios requested but no "
                     "channel network was built; skipping.\n";
    }

    // -------------------------------------------------------------------------
    // Read basin outlet CSV
    // -------------------------------------------------------------------------
    const std::string basin_csv_path =
        pp.add_read_path_if_needed(sf["basin_outlet_csv"]);

    std::cout << "\nReading basin outlet CSV: " << basin_csv_path << "\n";

    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();

    std::vector<lsdtt::lsdtt_input_point> requested_outlets;
    try
    {
        requested_outlets =
            lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
                basin_csv_path,
                meta
            );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Read " << requested_outlets.size()
              << " outlet point(s)\n";

    // -------------------------------------------------------------------------
    // Read optional gauge_id / gauge_name columns from the outlet CSV
    //
    // The point loader preserves row order, so gauge_id_vec[i] corresponds to
    // requested_outlets[i] (point_id == i+1).
    // -------------------------------------------------------------------------
    std::vector<std::string> gauge_id_vec;
    std::vector<std::string> gauge_name_vec;
    bool has_gauge_id   = false;
    bool has_gauge_name = false;

    {
        std::ifstream csv_in(basin_csv_path);
        if (csv_in)
        {

            std::string line;
            if (std::getline(csv_in, line))
            {
                const auto header = lsdtt::csv_split_row_trimmed(line);
                int gauge_id_col   = -1;
                int gauge_name_col = -1;

                for (int ci = 0; ci < static_cast<int>(header.size()); ++ci)
                {
                    if (header[ci] == "gauge_id")
                    {
                        gauge_id_col   = ci;
                        has_gauge_id   = true;
                    }
                    if (header[ci] == "gauge_name")
                    {
                        gauge_name_col = ci;
                        has_gauge_name = true;
                    }
                }

                while (std::getline(csv_in, line))
                {
                    if (line.empty() || line.find_first_not_of(" \t\r\n") ==
                        std::string::npos)
                    {
                        continue;
                    }

                    const auto fields = lsdtt::csv_split_row_trimmed(line);

                    gauge_id_vec.push_back(
                        (gauge_id_col >= 0 &&
                         gauge_id_col < static_cast<int>(fields.size()))
                            ? fields[gauge_id_col] : "");

                    gauge_name_vec.push_back(
                        (gauge_name_col >= 0 &&
                         gauge_name_col < static_cast<int>(fields.size()))
                            ? fields[gauge_name_col] : "");
                }
            }
        }
    }

    if (has_gauge_id)   std::cout << "  Found gauge_id column\n";
    if (has_gauge_name) std::cout << "  Found gauge_name column\n";

    // -------------------------------------------------------------------------
    // CRS transform for centroid lat/lon output
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_coordinate_transform to_wgs84;
    try
    {
        to_wgs84.initialise_from_epsg(meta.crs_wkt, 4326);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning: CRS transform initialisation failed: "
                  << e.what() << "\n"
                  << "Centroid lat/lon will be 0.\n";
    }

    // -------------------------------------------------------------------------
    // Basin statistics
    // -------------------------------------------------------------------------
    std::vector<basin_stats_record> basin_rows;
    basin_rows.reserve(requested_outlets.size());

    const lsdtt::lsdtt_array2d_f&  dem_data = dem.get_data_array();
    const lsdtt::lsdtt_array2d_u8& dem_mask = dem.get_mask_array();

    const double cell_size     = meta.cell_size;
    const double cell_area_km2 = cell_size * cell_size / 1.0e6;
    const double diag_len      = cell_size * std::sqrt(2.0);

    lsdtt::lsdtt_square_snap_options square_snap_options;
    square_snap_options.contributing_pixels_threshold =
        static_cast<std::uint64_t>(
            ifs["basin_outlet_contributing_pixels_threshold"]
        );
    square_snap_options.min_search_radius_pixels =
        ifs["basin_outlet_min_search_radius_pixels"];
    square_snap_options.max_search_radius_pixels =
        ifs["basin_outlet_max_search_radius_pixels"];

    lsdtt::lsdtt_downstream_snap_options downstream_snap_options;
    downstream_snap_options.contributing_pixels_threshold =
        static_cast<std::uint64_t>(
            ifs["basin_outlet_contributing_pixels_threshold"]
        );

    auto snap_status_to_string = [](lsdtt::lsdtt_snap_status status)
    {
        switch (status)
        {
            case lsdtt::lsdtt_snap_status::success:
                return std::string("success");
            case lsdtt::lsdtt_snap_status::input_out_of_bounds:
                return std::string("input_out_of_bounds");
            case lsdtt::lsdtt_snap_status::input_invalid_flow_cell:
                return std::string("input_invalid_flow_cell");
            case lsdtt::lsdtt_snap_status::threshold_not_found_in_search_window:
                return std::string("threshold_not_found_in_search_window");
            case lsdtt::lsdtt_snap_status::outlet_reached_before_threshold:
                return std::string("outlet_reached_before_threshold");
        }
        return std::string("unknown");
    };


    for (const auto& outlet : requested_outlets)
    {
        const double latitude  = outlet.original_y;
        const double longitude = outlet.original_x;

        // Gauge labels (aligned by 0-based index = point_id - 1)
        const int gauge_row = outlet.point_id - 1;
        const std::string g_id =
            (has_gauge_id && gauge_row >= 0 &&
             gauge_row < static_cast<int>(gauge_id_vec.size()))
                ? gauge_id_vec[gauge_row] : "";
        const std::string g_name =
            (has_gauge_name && gauge_row >= 0 &&
             gauge_row < static_cast<int>(gauge_name_vec.size()))
                ? gauge_name_vec[gauge_row] : "";

        std::cout << "\nOutlet " << outlet.point_id;
        if (!g_id.empty())   std::cout << "  gauge_id=" << g_id;
        if (!g_name.empty()) std::cout << "  gauge_name=" << g_name;
        std::cout << "  lat=" << latitude << " lon=" << longitude << "\n";

        try
        {
            if (outlet.status != lsdtt::lsdtt_point_input_status::success)
            {
                std::cerr << "Warning: outlet " << outlet.point_id
                          << " did not map to an in-bounds raster pixel; skipping.\n";
                continue;
            }

            int row = outlet.row;
            int col = outlet.col;
            std::int64_t flat = outlet.flat;

            lsdtt::lsdtt_xy xy{outlet.x, outlet.y};

            if (basin_outlet_snapping_method == "expanding_square")
            {
                const lsdtt::lsdtt_snap_result snapped =
                    lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(
                        *flow, row, col, square_snap_options
                    );

                if (snapped.status != lsdtt::lsdtt_snap_status::success)
                {
                    std::cerr << "Warning: outlet " << outlet.point_id
                              << " could not be snapped with expanding_square ("
                              << snap_status_to_string(snapped.status)
                              << "); skipping.\n";
                    continue;
                }

                row  = snapped.snapped_row;
                col  = snapped.snapped_col;
                flat = snapped.snapped_flat;
                xy   = lsdtt::cell_center_from_row_col(meta, row, col);

                std::cout << "  Snapped : row=" << row << " col=" << col
                          << " radius=" << snapped.search_radius_pixels
                          << " contributing=" << snapped.snapped_contributing_pixels << "\n";
            }
            else if (basin_outlet_snapping_method == "downstream")
            {
                const lsdtt::lsdtt_snap_result snapped =
                    lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(
                        *flow, row, col, downstream_snap_options
                    );

                if (snapped.status != lsdtt::lsdtt_snap_status::success)
                {
                    std::cerr << "Warning: outlet " << outlet.point_id
                              << " could not be snapped downstream ("
                              << snap_status_to_string(snapped.status)
                              << "); skipping.\n";
                    continue;
                }

                row  = snapped.snapped_row;
                col  = snapped.snapped_col;
                flat = snapped.snapped_flat;
                xy   = lsdtt::cell_center_from_row_col(meta, row, col);

                std::cout << "  Snapped : row=" << row << " col=" << col
                          << " steps=" << snapped.downstream_steps
                          << " contributing=" << snapped.snapped_contributing_pixels << "\n";
            }
            else
            {
                if (!flow->is_valid_cell(flat))
                {
                    std::cerr << "Warning: outlet " << outlet.point_id
                              << " maps to an invalid or nodata flow pixel; skipping.\n";
                    continue;
                }
            }

            const lsdtt::lsdtt_flat_index_slice upslope =
                flow->upslope_nodes(flat);

            // ------------------------------------------------------------------
            // Check whether any upslope pixel borders nodata or the raster edge
            // ------------------------------------------------------------------
            bool nodata_on_edge = false;
            for (const std::int64_t* it = upslope.begin();
                 it != upslope.end() && !nodata_on_edge;
                 ++it)
            {
                const int br = flow->row_from_flat_unchecked(*it);
                const int bc = flow->col_from_flat_unchecked(*it);

                for (int dr = -1; dr <= 1 && !nodata_on_edge; ++dr)
                {
                    for (int dc = -1; dc <= 1 && !nodata_on_edge; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;
                        const int nr = br + dr;
                        const int nc = bc + dc;
                        if (nr < 0 || nr >= meta.n_rows ||
                            nc < 0 || nc >= meta.n_cols)
                        {
                            nodata_on_edge = true;
                        }
                        else if (dem_mask(nr, nc) == 0)
                        {
                            nodata_on_edge = true;
                        }
                    }
                }
            }

            if (nodata_on_edge)
            {
                std::cout << "  Note: basin touches nodata / raster edge\n";
                if (bf["reject_nodata_edge_basins"])
                {
                    std::cout << "  Skipping (reject_nodata_edge_basins = true)\n";
                    continue;
                }
            }

            // ------------------------------------------------------------------
            // Flow distance at the outlet (absolute, from any baselevel outlet)
            // ------------------------------------------------------------------
            double outlet_fd = 0.0;
            if (fd_mask(row, col) != 0)
            {
                outlet_fd = fd_data(row, col);
            }

            // ------------------------------------------------------------------
            // Accumulate per-pixel values
            // ------------------------------------------------------------------
            std::vector<float> elev_values;
            std::vector<float> slope_values;
            std::vector<float> fdist_values;

            double sum_x = 0.0;
            double sum_y = 0.0;

            // Channel length by stream order (1-6); positive = accumulated
            std::map<int, double> chan_len;
            for (int o = 1; o <= 6; ++o) { chan_len[o] = 0.0; }

            elev_values.reserve(upslope.size());
            slope_values.reserve(upslope.size());
            fdist_values.reserve(upslope.size());

            for (const std::int64_t* it = upslope.begin();
                 it != upslope.end();
                 ++it)
            {
                const std::int64_t basin_flat = *it;
                if (basin_flat < 0) continue;

                const int basin_row = flow->row_from_flat_unchecked(basin_flat);
                const int basin_col = flow->col_from_flat_unchecked(basin_flat);

                if (!dem_mask.in_bounds(basin_row, basin_col)) continue;
                if (dem_mask(basin_row, basin_col) == 0) continue;

                // Elevation
                elev_values.push_back(dem_data(basin_row, basin_col));

                // Slope (gradient magnitude)
                if (grad_mask.in_bounds(basin_row, basin_col) &&
                    grad_mask(basin_row, basin_col) != 0)
                {
                    slope_values.push_back(grad_data(basin_row, basin_col));
                }

                // Flow distance from this outlet
                if (fd_mask.in_bounds(basin_row, basin_col) &&
                    fd_mask(basin_row, basin_col) != 0)
                {
                    const double rel_fd =
                        std::max(0.0,
                            static_cast<double>(fd_data(basin_row, basin_col))
                            - outlet_fd);
                    fdist_values.push_back(static_cast<float>(rel_fd));
                }

                // Centroid accumulation
                const lsdtt::lsdtt_xy pxy =
                    lsdtt::cell_center_from_row_col(meta, basin_row, basin_col);
                sum_x += pxy.x;
                sum_y += pxy.y;

                // Channel length by stream order
                if (has_channel_network &&
                    channel_network->is_channel(basin_flat))
                {
                    const std::int32_t junc =
                        channel_network->junction_upstream_of_flat(basin_flat);

                    if (junc >= 0)
                    {
                        const int order = static_cast<int>(
                            channel_network->stream_order_of_junction(junc)
                        );

                        if (order >= 1 && order <= 6)
                        {
                            const std::uint8_t dir =
                                flow->receiver_dir_of(basin_flat);
                            const bool diagonal =
                                (dir == 1 || dir == 3 ||
                                 dir == 5 || dir == 7);
                            chan_len[order] +=
                                diagonal ? diag_len : cell_size;
                        }
                    }
                }
            }

            if (elev_values.empty())
            {
                std::cerr << "Warning: basin for outlet " << outlet.point_id
                          << " contained no valid DEM pixels; skipping.\n";
                continue;
            }

            // ------------------------------------------------------------------
            // Descriptive statistics
            // ------------------------------------------------------------------
            const lsdtt::lsdtt_descriptive_stats elev_stats =
                lsdtt::compute_descriptive_stats(elev_values);

            // Slope percentiles
            double s25 = 0.0, s50 = 0.0, s75 = 0.0;
            if (!slope_values.empty())
            {
                s25 = lsdtt::compute_percentile(slope_values, 25.0);
                s50 = lsdtt::compute_percentile(slope_values, 50.0);
                s75 = lsdtt::compute_percentile(slope_values, 75.0);
            }

            // Flow distance percentiles
            double fd25 = 0.0, fd50 = 0.0, fd75 = 0.0;
            if (!fdist_values.empty())
            {
                fd25 = lsdtt::compute_percentile(fdist_values, 25.0);
                fd50 = lsdtt::compute_percentile(fdist_values, 50.0);
                fd75 = lsdtt::compute_percentile(fdist_values, 75.0);
            }

            // Centroid
            const std::size_t n_pix = elev_values.size();
            const double cx = sum_x / static_cast<double>(n_pix);
            const double cy = sum_y / static_cast<double>(n_pix);

            double cent_lat = 0.0, cent_lon = 0.0;
            if (to_wgs84.is_initialised())
            {
                const lsdtt::lsdtt_lonlat ll =
                    lsdtt::xy_to_lonlat(cx, cy, to_wgs84);
                cent_lat = ll.latitude;
                cent_lon = ll.longitude;
            }

            // PDF/CDF per basin
            if (bf["write_basin_pdf_cdf_csv"])
            {
                const lsdtt::lsdtt_histogram hist =
                    lsdtt::compute_histogram(
                        elev_values,
                        static_cast<std::size_t>(ifs["n_histogram_bins"])
                    );

                const std::string dist_csv =
                    basin_pdf_cdf_path(pp, outlet.point_id);
                std::cout << "  Writing PDF/CDF: " << dist_csv << "\n";
                write_basin_pdf_cdf_csv(dist_csv, outlet.point_id, hist);
            }

            // ------------------------------------------------------------------
            // Basin planform shape metrics
            //
            // Trace the basin outline polygon, convert the pixel-corner divide
            // ring to projected coordinates, and derive shape metrics.
            // See docs/tex_technical/basin_shape_metrics.tex for derivations and references.
            // ------------------------------------------------------------------
            lsdtt::lsdtt_basin_shape shape;
            bool shape_ok = false;
            if (do_shape_metrics)
            try
            {
                lsdtt::lsdtt_basin_cpu basin(flat, *flow);
                lsdtt::lsdtt_basin_outline_cpu outline = basin.get_basin_outline();

                const std::vector<lsdtt::lsdtt_pixel_rc>& div_rc =
                    outline.outline_divide_rc();
                const lsdtt::lsdtt_raster_meta& dgm = outline.divide_grid_meta();

                std::vector<lsdtt::lsdtt_xy> ring;
                ring.reserve(div_rc.size());
                for (const auto& rc : div_rc)
                {
                    ring.push_back(
                        lsdtt::cell_center_from_row_col(dgm, rc.row, rc.col));
                }

                if (ring.size() >= 3)
                {
                    double basin_len = 0.0;
                    if (basin_length_definition == "geometric")
                    {
                        const std::vector<lsdtt::lsdtt_xy> hull =
                            lsdtt::convex_hull(ring);
                        basin_len = lsdtt::feret_diameters(hull).d_max;
                    }
                    else  // flow_path
                    {
                        for (const float f : fdist_values)
                        {
                            basin_len = std::max(basin_len,
                                                 static_cast<double>(f));
                        }
                    }

                    const double simp_tol =
                        ff["shape_perimeter_simplify_cells"] * cell_size;

                    shape = lsdtt::compute_basin_shape(
                        ring,
                        basin_len,
                        elev_stats.min,
                        elev_stats.mean,
                        elev_stats.max,
                        simp_tol);
                    shape_ok = true;
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "  Warning: basin shape metrics failed for outlet "
                          << outlet.point_id << ": " << e.what() << "\n";
            }

            // ------------------------------------------------------------------
            // Assemble record
            // ------------------------------------------------------------------
            basin_stats_record row_out;
            row_out.basin_id           = outlet.point_id;
            row_out.gauge_id           = g_id;
            row_out.gauge_name         = g_name;
            row_out.latitude           = latitude;
            row_out.longitude          = longitude;
            row_out.outlet_x           = xy.x;
            row_out.outlet_y           = xy.y;
            row_out.outlet_row         = row;
            row_out.outlet_col         = col;
            row_out.outlet_flat        = flat;
            row_out.n_pixels           = n_pix;
            row_out.drainage_area_km2  = static_cast<double>(n_pix) * cell_area_km2;
            row_out.nodata_on_edge     = nodata_on_edge;
            row_out.centroid_x         = cx;
            row_out.centroid_y         = cy;
            row_out.centroid_latitude  = cent_lat;
            row_out.centroid_longitude = cent_lon;
            if (do_com_bearing)
            {
                row_out.com_bearing_deg =
                    lsdtt::bearing_from_north_degrees(xy.x, xy.y, cx, cy);
                row_out.com_distance_m = std::hypot(cx - xy.x, cy - xy.y);
            }
            row_out.stats              = elev_stats;
            row_out.slope_25th         = s25;
            row_out.slope_50th         = s50;
            row_out.slope_75th         = s75;
            row_out.flow_distance_25th = fd25;
            row_out.flow_distance_50th = fd50;
            row_out.flow_distance_75th = fd75;

            row_out.shape_valid           = shape_ok;
            row_out.basin_length_m        = shape.basin_length;
            row_out.elongation_ratio      = shape.elongation_ratio;
            row_out.form_factor           = shape.form_factor;
            row_out.circularity_ratio     = shape.circularity_ratio;
            row_out.gravelius_compactness = shape.gravelius_compactness;
            row_out.solidity              = shape.solidity;
            row_out.eccentricity          = shape.eccentricity;
            row_out.orientation_deg       = shape.orientation_deg;
            row_out.hypsometric_integral  = shape.hypsometric_integral;

            if (has_channel_network)
            {
                row_out.channel_length_by_order = chan_len;
            }

            basin_rows.push_back(row_out);

            std::cout << "  Basin pixels  : " << n_pix << "\n";
            std::cout << "  Drainage area : " << row_out.drainage_area_km2
                      << " km2\n";
            std::cout << "  Elev mean     : " << elev_stats.mean << " m\n";
            std::cout << "  Elev median   : " << elev_stats.median << " m\n";
            std::cout << "  Relief        : " << elev_stats.range << " m\n";
            std::cout << "  Slope 50th    : " << s50 << " m/m\n";
            std::cout << "  FlowDist 50th : " << fd50 << " m\n";
            if (shape_ok)
            {
                std::cout << "  Elongation Re : " << shape.elongation_ratio
                          << "  (L_b=" << shape.basin_length << " m)\n";
                std::cout << "  Circularity Rc: " << shape.circularity_ratio
                          << "  Solidity: " << shape.solidity << "\n";
            }
            if (has_channel_network)
            {
                std::cout << "  Channel len O1: " << chan_len[1] << " m\n";
                std::cout << "  Channel len O2: " << chan_len[2] << " m\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to process outlet "
                      << outlet.point_id << ": " << e.what() << "\n";
        }
    }


    if (basin_rows.empty())
    {
        std::cerr << "\nFATAL: no basin statistics were produced.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Write basin statistics CSV
    // -------------------------------------------------------------------------
    try
    {
        const std::string csv_out = pp.make_output_path( "_basin_stats.csv");
        std::cout << "\nWriting basin statistics CSV: " << csv_out << "\n";
        write_basin_stats_csv(csv_out, basin_rows, has_channel_network,
                              do_shape_metrics, do_com_bearing);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: failed to write basin statistics CSV: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Optional basin-polygon vector output (metrics attached to basin polygons)
    // -------------------------------------------------------------------------
    if (bf["write_basin_polygons"])
    {
        try
        {
            const lsdtt::vector_output_format poly_fmt =
                lsdtt::resolve_vector_output_format(sf["basin_polygon_format"]);
            const std::string ext = lsdtt::output_extension(poly_fmt);
            const std::string poly_out =
                pp.make_output_path("_basin_polygons" + ext);
            std::cout << "Writing basin polygons: " << poly_out << "\n";
            write_basin_polygons(poly_out, basin_rows, *flow, meta, poly_fmt,
                                 has_channel_network, do_shape_metrics,
                                 do_com_bearing, true);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to write basin polygons: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Echo parameters for reproducibility
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-basin-stats: done.\n";
    return 0;
}
