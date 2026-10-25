// lsdtt-basin-stats.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-basin-stats
//
// Reads a DEM and a CSV of basin outlet points, extracts basins, and computes
// a suite of hydromorphic statistics for each basin.
//
// Basin outlet CSV requirements:
//   - parameter: basin_outlet_fname
//   - must contain columns named: latitude, longitude
//   - optional columns: gauge_id, gauge_name
//   - coordinates are interpreted as WGS84 / EPSG:4326 decimal degrees
//
// Workflow:
//   1.  Load raster
//   2.  Check raster is projected in metres
//   3.  Remove seas by masking pixels at or below min_elevation to nodata
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
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>


#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_true_north.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
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
#include "../src/utils/lsdtt_string_utils.hpp"
#include "../src/cpu/lsdtt_basin_stats_cpu.hpp"


static const std::string k_status   = "stable";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-basin-stats";

// basin_stats_record, basin_stats_options, compute_basin_stats and
// basin_touches_nodata now live in src/cpu/lsdtt_basin_stats_cpu.hpp. The
// accumulation loop that used to sit in main() went with them; what stays here
// is the CSV and polygon schemas, the progress output, and the gauge labels,
// which come from optional columns of this tool's own outlet CSV.
using lsdtt::basin_stats_record;
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
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
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
// basin_stats_schema / basin_stats_values
//
// The per-basin statistics as typed fields, in table column order. Both the
// point layer and the CSV are written from these, so their columns cannot
// drift apart. Doubles reach the CSV through csv_number at 17 significant
// digits (exact round trip); the table used to be streamed at the default 6,
// which put latitude/longitude ~5 m from outlet_x/outlet_y.
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_vector_schema basin_stats_schema(
    bool has_channel_stats,
    bool has_shape_metrics,
    bool has_com_bearing)
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
    for (const char* f : {"latitude", "longitude", "outlet_x", "outlet_y"})
        add(f, FT::double_type);
    for (const char* f : {"outlet_row", "outlet_col", "outlet_flat", "n_pixels"})
        add(f, FT::int64_type);
    add("drainage_area_km2", FT::double_type);
    add("nodata_on_edge", FT::int64_type);
    for (const char* f : {"centroid_x", "centroid_y", "centroid_latitude", "centroid_longitude"})
        add(f, FT::double_type);
    if (has_com_bearing)
    {
        add("com_bearing_deg", FT::double_type);
        add("com_distance_m", FT::double_type);
    }
    for (const char* f : {"elevation_min", "elevation_max", "elevation_range", "elevation_mean",
                          "elevation_25th", "elevation_50th", "elevation_75th", "elevation_iqr",
                          "elevation_stddev", "elevation_mad", "slope_25th", "slope_50th",
                          "slope_75th", "flow_distance_25th", "flow_distance_50th",
                          "flow_distance_75th"})
        add(f, FT::double_type);
    if (has_shape_metrics)
    {
        for (const char* f : {"basin_length_m", "elongation_ratio", "form_factor",
                              "circularity_ratio", "gravelius_compactness", "solidity",
                              "eccentricity", "orientation_deg", "hypsometric_integral"})
            add(f, FT::double_type);
    }
    if (has_channel_stats)
    {
        for (int o = 1; o <= 6; ++o)
            add("channel_length_O" + std::to_string(o), FT::double_type);
    }
    return sch;
}

static std::vector<lsdtt::lsdtt_vector_field_value> basin_stats_values(
    const basin_stats_record& row,
    bool has_channel_stats,
    bool has_shape_metrics,
    bool has_com_bearing)
{
    std::vector<lsdtt::lsdtt_vector_field_value> v;
    v.emplace_back(static_cast<std::int64_t>(row.basin_id));
    v.emplace_back(row.gauge_id);
    v.emplace_back(row.gauge_name);
    for (double d : {row.latitude, row.longitude, row.outlet_x, row.outlet_y})
        v.emplace_back(d);
    for (std::int64_t i : {static_cast<std::int64_t>(row.outlet_row),
                           static_cast<std::int64_t>(row.outlet_col),
                           static_cast<std::int64_t>(row.outlet_flat),
                           static_cast<std::int64_t>(row.n_pixels)})
        v.emplace_back(i);
    v.emplace_back(row.drainage_area_km2);
    v.emplace_back(static_cast<std::int64_t>(row.nodata_on_edge ? 1 : 0));
    for (double d : {row.centroid_x, row.centroid_y, row.centroid_latitude,
                     row.centroid_longitude})
        v.emplace_back(d);
    if (has_com_bearing)
    {
        v.emplace_back(row.com_bearing_deg);
        v.emplace_back(row.com_distance_m);
    }
    for (double d : {row.stats.min, row.stats.max, row.stats.range, row.stats.mean,
                     row.stats.q25, row.stats.median, row.stats.q75, row.stats.iqr,
                     row.stats.stddev_sample, row.stats.median_absolute_deviation,
                     row.slope_25th, row.slope_50th, row.slope_75th,
                     row.flow_distance_25th, row.flow_distance_50th,
                     row.flow_distance_75th})
        v.emplace_back(d);
    if (has_shape_metrics)
    {
        for (double d : {row.basin_length_m, row.elongation_ratio, row.form_factor,
                         row.circularity_ratio, row.gravelius_compactness, row.solidity,
                         row.eccentricity, row.orientation_deg, row.hypsometric_integral})
            v.emplace_back(d);
    }
    if (has_channel_stats)
    {
        for (int o = 1; o <= 6; ++o)
        {
            auto it = row.channel_length_by_order.find(o);
            v.emplace_back(it != row.channel_length_by_order.end() ? it->second : 0.0);
        }
    }
    return v;
}

// -----------------------------------------------------------------------------
// write_basin_stats_table
//
// One record per basin: a point at the basin outlet (outlet_x, outlet_y, in
// the DEM's CRS) for FlatGeobuf / GeoPackage, or a CSV row.
// -----------------------------------------------------------------------------
static void write_basin_stats_table(
    const std::string&                     out_path,
    const std::vector<basin_stats_record>& rows,
    lsdtt::point_output_format             fmt,
    const std::string&                     crs_wkt,
    bool                                   has_channel_stats,
    bool                                   has_shape_metrics,
    bool                                   has_com_bearing)
{
    const lsdtt::lsdtt_vector_schema sch =
        basin_stats_schema(has_channel_stats, has_shape_metrics, has_com_bearing);
    lsdtt::lsdtt_data_loader writer;
    if (fmt == lsdtt::point_output_format::csv)
    {
        std::vector<std::vector<std::string>> cells;
        cells.reserve(rows.size());
        for (const auto& row : rows)
        {
            std::vector<std::string> line;
            for (const auto& value :
                 basin_stats_values(row, has_channel_stats, has_shape_metrics, has_com_bearing))
                line.push_back(lsdtt::vector_field_value_to_string(value));
            cells.push_back(std::move(line));
        }
        writer.write_csv(out_path, sch.field_names, cells);
        return;
    }
    // No spatial index: FlatGeobuf would reorder the features along a Hilbert
    // curve, and the table is kept in basin_id order, as the CSV is.
    writer.begin_vector_write(out_path, "basin_stats", sch, crs_wkt,
                              lsdtt::convert_point_to_vector_output_format(fmt), false);
    for (const auto& row : rows)
        writer.write_point(row.outlet_x, row.outlet_y,
                           basin_stats_values(row, has_channel_stats, has_shape_metrics,
                                              has_com_bearing));
    writer.finish_vector_write();
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
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

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
    string_default_map["basin_outlet_fname"] = "NULL";
    help_map["basin_outlet_fname"] = {"string","NULL",
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
        "Unlike extraction workflows, basin statistics defaults to none so an "
        "authoritative gauge/outlet pixel is not moved implicitly. Allowed values: "
        "none, expanding_square, downstream, largest_area. "
        "none: use the mapped raster pixel directly. "
        "expanding_square: search outward ring-by-ring until a qualifying pixel "
        "is found. downstream: follow the D8 receiver chain until a qualifying "
        "pixel is found. largest_area: take the pixel with the most "
        "contributing pixels in the square of half-width "
        "basin_outlet_snap_max_search_radius_pixels, if it is a qualifying "
        "pixel (a large radius can jump onto a bigger neighbouring river).",
        "Basin outlet snapping", "param", ""
    };

    int_default_map["basin_outlet_contributing_pixels_threshold"] =
        lsdtt::default_basin_outlet_contributing_pixels_threshold;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int",
        std::to_string(lsdtt::default_basin_outlet_contributing_pixels_threshold),
        "Contributing-pixel threshold used for optional basin outlet snapping.",
        "This is the total number of upslope contributing pixels from the "
        "flow-routing object, not the number of direct donors to a pixel.",
        "Basin outlet snapping", "param", ""
    };

    int_default_map["basin_outlet_snap_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_snap_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "0 means the mapped input pixel is tested first.",
        "Basin outlet snapping", "param", ""
    };

    int_default_map["basin_outlet_snap_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_snap_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square or "
        "largest_area. For expanding_square the search stops at the first radius "
        "containing one or more qualifying pixels; for largest_area it is the "
        "half-width of the square searched.",
        "Basin outlet snapping", "param", ""
    };

    // --- Sea removal ---
    float_default_map["min_elevation"] = 0.0f;
    help_map["min_elevation"] = {"float","0.0",
        "Pixels at or below this elevation are set to nodata. -9999 disables.",
        "Default of 0.0 removes sea and ocean pixels from many DEM products. "
        "Negative values apply as given: for a basin below sea level (Death "
        "Valley -86 m, Dead Sea -430 m) set it just below the lowest land. "
        "-9999 (or lower) disables masking.",
        "Sea removal", "param", ""};

    // --- Routing preparation ---
    lsdtt::declare_fill_before_routing_parameter(bool_default_map, help_map);

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope (m/m) imposed during pit filling or flat resolution.",
        "Used only when fill_before_routing is true. inline_gradient "
        "(default): imposed during the Priority-Flood fill. "
        "barnes2014b_gradient: the minimum gradient across flats. "
        "barnes2014b_float_step: ignored.",
        "Routing preparation", "param", ""};

    lsdtt::declare_fill_flat_resolution_parameter(string_default_map, help_map);

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
        "Sun azimuth in degrees clockwise from grid north (grid up).",
        "315 is the standard cartographic northwest illumination. It is "
        "relative to the raster grid, so in an oblique frame it lights from the "
        "frame's north-west; bearings_relative_to_true_north never changes it.",
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

    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

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
        "Format set by write_vector_format.",
        "Channel statistics", "switch", ""};

    string_default_map["basin_stats_write_vector_format"] = "flatgeobuf";
    help_map["basin_stats_write_vector_format"] = {"string","flatgeobuf",
        "Format of the basin statistics table: flatgeobuf, gpkg or csv.",
        "flatgeobuf (default) writes <write_fname>_basin_stats.fgb, one point per "
        "basin at its outlet in the DEM's CRS with every statistic as an "
        "attribute; gpkg the same as GeoPackage; csv writes "
        "<write_fname>_basin_stats.csv, where latitude/longitude and "
        "outlet_x/outlet_y locate each basin. Numbers are written at full "
        "double precision.",
        "Channel statistics", "param", ""};

    string_default_map["write_vector_format"] = "flatgeobuf";
    help_map["write_vector_format"] = {"string","flatgeobuf",
        "Vector format for the basin-polygon output.",
        "flatgeobuf (default) or gpkg.",
        "Channel statistics", "param", ""};

    bool_default_map["calculate_centre_of_mass_bearing"] = false;
    help_map["calculate_centre_of_mass_bearing"] = {"bool","false",
        "Compute the bearing from each basin outlet to its centre of mass.",
        "Opt-in. The basin centre of mass is the planimetric centroid of all "
        "upslope pixels (already reported as centroid_x/centroid_y). This adds "
        "com_bearing_deg (bearing clockwise from grid north, or from true north "
        "with bearings_relative_to_true_north, outlet -> centre of mass) and "
        "com_distance_m columns. The outlet is the user-"
        "supplied (optionally snapped) outlet pixel; the channel is not followed.",
        "Channel statistics", "switch", ""};

    // --- Basin shape metrics ---
    bool_default_map["calculate_shape_metrics"] = false;
    help_map["calculate_shape_metrics"] = {"bool","false",
        "Trace each basin outline polygon and compute planform shape metrics.",
        "Opt-in; off by default so the original curated hydromorphic stat set is "
        "reproduced unless explicitly requested. When true, adds "
        "basin_length_m, elongation_ratio, form_factor, circularity_ratio, "
        "gravelius_compactness, solidity, eccentricity, orientation_deg (principal "
        "axis, [0, 180) clockwise from grid north, or from true north with "
        "bearings_relative_to_true_north), and "
        "hypsometric_integral columns to the basin statistics table. When false, the outline "
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

    pp.set_deprecated_aliases({
        {"basin_outlet_csv", "basin_outlet_fname"},
        {"basin_outlet_min_search_radius_pixels", "basin_outlet_snap_min_search_radius_pixels"},
        {"basin_outlet_max_search_radius_pixels", "basin_outlet_snap_max_search_radius_pixels"},
        {"basin_polygon_format", "write_vector_format"},
        {"minimum_elevation", "min_elevation"} });
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const std::string basin_outlet_snapping_method =
        lsdtt::to_lower_ascii(sf["basin_outlet_snapping_method"]);

    if (basin_outlet_snapping_method != "none" &&
        basin_outlet_snapping_method != "expanding_square" &&
        basin_outlet_snapping_method != "downstream" &&
        basin_outlet_snapping_method != "largest_area")
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                  << "none, expanding_square, downstream, largest_area.\n";
        return 1;
    }

    if (ifs["basin_outlet_contributing_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: basin_outlet_contributing_pixels_threshold "
                  << "must be > 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_snap_min_search_radius_pixels"] < 0)
    {
        std::cerr << "FATAL: basin_outlet_snap_min_search_radius_pixels "
                  << "must be >= 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_snap_max_search_radius_pixels"] <
        ifs["basin_outlet_snap_min_search_radius_pixels"])
    {
        std::cerr << "FATAL: basin_outlet_snap_max_search_radius_pixels must be >= "
                  << "basin_outlet_snap_min_search_radius_pixels.\n";
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
        lsdtt::to_lower_ascii(sf["basin_length_definition"]);

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
                      << ifs["basin_outlet_snap_min_search_radius_pixels"] << "\n";
            std::cout << "Outlet snapping max search radius      : "
                      << ifs["basin_outlet_snap_max_search_radius_pixels"] << "\n";
        }
        else if (basin_outlet_snapping_method == "largest_area")
        {
            std::cout << "Outlet snapping search half-width      : "
                      << ifs["basin_outlet_snap_max_search_radius_pixels"] << "\n";
        }
    }

    if (sf["basin_outlet_fname"] == "NULL" || sf["basin_outlet_fname"].empty())
    {
        std::cerr << "FATAL: basin_outlet_fname is required.\n";
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

    // The basin statistics table: a point layer (default) or CSV.
    const std::string& stats_fmt_str = sf["basin_stats_write_vector_format"];
    if (stats_fmt_str != "flatgeobuf" && stats_fmt_str != "fgb" &&
        stats_fmt_str != "gpkg" && stats_fmt_str != "csv")
    {
        std::cerr << "FATAL: basin_stats_write_vector_format must be flatgeobuf, gpkg or csv "
                     "(got '" << stats_fmt_str << "').\n";
        return 1;
    }
    const lsdtt::point_output_format stats_fmt =
        lsdtt::resolve_point_output_format(stats_fmt_str);

    // How flats are made drainable (Barnes et al. 2014b options). Checked
    // here, before any raster work. Only a fill uses it: with
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
                flat_choice, ff["fill_min_slope"]);
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

    // com_bearing_deg and orientation_deg: grid-relative unless
    // bearings_relative_to_true_north.
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

    // -------------------------------------------------------------------------
    // Remove seas
    // -------------------------------------------------------------------------
    // -9999 (or lower) disables; any other value, negative included, applies.
    if (ff["min_elevation"] > -9998.0f)
    {
        std::cout << "\nMasking pixels at or below "
                  << ff["min_elevation"] << " m\n";

        dem = lsdtt::mask_to_nodata_using_threshold(
            dem,
            ff["min_elevation"],
            true
        );
    }

    // Interior nodata: kept as sinks or filled (interior_nodata_method).
    try
    {
        std::cout << lsdtt::apply_interior_nodata_method(
                         dem, sf.at("interior_nodata_method"),
                         sf.at("nodata_fill_method_choice"),
                         ifs.at("nodata_fill_max_iter"), ifs.at("n_threads"))
                  << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: interior nodata: " << e.what() << "\n";
        return 1;
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
            if (flat_choice == lsdtt::flat_resolution_choice::inline_gradient)
            {
                routing_dem = lsdtt::fill_pits(
                    dem,
                    ff["fill_min_slope"],
                    join_boundary_conditions(pp.get_boundary_conditions()),
                    ifs["n_threads"]
                );
            }
            else
            {
                routing_dem = lsdtt::fill_pits_reporting_flat_resolution(
                    dem,
                    ff["fill_min_slope"],
                    join_boundary_conditions(pp.get_boundary_conditions()),
                    flat_choice,
                    ifs["n_threads"]
                );
            }
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

        // PERIODIC BOUNDARIES ON A MEASURED DEM ARE ALMOST CERTAINLY A MISTAKE.
        //
        // They exist for numerical models, where the domain is meant to wrap.
        // A real DEM does not wrap: a basin that leaves one side of the raster
        // does not reappear on the other. Under a periodic boundary the flow
        // routing joins the two sides anyway, so upslope areas, basin extents
        // and the completeness test all silently describe a catchment that does
        // not exist on the ground.
        //
        // Warn rather than refuse: a synthetic DEM from lsdtt-lem or
        // lsdtt-synthetic-networks is a legitimate input to this tool, and for
        // those the periodic boundary is the point.
        if (fi_opts.boundary_conditions.find('p') != std::string::npos ||
            fi_opts.boundary_conditions.find('P') != std::string::npos)
        {
            std::cerr
                << "\nWARNING: boundary_conditions is '"
                << fi_opts.boundary_conditions
                << "', which makes at least one edge PERIODIC.\n"
                << "         Periodic boundaries are for numerical models. On a "
                   "measured DEM they are\n"
                << "         almost always a mistake: flow wraps from one side "
                   "of the raster to the other,\n"
                << "         so drainage areas, basin extents and the "
                   "basin-completeness test describe a\n"
                << "         catchment that does not exist on the ground. Use "
                   "'o' (open) or 'n' (no flux)\n"
                << "         unless this DEM came from a model that wraps.\n\n";
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
        pp.add_read_path_if_needed(sf["basin_outlet_fname"]);

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

    lsdtt::lsdtt_square_snap_options square_snap_options;
    square_snap_options.contributing_pixels_threshold =
        static_cast<std::uint64_t>(
            ifs["basin_outlet_contributing_pixels_threshold"]
        );
    square_snap_options.min_search_radius_pixels =
        ifs["basin_outlet_snap_min_search_radius_pixels"];
    square_snap_options.max_search_radius_pixels =
        ifs["basin_outlet_snap_max_search_radius_pixels"];

    lsdtt::lsdtt_downstream_snap_options downstream_snap_options;
    downstream_snap_options.contributing_pixels_threshold =
        static_cast<std::uint64_t>(
            ifs["basin_outlet_contributing_pixels_threshold"]
        );



    lsdtt::basin_stats_options stats_options;
    stats_options.snapping_method            = basin_outlet_snapping_method;
    stats_options.square_snap                = square_snap_options;
    stats_options.downstream_snap            = downstream_snap_options;
    stats_options.reject_nodata_edge_basins  = bf["reject_nodata_edge_basins"];
    stats_options.calculate_centre_of_mass_bearing = do_com_bearing;
    stats_options.true_north                 = true_north.get();
    stats_options.calculate_shape_metrics    = do_shape_metrics;
    stats_options.basin_length_definition    = basin_length_definition;
    stats_options.shape_perimeter_simplify_cells =
        ff["shape_perimeter_simplify_cells"];
    // The PDF/CDF writer needs the basin's elevation values, not the record.
    stats_options.collect_elevation_values   = bf["write_basin_pdf_cdf_csv"];

    const std::vector<lsdtt::basin_stats_outlet_result> outlet_results =
        lsdtt::compute_basin_stats(
            *flow, dem, gradient_raster, flow_dist_raster,
            channel_network.get(), to_wgs84,
            requested_outlets, stats_options);

    // Reporting is this CLI's job: compute_basin_stats does not print.
    for (std::size_t i = 0; i < outlet_results.size(); ++i)
    {
        const lsdtt::basin_stats_outlet_result& res = outlet_results[i];
        const lsdtt::lsdtt_input_point&      outlet = requested_outlets[i];

        // Gauge labels, aligned by 0-based index = point_id - 1. They come from
        // optional columns of this tool's outlet CSV, so they are filled in
        // here rather than by the compute layer.
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
        std::cout << "  lat=" << outlet.original_y
                  << " lon=" << outlet.original_x << "\n";

        if (res.snap_ran && res.snap.status == lsdtt::lsdtt_snap_status::success)
        {
            std::cout << "  Snapped : row=" << res.snap.snapped_row
                      << " col=" << res.snap.snapped_col;
            if (basin_outlet_snapping_method == "expanding_square" ||
                basin_outlet_snapping_method == "largest_area")
                std::cout << " radius=" << res.snap.search_radius_pixels;
            else
                std::cout << " steps=" << res.snap.downstream_steps;
            std::cout << " contributing="
                      << res.snap.snapped_contributing_pixels << "\n";
        }

        if (res.nodata_on_edge)
        {
            std::cout << "  Note: basin touches nodata / raster edge\n";
        }

        switch (res.outcome)
        {
            case lsdtt::basin_stats_outcome::point_not_mapped:
                std::cerr << "Warning: outlet " << outlet.point_id
                          << " did not map to an in-bounds raster pixel; skipping.\n";
                continue;
            case lsdtt::basin_stats_outcome::snap_failed:
                std::cerr << "Warning: outlet " << outlet.point_id
                          << " could not be snapped with "
                          << basin_outlet_snapping_method << " ("
                          << lsdtt::snap_status_to_string(res.snap.status)
                          << "); skipping.\n";
                continue;
            case lsdtt::basin_stats_outcome::invalid_flow_pixel:
                std::cerr << "Warning: outlet " << outlet.point_id
                          << " maps to an invalid or nodata flow pixel; skipping.\n";
                continue;
            case lsdtt::basin_stats_outcome::rejected_nodata_edge:
                std::cout << "  Skipping (reject_nodata_edge_basins = true)\n";
                continue;
            case lsdtt::basin_stats_outcome::no_valid_pixels:
                std::cerr << "Warning: basin for outlet " << outlet.point_id
                          << " contained no valid DEM pixels; skipping.\n";
                continue;
            case lsdtt::basin_stats_outcome::failed:
                std::cerr << "Warning: failed to process outlet "
                          << outlet.point_id << ": " << res.message << "\n";
                continue;
            case lsdtt::basin_stats_outcome::success:
                break;
        }

        if (res.shape_failed)
        {
            std::cerr << "  Warning: basin shape metrics failed for outlet "
                      << outlet.point_id << ": " << res.message << "\n";
        }

        basin_stats_record row_out = res.record;
        row_out.gauge_id   = g_id;
        row_out.gauge_name = g_name;

        if (bf["write_basin_pdf_cdf_csv"])
        {
            const lsdtt::lsdtt_histogram hist =
                lsdtt::compute_histogram(
                    res.elevation_values,
                    static_cast<std::size_t>(ifs["n_histogram_bins"]));
            const std::string dist_csv = basin_pdf_cdf_path(pp, outlet.point_id);
            std::cout << "  Writing PDF/CDF: " << dist_csv << "\n";
            write_basin_pdf_cdf_csv(dist_csv, outlet.point_id, hist);
        }

        basin_rows.push_back(row_out);

        std::cout << "  Basin pixels  : " << row_out.n_pixels << "\n";
        std::cout << "  Drainage area : " << row_out.drainage_area_km2 << " km2\n";
        std::cout << "  Elev mean     : " << row_out.stats.mean << " m\n";
        std::cout << "  Elev median   : " << row_out.stats.median << " m\n";
        std::cout << "  Relief        : " << row_out.stats.range << " m\n";
        std::cout << "  Slope 50th    : " << row_out.slope_50th << " m/m\n";
        std::cout << "  FlowDist 50th : " << row_out.flow_distance_50th << " m\n";
        if (row_out.shape_valid)
        {
            std::cout << "  Elongation Re : " << row_out.elongation_ratio
                      << "  (L_b=" << row_out.basin_length_m << " m)\n";
            std::cout << "  Circularity Rc: " << row_out.circularity_ratio
                      << "  Solidity: " << row_out.solidity << "\n";
        }
        if (channel_network != nullptr)
        {
            std::cout << "  Channel len O1: "
                      << row_out.channel_length_by_order.at(1) << " m\n";
            std::cout << "  Channel len O2: "
                      << row_out.channel_length_by_order.at(2) << " m\n";
        }
    }


    if (basin_rows.empty())
    {
        std::cerr << "\nFATAL: no basin statistics were produced.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Write the basin statistics table (point layer or CSV)
    // -------------------------------------------------------------------------
    try
    {
        const std::string stats_out =
            pp.make_output_path("_basin_stats" + lsdtt::output_extension(stats_fmt));
        std::cout << "\nWriting basin statistics: " << stats_out << "\n";
        write_basin_stats_table(stats_out, basin_rows, stats_fmt, meta.crs_wkt,
                                has_channel_network, do_shape_metrics, do_com_bearing);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: failed to write basin statistics: "
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
                lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
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
                  k_status, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-basin-stats: done.\n";
    return 0;
}
