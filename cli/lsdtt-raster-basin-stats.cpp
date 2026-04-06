// lsdtt-raster-basin-stats.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-raster-basin-stats
//
// Reads a DEM and a CSV of user-specified basin outlet points, then computes
// descriptive elevation statistics for the upslope basin draining to each
// outlet.
//
// Basin outlet CSV requirements:
//   - parameter: basin_outlet_csv
//   - must contain columns named: latitude, longitude
//   - coordinates are interpreted as WGS84 / EPSG:4326 decimal degrees
//
// Workflow:
//   1. Load raster
//   2. Check raster is projected in metres
//   3. Remove seas by masking pixels at or below minimum_elevation to nodata
//   4. Optionally write a hillshade of the sea-masked DEM
//   5. Fill pits before routing (default: true)
//   6. Build D8 flow routing
//   7. Read outlet CSV
//   8. Transform each outlet from WGS84 to raster CRS
//   9. Convert to raster pixel
//  10. Get all upslope nodes from flow routing
//  11. Extract DEM elevations for those basin pixels
//  12. Compute descriptive statistics on those valid values
//  13. Write basin statistics to CSV
//
// IMPORTANT:
//   - statistics functions assume the supplied flat data are already valid
//   - this CLI is responsible for excluding nodata before calling statistics
//
// Usage:
//   lsdtt-raster-basin-stats <path> <parameter_file>
//   lsdtt-raster-basin-stats <parameter_file>
//   lsdtt-raster-basin-stats -h
//   lsdtt-raster-basin-stats -v
//   lsdtt-raster-basin-stats -c
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/utils/lsdtt_statistics.hpp"
#include "../src/utils/lsdtt_histogram.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"


static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// Local record types
// -----------------------------------------------------------------------------

struct basin_stats_record
{
    int          basin_id     = -1;
    double       latitude     = 0.0;
    double       longitude    = 0.0;
    double       outlet_x     = 0.0;
    double       outlet_y     = 0.0;
    int          outlet_row   = -1;
    int          outlet_col   = -1;
    std::int64_t outlet_flat  = -1;
    std::size_t  n_pixels     = 0;

    lsdtt::lsdtt_descriptive_stats stats;
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
//
// Writes one row per bin for a single basin.
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
        << "|| Welcome to lsdtt3 raster basin statistics           ||\n"
        << "|| DEM loading, sea masking, routing, basin stats      ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// parse_args
// -----------------------------------------------------------------------------
static std::pair<std::string, std::string> parse_args(int argc, char* argv[])
{
    auto is_flag = [](const std::string& s,
                      std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags)
        {
            if (s == f) return true;
        }
        return false;
    };

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];

        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout
                << "Usage: lsdtt-raster-basin-stats <path> <parameter_file>\n"
                << "       lsdtt-raster-basin-stats <parameter_file>\n"
                << "       lsdtt-raster-basin-stats -h   (this help)\n"
                << "       lsdtt-raster-basin-stats -v   (version)\n"
                << "       lsdtt-raster-basin-stats -c   (citation)\n\n"
                << "Required parameter:\n"
                << "  basin_outlet_csv: CSV file with columns latitude,longitude\n";
            std::exit(0);
        }

        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-raster-basin-stats version "
                      << k_version << "\n";
            std::exit(0);
        }

        if (is_flag(a, {"-c", "--cite", "citation", "-citation",
                        "--citation"}))
        {
            std::cout << "Citation: " << k_citation << "\n";
            std::exit(0);
        }
    }

    std::string path = "./";
    std::string fname;

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        std::exit(1);
    }
    else if (argc == 2)
    {
        fname = argv[1];
    }
    else
    {
        path  = argv[1];
        fname = argv[2];
    }

    return {path, fname};
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
// out_path
// -----------------------------------------------------------------------------
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& suffix)
{
    return pp.get_full_write_path() + suffix;
}

// -----------------------------------------------------------------------------
// write_basin_stats_csv
// -----------------------------------------------------------------------------
static void write_basin_stats_csv(
    const std::string& csv_path,
    const std::vector<basin_stats_record>& rows)
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
        << ",latitude"
        << ",longitude"
        << ",outlet_x"
        << ",outlet_y"
        << ",outlet_row"
        << ",outlet_col"
        << ",outlet_flat"
        << ",n_pixels"
        << ",min"
        << ",max"
        << ",range"
        << ",sum"
        << ",mean"
        << ",median"
        << ",q25"
        << ",q75"
        << ",iqr"
        << ",variance_population"
        << ",variance_sample"
        << ",stddev_population"
        << ",stddev_sample"
        << ",mean_absolute_deviation"
        << ",median_absolute_deviation"
        << ",mad_normal_scaled"
        << "\n";

    for (const auto& row : rows)
    {
        out
            << row.basin_id << ","
            << row.latitude << ","
            << row.longitude << ","
            << row.outlet_x << ","
            << row.outlet_y << ","
            << row.outlet_row << ","
            << row.outlet_col << ","
            << row.outlet_flat << ","
            << row.n_pixels << ","
            << row.stats.min << ","
            << row.stats.max << ","
            << row.stats.range << ","
            << row.stats.sum << ","
            << row.stats.mean << ","
            << row.stats.median << ","
            << row.stats.q25 << ","
            << row.stats.q75 << ","
            << row.stats.iqr << ","
            << row.stats.variance_population << ","
            << row.stats.variance_sample << ","
            << row.stats.stddev_population << ","
            << row.stats.stddev_sample << ","
            << row.stats.mean_absolute_deviation << ","
            << row.stats.median_absolute_deviation << ","
            << row.stats.mad_normal_scaled
            << "\n";
    }
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = parse_args(argc, argv);

    // -------------------------------------------------------------------------
    // Load parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(path, fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = single-threaded."};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format."};

    // --- Basin outlets ---
    string_default_map["basin_outlet_csv"] = "NULL";
    help_map["basin_outlet_csv"] = {"string","NULL",
        "CSV file of requested basin outlet points.",
        "Required. CSV must contain columns named latitude and longitude. "
        "Coordinates are interpreted as WGS84 (EPSG:4326) decimal degrees."};

    // --- Basin outlet snapping ---
    string_default_map["basin_outlet_snapping_method"] = "none";
    help_map["basin_outlet_snapping_method"] = {"string","none",
        "Optional snapping of basin outlet CSV points to pixels that exceed "
        "a contributing-pixel threshold.",
        "Allowed values: none, expanding_square, downstream. "
        "none: use the mapped raster pixel directly. "
        "expanding_square: search outward ring-by-ring until a qualifying pixel "
        "is found. downstream: follow the D8 receiver chain until a qualifying "
        "pixel is found."
    };

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for optional basin outlet snapping.",
        "This is the total number of upslope contributing pixels from the "
        "flow-routing object, not the number of direct donors to a pixel."
    };

    int_default_map["basin_outlet_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "0 means the mapped input pixel is tested first."
    };

    int_default_map["basin_outlet_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "The search stops at the first radius containing one or more qualifying pixels."
    };



    // --- Sea removal ---
    float_default_map["minimum_elevation"] = 0.0f;
    help_map["minimum_elevation"] = {"float","0.0",
        "Pixels at or below this elevation are set to nodata.",
        "Default of 0.0 removes sea and ocean pixels from many DEM products."};

    // --- Routing preparation ---
    bool_default_map["fill_before_routing"] = true;
    help_map["fill_before_routing"] = {"bool","true",
        "Fill pits before constructing flow routing.",
        "Recommended for basin extraction from DEMs."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed during pit filling.",
        "Used when fill_before_routing is true."};

    // --- Hillshade ---
    bool_default_map["write_hillshade_raster"] = false;
    help_map["write_hillshade_raster"] = {"bool","false",
        "Write a hillshade raster after sea masking.",
        "Useful for quick visual QC of the masked DEM."};

    float_default_map["hillshade_azimuth"] = 315.0f;
    help_map["hillshade_azimuth"] = {"float","315",
        "Sun azimuth in degrees clockwise from north.",
        "315 is the standard cartographic northwest illumination."};

    float_default_map["hillshade_altitude"] = 45.0f;
    help_map["hillshade_altitude"] = {"float","45",
        "Sun altitude above horizon in degrees.",
        "45 degrees is a common default."};

    float_default_map["hillshade_z_factor"] = 1.0f;
    help_map["hillshade_z_factor"] = {"float","1.0",
        "Vertical exaggeration factor for hillshade.",
        "Values greater than 1 emphasise relief."};

    // --- Basin distribution output ---
    bool_default_map["write_basin_pdf_cdf_csv"] = false;
    help_map["write_basin_pdf_cdf_csv"] = {"bool","false",
        "Write PDF/CDF CSV output for each basin.",
        "If true, writes one CSV per basin containing histogram counts, "
        "probabilities, PDF, and CDF of basin elevations."};

    int_default_map["n_histogram_bins"] = 100;
    help_map["n_histogram_bins"] = {"int","100",
        "Number of bins used for basin elevation histogram/PDF/CDF output.",
        "Used when write_basin_pdf_cdf_csv is true. Default is 100 bins."};


    // -------------------------------------------------------------------------
    // Parse parameters
    // -------------------------------------------------------------------------
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
                      << fmt_str << "' — defaulting to cog.\n";
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
        std::cout << "Azimuth: "  << ff["hillshade_azimuth"]  << " deg\n";
        std::cout << "Altitude: " << ff["hillshade_altitude"] << " deg\n";
        std::cout << "Z-factor: " << ff["hillshade_z_factor"] << "\n";

        try
        {
            lsdtt::lsdtt_raster_cpu hillshade = lsdtt::compute_hillshade(
                dem,
                ff["hillshade_azimuth"],
                ff["hillshade_altitude"],
                ff["hillshade_z_factor"]
            );

            const std::string o = out_path(pp, "_hillshade.tif");
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
            routing_dem = lsdtt::fill_pits(
                dem,
                ff["fill_min_slope"],
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
    // Basin statistics
    // -------------------------------------------------------------------------
    std::vector<basin_stats_record> basin_rows;
    basin_rows.reserve(requested_outlets.size());

    const lsdtt::lsdtt_array2d_f& dem_data = dem.get_data_array();
    const lsdtt::lsdtt_array2d_u8& dem_mask = dem.get_mask_array();


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

        std::cout << "\nOutlet " << outlet.point_id
                  << "  lat=" << latitude
                  << " lon=" << longitude << "\n";

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
                        *flow,
                        row,
                        col,
                        square_snap_options
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

                std::cout << "  Snapped      : row=" << row
                          << " col=" << col
                          << " radius=" << snapped.search_radius_pixels
                          << " contributing_pixels="
                          << snapped.snapped_contributing_pixels << "\n";
            }
            else if (basin_outlet_snapping_method == "downstream")
            {
                const lsdtt::lsdtt_snap_result snapped =
                    lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(
                        *flow,
                        row,
                        col,
                        downstream_snap_options
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

                std::cout << "  Snapped      : row=" << row
                          << " col=" << col
                          << " steps=" << snapped.downstream_steps
                          << " contributing_pixels="
                          << snapped.snapped_contributing_pixels << "\n";
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

            std::vector<float> basin_values;
            basin_values.reserve(upslope.size());

            for (const std::int64_t* it = upslope.begin();
                 it != upslope.end();
                 ++it)
            {
                const std::int64_t basin_flat = *it;

                if (basin_flat < 0)
                {
                    continue;
                }

                const int basin_row = flow->row_from_flat_unchecked(basin_flat);
                const int basin_col = flow->col_from_flat_unchecked(basin_flat);

                if (!dem_mask.in_bounds(basin_row, basin_col))
                {
                    continue;
                }

                // Statistics expects valid flat data only.
                // Keep this defensive check even though flow routing only stores
                // valid cells.
                if (dem_mask(basin_row, basin_col) == 0)
                {
                    continue;
                }

                basin_values.push_back(dem_data(basin_row, basin_col));
            }

            if (basin_values.empty())
            {
                std::cerr << "Warning: basin for outlet " << outlet.point_id
                          << " contained no valid DEM pixels; skipping.\n";
                continue;
            }

            const lsdtt::lsdtt_descriptive_stats stats =
                lsdtt::compute_descriptive_stats(basin_values);

            if (bf["write_basin_pdf_cdf_csv"])
            {
                const lsdtt::lsdtt_histogram hist =
                    lsdtt::compute_histogram(
                        basin_values,
                        static_cast<std::size_t>(ifs["n_histogram_bins"])
                    );

                const std::string dist_csv =
                    basin_pdf_cdf_path(pp, outlet.point_id);

                std::cout << "  Writing PDF/CDF: " << dist_csv << "\n";

                write_basin_pdf_cdf_csv(
                    dist_csv,
                    outlet.point_id,
                    hist
                );
            }

            basin_stats_record row_out;
            row_out.basin_id    = outlet.point_id;
            row_out.latitude    = latitude;
            row_out.longitude   = longitude;
            row_out.outlet_x    = xy.x;
            row_out.outlet_y    = xy.y;
            row_out.outlet_row  = row;
            row_out.outlet_col  = col;
            row_out.outlet_flat = flat;
            row_out.n_pixels    = basin_values.size();
            row_out.stats       = stats;

            basin_rows.push_back(row_out);

            std::cout << "  Basin pixels : " << row_out.n_pixels << "\n";
            std::cout << "  Mean         : " << row_out.stats.mean << "\n";
            std::cout << "  Median       : " << row_out.stats.median << "\n";
            std::cout << "  Relief       : " << row_out.stats.range << "\n";
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
        const std::string csv_out = out_path(pp, "_basin_stats.csv");
        std::cout << "\nWriting basin statistics CSV: " << csv_out << "\n";
        write_basin_stats_csv(csv_out, basin_rows);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: failed to write basin statistics CSV: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Echo parameters for reproducibility
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-raster-basin-stats: done.\n";
    return 0;
}
