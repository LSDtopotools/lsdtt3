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

#include <ogr_spatialref.h>

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

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// Local record types
// -----------------------------------------------------------------------------
struct basin_outlet_record
{
    int    outlet_id = -1;
    double latitude  = 0.0;
    double longitude = 0.0;
};

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
// trim / string helpers
// -----------------------------------------------------------------------------
static std::string trim(const std::string& s)
{
    std::size_t begin = 0;
    while (begin < s.size() &&
           std::isspace(static_cast<unsigned char>(s[begin])))
    {
        ++begin;
    }

    std::size_t end = s.size();
    while (end > begin &&
           std::isspace(static_cast<unsigned char>(s[end - 1])))
    {
        --end;
    }

    return s.substr(begin, end - begin);
}

static std::string to_lower_copy(std::string s)
{
    for (char& ch : s)
    {
        ch = static_cast<char>(
            std::tolower(static_cast<unsigned char>(ch))
        );
    }
    return s;
}

static std::string strip_utf8_bom(const std::string& s)
{
    if (s.size() >= 3 &&
        static_cast<unsigned char>(s[0]) == 0xEF &&
        static_cast<unsigned char>(s[1]) == 0xBB &&
        static_cast<unsigned char>(s[2]) == 0xBF)
    {
        return s.substr(3);
    }

    return s;
}

// -----------------------------------------------------------------------------
// split_csv_line
//
// Simple CSV splitter for the expected outlet CSV.
// No quoted-field support is implemented here.
// -----------------------------------------------------------------------------
static std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string item;

    while (std::getline(ss, item, ','))
    {
        tokens.push_back(trim(item));
    }

    return tokens;
}

// -----------------------------------------------------------------------------
// load_basin_outlet_csv
//
// Expects header columns named latitude and longitude.
// -----------------------------------------------------------------------------
static std::vector<basin_outlet_record> load_basin_outlet_csv(
    const std::string& csv_path)
{
    std::ifstream infile(csv_path);
    if (!infile)
    {
        throw std::runtime_error(
            "Could not open basin_outlet_csv: " + csv_path
        );
    }

    std::string header_line;
    if (!std::getline(infile, header_line))
    {
        throw std::runtime_error(
            "Basin outlet CSV is empty: " + csv_path
        );
    }

    header_line = strip_utf8_bom(header_line);

    const std::vector<std::string> header_tokens = split_csv_line(header_line);

    int lat_col = -1;
    int lon_col = -1;

    for (std::size_t i = 0; i < header_tokens.size(); ++i)
    {
        const std::string key = to_lower_copy(header_tokens[i]);

        if (key == "latitude")
        {
            lat_col = static_cast<int>(i);
        }
        else if (key == "longitude")
        {
            lon_col = static_cast<int>(i);
        }
    }

    if (lat_col < 0 || lon_col < 0)
    {
        throw std::runtime_error(
            "Basin outlet CSV must contain columns named latitude and longitude"
        );
    }

    std::vector<basin_outlet_record> outlets;
    std::string line;
    int next_id = 1;

    while (std::getline(infile, line))
    {
        if (trim(line).empty())
        {
            continue;
        }

        const std::vector<std::string> tokens = split_csv_line(line);

        if (static_cast<int>(tokens.size()) <= std::max(lat_col, lon_col))
        {
            std::cerr << "Warning: skipping malformed CSV row: " << line << "\n";
            continue;
        }

        basin_outlet_record rec;
        rec.outlet_id = next_id++;

        try
        {
            rec.latitude  = std::stod(tokens[static_cast<std::size_t>(lat_col)]);
            rec.longitude = std::stod(tokens[static_cast<std::size_t>(lon_col)]);
        }
        catch (const std::exception&)
        {
            std::cerr << "Warning: skipping CSV row with invalid numeric "
                      << "latitude/longitude: " << line << "\n";
            continue;
        }

        outlets.push_back(rec);
    }

    if (outlets.empty())
    {
        throw std::runtime_error(
            "No valid basin outlet rows were read from: " + csv_path
        );
    }

    return outlets;
}

// -----------------------------------------------------------------------------
// lonlat_to_raster_xy
//
// Input longitude and latitude are WGS84 / EPSG:4326.
// Output x,y are in the raster CRS.
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_xy lonlat_to_raster_xy(
    double longitude,
    double latitude,
    const lsdtt::lsdtt_raster_meta& meta)
{
    if (meta.crs_wkt.empty())
    {
        throw std::runtime_error(
            "Raster CRS WKT is empty; cannot transform outlet coordinates"
        );
    }

    OGRSpatialReference src_srs;
    OGRSpatialReference dst_srs;

    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (src_srs.importFromEPSG(4326) != OGRERR_NONE)
    {
        throw std::runtime_error("Failed to initialise WGS84 spatial reference");
    }

    char* raster_wkt = const_cast<char*>(meta.crs_wkt.c_str());
    if (dst_srs.importFromWkt(&raster_wkt) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "Failed to parse raster CRS WKT for coordinate transformation"
        );
    }

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&src_srs, &dst_srs);

    if (ct == nullptr)
    {
        throw std::runtime_error(
            "Failed to create WGS84 -> raster CRS coordinate transformation"
        );
    }

    double x = longitude;
    double y = latitude;
    double z = 0.0;

    const int ok = ct->Transform(1, &x, &y, &z);
    OCTDestroyCoordinateTransformation(ct);

    if (!ok)
    {
        throw std::runtime_error(
            "Coordinate transformation from WGS84 to raster CRS failed"
        );
    }

    return lsdtt::lsdtt_xy{x, y};
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

    std::vector<basin_outlet_record> requested_outlets;
    try
    {
        requested_outlets = load_basin_outlet_csv(basin_csv_path);
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

    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const lsdtt::lsdtt_array2d_f& dem_data = dem.get_data_array();
    const lsdtt::lsdtt_array2d_u8& dem_mask = dem.get_mask_array();

    for (const auto& outlet : requested_outlets)
    {
        std::cout << "\nOutlet " << outlet.outlet_id
                  << "  lat=" << outlet.latitude
                  << " lon=" << outlet.longitude << "\n";

        try
        {
            // CSV is latitude,longitude.
            // Transform expects x=longitude, y=latitude.
            const lsdtt::lsdtt_xy xy =
                lonlat_to_raster_xy(outlet.longitude, outlet.latitude, meta);

            const int col = lsdtt::col_from_x(meta, xy.x);
            const int row = lsdtt::row_from_y(meta, xy.y);
            const std::int64_t flat = flow->flat_index(row, col);

            if (!flow->is_valid_cell(flat))
            {
                std::cerr << "Warning: outlet " << outlet.outlet_id
                          << " maps to an invalid or nodata pixel; skipping.\n";
                continue;
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
                std::cerr << "Warning: basin for outlet " << outlet.outlet_id
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
                    basin_pdf_cdf_path(pp, outlet.outlet_id);

                std::cout << "  Writing PDF/CDF: " << dist_csv << "\n";

                write_basin_pdf_cdf_csv(
                    dist_csv,
                    outlet.outlet_id,
                    hist
                );
            }



            basin_stats_record row_out;
            row_out.basin_id    = outlet.outlet_id;
            row_out.latitude    = outlet.latitude;
            row_out.longitude   = outlet.longitude;
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
                      << outlet.outlet_id << ": " << e.what() << "\n";
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
