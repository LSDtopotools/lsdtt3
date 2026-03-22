// lsdtt-test-flowinfo.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-test-flowinfo
//
// Builds D8 flow routing from a DEM and writes diagnostic rasters.
// Includes optional preprocessing (breach, fill) before routing.
//
// Usage:
//   lsdtt-test-flowinfo <path> <parameter_file>
//   lsdtt-test-flowinfo <parameter_file>          (path assumed to be ./)
//   lsdtt-test-flowinfo -h                        (print help and exit)
//   lsdtt-test-flowinfo -v                        (print version and exit)
//   lsdtt-test-flowinfo -c                        (print citation and exit)
//
// Parameter file format:
//   key: value
//   Lines beginning with # are ignored.
//   All keys use snake_case. No spaces in key names.
//
// Example parameter file:
//   read_path: /data/dems/
//   read_fname: my_dem.tif
//   write_path: /data/output/
//   write_fname: my_dem
//   boundary_conditions: oooo
//   breach_before_fill: true
//   fill_before_routing: true
//   write_n_contributing_pixels_raster: true
//   write_flow_direction_raster: true
//   write_donor_count_raster: true
//   write_outlet_raster: true
//
// Output filenames follow the convention:
//   <write_fname>_<descriptor>.tif
// e.g. for read_fname: my_dem.tif with the default write_fname:
//   my_dem_n_contributing_pixels.tif
//   my_dem_flow_direction.tif
//   my_dem_donor_count.tif
//   my_dem_outlet.tif
// -----------------------------------------------------------------------------

#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/utils/lsdtt_array2d.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt3 flowinfo test                     ||\n"
        << "|| D8 flow routing and diagnostic raster output        ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                    ||\n"
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
    auto is_flag = [](const std::string& s, std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags) if (s == f) return true;
        return false;
    };

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];

        if (is_flag(a, {"-h", "--help", "help", "-H"}))
        {
            std::cout
                << "Usage: lsdtt-test-flowinfo <path> <parameter_file>\n"
                << "       lsdtt-test-flowinfo <parameter_file>\n"
                << "       lsdtt-test-flowinfo -h   (this help)\n"
                << "       lsdtt-test-flowinfo -v   (version)\n"
                << "       lsdtt-test-flowinfo -c   (citation)\n\n"
                << "See docs/lsdtt_architecture.md for parameter file documentation.\n";
            exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-test-flowinfo version " << k_version << "\n";
            exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
        {
            std::cout << "Citation: " << k_citation << "\n";
            exit(0);
        }
    }

    std::string path  = "./";
    std::string fname;

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        exit(1);
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
// out_path — write_path + write_fname + suffix (suffix includes extension)
// -----------------------------------------------------------------------------
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& suffix)
{
    return pp.get_full_write_path() + suffix;
}

// -----------------------------------------------------------------------------
// resolve_output_format
// -----------------------------------------------------------------------------
static lsdtt::raster_output_format resolve_output_format(
    const std::string& fmt_str)
{
    if (fmt_str == "gtiff")         return lsdtt::raster_output_format::gtiff;
    if (fmt_str == "same_as_input") return lsdtt::raster_output_format::same_as_input;
    if (fmt_str != "cog")
    {
        std::cerr << "Warning: unknown write_output_format '"
                  << fmt_str << "' — defaulting to cog.\n";
    }
    return lsdtt::raster_output_format::cog;
}

// -----------------------------------------------------------------------------
// build_boundary_string
//
// The param parser stores boundary conditions as a vector<string>, one
// character per entry. lsdtt_flowinfo_options expects a 4-character string
// (north, east, south, west). Pad or trim to exactly 4 characters;
// missing entries default to 'o'.
// -----------------------------------------------------------------------------
static std::string build_boundary_string(
    const lsdtt::lsdtt_param_parser& pp)
{
    const auto& bc_vec = pp.get_boundary_conditions();

    std::string result = "oooo";
    for (std::size_t i = 0; i < 4 && i < bc_vec.size(); ++i)
    {
        if (!bc_vec[i].empty()) result[i] = bc_vec[i][0];
    }
    return result;
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = parse_args(argc, argv);

    // -------------------------------------------------------------------------
    // Parameter file
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
    // Default parameter maps
    // -------------------------------------------------------------------------
    const std::map<std::string, bool> bool_defaults =
    {
        // Preprocessing
        {"breach_before_fill",                  false},
        {"fill_remainder_after_breach",         true},
        {"fill_before_routing",                 true},

        // Diagnostic raster outputs
        {"write_n_contributing_pixels_raster",  false},
        {"write_flow_direction_raster",         false},
        {"write_donor_count_raster",            false},
        {"write_outlet_raster",                 false},

        // Pass-through preprocessing outputs (matches lsdtt-raster-preprocessing
        // conventions so parameter files can be shared)
        {"write_fill_raster",                   false},
        {"write_breached_raster",               false},
    };

    const std::map<std::string, float> float_defaults =
    {
        {"fill_min_slope",                      0.0001f},
        {"minimum_elevation",                   0.0f},
        {"maximum_elevation",                   30000.0f},
        {"breach_lc_max_cost",                  std::numeric_limits<float>::max()},
        {"breach_lc_flat_increment",            0.0f},
    };

    const std::map<std::string, int> int_defaults =
    {
        {"n_threads",                           0},
        {"breach_lc_max_dist",                  1000},
    };

    const std::map<std::string, std::string> string_defaults =
    {
        {"fill_algorithm_choice",               "barnes2014"},
        {"breach_algorithm_choice",             "lindsay2016"},
        {"write_output_format",                 "cog"},
    };

    // Help map — {type, default, description, guidance}
    const std::map<std::string, std::vector<std::string>> help_map =
    {
        {"breach_before_fill",
            {"bool","false","Breach depressions before pit filling.",
             "Recommended for lidar DEMs. Uses breach_algorithm_choice."}},
        {"fill_remainder_after_breach",
            {"bool","true","Fill remaining pits after breaching.",
             "Used when breach_before_fill is true."}},
        {"fill_before_routing",
            {"bool","true","Fill pits before building the flow graph.",
             "Highly recommended. The DEM must be pit-free for D8 routing "
             "to produce a valid acyclic graph."}},
        {"fill_algorithm_choice",
            {"string","barnes2014","Algorithm for pit filling.",
             "barnes2014: Improved Priority-Flood (recommended). "
             "wang_liu: original Wang and Liu (2006), routed to barnes2014."}},
        {"breach_algorithm_choice",
            {"string","lindsay2016","Algorithm for depression breaching.",
             "lindsay2016: Lindsay (2016) efficient breaching. "
             "least_cost: Lindsay (2020) pit-centric least-cost breaching."}},
        {"fill_min_slope",
            {"float","0.0001","Minimum slope between filled cells (m/m).",
             "Prevents flat regions in filled DEMs. Set to 0 for exact flats."}},
        {"minimum_elevation",
            {"float","0.0","Pixels below this elevation are set to nodata.",
             "Applied before preprocessing. Useful for removing ocean pixels."}},
        {"maximum_elevation",
            {"float","30000","Pixels above this elevation are set to nodata.",
             "Applied before preprocessing. Useful for removing artefacts."}},
        {"write_n_contributing_pixels_raster",
            {"bool","false","Write a raster of upstream contributing pixel count.",
             "Each pixel contains the number of pixels draining through it "
             "(including itself). Proportional to drainage area."}},
        {"write_flow_direction_raster",
            {"bool","false","Write the D8 flow direction raster.",
             "Direction codes: 0=N 1=NE 2=E 3=SE 4=S 5=SW 6=W 7=NW. "
             "Nodata cells are coded -1. Outlets are coded -2."}},
        {"write_donor_count_raster",
            {"bool","false","Write a raster of D8 donor (upstream neighbour) count.",
             "Values 0-8. Headwaters have donor count 0."}},
        {"write_outlet_raster",
            {"bool","false","Write a binary raster marking outlet cells.",
             "Outlet cells (value 1) drain off-grid or to a boundary. "
             "All other valid cells are 0. Nodata cells are nodata."}},
        {"write_fill_raster",
            {"bool","false","Write the pit-filled DEM.",
             "Written before flow routing if fill_before_routing is true."}},
        {"write_breached_raster",
            {"bool","false","Write the depression-breached DEM.",
             "Written before pit filling if breach_before_fill is true."}},
        {"write_output_format",
            {"string","cog","Output raster format.",
             "cog: Cloud Optimised GeoTIFF (default). "
             "gtiff: plain compressed GeoTIFF. "
             "same_as_input: match the input raster format."}},
        {"n_threads",
            {"int","0","Number of OpenMP threads. 0 = all available.",""}},
        {"breach_lc_max_dist",
            {"int","1000","Maximum search radius in pixels for least-cost breaching.",
             "Pits beyond this radius are left unbreached."}},
        {"breach_lc_max_cost",
            {"float","unlimited","Maximum total elevation lowering for a breach path.",
             "Paths more costly than this are abandoned."}},
    };

    pp.parse_all_parameters(float_defaults, int_defaults, bool_defaults, string_defaults);

    const auto bf  = pp.get_bool_parameters();
    const auto ff  = pp.get_float_parameters();
    const auto ifs = pp.get_int_parameters();
    const auto sf  = pp.get_string_parameters();

    const int n_threads = ifs.at("n_threads");
    const lsdtt::raster_output_format out_fmt =
        resolve_output_format(sf.at("write_output_format"));
    const std::string bc_string = build_boundary_string(pp);

    std::cout << "Boundary conditions: " << bc_string << "\n";
    std::cout << "Output format: " << sf.at("write_output_format") << "\n";

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    try
    {
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    const lsdtt::lsdtt_raster_meta meta = dem.get_meta();

    // -------------------------------------------------------------------------
    // Elevation threshold masking
    // -------------------------------------------------------------------------
    const float min_elev = ff.at("minimum_elevation");
    const float max_elev = ff.at("maximum_elevation");

    if (min_elev > 0.0f)
    {
        std::cout << "Masking pixels below " << min_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, min_elev, true);
    }
    if (max_elev < 30000.0f)
    {
        std::cout << "Masking pixels above " << max_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, max_elev, false);
    }

    // -------------------------------------------------------------------------
    // Preprocessing — breach then fill
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu working_dem = dem;

    if (bf.at("breach_before_fill") || bf.at("write_breached_raster"))
    {
        std::cout << "\n--- Depression breaching ---\n";
        const std::string algo = sf.at("breach_algorithm_choice");

        try
        {
            lsdtt::lsdtt_raster_cpu breached;

            if (algo == "least_cost")
            {
                std::cout << "Algorithm: Lindsay (2020) least-cost breaching\n";
                breached = lsdtt::breach_depressions_least_cost(
                    working_dem,
                    ifs.at("breach_lc_max_dist"),
                    ff.at("breach_lc_max_cost"),
                    true,
                    ff.at("breach_lc_flat_increment"),
                    bf.at("fill_remainder_after_breach"),
                    n_threads
                );
            }
            else
            {
                std::cout << "Algorithm: Lindsay (2016) efficient breaching\n";
                breached = lsdtt::breach_depressions(working_dem, n_threads);
            }

            if (bf.at("write_breached_raster"))
            {
                const std::string out = out_path(pp, "_breached.tif");
                std::cout << "Writing: " << out << "\n";
                loader.write_full(breached.get_data(), breached.get_meta(),
                                  out, out_fmt);
            }

            if (bf.at("breach_before_fill"))
            {
                working_dem = std::move(breached);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: breaching failed: " << e.what() << "\n";
        }
    }

    if (bf.at("fill_before_routing") || bf.at("write_fill_raster"))
    {
        std::cout << "\n--- Pit filling ---\n";

        try
        {
            lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(
                working_dem,
                ff.at("fill_min_slope"),
                n_threads
            );

            if (bf.at("write_fill_raster"))
            {
                const std::string out = out_path(pp, "_fill.tif");
                std::cout << "Writing: " << out << "\n";
                loader.write_full(filled.get_data(), filled.get_meta(),
                                  out, out_fmt);
            }

            if (bf.at("fill_before_routing"))
            {
                working_dem = std::move(filled);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: pit filling failed: " << e.what() << "\n";
            std::cerr << "Proceeding with unfilled DEM. Routing may fail "
                         "if pits remain.\n";
        }
    }

    // -------------------------------------------------------------------------
    // Build flowinfo
    // -------------------------------------------------------------------------
    std::cout << "\n--- Building D8 flow routing ---\n";
    std::cout << "Boundary conditions: " << bc_string << "\n";

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions      = bc_string;
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.require_acyclic_graph    = true;
    fi_opts.build_topological_order  = true;
    fi_opts.build_s_vector           = true;
    fi_opts.n_threads                = n_threads;

    const lsdtt::lsdtt_raster_meta& fi_meta  = working_dem.get_meta();
    const lsdtt::lsdtt_array2d_f&   dem_data = working_dem.get_data_array();
    const lsdtt::lsdtt_array2d_u8&  dem_mask = working_dem.get_mask_array();

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try
    {
        flow = lsdtt::make_flow_base(dem_data, dem_mask, fi_meta, fi_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Fatal: flow routing failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Flow routing complete.\n";
    std::cout << "  Valid cells  : " << flow->n_valid_cells()  << "\n";
    std::cout << "  Outlets      : " << flow->n_outlets()      << "\n";
    std::cout << "  Acyclic graph: " << (flow->is_acyclic() ? "yes" : "no") << "\n";

    // -------------------------------------------------------------------------
    // Diagnostic raster output helpers
    // -------------------------------------------------------------------------
    const int    nrows   = fi_meta.n_rows;
    const int    ncols   = fi_meta.n_cols;
    const float  nodata  = static_cast<float>(fi_meta.no_data_value);

    // Convenience: write a float raster assembled from a per-cell lambda.
    // The lambda receives a flat_index_type and returns float. Nodata and
    // invalid cells are written as nodata.
    auto write_diagnostic = [&](
        const std::string&                                          suffix,
        std::function<float(lsdtt::lsdtt_flow_base::flat_index_type)> value_fn)
    {
        const std::string out = out_path(pp, suffix);
        std::cout << "Writing: " << out << "\n";

        lsdtt::lsdtt_array2d_f raster(nrows, ncols, nodata);

        for (int r = 0; r < nrows; ++r)
        {
            for (int c = 0; c < ncols; ++c)
            {
                const auto flat = flow->flat_index_unchecked(r, c);
                if (!flow->is_valid_cell(flat)) continue;
                raster(r, c) = value_fn(flat);
            }
        }

        loader.write_full(raster.get_vector(), fi_meta, out, out_fmt);
    };

    // -------------------------------------------------------------------------
    // n_contributing_pixels
    //
    // The number of pixels draining through each cell (including itself) is
    // the size of that cell's upslope interval in the s_vector. This is
    // proportional to drainage area.
    // -------------------------------------------------------------------------
    if (bf.at("write_n_contributing_pixels_raster"))
    {
        std::cout << "\n--- n_contributing_pixels ---\n";
        write_diagnostic(
            "_n_contributing_pixels.tif",
            [&](lsdtt::lsdtt_flow_base::flat_index_type flat) -> float
            {
                return static_cast<float>(flow->s_vector_size_unchecked(flat));
            }
        );
    }

    // -------------------------------------------------------------------------
    // Flow direction
    //
    // Direction codes (receiver_dir_):
    //   0=N 1=NE 2=E 3=SE 4=S 5=SW 6=W 7=NW
    //   no_direction (255) for outlet cells → written as -2 to distinguish
    //   from nodata (-1 / no_data_value).
    // -------------------------------------------------------------------------
    if (bf.at("write_flow_direction_raster"))
    {
        std::cout << "\n--- Flow direction ---\n";
        write_diagnostic(
            "_flow_direction.tif",
            [&](lsdtt::lsdtt_flow_base::flat_index_type flat) -> float
            {
                const std::uint8_t dir = flow->receiver_dir_of_unchecked(flat);
                if (dir == lsdtt::lsdtt_flow_base::no_direction)
                {
                    return -2.0f;
                }
                return static_cast<float>(dir);
            }
        );
    }

    // -------------------------------------------------------------------------
    // Donor count
    // -------------------------------------------------------------------------
    if (bf.at("write_donor_count_raster"))
    {
        std::cout << "\n--- Donor count ---\n";
        write_diagnostic(
            "_donor_count.tif",
            [&](lsdtt::lsdtt_flow_base::flat_index_type flat) -> float
            {
                return static_cast<float>(flow->donor_count_of_unchecked(flat));
            }
        );
    }

    // -------------------------------------------------------------------------
    // Outlet mask
    // -------------------------------------------------------------------------
    if (bf.at("write_outlet_raster"))
    {
        std::cout << "\n--- Outlet raster ---\n";
        write_diagnostic(
            "_outlet.tif",
            [&](lsdtt::lsdtt_flow_base::flat_index_type flat) -> float
            {
                return flow->is_outlet_unchecked(flat) ? 1.0f : 0.0f;
            }
        );
    }

    // -------------------------------------------------------------------------
    // Help and parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map,
                  pp.get_write_fname() + "_README",
                  k_version,
                  k_citation);

    pp.print_parameters();

    std::cout << "\nlsdtt-test-flowinfo: done.\n";
    return 0;
}
