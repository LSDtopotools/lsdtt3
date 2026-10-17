// lsdtt-lithology-erodibility.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-lithology-erodibility
//
// Turns a lithology map (polygons with a lithology field, e.g. a Macrostrat
// download, or a categorical code raster) and a lithology -> K property table
// into a stream-power erodibility raster for lsdtt-lem (erodibility_fname),
// on the grid the LEM run will use: a DEM's, a synthetic grid's (with the same
// placements as lsdtt-lem: clips, relabelling, warping), or the map's own.
//
// Outputs: <write_fname>_erodibility.tif (K in the table's convention),
// <write_fname>_lithology_code.tif, <write_fname>_lithology_legend.csv and
// <write_fname>_lem_erodibility_parameters.txt (the lsdtt-lem lines stating
// the table's calibration).
//
// Usage:
//   lsdtt-lithology-erodibility <path> <parameter_file>
//   lsdtt-lithology-erodibility <parameter_file>   (path assumed to be ./)
//   lsdtt-lithology-erodibility -h                 (print help and exit)
//   lsdtt-lithology-erodibility -v                 (print version and exit)
//   lsdtt-lithology-erodibility -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_lem_input_grid_cpu.hpp"
#include "../src/cpu/lsdtt_lithology_map_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_stream_power_cpu.hpp"
#include "../src/cpu/lsdtt_vector_rasterise_cpu.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-lithology-erodibility";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-lithology-erodibility              ||\n"
        << "|| I turn a lithology map into an erodibility raster   ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
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

namespace
{

bool looks_like_raster(const std::string& path)
{
    std::string ext = std::filesystem::path(path).extension().string();
    for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    for (const char* e : {".tif", ".tiff", ".asc", ".bil", ".img", ".vrt", ".flt", ".nc"})
        if (ext == e) return true;
    return false;
}

// A synthetic, all-valid raster on `meta`, used as the rasterisation target.
lsdtt::lsdtt_raster_cpu blank_on(const lsdtt::lsdtt_raster_meta& meta)
{
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    return lsdtt::lsdtt_raster_cpu(std::vector<float>(n, 0.0f), std::vector<std::uint8_t>(n, 1),
                                   meta);
}

// Refuses a source CRS that cannot carry a metric clip point or grid.
void check_metric_crs(const lsdtt::lsdtt_raster_meta& m, const std::string& label)
{
    if (m.crs_wkt.empty())
        throw std::invalid_argument(label + " has no CRS");
    if (m.is_geographic())
        throw std::invalid_argument(label + " is in a geographic (degree) CRS; reproject it to a "
                                    "projected CRS in metres first (e.g. ogr2ogr -t_srs EPSG:<UTM>)");
    if (!m.is_in_metres())
        throw std::invalid_argument(label + " CRS unit is '" + m.linear_unit_name() +
                                    "', not metres");
    lsdtt::check_lem_projection_distortion(m, label);
}

} // namespace

int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    std::map<std::string, double>      double_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Lithology input ---
    string_default_map["lithology_fname"] = "";
    help_map["lithology_fname"] = {"string","",
        "Lithology map (in read_path): a polygon vector file or a categorical code raster.",
        "Vectors: any OGR polygon layer (GeoPackage, FlatGeobuf, GeoJSON, Shapefile), e.g. a "
        "Macrostrat download. Rasters: integer lithology codes. Must be in a projected CRS in "
        "metres for clips and for taking the grid from it.",
        "Lithology input", "param", ""};

    string_default_map["lithology_source"] = "auto";
    help_map["lithology_source"] = {"string","auto",
        "auto, vector_map or code_raster.",
        "auto decides from the file extension (.tif, .tiff, .asc, .bil, .img, .vrt, .flt, .nc "
        "are rasters).",
        "Lithology input", "param", ""};

    string_default_map["lithology_layer_name"] = "";
    help_map["lithology_layer_name"] = {"string","",
        "Vector layer to read. Empty (default) = the first layer.", "",
        "Lithology input", "param", "lithology_fname"};

    string_default_map["lithology_field_name"] = "lith_primary";
    help_map["lithology_field_name"] = {"string","lith_primary",
        "Attribute of the vector layer that holds the lithology.",
        "Values are matched to the property table by integer code, by table name "
        "(case-insensitive), or through lithology_crosswalk_fname. lith_primary is the first "
        "lithology of each unit in lsdtt-fetch-lithology output.",
        "Lithology input", "param", "lithology_fname"};

    string_default_map["lithology_crosswalk_fname"] = "";
    help_map["lithology_crosswalk_fname"] = {"string","",
        "Optional CSV (in read_path) with columns value,code mapping map values to table codes.",
        "When given, every map value present must be in it. Use it to map many source names "
        "(e.g. Macrostrat lithology terms) onto a few table units.",
        "Lithology input", "param", "lithology_fname"};

    string_default_map["lithology_property_table_fname"] = "";
    help_map["lithology_property_table_fname"] = {"string","",
        "Lithology property table (in read_path): code, name, erodibility per unit.",
        "CSV with '#' header lines stating erodibility_convention (map_area or specific_area), "
        "m_exponent and n_exponent (all required) and optionally "
        "calibration_cell_size_metres and provenance, then a header row code,name,"
        "erodibility[,reference]. K must be > 0. See the CLI reference for the format.",
        "Lithology input", "param", ""};

    int_default_map["lithology_unmapped_code"] = -1;
    help_map["lithology_unmapped_code"] = {"int","-1",
        "Table code given to model-domain pixels with no lithology. -1 (default) = refuse.",
        "Pixels a map does not cover (holes, the edge of a map) inside the domain are refused "
        "unless a table unit is named here. With model_grid_source: lithology the map's "
        "footprint is the domain, so uncovered pixels are nodata instead.",
        "Lithology input", "param", ""};

    // --- Grid ---
    string_default_map["model_grid_source"] = "auto";
    help_map["model_grid_source"] = {"string","auto",
        "What defines the output grid: auto, dem, parameters or lithology.",
        "dem: the grid of dem_fname (use the same DEM in lsdtt-lem). parameters: n_rows, n_cols, "
        "cell_size_metres on the synthetic georeference (the same grid lsdtt-lem builds from "
        "the same keys). lithology: the code raster's own grid, or a vector's extent at "
        "cell_size_metres (required), with the map's footprint as the domain; run lsdtt-lem "
        "with model_grid_source: erodibility_raster on the result. auto = dem when dem_fname is "
        "given, else lithology.",
        "Grid", "param", ""};

    string_default_map["dem_fname"] = "";
    help_map["dem_fname"] = {"string","",
        "DEM (in read_path) whose grid the K raster is made on (model_grid_source: dem).", "",
        "Grid", "param", ""};

    int_default_map["n_rows"] = 200;
    help_map["n_rows"] = {"int","200", "Rows of the synthetic grid (model_grid_source: parameters).",
        "", "Grid", "param", "model_grid_source"};
    int_default_map["n_cols"] = 400;
    help_map["n_cols"] = {"int","400", "Columns of the synthetic grid.", "", "Grid", "param",
        "model_grid_source"};
    float_default_map["cell_size_metres"] = 30.0f;
    help_map["cell_size_metres"] = {"float","30",
        "Cell size of the synthetic grid, or of the grid taken from a vector map.", "",
        "Grid", "param", "model_grid_source"};
    string_default_map["synthetic_crs_wkt"] = "";
    help_map["synthetic_crs_wkt"] = {"string","",
        "CRS of the synthetic grid. Empty = Watopia (EPSG:32758), as in lsdtt-lem.", "",
        "Grid", "param", "model_grid_source"};

    string_default_map["input_raster_placement"] = "as_georeferenced";
    help_map["input_raster_placement"] = {"string","as_georeferenced",
        "How the lithology map is placed on the grid (the lsdtt-lem placements).",
        "as_georeferenced: in place (vectors are rasterised in the grid's CRS; rasters must "
        "already be on the grid). warp_to_model_grid: rasters warped (near or mode). "
        "relabel_whole_raster: a code raster with the grid's rows, cols and cell size. "
        "relabel_centre_clip / relabel_lower_left_clip: an n_rows x n_cols window at "
        "cell_size_metres around / from the clip point in the map's CRS, relabelled onto the "
        "synthetic grid. relabel_* need model_grid_source: parameters.",
        "Grid", "param", ""};

    double_default_map["input_raster_clip_easting_metres"] =
        std::numeric_limits<double>::quiet_NaN();
    help_map["input_raster_clip_easting_metres"] = {"double","nan",
        "Easting of the clip centre or lower-left corner, in the map's CRS.", "",
        "Grid", "param", "input_raster_placement"};
    double_default_map["input_raster_clip_northing_metres"] =
        std::numeric_limits<double>::quiet_NaN();
    help_map["input_raster_clip_northing_metres"] = {"double","nan",
        "Northing of the clip centre or lower-left corner, in the map's CRS.", "",
        "Grid", "param", "input_raster_placement"};
    string_default_map["input_raster_resampling_method"] = "near";
    help_map["input_raster_resampling_method"] = {"string","near",
        "Resampling for a warped code raster: near or mode (clips always use near).", "",
        "Grid", "param", "input_raster_placement"};

    // --- Output ---
    bool_default_map["write_lithology_code_raster"] = true;
    help_map["write_lithology_code_raster"] = {"bool","true",
        "Write the lithology codes the K raster was made from.",
        "Output: <write_fname>_lithology_code.tif.", "Output", "switch", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Raster output format: cog or gtiff.", "", "Output", "io", ""};

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

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
    pp.parse_all_parameters(float_default_map, int_default_map, bool_default_map,
                            string_default_map, double_default_map);
    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();
    auto df  = pp.get_double_parameters();

    if (sf["write_output_format"] != "cog" && sf["write_output_format"] != "gtiff")
    {
        std::cerr << "FATAL: write_output_format must be cog or gtiff.\n";
        return 1;
    }
    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    try
    {
        // ---------------------------------------------------------------------
        // Inputs
        // ---------------------------------------------------------------------
        if (sf["lithology_fname"].empty())
            throw std::invalid_argument("lithology_fname is required");
        if (sf["lithology_property_table_fname"].empty())
            throw std::invalid_argument("lithology_property_table_fname is required");
        const std::string lith_path = pp.add_read_path_if_needed(sf["lithology_fname"]);
        if (!std::filesystem::exists(lith_path))
            throw std::invalid_argument("lithology_fname does not exist: " + lith_path);

        std::string source = sf["lithology_source"];
        if (source == "auto") source = looks_like_raster(lith_path) ? "code_raster" : "vector_map";
        if (source != "vector_map" && source != "code_raster")
            throw std::invalid_argument("lithology_source must be auto, vector_map or code_raster "
                                        "(got '" + sf["lithology_source"] + "')");
        const bool is_vector = source == "vector_map";

        const std::string table_path =
            pp.add_read_path_if_needed(sf["lithology_property_table_fname"]);
        const lsdtt::lsdtt_lithology_property_table table =
            lsdtt::load_lithology_property_table(table_path);
        std::map<std::string, int> crosswalk;
        const bool have_crosswalk = !sf["lithology_crosswalk_fname"].empty();
        if (have_crosswalk)
            crosswalk = lsdtt::load_lithology_crosswalk(
                pp.add_read_path_if_needed(sf["lithology_crosswalk_fname"]));

        const bool table_specific =
            table.calibration.area_convention == lsdtt::erodibility_area_convention::specific_area;
        std::cout << "\nLithology map   : " << lith_path << " ("
                  << (is_vector ? "vector, field " + sf["lithology_field_name"] : std::string("code raster"))
                  << ")\nProperty table  : " << table_path << " (" << table.units.size()
                  << " units)\n  convention    : "
                  << (table_specific ? "specific_area" : "map_area") << ", m = "
                  << table.calibration.m_exponent << ", n = " << table.calibration.n_exponent
                  << "\n  K units       : "
                  << lsdtt::describe_erodibility_units(table.calibration.m_exponent, table_specific)
                  << "\n";
        if (!table.provenance.empty()) std::cout << "  provenance    : " << table.provenance << "\n";
        for (const auto& n : table.notes) std::cout << "  " << n << "\n";
        for (const auto& w : table.warnings) std::cerr << "WARNING: " << w << "\n";

        // ---------------------------------------------------------------------
        // Grid and placement
        // ---------------------------------------------------------------------
        lsdtt::lsdtt_lem_input_raster_placement_options place;
        place.placement = lsdtt::parse_lem_input_raster_placement(sf["input_raster_placement"]);
        place.clip_easting = df["input_raster_clip_easting_metres"];
        place.clip_northing = df["input_raster_clip_northing_metres"];
        place.resampling_method = sf["input_raster_resampling_method"];
        place.categorical = true;
        using P = lsdtt::lem_input_raster_placement;
        const bool clip = place.placement == P::relabel_centre_clip ||
                          place.placement == P::relabel_lower_left_clip;
        const bool relabel = clip || place.placement == P::relabel_whole_raster;
        if (clip && !(std::isfinite(place.clip_easting) && std::isfinite(place.clip_northing)))
            throw std::invalid_argument("input_raster_placement: " + sf["input_raster_placement"] +
                                        " needs input_raster_clip_easting_metres and "
                                        "input_raster_clip_northing_metres");

        std::string grid_source = sf["model_grid_source"];
        if (grid_source == "auto") grid_source = sf["dem_fname"].empty() ? "lithology" : "dem";
        if (grid_source != "dem" && grid_source != "parameters" && grid_source != "lithology")
            throw std::invalid_argument("model_grid_source must be auto, dem, parameters or "
                                        "lithology (got '" + sf["model_grid_source"] + "')");
        if (relabel && grid_source != "parameters")
            throw std::invalid_argument("input_raster_placement: " + sf["input_raster_placement"] +
                                        " relabels the map onto a synthetic grid, so it needs "
                                        "model_grid_source: parameters");
        place.model_grid_is_synthetic = grid_source == "parameters";

        lsdtt::lsdtt_raster_cpu grid;          // mask = the model domain
        bool footprint_is_domain = false;
        lsdtt::polygon_layer_info layer;
        if (is_vector) layer = lsdtt::read_polygon_layer_info(lith_path, sf["lithology_layer_name"]);

        if (grid_source == "dem")
        {
            if (sf["dem_fname"].empty())
                throw std::invalid_argument("model_grid_source: dem needs dem_fname");
            const std::string dem_path = pp.add_read_path_if_needed(sf["dem_fname"]);
            lsdtt::lsdtt_data_loader l;
            l.open_raster(dem_path);
            l.get_meta().check_projected_metres();
            grid = lsdtt::lsdtt_raster_cpu(l);
            std::cout << "Grid            : the DEM " << dem_path << " ("
                      << grid.get_meta().n_rows << " x " << grid.get_meta().n_cols << " at "
                      << grid.get_meta().cell_size << " m)\n";
        }
        else if (grid_source == "parameters")
        {
            std::vector<std::string> warnings;
            grid = blank_on(lsdtt::make_lem_synthetic_grid(
                ifs["n_rows"], ifs["n_cols"], static_cast<double>(ff["cell_size_metres"]),
                sf["synthetic_crs_wkt"], &warnings));
            for (const auto& w : warnings) std::cerr << "WARNING: " << w << "\n";
            std::cout << "Grid            : synthetic " << ifs["n_rows"] << " x " << ifs["n_cols"]
                      << " at " << ff["cell_size_metres"] << " m ("
                      << grid.get_meta().crs_name() << ")\n";
        }
        else   // lithology: the map's own grid; its footprint is the domain
        {
            footprint_is_domain = true;
            if (is_vector)
            {
                if (!pp.parameter_supplied("cell_size_metres"))
                    throw std::invalid_argument("a vector map has no cell size: with "
                                                "model_grid_source: lithology give cell_size_metres");
                const double dx = ff["cell_size_metres"];
                lsdtt::lsdtt_raster_meta m;
                m.cell_size = dx;
                m.x_min = layer.x_min;
                m.y_min = layer.y_min;
                m.n_cols = std::max(1, static_cast<int>(std::ceil((layer.x_max - layer.x_min) / dx)));
                m.n_rows = std::max(1, static_cast<int>(std::ceil((layer.y_max - layer.y_min) / dx)));
                m.no_data_value = -9999.0;
                m.crs_wkt = layer.crs_wkt;
                check_metric_crs(m, "lithology vector map");
                grid = blank_on(m);
            }
            else
            {
                lsdtt::lsdtt_data_loader l;
                l.open_raster(lith_path);
                grid = blank_on(lsdtt::make_lem_model_grid_from_raster(
                    lsdtt::describe_lem_source_grid(l, "lithology code raster"),
                    "lithology code raster"));
            }
            std::cout << "Grid            : the lithology map's own ("
                      << grid.get_meta().n_rows << " x " << grid.get_meta().n_cols << " at "
                      << grid.get_meta().cell_size << " m, " << grid.get_meta().crs_name()
                      << "); its footprint is the model domain\n";
        }
        std::cout << "Placement       : " << sf["input_raster_placement"] << "\n";

        // ---------------------------------------------------------------------
        // Lithology values on the grid
        // ---------------------------------------------------------------------
        std::vector<std::int32_t> value_index;
        std::vector<std::string> values;
        const lsdtt::lsdtt_raster_meta& gm = grid.get_meta();
        const std::size_t npix = static_cast<std::size_t>(gm.n_rows) * gm.n_cols;

        if (is_vector)
        {
            if (place.placement == P::relabel_whole_raster)
                throw std::invalid_argument("relabel_whole_raster needs a raster: a vector map "
                                            "has no rows, cols or cell size (use a clip "
                                            "placement)");
            if (place.placement == P::warp_to_model_grid)
                std::cout << "  (a vector is reprojected exactly; warp_to_model_grid and "
                             "as_georeferenced are the same for it)\n";
            lsdtt::polygon_field_rasterise_result r;
            if (clip)
            {
                lsdtt::lsdtt_raster_meta src;
                src.crs_wkt = layer.crs_wkt;
                src.n_rows = gm.n_rows;
                src.n_cols = gm.n_cols;
                src.cell_size = gm.cell_size;
                const lsdtt::lsdtt_lem_clip_extent w = lsdtt::compute_lem_clip_extent(
                    place.placement, place.clip_easting, place.clip_northing, gm.n_rows,
                    gm.n_cols, gm.cell_size);
                src.x_min = w.x_min;
                src.y_min = w.y_min;
                src.no_data_value = -9999.0;
                check_metric_crs(src, "lithology vector map");
                lsdtt::check_lem_clip_extent_within(
                    w, {layer.x_min, layer.y_min, layer.x_max, layer.y_max}, "lithology vector map");
                r = lsdtt::rasterise_polygon_field_values(blank_on(src), lith_path,
                                                          sf["lithology_layer_name"],
                                                          sf["lithology_field_name"]);
                std::cout << "  clip window [" << std::setprecision(12) << w.x_min << ", "
                          << w.x_max << "] x [" << w.y_min << ", " << w.y_max
                          << "] in the map's CRS, relabelled onto the synthetic grid\n"
                          << std::setprecision(6);
            }
            else
                r = lsdtt::rasterise_polygon_field_values(grid, lith_path,
                                                          sf["lithology_layer_name"],
                                                          sf["lithology_field_name"]);
            value_index = std::move(r.value_index);
            values = std::move(r.values);
            std::cout << "  rasterised " << r.n_features_read << " features ("
                      << r.n_features_skipped << " skipped), " << values.size()
                      << " distinct values; " << r.n_gap_filled
                      << " sliver pixels filled by coverage; " << r.n_overlap_pixels
                      << " overlapping pixels\n";
        }
        else
        {
            lsdtt::lsdtt_lem_prepared_input_raster prepared;
            if (footprint_is_domain)
            {
                lsdtt::lsdtt_lem_input_raster_placement_options asis = place;
                asis.placement = P::as_georeferenced;
                prepared = lsdtt::prepare_lem_input_raster(lith_path, gm, asis,
                                                           "lithology code raster");
            }
            else
                prepared = lsdtt::prepare_lem_input_raster(lith_path, gm, place,
                                                           "lithology code raster");
            for (const auto& n : prepared.notes) std::cout << "  " << n << "\n";
            for (const auto& w : prepared.warnings) std::cerr << "WARNING: " << w << "\n";
            const auto& d = prepared.raster.get_data();
            const auto& m = prepared.raster.get_mask();
            std::map<long, std::int32_t> seen;
            std::vector<long> code_of_pixel(npix, -1);
            for (std::size_t i = 0; i < npix; ++i)
            {
                if (!m[i] || !grid.get_mask()[i]) continue;
                const float v = d[i];
                if (!(v >= 0.0f) || v != std::floor(v) || v > 65534.0f)
                {
                    std::ostringstream msg;
                    msg << "lithology code raster: value " << v << " at pixel (" << i / gm.n_cols
                        << ", " << i % gm.n_cols << ") is not an integer code in 0..65534";
                    throw std::invalid_argument(msg.str());
                }
                code_of_pixel[i] = static_cast<long>(v);
                seen[code_of_pixel[i]] = 0;
            }
            for (auto& [code, idx] : seen)
            {
                idx = static_cast<std::int32_t>(values.size());
                values.push_back(std::to_string(code));
            }
            value_index.assign(npix, -1);
            for (std::size_t i = 0; i < npix; ++i)
                if (code_of_pixel[i] >= 0) value_index[i] = seen[code_of_pixel[i]];
        }

        // Coverage on a synthetic grid placed in place is usually poor.
        {
            std::size_t mapped = 0, domain = 0;
            for (std::size_t i = 0; i < npix; ++i)
            {
                if (!grid.get_mask()[i]) continue;
                ++domain;
                if (value_index[i] >= 0) ++mapped;
            }
            if (mapped == 0)
                throw std::invalid_argument("the lithology map covers no pixel of the grid; check "
                                            "its CRS and extent, or use a clip placement");
            if (!footprint_is_domain && 2 * mapped < domain)
                std::cerr << "WARNING: the lithology map covers only "
                          << (100.0 * static_cast<double>(mapped) / domain)
                          << "% of the model domain\n";
        }

        std::vector<std::size_t> counts(values.size(), 0);
        for (std::int32_t v : value_index)
            if (v >= 0) ++counts[static_cast<std::size_t>(v)];
        const std::vector<int> codes = lsdtt::lithology_codes_for_values(
            values, counts, have_crosswalk ? &crosswalk : nullptr, table);

        lsdtt::lsdtt_lithology_erodibility_options eo;
        eo.unmapped_code = ifs["lithology_unmapped_code"];
        eo.unmapped_inside_domain = !footprint_is_domain;
        const lsdtt::lsdtt_lithology_erodibility_result res =
            lsdtt::make_lithology_erodibility(grid, value_index, values, codes, table, eo);
        for (const auto& n : res.notes) std::cout << "  " << n << "\n";
        for (const auto& w : res.warnings) std::cerr << "WARNING: " << w << "\n";

        // ---------------------------------------------------------------------
        // Outputs
        // ---------------------------------------------------------------------
        std::cout << "\nLithology units on the grid:\n";
        for (const auto& row : res.legend)
            std::cout << "  " << std::setw(6) << row.code << "  " << std::left << std::setw(24)
                      << row.name << std::right << "  K = " << row.erodibility << "  ("
                      << row.n_pixels << " pixels)\n";

        lsdtt::lsdtt_data_loader writer;
        const std::string k_path = pp.make_output_path("_erodibility", ".tif");
        writer.write_full(res.erodibility_raster.get_data(), res.erodibility_raster.get_meta(),
                          k_path, out_fmt);
        std::cout << "\nWrote " << k_path << "\n";
        if (bf["write_lithology_code_raster"])
        {
            const std::string c_path = pp.make_output_path("_lithology_code", ".tif");
            writer.write_full(res.code_raster.get_data(), res.code_raster.get_meta(), c_path,
                              out_fmt);
            std::cout << "Wrote " << c_path << "\n";
        }
        {
            const std::string l_path = pp.make_output_path("_lithology_legend", ".csv");
            std::ofstream csv(l_path);
            csv << "code,name,erodibility,n_pixels,map_values\n";
            csv << std::setprecision(7);
            for (const auto& row : res.legend)
            {
                std::string mv;
                for (const auto& v : row.map_values) mv += (mv.empty() ? "" : "; ") + v;
                std::string q = "\"";
                for (char c : mv) { if (c == '"') q += "\"\""; else q += c; }
                q += "\"";
                csv << row.code << ",\"" << row.name << "\"," << row.erodibility << ","
                    << row.n_pixels << "," << q << "\n";
            }
            std::cout << "Wrote " << l_path << "\n";
        }
        {
            const std::string p_path = pp.make_output_path("_lem_erodibility_parameters", ".txt");
            std::ofstream txt(p_path);
            txt << "# lsdtt-lem parameters for " << k_path << "\n"
                << "# K is in the property table's convention; lsdtt-lem converts it to the "
                   "run's.\n";
            if (grid_source == "lithology")
                txt << "# The grid is the lithology map's own: run lsdtt-lem with\n"
                       "#   model_grid_source: erodibility_raster\n";
            else if (grid_source == "parameters")
                txt << "# The grid is synthetic: use the same n_rows, n_cols, cell_size_metres "
                       "and synthetic_crs_wkt in lsdtt-lem.\n";
            else
                txt << "# The grid is the DEM's: run lsdtt-lem on the same DEM "
                       "(read_initial_dem: true).\n";
            txt << lsdtt::lithology_erodibility_lem_parameters(
                table, std::filesystem::path(k_path).filename().string());
            std::cout << "Wrote " << p_path << "\n\nFor lsdtt-lem:\n"
                      << lsdtt::lithology_erodibility_lem_parameters(
                             table, std::filesystem::path(k_path).filename().string());
            if (grid_source == "lithology") std::cout << "model_grid_source: erodibility_raster\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README", k_version, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-lithology-erodibility: done.\n";
    return 0;
}
