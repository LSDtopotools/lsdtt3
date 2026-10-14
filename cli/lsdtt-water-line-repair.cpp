// lsdtt-water-line-repair.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-water-line-repair
//
// Reconnects stream/water polylines that have been truncated where roads cross
// them (an OS VectorMap Local artefact). Reads a water LINESTRING layer and a
// road centreline layer from a vector dataset, bridges facing dangling water
// tips whose gap crosses a road, and writes a new line layer (originals + the
// inserted bridge segments, with typed provenance attributes).
//
// Pure vector: no DEM required. Input layers must share the same CRS.
//
// Usage:
//   lsdtt-water-line-repair <path> <parameter_file>
//   lsdtt-water-line-repair <parameter_file>
//   lsdtt-water-line-repair -h | -v | -c
// -----------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_vector_line_repair_cpu.hpp"
#include "../src/cpu/lsdtt_vector_ops_cpu.hpp"
#include "../src/utils/lsdtt_linear_algebra.hpp"

#include <cstdio>
#include <fstream>
#include <limits>
#include <set>

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-water-line-repair";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-water-line-repair                  ||\n"
        << "|| Reconnects water lines truncated by roads           ||\n"
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

int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    
    // -------------------------------------------------------------------------
    // Parameters (five-map pattern)
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    string_default_map["input_vector_fname"] = "";
    help_map["input_vector_fname"] = {"string","",
        "Input vector dataset (e.g. a GeoPackage) holding both layers.",
        "Required. A relative path is read from read_path. Both layers must "
        "share the same CRS (projected metres recommended).",
        "Input vectors", "io", ""};

    string_default_map["water_layer"] = "";
    help_map["water_layer"] = {"string","",
        "Name of the water LINESTRING layer to repair.",
        "Required. e.g. vml-nr73sw_Water_Line.",
        "Input vectors", "param", ""};

    string_default_map["road_layer"] = "";
    help_map["road_layer"] = {"string","",
        "Name of the road centreline LINESTRING layer used as the crossing test.",
        "Required. e.g. vml-nr73sw_Road_Centreline.",
        "Input vectors", "param", ""};

    string_default_map["road_vector_fname"] = "";
    help_map["road_vector_fname"] = {"string","",
        "Dataset holding the road layer, if different from input_vector_fname.",
        "Use when the water network is in one file (e.g. the centreline .fgb) "
        "and roads are in another (e.g. the source .gpkg). Empty = same as "
        "input_vector_fname. A relative path is read from read_path.",
        "Input vectors", "io", ""};

    string_default_map["road_name_field"] = "roadName";
    help_map["road_name_field"] = {"string","roadName",
        "Road attribute copied into the bridge crossed_road_name field.",
        "Falls back to 'roadNumber' when this field is empty/absent.",
        "Input vectors", "param", ""};

    string_default_map["write_vector_format"] = "flatgeobuf";
    help_map["write_vector_format"] = {"string","flatgeobuf",
        "Output VECTOR format: flatgeobuf or gpkg.",
        "Governs every vector layer this tool writes. Distinct from the raster "
        "key write_output_format used by the raster tools: this CLI writes "
        "vectors only. FlatGeobuf is preferred: it is compact, single-layer, "
        "carries its own georeferencing, and does not have GeoPackage's "
        "problems writing across drives (a common WSL2 situation). GeoPackage "
        "is an SQLite database, and SQLite cannot take the file locks it needs "
        "on a drive mounted into WSL2, so a gpkg write there fails.",
        "Outputs", "io", ""};

    float_default_map["snap_tolerance_metres"] = 0.5f;
    help_map["snap_tolerance_metres"] = {"float","0.5",
        "Endpoint coincidence tolerance (m) for degree counting.",
        "Endpoints within this distance are treated as the same node, so shared "
        "confluences are not mistaken for dangling tips.",
        "Gap bridging", "param", ""};

    float_default_map["max_gap_length_metres"] = 25.0f;
    help_map["max_gap_length_metres"] = {"float","25.0",
        "Maximum ROAD-crossing bridge length (m); longer road gaps are rejected.",
        "Roads are narrow, but raise for wide / dual-carriageway crossings or "
        "wide verges where a true crossing gap exceeds the default.",
        "Gap bridging", "param", ""};

    float_default_map["max_nonroad_gap_metres"] = 30.0f;
    help_map["max_nonroad_gap_metres"] = {"float","30.0",
        "Maximum bridge length (m) for gaps NOT associated with a road.",
        "Fills digitising breaks regardless of roads. The STRICT facing/heading "
        "gates still apply, so genuine headwaters and parallel streams are not "
        "joined. Set 0 to disable non-road bridging.",
        "Gap bridging", "param", ""};

    float_default_map["max_facing_angle_degrees"] = 45.0f;
    help_map["max_facing_angle_degrees"] = {"float","45.0",
        "Strict tier: max angle (deg) EACH tip's tangent may deviate from the "
        "line to its partner. Used for non-road gaps.",
        "Both tips must point toward each other across the gap.",
        "Gap bridging", "param", ""};

    float_default_map["max_heading_misalignment_degrees"] = 50.0f;
    help_map["max_heading_misalignment_degrees"] = {"float","50.0",
        "Strict tier: max deviation from anti-parallel (deg) between the two "
        "tips' tangents. Used for non-road gaps.",
        "A true continuation has the two outward tangents ~opposite (180 deg). "
        "Suppresses bridges between parallel unrelated streams.",
        "Gap bridging", "param", ""};

    float_default_map["road_max_facing_angle_degrees"] = 60.0f;
    help_map["road_max_facing_angle_degrees"] = {"float","60.0",
        "Road tier: a road-associated bridge only needs AT LEAST ONE tip to "
        "face within this angle (no heading constraint).",
        "The channel kinks at a road crossing so the dangling tip's own tangent "
        "is often unreliable; the partner across the road faces well.",
        "Gap bridging", "param", ""};

    float_default_map["tangent_fit_length_metres"] = 15.0f;
    help_map["tangent_fit_length_metres"] = {"float","15.0",
        "Inward arc length (m) used to estimate each tip's outward tangent.",
        "Longer values smooth over digitising wiggle near the tip.",
        "Gap bridging", "param", ""};

    float_default_map["road_intersection_tolerance_metres"] = 10.0f;
    help_map["road_intersection_tolerance_metres"] = {"float","10.0",
        "Road-association distance (m): a bridge is 'road-associated' if its "
        "connector comes within this distance of a road centreline.",
        "Covers staggered severances where the stub ends sit a few metres off "
        "the road. Uses segment-to-segment distance, not just a midpoint test.",
        "Gap bridging", "param", ""};

    bool_default_map["require_mutual_best"] = true;
    help_map["require_mutual_best"] = {"bool","true",
        "Require each tip to be the other's best partner before bridging.",
        "Strongly suppresses spurious bridges; disable only to be more inclusive.",
        "Gap bridging", "param", ""};

    bool_default_map["merge_chains"] = false;
    help_map["merge_chains"] = {"bool","false",
        "Dissolve connected chains into single features (NOT yet implemented).",
        "Leave false: keep originals + add bridge lines. Coincident endpoints "
        "already make the network topologically continuous downstream.",
        "Gap bridging", "param", ""};

    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool","true",
        "Build a spatial index when writing the output layer.",
        "Speeds GIS queries; buffers the layer in memory before writing.",
        "Outputs", "io", ""};

    bool_default_map["write_tip_diagnostics"] = false;
    help_map["write_tip_diagnostics"] = {"bool","false",
        "Write a diagnostic point layer of every dangling tip.",
        "Output: <write_fname>_repair_tips.fgb. Attributes: bridged (1 if the "
        "tip was joined), road_dist_m (distance to nearest road), line_index. "
        "Use it to see which tips near roads stayed unbridged and why.",
        "Outputs", "switch", ""};

    // -------------------------------------------------------------------------
    // Parse
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

    // write_output_format meant a VECTOR format here, but that name means a
    // RASTER format in every other CLI. Renamed to write_vector_format.
    pp.set_deprecated_aliases({ {"write_output_format", "write_vector_format"} });

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto sf = pp.get_string_parameters();

    if (sf["input_vector_fname"].empty())
    {
        std::cerr << "FATAL: input_vector_fname is required.\n";
        return 1;
    }
    if (sf["water_layer"].empty() || sf["road_layer"].empty())
    {
        std::cerr << "FATAL: water_layer and road_layer are required.\n";
        return 1;
    }

    const std::string input_path = pp.add_read_path_if_needed(sf["input_vector_fname"]);

    lsdtt::vector_output_format out_fmt =
        lsdtt::resolve_vector_output_format(sf["write_vector_format"]);

    lsdtt::water_line_repair_options opt;
    opt.snap_tolerance_metres            = static_cast<double>(ff["snap_tolerance_metres"]);
    opt.max_gap_length_metres            = static_cast<double>(ff["max_gap_length_metres"]);
    opt.max_nonroad_gap_metres           = static_cast<double>(ff["max_nonroad_gap_metres"]);
    opt.max_facing_angle_degrees         = static_cast<double>(ff["max_facing_angle_degrees"]);
    opt.max_heading_misalignment_degrees = static_cast<double>(ff["max_heading_misalignment_degrees"]);
    opt.road_max_facing_angle_degrees    = static_cast<double>(ff["road_max_facing_angle_degrees"]);
    opt.tangent_fit_length_metres        = static_cast<double>(ff["tangent_fit_length_metres"]);
    opt.road_intersection_tolerance_metres =
        static_cast<double>(ff["road_intersection_tolerance_metres"]);
    opt.require_mutual_best = bf["require_mutual_best"];
    opt.merge_chains        = bf["merge_chains"];
    opt.road_name_field     = sf["road_name_field"];

    // -------------------------------------------------------------------------
    // Load layers
    // -------------------------------------------------------------------------
    std::cout << "Opening vector: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_vector_layer water;
    lsdtt::lsdtt_vector_layer roads;
    std::string water_crs;
    std::string road_crs;
    const std::string road_path = sf["road_vector_fname"].empty()
        ? input_path
        : pp.add_read_path_if_needed(sf["road_vector_fname"]);

    try
    {
        // Water layer (from input_vector_fname).
        loader.open_vector(input_path);
        {
            const std::vector<std::string> layers = loader.get_layer_names();
            bool found = false;
            for (const auto& l : layers) if (l == sf["water_layer"]) found = true;
            if (!found)
                throw std::runtime_error("water_layer not found: " + sf["water_layer"]);
        }
        water     = loader.read_vector_layer(sf["water_layer"]);
        water_crs = loader.get_vector_layer_crs_wkt(sf["water_layer"]);
        loader.close_vector();

        // Road layer (from road_vector_fname, possibly a different dataset).
        lsdtt::lsdtt_data_loader road_loader;
        road_loader.open_vector(road_path);
        {
            const std::vector<std::string> layers = road_loader.get_layer_names();
            bool found = false;
            for (const auto& l : layers) if (l == sf["road_layer"]) found = true;
            if (!found)
                throw std::runtime_error("road_layer not found: " + sf["road_layer"]);
        }
        roads    = road_loader.read_vector_layer(sf["road_layer"]);
        road_crs = road_loader.get_vector_layer_crs_wkt(sf["road_layer"]);
        road_loader.close_vector();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not load layers: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Water lines: " << water.size()
              << "   Road lines: " << roads.size() << "\n";

    if (!lsdtt::crs_wkt_equivalent(water_crs, road_crs))
    {
        std::cerr << "FATAL: water_layer and road_layer have different CRS. "
                  << "Reproject them to a common CRS first (see "
                  << "lsdtt_coordinate_transform in lsdtt_vector_utils).\n";
        return 1;
    }
    if (water_crs.empty())
    {
        std::cerr << "Warning: input layers have no CRS WKT; output will be "
                  << "written without georeferencing.\n";
    }

    // -------------------------------------------------------------------------
    // Repair
    // -------------------------------------------------------------------------
    lsdtt::water_line_repair_result result;
    try
    {
        result = lsdtt::repair_water_lines(water, roads, opt);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: repair failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "\n--- Repair summary ---\n";
    std::cout << "Dangling tips        : " << result.n_dangling_tips << "\n";
    std::cout << "Road-crossed pairs   : " << result.n_candidate_pairs << "\n";
    std::cout << "Bridges accepted     : " << result.n_bridges_accepted << "\n";
    std::cout << "Output features      : " << result.output_layer.size() << "\n";

    // -------------------------------------------------------------------------
    // Write output
    // -------------------------------------------------------------------------
    const std::string ext = lsdtt::output_extension(out_fmt);
    const std::string out_path =
        lsdtt::lsdtt_data_loader::resolve_writable_path(pp.get_full_write_path() + "_water_line_repaired", ext);

    std::cout << "\nWriting: " << out_path << "\n";
    try
    {
        lsdtt::lsdtt_data_loader writer;
        writer.begin_vector_write(
            out_path, "water_line_repaired", result.output_layer.schema,
            water_crs, out_fmt, bf["calculate_spatial_index"]);

        for (const auto& rec : result.output_layer.records)
        {
            const auto* line =
                std::get_if<lsdtt::lsdtt_line_geom>(&rec.geometry);
            if (!line) continue;
            writer.write_line(line->points, rec.field_values);
        }
        writer.finish_vector_write();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not write output: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Optional tip diagnostics
    // -------------------------------------------------------------------------
    if (bf["write_tip_diagnostics"])
    {
        std::set<std::size_t> bridged;
        for (const auto& b : result.accepted_bridges)
        {
            bridged.insert(b.tip_index_a);
            bridged.insert(b.tip_index_b);
        }

        auto nearest_road_dist = [&](const lsdtt::lsdtt_xy& p) -> double
        {
            double best = std::numeric_limits<double>::infinity();
            for (const auto& rec : roads.records)
            {
                const auto* line =
                    std::get_if<lsdtt::lsdtt_line_geom>(&rec.geometry);
                if (!line || line->points.size() < 2) continue;
                for (std::size_t k = 1; k < line->points.size(); ++k)
                {
                    const double d = lsdtt::point_to_segment_distance(
                        p.x, p.y,
                        line->points[k - 1].x, line->points[k - 1].y,
                        line->points[k].x, line->points[k].y);
                    if (d < best) best = d;
                }
            }
            return best;
        };

        lsdtt::lsdtt_vector_schema schema;
        schema.field_names = {"bridged", "road_dist_m", "line_index",
                              "at_start", "tx", "ty"};
        schema.field_types = {
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::double_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::double_type,
            lsdtt::lsdtt_vector_field_type::double_type
        };

        const std::string tip_path =
            lsdtt::lsdtt_data_loader::resolve_writable_path(pp.get_full_write_path() + "_repair_tips", ".fgb");
        std::cout << "\nWriting tip diagnostics: " << tip_path << "\n";
        int n_unbridged_near_road = 0;
        try
        {
            lsdtt::lsdtt_data_loader tw;
            tw.begin_vector_write(tip_path, "repair_tips", schema, water_crs,
                                  lsdtt::vector_output_format::flatgeobuf, false);
            for (std::size_t i = 0; i < result.tips.size(); ++i)
            {
                const auto& t = result.tips[i];
                const double rd = nearest_road_dist(t.xy);
                const std::int64_t is_bridged = bridged.count(i) ? 1 : 0;
                if (!is_bridged && rd <= 30.0) ++n_unbridged_near_road;
                tw.write_point(t.xy.x, t.xy.y, {
                    is_bridged, rd,
                    static_cast<std::int64_t>(t.line_index),
                    static_cast<std::int64_t>(t.at_start ? 1 : 0),
                    t.tx, t.ty
                });
            }
            tw.finish_vector_write();
            std::cout << "  unbridged tips within 30 m of a road: "
                      << n_unbridged_near_road << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: could not write tip diagnostics: "
                      << e.what() << "\n";
        }
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-water-line-repair: done.\n";
    return 0;
}
