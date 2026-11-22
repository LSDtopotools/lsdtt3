// lsdtt-hillslope-channel-tagging.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-hillslope-channel-tagging
//
// Transposes channel information onto hillslope pixels. Every hillslope pixel is
// linked to the first channel pixel on its D8 path, whichever ridge method is
// used, and from that link the tool reports:
//
//   * HAND, the relief from the channel to a ridge, and relative height
//     (HAND / ridge relief);
//   * tags: the chi or k_sn of the channel pixel, or any imported channel
//     attribute such as sinuosity or planform curvature.
//
// Usage:
//   lsdtt-hillslope-channel-tagging <path> <parameter_file>
//   lsdtt-hillslope-channel-tagging <parameter_file>   (path assumed ./)
//   lsdtt-hillslope-channel-tagging -h                 (help)
//   lsdtt-hillslope-channel-tagging -v                 (version)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_channel_attribute_input.hpp"
#include "../src/core/lsdtt_channel_network_ingest.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_hillslope_channel_tagging_cpu.hpp"
#include "../src/cpu/lsdtt_outlet_resolution_cpu.hpp"
#include "../src/cpu/lsdtt_raster_breach_least_cost_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = lsdtt::lsdtt3_citation_url();
static const std::string k_tool_name = "lsdtt-hillslope-channel-tagging";

namespace
{

void print_banner()
{
    std::cout
        << "\n=========================================================\n"
        << "|| Welcome to the lsdtt3 hillslope-channel tagging tool  ||\n"
        << "|| I transpose channel information onto hillslopes.      ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| Docs     : " << lsdtt::lsdtt3_cli_documentation_url(k_tool_name) << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

std::vector<std::string> split_columns(const std::string& text)
{
    std::vector<std::string> out;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ','))
    {
        const auto first = token.find_first_not_of(" \t");
        if (first == std::string::npos) continue;
        const auto last = token.find_last_not_of(" \t");
        out.push_back(token.substr(first, last - first + 1));
    }
    return out;
}

// A column name becomes part of an output file name.
void check_column_name(const std::string& name)
{
    for (const unsigned char c : name)
        if (!std::isalnum(c) && c != '_' && c != '-')
            throw std::invalid_argument(
                "channel_attribute_columns: '" + name +
                "' can only contain letters, digits, '_' and '-', because it names an output file");
}

void declare_parameters(std::map<std::string, bool>&        bools,
                        std::map<std::string, float>&       floats,
                        std::map<std::string, int>&         ints,
                        std::map<std::string, std::string>& strings,
                        lsdtt::lsdtt_help_registry&         help)
{
    lsdtt::declare_routing_core_parameters(ints, strings, help);
    strings["boundary_conditions"] = "oooo";
    help["boundary_conditions"][1] = "oooo";

    lsdtt::declare_min_elevation_parameter(floats, help, -9999.0f, "Input");
    lsdtt::declare_interior_nodata_parameters(ints, strings, help);

    lsdtt::declare_breach_before_fill_parameters(bools, ints, floats, strings, help,
                                                 lsdtt::k_breach_lc_default_max_dist_outlet);
    lsdtt::declare_fill_before_routing_parameter(bools, help);
    lsdtt::declare_fill_min_slope_parameter(floats, help);

    lsdtt::declare_channel_initiation_parameters(ints, help);
    strings["channel_network_fname"] = "";
    help["channel_network_fname"] = {"string", "",
        "Optional channel network to ingest instead of using the threshold.",
        "A point file written by lsdtt-channel-extraction "
        "(write_channel_network_points), with row and col. It must have been "
        "extracted on this DEM's grid.",
        "Channel network", "param", ""};

    strings["hillslope_ridge_method"] = "farthest";
    help["hillslope_ridge_method"] = {"string", "farthest",
        "How each hillslope pixel's ridge is chosen: farthest, highest or aspect_ray.",
        "farthest: the pixel with the longest D8 path to the channel in the pixel's "
        "group. highest: the group's highest pixel. aspect_ray: a straight ray from "
        "the pixel, upslope along its aspect, stopping at the first drainage divide "
        "that is at or above the pixel and is a crest along the ray; a ray reaching "
        "a channel, nodata or the grid edge first uses the highest cell on the ray. "
        "The channel a pixel is linked to is the D8 one for all three. With farthest "
        "and a shared group, relative height can exceed 1: it is not clipped.",
        "Ridges", "param", ""};
    strings["hillslope_relief_grouping"] = "channel_pixel";
    help["hillslope_relief_grouping"] = {"string", "channel_pixel",
        "The pixels sharing a ridge, for farthest and highest.",
        "channel_pixel: all pixels draining to one channel pixel. entry_pixel: all "
        "pixels entering the channel through the same last hillslope pixel. pixel: "
        "the pixel and those draining through it; about half of the hillslope pixels "
        "of a DEM have no donors, so these are their own ridge with relative height 1. "
        "Not used by aspect_ray.",
        "Ridges", "param", "hillslope_ridge_method"};
    strings["hillslope_aspect_method"] = "polyfit";
    help["hillslope_aspect_method"] = {"string", "polyfit",
        "Aspect used for aspect_ray: polyfit or dinf.",
        "polyfit: from the fitted polynomial surface, at polyfit_lambda_metres. dinf: "
        "the D-infinity facet direction, which is noisier.",
        "Ridges", "param", "hillslope_ridge_method"};
    floats["polyfit_lambda_metres"] = 70.0f;
    help["polyfit_lambda_metres"] = {"float", "70.0",
        "Polyfit length scale (half-wavelength) for the aspect of aspect_ray (m).",
        "The fitting window follows from this half-wavelength, as in "
        "lsdtt-surface-metrics. 0 or less uses three cells.",
        "Ridges", "param", "hillslope_ridge_method"};
    bools["polyfit_on_routing_surface"] = false;
    help["polyfit_on_routing_surface"] = {"bool", "false",
        "Fit the aspect to the routing surface rather than the input DEM.",
        "The routing surface has carved and filled areas whose slopes come from the "
        "conditioning rather than the landscape.",
        "Ridges", "switch", "hillslope_ridge_method"};

    bools["tag_chi"] = false;
    help["tag_chi"] = {"bool", "false",
        "Tag hillslope pixels with the chi of their channel pixel.",
        "Needs basin_outlet_fname or chi_max_outlet_elevation_metres, because chi "
        "needs a true baselevel. Output: <write_fname>_tag_chi.",
        "Chi and k_sn tags", "switch", ""};
    bools["tag_k_sn"] = false;
    help["tag_k_sn"] = {"bool", "false",
        "Tag hillslope pixels with the k_sn of their channel pixel.",
        "Needs basin_outlet_fname or chi_max_outlet_elevation_metres, because chi "
        "needs a true baselevel. Output: <write_fname>_tag_k_sn.",
        "Chi and k_sn tags", "switch", ""};
    floats["m_over_n"] = 0.45f;
    help["m_over_n"] = {"float", "0.45", "Concavity index for chi.", "",
        "Chi and k_sn tags", "param", "tag_chi"};
    floats["chi_A0"] = 1.0f;
    help["chi_A0"] = {"float", "1.0", "Reference drainage area A0 for chi (m2).",
        "Leave at 1 so that the slope of a chi profile is k_sn.",
        "Chi and k_sn tags", "param", "tag_chi"};
    floats["fit_critical_divergence_metres"] = 5.0f;
    help["fit_critical_divergence_metres"] = {"float", "5.0",
        "Critical elevation divergence (m) used to split chi-elevation branches.", "",
        "Chi and k_sn tags", "param", "tag_chi"};
    lsdtt::declare_basin_outlet_parameters(bools, ints, floats, strings, help);
    floats["chi_max_outlet_elevation_metres"] = 0.0f;
    help["chi_max_outlet_elevation_metres"] = {"float", "not set",
        "Use every baselevel basin with its outlet at or below this elevation (m).",
        "An alternative to basin_outlet_fname; give one of them. On a DEM with the sea "
        "masked, a few metres keeps the basins that reach the coast. A basin that a "
        "clipped DEM cuts off has a high outlet and is left out. Basins touching the "
        "grid edge are always left out. Their hillslopes are nodata and counted.",
        "Chi and k_sn tags", "param", "tag_chi"};

    strings["channel_attributes_fname"] = "";
    help["channel_attributes_fname"] = {"string", "",
        "CSV or vector file of channel points or lines with numeric attributes.",
        "With row, col and flat columns the join is exact. Otherwise CSV needs "
        "latitude and longitude (EPSG:4326), and points or lines in a vector file are "
        "reprojected and snapped to the nearest channel pixel.",
        "Imported channel attributes", "param", ""};
    strings["channel_attribute_columns"] = "";
    help["channel_attribute_columns"] = {"string", "",
        "Comma-separated numeric columns to tag onto hillslopes.",
        "Each gives one raster, <write_fname>_tag_<column>. A missing or text column "
        "stops the run. Several records on a channel pixel give their median.",
        "Imported channel attributes", "param", "channel_attributes_fname"};
    floats["channel_attribute_max_snap_distance_metres"] = 60.0f;
    help["channel_attribute_max_snap_distance_metres"] = {"float", "60.0",
        "Furthest a record can be from the channel pixel it joins (m).",
        "Records further from every channel pixel are dropped and counted.",
        "Imported channel attributes", "param", "channel_attributes_fname"};
    floats["channel_attribute_max_gap_fill_distance_metres"] = 150.0f;
    help["channel_attribute_max_gap_fill_distance_metres"] = {"float", "150.0",
        "Furthest along the channel a record's value can fill a gap (m).",
        "A channel pixel with no record takes the nearest record on its own "
        "junction-to-junction link by D8 channel distance. Filling never crosses a "
        "junction, and equidistant records give their median. 0 turns filling off.",
        "Imported channel attributes", "param", "channel_attributes_fname"};
    floats["channel_attribute_nodata_value"] = -9999.0f;
    help["channel_attribute_nodata_value"] = {"float", "-9999",
        "Attribute value meaning 'no data'.",
        "Empty and NaN values are always nodata as well.",
        "Imported channel attributes", "param", "channel_attributes_fname"};

    lsdtt::declare_raster_output_format_parameter(strings, help);
    const auto output = [&](const std::string& key, bool value, const std::string& what,
                            const std::string& suffix) {
        bools[key] = value;
        help[key] = {"bool", value ? "true" : "false", "Write " + what + ".",
                     "Output: <write_fname>" + suffix + ".", "Output", "switch", ""};
    };
    output("write_hand_raster", false, "height above the nearest drainage (m)", "_HAND");
    output("write_ridge_relief_raster", true, "relief from channel to ridge (m)", "_ridge_relief");
    output("write_relative_height_raster", true, "HAND / ridge relief", "_relative_height");
    output("write_ridge_distance_raster", false,
           "ridge distance (m): the D8 path from ridge to channel, or for aspect_ray the "
           "straight line from pixel to ridge", "_ridge_distance");
    output("write_ridge_status_raster", false,
           "how each ridge was found (0 channel, 1 no channel, 2 D8, 3 crest on a divide, "
           "4 ray reached a channel, 5 ray reached nodata or the edge, 6 no aspect)",
           "_ridge_status");
    output("write_channel_pixel_summary", false,
           "a CSV with one row per channel pixel: hillslope pixels and median relative height",
           "_channel_pixel_summary.csv");
}

void report_tag(const std::string& name, const lsdtt::hillslope_tag_result& t)
{
    std::cout << "  " << name << ": " << t.n_hillslope_tagged << " hillslope pixels tagged, "
              << t.n_no_value << " on a channel pixel with no value, " << t.n_excluded
              << " in excluded basins, " << t.n_no_channel << " reaching no channel\n";
    const std::string key = "tag_" + name;
    lsdtt::lsdtt_run_record::add_provenance(key + "_n_hillslope_tagged", std::to_string(t.n_hillslope_tagged));
    lsdtt::lsdtt_run_record::add_provenance(key + "_n_no_value", std::to_string(t.n_no_value));
    lsdtt::lsdtt_run_record::add_provenance(key + "_n_excluded", std::to_string(t.n_excluded));
}

} // namespace

int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();
    const lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry         help_map;
    declare_parameters(bool_default_map, float_default_map, int_default_map,
                       string_default_map, help_map);

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(k_tool_name, k_status, k_citation,
                                                       help_map, cli_args);
        return 0;
    }

    lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);

    try
    {
        lsdtt::lsdtt_param_parser pp(cli_args.path, cli_args.fname);
        pp.parse_all_parameters(float_default_map, int_default_map, bool_default_map,
                                string_default_map);
        const auto bf = pp.get_bool_parameters();
        const auto ff = pp.get_float_parameters();
        const auto ip = pp.get_int_parameters();
        const auto sf = pp.get_string_parameters();

        // ---- validate before any raster work --------------------------------
        const std::string ridge = sf.at("hillslope_ridge_method");
        if (ridge != "farthest" && ridge != "highest" && ridge != "aspect_ray")
            throw std::invalid_argument("hillslope_ridge_method must be farthest, highest or aspect_ray");
        const std::string grouping = sf.at("hillslope_relief_grouping");
        if (grouping != "channel_pixel" && grouping != "entry_pixel" && grouping != "pixel")
            throw std::invalid_argument(
                "hillslope_relief_grouping must be channel_pixel, entry_pixel or pixel");
        if (ridge == "aspect_ray" && grouping != "channel_pixel")
            throw std::invalid_argument(
                "hillslope_relief_grouping does not apply to aspect_ray: each ray has its own ridge");
        const std::string aspect_method = sf.at("hillslope_aspect_method");
        if (aspect_method != "polyfit" && aspect_method != "dinf")
            throw std::invalid_argument("hillslope_aspect_method must be polyfit or dinf");
        if (ip.at("threshold_contributing_pixels") <= 0)
            throw std::invalid_argument("threshold_contributing_pixels must be > 0");

        const bool want_chi_tags = bf.at("tag_chi") || bf.at("tag_k_sn");
        const std::string outlet_file = sf.at("basin_outlet_fname");
        const bool has_outlet_file = !outlet_file.empty() && outlet_file != "NULL";
        const bool has_outlet_limit = pp.parameter_supplied("chi_max_outlet_elevation_metres");
        if (want_chi_tags && has_outlet_file == has_outlet_limit)
            throw std::invalid_argument(
                "tag_chi and tag_k_sn need exactly one of basin_outlet_fname and "
                "chi_max_outlet_elevation_metres: chi needs a true baselevel, and a DEM's "
                "nodata cannot separate a coast from a clipped edge");
        if (want_chi_tags && (ff.at("m_over_n") <= 0.0f || ff.at("chi_A0") <= 0.0f))
            throw std::invalid_argument("m_over_n and chi_A0 must be > 0");

        const std::string attributes_file = sf.at("channel_attributes_fname");
        const std::vector<std::string> attribute_columns = split_columns(sf.at("channel_attribute_columns"));
        if (attributes_file.empty() != attribute_columns.empty())
            throw std::invalid_argument(
                "channel_attributes_fname and channel_attribute_columns must be given together");
        for (const auto& name : attribute_columns)
        {
            check_column_name(name);
            if ((name == "chi" && bf.at("tag_chi")) || (name == "k_sn" && bf.at("tag_k_sn")))
                throw std::invalid_argument("channel_attribute_columns: '" + name +
                                            "' would overwrite the computed " + name + " tag");
        }
        const double max_snap = ff.at("channel_attribute_max_snap_distance_metres");
        const double max_gap = ff.at("channel_attribute_max_gap_fill_distance_metres");
        if (max_snap < 0.0 || max_gap < 0.0)
            throw std::invalid_argument(
                "channel_attribute_max_snap_distance_metres and "
                "channel_attribute_max_gap_fill_distance_metres must not be negative");

        const std::string breach_choice = sf.at("breach_algorithm_choice");
        if (breach_choice != "least_cost" && breach_choice != "lindsay2016")
            throw std::invalid_argument("breach_algorithm_choice must be least_cost or lindsay2016");
        lsdtt::breach_least_cost_options breach_options;
        if (bf.at("breach_before_fill") && breach_choice == "least_cost")
        {
            breach_options.max_cost = ff.at("breach_lc_max_cost") >= std::numeric_limits<float>::max()
                                          ? std::numeric_limits<double>::infinity()
                                          : static_cast<double>(ff.at("breach_lc_max_cost"));
            breach_options.minimize_dist  = bf.at("breach_lc_minimize_dist");
            breach_options.flat_increment = static_cast<double>(ff.at("breach_lc_flat_increment"));
            breach_options.fill_remainder = false;
            breach_options.edge_policy =
                lsdtt::parse_breach_lc_edge_policy(sf.at("breach_lc_edge_policy"), "breach_lc_edge_policy");
            breach_options.max_dist = pp.parameter_supplied("breach_lc_max_dist")
                                          ? ip.at("breach_lc_max_dist")
                                          : lsdtt::default_breach_lc_max_dist(breach_options.edge_policy);
            lsdtt::validate_breach_least_cost_options(breach_options);
        }

        const lsdtt::raster_output_format raster_format =
            lsdtt::resolve_raster_output_format(sf.at("write_output_format"));
        const std::string raster_ext = lsdtt::output_extension(raster_format, pp.get_read_ext());
        const int n_threads = ip.at("n_threads");

        // ---- DEM -------------------------------------------------------------
        lsdtt::lsdtt_data_loader loader;
        const std::string input_path =
            pp.get_read_path() + pp.get_read_fname() +
            (pp.get_read_ext().empty() ? std::string(".tif") : pp.get_read_ext());
        loader.open_raster(input_path);
        loader.print_raster_info();
        loader.get_meta().check_projected_metres();
        lsdtt::lsdtt_raster_cpu input_dem(loader);

        const std::size_t n_masked = lsdtt::apply_min_elevation(input_dem, ff.at("min_elevation"));
        if (n_masked > 0) std::cout << "Masked " << n_masked << " pixels at or below min_elevation\n";
        std::cout << lsdtt::apply_interior_nodata_method(
                         input_dem, sf.at("interior_nodata_method"), sf.at("nodata_fill_method_choice"),
                         ip.at("nodata_fill_max_iter"), n_threads)
                  << "\n";

        std::cout << "\n--- Routing surface ---\n";
        lsdtt::lsdtt_raster_cpu routing = input_dem;
        if (bf.at("breach_before_fill"))
        {
            if (breach_choice == "least_cost")
            {
                lsdtt::breach_least_cost_diagnostics diag;
                routing = lsdtt::breach_depressions_least_cost(routing, breach_options, &diag);
                std::cout << "Least-cost breaching: pits " << diag.pits_detected << ", solved "
                          << diag.pits_solved << "\n";
            }
            else
            {
                routing = lsdtt::breach_depressions(routing, n_threads);
            }
        }
        if (bf.at("fill_before_routing"))
            routing = lsdtt::fill_pits(routing, ff.at("fill_min_slope"), sf.at("boundary_conditions"), n_threads);

        const lsdtt::lsdtt_raster_meta& meta = routing.get_meta();
        const std::size_t n_cells = static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

        lsdtt::lsdtt_flowinfo_options flow_options;
        flow_options.boundary_conditions     = sf.at("boundary_conditions");
        flow_options.build_topological_order = true;
        flow_options.build_s_vector          = true;
        const lsdtt::lsdtt_flowinfo flow(routing.get_data_array(), routing.get_mask_array(), meta, flow_options);

        // ---- channels ---------------------------------------------------------
        std::cout << "\n--- Channel network ---\n";
        const int threshold = ip.at("threshold_contributing_pixels");
        std::vector<std::uint8_t> channel(n_cells, 0u);
        std::int64_t n_channel = 0;
        if (!sf.at("channel_network_fname").empty())
        {
            const std::string path = pp.add_read_path_if_needed(sf.at("channel_network_fname"));
            lsdtt::lsdtt_ingested_channel_network ingested = lsdtt::ingest_channel_network_points(path, flow);
            channel = std::move(ingested.mask);
            n_channel = ingested.n_channel;
            std::cout << "Ingested " << ingested.n_points << " points -> " << n_channel << " channel pixels";
            if (ingested.n_on_nodata > 0) std::cout << " (" << ingested.n_on_nodata << " on nodata, dropped)";
            std::cout << "\n";
        }
        else
        {
            const lsdtt::lsdtt_raster_result contributing = flow.get_n_contributing_pixels_raster();
            for (std::size_t i = 0; i < n_cells; ++i)
                if (contributing.mask[i] && contributing.data[i] >= static_cast<float>(threshold))
                { channel[i] = 1u; ++n_channel; }
            std::cout << n_channel << " channel pixels at threshold " << threshold << "\n";
        }
        if (n_channel == 0)
            throw std::runtime_error("no channel pixels: lower threshold_contributing_pixels");

        // ---- relief ------------------------------------------------------------
        std::cout << "\n--- Ridges and relative height (" << ridge << ") ---\n";
        lsdtt::hillslope_relief_options relief_options;
        relief_options.ridge     = ridge;
        relief_options.grouping  = grouping;
        relief_options.n_threads = n_threads;
        lsdtt::lsdtt_raster_result ray_aspect;
        if (ridge == "aspect_ray")
        {
            const lsdtt::lsdtt_raster_cpu& aspect_surface =
                bf.at("polyfit_on_routing_surface") || aspect_method == "dinf" ? routing : input_dem;
            ray_aspect = lsdtt::compute_hillslope_ray_aspect(
                aspect_surface, aspect_method, ff.at("polyfit_lambda_metres"),
                sf.at("boundary_conditions"), n_threads);
            relief_options.ray_aspect = &ray_aspect;
        }
        const lsdtt::hillslope_relief_result relief = lsdtt::compute_hillslope_relief(
            flow, routing.get_data(), routing.get_mask(), channel, relief_options);

        std::cout << "Hillslope pixels: " << relief.n_hillslope << ", reaching no channel: "
                  << relief.n_no_channel << ", zero relief: " << relief.n_zero_relief
                  << ", relative height above 1: " << relief.n_above_ridge << "\n";
        if (ridge == "aspect_ray")
            std::cout << "Rays: crest on a divide " << relief.n_crest_divide << ", reached a channel "
                      << relief.n_ray_channel << ", reached nodata or the edge " << relief.n_ray_boundary
                      << ", no aspect " << relief.n_no_aspect << "\n";
        lsdtt::lsdtt_run_record::add_provenance("n_channel_pixels", std::to_string(relief.n_channel));
        lsdtt::lsdtt_run_record::add_provenance("n_hillslope_pixels", std::to_string(relief.n_hillslope));
        lsdtt::lsdtt_run_record::add_provenance("n_no_channel_pixels", std::to_string(relief.n_no_channel));
        lsdtt::lsdtt_run_record::add_provenance("n_zero_relief_pixels", std::to_string(relief.n_zero_relief));
        lsdtt::lsdtt_run_record::add_provenance("n_above_ridge_pixels", std::to_string(relief.n_above_ridge));

        const auto write_raster = [&](const lsdtt::lsdtt_raster_result& raster, const std::string& suffix) {
            const std::string path = pp.make_output_path(suffix, raster_ext);
            loader.write_full(raster.data, meta, path, raster_format);
            std::cout << "Wrote " << path << "\n";
        };

        if (bf.at("write_hand_raster"))            write_raster(relief.hand, "_HAND");
        if (bf.at("write_ridge_relief_raster"))    write_raster(relief.ridge_relief, "_ridge_relief");
        if (bf.at("write_relative_height_raster")) write_raster(relief.relative_height, "_relative_height");
        if (bf.at("write_ridge_distance_raster"))  write_raster(relief.ridge_distance, "_ridge_distance");
        if (bf.at("write_ridge_status_raster"))
        {
            std::vector<float> status(n_cells, static_cast<float>(meta.no_data_value));
            const auto& valid = routing.get_mask();
            for (std::size_t i = 0; i < n_cells; ++i)
                if (valid[i]) status[i] = static_cast<float>(relief.ridge_status[i]);
            const std::string path = pp.make_output_path("_ridge_status", raster_ext);
            loader.write_full(status, meta, path, raster_format);
            std::cout << "Wrote " << path << "\n";
        }

        if (bf.at("write_channel_pixel_summary"))
        {
            const auto summary =
                lsdtt::summarise_hillslopes_by_channel_pixel(relief.channel_flat, relief.relative_height);
            const lsdtt::lsdtt_coordinate_transform to_lonlat =
                lsdtt::lsdtt_coordinate_transform::from_epsg(meta.crs_wkt, 4326);
            const std::string path = pp.make_output_path("_channel_pixel_summary", ".csv");
            std::ofstream file(path);
            if (!file) throw std::runtime_error("cannot write " + path);
            file << "flat,row,col,latitude,longitude,n_hillslope_pixels,n_with_relative_height,median_relative_height\n";
            file << std::setprecision(10);
            for (std::size_t k = 0; k < summary.flat.size(); ++k)
            {
                const int row = static_cast<int>(summary.flat[k] / meta.n_cols);
                const int col = static_cast<int>(summary.flat[k] % meta.n_cols);
                const lsdtt::lsdtt_lonlat ll = lsdtt::cell_center_to_lonlat(meta, row, col, to_lonlat);
                file << summary.flat[k] << ',' << row << ',' << col << ',' << ll.latitude << ','
                     << ll.longitude << ',' << summary.n_hillslope[k] << ',' << summary.n_with_value[k] << ',';
                if (!std::isnan(summary.median[k])) file << summary.median[k];
                file << '\n';
            }
            std::cout << "Wrote " << path << "\n";
        }

        // ---- chi and k_sn ------------------------------------------------------
        if (want_chi_tags)
        {
            std::cout << "\n--- Chi and k_sn tags ---\n";
            lsdtt::hillslope_tag_basin_options basin_options;
            basin_options.minimum_contributing_pixels = threshold;
            if (has_outlet_file)
            {
                lsdtt::outlet_resolution_options outlet_options;
                outlet_options.snapping_method = sf.at("basin_outlet_snapping_method");
                const auto snap_threshold =
                    static_cast<std::uint64_t>(ip.at("basin_outlet_contributing_pixels_threshold"));
                outlet_options.square_snap.contributing_pixels_threshold = snap_threshold;
                outlet_options.square_snap.min_search_radius_pixels = ip.at("basin_outlet_snap_min_search_radius_pixels");
                outlet_options.square_snap.max_search_radius_pixels = ip.at("basin_outlet_snap_max_search_radius_pixels");
                outlet_options.downstream_snap.contributing_pixels_threshold = snap_threshold;
                lsdtt::outlet_resolution_report report;
                basin_options.outlet_flats = lsdtt::resolve_outlet_flats_from_file(
                    routing, flow, pp.add_read_path_if_needed(outlet_file), outlet_options, report);
            }
            else
            {
                basin_options.max_outlet_elevation_m = ff.at("chi_max_outlet_elevation_metres");
            }
            const lsdtt::hillslope_tag_basins basins =
                lsdtt::select_hillslope_tag_basins(flow, routing.get_data(), basin_options);
            std::cout << "Basins: " << basins.n_candidates << " candidates, " << basins.outlets.size()
                      << " kept; left out: " << basins.n_excluded_elevation << " outlet too high, "
                      << basins.n_excluded_grid_edge << " on the grid edge, " << basins.n_excluded_too_small
                      << " below the channel threshold, " << basins.n_excluded_nested << " nested\n";
            lsdtt::lsdtt_run_record::add_provenance("n_chi_basins", std::to_string(basins.outlets.size()));

            lsdtt::hillslope_chi_tag_options chi_options;
            chi_options.chi.m_over_n = ff.at("m_over_n");
            chi_options.chi.A_0 = ff.at("chi_A0");
            chi_options.critical_divergence_m = ff.at("fit_critical_divergence_metres");
            chi_options.minimum_contributing_pixels = threshold;
            const lsdtt::hillslope_chi_tag_result tags =
                lsdtt::compute_hillslope_chi_tags(flow, routing, relief.channel_flat, basins, chi_options);
            std::cout << "Fitted " << tags.n_fitted_nodes << " channel nodes, " << tags.n_segments << " segments\n";
            if (bf.at("tag_chi"))  { report_tag("chi", tags.chi);   write_raster(tags.chi.tag, "_tag_chi"); }
            if (bf.at("tag_k_sn")) { report_tag("k_sn", tags.k_sn); write_raster(tags.k_sn.tag, "_tag_k_sn"); }
        }

        // ---- imported attributes -----------------------------------------------
        if (!attributes_file.empty())
        {
            std::cout << "\n--- Imported channel attributes ---\n";
            const std::string path = pp.add_read_path_if_needed(attributes_file);
            const lsdtt::channel_attribute_table table = lsdtt::read_channel_attribute_table(
                path, meta, attribute_columns, ff.at("channel_attribute_nodata_value"));
            const lsdtt::channel_attribute_join_result joined =
                lsdtt::join_channel_attributes(flow, routing.get_mask(), channel, table, max_snap, max_gap);
            std::cout << "Records: " << joined.n_records << " read, " << joined.n_records_matched
                      << " matched, " << joined.n_records_dropped << " dropped, "
                      << joined.n_records_snapped_beyond_one_cell << " moved more than one cell";
            if (joined.n_exact_samples_off_channel > 0)
                std::cout << ", " << joined.n_exact_samples_off_channel << " exact samples off the channel";
            std::cout << "\n";
            lsdtt::lsdtt_run_record::add_provenance("n_attribute_records", std::to_string(joined.n_records));
            lsdtt::lsdtt_run_record::add_provenance("n_attribute_records_dropped", std::to_string(joined.n_records_dropped));
            for (const auto& column : joined.columns)
            {
                std::cout << "  " << column.name << ": of " << joined.n_channel_pixels << " channel pixels, "
                          << column.n_direct << " have a record (" << column.n_multiple_records
                          << " with several), " << column.n_filled << " gap-filled, " << column.n_missing
                          << " empty\n";
                const lsdtt::hillslope_tag_result tag = lsdtt::tag_hillslopes_from_channel(
                    meta, relief.channel_flat, routing.get_mask(), column.values, column.mask);
                report_tag(column.name, tag);
                write_raster(tag.tag, "_tag_" + column.name);
            }
        }

        pp.print_help(help_map, pp.get_write_fname() + "_README", k_status, k_citation);
        pp.print_parameters();
        std::cout << "\n" << k_tool_name << ": done.\n";
        return run_record.succeed();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }
}
