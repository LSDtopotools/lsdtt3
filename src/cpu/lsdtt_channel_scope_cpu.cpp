// lsdtt_channel_scope_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_channel_scope_cpu
//
// Implementation of channel scope resolution. See header for the mode and
// semantics documentation.
// -----------------------------------------------------------------------------

#include "lsdtt_channel_scope_cpu.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

void validate_request(
    const lsdtt_flow_base&       flow,
    const channel_scope_request& request)
{
    if (flow.empty())
    {
        throw std::invalid_argument(
            "resolve_channel_scope: flow object is empty.");
    }
    if (!flow.has_s_vector())
    {
        throw std::invalid_argument(
            "resolve_channel_scope: flow object was built without an "
            "s_vector; contributing-pixel queries are unavailable.");
    }
    if (request.source_flats.empty() && !request.outlet_flats.empty())
    {
        throw std::invalid_argument(
            "resolve_channel_scope: outlet_flats supplied without "
            "source_flats. Outlets are only meaningful for user-source "
            "modes.");
    }
    if (request.auto_detect_roles && !request.outlet_flats.empty())
    {
        throw std::invalid_argument(
            "resolve_channel_scope: auto_detect_roles cannot be combined with "
            "explicit outlet_flats. Supply all points in source_flats, or give "
            "explicit source/outlet lists with auto_detect_roles disabled.");
    }
    if (request.threshold_contributing_pixels < 1)
    {
        throw std::invalid_argument(
            "resolve_channel_scope: threshold_contributing_pixels must be "
            ">= 1, got " +
            std::to_string(request.threshold_contributing_pixels) + ".");
    }
}

// Trace the D8 receiver chain from source_flat. Appends visited flats to
// node_list. When stop_at_claimed is true, tracing stops *after* appending
// the first pixel already claimed in source_id (so the traced link joins the
// existing channel). Otherwise the trace runs to the DEM outlet.
void trace_downstream(
    const lsdtt_flow_base&           flow,
    std::int64_t                     source_flat,
    bool                             stop_at_claimed,
    const std::vector<std::int32_t>& source_id,
    std::vector<std::int64_t>&       node_list)
{
    std::int64_t current = source_flat;
    while (true)
    {
        node_list.push_back(current);

        if (stop_at_claimed && source_id[current] >= 0) return;
        if (flow.is_outlet_unchecked(current)) return;

        const std::int64_t next = flow.receiver_of_unchecked(current);
        if (next == current || !flow.is_flat_in_range(next)) return;
        current = next;
    }
}

// Classify snapped points into sources and outlets by downstream membership.
// A point is an outlet (returns 1) when its snapped flat lies STRICTLY
// downstream of some other supplied point - i.e. another point's D8 receiver
// chain passes through it. Points downstream of no other supplied point are
// sources (returns 0). Entries with a failed snap (flat < 0) return 0 and are
// handled by the caller. Coincident points (same snapped cell) do not mark each
// other because the walk starts at each point's receiver, excluding its own
// cell. A DAG always has at least one source (the topologically highest point).
std::vector<char> classify_downstream_outlets(
    const lsdtt_flow_base&           flow,
    const std::vector<std::int64_t>& snapped)
{
    const std::size_t n = snapped.size();
    std::vector<char> is_outlet(n, 0);

    std::unordered_map<std::int64_t, std::vector<std::size_t>> at_flat;
    for (std::size_t i = 0; i < n; ++i)
        if (snapped[i] >= 0) at_flat[snapped[i]].push_back(i);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (snapped[i] < 0) continue;

        std::int64_t cur = snapped[i];
        if (flow.is_outlet_unchecked(cur)) continue;

        // Walk strictly downstream (start at the receiver, exclude cur itself).
        cur = flow.receiver_of_unchecked(cur);
        while (true)
        {
            auto it = at_flat.find(cur);
            if (it != at_flat.end())
                for (const std::size_t j : it->second)
                    if (j != i) is_outlet[j] = 1;

            if (flow.is_outlet_unchecked(cur)) break;
            const std::int64_t nxt = flow.receiver_of_unchecked(cur);
            if (nxt == cur || !flow.is_flat_in_range(nxt)) break;
            cur = nxt;
        }
    }
    return is_outlet;
}

} // namespace


channel_scope_result resolve_channel_scope(
    const lsdtt_flow_base&       flow,
    const channel_scope_request& request)
{
    validate_request(flow, request);

    channel_scope_result result;

    const std::size_t grid_size = flow.full_grid_size();
    result.channel_mask.assign(grid_size, 0);
    result.channel_source_id.assign(grid_size, -1);

    const bool mode_1 = request.source_flats.empty();

    // Outlet flats actually used for truncation: the explicit request outlets,
    // or (in auto-role mode) the points detected as lying downstream of others.
    std::vector<std::int64_t> eff_outlets;

    // -------------------------------------------------------------------------
    // Resolve sources
    // -------------------------------------------------------------------------
    if (mode_1)
    {
        result.resolved_sources = flow.get_sources_index_threshold(
            request.threshold_contributing_pixels);
    }
    else if (request.auto_detect_roles)
    {
        // Snap every supplied point, then classify by downstream membership.
        lsdtt_square_snap_options snap_opts;
        snap_opts.contributing_pixels_threshold =
            static_cast<std::uint64_t>(request.threshold_contributing_pixels);
        snap_opts.max_search_radius_pixels =
            request.source_snap_max_radius_pixels;

        result.source_snaps =
            snap_flats_to_contributing_pixels_threshold_expanding_square(
                flow, request.source_flats, snap_opts);

        std::vector<std::int64_t> snapped(result.source_snaps.size(), -1);
        for (std::size_t i = 0; i < result.source_snaps.size(); ++i)
            if (result.source_snaps[i].status == lsdtt_snap_status::success)
                snapped[i] = result.source_snaps[i].snapped_flat;

        const std::vector<char> is_outlet =
            classify_downstream_outlets(flow, snapped);

        result.input_point_role.assign(snapped.size(), -1);
        for (std::size_t i = 0; i < snapped.size(); ++i)
        {
            if (snapped[i] < 0) continue;               // -1: snap failed
            if (is_outlet[i])
            {
                result.input_point_role[i] = 1;         // outlet
                eff_outlets.push_back(snapped[i]);
            }
            else
            {
                result.input_point_role[i] = 0;         // source
                result.resolved_sources.push_back(snapped[i]);
            }
        }
    }
    else if (request.snap_sources)
    {
        lsdtt_square_snap_options snap_opts;
        snap_opts.contributing_pixels_threshold =
            static_cast<std::uint64_t>(request.threshold_contributing_pixels);
        snap_opts.max_search_radius_pixels =
            request.source_snap_max_radius_pixels;

        result.source_snaps =
            snap_flats_to_contributing_pixels_threshold_expanding_square(
                flow, request.source_flats, snap_opts);

        for (const lsdtt_snap_result& snap : result.source_snaps)
        {
            if (snap.status == lsdtt_snap_status::success)
                result.resolved_sources.push_back(snap.snapped_flat);
        }
        eff_outlets = request.outlet_flats;
    }
    else
    {
        for (const std::int64_t flat : request.source_flats)
        {
            if (flow.is_flat_in_range(flat) && flow.is_valid_cell(flat))
                result.resolved_sources.push_back(flat);
        }
        eff_outlets = request.outlet_flats;
    }

    const bool bounded = !eff_outlets.empty();
    result.mode = mode_1 ? channel_scope_mode::all_channels_threshold
                : bounded ? channel_scope_mode::user_sources_bounded
                          : channel_scope_mode::user_sources;

    const std::size_t n_sources = result.resolved_sources.size();
    result.resolved_outlets.assign(n_sources, -1);
    result.channel_node_lists.resize(n_sources);

    // -------------------------------------------------------------------------
    // Trace channels downstream
    //
    // Mode 1 stops each trace where it meets an already-claimed pixel; claims
    // are made during tracing so total work is O(n_pixels). User-source modes
    // trace full profiles; claiming happens after outlet truncation below.
    // -------------------------------------------------------------------------
    const bool user_modes = !mode_1;

    for (std::size_t s = 0; s < n_sources; ++s)
    {
        trace_downstream(flow, result.resolved_sources[s],
                         /*stop_at_claimed=*/mode_1,
                         result.channel_source_id,
                         result.channel_node_lists[s]);

        if (mode_1)
        {
            for (const std::int64_t flat : result.channel_node_lists[s])
            {
                if (result.channel_source_id[flat] < 0)
                {
                    result.channel_source_id[flat] =
                        static_cast<std::int32_t>(s);
                    result.channel_mask[flat] = 1;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Match outlets onto traced paths and truncate (bounded mode)
    //
    // Each outlet snaps to the nearest node across all traced paths. When a
    // path is matched by several outlets, the most upstream truncation wins.
    // -------------------------------------------------------------------------
    if (user_modes && bounded)
    {
        const double cell = flow.cell_size();
        const std::size_t n_outlets = eff_outlets.size();
        result.outlet_matched_source.assign(n_outlets, -1);
        result.outlet_snap_distance_metres.assign(
            n_outlets, std::numeric_limits<double>::quiet_NaN());

        // Truncation index per source (one past the outlet node); default =
        // full path.
        std::vector<std::size_t> truncate_at(n_sources);
        for (std::size_t s = 0; s < n_sources; ++s)
            truncate_at[s] = result.channel_node_lists[s].size();

        for (std::size_t o = 0; o < n_outlets; ++o)
        {
            const std::int64_t raw = eff_outlets[o];
            if (!flow.is_flat_in_range(raw)) continue;

            const double raw_row = flow.row_from_flat_unchecked(raw);
            const double raw_col = flow.col_from_flat_unchecked(raw);

            double      best_d2   = std::numeric_limits<double>::max();
            std::size_t best_s    = 0;
            std::size_t best_node = 0;
            bool        found     = false;

            for (std::size_t s = 0; s < n_sources; ++s)
            {
                const std::vector<std::int64_t>& nodes =
                    result.channel_node_lists[s];
                for (std::size_t n = 0; n < nodes.size(); ++n)
                {
                    const double dr =
                        flow.row_from_flat_unchecked(nodes[n]) - raw_row;
                    const double dc =
                        flow.col_from_flat_unchecked(nodes[n]) - raw_col;
                    const double d2 = dr * dr + dc * dc;
                    if (d2 < best_d2)
                    {
                        best_d2   = d2;
                        best_s    = s;
                        best_node = n;
                        found     = true;
                    }
                }
            }

            if (!found) continue;

            const double dist_m = std::sqrt(best_d2) * cell;
            if (request.outlet_snap_max_distance_metres > 0.0 &&
                dist_m > request.outlet_snap_max_distance_metres)
            {
                continue;
            }

            result.outlet_matched_source[o] =
                static_cast<std::int32_t>(best_s);
            result.outlet_snap_distance_metres[o] = dist_m;

            // Truncate every channel that passes through the matched flat (so an
            // outlet below a confluence cuts all channels merging above it).
            // Most upstream truncation wins.
            const std::int64_t match_flat =
                result.channel_node_lists[best_s][best_node];
            for (std::size_t s = 0; s < n_sources; ++s)
            {
                const std::vector<std::int64_t>& nodes =
                    result.channel_node_lists[s];
                for (std::size_t n = 0; n < nodes.size(); ++n)
                {
                    if (nodes[n] == match_flat)
                    {
                        if (n + 1 < truncate_at[s]) truncate_at[s] = n + 1;
                        break;
                    }
                }
            }
        }

        for (std::size_t s = 0; s < n_sources; ++s)
        {
            std::vector<std::int64_t>& nodes = result.channel_node_lists[s];
            if (truncate_at[s] < nodes.size())
            {
                nodes.resize(truncate_at[s]);
                result.resolved_outlets[s] = nodes.back();
            }
        }
    }

    // -------------------------------------------------------------------------
    // Claim mask/source-id pixels for user-source modes (post-truncation)
    // -------------------------------------------------------------------------
    if (user_modes)
    {
        for (std::size_t s = 0; s < n_sources; ++s)
        {
            for (const std::int64_t flat : result.channel_node_lists[s])
            {
                if (result.channel_source_id[flat] < 0)
                {
                    result.channel_source_id[flat] =
                        static_cast<std::int32_t>(s);
                    result.channel_mask[flat] = 1;
                }
            }
        }
    }

    return result;
}

// =============================================================================
// Shared CLI-side plumbing. See the header for why it is here and not in a CLI.
// =============================================================================

lsdtt_point_input_format resolve_point_input_format(
    const std::string& format_name)
{
    if (format_name == "auto_detect" || format_name == "auto")
        return lsdtt_point_input_format::auto_detect;
    if (format_name == "csv" ||
        format_name == "csv_latitude_longitude_epsg4326")
        return lsdtt_point_input_format::csv_latitude_longitude_epsg4326;
    if (format_name == "ogr_vector_layer" || format_name == "vector" ||
        format_name == "ogr")
        return lsdtt_point_input_format::ogr_vector_layer;

    throw std::runtime_error("Unknown point input format: " + format_name);
}

std::vector<std::int64_t> load_role_point_flats(
    const point_role_input&  role,
    const lsdtt_raster_meta& meta)
{
    if (role.path.empty()) return {};

    lsdtt_point_input_options popts;
    popts.input_format             = resolve_point_input_format(role.input_format);
    popts.vector_layer_name        = role.layer_name;
    popts.csv_latitude_field_name  = role.latitude_field_name;
    popts.csv_longitude_field_name = role.longitude_field_name;
    // Keep out-of-raster points so the snap diagnostics can report them as
    // unmapped rather than silently dropping them at load.
    popts.keep_points_outside_raster = true;

    std::cout << "Loading " << role.role_name << " points from: "
              << role.path << "\n";
    const std::vector<lsdtt_input_point> pts =
        load_points_in_raster_crs(role.path, meta, popts);
    std::vector<std::int64_t> flats =
        successfully_mapped_input_points_to_flats(pts);
    std::cout << "  " << pts.size() << " point(s) loaded, "
              << flats.size() << " mapped to valid cells\n";
    return flats;
}

// -----------------------------------------------------------------------------
// write_channel_scope_diagnostics
// -----------------------------------------------------------------------------
void write_channel_scope_diagnostics(
    const channel_scope_result&             scope,
    const lsdtt_flow_base&                  flow,
    const lsdtt_raster_cpu&                 dem,
    lsdtt_data_loader&                      loader,
    const std::string&                      output_stem,
    const channel_scope_diagnostic_request& request)
{
    const lsdtt_raster_meta& meta = dem.get_meta();

    auto xy_of_rc = [&](int r, int c)
    { return cell_center_from_row_col(meta, r, c); };

    if (request.write_scope_raster)
    {
        try
        {
            // source_id + 1, so 0 stays background and the raster reads the
            // same way as every other labelled mask in the project.
            std::vector<float> scope_vals(scope.channel_source_id.size(), 0.0f);
            for (std::size_t i = 0; i < scope_vals.size(); ++i)
            {
                if (scope.channel_source_id[i] >= 0)
                    scope_vals[i] =
                        static_cast<float>(scope.channel_source_id[i] + 1);
            }
            const std::string path = output_stem + "_channel_scope.tif";
            std::cout << "Writing channel scope raster (UInt32): " << path << "\n";
            loader.write_full(scope_vals, meta, path, request.raster_format,
                              GDT_UInt32);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing channel scope raster failed: "
                      << e.what() << "\n";
        }
    }

    if (request.write_scope_nodes_csv)
    {
        try
        {
            const std::string path = output_stem + "_channel_scope_nodes.csv";
            std::cout << "Writing channel scope nodes: " << path << "\n";
            const std::vector<float>& dem_d = dem.get_data();
            std::vector<std::vector<std::string>> rows;
            for (std::size_t s = 0; s < scope.channel_node_lists.size(); ++s)
            {
                const std::vector<std::int64_t>& nodes =
                    scope.channel_node_lists[s];
                for (std::size_t n = 0; n < nodes.size(); ++n)
                {
                    const int r = flow.row_from_flat_unchecked(nodes[n]);
                    const int c = flow.col_from_flat_unchecked(nodes[n]);
                    const lsdtt_xy xy = xy_of_rc(r, c);
                    rows.push_back({
                        std::to_string(s),
                        std::to_string(n),
                        std::to_string(xy.x),
                        std::to_string(xy.y),
                        std::to_string(dem_d[static_cast<std::size_t>(nodes[n])]),
                        std::to_string(flow.s_vector_size(nodes[n])) });
                }
            }
            loader.write_csv(path,
                {"source_id","node_seq","x","y","elevation",
                 "contributing_pixels"},
                rows);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing channel scope nodes failed: "
                      << e.what() << "\n";
        }
    }

    if (!request.write_point_snap_csv) return;

    try
    {
        auto status_name = [](lsdtt_snap_status st)
        {
            switch (st)
            {
                case lsdtt_snap_status::success:
                    return "success";
                case lsdtt_snap_status::input_out_of_bounds:
                    return "input_out_of_bounds";
                case lsdtt_snap_status::input_invalid_flow_cell:
                    return "input_invalid_flow_cell";
                case lsdtt_snap_status::threshold_not_found_in_search_window:
                    return "threshold_not_found";
                default:
                    return "outlet_reached_before_threshold";
            }
        };

        struct scope_pt { double x, y; std::string role; int id; };
        std::vector<scope_pt> raw_pts, snapped_pts;
        std::vector<std::vector<std::string>> rows;

        // In auto-role mode the role comes from input_point_role: one supplied
        // point may be a source or an auto-detected outlet.
        const bool auto_roles = !scope.input_point_role.empty();

        for (std::size_t i = 0; i < scope.source_snaps.size(); ++i)
        {
            const lsdtt_snap_result& sr = scope.source_snaps[i];
            const std::string role =
                (auto_roles && i < scope.input_point_role.size() &&
                 scope.input_point_role[i] == 1) ? "outlet" : "source";
            const lsdtt_xy raw = xy_of_rc(sr.input_row, sr.input_col);
            const bool ok = sr.status == lsdtt_snap_status::success;
            lsdtt_xy snp{0.0, 0.0};
            double snap_dist = -1.0;
            if (ok)
            {
                snp = xy_of_rc(sr.snapped_row, sr.snapped_col);
                const double dx = snp.x - raw.x, dy = snp.y - raw.y;
                snap_dist = std::sqrt(dx * dx + dy * dy);
            }
            rows.push_back({
                role, std::to_string(i), status_name(sr.status),
                std::to_string(raw.x), std::to_string(raw.y),
                ok ? std::to_string(snp.x) : "",
                ok ? std::to_string(snp.y) : "",
                std::to_string(snap_dist),
                std::to_string(sr.snapped_contributing_pixels) });
            raw_pts.push_back({raw.x, raw.y, role, static_cast<int>(i)});
            if (ok)
                snapped_pts.push_back({snp.x, snp.y, role, static_cast<int>(i)});
        }

        // Explicit-outlet mode only: the separate outlet file's diagnostics.
        const std::size_t n_outlet_rows =
            auto_roles ? 0 : scope.outlet_matched_source.size();
        for (std::size_t o = 0; o < n_outlet_rows; ++o)
        {
            const std::int32_t ms = scope.outlet_matched_source[o];
            const bool ok = ms >= 0;
            std::string sx, sy;
            if (ok && static_cast<std::size_t>(ms) < scope.resolved_outlets.size()
                && scope.resolved_outlets[static_cast<std::size_t>(ms)] >= 0)
            {
                const std::int64_t of =
                    scope.resolved_outlets[static_cast<std::size_t>(ms)];
                const lsdtt_xy oxy = xy_of_rc(
                    flow.row_from_flat_unchecked(of),
                    flow.col_from_flat_unchecked(of));
                sx = std::to_string(oxy.x);
                sy = std::to_string(oxy.y);
                snapped_pts.push_back({oxy.x, oxy.y, "outlet",
                                       static_cast<int>(o)});
            }
            rows.push_back({
                "outlet", std::to_string(o), ok ? "matched" : "unmatched",
                "", "", sx, sy,
                std::to_string(scope.outlet_snap_distance_metres.empty()
                    ? -1.0 : scope.outlet_snap_distance_metres[o]),
                "" });
        }

        const std::string path = output_stem + "_point_snaps.csv";
        std::cout << "Writing point snap diagnostics: " << path << "\n";
        loader.write_csv(path,
            {"role","point_id","status","raw_x","raw_y",
             "snapped_x","snapped_y","snap_distance_m","contributing_pixels"},
            rows);

        // Raw and snapped point layers as well as the CSV, so both styles can
        // be overlaid in one QC figure.
        auto write_pt_layer = [&](const std::vector<scope_pt>& pts,
                                  const std::string& suffix,
                                  const std::string& layer)
        {
            if (pts.empty()) return;
            const std::string vp = output_stem + suffix;
            std::cout << "Writing " << pts.size() << " point(s): " << vp << "\n";
            lsdtt_vector_schema schema;
            schema.field_names = {"role", "point_id"};
            schema.field_types = {
                lsdtt_vector_field_type::string_type,
                lsdtt_vector_field_type::int64_type };
            loader.begin_vector_write(
                vp, layer, schema, meta.crs_wkt,
                vector_output_format::flatgeobuf, false);
            for (const scope_pt& pt : pts)
                loader.write_point(pt.x, pt.y,
                    { pt.role, static_cast<std::int64_t>(pt.id) });
            loader.finish_vector_write();
        };

        write_pt_layer(raw_pts, "_scope_points_raw.fgb", "scope_points_raw");
        write_pt_layer(snapped_pts, "_scope_points_snapped.fgb",
                       "scope_points_snapped");
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning: writing point snap diagnostics failed: "
                  << e.what() << "\n";
    }
}

// -----------------------------------------------------------------------------
// write_channel_network_layers
// -----------------------------------------------------------------------------
void write_channel_network_layers(
    const lsdtt_junction_network_base&   jnet,
    const std::string&                   output_stem,
    const channel_network_write_request& request)
{
    try
    {
        const std::string ext = output_extension(request.format);

        if (request.write_lines)
        {
            const std::string path = output_stem + "_channel_network_lines" + ext;
            std::cout << "Writing: " << path << "\n";
            jnet.write_channel_network_lines(
                path, request.format, "channel_network_lines", true,
                request.channel_source_id);
        }

        // Point layers have their own format enum; GPKG is the only case where
        // the two differ, everything else lands in FlatGeobuf.
        const point_output_format pt_fmt =
            (request.format == vector_output_format::gpkg)
            ? point_output_format::gpkg
            : point_output_format::flatgeobuf;

        if (request.write_channel_points)
        {
            const std::string path =
                output_stem + "_channel_network_points" + ext;
            std::cout << "Writing: " << path << "\n";
            jnet.write_channel_network_points(path, pt_fmt);
        }

        if (request.write_junction_points)
        {
            const std::string path =
                output_stem + "_junction_network_points" + ext;
            std::cout << "Writing: " << path << "\n";
            jnet.write_junction_network_points(path, pt_fmt);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning: writing channel network failed: "
                  << e.what() << "\n";
    }
}

} // namespace lsdtt
