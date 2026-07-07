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
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

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
// downstream of some other supplied point — i.e. another point's D8 receiver
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

} // namespace lsdtt
