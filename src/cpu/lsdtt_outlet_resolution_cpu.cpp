// lsdtt_outlet_resolution_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_outlet_resolution_cpu.hpp"

#include <algorithm>
#include <stdexcept>

#include "lsdtt_flowinfo_cpu.hpp"

namespace lsdtt
{

std::vector<std::int64_t> resolve_outlet_flats_auto(
    const lsdtt_flow_base&           flow,
    const outlet_resolution_options& options,
    outlet_resolution_report&        report
)
{
    report = outlet_resolution_report{};
    report.from_file = false;

    std::vector<std::int64_t> all_outlets =
        get_all_basin_outlet_flats(flow, options.source_pixels_threshold);

    report.n_detected = all_outlets.size();

    if (all_outlets.empty())
    {
        throw std::runtime_error("No basin outlets found in the flow object.");
    }

    std::vector<std::int64_t> resolved;
    resolved.reserve(all_outlets.size());

    for (const std::int64_t o : all_outlets)
    {
        std::int64_t outlet = o;

        // An INTERIOR outlet drains to nodata rather than off the raster edge --
        // an island coast, for instance -- so it sits on the nodata fringe where
        // chi and drainage area are unreliable. Walking up the main stem lifts
        // it clear. Edge outlets are already fine and are left alone.
        if (options.n_pixels_move_outlet_upstream > 0 &&
            !is_flat_on_dem_edge(flow, o))
        {
            outlet = move_outlet_upstream(
                flow, o, options.n_pixels_move_outlet_upstream);
            if (outlet != o) ++report.n_moved_upstream;
        }

        // Re-check: moving upstream reduces contributing area, so an outlet that
        // qualified at its original position may not qualify here.
        if (flow.s_vector_size_unchecked(outlet) <
            static_cast<std::uint64_t>(options.source_pixels_threshold))
        {
            ++report.n_below_threshold;
            continue;
        }

        resolved.push_back(outlet);
    }

    // Two outlets can walk onto the same pixel, so de-duplicate.
    const std::size_t before = resolved.size();
    std::sort(resolved.begin(), resolved.end());
    resolved.erase(std::unique(resolved.begin(), resolved.end()), resolved.end());
    report.n_duplicates = before - resolved.size();

    if (resolved.empty())
    {
        throw std::runtime_error(
            "No basin outlets remain after the upstream move / threshold filter.");
    }

    return resolved;
}


std::vector<std::int64_t> resolve_outlet_flats_from_file(
    const lsdtt_raster_cpu&          dem,
    const lsdtt_flow_base&           flow,
    const std::string&               outlet_file_path,
    const outlet_resolution_options& options,
    outlet_resolution_report&        report
)
{
    report = outlet_resolution_report{};
    report.from_file = true;

    std::vector<lsdtt_input_point> requested =
        load_points_in_raster_crs(outlet_file_path, dem.get_meta());

    if (requested.empty())
    {
        throw std::runtime_error("Outlet file contained no usable points.");
    }

    for (const lsdtt_input_point& pt : requested)
    {
        switch (pt.status)
        {
            case lsdtt_point_input_status::coordinate_transform_failed:
                ++report.input_summary.transform_failed;   break;
            case lsdtt_point_input_status::unsupported_geometry_type:
                ++report.input_summary.unsupported_geom;   break;
            case lsdtt_point_input_status::outside_raster_extent:
                ++report.input_summary.outside_raster;     break;
            case lsdtt_point_input_status::empty_geometry:
                ++report.input_summary.empty_geometry;     break;
            case lsdtt_point_input_status::missing_source_crs:
                ++report.input_summary.missing_source_crs; break;
            default: break;
        }
    }

    std::vector<std::int64_t> outlet_flats;
    outlet_flats.reserve(requested.size());

    for (const lsdtt_input_point& pt : requested)
    {
        if (pt.status != lsdtt_point_input_status::success)
        {
            report.unsnappable_point_ids.push_back(pt.point_id);
            continue;
        }

        std::int64_t flat = pt.flat;

        if (options.snapping_method == "expanding_square" ||
            options.snapping_method == "downstream")
        {
            const lsdtt_snap_result snapped = snap_pixel_by_method(
                flow, pt.row, pt.col, options.snapping_method,
                options.square_snap, options.downstream_snap);

            if (snapped.status != lsdtt_snap_status::success)
            {
                report.unsnappable_point_ids.push_back(pt.point_id);
                continue;
            }
            flat = snapped.snapped_flat;
        }
        else if (!flow.is_valid_cell(flat))
        {
            report.unsnappable_point_ids.push_back(pt.point_id);
            continue;
        }

        outlet_flats.push_back(flat);
    }

    if (outlet_flats.empty())
    {
        throw std::runtime_error(
            "No valid outlet flats were produced from the outlet file.");
    }

    return outlet_flats;
}

} // namespace lsdtt
