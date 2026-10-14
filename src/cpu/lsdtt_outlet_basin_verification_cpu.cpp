// lsdtt_outlet_basin_verification_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Boris Gailleton
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_outlet_basin_verification_cpu.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <memory>

#include <ogr_spatialref.h>

#include "lsdtt_flow_snapping_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {
namespace {

std::string lower_copy(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

} // namespace

lsdtt_outlet_basin_verification_result verify_outlet_group_basins(
    const lsdtt_raster_cpu& dem,
    const lsdtt_dem_group& group,
    const std::vector<lsdtt_outlet_bounds_record>& outlets,
    const lsdtt_outlet_basin_verification_options& options)
{
    lsdtt_outlet_basin_verification_result out;
    out.min_edge_margin_pixels = std::numeric_limits<double>::infinity();
    out.min_nodata_margin_pixels = std::numeric_limits<double>::infinity();
    if (options.required_edge_buffer_pixels < 0 ||
        options.required_nodata_buffer_pixels < 0 || options.tight_bbox_buffer_pixels < 0)
    {
        out.message = "basin verification buffers must be non-negative";
        return out;
    }

    const lsdtt_raster_meta& meta = dem.get_meta();
    const lsdtt_array2d_u8& mask = dem.get_mask_array();
    lsdtt_flowinfo_options flow_options;
    flow_options.boundary_conditions = "oooo";
    flow_options.interior_nodata_is_outlet = true;
    flow_options.allow_pits_as_outlets = true;

    lsdtt_flowinfo flow;
    try { flow.create(dem.get_data_array(), mask, meta, flow_options); }
    catch (const std::exception& e)
    {
        out.message = std::string("could not build flow routing for verification: ") + e.what();
        return out;
    }

    OGRSpatialReference wgs84;
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    wgs84.importFromEPSG(4326);
    OGRSpatialReference dem_srs;
    dem_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (meta.crs_wkt.empty() || dem_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE)
    {
        out.message = "could not parse DEM projection for verification";
        return out;
    }
    std::unique_ptr<OGRCoordinateTransformation, decltype(&OCTDestroyCoordinateTransformation)> to_dem(
        OGRCreateCoordinateTransformation(&wgs84, &dem_srs), &OCTDestroyCoordinateTransformation);
    std::unique_ptr<OGRCoordinateTransformation, decltype(&OCTDestroyCoordinateTransformation)> to_wgs84(
        OGRCreateCoordinateTransformation(&dem_srs, &wgs84), &OCTDestroyCoordinateTransformation);
    if (!to_dem || !to_wgs84)
    {
        out.message = "could not create DEM coordinate transformations for verification";
        return out;
    }

    const std::string snapping = lower_copy(options.snapping_method);
    if (snapping != "none" && snapping != "expanding_square" && snapping != "downstream")
    {
        out.message = "snapping_method must be none, expanding_square, or downstream";
        return out;
    }
    lsdtt_square_snap_options square_options;
    square_options.contributing_pixels_threshold = options.contributing_pixels_threshold;
    square_options.min_search_radius_pixels = options.min_search_radius_pixels;
    square_options.max_search_radius_pixels = options.max_search_radius_pixels;
    lsdtt_downstream_snap_options downstream_options;
    downstream_options.contributing_pixels_threshold = options.contributing_pixels_threshold;

    const int rows = meta.n_rows;
    const int cols = meta.n_cols;
    out.basin_labels.assign(static_cast<std::size_t>(rows) * cols, -1);
    int row_min = rows, row_max = -1, col_min = cols, col_max = -1;
    bool has_pixels = false;
    bool all_ok = true;
    std::string failure;
    for (const std::size_t index : group.record_indices)
    {
        if (index >= outlets.size()) { all_ok = false; failure = "group references an unknown outlet"; break; }
        double x = outlets[index].longitude;
        double y = outlets[index].latitude;
        if (!to_dem->Transform(1, &x, &y)) { all_ok = false; failure = "outlet transform failed on DEM"; break; }
        int col = -1, row = -1;
        try { col = col_from_x(meta, x); row = row_from_y(meta, y); }
        catch (const std::exception&) { all_ok = false; failure = "outlet lies outside DEM"; break; }
        if (row < 0 || row >= rows || col < 0 || col >= cols)
        {
            all_ok = false; failure = "outlet lies outside DEM"; break;
        }

        const lsdtt_snap_result result = snap_pixel_by_method(
            flow, row, col, snapping, square_options, downstream_options);
        if (result.status != lsdtt_snap_status::success)
        {
            all_ok = false;
            failure = snapping == "none"
                ? "outlet pixel is invalid on DEM"
                : "could not snap outlet to channel cell (" + snapping + ")";
            break;
        }
        const std::int64_t snapped = result.snapped_flat;

        out.snapped_outlets.push_back({outlets[index].point_id, snapped,
                                       flow.row_from_flat_unchecked(snapped),
                                       flow.col_from_flat_unchecked(snapped)});

        for (const std::int64_t flat : flow.upslope_nodes(snapped))
        {
            const int basin_row = flow.row_from_flat_unchecked(flat);
            const int basin_col = flow.col_from_flat_unchecked(flat);
            has_pixels = true;
            const std::size_t flat_index = static_cast<std::size_t>(flat);
            if (out.basin_labels[flat_index] < 0)
                out.basin_labels[flat_index] = outlets[index].point_id;
            row_min = std::min(row_min, basin_row); row_max = std::max(row_max, basin_row);
            col_min = std::min(col_min, basin_col); col_max = std::max(col_max, basin_col);
            const int edge = std::min(std::min(basin_row, rows - 1 - basin_row),
                                      std::min(basin_col, cols - 1 - basin_col));
            out.min_edge_margin_pixels = std::min(out.min_edge_margin_pixels, static_cast<double>(edge));
            if (edge < options.required_edge_buffer_pixels)
            {
                all_ok = false; failure = "basin pixel buffer to raster edge is too small"; break;
            }

            bool found_nodata = false;
            int nodata_margin = options.required_nodata_buffer_pixels + 1;
            for (int r = std::max(0, basin_row - options.required_nodata_buffer_pixels);
                 r <= std::min(rows - 1, basin_row + options.required_nodata_buffer_pixels) && !found_nodata; ++r)
                for (int c = std::max(0, basin_col - options.required_nodata_buffer_pixels);
                     c <= std::min(cols - 1, basin_col + options.required_nodata_buffer_pixels); ++c)
                    if (mask(r, c) == 0)
                    {
                        nodata_margin = std::max(std::abs(r - basin_row), std::abs(c - basin_col));
                        found_nodata = true;
                        break;
                    }
            out.min_nodata_margin_pixels = std::min(
                out.min_nodata_margin_pixels, static_cast<double>(nodata_margin));
            if (found_nodata && nodata_margin <= options.required_nodata_buffer_pixels)
            {
                all_ok = false; failure = "basin pixel buffer to nodata is too small"; break;
            }
        }
        if (!all_ok) break;
    }

    if (has_pixels)
    {
        const int r0 = std::max(0, row_min - options.tight_bbox_buffer_pixels);
        const int r1 = std::min(rows - 1, row_max + options.tight_bbox_buffer_pixels);
        const int c0 = std::max(0, col_min - options.tight_bbox_buffer_pixels);
        const int c1 = std::min(cols - 1, col_max + options.tight_bbox_buffer_pixels);
        double x[4] = {meta.x_min + c0 * meta.cell_size, meta.x_min + c0 * meta.cell_size,
                       meta.x_min + (c1 + 1) * meta.cell_size, meta.x_min + (c1 + 1) * meta.cell_size};
        double y[4] = {meta.y_max() - (r1 + 1) * meta.cell_size, meta.y_max() - r0 * meta.cell_size,
                       meta.y_max() - (r1 + 1) * meta.cell_size, meta.y_max() - r0 * meta.cell_size};
        if (to_wgs84->Transform(4, x, y))
        {
            out.tight_west = std::max(-180.0, *std::min_element(x, x + 4));
            out.tight_east = std::min(180.0, *std::max_element(x, x + 4));
            out.tight_south = std::max(-90.0, *std::min_element(y, y + 4));
            out.tight_north = std::min(90.0, *std::max_element(y, y + 4));
            out.tight_bbox_valid = out.tight_west < out.tight_east && out.tight_south < out.tight_north;
        }
    }
    if (!has_pixels && all_ok) { all_ok = false; failure = "no basin pixels found for outlet group"; }
    if (!all_ok)
    {
        out.snapped_outlets.clear();
        out.basin_labels.clear();
    }
    if (!std::isfinite(out.min_edge_margin_pixels)) out.min_edge_margin_pixels = -1.0;
    if (!std::isfinite(out.min_nodata_margin_pixels)) out.min_nodata_margin_pixels = -1.0;
    out.ok = all_ok;
    out.message = all_ok ? "verification passed" : failure;
    return out;
}

} // namespace lsdtt
