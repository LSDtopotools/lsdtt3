// lsdtt_basin_stats_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include "lsdtt_basin_stats_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <exception>

#include "lsdtt_basin_cpu.hpp"
#include "lsdtt_channel_network_analysis_cpu.hpp"

namespace lsdtt
{

std::vector<basin_stats_outlet_result> compute_basin_stats(
    const lsdtt_flow_base&                flow,
    const lsdtt_raster_cpu&               dem,
    const lsdtt_raster_cpu&               gradient,
    const lsdtt_raster_cpu&               flow_distance,
    const lsdtt_junction_network_base*    channel_network,
    const lsdtt_coordinate_transform&     to_wgs84,
    const std::vector<lsdtt_input_point>& outlets,
    const basin_stats_options&            options
)
{
    const lsdtt_raster_meta& meta = dem.get_meta();

    const lsdtt_array2d_f&  dem_data  = dem.get_data_array();
    const lsdtt_array2d_u8& dem_mask  = dem.get_mask_array();
    const lsdtt_array2d_f&  grad_data = gradient.get_data_array();
    const lsdtt_array2d_u8& grad_mask = gradient.get_mask_array();
    const lsdtt_array2d_f&  fd_data   = flow_distance.get_data_array();
    const lsdtt_array2d_u8& fd_mask   = flow_distance.get_mask_array();

    const double cell_size    = meta.cell_size;
    const double cell_area_km2 = cell_size * cell_size / 1.0e6;
    const double diag_len      = cell_size * std::sqrt(2.0);

    const bool has_channel_network = (channel_network != nullptr);

    std::vector<basin_stats_outlet_result> results;
    results.reserve(outlets.size());

    for (const lsdtt_input_point& outlet : outlets)
    {
        basin_stats_outlet_result out;
        out.point_id = outlet.point_id;

        try
        {
            if (outlet.status != lsdtt_point_input_status::success)
            {
                out.outcome = basin_stats_outcome::point_not_mapped;
                results.push_back(std::move(out));
                continue;
            }

            int          row  = outlet.row;
            int          col  = outlet.col;
            std::int64_t flat = outlet.flat;
            lsdtt_xy     xy{outlet.x, outlet.y};

            if (options.snapping_method == "expanding_square" ||
                options.snapping_method == "downstream")
            {
                const lsdtt_snap_result snapped =
                    (options.snapping_method == "expanding_square")
                        ? snap_pixel_to_contributing_pixels_threshold_expanding_square(
                              flow, row, col, options.square_snap)
                        : snap_pixel_to_contributing_pixels_threshold_downstream(
                              flow, row, col, options.downstream_snap);

                out.snap_ran = true;
                out.snap     = snapped;

                if (snapped.status != lsdtt_snap_status::success)
                {
                    out.outcome = basin_stats_outcome::snap_failed;
                    results.push_back(std::move(out));
                    continue;
                }

                row  = snapped.snapped_row;
                col  = snapped.snapped_col;
                flat = snapped.snapped_flat;
                xy   = cell_center_from_row_col(meta, row, col);
            }
            else if (!flow.is_valid_cell(flat))
            {
                out.outcome = basin_stats_outcome::invalid_flow_pixel;
                results.push_back(std::move(out));
                continue;
            }

            // ONE basin object: it owns the upslope slice, and it answers the
            // completeness question. Resolving the slice separately and then
            // constructing a basin for the shape metrics did both twice.
            const lsdtt_basin_cpu basin(flat, flow);
            const lsdtt_flat_index_slice& upslope = basin.get_slice();

            // Boundary-condition aware: see the note in the header. A periodic
            // boundary is not truncation; an open one is.
            out.nodata_on_edge = !basin.check_basin_completeness();
            if (out.nodata_on_edge && options.reject_nodata_edge_basins)
            {
                out.outcome = basin_stats_outcome::rejected_nodata_edge;
                results.push_back(std::move(out));
                continue;
            }

            // Flow distance at the outlet, so the per-pixel values below can be
            // made relative to it. Absolute flow distance is measured from
            // whatever baselevel the routing found, not from this outlet.
            double outlet_fd = 0.0;
            if (fd_mask(row, col) != 0)
            {
                outlet_fd = fd_data(row, col);
            }

            std::vector<float> elev_values;
            std::vector<float> slope_values;
            std::vector<float> fdist_values;
            double sum_x = 0.0;
            double sum_y = 0.0;

            std::map<int, double> chan_len;
            for (int o = 1; o <= 6; ++o) { chan_len[o] = 0.0; }

            elev_values.reserve(upslope.size());
            slope_values.reserve(upslope.size());
            fdist_values.reserve(upslope.size());

            // ONE pass, deliberately. lsdtt_basin_cpu::get_basin_values(raster)
            // would give the elevation vector on its own, but this loop also
            // accumulates slope, flow distance relative to THIS outlet, the
            // centroid sums and channel length by Strahler order. Three or four
            // get_basin_values calls would be three or four passes over the
            // basin for values that are all available here at once.
            for (const std::int64_t* it = upslope.begin();
                 it != upslope.end(); ++it)
            {
                const std::int64_t basin_flat = *it;
                if (basin_flat < 0) continue;

                const int basin_row = flow.row_from_flat_unchecked(basin_flat);
                const int basin_col = flow.col_from_flat_unchecked(basin_flat);

                if (!dem_mask.in_bounds(basin_row, basin_col)) continue;
                if (dem_mask(basin_row, basin_col) == 0)       continue;

                elev_values.push_back(dem_data(basin_row, basin_col));

                if (grad_mask.in_bounds(basin_row, basin_col) &&
                    grad_mask(basin_row, basin_col) != 0)
                {
                    slope_values.push_back(grad_data(basin_row, basin_col));
                }

                if (fd_mask.in_bounds(basin_row, basin_col) &&
                    fd_mask(basin_row, basin_col) != 0)
                {
                    const double rel_fd =
                        std::max(0.0,
                                 static_cast<double>(fd_data(basin_row, basin_col))
                                 - outlet_fd);
                    fdist_values.push_back(static_cast<float>(rel_fd));
                }

                const lsdtt_xy pxy =
                    cell_center_from_row_col(meta, basin_row, basin_col);
                sum_x += pxy.x;
                sum_y += pxy.y;

                if (has_channel_network &&
                    channel_network->is_channel(basin_flat))
                {
                    const std::int32_t junc =
                        channel_network->junction_upstream_of_flat(basin_flat);
                    if (junc >= 0)
                    {
                        const int order = static_cast<int>(
                            channel_network->stream_order_of_junction(junc));
                        if (order >= 1 && order <= 6)
                        {
                            const std::uint8_t dir = flow.receiver_dir_of(basin_flat);
                            const bool diagonal =
                                (dir == 1 || dir == 3 || dir == 5 || dir == 7);
                            chan_len[order] += diagonal ? diag_len : cell_size;
                        }
                    }
                }
            }

            if (elev_values.empty())
            {
                out.outcome = basin_stats_outcome::no_valid_pixels;
                results.push_back(std::move(out));
                continue;
            }

            const lsdtt_descriptive_stats elev_stats =
                compute_descriptive_stats(elev_values);

            double s25 = 0.0, s50 = 0.0, s75 = 0.0;
            if (!slope_values.empty())
            {
                s25 = compute_percentile(slope_values, 25.0);
                s50 = compute_percentile(slope_values, 50.0);
                s75 = compute_percentile(slope_values, 75.0);
            }

            double fd25 = 0.0, fd50 = 0.0, fd75 = 0.0;
            if (!fdist_values.empty())
            {
                fd25 = compute_percentile(fdist_values, 25.0);
                fd50 = compute_percentile(fdist_values, 50.0);
                fd75 = compute_percentile(fdist_values, 75.0);
            }

            const std::size_t n_pix = elev_values.size();
            const double cx = sum_x / static_cast<double>(n_pix);
            const double cy = sum_y / static_cast<double>(n_pix);

            double cent_lat = 0.0, cent_lon = 0.0;
            if (to_wgs84.is_initialised())
            {
                const lsdtt_lonlat ll = xy_to_lonlat(cx, cy, to_wgs84);
                cent_lat = ll.latitude;
                cent_lon = ll.longitude;
            }

            // Basin planform shape metrics. Failure here is reported and the
            // record is still returned, as it was before: a basin whose outline
            // could not be traced still has valid elevation and slope figures.
            lsdtt_basin_shape shape;
            bool shape_ok = false;
            if (options.calculate_shape_metrics)
            {
                try
                {
                    lsdtt_basin_outline_cpu outline = basin.get_basin_outline();

                    const std::vector<lsdtt_pixel_rc>& div_rc =
                        outline.outline_divide_rc();
                    const lsdtt_raster_meta& dgm = outline.divide_grid_meta();

                    std::vector<lsdtt_xy> ring;
                    ring.reserve(div_rc.size());
                    for (const lsdtt_pixel_rc& rc : div_rc)
                    {
                        ring.push_back(cell_center_from_row_col(dgm, rc.row, rc.col));
                    }

                    if (ring.size() >= 3)
                    {
                        double basin_len = 0.0;
                        if (options.basin_length_definition == "geometric")
                        {
                            const std::vector<lsdtt_xy> hull = convex_hull(ring);
                            basin_len = feret_diameters(hull).d_max;
                        }
                        else  // flow_path
                        {
                            for (const float f : fdist_values)
                            {
                                basin_len = std::max(basin_len,
                                                     static_cast<double>(f));
                            }
                        }

                        const double simp_tol =
                            options.shape_perimeter_simplify_cells * cell_size;

                        shape = compute_basin_shape(
                            ring, basin_len,
                            elev_stats.min, elev_stats.mean, elev_stats.max,
                            simp_tol);
                        shape_ok = true;
                    }
                }
                catch (const std::exception& e)
                {
                    out.shape_failed = true;
                    out.message      = e.what();
                }
            }

            basin_stats_record& rec = out.record;
            rec.basin_id           = outlet.point_id;
            rec.latitude           = outlet.original_y;
            rec.longitude          = outlet.original_x;
            rec.outlet_x           = xy.x;
            rec.outlet_y           = xy.y;
            rec.outlet_row         = row;
            rec.outlet_col         = col;
            rec.outlet_flat        = flat;
            rec.n_pixels           = n_pix;
            rec.drainage_area_km2  = static_cast<double>(n_pix) * cell_area_km2;
            rec.nodata_on_edge     = out.nodata_on_edge;
            rec.centroid_x         = cx;
            rec.centroid_y         = cy;
            rec.centroid_latitude  = cent_lat;
            rec.centroid_longitude = cent_lon;

            if (options.calculate_centre_of_mass_bearing)
            {
                rec.com_bearing_deg = bearing_from_north_degrees(xy.x, xy.y, cx, cy);
                rec.com_distance_m  = std::hypot(cx - xy.x, cy - xy.y);
            }

            rec.stats              = elev_stats;
            rec.slope_25th         = s25;
            rec.slope_50th         = s50;
            rec.slope_75th         = s75;
            rec.flow_distance_25th = fd25;
            rec.flow_distance_50th = fd50;
            rec.flow_distance_75th = fd75;

            rec.shape_valid           = shape_ok;
            rec.basin_length_m        = shape.basin_length;
            rec.elongation_ratio      = shape.elongation_ratio;
            rec.form_factor           = shape.form_factor;
            rec.circularity_ratio     = shape.circularity_ratio;
            rec.gravelius_compactness = shape.gravelius_compactness;
            rec.solidity              = shape.solidity;
            rec.eccentricity          = shape.eccentricity;
            rec.orientation_deg       = shape.orientation_deg;
            rec.hypsometric_integral  = shape.hypsometric_integral;

            if (has_channel_network)
            {
                rec.channel_length_by_order = chan_len;
            }

            if (options.collect_elevation_values)
            {
                out.elevation_values = std::move(elev_values);
            }

            out.outcome = basin_stats_outcome::success;
            results.push_back(std::move(out));
        }
        catch (const std::exception& e)
        {
            out.outcome = basin_stats_outcome::failed;
            out.message = e.what();
            results.push_back(std::move(out));
        }
    }

    return results;
}

} // namespace lsdtt
