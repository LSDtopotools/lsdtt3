// lsdtt_wetness_index_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_wetness_index_cpu.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

lsdtt_raster_cpu compute_twi_from_area_and_slope(
    const lsdtt_raster_cpu&  drainage_area,
    const lsdtt_raster_cpu&  slope,
    const lsdtt_raster_meta& meta)
{
    const int    n_rows     = meta.n_rows;
    const int    n_cols     = meta.n_cols;
    const double pixel_area = meta.cell_size * meta.cell_size;
    const float  nodata_val = static_cast<float>(meta.no_data_value);

    const lsdtt_array2d_f&  area_arr  = drainage_area.get_data_array();
    const lsdtt_array2d_u8& area_msk  = drainage_area.get_mask_array();
    const lsdtt_array2d_f&  slope_arr = slope.get_data_array();
    const lsdtt_array2d_u8& slope_msk = slope.get_mask_array();

    const lsdtt_raster_meta& area_meta  = drainage_area.get_meta();
    const lsdtt_raster_meta& slope_meta = slope.get_meta();

    if (area_meta.n_rows != n_rows || area_meta.n_cols != n_cols ||
        slope_meta.n_rows != n_rows || slope_meta.n_cols != n_cols)
    {
        throw std::runtime_error(
            "compute_twi_from_area_and_slope: raster dimensions disagree "
            "(area " + std::to_string(area_meta.n_rows) + "x" +
            std::to_string(area_meta.n_cols) + ", slope " +
            std::to_string(slope_meta.n_rows) + "x" +
            std::to_string(slope_meta.n_cols) + ", requested " +
            std::to_string(n_rows) + "x" + std::to_string(n_cols) + ")");
    }

    lsdtt_array2d_f  twi_data(n_rows, n_cols, nodata_val);
    lsdtt_array2d_u8 twi_mask(n_rows, n_cols, uint8_t(0));

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            if (area_msk(row, col) == 0) continue;
            if (slope_msk(row, col) == 0) continue;

            const float area_val  = area_arr(row, col);
            const float slope_val = slope_arr(row, col);

            if (area_val > 0.0f && slope_val > 0.0f)
            {
                // Float arithmetic throughout -- see the precision contract in
                // the header. Do not promote this to double.
                twi_data(row, col) = std::log(
                    (area_val / static_cast<float>(pixel_area)) / slope_val
                );
                twi_mask(row, col) = uint8_t(1);
            }
        }
    }

    return lsdtt_raster_cpu(std::move(twi_data), std::move(twi_mask), meta);
}

lsdtt_raster_cpu compute_twi(
    const lsdtt_raster_cpu&    routed_dem,
    const lsdtt_flowgraph_cpu& graph,
    const twi_options&         opts)
{
    poly_fit_coeffs coeffs;
    try
    {
        coeffs = direct_polyfit(
            routed_dem,
            opts.slope_window_radius_metres,
            opts.n_threads
        );
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error(
            "compute_twi: polyfit for slope failed: " + std::string(e.what()));
    }

    polyfit_metric_request slope_req;
    slope_req.metrics   = { polyfit_metric::slope };
    slope_req.n_threads = opts.n_threads;

    polyfit_metric_result slope_res;
    try
    {
        slope_res = compute_polyfit_metrics(coeffs, slope_req);
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error(
            "compute_twi: slope metric computation failed: " +
            std::string(e.what()));
    }

    lsdtt_raster_result area_result;
    try
    {
        area_result = graph.accumulate_raster();
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error(
            "compute_twi: drainage area accumulation failed: " +
            std::string(e.what()));
    }

    const lsdtt_raster_cpu area_raster(std::move(area_result));

    return compute_twi_from_area_and_slope(
        area_raster,
        slope_res.at(polyfit_metric::slope),
        routed_dem.get_meta()
    );
}

} // namespace lsdtt
