// lsdtt_channel_mask_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_channel_mask_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_vector_rasterise_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

namespace
{

void require_positive(double v, const char* name, const char* who)
{
    if (!(std::isfinite(v) && v > 0.0))
        throw std::invalid_argument(
            std::string(who) + ": " + name + " must be finite and > 0 (got " +
            std::to_string(v) + ")");
}

void require_same_grid(const lsdtt_raster_cpu& r, const lsdtt_raster_meta& meta,
                       const char* name, const char* who)
{
    const lsdtt_raster_meta& m = r.get_meta();
    if (m.n_rows != meta.n_rows || m.n_cols != meta.n_cols)
        throw std::invalid_argument(
            std::string(who) + ": " + name + " is " + std::to_string(m.n_rows) +
            "x" + std::to_string(m.n_cols) + " but the grid is " +
            std::to_string(meta.n_rows) + "x" + std::to_string(meta.n_cols));
}

} // namespace

double emergent_channel_area_m2(const channel_mask_options& opts)
{
    const char* who = "emergent_channel_area_m2";
    require_positive(opts.hydraulic_radius_coeff, "hydraulic_radius_coeff", who);
    require_positive(opts.hydraulic_radius_exponent, "hydraulic_radius_exponent", who);
    require_positive(opts.reference_area_m2, "reference_area_m2", who);
    require_positive(opts.sheet_flow_min_depth_metres,
                     "sheet_flow_min_depth_metres", who);

    return opts.reference_area_m2 *
           std::pow(opts.sheet_flow_min_depth_metres / opts.hydraulic_radius_coeff,
                    1.0 / opts.hydraulic_radius_exponent);
}

channel_mask_result compute_channel_mask(
    const lsdtt_raster_cpu&     drainage_area,
    const lsdtt_raster_cpu*     slope,
    const lsdtt_raster_cpu*     supplied_mask,
    const lsdtt_raster_meta&    meta,
    const channel_mask_options& opts)
{
    const char* who = "compute_channel_mask";

    require_same_grid(drainage_area, meta, "drainage_area", who);
    if (slope)        require_same_grid(*slope, meta, "slope", who);
    if (supplied_mask) require_same_grid(*supplied_mask, meta, "supplied_mask", who);

    switch (opts.mode)
    {
        case channel_mask_mode::none:
            // Validated inside emergent_channel_area_m2().
            break;
        case channel_mask_mode::area_threshold:
            require_positive(opts.channel_area_threshold_m2,
                             "channel_area_threshold_m2", who);
            break;
        case channel_mask_mode::area_slope_threshold:
            if (slope == nullptr)
                throw std::invalid_argument(
                    std::string(who) + ": slope is required in "
                    "area_slope_threshold mode");
            require_positive(opts.area_slope_threshold, "area_slope_threshold", who);
            require_positive(opts.min_slope, "min_slope", who);
            break;
        case channel_mask_mode::extraction:
            if (supplied_mask == nullptr)
                throw std::invalid_argument(
                    std::string(who) + ": supplied_mask is required in "
                    "extraction mode");
            break;
        case channel_mask_mode::vector:
            if (opts.channel_vector_fname.empty())
                throw std::invalid_argument(
                    std::string(who) + ": channel_vector_fname is required in "
                    "vector mode");
            if (opts.snap_search_radius_pixels < 0)
                throw std::invalid_argument(
                    std::string(who) + ": snap_search_radius_pixels must be >= 0");
            break;
    }

    const int   n_rows = meta.n_rows;
    const int   n_cols = meta.n_cols;
    const float ndv    = static_cast<float>(meta.no_data_value);

    const lsdtt_array2d_f&  a_d = drainage_area.get_data_array();
    const lsdtt_array2d_u8& a_m = drainage_area.get_mask_array();

    lsdtt_array2d_f  m_data(n_rows, n_cols, ndv);
    lsdtt_array2d_u8 m_mask(n_rows, n_cols, std::uint8_t(0));

    channel_mask_result out;

    const double emergent_area =
        (opts.mode == channel_mask_mode::none) ? emergent_channel_area_m2(opts) : 0.0;

    // --- vector mode: burn the lines first, then snap ---
    std::vector<std::uint8_t> burned;
    if (opts.mode == channel_mask_mode::vector)
    {
        lsdtt_raster_cpu lines;
        try
        {
            lines = rasterise_polylines(drainage_area, opts.channel_vector_fname);
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error(
                std::string(who) + ": could not rasterise '" +
                opts.channel_vector_fname + "': " + e.what());
        }
        require_same_grid(lines, meta, "rasterised channel vector", who);

        burned.assign(static_cast<std::size_t>(n_rows) * n_cols, 0);
        const auto& l_d = lines.get_data_array();
        const auto& l_m = lines.get_mask_array();

        std::vector<std::uint8_t> snapped(burned.size(), 0);

        for (int row = 0; row < n_rows; ++row)
            for (int col = 0; col < n_cols; ++col)
            {
                if (l_m(row, col) == 0) continue;
                if (l_d(row, col) == 0.0f) continue;
                ++out.n_burned;

                // Snap to the largest drainage area within the search radius.
                // The survey and the DEM disagree by a pixel or two; without
                // this, a line that sits one cell off the modelled flow path
                // paints a hillslope as channel.
                int best_r = row, best_c = col;
                double best_a = (a_m(row, col) != 0)
                    ? static_cast<double>(a_d(row, col)) : -1.0;

                const int rad = opts.snap_search_radius_pixels;
                for (int dr = -rad; dr <= rad; ++dr)
                    for (int dc = -rad; dc <= rad; ++dc)
                    {
                        const int rr = row + dr, cc = col + dc;
                        if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) continue;
                        if (a_m(rr, cc) == 0) continue;
                        const double a = static_cast<double>(a_d(rr, cc));
                        if (a > best_a) { best_a = a; best_r = rr; best_c = cc; }
                    }

                if (best_a < 0.0) continue;                   // nothing valid nearby
                if (opts.snap_min_area_m2 > 0.0 && best_a < opts.snap_min_area_m2)
                {
                    ++out.n_rejected_by_area;
                    continue;
                }
                if (best_r != row || best_c != col) ++out.n_snapped;

                burned[static_cast<std::size_t>(best_r) * n_cols + best_c] = 1;
            }
    }

    for (int row = 0; row < n_rows; ++row)
        for (int col = 0; col < n_cols; ++col)
        {
            if (a_m(row, col) == 0) continue;
            const double A = static_cast<double>(a_d(row, col));
            if (!(A > 0.0)) continue;

            ++out.n_valid;
            bool is_channel = false;

            switch (opts.mode)
            {
                case channel_mask_mode::none:
                    // Where hydraulic-geometry depth has overtaken sheet-flow
                    // depth. This is a DIAGNOSTIC of the geometry, not an input
                    // to it -- the velocity closure already transitions here
                    // smoothly whether or not anyone asks for the mask.
                    is_channel = (A >= emergent_area);
                    break;

                case channel_mask_mode::area_threshold:
                    is_channel = (A >= opts.channel_area_threshold_m2);
                    break;

                case channel_mask_mode::area_slope_threshold:
                {
                    if (slope->get_mask_array()(row, col) == 0) break;
                    const double S = std::max(
                        static_cast<double>(slope->get_data_array()(row, col)),
                        opts.min_slope);
                    is_channel = (A * std::pow(S, opts.area_slope_theta) >=
                                  opts.area_slope_threshold);
                    break;
                }

                case channel_mask_mode::extraction:
                {
                    if (supplied_mask->get_mask_array()(row, col) == 0) break;
                    is_channel =
                        (supplied_mask->get_data_array()(row, col) != 0.0f);
                    break;
                }

                case channel_mask_mode::vector:
                    is_channel =
                        burned[static_cast<std::size_t>(row) * n_cols + col] != 0;
                    break;
            }

            m_data(row, col) = is_channel ? 1.0f : 0.0f;
            m_mask(row, col) = std::uint8_t(1);
            if (is_channel) ++out.n_channel;
        }

    out.mask = lsdtt_raster_cpu(std::move(m_data), std::move(m_mask), meta);
    return out;
}

} // namespace lsdtt
