// lsdtt_raster_binary_mask_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html


#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <deque>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "lsdtt_proximity_cpu.hpp"  // euclidean_distance_transform

namespace lsdtt
{
namespace
{




// -----------------------------------------------------------------------------
// filter_binary_components_by_shape
//
// Input must be a binary raster with values 0/1 on valid pixels.
// Connected components are 8-neighbourhood foreground regions.
// Components are retained if they pass all of:
//   - minimum pixel count
//   - minimum major-axis length (pixels)
//   - minimum elongation = major / minor
//
// The major/minor axes are estimated from the eigenvalues of the 2D covariance
// matrix of pixel-centre coordinates within each component.
// -----------------------------------------------------------------------------
[[maybe_unused]] static lsdtt_raster_cpu filter_binary_components_by_shape(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels,
    double                  min_major_axis_pixels,
    double                  min_elongation)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;


    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});
    std::vector<uint8_t> visited(npix, uint8_t{0});

    auto is_valid = [&](std::size_t idx) -> bool
    {
        return in_mask[idx] != 0 &&
               !meta.is_nodata(in_data[idx]) &&
               std::isfinite(in_data[idx]);
    };

    auto is_foreground = [&](std::size_t idx) -> bool
    {
        return is_valid(idx) && in_data[idx] > 0.5f;
    };

    for (std::size_t idx = 0; idx < npix; ++idx)
    {
        if (is_valid(idx))
        {
            out_mask[idx] = uint8_t{1};
        }
    }

    for (int row0 = 0; row0 < nrows; ++row0)
    {
        for (int col0 = 0; col0 < ncols; ++col0)
        {
            const std::size_t seed =
                static_cast<std::size_t>(row0) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(col0);

            if (visited[seed] != 0 || !is_foreground(seed))
            {
                continue;
            }

            std::vector<std::size_t> stack;
            std::vector<std::size_t> component;

            stack.push_back(seed);
            visited[seed] = uint8_t{1};

            double sum_r  = 0.0;
            double sum_c  = 0.0;
            double sum_rr = 0.0;
            double sum_cc = 0.0;
            double sum_rc = 0.0;

            while (!stack.empty())
            {
                const std::size_t idx = stack.back();
                stack.pop_back();

                component.push_back(idx);

                const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
                const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));

                const double r = static_cast<double>(row);
                const double c = static_cast<double>(col);

                sum_r  += r;
                sum_c  += c;
                sum_rr += r * r;
                sum_cc += c * c;
                sum_rc += r * c;

                for (int dr = -1; dr <= 1; ++dr)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;

                        const int rr = row + dr;
                        const int cc = col + dc;

                        if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols)
                        {
                            continue;
                        }

                        const std::size_t nidx =
                            static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) +
                            static_cast<std::size_t>(cc);

                        if (visited[nidx] == 0 && is_foreground(nidx))
                        {
                            visited[nidx] = uint8_t{1};
                            stack.push_back(nidx);
                        }
                    }
                }
            }

            const double n = static_cast<double>(component.size());

            double major_axis_pixels = 0.0;
            double minor_axis_pixels = 0.0;
            double elongation        = 1.0;

            if (n >= 2.0)
            {
                const double mean_r = sum_r / n;
                const double mean_c = sum_c / n;

                const double cov_rr = std::max(0.0, sum_rr / n - mean_r * mean_r);
                const double cov_cc = std::max(0.0, sum_cc / n - mean_c * mean_c);
                const double cov_rc =               sum_rc / n - mean_r * mean_c;

                const double trace = cov_rr + cov_cc;
                const double disc = std::sqrt(
                    std::max(0.0, (cov_rr - cov_cc) * (cov_rr - cov_cc) +
                                  4.0 * cov_rc * cov_rc)
                );

                const double lambda_max = std::max(0.0, 0.5 * (trace + disc));
                const double lambda_min = std::max(0.0, 0.5 * (trace - disc));

                major_axis_pixels = 4.0 * std::sqrt(lambda_max);
                minor_axis_pixels = 4.0 * std::sqrt(lambda_min);

                if (minor_axis_pixels > 1e-12)
                {
                    elongation = major_axis_pixels / minor_axis_pixels;
                }
                else if (major_axis_pixels > 0.0)
                {
                    elongation = std::numeric_limits<double>::infinity();
                }
            }

            const bool keep =
                static_cast<int>(component.size()) >= min_component_pixels &&
                major_axis_pixels >= min_major_axis_pixels &&
                elongation >= min_elongation;

            if (keep)
            {
                for (std::size_t idx : component)
                {
                    out_data[idx] = 1.0f;
                }
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_directional_closing - internal helpers
// =============================================================================

// Build the pixel offsets for a symmetric line SE of the given half-length
// at angle_degrees from horizontal (measured counterclockwise).
// In raster coords row increases downward, so the vertical component is negated.
// Duplicate offsets (produced by rounding at certain angles) are removed.
static std::vector<std::pair<int,int>> make_line_se_offsets(
    int half_length, double angle_degrees)
{
    const double theta  = angle_degrees * (M_PI / 180.0);
    const double cos_t  = std::cos(theta);
    const double sin_t  = std::sin(theta);

    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(static_cast<std::size_t>(2 * half_length + 1));

    for (int step = -half_length; step <= half_length; ++step)
    {
        const double t = static_cast<double>(step);
        const int dc   = static_cast<int>(std::round( t * cos_t));
        const int dr   = static_cast<int>(std::round(-t * sin_t));
        offsets.emplace_back(dr, dc);
    }

    // Remove duplicates that arise from rounding (common at multiples of 45 deg).
    std::sort(offsets.begin(), offsets.end());
    offsets.erase(std::unique(offsets.begin(), offsets.end()), offsets.end());

    return offsets;
}

// Binary dilation: pixel is foreground if ANY valid SE-neighbour in in_data is.
// Returns a uint8_t plane (1 = fg, 0 = bg/nodata).
static std::vector<uint8_t> binary_dilate_se(
    const std::vector<float>&              in_data,
    const std::vector<uint8_t>&            in_mask,
    const std::vector<std::pair<int,int>>& offsets,
    int nrows, int ncols)
{
    const std::size_t npix = static_cast<std::size_t>(nrows) *
                             static_cast<std::size_t>(ncols);
    std::vector<uint8_t> dilated(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(ncols) +
                                    static_cast<std::size_t>(col);
            if (in_mask[idx] == 0) continue;

            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;
                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;
                const std::size_t nidx = static_cast<std::size_t>(rr) *
                                         static_cast<std::size_t>(ncols) +
                                         static_cast<std::size_t>(cc);
                if (in_mask[nidx] != 0 && in_data[nidx] > 0.5f)
                {
                    dilated[idx] = uint8_t{1};
                    break;
                }
            }
        }
    }
    return dilated;
}

// Binary erosion on a uint8_t plane produced by binary_dilate_se.
// A pixel stays foreground only if ALL valid-mask SE-neighbours are also
// foreground.  Nodata neighbours (in_mask==0) are skipped so that valid pixels
// near the DEM border are not incorrectly eroded.
static std::vector<uint8_t> binary_erode_se(
    const std::vector<uint8_t>&            dilated,
    const std::vector<uint8_t>&            in_mask,
    const std::vector<std::pair<int,int>>& offsets,
    int nrows, int ncols)
{
    const std::size_t npix = static_cast<std::size_t>(nrows) *
                             static_cast<std::size_t>(ncols);
    std::vector<uint8_t> eroded(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(ncols) +
                                    static_cast<std::size_t>(col);
            if (in_mask[idx] == 0 || dilated[idx] == 0) continue;

            bool all_fg = true;
            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;
                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;
                const std::size_t nidx = static_cast<std::size_t>(rr) *
                                         static_cast<std::size_t>(ncols) +
                                         static_cast<std::size_t>(cc);
                if (in_mask[nidx] == 0) continue;   // ignore nodata neighbours
                if (dilated[nidx] == 0) { all_fg = false; break; }
            }

            if (all_fg) eroded[idx] = uint8_t{1};
        }
    }
    return eroded;
}

} // unnamed namespace






// .##...##...####...#####...#####...##..##...####...##.......####....####...##..##.
// .###.###..##..##..##..##..##..##..##..##..##..##..##......##..##..##.......####..
// .##.#.##..##..##..#####...#####...######..##..##..##......##..##..##.###....##...
// .##...##..##..##..##..##..##......##..##..##..##..##......##..##..##..##....##...
// .##...##...####...##..##..##......##..##...####...######...####....####.....##...

// -----------------------------------------------------------------------------
// binary_morphological_closing
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_morphological_closing(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::binary_morphological_closing: input raster metadata is invalid"
        );
    }

    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::binary_morphological_closing: radius_metres must be finite and > 0"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const double cell_size = meta.cell_size;

    const int pixel_radius = std::max(
        1, static_cast<int>(std::ceil(radius_metres / cell_size))
    );

    // Pre-build circular kernel offsets
    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(
        static_cast<std::size_t>(2 * pixel_radius + 1) *
        static_cast<std::size_t>(2 * pixel_radius + 1)
    );
    for (int dr = -pixel_radius; dr <= pixel_radius; ++dr)
    {
        for (int dc = -pixel_radius; dc <= pixel_radius; ++dc)
        {
            const double dx   = static_cast<double>(dc) * cell_size;
            const double dy   = static_cast<double>(dr) * cell_size;
            const double dist = std::sqrt(dx * dx + dy * dy);
            if (dist <= radius_metres + 1e-12)
            {
                offsets.emplace_back(dr, dc);
            }
        }
    }

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    // -------------------------------------------------------------------------
    // Dilation pass: pixel is foreground if ANY valid kernel neighbour is
    // -------------------------------------------------------------------------
    std::vector<uint8_t> dilated(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
                + static_cast<std::size_t>(col);

            if (in_mask[idx] == 0) continue;

            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;

                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;

                const std::size_t nidx =
                    static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols)
                    + static_cast<std::size_t>(cc);

                if (in_mask[nidx] != 0 && in_data[nidx] > 0.5f)
                {
                    dilated[idx] = uint8_t{1};
                    break;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Erosion pass: pixel stays foreground only if ALL valid kernel neighbours
    // are also foreground in the dilated result.  Nodata neighbours are skipped
    // so that valid pixels near the DEM border are not incorrectly eroded.
    // -------------------------------------------------------------------------
    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
                + static_cast<std::size_t>(col);

            if (in_mask[idx] == 0) continue;
            out_mask[idx] = uint8_t{1};

            if (dilated[idx] == 0)
            {
                // Was background even after dilation - stays 0
                continue;
            }

            bool all_valid_on = true;
            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;

                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;

                const std::size_t nidx =
                    static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols)
                    + static_cast<std::size_t>(cc);

                if (in_mask[nidx] == 0) continue;  // ignore nodata neighbours

                if (dilated[nidx] == 0)
                {
                    all_valid_on = false;
                    break;
                }
            }

            if (all_valid_on)
            {
                out_data[idx] = 1.0f;
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_morphological_opening
// =============================================================================
lsdtt_raster_cpu binary_morphological_opening(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_morphological_opening: input raster metadata is invalid");
    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
        throw std::runtime_error(
            "lsdtt::binary_morphological_opening: radius_metres must be finite and > 0");

    const lsdtt_raster_cpu eroded =
        binary_morphological_erosion(binary, radius_metres, shape);

    return binary_morphological_dilation(eroded, radius_metres, shape);
}




// =============================================================================
// binary_morphological_dilation
// =============================================================================
lsdtt_raster_cpu binary_morphological_dilation(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_morphological_dilation: input raster metadata is invalid");
    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
        throw std::runtime_error(
            "lsdtt::binary_morphological_dilation: radius_metres must be finite and > 0");

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const double cs    = meta.cell_size;
    const int    r_px  = std::max(1, static_cast<int>(std::ceil(radius_metres / cs)));

    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(static_cast<std::size_t>(2*r_px+1) * static_cast<std::size_t>(2*r_px+1));
    for (int dr = -r_px; dr <= r_px; ++dr)
    {
        for (int dc = -r_px; dc <= r_px; ++dc)
        {
            bool in_se = (shape == window_shape::square)
                ? (std::abs(dr) <= r_px && std::abs(dc) <= r_px)
                : (std::sqrt(static_cast<double>(dr)*dr + static_cast<double>(dc)*dc) * cs
                   <= radius_metres + 1e-12);
            if (in_se) offsets.emplace_back(dr, dc);
        }
    }

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) + col;
            if (in_mask[idx] == 0) continue;
            out_mask[idx] = uint8_t{1};
            out_data[idx] = 0.0f;

            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;
                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;
                const std::size_t nidx =
                    static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) + cc;
                if (in_mask[nidx] != 0 && in_data[nidx] > 0.5f)
                {
                    out_data[idx] = 1.0f;
                    break;
                }
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_morphological_erosion
// =============================================================================
lsdtt_raster_cpu binary_morphological_erosion(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_morphological_erosion: input raster metadata is invalid");
    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
        throw std::runtime_error(
            "lsdtt::binary_morphological_erosion: radius_metres must be finite and > 0");

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const double cs    = meta.cell_size;
    const int    r_px  = std::max(1, static_cast<int>(std::ceil(radius_metres / cs)));

    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(static_cast<std::size_t>(2*r_px+1) * static_cast<std::size_t>(2*r_px+1));
    for (int dr = -r_px; dr <= r_px; ++dr)
    {
        for (int dc = -r_px; dc <= r_px; ++dc)
        {
            bool in_se = (shape == window_shape::square)
                ? (std::abs(dr) <= r_px && std::abs(dc) <= r_px)
                : (std::sqrt(static_cast<double>(dr)*dr + static_cast<double>(dc)*dc) * cs
                   <= radius_metres + 1e-12);
            if (in_se) offsets.emplace_back(dr, dc);
        }
    }

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) + col;
            if (in_mask[idx] == 0) continue;
            out_mask[idx] = uint8_t{1};
            out_data[idx] = 0.0f;
            if (in_data[idx] <= 0.5f) continue;  // background stays background

            bool all_fg = true;
            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;
                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;
                const std::size_t nidx =
                    static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) + cc;
                if (in_mask[nidx] == 0) continue;  // skip nodata neighbours
                if (in_data[nidx] <= 0.5f) { all_fg = false; break; }
            }
            if (all_fg) out_data[idx] = 1.0f;
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_opening_by_reconstruction
// =============================================================================
lsdtt_raster_cpu binary_opening_by_reconstruction(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape)
{
    // Step 1: erode to get seeds - pixels that survive the radius filter.
    const lsdtt_raster_cpu eroded =
        binary_morphological_erosion(binary, radius_metres, shape);

    const lsdtt_raster_meta& meta  = binary.get_meta();
    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   orig_data = binary.get_data();
    const std::vector<uint8_t>& orig_mask = binary.get_mask();
    const std::vector<float>&   seed_data = eroded.get_data();
    const std::vector<uint8_t>& seed_mask = eroded.get_mask();

    const float nd = static_cast<float>(meta.no_data_value);

    // Output inherits the original nodata structure; valid pixels start as
    // background (0) and are set to foreground (1) only by the BFS.
    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask = orig_mask;
    for (std::size_t i = 0; i < npix; ++i)
        if (orig_mask[i] != 0) out_data[i] = 0.0f;

    // Step 2: BFS from eroded-foreground seeds, constrained to original
    // foreground.  Uses 8-connectivity so diagonal necks are also restored.
    static constexpr const auto& DR8 = grid::d8_rowmajor_row;
    static constexpr const auto& DC8 = grid::d8_rowmajor_col;

    std::vector<bool> visited(npix, false);
    std::queue<int>   q;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (seed_mask[i] == 0 || seed_data[i] <= 0.5f) continue;
        if (orig_mask[i] == 0 || orig_data[i] <= 0.5f) continue; // guard
        visited[i]  = true;
        out_data[i] = 1.0f;
        q.push(static_cast<int>(i));
    }

    while (!q.empty())
    {
        const int idx = q.front(); q.pop();
        const int row = idx / ncols;
        const int col = idx % ncols;
        for (int d = 0; d < 8; ++d)
        {
            const int nr = row + DR8[d];
            const int nc = col + DC8[d];
            if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
            const int nidx = nr * ncols + nc;
            if (visited[nidx])             continue;
            if (orig_mask[nidx] == 0)      continue; // nodata boundary
            if (orig_data[nidx] <= 0.5f)   continue; // not in domain
            visited[nidx]    = true;
            out_data[nidx]   = 1.0f;
            q.push(nidx);
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_closing_by_reconstruction
// =============================================================================
lsdtt_raster_cpu binary_closing_by_reconstruction(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_closing_by_reconstruction: input raster metadata is invalid");
    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
        throw std::runtime_error(
            "lsdtt::binary_closing_by_reconstruction: radius_metres must be finite and > 0");

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const float nd = static_cast<float>(meta.no_data_value);

    // -------------------------------------------------------------------------
    // Step 1: complement over the DEM-valid domain
    // -------------------------------------------------------------------------
    std::vector<float>   comp_data(npix, nd);
    std::vector<uint8_t> comp_mask = in_mask;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        comp_data[i] = (in_data[i] > 0.5f) ? 0.0f : 1.0f;
    }

    const lsdtt_raster_cpu complement(
        std::move(comp_data), std::move(comp_mask), meta);

    // -------------------------------------------------------------------------
    // Step 2: opening-by-reconstruction on the complement
    // -------------------------------------------------------------------------
    const lsdtt_raster_cpu opened_complement =
        binary_opening_by_reconstruction(complement, radius_metres, shape);

    // -------------------------------------------------------------------------
    // Step 3: complement back
    // -------------------------------------------------------------------------
    const std::vector<float>&   oc_data = opened_complement.get_data();
    const std::vector<uint8_t>& oc_mask = opened_complement.get_mask();

    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask = oc_mask;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (oc_mask[i] == 0) continue;
        out_data[i] = (oc_data[i] > 0.5f) ? 0.0f : 1.0f;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}



// -----------------------------------------------------------------------------
//.##...##...####....####...##..##...........####...##......######...####...##..##..##..##..#####..
//.###.###..##..##..##......##.##...........##..##..##......##......##..##..###.##..##..##..##..##.
//.##.#.##..######...####...####............##......##......####....######..##.###..##..##..#####..
//.##...##..##..##......##..##.##...........##..##..##......##......##..##..##..##..##..##..##.....
//.##...##..##..##...####...##..##...........####...######..######..##..##..##..##...####...##.....
// -----------------------------------------------------------------------------
// =============================================================================
// binary_fill_holes
// =============================================================================
lsdtt_raster_cpu binary_fill_holes(
    const lsdtt_raster_cpu& binary)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_fill_holes: input raster metadata is invalid");

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    // Mark DEM-valid background pixels reachable from the DEM border via
    // 4-connected flood fill.  Unreachable valid background pixels are holes.
    std::vector<bool> reachable(npix, false);
    std::queue<int>   q;

    auto try_enqueue = [&](int r, int c)
    {
        if (r < 0 || r >= nrows || c < 0 || c >= ncols) return;
        const int i = r * ncols + c;
        if (in_mask[i] == 0 || reachable[i]) return;
        if (in_data[i] > 0.5f) return;   // foreground - not background
        reachable[i] = true;
        q.push(i);
    };

    for (int c = 0; c < ncols; ++c) { try_enqueue(0, c); try_enqueue(nrows-1, c); }
    for (int r = 0; r < nrows;  ++r) { try_enqueue(r, 0); try_enqueue(r, ncols-1); }

    static constexpr const auto& DR = grid::d4_nswe_row;
    static constexpr const auto& DC = grid::d4_nswe_col;

    // Also seed valid background pixels adjacent to any nodata pixel - they sit
    // at the effective exterior boundary of the valid DEM region and are outside
    // any real hole even when the grid border is nodata.
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const int i = r * ncols + c;
            if (in_mask[i] == 0 || reachable[i] || in_data[i] > 0.5f) continue;
            for (int d = 0; d < 4; ++d)
            {
                const int nr = r + DR[d], nc = c + DC[d];
                if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                if (in_mask[nr * ncols + nc] == 0)
                { try_enqueue(r, c); break; }
            }
        }
    }
    while (!q.empty())
    {
        const int i   = q.front(); q.pop();
        const int row = i / ncols;
        const int col = i % ncols;
        for (int d = 0; d < 4; ++d) try_enqueue(row + DR[d], col + DC[d]);
    }

    // Unreachable valid background pixels -> foreground.
    // Nodata pixels (in_mask==0) get the nodata sentinel in the data array.
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask = in_mask;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (in_data[i] <= 0.5f && !reachable[i])
            out_data[i] = 1.0f;   // hole -> foreground
        else
            out_data[i] = in_data[i];
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}




// =============================================================================
// binary_directional_closing
// =============================================================================
lsdtt_raster_cpu binary_directional_closing(
    const lsdtt_raster_cpu& binary,
    int                     line_length_pixels,
    int                     n_angles)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_directional_closing: input raster metadata is invalid"
        );
    if (line_length_pixels < 3)
        throw std::invalid_argument(
            "lsdtt::binary_directional_closing: line_length_pixels must be >= 3"
        );
    if (n_angles < 1)
        throw std::invalid_argument(
            "lsdtt::binary_directional_closing: n_angles must be >= 1"
        );

    // Force odd length so the SE is symmetric about the origin pixel.
    if (line_length_pixels % 2 == 0) ++line_length_pixels;
    const int half_len = (line_length_pixels - 1) / 2;

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix = static_cast<std::size_t>(nrows) *
                             static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    // Accumulate OR of all directional closings.
    std::vector<uint8_t> acc(npix, uint8_t{0});

    for (int ai = 0; ai < n_angles; ++ai)
    {
        // Evenly spaced orientations in [0 deg, 180 deg).
        // A line at theta and theta+180 deg is the same, so we only need half the circle.
        const double angle = static_cast<double>(ai) * 180.0 /
                             static_cast<double>(n_angles);

        const auto se      = make_line_se_offsets(half_len, angle);
        const auto dilated = binary_dilate_se(in_data, in_mask, se, nrows, ncols);
        const auto closed  = binary_erode_se (dilated,  in_mask, se, nrows, ncols);

        for (std::size_t i = 0; i < npix; ++i)
            if (closed[i] != 0) acc[i] = uint8_t{1};
    }

    // Convert accumulator to float output raster.
    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        out_mask[i] = uint8_t{1};
        if (acc[i] != 0) out_data[i] = 1.0f;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}




// ..####....####...##..##..##..##..######...####...######..######..#####.......... 
// .##..##..##..##..###.##..###.##..##......##..##....##....##......##..##......... 
// .##......##..##..##.###..##.###..####....##........##....####....##..##......... 
// .##..##..##..##..##..##..##..##..##......##..##....##....##......##..##......... 
// ..####....####...##..##..##..##..######...####.....##....######..#####.......... 
// ................................................................................ 
// ..####....####...##...##..#####....####...##..##..######..##..##..######...####..
// .##..##..##..##..###.###..##..##..##..##..###.##..##......###.##....##....##.....
// .##......##..##..##.#.##..#####...##..##..##.###..####....##.###....##.....####..
// .##..##..##..##..##...##..##......##..##..##..##..##......##..##....##........##.
// ..####....####...##...##..##.......####...##..##..######..##..##....##.....####..

// -----------------------------------------------------------------------------
// filter_binary_components_by_size
// -----------------------------------------------------------------------------
lsdtt_raster_cpu filter_binary_components_by_size(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels,
    pixel_connectivity      connectivity)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::filter_binary_components_by_size: input raster metadata is invalid"
        );
    }

    if (min_component_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::filter_binary_components_by_size: min_component_pixels must be >= 1"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});
    std::vector<uint8_t> visited(npix, uint8_t{0});

    auto is_valid = [&](std::size_t idx) -> bool
    {
        return in_mask[idx] != 0 &&
               !meta.is_nodata(in_data[idx]) &&
               std::isfinite(in_data[idx]);
    };

    auto is_foreground = [&](std::size_t idx) -> bool
    {
        return is_valid(idx) && in_data[idx] > 0.5f;
    };

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (is_valid(i))
        {
            out_mask[i] = uint8_t{1};
        }
    }

    for (int row0 = 0; row0 < nrows; ++row0)
    {
        for (int col0 = 0; col0 < ncols; ++col0)
        {
            const std::size_t seed =
                static_cast<std::size_t>(row0) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(col0);

            if (visited[seed] != 0 || !is_foreground(seed))
            {
                continue;
            }

            std::vector<std::size_t> stack;
            std::vector<std::size_t> component;

            stack.push_back(seed);
            visited[seed] = uint8_t{1};

            while (!stack.empty())
            {
                const std::size_t idx = stack.back();
                stack.pop_back();
                component.push_back(idx);

                const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
                const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));

                for (int dr = -1; dr <= 1; ++dr)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;
                        if (connectivity == pixel_connectivity::d4 &&
                            dr != 0 && dc != 0) continue;

                        const int rr = row + dr;
                        const int cc = col + dc;

                        if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols)
                        {
                            continue;
                        }

                        const std::size_t nidx =
                            static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) +
                            static_cast<std::size_t>(cc);

                        if (visited[nidx] == 0 && is_foreground(nidx))
                        {
                            visited[nidx] = uint8_t{1};
                            stack.push_back(nidx);
                        }
                    }
                }
            }

            if (static_cast<int>(component.size()) >= min_component_pixels)
            {
                for (std::size_t idx : component)
                {
                    out_data[idx] = 1.0f;
                }
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// -----------------------------------------------------------------------------
// label_binary_components_by_size
// -----------------------------------------------------------------------------
lsdtt_raster_cpu label_binary_components_by_size(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels,
    pixel_connectivity      connectivity)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::label_binary_components_by_size: input raster metadata is invalid"
        );
    }

    if (min_component_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::label_binary_components_by_size: min_component_pixels must be >= 1"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    // Output is UInt32: labels are positive integers, nodata = 0.
    lsdtt_raster_meta out_meta   = meta;
    out_meta.data_type           = lsdtt::lsdtt_pixel_type::uint32;
    out_meta.no_data_value       = 0.0;
    const float nodata           = 0.0f;

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});
    std::vector<uint8_t> visited(npix, uint8_t{0});

    auto is_valid = [&](std::size_t idx) -> bool
    {
        return in_mask[idx] != 0 &&
               !meta.is_nodata(in_data[idx]) &&
               std::isfinite(in_data[idx]);
    };

    auto is_foreground = [&](std::size_t idx) -> bool
    {
        return is_valid(idx) && in_data[idx] > 0.5f;
    };

    int next_label = 1;

    for (int row0 = 0; row0 < nrows; ++row0)
    {
        for (int col0 = 0; col0 < ncols; ++col0)
        {
            const std::size_t seed =
                static_cast<std::size_t>(row0) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(col0);

            if (visited[seed] != 0 || !is_foreground(seed))
            {
                continue;
            }

            std::vector<std::size_t> stack;
            std::vector<std::size_t> component;

            stack.push_back(seed);
            visited[seed] = uint8_t{1};

            while (!stack.empty())
            {
                const std::size_t idx = stack.back();
                stack.pop_back();
                component.push_back(idx);

                const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
                const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));

                for (int dr = -1; dr <= 1; ++dr)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;
                        if (connectivity == pixel_connectivity::d4 &&
                            dr != 0 && dc != 0) continue;

                        const int rr = row + dr;
                        const int cc = col + dc;

                        if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols)
                        {
                            continue;
                        }

                        const std::size_t nidx =
                            static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) +
                            static_cast<std::size_t>(cc);

                        if (visited[nidx] == 0 && is_foreground(nidx))
                        {
                            visited[nidx] = uint8_t{1};
                            stack.push_back(nidx);
                        }
                    }
                }
            }

            if (static_cast<int>(component.size()) >= min_component_pixels)
            {
                const float label = static_cast<float>(next_label);
                ++next_label;

                for (std::size_t idx : component)
                {
                    out_data[idx] = label;
                    out_mask[idx] = uint8_t{1};
                }
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), out_meta);
}

// Set / unset sentinels for the byte masks this operation works on. Local to
// the file: the public interface is rasters, where > 0.5 means set.
namespace
{
constexpr std::uint8_t k_off = 0u;
constexpr std::uint8_t k_on  = 1u;

// 8-neighbour offsets, row then column.
static constexpr const auto& k_dr = grid::d8_rowmajor_row;
static constexpr const auto& k_dc = grid::d8_rowmajor_col;
}

// -----------------------------------------------------------------------------
// filter_components_by_seed_reachability_in_place
//
// A morphological opening CONSTRAINED BY DISTANCE TO A SEED SET. Generic: the
// seed is whatever the caller says it is (a channel network, a coastline, a
// road centreline), and nothing here knows what the mask represents.
//
//   1. erode by n_dilate iterations (the seed and its immediate neighbours are
//      pinned on, so the mask never detaches from what anchors it);
//   2. label 8-connected patches and drop those with no pixel on or 8-adjacent
//      to a seed, and those whose along-ground path to a seed exceeds
//      max_detour_ratio times the straight-line distance - a patch that is NEAR
//      a seed but only reachable the long way round is a different object;
//   3. dilate 2 * n_dilate times (regrow the eroded rim and fill interior
//      holes), retaining a dilated pixel only when its distance to a seed does
//      not exceed the farthest pixel of the patch it grows from, plus a margin
//      covering the erosion depth - so patches regrow and holes fill, but
//      nothing extends outward away from the seed.
//
// Step 3's cap is what stops the regrow creeping outward: a distance cap alone
// does not bind near the seed, where the rim would otherwise climb whatever
// rises beside it.
//
// Operates in place on a 0/1 byte mask, because the callers that need it in the
// inner loop already hold one. See the raster overload for the ordinary case.
// -----------------------------------------------------------------------------
void filter_components_by_seed_reachability_in_place(
    std::vector<std::uint8_t>&       mask,
    const std::vector<float>&        distance_to_seed,
    const std::vector<std::uint8_t>& distance_mask,
    const std::vector<std::int64_t>& seed_flats,
    int n_dilate, bool require_seed_adjacent, double max_detour_ratio,
    bool require_seed_line_of_sight,
    double cell, int n_rows, int n_cols)
{
    const std::size_t n = mask.size();
    const int n_erode = n_dilate;

    // Growth domain for the final dilation: the ORIGINAL (pre-erosion)
    // footprint plus its enclosed holes. The dilation must never expand the
    // mask's outer boundary - a distance cap alone does not bind near the
    // seed and the regrow rim creeps up the risers. Enclosed holes =
    // background not 4-connected to the tile border across non-mask
    // ground.
    std::vector<std::uint8_t> allowed(mask);
    {
        std::vector<std::uint8_t> outside(n, k_off);
        std::deque<std::int64_t> oq;
        auto push_bg = [&](int r, int c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * n_cols + c;
            if (!mask[i] && !outside[i])
            {
                outside[i] = k_on;
                oq.push_back(static_cast<std::int64_t>(i));
            }
        };
        for (int c = 0; c < n_cols; ++c)
        {
            push_bg(0, c);
            push_bg(n_rows - 1, c);
        }
        for (int r = 0; r < n_rows; ++r)
        {
            push_bg(r, 0);
            push_bg(r, n_cols - 1);
        }
        while (!oq.empty())
        {
            const std::int64_t flat = oq.front();
            oq.pop_front();
            const int r = static_cast<int>(flat / n_cols);
            const int c = static_cast<int>(flat % n_cols);
            if (r > 0)          push_bg(r - 1, c);
            if (r + 1 < n_rows) push_bg(r + 1, c);
            if (c > 0)          push_bg(r, c - 1);
            if (c + 1 < n_cols) push_bg(r, c + 1);
        }
        for (std::size_t i = 0; i < n; ++i)
            if (!mask[i] && !outside[i]) allowed[i] = k_on;
    }

    // 0. Erode (seed + neighbours protected) to sever connecting
    //    corridors narrower than ~2 * n_erode + 1 pixels.
    std::vector<std::uint8_t> protect(n, k_off);
    for (const std::int64_t cf : seed_flats)
    {
        if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
        protect[static_cast<std::size_t>(cf)] = k_on;
        const int r = static_cast<int>(cf / n_cols);
        const int c = static_cast<int>(cf % n_cols);
        for (int k = 0; k < 8; ++k)
        {
            const int rr = r + k_dr[k], cc = c + k_dc[k];
            if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) continue;
            protect[static_cast<std::size_t>(rr) * n_cols + cc] = k_on;
        }
    }
    for (int it = 0; it < n_erode; ++it)
    {
        std::vector<std::uint8_t> eroded = mask;
        bool changed = false;
        for (int r = 0; r < n_rows; ++r)
        {
            for (int c = 0; c < n_cols; ++c)
            {
                const std::size_t i =
                    static_cast<std::size_t>(r) * n_cols + c;
                if (!mask[i] || protect[i]) continue;
                for (int k = 0; k < 8; ++k)
                {
                    const int rr = r + k_dr[k], cc = c + k_dc[k];
                    if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols ||
                        !mask[static_cast<std::size_t>(rr) * n_cols + cc])
                    {
                        eroded[i] = k_off;
                        changed = true;
                        break;
                    }
                }
            }
        }
        mask.swap(eroded);
        if (!changed) break;
    }
    n_dilate = 2 * n_dilate;   // regrow the rim + fill holes

    // 1. Label patches.
    std::vector<std::int32_t> lbl(n, -1);
    std::vector<double> patch_max_dist;
    std::deque<std::int64_t> queue;
    for (std::size_t seed = 0; seed < n; ++seed)
    {
        if (!mask[seed] || lbl[seed] >= 0) continue;
        const std::int32_t id =
            static_cast<std::int32_t>(patch_max_dist.size());
        patch_max_dist.push_back(0.0);
        lbl[seed] = id;
        queue.push_back(static_cast<std::int64_t>(seed));
        while (!queue.empty())
        {
            const std::int64_t flat = queue.front();
            queue.pop_front();
            const std::size_t fi = static_cast<std::size_t>(flat);
            if (distance_mask[fi] && std::isfinite(distance_to_seed[fi]))
                patch_max_dist[id] = std::max(
                    patch_max_dist[id],
                    static_cast<double>(distance_to_seed[fi]));
            const int r = static_cast<int>(flat / n_cols);
            const int c = static_cast<int>(flat % n_cols);
            for (int k = 0; k < 8; ++k)
            {
                const int rr = r + k_dr[k], cc = c + k_dc[k];
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                    continue;
                const std::size_t nf =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (mask[nf] && lbl[nf] < 0)
                {
                    lbl[nf] = id;
                    queue.push_back(static_cast<std::int64_t>(nf));
                }
            }
        }
    }

    // 2. Remove patches not adjacent to a seed (before regrowth, so a
    //    severed far-off patch cannot re-bridge during dilation).
    if (require_seed_adjacent && !patch_max_dist.empty())
    {
        std::vector<char> adjacent(patch_max_dist.size(), 0);
        for (const std::int64_t cf : seed_flats)
        {
            const int r = static_cast<int>(cf / n_cols);
            const int c = static_cast<int>(cf % n_cols);
            for (int k = -1; k < 8; ++k)
            {
                const int rr = (k < 0) ? r : r + k_dr[k];
                const int cc = (k < 0) ? c : c + k_dc[k];
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                    continue;
                const std::size_t nf =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (mask[nf] && lbl[nf] >= 0) adjacent[lbl[nf]] = 1;
            }
        }
        for (std::size_t i = 0; i < n; ++i)
            if (mask[i] && (lbl[i] < 0 || !adjacent[lbl[i]]))
            {
                mask[i] = k_off;
                lbl[i] = -1;
            }
    }

    // 2b. Detour-ratio pruning: geodesic distance from a seed THROUGH
    //     the mask mask (Dijkstra, 8-connected) must not exceed
    //     max_detour_ratio times the Euclidean distance to a seed.
    //     Genuine floodplain is directly reachable (ratio ~1); a
    //     neighbouring seed's floodplain hangs on a long corridor detour
    //     and prunes away. A ten-cell floor keeps near-seed pixels behind
    //     small obstructions.
    if (max_detour_ratio > 0.0)
    {
        constexpr float k_unreached = std::numeric_limits<float>::infinity();
        std::vector<float> geo(n, k_unreached);
        using qnode = std::pair<float, std::int64_t>;
        std::priority_queue<qnode, std::vector<qnode>,
                            std::greater<qnode>> pq;
        for (const std::int64_t cf : seed_flats)
        {
            if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
            const std::size_t fi = static_cast<std::size_t>(cf);
            if (!mask[fi]) continue;
            geo[fi] = 0.0f;
            pq.emplace(0.0f, cf);
        }
        const float straight = static_cast<float>(cell);
        const float diag = static_cast<float>(cell * grid::k_sqrt2_d);
        while (!pq.empty())
        {
            const auto [d, flat] = pq.top();
            pq.pop();
            const std::size_t fi = static_cast<std::size_t>(flat);
            if (d > geo[fi]) continue;
            const int r = static_cast<int>(flat / n_cols);
            const int c = static_cast<int>(flat % n_cols);
            for (int k = 0; k < 8; ++k)
            {
                const int rr = r + k_dr[k], cc = c + k_dc[k];
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                    continue;
                const std::size_t nf =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (!mask[nf]) continue;
                const float nd =
                    d + ((k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight);
                if (nd < geo[nf])
                {
                    geo[nf] = nd;
                    pq.emplace(nd, static_cast<std::int64_t>(nf));
                }
            }
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!mask[i]) continue;
            if (!std::isfinite(geo[i])) { mask[i] = k_off; continue; }
            const double eu =
                (distance_mask[i] && std::isfinite(distance_to_seed[i]))
                    ? static_cast<double>(distance_to_seed[i]) : 0.0;
            const double allowed_detour =
                std::max(max_detour_ratio * eu, 10.0 * cell);
            if (static_cast<double>(geo[i]) > allowed_detour) mask[i] = k_off;
        }
    }

    // 2c. Line-of-sight pruning: a genuine mask pixel sees its nearest
    //     seed pixel across low ground - the straight line between them
    //     stays inside the mask mask (small gaps tolerated: roads,
    //     embankments). A neighbouring mask's flat that is merely
    //     CONNECTED to this mask (e.g. around a broad glacial col) has
    //     intervening high ground on the direct line and prunes away.
    //
    // Gated: this is a strong filter, and on a broad meandering floodplain it
    // prunes genuine valley floor that simply has no straight line to the
    // channel. See require_seed_line_of_sight in the header.
    if (require_seed_line_of_sight)
    {
        // Nearest seed pixel per grid cell: multi-source Dijkstra over
        // the full valid grid (grid distance ~ Euclidean), carrying the
        // seed's flat index.
        std::vector<float>        gdist(n,
            std::numeric_limits<float>::infinity());
        std::vector<std::int64_t> gsrc(n, -1);
        using qnode = std::pair<float, std::int64_t>;
        std::priority_queue<qnode, std::vector<qnode>,
                            std::greater<qnode>> pq;
        for (const std::int64_t cf : seed_flats)
        {
            if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
            gdist[static_cast<std::size_t>(cf)] = 0.0f;
            gsrc[static_cast<std::size_t>(cf)]  = cf;
            pq.emplace(0.0f, cf);
        }
        const float straight = static_cast<float>(cell);
        const float diag = static_cast<float>(cell * grid::k_sqrt2_d);
        while (!pq.empty())
        {
            const auto [d, flat] = pq.top();
            pq.pop();
            const std::size_t fi = static_cast<std::size_t>(flat);
            if (d > gdist[fi]) continue;
            const int r = static_cast<int>(flat / n_cols);
            const int c = static_cast<int>(flat % n_cols);
            for (int k = 0; k < 8; ++k)
            {
                const int rr = r + k_dr[k], cc = c + k_dc[k];
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                    continue;
                const std::size_t nf =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (!distance_mask[nf]) continue;
                const float nd =
                    d + ((k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight);
                if (nd < gdist[nf])
                {
                    gdist[nf] = nd;
                    gsrc[nf]  = gsrc[fi];
                    pq.emplace(nd, static_cast<std::int64_t>(nf));
                }
            }
        }

        for (std::size_t i = 0; i < n; ++i)
        {
            if (!mask[i] || gsrc[i] < 0) continue;
            const std::int64_t tgt = gsrc[i];
            const int r0 = static_cast<int>(i / n_cols);
            const int c0 = static_cast<int>(i % n_cols);
            const int r1 = static_cast<int>(tgt / n_cols);
            const int c1 = static_cast<int>(tgt % n_cols);
            const int steps = std::max(std::abs(r1 - r0),
                                       std::abs(c1 - c0));
            if (steps <= 2) continue;
            int miss = 0;
            for (int t = 1; t < steps; ++t)
            {
                const int rr = r0 + (r1 - r0) * t / steps;
                const int cc = c0 + (c1 - c0) * t / steps;
                if (!mask[static_cast<std::size_t>(rr) * n_cols + cc])
                    ++miss;
            }
            // Tolerate small gaps: two cells or 10% of the line.
            if (miss > std::max(2, steps / 10)) mask[i] = k_off;
        }
    }

    // 3. Constrained dilation (frontier-based, label-propagating), with a
    //    margin on the cap covering the erosion depth so the rim can regrow.
    const double cap_margin = (n_erode + 0.5) * cell;
    for (int it = 0; it < n_dilate; ++it)
    {
        std::vector<std::pair<std::size_t, std::int32_t>> grown;
        for (int r = 0; r < n_rows; ++r)
        {
            for (int c = 0; c < n_cols; ++c)
            {
                const std::size_t i =
                    static_cast<std::size_t>(r) * n_cols + c;
                if (mask[i] || !allowed[i]) continue;
                if (!distance_mask[i] || !std::isfinite(distance_to_seed[i]))
                    continue;
                // Grow from the neighbouring patch with the largest reach.
                std::int32_t best = -1;
                double best_max = -1.0;
                for (int k = 0; k < 8; ++k)
                {
                    const int rr = r + k_dr[k], cc = c + k_dc[k];
                    if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                        continue;
                    const std::size_t nf =
                        static_cast<std::size_t>(rr) * n_cols + cc;
                    if (mask[nf] && lbl[nf] >= 0 &&
                        patch_max_dist[lbl[nf]] > best_max)
                    {
                        best_max = patch_max_dist[lbl[nf]];
                        best = lbl[nf];
                    }
                }
                if (best >= 0 &&
                    static_cast<double>(distance_to_seed[i]) <=
                        best_max + cap_margin)
                    grown.emplace_back(i, best);
            }
        }
        if (grown.empty()) break;
        for (const auto& [i, id] : grown)
        {
            mask[i] = k_on;
            lbl[i] = id;
        }
    }
}
// -----------------------------------------------------------------------------
// filter_binary_components_by_seed_reachability (public - see header)
// -----------------------------------------------------------------------------
lsdtt_raster_cpu filter_binary_components_by_seed_reachability(
    const lsdtt_raster_cpu&          mask,
    const lsdtt_raster_cpu&          seed_mask,
    const seed_reachability_options& options)
{
    static const char* fn = "lsdtt::filter_binary_components_by_seed_reachability";

    const lsdtt_raster_meta& meta = mask.get_meta();
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    if (seed_mask.get_data().size() != n)
    {
        throw std::invalid_argument(
            std::string(fn) + ": seed_mask does not match mask");
    }

    const std::vector<float>&   mv = mask.get_data();
    const std::vector<uint8_t>& mm = mask.get_mask();
    const std::vector<float>&   sv = seed_mask.get_data();
    const std::vector<uint8_t>& sm = seed_mask.get_mask();

    std::vector<std::uint8_t> work(n, 0u);
    for (std::size_t i = 0; i < n; ++i)
        if (mm[i] && mv[i] > 0.5f) work[i] = 1u;

    std::vector<std::int64_t> seed_flats;
    for (std::size_t i = 0; i < n; ++i)
        if (sm[i] && sv[i] > 0.5f) seed_flats.push_back(static_cast<std::int64_t>(i));

    // Nothing to be reachable from. Returning the input unchanged is the honest
    // answer; emptying it would look like the filter had rejected everything.
    if (seed_flats.empty()) return mask;

    const lsdtt_raster_cpu dist =
        euclidean_distance_transform(seed_mask, options.n_threads);

    filter_components_by_seed_reachability_in_place(
        work, dist.get_data(), dist.get_mask(), seed_flats,
        options.dilation_pixels, options.require_seed_adjacent,
        options.max_detour_ratio, options.require_seed_line_of_sight,
        meta.cell_size, meta.n_rows, meta.n_cols);

    std::vector<float>   od(n, 0.0f);
    std::vector<uint8_t> om(n, 1u);
    for (std::size_t i = 0; i < n; ++i)
        if (work[i]) od[i] = 1.0f;

    return lsdtt_raster_cpu(std::move(od), std::move(om), meta);
}

// -----------------------------------------------------------------------------
// measure_mask_width_across_polyline (public - see header)
// -----------------------------------------------------------------------------
mask_width_result measure_mask_width_across_polyline(
    const lsdtt_raster_cpu&          mask,
    const std::vector<std::int64_t>& nodes,
    double section_spacing_metres,
    double max_half_length_metres,
    int    tangent_span_nodes)
{
    mask_width_result out;
    std::vector<double> envelopes;

    const lsdtt_raster_meta& meta = mask.get_meta();
    const int n_rows = meta.n_rows, n_cols = meta.n_cols;
    const double cell = meta.cell_size;
    const std::vector<float>& mv = mask.get_data();
    if (nodes.size() < 2 || !(cell > 0.0)) return out;

    const int last = static_cast<int>(nodes.size()) - 1;
    const int span = std::max(1, tangent_span_nodes);
    const int max_steps = std::max(
        1, static_cast<int>(max_half_length_metres / cell));

    // Stations every section_spacing_metres of along-line distance.
    double since = section_spacing_metres;   // force one at the start
    for (int k = 0; k <= last; ++k)
    {
        if (k > 0)
        {
            const std::int64_t a = nodes[static_cast<std::size_t>(k - 1)];
            const std::int64_t b = nodes[static_cast<std::size_t>(k)];
            const double dr = static_cast<double>(a / n_cols - b / n_cols);
            const double dc = static_cast<double>(a % n_cols - b % n_cols);
            since += cell * std::sqrt(dr * dr + dc * dc);
        }
        if (since < section_spacing_metres) continue;
        since = 0.0;

        // Local direction over +/- span nodes, then the perpendicular.
        const int ka = std::max(k - span, 0);
        const int kb = std::min(k + span, last);
        const std::int64_t fa = nodes[static_cast<std::size_t>(ka)];
        const std::int64_t fb = nodes[static_cast<std::size_t>(kb)];
        double dx = static_cast<double>((fb % n_cols) - (fa % n_cols));
        double dy = static_cast<double>((fb / n_cols) - (fa / n_cols));
        const double len = std::hypot(dx, dy);
        if (!(len > 0.0)) continue;
        dx /= len; dy /= len;
        const double px = -dy, py = dx;      // perpendicular

        const std::int64_t fk = nodes[static_cast<std::size_t>(k)];
        const double r0 = static_cast<double>(fk / n_cols);
        const double c0 = static_cast<double>(fk % n_cols);

        // Walk both ways, summing SET cells only. A section that runs to the
        // limit on both sides is measuring the limit, not the mask.
        std::int64_t n_set = 0;
        int far_set[2] = {0, 0};     // furthest set cell each way, in steps
        for (int side = -1; side <= 1; side += 2)
        {
            for (int s = 1; s <= max_steps; ++s)
            {
                const int rr = static_cast<int>(std::llround(r0 + side * py * s));
                const int cc = static_cast<int>(std::llround(c0 + side * px * s));
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) break;
                const std::size_t i =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (mv[i] > 0.5f)
                {
                    ++n_set;
                    far_set[(side < 0) ? 0 : 1] = s;
                }
            }
        }
        // A section still inside the mask at the search limit on BOTH sides is
        // measuring max_half_length, not the valley, so it is no measurement.
        // The test is on the furthest SET cell, not on how far the walk got:
        // walking the full half-length over ground that left the mask 100 m out
        // is exactly the ordinary case, and rejecting those leaves only the few
        // sections whose walk was cut short by the grid edge.
        if (far_set[0] >= max_steps && far_set[1] >= max_steps) continue;

        // The node's own cell counts once when it is inside the mask.
        const std::size_t fi = static_cast<std::size_t>(fk);
        if (fi < mv.size() && mv[fi] > 0.5f) ++n_set;

        // Both recorded only for sections that actually cross the mask; an
        // empty section is not a zero-width valley, it is no measurement at
        // all, and averaging it in drags the median toward the cell size.
        if (n_set > 0)
        {
            out.section_width_metres.push_back(
                static_cast<double>(n_set) * cell);
            envelopes.push_back(
                static_cast<double>(far_set[0] + far_set[1] + 1) * cell);
        }
    }

    out.n_sections = static_cast<std::int64_t>(out.section_width_metres.size());
    if (out.n_sections == 0) return out;

    std::vector<double> w = out.section_width_metres;
    std::sort(w.begin(), w.end());
    auto pct = [&](double p) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(
            std::clamp(p * static_cast<double>(w.size() - 1), 0.0,
                       static_cast<double>(w.size() - 1)));
        return w[idx];
    };
    out.median_width_metres = pct(0.5);
    out.p25_width_metres    = pct(0.25);
    out.p75_width_metres    = pct(0.75);

    if (!envelopes.empty())
    {
        std::sort(envelopes.begin(), envelopes.end());
        out.median_envelope_metres = envelopes[envelopes.size() / 2];
    }
    return out;
}

} // end namespace lsdtt