// lsdtt_raster_binary_mask_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html


#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include "lsdtt_raster_binary_mask_cpu.hpp"

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
// binary_directional_closing — internal helpers
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

    // Remove duplicates that arise from rounding (common at multiples of 45°).
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
                // Was background even after dilation — stays 0
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
    // Step 1: erode to get seeds — pixels that survive the radius filter.
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
    constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

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
        if (in_data[i] > 0.5f) return;   // foreground — not background
        reachable[i] = true;
        q.push(i);
    };

    for (int c = 0; c < ncols; ++c) { try_enqueue(0, c); try_enqueue(nrows-1, c); }
    for (int r = 0; r < nrows;  ++r) { try_enqueue(r, 0); try_enqueue(r, ncols-1); }

    constexpr int DR[4] = {-1, 1,  0, 0};
    constexpr int DC[4] = { 0, 0, -1, 1};

    // Also seed valid background pixels adjacent to any nodata pixel — they sit
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

    // Unreachable valid background pixels → foreground.
    // Nodata pixels (in_mask==0) get the nodata sentinel in the data array.
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask = in_mask;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (in_data[i] <= 0.5f && !reachable[i])
            out_data[i] = 1.0f;   // hole → foreground
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
        // Evenly spaced orientations in [0°, 180°).
        // A line at θ and θ+180° is the same, so we only need half the circle.
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
    out_meta.data_type           = GDT_UInt32;
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

} // end namespace lsdtt