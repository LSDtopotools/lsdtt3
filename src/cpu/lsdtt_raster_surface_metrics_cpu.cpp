// lsdtt_surface_metrics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"    // need this for the polyfit struct
#include "../cpu/lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"

#include <map>
#include <array>
#include <set>
#include <string>
#include <thread>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cctype>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// -----------------------------------------------------------------------------
// Internal constants
// -----------------------------------------------------------------------------
static constexpr double k_pi      = M_PI;
static constexpr double k_deg2rad = k_pi / 180.0;

// -----------------------------------------------------------------------------
// Internal helper — build lsdtt_array2d views from a raster base object.
//
// These are non-owning views: the array2d objects share memory with the
// raster's internal vectors. They are used only within this translation
// unit. The raster must outlive the views.
// -----------------------------------------------------------------------------
static lsdtt_array2d<float> data_view(const lsdtt_raster_cpu&raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    return lsdtt_array2d<float>(
        meta.n_rows, meta.n_cols,
        lsdtt_array2d<float>::from_ptr,
        raster.get_data().data()
    );
}

static lsdtt_array2d<uint8_t> mask_view(const lsdtt_raster_cpu&raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    return lsdtt_array2d<uint8_t>(
        meta.n_rows, meta.n_cols,
        lsdtt_array2d<uint8_t>::from_ptr,
        raster.get_mask().data()
    );
}

// -----------------------------------------------------------------------------
// Internal helper — check whether all pixels in a 3x3 window are valid.
// Used by algorithms that require a full 3x3 neighbourhood.
// -----------------------------------------------------------------------------
static bool window3x3_has_nodata(
    const lsdtt_array2d<float>&   dem,
    const lsdtt_array2d<uint8_t>& mask,
    const lsdtt_raster_meta&      meta,
    int row, int col
)
{
    for (int dr = -1; dr <= 1; ++dr)
    {
        for (int dc = -1; dc <= 1; ++dc)
        {
            if (mask(row + dr, col + dc) == 0) return true;
            if (meta.is_nodata(dem(row + dr, col + dc))) return true;
        }
    }
    return false;
}


// -----------------------------------------------------------------------------
// window_offsets
// -----------------------------------------------------------------------------
struct window_offsets
{
    int pixel_radius = 0;
    std::vector<std::pair<int,int>> offsets;
};

// -----------------------------------------------------------------------------
// make_window_offsets
//
// Convert a physical window radius in metres to a list of row/col offsets.
// For circular windows, offsets are included only if their cell-centre
// distance from the origin is <= radius_metres.
// -----------------------------------------------------------------------------
static window_offsets make_window_offsets(
    double       radius_metres,
    double       cell_size,
    window_shape shape
)
{
    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt surface metrics: window radius must be finite and > 0 metres"
        );
    }

    if (!std::isfinite(cell_size) || cell_size <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt surface metrics: cell_size must be finite and > 0"
        );
    }

    const int pixel_radius = std::max(
        1,
        static_cast<int>(std::ceil(radius_metres / cell_size))
    );

    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(
        static_cast<std::size_t>(2 * pixel_radius + 1) *
        static_cast<std::size_t>(2 * pixel_radius + 1)
    );

    for (int dr = -pixel_radius; dr <= pixel_radius; ++dr)
    {
        for (int dc = -pixel_radius; dc <= pixel_radius; ++dc)
        {
            if (shape == window_shape::square)
            {
                offsets.emplace_back(dr, dc);
            }
            else
            {
                const double dx = static_cast<double>(dc) * cell_size;
                const double dy = static_cast<double>(dr) * cell_size;
                const double dist = std::sqrt(dx * dx + dy * dy);

                if (dist <= radius_metres + 1e-12)
                {
                    offsets.emplace_back(dr, dc);
                }
            }
        }
    }

    if (offsets.empty())
    {
        offsets.emplace_back(0, 0);
    }

    return window_offsets{pixel_radius, std::move(offsets)};
}

static inline std::size_t linear_index(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
         + static_cast<std::size_t>(col);
}

static int resolve_windowed_metric_threads(int n_threads)
{
    if (n_threads > 0) return n_threads;
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    const unsigned hc = std::thread::hardware_concurrency();
    return (hc == 0u) ? 1 : static_cast<int>(hc);
#endif
}

static double median_in_place(std::vector<float>& values)
{
    const std::size_t n = values.size();
    const std::size_t mid = n / 2;

    std::nth_element(
        values.begin(),
        values.begin() + static_cast<std::ptrdiff_t>(mid),
        values.end()
    );
    const double upper = static_cast<double>(values[mid]);

    if ((n % 2) == 1)
    {
        return upper;
    }

    std::nth_element(
        values.begin(),
        values.begin() + static_cast<std::ptrdiff_t>(mid - 1),
        values.begin() + static_cast<std::ptrdiff_t>(mid)
    );
    const double lower = static_cast<double>(values[mid - 1]);

    return 0.5 * (lower + upper);
}



// -----------------------------------------------------------------------------
// compute_hillshade
// -----------------------------------------------------------------------------

lsdtt_raster_cpu compute_hillshade(
    const lsdtt_raster_cpu&raster,
    double azimuth,
    double altitude,
    double z_factor
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_hillshade: input raster metadata is invalid"
        );
    }

    if (z_factor <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_hillshade: z_factor must be > 0"
        );
    }

    const int    nrows   = meta.n_rows;
    const int    ncols   = meta.n_cols;
    const double L       = meta.cell_size;
    const float  nodata  = static_cast<float>(meta.no_data_value);

    // Build 2D views of the input data and mask
    lsdtt_array2d<float>   dem  = data_view(raster);
    lsdtt_array2d<uint8_t> mask = mask_view(raster);

    // Precompute sun vector components
    const double alt_rad = altitude * k_deg2rad;
    const double az_rad  = azimuth  * k_deg2rad;
    const double sin_alt = std::sin(alt_rad);
    const double cos_alt = std::cos(alt_rad);

    // Output arrays — initialised to nodata / invalid
    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

    // Interior pixels only — edge pixels remain nodata because the
    // 3x3 window cannot be fully populated there.
#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            if (window3x3_has_nodata(dem, mask, meta, row, col)) continue;

            // Cardinal neighbours (Zevenbergen and Thorne eq. 13)
            // north = row-1, south = row+1, east = col+1, west = col-1
            const double z_north = static_cast<double>(dem(row - 1, col)) * z_factor;
            const double z_south = static_cast<double>(dem(row + 1, col)) * z_factor;
            const double z_east  = static_cast<double>(dem(row, col + 1)) * z_factor;
            const double z_west  = static_cast<double>(dem(row, col - 1)) * z_factor;

            // First derivatives
            const double p = (z_east  - z_west)  / (2.0 * L);   // dz/dx
            const double q = (z_north - z_south) / (2.0 * L);   // dz/dy

            // Slope (radians)
            const double slope_rad = std::atan(std::sqrt(p * p + q * q));

            // Aspect (radians, clockwise from north)
            double aspect_rad;
            if (std::abs(p) < 1e-10 && std::abs(q) < 1e-10)
            {
                aspect_rad = 0.0;  // flat pixel — aspect undefined
            }
            else
            {
                aspect_rad = std::atan2(p, q);
                if (aspect_rad < 0.0) aspect_rad += 2.0 * k_pi;
            }

            // Lambertian hillshade
            double hs = 1.0 + 254.0 * (
                  sin_alt * std::cos(slope_rad)
                + cos_alt * std::sin(slope_rad)
                          * std::cos(az_rad - (k_pi / 2.0) - aspect_rad)
            );

            // Clamp to [0, 255]
            hs = std::max(0.0, std::min(255.0, hs));

            out(row, col)      = static_cast<float>(hs);
            out_mask(row, col) = 1;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// moving_window_mean
// -----------------------------------------------------------------------------

lsdtt_raster_cpu moving_window_mean(
    const lsdtt_raster_cpu&raster,
    double       radius_metres,
    window_shape shape
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_mean: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const window_offsets win = make_window_offsets(
        radius_metres,
        meta.cell_size,
        shape
    );

    // Build 2D views
    lsdtt_array2d<float>   dem  = data_view(raster);
    lsdtt_array2d<uint8_t> mask = mask_view(raster);

    // Output arrays
    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            double sum   = 0.0;
            int    count = 0;

            for (const auto& [dr, dc] : win.offsets)
            {
                const int wr = row + dr;
                const int wc = col + dc;

                if (!dem.in_bounds(wr, wc)) continue;
                if (mask(wr, wc) == 0) continue;
                if (meta.is_nodata(dem(wr, wc))) continue;

                sum += static_cast<double>(dem(wr, wc));
                ++count;
            }

            if (count > 0)
            {
                out(row, col)      = static_cast<float>(sum / count);
                out_mask(row, col) = 1;
            }
            // else: remains nodata / mask = 0
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}




// -----------------------------------------------------------------------------
// moving_window_local_relief
// -----------------------------------------------------------------------------

lsdtt_raster_cpu moving_window_local_relief(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres,
    window_shape            shape
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_local_relief: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const window_offsets win = make_window_offsets(
        radius_metres,
        meta.cell_size,
        shape
    );

    lsdtt_array2d<float>   dem  = data_view(raster);
    lsdtt_array2d<uint8_t> mask = mask_view(raster);

    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            double local_min =  std::numeric_limits<double>::infinity();
            double local_max = -std::numeric_limits<double>::infinity();
            int    count     = 0;

            for (const auto& [dr, dc] : win.offsets)
            {
                const int wr = row + dr;
                const int wc = col + dc;

                if (!dem.in_bounds(wr, wc)) continue;
                if (mask(wr, wc) == 0) continue;
                if (meta.is_nodata(dem(wr, wc))) continue;

                const double z = static_cast<double>(dem(wr, wc));
                if (z < local_min) local_min = z;
                if (z > local_max) local_max = z;
                ++count;
            }

            if (count > 0)
            {
                out(row, col)      = static_cast<float>(local_max - local_min);
                out_mask(row, col) = 1;
            }
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}




// -----------------------------------------------------------------------------
// calculate_madk2
//
// Computes the MADk2 geostatistical roughness metric following
// Trevisani, Teza and Guth (2023) using second-order directional differences
// within a moving circular window.
// -----------------------------------------------------------------------------
madk2_result calculate_madk2(
    const lsdtt_raster_cpu& raster,
    double                  window_radius_metres,
    int                     lag_pixels,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::calculate_madk2: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    if (!std::isfinite(window_radius_metres) || window_radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::calculate_madk2: window_radius_metres must be finite and > 0"
        );
    }

    if (lag_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::calculate_madk2: lag_pixels must be >= 1"
        );
    }

    const double effective_radius_metres =
        std::max(window_radius_metres, meta.cell_size);

    const window_offsets win = make_window_offsets(
        effective_radius_metres,
        meta.cell_size,
        window_shape::circular
    );

    const int pixel_radius = win.pixel_radius;
    const int nrows        = meta.n_rows;
    const int ncols        = meta.n_cols;
    const int threads      = resolve_windowed_metric_threads(n_threads);
    const float nodata     = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<float>   or_data (npix, nodata);
    std::vector<uint8_t> or_mask (npix, uint8_t{0});
    std::vector<float>   dmc_data(npix, nodata);
    std::vector<uint8_t> dmc_mask(npix, uint8_t{0});
    std::vector<float>   ra_data (npix, nodata);
    std::vector<uint8_t> ra_mask (npix, uint8_t{0});

    struct madk2_direction
    {
        int    dr;
        int    dc;
        double axis_angle_rad;
        double dist_sq;
    };

    const int lag_sq = lag_pixels * lag_pixels;
    const std::array<madk2_direction, 4> dirs =
    {{
        { 0,          lag_pixels, 0.0,              static_cast<double>(lag_sq)       },
        { lag_pixels, 0,          k_pi / 2.0,       static_cast<double>(lag_sq)       },
        { lag_pixels, lag_pixels, k_pi / 4.0,       2.0 * static_cast<double>(lag_sq) },
        {-lag_pixels, lag_pixels, 3.0 * k_pi / 4.0, 2.0 * static_cast<double>(lag_sq) }
    }};

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic) num_threads(threads)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = linear_index(row, col, ncols);

            if (mask[idx] == 0 || meta.is_nodata(data[idx]) || !std::isfinite(data[idx]))
            {
                continue;
            }

            if (row - pixel_radius < 0 || row + pixel_radius + 1 > nrows ||
                col - pixel_radius < 0 || col + pixel_radius + 1 > ncols)
            {
                continue;
            }

            bool edge_lag = false;
            for (const auto& dir : dirs)
            {
                const int ri_neg = row - dir.dr;
                const int rj_neg = col - dir.dc;
                const int ri_pos = row + dir.dr;
                const int rj_pos = col + dir.dc;

                if (ri_neg < 0 || ri_neg >= nrows || rj_neg < 0 || rj_neg >= ncols ||
                    ri_pos < 0 || ri_pos >= nrows || rj_pos < 0 || rj_pos >= ncols)
                {
                    edge_lag = true;
                    break;
                }
            }
            if (edge_lag)
            {
                continue;
            }

            std::array<double, 4> mad_dir{{0.0, 0.0, 0.0, 0.0}};
            bool valid = true;

            for (std::size_t d = 0; d < dirs.size(); ++d)
            {
                std::vector<float> k2_diffs;
                k2_diffs.reserve(win.offsets.size());

                for (const auto& [dr_win, dc_win] : win.offsets)
                {
                    const int ci = row + dr_win;
                    const int cj = col + dc_win;

                    const int ri_pos = ci + dirs[d].dr;
                    const int rj_pos = cj + dirs[d].dc;
                    const int ri_neg = ci - dirs[d].dr;
                    const int rj_neg = cj - dirs[d].dc;

                    if (ri_pos < 0 || ri_pos >= nrows || rj_pos < 0 || rj_pos >= ncols ||
                        ri_neg < 0 || ri_neg >= nrows || rj_neg < 0 || rj_neg >= ncols)
                    {
                        continue;
                    }

                    const std::size_t idx_pos = linear_index(ri_pos, rj_pos, ncols);
                    const std::size_t idx_cen = linear_index(ci,     cj,     ncols);
                    const std::size_t idx_neg = linear_index(ri_neg, rj_neg, ncols);

                    if (mask[idx_pos] == 0 || mask[idx_cen] == 0 || mask[idx_neg] == 0)
                    {
                        continue;
                    }

                    const float z_pos = data[idx_pos];
                    const float z_cen = data[idx_cen];
                    const float z_neg = data[idx_neg];

                    if (meta.is_nodata(z_pos) || meta.is_nodata(z_cen) || meta.is_nodata(z_neg) ||
                        !std::isfinite(z_pos) || !std::isfinite(z_cen) || !std::isfinite(z_neg))
                    {
                        continue;
                    }

                    const double k2 =
                        static_cast<double>(z_neg) -
                        2.0 * static_cast<double>(z_cen) +
                        static_cast<double>(z_pos);

                    const double k2_norm = std::abs(k2) / dirs[d].dist_sq;
                    k2_diffs.push_back(static_cast<float>(k2_norm));
                }

                if (k2_diffs.size() < 3)
                {
                    valid = false;
                    break;
                }

                mad_dir[d] = median_in_place(k2_diffs);
            }

            if (!valid)
            {
                continue;
            }

            double or_val = 0.0;
            for (double v : mad_dir)
            {
                or_val += v;
            }
            or_val /= static_cast<double>(mad_dir.size());

            double rx = 0.0;
            double ry = 0.0;
            double mad_sum = 0.0;

            for (std::size_t d = 0; d < dirs.size(); ++d)
            {
                rx      += mad_dir[d] * std::cos(dirs[d].axis_angle_rad);
                ry      += mad_dir[d] * std::sin(dirs[d].axis_angle_rad);
                mad_sum += mad_dir[d];
            }

            const double r_mag  = std::sqrt(rx * rx + ry * ry);
            const double ra_val = (mad_sum > 0.0) ? (r_mag / mad_sum) : 0.0;

            const double max_rough_angle_rad = std::atan2(ry, rx);
            double dmc_deg = (max_rough_angle_rad + k_pi / 2.0) * (180.0 / k_pi);

            while (dmc_deg < 0.0)
            {
                dmc_deg += 180.0;
            }
            while (dmc_deg >= 180.0)
            {
                dmc_deg -= 180.0;
            }

            or_data[idx]  = static_cast<float>(or_val);
            or_mask[idx]  = uint8_t{1};
            dmc_data[idx] = static_cast<float>(dmc_deg);
            dmc_mask[idx] = uint8_t{1};
            ra_data[idx]  = static_cast<float>(ra_val);
            ra_mask[idx]  = uint8_t{1};
        }
    }

    return madk2_result
    {
        lsdtt_raster_cpu(std::move(or_data),  std::move(or_mask),  meta),
        lsdtt_raster_cpu(std::move(dmc_data), std::move(dmc_mask), meta),
        lsdtt_raster_cpu(std::move(ra_data),  std::move(ra_mask),  meta)
    };
}







// -----------------------------------------------------------------------------
// Polyfit-derived surface metrics
//
// Derivative notation used throughout this translation unit:
//
//   p = z_x
//   q = z_y
//   r = z_xx
//   s = z_xy
//   t = z_yy
//
// For the fitted quadratic
//
//   z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
//
// evaluated at the pixel centre:
//
//   p = d
//   q = e
//   r = 2*a
//   s = c
//   t = 2*b
//
// Coordinate convention:
//   x  — positive east
//   y  — positive north
//
// The derivative naming convention follows:
//
//   Florinsky, I. V. (2017). An illustrated introduction to general
//   geomorphometry. Progress in Physical Geography: Earth and Environment,
//   41(6), 723-752. doi:10.1177/0309133317733667
//
// Original polyfit implementation:
//
//   Hurst, M. D., S. M. Mudd, R. Walcott, M. Attal, and K. Yoo (2012),
//   Using hilltop curvature to derive the spatial distribution of erosion
//   rates, Journal of Geophysical Research, 117, F02017,
//   doi:10.1029/2011JF002057
//
// Revised in lsdtt3 to use wavelet fitting for computational efficiency.
// -----------------------------------------------------------------------------

namespace
{

static constexpr double k_rad2deg = 180.0 / k_pi;

// -----------------------------------------------------------------------------
// flat_index
// -----------------------------------------------------------------------------
static inline std::size_t flat_index(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
         + static_cast<std::size_t>(col);
}

// -----------------------------------------------------------------------------
// resolve_thread_count_surface_metrics
// -----------------------------------------------------------------------------
static int resolve_thread_count_surface_metrics(int n_threads)
{
    int threads = n_threads;

    if (threads <= 0)
    {
#ifdef _OPENMP
        threads = omp_get_max_threads();
#else
        threads = static_cast<int>(std::thread::hardware_concurrency());
#endif
    }

    if (threads < 1) threads = 1;
    return threads;
}

// -----------------------------------------------------------------------------
// assert_same_meta
// -----------------------------------------------------------------------------
static void assert_same_meta(
    const lsdtt_raster_meta& a,
    const lsdtt_raster_meta& b,
    const std::string&       a_name,
    const std::string&       b_name)
{
    if (a.n_rows != b.n_rows ||
        a.n_cols != b.n_cols ||
        a.cell_size != b.cell_size ||
        a.x_min != b.x_min ||
        a.y_min != b.y_min ||
        a.no_data_value != b.no_data_value)
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: incompatible coefficient rasters: "
            + a_name + " and " + b_name
        );
    }
}

// -----------------------------------------------------------------------------
// assert_polyfit_compatible
// -----------------------------------------------------------------------------
static void assert_polyfit_compatible(const poly_fit_coeffs& polycoeff)
{
    const lsdtt_raster_meta& meta = polycoeff.a.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: invalid coefficient raster metadata"
        );
    }

    meta.check_projected_metres();

    assert_same_meta(meta, polycoeff.b.get_meta(), "a", "b");
    assert_same_meta(meta, polycoeff.c.get_meta(), "a", "c");
    assert_same_meta(meta, polycoeff.d.get_meta(), "a", "d");
    assert_same_meta(meta, polycoeff.e.get_meta(), "a", "e");
    assert_same_meta(meta, polycoeff.f.get_meta(), "a", "f");

    const std::size_t n_pixels =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    if (polycoeff.a.get_data().size() != n_pixels || polycoeff.a.get_mask().size() != n_pixels ||
        polycoeff.b.get_data().size() != n_pixels || polycoeff.b.get_mask().size() != n_pixels ||
        polycoeff.c.get_data().size() != n_pixels || polycoeff.c.get_mask().size() != n_pixels ||
        polycoeff.d.get_data().size() != n_pixels || polycoeff.d.get_mask().size() != n_pixels ||
        polycoeff.e.get_data().size() != n_pixels || polycoeff.e.get_mask().size() != n_pixels ||
        polycoeff.f.get_data().size() != n_pixels || polycoeff.f.get_mask().size() != n_pixels)
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: coefficient raster data/mask size mismatch"
        );
    }
}

// -----------------------------------------------------------------------------
// output_slot
// -----------------------------------------------------------------------------
struct output_slot
{
    std::vector<float>*   data = nullptr;
    std::vector<uint8_t>* mask = nullptr;

    bool requested() const
    {
        return data != nullptr && mask != nullptr;
    }

    void write(std::size_t idx, float value)
    {
        (*data)[idx] = value;
        (*mask)[idx] = uint8_t{1};
    }
};

// -----------------------------------------------------------------------------
// output_arrays
// -----------------------------------------------------------------------------
struct output_arrays
{
    std::vector<float>   data;
    std::vector<uint8_t> mask;
};

} // unnamed namespace

// -----------------------------------------------------------------------------
// compute_polyfit_metrics
// -----------------------------------------------------------------------------
polyfit_metric_result compute_polyfit_metrics(
    const poly_fit_coeffs&        polycoeff,
    const polyfit_metric_request& request)
{
    assert_polyfit_compatible(polycoeff);

    if (!std::isfinite(request.z_factor) || request.z_factor <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: z_factor must be finite and > 0"
        );
    }

    if (!std::isfinite(request.flat_slope_tolerance) ||
        request.flat_slope_tolerance <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: flat_slope_tolerance must be finite and > 0"
        );
    }

    std::set<polyfit_metric> requested(
        request.metrics.begin(),
        request.metrics.end()
    );

    if (requested.empty())
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: no metrics were requested"
        );
    }

    const lsdtt_raster_meta& meta   = polycoeff.a.get_meta();
    const int                nrows  = meta.n_rows;
    const int                ncols  = meta.n_cols;
    const std::size_t        npix   =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const float              nodata = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   a_data = polycoeff.a.get_data();
    const std::vector<float>&   b_data = polycoeff.b.get_data();
    const std::vector<float>&   c_data = polycoeff.c.get_data();
    const std::vector<float>&   d_data = polycoeff.d.get_data();
    const std::vector<float>&   e_data = polycoeff.e.get_data();
    const std::vector<float>&   f_data = polycoeff.f.get_data();

    const std::vector<uint8_t>& a_mask = polycoeff.a.get_mask();
    const std::vector<uint8_t>& b_mask = polycoeff.b.get_mask();
    const std::vector<uint8_t>& c_mask = polycoeff.c.get_mask();
    const std::vector<uint8_t>& d_mask = polycoeff.d.get_mask();
    const std::vector<uint8_t>& e_mask = polycoeff.e.get_mask();
    const std::vector<uint8_t>& f_mask = polycoeff.f.get_mask();

    std::map<polyfit_metric, output_arrays> outputs;
    for (polyfit_metric metric : requested)
    {
        outputs[metric] = output_arrays{
            std::vector<float>(npix, nodata),
            std::vector<uint8_t>(npix, uint8_t{0})
        };
    }

    auto make_slot = [&](polyfit_metric metric) -> output_slot
    {
        auto it = outputs.find(metric);
        if (it == outputs.end()) return output_slot{};
        return output_slot{ &it->second.data, &it->second.mask };
    };

    output_slot out_fitted_surface              = make_slot(polyfit_metric::fitted_surface);
    output_slot out_dzdx                        = make_slot(polyfit_metric::dzdx);
    output_slot out_dzdy                        = make_slot(polyfit_metric::dzdy);
    output_slot out_d2zdx2                      = make_slot(polyfit_metric::d2zdx2);
    output_slot out_d2zdxdy                     = make_slot(polyfit_metric::d2zdxdy);
    output_slot out_d2zdy2                      = make_slot(polyfit_metric::d2zdy2);
    output_slot out_slope                       = make_slot(polyfit_metric::slope);
    output_slot out_aspect_radians              = make_slot(polyfit_metric::aspect_radians);
    output_slot out_aspect_degrees              = make_slot(polyfit_metric::aspect_degrees);
    output_slot out_laplacian                   = make_slot(polyfit_metric::laplacian);
    output_slot out_hessian_determinant         = make_slot(polyfit_metric::hessian_determinant);
    output_slot out_mean_curvature              = make_slot(polyfit_metric::mean_curvature);
    output_slot out_gaussian_curvature          = make_slot(polyfit_metric::gaussian_curvature);
    output_slot out_maximum_curvature           = make_slot(polyfit_metric::maximum_curvature);
    output_slot out_minimum_curvature           = make_slot(polyfit_metric::minimum_curvature);
    output_slot out_curvedness                  = make_slot(polyfit_metric::curvedness);
    output_slot out_shape_index                 = make_slot(polyfit_metric::shape_index);
    output_slot out_profile_curvature           = make_slot(polyfit_metric::profile_curvature);
    output_slot out_plan_curvature              = make_slot(polyfit_metric::plan_curvature);
    output_slot out_tangential_curvature        = make_slot(polyfit_metric::tangential_curvature);
    output_slot out_difference_curvature        = make_slot(polyfit_metric::difference_curvature);
    output_slot out_horizontal_excess_curvature = make_slot(polyfit_metric::horizontal_excess_curvature);
    output_slot out_vertical_excess_curvature   = make_slot(polyfit_metric::vertical_excess_curvature);
    output_slot out_accumulation_curvature      = make_slot(polyfit_metric::accumulation_curvature);
    output_slot out_ring_curvature              = make_slot(polyfit_metric::ring_curvature);
    output_slot out_rotor                       = make_slot(polyfit_metric::rotor);

    const bool need_first_derivs =
        out_dzdx.requested() ||
        out_dzdy.requested() ||
        out_slope.requested() ||
        out_aspect_radians.requested() ||
        out_aspect_degrees.requested() ||
        out_mean_curvature.requested() ||
        out_gaussian_curvature.requested() ||
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_profile_curvature.requested() ||
        out_plan_curvature.requested() ||
        out_tangential_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested() ||
        out_rotor.requested();

    const bool need_second_derivs =
        out_d2zdx2.requested() ||
        out_d2zdxdy.requested() ||
        out_d2zdy2.requested() ||
        out_laplacian.requested() ||
        out_hessian_determinant.requested() ||
        out_mean_curvature.requested() ||
        out_gaussian_curvature.requested() ||
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_profile_curvature.requested() ||
        out_plan_curvature.requested() ||
        out_tangential_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested() ||
        out_rotor.requested();

    const bool need_gradient =
        out_slope.requested() ||
        out_aspect_radians.requested() ||
        out_aspect_degrees.requested() ||
        out_mean_curvature.requested() ||
        out_gaussian_curvature.requested() ||
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_profile_curvature.requested() ||
        out_plan_curvature.requested() ||
        out_tangential_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested() ||
        out_rotor.requested();

    const bool need_first_fundamental =
        out_mean_curvature.requested() ||
        out_gaussian_curvature.requested() ||
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_profile_curvature.requested() ||
        out_tangential_curvature.requested();

    const bool need_principal_curvatures =
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_ring_curvature.requested();

    const bool need_vertical_curvature =
        out_profile_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested();

    const bool need_tangential_curvature =
        out_tangential_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested();

    const bool need_plan_curvature =
        out_plan_curvature.requested();

    const bool need_rotor =
        out_rotor.requested();

    const double z_factor         = request.z_factor;
    const double flat_tol         = request.flat_slope_tolerance;
    const double flat_tol_squared = flat_tol * flat_tol;
    const int    threads          = resolve_thread_count_surface_metrics(request.n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = flat_index(row, col, ncols);

            const bool coeff_valid =
                a_mask[idx] != 0 && b_mask[idx] != 0 && c_mask[idx] != 0 &&
                d_mask[idx] != 0 && e_mask[idx] != 0 && f_mask[idx] != 0 &&
                !meta.is_nodata(a_data[idx]) &&
                !meta.is_nodata(b_data[idx]) &&
                !meta.is_nodata(c_data[idx]) &&
                !meta.is_nodata(d_data[idx]) &&
                !meta.is_nodata(e_data[idx]) &&
                !meta.is_nodata(f_data[idx]) &&
                std::isfinite(a_data[idx]) &&
                std::isfinite(b_data[idx]) &&
                std::isfinite(c_data[idx]) &&
                std::isfinite(d_data[idx]) &&
                std::isfinite(e_data[idx]) &&
                std::isfinite(f_data[idx]);

            if (!coeff_valid) continue;

            const double a = z_factor * static_cast<double>(a_data[idx]);
            const double b = z_factor * static_cast<double>(b_data[idx]);
            const double c = z_factor * static_cast<double>(c_data[idx]);
            const double d = z_factor * static_cast<double>(d_data[idx]);
            const double e = z_factor * static_cast<double>(e_data[idx]);
            const double f = z_factor * static_cast<double>(f_data[idx]);

            double p = 0.0;
            double q = 0.0;
            double r = 0.0;
            double s = 0.0;
            double t = 0.0;

            if (need_first_derivs)
            {
                p = d;
                q = e;
            }

            if (need_second_derivs)
            {
                r = 2.0 * a;
                s = c;
                t = 2.0 * b;
            }

            double p2 = 0.0;
            double q2 = 0.0;
            double g2 = 0.0;
            double g  = 0.0;

            if (need_gradient)
            {
                p2 = p * p;
                q2 = q * q;
                g2 = p2 + q2;
                g  = std::sqrt(g2);
            }

            double w2 = 0.0;
            double w  = 0.0;
            double w3 = 0.0;
            double w4 = 0.0;

            if (need_first_fundamental)
            {
                w2 = 1.0 + g2;
                w  = std::sqrt(w2);
                w3 = w2 * w;
                w4 = w2 * w2;
            }

            double H = std::numeric_limits<double>::quiet_NaN();
            double K = std::numeric_limits<double>::quiet_NaN();

            if (out_mean_curvature.requested() ||
                out_gaussian_curvature.requested() ||
                need_principal_curvatures)
            {
                // mean_curvature =
                // ((1 + q^2)*r - 2*p*q*s + (1 + p^2)*t) /
                // (2*(1 + p^2 + q^2)^(3/2))
                if (w3 > 0.0)
                {
                    H = ((1.0 + q2) * r - 2.0 * p * q * s + (1.0 + p2) * t)
                      / (2.0 * w3);
                }

                // gaussian_curvature =
                // (r*t - s^2) / (1 + p^2 + q^2)^2
                if (w4 > 0.0)
                {
                    K = (r * t - s * s) / w4;
                }
            }

            const double unset = std::numeric_limits<double>::quiet_NaN();

            double principal_numeric_max = unset;   // k_max
            double principal_numeric_min = unset;   // k_min

            double plan_curvature_value               = unset;
            double tangential_curvature_value         = unset;
            double vertical_curvature_value           = unset;
            double difference_curvature_value         = unset;
            double horizontal_excess_curvature_value  = unset;
            double vertical_excess_curvature_value    = unset;
            double accumulation_curvature_value       = unset;
            double ring_curvature_value               = unset;
            double rotor_value                        = unset;

            if (need_principal_curvatures)
            {
                const double disc = std::max(0.0, H * H - K);
                const double root = std::sqrt(disc);

                // k_max = H + sqrt(H^2 - K)
                // k_min = H - sqrt(H^2 - K)
                principal_numeric_max = H + root;
                principal_numeric_min = H - root;
            }

            // fitted_surface = f
            if (out_fitted_surface.requested())
            {
                out_fitted_surface.write(idx, static_cast<float>(f));
            }

            // dzdx = p
            if (out_dzdx.requested())
            {
                out_dzdx.write(idx, static_cast<float>(p));
            }

            // dzdy = q
            if (out_dzdy.requested())
            {
                out_dzdy.write(idx, static_cast<float>(q));
            }

            // d2zdx2 = r
            if (out_d2zdx2.requested())
            {
                out_d2zdx2.write(idx, static_cast<float>(r));
            }

            // d2zdxdy = s
            if (out_d2zdxdy.requested())
            {
                out_d2zdxdy.write(idx, static_cast<float>(s));
            }

            // d2zdy2 = t
            if (out_d2zdy2.requested())
            {
                out_d2zdy2.write(idx, static_cast<float>(t));
            }

            // slope = sqrt(p^2 + q^2)
            if (out_slope.requested())
            {
                out_slope.write(idx, static_cast<float>(g));
            }

            // aspect_radians = atan2(p, q), clockwise from north on [0, 2*pi)
            if (out_aspect_radians.requested())
            {
                if (g2 > flat_tol_squared)
                {
                    double aspect = std::atan2(p, q);
                    if (aspect < 0.0) aspect += 2.0 * k_pi;
                    out_aspect_radians.write(idx, static_cast<float>(aspect));
                }
            }

            // aspect_degrees = aspect_radians * 180 / pi
            if (out_aspect_degrees.requested())
            {
                if (g2 > flat_tol_squared)
                {
                    double aspect = std::atan2(p, q);
                    if (aspect < 0.0) aspect += 2.0 * k_pi;
                    out_aspect_degrees.write(
                        idx,
                        static_cast<float>(aspect * k_rad2deg)
                    );
                }
            }

            // laplacian = r + t
            if (out_laplacian.requested())
            {
                out_laplacian.write(idx, static_cast<float>(r + t));
            }

            // hessian_determinant = r*t - s^2
            if (out_hessian_determinant.requested())
            {
                out_hessian_determinant.write(
                    idx,
                    static_cast<float>(r * t - s * s)
                );
            }

            // mean_curvature =
            // ((1 + q^2)*r - 2*p*q*s + (1 + p^2)*t) /
            // (2*(1 + p^2 + q^2)^(3/2))
            if (out_mean_curvature.requested())
            {
                if (std::isfinite(H))
                {
                    out_mean_curvature.write(idx, static_cast<float>(H));
                }
            }

            // gaussian_curvature =
            // (r*t - s^2) / (1 + p^2 + q^2)^2
            if (out_gaussian_curvature.requested())
            {
                if (std::isfinite(K))
                {
                    out_gaussian_curvature.write(idx, static_cast<float>(K));
                }
            }

            if (need_plan_curvature)
            {
                // plan_curvature =
                // - (q^2*r - 2*p*q*s + p^2*t) / sqrt((p^2 + q^2)^3)
                //
                // Equivalent denominator used in code:
                // (p^2 + q^2) * sqrt(p^2 + q^2)
                const double denom = g2 * std::sqrt(g2);

                if (denom > 0.0 && g2 > flat_tol_squared)
                {
                    plan_curvature_value =
                        - (q2 * r - 2.0 * p * q * s + p2 * t) / denom;
                }
            }

            if (out_plan_curvature.requested())
            {
                if (std::isfinite(plan_curvature_value))
                {
                    out_plan_curvature.write(
                        idx,
                        static_cast<float>(plan_curvature_value)
                    );
                }
            }

            if (need_tangential_curvature)
            {
                // tangential_curvature =
                // - (q^2*r - 2*p*q*s + p^2*t) /
                //   ((p^2 + q^2) * sqrt(1 + p^2 + q^2))
                const double denom = g2 * w;

                if (denom > 0.0 && g2 > flat_tol_squared)
                {
                    tangential_curvature_value =
                        - (q2 * r - 2.0 * p * q * s + p2 * t) / denom;
                }
            }

            if (out_tangential_curvature.requested())
            {
                if (std::isfinite(tangential_curvature_value))
                {
                    out_tangential_curvature.write(
                        idx,
                        static_cast<float>(tangential_curvature_value)
                    );
                }
            }

            if (need_vertical_curvature)
            {
                // profile_curvature = vertical_curvature =
                // - (p^2*r + 2*p*q*s + q^2*t) /
                //   ((p^2 + q^2) * (1 + p^2 + q^2)^(3/2))
                const double denom = g2 * w3;

                if (denom > 0.0 && g2 > flat_tol_squared)
                {
                    vertical_curvature_value =
                        - (p2 * r + 2.0 * p * q * s + q2 * t) / denom;
                }
            }

            if (out_profile_curvature.requested())
            {
                if (std::isfinite(vertical_curvature_value))
                {
                    out_profile_curvature.write(
                        idx,
                        static_cast<float>(vertical_curvature_value)
                    );
                }
            }

            if (need_rotor)
            {
                // rotor =
                // ((p^2 - q^2)*s - p*q*(r - t)) / sqrt((p^2 + q^2)^3)
                //
                // Equivalent denominator used in code:
                // (p^2 + q^2) * sqrt(p^2 + q^2)
                const double denom = g2 * std::sqrt(g2);

                if (denom > 0.0 && g2 > flat_tol_squared)
                {
                    rotor_value =
                        ((p2 - q2) * s - p * q * (r - t)) / denom;
                }
            }

            if (out_rotor.requested())
            {
                if (std::isfinite(rotor_value))
                {
                    out_rotor.write(
                        idx,
                        static_cast<float>(rotor_value)
                    );
                }
            }

            // maximum_curvature = k_max = H + sqrt(H^2 - K)
            if (out_maximum_curvature.requested())
            {
                if (std::isfinite(principal_numeric_max))
                {
                    out_maximum_curvature.write(
                        idx,
                        static_cast<float>(principal_numeric_max)
                    );
                }
            }

            // minimum_curvature = k_min = H - sqrt(H^2 - K)
            if (out_minimum_curvature.requested())
            {
                if (std::isfinite(principal_numeric_min))
                {
                    out_minimum_curvature.write(
                        idx,
                        static_cast<float>(principal_numeric_min)
                    );
                }
            }

            // curvedness = sqrt((k_max^2 + k_min^2) / 2)
            if (out_curvedness.requested())
            {
                if (std::isfinite(principal_numeric_max) &&
                    std::isfinite(principal_numeric_min))
                {
                    const double curvedness_value =
                        std::sqrt(
                            0.5 * (
                                principal_numeric_max * principal_numeric_max +
                                principal_numeric_min * principal_numeric_min
                            )
                        );

                    if (std::isfinite(curvedness_value))
                    {
                        out_curvedness.write(
                            idx,
                            static_cast<float>(curvedness_value)
                        );
                    }
                }
            }

            // shape_index =
            // (2 / pi) * atan((k_min + k_max) / (k_min - k_max))
            if (out_shape_index.requested())
            {
                if (std::isfinite(principal_numeric_max) &&
                    std::isfinite(principal_numeric_min))
                {
                    const double denom =
                        principal_numeric_min - principal_numeric_max;

                    if (std::abs(denom) > 1e-15)
                    {
                        const double shape_index_value =
                            (2.0 / k_pi) * std::atan(
                                (principal_numeric_min + principal_numeric_max) /
                                denom
                            );

                        if (std::isfinite(shape_index_value))
                        {
                            out_shape_index.write(
                                idx,
                                static_cast<float>(shape_index_value)
                            );
                        }
                    }
                }
            }

            if (std::isfinite(vertical_curvature_value) &&
                std::isfinite(tangential_curvature_value))
            {
                // difference_curvature = 0.5 * (k_v - k_h)
                difference_curvature_value =
                    0.5 * (
                        vertical_curvature_value -
                        tangential_curvature_value
                    );
            }

            if (out_difference_curvature.requested())
            {
                if (std::isfinite(difference_curvature_value))
                {
                    out_difference_curvature.write(
                        idx,
                        static_cast<float>(difference_curvature_value)
                    );
                }
            }

            if (std::isfinite(tangential_curvature_value) &&
                std::isfinite(principal_numeric_min))
            {
                // horizontal_excess_curvature = k_h - k_min
                horizontal_excess_curvature_value =
                    tangential_curvature_value - principal_numeric_min;
            }

            if (out_horizontal_excess_curvature.requested())
            {
                if (std::isfinite(horizontal_excess_curvature_value))
                {
                    out_horizontal_excess_curvature.write(
                        idx,
                        static_cast<float>(horizontal_excess_curvature_value)
                    );
                }
            }

            if (std::isfinite(vertical_curvature_value) &&
                std::isfinite(principal_numeric_min))
            {
                // vertical_excess_curvature = k_v - k_min
                vertical_excess_curvature_value =
                    vertical_curvature_value - principal_numeric_min;
            }

            if (out_vertical_excess_curvature.requested())
            {
                if (std::isfinite(vertical_excess_curvature_value))
                {
                    out_vertical_excess_curvature.write(
                        idx,
                        static_cast<float>(vertical_excess_curvature_value)
                    );
                }
            }

            if (std::isfinite(tangential_curvature_value) &&
                std::isfinite(vertical_curvature_value))
            {
                // accumulation_curvature = k_h * k_v
                // Units are 1/m^2, so this is not a true curvature.
                accumulation_curvature_value =
                    tangential_curvature_value * vertical_curvature_value;
            }

            if (out_accumulation_curvature.requested())
            {
                if (std::isfinite(accumulation_curvature_value))
                {
                    out_accumulation_curvature.write(
                        idx,
                        static_cast<float>(accumulation_curvature_value)
                    );
                }
            }

            if (std::isfinite(horizontal_excess_curvature_value) &&
                std::isfinite(vertical_excess_curvature_value))
            {
                // ring_curvature = k_he * k_ve
                // Units are 1/m^2, so this is not a true curvature.
                ring_curvature_value =
                    horizontal_excess_curvature_value *
                    vertical_excess_curvature_value;
            }

            if (out_ring_curvature.requested())
            {
                if (std::isfinite(ring_curvature_value))
                {
                    out_ring_curvature.write(
                        idx,
                        static_cast<float>(ring_curvature_value)
                    );
                }
            }
        }
    }

    polyfit_metric_result result;

    for (auto& kv : outputs)
    {
        result.rasters.emplace(
            kv.first,
            lsdtt_raster_cpu(
                std::move(kv.second.data),
                std::move(kv.second.mask),
                meta
            )
        );
    }

    return result;
}

} // namespace lsdtt
