// lsdtt_raster_linear_features_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lsdtt_raster_linear_features_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"   // moving_window_mean, window_shape

namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// resolve_thread_count
//
// 0 (or negative) -> all available hardware threads; otherwise the request.
// -----------------------------------------------------------------------------
int resolve_thread_count(int requested)
{
#ifdef _OPENMP
    if (requested <= 0) return omp_get_max_threads();
    return requested;
#else
    (void)requested;
    return 1;
#endif
}

// -----------------------------------------------------------------------------
// family_longest_path
//
// For a binary set b (1 = foreground), compute the longest near-straight path
// length passing through each foreground pixel within one orientation family,
// using two topological sweeps over the implied DAG:
//
//   up(p)   = longest path ENDING at p     (number of pixels)
//   down(p) = longest path STARTING at p
//   length through p = up(p) + down(p) - 1
//
// advance_columns == true  -> "horizontalish": each step advances one column
//                             and may deviate by -1/0/+1 row (+/-45 deg cone).
// advance_columns == false -> "verticalish":   each step advances one row
//                             and may deviate by -1/0/+1 column.
//
// up and down are overwritten (must be sized nrows*ncols).
// -----------------------------------------------------------------------------
void family_longest_path(
    const std::vector<std::uint8_t>& b,
    int                              nrows,
    int                              ncols,
    bool                             advance_columns,
    std::vector<int>&                up,
    std::vector<int>&                down)
{
    auto at = [ncols](int r, int c) -> std::size_t
    {
        return static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) + c;
    };

    if (advance_columns)
    {
        // up: predecessors are in column c-1 at rows r-1, r, r+1
        for (int c = 0; c < ncols; ++c)
        {
            for (int r = 0; r < nrows; ++r)
            {
                const std::size_t idx = at(r, c);
                if (!b[idx]) { up[idx] = 0; continue; }
                int best = 0;
                if (c > 0)
                {
                    for (int dr = -1; dr <= 1; ++dr)
                    {
                        const int rr = r + dr;
                        if (rr >= 0 && rr < nrows)
                            best = std::max(best, up[at(rr, c - 1)]);
                    }
                }
                up[idx] = 1 + best;
            }
        }
        // down: successors are in column c+1 at rows r-1, r, r+1
        for (int c = ncols - 1; c >= 0; --c)
        {
            for (int r = 0; r < nrows; ++r)
            {
                const std::size_t idx = at(r, c);
                if (!b[idx]) { down[idx] = 0; continue; }
                int best = 0;
                if (c + 1 < ncols)
                {
                    for (int dr = -1; dr <= 1; ++dr)
                    {
                        const int rr = r + dr;
                        if (rr >= 0 && rr < nrows)
                            best = std::max(best, down[at(rr, c + 1)]);
                    }
                }
                down[idx] = 1 + best;
            }
        }
    }
    else
    {
        // up: predecessors are in row r-1 at columns c-1, c, c+1
        for (int r = 0; r < nrows; ++r)
        {
            for (int c = 0; c < ncols; ++c)
            {
                const std::size_t idx = at(r, c);
                if (!b[idx]) { up[idx] = 0; continue; }
                int best = 0;
                if (r > 0)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        const int cc = c + dc;
                        if (cc >= 0 && cc < ncols)
                            best = std::max(best, up[at(r - 1, cc)]);
                    }
                }
                up[idx] = 1 + best;
            }
        }
        // down: successors are in row r+1 at columns c-1, c, c+1
        for (int r = nrows - 1; r >= 0; --r)
        {
            for (int c = 0; c < ncols; ++c)
            {
                const std::size_t idx = at(r, c);
                if (!b[idx]) { down[idx] = 0; continue; }
                int best = 0;
                if (r + 1 < nrows)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        const int cc = c + dc;
                        if (cc >= 0 && cc < ncols)
                            best = std::max(best, down[at(r + 1, cc)]);
                    }
                }
                down[idx] = 1 + best;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// path_opening_impl
//
// Shared implementation of the complete grayscale path opening. Returns the
// opened raster. The caller handles the opening/closing duality by negating.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu path_opening_impl(
    const lsdtt_raster_cpu& response,
    int                     min_path_length_pixels,
    int                     n_quantization_levels,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = response.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "grayscale_path_opening: input raster metadata is invalid.");
    }
    if (min_path_length_pixels < 2)
    {
        throw std::invalid_argument(
            "grayscale_path_opening: min_path_length_pixels must be >= 2.");
    }
    if (n_quantization_levels < 2)
    {
        throw std::invalid_argument(
            "grayscale_path_opening: n_quantization_levels must be >= 2.");
    }

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const int    L     = min_path_length_pixels;
    const int    K     = n_quantization_levels;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&        in_data = response.get_data();
    const std::vector<std::uint8_t>& in_mask = response.get_mask();

    auto is_valid = [&](std::size_t idx) -> bool
    {
        return in_mask[idx] != 0 &&
               !meta.is_nodata(in_data[idx]) &&
               std::isfinite(in_data[idx]);
    };

    // --- value range over valid pixels --------------------------------------
    float vmin = std::numeric_limits<float>::infinity();
    float vmax = -std::numeric_limits<float>::infinity();
    for (std::size_t idx = 0; idx < npix; ++idx)
    {
        if (!is_valid(idx)) continue;
        vmin = std::min(vmin, in_data[idx]);
        vmax = std::max(vmax, in_data[idx]);
    }

    std::vector<float>        out_data(npix,
        static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> out_mask(npix, std::uint8_t{0});

    // Degenerate cases: no valid pixels, or a constant field. A constant field
    // is its own opening, so copy the valid pixels straight through.
    const bool degenerate =
        !std::isfinite(vmin) || !std::isfinite(vmax) || (vmax <= vmin);
    if (degenerate)
    {
        for (std::size_t idx = 0; idx < npix; ++idx)
        {
            if (is_valid(idx))
            {
                out_data[idx] = in_data[idx];
                out_mask[idx] = std::uint8_t{1};
            }
        }
        return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
    }

    const double span    = static_cast<double>(vmax) - static_cast<double>(vmin);
    const double inv_span = 1.0 / span;

    // Quantise each valid pixel to a level index in [0, K-1] (nearest level).
    std::vector<int> q(npix, -1);
    for (std::size_t idx = 0; idx < npix; ++idx)
    {
        if (!is_valid(idx)) continue;
        const double frac =
            (static_cast<double>(in_data[idx]) - static_cast<double>(vmin))
            * inv_span;
        int lvl = static_cast<int>(std::lround(frac * (K - 1)));
        lvl = std::max(0, std::min(K - 1, lvl));
        q[idx] = lvl;
    }

    // out_level[p] = highest threshold k at which p still lies on a path of
    // length >= L (union of the two orientation families). Baseline 0 for valid
    // pixels = the global minimum level; -1 for invalid pixels.
    std::vector<int> out_level(npix, -1);
    for (std::size_t idx = 0; idx < npix; ++idx)
        if (is_valid(idx)) out_level[idx] = 0;

    const int nthreads = resolve_thread_count(n_threads);

    // Per-thread accumulator of the highest passing level, reduced by max.
    std::vector<std::vector<int>> local(
        static_cast<std::size_t>(nthreads),
        std::vector<int>(npix, 0));

#ifdef _OPENMP
    #pragma omp parallel num_threads(nthreads)
#endif
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        std::vector<int>& loc = local[static_cast<std::size_t>(tid)];

        std::vector<std::uint8_t> b(npix);
        std::vector<int>          up(npix);
        std::vector<int>          down(npix);

        // Thresholds k = 1 .. K-1 (k = 0 maps to vmin = background baseline).
#ifdef _OPENMP
        #pragma omp for schedule(dynamic)
#endif
        for (int k = 1; k < K; ++k)
        {
            for (std::size_t idx = 0; idx < npix; ++idx)
                b[idx] = (q[idx] >= k) ? std::uint8_t{1} : std::uint8_t{0};

            // Family 1: horizontalish
            family_longest_path(b, nrows, ncols, /*advance_columns=*/true,
                                up, down);
            for (std::size_t idx = 0; idx < npix; ++idx)
            {
                if (b[idx] && (up[idx] + down[idx] - 1) >= L)
                    loc[idx] = std::max(loc[idx], k);
            }

            // Family 2: verticalish
            family_longest_path(b, nrows, ncols, /*advance_columns=*/false,
                                up, down);
            for (std::size_t idx = 0; idx < npix; ++idx)
            {
                if (b[idx] && (up[idx] + down[idx] - 1) >= L)
                    loc[idx] = std::max(loc[idx], k);
            }
        }
    }

    // Reduce per-thread maxima.
    for (int t = 0; t < nthreads; ++t)
    {
        const std::vector<int>& loc = local[static_cast<std::size_t>(t)];
        for (std::size_t idx = 0; idx < npix; ++idx)
            if (out_level[idx] >= 0)
                out_level[idx] = std::max(out_level[idx], loc[idx]);
    }

    // Map surviving level back to a value, clamped to be anti-extensive.
    for (std::size_t idx = 0; idx < npix; ++idx)
    {
        if (out_level[idx] < 0) continue;   // invalid pixel
        const double value =
            static_cast<double>(vmin)
            + span * (static_cast<double>(out_level[idx])
                      / static_cast<double>(K - 1));
        out_data[idx] =
            std::min(static_cast<float>(value), in_data[idx]);
        out_mask[idx] = std::uint8_t{1};
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// grayscale_path_opening
// -----------------------------------------------------------------------------
lsdtt_raster_cpu grayscale_path_opening(
    const lsdtt_raster_cpu& response,
    int                     min_path_length_pixels,
    int                     n_quantization_levels,
    int                     n_threads)
{
    return path_opening_impl(
        response, min_path_length_pixels, n_quantization_levels, n_threads);
}

// -----------------------------------------------------------------------------
// grayscale_path_closing  =  -opening(-f)
// -----------------------------------------------------------------------------
lsdtt_raster_cpu grayscale_path_closing(
    const lsdtt_raster_cpu& response,
    int                     min_path_length_pixels,
    int                     n_quantization_levels,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = response.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "grayscale_path_closing: input raster metadata is invalid.");
    }

    const std::vector<float>&        in_data = response.get_data();
    const std::vector<std::uint8_t>& in_mask = response.get_mask();
    const std::size_t npix = in_data.size();

    // Negate valid pixels, leave nodata untouched.
    std::vector<float>        neg_data(in_data);
    std::vector<std::uint8_t> neg_mask(in_mask);
    for (std::size_t idx = 0; idx < npix; ++idx)
    {
        if (in_mask[idx] != 0 &&
            !meta.is_nodata(in_data[idx]) &&
            std::isfinite(in_data[idx]))
        {
            neg_data[idx] = -in_data[idx];
        }
    }
    lsdtt_raster_cpu neg(std::move(neg_data), std::move(neg_mask), meta);

    lsdtt_raster_cpu opened = path_opening_impl(
        neg, min_path_length_pixels, n_quantization_levels, n_threads);

    // Negate back.
    std::vector<float>&       od = opened.get_data();
    const std::vector<std::uint8_t>& om = opened.get_mask();
    for (std::size_t idx = 0; idx < od.size(); ++idx)
    {
        if (om[idx] != 0 && !meta.is_nodata(od[idx]) && std::isfinite(od[idx]))
            od[idx] = -od[idx];
    }
    return opened;
}

// -----------------------------------------------------------------------------
// compute_structure_tensor_coherence
// -----------------------------------------------------------------------------
structure_tensor_result compute_structure_tensor_coherence(
    const lsdtt_raster_cpu& field,
    double                  derivative_radius_metres,
    double                  integration_radius_metres,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = field.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "compute_structure_tensor_coherence: input metadata is invalid.");
    }
    if (integration_radius_metres <= 0.0)
    {
        throw std::invalid_argument(
            "compute_structure_tensor_coherence: integration_radius_metres "
            "must be > 0.");
    }
    if (derivative_radius_metres < 0.0)
    {
        throw std::invalid_argument(
            "compute_structure_tensor_coherence: derivative_radius_metres "
            "must be >= 0.");
    }

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const double cell  = meta.cell_size;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const float  nd    = static_cast<float>(meta.no_data_value);
    const int    nthreads = resolve_thread_count(n_threads);

    auto at = [ncols](int r, int c) -> std::size_t
    {
        return static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) + c;
    };

    // 1. Optional pre-smoothing to set the gradient noise scale.
    const lsdtt_raster_cpu  smoothed = (derivative_radius_metres > 0.0)
        ? moving_window_mean(field, derivative_radius_metres,
                             window_shape::circular)
        : field;
    const std::vector<float>&        f  = smoothed.get_data();
    const std::vector<std::uint8_t>& fm = smoothed.get_mask();

    auto fvalid = [&](std::size_t idx) -> bool
    {
        return fm[idx] != 0 && !meta.is_nodata(f[idx]) && std::isfinite(f[idx]);
    };

    // 2. Central-difference gradient; build outer-product components.
    std::vector<float>        jxx(npix, nd), jxy(npix, nd), jyy(npix, nd);
    std::vector<std::uint8_t> jmask(npix, std::uint8_t{0});
    const double inv2cell = 1.0 / (2.0 * cell);

#ifdef _OPENMP
    #pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t idx = at(r, c);
            if (r == 0 || r == nrows - 1 || c == 0 || c == ncols - 1)
                continue;
            if (!fvalid(idx) ||
                !fvalid(at(r, c - 1)) || !fvalid(at(r, c + 1)) ||
                !fvalid(at(r - 1, c)) || !fvalid(at(r + 1, c)))
                continue;

            const double gx =
                (static_cast<double>(f[at(r, c + 1)]) -
                 static_cast<double>(f[at(r, c - 1)])) * inv2cell;
            // row index increases southward
            const double gy =
                (static_cast<double>(f[at(r + 1, c)]) -
                 static_cast<double>(f[at(r - 1, c)])) * inv2cell;

            jxx[idx]   = static_cast<float>(gx * gx);
            jxy[idx]   = static_cast<float>(gx * gy);
            jyy[idx]   = static_cast<float>(gy * gy);
            jmask[idx] = std::uint8_t{1};
        }
    }

    // 3. Average each tensor component over the integration window.
    lsdtt_raster_cpu jxx_r(std::move(jxx), std::vector<std::uint8_t>(jmask),
                           meta);
    lsdtt_raster_cpu jxy_r(std::move(jxy), std::vector<std::uint8_t>(jmask),
                           meta);
    lsdtt_raster_cpu jyy_r(std::move(jyy), std::move(jmask), meta);

    const lsdtt_raster_cpu sxx =
        moving_window_mean(jxx_r, integration_radius_metres,
                           window_shape::circular);
    const lsdtt_raster_cpu sxy =
        moving_window_mean(jxy_r, integration_radius_metres,
                           window_shape::circular);
    const lsdtt_raster_cpu syy =
        moving_window_mean(jyy_r, integration_radius_metres,
                           window_shape::circular);

    const std::vector<float>&        axx = sxx.get_data();
    const std::vector<float>&        axy = sxy.get_data();
    const std::vector<float>&        ayy = syy.get_data();
    const std::vector<std::uint8_t>& amk = sxx.get_mask();

    // 4. Analytic 2x2 eigen-decomposition -> coherence / orientation / energy.
    std::vector<float>        coh (npix, nd), ori (npix, nd), eng (npix, nd);
    std::vector<std::uint8_t> omask(npix, std::uint8_t{0});

    constexpr double kRad2Deg = 57.29577951308232;

#ifdef _OPENMP
    #pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
    for (std::size_t idx = 0; idx < npix; ++idx)
    {
        if (amk[idx] == 0 || meta.is_nodata(axx[idx])) continue;

        const double a = static_cast<double>(axx[idx]);
        const double b = static_cast<double>(axy[idx]);
        const double cc= static_cast<double>(ayy[idx]);

        const double trace = a + cc;
        const double disc  = std::sqrt(std::max(0.0,
            0.25 * (a - cc) * (a - cc) + b * b));
        const double lam1  = 0.5 * trace + disc;   // >= lam2
        const double lam2  = 0.5 * trace - disc;

        eng[idx]   = static_cast<float>(lam1 + lam2);
        coh[idx]   = (trace > 0.0)
            ? static_cast<float>((lam1 - lam2) / trace)   // = 2*disc/trace
            : 0.0f;

        // Dominant gradient orientation (from +x, toward +y=south).
        const double theta_grad = 0.5 * std::atan2(2.0 * b, a - cc);
        // Feature tangent is perpendicular to the gradient.
        const double beta = theta_grad + 1.5707963267948966;  // +90 deg
        // Tangent vector in image coords (x east, y south).
        const double tx =  std::cos(beta);
        const double ty =  std::sin(beta);
        // Compass azimuth (clockwise from north): east = +x, north = -y.
        double az = std::atan2(tx, -ty) * kRad2Deg;   // degrees
        // Fold undirected line orientation into [0, 180).
        az = std::fmod(az, 180.0);
        if (az < 0.0) az += 180.0;
        ori[idx]   = static_cast<float>(az);

        omask[idx] = std::uint8_t{1};
    }

    structure_tensor_result out;
    out.coherence      =
        lsdtt_raster_cpu(std::move(coh), std::vector<std::uint8_t>(omask), meta);
    out.orientation_deg =
        lsdtt_raster_cpu(std::move(ori), std::vector<std::uint8_t>(omask), meta);
    out.energy         =
        lsdtt_raster_cpu(std::move(eng), std::move(omask), meta);
    return out;
}

// -----------------------------------------------------------------------------
// oriented_nms_hysteresis
// -----------------------------------------------------------------------------
oriented_ridge_result oriented_nms_hysteresis(
    const lsdtt_raster_cpu& response,
    const lsdtt_raster_cpu& orientation_deg,
    double                  high_threshold,
    double                  low_threshold,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = response.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "oriented_nms_hysteresis: input metadata is invalid.");
    }
    const lsdtt_raster_meta& ometa = orientation_deg.get_meta();
    if (ometa.n_rows != meta.n_rows || ometa.n_cols != meta.n_cols)
    {
        throw std::runtime_error(
            "oriented_nms_hysteresis: response and orientation rasters differ "
            "in dimensions.");
    }
    if (low_threshold > high_threshold)
    {
        throw std::invalid_argument(
            "oriented_nms_hysteresis: low_threshold must be <= high_threshold.");
    }

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const float  nd       = static_cast<float>(meta.no_data_value);
    const int    nthreads = resolve_thread_count(n_threads);
    const double deg2rad  = 0.017453292519943295;

    const std::vector<float>&        resp = response.get_data();
    const std::vector<std::uint8_t>& rmsk = response.get_mask();
    const std::vector<float>&        ori  = orientation_deg.get_data();
    const std::vector<std::uint8_t>& omsk = orientation_deg.get_mask();

    auto at = [ncols](int r, int c) -> std::size_t
    {
        return static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) + c;
    };
    auto rvalid = [&](std::size_t idx) -> bool
    {
        return rmsk[idx] != 0 && !meta.is_nodata(resp[idx]) &&
               std::isfinite(resp[idx]);
    };

    // Bilinear sample of the response at fractional (fr, fc); returns -inf if
    // any of the four corners is off-grid or invalid (so the centre wins NMS).
    auto bilerp = [&](double fr, double fc) -> double
    {
        const double NEG = -std::numeric_limits<double>::infinity();
        const int r0 = static_cast<int>(std::floor(fr));
        const int c0 = static_cast<int>(std::floor(fc));
        if (r0 < 0 || r0 + 1 >= nrows || c0 < 0 || c0 + 1 >= ncols) return NEG;
        const std::size_t i00 = at(r0, c0),     i01 = at(r0, c0 + 1);
        const std::size_t i10 = at(r0 + 1, c0), i11 = at(r0 + 1, c0 + 1);
        if (!rvalid(i00) || !rvalid(i01) || !rvalid(i10) || !rvalid(i11))
            return NEG;
        const double tr = fr - r0, tc = fc - c0;
        const double top = (1.0 - tc) * resp[i00] + tc * resp[i01];
        const double bot = (1.0 - tc) * resp[i10] + tc * resp[i11];
        return (1.0 - tr) * top + tr * bot;
    };

    // --- Oriented non-maximum suppression ---
    std::vector<float>        thin (npix, nd);
    std::vector<std::uint8_t> tmask(npix, std::uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t idx = at(r, c);
            if (!rvalid(idx) || omsk[idx] == 0 || meta.is_nodata(ori[idx]))
                continue;

            const double az = static_cast<double>(ori[idx]) * deg2rad;
            // Feature normal (across the ridge) in grid steps: (drow, dcol).
            const double dcol = std::cos(az);
            const double drow = std::sin(az);

            const double here = static_cast<double>(resp[idx]);
            const double plus  = bilerp(r + drow, c + dcol);
            const double minus = bilerp(r - drow, c - dcol);

            tmask[idx] = std::uint8_t{1};   // valid output pixel either way
            thin[idx]  = (here >= plus && here >= minus)
                ? resp[idx] : 0.0f;
        }
    }

    // --- Hysteresis linking (8-connected flood-fill from strong seeds) ---
    std::vector<float>        mdata(npix, nd);
    std::vector<std::uint8_t> mmask(npix, std::uint8_t{0});
    for (std::size_t idx = 0; idx < npix; ++idx)
        if (tmask[idx]) { mdata[idx] = 0.0f; mmask[idx] = std::uint8_t{1}; }

    std::vector<std::uint8_t> kept(npix, std::uint8_t{0});
    std::vector<std::size_t>  stack;
    for (std::size_t idx = 0; idx < npix; ++idx)
    {
        if (tmask[idx] && static_cast<double>(thin[idx]) >= high_threshold)
        {
            kept[idx] = std::uint8_t{1};
            stack.push_back(idx);
        }
    }

    while (!stack.empty())
    {
        const std::size_t idx = stack.back();
        stack.pop_back();
        const int r = static_cast<int>(idx / ncols);
        const int c = static_cast<int>(idx % ncols);
        for (int dr = -1; dr <= 1; ++dr)
        {
            for (int dc = -1; dc <= 1; ++dc)
            {
                if (dr == 0 && dc == 0) continue;
                const int rr = r + dr, ccx = c + dc;
                if (rr < 0 || rr >= nrows || ccx < 0 || ccx >= ncols) continue;
                const std::size_t nidx = at(rr, ccx);
                if (kept[nidx] || !tmask[nidx]) continue;
                if (static_cast<double>(thin[nidx]) >= low_threshold)
                {
                    kept[nidx] = std::uint8_t{1};
                    stack.push_back(nidx);
                }
            }
        }
    }

    for (std::size_t idx = 0; idx < npix; ++idx)
        if (kept[idx]) mdata[idx] = 1.0f;

    oriented_ridge_result out;
    out.thinned = lsdtt_raster_cpu(std::move(thin),  std::move(tmask), meta);
    out.mask    = lsdtt_raster_cpu(std::move(mdata), std::move(mmask), meta);
    return out;
}

// -----------------------------------------------------------------------------
// coherence_enhancing_diffusion
// -----------------------------------------------------------------------------
lsdtt_raster_cpu coherence_enhancing_diffusion(
    const lsdtt_raster_cpu& field,
    int                     n_iterations,
    double                  time_step,
    double                  derivative_radius_metres,
    double                  integration_radius_metres,
    double                  contrast,
    double                  min_diffusivity,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = field.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "coherence_enhancing_diffusion: input metadata is invalid.");
    }
    if (n_iterations < 1)
        throw std::invalid_argument(
            "coherence_enhancing_diffusion: n_iterations must be >= 1.");
    if (time_step <= 0.0 || time_step > 0.25)
        throw std::invalid_argument(
            "coherence_enhancing_diffusion: time_step must be in (0, 0.25].");
    if (integration_radius_metres <= 0.0)
        throw std::invalid_argument(
            "coherence_enhancing_diffusion: integration_radius_metres must be > 0.");
    if (contrast <= 0.0)
        throw std::invalid_argument(
            "coherence_enhancing_diffusion: contrast must be > 0.");
    if (min_diffusivity <= 0.0 || min_diffusivity >= 1.0)
        throw std::invalid_argument(
            "coherence_enhancing_diffusion: min_diffusivity must be in (0, 1).");

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const int    nthreads = resolve_thread_count(n_threads);
    const double alpha    = min_diffusivity;
    const double deg2rad  = 0.017453292519943295;

    auto at = [ncols](int r, int c) -> std::size_t
    {
        return static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) + c;
    };

    // Working copy and a fixed validity mask (nodata stays nodata throughout).
    lsdtt_raster_cpu work = field;
    std::vector<float>&              u  = work.get_data();
    const std::vector<std::uint8_t>& wm = work.get_mask();
    std::vector<std::uint8_t> valid(npix);
    for (std::size_t idx = 0; idx < npix; ++idx)
        valid[idx] = (wm[idx] != 0 && !meta.is_nodata(u[idx]) &&
                      std::isfinite(u[idx])) ? std::uint8_t{1} : std::uint8_t{0};

    // No-flux sampler: substitute the centre value for any off-grid or invalid
    // neighbour (gives a zero gradient across that boundary -> no flux).
    auto sample = [&](const std::vector<float>& a, int r, int c, float ctr) -> double
    {
        if (r < 0 || r >= nrows || c < 0 || c >= ncols) return ctr;
        const std::size_t idx = at(r, c);
        if (!valid[idx]) return ctr;
        return static_cast<double>(a[idx]);
    };

    std::vector<float> dxx(npix), dxy(npix), dyy(npix);
    std::vector<float> jx(npix, 0.0f), jy(npix, 0.0f);
    std::vector<float> unew(npix);

    for (int iter = 0; iter < n_iterations; ++iter)
    {
        // 1. Structure tensor of the current field -> coherence + orientation.
        const structure_tensor_result st = compute_structure_tensor_coherence(
            work, derivative_radius_metres, integration_radius_metres, nthreads);
        const std::vector<float>&        kap = st.coherence.get_data();
        const std::vector<float>&        ori = st.orientation_deg.get_data();
        const std::vector<std::uint8_t>& stm = st.coherence.get_mask();

        // 2. Diffusion tensor D per pixel.
#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
        for (std::size_t idx = 0; idx < npix; ++idx)
        {
            double a = alpha, b = 0.0, c = alpha;   // isotropic floor default
            if (stm[idx] != 0 && !meta.is_nodata(kap[idx]))
            {
                const double az  = static_cast<double>(ori[idx]) * deg2rad;
                // tangent (along feature) and normal (across feature) in image
                // coords (x east, y south); azimuth is clockwise from north.
                const double v2x = std::sin(az),  v2y = -std::cos(az);  // along
                const double v1x = std::cos(az),  v1y =  std::sin(az);  // across
                const double k2  = static_cast<double>(kap[idx]) *
                                   static_cast<double>(kap[idx]);
                const double lam_perp = alpha;
                const double lam_par  = (k2 > 0.0)
                    ? alpha + (1.0 - alpha) * std::exp(-contrast / k2)
                    : alpha;
                a = lam_perp * v1x * v1x + lam_par * v2x * v2x;
                b = lam_perp * v1x * v1y + lam_par * v2x * v2y;
                c = lam_perp * v1y * v1y + lam_par * v2y * v2y;
            }
            dxx[idx] = static_cast<float>(a);
            dxy[idx] = static_cast<float>(b);
            dyy[idx] = static_cast<float>(c);
        }

        // 3a. Fluxes j = D grad u  (pixel-unit central-difference gradient).
#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
        for (int r = 0; r < nrows; ++r)
        {
            for (int c = 0; c < ncols; ++c)
            {
                const std::size_t idx = at(r, c);
                if (!valid[idx]) { jx[idx] = 0.0f; jy[idx] = 0.0f; continue; }
                const float ctr = u[idx];
                const double ux = 0.5 * (sample(u, r, c + 1, ctr) -
                                         sample(u, r, c - 1, ctr));
                const double uy = 0.5 * (sample(u, r + 1, c, ctr) -
                                         sample(u, r - 1, c, ctr));
                jx[idx] = static_cast<float>(dxx[idx] * ux + dxy[idx] * uy);
                jy[idx] = static_cast<float>(dxy[idx] * ux + dyy[idx] * uy);
            }
        }

        // 3b. Explicit update  u += tau * div(j).
#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
        for (int r = 0; r < nrows; ++r)
        {
            for (int c = 0; c < ncols; ++c)
            {
                const std::size_t idx = at(r, c);
                if (!valid[idx]) { unew[idx] = u[idx]; continue; }
                const double djx = 0.5 * (sample(jx, r, c + 1, jx[idx]) -
                                          sample(jx, r, c - 1, jx[idx]));
                const double djy = 0.5 * (sample(jy, r + 1, c, jy[idx]) -
                                          sample(jy, r - 1, c, jy[idx]));
                unew[idx] = static_cast<float>(
                    static_cast<double>(u[idx]) + time_step * (djx + djy));
            }
        }

        u.swap(unew);
    }

    return work;
}

// -----------------------------------------------------------------------------
// trim_valid_border
// -----------------------------------------------------------------------------
lsdtt_raster_cpu trim_valid_border(
    const lsdtt_raster_cpu& raster,
    int                     border_pixels)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "trim_valid_border: input raster metadata is invalid.");
    }
    if (border_pixels < 1)
    {
        throw std::invalid_argument(
            "trim_valid_border: border_pixels must be >= 1.");
    }

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const int    N     = border_pixels;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&        in_data = raster.get_data();
    const std::vector<std::uint8_t>& in_mask = raster.get_mask();

    auto at = [ncols](int r, int c) -> std::size_t
    {
        return static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) + c;
    };

    std::vector<std::uint8_t> valid(npix);
    for (std::size_t idx = 0; idx < npix; ++idx)
    {
        valid[idx] = (in_mask[idx] != 0 &&
                      !meta.is_nodata(in_data[idx]) &&
                      std::isfinite(in_data[idx]))
                     ? std::uint8_t{1} : std::uint8_t{0};
    }

    // Horizontal erosion: survive only if all columns in [c-N, c+N] (within the
    // grid) are valid. Out-of-grid columns count as invalid.
    std::vector<std::uint8_t> htmp(npix, std::uint8_t{0});
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (c - N < 0 || c + N >= ncols) continue;   // grid-edge margin
            bool ok = true;
            for (int dc = -N; dc <= N && ok; ++dc)
                if (!valid[at(r, c + dc)]) ok = false;
            if (ok) htmp[at(r, c)] = std::uint8_t{1};
        }
    }

    // Vertical erosion over the horizontal result.
    std::vector<float>        out_data(in_data);
    std::vector<std::uint8_t> out_mask(in_mask);
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            bool keep = true;
            if (r - N < 0 || r + N >= nrows) keep = false;   // grid-edge margin
            for (int dr = -N; dr <= N && keep; ++dr)
                if (!htmp[at(r + dr, c)]) keep = false;
            if (!keep)
            {
                out_data[at(r, c)] = static_cast<float>(meta.no_data_value);
                out_mask[at(r, c)] = std::uint8_t{0};
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

} // namespace lsdtt
