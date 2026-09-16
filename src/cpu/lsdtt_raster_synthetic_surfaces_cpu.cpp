// lsdtt_raster_synthetic_surfaces.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"

#include "../utils/lsdtt_thread_utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>
#include <cstdint>

#ifdef LSDTT_HAS_FFTW
#include <fftw3.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------
static constexpr double k_pi = 3.14159265358979323846264338327950288;

// -----------------------------------------------------------------------------
// flat_index
// -----------------------------------------------------------------------------
static inline std::size_t flat_index(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
         + static_cast<std::size_t>(col);
}

// -----------------------------------------------------------------------------
// make_raster
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu make_raster(
    std::vector<float>&&     data,
    const lsdtt_raster_meta& meta)
{
    const std::size_t n_pixels =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    std::vector<uint8_t> mask(n_pixels, uint8_t{1});

    return lsdtt_raster_cpu(
        std::move(data),
        std::move(mask),
        meta
    );
}

// -----------------------------------------------------------------------------
// build_meta
// -----------------------------------------------------------------------------
static lsdtt_raster_meta build_meta(
    int    n_rows,
    int    n_cols,
    double cell_size)
{
    if (n_rows <= 0 || n_cols <= 0)
    {
        throw std::runtime_error(
            "lsdtt synthetic surface generator: n_rows and n_cols must be > 0"
        );
    }

    if (!std::isfinite(cell_size) || cell_size <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt synthetic surface generator: cell_size must be finite and > 0"
        );
    }

    lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell_size;
    meta.no_data_value = -9999.0;
    meta.data_type     = lsdtt_pixel_type::float32;

    // Synthetic surfaces have no input raster, so default to a sensible
    // projected georeference centred on "Watopia" (Teanu Island, UTM zone 58S).
    meta.set_default_synthetic_georeference();

    return meta;
}



// -----------------------------------------------------------------------------
// Hash-based deterministic RNG utilities
//
// We use a SplitMix64-style mixer keyed by:
//   seed, tag, row, col, iteration...
//
// This makes the stochastic surfaces deterministic for a fixed seed and
// independent of thread scheduling.
// -----------------------------------------------------------------------------
static inline std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    x =  x ^ (x >> 31);
    return x;
}

static inline std::uint64_t hash_combine_u64(std::uint64_t h, std::uint64_t v)
{
    return splitmix64(h ^ (v + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2)));
}

static std::uint64_t make_hash_key(
    int           seed,
    std::uint64_t tag,
    std::uint64_t a = 0,
    std::uint64_t b = 0,
    std::uint64_t c = 0,
    std::uint64_t d = 0)
{
    std::uint64_t h = splitmix64(static_cast<std::uint64_t>(
        static_cast<std::uint32_t>(seed)
    ));
    h = hash_combine_u64(h, tag);
    h = hash_combine_u64(h, a);
    h = hash_combine_u64(h, b);
    h = hash_combine_u64(h, c);
    h = hash_combine_u64(h, d);
    return splitmix64(h);
}

static inline double unit_uniform_01(std::uint64_t key)
{
    // Use the top 53 bits to form a double in [0, 1)
    return static_cast<double>(key >> 11) *
           (1.0 / 9007199254740992.0);  // 2^53
}

static inline double symmetric_uniform(std::uint64_t key)
{
    return 2.0 * unit_uniform_01(key) - 1.0;
}

// -----------------------------------------------------------------------------
// normalise_to_mean_std
//
// Serial on purpose for strict determinism across thread counts.
// -----------------------------------------------------------------------------
static void normalise_to_mean_std(
    std::vector<float>& data,
    double              target_mean,
    double              target_stddev)
{
    if (!std::isfinite(target_mean) ||
        !std::isfinite(target_stddev) ||
        target_stddev < 0.0)
    {
        throw std::runtime_error(
            "lsdtt synthetic noise generator: invalid target mean/stddev"
        );
    }

    const std::size_t n = data.size();
    if (n == 0)
    {
        throw std::runtime_error(
            "lsdtt synthetic noise generator: empty data array"
        );
    }

    double mean = 0.0;
    for (float v : data)
    {
        mean += static_cast<double>(v);
    }
    mean /= static_cast<double>(n);

    double var = 0.0;
    for (float v : data)
    {
        const double dv = static_cast<double>(v) - mean;
        var += dv * dv;
    }
    var /= static_cast<double>(n);

    const double stddev = std::sqrt(std::max(0.0, var));

    if (stddev < 1e-20)
    {
        for (float& v : data)
        {
            v = static_cast<float>(target_mean);
        }
        return;
    }

    const double scale = target_stddev / stddev;

    for (float& v : data)
    {
        const double out = target_mean + (static_cast<double>(v) - mean) * scale;
        v = static_cast<float>(out);
    }
}

// -----------------------------------------------------------------------------
// wrap_index
// -----------------------------------------------------------------------------
static inline int wrap_index(int i, int n)
{
    i %= n;
    if (i < 0) i += n;
    return i;
}

// -----------------------------------------------------------------------------
// binomial_blur_5x5_periodic
//
// Separable 5-tap periodic blur:
//   kernel = [1, 4, 6, 4, 1] / 16
//
// Deterministic under multithreading because each output pixel is computed
// independently with a fixed evaluation order.
// -----------------------------------------------------------------------------
static void binomial_blur_5x5_periodic(
    const std::vector<float>& src,
    std::vector<float>&       tmp,
    std::vector<float>&       dst,
    int                       n_rows,
    int                       n_cols,
    int                       n_threads)
{
    const int threads = resolve_thread_count(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const float v =
                  1.0f * src[flat_index(row, wrap_index(col - 2, n_cols), n_cols)]
                + 4.0f * src[flat_index(row, wrap_index(col - 1, n_cols), n_cols)]
                + 6.0f * src[flat_index(row, wrap_index(col    , n_cols), n_cols)]
                + 4.0f * src[flat_index(row, wrap_index(col + 1, n_cols), n_cols)]
                + 1.0f * src[flat_index(row, wrap_index(col + 2, n_cols), n_cols)];

            tmp[flat_index(row, col, n_cols)] = v * (1.0f / 16.0f);
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const float v =
                  1.0f * tmp[flat_index(wrap_index(row - 2, n_rows), col, n_cols)]
                + 4.0f * tmp[flat_index(wrap_index(row - 1, n_rows), col, n_cols)]
                + 6.0f * tmp[flat_index(wrap_index(row    , n_rows), col, n_cols)]
                + 4.0f * tmp[flat_index(wrap_index(row + 1, n_rows), col, n_cols)]
                + 1.0f * tmp[flat_index(wrap_index(row + 2, n_rows), col, n_cols)];

            dst[flat_index(row, col, n_cols)] = v * (1.0f / 16.0f);
        }
    }
}

// -----------------------------------------------------------------------------
// generate_white_noise_core
//
// Generates symmetric white noise in [-1, 1] using deterministic per-pixel hash
// keys. Intended as a base field for red/blue noise.
// -----------------------------------------------------------------------------
static std::vector<float> generate_white_noise_core(
    const lsdtt_raster_meta& meta,
    int                      seed,
    std::uint64_t            tag,
    int                      n_threads)
{
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;

    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    std::vector<float> data(n_pixels, 0.0f);

    const int threads = resolve_thread_count(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::uint64_t key = make_hash_key(
                seed,
                tag,
                static_cast<std::uint64_t>(row),
                static_cast<std::uint64_t>(col)
            );

            data[flat_index(row, col, n_cols)] =
                static_cast<float>(symmetric_uniform(key));
        }
    }

    return data;
}

// -----------------------------------------------------------------------------
// superimpose_ns_parabola
// -----------------------------------------------------------------------------
lsdtt_raster_cpu superimpose_ns_parabola(
    const lsdtt_raster_cpu& dem,
    float peak_elevation)
{
    const lsdtt_raster_meta& meta  = dem.get_meta();
    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float centre = static_cast<float>(nrows - 1) * 0.5f;

    std::vector<float> data = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    for (int row = 0; row < nrows; ++row)
    {
        const float t = (static_cast<float>(row) - centre) / centre;
        const float dz = peak_elevation * (1.0f - t * t);
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = flat_index(row, col, ncols);
            if (mask[idx])
                data[idx] += dz;
        }
    }

    return lsdtt_raster_cpu(std::move(data), std::vector<uint8_t>(mask), meta);
}

// -----------------------------------------------------------------------------
// superimpose_open_boundary_gradient
// -----------------------------------------------------------------------------
lsdtt_raster_cpu superimpose_open_boundary_gradient(
    const lsdtt_raster_cpu& dem,
    const std::string& boundary_conditions,
    float slope)
{
    if (slope < 0.0f)
        throw std::invalid_argument("open-boundary gradient slope must be non-negative");
    if (boundary_conditions.size() != 4)
        throw std::invalid_argument("boundary_conditions must contain four NESW codes");

    auto is_open = [](char bc) { return bc == 'o' || bc == 'b'; };
    const bool open_n = is_open(boundary_conditions[0]);
    const bool open_e = is_open(boundary_conditions[1]);
    const bool open_s = is_open(boundary_conditions[2]);
    const bool open_w = is_open(boundary_conditions[3]);
    if (slope == 0.0f || !(open_n || open_e || open_s || open_w)) return dem;

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const float dx = static_cast<float>(meta.cell_size);
    std::vector<float> data = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    for (int row = 0; row < nrows; ++row)
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = flat_index(row, col, ncols);
            if (!mask[idx]) continue;

            int distance_cells = nrows + ncols;
            if (open_n) distance_cells = std::min(distance_cells, row);
            if (open_e) distance_cells = std::min(distance_cells, ncols - 1 - col);
            if (open_s) distance_cells = std::min(distance_cells, nrows - 1 - row);
            if (open_w) distance_cells = std::min(distance_cells, col);
            data[idx] += slope * dx * static_cast<float>(distance_cells);
        }

    return lsdtt_raster_cpu(std::move(data), std::vector<uint8_t>(mask), meta);
}

// -----------------------------------------------------------------------------
// generate_uniform_random_noise_surface
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_uniform_random_noise_surface(
    const lsdtt_raster_meta& meta,
    int    seed,
    double min_value,
    double max_value,
    int    n_threads)
{
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::generate_uniform_random_noise_surface: invalid raster metadata"
        );
    }

    if (!std::isfinite(min_value) ||
        !std::isfinite(max_value) ||
        min_value > max_value)
    {
        throw std::runtime_error(
            "lsdtt::generate_uniform_random_noise_surface: invalid value range"
        );
    }

    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;

    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    std::vector<float> data(n_pixels, 0.0f);

    const double span = max_value - min_value;
    const int threads = resolve_thread_count(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::uint64_t key = make_hash_key(
                seed,
                0x554E49464F524DULL, // "UNIFORM" tag-ish
                static_cast<std::uint64_t>(row),
                static_cast<std::uint64_t>(col)
            );

            const double u = unit_uniform_01(key);
            const double z = min_value + span * u;

            data[flat_index(row, col, n_cols)] = static_cast<float>(z);
        }
    }

    return make_raster(std::move(data), meta);
}

lsdtt_raster_cpu generate_uniform_random_noise_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed,
    double min_value,
    double max_value,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);

    return generate_uniform_random_noise_surface(
        meta,
        seed,
        min_value,
        max_value,
        n_threads
    );
}

// -----------------------------------------------------------------------------
// generate_red_noise_surface
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_red_noise_surface(
    const lsdtt_raster_meta& meta,
    int    seed,
    double target_mean,
    double target_stddev,
    int    blur_passes,
    int    n_threads)
{
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::generate_red_noise_surface: invalid raster metadata"
        );
    }

    if (blur_passes < 1)
    {
        throw std::runtime_error(
            "lsdtt::generate_red_noise_surface: blur_passes must be >= 1"
        );
    }

    std::vector<float> cur = generate_white_noise_core(
        meta, seed, 0x5245444E4F495345ULL, n_threads // "REDNOISE"
    );

    std::vector<float> tmp(cur.size(), 0.0f);
    std::vector<float> nxt(cur.size(), 0.0f);

    for (int pass = 0; pass < blur_passes; ++pass)
    {
        binomial_blur_5x5_periodic(
            cur, tmp, nxt,
            meta.n_rows, meta.n_cols,
            n_threads
        );
        cur.swap(nxt);
    }

    normalise_to_mean_std(cur, target_mean, target_stddev);

    return make_raster(std::move(cur), meta);
}

lsdtt_raster_cpu generate_red_noise_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed,
    double target_mean,
    double target_stddev,
    int    blur_passes,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);

    return generate_red_noise_surface(
        meta,
        seed,
        target_mean,
        target_stddev,
        blur_passes,
        n_threads
    );
}

// -----------------------------------------------------------------------------
// generate_blue_noise_surface
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_blue_noise_surface(
    const lsdtt_raster_meta& meta,
    int    seed,
    double target_mean,
    double target_stddev,
    int    blur_passes,
    int    n_threads)
{
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::generate_blue_noise_surface: invalid raster metadata"
        );
    }

    if (blur_passes < 1)
    {
        throw std::runtime_error(
            "lsdtt::generate_blue_noise_surface: blur_passes must be >= 1"
        );
    }

    const std::vector<float> white = generate_white_noise_core(
        meta, seed, 0x424C55454E4F4953ULL, n_threads // "BLUENOIS"
    );

    std::vector<float> low = white;
    std::vector<float> tmp(low.size(), 0.0f);
    std::vector<float> nxt(low.size(), 0.0f);

    for (int pass = 0; pass < blur_passes; ++pass)
    {
        binomial_blur_5x5_periodic(
            low, tmp, nxt,
            meta.n_rows, meta.n_cols,
            n_threads
        );
        low.swap(nxt);
    }

    std::vector<float> blue(white.size(), 0.0f);

    const int threads = resolve_thread_count(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (std::size_t i = 0; i < blue.size(); ++i)
    {
        blue[i] = white[i] - low[i];
    }

    normalise_to_mean_std(blue, target_mean, target_stddev);

    return make_raster(std::move(blue), meta);
}

lsdtt_raster_cpu generate_blue_noise_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed,
    double target_mean,
    double target_stddev,
    int    blur_passes,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);

    return generate_blue_noise_surface(
        meta,
        seed,
        target_mean,
        target_stddev,
        blur_passes,
        n_threads
    );
}

// -----------------------------------------------------------------------------
// generate_diamond_square_surface
//
// Generates on the smallest enclosing square grid of size (2^k + 1), then
// crops the upper-left n_rows x n_cols block.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_diamond_square_surface(
    const lsdtt_raster_meta& meta,
    int    seed,
    double roughness,
    double initial_scale,
    double target_mean,
    double target_stddev,
    int    n_threads)
{
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::generate_diamond_square_surface: invalid raster metadata"
        );
    }

    if (!std::isfinite(roughness) || roughness <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::generate_diamond_square_surface: roughness must be > 0"
        );
    }

    if (!std::isfinite(initial_scale) || initial_scale < 0.0)
    {
        throw std::runtime_error(
            "lsdtt::generate_diamond_square_surface: initial_scale must be >= 0"
        );
    }

    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;

    const int span_target = std::max(n_rows - 1, n_cols - 1);

    int span = 1;
    while (span < span_target)
    {
        if (span > std::numeric_limits<int>::max() / 2)
        {
            throw std::runtime_error(
                "lsdtt::generate_diamond_square_surface: grid size overflow"
            );
        }
        span *= 2;
    }

    const int side = span + 1;

    std::vector<double> grid(
        static_cast<std::size_t>(side) * static_cast<std::size_t>(side),
        0.0
    );

    auto gidx = [side](int r, int c) -> std::size_t
    {
        return static_cast<std::size_t>(r) * static_cast<std::size_t>(side)
             + static_cast<std::size_t>(c);
    };

    // Seed the four corners deterministically
    grid[gidx(0,       0      )] = initial_scale * symmetric_uniform(
        make_hash_key(seed, 0x44495351434F524EULL, 0, 0)
    );
    grid[gidx(0,       side-1 )] = initial_scale * symmetric_uniform(
        make_hash_key(seed, 0x44495351434F524EULL, 0, 1)
    );
    grid[gidx(side-1,  0      )] = initial_scale * symmetric_uniform(
        make_hash_key(seed, 0x44495351434F524EULL, 1, 0)
    );
    grid[gidx(side-1,  side-1 )] = initial_scale * symmetric_uniform(
        make_hash_key(seed, 0x44495351434F524EULL, 1, 1)
    );

    int step   = span;
    double amp = initial_scale;
    const int threads = resolve_thread_count(n_threads);

    while (step > 1)
    {
        const int half = step / 2;

        // Diamond step
#ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
        for (int row = half; row < side - 1; row += step)
        {
            for (int col = half; col < side - 1; col += step)
            {
                const double z00 = grid[gidx(row - half, col - half)];
                const double z01 = grid[gidx(row - half, col + half)];
                const double z10 = grid[gidx(row + half, col - half)];
                const double z11 = grid[gidx(row + half, col + half)];

                const double avg = 0.25 * (z00 + z01 + z10 + z11);
                const double dz  = amp * symmetric_uniform(
                    make_hash_key(
                        seed,
                        0x444953514449414DULL, // "DISQDIAM"
                        static_cast<std::uint64_t>(step),
                        static_cast<std::uint64_t>(row),
                        static_cast<std::uint64_t>(col)
                    )
                );

                grid[gidx(row, col)] = avg + dz;
            }
        }

        // Square step
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
#endif
        for (int row = 0; row < side; row += half)
        {
            const int col_start = ((row / half) % 2 == 0) ? half : 0;

            for (int col = col_start; col < side; col += step)
            {
                double sum = 0.0;
                int    n   = 0;

                if (row - half >= 0)
                {
                    sum += grid[gidx(row - half, col)];
                    ++n;
                }
                if (row + half < side)
                {
                    sum += grid[gidx(row + half, col)];
                    ++n;
                }
                if (col - half >= 0)
                {
                    sum += grid[gidx(row, col - half)];
                    ++n;
                }
                if (col + half < side)
                {
                    sum += grid[gidx(row, col + half)];
                    ++n;
                }

                const double avg = (n > 0) ? (sum / static_cast<double>(n)) : 0.0;
                const double dz  = amp * symmetric_uniform(
                    make_hash_key(
                        seed,
                        0x4449535153515541ULL, // "DISQSQUA"
                        static_cast<std::uint64_t>(step),
                        static_cast<std::uint64_t>(row),
                        static_cast<std::uint64_t>(col)
                    )
                );

                grid[gidx(row, col)] = avg + dz;
            }
        }

        step = half;
        amp *= roughness;
    }

    // Crop to requested output size and convert to float
    std::vector<float> data(
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols),
        0.0f
    );

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            data[flat_index(row, col, n_cols)] =
                static_cast<float>(grid[gidx(row, col)]);
        }
    }

    normalise_to_mean_std(data, target_mean, target_stddev);

    return make_raster(std::move(data), meta);
}

lsdtt_raster_cpu generate_diamond_square_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed,
    double roughness,
    double initial_scale,
    double target_mean,
    double target_stddev,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);

    return generate_diamond_square_surface(
        meta,
        seed,
        roughness,
        initial_scale,
        target_mean,
        target_stddev,
        n_threads
    );
}


// -----------------------------------------------------------------------------
// generate_power_law_fractal_surface_fft
//
// FFT-based synthetic fractal DEM generator for spectral background modelling.
//
// Method:
//   1. Generate deterministic white noise on a padded power-of-2 grid
//   2. Forward complex FFT
//   3. Apply amplitude scaling proportional to f^(-beta/2)
//   4. Inverse FFT
//   5. Crop back to the requested raster size
//   6. Normalise to target mean and target standard deviation
//
// Notes:
//   - beta is the PSD slope exponent, so amplitude is scaled by beta/2
//   - same seed => identical output
//   - if built without FFTW, throws at runtime with a clear message
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_power_law_fractal_surface_fft(
    const lsdtt_raster_meta& meta,
    int    seed,
    double beta,
    double target_mean,
    double target_stddev,
    int    n_threads)
{
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::generate_power_law_fractal_surface_fft: invalid raster metadata"
        );
    }

    meta.check_projected_metres();

    if (!std::isfinite(beta))
    {
        throw std::runtime_error(
            "lsdtt::generate_power_law_fractal_surface_fft: beta must be finite"
        );
    }

    if (!std::isfinite(target_mean) ||
        !std::isfinite(target_stddev) ||
        target_stddev < 0.0)
    {
        throw std::runtime_error(
            "lsdtt::generate_power_law_fractal_surface_fft: invalid target mean/stddev"
        );
    }

#ifndef LSDTT_HAS_FFTW
    (void)seed;
    (void)n_threads;

    throw std::runtime_error(
        "lsdtt::generate_power_law_fractal_surface_fft: lsdtt3 was built "
        "without FFTW support. Install FFTW development files and rebuild."
    );
#else
    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;

    auto next_pow2_local = [](int n) -> int
    {
        if (n <= 1) return 1;
        int p = 1;
        while (p < n) p <<= 1;
        return p;
    };

    const int padded_rows = next_pow2_local(n_rows);
    const int padded_cols = next_pow2_local(n_cols);

    const double dfx = 1.0 / (cell_size * static_cast<double>(padded_cols));
    const double dfy = 1.0 / (cell_size * static_cast<double>(padded_rows));

    const std::size_t n_fft =
        static_cast<std::size_t>(padded_rows) *
        static_cast<std::size_t>(padded_cols);

    fftw_complex* in  = nullptr;
    fftw_complex* out = nullptr;
    fftw_plan forward  = nullptr;
    fftw_plan backward = nullptr;

    try
    {
        in = reinterpret_cast<fftw_complex*>(
            fftw_malloc(sizeof(fftw_complex) * n_fft)
        );
        out = reinterpret_cast<fftw_complex*>(
            fftw_malloc(sizeof(fftw_complex) * n_fft)
        );

        if (!in || !out)
        {
            throw std::runtime_error(
                "lsdtt::generate_power_law_fractal_surface_fft: FFTW allocation failed"
            );
        }

        forward = fftw_plan_dft_2d(
            padded_rows, padded_cols,
            in, out,
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );

        backward = fftw_plan_dft_2d(
            padded_rows, padded_cols,
            out, in,
            FFTW_BACKWARD,
            FFTW_ESTIMATE
        );

        if (!forward || !backward)
        {
            throw std::runtime_error(
                "lsdtt::generate_power_law_fractal_surface_fft: FFTW plan creation failed"
            );
        }

        const int threads = resolve_thread_count(n_threads);

        // ---------------------------------------------------------------------
        // 1. Deterministic white noise on the full padded grid
        // ---------------------------------------------------------------------
#ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
        for (int row = 0; row < padded_rows; ++row)
        {
            for (int col = 0; col < padded_cols; ++col)
            {
                const std::size_t idx = flat_index(row, col, padded_cols);

                const std::uint64_t key = make_hash_key(
                    seed,
                    0x504F574C41574646ULL, // "POWLAWFF"
                    static_cast<std::uint64_t>(row),
                    static_cast<std::uint64_t>(col)
                );

                const double z = symmetric_uniform(key);

                in[idx][0] = z;
                in[idx][1] = 0.0;
            }
        }

        // ---------------------------------------------------------------------
        // 2. Forward FFT
        // ---------------------------------------------------------------------
        fftw_execute(forward);

        // ---------------------------------------------------------------------
        // 3. Apply spectral amplitude scaling proportional_to f^(-beta/2)
        //
        // Frequency convention is unshifted FFT order:
        //   fx = col * dfx                    for col <= padded_cols/2
        //      = (col - padded_cols) * dfx    otherwise
        //
        //   fy = row * dfy                    for row <= padded_rows/2
        //      = (row - padded_rows) * dfy    otherwise
        //
        // DC term is set to 0 to avoid singular scaling at f = 0 and to keep
        // the generated field mean-free before final normalisation.
        // ---------------------------------------------------------------------
        const double amplitude_exponent = 0.5 * beta;

#ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
        for (int row = 0; row < padded_rows; ++row)
        {
            for (int col = 0; col < padded_cols; ++col)
            {
                const double fy =
                    (row <= padded_rows / 2)
                        ? static_cast<double>(row) * dfy
                        : static_cast<double>(row - padded_rows) * dfy;

                const double fx =
                    (col <= padded_cols / 2)
                        ? static_cast<double>(col) * dfx
                        : static_cast<double>(col - padded_cols) * dfx;

                const double f = std::sqrt(fx * fx + fy * fy);

                double scale = 0.0;
                if (f > 0.0)
                {
                    scale = std::pow(f, -amplitude_exponent);
                }

                const std::size_t idx = flat_index(row, col, padded_cols);

                out[idx][0] *= scale;
                out[idx][1] *= scale;
            }
        }

        // ---------------------------------------------------------------------
        // 4. Inverse FFT
        // ---------------------------------------------------------------------
        fftw_execute(backward);

        // ---------------------------------------------------------------------
        // 5. Crop back to requested raster size and normalise FFTW inverse by N
        // ---------------------------------------------------------------------
        std::vector<float> cropped(
            static_cast<std::size_t>(n_rows) *
            static_cast<std::size_t>(n_cols),
            0.0f
        );

        const double inv_n = 1.0 / static_cast<double>(n_fft);

#ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
        for (int row = 0; row < n_rows; ++row)
        {
            for (int col = 0; col < n_cols; ++col)
            {
                const std::size_t src_idx = flat_index(row, col, padded_cols);
                const std::size_t dst_idx = flat_index(row, col, n_cols);

                cropped[dst_idx] = static_cast<float>(in[src_idx][0] * inv_n);
            }
        }

        // ---------------------------------------------------------------------
        // 6. Match requested mean and standard deviation
        // ---------------------------------------------------------------------
        normalise_to_mean_std(cropped, target_mean, target_stddev);

        if (forward)  fftw_destroy_plan(forward);
        if (backward) fftw_destroy_plan(backward);
        if (in)       fftw_free(in);
        if (out)      fftw_free(out);

        return make_raster(std::move(cropped), meta);
    }
    catch (...)
    {
        if (forward)  fftw_destroy_plan(forward);
        if (backward) fftw_destroy_plan(backward);
        if (in)       fftw_free(in);
        if (out)      fftw_free(out);
        throw;
    }
#endif
}

lsdtt_raster_cpu generate_power_law_fractal_surface_fft(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed,
    double beta,
    double target_mean,
    double target_stddev,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);

    return generate_power_law_fractal_surface_fft(
        meta,
        seed,
        beta,
        target_mean,
        target_stddev,
        n_threads
    );
}




// -----------------------------------------------------------------------------
// validate_periodic_inputs
// -----------------------------------------------------------------------------
static void validate_periodic_inputs(
    const lsdtt_raster_meta& meta,
    int    n_components,
    double amplitude_min,
    double amplitude_max,
    int    mode_min,
    int    mode_max)
{
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::generate_periodic_wave_test_surface: invalid raster metadata"
        );
    }

    if (n_components < 1)
    {
        throw std::runtime_error(
            "lsdtt::generate_periodic_wave_test_surface: n_components must be >= 1"
        );
    }

    if (!std::isfinite(amplitude_min) ||
        !std::isfinite(amplitude_max) ||
        amplitude_min < 0.0 ||
        amplitude_max <= 0.0 ||
        amplitude_min > amplitude_max)
    {
        throw std::runtime_error(
            "lsdtt::generate_periodic_wave_test_surface: invalid amplitude range"
        );
    }

    if (mode_min < 1 || mode_max < mode_min)
    {
        throw std::runtime_error(
            "lsdtt::generate_periodic_wave_test_surface: mode_min must be >= 1 "
            "and mode_max must be >= mode_min"
        );
    }
}

// -----------------------------------------------------------------------------
// validate_quadratic_inputs
// -----------------------------------------------------------------------------
static void validate_quadratic_inputs(
    const lsdtt_raster_meta& meta,
    const lsdtt_global_quadratic_coefficients& coeffs)
{
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::generate_global_quadratic_test_surface: invalid raster metadata"
        );
    }

    const double vals[6] = {
        coeffs.a, coeffs.b, coeffs.c, coeffs.d, coeffs.e, coeffs.f
    };

    for (double v : vals)
    {
        if (!std::isfinite(v))
        {
            throw std::runtime_error(
                "lsdtt::generate_global_quadratic_test_surface: "
                "all coefficients must be finite"
            );
        }
    }
}

// -----------------------------------------------------------------------------
// draw_mode_component
// -----------------------------------------------------------------------------
static int draw_mode_component(
    std::mt19937& rng,
    int           mode_min,
    int           mode_max)
{
    std::uniform_int_distribution<int> dist(-mode_max, mode_max);

    while (true)
    {
        const int v = dist(rng);

        if (v == 0) return 0;
        if (std::abs(v) >= mode_min) return v;
    }
}

// -----------------------------------------------------------------------------
// draw_symmetric_uniform
// -----------------------------------------------------------------------------
static double draw_symmetric_uniform(
    std::mt19937& rng,
    double        abs_max,
    bool          random_sign)
{
    if (!std::isfinite(abs_max) || abs_max < 0.0)
    {
        throw std::runtime_error(
            "lsdtt synthetic surface generator: invalid coefficient range"
        );
    }

    std::uniform_real_distribution<double> mag_dist(0.0, abs_max);

    double v = mag_dist(rng);

    if (random_sign)
    {
        std::uniform_int_distribution<int> sign_dist(0, 1);
        if (sign_dist(rng) == 0) v = -v;
    }

    return v;
}

// -----------------------------------------------------------------------------
// generate_random_components
// -----------------------------------------------------------------------------
static std::vector<lsdtt_periodic_wave_component> generate_random_components(
    int    n_components,
    int    seed,
    double amplitude_min,
    double amplitude_max,
    int    mode_min,
    int    mode_max,
    bool   random_sign)
{
    std::mt19937 rng(static_cast<std::mt19937::result_type>(seed));

    std::uniform_real_distribution<double> amp_dist(amplitude_min, amplitude_max);
    std::uniform_real_distribution<double> phase_dist(0.0, 2.0 * k_pi);
    std::uniform_int_distribution<int>     sign_dist(0, 1);

    std::vector<lsdtt_periodic_wave_component> components;
    components.reserve(static_cast<std::size_t>(n_components));

    for (int i = 0; i < n_components; ++i)
    {
        int mx = 0;
        int my = 0;

        while (mx == 0 && my == 0)
        {
            mx = draw_mode_component(rng, mode_min, mode_max);
            my = draw_mode_component(rng, mode_min, mode_max);
        }

        double amplitude = amp_dist(rng);
        if (random_sign && sign_dist(rng) == 0)
        {
            amplitude = -amplitude;
        }

        lsdtt_periodic_wave_component c;
        c.mode_x    = mx;
        c.mode_y    = my;
        c.amplitude = amplitude;
        c.phase     = phase_dist(rng);

        components.push_back(c);
    }

    return components;
}

// -----------------------------------------------------------------------------
// generate_periodic_wave_test_surface
// -----------------------------------------------------------------------------
lsdtt_periodic_wave_surface generate_periodic_wave_test_surface(
    const lsdtt_raster_meta& meta,
    int    n_components,
    int    seed,
    double amplitude_min,
    double amplitude_max,
    int    mode_min,
    int    mode_max,
    bool   random_sign,
    int    n_threads)
{
    validate_periodic_inputs(
        meta,
        n_components,
        amplitude_min,
        amplitude_max,
        mode_min,
        mode_max
    );

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;

    const double Lx = static_cast<double>(n_cols) * cell_size;
    const double Ly = static_cast<double>(n_rows) * cell_size;

    if (!std::isfinite(Lx) || !std::isfinite(Ly) || Lx <= 0.0 || Ly <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::generate_periodic_wave_test_surface: invalid raster extents"
        );
    }

    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    std::vector<float> z_data(n_pixels, 0.0f);
    std::vector<float> a_data(n_pixels, 0.0f);
    std::vector<float> b_data(n_pixels, 0.0f);
    std::vector<float> c_data(n_pixels, 0.0f);
    std::vector<float> d_data(n_pixels, 0.0f);
    std::vector<float> e_data(n_pixels, 0.0f);
    std::vector<float> f_data(n_pixels, 0.0f);

    std::vector<lsdtt_periodic_wave_component> components =
        generate_random_components(
            n_components,
            seed,
            amplitude_min,
            amplitude_max,
            mode_min,
            mode_max,
            random_sign
        );

    std::mt19937 offset_rng(static_cast<std::mt19937::result_type>(seed + 7919));
    std::uniform_real_distribution<double> offset_dist(-amplitude_max, amplitude_max);
    const double offset = offset_dist(offset_rng);

    const int threads = resolve_thread_count(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx = flat_index(row, col, n_cols);

            // Pixel-centre coordinates measured from the lower-left corner
            // of the raster domain. This makes the constructed surface periodic
            // over [0, Lx] x [0, Ly].
            const double x =
                (static_cast<double>(col) + 0.5) * cell_size;

            const double y =
                (static_cast<double>(n_rows - row) - 0.5) * cell_size;

            double z    = offset;
            double z_x  = 0.0;
            double z_y  = 0.0;
            double z_xx = 0.0;
            double z_yy = 0.0;
            double z_xy = 0.0;

            for (const auto& comp : components)
            {
                const double alpha =
                    2.0 * k_pi * static_cast<double>(comp.mode_x) / Lx;

                const double beta =
                    2.0 * k_pi * static_cast<double>(comp.mode_y) / Ly;

                const double theta = alpha * x + beta * y + comp.phase;

                const double cos_t = std::cos(theta);
                const double sin_t = std::sin(theta);

                const double term  = comp.amplitude * cos_t;
                const double d1    = -comp.amplitude * sin_t;
                const double d2    = -comp.amplitude * cos_t;

                z    += term;
                z_x  += d1 * alpha;
                z_y  += d1 * beta;
                z_xx += d2 * alpha * alpha;
                z_yy += d2 * beta  * beta;
                z_xy += d2 * alpha * beta;
            }

            const double a = 0.5 * z_xx;
            const double b = 0.5 * z_yy;
            const double c = z_xy;
            const double d = z_x;
            const double e = z_y;
            const double f = z;

            z_data[idx] = static_cast<float>(z);
            a_data[idx] = static_cast<float>(a);
            b_data[idx] = static_cast<float>(b);
            c_data[idx] = static_cast<float>(c);
            d_data[idx] = static_cast<float>(d);
            e_data[idx] = static_cast<float>(e);
            f_data[idx] = static_cast<float>(f);
        }
    }

    lsdtt_periodic_wave_surface out;
    out.dem = make_raster(std::move(z_data), meta);
    out.analytic_taylor = poly_fit_coeffs
    {
        make_raster(std::move(a_data), meta),
        make_raster(std::move(b_data), meta),
        make_raster(std::move(c_data), meta),
        make_raster(std::move(d_data), meta),
        make_raster(std::move(e_data), meta),
        make_raster(std::move(f_data), meta)
    };
    out.components = std::move(components);
    out.offset     = offset;

    return out;
}

// -----------------------------------------------------------------------------
// generate_periodic_wave_test_surface convenience overload
// -----------------------------------------------------------------------------
lsdtt_periodic_wave_surface generate_periodic_wave_test_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    n_components,
    int    seed,
    double amplitude_min,
    double amplitude_max,
    int    mode_min,
    int    mode_max,
    bool   random_sign,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);

    return generate_periodic_wave_test_surface(
        meta,
        n_components,
        seed,
        amplitude_min,
        amplitude_max,
        mode_min,
        mode_max,
        random_sign,
        n_threads
    );
}

// -----------------------------------------------------------------------------
// generate_global_quadratic_test_surface
// -----------------------------------------------------------------------------
lsdtt_global_quadratic_surface generate_global_quadratic_test_surface(
    const lsdtt_raster_meta&                  meta,
    const lsdtt_global_quadratic_coefficients& coeffs,
    int                                       n_threads)
{
    validate_quadratic_inputs(meta, coeffs);

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;

    const double Lx = static_cast<double>(n_cols) * cell_size;
    const double Ly = static_cast<double>(n_rows) * cell_size;

    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    std::vector<float> z_data(n_pixels, 0.0f);
    std::vector<float> a_data(n_pixels, static_cast<float>(coeffs.a));
    std::vector<float> b_data(n_pixels, static_cast<float>(coeffs.b));
    std::vector<float> c_data(n_pixels, static_cast<float>(coeffs.c));
    std::vector<float> d_data(n_pixels, 0.0f);
    std::vector<float> e_data(n_pixels, 0.0f);
    std::vector<float> f_data(n_pixels, 0.0f);

    const int threads = resolve_thread_count(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx = flat_index(row, col, n_cols);

            // Coordinates centred on the raster centre to keep values moderate.
            const double x =
                (static_cast<double>(col) + 0.5) * cell_size - 0.5 * Lx;

            const double y =
                (static_cast<double>(n_rows - row) - 0.5) * cell_size - 0.5 * Ly;

            const double z =
                  coeffs.a * x * x
                + coeffs.b * y * y
                + coeffs.c * x * y
                + coeffs.d * x
                + coeffs.e * y
                + coeffs.f;

            const double z_x = 2.0 * coeffs.a * x + coeffs.c * y + coeffs.d;
            const double z_y = coeffs.c * x + 2.0 * coeffs.b * y + coeffs.e;

            z_data[idx] = static_cast<float>(z);
            d_data[idx] = static_cast<float>(z_x);
            e_data[idx] = static_cast<float>(z_y);
            f_data[idx] = static_cast<float>(z);
        }
    }

    lsdtt_global_quadratic_surface out;
    out.dem = make_raster(std::move(z_data), meta);
    out.analytic_taylor = poly_fit_coeffs
    {
        make_raster(std::move(a_data), meta),
        make_raster(std::move(b_data), meta),
        make_raster(std::move(c_data), meta),
        make_raster(std::move(d_data), meta),
        make_raster(std::move(e_data), meta),
        make_raster(std::move(f_data), meta)
    };
    out.global_coeffs = coeffs;

    return out;
}

// -----------------------------------------------------------------------------
// generate_global_quadratic_test_surface convenience overload
// -----------------------------------------------------------------------------
lsdtt_global_quadratic_surface generate_global_quadratic_test_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    const lsdtt_global_quadratic_coefficients& coeffs,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);

    return generate_global_quadratic_test_surface(
        meta,
        coeffs,
        n_threads
    );
}

// -----------------------------------------------------------------------------
// generate_random_quadratic_test_surface
// -----------------------------------------------------------------------------
lsdtt_global_quadratic_surface generate_random_quadratic_test_surface(
    const lsdtt_raster_meta& meta,
    int    seed,
    double quadratic_abs_max,
    double linear_abs_max,
    double constant_abs_max,
    bool   random_sign,
    int    n_threads)
{
    if (!std::isfinite(quadratic_abs_max) || quadratic_abs_max < 0.0 ||
        !std::isfinite(linear_abs_max)    || linear_abs_max    < 0.0 ||
        !std::isfinite(constant_abs_max)  || constant_abs_max  < 0.0)
    {
        throw std::runtime_error(
            "lsdtt::generate_random_quadratic_test_surface: invalid coefficient ranges"
        );
    }

    std::mt19937 rng(static_cast<std::mt19937::result_type>(seed));

    lsdtt_global_quadratic_coefficients coeffs;
    coeffs.a = draw_symmetric_uniform(rng, quadratic_abs_max, random_sign);
    coeffs.b = draw_symmetric_uniform(rng, quadratic_abs_max, random_sign);
    coeffs.c = draw_symmetric_uniform(rng, quadratic_abs_max, random_sign);
    coeffs.d = draw_symmetric_uniform(rng, linear_abs_max,    random_sign);
    coeffs.e = draw_symmetric_uniform(rng, linear_abs_max,    random_sign);
    coeffs.f = draw_symmetric_uniform(rng, constant_abs_max,  random_sign);

    return generate_global_quadratic_test_surface(
        meta,
        coeffs,
        n_threads
    );
}

// -----------------------------------------------------------------------------
// generate_random_quadratic_test_surface convenience overload
// -----------------------------------------------------------------------------
lsdtt_global_quadratic_surface generate_random_quadratic_test_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed,
    double quadratic_abs_max,
    double linear_abs_max,
    double constant_abs_max,
    bool   random_sign,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);

    return generate_random_quadratic_test_surface(
        meta,
        seed,
        quadratic_abs_max,
        linear_abs_max,
        constant_abs_max,
        random_sign,
        n_threads
    );
}


// =============================================================================
// Mountain-range synthetic surface generators
// =============================================================================

// Pixels whose alpha >= this threshold are labelled mountain in the mask.
// exp(-0.5) places the front at the 1-sigma half-width of the ridge Gaussian,
// giving a clean, well-defined boundary that is still gradual enough in elevation
// to make the segmentation problem non-trivial.
static constexpr double k_mountain_mask_threshold = 0.6065306597633104; // exp(-0.5)

// -----------------------------------------------------------------------------
// build_range_surface_from_alpha
//
// Common kernel for all four range generators.  Takes a pre-computed alpha
// field (0..1) and produces the DEM and ground-truth mask.
//
//   z(pixel) = piedmont_elevation
//            + (ridge_elevation - piedmont_elevation) * alpha
//            + diamond_square_noise * noise_amplitude(alpha)
//
//   noise_amplitude(alpha) = piedmont_noise + (mountain_noise - piedmont_noise) * alpha
//
//   mask(pixel) = (alpha >= k_mountain_mask_threshold) ? 1.0 : 0.0
// -----------------------------------------------------------------------------
static lsdtt_mountain_range_surface build_range_surface_from_alpha(
    const lsdtt_raster_meta&  meta,
    const std::vector<float>& alpha_data,
    double ridge_elevation,
    double piedmont_elevation,
    double mountain_noise_amplitude,
    double piedmont_noise_amplitude,
    double noise_roughness,
    int    noise_seed,
    int    n_threads)
{
    const int          n_rows   = meta.n_rows;
    const int          n_cols   = meta.n_cols;
    const std::size_t  n_pixels = static_cast<std::size_t>(n_rows)
                                * static_cast<std::size_t>(n_cols);
    const int          threads  = resolve_thread_count(n_threads);

    // Noise raster normalised to mean=0, stddev=1.
    // Per-pixel amplitude is applied below so mountain and piedmont zones have
    // independently controllable roughness.
    const lsdtt_raster_cpu noise_raster = generate_diamond_square_surface(
        meta,
        noise_seed,
        noise_roughness,
        1.0,   // initial_scale
        0.0,   // target_mean
        1.0,   // target_stddev
        n_threads
    );
    const std::vector<float>& noise = noise_raster.get_data();

    std::vector<float> dem_data (n_pixels, 0.0f);
    std::vector<float> mask_data(n_pixels, 0.0f);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (std::size_t i = 0; i < n_pixels; ++i)
    {
        const double alpha = static_cast<double>(alpha_data[i]);

        const double z_base    = piedmont_elevation
                               + (ridge_elevation - piedmont_elevation) * alpha;

        const double noise_amp = piedmont_noise_amplitude
                               + (mountain_noise_amplitude - piedmont_noise_amplitude) * alpha;

        dem_data[i]  = static_cast<float>(z_base + static_cast<double>(noise[i]) * noise_amp);
        mask_data[i] = (alpha_data[i] >= static_cast<float>(k_mountain_mask_threshold))
                       ? 1.0f : 0.0f;
    }

    lsdtt_mountain_range_surface out;
    out.dem        = make_raster(std::move(dem_data),  meta);
    out.front_mask = make_raster(std::move(mask_data), meta);
    return out;
}

// -----------------------------------------------------------------------------
// validate_range_noise_params  - shared parameter checks
// -----------------------------------------------------------------------------
static void validate_range_noise_params(
    const char*  fn_name,
    double       mountain_noise_amplitude,
    double       piedmont_noise_amplitude,
    double       noise_roughness)
{
    if (mountain_noise_amplitude < 0.0)
    {
        throw std::runtime_error(
            std::string(fn_name) + ": mountain_noise_amplitude must be >= 0"
        );
    }
    if (piedmont_noise_amplitude < 0.0)
    {
        throw std::runtime_error(
            std::string(fn_name) + ": piedmont_noise_amplitude must be >= 0"
        );
    }
    if (!std::isfinite(noise_roughness) || noise_roughness <= 0.0)
    {
        throw std::runtime_error(
            std::string(fn_name) + ": noise_roughness must be finite and > 0"
        );
    }
}

// -----------------------------------------------------------------------------
// generate_linear_range_surface
// -----------------------------------------------------------------------------
lsdtt_mountain_range_surface generate_linear_range_surface(
    const lsdtt_raster_meta& meta,
    double ridge_elevation,
    double piedmont_elevation,
    double ridge_position_fraction,
    double ridge_half_width_fraction,
    double mountain_noise_amplitude,
    double piedmont_noise_amplitude,
    bool   ridge_runs_north_south,
    double noise_roughness,
    int    noise_seed,
    int    n_threads)
{
    static const char* fn = "lsdtt::generate_linear_range_surface";

    if (!meta.is_valid())
    {
        throw std::runtime_error(std::string(fn) + ": invalid raster metadata");
    }
    if (ridge_position_fraction < 0.0 || ridge_position_fraction > 1.0)
    {
        throw std::runtime_error(
            std::string(fn) + ": ridge_position_fraction must be in [0, 1]"
        );
    }
    if (ridge_half_width_fraction <= 0.0 || ridge_half_width_fraction > 0.5)
    {
        throw std::runtime_error(
            std::string(fn) + ": ridge_half_width_fraction must be in (0, 0.5]"
        );
    }
    validate_range_noise_params(fn,
        mountain_noise_amplitude, piedmont_noise_amplitude, noise_roughness);

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;
    const double Lx        = static_cast<double>(n_cols) * cell_size;
    const double Ly        = static_cast<double>(n_rows) * cell_size;
    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    // Cross-axis extent, ridge position, and half-width in physical coordinates
    const double L_cross  = ridge_runs_north_south ? Lx : Ly;
    const double ridge_pos = L_cross * ridge_position_fraction;
    const double ridge_hw  = L_cross * ridge_half_width_fraction;

    const int threads = resolve_thread_count(n_threads);

    std::vector<float> alpha(n_pixels);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double x = (static_cast<double>(col) + 0.5) * cell_size;
            const double y = (static_cast<double>(n_rows - row) - 0.5) * cell_size;

            const double coord = ridge_runs_north_south ? x : y;
            const double t     = (coord - ridge_pos) / ridge_hw;

            alpha[flat_index(row, col, n_cols)] =
                static_cast<float>(std::exp(-0.5 * t * t));
        }
    }

    return build_range_surface_from_alpha(
        meta, alpha,
        ridge_elevation, piedmont_elevation,
        mountain_noise_amplitude, piedmont_noise_amplitude,
        noise_roughness, noise_seed, n_threads
    );
}

lsdtt_mountain_range_surface generate_linear_range_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    double ridge_elevation,
    double piedmont_elevation,
    double ridge_position_fraction,
    double ridge_half_width_fraction,
    double mountain_noise_amplitude,
    double piedmont_noise_amplitude,
    bool   ridge_runs_north_south,
    double noise_roughness,
    int    noise_seed,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);
    return generate_linear_range_surface(
        meta,
        ridge_elevation, piedmont_elevation,
        ridge_position_fraction, ridge_half_width_fraction,
        mountain_noise_amplitude, piedmont_noise_amplitude,
        ridge_runs_north_south,
        noise_roughness, noise_seed, n_threads
    );
}

// -----------------------------------------------------------------------------
// generate_halfspace_range_surface
// -----------------------------------------------------------------------------
lsdtt_mountain_range_surface generate_halfspace_range_surface(
    const lsdtt_raster_meta& meta,
    lsdtt_range_edge edge,
    double ridge_elevation,
    double piedmont_elevation,
    double front_fraction,
    double mountain_noise_amplitude,
    double piedmont_noise_amplitude,
    double noise_roughness,
    int    noise_seed,
    int    n_threads)
{
    static const char* fn = "lsdtt::generate_halfspace_range_surface";

    if (!meta.is_valid())
    {
        throw std::runtime_error(std::string(fn) + ": invalid raster metadata");
    }
    if (front_fraction <= 0.0 || front_fraction >= 1.0)
    {
        throw std::runtime_error(
            std::string(fn) + ": front_fraction must be in (0, 1)"
        );
    }
    validate_range_noise_params(fn,
        mountain_noise_amplitude, piedmont_noise_amplitude, noise_roughness);

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;
    const double Lx        = static_cast<double>(n_cols) * cell_size;
    const double Ly        = static_cast<double>(n_rows) * cell_size;
    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    // Decay length: the exp(-0.5) contour sits at d = decay_length from the edge
    const double L_perp    = (edge == lsdtt_range_edge::north ||
                              edge == lsdtt_range_edge::south) ? Ly : Lx;
    const double decay_len = L_perp * front_fraction;

    const int threads = resolve_thread_count(n_threads);

    std::vector<float> alpha(n_pixels);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double x = (static_cast<double>(col) + 0.5) * cell_size;
            const double y = (static_cast<double>(n_rows - row) - 0.5) * cell_size;

            double d = 0.0;
            switch (edge)
            {
            case lsdtt_range_edge::north: d = Ly - y; break;
            case lsdtt_range_edge::south: d = y;      break;
            case lsdtt_range_edge::east:  d = Lx - x; break;
            case lsdtt_range_edge::west:  d = x;      break;
            }

            const double t = d / decay_len;
            alpha[flat_index(row, col, n_cols)] =
                static_cast<float>(std::exp(-0.5 * t * t));
        }
    }

    return build_range_surface_from_alpha(
        meta, alpha,
        ridge_elevation, piedmont_elevation,
        mountain_noise_amplitude, piedmont_noise_amplitude,
        noise_roughness, noise_seed, n_threads
    );
}

lsdtt_mountain_range_surface generate_halfspace_range_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    lsdtt_range_edge edge,
    double ridge_elevation,
    double piedmont_elevation,
    double front_fraction,
    double mountain_noise_amplitude,
    double piedmont_noise_amplitude,
    double noise_roughness,
    int    noise_seed,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);
    return generate_halfspace_range_surface(
        meta, edge,
        ridge_elevation, piedmont_elevation,
        front_fraction,
        mountain_noise_amplitude, piedmont_noise_amplitude,
        noise_roughness, noise_seed, n_threads
    );
}

// -----------------------------------------------------------------------------
// generate_circular_range_surface
// -----------------------------------------------------------------------------
lsdtt_mountain_range_surface generate_circular_range_surface(
    const lsdtt_raster_meta& meta,
    double ridge_elevation,
    double piedmont_elevation,
    double centre_x_fraction,
    double centre_y_fraction,
    double radius_x_fraction,
    double radius_y_fraction,
    double mountain_noise_amplitude,
    double piedmont_noise_amplitude,
    double noise_roughness,
    int    noise_seed,
    int    n_threads)
{
    static const char* fn = "lsdtt::generate_circular_range_surface";

    if (!meta.is_valid())
    {
        throw std::runtime_error(std::string(fn) + ": invalid raster metadata");
    }
    if (radius_x_fraction <= 0.0 || radius_y_fraction <= 0.0)
    {
        throw std::runtime_error(
            std::string(fn) + ": radius fractions must be > 0"
        );
    }
    validate_range_noise_params(fn,
        mountain_noise_amplitude, piedmont_noise_amplitude, noise_roughness);

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;
    const double Lx        = static_cast<double>(n_cols) * cell_size;
    const double Ly        = static_cast<double>(n_rows) * cell_size;
    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    const double cx = Lx * centre_x_fraction;
    const double cy = Ly * centre_y_fraction;
    const double rx = Lx * radius_x_fraction;
    const double ry = Ly * radius_y_fraction;

    const int threads = resolve_thread_count(n_threads);

    std::vector<float> alpha(n_pixels);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double x  = (static_cast<double>(col) + 0.5) * cell_size;
            const double y  = (static_cast<double>(n_rows - row) - 0.5) * cell_size;
            const double dx = (x - cx) / rx;
            const double dy = (y - cy) / ry;
            const double r2 = dx * dx + dy * dy;

            alpha[flat_index(row, col, n_cols)] =
                static_cast<float>(std::exp(-0.5 * r2));
        }
    }

    return build_range_surface_from_alpha(
        meta, alpha,
        ridge_elevation, piedmont_elevation,
        mountain_noise_amplitude, piedmont_noise_amplitude,
        noise_roughness, noise_seed, n_threads
    );
}

lsdtt_mountain_range_surface generate_circular_range_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    double ridge_elevation,
    double piedmont_elevation,
    double centre_x_fraction,
    double centre_y_fraction,
    double radius_x_fraction,
    double radius_y_fraction,
    double mountain_noise_amplitude,
    double piedmont_noise_amplitude,
    double noise_roughness,
    int    noise_seed,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);
    return generate_circular_range_surface(
        meta,
        ridge_elevation, piedmont_elevation,
        centre_x_fraction, centre_y_fraction,
        radius_x_fraction, radius_y_fraction,
        mountain_noise_amplitude, piedmont_noise_amplitude,
        noise_roughness, noise_seed, n_threads
    );
}

// -----------------------------------------------------------------------------
// generate_arc_range_surface
//
// Alpha is computed using a distance-from-edge Gaussian whose effective decay
// length is modulated by a parabolic lateral weight:
//
//   u            = (parallel_coord - arc_centre) / arc_half_width
//   lw(u)        = max(0, 1 - u^2)          [parabolic lateral weight]
//   eff_decay(u) = decay_length * max(lw, epsilon)
//   alpha        = lw > 0 ? exp(-0.5 * (d / eff_decay)^2) : 0
//
// The mask threshold contour is then the parabola
//   d_front(u) = decay_length * lw(u)
// which gives depth = decay_length at the centre, tapering to 0 at u = +/-1.
// The mountain zone is deepest at the arc centre and shallows toward the ends,
// producing a convex front as seen from the piedmont.
// -----------------------------------------------------------------------------
lsdtt_mountain_range_surface generate_arc_range_surface(
    const lsdtt_raster_meta& meta,
    lsdtt_range_edge edge,
    double ridge_elevation,
    double piedmont_elevation,
    double front_fraction,
    double arc_centre_fraction,
    double arc_half_width_fraction,
    double mountain_noise_amplitude,
    double piedmont_noise_amplitude,
    double noise_roughness,
    int    noise_seed,
    int    n_threads)
{
    static const char* fn = "lsdtt::generate_arc_range_surface";

    if (!meta.is_valid())
    {
        throw std::runtime_error(std::string(fn) + ": invalid raster metadata");
    }
    if (front_fraction <= 0.0 || front_fraction >= 1.0)
    {
        throw std::runtime_error(
            std::string(fn) + ": front_fraction must be in (0, 1)"
        );
    }
    if (arc_half_width_fraction <= 0.0)
    {
        throw std::runtime_error(
            std::string(fn) + ": arc_half_width_fraction must be > 0"
        );
    }
    validate_range_noise_params(fn,
        mountain_noise_amplitude, piedmont_noise_amplitude, noise_roughness);

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;
    const double Lx        = static_cast<double>(n_cols) * cell_size;
    const double Ly        = static_cast<double>(n_rows) * cell_size;
    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    const bool ns_edge = (edge == lsdtt_range_edge::north ||
                          edge == lsdtt_range_edge::south);

    // Perpendicular extent (depth axis) and parallel extent (along-range axis)
    const double L_perp     = ns_edge ? Ly : Lx;
    const double L_par      = ns_edge ? Lx : Ly;
    const double decay_len  = L_perp * front_fraction;
    const double arc_centre = L_par  * arc_centre_fraction;
    const double arc_hw     = L_par  * arc_half_width_fraction;

    // Guard against numerical issues when lw is very small
    static constexpr double k_lw_epsilon = 1e-6;

    const int threads = resolve_thread_count(n_threads);

    std::vector<float> alpha(n_pixels);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double x = (static_cast<double>(col) + 0.5) * cell_size;
            const double y = (static_cast<double>(n_rows - row) - 0.5) * cell_size;

            // Distance from the edge (depth into the DEM)
            double d        = 0.0;
            double par_coord = 0.0;
            switch (edge)
            {
            case lsdtt_range_edge::north: d = Ly - y; par_coord = x; break;
            case lsdtt_range_edge::south: d = y;      par_coord = x; break;
            case lsdtt_range_edge::east:  d = Lx - x; par_coord = y; break;
            case lsdtt_range_edge::west:  d = x;      par_coord = y; break;
            }

            const double u  = (par_coord - arc_centre) / arc_hw;
            const double lw = std::max(0.0, 1.0 - u * u);

            float a = 0.0f;
            if (lw > k_lw_epsilon)
            {
                const double eff_decay = decay_len * lw;
                const double t         = d / eff_decay;
                a = static_cast<float>(lw * std::exp(-0.5 * t * t));
            }

            alpha[flat_index(row, col, n_cols)] = a;
        }
    }

    return build_range_surface_from_alpha(
        meta, alpha,
        ridge_elevation, piedmont_elevation,
        mountain_noise_amplitude, piedmont_noise_amplitude,
        noise_roughness, noise_seed, n_threads
    );
}

lsdtt_mountain_range_surface generate_arc_range_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    lsdtt_range_edge edge,
    double ridge_elevation,
    double piedmont_elevation,
    double front_fraction,
    double arc_centre_fraction,
    double arc_half_width_fraction,
    double mountain_noise_amplitude,
    double piedmont_noise_amplitude,
    double noise_roughness,
    int    noise_seed,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);
    return generate_arc_range_surface(
        meta, edge,
        ridge_elevation, piedmont_elevation,
        front_fraction,
        arc_centre_fraction, arc_half_width_fraction,
        mountain_noise_amplitude, piedmont_noise_amplitude,
        noise_roughness, noise_seed, n_threads
    );
}

// -----------------------------------------------------------------------------
// strip_nonlinear_slope
//
// Roering (2007) steady-state slope magnitude S at distance x from the divide,
// the positive root of  D*S/(1 - (S/Sc)^2) = U*x, i.e.
//   (U*x/Sc^2) S^2 + D S - U*x = 0.
// Reduces to S = U*x/D as x -> 0 (the linear limit) and approaches Sc from
// below as x grows.
// -----------------------------------------------------------------------------
static double strip_nonlinear_slope(double x, double U, double D, double Sc)
{
    const double u = U * x;
    if (u <= 0.0)
    {
        return 0.0;
    }
    const double a = u / (Sc * Sc);
    // S = (-D + sqrt(D^2 + 4 a u)) / (2 a)
    return (-D + std::sqrt(D * D + 4.0 * a * u)) / (2.0 * a);
}

// -----------------------------------------------------------------------------
// generate_hillslope_strip_surface
// -----------------------------------------------------------------------------
lsdtt_hillslope_strip_surface generate_hillslope_strip_surface(
    const lsdtt_raster_meta&   meta,
    lsdtt_strip_channel_config channel_config,
    lsdtt_strip_init_profile   init_profile,
    bool                       cross_axis_is_x,
    double                     baselevel_elevation,
    double                     half_width,
    double                     uplift_U,
    double                     diffusivity_D,
    double                     critical_slope_Sc,
    bool                       nonlinear_profile,
    double                     init_amplitude,
    int                        perturb_seed,
    int                        n_threads)
{
    static const char* fn = "lsdtt::generate_hillslope_strip_surface";

    if (!meta.is_valid())
    {
        throw std::runtime_error(std::string(fn) + ": invalid raster metadata");
    }
    if (diffusivity_D <= 0.0)
    {
        throw std::runtime_error(std::string(fn) + ": diffusivity_D must be > 0");
    }
    if (critical_slope_Sc <= 0.0)
    {
        throw std::runtime_error(std::string(fn) + ": critical_slope_Sc must be > 0");
    }

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;

    // Number of cells along the cross-axis (the direction elevation varies in).
    const int n_cross = cross_axis_is_x ? n_cols : n_rows;
    if (n_cross < 3)
    {
        throw std::runtime_error(
            std::string(fn) + ": cross-axis must span at least 3 cells");
    }

    // -------------------------------------------------------------------------
    // Per cross-index: distance-from-divide xd, channel flag, and the
    // auto-derived half-width L (divide-to-channel distance).
    // -------------------------------------------------------------------------
    std::vector<double>  xd(static_cast<std::size_t>(n_cross), 0.0);
    std::vector<uint8_t> chan(static_cast<std::size_t>(n_cross), uint8_t{0});
    double L_auto = 0.0;

    switch (channel_config)
    {
        case lsdtt_strip_channel_config::both_ends:
        {
            // Divide at the centre; channels at both extreme cells.
            const double centre = 0.5 * static_cast<double>(n_cross - 1);
            L_auto = centre * cell_size;
            for (int i = 0; i < n_cross; ++i)
            {
                xd[static_cast<std::size_t>(i)] =
                    std::abs(static_cast<double>(i) - centre) * cell_size;
            }
            chan[0]                                          = uint8_t{1};
            chan[static_cast<std::size_t>(n_cross - 1)]      = uint8_t{1};
            break;
        }
        case lsdtt_strip_channel_config::single_end_low:
        {
            // Channel at index 0; divide (no-flux) at the high-index end.
            L_auto = static_cast<double>(n_cross - 1) * cell_size;
            for (int i = 0; i < n_cross; ++i)
            {
                xd[static_cast<std::size_t>(i)] =
                    static_cast<double>((n_cross - 1) - i) * cell_size;
            }
            chan[0] = uint8_t{1};
            break;
        }
        case lsdtt_strip_channel_config::single_end_high:
        {
            // Channel at the high-index end; divide (no-flux) at index 0.
            L_auto = static_cast<double>(n_cross - 1) * cell_size;
            for (int i = 0; i < n_cross; ++i)
            {
                xd[static_cast<std::size_t>(i)] =
                    static_cast<double>(i) * cell_size;
            }
            chan[static_cast<std::size_t>(n_cross - 1)] = uint8_t{1};
            break;
        }
    }

    // Half-width used by the analytic profiles: auto-derived unless overridden.
    const double L = (half_width > 0.0) ? half_width : L_auto;

    // -------------------------------------------------------------------------
    // Build the 1-D profile along the cross-axis (one value per cross-index).
    // Channel cells are forced to baselevel below.
    // -------------------------------------------------------------------------
    std::vector<double> profile(static_cast<std::size_t>(n_cross), baselevel_elevation);

    switch (init_profile)
    {
        case lsdtt_strip_init_profile::flat:
        {
            // already baselevel everywhere
            break;
        }
        case lsdtt_strip_init_profile::analytic_parabola:
        {
            if (!nonlinear_profile)
            {
                // Linear steady parabola: z = z_ch + (U/2D)(L^2 - xd^2)
                const double k = uplift_U / (2.0 * diffusivity_D);
                for (int i = 0; i < n_cross; ++i)
                {
                    const double x = xd[static_cast<std::size_t>(i)];
                    profile[static_cast<std::size_t>(i)] =
                        baselevel_elevation + k * (L * L - x * x);
                }
            }
            else
            {
                // Nonlinear Roering steady profile: integrate slope from the
                // channel (xd = L) back to the divide (xd = 0).  Trapezoidal
                // integration on a fine sub-grid keyed by distance-from-divide.
                const int n_sub = 4096;
                const double dx_sub = L / static_cast<double>(n_sub);
                std::vector<double> z_of_s(static_cast<std::size_t>(n_sub + 1), 0.0);
                // z at the channel (s = L) is baselevel; accumulate upslope.
                z_of_s[static_cast<std::size_t>(n_sub)] = baselevel_elevation;
                for (int s = n_sub - 1; s >= 0; --s)
                {
                    const double x_hi = static_cast<double>(s + 1) * dx_sub;
                    const double x_lo = static_cast<double>(s)     * dx_sub;
                    const double S_hi = strip_nonlinear_slope(
                        x_hi, uplift_U, diffusivity_D, critical_slope_Sc);
                    const double S_lo = strip_nonlinear_slope(
                        x_lo, uplift_U, diffusivity_D, critical_slope_Sc);
                    // dz = +S dx going upslope (decreasing x toward divide)
                    z_of_s[static_cast<std::size_t>(s)] =
                        z_of_s[static_cast<std::size_t>(s + 1)]
                        + 0.5 * (S_hi + S_lo) * dx_sub;
                }
                for (int i = 0; i < n_cross; ++i)
                {
                    double x = xd[static_cast<std::size_t>(i)];
                    if (x > L) x = L;
                    const double f   = x / dx_sub;
                    int    s0 = static_cast<int>(f);
                    if (s0 < 0)      s0 = 0;
                    if (s0 > n_sub)  s0 = n_sub;
                    int s1 = (s0 < n_sub) ? s0 + 1 : n_sub;
                    const double w = f - static_cast<double>(s0);
                    profile[static_cast<std::size_t>(i)] =
                        (1.0 - w) * z_of_s[static_cast<std::size_t>(s0)]
                        +       w * z_of_s[static_cast<std::size_t>(s1)];
                }
            }
            break;
        }
        case lsdtt_strip_init_profile::raised_plateau:
        {
            for (int i = 0; i < n_cross; ++i)
            {
                profile[static_cast<std::size_t>(i)] =
                    baselevel_elevation + init_amplitude;
            }
            break;
        }
        case lsdtt_strip_init_profile::perturbed:
        {
            for (int i = 0; i < n_cross; ++i)
            {
                const std::uint64_t key = make_hash_key(
                    perturb_seed, 0xB17, static_cast<std::uint64_t>(i));
                profile[static_cast<std::size_t>(i)] =
                    baselevel_elevation + init_amplitude * symmetric_uniform(key);
            }
            break;
        }
        case lsdtt_strip_init_profile::linear_ramp:
        {
            for (int i = 0; i < n_cross; ++i)
            {
                const double x = xd[static_cast<std::size_t>(i)];
                const double frac = (L > 0.0) ? (1.0 - x / L) : 0.0;
                profile[static_cast<std::size_t>(i)] =
                    baselevel_elevation + init_amplitude * frac;
            }
            break;
        }
    }

    // -------------------------------------------------------------------------
    // Expand the 1-D profile across the uniform axis and stamp the channel mask.
    // -------------------------------------------------------------------------
    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    std::vector<float> dem_data(n_pixels, 0.0f);
    std::vector<float> mask_data(n_pixels, 0.0f);

    const int threads = resolve_thread_count(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const int cross_i = cross_axis_is_x ? col : row;
            const std::size_t sc = static_cast<std::size_t>(cross_i);
            const std::size_t flat = flat_index(row, col, n_cols);

            const bool is_chan = (chan[sc] != 0);
            dem_data[flat] = static_cast<float>(
                is_chan ? baselevel_elevation : profile[sc]);
            mask_data[flat] = is_chan ? 1.0f : 0.0f;
        }
    }

    lsdtt_hillslope_strip_surface out;
    out.dem          = make_raster(std::move(dem_data),  meta);
    out.channel_mask = make_raster(std::move(mask_data), meta);
    return out;
}

lsdtt_hillslope_strip_surface generate_hillslope_strip_surface(
    int                        n_rows,
    int                        n_cols,
    double                     cell_size,
    lsdtt_strip_channel_config channel_config,
    lsdtt_strip_init_profile   init_profile,
    bool                       cross_axis_is_x,
    double                     baselevel_elevation,
    double                     half_width,
    double                     uplift_U,
    double                     diffusivity_D,
    double                     critical_slope_Sc,
    bool                       nonlinear_profile,
    double                     init_amplitude,
    int                        perturb_seed,
    int                        n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);
    return generate_hillslope_strip_surface(
        meta, channel_config, init_profile, cross_axis_is_x,
        baselevel_elevation, half_width,
        uplift_U, diffusivity_D, critical_slope_Sc,
        nonlinear_profile, init_amplitude, perturb_seed, n_threads);
}

// =============================================================================
// Crater field
//
// A tilted plane with non-overlapping closed depressions cut into it. Every
// reported quantity is derived from the crater specs, so the surface can be
// used as an oracle for depression detection and for sink-terminated filling.
// See the header for the surface definition and for the tolerances a caller
// should expect against a raster measurement.
// =============================================================================

// -----------------------------------------------------------------------------
// crater_geometry
//
// The per-crater scalars that both the surface evaluation and the analytic
// integrals need, derived once from the spec.
// -----------------------------------------------------------------------------
struct crater_geometry
{
    double cx     = 0.0;   // centre easting
    double cy     = 0.0;   // centre northing
    double R      = 0.0;   // rim radius
    double F      = 0.0;   // floor radius
    double D      = 0.0;   // R - F, the width of the sloping wall
    double depth  = 0.0;
    double s      = 0.0;   // background slope, m/m, rising north
    double z_c    = 0.0;   // background elevation at the centre
    double z_pit  = 0.0;   // z_c - depth
    double z_spill= 0.0;   // z_c - s * R, the lowest point of the rim
};

// -----------------------------------------------------------------------------
// crater_boundary_radius
//
// Distance from the crater centre, along the ray whose northward direction
// cosine is uy, at which the surface reaches the spill elevation. This is the
// boundary of the closed depression: inside it water is retained, outside it
// water has already escaped over the rim.
//
// The surface on the wall is
//   z(d) = z_c + a*d*(d - F)/D - depth*(R - d)/D,   a = s * uy
// so setting z(d) = z_spill = z_c - s*R gives
//   a*d^2 + (depth - a*F)*d + (s*R*D - depth*R) = 0
// which has exactly one root in [F, R] because z is strictly increasing there
// (guaranteed by the depth constraint validated in the generator).
//
// With a == 0 -- an untilted background, or the two rays running due east and
// due west -- the quadratic degenerates to a linear equation, so the stable
// citardauq form is used rather than the schoolbook one.
// -----------------------------------------------------------------------------
static double crater_boundary_radius(const crater_geometry& g, double uy)
{
    const double a = g.s * uy;
    const double B = g.depth - a * g.F;
    const double C = g.s * g.R * g.D - g.depth * g.R;

    double d = 0.0;

    if (std::abs(a) < 1e-300)
    {
        d = -C / B;
    }
    else
    {
        const double disc = B * B - 4.0 * a * C;
        const double root = std::sqrt(disc > 0.0 ? disc : 0.0);

        // q is the larger-magnitude of the two numerators, so q/a and C/q are
        // the two roots with no cancellation in either.
        const double q  = -0.5 * (B + (B >= 0.0 ? root : -root));
        const double r1 = q / a;
        const double r2 = (std::abs(q) > 0.0) ? (C / q) : r1;

        const bool r1_ok = (r1 >= g.F - 1e-9) && (r1 <= g.R + 1e-9);
        d = r1_ok ? r1 : r2;
    }

    if (d < g.F) d = g.F;
    if (d > g.R) d = g.R;
    return d;
}

// -----------------------------------------------------------------------------
// crater_wall_volume_antiderivative
//
// Antiderivative in d of h(d) * d over the sloping wall, where
//   h(d) = z_spill - z(d) = -s*R - a*d*(d - F)/D + depth*(R - d)/D
// is the depth of standing water. Multiplying by d is the polar area element,
// so integrating this radially and then over azimuth gives the fill volume.
// -----------------------------------------------------------------------------
static double crater_wall_volume_antiderivative(
    const crater_geometry& g,
    double a,
    double d)
{
    const double d2 = d * d;
    const double d3 = d2 * d;
    const double d4 = d2 * d2;

    return -g.s * g.R * d2 * 0.5
         - (a / g.D) * (0.25 * d4 - g.F * d3 / 3.0)
         + (g.depth / g.D) * (0.5 * g.R * d2 - d3 / 3.0);
}

// -----------------------------------------------------------------------------
// crater_area_and_volume
//
// Planform area and fill volume of the closed depression, both exact.
//
// With no background tilt the depression boundary is the rim circle and both
// integrals close in elementary form. With a tilt the boundary is no longer a
// circle -- the depression is cut short on the down-tilt side -- so the radial
// integral (which is still closed form) is swept over azimuth with a uniform
// trapezoid rule. On a smooth periodic integrand the trapezoid rule converges
// far faster than its formal second order, and the sample count below puts the
// result well inside 1e-12 relative. Nothing here touches the raster, so the
// answer is independent of cell size.
// -----------------------------------------------------------------------------
static void crater_area_and_volume(
    const crater_geometry& g,
    double& area_m2,
    double& volume_m3)
{
    const double H = g.z_spill - g.z_pit;   // water depth over the flat floor

    if (g.s == 0.0)
    {
        area_m2 = k_pi * g.R * g.R;
        volume_m3 = 2.0 * k_pi
            * ( 0.5 * H * g.F * g.F
              + crater_wall_volume_antiderivative(g, 0.0, g.R)
              - crater_wall_volume_antiderivative(g, 0.0, g.F) );
        return;
    }

    constexpr int n_theta = 131072;

    double area_acc   = 0.0;
    double volume_acc = 0.0;

    for (int i = 0; i < n_theta; ++i)
    {
        const double theta = 2.0 * k_pi * static_cast<double>(i)
                           / static_cast<double>(n_theta);
        const double uy = std::sin(theta);
        const double a  = g.s * uy;
        const double r  = crater_boundary_radius(g, uy);

        area_acc += 0.5 * r * r;
        volume_acc += 0.5 * H * g.F * g.F
                    + crater_wall_volume_antiderivative(g, a, r)
                    - crater_wall_volume_antiderivative(g, a, g.F);
    }

    const double dtheta = 2.0 * k_pi / static_cast<double>(n_theta);
    area_m2   = area_acc   * dtheta;
    volume_m3 = volume_acc * dtheta;
}

// -----------------------------------------------------------------------------
// generate_crater_field_surface
// -----------------------------------------------------------------------------
lsdtt_crater_field_surface generate_crater_field_surface(
    const lsdtt_raster_meta&              meta,
    const std::vector<lsdtt_crater_spec>& craters,
    double background_slope,
    double background_elevation)
{
    static const char* fn = "lsdtt::generate_crater_field_surface";

    if (!meta.is_valid())
    {
        throw std::runtime_error(std::string(fn) + ": invalid raster metadata");
    }
    if (!std::isfinite(background_slope) || background_slope < 0.0)
    {
        throw std::invalid_argument(
            std::string(fn) + ": background_slope must be finite and >= 0"
        );
    }
    if (!std::isfinite(background_elevation))
    {
        throw std::invalid_argument(
            std::string(fn) + ": background_elevation must be finite"
        );
    }

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;
    const double x_min     = meta.x_min;
    const double y_min     = meta.y_min;
    const double x_max     = meta.x_max();
    const double y_max     = meta.y_max();
    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    // -------------------------------------------------------------------------
    // Validate every spec before touching the raster, so a bad field never
    // produces a half-built surface.
    // -------------------------------------------------------------------------
    std::vector<crater_geometry> geom;
    geom.reserve(craters.size());

    for (std::size_t i = 0; i < craters.size(); ++i)
    {
        const lsdtt_crater_spec& cs = craters[i];
        const std::string tag = " (crater " + std::to_string(i) + ")";

        if (!std::isfinite(cs.centre_x) || !std::isfinite(cs.centre_y))
        {
            throw std::invalid_argument(
                std::string(fn) + ": centre_x and centre_y must be finite" + tag
            );
        }
        if (!std::isfinite(cs.radius_m) || cs.radius_m <= 0.0)
        {
            throw std::invalid_argument(
                std::string(fn) + ": radius_m must be finite and > 0" + tag
            );
        }
        if (!std::isfinite(cs.depth_m) || cs.depth_m <= 0.0)
        {
            throw std::invalid_argument(
                std::string(fn) + ": depth_m must be finite and > 0" + tag
            );
        }
        if (!std::isfinite(cs.floor_fraction) ||
            cs.floor_fraction < 0.0 || cs.floor_fraction >= 1.0)
        {
            throw std::invalid_argument(
                std::string(fn) + ": floor_fraction must be in [0, 1)" + tag
            );
        }

        // The rim must sit inside the grid with a cell of margin, otherwise the
        // depression opens onto the DEM edge and is not closed at all.
        if (cs.centre_x - cs.radius_m < x_min + cell_size ||
            cs.centre_x + cs.radius_m > x_max - cell_size ||
            cs.centre_y - cs.radius_m < y_min + cell_size ||
            cs.centre_y + cs.radius_m > y_max - cell_size)
        {
            throw std::invalid_argument(
                std::string(fn) + ": the rim runs off the raster, so the "
                "depression would not be closed" + tag
            );
        }

        // The wall must outrun the regional tilt along every ray, including the
        // down-tilt one, or there is no closed depression to find.
        const double tilt_budget =
            background_slope * cs.radius_m * (2.0 - cs.floor_fraction);
        if (cs.depth_m <= tilt_budget)
        {
            throw std::invalid_argument(
                std::string(fn) + ": depth_m (" + std::to_string(cs.depth_m) +
                ") must exceed background_slope * radius_m * "
                "(2 - floor_fraction) = " + std::to_string(tilt_budget) +
                ", otherwise the regional tilt outruns the crater wall on the "
                "down-tilt side and the depression is not closed" + tag
            );
        }

        crater_geometry g;
        g.cx      = cs.centre_x;
        g.cy      = cs.centre_y;
        g.R       = cs.radius_m;
        g.F       = cs.floor_fraction * cs.radius_m;
        g.D       = g.R - g.F;
        g.depth   = cs.depth_m;
        g.s       = background_slope;
        g.z_c     = background_elevation + background_slope * (g.cy - y_min);
        g.z_pit   = g.z_c - g.depth;
        g.z_spill = g.z_c - background_slope * g.R;
        geom.push_back(g);
    }

    for (std::size_t i = 0; i < geom.size(); ++i)
    {
        for (std::size_t j = i + 1; j < geom.size(); ++j)
        {
            const double dx = geom[i].cx - geom[j].cx;
            const double dy = geom[i].cy - geom[j].cy;
            const double sep = std::sqrt(dx * dx + dy * dy);
            if (sep <= geom[i].R + geom[j].R)
            {
                throw std::invalid_argument(
                    std::string(fn) + ": craters " + std::to_string(i) +
                    " and " + std::to_string(j) + " overlap (centres " +
                    std::to_string(sep) + " m apart, rims sum to " +
                    std::to_string(geom[i].R + geom[j].R) +
                    " m). Overlapping rims merge into one depression whose pit "
                    "and spill point neither spec describes."
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // Background plane, rising to the north.
    // -------------------------------------------------------------------------
    std::vector<float> dem(n_pixels, 0.0f);

    for (int row = 0; row < n_rows; ++row)
    {
        const double y = y_max - (static_cast<double>(row) + 0.5) * cell_size;
        const float  z = static_cast<float>(
            background_elevation + background_slope * (y - y_min));

        for (int col = 0; col < n_cols; ++col)
        {
            dem[flat_index(row, col, n_cols)] = z;
        }
    }

    // -------------------------------------------------------------------------
    // Cut each crater into its own bounding box. Rims cannot overlap, so no
    // cell is written twice.
    // -------------------------------------------------------------------------
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   id_data(n_pixels, nodata);
    std::vector<uint8_t> id_mask(n_pixels, uint8_t{0});
    std::vector<float>   sink_data(n_pixels, 0.0f);

    lsdtt_crater_field_surface out{
        lsdtt_raster_cpu(), lsdtt_raster_cpu(), lsdtt_raster_cpu(),
        {}, {}, {}, {}, {}, {}
    };

    for (std::size_t i = 0; i < geom.size(); ++i)
    {
        const crater_geometry& g = geom[i];

        const int col_lo = std::max(0,
            static_cast<int>(std::floor((g.cx - g.R - x_min) / cell_size)));
        const int col_hi = std::min(n_cols - 1,
            static_cast<int>(std::ceil((g.cx + g.R - x_min) / cell_size)));
        const int row_lo = std::max(0,
            static_cast<int>(std::floor((y_max - (g.cy + g.R)) / cell_size)));
        const int row_hi = std::min(n_rows - 1,
            static_cast<int>(std::ceil((y_max - (g.cy - g.R)) / cell_size)));

        for (int row = row_lo; row <= row_hi; ++row)
        {
            const double y = y_max - (static_cast<double>(row) + 0.5) * cell_size;

            for (int col = col_lo; col <= col_hi; ++col)
            {
                const double x =
                    x_min + (static_cast<double>(col) + 0.5) * cell_size;

                const double dx = x - g.cx;
                const double dy = y - g.cy;
                const double d  = std::sqrt(dx * dx + dy * dy);
                if (d >= g.R) continue;

                const double w = (d <= g.F) ? 0.0 : (d - g.F) / g.D;
                const double z_bg =
                    background_elevation + background_slope * (y - y_min);

                const std::size_t k = flat_index(row, col, n_cols);
                dem[k] = static_cast<float>(
                    g.z_c + (z_bg - g.z_c) * w - g.depth * (1.0 - w));

                id_data[k] = static_cast<float>(i);
                id_mask[k] = uint8_t{1};
            }
        }

        // The pit is the cell containing the analytic apex. On a square grid
        // the containing cell is also the cell whose centre is nearest the
        // apex, so this is the raster minimum for any centre at or near a
        // pixel centre -- see the header note.
        int pit_row = static_cast<int>(std::floor((y_max - g.cy) / cell_size));
        int pit_col = static_cast<int>(std::floor((g.cx - x_min) / cell_size));
        if (pit_row < 0) pit_row = 0;
        if (pit_col < 0) pit_col = 0;
        if (pit_row > n_rows - 1) pit_row = n_rows - 1;
        if (pit_col > n_cols - 1) pit_col = n_cols - 1;

        const std::size_t pit_k = flat_index(pit_row, pit_col, n_cols);
        sink_data[pit_k] = 1.0f;

        double area_m2   = 0.0;
        double volume_m3 = 0.0;
        crater_area_and_volume(g, area_m2, volume_m3);

        out.pit_rc.push_back(lsdtt_pixel_rc{ pit_row, pit_col });
        out.pit_flats.push_back(static_cast<std::int64_t>(pit_k));
        out.pit_elevation.push_back(g.z_pit);
        out.spill_elevation.push_back(g.z_spill);
        out.area_m2.push_back(area_m2);
        out.volume_m3.push_back(volume_m3);
    }

    out.dem       = make_raster(std::move(dem), meta);
    out.sink_mask = make_raster(std::move(sink_data), meta);
    out.crater_id = lsdtt_raster_cpu(
        std::move(id_data), std::move(id_mask), meta);

    return out;
}

lsdtt_crater_field_surface generate_crater_field_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    const std::vector<lsdtt_crater_spec>& craters,
    double background_slope,
    double background_elevation)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);
    return generate_crater_field_surface(
        meta, craters, background_slope, background_elevation);
}

// =============================================================================
// generate_terrace_staircase_surface
// =============================================================================
lsdtt_terrace_staircase_surface generate_terrace_staircase_surface(
    const lsdtt_raster_meta& meta,
    int    n_terraces,
    double tread_width_m,
    double riser_height_m,
    double valley_gradient,
    double riser_height_downstream_taper,
    double tread_cross_gradient,
    int    n_riser_gaps,
    double riser_gap_width_m,
    double noise_amplitude_m,
    double channel_elevation,
    double valley_wall_gradient,
    int    seed,
    int    n_threads)
{
    static const char* fn = "lsdtt::generate_terrace_staircase_surface";

    if (!meta.is_valid())
    {
        throw std::runtime_error(std::string(fn) + ": invalid raster metadata");
    }
    if (n_terraces < 1)
    {
        throw std::runtime_error(std::string(fn) + ": n_terraces must be >= 1");
    }
    if (!(tread_width_m > 0.0))
    {
        throw std::runtime_error(std::string(fn) + ": tread_width_m must be > 0");
    }
    if (!(riser_height_m > 0.0))
    {
        throw std::runtime_error(std::string(fn) + ": riser_height_m must be > 0");
    }
    if (riser_height_downstream_taper < 0.0 || riser_height_downstream_taper >= 1.0)
    {
        throw std::runtime_error(
            std::string(fn) + ": riser_height_downstream_taper must be in [0, 1)");
    }
    if (n_riser_gaps < 0)
    {
        throw std::runtime_error(std::string(fn) + ": n_riser_gaps must be >= 0");
    }
    if (n_riser_gaps > 0 && !(riser_gap_width_m > 0.0))
    {
        throw std::runtime_error(
            std::string(fn) + ": riser_gap_width_m must be > 0 when n_riser_gaps > 0");
    }

    // The cross-fall must not swallow the riser, or the staircase stops being
    // monotonic and "terrace j" no longer means a height above the channel.
    // Checked against the SMALLEST riser height, which is the downstream end
    // when the taper is in use.
    const double min_riser_height =
        riser_height_m * (1.0 - riser_height_downstream_taper);
    if (tread_cross_gradient < 0.0 ||
        tread_cross_gradient * tread_width_m >= min_riser_height)
    {
        throw std::runtime_error(
            std::string(fn) + ": tread_cross_gradient * tread_width_m must be "
            ">= 0 and < the smallest riser height (" +
            std::to_string(min_riser_height) + " m); otherwise the tread "
            "cross-fall exceeds the riser it steps over");
    }

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;

    if (n_rows < 3 || n_cols < 3)
    {
        throw std::runtime_error(
            std::string(fn) + ": grid must span at least 3x3 cells");
    }

    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    std::vector<float> dem_data(n_pixels, 0.0f);
    std::vector<float> tread_data(n_pixels, 0.0f);
    std::vector<float> chan_data(n_pixels, 0.0f);

    // The channel occupies ONE WHOLE COLUMN, so the centre must be an integer
    // column index. A fractional centre (0.5 * (n_cols - 1)) puts the axis
    // between two cells on an even-width grid, and then no cell has x = 0 and
    // the channel mask comes out empty -- taking every HAND-based test with it.
    // On an even-width grid this leaves one more column on the east flank,
    // which is immaterial to a symmetric flight.
    const int    centre_col_i = (n_cols - 1) / 2;
    const double centre_col   = static_cast<double>(centre_col_i);

    // Gap band centres, evenly spaced down the valley and kept off both ends so
    // a gap is never clipped by the grid edge.
    std::vector<double> gap_centres_m;
    gap_centres_m.reserve(static_cast<std::size_t>(n_riser_gaps));
    const double valley_length_m = static_cast<double>(n_rows - 1) * cell_size;
    for (int g = 0; g < n_riser_gaps; ++g)
    {
        const double f = (static_cast<double>(g) + 0.5) /
                          static_cast<double>(n_riser_gaps);
        gap_centres_m.push_back(f * valley_length_m);
    }

    const int threads = resolve_thread_count(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx = flat_index(row, col, n_cols);

            const double downstream_m = static_cast<double>(row) * cell_size;
            const double f_down       = (n_rows > 1)
                ? downstream_m / valley_length_m
                : 0.0;

            const double z_channel =
                channel_elevation - valley_gradient * downstream_m;

            const double h_local =
                riser_height_m * (1.0 - riser_height_downstream_taper * f_down);

            // Cross-valley distance from the channel centre.
            const double x = std::abs(static_cast<double>(col) - centre_col)
                             * cell_size;

            // Is this row inside a riser gap?
            bool in_gap = false;
            for (const double gc : gap_centres_m)
            {
                if (std::abs(downstream_m - gc) <= 0.5 * riser_gap_width_m)
                {
                    in_gap = true;
                    break;
                }
            }

            const double flight_width_m =
                static_cast<double>(n_terraces) * tread_width_m;

            double z          = 0.0;
            int    tread_here = 0;
            bool   is_channel = false;

            if (x < 0.5 * cell_size)
            {
                // The channel itself: one cell wide at the centre column.
                z          = z_channel;
                is_channel = true;
            }
            else if (x >= flight_width_m)
            {
                // Valley wall beyond the outermost tread. Anchored to the top
                // of the flight so the wall never dips below the last tread.
                z = z_channel
                  + static_cast<double>(n_terraces) * h_local
                  + valley_wall_gradient * (x - flight_width_m);
            }
            else if (in_gap)
            {
                // A graded ramp instead of the stepped cross-section. The
                // terrace index is interpolated rather than floored, so
                // elevation climbs continuously from the channel to the top of
                // the flight at gradient h_local / tread_width_m.
                const double j_continuous = x / tread_width_m;
                z = z_channel + j_continuous * h_local;
            }
            else
            {
                const int j = std::min(
                    n_terraces,
                    static_cast<int>(std::floor(x / tread_width_m)) + 1);

                // Falls toward the channel across the tread, measured from the
                // tread's outer edge so the outer edge sits exactly j risers up.
                const double outer_edge_m =
                    static_cast<double>(j) * tread_width_m;

                z = z_channel
                  + static_cast<double>(j) * h_local
                  - tread_cross_gradient * (outer_edge_m - x);

                tread_here = j;
            }

            if (noise_amplitude_m > 0.0)
            {
                const std::uint64_t key = make_hash_key(
                    seed,
                    0x5445525241434553ULL,   // "TERRACES"
                    static_cast<std::uint64_t>(row),
                    static_cast<std::uint64_t>(col));
                z += noise_amplitude_m * symmetric_uniform(key);
            }

            dem_data[idx]   = static_cast<float>(z);
            tread_data[idx] = static_cast<float>(tread_here);
            chan_data[idx]  = is_channel ? 1.0f : 0.0f;
        }
    }

    lsdtt_terrace_staircase_surface out;
    out.dem          = make_raster(std::move(dem_data),   meta);
    out.tread_id     = make_raster(std::move(tread_data), meta);
    out.channel_mask = make_raster(std::move(chan_data),  meta);
    return out;
}

lsdtt_terrace_staircase_surface generate_terrace_staircase_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    n_terraces,
    double tread_width_m,
    double riser_height_m,
    double valley_gradient,
    double riser_height_downstream_taper,
    double tread_cross_gradient,
    int    n_riser_gaps,
    double riser_gap_width_m,
    double noise_amplitude_m,
    double channel_elevation,
    double valley_wall_gradient,
    int    seed,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);
    return generate_terrace_staircase_surface(
        meta, n_terraces, tread_width_m, riser_height_m, valley_gradient,
        riser_height_downstream_taper, tread_cross_gradient,
        n_riser_gaps, riser_gap_width_m, noise_amplitude_m,
        channel_elevation, valley_wall_gradient, seed, n_threads);
}

// =============================================================================
// generate_alluvial_fan_surface
// =============================================================================
lsdtt_alluvial_fan_surface generate_alluvial_fan_surface(
    const lsdtt_raster_meta& meta,
    double apex_x_fraction,
    double apex_y_fraction,
    double apex_elevation,
    double fan_gradient,
    double fan_concavity,
    int    n_inset_lobes,
    double lobe_inset_m,
    double lobe_half_angle_deg,
    double noise_amplitude_m,
    int    seed,
    int    n_threads)
{
    static const char* fn = "lsdtt::generate_alluvial_fan_surface";

    if (!meta.is_valid())
    {
        throw std::runtime_error(std::string(fn) + ": invalid raster metadata");
    }
    if (!(fan_gradient > 0.0))
    {
        throw std::runtime_error(std::string(fn) + ": fan_gradient must be > 0");
    }
    if (fan_concavity < 0.0 || fan_concavity >= fan_gradient)
    {
        throw std::runtime_error(
            std::string(fn) + ": fan_concavity must be in [0, fan_gradient); "
            "at or above fan_gradient the surface stops falling with distance "
            "from the apex and is no longer a fan");
    }
    if (n_inset_lobes < 0)
    {
        throw std::runtime_error(std::string(fn) + ": n_inset_lobes must be >= 0");
    }
    if (n_inset_lobes > 0 &&
        (!(lobe_half_angle_deg > 0.0) || lobe_half_angle_deg >= 90.0))
    {
        throw std::runtime_error(
            std::string(fn) + ": lobe_half_angle_deg must be in (0, 90)");
    }

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;

    if (n_rows < 3 || n_cols < 3)
    {
        throw std::runtime_error(
            std::string(fn) + ": grid must span at least 3x3 cells");
    }

    const std::size_t n_pixels =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    std::vector<float> dem_data(n_pixels, 0.0f);
    std::vector<float> lobe_data(n_pixels, 0.0f);

    // Apex in cell units. y increases southward, matching row order.
    const double apex_col = apex_x_fraction * static_cast<double>(n_cols - 1);
    const double apex_row = apex_y_fraction * static_cast<double>(n_rows - 1);

    const double width_m  = static_cast<double>(n_cols - 1) * cell_size;
    const double height_m = static_cast<double>(n_rows - 1) * cell_size;
    const double r_ref    = std::sqrt(width_m * width_m + height_m * height_m);

    // Lobe axes are spread evenly across the downslope half-plane. With the
    // apex near the south edge the fan opens northward, so the half-plane is
    // centred on -90 degrees in atan2(dy_north, dx_east) terms.
    const double half_angle_rad = lobe_half_angle_deg * M_PI / 180.0;
    std::vector<double> lobe_axis_rad;
    lobe_axis_rad.reserve(static_cast<std::size_t>(n_inset_lobes));
    for (int l = 0; l < n_inset_lobes; ++l)
    {
        const double f = (static_cast<double>(l) + 0.5) /
                          static_cast<double>(n_inset_lobes);
        // Spread across the 180 degree downslope fan, inset by the half angle
        // so no lobe is clipped by the edge of the half-plane.
        lobe_axis_rad.push_back(
            -M_PI + half_angle_rad + f * (M_PI - 2.0 * half_angle_rad));
    }

    const int threads = resolve_thread_count(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx = flat_index(row, col, n_cols);

            const double dx = (static_cast<double>(col) - apex_col) * cell_size;
            const double dy = (static_cast<double>(row) - apex_row) * cell_size;
            const double r  = std::sqrt(dx * dx + dy * dy);

            // The older fan surface: radial fall, slackening with distance when
            // fan_concavity > 0.
            double z = apex_elevation
                     - fan_gradient * r
                     + fan_concavity * r * r / (2.0 * r_ref);

            int lobe_here = 0;

            if (n_inset_lobes > 0 && r > 0.0)
            {
                const double theta = std::atan2(dy, dx);
                for (int l = 0; l < n_inset_lobes; ++l)
                {
                    // Smallest signed angle between this cell and the lobe axis.
                    double d_theta = theta - lobe_axis_rad[
                        static_cast<std::size_t>(l)];
                    while (d_theta >  M_PI) d_theta -= 2.0 * M_PI;
                    while (d_theta < -M_PI) d_theta += 2.0 * M_PI;

                    if (std::abs(d_theta) <= half_angle_rad)
                    {
                        // Parallel to the older surface, not a radial step:
                        // the residual to a fitted fan is then flat within the
                        // lobe and offset from the older surface by exactly
                        // lobe_inset_m.
                        z -= lobe_inset_m;
                        lobe_here = l + 1;
                        break;
                    }
                }
            }

            if (noise_amplitude_m > 0.0)
            {
                const std::uint64_t key = make_hash_key(
                    seed,
                    0x414C4C55564641ULL,   // "ALLUVFA"
                    static_cast<std::uint64_t>(row),
                    static_cast<std::uint64_t>(col));
                z += noise_amplitude_m * symmetric_uniform(key);
            }

            dem_data[idx]  = static_cast<float>(z);
            lobe_data[idx] = static_cast<float>(lobe_here);
        }
    }

    lsdtt_alluvial_fan_surface out;
    out.dem     = make_raster(std::move(dem_data),  meta);
    out.lobe_id = make_raster(std::move(lobe_data), meta);
    return out;
}

lsdtt_alluvial_fan_surface generate_alluvial_fan_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    double apex_x_fraction,
    double apex_y_fraction,
    double apex_elevation,
    double fan_gradient,
    double fan_concavity,
    int    n_inset_lobes,
    double lobe_inset_m,
    double lobe_half_angle_deg,
    double noise_amplitude_m,
    int    seed,
    int    n_threads)
{
    const lsdtt_raster_meta meta = build_meta(n_rows, n_cols, cell_size);
    return generate_alluvial_fan_surface(
        meta, apex_x_fraction, apex_y_fraction, apex_elevation,
        fan_gradient, fan_concavity, n_inset_lobes, lobe_inset_m,
        lobe_half_angle_deg, noise_amplitude_m, seed, n_threads);
}

} // namespace lsdtt
