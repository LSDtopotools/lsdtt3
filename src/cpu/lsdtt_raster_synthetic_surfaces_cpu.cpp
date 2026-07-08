// lsdtt_raster_synthetic_surfaces.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"

#include "../utils/lsdtt_thread_utils.hpp"

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
    meta.data_type     = GDT_Float32;

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
        // 3. Apply spectral amplitude scaling ∝ f^(-beta/2)
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
// validate_range_noise_params  — shared parameter checks
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
// which gives depth = decay_length at the centre, tapering to 0 at u = ±1.
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

} // namespace lsdtt
