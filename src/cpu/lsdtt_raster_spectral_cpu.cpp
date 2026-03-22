// lsdtt_raster_spectral.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Spectral raster analysis — FFTW-based free functions operating on
// lsdtt_raster_base objects.
//
// Current contents:
//   wavelet_quadratic_fit
//
// Notes:
//   - Gaussian kernel only.
//   - Uses the double-precision FFTW API (fftw_*).
//   - Does not use FFTW thread-planning routines.
//   - Conservative nodata handling: if the support square overlaps nodata or
//     the raster edge, output coefficients are set to nodata.
// -----------------------------------------------------------------------------

#include "lsdtt_raster_spectral_cpu.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef LSDTT_HAS_FFTW
#include <fftw3.h>
#endif

namespace lsdtt
{

#ifndef LSDTT_HAS_FFTW

poly_fit_coeffs wavelet_quadratic_fit(
    const lsdtt_raster_base&,
    double,
    int)
{
    throw std::runtime_error(
        "lsdtt::wavelet_quadratic_fit: lsdtt3 was built without FFTW support. "
        "Install FFTW development files and rebuild."
    );
}

#else

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------
static constexpr double k_pi           = 3.14159265358979323846264338327950288;
static constexpr double k_lu_rel_tol   = 1e-12;
static constexpr double k_residual_tol = 1e-8;

// -----------------------------------------------------------------------------
// checked_add / checked_mul
// -----------------------------------------------------------------------------
static std::size_t checked_add(
    std::size_t       a,
    std::size_t       b,
    const std::string& what)
{
    if (b > std::numeric_limits<std::size_t>::max() - a)
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: overflow while computing " + what
        );
    }
    return a + b;
}

static std::size_t checked_mul(
    std::size_t       a,
    std::size_t       b,
    const std::string& what)
{
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a)
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: overflow while computing " + what
        );
    }
    return a * b;
}

// -----------------------------------------------------------------------------
// flat_index
// -----------------------------------------------------------------------------
static inline std::size_t flat_index(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
         + static_cast<std::size_t>(col);
}

// -----------------------------------------------------------------------------
// FFTW RAII wrappers
// -----------------------------------------------------------------------------
class fftw_complex_buffer
{
public:
    fftw_complex_buffer(
        std::size_t       n_complex,
        const std::string& name)
    {
        const std::size_t n_bytes = checked_mul(
            n_complex,
            sizeof(fftw_complex),
            "FFTW buffer bytes for " + name
        );

        ptr_ = reinterpret_cast<fftw_complex*>(fftw_malloc(n_bytes));
        if (!ptr_)
        {
            throw std::runtime_error(
                "lsdtt::wavelet_quadratic_fit: FFTW allocation failed for " + name
            );
        }

        n_bytes_ = n_bytes;
    }

    ~fftw_complex_buffer()
    {
        if (ptr_) fftw_free(ptr_);
    }

    fftw_complex_buffer(const fftw_complex_buffer&)            = delete;
    fftw_complex_buffer& operator=(const fftw_complex_buffer&) = delete;

    fftw_complex* get()
    {
        return ptr_;
    }

    const fftw_complex* get() const
    {
        return ptr_;
    }

    std::size_t bytes() const
    {
        return n_bytes_;
    }

private:
    fftw_complex* ptr_     = nullptr;
    std::size_t   n_bytes_ = 0;
};

class fftw_plan_handle
{
public:
    fftw_plan_handle(
        int               nrows,
        int               ncols,
        fftw_complex*     in,
        fftw_complex*     out,
        int               sign,
        unsigned          flags,
        const std::string& name)
    {
        plan_ = fftw_plan_dft_2d(nrows, ncols, in, out, sign, flags);
        if (!plan_)
        {
            throw std::runtime_error(
                "lsdtt::wavelet_quadratic_fit: FFTW plan creation failed for "
                + name
            );
        }
    }

    ~fftw_plan_handle()
    {
        if (plan_) fftw_destroy_plan(plan_);
    }

    fftw_plan_handle(const fftw_plan_handle&)            = delete;
    fftw_plan_handle& operator=(const fftw_plan_handle&) = delete;

    fftw_plan get() const
    {
        return plan_;
    }

private:
    fftw_plan plan_ = nullptr;
};

// -----------------------------------------------------------------------------
// support_sample
//
// Stores one truncated-support sample used to build both the kernel and the
// Gram matrix.
//
// The kernel is stored for FFT convolution. Since DFT multiplication produces
// convolution rather than correlation, the stored kernel at offset delta must
// use the polynomial basis evaluated at -delta. This matters for the odd terms
// (x and y), and is the reason x/y below use the negated stored offset.
// -----------------------------------------------------------------------------
struct support_sample
{
    std::size_t           wrapped_flat = 0;
    double                weight       = 0.0;
    std::array<double, 6> basis{};
};

// -----------------------------------------------------------------------------
// build_support_samples
//
// Important details:
//   1. Gaussian weights are computed in pixel coordinates because the scale
//      parameter s_pixels comes from lambda / cell_size.
//   2. Polynomial coordinates are built in physical units (metres), then
//      normalised by normalisation_length for numerical stability.
//   3. Because FFT multiplication implements convolution, basis functions are
//      evaluated at the negated stored offset so that the resulting moments
//      correspond to a local weighted fit with the expected sign convention.
//
// Coordinate convention:
//   x  — positive east
//   y  — positive north
// -----------------------------------------------------------------------------
static std::vector<support_sample> build_support_samples(
    int    pad_size,
    int    padded_rows,
    int    padded_cols,
    double cell_size,
    double s_pixels,
    double normalisation_length)
{
    if (!std::isfinite(s_pixels) || s_pixels <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: invalid Gaussian scale parameter"
        );
    }

    const double two_s2 = 2.0 * s_pixels * s_pixels;

    const std::size_t kernel_width = static_cast<std::size_t>(2 * pad_size + 1);

    std::vector<support_sample> support;
    support.reserve(checked_mul(kernel_width, kernel_width, "support sample count"));

    for (int dr = -pad_size; dr <= pad_size; ++dr)
    {
        for (int dc = -pad_size; dc <= pad_size; ++dc)
        {
            // Stored offset delta = (dr, dc) in raster index coordinates
            // Physical coordinate for the source sample contributing at output
            // is -delta:
            //   x = -dc * cell_size
            //   y = +dr * cell_size
            const double x = -static_cast<double>(dc) * cell_size;
            const double y =  static_cast<double>(dr) * cell_size;

            const double u = x / normalisation_length;
            const double v = y / normalisation_length;

            const double r2_pixels =
                static_cast<double>(dr) * static_cast<double>(dr) +
                static_cast<double>(dc) * static_cast<double>(dc);

            const double w = std::exp(-r2_pixels / two_s2);

            // Wrap to top-left origin for circular convolution
            const int wr = (dr >= 0) ? dr : padded_rows + dr;
            const int wc = (dc >= 0) ? dc : padded_cols + dc;

            support_sample sample;
            sample.wrapped_flat = flat_index(wr, wc, padded_cols);
            sample.weight       = w;
            sample.basis        = { u * u, v * v, u * v, u, v, 1.0 };

            support.push_back(sample);
        }
    }

    return support;
}

// -----------------------------------------------------------------------------
// build_gram_matrix
// -----------------------------------------------------------------------------
static void build_gram_matrix(
    const std::vector<support_sample>& support,
    double                             cell_area,
    double                             gram[6][6])
{
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            gram[i][j] = 0.0;
        }
    }

    for (const auto& sample : support)
    {
        const double w = sample.weight * cell_area;

        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                gram[i][j] += sample.basis[i] * sample.basis[j] * w;
            }
        }
    }

    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            if (!std::isfinite(gram[i][j]))
            {
                throw std::runtime_error(
                    "lsdtt::wavelet_quadratic_fit: non-finite Gram matrix entry"
                );
            }
        }
    }
}

// -----------------------------------------------------------------------------
// lu6_factor
//
// In-place LU factorisation with partial pivoting.
// Uses a relative pivot tolerance scaled to the largest absolute matrix entry.
// -----------------------------------------------------------------------------
static bool lu6_factor(
    double mat[6][6],
    int    piv[6],
    double rel_tol)
{
    double max_abs = 0.0;

    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            if (!std::isfinite(mat[i][j])) return false;
            max_abs = std::max(max_abs, std::abs(mat[i][j]));
        }
    }

    if (max_abs <= 0.0) return false;

    const double pivot_tol = rel_tol * max_abs;

    for (int k = 0; k < 6; ++k)
    {
        int    max_row = k;
        double max_val = std::abs(mat[k][k]);

        for (int i = k + 1; i < 6; ++i)
        {
            const double val = std::abs(mat[i][k]);
            if (val > max_val)
            {
                max_val = val;
                max_row = i;
            }
        }

        if (!std::isfinite(max_val) || max_val <= pivot_tol)
        {
            return false;
        }

        piv[k] = max_row;

        if (max_row != k)
        {
            for (int j = 0; j < 6; ++j)
            {
                std::swap(mat[k][j], mat[max_row][j]);
            }
        }

        for (int i = k + 1; i < 6; ++i)
        {
            mat[i][k] /= mat[k][k];
            if (!std::isfinite(mat[i][k])) return false;

            for (int j = k + 1; j < 6; ++j)
            {
                mat[i][j] -= mat[i][k] * mat[k][j];
                if (!std::isfinite(mat[i][j])) return false;
            }
        }
    }

    return true;
}

// -----------------------------------------------------------------------------
// lu6_solve
//
// Solves LU*x = rhs using the factorisation produced by lu6_factor().
// Returns false if a non-finite or zero diagonal is encountered.
// -----------------------------------------------------------------------------
static bool lu6_solve(
    const double lu[6][6],
    const int    piv[6],
    const double rhs[6],
    double       x[6])
{
    double y[6];
    for (int i = 0; i < 6; ++i) y[i] = rhs[i];

    // Apply the same row swaps used during factorisation
    for (int k = 0; k < 6; ++k)
    {
        if (piv[k] != k)
        {
            std::swap(y[k], y[piv[k]]);
        }
    }

    // Forward substitution: L * z = P * rhs
    for (int i = 1; i < 6; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            y[i] -= lu[i][j] * y[j];
        }
        if (!std::isfinite(y[i])) return false;
    }

    // Back substitution: U * x = z
    for (int i = 5; i >= 0; --i)
    {
        double sum = y[i];

        for (int j = i + 1; j < 6; ++j)
        {
            sum -= lu[i][j] * x[j];
        }

        const double diag = lu[i][i];
        if (!std::isfinite(diag) || std::abs(diag) <= 0.0)
        {
            return false;
        }

        x[i] = sum / diag;
        if (!std::isfinite(x[i])) return false;
    }

    return true;
}

// -----------------------------------------------------------------------------
// relative_residual_6x6
// -----------------------------------------------------------------------------
static double relative_residual_6x6(
    const double a[6][6],
    const double x[6],
    const double b[6])
{
    double r2 = 0.0;
    double b2 = 0.0;

    for (int i = 0; i < 6; ++i)
    {
        double ri = -b[i];
        for (int j = 0; j < 6; ++j)
        {
            ri += a[i][j] * x[j];
        }

        r2 += ri * ri;
        b2 += b[i] * b[i];
    }

    const double denom = std::max(b2, 1e-30);
    return std::sqrt(r2 / denom);
}

// -----------------------------------------------------------------------------
// build_invalid_prefix
//
// Builds a summed-area table over invalid pixels:
//   invalid = mask == 0 OR value is nodata
// -----------------------------------------------------------------------------
static std::vector<std::uint64_t> build_invalid_prefix(
    const std::vector<float>&   data,
    const std::vector<uint8_t>& mask,
    const lsdtt_raster_meta&    meta)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const std::size_t prefix_rows = static_cast<std::size_t>(nrows) + 1;
    const std::size_t prefix_cols = static_cast<std::size_t>(ncols) + 1;

    std::vector<std::uint64_t> prefix(
        checked_mul(prefix_rows, prefix_cols, "invalid prefix array size"),
        0
    );

    for (int r = 0; r < nrows; ++r)
    {
        std::uint64_t row_sum = 0;

        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t idx = flat_index(r, c, ncols);
            const bool invalid = (mask[idx] == 0) || meta.is_nodata(data[idx]);

            row_sum += invalid ? 1ULL : 0ULL;

            prefix[static_cast<std::size_t>(r + 1) * prefix_cols
                 + static_cast<std::size_t>(c + 1)] =
                prefix[static_cast<std::size_t>(r) * prefix_cols
                     + static_cast<std::size_t>(c + 1)] + row_sum;
        }
    }

    return prefix;
}

// -----------------------------------------------------------------------------
// prefix_sum_rect
//
// Rectangle sum from summed-area table using half-open bounds:
//   rows [r0, r1)
//   cols [c0, c1)
// -----------------------------------------------------------------------------
static std::uint64_t prefix_sum_rect(
    const std::vector<std::uint64_t>& prefix,
    int                               ncols,
    int                               r0,
    int                               c0,
    int                               r1,
    int                               c1)
{
    const std::size_t stride = static_cast<std::size_t>(ncols) + 1;

    const std::size_t i00 = static_cast<std::size_t>(r0) * stride + static_cast<std::size_t>(c0);
    const std::size_t i01 = static_cast<std::size_t>(r0) * stride + static_cast<std::size_t>(c1);
    const std::size_t i10 = static_cast<std::size_t>(r1) * stride + static_cast<std::size_t>(c0);
    const std::size_t i11 = static_cast<std::size_t>(r1) * stride + static_cast<std::size_t>(c1);

    return prefix[i11] - prefix[i01] - prefix[i10] + prefix[i00];
}

// -----------------------------------------------------------------------------
// build_support_valid_mask
//
// Conservative support validity:
//   valid only if the full square support [r-pad, r+pad] x [c-pad, c+pad]
//   lies within the raster AND contains no invalid pixels.
// -----------------------------------------------------------------------------
static std::vector<uint8_t> build_support_valid_mask(
    const std::vector<std::uint64_t>& invalid_prefix,
    const lsdtt_raster_meta&          meta,
    int                               pad_size)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    std::vector<uint8_t> support_valid(
        checked_mul(static_cast<std::size_t>(nrows),
                    static_cast<std::size_t>(ncols),
                    "support-valid mask size"),
        uint8_t{0}
    );

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (r < pad_size || r >= nrows - pad_size ||
                c < pad_size || c >= ncols - pad_size)
            {
                continue;
            }

            const std::uint64_t n_invalid = prefix_sum_rect(
                invalid_prefix,
                ncols,
                r - pad_size,
                c - pad_size,
                r + pad_size + 1,
                c + pad_size + 1
            );

            if (n_invalid == 0)
            {
                support_valid[flat_index(r, c, ncols)] = uint8_t{1};
            }
        }
    }

    return support_valid;
}

// -----------------------------------------------------------------------------
// resolve_thread_count
//
// Used only for OpenMP loops in this file.
// FFTW itself is used in serial mode here.
// -----------------------------------------------------------------------------
static int resolve_thread_count(int n_threads)
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
// wavelet_quadratic_fit
// -----------------------------------------------------------------------------
poly_fit_coeffs wavelet_quadratic_fit(
    const lsdtt_raster_cpu&raster,
    double                   lambda,
    int                      n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: invalid raster metadata"
        );
    }

    meta.check_projected_metres();

    if (!std::isfinite(lambda) || lambda <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: lambda must be finite and > 0"
        );
    }

    const int    nrows     = meta.n_rows;
    const int    ncols     = meta.n_cols;
    const double cell_size = meta.cell_size;
    const float  nodata    = static_cast<float>(meta.no_data_value);

    if (lambda < cell_size)
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: lambda must be >= cell_size"
        );
    }

    const std::size_t n_pixels = checked_mul(
        static_cast<std::size_t>(nrows),
        static_cast<std::size_t>(ncols),
        "raster pixel count"
    );

    const std::vector<float>&   data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    if (data.size() != n_pixels || mask.size() != n_pixels)
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: data / mask size does not match metadata"
        );
    }

    // Torrence and Compo mapping from wavelength to Gaussian scale, in pixels
    const double s_pixels =
        ((lambda / cell_size) * std::sqrt(5.0 / 2.0)) / (2.0 * k_pi);

    if (!std::isfinite(s_pixels) || s_pixels <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: invalid Gaussian scale parameter"
        );
    }

    const int pad_size = std::max(1, static_cast<int>(std::ceil(8.0 * s_pixels)));

    if (2 * pad_size >= nrows || 2 * pad_size >= ncols)
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: support radius is too large for the raster "
            "dimensions. Reduce lambda or use a larger DEM."
        );
    }

    const std::size_t halo_size = checked_mul(
        static_cast<std::size_t>(2),
        static_cast<std::size_t>(pad_size),
        "padded raster halo"
    );

    const std::size_t padded_rows_sz = checked_add(
        static_cast<std::size_t>(nrows),
        halo_size,
        "padded row count"
    );

    const std::size_t padded_cols_sz = checked_add(
        static_cast<std::size_t>(ncols),
        halo_size,
        "padded column count"
    );

    if (padded_rows_sz > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        padded_cols_sz > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: padded FFT dimensions exceed FFTW int limits"
        );
    }

    const int padded_rows = static_cast<int>(padded_rows_sz);
    const int padded_cols = static_cast<int>(padded_cols_sz);

    const std::size_t n_complex = checked_mul(
        padded_rows_sz,
        padded_cols_sz,
        "padded FFT element count"
    );

    const double normalisation_length = std::max(lambda, cell_size);

    const std::vector<support_sample> support = build_support_samples(
        pad_size,
        padded_rows,
        padded_cols,
        cell_size,
        s_pixels,
        normalisation_length
    );

    double gram[6][6];
    build_gram_matrix(support, cell_size * cell_size, gram);

    double lu[6][6];
    std::memcpy(lu, gram, sizeof(gram));

    int piv[6];
    if (!lu6_factor(lu, piv, k_lu_rel_tol))
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: Gram matrix factorisation failed. "
            "The system is singular or badly conditioned."
        );
    }

    const std::vector<std::uint64_t> invalid_prefix =
        build_invalid_prefix(data, mask, meta);

    const std::vector<uint8_t> support_valid =
        build_support_valid_mask(invalid_prefix, meta, pad_size);

    const std::size_t n_valid_support = static_cast<std::size_t>(
        std::count_if(
            support_valid.begin(),
            support_valid.end(),
            [](uint8_t v) { return v != 0; }
        )
    );

    if (n_valid_support == 0)
    {
        throw std::runtime_error(
            "lsdtt::wavelet_quadratic_fit: no valid pixels remain after "
            "conservative support masking. Reduce lambda or remove/fill nodata."
        );
    }

    const int threads = resolve_thread_count(n_threads);

    fftw_complex_buffer dem_in   (n_complex, "dem_in");
    fftw_complex_buffer dem_freq (n_complex, "dem_freq");
    fftw_complex_buffer ker_in   (n_complex, "ker_in");
    fftw_complex_buffer ker_freq (n_complex, "ker_freq");
    fftw_complex_buffer prod_freq(n_complex, "prod_freq");
    fftw_complex_buffer out_space(n_complex, "out_space");

    fftw_plan_handle plan_dem(
        padded_rows, padded_cols,
        dem_in.get(), dem_freq.get(),
        FFTW_FORWARD, FFTW_MEASURE,
        "DEM forward transform"
    );

    fftw_plan_handle plan_ker(
        padded_rows, padded_cols,
        ker_in.get(), ker_freq.get(),
        FFTW_FORWARD, FFTW_MEASURE,
        "kernel forward transform"
    );

    fftw_plan_handle plan_inv(
        padded_rows, padded_cols,
        prod_freq.get(), out_space.get(),
        FFTW_BACKWARD, FFTW_MEASURE,
        "inverse transform"
    );

    std::memset(dem_in.get(), 0, dem_in.bytes());

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t src_idx = flat_index(r, c, ncols);

            if (mask[src_idx] == 0) continue;
            if (meta.is_nodata(data[src_idx])) continue;

            const std::size_t dst_idx =
                flat_index(r + pad_size, c + pad_size, padded_cols);

            dem_in.get()[dst_idx][0] = static_cast<double>(data[src_idx]);
            dem_in.get()[dst_idx][1] = 0.0;
        }
    }

    fftw_execute(plan_dem.get());

    const double unset = std::numeric_limits<double>::quiet_NaN();

    std::array<std::vector<double>, 6> moments;
    for (auto& m : moments)
    {
        m.assign(n_pixels, unset);
    }

    const double conv_scale =
        (cell_size * cell_size) / static_cast<double>(n_complex);

    for (int k = 0; k < 6; ++k)
    {
        std::memset(ker_in.get(), 0, ker_in.bytes());

        for (const auto& sample : support)
        {
            ker_in.get()[sample.wrapped_flat][0] =
                sample.basis[k] * sample.weight;
            ker_in.get()[sample.wrapped_flat][1] = 0.0;
        }

        fftw_execute(plan_ker.get());

        for (std::size_t i = 0; i < n_complex; ++i)
        {
            const double ar = dem_freq.get()[i][0];
            const double ai = dem_freq.get()[i][1];
            const double br = ker_freq.get()[i][0];
            const double bi = ker_freq.get()[i][1];

            prod_freq.get()[i][0] = ar * br - ai * bi;
            prod_freq.get()[i][1] = ar * bi + ai * br;
        }

        fftw_execute(plan_inv.get());

        for (int r = 0; r < nrows; ++r)
        {
            for (int c = 0; c < ncols; ++c)
            {
                const std::size_t out_idx = flat_index(r, c, ncols);

                if (support_valid[out_idx] == 0) continue;

                const std::size_t padded_idx =
                    flat_index(r + pad_size, c + pad_size, padded_cols);

                const double m = out_space.get()[padded_idx][0] * conv_scale;

                if (std::isfinite(m))
                {
                    moments[k][out_idx] = m;
                }
            }
        }
    }

    std::array<std::vector<float>, 6>   coeff_data;
    std::array<std::vector<uint8_t>, 6> coeff_mask;

    for (auto& v : coeff_data)
    {
        v.assign(n_pixels, nodata);
    }

    for (auto& m : coeff_mask)
    {
        m.assign(n_pixels, uint8_t{0});
    }

    const double inv_l  = 1.0 / normalisation_length;
    const double inv_l2 = inv_l * inv_l;

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t idx = flat_index(r, c, ncols);

            if (support_valid[idx] == 0) continue;

            double rhs[6];
            bool   valid = true;

            for (int k = 0; k < 6; ++k)
            {
                const double m = moments[k][idx];
                if (!std::isfinite(m))
                {
                    valid = false;
                    break;
                }
                rhs[k] = m;
            }

            if (!valid) continue;

            double x[6];
            if (!lu6_solve(lu, piv, rhs, x)) continue;

            const double rel_resid = relative_residual_6x6(gram, x, rhs);
            if (!std::isfinite(rel_resid) || rel_resid > k_residual_tol) continue;

            // Convert from scaled local coordinates (u, v) back to physical
            // coordinates (x, y) in metres.
            //
            // If z = a_u u^2 + b_u v^2 + c_u u v + d_u u + e_u v + f_u
            // and u = x/L, v = y/L, then:
            //   a = a_u / L^2
            //   b = b_u / L^2
            //   c = c_u / L^2
            //   d = d_u / L
            //   e = e_u / L
            //   f = f_u
            const double phys[6] =
            {
                x[0] * inv_l2,
                x[1] * inv_l2,
                x[2] * inv_l2,
                x[3] * inv_l,
                x[4] * inv_l,
                x[5]
            };

            for (int k = 0; k < 6; ++k)
            {
                if (!std::isfinite(phys[k]))
                {
                    valid = false;
                    break;
                }
            }

            if (!valid) continue;

            for (int k = 0; k < 6; ++k)
            {
                coeff_data[k][idx] = static_cast<float>(phys[k]);
                coeff_mask[k][idx] = uint8_t{1};
            }
        }
    }

    auto make_raster = [&](int k) -> lsdtt_raster_cpu
    {
        return lsdtt_raster_cpu(
            std::move(coeff_data[k]),
            std::move(coeff_mask[k]),
            meta
        );
    };

    return poly_fit_coeffs
    {
        make_raster(0),
        make_raster(1),
        make_raster(2),
        make_raster(3),
        make_raster(4),
        make_raster(5)
    };
}

#endif  // LSDTT_HAS_FFTW

} // namespace lsdtt
