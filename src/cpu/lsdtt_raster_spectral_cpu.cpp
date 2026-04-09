// lsdtt_raster_spectral.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Spectral raster analysis — FFTW-based free functions operating on
// lsdtt_raster_base objects.
//
// Current contents:
//   wavelet_quadratic_fit
//   various specral filters
//
// Notes:
//   - Gaussian kernel only.
//   - Uses the double-precision FFTW API (fftw_*).
//   - Does not use FFTW thread-planning routines.
//   - Conservative nodata handling: if the support square overlaps nodata or
//     the raster edge, output coefficients are set to nodata.
// -----------------------------------------------------------------------------

#include "lsdtt_raster_spectral_cpu.hpp"
#include "lsdtt_raster_nodata_filling_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"

#include "../utils/lsdtt_histogram.hpp"
#include "../utils/lsdtt_regression.hpp"
#include "../utils/lsdtt_thread_utils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
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
    const lsdtt_raster_cpu&,
    double,
    int)
{
    throw std::runtime_error(
        "lsdtt::wavelet_quadratic_fit: lsdtt3 was built without FFTW support. "
        "Install FFTW development files and rebuild."
    );
}

static const char* k_no_fftw_msg =
    "lsdtt3 was built without FFTW support. "
    "Install FFTW development files and rebuild.";

lsdtt_raster_cpu fftw2D_lowpass_filter(
    const lsdtt_raster_cpu&, double, spectral_window_type, bool, int)
{
    throw std::runtime_error(std::string("lsdtt::fftw2D_lowpass_filter: ") + k_no_fftw_msg);
}

lsdtt_raster_cpu fftw2D_highpass_filter(
    const lsdtt_raster_cpu&, double, spectral_window_type, bool, int)
{
    throw std::runtime_error(std::string("lsdtt::fftw2D_highpass_filter: ") + k_no_fftw_msg);
}

lsdtt_raster_cpu fftw2D_bandpass_filter(
    const lsdtt_raster_cpu&, double, double, spectral_window_type, bool, int)
{
    throw std::runtime_error(std::string("lsdtt::fftw2D_bandpass_filter: ") + k_no_fftw_msg);
}

lsdtt_raster_cpu fftw2D_wiener_filter(
    const lsdtt_raster_cpu&, double, double, spectral_window_type, bool, int)
{
    throw std::runtime_error(std::string("lsdtt::fftw2D_wiener_filter: ") + k_no_fftw_msg);
}

lsdtt_spectral_analysis_result full_spectral_analysis(
    const lsdtt_raster_cpu&,
    const lsdtt_spectral_analysis_options&)
{
    throw std::runtime_error(
        "lsdtt::full_spectral_analysis: lsdtt3 was built without FFTW support. "
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
// fftw_double_buffer — RAII wrapper for a real double array allocated via FFTW
// -----------------------------------------------------------------------------
class fftw_double_buffer
{
public:
    fftw_double_buffer(std::size_t n, const std::string& name)
        : n_(n)
    {
        ptr_ = reinterpret_cast<double*>(fftw_malloc(sizeof(double) * n));
        if (!ptr_)
        {
            throw std::runtime_error(
                "lsdtt spectral filter: FFTW double buffer allocation failed for " + name
            );
        }
    }

    ~fftw_double_buffer() { if (ptr_) fftw_free(ptr_); }

    fftw_double_buffer(const fftw_double_buffer&)            = delete;
    fftw_double_buffer& operator=(const fftw_double_buffer&) = delete;

    double*     get()  { return ptr_; }
    std::size_t size() const { return n_; }

private:
    double*     ptr_ = nullptr;
    std::size_t n_   = 0;
};

// -----------------------------------------------------------------------------
// fftw_halfcomplex_buffer
//
// Holds the r2c output: Ly rows × (Lx/2+1) complex columns.
// The non-redundant half of the spectrum for a real input of size Ly × Lx.
// -----------------------------------------------------------------------------
class fftw_halfcomplex_buffer
{
public:
    fftw_halfcomplex_buffer(std::size_t n, const std::string& name)
        : n_(n)
    {
        ptr_ = reinterpret_cast<fftw_complex*>(
            fftw_malloc(sizeof(fftw_complex) * n)
        );
        if (!ptr_)
        {
            throw std::runtime_error(
                "lsdtt spectral filter: FFTW halfcomplex buffer allocation failed for " + name
            );
        }
    }

    ~fftw_halfcomplex_buffer() { if (ptr_) fftw_free(ptr_); }

    fftw_halfcomplex_buffer(const fftw_halfcomplex_buffer&)            = delete;
    fftw_halfcomplex_buffer& operator=(const fftw_halfcomplex_buffer&) = delete;

    fftw_complex* get()  { return ptr_; }
    std::size_t   size() const { return n_; }

private:
    fftw_complex* ptr_ = nullptr;
    std::size_t   n_   = 0;
};

// -----------------------------------------------------------------------------
// fftw_r2c_plan_handle / fftw_c2r_plan_handle
//
// RAII plan handles for real-to-complex and complex-to-real 2D transforms.
// Always uses FFTW_ESTIMATE so planning never touches the data arrays.
// The c2r transform may destroy its complex input on execution — this is
// expected and accounted for in spectral_filter_pipeline.
// -----------------------------------------------------------------------------
class fftw_r2c_plan_handle
{
public:
    fftw_r2c_plan_handle(
        int            nrows,
        int            ncols,
        double*        in,
        fftw_complex*  out,
        const std::string& name)
    {
        plan_ = fftw_plan_dft_r2c_2d(nrows, ncols, in, out, FFTW_ESTIMATE);
        if (!plan_)
        {
            throw std::runtime_error(
                "lsdtt spectral filter: r2c plan creation failed for " + name
            );
        }
    }

    ~fftw_r2c_plan_handle() { if (plan_) fftw_destroy_plan(plan_); }

    fftw_r2c_plan_handle(const fftw_r2c_plan_handle&)            = delete;
    fftw_r2c_plan_handle& operator=(const fftw_r2c_plan_handle&) = delete;

    fftw_plan get() const { return plan_; }

private:
    fftw_plan plan_ = nullptr;
};

class fftw_c2r_plan_handle
{
public:
    fftw_c2r_plan_handle(
        int            nrows,
        int            ncols,
        fftw_complex*  in,
        double*        out,
        const std::string& name)
    {
        plan_ = fftw_plan_dft_c2r_2d(nrows, ncols, in, out, FFTW_ESTIMATE);
        if (!plan_)
        {
            throw std::runtime_error(
                "lsdtt spectral filter: c2r plan creation failed for " + name
            );
        }
    }

    ~fftw_c2r_plan_handle() { if (plan_) fftw_destroy_plan(plan_); }

    fftw_c2r_plan_handle(const fftw_c2r_plan_handle&)            = delete;
    fftw_c2r_plan_handle& operator=(const fftw_c2r_plan_handle&) = delete;

    fftw_plan get() const { return plan_; }

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

// =============================================================================
// Spectral filter helpers
// =============================================================================

// Next power of 2 >= n
static int next_pow2(int n)
{
    if (n <= 1) return 1;
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

// -----------------------------------------------------------------------------
// fit_detrend_plane
//
// Fits z = a*col + b*row + c to all valid pixels using 3×3 normal equations
// with partial pivoting.  Row/col are used as x/y (pixel index coordinates)
// since the trend is only used as a DC/tilt removal step before FFT; physical
// coordinates are not needed.
//
// Returns true on success.  On failure (< 3 valid pixels or singular system)
// returns false and sets all coefficients to 0 (no detrend applied).
// -----------------------------------------------------------------------------
static bool fit_detrend_plane(
    const std::vector<float>&   data,
    const std::vector<uint8_t>& mask,
    const lsdtt_raster_meta&    meta,
    double& a_out,
    double& b_out,
    double& c_out)
{
    a_out = b_out = c_out = 0.0;

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    double Sxx = 0, Sxy = 0, Sx = 0;
    double             Syy = 0, Sy = 0;
    double                         Sn = 0;
    double Sxz = 0, Syz = 0, Sz = 0;

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t i =
                static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) + c;
            if (mask[i] == 0 || meta.is_nodata(data[i])) continue;

            const double x = static_cast<double>(c);
            const double y = static_cast<double>(r);
            const double z = static_cast<double>(data[i]);

            Sxx += x*x;  Sxy += x*y;  Sx += x;
            Syy += y*y;  Sy  += y;    Sn += 1.0;
            Sxz += x*z;  Syz += y*z;  Sz += z;
        }
    }

    if (Sn < 3.0) return false;

    // Augmented matrix [A | b] for Gaussian elimination
    double A[3][4] =
    {
        { Sxx, Sxy, Sx, Sxz },
        { Sxy, Syy, Sy, Syz },
        { Sx,  Sy,  Sn, Sz  }
    };

    for (int col = 0; col < 3; ++col)
    {
        int    prow = col;
        double pval = std::abs(A[col][col]);
        for (int row = col + 1; row < 3; ++row)
        {
            if (std::abs(A[row][col]) > pval)
            {
                pval = std::abs(A[row][col]);
                prow = row;
            }
        }
        if (pval < 1e-14) return false;

        if (prow != col)
            for (int j = 0; j <= 3; ++j) std::swap(A[col][j], A[prow][j]);

        for (int row = col + 1; row < 3; ++row)
        {
            const double f = A[row][col] / A[col][col];
            for (int j = col; j <= 3; ++j) A[row][j] -= f * A[col][j];
        }
    }

    double x[3];
    for (int i = 2; i >= 0; --i)
    {
        x[i] = A[i][3];
        for (int j = i + 1; j < 3; ++j) x[i] -= A[i][j] * x[j];
        if (std::abs(A[i][i]) < 1e-30) return false;
        x[i] /= A[i][i];
        if (!std::isfinite(x[i])) return false;
    }

    a_out = x[0];
    b_out = x[1];
    c_out = x[2];
    return true;
}

// -----------------------------------------------------------------------------
// build_2d_window
//
// Returns a flat (nrows × ncols, row-major) array of window weights in [0, 1].
//
// For hann:    w(r) = 0.5 * (1 + cos(pi * r / r'))
// For hamming: w(r) = 0.54 + 0.46 * cos(pi * r / r')
//
// where r is the Euclidean distance from the raster centre and r' is the
// distance to the edge of the inscribed ellipse at the same angle.
// Both are 1.0 at the centre, taper to near-zero at the ellipse boundary,
// and are 0 outside.  The Hamming window retains a small non-zero base level
// (~0.08) that further reduces spectral leakage.
// -----------------------------------------------------------------------------
static std::vector<double> build_2d_window(
    int                  nrows,
    int                  ncols,
    spectral_window_type win_type)
{
    const std::size_t npix = static_cast<std::size_t>(nrows) * ncols;
    std::vector<double> w(npix, 1.0);

    if (win_type == spectral_window_type::none) return w;

    const double a = (ncols - 1) * 0.5;   // x semi-axis (cols)
    const double b = (nrows - 1) * 0.5;   // y semi-axis (rows)

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const double x    = c - a;
            const double y    = r - b;
            const double dist = std::sqrt(x*x + y*y);

            double r_prime;
            if (dist < 1e-12)
            {
                r_prime = std::min(a, b);
            }
            else
            {
                const double cos_t = x / dist;
                const double sin_t = y / dist;
                r_prime = (a * b) /
                          std::sqrt(b*b*cos_t*cos_t + a*a*sin_t*sin_t);
            }

            double coeff = 0.0;
            if (dist < r_prime)
            {
                const double ratio = k_pi * dist / r_prime;
                if (win_type == spectral_window_type::hann)
                    coeff = 0.5 * (1.0 + std::cos(ratio));
                else   // hamming — corrected formula, equals 1 at centre
                    coeff = 0.54 + 0.46 * std::cos(ratio);
            }

            w[static_cast<std::size_t>(r) * ncols + c] = coeff;
        }
    }

    return w;
}

// -----------------------------------------------------------------------------
// radial_freq_r2c
//
// Radial frequency (cycles/metre) of element (i, j) in the UNSHIFTED r2c
// output of a transform with padded dimensions Ly × Lx:
//
//   f_x = j / (Lx * cell_size)                        j in [0, Lx/2]
//   f_y = i / (Ly * cell_size)      for i in [0, Ly/2]
//   f_y = (i - Ly) / (Ly * cell_size)  for i in [Ly/2+1, Ly-1]
// -----------------------------------------------------------------------------
static inline double radial_freq_r2c(
    int i, int j, int Ly, double dfx, double dfy)
{
    const double fy = (i <= Ly / 2) ? i * dfy : (i - Ly) * dfy;
    const double fx = j * dfx;
    return std::sqrt(fx*fx + fy*fy);
}

// -----------------------------------------------------------------------------
// spectral_filter_pipeline
//
// Shared forward / inverse FFT pipeline.  The caller provides apply_kernel,
// which receives a pointer to the r2c spectrum (Ly × Lx_half complex elements,
// row-major, unshifted) and modifies it in place.
//
// Parameters passed to apply_kernel:
//   spec     — pointer to Ly × Lx_half fftw_complex array
//   Ly, Lx   — padded dimensions
//   Lx_half  — Lx/2 + 1
//   dfx, dfy — frequency increments in cycles/metre
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu spectral_filter_pipeline(
    const lsdtt_raster_cpu& dem,
    spectral_window_type     win_type,
    const std::function<void(fftw_complex*, int, int, int, double, double)>&
        apply_kernel)
{
    const lsdtt_raster_meta& meta      = dem.get_meta();
    const int                NRows     = meta.n_rows;
    const int                NCols     = meta.n_cols;
    const double             cell_size = meta.cell_size;
    const float              nodata_f  = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   data  = dem.get_data();
    const std::vector<uint8_t>& valid = dem.get_mask();
    const std::size_t npix = static_cast<std::size_t>(NRows) * NCols;

    // -------------------------------------------------------------------------
    // 1. Detrend: fit plane z = a*col + b*row + c, subtract from valid pixels.
    //    Set nodata and invalid pixels to 0 (they become zero-padding).
    // -------------------------------------------------------------------------
    double a_tr = 0.0, b_tr = 0.0, c_tr = 0.0;
    fit_detrend_plane(data, valid, meta, a_tr, b_tr, c_tr);

    std::vector<double> detrended(npix, 0.0);
    std::vector<double> trend_vals(npix, 0.0);

    for (int r = 0; r < NRows; ++r)
    {
        for (int c = 0; c < NCols; ++c)
        {
            const std::size_t i =
                static_cast<std::size_t>(r) * NCols + c;
            const double t = a_tr * c + b_tr * r + c_tr;
            trend_vals[i] = t;

            if (valid[i] != 0 && !meta.is_nodata(data[i]))
                detrended[i] = static_cast<double>(data[i]) - t;
        }
    }

    // -------------------------------------------------------------------------
    // 2. Window: multiply detrended values by 2D taper
    // -------------------------------------------------------------------------
    const std::vector<double> win_w = build_2d_window(NRows, NCols, win_type);

    // -------------------------------------------------------------------------
    // 3. Pad to next power-of-2 dimensions
    // -------------------------------------------------------------------------
    const int Ly      = next_pow2(NRows);
    const int Lx      = next_pow2(NCols);
    const int Lx_half = Lx / 2 + 1;

    const double dfx = 1.0 / (cell_size * Lx);
    const double dfy = 1.0 / (cell_size * Ly);

    const std::size_t n_real    = static_cast<std::size_t>(Ly)  * Lx;
    const std::size_t n_complex = static_cast<std::size_t>(Ly)  * Lx_half;

    fftw_double_buffer      real_buf (n_real,    "real");
    fftw_halfcomplex_buffer cplx_buf (n_complex, "spectrum");

    // Create both plans before writing any data (FFTW_ESTIMATE never touches
    // the arrays during planning; safe to plan before filling).
    fftw_r2c_plan_handle r2c_plan(Ly, Lx, real_buf.get(), cplx_buf.get(), "r2c");
    fftw_c2r_plan_handle c2r_plan(Ly, Lx, cplx_buf.get(), real_buf.get(), "c2r");

    // Fill real buffer: windowed detrended signal, zero-padded
    std::fill(real_buf.get(), real_buf.get() + n_real, 0.0);
    for (int r = 0; r < NRows; ++r)
    {
        for (int c = 0; c < NCols; ++c)
        {
            const std::size_t i_dem = static_cast<std::size_t>(r) * NCols + c;
            real_buf.get()[static_cast<std::size_t>(r) * Lx + c] =
                detrended[i_dem] * win_w[i_dem];
        }
    }

    // -------------------------------------------------------------------------
    // 4. Forward FFT (r2c)
    // -------------------------------------------------------------------------
    fftw_execute(r2c_plan.get());

    // -------------------------------------------------------------------------
    // 5. Apply filter kernel (caller-provided)
    // -------------------------------------------------------------------------
    apply_kernel(cplx_buf.get(), Ly, Lx, Lx_half, dfx, dfy);

    // -------------------------------------------------------------------------
    // 6. Inverse FFT (c2r), normalise, restore trend, re-apply nodata mask
    //    c2r may destroy cplx_buf — that is expected.
    // -------------------------------------------------------------------------
    fftw_execute(c2r_plan.get());

    const double norm = 1.0 / static_cast<double>(n_real);  // 1 / (Ly * Lx)

    std::vector<float>   out_data(npix, nodata_f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (int r = 0; r < NRows; ++r)
    {
        for (int c = 0; c < NCols; ++c)
        {
            const std::size_t i_dem = static_cast<std::size_t>(r) * NCols + c;
            if (valid[i_dem] == 0 || meta.is_nodata(data[i_dem])) continue;

            const double z = real_buf.get()[static_cast<std::size_t>(r) * Lx + c]
                           * norm
                           + trend_vals[i_dem];

            out_data[i_dem] = static_cast<float>(z);
            out_mask[i_dem] = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}



namespace
{




// -----------------------------------------------------------------------------
// spectral_psd_product
//
// Internal product from a full DEM periodogram calculation.
// Stores the 2D PSD in a flat row-major vector and the derived radial PSD.
// -----------------------------------------------------------------------------
struct spectral_psd_product
{
    int padded_rows = 0;
    int padded_cols = 0;

    double dfx = 0.0;
    double dfy = 0.0;
    double nyquist_frequency = 0.0;

    double window_sum_squared = 0.0;
    double detrended_variance = 0.0;

    std::vector<double> psd_2d_flat;
    std::vector<double> radial_frequency;
    std::vector<double> radial_psd;
};

// -----------------------------------------------------------------------------
// background_spectrum_product
//
// Internal averaged background spectrum built from many synthetic DEMs.
// -----------------------------------------------------------------------------
struct background_spectrum_product
{
    int padded_rows = 0;
    int padded_cols = 0;

    double dfx = 0.0;
    double dfy = 0.0;
    double nyquist_frequency = 0.0;

    std::vector<double> background_psd_2d_flat;
    std::vector<double> background_radial_frequency;
    std::vector<double> background_radial_psd;
};






// -----------------------------------------------------------------------------
// compute_window_sum_squared
//
// Returns sum(w^2) for the chosen window weights.
// -----------------------------------------------------------------------------
static double compute_window_sum_squared(
    const std::vector<double>& window_weights)
{
    double wss = 0.0;
    for (double w : window_weights)
    {
        wss += w * w;
    }
    return wss;
}

// -----------------------------------------------------------------------------
// compute_valid_mean_and_variance
//
// Computes mean and population variance over valid raster pixels only.
// Valid means:
//   mask != 0
//   value is not nodata
//
// Throws if there are no valid pixels.
// -----------------------------------------------------------------------------
static std::pair<double, double> compute_valid_mean_and_variance(
    const std::vector<float>&   data,
    const std::vector<uint8_t>& mask,
    const lsdtt_raster_meta&    meta)
{
    double sum = 0.0;
    double sumsq = 0.0;
    std::size_t count = 0;

    const std::size_t n = data.size();

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(data[i])) continue;

        const double z = static_cast<double>(data[i]);
        if (!std::isfinite(z)) continue;

        sum += z;
        sumsq += z * z;
        ++count;
    }

    if (count == 0)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: no valid pixels available for mean/variance"
        );
    }

    const double mean = sum / static_cast<double>(count);
    const double var =
        std::max(0.0, sumsq / static_cast<double>(count) - mean * mean);

    return std::make_pair(mean, var);
}


// -----------------------------------------------------------------------------
// sort_xy_pairs_by_x
//
// Sorts paired vectors (x, y) by ascending x.
// Throws if the input sizes differ.
// -----------------------------------------------------------------------------
static void sort_xy_pairs_by_x(
    std::vector<double>& x,
    std::vector<double>& y)
{
    if (x.size() != y.size())
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: sort_xy_pairs_by_x size mismatch"
        );
    }

    std::vector<std::size_t> index(x.size(), 0);
    for (std::size_t i = 0; i < index.size(); ++i)
    {
        index[i] = i;
    }

    std::stable_sort(
        index.begin(),
        index.end(),
        [&x](std::size_t a, std::size_t b)
        {
            return x[a] < x[b];
        }
    );

    std::vector<double> x_sorted(x.size(), 0.0);
    std::vector<double> y_sorted(y.size(), 0.0);

    for (std::size_t i = 0; i < index.size(); ++i)
    {
        x_sorted[i] = x[index[i]];
        y_sorted[i] = y[index[i]];
    }

    x.swap(x_sorted);
    y.swap(y_sorted);
}

// -----------------------------------------------------------------------------
// collapse_shifted_psd_to_radial
//
// Converts a shifted 2D PSD into a radial PSD.
//
// Assumptions:
//   - psd_2d_flat is row-major with size padded_rows * padded_cols
//   - the PSD has already been fftshift-ed so that zero frequency is at:
//         row = padded_rows / 2
//         col = padded_cols / 2
//
// Legacy-parity note:
//   This follows the old LSDRasterSpectral logic and uses only the left half
//   of the shifted spectrum (cols 0..Lx/2 inclusive), multiplying every kept
//   PSD value by 2.0 to account for Hermitian symmetry of the FFT of a real
//   input field.
// -----------------------------------------------------------------------------
static std::pair<std::vector<double>, std::vector<double>>
collapse_shifted_psd_to_radial(
    const std::vector<double>& psd_2d_flat,
    int                        padded_rows,
    int                        padded_cols,
    double                     dfx,
    double                     dfy,
    double                     nyquist_frequency)
{
    if (padded_rows <= 0 || padded_cols <= 0)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: invalid padded dimensions in "
            "collapse_shifted_psd_to_radial"
        );
    }

    const std::size_t expected_size =
        static_cast<std::size_t>(padded_rows) *
        static_cast<std::size_t>(padded_cols);

    if (psd_2d_flat.size() != expected_size)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: 2D PSD size mismatch in "
            "collapse_shifted_psd_to_radial"
        );
    }

    if (!std::isfinite(dfx) || dfx <= 0.0 ||
        !std::isfinite(dfy) || dfy <= 0.0 ||
        !std::isfinite(nyquist_frequency) || nyquist_frequency <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: invalid frequency metadata in "
            "collapse_shifted_psd_to_radial"
        );
    }

    const int centre_row = padded_rows / 2;
    const int centre_col = padded_cols / 2;
    const int half_cols  = padded_cols / 2;

    std::vector<double> radial_frequency;
    std::vector<double> radial_psd;

    radial_frequency.reserve(
        static_cast<std::size_t>(padded_rows) *
        static_cast<std::size_t>(half_cols + 1)
    );
    radial_psd.reserve(
        static_cast<std::size_t>(padded_rows) *
        static_cast<std::size_t>(half_cols + 1)
    );

    for (int row = 0; row < padded_rows; ++row)
    {
        for (int col = 0; col <= half_cols; ++col)
        {
            const double y = static_cast<double>(row - centre_row);
            const double x = static_cast<double>(col - centre_col);

            const double radial_freq =
                std::sqrt(y * y * dfy * dfy + x * x * dfx * dfx);

            if (radial_freq > nyquist_frequency)
            {
                continue;
            }

            const std::size_t idx = flat_index(row, col, padded_cols);
            const double p = psd_2d_flat[idx];

            if (!std::isfinite(p))
            {
                continue;
            }

            radial_frequency.push_back(radial_freq);
            radial_psd.push_back(2.0 * p);
        }
    }

    sort_xy_pairs_by_x(radial_frequency, radial_psd);

    return std::make_pair(
        std::move(radial_frequency),
        std::move(radial_psd)
    );
}



// -----------------------------------------------------------------------------
// Forward declarations for helpers we will add next
// -----------------------------------------------------------------------------
static spectral_psd_product compute_periodogram(
    const lsdtt_raster_cpu& dem,
    spectral_window_type    window,
    int                     n_threads
);

static background_spectrum_product construct_background_spectrum(
    const lsdtt_raster_cpu&   dem,
    const spectral_psd_product& observed_psd,
    double                    sub_rollover_beta,
    int                       n_iterations,
    spectral_window_type      window,
    int                       n_threads
);


// -----------------------------------------------------------------------------
// compute_periodogram
//
// Computes the shifted 2D periodogram and the derived radial PSD for a DEM.
//
// Workflow:
//   1. Fit and subtract planar trend from valid pixels
//   2. Apply selected 2D window over the raster footprint
//   3. Zero-pad to next power-of-2 dimensions
//   4. Forward complex FFT
//   5. Compute shifted 2D PSD using legacy-compatible scaling:
//          PSD = |F|^2 / (Ly * Lx * WSS)
//   6. Collapse shifted PSD to radial PSD
//
// Notes:
//   - Invalid / nodata pixels are set to zero before FFT.
//   - detrended_variance is computed from the detrended valid pixels
//     before windowing, matching the old background-spectrum workflow.
// -----------------------------------------------------------------------------
static spectral_psd_product compute_periodogram(
    const lsdtt_raster_cpu& dem,
    spectral_window_type    window,
    int                     n_threads)
{
    (void)n_threads;

    const lsdtt_raster_meta& meta = dem.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: invalid raster metadata in compute_periodogram"
        );
    }

    meta.check_projected_metres();

    const int    n_rows    = meta.n_rows;
    const int    n_cols    = meta.n_cols;
    const double cell_size = meta.cell_size;
    const float  nodata    = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   data = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    const std::size_t n_pixels =
        checked_mul(
            static_cast<std::size_t>(n_rows),
            static_cast<std::size_t>(n_cols),
            "compute_periodogram input pixel count"
        );

    if (data.size() != n_pixels || mask.size() != n_pixels)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: data/mask size mismatch in compute_periodogram"
        );
    }

    // -------------------------------------------------------------------------
    // 1. Detrend valid pixels using fitted plane z = a*col + b*row + c
    // -------------------------------------------------------------------------
    double a_tr = 0.0;
    double b_tr = 0.0;
    double c_tr = 0.0;

    fit_detrend_plane(data, mask, meta, a_tr, b_tr, c_tr);

    std::vector<float> detrended(n_pixels, nodata);

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx = flat_index(row, col, n_cols);

            if (mask[idx] == 0) continue;
            if (meta.is_nodata(data[idx])) continue;

            const double trend =
                a_tr * static_cast<double>(col) +
                b_tr * static_cast<double>(row) +
                c_tr;

            const double z_dt = static_cast<double>(data[idx]) - trend;
            detrended[idx] = static_cast<float>(z_dt);
        }
    }

    const std::pair<double, double> mean_var =
        compute_valid_mean_and_variance(detrended, mask, meta);

    const double detrended_variance = mean_var.second;

    // -------------------------------------------------------------------------
    // 2. Build 2D window and WSS over the original raster footprint
    // -------------------------------------------------------------------------
    const std::vector<double> window_weights =
        build_2d_window(n_rows, n_cols, window);

    if (window_weights.size() != n_pixels)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: window size mismatch in compute_periodogram"
        );
    }

    const double window_sum_squared =
        compute_window_sum_squared(window_weights);

    if (!std::isfinite(window_sum_squared) || window_sum_squared <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: invalid WSS in compute_periodogram"
        );
    }

    // -------------------------------------------------------------------------
    // 3. Pad to next power-of-2 dimensions
    // -------------------------------------------------------------------------
    const int padded_rows = next_pow2(n_rows);
    const int padded_cols = next_pow2(n_cols);

    const double dfx = 1.0 / (cell_size * static_cast<double>(padded_cols));
    const double dfy = 1.0 / (cell_size * static_cast<double>(padded_rows));
    const double nyquist_frequency = 1.0 / (2.0 * cell_size);

    const std::size_t n_fft =
        checked_mul(
            static_cast<std::size_t>(padded_rows),
            static_cast<std::size_t>(padded_cols),
            "compute_periodogram FFT element count"
        );

    // -------------------------------------------------------------------------
    // 4. Full complex FFTW buffers and forward plan
    //
    // Important:
    //   create the plan before filling the input buffer because FFTW_MEASURE
    //   may overwrite working arrays during planning.
    // -------------------------------------------------------------------------
    fftw_complex_buffer fft_in (n_fft, "periodogram fft_in");
    fftw_complex_buffer fft_out(n_fft, "periodogram fft_out");


    fftw_plan_handle plan_fwd(
        padded_rows,
        padded_cols,
        fft_in.get(),
        fft_out.get(),
        FFTW_FORWARD,
        FFTW_MEASURE,
        "periodogram forward transform"
    );

    std::memset(fft_in.get(), 0, fft_in.bytes());

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t src_idx = flat_index(row, col, n_cols);

            if (mask[src_idx] == 0) continue;
            if (meta.is_nodata(detrended[src_idx])) continue;

            const std::size_t dst_idx = flat_index(row, col, padded_cols);

            fft_in.get()[dst_idx][0] =
                static_cast<double>(detrended[src_idx]) * window_weights[src_idx];
            fft_in.get()[dst_idx][1] = 0.0;
        }
    }

    fftw_execute(plan_fwd.get());

    // -------------------------------------------------------------------------
    // 5. Compute shifted 2D PSD
    //
    // Legacy-compatible scaling:
    //   P_DFT = |F|^2 / (Ly * Lx * WSS)
    //
    // We build the PSD directly into shifted order so that zero frequency lies
    // at (padded_rows/2, padded_cols/2), matching the old LSDRasterSpectral
    // workflow.
    // -------------------------------------------------------------------------
    std::vector<double> psd_2d_flat(n_fft, 0.0);

    const int half_rows = padded_rows / 2;
    const int half_cols = padded_cols / 2;
    const double denom =
        static_cast<double>(padded_rows) *
        static_cast<double>(padded_cols) *
        window_sum_squared;

    for (int row = 0; row < padded_rows; ++row)
    {
        for (int col = 0; col < padded_cols; ++col)
        {
            const std::size_t src_idx = flat_index(row, col, padded_cols);

            const double re = fft_out.get()[src_idx][0];
            const double im = fft_out.get()[src_idx][1];
            const double power = (re * re + im * im) / denom;

            const int shifted_row = (row + half_rows) % padded_rows;
            const int shifted_col = (col + half_cols) % padded_cols;
            const std::size_t dst_idx =
                flat_index(shifted_row, shifted_col, padded_cols);

            psd_2d_flat[dst_idx] = power;
        }
    }

    // -------------------------------------------------------------------------
    // 6. Collapse shifted PSD to radial PSD
    // -------------------------------------------------------------------------
    const std::pair<std::vector<double>, std::vector<double>> radial =
        collapse_shifted_psd_to_radial(
            psd_2d_flat,
            padded_rows,
            padded_cols,
            dfx,
            dfy,
            nyquist_frequency
        );

    spectral_psd_product out;
    out.padded_rows = padded_rows;
    out.padded_cols = padded_cols;
    out.dfx = dfx;
    out.dfy = dfy;
    out.nyquist_frequency = nyquist_frequency;
    out.window_sum_squared = window_sum_squared;
    out.detrended_variance = detrended_variance;
    out.psd_2d_flat = std::move(psd_2d_flat);
    out.radial_frequency = std::move(radial.first);
    out.radial_psd = std::move(radial.second);

    return out;
}

// -----------------------------------------------------------------------------
// construct_background_spectrum
//
// Builds an averaged background spectrum by:
//
//   1. Generating synthetic power-law fractal DEMs
//   2. Matching their detrended standard deviation to the observed DEM
//   3. Computing the same periodogram product for each synthetic DEM
//   4. Averaging the 2D PSD and radial PSD across iterations
//
// Notes:
//   - The synthetic DEM generator is intentionally delegated to
//     lsdtt_raster_synthetic_surfaces_cpu.
//   - This helper assumes the observed periodogram has already been computed.
//   - Radial frequencies should be identical between observed and synthetic
//     runs because raster dimensions, padding, and PSD collapse rules match.
// -----------------------------------------------------------------------------
static background_spectrum_product construct_background_spectrum(
    const lsdtt_raster_cpu&     dem,
    const spectral_psd_product& observed_psd,
    double                      sub_rollover_beta,
    int                         n_iterations,
    spectral_window_type        window,
    int                         n_threads)
{
    if (n_iterations < 1)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: n_iterations must be >= 1 in "
            "construct_background_spectrum"
        );
    }

    if (!std::isfinite(sub_rollover_beta) || sub_rollover_beta <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: sub_rollover_beta must be finite and > 0 "
            "in construct_background_spectrum"
        );
    }

    const lsdtt_raster_meta& meta = dem.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: invalid raster metadata in "
            "construct_background_spectrum"
        );
    }

    meta.check_projected_metres();

    if (observed_psd.padded_rows <= 0 || observed_psd.padded_cols <= 0)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: invalid observed padded dimensions in "
            "construct_background_spectrum"
        );
    }

    const std::size_t n_psd_2d =
        checked_mul(
            static_cast<std::size_t>(observed_psd.padded_rows),
            static_cast<std::size_t>(observed_psd.padded_cols),
            "construct_background_spectrum background 2D PSD size"
        );

    if (observed_psd.psd_2d_flat.size() != n_psd_2d)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: observed 2D PSD size mismatch in "
            "construct_background_spectrum"
        );
    }

    const std::size_t n_radial = observed_psd.radial_psd.size();

    if (observed_psd.radial_frequency.size() != n_radial)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: observed radial PSD/frequency size mismatch "
            "in construct_background_spectrum"
        );
    }

    if (!std::isfinite(observed_psd.detrended_variance) ||
        observed_psd.detrended_variance < 0.0)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: invalid observed detrended variance in "
            "construct_background_spectrum"
        );
    }

    std::vector<double> background_psd_2d_accum(n_psd_2d, 0.0);
    std::vector<double> background_radial_psd_accum(n_radial, 0.0);
    std::vector<double> background_radial_frequency = observed_psd.radial_frequency;

    const double target_stddev = std::sqrt(observed_psd.detrended_variance);

    for (int iter = 0; iter < n_iterations; ++iter)
    {
        // Use a simple deterministic seed mapping for repeatability.
        const int seed = 1000 + iter;

        lsdtt_raster_cpu synthetic_dem =
            generate_power_law_fractal_surface_fft(
                meta,
                seed,
                sub_rollover_beta,
                0.0,             // target_mean
                target_stddev,   // target_stddev
                n_threads
            );

        spectral_psd_product synthetic_psd =
            compute_periodogram(
                synthetic_dem,
                window,
                n_threads
            );

        if (synthetic_psd.padded_rows != observed_psd.padded_rows ||
            synthetic_psd.padded_cols != observed_psd.padded_cols)
        {
            throw std::runtime_error(
                "lsdtt spectral analysis: synthetic padded dimensions do not match "
                "observed padded dimensions in construct_background_spectrum"
            );
        }

        if (synthetic_psd.psd_2d_flat.size() != n_psd_2d)
        {
            throw std::runtime_error(
                "lsdtt spectral analysis: synthetic 2D PSD size mismatch in "
                "construct_background_spectrum"
            );
        }

        if (synthetic_psd.radial_psd.size() != n_radial ||
            synthetic_psd.radial_frequency.size() != n_radial)
        {
            throw std::runtime_error(
                "lsdtt spectral analysis: synthetic radial PSD size mismatch in "
                "construct_background_spectrum"
            );
        }

        for (std::size_t i = 0; i < n_radial; ++i)
        {
            const double f_obs = observed_psd.radial_frequency[i];
            const double f_syn = synthetic_psd.radial_frequency[i];

            if (!std::isfinite(f_obs) || !std::isfinite(f_syn) ||
                std::abs(f_obs - f_syn) > 1e-12)
            {
                throw std::runtime_error(
                    "lsdtt spectral analysis: synthetic radial frequency mismatch "
                    "in construct_background_spectrum"
                );
            }
        }

        for (std::size_t i = 0; i < n_psd_2d; ++i)
        {
            background_psd_2d_accum[i] += synthetic_psd.psd_2d_flat[i];
        }

        for (std::size_t i = 0; i < n_radial; ++i)
        {
            background_radial_psd_accum[i] += synthetic_psd.radial_psd[i];
        }
    }

    const double inv_n = 1.0 / static_cast<double>(n_iterations);

    for (double& v : background_psd_2d_accum)
    {
        v *= inv_n;
    }

    for (double& v : background_radial_psd_accum)
    {
        v *= inv_n;
    }

    background_spectrum_product out;
    out.padded_rows = observed_psd.padded_rows;
    out.padded_cols = observed_psd.padded_cols;
    out.dfx = observed_psd.dfx;
    out.dfy = observed_psd.dfy;
    out.nyquist_frequency = observed_psd.nyquist_frequency;
    out.background_psd_2d_flat = std::move(background_psd_2d_accum);
    out.background_radial_frequency = std::move(background_radial_frequency);
    out.background_radial_psd = std::move(background_radial_psd_accum);

    return out;
}



// -----------------------------------------------------------------------------
// rollover_fit_product
//
// Internal result from binned-spectrum rollover estimation.
// -----------------------------------------------------------------------------
struct rollover_fit_product
{
    double rollover_frequency = 0.0;
    double rollover_beta = 0.0;
    double sub_rollover_beta = 0.0;

    std::vector<double> binned_frequency;
    std::vector<double> binned_psd;
    std::vector<double> beta_by_bin;
    std::vector<double> r_sq_by_bin;
};

// -----------------------------------------------------------------------------
// estimate_rollover_from_binned_spectrum
//
// Legacy-inspired rollover finder:
//
//   1. Log-bin the radial PSD
//   2. Repeatedly fit linear regressions to suffixes in log-log space
//   3. Track beta = -slope and R^2 at each suffix start
//   4. Select rollover using max beta, with R^2 as fallback
//
// Returns:
//   - rollover_frequency
//   - rollover_beta
//   - sub_rollover_beta
//   - binned spectrum
//   - beta and R^2 traces by bin
// -----------------------------------------------------------------------------
static rollover_fit_product estimate_rollover_from_binned_spectrum(
    const std::vector<double>& radial_frequency,
    const std::vector<double>& radial_psd,
    double                     log_bin_width)
{
    rollover_fit_product out;

    const std::pair<std::vector<double>, std::vector<double>> binned =
        log_bin_mean_xy(radial_frequency, radial_psd, log_bin_width);

    out.binned_frequency = binned.first;
    out.binned_psd = binned.second;

    const std::size_t n_bins = out.binned_frequency.size();

    if (n_bins != out.binned_psd.size() || n_bins < 3)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: insufficient binned data for rollover fit"
        );
    }

    std::vector<double> log_f(n_bins, 0.0);
    std::vector<double> log_p(n_bins, 0.0);

    for (std::size_t i = 0; i < n_bins; ++i)
    {
        if (!std::isfinite(out.binned_frequency[i]) ||
            !std::isfinite(out.binned_psd[i]) ||
            out.binned_frequency[i] <= 0.0 ||
            out.binned_psd[i] <= 0.0)
        {
            throw std::runtime_error(
                "lsdtt spectral analysis: invalid binned spectrum values in "
                "estimate_rollover_from_binned_spectrum"
            );
        }

        log_f[i] = std::log10(out.binned_frequency[i]);
        log_p[i] = std::log10(out.binned_psd[i]);
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    out.beta_by_bin.assign(n_bins, nan);
    out.r_sq_by_bin.assign(n_bins, nan);

    std::size_t cutoff = n_bins - n_bins / 5;
    if (cutoff < 1) cutoff = 1;
    if (cutoff > n_bins) cutoff = n_bins;

    double best_beta = -std::numeric_limits<double>::infinity();
    double best_rsq  = -std::numeric_limits<double>::infinity();

    std::size_t best_beta_index = 0;
    std::size_t best_rsq_index  = 0;

    bool any_fit = false;

    for (std::size_t start = 0; start < cutoff; ++start)
    {
        std::vector<double> x(log_f.begin() + static_cast<std::ptrdiff_t>(start), log_f.end());
        std::vector<double> y(log_p.begin() + static_cast<std::ptrdiff_t>(start), log_p.end());

        const lsdtt_linear_regression_result fit = compute_linear_regression(x, y);

        if (!fit.has_finite_slope) continue;

        any_fit = true;

        const double beta_here = -fit.slope;

        out.beta_by_bin[start] = beta_here;
        out.r_sq_by_bin[start] = fit.r_squared;

        if (std::isfinite(beta_here) && beta_here > best_beta)
        {
            best_beta = beta_here;
            best_beta_index = start;
        }

        if (std::isfinite(fit.r_squared) && fit.r_squared > best_rsq)
        {
            best_rsq = fit.r_squared;
            best_rsq_index = start;
        }
    }

    if (!any_fit)
    {
        throw std::runtime_error(
            "lsdtt spectral analysis: no valid regression fits in rollover estimation"
        );
    }

    std::size_t rollover_index = best_beta_index;

    if (rollover_index == 0 || rollover_index >= n_bins - 1)
    {
        rollover_index = best_rsq_index;
    }

    if (rollover_index == 0 && n_bins > 1)
    {
        rollover_index = n_bins - 1;
    }

    out.rollover_frequency = out.binned_frequency[rollover_index];

    if (std::isfinite(out.beta_by_bin[rollover_index]))
    {
        out.rollover_beta = out.beta_by_bin[rollover_index];
    }
    else if (std::isfinite(best_beta))
    {
        out.rollover_beta = best_beta;
    }

    // -------------------------------------------------------------------------
    // Sub-rollover fit:
    // mimic the old code roughly by skipping the first two lowest-frequency bins
    // when enough bins are available.
    // -------------------------------------------------------------------------
    std::size_t sub_begin = 0;
    if (rollover_index > 2)
    {
        sub_begin = 2;
    }

    std::size_t sub_end = rollover_index + 1;

    if (sub_end > sub_begin && (sub_end - sub_begin) >= 2)
    {
        std::vector<double> x(
            log_f.begin() + static_cast<std::ptrdiff_t>(sub_begin),
            log_f.begin() + static_cast<std::ptrdiff_t>(sub_end)
        );
        std::vector<double> y(
            log_p.begin() + static_cast<std::ptrdiff_t>(sub_begin),
            log_p.begin() + static_cast<std::ptrdiff_t>(sub_end)
        );

        const lsdtt_linear_regression_result fit = compute_linear_regression(x, y);

        if (fit.has_finite_slope)
        {
            out.sub_rollover_beta = -fit.slope;
        }
    }

    return out;
}





} // unnamed namespace





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




// =============================================================================
// Public spectral filter functions
// =============================================================================

// -----------------------------------------------------------------------------
// fftw2D_lowpass_filter
//
// H(f) = exp( -ln(2) * (f / f_c)^2 )   with f_c = 1 / lambda_cutoff_metres
//
// H(0)   = 1.0   (DC and low frequencies pass unchanged)
// H(f_c) = 0.5   (-3 dB at the cutoff)
// H >> f_c → 0  (high frequencies suppressed)
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fftw2D_lowpass_filter(
    const lsdtt_raster_cpu& dem,
    double               lambda_cutoff_metres,
    spectral_window_type window,
    bool                 clip_to_valid_rectangle,
    int                  /*n_threads*/)
{
    if (!std::isfinite(lambda_cutoff_metres) || lambda_cutoff_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::fftw2D_lowpass_filter: lambda_cutoff_metres must be > 0"
        );
    }

    dem.get_meta().check_projected_metres();

    // Perform rectangle clip once, before the pipeline, if requested.
    const lsdtt_raster_cpu clipped_storage =
        clip_to_valid_rectangle ? clip_to_maximum_rectangle(dem) : lsdtt_raster_cpu{};
    const lsdtt_raster_cpu& work = clip_to_valid_rectangle ? clipped_storage : dem;

    const double f_c    = 1.0 / lambda_cutoff_metres;
    const double ln2    = std::log(2.0);
    const double inv_fc2 = 1.0 / (f_c * f_c);

    return spectral_filter_pipeline(
        work, window,
        [ln2, inv_fc2](
            fftw_complex* spec,
            int Ly, int /*Lx*/, int Lx_half,
            double dfx, double dfy)
        {
            for (int i = 0; i < Ly; ++i)
            {
                for (int j = 0; j < Lx_half; ++j)
                {
                    const double f = radial_freq_r2c(i, j, Ly, dfx, dfy);
                    const double H = (f == 0.0)
                        ? 1.0
                        : std::exp(-ln2 * f * f * inv_fc2);

                    const std::size_t k =
                        static_cast<std::size_t>(i) * Lx_half + j;
                    spec[k][0] *= H;
                    spec[k][1] *= H;
                }
            }
        }
    );
}

// -----------------------------------------------------------------------------
// fftw2D_highpass_filter
//
// H(f) = 1 - exp( -ln(2) * (f / f_c)^2 )
//
// H(0) = 0   (DC removed)
// H(f_c) = 0.5   (-3 dB at cutoff)
// H >> f_c → 1   (high frequencies pass unchanged)
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fftw2D_highpass_filter(
    const lsdtt_raster_cpu& dem,
    double               lambda_cutoff_metres,
    spectral_window_type window,
    bool                 clip_to_valid_rectangle,
    int                  /*n_threads*/)
{
    if (!std::isfinite(lambda_cutoff_metres) || lambda_cutoff_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::fftw2D_highpass_filter: lambda_cutoff_metres must be > 0"
        );
    }

    dem.get_meta().check_projected_metres();

    const lsdtt_raster_cpu clipped_storage =
        clip_to_valid_rectangle ? clip_to_maximum_rectangle(dem) : lsdtt_raster_cpu{};
    const lsdtt_raster_cpu& work = clip_to_valid_rectangle ? clipped_storage : dem;

    const double f_c     = 1.0 / lambda_cutoff_metres;
    const double ln2     = std::log(2.0);
    const double inv_fc2 = 1.0 / (f_c * f_c);

    return spectral_filter_pipeline(
        work, window,
        [ln2, inv_fc2](
            fftw_complex* spec,
            int Ly, int /*Lx*/, int Lx_half,
            double dfx, double dfy)
        {
            for (int i = 0; i < Ly; ++i)
            {
                for (int j = 0; j < Lx_half; ++j)
                {
                    const double f = radial_freq_r2c(i, j, Ly, dfx, dfy);
                    const double H = (f == 0.0)
                        ? 0.0
                        : 1.0 - std::exp(-ln2 * f * f * inv_fc2);

                    const std::size_t k =
                        static_cast<std::size_t>(i) * Lx_half + j;
                    spec[k][0] *= H;
                    spec[k][1] *= H;
                }
            }
        }
    );
}

// -----------------------------------------------------------------------------
// fftw2D_bandpass_filter
//
// H(f) = H_high(f, f_low) * H_low(f, f_high)
//
// where:
//   f_low  = 1 / lambda_long_metres   (low-frequency edge, long-wavelength cutoff)
//   f_high = 1 / lambda_short_metres  (high-frequency edge, short-wavelength cutoff)
//
// H_high(f, f_low)  = 1 - exp(-ln2 * (f/f_low)^2)   — passes f >> f_low
// H_low(f, f_high)  =     exp(-ln2 * (f/f_high)^2)   — passes f << f_high
//
// The combined response is 0.25 (-6 dB) at each edge and peaks near the
// geometric mean of f_low and f_high.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fftw2D_bandpass_filter(
    const lsdtt_raster_cpu& dem,
    double               lambda_short_metres,
    double               lambda_long_metres,
    spectral_window_type window,
    bool                 clip_to_valid_rectangle,
    int                  /*n_threads*/)
{
    if (!std::isfinite(lambda_short_metres) || lambda_short_metres <= 0.0 ||
        !std::isfinite(lambda_long_metres)  || lambda_long_metres  <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::fftw2D_bandpass_filter: wavelengths must be > 0"
        );
    }
    if (lambda_short_metres >= lambda_long_metres)
    {
        throw std::runtime_error(
            "lsdtt::fftw2D_bandpass_filter: lambda_short_metres must be < "
            "lambda_long_metres"
        );
    }

    dem.get_meta().check_projected_metres();

    const lsdtt_raster_cpu clipped_storage =
        clip_to_valid_rectangle ? clip_to_maximum_rectangle(dem) : lsdtt_raster_cpu{};
    const lsdtt_raster_cpu& work = clip_to_valid_rectangle ? clipped_storage : dem;

    // f_low  = 1/lambda_long  (lower frequency edge — long-wavelength cutoff)
    // f_high = 1/lambda_short (higher frequency edge — short-wavelength cutoff)
    const double f_low   = 1.0 / lambda_long_metres;
    const double f_high  = 1.0 / lambda_short_metres;
    const double ln2     = std::log(2.0);
    const double inv_fl2 = 1.0 / (f_low  * f_low);
    const double inv_fh2 = 1.0 / (f_high * f_high);

    return spectral_filter_pipeline(
        work, window,
        [ln2, f_low, inv_fl2, inv_fh2](
            fftw_complex* spec,
            int Ly, int /*Lx*/, int Lx_half,
            double dfx, double dfy)
        {
            for (int i = 0; i < Ly; ++i)
            {
                for (int j = 0; j < Lx_half; ++j)
                {
                    const double f = radial_freq_r2c(i, j, Ly, dfx, dfy);
                    double H;
                    if (f == 0.0)
                    {
                        H = 0.0;
                    }
                    else
                    {
                        const double H_hp =
                            1.0 - std::exp(-ln2 * f * f * inv_fl2);
                        const double H_lp =
                            std::exp(-ln2 * f * f * inv_fh2);
                        H = H_hp * H_lp;
                    }

                    const std::size_t k =
                        static_cast<std::size_t>(i) * Lx_half + j;
                    spec[k][0] *= H;
                    spec[k][1] *= H;
                }
            }
            (void)f_low;  // captured but only used via inv_fl2
        }
    );
}

// -----------------------------------------------------------------------------
// fftw2D_wiener_filter
//
// Adaptive noise-removal filter.  Algorithm:
//
//   1. Forward FFT (inside spectral_filter_pipeline).
//   2. Compute radial PSD: P(f) = |F(f)|^2 / (Ly * Lx) for all r2c bins.
//   3. Collect points in the signal band [f_sig_low, f_sig_high] in log-log
//      space and fit a power law log(P) = m*log(f) + log(c).
//   4. Estimate white noise N from the mean PSD at f > 0.8 * f_nyquist
//      (the near-Nyquist high-frequency region where noise dominates).
//   5. Build Wiener weight phi(f) = S(f) / (S(f) + N)
//      where S(f) = c * f^m is the signal model.
//   6. Apply phi(f) to spectrum.
//
// The Wiener weight approaches 1 at low frequencies (signal >> noise) and 0
// at high frequencies (noise >> signal).  The transition band is determined
// automatically from the spectrum of the DEM itself, so no explicit cutoff
// needs to be specified.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fftw2D_wiener_filter(
    const lsdtt_raster_cpu& dem,
    double               signal_wavelength_long_m,
    double               signal_wavelength_short_m,
    spectral_window_type window,
    bool                 clip_to_valid_rectangle,
    int                  /*n_threads*/)
{
    if (!std::isfinite(signal_wavelength_long_m) || signal_wavelength_long_m <= 0.0 ||
        !std::isfinite(signal_wavelength_short_m)|| signal_wavelength_short_m <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::fftw2D_wiener_filter: signal wavelengths must be > 0"
        );
    }
    if (signal_wavelength_short_m >= signal_wavelength_long_m)
    {
        throw std::runtime_error(
            "lsdtt::fftw2D_wiener_filter: signal_wavelength_short_m must be < "
            "signal_wavelength_long_m"
        );
    }

    dem.get_meta().check_projected_metres();

    const lsdtt_raster_cpu clipped_storage =
        clip_to_valid_rectangle ? clip_to_maximum_rectangle(dem) : lsdtt_raster_cpu{};
    const lsdtt_raster_cpu& work = clip_to_valid_rectangle ? clipped_storage : dem;

    const double cell_size  = work.get_meta().cell_size;
    const double f_nyquist  = 1.0 / (2.0 * cell_size);
    const double f_sig_low  = 1.0 / signal_wavelength_long_m;   // lower freq of signal band
    const double f_sig_high = 1.0 / signal_wavelength_short_m;  // upper freq of signal band
    const double f_noise_lo = 0.8 * f_nyquist;                  // near-Nyquist noise region

    return spectral_filter_pipeline(
        work, window,
        [f_sig_low, f_sig_high, f_noise_lo, f_nyquist](
            fftw_complex* spec,
            int Ly, int Lx, int Lx_half,
            double dfx, double dfy)
        {
            const double inv_n = 1.0 / static_cast<double>(
                static_cast<std::size_t>(Ly) * Lx);

            // --- Step 1: compute radial PSD and identify signal / noise bins ---
            double sum_log_f   = 0.0, sum_log_p   = 0.0;
            double sum_log_f2  = 0.0, sum_log_fp  = 0.0;
            int    n_sig       = 0;

            double noise_sum   = 0.0;
            int    n_noise     = 0;

            for (int i = 0; i < Ly; ++i)
            {
                for (int j = 0; j < Lx_half; ++j)
                {
                    const double f = radial_freq_r2c(i, j, Ly, dfx, dfy);
                    if (f == 0.0 || f > f_nyquist) continue;

                    const std::size_t k =
                        static_cast<std::size_t>(i) * Lx_half + j;
                    const double re  = spec[k][0];
                    const double im  = spec[k][1];
                    // Factor-of-2 for the folded half-spectrum (except DC and Nyquist)
                    const double deg = (j == 0 || j == Lx_half - 1) ? 1.0 : 2.0;
                    const double psd = deg * (re*re + im*im) * inv_n;

                    // Signal band — log-log regression
                    if (f >= f_sig_low && f <= f_sig_high && psd > 0.0)
                    {
                        const double lf = std::log10(f);
                        const double lp = std::log10(psd);
                        sum_log_f  += lf;
                        sum_log_p  += lp;
                        sum_log_f2 += lf * lf;
                        sum_log_fp += lf * lp;
                        ++n_sig;
                    }

                    // Noise region — near Nyquist
                    if (f >= f_noise_lo)
                    {
                        noise_sum += psd;
                        ++n_noise;
                    }
                }
            }

            // --- Step 2: fit power law log(P) = m*log(f) + log(c) ---
            double c_model = 1.0;
            double m_model = -2.0;   // typical topographic spectral slope

            if (n_sig >= 2)
            {
                const double fn   = static_cast<double>(n_sig);
                const double denom = fn * sum_log_f2 - sum_log_f * sum_log_f;
                if (std::abs(denom) > 1e-30)
                {
                    m_model  = (fn * sum_log_fp - sum_log_f * sum_log_p) / denom;
                    const double log_c = (sum_log_p - m_model * sum_log_f) / fn;
                    c_model  = std::pow(10.0, log_c);
                }
            }

            // --- Step 3: estimate white noise level ---
            double noise_level = (n_noise > 0)
                ? noise_sum / static_cast<double>(n_noise)
                : 0.0;

            // Guard: if noise estimate is zero or very small relative to
            // signal, set it to a tiny fraction of the DC component so the
            // Wiener filter doesn't become a unity passthrough everywhere.
            if (noise_level <= 0.0)
            {
                const double psd_dc =
                    (spec[0][0]*spec[0][0] + spec[0][1]*spec[0][1]) * inv_n;
                noise_level = std::max(psd_dc * 1e-6, 1e-30);
            }

            // --- Step 4: apply Wiener weights ---
            for (int i = 0; i < Ly; ++i)
            {
                for (int j = 0; j < Lx_half; ++j)
                {
                    const double f = radial_freq_r2c(i, j, Ly, dfx, dfy);

                    double phi;
                    if (f == 0.0)
                    {
                        phi = 1.0;   // keep DC
                    }
                    else
                    {
                        const double S = c_model * std::pow(f, m_model);
                        phi = S / (S + noise_level);
                        phi = std::max(0.0, std::min(1.0, phi));
                    }

                    const std::size_t k =
                        static_cast<std::size_t>(i) * Lx_half + j;
                    spec[k][0] *= phi;
                    spec[k][1] *= phi;
                }
            }
        }
    );
}


// -----------------------------------------------------------------------------
// full_spectral_analysis
// -----------------------------------------------------------------------------
lsdtt_spectral_analysis_result full_spectral_analysis(
    const lsdtt_raster_cpu&                dem,
    const lsdtt_spectral_analysis_options& opts)
{
    if (!dem.get_meta().is_valid())
    {
        throw std::runtime_error(
            "lsdtt::full_spectral_analysis: invalid raster metadata"
        );
    }

    dem.get_meta().check_projected_metres();

    if (!std::isfinite(opts.log_bin_width) || opts.log_bin_width <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::full_spectral_analysis: log_bin_width must be finite and > 0"
        );
    }

    if (opts.n_background_iterations < 1)
    {
        throw std::runtime_error(
            "lsdtt::full_spectral_analysis: n_background_iterations must be >= 1"
        );
    }

    // Optional crop to largest valid rectangle before any spectral work.
    const lsdtt_raster_cpu clipped_storage =
        opts.clip_to_valid_rectangle ? clip_to_maximum_rectangle(dem)
                                     : lsdtt_raster_cpu{};
    const lsdtt_raster_cpu& work =
        opts.clip_to_valid_rectangle ? clipped_storage : dem;

    // -------------------------------------------------------------------------
    // 1. Observed periodogram
    // -------------------------------------------------------------------------
    const spectral_psd_product observed_psd =
        compute_periodogram(
            work,
            opts.window,
            opts.n_threads
        );

    // -------------------------------------------------------------------------
    // 2. Rollover estimation from binned radial spectrum
    // -------------------------------------------------------------------------
    const rollover_fit_product rollover =
        estimate_rollover_from_binned_spectrum(
            observed_psd.radial_frequency,
            observed_psd.radial_psd,
            opts.log_bin_width
        );

    // Use the sub-rollover beta if available; otherwise fall back to the main
    // rollover beta, and finally to a conservative default.
    double background_beta = rollover.sub_rollover_beta;
    if (!std::isfinite(background_beta) || background_beta <= 0.0)
    {
        background_beta = rollover.rollover_beta;
    }
    if (!std::isfinite(background_beta) || background_beta <= 0.0)
    {
        background_beta = 2.0;
    }

    // -------------------------------------------------------------------------
    // 3. Monte Carlo background spectrum
    // -------------------------------------------------------------------------
    const background_spectrum_product background =
        construct_background_spectrum(
            work,
            observed_psd,
            background_beta,
            opts.n_background_iterations,
            opts.window,
            opts.n_threads
        );

    // -------------------------------------------------------------------------
    // 4. Validate observed/background frequency compatibility
    // -------------------------------------------------------------------------
    const std::size_t n_radial = observed_psd.radial_frequency.size();

    if (observed_psd.radial_psd.size() != n_radial ||
        background.background_radial_frequency.size() != n_radial ||
        background.background_radial_psd.size() != n_radial)
    {
        throw std::runtime_error(
            "lsdtt::full_spectral_analysis: radial spectrum size mismatch"
        );
    }

    for (std::size_t i = 0; i < n_radial; ++i)
    {
        const double f_obs = observed_psd.radial_frequency[i];
        const double f_bg  = background.background_radial_frequency[i];

        if (!std::isfinite(f_obs) || !std::isfinite(f_bg) ||
            std::abs(f_obs - f_bg) > 1e-12)
        {
            throw std::runtime_error(
                "lsdtt::full_spectral_analysis: radial frequency mismatch between "
                "observed and background spectra"
            );
        }
    }

    // -------------------------------------------------------------------------
    // 5. Populate result
    // -------------------------------------------------------------------------
    lsdtt_spectral_analysis_result result;

    result.padded_rows = observed_psd.padded_rows;
    result.padded_cols = observed_psd.padded_cols;
    result.dfx = observed_psd.dfx;
    result.dfy = observed_psd.dfy;
    result.nyquist_frequency = observed_psd.nyquist_frequency;
    result.window_sum_squared = observed_psd.window_sum_squared;
    result.detrended_variance = observed_psd.detrended_variance;

    result.rollover_frequency = rollover.rollover_frequency;
    result.rollover_beta = rollover.rollover_beta;
    result.sub_rollover_beta = background_beta;

    // Convert flat 2D PSD vectors into array2d outputs
    result.psd_2d = lsdtt_array2d_f(
        observed_psd.padded_rows,
        observed_psd.padded_cols,
        0.0f
    );

    result.background_psd_2d = lsdtt_array2d_f(
        background.padded_rows,
        background.padded_cols,
        0.0f
    );

    for (int row = 0; row < observed_psd.padded_rows; ++row)
    {
        for (int col = 0; col < observed_psd.padded_cols; ++col)
        {
            const std::size_t idx =
                flat_index(row, col, observed_psd.padded_cols);

            result.psd_2d(row, col) =
                static_cast<float>(observed_psd.psd_2d_flat[idx]);

            result.background_psd_2d(row, col) =
                static_cast<float>(background.background_psd_2d_flat[idx]);
        }
    }

    result.radial_frequency = observed_psd.radial_frequency;
    result.radial_psd = observed_psd.radial_psd;

    result.radial_wavelength.resize(n_radial, 0.0);
    result.background_psd.resize(n_radial, 0.0);
    result.normalised_psd.resize(n_radial, 0.0);
    result.ci95.resize(n_radial, 0.0);
    result.norm_ci95.resize(n_radial, 0.0);
    result.norm_ci99.resize(n_radial, 0.0);

    const double chi_squared_95 = 5.991;
    const double chi_squared_99 = 9.210;
    const double norm_ci95 = 0.5 * chi_squared_95;
    const double norm_ci99 = 0.5 * chi_squared_99;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (std::size_t i = 0; i < n_radial; ++i)
    {
        const double f   = observed_psd.radial_frequency[i];
        const double psd = observed_psd.radial_psd[i];
        const double bg  = background.background_radial_psd[i];

        result.radial_wavelength[i] =
            (f > 0.0) ? (1.0 / f) : std::numeric_limits<double>::infinity();

        result.background_psd[i] = bg;

        if (std::isfinite(psd) && std::isfinite(bg) && bg > 0.0)
        {
            result.normalised_psd[i] = psd / bg;
            result.ci95[i] = norm_ci95 * bg;
        }
        else
        {
            result.normalised_psd[i] = nan;
            result.ci95[i] = nan;
        }

        result.norm_ci95[i] = norm_ci95;
        result.norm_ci99[i] = norm_ci99;
    }

    // -------------------------------------------------------------------------
    // 6. Shared-valid binned outputs
    //
    // To keep the binned frequency grid consistent across PSD, background PSD,
    // and normalised PSD, only retain raw points that are valid and positive in
    // all three series before log-binning.
    // -------------------------------------------------------------------------
    std::vector<double> shared_freq;
    std::vector<double> shared_obs_psd;
    std::vector<double> shared_bg_psd;
    std::vector<double> shared_norm_psd;

    shared_freq.reserve(n_radial);
    shared_obs_psd.reserve(n_radial);
    shared_bg_psd.reserve(n_radial);
    shared_norm_psd.reserve(n_radial);

    for (std::size_t i = 0; i < n_radial; ++i)
    {
        const double f  = result.radial_frequency[i];
        const double p  = result.radial_psd[i];
        const double bg = result.background_psd[i];
        const double np = result.normalised_psd[i];

        if (!std::isfinite(f)  || f  <= 0.0) continue;
        if (!std::isfinite(p)  || p  <= 0.0) continue;
        if (!std::isfinite(bg) || bg <= 0.0) continue;
        if (!std::isfinite(np) || np <= 0.0) continue;

        shared_freq.push_back(f);
        shared_obs_psd.push_back(p);
        shared_bg_psd.push_back(bg);
        shared_norm_psd.push_back(np);
    }

    const std::pair<std::vector<double>, std::vector<double>> binned_obs =
        log_bin_mean_xy(shared_freq, shared_obs_psd, opts.log_bin_width);

    const std::pair<std::vector<double>, std::vector<double>> binned_bg =
        log_bin_mean_xy(shared_freq, shared_bg_psd, opts.log_bin_width);

    const std::pair<std::vector<double>, std::vector<double>> binned_norm =
        log_bin_mean_xy(shared_freq, shared_norm_psd, opts.log_bin_width);

    if (binned_obs.first.size() != binned_bg.first.size() ||
        binned_obs.first.size() != binned_norm.first.size() ||
        binned_obs.second.size() != binned_bg.second.size() ||
        binned_obs.second.size() != binned_norm.second.size())
    {
        throw std::runtime_error(
            "lsdtt::full_spectral_analysis: binned spectrum size mismatch"
        );
    }

    result.binned.frequency = binned_obs.first;
    result.binned.psd = binned_obs.second;
    result.binned.background_psd = binned_bg.second;
    result.binned.normalised_psd = binned_norm.second;
    result.binned.beta_by_bin = rollover.beta_by_bin;
    result.binned.r_sq_by_bin = rollover.r_sq_by_bin;

    result.binned.wavelength.resize(result.binned.frequency.size(), 0.0);
    for (std::size_t i = 0; i < result.binned.frequency.size(); ++i)
    {
        const double f = result.binned.frequency[i];
        result.binned.wavelength[i] =
            (f > 0.0) ? (1.0 / f) : std::numeric_limits<double>::infinity();
    }

    return result;
}




#endif  // LSDTT_HAS_FFTW

} // namespace lsdtt
