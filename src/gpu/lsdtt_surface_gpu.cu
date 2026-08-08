// lsdtt_surface_gpu.cu
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// CUDA implementation of the GPU surface-fitting kernels declared in
// lsdtt_surface_gpu.hpp. Compiled into lsdtt3_gpu only when -DLSDTT_GPU=ON.
//
// NOTE (design-ahead): this file was written and reviewed against the CPU
// reference (lsdtt::direct_polyfit) but has not yet been compiled or run on a
// CUDA device. It must be parity-tested against the CPU path on real hardware
// before being relied upon. The default LSDTT_GPU=OFF build never compiles it.
// -----------------------------------------------------------------------------

#include "lsdtt_surface_gpu.hpp"

#include <cuda_runtime.h>
#include <cstddef>
#include <cmath>

namespace lsdtt {
namespace gpu {

namespace {

// Small RAII-free helper: free a list of device pointers (nullptr-safe).
inline void free_all(void* a, void* b, void* c, void* d,
                     void* e, void* f, void* g, void* h)
{
    if (a) cudaFree(a);
    if (b) cudaFree(b);
    if (c) cudaFree(c);
    if (d) cudaFree(d);
    if (e) cudaFree(e);
    if (f) cudaFree(f);
    if (g) cudaFree(g);
    if (h) cudaFree(h);
}

// Closed-form inverse of a real 3x3 matrix (row-major, 9 doubles).
// Returns false if (near-)singular. Used host-side to precompute the constant
// normal-matrix inverse before launch.
bool invert3x3_rowmajor(const double m[9], double inv[9])
{
    const double a = m[0], b = m[1], c = m[2];
    const double d = m[3], e = m[4], f = m[5];
    const double g = m[6], h = m[7], i = m[8];

    const double A =  (e * i - f * h);
    const double B = -(d * i - f * g);
    const double C =  (d * h - e * g);
    const double det = a * A + b * B + c * C;

    if (!(std::fabs(det) > 1e-300))
        return false;

    const double inv_det = 1.0 / det;

    inv[0] =  A * inv_det;
    inv[1] = -(b * i - c * h) * inv_det;
    inv[2] =  (b * f - c * e) * inv_det;
    inv[3] =  B * inv_det;
    inv[4] =  (a * i - c * g) * inv_det;
    inv[5] = -(a * f - c * d) * inv_det;
    inv[6] =  C * inv_det;
    inv[7] = -(a * h - b * g) * inv_det;
    inv[8] =  (a * e - b * d) * inv_det;
    return true;
}

// Constants needed inside the kernel. All precomputed on the host.
struct polyfit_constants
{
    int    nrows;
    int    ncols;
    int    W;
    float  cell_size;     // window x/y use double via promotion; float is plenty
    float  nodata;        // sentinel for the equality check
    double inv_Sx2;       // 1 / Sx2   (denominator for D and E)
    double inv_Sx2y2;     // 1 / Sx2y2 (denominator for C)
    double invM3[9];      // inverse of the 3x3 [a,b,F] normal matrix (row-major)
};

__device__ inline bool pixel_invalid(float v, unsigned char m, float nodata)
{
    return (m == 0) || (v == nodata) || isnan(v);
}

// Mark an output pixel as nodata across all six coefficient rasters + mask.
// The whole device output buffer is copied back to the host, so skipped pixels
// must be written explicitly rather than left as uninitialised device memory.
__device__ inline void write_nodata(
    long idx, float nodata,
    float* a, float* b, float* c, float* d, float* e, float* f,
    unsigned char* m)
{
    a[idx] = nodata; b[idx] = nodata; c[idx] = nodata;
    d[idx] = nodata; e[idx] = nodata; f[idx] = nodata;
    m[idx] = 0;
}

// One thread per output pixel. Brute-force window gather; double accumulation
// for parity with the CPU prefix-sum implementation.
__global__ void direct_polyfit_kernel_bruteforce(
    const float*         z,
    const unsigned char* mask,
    polyfit_constants    k,
    float*               out_a,
    float*               out_b,
    float*               out_c,
    float*               out_d,
    float*               out_e,
    float*               out_f,
    unsigned char*       out_mask)
{
    const long idx = static_cast<long>(blockIdx.x) * blockDim.x + threadIdx.x;
    const long npix = static_cast<long>(k.nrows) * k.ncols;
    if (idx >= npix) return;

    const int r0 = static_cast<int>(idx / k.ncols);
    const int c0 = static_cast<int>(idx % k.ncols);

    // Interior-only: full window must lie inside the raster.
    if (r0 < k.W || r0 >= k.nrows - k.W ||
        c0 < k.W || c0 >= k.ncols - k.W)
    {
        write_nodata(idx, k.nodata,
                     out_a, out_b, out_c, out_d, out_e, out_f, out_mask);
        return;
    }

    const double cell = static_cast<double>(k.cell_size);

    double Sz = 0.0, Sxz = 0.0, Syz = 0.0;
    double Sxyz = 0.0, Sx2z = 0.0, Sy2z = 0.0;

    for (int dr = -k.W; dr <= k.W; ++dr)
    {
        const long base = static_cast<long>(r0 + dr) * k.ncols + c0;
        const double y = -static_cast<double>(dr) * cell;   // north positive

        for (int dc = -k.W; dc <= k.W; ++dc)
        {
            const long widx = base + dc;
            const float  vv = z[widx];
            if (pixel_invalid(vv, mask[widx], k.nodata))
            {
                // Any nodata in the window -> whole pixel is nodata.
                write_nodata(idx, k.nodata,
                             out_a, out_b, out_c, out_d, out_e, out_f, out_mask);
                return;
            }
            const double zz = static_cast<double>(vv);
            const double x  = static_cast<double>(dc) * cell;

            Sz   += zz;
            Sxz  += x * zz;
            Syz  += y * zz;
            Sxyz += x * y * zz;
            Sx2z += x * x * zz;
            Sy2z += y * y * zz;
        }
    }

    const double C = Sxyz * k.inv_Sx2y2;
    const double D = Sxz  * k.inv_Sx2;
    const double E = Syz  * k.inv_Sx2;

    const double a = k.invM3[0] * Sx2z + k.invM3[1] * Sy2z + k.invM3[2] * Sz;
    const double b = k.invM3[3] * Sx2z + k.invM3[4] * Sy2z + k.invM3[5] * Sz;
    const double F = k.invM3[6] * Sx2z + k.invM3[7] * Sy2z + k.invM3[8] * Sz;

    out_a[idx] = static_cast<float>(a);
    out_b[idx] = static_cast<float>(b);
    out_c[idx] = static_cast<float>(C);
    out_d[idx] = static_cast<float>(D);
    out_e[idx] = static_cast<float>(E);
    out_f[idx] = static_cast<float>(F);
    out_mask[idx] = 1;
}

// Pass 1 (rolling): horizontal moments per row/centre column.
//
// For each row r and centre c we build:
//   H0 = sum z
//   H1 = sum (dc)   z
//   H2 = sum (dc^2) z
//   Hn = count of invalid samples in [c-W, c+W]
// where dc = (j - c) in pixel units.
//
// Implementation detail: maintain rolling absolute-index sums
//   J0 = sum z, J1 = sum j z, J2 = sum j^2 z
// and convert to centred moments by
//   H1 = J1 - c J0
//   H2 = J2 - 2c J1 + c^2 J0
// This avoids per-window reweight loops.
__global__ void horizontal_moments_kernel(
    const float*         z,
    const unsigned char* mask,
    polyfit_constants    k,
    double*              H0,
    double*              H1,
    double*              H2,
    int*                 Hn)
{
    const int r = static_cast<int>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (r >= k.nrows) return;

    const int ncols = k.ncols;
    const int W = k.W;

    double J0 = 0.0;
    double J1 = 0.0;
    double J2 = 0.0;
    int invalid = 0;

    // Initialise the c=0 window [ -W, +W ].
    for (int j = -W; j <= W; ++j)
    {
        if (j < 0 || j >= ncols)
        {
            ++invalid;
            continue;
        }

        const long idx = static_cast<long>(r) * ncols + j;
        const float v = z[idx];
        if (pixel_invalid(v, mask[idx], k.nodata))
        {
            ++invalid;
            continue;
        }

        const double jd = static_cast<double>(j);
        const double zd = static_cast<double>(v);
        J0 += zd;
        J1 += jd * zd;
        J2 += jd * jd * zd;
    }

    for (int c = 0; c < ncols; ++c)
    {
        const long out = static_cast<long>(r) * ncols + c;
        const double cd = static_cast<double>(c);

        H0[out] = J0;
        H1[out] = J1 - cd * J0;
        H2[out] = J2 - 2.0 * cd * J1 + cd * cd * J0;
        Hn[out] = invalid;

        if (c == ncols - 1) break;

        // Slide centre c -> c+1: remove j=c-W, add j=c+W+1.
        const int j_rem = c - W;
        if (j_rem < 0 || j_rem >= ncols)
        {
            --invalid;
        }
        else
        {
            const long idx = static_cast<long>(r) * ncols + j_rem;
            const float v = z[idx];
            if (pixel_invalid(v, mask[idx], k.nodata))
            {
                --invalid;
            }
            else
            {
                const double jd = static_cast<double>(j_rem);
                const double zd = static_cast<double>(v);
                J0 -= zd;
                J1 -= jd * zd;
                J2 -= jd * jd * zd;
            }
        }

        const int j_add = c + W + 1;
        if (j_add < 0 || j_add >= ncols)
        {
            ++invalid;
        }
        else
        {
            const long idx = static_cast<long>(r) * ncols + j_add;
            const float v = z[idx];
            if (pixel_invalid(v, mask[idx], k.nodata))
            {
                ++invalid;
            }
            else
            {
                const double jd = static_cast<double>(j_add);
                const double zd = static_cast<double>(v);
                J0 += zd;
                J1 += jd * zd;
                J2 += jd * jd * zd;
            }
        }
    }
}

// Pass 2 (rolling): vertical combine to full 2D moments and solve coefficients.
// Complexity per pixel is O(W), not O(W^2).
__global__ void direct_polyfit_kernel_rolling(
    const double*        H0,
    const double*        H1,
    const double*        H2,
    const int*           Hn,
    polyfit_constants    k,
    float*               out_a,
    float*               out_b,
    float*               out_c,
    float*               out_d,
    float*               out_e,
    float*               out_f,
    unsigned char*       out_mask)
{
    const long idx = static_cast<long>(blockIdx.x) * blockDim.x + threadIdx.x;
    const long npix = static_cast<long>(k.nrows) * k.ncols;
    if (idx >= npix) return;

    const int r0 = static_cast<int>(idx / k.ncols);
    const int c0 = static_cast<int>(idx % k.ncols);

    if (r0 < k.W || r0 >= k.nrows - k.W ||
        c0 < k.W || c0 >= k.ncols - k.W)
    {
        write_nodata(idx, k.nodata,
                     out_a, out_b, out_c, out_d, out_e, out_f, out_mask);
        return;
    }

    double Sz = 0.0;
    double Sxzp = 0.0;
    double Sx2zp = 0.0;
    double Syzp = 0.0;
    double Sy2zp = 0.0;
    double Sxyzp = 0.0;
    int invalid = 0;

    for (int dr = -k.W; dr <= k.W; ++dr)
    {
        const int rr = r0 + dr;
        const long ridx = static_cast<long>(rr) * k.ncols + c0;

        const double h0 = H0[ridx];
        const double h1 = H1[ridx];
        const double h2 = H2[ridx];
        invalid += Hn[ridx];

        const double ddr = static_cast<double>(dr);
        Sz    += h0;
        Sxzp  += h1;
        Sx2zp += h2;
        Syzp  += (-ddr) * h0;
        Sy2zp += (ddr * ddr) * h0;
        Sxyzp += (-ddr) * h1;
    }

    if (invalid != 0)
    {
        write_nodata(idx, k.nodata,
                     out_a, out_b, out_c, out_d, out_e, out_f, out_mask);
        return;
    }

    const double cell = static_cast<double>(k.cell_size);
    const double d2 = cell * cell;

    const double Sxz  = cell * Sxzp;
    const double Syz  = cell * Syzp;
    const double Sx2z = d2 * Sx2zp;
    const double Sy2z = d2 * Sy2zp;
    const double Sxyz = d2 * Sxyzp;

    const double C = Sxyz * k.inv_Sx2y2;
    const double D = Sxz  * k.inv_Sx2;
    const double E = Syz  * k.inv_Sx2;

    const double a = k.invM3[0] * Sx2z + k.invM3[1] * Sy2z + k.invM3[2] * Sz;
    const double b = k.invM3[3] * Sx2z + k.invM3[4] * Sy2z + k.invM3[5] * Sz;
    const double F = k.invM3[6] * Sx2z + k.invM3[7] * Sy2z + k.invM3[8] * Sz;

    out_a[idx] = static_cast<float>(a);
    out_b[idx] = static_cast<float>(b);
    out_c[idx] = static_cast<float>(C);
    out_d[idx] = static_cast<float>(D);
    out_e[idx] = static_cast<float>(E);
    out_f[idx] = static_cast<float>(F);
    out_mask[idx] = 1;
}

// IDs must match lsdtt::polyfit_metric enum order in
// src/cpu/lsdtt_raster_surface_metrics_cpu.hpp.
enum polyfit_metric_id : int
{
    PM_FITTED_SURFACE = 0,
    PM_DZDX,
    PM_DZDY,
    PM_D2ZDX2,
    PM_D2ZDXDY,
    PM_D2ZDY2,
    PM_SLOPE,
    PM_SLOPE_DEGREES,
    PM_ASPECT_RADIANS,
    PM_ASPECT_DEGREES,
    PM_LAPLACIAN,
    PM_HESSIAN_DETERMINANT,
    PM_MEAN_CURVATURE,
    PM_GAUSSIAN_CURVATURE,
    PM_MAXIMUM_CURVATURE,
    PM_MINIMUM_CURVATURE,
    PM_CURVEDNESS,
    PM_SHAPE_INDEX,
    PM_PROFILE_CURVATURE,
    PM_PLAN_CURVATURE,
    PM_TANGENTIAL_CURVATURE,
    PM_DIFFERENCE_CURVATURE,
    PM_HORIZONTAL_EXCESS_CURVATURE,
    PM_VERTICAL_EXCESS_CURVATURE,
    PM_ACCUMULATION_CURVATURE,
    PM_RING_CURVATURE,
    PM_ROTOR
};

constexpr double k_pi_gpu = 3.14159265358979323846;

__global__ void compute_polyfit_metrics_kernel(
    const float* a_data,
    const float* b_data,
    const float* c_data,
    const float* d_data,
    const float* e_data,
    const float* f_data,
    const unsigned char* a_mask,
    const unsigned char* b_mask,
    const unsigned char* c_mask,
    const unsigned char* d_mask,
    const unsigned char* e_mask,
    const unsigned char* f_mask,
    int nrows,
    int ncols,
    float nodata,
    double z_factor,
    double flat_tol,
    const int* metric_ids,
    int n_metrics,
    float* out_data,
    unsigned char* out_mask)
{
    const long idx = static_cast<long>(blockIdx.x) * blockDim.x + threadIdx.x;
    const long npix = static_cast<long>(nrows) * ncols;
    if (idx >= npix) return;

    // Initialize all requested outputs at this pixel to nodata/mask0.
    for (int m = 0; m < n_metrics; ++m)
    {
        const long o = static_cast<long>(m) * npix + idx;
        out_data[o] = nodata;
        out_mask[o] = 0;
    }

    const bool coeff_valid =
        a_mask[idx] != 0 && b_mask[idx] != 0 && c_mask[idx] != 0 &&
        d_mask[idx] != 0 && e_mask[idx] != 0 && f_mask[idx] != 0 &&
        (a_data[idx] != nodata) && (b_data[idx] != nodata) &&
        (c_data[idx] != nodata) && (d_data[idx] != nodata) &&
        (e_data[idx] != nodata) && (f_data[idx] != nodata) &&
        isfinite(a_data[idx]) && isfinite(b_data[idx]) &&
        isfinite(c_data[idx]) && isfinite(d_data[idx]) &&
        isfinite(e_data[idx]) && isfinite(f_data[idx]);

    if (!coeff_valid) return;

    const double a = z_factor * static_cast<double>(a_data[idx]);
    const double b = z_factor * static_cast<double>(b_data[idx]);
    const double c = z_factor * static_cast<double>(c_data[idx]);
    const double d = z_factor * static_cast<double>(d_data[idx]);
    const double e = z_factor * static_cast<double>(e_data[idx]);
    const double f = z_factor * static_cast<double>(f_data[idx]);

    const double p = d;
    const double q = e;
    const double r = 2.0 * a;
    const double s = c;
    const double t = 2.0 * b;

    const double p2 = p * p;
    const double q2 = q * q;
    const double g2 = p2 + q2;
    const double g  = sqrt(g2);
    const double w2 = 1.0 + g2;
    const double w  = sqrt(w2);
    const double w3 = w2 * w;
    const double w4 = w2 * w2;

    double H = NAN;
    double K = NAN;
    if (w3 > 0.0)
        H = ((1.0 + q2) * r - 2.0 * p * q * s + (1.0 + p2) * t) / (2.0 * w3);
    if (w4 > 0.0)
        K = (r * t - s * s) / w4;

    const double disc = fmax(0.0, H * H - K);
    const double root = sqrt(disc);
    const double kmax = H + root;
    const double kmin = H - root;

    const double flat_tol_sq = flat_tol * flat_tol;
    const double denom_g3 = g2 * sqrt(g2);

    double plan_curv = NAN;
    if (denom_g3 > 0.0 && g2 > flat_tol_sq)
        plan_curv = (q2 * r - 2.0 * p * q * s + p2 * t) / denom_g3;

    double tan_curv = NAN;
    {
        const double denom = g2 * w;
        if (denom > 0.0 && g2 > flat_tol_sq)
            tan_curv = (q2 * r - 2.0 * p * q * s + p2 * t) / denom;
    }

    double vert_curv = NAN;
    {
        const double denom = g2 * w3;
        if (denom > 0.0 && g2 > flat_tol_sq)
            vert_curv = (p2 * r + 2.0 * p * q * s + q2 * t) / denom;
    }

    double rotor = NAN;
    if (denom_g3 > 0.0 && g2 > flat_tol_sq)
        rotor = ((p2 - q2) * s - p * q * (r - t)) / denom_g3;

    const double diff_curv = (isfinite(vert_curv) && isfinite(tan_curv))
        ? 0.5 * (vert_curv - tan_curv) : NAN;
    const double h_excess = (isfinite(tan_curv) && isfinite(kmin))
        ? (tan_curv - kmin) : NAN;
    const double v_excess = (isfinite(vert_curv) && isfinite(kmin))
        ? (vert_curv - kmin) : NAN;
    const double accum_curv = (isfinite(tan_curv) && isfinite(vert_curv))
        ? (tan_curv * vert_curv) : NAN;
    const double ring_curv = (isfinite(h_excess) && isfinite(v_excess))
        ? (h_excess * v_excess) : NAN;

    for (int m = 0; m < n_metrics; ++m)
    {
        const int metric = metric_ids[m];
        const long o = static_cast<long>(m) * npix + idx;

        double v = NAN;
        bool valid = false;

        switch (metric)
        {
            case PM_FITTED_SURFACE:      v = f; valid = true; break;
            case PM_DZDX:                v = p; valid = true; break;
            case PM_DZDY:                v = q; valid = true; break;
            case PM_D2ZDX2:              v = r; valid = true; break;
            case PM_D2ZDXDY:             v = s; valid = true; break;
            case PM_D2ZDY2:              v = t; valid = true; break;
            case PM_SLOPE:               v = g; valid = true; break;
            case PM_SLOPE_DEGREES:       v = atan(g) * (180.0 / k_pi_gpu); valid = true; break;
            case PM_ASPECT_RADIANS:
                if (g2 > flat_tol_sq)
                {
                    v = atan2(-p, -q);
                    if (v < 0.0) v += 2.0 * k_pi_gpu;
                    valid = true;
                }
                break;
            case PM_ASPECT_DEGREES:
                if (g2 > flat_tol_sq)
                {
                    v = atan2(-p, -q);
                    if (v < 0.0) v += 2.0 * k_pi_gpu;
                    v *= (180.0 / k_pi_gpu);
                    valid = true;
                }
                break;
            case PM_LAPLACIAN:           v = r + t; valid = true; break;
            case PM_HESSIAN_DETERMINANT: v = r * t - s * s; valid = true; break;
            case PM_MEAN_CURVATURE:      v = H; valid = isfinite(v); break;
            case PM_GAUSSIAN_CURVATURE:  v = K; valid = isfinite(v); break;
            case PM_MAXIMUM_CURVATURE:   v = kmax; valid = isfinite(v); break;
            case PM_MINIMUM_CURVATURE:   v = kmin; valid = isfinite(v); break;
            case PM_CURVEDNESS:
                if (isfinite(kmax) && isfinite(kmin))
                {
                    v = sqrt(0.5 * (kmax * kmax + kmin * kmin));
                    valid = isfinite(v);
                }
                break;
            case PM_SHAPE_INDEX:
                if (isfinite(kmax) && isfinite(kmin))
                {
                    const double denom = kmin - kmax;
                    if (fabs(denom) > 1e-15)
                    {
                        v = (2.0 / k_pi_gpu) * atan((kmin + kmax) / denom);
                        valid = isfinite(v);
                    }
                }
                break;
            case PM_PROFILE_CURVATURE:           v = vert_curv; valid = isfinite(v); break;
            case PM_PLAN_CURVATURE:              v = plan_curv; valid = isfinite(v); break;
            case PM_TANGENTIAL_CURVATURE:        v = tan_curv; valid = isfinite(v); break;
            case PM_DIFFERENCE_CURVATURE:        v = diff_curv; valid = isfinite(v); break;
            case PM_HORIZONTAL_EXCESS_CURVATURE: v = h_excess; valid = isfinite(v); break;
            case PM_VERTICAL_EXCESS_CURVATURE:   v = v_excess; valid = isfinite(v); break;
            case PM_ACCUMULATION_CURVATURE:      v = accum_curv; valid = isfinite(v); break;
            case PM_RING_CURVATURE:              v = ring_curv; valid = isfinite(v); break;
            case PM_ROTOR:                       v = rotor; valid = isfinite(v); break;
            default: break;
        }

        if (valid)
        {
            out_data[o] = static_cast<float>(v);
            out_mask[o] = 1;
        }
    }
}

} // anonymous namespace

bool gpu_available()
{
    int n = 0;
    cudaError_t err = cudaGetDeviceCount(&n);
    return (err == cudaSuccess && n > 0);
}

namespace {

bool build_polyfit_constants(
    int                 nrows,
    int                 ncols,
    double              cell_size,
    double              nodata_value,
    int                 W,
    polyfit_constants&  k)
{
    if (nrows <= 0 || ncols <= 0 || W < 1)
        return false;

    const double Wd  = static_cast<double>(W);
    const double W3d = static_cast<double>(2 * W + 1);
    const double sk2 = Wd * (Wd + 1.0) * W3d / 3.0;
    const double sk4 = 2.0 * Wd * (Wd + 1.0) * W3d *
                       (3.0 * Wd * Wd + 3.0 * Wd - 1.0) / 30.0;

    const double N      = W3d * W3d;
    const double d2     = cell_size * cell_size;
    const double d4     = d2 * d2;
    const double Sx2    = W3d * d2 * sk2;   // == Sy2
    const double Sx4    = W3d * d4 * sk4;   // == Sy4
    const double Sx2y2  = d4 * sk2 * sk2;

    if (!(Sx2y2 > 1e-40) || !(Sx2 > 1e-40))
        return false;

    const double M3[9] = {
        Sx4,   Sx2y2, Sx2,
        Sx2y2, Sx4,   Sx2,
        Sx2,   Sx2,   N
    };
    if (!invert3x3_rowmajor(M3, k.invM3))
        return false;

    k.nrows     = nrows;
    k.ncols     = ncols;
    k.W         = W;
    k.cell_size = static_cast<float>(cell_size);
    k.nodata    = static_cast<float>(nodata_value);
    k.inv_Sx2   = 1.0 / Sx2;
    k.inv_Sx2y2 = 1.0 / Sx2y2;
    return true;
}

} // anonymous namespace

// Retained as a fallback/reference implementation.
bool direct_polyfit_gpu_bruteforce(
    const float*        z,
    const std::uint8_t* mask,
    int                 nrows,
    int                 ncols,
    double              cell_size,
    double              nodata_value,
    int                 W,
    float*              out_a,
    float*              out_b,
    float*              out_c,
    float*              out_d,
    float*              out_e,
    float*              out_f,
    std::uint8_t*       out_mask)
{
    if (!gpu_available())
        return false;

    polyfit_constants k;
    if (!build_polyfit_constants(nrows, ncols, cell_size, nodata_value, W, k))
        return false;

    const std::size_t npix  = static_cast<std::size_t>(nrows) * ncols;
    const std::size_t fbytes = npix * sizeof(float);
    const std::size_t mbytes = npix * sizeof(unsigned char);

    // -------------------------------------------------------------------------
    // Device allocation
    // -------------------------------------------------------------------------
    float*         d_z = nullptr;
    unsigned char* d_mask = nullptr;
    float*         d_a = nullptr; float* d_b = nullptr; float* d_c = nullptr;
    float*         d_d = nullptr; float* d_e = nullptr; float* d_f = nullptr;
    unsigned char* d_omask = nullptr;

    cudaError_t err = cudaSuccess;
    bool alloc_ok =
        (cudaMalloc(&d_z,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_mask,  mbytes) == cudaSuccess) &&
        (cudaMalloc(&d_a,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_b,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_c,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_d,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_e,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_f,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_omask, mbytes) == cudaSuccess);

    if (!alloc_ok)
    {
        free_all(d_z, d_mask, d_a, d_b, d_c, d_d, d_e, d_f);
        if (d_omask) cudaFree(d_omask);
        return false;
    }

    // Copy inputs to device.
    if (cudaMemcpy(d_z,    z,    fbytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_mask, mask, mbytes, cudaMemcpyHostToDevice) != cudaSuccess)
    {
        free_all(d_z, d_mask, d_a, d_b, d_c, d_d, d_e, d_f);
        cudaFree(d_omask);
        return false;
    }

    // -------------------------------------------------------------------------
    // Launch: 1D grid over all pixels (boundary pixels short-circuit).
    // -------------------------------------------------------------------------
    const int threads = 256;
    const long blocks = (static_cast<long>(npix) + threads - 1) / threads;

    direct_polyfit_kernel_bruteforce<<<static_cast<unsigned int>(blocks), threads>>>(
        d_z, d_mask, k, d_a, d_b, d_c, d_d, d_e, d_f, d_omask);

    err = cudaGetLastError();
    if (err == cudaSuccess)
        err = cudaDeviceSynchronize();

    if (err != cudaSuccess)
    {
        free_all(d_z, d_mask, d_a, d_b, d_c, d_d, d_e, d_f);
        cudaFree(d_omask);
        return false;
    }

    // Copy results back.
    bool copy_ok =
        (cudaMemcpy(out_a,    d_a,     fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_b,    d_b,     fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_c,    d_c,     fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_d,    d_d,     fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_e,    d_e,     fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_f,    d_f,     fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_mask, d_omask, mbytes, cudaMemcpyDeviceToHost) == cudaSuccess);

    free_all(d_z, d_mask, d_a, d_b, d_c, d_d, d_e, d_f);
    cudaFree(d_omask);

    return copy_ok;
}

struct polyfit_device_bundle
{
    int nrows = 0;
    int ncols = 0;
    float nodata = 0.0f;
    float* d_a = nullptr;
    float* d_b = nullptr;
    float* d_c = nullptr;
    float* d_d = nullptr;
    float* d_e = nullptr;
    float* d_f = nullptr;
    unsigned char* d_mask = nullptr;
    float* d_metric_out_data = nullptr;
    unsigned char* d_metric_out_mask = nullptr;
    int* d_metric_id = nullptr;
};

namespace {

void free_polyfit_device_bundle_arrays(polyfit_device_bundle* bundle)
{
    if (!bundle) return;
    if (bundle->d_a) cudaFree(bundle->d_a);
    if (bundle->d_b) cudaFree(bundle->d_b);
    if (bundle->d_c) cudaFree(bundle->d_c);
    if (bundle->d_d) cudaFree(bundle->d_d);
    if (bundle->d_e) cudaFree(bundle->d_e);
    if (bundle->d_f) cudaFree(bundle->d_f);
    if (bundle->d_mask) cudaFree(bundle->d_mask);
    if (bundle->d_metric_out_data) cudaFree(bundle->d_metric_out_data);
    if (bundle->d_metric_out_mask) cudaFree(bundle->d_metric_out_mask);
    if (bundle->d_metric_id) cudaFree(bundle->d_metric_id);
    bundle->d_a = nullptr;
    bundle->d_b = nullptr;
    bundle->d_c = nullptr;
    bundle->d_d = nullptr;
    bundle->d_e = nullptr;
    bundle->d_f = nullptr;
    bundle->d_mask = nullptr;
    bundle->d_metric_out_data = nullptr;
    bundle->d_metric_out_mask = nullptr;
    bundle->d_metric_id = nullptr;
}

} // anonymous namespace

bool direct_polyfit_gpu_bundle_create(
    const float*        z,
    const std::uint8_t* mask,
    int                 nrows,
    int                 ncols,
    double              cell_size,
    double              nodata_value,
    int                 W,
    polyfit_device_bundle** out_bundle)
{
    if (!gpu_available())
        return false;
    if (out_bundle == nullptr)
        return false;
    *out_bundle = nullptr;

    polyfit_constants k;
    if (!build_polyfit_constants(nrows, ncols, cell_size, nodata_value, W, k))
        return false;

    const std::size_t npix   = static_cast<std::size_t>(nrows) * ncols;
    const std::size_t fbytes = npix * sizeof(float);
    const std::size_t mbytes = npix * sizeof(unsigned char);
    const std::size_t dbytes = npix * sizeof(double);
    const std::size_t ibytes = npix * sizeof(int);

    float*         d_z = nullptr;
    unsigned char* d_mask = nullptr;
    polyfit_device_bundle* bundle = new polyfit_device_bundle();
    bundle->nrows = nrows;
    bundle->ncols = ncols;
    bundle->nodata = static_cast<float>(nodata_value);

    double*        d_H0 = nullptr; double* d_H1 = nullptr; double* d_H2 = nullptr;
    int*           d_Hn = nullptr;

    auto cleanup = [&]() {
        if (d_z)     cudaFree(d_z);
        if (d_mask)  cudaFree(d_mask);
        if (d_H0)    cudaFree(d_H0);
        if (d_H1)    cudaFree(d_H1);
        if (d_H2)    cudaFree(d_H2);
        if (d_Hn)    cudaFree(d_Hn);
        free_polyfit_device_bundle_arrays(bundle);
        delete bundle;
    };

    bool alloc_ok =
        (cudaMalloc(&d_z,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_mask,  mbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_a,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_b,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_c,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_d,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_e,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_f,     fbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_mask, mbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_metric_out_data, fbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_metric_out_mask, mbytes) == cudaSuccess) &&
        (cudaMalloc(&bundle->d_metric_id, sizeof(int)) == cudaSuccess) &&
        (cudaMalloc(&d_H0,    dbytes) == cudaSuccess) &&
        (cudaMalloc(&d_H1,    dbytes) == cudaSuccess) &&
        (cudaMalloc(&d_H2,    dbytes) == cudaSuccess) &&
        (cudaMalloc(&d_Hn,    ibytes) == cudaSuccess);

    if (!alloc_ok)
    {
        cleanup();
        return false;
    }

    if (cudaMemcpy(d_z, z, fbytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_mask, mask, mbytes, cudaMemcpyHostToDevice) != cudaSuccess)
    {
        cleanup();
        return false;
    }

    const int h_threads = 128;
    const unsigned int h_blocks = static_cast<unsigned int>((nrows + h_threads - 1) / h_threads);
    horizontal_moments_kernel<<<h_blocks, h_threads>>>(d_z, d_mask, k, d_H0, d_H1, d_H2, d_Hn);

    cudaError_t err = cudaGetLastError();
    if (err == cudaSuccess)
        err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        cleanup();
        return false;
    }

    const int v_threads = 256;
    const unsigned int v_blocks = static_cast<unsigned int>((static_cast<long>(npix) + v_threads - 1) / v_threads);
    direct_polyfit_kernel_rolling<<<v_blocks, v_threads>>>(
        d_H0, d_H1, d_H2, d_Hn, k,
        bundle->d_a, bundle->d_b, bundle->d_c,
        bundle->d_d, bundle->d_e, bundle->d_f, bundle->d_mask);

    err = cudaGetLastError();
    if (err == cudaSuccess)
        err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        cleanup();
        return false;
    }

    if (d_z) cudaFree(d_z);
    if (d_mask) cudaFree(d_mask);
    if (d_H0) cudaFree(d_H0);
    if (d_H1) cudaFree(d_H1);
    if (d_H2) cudaFree(d_H2);
    if (d_Hn) cudaFree(d_Hn);

    *out_bundle = bundle;
    return true;
}

bool direct_polyfit_gpu_bundle_download(
    const polyfit_device_bundle* bundle,
    float*                       out_a,
    float*                       out_b,
    float*                       out_c,
    float*                       out_d,
    float*                       out_e,
    float*                       out_f,
    std::uint8_t*                out_mask)
{
    if (!bundle || !out_a || !out_b || !out_c || !out_d || !out_e || !out_f || !out_mask)
        return false;

    const std::size_t npix = static_cast<std::size_t>(bundle->nrows) * bundle->ncols;
    const std::size_t fbytes = npix * sizeof(float);
    const std::size_t mbytes = npix * sizeof(unsigned char);

    return
        (cudaMemcpy(out_a,    bundle->d_a,    fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_b,    bundle->d_b,    fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_c,    bundle->d_c,    fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_d,    bundle->d_d,    fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_e,    bundle->d_e,    fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_f,    bundle->d_f,    fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_mask, bundle->d_mask, mbytes, cudaMemcpyDeviceToHost) == cudaSuccess);
}

void destroy_polyfit_device_bundle(polyfit_device_bundle* bundle)
{
    if (!bundle) return;
    free_polyfit_device_bundle_arrays(bundle);
    delete bundle;
}

bool direct_polyfit_gpu_rolling(
    const float*        z,
    const std::uint8_t* mask,
    int                 nrows,
    int                 ncols,
    double              cell_size,
    double              nodata_value,
    int                 W,
    float*              out_a,
    float*              out_b,
    float*              out_c,
    float*              out_d,
    float*              out_e,
    float*              out_f,
    std::uint8_t*       out_mask)
{
    polyfit_device_bundle* bundle = nullptr;
    if (!direct_polyfit_gpu_bundle_create(
            z, mask, nrows, ncols, cell_size, nodata_value, W, &bundle))
    {
        return false;
    }

    const bool ok = direct_polyfit_gpu_bundle_download(
        bundle, out_a, out_b, out_c, out_d, out_e, out_f, out_mask);
    destroy_polyfit_device_bundle(bundle);
    return ok;
}

bool direct_polyfit_gpu(
    const float*        z,
    const std::uint8_t* mask,
    int                 nrows,
    int                 ncols,
    double              cell_size,
    double              nodata_value,
    int                 W,
    float*              out_a,
    float*              out_b,
    float*              out_c,
    float*              out_d,
    float*              out_e,
    float*              out_f,
    std::uint8_t*       out_mask)
{
    if (direct_polyfit_gpu_rolling(
            z, mask, nrows, ncols, cell_size, nodata_value, W,
            out_a, out_b, out_c, out_d, out_e, out_f, out_mask))
    {
        return true;
    }

    return direct_polyfit_gpu_bruteforce(
        z, mask, nrows, ncols, cell_size, nodata_value, W,
        out_a, out_b, out_c, out_d, out_e, out_f, out_mask);
}

bool compute_polyfit_metrics_gpu(
    const float*        a_data,
    const float*        b_data,
    const float*        c_data,
    const float*        d_data,
    const float*        e_data,
    const float*        f_data,
    const std::uint8_t* a_mask,
    const std::uint8_t* b_mask,
    const std::uint8_t* c_mask,
    const std::uint8_t* d_mask,
    const std::uint8_t* e_mask,
    const std::uint8_t* f_mask,
    int                 nrows,
    int                 ncols,
    double              nodata_value,
    double              z_factor,
    double              flat_slope_tolerance,
    const int*          metric_ids,
    int                 n_metrics,
    float*              out_data,
    std::uint8_t*       out_mask)
{
    if (!gpu_available())
        return false;
    if (nrows <= 0 || ncols <= 0 || n_metrics <= 0 || metric_ids == nullptr)
        return false;

    const std::size_t npix = static_cast<std::size_t>(nrows) * ncols;
    const std::size_t fbytes = npix * sizeof(float);
    const std::size_t mbytes = npix * sizeof(unsigned char);
    const std::size_t out_fbytes = static_cast<std::size_t>(n_metrics) * fbytes;
    const std::size_t out_mbytes = static_cast<std::size_t>(n_metrics) * mbytes;
    const std::size_t id_bytes = static_cast<std::size_t>(n_metrics) * sizeof(int);

    float *d_a = nullptr, *d_b = nullptr, *d_c = nullptr;
    float *d_d = nullptr, *d_e = nullptr, *d_f = nullptr;
    unsigned char *d_am = nullptr, *d_bm = nullptr, *d_cm = nullptr;
    unsigned char *d_dm = nullptr, *d_em = nullptr, *d_fm = nullptr;
    int* d_ids = nullptr;
    float* d_out_data = nullptr;
    unsigned char* d_out_mask = nullptr;

    auto cleanup = [&]() {
        if (d_a) cudaFree(d_a);
        if (d_b) cudaFree(d_b);
        if (d_c) cudaFree(d_c);
        if (d_d) cudaFree(d_d);
        if (d_e) cudaFree(d_e);
        if (d_f) cudaFree(d_f);
        if (d_am) cudaFree(d_am);
        if (d_bm) cudaFree(d_bm);
        if (d_cm) cudaFree(d_cm);
        if (d_dm) cudaFree(d_dm);
        if (d_em) cudaFree(d_em);
        if (d_fm) cudaFree(d_fm);
        if (d_ids) cudaFree(d_ids);
        if (d_out_data) cudaFree(d_out_data);
        if (d_out_mask) cudaFree(d_out_mask);
    };

    bool alloc_ok =
        (cudaMalloc(&d_a, fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_b, fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_c, fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_d, fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_e, fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_f, fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_am, mbytes) == cudaSuccess) &&
        (cudaMalloc(&d_bm, mbytes) == cudaSuccess) &&
        (cudaMalloc(&d_cm, mbytes) == cudaSuccess) &&
        (cudaMalloc(&d_dm, mbytes) == cudaSuccess) &&
        (cudaMalloc(&d_em, mbytes) == cudaSuccess) &&
        (cudaMalloc(&d_fm, mbytes) == cudaSuccess) &&
        (cudaMalloc(&d_ids, id_bytes) == cudaSuccess) &&
        (cudaMalloc(&d_out_data, out_fbytes) == cudaSuccess) &&
        (cudaMalloc(&d_out_mask, out_mbytes) == cudaSuccess);

    if (!alloc_ok)
    {
        cleanup();
        return false;
    }

    bool copy_in_ok =
        (cudaMemcpy(d_a, a_data, fbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_b, b_data, fbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_c, c_data, fbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_d, d_data, fbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_e, e_data, fbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_f, f_data, fbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_am, a_mask, mbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_bm, b_mask, mbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_cm, c_mask, mbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_dm, d_mask, mbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_em, e_mask, mbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_fm, f_mask, mbytes, cudaMemcpyHostToDevice) == cudaSuccess) &&
        (cudaMemcpy(d_ids, metric_ids, id_bytes, cudaMemcpyHostToDevice) == cudaSuccess);

    if (!copy_in_ok)
    {
        cleanup();
        return false;
    }

    const int threads = 256;
    const unsigned int blocks = static_cast<unsigned int>((static_cast<long>(npix) + threads - 1) / threads);
    compute_polyfit_metrics_kernel<<<blocks, threads>>>(
        d_a, d_b, d_c, d_d, d_e, d_f,
        d_am, d_bm, d_cm, d_dm, d_em, d_fm,
        nrows, ncols,
        static_cast<float>(nodata_value), z_factor, flat_slope_tolerance,
        d_ids, n_metrics,
        d_out_data, d_out_mask);

    cudaError_t err = cudaGetLastError();
    if (err == cudaSuccess)
        err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        cleanup();
        return false;
    }

    const bool copy_out_ok =
        (cudaMemcpy(out_data, d_out_data, out_fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_mask, d_out_mask, out_mbytes, cudaMemcpyDeviceToHost) == cudaSuccess);

    cleanup();
    return copy_out_ok;
}

bool compute_polyfit_metric_from_bundle_gpu(
    const polyfit_device_bundle* bundle,
    double                       z_factor,
    double                       flat_slope_tolerance,
    int                          metric_id,
    float*                       out_data,
    std::uint8_t*                out_mask)
{
    if (!bundle || !out_data || !out_mask)
        return false;
    if (!bundle->d_metric_out_data || !bundle->d_metric_out_mask || !bundle->d_metric_id)
        return false;

    const std::size_t npix = static_cast<std::size_t>(bundle->nrows) * bundle->ncols;
    const std::size_t fbytes = npix * sizeof(float);
    const std::size_t mbytes = npix * sizeof(unsigned char);

    if (cudaMemcpy(bundle->d_metric_id, &metric_id, sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess)
        return false;

    const int threads = 256;
    const unsigned int blocks = static_cast<unsigned int>((static_cast<long>(npix) + threads - 1) / threads);
    compute_polyfit_metrics_kernel<<<blocks, threads>>>(
        bundle->d_a, bundle->d_b, bundle->d_c,
        bundle->d_d, bundle->d_e, bundle->d_f,
        bundle->d_mask, bundle->d_mask, bundle->d_mask,
        bundle->d_mask, bundle->d_mask, bundle->d_mask,
        bundle->nrows, bundle->ncols,
        bundle->nodata, z_factor, flat_slope_tolerance,
        bundle->d_metric_id, 1,
        bundle->d_metric_out_data, bundle->d_metric_out_mask);

    cudaError_t err = cudaGetLastError();
    if (err == cudaSuccess)
        err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
        return false;

    const bool ok =
        (cudaMemcpy(out_data, bundle->d_metric_out_data, fbytes, cudaMemcpyDeviceToHost) == cudaSuccess) &&
        (cudaMemcpy(out_mask, bundle->d_metric_out_mask, mbytes, cudaMemcpyDeviceToHost) == cudaSuccess);

    return ok;
}

} // namespace gpu
} // namespace lsdtt
