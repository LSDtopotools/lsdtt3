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
__global__ void direct_polyfit_kernel(
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

} // anonymous namespace

bool gpu_available()
{
    int n = 0;
    cudaError_t err = cudaGetDeviceCount(&n);
    return (err == cudaSuccess && n > 0);
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
    if (!gpu_available())
        return false;
    if (nrows <= 0 || ncols <= 0 || W < 1)
        return false;

    // -------------------------------------------------------------------------
    // Precompute the constant normal-matrix quantities (same algebra as the CPU
    // direct_polyfit). The symmetric centred window decouples the 6x6 normal
    // matrix into three scalar equations (C, D, E) and one 3x3 system (a,b,F).
    // -------------------------------------------------------------------------
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
    polyfit_constants k;
    if (!invert3x3_rowmajor(M3, k.invM3))
        return false;

    k.nrows     = nrows;
    k.ncols     = ncols;
    k.W         = W;
    k.cell_size = static_cast<float>(cell_size);
    k.nodata    = static_cast<float>(nodata_value);
    k.inv_Sx2   = 1.0 / Sx2;
    k.inv_Sx2y2 = 1.0 / Sx2y2;

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

    direct_polyfit_kernel<<<static_cast<unsigned int>(blocks), threads>>>(
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

} // namespace gpu
} // namespace lsdtt
