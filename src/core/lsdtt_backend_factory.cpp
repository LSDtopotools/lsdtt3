// lsdtt_backend_factory.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_backend_factory.hpp"
#include "../cpu/lsdtt_raster_cpu.hpp"

// GPU backend is only included when compiled with -DLSDTT_GPU=ON
#ifdef LSDTT_GPU_ENABLED
#include "../gpu/lsdtt_raster_gpu.hpp"
#include <cuda_runtime.h>
#endif

#include <stdexcept>

namespace lsdtt {

// -----------------------------------------------------------------------------
// gpu_available
// -----------------------------------------------------------------------------

bool gpu_available()
{
#ifdef LSDTT_GPU_ENABLED
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    return (err == cudaSuccess && device_count > 0);
#else
    return false;
#endif
}

// -----------------------------------------------------------------------------
// Global compute-backend preference
//
// Pure host-side state. No CUDA here: whether a GPU is actually used is decided
// at each dispatch site (in lsdtt3_cpu) by combining this preference with the
// runtime device check lsdtt::gpu::gpu_available() from the GPU kernel library.
// -----------------------------------------------------------------------------
namespace {
backend_type g_compute_backend = backend_type::auto_detect;
}

void set_compute_backend(backend_type backend)
{
    g_compute_backend = backend;
}

backend_type get_compute_backend()
{
    return g_compute_backend;
}

// -----------------------------------------------------------------------------
// backend_name
// -----------------------------------------------------------------------------

std::string backend_name(backend_type backend)
{
    switch (backend)
    {
        case backend_type::cpu:         return "cpu";
        case backend_type::gpu:         return "gpu";
        case backend_type::auto_detect: return "auto_detect";
        default:                        return "unknown";
    }
}

// -----------------------------------------------------------------------------
// make_raster_engine
// -----------------------------------------------------------------------------

std::unique_ptr<lsdtt_raster_base> make_raster_engine(
    const lsdtt_data_loader& loader,
    backend_type             backend
)
{
    if (!loader.is_open())
    {
        throw std::runtime_error(
            "lsdtt::make_raster_engine: loader has no raster open"
        );
    }

    // Resolve auto_detect to a concrete backend
    if (backend == backend_type::auto_detect)
    {
        backend = gpu_available() ? backend_type::gpu : backend_type::cpu;
    }

    // Construct the chosen backend
    if (backend == backend_type::cpu)
    {
        return std::make_unique<lsdtt_raster_cpu>(loader);
    }

    if (backend == backend_type::gpu)
    {
#ifdef LSDTT_GPU_ENABLED
        if (!gpu_available())
        {
            throw std::runtime_error(
                "lsdtt::make_raster_engine: GPU backend requested but no "
                "usable GPU was detected at runtime"
            );
        }
        return std::make_unique<lsdtt_raster_gpu>(loader);
#else
        throw std::runtime_error(
            "lsdtt::make_raster_engine: GPU backend requested but lsdtt3 "
            "was not compiled with GPU support. "
            "Rebuild with -DLSDTT_GPU=ON."
        );
#endif
    }

    // Should never reach here
    throw std::runtime_error(
        "lsdtt::make_raster_engine: unknown backend_type"
    );
}

} // namespace lsdtt
