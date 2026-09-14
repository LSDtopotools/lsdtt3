// lsdtt_backend_factory.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_backend_factory
//
// Factory function for creating raster compute engines.
// Returns an lsdtt_raster_cpu or lsdtt_raster_gpu wrapped behind the
// lsdtt_raster_base interface. Calling code never needs to know which
// backend it received.
//
// Usage:
//
//   lsdtt_data_loader loader;
//   loader.open_raster("dem.tif");
//
//   // Let the factory decide (GPU if available, CPU otherwise)
//   auto engine = lsdtt::make_raster_engine(loader);
//
//   // Force CPU
//   auto engine = lsdtt::make_raster_engine(loader, lsdtt::backend_type::cpu);
//
//   // Force GPU (throws if no GPU present)
//   auto engine = lsdtt::make_raster_engine(loader, lsdtt::backend_type::gpu);
//
// GPU support is only compiled in when the build is configured with
// -DLSDTT_GPU=ON. If GPU is requested at runtime but support was not compiled
// in, make_raster_engine() throws std::runtime_error.
// -----------------------------------------------------------------------------

#pragma once

#include <memory>
#include <string>

#include "lsdtt_raster_base.hpp"
#include "lsdtt_data_loader.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// backend_type
// -----------------------------------------------------------------------------
enum class backend_type
{
    cpu,         // always use CPU backend
    gpu,         // always use GPU backend (throws if unavailable)
    auto_detect  // use GPU if available, fall back to CPU silently
};

// -----------------------------------------------------------------------------
// make_raster_engine
//
// Creates and returns a raster compute engine loaded with data from loader.
//
// Parameters:
//   loader  - an lsdtt_data_loader with an open raster
//   backend - which backend to use (default: auto_detect)
//
// Returns:
//   std::unique_ptr<lsdtt_raster_base> owning the engine
//
// Throws:
//   std::runtime_error if backend == gpu and no GPU is available or GPU
//   support was not compiled in
//   std::runtime_error if loader has no raster open
// -----------------------------------------------------------------------------
std::unique_ptr<lsdtt_raster_base> make_raster_engine(
    const lsdtt_data_loader& loader,
    backend_type backend = backend_type::auto_detect
);

// -----------------------------------------------------------------------------
// gpu_available
//
// Returns true if at least one usable GPU is detected at runtime.
// Always returns false if the build was configured without -DLSDTT_GPU=ON.
// -----------------------------------------------------------------------------
bool gpu_available();

// -----------------------------------------------------------------------------
// Global compute-backend preference (selective-dispatch compute kernels)
//
// This is a *preference*, not a hard mode switch. It is consulted by the small
// set of heavy per-pixel compute functions that have a GPU implementation
// (currently direct_polyfit, with more roughness/polyfit kernels to follow).
// Each such function dispatches as:
//
//     if (preference != cpu) AND (a GPU kernel exists for this function)
//        AND (a device is present)  ->  run the GPU kernel
//     else                          ->  run the CPU implementation
//
// Functions with no GPU implementation ignore the preference entirely and run
// on the CPU. File I/O (lsdtt_data_loader) is never affected. A typical "GPU"
// run is therefore a mix: CPU does all I/O and orchestration, GPU runs the
// ported kernels, and anything unported stays on the CPU automatically.
//
// This toggle is independent of the GPU *raster-engine backend* selected by
// make_raster_engine(): that is a separate (future) axis. set_compute_backend()
// has no effect when lsdtt3 was built without GPU kernel support.
//
// Note: the preference is process-wide global state. Set it once at start-up
// (e.g. from a CLI parameter); do not flip it from multiple threads while a
// dispatched compute function is running.
//
// Default: backend_type::auto_detect (use the GPU when kernels and a device
// are both available, otherwise CPU).
// -----------------------------------------------------------------------------
void         set_compute_backend(backend_type backend);
backend_type get_compute_backend();

// -----------------------------------------------------------------------------
// backend_name
//
// Returns a human-readable string for a backend_type value.
// Useful for logging and CLI output.
// -----------------------------------------------------------------------------
std::string backend_name(backend_type backend);

} // namespace lsdtt
