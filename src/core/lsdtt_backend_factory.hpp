// lsdtt_backend_factory.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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
//   loader  — an lsdtt_data_loader with an open raster
//   backend — which backend to use (default: auto_detect)
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
// backend_name
//
// Returns a human-readable string for a backend_type value.
// Useful for logging and CLI output.
// -----------------------------------------------------------------------------
std::string backend_name(backend_type backend);

} // namespace lsdtt
