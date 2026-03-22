// lsdtt_raster_base.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_base
//
// Abstract base class for all raster compute objects in lsdtt3.
// Defines the dispatch interface that lsdtt_raster_cpu and lsdtt_raster_gpu
// must satisfy.
//
// Design:
//   This class is a pure dispatch interface. It does NOT expose data accessors
//   (get_data / get_mask) because those are CPU-specific concepts — GPU memory
//   lives in device memory, not in std::vector. Analysis free functions in
//   src/cpu/ take lsdtt_raster_cpu directly and access data through its own
//   member functions.
//
//   get_meta() is retained here because raster metadata (dimensions, CRS,
//   cell size, nodata value) is always host-side regardless of backend.
//
// Virtual method categories:
//
//   Pure virtual  — must be implemented by every backend. Operations where
//                   GPU coverage is guaranteed from the start, or where the
//                   class cannot function without them.
//
//   CPU fallback  — base provides a working CPU implementation. The GPU
//                   subclass inherits and uses it until a GPU kernel is
//                   written, at which point the subclass overrides.
//
// Reference:
//   Zevenbergen, L.W. and Thorne, C.R. (1987) Quantitative analysis of land
//   surface topography. Earth Surface Processes and Landforms, 12, 47-56.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>

#include "lsdtt_raster_meta.hpp"

// Forward declaration
namespace lsdtt
{
    class lsdtt_data_loader;
}

namespace lsdtt {

class lsdtt_raster_base
{
public:

    // -------------------------------------------------------------------------
    // Metadata accessor
    //
    // Metadata is always host-side regardless of backend. This is the only
    // data accessor on the base class. For data access (float*, vector<float>)
    // use lsdtt_raster_cpu directly — those are CPU-specific types.
    // -------------------------------------------------------------------------

    virtual const lsdtt_raster_meta& get_meta() const = 0;

    // -------------------------------------------------------------------------
    // Pure virtual — every backend must implement these
    // -------------------------------------------------------------------------

    // Data exchange with the loader
    virtual void write_to_loader(lsdtt_data_loader& loader)        = 0;
    virtual void read_from_loader(const lsdtt_data_loader& loader) = 0;

    // Core filtering — OpenCV operations where GPU coverage is straightforward
    virtual void smooth_gaussian(double sigma)                              = 0;
    virtual void smooth_bilateral(double sigma_space, double sigma_colour)  = 0;

    // -------------------------------------------------------------------------
    // CPU fallback — base provides working CPU implementation
    // GPU subclass overrides when GPU kernel is ready
    // -------------------------------------------------------------------------

    // Flow routing — complex GPU implementation, developed later
    virtual void compute_flow_direction();
    virtual void compute_flow_accumulation();

    // Spectral — GPU version added incrementally
    virtual void compute_fft();
    virtual void apply_wavelet_filter();

    // -------------------------------------------------------------------------
    // Destructor
    // -------------------------------------------------------------------------

    virtual ~lsdtt_raster_base() = default;

protected:

    // CPU implementations used by the fallback methods above.
    // Also callable from lsdtt_raster_cpu for direct use.
    void cpu_compute_flow_direction();
    void cpu_compute_flow_accumulation();
    void cpu_compute_fft();
    void cpu_apply_wavelet_filter();
};

} // namespace lsdtt
