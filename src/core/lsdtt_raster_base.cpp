// lsdtt_raster_base.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// CPU fallback implementations for lsdtt_raster_base virtual methods.
//
// These are placeholder implementations that will be replaced as the CPU
// algorithms are developed. The GPU subclass inherits these and silently
// uses them until proper GPU kernels are written.
// -----------------------------------------------------------------------------

#include "lsdtt_raster_base.hpp"

#include <stdexcept>

namespace lsdtt {

// -----------------------------------------------------------------------------
// CPU fallback — flow routing
// These will be implemented properly in a later development phase.
// -----------------------------------------------------------------------------

void lsdtt_raster_base::compute_flow_direction()
{
    cpu_compute_flow_direction();
}

void lsdtt_raster_base::compute_flow_accumulation()
{
    cpu_compute_flow_accumulation();
}

// -----------------------------------------------------------------------------
// CPU fallback — spectral
// -----------------------------------------------------------------------------

void lsdtt_raster_base::compute_fft()
{
    cpu_compute_fft();
}

void lsdtt_raster_base::apply_wavelet_filter()
{
    cpu_apply_wavelet_filter();
}

// -----------------------------------------------------------------------------
// Protected CPU implementations
// Called by the fallbacks above; also callable directly from lsdtt_raster_cpu.
// Implementations will be filled in when each algorithm is developed.
// -----------------------------------------------------------------------------

void lsdtt_raster_base::cpu_compute_flow_direction()
{
    throw std::runtime_error(
        "lsdtt_raster_base: cpu_compute_flow_direction() not yet implemented"
    );
}

void lsdtt_raster_base::cpu_compute_flow_accumulation()
{
    throw std::runtime_error(
        "lsdtt_raster_base: cpu_compute_flow_accumulation() not yet implemented"
    );
}

void lsdtt_raster_base::cpu_compute_fft()
{
    throw std::runtime_error(
        "lsdtt_raster_base: cpu_compute_fft() not yet implemented"
    );
}

void lsdtt_raster_base::cpu_apply_wavelet_filter()
{
    throw std::runtime_error(
        "lsdtt_raster_base: cpu_apply_wavelet_filter() not yet implemented"
    );
}

} // namespace lsdtt
