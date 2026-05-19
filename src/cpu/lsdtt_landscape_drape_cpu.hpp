// lsdtt_landscape_drape_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Landscape drape — computes an upper envelope surface that drapes over the
// topography by resting on the highest points (ridges), and derives a raster
// of incision depth (envelope − DEM).
//
// Algorithm: cloth simulation (Verlet integration)
//
//   A virtual cloth is initialised as a flat horizontal sheet at the maximum
//   DEM elevation.  Under gravity the cloth falls and comes to rest on the
//   ridges and summits of the terrain; spring forces between neighbouring cloth
//   particles keep the surface smooth between contact points.
//
//   This is mathematically equivalent to the Zhang et al. (2016) Cloth
//   Simulation Filter (CSF) used for LiDAR ground classification, applied
//   directly (without DEM inversion) to identify the upper landscape envelope.
//
//   At equilibrium, the cloth satisfies the discrete Laplace equation
//   (∇²f = 0) between its terrain contact points — identical to harmonic
//   (bilinear spline) interpolation through the ridge points.  The explicit
//   Verlet time-stepping naturally enforces the above-terrain constraint that
//   a direct solver cannot.
//
// Key parameters:
//   rigidity        — spring stiffness; controls how smoothly the cloth bridges
//                     between ridges.  Higher = smoother, longer-wavelength
//                     envelope.  Analogous to the Tension parameter in the
//                     Carretier & Valenti (GRASS) EnvelopeSurface script and to
//                     spring stiffness in CSF.
//   class_threshold — contact tolerance (m).  A cloth particle is pinned when
//                     it comes within this distance of the terrain surface.
//                     Increase for noisy or pitted DEMs to avoid premature
//                     pinning in artefact pits.
//
// Reference:
//   Carretier S. & Valenti C. — EnvelopSurface.bash (IRD/GET, France)
//   Zhang W. et al. (2016) An Easy-to-Use Airborne LiDAR Data Filtering Method
//     Based on Cloth Simulation. Remote Sensing 8(6):501.
//     DOI: 10.3390/rs8060501
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// lsdtt_cloth_drape_options
// -----------------------------------------------------------------------------
struct lsdtt_cloth_drape_options
{
    float gravity           = 0.3f;    // downward force per step (m/step²)
    float rigidity          = 1.0f;    // spring stiffness — higher = smoother envelope
    float time_step         = 0.65f;   // Verlet dt
    float class_threshold   = 0.5f;    // contact tolerance (m)
    int   max_iterations    = 500;
    float convergence_tol   = 0.005f;  // max displacement (m) to declare settled
    int   n_threads         = 0;       // 0 = OMP default / hardware maximum
};

// -----------------------------------------------------------------------------
// lsdtt_cloth_drape_result
// -----------------------------------------------------------------------------
struct lsdtt_cloth_drape_result
{
    lsdtt_raster_cpu envelope;              // cloth resting position = upper landscape envelope
    lsdtt_raster_cpu incision;             // envelope − DEM  (negative where DEM > envelope)
    lsdtt_raster_cpu crest_mask;           // 1 = pinned contact pixel (ridge / summit)
    int              iterations_taken      = 0;
    float            final_max_displacement = 0.0f;
};

// -----------------------------------------------------------------------------
// compute_cloth_drape
//
// Runs the Verlet cloth simulation on dem and returns the envelope, incision,
// and crest mask rasters together with convergence diagnostics.
//
// The DEM must already be in a projected CRS with units of metres.
// Preprocessing (pit fill / breach) is the caller's responsibility.
// -----------------------------------------------------------------------------
lsdtt_cloth_drape_result compute_cloth_drape(
    const lsdtt_raster_cpu&          dem,
    const lsdtt_cloth_drape_options& opts = {}
);

} // namespace lsdtt
