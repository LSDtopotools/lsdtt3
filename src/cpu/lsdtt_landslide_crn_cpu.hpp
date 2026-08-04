// lsdtt_landslide_crn_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Landslides -> cosmogenic nuclide concentrations.
//
// This is the whole point of the landslide module: what does an imposed
// landslide flux do to the nuclide concentration a basin delivers? Niemi et al.
// (2005) and Yanites et al. (2009) both ask exactly this question, and
// reproducing them is the validation target.
//
// *** NO NEW CRN PHYSICS IS DEFINED HERE. ***
//
// The cosmogenic numeric core is locked. Every number below comes from existing,
// validated calls in lsdtt_crn_column:
//
//   1. INHERITANCE. A landslide that strips effective depth d exposes rock that
//      was sitting at depth d, carrying the concentration it had there. That is
//      lsdtt_crn_column::steady_state(d, eps, ...) under the background rate.
//
//   2. RE-ACCUMULATION. Since failure the fresh surface has been exposed and
//      eroding. That is lsdtt_crn_column::run_transient over an erosion history
//      containing an instantaneous REMOVAL event -- a facility the transient
//      integrator already has (rate = 0, removal = depth in metres), and which
//      is precisely a landslide.
//
// So the composition is: background steady state, one removal event at the
// failure time, background erosion since. No approximation of "decay towards a
// new steady state" is written here, because writing one would be new physics
// sitting next to a locked core.
//
// HISTORY ROW ORDER: row 0 is the MOST RECENT step and yr BP increases with
// index. A reversed record is rejected fatally by the integrator -- it used to
// return a confidently wrong number instead, which is why the check exists. The
// history built here is ordered correctly by construction.
//
// CACHING: every pixel within one scar shares a failure time AND an excavation
// depth, so the transient solve is done once per (depth, time) pair rather than
// once per pixel. On a typical event that is a few hundred solves instead of a
// few thousand.
// -----------------------------------------------------------------------------

#ifndef LSDTT_LANDSLIDE_CRN_CPU_HPP
#define LSDTT_LANDSLIDE_CRN_CPU_HPP

#include <string>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct LandslideCrnConfig
{
    // Background effective erosion rate, applied before and after failure.
    double background_eff_erate_g_cm2_yr = 0.001;

    // Site production scaling (Stone 2000; 1.0 at SLHL) and shielding. These go
    // to the library SEPARATELY -- scaling as a flat P0 multiply, shielding
    // through the Vermeesch redistribution -- because merging them is wrong and
    // has already cost this codebase an 11.4 % disagreement once.
    double production_scaling = 1.0;
    double total_shielding    = 1.0;

    std::string muon_scheme        = "BraucherBorchers";
    double      rock_density_kg_m3 = 2650.0;

    // Present time on the same clock as time_of_last_failure.
    double t_now_yr = 0.0;

    // Nuclide selector for the output raster: "Be10", "Al26" or "C14".
    std::string nuclide = "Be10";
};


// Surface concentration (atoms/g) at every valid pixel.
//
// `excavation_depth_m` and `time_of_last_failure` are the two rasters
// trigger_landslide_event writes. A pixel that has never failed
// (time_of_last_failure == its no-data value) gets the ordinary steady-state
// surface concentration, so a landscape with no landslides returns a uniform
// field and is the natural control.
//
// Throws std::invalid_argument if the rasters disagree in geometry or the
// nuclide is unknown.
lsdtt_raster_cpu compute_landslide_surface_concentration(
    const lsdtt_raster_cpu&   excavation_depth_m,
    const lsdtt_raster_cpu&   time_of_last_failure,
    const LandslideCrnConfig& config);


// The undisturbed surface concentration for this configuration -- the control
// value every unfailed pixel takes. Exposed so a caller can report the
// reduction a landslide field causes without recomputing it.
double landslide_crn_background_concentration(const LandslideCrnConfig& config);

}  // namespace lsdtt

#endif  // LSDTT_LANDSLIDE_CRN_CPU_HPP
