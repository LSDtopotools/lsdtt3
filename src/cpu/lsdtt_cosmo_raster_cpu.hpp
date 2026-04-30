#pragma once
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_crn_parameters.hpp"
#include <string>
#include <vector>

namespace lsdtt
{

// ─────────────────────────────────────────────────────────────────────────────
// lsdtt_cosmo_raster_cpu
//
// Raster-level operations for cosmogenic nuclide analysis.
//
// The "new erate" inversion pipeline is:
//   1. calculate_production_raster     — Stone 2000 scaling per pixel
//   2. calculate_crn_concentration_raster — steady-state concentration per pixel
//   3. calculate_accumulated_crn_concentration — erosion-weighted average to outlet
//   4. calculate_eff_erate_from_conc   — Newton-Raphson + bisection inversion
//
// All rasters must be in a projected CRS (metres), and must cover the same
// spatial extent as the DEM used to build the flow routing object.
//
// Units:
//   Concentrations   : atoms g⁻¹
//   Effective e-rate : g cm⁻² yr⁻¹
//   To convert to mm/kyr: multiply by 1e7 / rock_density_kg_m3
// ─────────────────────────────────────────────────────────────────────────────

// ── Production raster ────────────────────────────────────────────────────────
// Computes the Stone (2000) spallation scaling factor at every valid pixel.
// Requires a projected DEM (UTM/polar). Atmospheric pressure is estimated
// from elevation using the standard-atmosphere barometric formula.
// Returns a raster with the same geometry as dem; nodata pixels are preserved.
lsdtt_raster_cpu calculate_production_raster(const lsdtt_raster_cpu& dem);

// ── CRN concentration raster ─────────────────────────────────────────────────
// Computes steady-state cosmogenic concentration at every pixel upslope of
// outlet_node (all other pixels remain nodata).
//
// nuclide       : "Be10", "Al26", or "C14"
// muon_scheme   : parameter scheme name (see lsdtt_crn_parameters::reset_scaling)
// erosion_raster: effective erosion rate (g cm⁻² yr⁻¹), uniform is fine
// prod_raster   : Stone 2000 scaling raster (from calculate_production_raster)
// topo_shield   : topographic shielding raster (dimensionless, 0–1; 1 = unshielded)
// snow_shield   : snow shielding as effective depth (g cm⁻²); 0 = no snow
// self_shield   : self/sample shielding as effective depth (g cm⁻²); 0 = surface
// flow          : flow routing object (for upslope traversal)
// outlet_node   : flat index of the basin outlet node
lsdtt_raster_cpu calculate_crn_concentration_raster(
    const std::string&      nuclide,
    const std::string&      muon_scheme,
    const lsdtt_raster_cpu& erosion_raster,
    const lsdtt_raster_cpu& prod_raster,
    const lsdtt_raster_cpu& topo_shield,
    const lsdtt_raster_cpu& snow_shield,
    const lsdtt_raster_cpu& self_shield,
    const lsdtt_flowinfo& flow,
    std::int64_t            outlet_node);

// ── Accumulated concentration ─────────────────────────────────────────────────
// Returns the erosion-rate-weighted average concentration across all upslope
// pixels (including outlet_node):
//   C_outlet = Σ (C_node × ε_node) / Σ ε_node
// This is the predicted concentration at the outlet for a given erosion field.
double calculate_accumulated_crn_concentration(
    const lsdtt_raster_cpu& conc_raster,
    const lsdtt_raster_cpu& erosion_raster,
    const lsdtt_flowinfo& flow,
    std::int64_t outlet_node);

// ── Erosion-rate inversion ────────────────────────────────────────────────────
// Inverts a measured nuclide concentration to an effective erosion rate using
// Newton-Raphson followed by bisection if N-R does not converge.
//
// measured_conc    : measured nuclide concentration (atoms g⁻¹)
// nuclide / muon_scheme : as in calculate_crn_concentration_raster
// prod_raster / topo_shield / snow_shield / self_shield : shielding rasters
// flow / outlet_node : flow routing and basin outlet
// predicted_conc   : OUTPUT — the predicted concentration at convergence
//
// Returns the best-estimate effective erosion rate (g cm⁻² yr⁻¹).
double calculate_eff_erate_from_conc(
    double                  measured_conc,
    const std::string&      nuclide,
    const std::string&      muon_scheme,
    const lsdtt_raster_cpu& prod_raster,
    const lsdtt_raster_cpu& topo_shield,
    const lsdtt_raster_cpu& snow_shield,
    const lsdtt_raster_cpu& self_shield,
    const lsdtt_flowinfo& flow,
    std::int64_t            outlet_node,
    double&                 predicted_conc);

} // namespace lsdtt
