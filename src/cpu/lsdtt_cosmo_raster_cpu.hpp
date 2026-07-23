#pragma once
// =============================================================================
//  DO NOT CHANGE THE NUMBERS THIS CODE PRODUCES — ask Simon Mudd first.
//
//  This applies to the concentration and erosion-rate (inversion) calculations
//  in this file. The raster reading/writing and flow-routing helpers here are
//  ordinary code and can be changed normally.
//
//  Those cosmogenic calculations were checked by hand over several years. To
//  stop anyone changing their results by accident, the test suite keeps a stored
//  list of known-correct output values in
//  tests/expected/cosmo/numerics_sentinel.csv. The tests in
//  tests/cpp/test_cosmo_physics.cpp re-run these calculations and compare the
//  results against that stored list to 12 decimal places.
//
//  So if you edit those calculations and any output value changes — even from a
//  refactor or an algebra rewrite you believe is harmless — those tests will
//  fail. When that happens, do NOT edit the stored list to make the tests pass
//  again. A changed value has to be reviewed and approved by Simon Mudd before
//  it is committed.
// =============================================================================
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_crn_parameters.hpp"
#include <string>
#include <vector>

namespace lsdtt
{

// 
// lsdtt_cosmo_raster_cpu
//
// Raster-level operations for cosmogenic nuclide analysis.
//
// The "new erate" inversion pipeline is:
//   1. calculate_production_raster     - Stone 2000 scaling per pixel
//   2. calculate_crn_concentration_raster - steady-state concentration per pixel
//   3. calculate_accumulated_crn_concentration - erosion-weighted average to outlet
//   4. calculate_eff_erate_from_conc   - Newton-Raphson + bisection inversion
//
// All rasters must be in a projected CRS (metres), and must cover the same
// spatial extent as the DEM used to build the flow routing object.
//
// Units:
//   Concentrations   : atoms g^-^1
//   Effective e-rate : g cm^-^2 yr^-^1
//   To convert to mm/kyr: multiply by 1e7 / rock_density_kg_m3
// 

//  Production raster 
// Computes the Stone (2000) spallation scaling factor at every valid pixel.
// Requires a projected DEM (UTM/polar).
//
// nc      : ERA40 / paleomagnetic data loaded via lsdtt_crn_parameters::load_nc_data().
//           When nc.loaded == false, falls back to the standard-atmosphere formula and
//           no paleomagnetic correction (equivalent to the single-argument overload).
// age_yr  : Sample exposure age in years.  When > 0 and nc.loaded == true, applies
//           the Muscheler VDM time-integrated paleomagnetic correction to the Stone
//           scaling factor.  When 0 (default), modern field only.
//
// Returns a raster with the same geometry as dem; nodata pixels are preserved.
lsdtt_raster_cpu calculate_production_raster(const lsdtt_raster_cpu& dem,
                                              const CrnNcData&         nc,
                                              double                   age_yr = 0.0);

// Backward-compatible overload: standard atmosphere, no paleomagnetic correction.
inline lsdtt_raster_cpu calculate_production_raster(const lsdtt_raster_cpu& dem)
{
    return calculate_production_raster(dem, CrnNcData{}, 0.0);
}

//  CRN concentration raster 
// Computes steady-state cosmogenic concentration at every pixel upslope of
// outlet_node (all other pixels remain nodata).
//
// nuclide       : "Be10", "Al26", or "C14"
// muon_scheme   : parameter scheme name (see lsdtt_crn_parameters::reset_scaling)
// erosion_raster: effective erosion rate (g cm^-^2 yr^-^1), uniform is fine
// prod_raster   : Stone 2000 scaling raster (from calculate_production_raster)
// topo_shield   : topographic shielding raster (dimensionless, 0-1; 1 = unshielded)
// snow_shield   : snow shielding as effective depth (g cm^-^2); 0 = no snow
// self_shield   : self/sample shielding as effective depth (g cm^-^2); 0 = surface
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

//  Accumulated concentration 
// Returns the erosion-rate-weighted average concentration across all upslope
// pixels (including outlet_node):
//   C_outlet = Sigma (C_node x epsilon_node) / Sigma epsilon_node
// This is the predicted concentration at the outlet for a given erosion field.
double calculate_accumulated_crn_concentration(
    const lsdtt_raster_cpu& conc_raster,
    const lsdtt_raster_cpu& erosion_raster,
    const lsdtt_flowinfo& flow,
    std::int64_t outlet_node);

//  Whole-raster steady-state CRN concentration 
// Same physics as calculate_crn_concentration_raster but operates on every
// valid pixel (no flow routing / outlet restriction needed).
// nuclide       : "Be10", "Al26", or "C14"
// muon_scheme   : see lsdtt_crn_parameters::reset_scaling
// erosion_raster: effective erosion rate (g cm^-^2 yr^-^1) per pixel
// prod_raster   : Stone 2000 scaling raster (from calculate_production_raster)
// topo_shield   : topographic shielding (dimensionless 0-1; 1 = unshielded)
// snow_shield   : snow shielding as effective depth (g cm^-^2); 0 = no snow
// self_shield   : self shielding as effective depth (g cm^-^2); 0 = surface
lsdtt_raster_cpu calculate_crn_concentration_raster_whole(
    const std::string&      nuclide,
    const std::string&      muon_scheme,
    const lsdtt_raster_cpu& erosion_raster,
    const lsdtt_raster_cpu& prod_raster,
    const lsdtt_raster_cpu& topo_shield,
    const lsdtt_raster_cpu& snow_shield,
    const lsdtt_raster_cpu& self_shield);

//  Step-change CRN concentration (whole raster) 
// Computes the transient concentration after an instantaneous step change in
// erosion rate.  Surface sample (effective depth = 0).
//
// At pixels where time_since_change > 0:
//   1. start depth = eff_erate_new x time_since_change
//   2. initialise SS at eff_erate_old at that starting depth
//   3. evolve for time_since_change yr at eff_erate_new -> surface concentration
// At pixels where time_since_change <= 0: returns SS at eff_erate_new.
//
// eff_erate_old_g_cm2_yr / eff_erate_new_g_cm2_yr : uniform scalar erosion rates
// time_since_change : raster of years since the step change occurred (0 = no change)
lsdtt_raster_cpu calculate_crn_concentration_raster_step_change(
    const std::string&      nuclide,
    const std::string&      muon_scheme,
    const lsdtt_raster_cpu& time_since_change,
    const lsdtt_raster_cpu& prod_raster,
    const lsdtt_raster_cpu& topo_shield,
    double                  eff_erate_old_g_cm2_yr,
    double                  eff_erate_new_g_cm2_yr);

//  Whole-raster CRN accumulation 
// Computes the erosion-weighted mean concentration at every pixel using D8
// flow accumulation:
//   C_acc[i] = Sigma_upslope(C x epsilon) / Sigma_upslope(epsilon)
// Returns nodata where the accumulated erosion is zero.
// Requires flow routing built with build_topological_order = true.
lsdtt_raster_cpu calculate_accumulated_crn_raster(
    const lsdtt_raster_cpu& conc_raster,
    const lsdtt_raster_cpu& erosion_raster,
    const lsdtt_flowinfo&   flow);

//  Erosion-rate inversion 
// Inverts a measured nuclide concentration to an effective erosion rate using
// Newton-Raphson followed by bisection if N-R does not converge.
//
// measured_conc    : measured nuclide concentration (atoms g^-^1)
// nuclide / muon_scheme : as in calculate_crn_concentration_raster
// prod_raster / topo_shield / snow_shield / self_shield : shielding rasters
// flow / outlet_node : flow routing and basin outlet
// predicted_conc   : OUTPUT - the predicted concentration at convergence
//
// Returns the best-estimate effective erosion rate (g cm^-^2 yr^-^1).
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
