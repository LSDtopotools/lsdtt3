// lsdtt_niemi_ensemble_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Basin sampling, apparent-versus-true erosion rate, and the area statistics
// that reproduce Niemi et al. (2005).
//
// WHAT NIEMI'S RESULT ACTUALLY IS. It is not a number, it is a DISTRIBUTION:
// basin-averaged cosmogenic concentrations in a landsliding landscape are biased
// and highly variable, the variability falls as basin area grows, and only above
// some area does a single sample recover the true erosion rate reliably. His own
// conclusions are (i) that a single sample in a small landsliding basin is
// uninterpretable, and (ii) that the MEDIAN over several similar basins is
// nevertheless representative. Reproducing that needs an ensemble, not a run.
//
// THE TWO RATES BEING COMPARED, precisely:
//
//   TRUE      the volumetric erosion rate -- the mass actually leaving each
//             pixel per unit area per year. That is exactly the DELIVERY-MASS
//             field the landslide coupling already produces (background erosion
//             everywhere, plus the excavated slab spread over the flush window
//             where a scar is delivering). Averaged over the basin by area.
//
//   APPARENT  what an analyst gets: invert the flux-weighted basin-averaged
//             concentration through the ordinary CAIRN inversion, as though the
//             landscape were eroding steadily. This is calculate_eff_erate_from_conc,
//             the same routine lsdtt-cosmo-tool uses, so the apparent rate here
//             is the number the tool would report for that basin.
//
// Their ratio is Niemi's E_c/E_v. A control run with no landslides must return
// 1 to within the inversion's own tolerance; anything else is a defect in the
// forward/inverse pair rather than a result about landslides.
//
// WHY THE SAMPLING IS STRATIFIED BY AREA. The claim under test is about the
// AREA dependence, so basins must span decades of area. Drawing outlets
// uniformly from the network would swamp the sample with small headwater basins
// -- contributing area is power-law distributed -- and leave the large end,
// where the interesting convergence happens, with almost no points.
// -----------------------------------------------------------------------------

#ifndef LSDTT_NIEMI_ENSEMBLE_CPU_HPP
#define LSDTT_NIEMI_ENSEMBLE_CPU_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// One sampled basin in one realisation.
// -----------------------------------------------------------------------------
struct NiemiBasinSample
{
    int          realisation      = 0;
    std::int64_t outlet_flat      = -1;
    std::int64_t n_pixels         = 0;
    double       area_m2          = 0.0;
    double       area_km2         = 0.0;

    std::int64_t n_failed_pixels  = 0;
    double       fraction_failed  = 0.0;

    double accumulated_conc       = 0.0;   // atoms/g, flux-weighted basin mean
    double apparent_erate_g_cm2_yr = 0.0;  // inverted as an analyst would
    double true_erate_g_cm2_yr     = 0.0;  // area mean of the delivery mass
    double ratio                   = 0.0;  // apparent / true

    bool   valid = false;                  // false if the inversion failed
};

// -----------------------------------------------------------------------------
// Statistics of apparent/true within one band of basin area.
// -----------------------------------------------------------------------------
struct NiemiAreaBin
{
    double area_lo_km2 = 0.0;
    double area_hi_km2 = 0.0;
    int    n           = 0;

    double median = 0.0;
    double p25 = 0.0, p75 = 0.0;
    double p05 = 0.0, p95 = 0.0;
    double iqr = 0.0;                 // the spread that must fall with area
    double frac_within_10pc = 0.0;
    double frac_within_20pc = 0.0;
};

// -----------------------------------------------------------------------------
// sample_basin_outlets_stratified
//
// Outlets spanning the available range of contributing area, stratified by
// decade so the large end is represented. Candidates are every cell whose
// upslope count is in [min_pixels, max_pixels]; within each log-area band up to
// per_band are taken, spread through the band rather than clustered.
//
// Deterministic: the same flow object and arguments always give the same list,
// so an ensemble varies only in its landslide seed.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::int64_t> sample_basin_outlets_stratified(
    const lsdtt_flow_base& flow,
    std::int64_t           min_pixels,
    std::int64_t           max_pixels,
    int                    per_band,
    int                    bands_per_decade = 2);

// -----------------------------------------------------------------------------
// evaluate_niemi_basins
//
// For each outlet: form the flux-weighted basin-averaged concentration, invert
// it to an apparent erosion rate, and compare with the area-averaged true rate.
//
// weight_raster is what the mixture is weighted by -- the delivery-mass or
// quartz-weight field, NOT the erosion rate, whenever landslides are active.
// true_rate_raster is in g/cm2/yr.
//
// Rasters must all share the flow object's grid.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<NiemiBasinSample> evaluate_niemi_basins(
    const lsdtt_flowinfo&            flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu&          concentration,
    const lsdtt_raster_cpu&          weight_raster,
    const lsdtt_raster_cpu&          true_rate_g_cm2_yr,
    const lsdtt_raster_cpu&          excavation_depth_m,
    const lsdtt_raster_cpu&          prod_raster,
    const lsdtt_raster_cpu&          topo_shield,
    const lsdtt_raster_cpu&          snow_shield,
    const lsdtt_raster_cpu&          self_shield,
    const std::string&               nuclide,
    const std::string&               muon_scheme,
    int                              realisation = 0);

// -----------------------------------------------------------------------------
// bin_niemi_by_area
//
// Bands are logarithmic in area, because the claim is about decades.
// Bins with fewer than min_per_bin samples are dropped rather than reported on
// two points.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<NiemiAreaBin> bin_niemi_by_area(
    const std::vector<NiemiBasinSample>& samples,
    int                                  bands_per_decade = 2,
    int                                  min_per_bin      = 5);

// -----------------------------------------------------------------------------
// niemi_averaging_area_km2
//
// The area above which a stated fraction of samples fall within a stated
// tolerance of the true rate -- our equivalent of Niemi's Eq. 7 rule of thumb,
//
//     A_avg [km2] = 100 [km3/Myr] / E [km/Myr]
//
// Returns the lower edge of the smallest area band at and above which every
// band satisfies the criterion; negative if no such band exists, which is a
// meaningful answer (the landscape never averages out within the range sampled).
// -----------------------------------------------------------------------------
[[nodiscard]] double niemi_averaging_area_km2(
    const std::vector<NiemiAreaBin>& bins,
    double                           tolerance = 0.20,
    double                           fraction  = 0.90);

// -----------------------------------------------------------------------------
// CSV writers. Orchestration only -- no statistics are computed here that are
// not already in the structs above.
// -----------------------------------------------------------------------------
void write_niemi_basins_csv(const std::string& path,
                            const std::vector<NiemiBasinSample>& samples);

void write_niemi_area_bins_csv(const std::string& path,
                               const std::vector<NiemiAreaBin>& bins);

} // namespace lsdtt

#endif // LSDTT_NIEMI_ENSEMBLE_CPU_HPP
