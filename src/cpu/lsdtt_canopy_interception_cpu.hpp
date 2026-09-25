// lsdtt_canopy_interception_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Canopy interception, store [1] of the column
//
// Phase 4 of the distributed-hydrology work. A Rutter-type bucket sitting above
// the root zone: rainfall partly misses the canopy, the rest fills a store that
// drains as throughfall once full, and the stored water evaporates.
//
// INTERCEPTION IS NOT A BOLT-ON -- IT RE-PARTITIONS THE EVAPORATIVE DEMAND
// ------------------------------------------------------------------------
// This is the single easiest thing to get wrong here, so it is stated first.
// Adding an interception store does NOT leave the rest of the column's
// behaviour intact with a bit of water skimmed off the top. Evaporation from a
// WET CANOPY happens at the potential rate and suppresses transpiration while
// it lasts, so only the RESIDUAL demand reaches the root zone:
//
//     E_canopy   = min( S_canopy , PET * evap_factor )
//     PET_residual = max( 0 , PET - E_canopy / evap_factor )
//
// The division by evap_factor is the part that is easy to drop. `evap_factor`
// exists because wet-canopy evaporation can EXCEED grass-reference PET for tall
// aerodynamically rough vegetation (Rutter; Calder) -- the canopy is better
// ventilated than the reference surface. With evap_factor = 2, one millimetre
// of canopy evaporation only consumes half a millimetre of transpiration
// demand, so total evapotranspiration legitimately exceeds PET. Subtracting
// E_canopy directly instead would silently cap total ET at PET and understate
// forest water use.
//
// LAI DRIVES BOTH CANOPY PARAMETERS
// ---------------------------------
//     C = k_S * LAI                canopy storage capacity, k_S ~ 0.2 mm
//     p = exp( -kappa * LAI )      free throughfall fraction, kappa ~ 0.5-0.7
//
// C from Dickinson (1984) and the S ~ S_L LAI form of van Dijk & Bruijnzeel
// (2001); the Beer-Lambert gap fraction for p from Aston (1979).
//
// SEASONALITY IS DECISIVE, NOT DECORATIVE
// ---------------------------------------
// Winter deciduous LAI is ~0.5-1 against ~5 in summer -- an order of magnitude
// in interception capacity. UK flood events are winter-dominated, so newly
// planted deciduous woodland intercepts LEAST exactly when it matters most. A
// study that runs winter storms with annual-mean or summer LAI will overstate
// the benefit of planting.
//
// `lai_seasonal_factor` therefore multiplies LAI at construction, and it is a
// parameter a caller must think about rather than a default that quietly
// applies summer conditions to a January storm.
//
// References
// ----------
//   Rutter A.J. et al. (1971) Agric. Meteorol. 9:367-384 -- the store model
//   Gash J.H.C. (1979) QJRMS 105:43-55 -- the analytical alternative
//   Aston A.R. (1979) J. Hydrol. 42:383-396 -- gap fraction
//   Dickinson R.E. (1984) -- C = 0.2 LAI
//   van Dijk A.I.J.M. & Bruijnzeel L.A. (2001) J. Hydrol. 247:230-238
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace lsdtt
{

// ---------------------------------------------------------------------------
// Closed forms, exposed for testing and per-cell use.
// ---------------------------------------------------------------------------

// Canopy storage capacity [m] from LAI. `mm_per_lai` is k_S in millimetres, so
// the conversion to metres happens here rather than at every call site.
[[nodiscard]] double canopy_capacity_from_lai(double lai,
                                              double mm_per_lai) noexcept;

// Free throughfall fraction [-]: the share of rain that reaches the ground
// without touching the canopy. exp(-kappa LAI), so bare ground passes
// everything and a dense canopy passes almost nothing.
[[nodiscard]] double free_throughfall_from_lai(double lai,
                                               double gap_coefficient) noexcept;

// ---------------------------------------------------------------------------
struct canopy_parameters
{
    // C = canopy_storage_mm_per_lai * LAI   [mm per unit LAI]
    double canopy_storage_mm_per_lai = 0.2;

    // kappa in p = exp(-kappa LAI)
    double lai_gap_coefficient = 0.6;

    // Wet-canopy evaporation as a multiple of PET. >= 1 for tall rough canopies.
    double interception_evap_factor = 1.0;

    // Multiplies every LAI at construction. SEE THE HEADER: leaving this at 1.0
    // for a winter storm on deciduous woodland overstates interception by up to
    // an order of magnitude.
    double lai_seasonal_factor = 1.0;
};

struct canopy_step_result
{
    // All depths are domain-mean metres for the step, so they close as a
    // balance directly.
    double rainfall           = 0.0;
    double throughfall        = 0.0;   // reaches the root zone
    double canopy_evaporation = 0.0;   // leaves to the atmosphere
    double potential_et       = 0.0;
    double residual_pet       = 0.0;   // what the root zone may still draw on
    double canopy_storage     = 0.0;   // domain-mean store at END of step
};

struct canopy_water_balance
{
    double total_rainfall           = 0.0;
    double total_throughfall        = 0.0;
    double total_canopy_evaporation = 0.0;
    double initial_storage          = 0.0;
    double final_storage            = 0.0;

    [[nodiscard]] double closure_error() const noexcept
    {
        return (final_storage - initial_storage)
             - (total_rainfall - total_throughfall - total_canopy_evaporation);
    }
};

// ---------------------------------------------------------------------------
// The canopy engine.
//
// Takes LAI and area fraction per unit, using the SAME unit list as
// lsdtt_topmodel_index_cpu. That alignment is deliberate: two independently
// ordered unit lists would compose into a silently wrong model, and there is no
// diagnostic that would catch it.
// ---------------------------------------------------------------------------
class lsdtt_canopy_interception_cpu
{
public:
    // Throws std::invalid_argument if the vectors differ in length, are empty,
    // contain a negative or non-finite LAI, have area fractions not summing to
    // 1 within 1e-6, or carry a non-positive parameter.
    lsdtt_canopy_interception_cpu(std::vector<double> lai_per_unit,
                                  std::vector<double> area_fraction,
                                  const canopy_parameters& params);

    // Advance one step with uniform forcing. Per-unit throughfall and residual
    // PET are written to the output vectors, which are resized as needed; both
    // may be the same object only if the caller does not need both.
    canopy_step_result step(double rainfall_metres,
                            double potential_et_metres,
                            std::vector<double>& throughfall_per_unit,
                            std::vector<double>& residual_pet_per_unit);

    // Convenience for the lumped case, where only the aggregates are wanted.
    canopy_step_result step(double rainfall_metres, double potential_et_metres);

    [[nodiscard]] const std::vector<double>& canopy_storage() const noexcept
    { return storage_; }
    [[nodiscard]] const std::vector<double>& capacity() const noexcept
    { return capacity_; }
    [[nodiscard]] const std::vector<double>& free_throughfall() const noexcept
    { return free_tf_; }
    [[nodiscard]] const canopy_water_balance& water_balance() const noexcept
    { return balance_; }
    [[nodiscard]] std::size_t n_units() const noexcept { return capacity_.size(); }

private:
    std::vector<double> area_;
    std::vector<double> capacity_;
    std::vector<double> free_tf_;
    std::vector<double> storage_;
    canopy_parameters   params_;
    canopy_water_balance balance_;
};

// ---------------------------------------------------------------------------
// LAI per unit from a land-cover class raster plus a class -> LAI table.
//
// This is the initial ingest path. A direct LAI raster is the later one and
// needs no helper: its values ARE the per-unit LAI.
//
// `class_values` is one land-cover class code per unit. Codes absent from the
// table take `default_lai` and are counted in `n_unmatched` -- a silent zero
// there would remove interception from exactly the cells nobody checked.
// ---------------------------------------------------------------------------
struct lai_from_landcover_result
{
    std::vector<double> lai;
    std::int64_t        n_unmatched = 0;
};

[[nodiscard]]
lai_from_landcover_result lai_from_landcover_classes(
    const std::vector<int>&                    class_values,
    const std::unordered_map<int, double>&     lai_by_class,
    double                                     default_lai = 0.0
);

} // namespace lsdtt
