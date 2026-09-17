// lsdtt_topmodel_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// TOPMODEL column physics, index mode
//
// Phase 2 of the distributed-hydrology work. This is the vertical column of
// Beven & Kirkby (1979) as implemented in TOPMODEL 95.02, generalised so that
// the computational unit is either a ln(a/tanB) HISTOGRAM BIN (which reproduces
// the lumped model) or a RASTER CELL (which is the distributed model). Both run
// through one engine, which is what makes the two comparable at all.
//
// Nothing here routes water laterally and nothing here has a velocity. The
// output is a runoff depth per unit per timestep; Phase 5 routes it.
//
// THE STORES
// ----------
//   S_rz   root zone storage DEFICIT   [m]   0 at field capacity, grows as it dries
//   S_uz   unsaturated zone STORAGE    [m]   0 at field capacity, grows as it wets
//   S_bar  catchment mean saturation DEFICIT [m], one value for the whole domain
//
// The sign conventions are the awkward part and are worth stating plainly:
// S_rz and S_bar are deficits, S_uz is a storage. That is 95.02's convention and
// it is kept, because every published TOPMODEL parameter set assumes it.
//
// INDEX MODE, AND WHAT IT ASSUMES
// -------------------------------
// The local saturation deficit is not a state. It is derived every timestep
// from the catchment mean and the topographic index:
//
//     D_i = S_bar + m * (lambda - TWI_i)          clipped at 0
//
// where lambda is the areal mean of TWI. This is the quasi-steady-state
// assumption: it says the water table is everywhere parallel to a steady
// recharge solution, so one number plus topography fixes the whole pattern.
// It is cheap and it is what the published parameter sets were fitted to, but
// it cannot represent a cell drying at a different rate from its neighbours.
// `dynamic` mode (Phase 7) makes D_i prognostic per cell and drops the
// assumption.
//
// Time stepping is EXPLICIT and uses start-of-step values, again as 95.02 does:
// S_bar at the start of the step sets the contributing area, and the baseflow
// computed at the start is what updates S_bar at the end. With long (daily)
// steps the contributing area therefore depends on the initial value plus any
// volume-filling within the step.
//
// SATURATION EXCESS: TWO RULES, AND WHY BOTH EXIST
// ------------------------------------------------
// A unit whose S_uz exceeds its local deficit spills the difference. How that
// spill is totalled differs between the two representations, and conflating
// them is the single easiest way to get a wrong answer that looks plausible:
//
//   per_unit     each unit contributes its own excess times its own area.
//                Correct for raster cells, where a cell is simply saturated or
//                not.
//
//   trapezoid    95.02's rule for a histogram. Within a bin that is only
//                partly saturated it interpolates the saturated fraction
//                linearly between the bin edges, because a bin is a RANGE of
//                index values, not a single value.
//
// The trapezoid rule is a quadrature approximation of the same integral the
// per-unit rule evaluates exactly. They therefore do NOT agree exactly at a
// finite bin count -- they converge as the bins narrow. Any test asserting
// "binned == raster" must be a convergence test, not an equality test.
//
// References
// ----------
//   Beven K.J. & Kirkby M.J. (1979) Hydrol. Sci. Bull. 24(1):43-69
//   Beven K.J. et al. (1995) TOPMODEL 95.02
//   Beven K.J. et al. (2021) HESS 25:527-549
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace lsdtt
{

// ---------------------------------------------------------------------------
// Transmissivity, and the units trap that comes with it
//
// lsdtt_sti_cpu's `ln_kd` field IS ln(T0): build_sti_lnkd_table_from_csv
// accumulates `ks_scaled * dz_cm` over a profile and takes its log, which is
// conductivity integrated over depth. That makes the measured soil dataset a
// direct source for TOPMODEL's `ln_t0`, per soil series and land use, instead
// of one calibrated number for a whole catchment.
//
// MEASURED 2026-09-17: that accumulation mixes units. `ks_scaled` is metres per
// day (via sti_ks_unit_scale_to_metres_per_day) but `dz_cm` is CENTIMETRES, so
// the product is m.cm/day and not m^2/day -- a factor of 100, which is an
// offset of ln(100) = 4.60517 in log space. It is harmless inside STI, where
// the index is used comparatively and a constant offset shifts every pixel
// equally, but it is NOT harmless here: T0 enters as exp(), so a factor of 100
// in transmissivity is a factor of 100 in baseflow.
//
// These two functions make the conversion explicit and are pinned by
// unit_transmissivity_units_from_sti. Never feed a raw `ln_kd` to
// topmodel_parameters::ln_t0.
// ---------------------------------------------------------------------------

// ln(T0 [m^2/day]) from lsdtt_sti_cpu's ln_kd. Subtracts ln(100).
[[nodiscard]] double sti_ln_kd_to_ln_t0_m2_per_day(double sti_ln_kd) noexcept;

// ln(T0 [m^2/h]) from lsdtt_sti_cpu's ln_kd. Subtracts ln(100) and ln(24).
// This is the form TOPMODEL wants when the timestep is expressed in hours.
[[nodiscard]] double sti_ln_kd_to_ln_t0_m2_per_hour(double sti_ln_kd) noexcept;

// ---------------------------------------------------------------------------
// Parameters. Names follow docs/developer_notes/cli_parameter_conventions.md.
// ---------------------------------------------------------------------------
struct topmodel_parameters
{
    // The famous m: the exponential store depth controlling how fast
    // transmissivity falls as the deficit grows. [m]
    double szm_metres = 0.032;

    // ln of saturated transmissivity at zero deficit, T0 in m^2/h.
    // See sti_ln_kd_to_ln_t0_m2_per_hour() before sourcing this from STI.
    double ln_t0 = 5.0;

    // Unsaturated-zone time delay per unit deficit. [h/m]
    double td_hours_per_metre = 50.0;

    // Maximum root zone deficit, i.e. plant-available storage. [m]
    double srmax_metres = 0.05;

    // Initial subsurface flow per unit area, used to initialise S_bar. [m/step]
    double q0_metres_per_step = 3.28e-5;

    // Initial root zone deficit. [m]
    double sr0_metres = 0.002;

    // Timestep. [h]
    double dt_hours = 1.0;

    // Areal mean topographic index (lambda), if it must be supplied rather than
    // derived from the units.
    //
    // This exists because 95.02 forms the two quantities from DIFFERENT
    // weights, which is easy to conflate and was conflated here first time:
    // TREAD computes TL trapezoidally from the PLAIN normalised bin areas,
    // before acf is formed at all, while the timestep loop weights fluxes by
    // acf = 0.5(a_i + a_{i+1}). Deriving lambda from acf-weighted units gives
    // 5.608 where 95.02 gives 5.455, which propagates straight into SZQ and
    // every deficit.
    //
    // Leave it NaN to derive lambda from the units. Set it -- via
    // topmodel_lambda_from_histogram() -- when reproducing 95.02.
    double lambda_override = std::numeric_limits<double>::quiet_NaN();
};

// ---------------------------------------------------------------------------
// The areal mean topographic index the way 95.02's TREAD forms it: a
// trapezoidal integral over the PLAIN normalised bin areas,
//
//     lambda = SUM_{j>=1} a_j * (idx_j + idx_{j-1}) / 2
//
// Note the sum starts at the second bin, so the first bin's area contributes
// to the normalisation but not to lambda. That is 95.02's indexing, not an
// oversight.
// ---------------------------------------------------------------------------
[[nodiscard]]
double topmodel_lambda_from_histogram(
    const std::vector<double>& area_fractions,
    const std::vector<double>& topographic_index_values
);

// ---------------------------------------------------------------------------
// One computational unit: a TWI bin, or a raster cell.
// ---------------------------------------------------------------------------
struct topmodel_unit
{
    // ln(a/tanB). For a bin this is the value at the bin EDGE, which is what
    // the trapezoid rule interpolates between.
    double topographic_index = 0.0;

    // Fractional area. Must sum to 1 across all units.
    double area_fraction = 0.0;
};

enum class topmodel_excess_rule
{
    per_unit,    // raster cells
    trapezoid    // 95.02 histogram bins
};

// ---------------------------------------------------------------------------
// Per-timestep diagnostics. All depths are metres per timestep, averaged over
// the whole domain, so they add up as a water balance directly.
// ---------------------------------------------------------------------------
struct topmodel_step_result
{
    double rainfall            = 0.0;  // p
    double potential_et        = 0.0;  // ep
    double actual_et           = 0.0;  // sumae contribution
    double saturation_excess   = 0.0;  // qof
    double unsaturated_drainage= 0.0;  // quz
    double baseflow            = 0.0;  // qb
    double total_runoff        = 0.0;  // qof + qb
    double mean_deficit        = 0.0;  // S_bar at END of step
    double contributing_area   = 0.0;  // saturated fraction [-]
};

// ---------------------------------------------------------------------------
// The water balance. Positive is storage gained.
//
// closure_error() is the number to watch: it is a small difference of large
// numbers, so it is accumulated in double and it is the single best tripwire
// this model has for a structural mistake.
// ---------------------------------------------------------------------------
struct topmodel_water_balance
{
    double total_rainfall      = 0.0;
    double total_actual_et     = 0.0;
    double total_runoff        = 0.0;
    double initial_storage     = 0.0;
    double final_storage       = 0.0;

    // Water discarded by the unsaturated store's small-value truncation.
    //
    // TOPMODEL 95.02 contains `if(stor_unsat_zone[ia] < 0.0000001)
    // stor_unsat_zone[ia] = 0.0;`, which does not move that water anywhere --
    // it deletes it. The line is kept here because the published parameter sets
    // were fitted with it, but the deleted volume is ACCOUNTED rather than left
    // to leak: it is accumulated here and enters closure_error() as an explicit
    // sink.
    //
    // That turns an invisible 1e-6 relative imbalance into a number you can
    // look at. If it is ever large relative to runoff, the truncation threshold
    // is interacting with the timestep and should be questioned.
    double truncation_loss     = 0.0;

    [[nodiscard]] double closure_error() const noexcept
    {
        return (final_storage - initial_storage)
             - (total_rainfall - total_actual_et - total_runoff - truncation_loss);
    }
};

// ---------------------------------------------------------------------------
// The engine.
//
// Construct with the units and parameters, then advance one timestep at a time.
// The same object serves the lumped and the distributed model; only the units
// and the excess rule differ.
// ---------------------------------------------------------------------------
class lsdtt_topmodel_index_cpu
{
public:
    // Throws std::invalid_argument for an empty unit list, area fractions that
    // do not sum to 1 within 1e-6, a non-positive szm/srmax/td/dt, or a
    // non-positive q0 (which would make the initial S_bar undefined, since it
    // is found from log(q0/szq)).
    lsdtt_topmodel_index_cpu(std::vector<topmodel_unit> units,
                             const topmodel_parameters& params,
                             topmodel_excess_rule rule = topmodel_excess_rule::per_unit);

    // Advance one step. rainfall and potential_et are depths for the step [m].
    topmodel_step_result step(double rainfall_metres, double potential_et_metres);

    [[nodiscard]] double mean_deficit() const noexcept { return sbar_; }
    [[nodiscard]] double lambda() const noexcept { return lambda_; }
    [[nodiscard]] double szq() const noexcept { return szq_; }

    // Domain-mean storage, for the balance. Deficits enter negative.
    [[nodiscard]] double total_storage() const noexcept;

    [[nodiscard]] const topmodel_water_balance& water_balance() const noexcept
    { return balance_; }

    // Local saturation deficit of each unit under the CURRENT S_bar [m].
    [[nodiscard]] std::vector<double> local_deficits() const;

    // The per-unit stores. Exposed because their SIGN is a correctness
    // statement that the water balance cannot make: a negative unsaturated
    // storage still conserves mass (the deficit simply carries the other way)
    // while being physically impossible. A mutation removing the drainage cap
    // survived every other check in this module for exactly that reason.
    [[nodiscard]] const std::vector<double>& unsaturated_storage() const noexcept
    { return s_uz_; }
    [[nodiscard]] const std::vector<double>& root_zone_deficit() const noexcept
    { return s_rz_; }

    [[nodiscard]] std::size_t n_units() const noexcept { return units_.size(); }

private:
    std::vector<topmodel_unit> units_;
    topmodel_parameters        params_;
    topmodel_excess_rule       rule_;

    std::vector<double> s_rz_;   // root zone deficit per unit
    std::vector<double> s_uz_;   // unsaturated storage per unit

    double lambda_ = 0.0;        // areal mean topographic index
    double szq_    = 0.0;        // exp(ln_t0 + ln(dt) - lambda)
    double sbar_   = 0.0;

    topmodel_water_balance balance_;
};

// ---------------------------------------------------------------------------
// Build histogram-bin units from a ln(a/tanB) distribution, as a subcat.dat
// file supplies it: area fraction per bin, index value at the bin edge.
//
// The areal mean lambda is formed the way 95.02 forms it -- trapezoidally
// between adjacent bin edges -- so a model built this way reproduces the lumped
// one. Normalises the area fractions if they do not already sum to 1.
//
// Throws std::invalid_argument if the two vectors differ in length, if there
// are fewer than two bins, or if the areas sum to zero.
// ---------------------------------------------------------------------------
[[nodiscard]]
std::vector<topmodel_unit> make_units_from_histogram(
    const std::vector<double>& area_fractions,
    const std::vector<double>& topographic_index_values
);

// ---------------------------------------------------------------------------
// Histogram bins using TOPMODEL 95.02's own area convention.
//
// 95.02 does NOT weight bin i by its own area. TREAD/TOPMOD weight it by
//
//     acf_i = 0.5 * ( a_i + a_{i+1} )        with a_{n+1} = 0
//
// i.e. the mean of this bin's area and the next one's, because a bin edge sits
// between two areas. The weights therefore sum to 1 - 0.5*a_1 - 0.5*a_n rather
// than to 1, and 95.02 does not renormalise. For the published demo catchment
// that shortfall is 5e-7, which is why the difference is easy to miss and why
// it is spelled out here.
//
// Use this when reproducing 95.02. Use make_units_from_histogram() otherwise.
//
// `renormalise` divides through so the weights sum to 1 exactly. That is the
// defensible choice for new work and the wrong one for bit-matching the
// original, so it is a parameter rather than a decision made here.
// ---------------------------------------------------------------------------
[[nodiscard]]
std::vector<topmodel_unit> make_units_from_histogram_95_02(
    const std::vector<double>& area_fractions,
    const std::vector<double>& topographic_index_values,
    bool renormalise = true
);

// ---------------------------------------------------------------------------
// Build one unit per valid cell from a TWI raster, each weighted equally.
// This is the distributed model.
//
// Throws std::invalid_argument if the raster has no valid cells.
// ---------------------------------------------------------------------------
[[nodiscard]]
std::vector<topmodel_unit> make_units_from_twi(
    const std::vector<float>&        twi_values,
    const std::vector<std::uint8_t>& twi_mask
);

} // namespace lsdtt
