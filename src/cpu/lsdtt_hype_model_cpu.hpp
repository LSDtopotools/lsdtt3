// lsdtt_hype_model_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The HYPE model of lsdtt-distributed-hydrology: sub-basins x classes, routed
//
// Assembles the HYPE land column (lsdtt_hype_column_cpu) and HYPE river
// routing (lsdtt_hype_routing_cpu) over a set-up (lsdtt_hype_parameters):
//
//   * UNITS are the (sub-basin, class) pairs with a non-zero area fraction,
//     stored in one flat array ordered by sub-basin then class. Sub-basin s
//     owns units unit_offset[s] .. unit_offset[s+1]-1.
//     Aggregating class runoff to a sub-basin is a loop over that contiguous
//     range in fixed order, so the result does not depend on how units are
//     stepped (the basis of Stage B's thread-count invariance).
//   * Each step: every unit runs the column; sub-basin runoff is the
//     area-fraction-weighted class runoff (model_hype.f90:2361); the network
//     is routed upstream-first (hype_route_network_step).
//   * The domain water balance accounts for snow, soil water, river damping
//     boxes and water in transit in the translation queues (HYPE's own river
//     volume, model_hype.f90:1424-1425). It is optional (track_balance): the
//     storage term is a full pass over units and river queues every step.
//   * The land loop runs in parallel over units (OpenMP, schedule(static)).
//     Each unit touches only its own state, and the per-sub-basin sums are a
//     separate serial loop in fixed order (no reduction clause), so results
//     are bitwise identical at any thread count.
//
// Sub-basin order is the set-up's order; ids are kept for output. Forcing is
// one precipitation and one temperature value per step for the whole domain
// (first version), optionally with observed PET. Daily
// steps only.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#include "lsdtt_hype_column_cpu.hpp"
#include "lsdtt_hype_parameters.hpp"
#include "lsdtt_hype_routing_cpu.hpp"

namespace lsdtt
{

struct hype_model_unit
{
    std::size_t subbasin = 0;    // index into the set-up's sub-basins
    int         class_id = 0;
    double      fraction = 0.0;  // area fraction of the sub-basin
    hype_unit_constants constants;
    hype_unit_state     state;
    // The full per-step flux record is kept outside the unit
    // (lsdtt_hype_model_cpu::unit_fluxes()) so the unit stays small.
};

// Domain totals of the last step and the running balance [m3].
//
// With balance tracking off (lsdtt_hype_model_cpu, track_balance = false) the
// cumulative precipitation, evaporation and outflow are still kept, since they
// are by-products of the step. The storage terms are not computed: both are
// NaN and storage_tracked is false, so closure_error_m3() is NaN.
struct hype_model_balance
{
    bool   storage_tracked    = true;  // false: storage terms are NaN
    double storage_start_m3   = 0.0;   // storage before the first step
    double storage_m3         = 0.0;   // storage now
    double precipitation_m3   = 0.0;   // cumulative, after interception
    double evaporation_m3     = 0.0;   // cumulative
    double outflow_m3         = 0.0;   // cumulative, at domain outlets

    // (storage now - start) - (precipitation - evaporation - outflow);
    // NaN when storage is not tracked.
    [[nodiscard]] double closure_error_m3() const noexcept
    {
        return (storage_m3 - storage_start_m3) -
               (precipitation_m3 - evaporation_m3 - outflow_m3);
    }
};

class lsdtt_hype_model_cpu
{
public:
    // Validates the set-up (validate_hype_setup) and the first-version
    // options. dt_seconds must be 86400 (daily); anything else throws.
    // track_balance = false skips the domain storage pass on every step (see
    // hype_model_balance); flows and states are identical either way.
    explicit lsdtt_hype_model_cpu(const hype_setup& setup, double dt_seconds = 86400.0,
                                  bool track_balance = true);

    // Advance one step. observed_pet_mm NaN = use HYPE's calculated PET.
    // If a unit's column fails (hype_unit_step throws), the exception of the
    // lowest-index failing unit is rethrown after the land loop, whatever the
    // thread count; the model state is then partly stepped and not usable.
    void step(double precipitation_mm, double temperature_c, int day_of_year,
              double observed_pet_mm = std::numeric_limits<double>::quiet_NaN());

    [[nodiscard]] std::size_t n_subbasins() const noexcept { return subids_.size(); }
    [[nodiscard]] const std::vector<int>& subids() const noexcept { return subids_; }

    // Per sub-basin, last step: runoff [mm] (HYPE crun), local-river outflow
    // and sub-basin outflow [m3/s] (HYPE cloc, cout).
    [[nodiscard]] const std::vector<double>& runoff_mm() const noexcept { return runoff_mm_; }
    [[nodiscard]] const std::vector<double>& local_outflow_m3s() const noexcept { return flows_.local_outflow_m3s; }
    [[nodiscard]] const std::vector<double>& outflow_m3s() const noexcept { return flows_.outflow_m3s; }

    // Units in sub-basin-then-class order, and the offsets into them.
    [[nodiscard]] const std::vector<hype_model_unit>& units() const noexcept { return units_; }
    [[nodiscard]] const std::vector<std::size_t>& unit_offset() const noexcept { return unit_offset_; }

    // Per-unit flux records of the last step, in units() order.
    // Kept by default; set_store_unit_fluxes(false) stops writing them, which
    // cuts per-step memory traffic when only sub-basin results are needed.
    // When off, unit_fluxes() is empty.
    void set_store_unit_fluxes(bool store);
    [[nodiscard]] bool store_unit_fluxes() const noexcept { return store_unit_fluxes_; }
    [[nodiscard]] const std::vector<hype_unit_fluxes>& unit_fluxes() const noexcept { return unit_fluxes_; }

    [[nodiscard]] const hype_model_balance& balance() const noexcept { return balance_; }
    [[nodiscard]] bool balance_tracked() const noexcept { return track_balance_; }

    // FOR TESTS ONLY. From the next step on, every unit with index >= unit_index
    // throws std::runtime_error("hype model: injected test failure in unit N")
    // inside the parallel land loop. Validation leaves no input that makes the
    // column fail, so this is the only way to exercise exception propagation
    // out of the OpenMP region. npos (the default) disables it.
    void inject_failure_from_unit_for_testing(std::size_t unit_index) noexcept { inject_from_unit_ = unit_index; }

    // Domain storage [m3]: snow and soil water of every unit, plus river
    // damping boxes and water in transit in the translation queues. Always
    // computed on request, whether or not the balance is tracked.
    [[nodiscard]] double storage_m3() const;

private:
    std::size_t inject_from_unit_ = static_cast<std::size_t>(-1);   // test hook; npos = off
    hype_general_parameters  general_;
    double                   dt_seconds_ = 86400.0;
    bool                     track_balance_ = true;
    bool                     store_unit_fluxes_ = true;
    std::vector<hype_unit_fluxes> unit_fluxes_;      // per unit, when stored
    std::vector<double>      unit_runoff_mm_;        // per unit, always: what aggregation needs
    std::vector<double>      unit_prec_mm_;
    std::vector<double>      unit_evap_mm_;
    std::vector<int>         subids_;
    std::vector<hype_model_unit> units_;
    std::vector<std::size_t> unit_offset_;
    hype_river_network       net_;
    std::vector<double>      runoff_mm_;
    hype_network_flows       flows_;
    hype_model_balance       balance_;
};

// Water in a river [m3]: damping box plus water in transit in the queue
// (model_hype.f90:1424).
[[nodiscard]] double hype_river_volume_m3(const hype_river& r) noexcept;

} // namespace lsdtt
