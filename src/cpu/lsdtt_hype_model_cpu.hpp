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
//     owns units unit_offset[s] .. unit_offset[s+1]-1 (plan section 7.3).
//     Aggregating class runoff to a sub-basin is a loop over that contiguous
//     range in fixed order, so the result does not depend on how units are
//     stepped (the basis of Stage B's thread-count invariance).
//   * Each step: every unit runs the column; sub-basin runoff is the
//     area-fraction-weighted class runoff (model_hype.f90:2361); the network
//     is routed upstream-first (hype_route_network_step).
//   * The domain water balance accounts for snow, soil water, river damping
//     boxes and water in transit in the translation queues (HYPE's own river
//     volume, model_hype.f90:1424-1425).
//
// Sub-basin order is the set-up's order; ids are kept for output. Forcing is
// one precipitation and one temperature value per step for the whole domain
// (first version, plan section 7.1), optionally with observed PET. Daily
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
    hype_unit_fluxes    fluxes;  // of the last step
};

// Domain totals of the last step and the running balance [m3].
struct hype_model_balance
{
    double storage_start_m3   = 0.0;   // storage before the first step
    double storage_m3         = 0.0;   // storage now
    double precipitation_m3   = 0.0;   // cumulative, after interception
    double evaporation_m3     = 0.0;   // cumulative
    double outflow_m3         = 0.0;   // cumulative, at domain outlets

    // (storage now - start) - (precipitation - evaporation - outflow)
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
    explicit lsdtt_hype_model_cpu(const hype_setup& setup, double dt_seconds = 86400.0);

    // Advance one step. observed_pet_mm NaN = use HYPE's calculated PET.
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

    [[nodiscard]] const hype_model_balance& balance() const noexcept { return balance_; }

    // Domain storage [m3]: snow and soil water of every unit, plus river
    // damping boxes and water in transit in the translation queues.
    [[nodiscard]] double storage_m3() const;

private:
    hype_general_parameters  general_;
    double                   dt_seconds_ = 86400.0;
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
