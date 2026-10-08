// lsdtt_canopy_interception_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_canopy_interception_cpu.hpp"
#include "lsdtt_distributed_hydrology_argument_checks_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace lsdtt
{

double canopy_capacity_from_lai(double lai, double mm_per_lai) noexcept
{
    if (!(lai > 0.0) || !(mm_per_lai > 0.0)) return 0.0;
    return mm_per_lai * lai * 1.0e-3;   // mm -> m
}

double free_throughfall_from_lai(double lai, double gap_coefficient) noexcept
{
    if (!(lai > 0.0)) return 1.0;       // bare ground intercepts nothing
    return std::exp(-gap_coefficient * lai);
}

lsdtt_canopy_interception_cpu::lsdtt_canopy_interception_cpu(
    std::vector<double> lai_per_unit,
    std::vector<double> area_fraction,
    const canopy_parameters& params)
    : area_(std::move(area_fraction)), params_(params)
{
    const char* who = "lsdtt_canopy_interception_cpu";

    if (lai_per_unit.empty())
        throw std::invalid_argument(std::string(who) + ": LAI vector is empty");
    if (lai_per_unit.size() != area_.size())
        throw std::invalid_argument(
            std::string(who) + ": LAI has " + std::to_string(lai_per_unit.size()) +
            " entries but area_fraction has " + std::to_string(area_.size()));

    require_positive_finite(params_.canopy_storage_mm_per_lai,
                     "canopy_storage_mm_per_lai", who);
    require_positive_finite(params_.interception_evap_factor,
                     "interception_evap_factor", who);
    require_positive_finite(params_.lai_seasonal_factor, "lai_seasonal_factor", who);
    if (!(std::isfinite(params_.lai_gap_coefficient) &&
          params_.lai_gap_coefficient >= 0.0))
        throw std::invalid_argument(
            std::string(who) + ": lai_gap_coefficient must be finite and >= 0");

    double area_sum = 0.0;
    for (std::size_t i = 0; i < area_.size(); ++i)
    {
        if (!std::isfinite(lai_per_unit[i]) || lai_per_unit[i] < 0.0)
            throw std::invalid_argument(
                std::string(who) + ": a unit has a negative or non-finite LAI");
        if (!std::isfinite(area_[i]) || area_[i] < 0.0)
            throw std::invalid_argument(
                std::string(who) + ": a unit has a negative or non-finite area");
        area_sum += area_[i];
    }
    if (std::fabs(area_sum - 1.0) > 1.0e-6)
        throw std::invalid_argument(
            std::string(who) + ": area fractions sum to " +
            std::to_string(area_sum) + ", not 1");

    capacity_.resize(lai_per_unit.size());
    free_tf_.resize(lai_per_unit.size());
    storage_.assign(lai_per_unit.size(), 0.0);

    for (std::size_t i = 0; i < lai_per_unit.size(); ++i)
    {
        // Seasonality is applied ONCE, here, so every derived quantity sees the
        // same LAI. Scaling capacity but not the gap fraction -- or the other
        // way round -- would be internally inconsistent.
        const double lai = lai_per_unit[i] * params_.lai_seasonal_factor;
        capacity_[i] = canopy_capacity_from_lai(lai, params_.canopy_storage_mm_per_lai);
        free_tf_[i]  = free_throughfall_from_lai(lai, params_.lai_gap_coefficient);
    }

    balance_.initial_storage = 0.0;
}

canopy_step_result lsdtt_canopy_interception_cpu::step(
    double rainfall_metres,
    double potential_et_metres,
    std::vector<double>& throughfall_per_unit,
    std::vector<double>& residual_pet_per_unit)
{
    const double p  = std::max(0.0, rainfall_metres);
    const double ep = std::max(0.0, potential_et_metres);

    const std::size_t n = capacity_.size();
    throughfall_per_unit.assign(n, 0.0);
    residual_pet_per_unit.assign(n, 0.0);

    canopy_step_result out;
    out.rainfall     = p;
    out.potential_et = ep;

    double tf_sum = 0.0, ec_sum = 0.0, rp_sum = 0.0, st_sum = 0.0;

    for (std::size_t i = 0; i < n; ++i)
    {
        // 1. Free throughfall: rain that never touches the canopy.
        const double direct = free_tf_[i] * p;
        const double to_canopy = p - direct;

        // 2. Fill the store; the overflow drains as canopy throughfall.
        storage_[i] += to_canopy;
        double drainage = 0.0;
        if (storage_[i] > capacity_[i])
        {
            drainage    = storage_[i] - capacity_[i];
            storage_[i] = capacity_[i];
        }

        const double throughfall = direct + drainage;

        // 3. Evaporate from the wet canopy at the potential rate, scaled.
        const double demand = ep * params_.interception_evap_factor;
        const double e_can  = std::min(storage_[i], demand);
        storage_[i] -= e_can;
        if (storage_[i] < 0.0) storage_[i] = 0.0;

        // 4. What is left of the transpiration demand. Dividing by the factor
        //    is essential -- see the header. A rough canopy evaporating at 2x
        //    PET consumes only half as much transpiration demand per unit
        //    evaporated, which is how total ET can exceed PET.
        const double consumed = (params_.interception_evap_factor > 0.0)
            ? e_can / params_.interception_evap_factor : e_can;
        const double residual = std::max(0.0, ep - consumed);

        throughfall_per_unit[i]  = throughfall;
        residual_pet_per_unit[i] = residual;

        const double a = area_[i];
        tf_sum += throughfall * a;
        ec_sum += e_can * a;
        rp_sum += residual * a;
        st_sum += storage_[i] * a;
    }

    out.throughfall        = tf_sum;
    out.canopy_evaporation = ec_sum;
    out.residual_pet       = rp_sum;
    out.canopy_storage     = st_sum;

    balance_.total_rainfall           += p;
    balance_.total_throughfall        += tf_sum;
    balance_.total_canopy_evaporation += ec_sum;
    balance_.final_storage             = st_sum;

    return out;
}

canopy_step_result lsdtt_canopy_interception_cpu::step(double rainfall_metres,
                                                       double potential_et_metres)
{
    std::vector<double> tf, rp;
    return step(rainfall_metres, potential_et_metres, tf, rp);
}

lai_from_landcover_result lai_from_landcover_classes(
    const std::vector<int>&                class_values,
    const std::unordered_map<int, double>& lai_by_class,
    double                                 default_lai)
{
    const char* who = "lai_from_landcover_classes";

    if (!(std::isfinite(default_lai) && default_lai >= 0.0))
        throw std::invalid_argument(
            std::string(who) + ": default_lai must be finite and >= 0");

    for (const auto& kv : lai_by_class)
        if (!(std::isfinite(kv.second) && kv.second >= 0.0))
            throw std::invalid_argument(
                std::string(who) + ": class " + std::to_string(kv.first) +
                " has a negative or non-finite LAI");

    lai_from_landcover_result out;
    out.lai.resize(class_values.size());

    for (std::size_t i = 0; i < class_values.size(); ++i)
    {
        const auto it = lai_by_class.find(class_values[i]);
        if (it == lai_by_class.end())
        {
            out.lai[i] = default_lai;
            ++out.n_unmatched;
        }
        else
        {
            out.lai[i] = it->second;
        }
    }
    return out;
}

} // namespace lsdtt
