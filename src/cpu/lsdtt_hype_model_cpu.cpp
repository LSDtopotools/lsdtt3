// lsdtt_hype_model_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// The HYPE model: assembly of column and routing. See the header.

#include "lsdtt_hype_model_cpu.hpp"

#include <cstdint>
#include <exception>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

namespace lsdtt
{

double hype_river_volume_m3(const hype_river& r) noexcept
{
    // model_hype.f90:1424: damping box + water in transit in the queue
    return r.water_m3 + hype_river_in_transit_m3(r);
}

lsdtt_hype_model_cpu::lsdtt_hype_model_cpu(const hype_setup& setup, double dt_seconds,
                                           bool track_balance)
    : general_(setup.general), dt_seconds_(dt_seconds), track_balance_(track_balance)
{
    if (dt_seconds != 86400.0)
    {
        throw std::invalid_argument("hype model: dt_seconds must be 86400 (daily); HYPE's rate "
                                    "parameters are per timestep and the first version is daily only");
    }
    validate_hype_setup(setup);
    validate_hype_column_options(setup.general);

    std::map<int, std::size_t> sub_index;
    for (std::size_t s = 0; s < setup.subbasins.size(); ++s)
    {
        subids_.push_back(setup.subbasins[s].subid);
        sub_index[setup.subbasins[s].subid] = s;
    }
    std::map<int, const hype_class_definition*> classes;
    for (const auto& c : setup.classes) classes[c.class_id] = &c;
    std::map<int, const hype_landuse_parameters*> landuse;
    for (const auto& l : setup.landuse) landuse[l.code] = &l;
    std::map<int, const hype_soil_parameters*> soils;
    for (const auto& x : setup.soils) soils[x.code] = &x;

    // Units: sub-basin order, then class id order within a sub-basin.
    std::vector<std::map<int, const hype_subbasin_class*>> per_sub(setup.subbasins.size());
    for (const auto& sc : setup.subbasin_classes)
    {
        if (sc.fraction > 0.0) per_sub[sub_index.at(sc.subid)][sc.class_id] = &sc;
    }
    unit_offset_.push_back(0);
    for (std::size_t s = 0; s < setup.subbasins.size(); ++s)
    {
        for (const auto& [cid, sc] : per_sub[s])
        {
            const hype_class_definition& cd = *classes.at(cid);
            hype_model_unit u;
            u.subbasin  = s;
            u.class_id  = cid;
            u.fraction  = sc->fraction;
            u.constants = make_hype_unit_constants(setup.general, *landuse.at(cd.landuse),
                                                   *soils.at(cd.soil), cd,
                                                   setup.subbasins[s].slope_mean_percent, sc->dh_m);
            u.state     = initial_hype_unit_state(u.constants);
            units_.push_back(u);
        }
        unit_offset_.push_back(units_.size());
    }

    // Rivers and network.
    for (const auto& sb : setup.subbasins)
    {
        net_.local.push_back(make_hype_river(sb.loc_rivlen_m, general_.rivvel, general_.damp, dt_seconds_));
        net_.main.push_back(make_hype_river(sb.rivlen_m, general_.rivvel, general_.damp, dt_seconds_));
        net_.downstream.push_back(sb.maindown == 0 ? -1 : static_cast<int>(sub_index.at(sb.maindown)));
        net_.area_m2.push_back(sb.area_m2);
    }
    net_.order = hype_subbasin_downstream_order(setup.subbasins);
    unit_fluxes_.assign(units_.size(), hype_unit_fluxes{});   // stored by default
    unit_runoff_mm_.assign(units_.size(), 0.0);
    unit_prec_mm_.assign(units_.size(), 0.0);
    unit_evap_mm_.assign(units_.size(), 0.0);
    runoff_mm_.assign(setup.subbasins.size(), 0.0);
    flows_.local_outflow_m3s.assign(setup.subbasins.size(), 0.0);
    flows_.outflow_m3s.assign(setup.subbasins.size(), 0.0);

    if (track_balance_)
    {
        balance_.storage_start_m3 = storage_m3();
        balance_.storage_m3       = balance_.storage_start_m3;
    }
    else
    {
        balance_.storage_tracked  = false;
        balance_.storage_start_m3 = std::numeric_limits<double>::quiet_NaN();
        balance_.storage_m3       = std::numeric_limits<double>::quiet_NaN();
    }
}

double lsdtt_hype_model_cpu::storage_m3() const
{
    double v = 0.0;
    for (const auto& u : units_)
    {
        v += hype_unit_storage_mm(u.state) * 1.0e-3 * u.fraction * net_.area_m2[u.subbasin];
    }
    for (std::size_t s = 0; s < net_.local.size(); ++s)
    {
        v += hype_river_volume_m3(net_.local[s]) + hype_river_volume_m3(net_.main[s]);
    }
    return v;
}

void lsdtt_hype_model_cpu::set_store_unit_fluxes(bool store)
{
    store_unit_fluxes_ = store;
    if (store) unit_fluxes_.assign(units_.size(), hype_unit_fluxes{});
    else       std::vector<hype_unit_fluxes>().swap(unit_fluxes_);
}

void lsdtt_hype_model_cpu::step(double precipitation_mm, double temperature_c, int day_of_year,
                                double observed_pet_mm)
{
    hype_unit_forcing forcing;
    forcing.precipitation_mm = precipitation_mm;
    forcing.temperature_c    = temperature_c;
    forcing.day_of_year      = day_of_year;
    forcing.observed_pet_mm  = observed_pet_mm;
    // WP7 B3: the PET seasonal factor once per step, not once per unit.
    forcing.pet_seasonal_factor = hype_pet_seasonal_factor(general_.cevpam, general_.cevpph, day_of_year);

    // Land: every unit in parallel (each touches only its own state and
    // fluxes), then fixed-order aggregation per sub-basin in a separate serial
    // loop, so the sums do not depend on the thread count. An exception must
    // not leave the parallel region: record the one from the lowest-index
    // failing unit (the same unit at any thread count) and rethrow after.
    const std::int64_t n_units = static_cast<std::int64_t>(units_.size());
    std::exception_ptr first_error;
    std::int64_t       first_error_unit = n_units;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (std::int64_t k = 0; k < n_units; ++k)
    {
        hype_model_unit& u = units_[static_cast<std::size_t>(k)];
        try
        {
            if (static_cast<std::size_t>(k) >= inject_from_unit_)
            {
                throw std::runtime_error("hype model: injected test failure in unit " + std::to_string(k));
            }
            hype_unit_fluxes  local;
            hype_unit_fluxes& f = store_unit_fluxes_ ? unit_fluxes_[static_cast<std::size_t>(k)] : local;
            hype_unit_step(general_, u.constants, forcing, u.state, f);
            unit_runoff_mm_[static_cast<std::size_t>(k)] = f.runoff_mm;
            unit_prec_mm_[static_cast<std::size_t>(k)]   = f.class_precipitation_mm;
            unit_evap_mm_[static_cast<std::size_t>(k)]   = f.evap_mm;
        }
        catch (...)
        {
#ifdef _OPENMP
            #pragma omp critical (lsdtt_hype_model_unit_failure)
#endif
            {
                if (k < first_error_unit)
                {
                    first_error_unit = k;
                    first_error      = std::current_exception();
                }
            }
        }
    }
    if (first_error) std::rethrow_exception(first_error);

    double prec_m3 = 0.0, evap_m3 = 0.0;
    for (std::size_t s = 0; s + 1 < unit_offset_.size(); ++s)
    {
        double runoff = 0.0;
        for (std::size_t k = unit_offset_[s]; k < unit_offset_[s + 1]; ++k)
        {
            const auto& u = units_[k];
            runoff  += unit_runoff_mm_[k] * u.fraction;                       // model_hype.f90:2361
            const double area_m2 = u.fraction * net_.area_m2[s];
            prec_m3 += unit_prec_mm_[k] * 1.0e-3 * area_m2;
            evap_m3 += unit_evap_mm_[k] * 1.0e-3 * area_m2;
        }
        runoff_mm_[s] = runoff;
    }

    // Rivers.
    hype_route_network_step(net_, runoff_mm_, dt_seconds_, flows_);

    double out_m3 = 0.0;
    for (std::size_t s = 0; s < net_.downstream.size(); ++s)
    {
        if (net_.downstream[s] < 0) out_m3 += flows_.outflow_m3s[s] * dt_seconds_;
    }
    balance_.precipitation_m3 += prec_m3;
    balance_.evaporation_m3   += evap_m3;
    balance_.outflow_m3       += out_m3;
    if (track_balance_) balance_.storage_m3 = storage_m3();
}

} // namespace lsdtt
