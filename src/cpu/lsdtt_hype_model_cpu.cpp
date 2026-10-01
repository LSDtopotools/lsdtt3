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

#include <map>
#include <stdexcept>
#include <string>

namespace lsdtt
{

double hype_river_volume_m3(const hype_river& r) noexcept
{
    // model_hype.f90:1424: water + SUM(qqueue(1:ttstep)) + qqueue(ttstep+1)*ttpart
    double v = r.water_m3;
    if (r.length_m > 0.0 && r.transtime > 0.0)
    {
        const std::size_t y = static_cast<std::size_t>(r.ttstep);
        for (std::size_t k = 1; k <= y; ++k) v += r.queue[k];
        v += r.queue[y + 1] * r.ttpart;
    }
    return v;
}

lsdtt_hype_model_cpu::lsdtt_hype_model_cpu(const hype_setup& setup, double dt_seconds)
    : general_(setup.general), dt_seconds_(dt_seconds)
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
    runoff_mm_.assign(setup.subbasins.size(), 0.0);
    flows_.local_outflow_m3s.assign(setup.subbasins.size(), 0.0);
    flows_.outflow_m3s.assign(setup.subbasins.size(), 0.0);

    balance_.storage_start_m3 = storage_m3();
    balance_.storage_m3       = balance_.storage_start_m3;
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

void lsdtt_hype_model_cpu::step(double precipitation_mm, double temperature_c, int day_of_year,
                                double observed_pet_mm)
{
    hype_unit_forcing forcing;
    forcing.precipitation_mm = precipitation_mm;
    forcing.temperature_c    = temperature_c;
    forcing.day_of_year      = day_of_year;
    forcing.observed_pet_mm  = observed_pet_mm;

    // Land: every unit, then fixed-order aggregation per sub-basin.
    for (auto& u : units_) hype_unit_step(general_, u.constants, forcing, u.state, u.fluxes);

    double prec_m3 = 0.0, evap_m3 = 0.0;
    for (std::size_t s = 0; s + 1 < unit_offset_.size(); ++s)
    {
        double runoff = 0.0;
        for (std::size_t k = unit_offset_[s]; k < unit_offset_[s + 1]; ++k)
        {
            const auto& u = units_[k];
            runoff  += u.fluxes.runoff_mm * u.fraction;                       // model_hype.f90:2361
            const double area_m2 = u.fraction * net_.area_m2[s];
            prec_m3 += u.fluxes.class_precipitation_mm * 1.0e-3 * area_m2;
            evap_m3 += u.fluxes.evap_mm * 1.0e-3 * area_m2;
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
    balance_.storage_m3        = storage_m3();
}

} // namespace lsdtt
