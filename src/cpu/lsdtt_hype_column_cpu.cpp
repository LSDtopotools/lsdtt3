// lsdtt_hype_column_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// The HYPE land column. See the header for scope, units and sources. Every
// "file:line" reference is to HYPE 5.36.0 (SourceForge release_hype_5_36_0).

#include "lsdtt_hype_column_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace lsdtt
{

namespace
{

// HYPE's "zero" for real comparisons (general_wc.f90:29).
constexpr double k_realzero = 1.0e-37;

// Tolerance HYPE uses when deciding that a layer is saturated
// (soil_proc.f90:1628 and :2273).
constexpr double k_mindiff = 0.000002;

// HYPE's pi constant (modvar.f90:100). Kept, rather than M_PI, so the
// seasonal PET factor uses the same number HYPE does.
constexpr double k_hype_pi = 3.1415927;

// remove_water (general_wc.f90:141): subtract, snap a tiny negative to zero,
// and fail where HYPE would stop the run.
void remove_water(double& vol, double q, const char* where)
{
    if (vol - q > 0.0)
    {
        vol = vol - q;
    }
    else if (vol - q > -k_realzero)
    {
        vol = 0.0;
    }
    else
    {
        throw std::runtime_error(std::string("hype column: removing more water than present in ") +
                                 where + " (HYPE would stop: remove water routine Q>VOL)");
    }
}

// add_water (general_wc.f90:181) without substances.
void add_water(double& vol, double q) noexcept
{
    const double newvol = vol + q;
    vol = (newvol > k_realzero) ? newvol : 0.0;
}

// water_is_above_field_capacity / water_above_field_capacity /
// water_above_wilting_point / water_relative_porevolume
// (soil_proc.f90:2800-2846).
bool above_fc(const hype_unit_constants& c, int k, double w) noexcept
{
    return w - c.wp_mm[k] - c.fc_mm[k] > 0.0;
}
double water_above_fc(const hype_unit_constants& c, int k, double w) noexcept
{
    return above_fc(c, k, w) ? w - c.wp_mm[k] - c.fc_mm[k] : 0.0;
}
double water_above_wp(const hype_unit_constants& c, int k, double w) noexcept
{
    return (w - c.wp_mm[k] > 0.0) ? w - c.wp_mm[k] : 0.0;
}
double water_relative_porevolume(const hype_unit_constants& c, int k, double w) noexcept
{
    return (c.thick_m[k] > 0.0) ? w - c.pw_mm[k] : 0.0;
}

// ---------------------------------------------------------------------------
// Snowpack (calculate_snow, soil_proc.f90:486, with calculate_snowmelt
// option 0, soil_proc.f90:759). Options fixed at the first version's scope:
// no snow heat, no sublimation (p_snowevap = 0), no pack reset, fsceff = 0.
// ---------------------------------------------------------------------------
void snowpack_step(const hype_general_parameters& g, const hype_unit_constants& c,
                   double snowfall, double temp, hype_unit_state& s, hype_unit_fluxes& f)
{
    double snow    = s.snow_mm;
    double snowliq = s.snowliq_mm;
    const double whc = g.whcsnow;

    // effcov = 1 - fsceff*(1 - snowcover) (:572). fsceff = 0 here, so 1.
    const double effcov = 1.0;

    double snowice = std::max(0.0, snow - snowliq);   // (soil_proc.f90:588)

    // Potential melt, temperature index (soil_proc.f90:826). cmltcorr is a
    // regional parameter, absent here, so the correction factor is 1.
    double melt = 0.0;
    if (snow > 0.0 && temp >= c.landuse.ttmp)
    {
        melt = c.landuse.cmlt * (temp - c.landuse.ttmp);
        melt = melt * effcov;
    }
    if (melt > snowice) melt = snowice;               // (:602)

    // Snow evaporation off (p_snowevap = 0): evap = 0 (:620).
    double evap = 0.0;
    if (evap > snow) evap = snow;

    const double abla0 = melt + evap;                 // (:625)
    if (abla0 > 0.0)
    {
        if (snow < abla0)
        {
            melt = melt * (snow / abla0);
            evap = snow - melt;
        }
    }
    else
    {
        melt = 0.0;
        evap = 0.0;
    }

    // Snowfall joins the pack; rain on snow bypasses it (:636-642).
    if (snowfall > 0.0)
    {
        add_water(snow, snowfall);
        snowice = snowice + snowfall;
    }

    // Liquid holding capacity of the remaining ice (:647).
    double maxliq = 0.0;
    if ((1.0 - whc) > 0.0)
    {
        maxliq = whc * std::max(0.0, snowice - melt - evap) / (1.0 - whc);
    }

    double srunoff = 0.0;
    if (melt > 0.0 || snowliq > maxliq)               // (:653)
    {
        snowice = snowice - melt;
        snowliq = snowliq + melt;
        srunoff = std::max(0.0, snowliq - maxliq);
        snowliq = std::max(0.0, snowliq - srunoff);
        if (srunoff < snow && snow > 0.0)
        {
            remove_water(snow, srunoff, "snowpack");
        }
        else
        {
            snow    = 0.0;
            snowliq = 0.0;
        }
    }
    // Sublimation branch (:676-695) is inert: evap is 0.

    if (snow <= 0.0 && snowliq > 0.0) snowliq = 0.0;  // (:733)

    s.snow_mm    = snow;
    s.snowliq_mm = snowliq;
    f.melt_mm        = melt;
    f.snow_runoff_mm = srunoff;
}

// ---------------------------------------------------------------------------
// Diversion of gross infiltration (calculate_infiltration_flow_diversion,
// soil_proc.f90:1725, p_surfacerunoff = 0, no ice lenses).
// ---------------------------------------------------------------------------
void diversion(const hype_unit_constants& c, double ginfilt, double water1,
               hype_unit_fluxes& f)
{
    double infilt      = ginfilt;
    double surfaceflow = 0.0;
    double macroflow   = 0.0;
    if (ginfilt != 0.0)                                 // (:1784)
    {
        double macrate          = c.soil.macrate;
        const double inflowtres = c.soil.mactrinf;
        const double smtresmm   = (c.wp_mm[0] + c.fc_mm[0]) * c.soil.mactrsm;
        double srrate           = c.soil.srrate;
        if (macrate + srrate > 1.0)                     // (:1794)
        {
            const double help = macrate + srrate;
            macrate = macrate / help;
            srrate  = srrate / help;
        }
        const double avail = ginfilt - inflowtres;      // (:1801)
        if (avail > 0.0 && water1 > smtresmm)
        {
            macroflow   = macrate * avail;
            surfaceflow = srrate * avail;
        }
        infilt = ginfilt - macroflow - surfaceflow;     // (:1808)
    }
    f.infiltration_mm        = infilt;
    f.excess_infiltration_mm = surfaceflow;
    f.macropore_mm           = macroflow;
}

// ---------------------------------------------------------------------------
// Actual evapotranspiration (calculate_actual_soil_evapotranspiration,
// soil_proc.f90:1385). No soil-temperature reduction (tredA = 0), all water
// liquid, bare fraction 1 (snow sublimation off makes effcov 0,
// soilmodel0.f90:329).
// ---------------------------------------------------------------------------
void evapotranspiration(const hype_general_parameters& g, const hype_unit_constants& c,
                        double temp, double epot, const std::array<double, 3>& water_snapshot,
                        std::array<double, 3>& water, hype_unit_fluxes& f)
{
    f.soil_evap_mm = {0.0, 0.0};
    if (epot <= 0.0) return;                            // (:1433)
    if (!(temp > c.landuse.ttmp)) return;               // (:1451)
    const double barefrac = 1.0;
    for (int k = 0; k < k_hype_evap_layers; ++k)
    {
        if (!(c.thick_m[k] > 0.0)) continue;
        const double efwater = water_snapshot[k] - c.wp_mm[k];
        if (efwater <= 0.0) continue;
        double evapsl = epot * c.epot_frac[k];          // soiltemp_factor = 1
        if (efwater < g.lp * c.fc_mm[k]) evapsl = evapsl * efwater / (g.lp * c.fc_mm[k]);
        if (evapsl + k_realzero > water[k] - c.wp_mm[k]) evapsl = water[k] - c.wp_mm[k];
        evapsl = evapsl * std::min(1.0, barefrac);
        evapsl = evapsl * 1.0;                          // liquid fraction
        if (evapsl > 0.0)
        {
            if (evapsl + k_realzero < water[k]) remove_water(water[k], evapsl, "soil layer (evapotranspiration)");
            else water[k] = 0.0;
        }
        f.soil_evap_mm[k] = evapsl;
    }
}

// ---------------------------------------------------------------------------
// Soil-layer runoff (calculate_soil_runoff, soil_proc.f90:1594, and
// calculate_pressurelevel_soillayer_runoff, :1686). All water liquid.
// `water_snapshot` is the soil water before evapotranspiration
// (soilmodel0.f90:327); `water` the current state.
// ---------------------------------------------------------------------------
void soil_runoff(const hype_unit_constants& c, const std::array<double, 3>& water_snapshot,
                 std::array<double, 3>& water, hype_unit_fluxes& f)
{
    f.soil_runoff_mm = {0.0, 0.0, 0.0};
    const double ddepth = c.stream_depth_m;
    if (ddepth <= 0.0) return;                          // (:1635)

    std::array<double, 3> delta{}, avail{}, epwater{}, ep{};
    for (int k = 0; k < k_hype_max_soil_layers; ++k)    // (:1645-1648)
    {
        ep[k]      = c.ep_mm[k];
        epwater[k] = water_above_fc(c, k, water_snapshot[k]);
        avail[k]   = std::min(water_above_fc(c, k, water[k]), water_above_wp(c, k, water[k]));
        if (ep[k] > 0.0) delta[k] = epwater[k] / ep[k] * c.thick_m[k];
    }
    for (int k = 0; k < k_hype_max_soil_layers; ++k)
    {
        if (!(c.thick_m[k] > 0.0)) continue;
        if (ddepth < c.bottom_m[k] - c.thick_m[k]) continue;   // (:1654)
        double deltah = delta[k];
        bool includedd = false;
        if (ddepth <= c.bottom_m[k]) includedd = true;
        if (k < k_hype_max_soil_layers - 1)
        {
            if (c.thick_m[k + 1] == 0.0) includedd = true;
        }
        else
        {
            includedd = true;
        }
        if (includedd) deltah = deltah + (ddepth - c.bottom_m[k]);   // (:1664)
        for (int k2 = k - 1; k2 >= 0; --k2)             // (:1666)
        {
            if (epwater[k2 + 1] - ep[k2 + 1] >= 0.0 - k_mindiff * ep[k2 + 1]) deltah = deltah + delta[k2];
            else break;
        }
        // calculate_pressurelevel_soillayer_runoff (:1700-1711)
        double runoff = 0.0;
        if (avail[k] > 0.0 && deltah > 0.0 && c.thick_m[k] > 0.0)
        {
            runoff = c.recession[k] * deltah / c.thick_m[k] * ep[k];
            if (runoff >= avail[k]) runoff = avail[k];
            remove_water(water[k], runoff, "soil layer (soil runoff)");
        }
        f.soil_runoff_mm[k] = runoff;
    }
}

// ---------------------------------------------------------------------------
// Macropore water to the layer holding the groundwater table
// (add_macropore_flow, soil_proc.f90:2157).
// ---------------------------------------------------------------------------
void add_macropore_flow(const hype_unit_constants& c, double macroflow,
                        std::array<double, 3>& water, std::array<double, 3>& to_layer)
{
    to_layer = {0.0, 0.0, 0.0};
    if (macroflow <= 0.0) return;
    const double gwat = hype_groundwater_table(c, water);
    if (-gwat > c.bottom_m[1])                          // table in layer 3 (:2188)
    {
        double fill = 0.0;
        const double newsoil = water[2] + macroflow;
        if (newsoil > c.pw_mm[2])
        {
            fill = newsoil - c.pw_mm[2];
            water[2] = c.pw_mm[2];
            to_layer[2] = macroflow - fill;
        }
        else
        {
            water[2] = newsoil;
            to_layer[2] = macroflow;
        }
        if (fill > 0.0)
        {
            double fill2 = 0.0;
            const double newsoil2 = water[1] + fill;
            if (newsoil2 > c.pw_mm[1])
            {
                fill2 = newsoil2 - c.pw_mm[1];
                water[1] = c.pw_mm[1];
                to_layer[1] = fill - fill2;
            }
            else
            {
                water[1] = newsoil2;
                to_layer[1] = fill;
            }
            if (fill2 > 0.0)
            {
                water[0] = water[0] + fill2;
                to_layer[0] = fill2;
            }
        }
    }
    else if (-gwat > c.bottom_m[0])                     // table in layer 2 (:2226)
    {
        double fill = 0.0;
        const double newsoil = water[1] + macroflow;
        if (newsoil > c.pw_mm[1])
        {
            fill = newsoil - c.pw_mm[1];
            water[1] = c.pw_mm[1];
            to_layer[1] = macroflow - fill;
        }
        else
        {
            water[1] = newsoil;
            to_layer[1] = macroflow;
        }
        if (fill > 0.0)
        {
            add_water(water[0], fill);
            to_layer[0] = fill;
        }
    }
    else                                                // table in layer 1 (:2248)
    {
        add_water(water[0], macroflow);
        to_layer[0] = macroflow;
    }
}

} // namespace

// ---------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------

void validate_hype_column_options(const hype_general_parameters& general)
{
    if (general.fsceff != 0.0)
    {
        throw std::invalid_argument(
            "hype column: fsceff = " + std::to_string(general.fsceff) +
            " but the first version supports only fsceff = 0 (no fractional snow cover)");
    }
}

hype_unit_constants make_hype_unit_constants(
    const hype_general_parameters& general,
    const hype_landuse_parameters& landuse,
    const hype_soil_parameters&    soil,
    const hype_class_definition&   class_def,
    double                         slope_mean_percent,
    double                         dh_m)
{
    validate_hype_column_options(general);
    if (class_def.n_layers < 1 || class_def.n_layers > k_hype_max_soil_layers)
    {
        throw std::invalid_argument("hype column: class " + std::to_string(class_def.class_id) +
                                    " has n_layers outside 1..3");
    }

    hype_unit_constants c;
    c.n_layers       = class_def.n_layers;
    c.landuse.ttmp   = landuse.ttmp;
    c.landuse.cmlt   = landuse.cmlt;
    c.landuse.cevp   = landuse.cevp;
    c.landuse.srrcs  = landuse.srrcs;
    c.landuse.pcluse = landuse.pcluse;
    c.soil.mperc1    = soil.mperc1;
    c.soil.mperc2    = soil.mperc2;
    c.soil.macrate   = soil.macrate;
    c.soil.srrate    = soil.srrate;
    c.soil.mactrinf  = soil.mactrinf;
    c.soil.mactrsm   = soil.mactrsm;
    c.stream_depth_m = class_def.stream_depth_m;
    c.dh_m           = dh_m;

    // Layer geometry (set_class_current_soildepth, modvar.f90:1555-1566).
    for (int k = 0; k < c.n_layers; ++k)
    {
        c.bottom_m[k] = class_def.layer_bottom_m[k];
        c.thick_m[k]  = (k == 0) ? c.bottom_m[0] : c.bottom_m[k] - c.bottom_m[k - 1];
    }
    for (int k = c.n_layers; k < k_hype_max_soil_layers; ++k)
    {
        c.bottom_m[k] = c.bottom_m[c.n_layers - 1];
        c.thick_m[k]  = 0.0;
    }

    // Pore volumes (calculate_porevolumes, soil_proc.f90:284). Per-layer
    // overrides (wcfc1..3 etc.) are not supported in the first version.
    for (int k = 0; k < k_hype_max_soil_layers; ++k)
    {
        const double thick1000 = c.thick_m[k] * 1000.0;
        c.wp_mm[k] = soil.wcwp * thick1000;
        c.fc_mm[k] = soil.wcfc * thick1000;
        c.ep_mm[k] = soil.wcep * thick1000;
        c.pw_mm[k] = c.wp_mm[k] + c.fc_mm[k] + c.ep_mm[k];
    }

    // PET split between the evaporating layers (soil_proc.f90:218-228).
    {
        const double coeff = general.epotdist;
        c.epot_frac[0] = c.thick_m[0] * std::exp(-coeff * c.bottom_m[0] / 2.0);
        double sums = c.epot_frac[0];
        for (int k = 1; k < k_hype_evap_layers; ++k)
        {
            c.epot_frac[k] = c.thick_m[k] *
                std::exp(-coeff * (c.bottom_m[k - 1] + (c.bottom_m[k] - c.bottom_m[k - 1]) / 2.0));
            sums = sums + c.epot_frac[k];
        }
        for (int k = 0; k < k_hype_evap_layers; ++k) c.epot_frac[k] = c.epot_frac[k] / sums;
    }

    // Soil-runoff recession coefficients (soil_proc.f90:337). The sub-basin
    // correction basincorr is a regional parameter, absent here, so 1.
    {
        const double rc2 = general.rrcs3;
        double rc0 = soil.rrcs1;
        if (rc0 > 1.0) rc0 = 1.0;
        rc0 = rc0 + rc2 * slope_mean_percent;            // (:373)
        if (rc0 > 1.0) rc0 = 1.0;
        double rc1 = soil.rrcs2;
        if (rc1 > 1.0) rc1 = 1.0;
        if (rc1 == 0.0) rc1 = rc0;
        const double b = std::log(rc0 / rc1) /
            ((c.bottom_m[2] - c.thick_m[2] / 2.0) - c.thick_m[0] / 2.0);   // (:378)
        c.recession[0] = rc0;
        c.recession[2] = rc1;
        double depth = 0.0;
        for (int k = 1; k < k_hype_max_soil_layers - 1; ++k)
        {
            depth = depth + c.thick_m[k - 1] / 2.0 + c.thick_m[k] / 2.0;
            c.recession[k] = rc0 * std::exp(-b * depth);   // (:385)
        }
    }
    return c;
}

hype_unit_state initial_hype_unit_state(const hype_unit_constants& c)
{
    hype_unit_state s;
    for (int k = 0; k < k_hype_max_soil_layers; ++k) s.water_mm[k] = c.wp_mm[k] + c.fc_mm[k];
    return s;
}

double hype_unit_storage_mm(const hype_unit_state& s) noexcept
{
    return s.snow_mm + s.water_mm[0] + s.water_mm[1] + s.water_mm[2];
}

double hype_rain_fraction(double precipitation_mm, double temperature_c,
                          double ttmp, double ttpd, double ttpi) noexcept
{
    // get_rainfall_fraction / get_snowfall_fraction (atm_proc.f90:476, :525)
    const double dtt = ttpi;
    const double tt  = ttmp + ttpd;
    double arain = 0.0;
    if (precipitation_mm > 0.0)
    {
        if (temperature_c >= tt + dtt)
        {
            arain = 1.0;
        }
        else if (temperature_c > tt - dtt && temperature_c < tt + dtt)
        {
            arain = (temperature_c - (tt - dtt)) / (2 * dtt);
        }
    }
    return arain;
}

double hype_pet_seasonal_factor(double cevpam, double cevpph, int day_of_year) noexcept
{
    // atm_proc.f90:2068. Daily step: tsofday = timesteps_per_day = 1, so the
    // day term is (dayno - 1 + 1/1 - cevpph).
    return 1.0 + cevpam * std::sin(2.0 * k_hype_pi *
        (static_cast<double>(day_of_year) - 1.0 + 1.0 - cevpph) / 365.0);
}

double hype_potential_evaporation_with_factor(double temperature_c, double ttmp, double cevp,
                                              double seasonal_factor) noexcept
{
    const double ce = cevp * seasonal_factor;
    // Option 0/1 (:2094/:2100); cevpcorr is regional and absent, so x 1.
    return (temperature_c > ttmp) ? ce * (temperature_c - ttmp) : 0.0;
}

double hype_potential_evaporation(double temperature_c, double ttmp, double cevp,
                                  double cevpam, double cevpph, int day_of_year) noexcept
{
    return hype_potential_evaporation_with_factor(
        temperature_c, ttmp, cevp, hype_pet_seasonal_factor(cevpam, cevpph, day_of_year));
}

void hype_percolation(const hype_unit_constants& c, std::array<double, 3>& water,
                      std::array<double, 2>& percflow)
{
    // percolation, soil_proc.f90:2064. All water liquid.
    if (!(c.thick_m[1] > 0.0)) return;                  // (:2103)
    const std::array<double, 2> first = percflow;
    const double maxperc1 = std::min(water_above_fc(c, 0, water[0]), c.soil.mperc1 - first[0]);
    const double maxperc2 = std::min(-water_relative_porevolume(c, 2, water[2]), c.soil.mperc2 - first[1]);
    const double water2help  = water[1] - c.wp_mm[1] - c.fc_mm[1];
    const double water2help2 = water_relative_porevolume(c, 1, water[1]);
    const double perc2 = std::min(std::max(water2help + maxperc1, 0.0), maxperc2);   // (:2115)
    const double perc1 = std::min(maxperc1, std::max(perc2 - water2help2, 0.0));       // (:2116)
    if (perc1 > 0.0)
    {
        remove_water(water[0], perc1, "soil layer 1 (percolation)");
        add_water(water[1], perc1);
        percflow[0] = perc1 + first[0];
    }
    if (perc2 > 0.0)
    {
        remove_water(water[1], perc2, "soil layer 2 (percolation)");
        add_water(water[2], perc2);
        percflow[1] = perc2 + first[1];
    }
}

double hype_groundwater_table(const hype_unit_constants& c, const std::array<double, 3>& water)
{
    // calculate_groundwater_table, soil_proc.f90:2260
    double gwatk = 0.0;
    for (int k = k_hype_max_soil_layers - 1; k >= 0; --k)
    {
        if (!(c.thick_m[k] > 0.0)) continue;
        if (above_fc(c, k, water[k]))
        {
            gwatk = water_relative_porevolume(c, k, water[k]) / c.ep_mm[k] * c.thick_m[k];
            if (k == 0)
            {
                if (gwatk > 0.0) gwatk = water_relative_porevolume(c, k, water[k]) * 0.001;
            }
            else
            {
                gwatk = gwatk - c.bottom_m[k - 1];
            }
            if (water_relative_porevolume(c, k, water[k]) >= 0.0 - k_mindiff * c.ep_mm[k]) continue;
            break;
        }
        gwatk = -c.bottom_m[k];
        break;
    }
    return gwatk;
}

void hype_unit_step(const hype_general_parameters& g, const hype_unit_constants& c,
                    const hype_unit_forcing& forcing, hype_unit_state& s,
                    hype_unit_fluxes& f)
{
    f = hype_unit_fluxes{};

    // --- Class forcing (calculate_class_atmospheric_forcing, atm_proc.f90:192).
    // Temperature: lapse rate on the class elevation deviation (:332). The
    // sub-basin corrections (tcobselev, tcelevadd, mlapse) are absent: identity.
    double temp = forcing.temperature_c;
    if (g.tcalt != 0.0) temp = temp - g.tcalt * c.dh_m * 0.01;
    // Precipitation: the elevation correction basinpreccorr is 1 with its
    // parameters absent (atm_proc.f90:164); interception pcluse (:277-279).
    double prec = forcing.precipitation_mm * 1.0;
    f.interception_mm = prec * c.landuse.pcluse;
    prec = prec * (1.0 - c.landuse.pcluse);
    f.class_temperature_c    = temp;
    f.class_precipitation_mm = prec;

    // --- Potential evaporation (soilmodel0.f90:281). Snow PET (fepotsnow) is
    // weighted by effcov, which is 0 with sublimation off, so epot is unchanged
    // (soilmodel0.f90:338).
    // Observed PET replaces the calculated value (option 0, atm_proc.f90:2091).
    const double seasonal = std::isfinite(forcing.pet_seasonal_factor)
        ? forcing.pet_seasonal_factor
        : hype_pet_seasonal_factor(g.cevpam, g.cevpph, forcing.day_of_year);
    const double epot = std::isfinite(forcing.observed_pet_mm)
        ? forcing.observed_pet_mm
        : hype_potential_evaporation_with_factor(temp, c.landuse.ttmp, c.landuse.cevp, seasonal);
    f.epot_mm = epot;

    // --- Rain/snow split (:284) and snowpack (:285).
    const double arain = hype_rain_fraction(prec, temp, c.landuse.ttmp, g.ttpd, g.ttpi);
    const double asnow = 1.0 - arain;
    f.rainfall_mm = prec * (1.0 - asnow);
    f.snowfall_mm = prec * asnow;
    snowpack_step(g, c, f.snowfall_mm, temp, s, f);

    // --- Gross infiltration (:290).
    const double ginfilt = f.rainfall_mm + f.snow_runoff_mm;
    f.gross_infiltration_mm = ginfilt;

    // --- Diversion (:311) and infiltration into layer 1 (:318).
    diversion(c, ginfilt, s.water_mm[0], f);
    if (f.infiltration_mm > 0.0) add_water(s.water_mm[0], f.infiltration_mm);   // add_infiltration (soil_proc.f90:1956)

    // --- Percolation, first call (:322).
    f.percolation_mm = {0.0, 0.0};
    hype_percolation(c, s.water_mm, f.percolation_mm);

    // --- Evapotranspiration (:328), from a snapshot that soil runoff also uses.
    const std::array<double, 3> soilwater = s.water_mm;   // (:327)
    evapotranspiration(g, c, temp, epot, soilwater, s.water_mm, f);
    f.evap_mm = f.soil_evap_mm[0] + f.soil_evap_mm[1] + 0.0;   // + snow evaporation (0)

    // --- Soil-layer runoff (:341). Tile drainage (:346) is off.
    soil_runoff(c, soilwater, s.water_mm, f);

    // --- Saturation-excess surface runoff (calculate_surface_runoff,
    // soil_proc.f90:1998; rrcscorr is regional and absent, so 1).
    {
        double sc = c.landuse.srrcs * 1.0;
        if (sc > 1.0) sc = 1.0;
        double satoverflow = std::max(sc * (s.water_mm[0] - c.pw_mm[0]), 0.0);
        if (satoverflow > 0.0) remove_water(s.water_mm[0], satoverflow, "soil layer 1 (surface runoff)");
        f.saturation_excess_mm = satoverflow;
    }

    // --- Macropore water (:431) and the second percolation (:437).
    add_macropore_flow(c, f.macropore_mm, s.water_mm, f.macropore_to_layer_mm);
    hype_percolation(c, s.water_mm, f.percolation_mm);

    // --- Class runoff and HYPE's output conventions (model_hype.f90).
    const double totalsurfaceflow = f.saturation_excess_mm + f.excess_infiltration_mm;
    f.runoff_mm = f.soil_runoff_mm[0] + f.soil_runoff_mm[1] + f.soil_runoff_mm[2] +
                  totalsurfaceflow + 0.0;                                       // (:2305), tile 0
    const double macro_sum = f.macropore_to_layer_mm[0] + f.macropore_to_layer_mm[1] +
                             f.macropore_to_layer_mm[2];
    const double snowmeltfrac = (ginfilt > 0.0) ? f.snow_runoff_mm / ginfilt : 0.0;   // soilmodel0.f90:557
    f.out_melt_mm = snowmeltfrac * (f.infiltration_mm + f.excess_infiltration_mm + macro_sum);  // (:2191)
    f.out_infi_mm = f.infiltration_mm + macro_sum;                                    // (:2244)
    const double standingwater = (s.water_mm[0] > c.pw_mm[0]) ? s.water_mm[0] - c.pw_mm[0] : 0.0;  // (:2265)
    f.out_sml_mm  = {s.water_mm[0] - standingwater, s.water_mm[1], s.water_mm[2]};
    f.out_snwc_mm = s.snowliq_mm;
}

} // namespace lsdtt
