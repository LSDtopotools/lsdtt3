// lsdtt_hype_column_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The HYPE land column for one (sub-basin, class) unit, one daily timestep
//
// A reimplementation of the land processes of HYPE (SMHI; Lindström et al.
// 2010, Hydrology Research 41:295-319) with HYPE's DEFAULT model options, for
// the HYPE column of lsdtt-distributed-hydrology. Written from the equations
// in HYPE 5.36.0's Fortran (LGPLv3); no Fortran is translated line for line,
// but the ORDER of operations, the caps and the clips are reproduced exactly,
// because this is Stage A of the port: it must reproduce HYPE, and is checked
// against HYPE's own output (tests/expected/hype, unit_hype_column).
//
// Every process function below cites the HYPE routine it reproduces as
// file:line in HYPE 5.36.0. The sequence in hype_unit_step() is
// soilmodel_0 (soilmodel0.f90:104), preceded by the class forcing of
// calculate_class_atmospheric_forcing (atm_proc.f90:192):
//
//   class temperature (lapse rate) and precipitation (pcluse interception)
//   potential evaporation                       (soilmodel0.f90:281)
//   rain/snow split, snowpack                   (:284, :285)
//   gross infiltration = rain + snow runoff     (:290)
//   diversion to surface runoff / macropores    (:311), infiltration (:318)
//   percolation, first call                     (:322)
//   actual evapotranspiration                   (:328)
//   soil-layer runoff                           (:341)
//   [tile drainage: off in the first version]   (:346)
//   saturation-excess surface runoff            (:397)
//   macropore water to the groundwater layer    (:431)
//   percolation, second call                    (:437)
//
// UNITS ARE HYPE'S: mm and mm per timestep for water, degrees C, metres for
// depths. The conversion from lsdtt3's metres happens at the model boundary,
// not here.
//
// FIRST-VERSION RESTRICTIONS, enforced by
// validate_hype_column_options(): daily timestep; fsceff = 0 (no fractional
// snow cover); no snow sublimation; no tile drainage; no frozen soil; no
// regional parameter corrections. Under those options several HYPE terms are
// identities and are written as such, with the reason in a comment.
//
// SURPRISES IN HYPE THAT ARE REPRODUCED ON PURPOSE
//   * Rain on snow bypasses the snowpack (soil_proc.f90:636).
//   * macrate and srrate are rescaled to sum to 1 when their sum exceeds 1
//     (soil_proc.f90:1794).
//   * Percolation runs twice per step; the second call is capped by what the
//     first already moved ("limited to same maxperc", soilmodel0.f90:436).
//   * Soil runoff takes its pressure level from the soil water BEFORE
//     evapotranspiration but its available water from the state AFTER it
//     (soilmodel0.f90:327, soil_proc.f90:1647 and :1646).
//   * HYPE's class output "melt" is the liquid water LEAVING the snowpack, not
//     the melt rate (model_hype.f90:2191); "sml1" excludes water above the
//     pore volume ("standing water", model_hype.f90:2265).
//
// Precision: double throughout (HYPE runs in single precision; the oracle's
// README records the resulting floor).
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstddef>
#include <limits>

#include "lsdtt_hype_parameters.hpp"

namespace lsdtt
{

// Number of soil layers HYPE allows, and the number that evaporate under the
// default soil evaporation option (p_evapmodel = 0, soil_proc.f90:212).
inline constexpr int k_hype_max_soil_layers = 3;
inline constexpr int k_hype_evap_layers     = 2;

// ---------------------------------------------------------------------------
// Constants of one unit, computed once (HYPE computes them at initialisation).
// ---------------------------------------------------------------------------
struct hype_unit_constants
{
    int n_layers = 0;

    // Layer geometry [m]. Unused layers have zero thickness and a bottom equal
    // to the last used layer's bottom (modvar.f90:1564-1566).
    std::array<double, 3> thick_m{};
    std::array<double, 3> bottom_m{};

    // Stacked pore volumes [mm] (calculate_porevolumes, soil_proc.f90:284).
    std::array<double, 3> wp_mm{};   // below wilting point
    std::array<double, 3> fc_mm{};   // wilting point to field capacity
    std::array<double, 3> ep_mm{};   // above field capacity (drainable)
    std::array<double, 3> pw_mm{};   // wp + fc + ep

    // Split of potential evaporation between the evaporating layers
    // (soil_proc.f90:218-228).
    std::array<double, 2> epot_frac{};

    // Soil-runoff recession coefficient per layer [1/timestep]
    // (calculate_soil_runoff_recession_coefficients, soil_proc.f90:337).
    std::array<double, 3> recession{};

    double stream_depth_m = 0.0;
    double dh_m           = 0.0;   // class elevation deviation from sub-basin mean

    // The land-use and soil parameters the column reads every step, and only
    // those: the full parameter rows are not copied into every unit, to keep
    // units small.
    struct landuse_part { double ttmp = 0, cmlt = 0, cevp = 0, srrcs = 0, pcluse = 0; };
    struct soil_part    { double mperc1 = 0, mperc2 = 0, macrate = 0, srrate = 0, mactrinf = 0, mactrsm = 0; };
    landuse_part landuse;
    soil_part    soil;
};

// slope_mean_percent is the sub-basin's mean slope in percent (100 x gradient).
[[nodiscard]] hype_unit_constants make_hype_unit_constants(
    const hype_general_parameters& general,
    const hype_landuse_parameters& landuse,
    const hype_soil_parameters&    soil,
    const hype_class_definition&   class_def,
    double                         slope_mean_percent,
    double                         dh_m);

// Throws std::invalid_argument for general parameters outside the first
// version's scope: fsceff != 0.
void validate_hype_column_options(const hype_general_parameters& general);

// ---------------------------------------------------------------------------
// State, forcing and per-step fluxes of one unit.
// ---------------------------------------------------------------------------
struct hype_unit_state
{
    double snow_mm     = 0.0;   // snowpack, ice + liquid
    double snowliq_mm  = 0.0;   // liquid water held in the snowpack
    std::array<double, 3> water_mm{};   // soil water per layer (total, incl. wp)
};

// HYPE's initial state: no snow, every layer at field capacity (wp + fc)
// (initiate_soil_water_state, soil_proc.f90:97; soiliniwet off).
[[nodiscard]] hype_unit_state initial_hype_unit_state(const hype_unit_constants& c);

struct hype_unit_forcing
{
    double precipitation_mm = 0.0;  // sub-basin precipitation for the step
    double temperature_c    = 0.0;  // sub-basin air temperature
    int    day_of_year      = 1;    // 1..366, for the PET seasonal factor
    // Seasonal PET factor 1 + cevpam*sin(...) for this step. It
    // depends only on the day and two general parameters, so a caller
    // stepping many units computes it once with hype_pet_seasonal_factor()
    // and passes it in. NaN means "compute it here".
    double pet_seasonal_factor = std::numeric_limits<double>::quiet_NaN();
    // Observed potential evaporation [mm/step]. When supplied (finite), it
    // replaces the calculated value, as HYPE's PET option 0 does with an
    // observed series (atm_proc.f90:2091-2092). NaN means "not supplied".
    double observed_pet_mm  = std::numeric_limits<double>::quiet_NaN();
};

struct hype_unit_fluxes
{
    // Class forcing
    double class_temperature_c   = 0.0;   // HYPE ctmp
    double class_precipitation_mm= 0.0;   // HYPE cprc (after interception)
    double interception_mm       = 0.0;   // pcluse loss

    // Snow
    double rainfall_mm   = 0.0;
    double snowfall_mm   = 0.0;
    double melt_mm       = 0.0;   // melt of the pack (internal)
    double snow_runoff_mm= 0.0;   // liquid leaving the pack

    // Evaporation
    double epot_mm = 0.0;                    // HYPE epot
    std::array<double, 2> soil_evap_mm{};    // per evaporating layer
    double evap_mm = 0.0;                    // HYPE evap (soil + snow; snow is 0 here)

    // Infiltration and diversion
    double gross_infiltration_mm = 0.0;      // rain + snow runoff
    double infiltration_mm       = 0.0;      // into layer 1 after diversion
    double excess_infiltration_mm= 0.0;      // HYPE ros2
    double macropore_mm          = 0.0;      // diverted to macropores
    std::array<double, 3> macropore_to_layer_mm{};

    // Percolation, cumulative over both calls: 1->2 and 2->3
    std::array<double, 2> percolation_mm{};

    // Runoff
    std::array<double, 3> soil_runoff_mm{};  // HYPE cro1..cro3
    double saturation_excess_mm = 0.0;       // HYPE ros1
    double runoff_mm            = 0.0;       // HYPE crun

    // HYPE class output conventions (see the header comment)
    double out_melt_mm = 0.0;                // HYPE melt
    double out_infi_mm = 0.0;                // HYPE infi
    std::array<double, 3> out_sml_mm{};      // HYPE sml1..sml3
    double out_snwc_mm = 0.0;                // HYPE snwc
};

// Advance one unit by one daily step. `state` is updated in place.
void hype_unit_step(const hype_general_parameters& general,
                    const hype_unit_constants&     c,
                    const hype_unit_forcing&       forcing,
                    hype_unit_state&               state,
                    hype_unit_fluxes&              fluxes);

// Storage of a unit [mm]: snow plus all soil water. For water-balance checks.
[[nodiscard]] double hype_unit_storage_mm(const hype_unit_state& s) noexcept;

// ---------------------------------------------------------------------------
// The individual processes, exposed for the analytic tests.
// ---------------------------------------------------------------------------

// Fraction of precipitation falling as rain (get_rainfall_fraction,
// atm_proc.f90:476, default option).
[[nodiscard]] double hype_rain_fraction(double precipitation_mm, double temperature_c,
                                        double ttmp, double ttpd, double ttpi) noexcept;

// Seasonal factor of potential evaporation, 1 + cevpam*sin(...)
// (atm_proc.f90:2068), daily steps.
[[nodiscard]] double hype_pet_seasonal_factor(double cevpam, double cevpph, int day_of_year) noexcept;

// Potential evaporation, option 0/1 (atm_proc.f90:2068, :2094).
[[nodiscard]] double hype_potential_evaporation(double temperature_c, double ttmp,
                                                double cevp, double cevpam, double cevpph,
                                                int day_of_year) noexcept;

// As above, with the seasonal factor already computed. Gives
// bitwise the same value as hype_potential_evaporation.
[[nodiscard]] double hype_potential_evaporation_with_factor(double temperature_c, double ttmp,
                                                            double cevp, double seasonal_factor) noexcept;

// Percolation 1->2 and 2->3 (soil_proc.f90:2064). `percolation_mm` carries
// what earlier calls in the same step moved, and is increased.
void hype_percolation(const hype_unit_constants& c, std::array<double, 3>& water_mm,
                      std::array<double, 2>& percolation_mm);

// Groundwater table [m, negative below the surface]
// (calculate_groundwater_table, soil_proc.f90:2260).
[[nodiscard]] double hype_groundwater_table(const hype_unit_constants& c,
                                            const std::array<double, 3>& water_mm);

} // namespace lsdtt
