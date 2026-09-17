// lsdtt_hydro_roughness_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The roughness field: where n comes from, and how it varies in channels
//
// Phase 3 of the distributed-hydrology work. Two jobs:
//
//   1. Assemble the roughness field from the RIGHT SOURCE for each regime.
//   2. Give channels a roughness that varies the way bed material does.
//
// ROUGHNESS HAS THREE REGIMES AND TWO SOURCES
// -------------------------------------------
//   hillslope    n from land cover -- vegetation and surface treatment control
//                sheet-flow resistance
//   floodplain   n from land cover -- same, and this is where out-of-bank NFM
//                planting acts
//   channel      n from GRAIN SIZE -- in-channel resistance is set by the bed
//                material, not by what grows on the bank
//
// Note honestly that this is a three-way partition with only TWO distinct
// sources today: floodplain and hillslope both draw from land cover. The
// distinction is carried anyway because it becomes load-bearing the moment
// compound-channel conveyance arrives (channel and floodplain conveying
// simultaneously at different roughnesses once flow goes out of bank), which is
// deliberately out of scope here. Until then, `surface_class` changes the
// reported statistics and nothing else, and the partition is effectively
// channel versus not-channel.
//
// WHY CHANNEL n IS NOT A LAND-COVER PROPERTY
// ------------------------------------------
// Because it is a grain-size property, and grain size fines downstream
// (Sternberg 1875). With Strickler's n ~ D50^(1/6) and D50 ~ A^-b over
// b ~ 0.3-0.7, channel n FALLS downstream as A^(-0.05..-0.12). That reinforces
// the "low-order channels are slower" behaviour of Phase 1 through a second,
// independent mechanism: headwater channels are both shallower AND rougher.
//
// It also means a land-cover scenario -- planting a field, aerating a soil --
// must leave channel n untouched. Changing what grows on the bank does not
// change the bed material. rigorous_roughness_partition asserts exactly that.
//
// WHERE MANNING STOPS WORKING
// ---------------------------
// At low relative submergence Manning is not merely imprecise, it is the wrong
// shape. In headwater channels h/D84 is often 1-3, where resistance rises
// steeply and no constant n describes it. Two modes address this:
//
//   jarrett        n = 0.39 S^0.38 R^(-0.16), fitted to high-gradient streams.
//                  Depth AND slope dependent, so n is not a constant field.
//   ferguson_vpe   abandons Manning's form entirely for the variable power
//                  equation, then reports the EQUIVALENT n that would give the
//                  same velocity. That keeps one currency -- everything
//                  downstream still consumes a roughness raster -- while the
//                  physics underneath is not Manning's.
//
// Both are depth-dependent, which is the same point Phase 1 makes about
// velocity: roughness alone is not a resistance model.
//
// References
// ----------
//   Strickler A. (1923); Sternberg H. (1875)
//   Jarrett R.D. (1984) ASCE J. Hydraul. Eng. 110(11):1519-1539
//   Ferguson R. (2007) WRR 43:W05427 -- the variable power equation
//   Rickenmann D. & Recking A. (2011) WRR 47:W07538 -- VPE validated on ~2500 measurements
//   Leopold L.B. & Maddock T. (1953) USGS Prof. Paper 252
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// How channel roughness varies. See the header note on each.
// ---------------------------------------------------------------------------
enum class channel_roughness_mode
{
    constant,        // one value everywhere in channel
    drainage_area,   // n = n_ref (A/A_ref)^(-e_A)        -- downstream fining
    ksn,             // n = n_ref (ksn/ksn_ref)^(e_k)     -- steeper is coarser
    slope_area,      // as ksn, with ksn formed as S A^theta rather than supplied
    grain_size,      // Strickler, n = a_str D50^(1/6)    -- measured bed material
    jarrett,         // n = 0.39 S^0.38 R^(-0.16)         -- high-gradient streams
    ferguson_vpe     // equivalent n from the variable power equation
};

struct channel_roughness_options
{
    channel_roughness_mode mode = channel_roughness_mode::constant;

    // The universal value, and n_ref for every scaling mode.
    double channel_roughness = 0.035;

    // n ~ D50^(1/6) with D50 ~ A^-b, b ~ 0.3-0.7, gives 0.05-0.12.
    double area_exponent = 0.08;

    // D50 ~ S at fixed area and n ~ D50^(1/6) gives ~0.17.
    double ksn_exponent = 0.17;

    double reference_area_m2 = 1.0e6;

    // 0 -> use the median ksn over channel cells, so the field is centred on
    // n_ref rather than on an arbitrary absolute value.
    double reference_ksn = 0.0;

    // Concavity used to form ksn = S * A^theta in slope_area mode.
    double chi_m_over_n = 0.45;

    // Strickler coefficient, D50 in metres.
    double strickler_coefficient = 0.047;

    // Jarrett (1984) was fitted over S 0.002-0.04 and R 0.15-2.1 m. Outside
    // that the equation is extrapolation; cells outside are COUNTED and
    // reported rather than silently accepted.
    double jarrett_min_slope  = 0.002;
    double jarrett_max_slope  = 0.04;
    double jarrett_min_radius = 0.15;
    double jarrett_max_radius = 2.1;

    // Ferguson (2007) variable power equation coefficients.
    double ferguson_a1 = 6.5;
    double ferguson_a2 = 2.5;

    double min_slope = 1.0e-5;
    double n_min     = 0.005;
    double n_max     = 1.0;
};

// ---------------------------------------------------------------------------
// Rasters the roughness field may consume. Only what the chosen mode needs is
// required; compute_channel_roughness() says which in its exception text.
// ---------------------------------------------------------------------------
struct roughness_field_inputs
{
    const lsdtt_raster_cpu* drainage_area    = nullptr;  // m^2
    const lsdtt_raster_cpu* slope            = nullptr;  // m/m
    const lsdtt_raster_cpu* hydraulic_radius = nullptr;  // m; jarrett, ferguson_vpe
    const lsdtt_raster_cpu* d50              = nullptr;  // m; grain_size
    const lsdtt_raster_cpu* d84              = nullptr;  // m; ferguson_vpe
    const lsdtt_raster_cpu* ksn              = nullptr;  // ksn mode

    // Manning n from the land-cover table. Used on hillslope AND floodplain.
    const lsdtt_raster_cpu* landcover_roughness = nullptr;

    // Non-zero marks a channel cell.
    const lsdtt_raster_cpu* channel_mask = nullptr;

    // Optional; from lsdtt-fluvial-surfaces. Cells equal to
    // roughness_partition_options::floodplain_class are counted as floodplain.
    const lsdtt_raster_cpu* surface_class = nullptr;
};

// ---------------------------------------------------------------------------
// Closed forms, exposed for testing and for callers that work per cell.
// ---------------------------------------------------------------------------

// Strickler: n = a_str * D50^(1/6), D50 in metres.
[[nodiscard]] double strickler_roughness(double d50_metres,
                                         double strickler_coefficient) noexcept;

// Jarrett (1984): n = 0.39 * S^0.38 * R^(-0.16).
[[nodiscard]] double jarrett_roughness(double slope,
                                       double hydraulic_radius_metres) noexcept;

// Ferguson (2007) VPE. Returns (8/f)^(1/2):
//
//     (8/f)^(1/2) = a1 a2 (R/D84) / sqrt( a1^2 + a2^2 (R/D84)^(5/3) )
//
// At large R/D84 this tends to a1 (R/D84)^(1/6), i.e. Manning-Strickler; at
// small R/D84 it tends to a2 (R/D84), the roughness-layer regime where a
// constant Manning n is simply the wrong shape.
[[nodiscard]] double ferguson_sqrt_8_over_f(double hydraulic_radius_metres,
                                            double d84_metres,
                                            double a1, double a2) noexcept;

// VPE velocity: v = (8/f)^(1/2) * sqrt(g R S).
[[nodiscard]] double ferguson_vpe_velocity(double hydraulic_radius_metres,
                                           double d84_metres,
                                           double slope,
                                           double a1, double a2) noexcept;

// The Manning n that would reproduce the VPE velocity at this R and S:
//
//     n_eq = R^(2/3) S^(1/2) / v_vpe
//
// This is how a non-Manning law is carried in a Manning-shaped pipeline. It is
// exact by construction at the (R, S) it was evaluated at, and it is NOT a
// constant -- reusing it at a different depth reintroduces the error it exists
// to avoid.
[[nodiscard]] double ferguson_equivalent_manning_n(double hydraulic_radius_metres,
                                                   double d84_metres,
                                                   double slope,
                                                   double a1, double a2) noexcept;

// ---------------------------------------------------------------------------
// Channel roughness over the whole grid. Cells outside the channel mask are
// nodata; supply no mask to treat every valid cell as channel.
//
// Throws std::invalid_argument for a missing required input, a grid mismatch,
// or an out-of-range option.
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_raster_cpu compute_channel_roughness(
    const roughness_field_inputs&    inputs,
    const lsdtt_raster_meta&         meta,
    const channel_roughness_options& opts
);

// ---------------------------------------------------------------------------
// The assembled field.
// ---------------------------------------------------------------------------
struct roughness_partition_options
{
    channel_roughness_options channel;

    // Value in `surface_class` that means floodplain.
    double floodplain_class = 1.0;

    // Used where no land-cover raster is supplied, or where it is nodata.
    double landcover_roughness_default = 0.06;
};

struct roughness_partition_stats
{
    std::int64_t n_channel      = 0;
    std::int64_t n_floodplain   = 0;
    std::int64_t n_hillslope    = 0;
    std::int64_t n_unclassified = 0;   // valid in no input; left nodata

    // Cells outside Jarrett's calibration envelope, when that mode is used.
    std::int64_t n_outside_jarrett_range = 0;
};

struct roughness_field_result
{
    lsdtt_raster_cpu          roughness;
    roughness_partition_stats stats;
};

// Assemble: channel cells take the channel law, everything else takes land
// cover. Throws as compute_channel_roughness() does.
[[nodiscard]]
roughness_field_result compute_roughness_field(
    const roughness_field_inputs&      inputs,
    const lsdtt_raster_meta&           meta,
    const roughness_partition_options& opts
);

} // namespace lsdtt
