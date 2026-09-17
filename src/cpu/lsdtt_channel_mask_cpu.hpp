// lsdtt_channel_mask_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Where the channels are -- and why that is a refinement, not a prerequisite
//
// Phase 6. Five ways to say which cells are channel, from "do not say at all"
// to "here is the Ordnance Survey river network".
//
// WHAT THE MASK IS ACTUALLY FOR
// -----------------------------
// Not for deciding where water goes. Saturation excess is generated on ANY
// cell, the saturated area expands and contracts through a storm, and the
// water routes over everything. The mask exists to choose WHICH RESISTANCE
// CLOSURE applies: sheet flow (depth in millimetres, n ~ 0.1-0.5) or channel
// flow (depth in metres, n ~ 0.035). That is a 2-3 order-of-magnitude contrast
// in velocity.
//
// WHICH IS EXACTLY WHY A HARD THRESHOLD IS THE WRONG INSTRUMENT
// -------------------------------------------------------------
// A three-order velocity step across one pixel boundary is numerically nasty
// and physically false: rills, ephemeral swales and zero-order hollows are a
// continuum, not a step. So the work is split in two:
//
//   * GEOMETRY transitions smoothly, with no threshold anywhere. The velocity
//     closure uses R = max(h_sheet, a_R (A/A_ref)^c_R), and a "channel" is
//     simply where hydraulic-geometry depth has overtaken sheet-flow depth.
//     Continuous, differentiable, no parameter to choose.
//
//   * The MASK only overrides ROUGHNESS (and, later, imposes a minimum depth).
//
// So the model runs correctly with NO mask at all and gets better when one is
// supplied. `mode = none` is the default and returns the emergent mask, which
// is a diagnostic of where the geometry has already put the channel head -- not
// an input to it.
//
// THE LADDER
// ----------
//   none                   emergent from the geometry. No parameter.
//   area_threshold         A >= threshold. One parameter, the simple option.
//   area_slope_threshold   A S^theta >= threshold (Montgomery & Dietrich).
//                          Better: channel heads are steeper at smaller area.
//   extraction             a channel raster from lsdtt-channel-extraction.
//   vector                 a mapped network -- OS Open Rivers, MasterMap Water
//                          Network -- burned in and snapped to the flow field.
//
// ON SNAPPING
// -----------
// A burned vector line rarely lands exactly on the modelled flow path: the
// survey and the DEM disagree by a pixel or two. Snapping moves each burned
// pixel to the largest-drainage-area cell within a small radius.
//
// This does NOT use lsdtt_flow_snapping_cpu, deliberately. Those routines snap
// on `flow.s_vector_size(flat)`, and lsdtt_flowgraph_cpu does not build an
// s-vector at all -- with D-infinity or MFD the upslope set is not a tree, so
// there is nothing for an s-vector to index. Snapping here therefore works on
// the drainage-area raster, which every flow method produces.
//
// References
// ----------
//   Montgomery D.R. & Dietrich W.E. (1988) Nature 336:232-234
//   Montgomery D.R. & Dietrich W.E. (1992) Science 255:826-830
//   Montgomery D.R. & Foufoula-Georgiou E. (1993) WRR 29:3925-3934
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

enum class channel_mask_mode
{
    none,
    area_threshold,
    area_slope_threshold,
    extraction,
    vector
};

struct channel_mask_options
{
    channel_mask_mode mode = channel_mask_mode::none;

    // --- mode none: the SAME geometry the velocity closure uses ---
    // These must match hydro_velocity_options or the emergent mask describes a
    // different channel head from the one the model actually has.
    double hydraulic_radius_coeff      = 0.1;
    double hydraulic_radius_exponent   = 0.4;
    double reference_area_m2           = 1.0e6;
    double sheet_flow_min_depth_metres = 0.001;

    // --- mode area_threshold ---
    double channel_area_threshold_m2 = 1.0e5;

    // --- mode area_slope_threshold: channel where A * S^theta >= threshold ---
    double area_slope_theta     = 2.0;
    double area_slope_threshold = 1.0e2;
    double min_slope            = 1.0e-5;

    // --- mode vector ---
    std::string channel_vector_fname;

    // Snap each burned pixel to the largest drainage area within this radius.
    // 0 disables snapping and takes the burn as-is.
    int snap_search_radius_pixels = 2;

    // A snapped cell must reach this area, or the burned pixel is dropped and
    // counted. Prevents a mis-registered line from painting a hillslope as
    // channel. 0 disables the check.
    double snap_min_area_m2 = 0.0;
};

struct channel_mask_result
{
    lsdtt_raster_cpu mask;            // 1 = channel, 0 = not, nodata = no data

    std::int64_t n_channel = 0;
    std::int64_t n_valid   = 0;

    // vector mode only
    std::int64_t n_burned          = 0;   // pixels the lines touched
    std::int64_t n_snapped         = 0;   // burned pixels moved by snapping
    std::int64_t n_rejected_by_area = 0;  // dropped for failing snap_min_area_m2

    [[nodiscard]] double channel_fraction() const noexcept
    { return n_valid > 0 ? static_cast<double>(n_channel) /
                           static_cast<double>(n_valid) : 0.0; }
};

// ---------------------------------------------------------------------------
// The drainage area at which the geometry ITSELF puts the channel head:
//
//     a_R (A/A_ref)^c_R = h_sheet   =>   A = A_ref (h_sheet/a_R)^(1/c_R)
//
// Exposed because it is the number a user actually wants when asking "what
// threshold is `none` implying?" -- and because it makes the point that mode
// `none` has an implied threshold that is DERIVED rather than chosen.
//
// Throws std::invalid_argument for a non-positive coefficient or exponent.
// ---------------------------------------------------------------------------
[[nodiscard]] double emergent_channel_area_m2(const channel_mask_options& opts);

// ---------------------------------------------------------------------------
// Build the mask.
//
// drainage_area is always required. slope is required by area_slope_threshold.
// supplied_mask is required by `extraction` (any non-zero value counts).
//
// Throws std::invalid_argument for a missing required input, a grid mismatch,
// or an out-of-range option; std::runtime_error if a vector file cannot be
// read.
// ---------------------------------------------------------------------------
[[nodiscard]]
channel_mask_result compute_channel_mask(
    const lsdtt_raster_cpu&      drainage_area,
    const lsdtt_raster_cpu*      slope,
    const lsdtt_raster_cpu*      supplied_mask,
    const lsdtt_raster_meta&     meta,
    const channel_mask_options&  opts
);

} // namespace lsdtt
