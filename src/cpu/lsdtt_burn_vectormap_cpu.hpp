// lsdtt_burn_vectormap_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_burn_vectormap
//
// Stream-burning strategy tailored to an ACCURATE vector river layer (OS
// VectorMap Local): follow the LiDAR flow paths as closely as possible and burn
// only where an obstruction (mostly a road) breaks the channel.
//
// Algorithm
// ---------
//   1. Compute a flow-accumulation raster on the filled DEM (contributing
//      pixels), and define a D8 channel mask: cells with accumulation above a
//      threshold.
//   2. Sample each vector line. At each sample point, look for an above-
//      threshold channel cell within a snap radius:
//        - found  -> "on-channel": the LiDAR already carries this reach; leave
//          the DEM untouched (the line is effectively snapped to the channel).
//        - none   -> "gap": no channel here. Just downstream of a road the
//          accumulation falls to ~0, so road crossings appear as gaps.
//   3. For each gap bracketed by on-channel reaches, burn a monotone-descending
//      trench directly under the vector line, anchored to the channel cells on
//      each side, so D8 routes across the obstruction.
//
// The result preserves the LiDAR thalweg everywhere a real channel exists and
// carves only the obstructed spans. There is no source seeding, so a vector
// source that happens to sit on a rough hillslope cannot divert the result.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct vectormap_burn_options
{
    // A DEM cell belongs to the D8 channel mask when its number of contributing
    // pixels (include-self) is at least this value.
    std::int64_t accumulation_threshold_pixels = 1000;

    // Lateral search radius (metres) for snapping a vector sample point to an
    // above-threshold channel cell. Points with no channel within this radius
    // are flagged as gap (obstruction) points.
    double snap_radius_m = 30.0;

    // Gaps shorter than this (metres) are ignored (treated as on-channel noise).
    double min_gap_length_m = 0.0;

    // Gaps longer than this (metres) are reported but NOT burned: a span this
    // long with no nearby channel is unlikely to be a single obstruction.
    // Set <= 0 to disable the cap.
    double max_gap_length_m = 250.0;

    // Enforced descending gradient (m/m) along a burned gap trench.
    double burn_min_slope = 0.0001;

    // Sampling step (metres) along vector lines. If <= 0, defaults to half the
    // DEM cell size.
    double sample_step_m = -1.0;

    int n_threads = 1;
};

// One burned obstruction gap, retained for diagnostic output.
struct vectormap_burn_gap
{
    std::vector<lsdtt_xy> path_xy;   // anchor -> line-across-gap -> anchor
    double length_m = 0.0;
    bool   burned   = false;         // false if skipped (too long / unanchored)
    std::string skip_reason;
};

struct vectormap_burn_result
{
    int          n_lines              = 0;
    int          n_gaps_detected      = 0;
    int          n_gaps_burned        = 0;
    int          n_gaps_skipped       = 0;
    double       total_burn_length_m  = 0.0;
    std::int64_t channel_mask_cells   = 0;

    std::vector<vectormap_burn_gap> gaps;
};

// Burn obstruction gaps under the vector lines into `dem` (modified in place).
// `flow` must have been built on the current state of `dem` (the filled DEM);
// its accumulation raster defines the channel mask. The vector file is assumed
// to be in the DEM CRS.
[[nodiscard]] vectormap_burn_result burn_vectormap_road_gaps(
    lsdtt_raster_cpu&              dem,
    const lsdtt_flow_base&        flow,
    const std::string&            vector_path,
    const vectormap_burn_options& options = vectormap_burn_options{}
);

// Channel-head seeds for network extraction by accumulation threshold: cells
// whose contributing-pixel count is >= threshold and which have no upstream
// donor cell that also meets the threshold.
[[nodiscard]] std::vector<std::int64_t> channel_heads_from_accumulation(
    const lsdtt_flow_base& flow,
    std::int64_t           accumulation_threshold_pixels
);

} // namespace lsdtt
