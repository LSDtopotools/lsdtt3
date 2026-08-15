// lsdtt_paired_flowpaths_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_paired_flowpaths_cpu
//
// Lengthwise asymmetry of paired flow paths -- the planform-complexity measure
// of Goren, L. and Shelef, E. (2024), "Channel concavity controls planform
// complexity of branching drainage networks", Earth Surf. Dynam. 12, 1347-1369,
// doi:10.5194/esurf-12-1347-2024, building on Shelef, E. and Hilley, G. E.
// (2014), Geophys. Res. Lett. 41, 3485-3493, doi:10.1002/2014GL059816.
//
// For two flow paths that DIVERGE at a divide and REJOIN at a common junction:
//
//     dL   = 2 |L_ij - L_ji| / (L_ij + L_ji)      ("network complexity")
//     dchi = 2 |chi_ij - chi_ji| / (chi_ij + chi_ji)
//
// where L_ij is the along-channel distance from where path i enters the channel
// network down to the junction the two paths share. Both are normalised to
// [0, 2], so they are dimensionless and scale-free.
//
// Published end-members, and they are the contrast this module exists to
// measure: a simple network of sub-parallel trunk streams running down the
// regional gradient has median dL = 0.19; a complex, tortuous network has
// median dL = 1.00.
//
// -----------------------------------------------------------------------------
// WHICH "DIVIDE" THIS USES, AND WHY IT IS NOT THE OTHER TWO
//
// lsdtt3 has three divide-related modules answering three different questions:
//
//   lsdtt_divide_migration_cpu    is this divide MOVING, and which way?
//                                 (cross-divide gradient / chi / relief)
//   lsdtt_ridgeline_cpu           what is this ridgeline LIKE?
//                                 (the divide as a geometric object)
//   lsdtt_divide_diagnostics_cpu  is the surface SELF-CONSISTENT across the
//                                 divide? (elevation step, an equilibration
//                                 convergence diagnostic)
//
// This module belongs with the SECOND. Goren & Shelef iterate over edges of a
// divide network to characterise its geometry, and their junction-proximity
// filter is an along-divide distance -- both are ridgeline properties. It is
// not a migration metric (nothing here is about disequilibrium) and not a
// divide-step metric (nothing here is about elevation discontinuity).
//
// The dchi computed here is likewise NOT the cross-divide dchi that
// lsdtt_divide_migration_cpu uses as a migration indicator. This one is the chi
// asymmetry BETWEEN PAIRED FLOW PATHS, a planform measure. Same symbol,
// different construction, different meaning.
//
// -----------------------------------------------------------------------------
// WHERE THE DIVIDE EDGES COME FROM
//
// The pairs must flank an INTERNAL divide -- one where both sides ultimately
// rejoin -- because a pair whose paths never meet has no common junction and no
// dL. Goren & Shelef enforce this by requiring both flanking pixels to lie in
// the same basin.
//
// So the input is a set of INTERFLUVIAL divide ridgelines, from
// lsdtt_basin_outline_cpu::get_interfluvial_divide_as_ridgeline(a, b, network),
// which returns the boundary shared between two basins. Basin OUTLINE
// ridgelines are the wrong object: an outline separates a basin from its
// neighbours, so its flanking pixels are in different basins and every pair
// would be discarded.
//
// Ridgeline nodes are on the pixel-CORNER divide grid ((n_rows+1) x (n_cols+1),
// shifted half a cell; periodic axes keep n). Consecutive nodes therefore span
// one pixel edge, and the two DEM pixels flanking that edge are the pair.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_ridgeline_cpu.hpp"

namespace lsdtt
{

struct lsdtt_paired_flowpath_options
{
    // Reject divide edges within this distance, ALONG THE DIVIDE, of either end
    // of the ridgeline. Ends are where the divide meets a junction, and pairs
    // there are trivially symmetric -- they have barely separated before they
    // rejoin -- so including them drags the median down. Goren & Shelef use
    // 1000 m. Set 0 to disable.
    double threshold_divide_distance_m = 1000.0;

    // Goren & Shelef iterate every divide edge and take data from the nearest
    // downstream channel, so two edges can yield the SAME channel pair and be
    // counted twice. false reproduces them; true keeps one record per distinct
    // channel pair.
    bool unique_pairs = false;

    // Goren & Shelef do NOT restrict pairs to channel heads: a paired path may
    // start further down a channel. true keeps only pairs where both nearest
    // channels are sources.
    bool only_take_channel_heads = false;
};

struct lsdtt_paired_flowpath
{
    std::int64_t pixel_a = -1;      // flat index of the flanking DEM pixel, side A
    std::int64_t pixel_b = -1;
    std::int64_t channel_a = -1;    // nearest downstream channel of each
    std::int64_t channel_b = -1;
    std::int64_t common_junction = -1;   // flat index where the paths rejoin

    double divide_distance_m = std::numeric_limits<double>::quiet_NaN();

    double L_a = std::numeric_limits<double>::quiet_NaN();
    double L_b = std::numeric_limits<double>::quiet_NaN();
    double chi_a = std::numeric_limits<double>::quiet_NaN();
    double chi_b = std::numeric_limits<double>::quiet_NaN();

    double delta_L   = std::numeric_limits<double>::quiet_NaN();
    double delta_chi = std::numeric_limits<double>::quiet_NaN();
};

struct lsdtt_paired_flowpath_result
{
    std::vector<lsdtt_paired_flowpath> pairs;

    std::size_t n_pairs           = 0;
    std::size_t n_edges_examined  = 0;
    std::size_t n_rejected_distance = 0;   // inside threshold_divide_distance_m
    std::size_t n_rejected_no_channel = 0; // a side had no downstream channel
    std::size_t n_rejected_no_junction = 0;// the paths never rejoin
    std::size_t n_rejected_duplicate = 0;
    std::size_t n_rejected_not_head  = 0;

    double median_delta_L   = std::numeric_limits<double>::quiet_NaN();
    double mean_delta_L     = std::numeric_limits<double>::quiet_NaN();
    double median_delta_chi = std::numeric_limits<double>::quiet_NaN();
};

// -----------------------------------------------------------------------------
// find_common_downstream_node
//
// The first node reached by BOTH flow paths -- the junction where two paths
// rejoin. Returns -1 if they never meet (different basins, or one leaves the
// grid first).
//
// Walks a to its outlet marking every node, then walks b until it hits a mark.
// -----------------------------------------------------------------------------
[[nodiscard]] std::int64_t find_common_downstream_node(
    const lsdtt_flow_base& flow,
    std::int64_t           a,
    std::int64_t           b
);

// -----------------------------------------------------------------------------
// compute_paired_flowpath_statistics
//
// `fd_at_channel` and `chi_at_channel` are the flow distance and chi OF the
// nearest downstream channel, per pixel -- i.e. the output of
// sample_at_nearest_channel(), not the distance TO the channel. That
// distinction is the whole reason sample_at_nearest_channel exists.
//
// `fd_from_outlet` and `chi_from_outlet` are the plain flow-distance and chi
// rasters, sampled at the common junction so that BOTH L and chi are
// differences against a shared datum:
//
//     L_i   = fd(channel_i)  - fd(junction)
//     chi_i = chi(channel_i) - chi(junction)
//
// THE INHERITED "ELEVATION CORRECTION" GAP. The reference implementation
// (LSDTT_Development, LSDDivides.cpp:2632) carries the comment
// "DOES NOT DO THE ELEVATION CORRECTION YET!" on dchi. Reading the code around
// it, the gap is exactly the datum asymmetry described above: there L_ij is
// differenced against the common junction but chi_ij is the raw chi, taken
// straight from the nearest-channel raster with no junction term. chi is an
// elevation-like coordinate, so referring it to the junction IS the elevation
// correction, and subtracting chi at the junction here closes the gap rather
// than inheriting it. (If the original author meant something else by the
// phrase -- a correction using true elevation rather than chi -- that would be
// an additional step, not a substitute for this one: L and chi have to share a
// datum either way.)
//
// The subtraction is not cosmetic. Two paths that rejoin share their entire
// downstream length, and from the outlet datum that shared part is common to
// both terms -- it cancels in the numerator of dL but not the denominator, so
// leaving it in drives the ratio towards zero and makes every network look
// symmetric. The junction datum is what makes the measure local to the pair.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_paired_flowpath_result compute_paired_flowpath_statistics(
    const lsdtt_flow_base&                  flow,
    const std::vector<lsdtt_ridgeline_cpu>& divide_ridgelines,
    const std::vector<float>&               fd_at_channel,
    const std::vector<float>&               chi_at_channel,
    const std::vector<float>&               fd_from_outlet,
    const std::vector<float>&               chi_from_outlet,
    const std::vector<std::uint8_t>&        channel_mask,
    const lsdtt_paired_flowpath_options&    opts = {}
);

// Per-pair CSV (one row per paired flow path).
void write_paired_flowpath_data(
    const lsdtt_paired_flowpath_result& result,
    const lsdtt_raster_meta&            meta,
    const std::string&                  out_path
);

// Single-row summary CSV (median dL and dchi, plus the rejection counts, which
// are how you tell "no complexity" from "no pairs survived the filters").
void write_paired_flowpath_summary(
    const lsdtt_paired_flowpath_result& result,
    const std::string&                  out_path
);

} // namespace lsdtt
