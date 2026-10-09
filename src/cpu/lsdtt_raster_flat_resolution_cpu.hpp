// lsdtt_raster_flat_resolution_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Flat resolution on a depression-filled DEM.
//
// Algorithm:
//   Barnes, R., Lehman, C., Mulla, D. (2014) An efficient assignment of
//   drainage direction over flat surfaces in raster digital elevation models.
//   Computers & Geosciences, 62, 128-135.
//   DOI: 10.1016/j.cageo.2013.01.009
//
// The float-step variant raises each flat cell by a whole number of float32
// ULPs (the Barnes mask value), which is what RichDEM's ResolveFlats does on a
// float32 raster. The mask reaches thousands of steps on large flats, so the
// step is done in O(1) by step_up_float32 rather than by a nextafter loop.
//
// Naming: in lsdtt3 "flat" on its own means a row-major cell index
// (flat_index, sink_flats). Equal-elevation regions are therefore called
// "flat regions" here (flat_region_label, n_flat_regions).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// step_up_float32
//
// Returns the float32 value reached by k successive std::nextafter(z, +inf)
// calls, bitwise identical to that loop (including the sign of zero), in O(1).
//
// Method: the float32 bit pattern is mapped to a monotone signed integer key
// (positive floats -> their bits, negative floats -> minus their magnitude
// bits, so +0.0 and -0.0 both map to 0), k is added, and the key is mapped back.
// A result key of 0 reached with k > 0 can only come from the negative side,
// and nextafter(-denorm_min, +inf) is -0.0, so that case returns -0.0
// (bit pattern 0x80000000), not +0.0.
//
// Precondition: k >= 0. (A negative k is not supported.)
//
// Behaviour:
//   * k == 0                    returns z unchanged, bitwise (-0.0 stays -0.0)
//   * stepping past FLT_MAX     returns +inf, as nextafter does
//   * z == +inf                 returns +inf
//   * z == -inf                 steps up from -FLT_MAX, as nextafter does
//   * z is NaN                  returns z unchanged (same bit pattern)
// The flat-resolution caller never reaches the infinite or NaN cases: it only
// steps finite, valid elevations by a bounded mask.
float step_up_float32(float z, std::int32_t k) noexcept;


// barnes2014b_flat_mask / compute_barnes2014b_flat_mask
//
// The Barnes et al. (2014b) flat mask: for every cell of every drainable flat
// region, the number of increments that makes the region drain towards its
// outlets and away from the higher terrain around it. This is the shared core
// of both Barnes flat-resolution options; it does not modify any elevation.
//
// Precondition: `filled` is depression-filled with true flats, i.e.
// fill_pits_reconstruction, or fill_pits with min_slope 0. Every flat must have
// an outlet. If any does not, std::invalid_argument is thrown with the message
//   lsdtt::compute_barnes2014b_flat_mask: <n> flat cell(s) in <m> flat
//   region(s) have no outlet; the input must be depression-filled first
//   (fill_pits_reconstruction or fill_pits with min_slope 0)
// (one line), where n counts the NO_FLOW cells left unlabelled and m their
// 8-connected equal-elevation components.
//
// Method (all comparisons are exact float comparisons, as in Barnes and
// RichDEM; neighbours are 8-connected, in grid::d8_clockwise order):
//   0. Validity comes from the raster's mask array, never from comparing to
//      no_data_value, so NaN nodata works.
//   1. A valid cell DRAINS if it has a strictly lower valid neighbour, or if it
//      is an exit: on the raster border, or touching nodata. This is the exit
//      rule of the fills. Every other valid cell is NO_FLOW. (OpenMP-parallel.)
//   2. Low edge: a draining cell with a NO_FLOW neighbour of equal elevation.
//      High edge: a NO_FLOW cell with a strictly higher neighbour.
//   3. Labels: an equal-elevation flood fill from each unlabelled low edge, in
//      raster scan order. A label covers every connected equal-elevation
//      valid cell, draining or not.
//   4. Away distance a: BFS from the high edges (a = 1) through same-label
//      NO_FLOW cells. H[label] = max a over the label (0 if none).
//   5. Towards distance t: BFS from the low edges (t = 1) through same-label
//      NO_FLOW cells. mask = 2t + H[label] - a where a > 0, else 2t.
// t, a and H are BFS distances or maxima over a fixed graph, so the result is
// independent of scan, queue and thread order.
//
// What mask means: the number of increments to add to the cell's elevation.
// mask is 0 for every cell outside a flat region and for labelled draining
// cells that are not low edges. Low edges get 2, which is what RichDEM does.
// Border rule: the mask here is the pure value above for every labelled
// cell, INCLUDING low edges on the raster border (they get 2). RichDEM never
// modifies raster-border cells; that rule is applied by the resolver, not here.
//
// Descent guarantee: every NO_FLOW cell c has a same-label neighbour n with
// t(n) = t(c) - 1. If n is NO_FLOW, mask(c) - mask(n) = 2 + a(n) - a(c), which
// is in {1, 2, 3}. If n is a low edge (mask 2), mask(c) - 2 is in [2, H + 1].
// So every NO_FLOW cell has a neighbour at least one increment lower.
//
// Memory (dense): int32 label and int32 mask (both returned), a uint8 drain
// status and a transient uint8 edge flag: at most 10 B per cell. The away
// distance is held in the array that becomes the mask, and the status is
// freed before the towards pass. On top of that: the edge-cell lists and the
// BFS queues, up to 8 B per flat-region cell each.
//
// n_threads: OpenMP threads for the drain-status and edge scans (0 = all).
// The flood fill and both BFS passes are serial. The result does not depend
// on n_threads.
struct barnes2014b_flat_mask
{
    std::vector<std::int32_t> flat_region_label;  // row-major; 0 = not in a flat region
    std::vector<std::int32_t> mask;               // row-major; 0 = not in a flat region

    // Diagnostics
    std::int64_t n_flat_regions    = 0;  // number of labels
    std::int64_t n_flat_cells      = 0;  // cells with label > 0
    std::int64_t n_no_flow         = 0;  // valid cells with no lower neighbour, not exits
    std::int64_t n_low_edge        = 0;
    std::int64_t n_high_edge       = 0;
    std::int64_t n_no_outlet_cells = 0;  // always 0 on return (non-zero throws)
    std::int32_t max_mask          = 0;
    std::int32_t max_H             = 0;  // largest away distance over all regions
    std::int32_t max_T             = 0;  // largest towards distance over all regions
};

barnes2014b_flat_mask compute_barnes2014b_flat_mask(const lsdtt_raster_cpu& filled,
                                                    int n_threads = 0);

}  // namespace lsdtt
