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


// flat_resolution_choice
//
// How flats left by a min_slope 0 fill are made drainable.
//   inline_gradient         the default: the Priority-Flood fill imposes
//                           min_slope during the flood (fill_pits). It is not
//                           a resolver and resolve_flats_barnes2014b refuses it.
//   barnes2014b_gradient    Barnes mask scaled to a physical gradient
//                           (not implemented yet).
//   barnes2014b_float_step  Barnes mask applied as float32 steps (RichDEM
//                           ResolveFlats on float32).
enum class flat_resolution_choice : std::uint8_t
{
    inline_gradient,
    barnes2014b_gradient,
    barnes2014b_float_step
};

struct flat_resolution_options
{
    flat_resolution_choice method = flat_resolution_choice::barnes2014b_float_step;
    float min_slope = 1.0e-4f;                // m/m; barnes2014b_gradient only
    // Near-zero guard (see near_zero_guard_step). Always on in every caller;
    // it is an option only so that the tests can show what it prevents.
    bool  float_step_near_zero_guard = true;
    int   n_threads = 0;                      // 0 = all
};

// Result of resolve_flats_barnes2014b. mask_diagnostics carries the counts of
// compute_barnes2014b_flat_mask; its flat_region_label and mask vectors are
// released (empty) to save 8 B per cell. Call compute_barnes2014b_flat_mask
// directly if the per-cell arrays are needed.
struct flat_resolution_result
{
    lsdtt_raster_cpu      resolved;
    barnes2014b_flat_mask mask_diagnostics;

    std::int64_t n_flats_near_zero_guarded = 0;  // flat REGIONS raised by the guard
    double       max_raise_m               = 0.0;  // max of double(z') - double(z)

    // Ordered 8-neighbour pairs (c, n), both valid, with z(n) > z(c) before
    // and z'(c) >= z'(n) after. RichDEM prints this count as "Cells
    // inappropriately raised above surrounding terrain" (it counts pairs).
    std::int64_t n_pairs_raised_above_surroundings = 0;
    // Cells c that are the first member of at least one such pair.
    std::int64_t n_cells_raised_above_surroundings = 0;
    // Valid cells that are not exits (not on the raster border, no nodata
    // neighbour) and have no strictly lower valid neighbour in z'. These are
    // the cells the D8 router will refuse ("could not assign a receiver").
    std::int64_t n_rim_cells_without_descent = 0;
};


// near_zero_guard_step
//
// The near-zero guard of the float-step option: for a flat region at
// elevation |z| < 1 m, each mask increment is a fixed 2^-23 m (the ULP of 1.0f)
// instead of one float32 ULP:
//   z' = float(double(z) + mask * 2^-23)
//
// Why: one float32 ULP at 0 m is 1.4e-45 m. lsdtt3's D8 router compares in
// double and still routes such a flat, but any consumer that forms the slope
// in float32, (z - z_r) / dist, gets exactly 0 (and TWI gets +inf). With the
// guard every D8 step on the flat is at least 2^-23 m, i.e. a float32 slope of
// at least 2^-23 / (sqrt(2) L). At |z| = 1 the native ULP is 2^-23, so the
// step size is continuous across the switch.
//
// Bound: consecutive values stay distinct only while z' < 2 m, where the ULP
// is <= 2^-23. Beyond that, two masks can round to the same float. So
// std::logic_error is thrown if the result is >= 2.0f, i.e. for mask >=
// (2 - z) * 2^23 (>= 8.4e6 for any z < 1). Real flats do not come near it.
//
// Preconditions (std::logic_error if violated): |z| < 1, mask >= 0.
float near_zero_guard_step(float z, std::int32_t mask);


// resolve_flats_barnes2014b
//
// Makes every flat of a depression-filled DEM drain, using the Barnes et al.
// (2014b) mask (compute_barnes2014b_flat_mask): flow converges across flats
// towards their outlets and away from the higher ground around them.
//
// method:
//   barnes2014b_float_step  implemented; see below.
//   barnes2014b_gradient    throws std::logic_error (not implemented yet).
//   inline_gradient         throws std::invalid_argument with the message
//     lsdtt::resolve_flats_barnes2014b: inline_gradient is the Priority-Flood
//     fill (fill_pits with min_slope > 0), not a flat resolver; choose
//     barnes2014b_float_step or barnes2014b_gradient
//
// Precondition: as compute_barnes2014b_flat_mask. Its no-outlet
// std::invalid_argument (prefix "lsdtt::compute_barnes2014b_flat_mask: ")
// propagates unchanged.
//
// Float step. Every cell with a flat-region label and mask > 0 that is NOT on
// the raster border is raised:
//   z' = step_up_float32(z, mask)            (mask float32 ULPs)
//   z' = near_zero_guard_step(z, mask)       if |z| < 1 and the guard is on
// Every cell of a flat region has the same z, so the choice is per region.
// Raster-border cells are never raised, whatever their mask (RichDEM never
// modifies them; they are exits anyway). Cells next to interior nodata are
// raised. Other cells are copied bitwise.
//
// RichDEM equivalence: on flats with |z| >= 1 m the output is bitwise equal to
// RichDEM ResolveFlats on a float32 raster, border included (checked on hand
// fixtures cell by cell). On guarded flats (|z| < 1 m) it differs by design.
//
// Rim violations (decision D6): the raise is not capped. Where a flat sits
// within mask-count ULPs of a higher neighbour, the flat can reach or pass it.
// This is reported in n_pairs/n_cells_raised_above_surroundings, not
// prevented. If it leaves a non-exit cell with no lower neighbour
// (n_rim_cells_without_descent > 0), the D8 router's existing "could not
// assign a receiver" error fires on that cell.
//
// What the slopes mean: a float-step drop is one or a few ULP(z), so the
// gradient on a resolved flat is about ULP(z) / L. It depends on elevation
// (it doubles at every power of two: 512, 1024, 2048, 4096 m) and has no
// physical meaning; it exists only to give the router a direction. Most
// physics consumers floor slope at 1e-5 to 1e-3, so to them a float-step flat
// is still flat. On large high flats the total raise (mask x ULP) can reach
// decimetres to metres; see max_raise_m.
//
// n_threads: OpenMP threads for the mask scans, the raise and the diagnostics
// (0 = all). The output and every diagnostic are bitwise independent of it.
flat_resolution_result resolve_flats_barnes2014b(const lsdtt_raster_cpu& filled,
                                                 const flat_resolution_options& options);

}  // namespace lsdtt
