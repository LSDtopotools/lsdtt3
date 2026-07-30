// lsdtt_raster_flow_preprocessing_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Flow preprocessing - CPU free functions that condition a DEM prior to flow
// routing. All functions are stateless: they take a raster and return a new
// raster. The input is never modified.
//
// Fill algorithm:
//   Barnes, R., Lehman, C., Mulla, D. (2014) Priority-flood: An optimal
//   depression-filling and watershed-labeling algorithm for digital elevation
//   models. Computers & Geosciences, 62, 117-127.
//   DOI: 10.1016/j.cageo.2013.04.024
//
// -----------------------------------------------------------------------------
// FLAT HANDLING: lsdtt3 and RichDEM do DIFFERENT things
//
// RichDEM is used as the external reference in the flow-routing benchmark
// (docs/testing_and_benchmarking.md), so it is worth being explicit about
// where the two agree and where they deliberately do not.
//
//   min_slope == 0  -- equivalent to RichDEM's FillDepressions(epsilon=false).
//                      Verified bit-exact against it on both a high-relief and
//                      a low-relief real DEM.
//
//   min_slope > 0   -- NOT equivalent to RichDEM's FillDepressions(epsilon=true),
//                      despite both being described as "epsilon fill".
//
//                      lsdtt3 imposes a true SLOPE: each cell in a flat is
//                      raised by min_slope * cell_size (times sqrt(2) on the
//                      diagonal). The imposed gradient is uniform across the
//                      raster, scales correctly with resolution, and is
//                      independent of elevation and of the storage type.
//
//                      RichDEM instead raises each cell to the next
//                      representable float above its neighbour
//                      (std::nextafter, i.e. one ULP). That step depends on
//                      elevation magnitude and data type: ~1.2e-4 m at 2000 m
//                      in float32, ~7.6e-6 m at 100 m, and a denormal 1.4e-45 m
//                      near zero elevation. RichDEM therefore refuses integer
//                      rasters outright.
//
// Consequence for anyone comparing outputs: the two agree closely in high
// relief, where the flat-gradient term is a negligible fraction of the total
// fill volume (measured 0.01 % on an Alpine DEM), and diverge substantially in
// low relief, where it is not (29.5 % on a Finnish Lakeland DEM, giving a 40x
// difference in material moved). Neither is wrong; they are answering the
// question differently. The slope-based formulation here is deliberate -- 
// the ULP approach was evaluated and rejected for lsdtt3 because of potential 
// errors in data type and rasters that cross 0 elevation.
// -----------------------------------------------------------------------------
//
// Breach algorithm:
//   Lindsay, J.B. (2016) Efficient hybrid breaching-filling sink removal
//   methods for flow path enforcement in digital elevation models.
//   Hydrological Processes, 30(6), 846-857.
//   DOI: https://doi.org/10.1002/hyp.10648
//   Adapted from RichDEM (r-barnes/richdem, MIT licence).
//
// Least-cost breach algorithm:
//   Lindsay, J.B. (2020) Pit-centric depression removal methods.
//   Geomorphometry 2020 Conference Proceedings.
//   Implementation reference: WhiteboxTools BreachDepressionsLeastCost.
// -----------------------------------------------------------------------------

#pragma once

#include <limits>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// mask_to_nodata_using_threshold
// -----------------------------------------------------------------------------
lsdtt_raster_cpu mask_to_nodata_using_threshold(
    const lsdtt_raster_cpu& raster,
    float threshold,
    bool  is_below_threshold_nodata
);

// -----------------------------------------------------------------------------
// fill_pits
//
// Depression filling. Dispatches on min_slope:
//
//   min_slope <= 0  -> fill_pits_reconstruction() (Vincent 1993), which is
//                      2-3x faster and produces a BIT-IDENTICAL result. Verified
//                      exactly equal to the priority-flood on every fixture and
//                      edge case in tests/cpp/test_flow_properties.cpp (P9),
//                      including perfectly flat surfaces and interior nodata.
//
//   min_slope >  0  -> fill_pits_priority_flood() (Barnes et al. 2014), which
//                      imposes the gradient across flats as it floods.
//                      Reconstruction cannot do this, so there is no faster
//                      option for this case.
//
// The dispatch is an optimisation only: for a given min_slope the output is
// unchanged. Call the two backends directly if you need to pin one.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits(
    const lsdtt_raster_cpu& raster,
    float min_slope = 0.0001f,
    int   n_threads = 0
);

// -----------------------------------------------------------------------------
// fill_pits_priority_flood
//
// The Barnes et al. (2014) Improved Priority-Flood backend, always, regardless
// of min_slope. This is what fill_pits() used to be unconditionally; it is
// exposed so that the reconstruction backend can be compared against it without
// the comparison becoming circular once fill_pits() started dispatching.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits_priority_flood(
    const lsdtt_raster_cpu& raster,
    float min_slope = 0.0001f,
    int   n_threads = 0
);

// -----------------------------------------------------------------------------
// fill_pits_reconstruction
//
// Alternative depression-filling backend: grayscale morphological
// reconstruction by erosion (Vincent 1993), using the hybrid algorithm --
// a forward raster scan, a backward raster scan, then a FIFO queue that
// propagates the few remaining changes.
//
// Same result, different cost class. fill_pits() uses the Barnes et al. (2014)
// Improved Priority-Flood, which is O(N log N) because nearly every cell passes
// through a binary heap. This routine is O(N) with a small constant and almost
// entirely sequential memory access, which is why TopoToolbox's equivalent is
// roughly twice as fast as ours (see docs/testing_and_benchmarking.md).
//
// *** IT DOES NOT SUPPORT min_slope > 0 ***
//
// Reconstruction fills depressions to exactly their spill level, leaving true
// flats. It has no mechanism for imposing a gradient across those flats, which
// is what fill_pits(min_slope > 0) does inline during the flood. The two agree
// bit-for-bit at min_slope == 0 and cannot agree otherwise -- so this is a fast
// path for the no-gradient case, not a drop-in replacement.
//
// Callers that need drainable flats (lsdtt3's own D8 router throws on a truly
// flat surface) must use fill_pits() with a positive min_slope.
//
// Throws std::invalid_argument rather than silently producing flats that the
// caller did not ask for.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits_reconstruction(
    const lsdtt_raster_cpu& raster,
    int n_threads = 0
);

// -----------------------------------------------------------------------------
// resolve_flats_gradient
//
// Imposes a drainage gradient across the flat regions of an already-filled DEM,
// so the surface can be routed.
//
// *** THIS IS NOT A PERFORMANCE OPTIMISATION -- IT IS SLOWER ***
//
// It was written to let the O(N) reconstruction backend cover the min_slope > 0
// case, and MEASURABLY FAILS at that. Reconstruction + this step runs 8-88 %
// slower than fill_pits_priority_flood() on the same input, across every size
// tested. Two independent implementations (fast-sweeping, and the work queue
// used here) both came out behind.
//
// The reason is structural, not a tuning problem. Imposing the gradient is a
// shortest-path problem over the surface: each cell must end one step above
// whichever neighbour it drains to, which is precisely the ordered propagation
// priority-flood already performs inline as it floods. Splitting it into a
// second pass means paying for that propagation separately, on top of the fill,
// and the heap it removes is cheaper than the second pass costs. Solving the
// second pass "properly" with a priority queue would simply be priority-flood
// again.
//
// It is retained because it is correct, tested, and defines the gradient
// differently from priority-flood -- by true propagated distance rather than
// FIFO first arrival -- which may be independently useful. It is not wired into
// any CLI and is not the default anywhere.
//
// Method: every cell that has a strictly lower valid neighbour can already
// drain, and is a seed at distance zero. Every other valid cell sits in a flat.
// A geodesic chamfer distance transform (1 per cardinal step, sqrt(2) per
// diagonal, propagated only through flat cells) gives each flat cell its
// distance to the nearest draining cell, and the cell is raised by
//
//     min_slope * cell_size * distance
//
// so elevation increases away from the outlet and flow runs back towards it.
//
// This reproduces what priority-flood's inline increment achieves, but not
// always the identical surface: priority-flood assigns each flat cell on FIRST
// ARRIVAL through a FIFO, which is not necessarily the shortest chamfer path,
// whereas this uses the true geodesic distance. Where they differ, this one is
// the better-defined answer -- but it is a DIFFERENT answer, so it is offered
// as a choice rather than substituted silently.
//
// Input must already be depression-free. Returns a new raster.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu resolve_flats_gradient(
    const lsdtt_raster_cpu& filled_raster,
    float min_slope,
    int   n_threads = 0
);

// -----------------------------------------------------------------------------
// fill_pits_reconstruction_with_gradient
//
// Reconstruction fill followed by resolve_flats_gradient(): the fast O(N)
// equivalent of fill_pits() for min_slope > 0.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits_reconstruction_with_gradient(
    const lsdtt_raster_cpu& raster,
    float min_slope,
    int   n_threads = 0
);

// -----------------------------------------------------------------------------
// fill_pits_with_boundary
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits_with_boundary(
    const lsdtt_raster_cpu&              raster,
    float                                min_slope,
    const std::vector<std::vector<int>>& boundary_rows,
    const std::vector<std::vector<int>>& boundary_cols,
    const std::vector<std::string>&      boundary_conditions,
    const std::string&                   outlet_label = "outlet",
    int                                  n_threads    = 0
);

// -----------------------------------------------------------------------------
// breach_depressions
// -----------------------------------------------------------------------------
lsdtt_raster_cpu breach_depressions(
    const lsdtt_raster_cpu& raster,
    int n_threads = 0
);

// -----------------------------------------------------------------------------
// breach_depressions_least_cost
// -----------------------------------------------------------------------------
lsdtt_raster_cpu breach_depressions_least_cost(
    const lsdtt_raster_cpu& raster,
    int   max_dist       = 1000,
    float max_cost       = (std::numeric_limits<float>::max)(),
    bool  min_dist       = true,
    float flat_increment = 0.0f,
    bool  fill_remainder = true,
    int   n_threads      = 0
);

// -----------------------------------------------------------------------------
// depression_polygons_result
// -----------------------------------------------------------------------------
struct depression_polygons_result
{
    lsdtt_vector_schema              schema;
    std::vector<lsdtt_vector_record> records;
};

// -----------------------------------------------------------------------------
// depressions_to_polygons
//
// Derives a vector polygon layer from a fill-difference raster.
// Each polygon represents one topographic depression with attributes:
//   depression_id  - sequential integer label (int64)
//   area_m2        - planimetric area of the depression in m^2  (double)
//   volume_m3      - total fill volume in m^3, computed as
//                    sum(fill_diff[i]) * cell_size^2             (double)
//
// Parameters
//   fill_diff      - filled minus original DEM; positive values are where
//                    filling raised the surface. Typically the output of
//                    raster_subtract(filled_dem, original_dem).
//   min_pixels     - components with fewer pixels than this are discarded.
//                    Must be >= 1.  Default 1 keeps everything.
//   min_depth      - cells where fill_diff <= min_depth are treated as
//                    background before component labelling.  Default 0.0f
//                    keeps any cell with a positive fill value.
//
// Returns
//   depression_polygons_result with schema and one record per depression.
//   Pass directly to lsdtt_data_loader::begin_vector_write /
//   write_polygon / finish_vector_write.
// -----------------------------------------------------------------------------
[[nodiscard]] depression_polygons_result
depressions_to_polygons(
    const lsdtt_raster_cpu& fill_diff,
    int   min_pixels = 1,
    float min_depth  = 0.0f
);

} // namespace lsdtt
