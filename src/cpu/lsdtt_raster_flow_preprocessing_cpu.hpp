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

#include <cstdint>
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
// fill_pits (periodic-aware overload)
//
// Same filling, but honouring a NESW boundary_conditions string so that a
// PERIODIC edge is treated as what it is: not a boundary. Two things change on
// a periodic axis --
//
//   1. its edges are NOT seeded as spill points, because water does not leave
//      the grid there; and
//   2. the neighbour walk WRAPS across the seam, so a depression that spills
//      round the edge is filled to the level it actually spills at.
//
// Without this, filling with open edges and then routing with 'p' leaves cells
// that have no downhill neighbour once the router wraps, and routing fails with
// "could not assign a receiver".
//
// Periodic in BOTH directions is refused: the grid would then have no drainage
// exit at all and every cell would fill to the global maximum.
//
// Always uses the priority-flood backend, even at min_slope <= 0: the
// reconstruction fast path has its own neighbour walk and does not wrap.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    const std::string&      boundary_conditions,
    int                     n_threads = 0
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
// fill_pits_with_sinks
//
// Depression filling that RETAINS a designated set of sinks.
//
// This is fill_pits_priority_flood() -- the Barnes et al. (2014) Improved
// Priority-Flood -- with one change: each cell listed in `sink_flats` joins the
// seed set alongside the raster border and the cells that touch an interior
// nodata void. A sink is therefore treated as a drainage exit, exactly as an
// interior nodata hole already is, but WITHOUT punching a hole in the DEM.
//
// The consequence is the point of the function. Plain filling raises a closed
// depression until it spills over its lowest saddle, which destroys it: flow
// routing then carries the depression's water out over that saddle. Seeding the
// flood from the sink instead means the flood works OUTWARDS from the
// depression floor, the depression is never raised, and every cell inside it
// keeps a descending path down to the sink. This is what internal drainage --
// moulins, playas, sabkhas, sinkholes, quarries -- requires.
//
// A sink cell is seeded at its OWN original elevation and is never raised.
//
// sink_flats
//   Flat (row-major) indices, i.e. row * n_cols + col. Order is irrelevant and
//   duplicates are harmless. An index outside the raster, or one landing on a
//   nodata or masked cell, throws std::invalid_argument -- a sink that is not a
//   real cell cannot act as an exit and silently dropping it would leave the
//   depression filled.
//
// min_slope
//   Gradient in m/m (rise over run) imposed across the flats produced by
//   filling; a cell is raised by min_slope * cell_size above its neighbour, or
//   min_slope * cell_size * sqrt(2) on a diagonal. min_slope <= 0 leaves true
//   flats, which lsdtt3's D8 router will not accept. See the FLAT HANDLING note
//   at the top of this header.
//
// An empty `sink_flats` returns a raster identical to
// fill_pits_priority_flood(raster, min_slope) -- there are no extra seeds, so
// it is the same flood. Note that this is the PRIORITY-FLOOD backend, not
// fill_pits(), which dispatches to the reconstruction backend when
// min_slope <= 0.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits_with_sinks(
    const lsdtt_raster_cpu&          raster,
    const std::vector<std::int64_t>& sink_flats,
    float min_slope = 0.0001f,
    int   n_threads = 0
);

// -----------------------------------------------------------------------------
// hydraulic_surface_fill_result / fill_hydraulic_surface
//
// A priority-flood projection of a hydraulic water surface.  `node_state`
// uses 0 = nodata/outside, 1 = active, and 3 = selected outlet. Code-3 cells
// are prescribed sinks wherever they occur. With check_valid_outlet=true,
// every connected active component needs one. `n` edges are closed and `p`
// edges wrap periodically.
//
// Unlike the historic terrain-fill helpers, nodata does NOT create drainage
// exits here: it is an impermeable excluded cell.  This distinction is needed
// by hydraulic models with holes, walls, or deliberately selected outlets.
//
// The returned added-depth field is `filled_surface - surface` on active
// cells.  Calling code that writes this into water depth must account for the
// resulting projection volume explicitly; the operation is not conservative.
// -----------------------------------------------------------------------------
struct hydraulic_surface_fill_result
{
    lsdtt_raster_cpu filled_surface;
    lsdtt_raster_cpu added_depth;
    double           added_volume_m3 = 0.0;
};

[[nodiscard]] hydraulic_surface_fill_result fill_hydraulic_surface(
    const lsdtt_raster_cpu& surface,
    const lsdtt_array2d_u8& node_state,
    const std::string&      boundary_conditions,
    float                   fill_slope_epsilon = 1.0e-8f,
    bool                    check_valid_outlet = true
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

// =============================================================================
// PERONA-MALIK ANISOTROPIC DIFFUSION
//
// Perona, P. & Malik, J. (1990) Scale-space and edge detection using anisotropic
// diffusion. IEEE PAMI 12(7), 629-639. Applied to topography by Passalacqua et
// al. (2010) WRR 46:W11535 (GeoNet) and Clubb et al. (2017) ESurf 5:369-385.
//
//     dz/dt = div( c(|grad z|) grad z )
//
//     c(s) = exp(-(s/K)^2)      [P&M eq. 3, `exponential`]
//     c(s) = 1 / (1 + (s/K)^2)  [P&M eq. 4, `quadratic`]
//
// Discretised on the 4-neighbour stencil (P&M eq. 8):
//
//     z(i,j)^{t+1} = z(i,j)^t + lambda * SUM_d c(|nabla_d z|) * nabla_d z
//
// for d in {N, S, E, W}, where nabla_d z is the neighbour difference.
//
// WHY THIS AND NOT A LOW-PASS FILTER. A mean or Gaussian filter removes DEM
// noise and the tread/riser break in equal measure. For terrace and valley-floor
// work the riser IS the signal: it is what a HAND break sits on, and blurring it
// moves every break the ladder finds. Anisotropic diffusion suppresses
// within-facet noise while leaving a high-gradient riser essentially untouched,
// because c(s) collapses toward zero wherever the local difference exceeds K.
//
// WHAT K MEANS AND WHY IT IS NOT GUESSED. K is the edge threshold in the SAME
// units as the stencil differences - metres of rise between adjacent cells. A
// difference above K is preserved, one below it is smoothed. A number carried
// over from another landscape is meaningless, so following Black et al. (1998)
// it is set to a high quantile of the DEM's own distribution of those
// differences - see stencil_difference_magnitude, and note that a
// central-difference gradient is NOT the right quantity to read it from. The
// value actually used comes back in the result, because it is the one number a
// reader needs to reproduce the conditioning.
//
// STABILITY. lambda <= 0.25 on this stencil. Above that the explicit update
// diverges, so it is REJECTED rather than warned about: an unstable run does not
// produce a worse DEM, it produces a meaningless one.
//
// This is a CONDITIONING step. Run it before filling or breaching, keep the
// smoothed raster as its own file, and let everything downstream read that -
// so which DEM an analysis was run on is visible rather than inferred.
// =============================================================================

// Which conductance function c(s) is used. See the block comment above.
enum class pm_conductance_choice
{
    exponential,   // exp(-(s/K)^2)      - P&M eq. 3
    quadratic      // 1 / (1 + (s/K)^2)  - P&M eq. 4
};

// -----------------------------------------------------------------------------
// perona_malik_options
//
// n_iterations        - explicit timesteps. More smooths further; the riser is
//                       preserved at every count, so this trades runtime for
//                       how much within-facet noise goes.
// lambda              - timestep. Must be in (0, 0.25]; see STABILITY above.
// conductance         - which c(s).
// k_gradient_quantile - quantile of the stencil-difference distribution (metres,
//                       see stencil_difference_magnitude) used to set K. 0 or
//                       negative means "use k_metres_per_cell instead".
// k_metres_per_cell   - K stated directly. Only read when k_gradient_quantile
//                       is not positive. Must then be > 0.
// n_threads           - OpenMP threads; 0 = runtime default.
// -----------------------------------------------------------------------------
struct perona_malik_options
{
    int    n_iterations        = 10;
    double lambda              = 0.20;
    pm_conductance_choice conductance = pm_conductance_choice::exponential;
    double k_gradient_quantile = 0.90;
    double k_metres_per_cell   = 0.0;
    int    n_threads           = 0;
};

// -----------------------------------------------------------------------------
// perona_malik_result
//
// smoothed          - the conditioned raster. Nodata cells are untouched and
//                     never diffused into.
// k_metres_per_cell - the K actually used, whether stated or estimated.
// n_iterations      - iterations actually run.
// -----------------------------------------------------------------------------
struct perona_malik_result
{
    lsdtt_raster_cpu smoothed;
    double k_metres_per_cell = 0.0;
    int    n_iterations      = 0;
};

// -----------------------------------------------------------------------------
// stencil_difference_magnitude
//
// Per cell, the LARGEST absolute elevation difference to a valid 4-neighbour,
// in metres. Nodata where the cell itself is nodata or has no valid neighbour.
//
// This is the exact quantity the conductance function is evaluated on, which is
// the point: K has to be a threshold on what the stencil actually sees. A
// central-difference gradient is the intuitive choice and is WRONG here,
// because it averages across two cells and so cancels one-cell noise the
// diffusion does see - on a checkerboard DEM it reports zero gradient
// everywhere and the estimated K collapses to nothing.
//
// Exposed because a K quantile is only interpretable against the distribution
// it was taken from.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu stencil_difference_magnitude(
    const lsdtt_raster_cpu& raster
);

// -----------------------------------------------------------------------------
// perona_malik_diffusion
//
// Throws:
//   std::runtime_error if n_iterations < 1, if lambda is outside (0, 0.25], if
//   neither K source is usable, or if the raster metadata are invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] perona_malik_result perona_malik_diffusion(
    const lsdtt_raster_cpu&      raster,
    const perona_malik_options&  options
);

} // namespace lsdtt
