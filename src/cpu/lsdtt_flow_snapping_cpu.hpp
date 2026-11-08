// lsdtt_flow_snapping_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_flow_snapping_cpu
//
// Flow-based pixel snapping helpers.
//
// These functions snap user-supplied raster pixels to cells that satisfy a
// contributing-pixel threshold derived from the flow-routing object:
//
//   contributing_pixels = flow.s_vector_size(flat)
//
// This is an upslope contributing-pixel count, not a direct-donor count.
//
// Four snapping behaviours are currently supported:
//
//   1. expanding_square
//      Search outward ring-by-ring in pixel space. At the first radius that
//      contains one or more qualifying cells, choose the candidate with:
//        a. smallest Euclidean distance in pixel units
//        b. if tied, largest contributing-pixel count
//        c. if still tied, smallest flat index
//
//   2. downstream
//      Follow the D8 receiver chain until the first cell is reached whose
//      contributing-pixel count meets or exceeds the threshold.
//
//   3. largest_area (snap_pixel_to_contributing_pixels_threshold_largest_area)
//      Take the valid cell with the most contributing pixels in the square of
//      half-width R around the input pixel (ties: smaller pixel distance, then
//      smaller flat index), provided its count meets or exceeds the threshold.
//      This is the GIS "snap pour point" rule: a point slightly off a river is
//      moved onto it, but a large R can jump onto a bigger neighbouring river.
//
//   4. reported_area (snap_pixel_to_reported_drainage_area_expanding_square)
//      For a point that comes with a reported drainage area (e.g. a sampled
//      basin). Among every channel cell (contributing pixels >= threshold) in
//      the square of radius R, take the NEAREST whose area is within a ratio
//      tolerance of the reported one; if none is, the one closest in area.
//      This stops a sample on a small tributary snapping to the trunk next to
//      it, which nearest-cell snapping cannot detect.
//   5. nearest (snap_point_to_contributing_pixels_threshold_nearest)
//      For a POINT (x, y in the raster CRS), not a pixel: the qualifying cell
//      whose centre is nearest the point itself, by straight-line distance,
//      whatever the D8 path from the point does. Unlike expanding_square,
//      distances are measured from the point, not from the centre of the
//      pixel it falls in, so a point near the edge of its pixel takes the
//      truly nearest cell. For field measurements made in a channel.
//
// Notes
// -----
// - These routines require a flow object built with s_vector enabled because
//   the threshold is based on upslope contributing pixels.
// - The expanding-square search can start from any in-bounds raster pixel.
// - The downstream search requires the starting pixel to be a valid flow cell.
// - Batch functions preserve input order.
// - Function names intentionally encode both the threshold kind and the search
//   method so that future snapping families can follow the same pattern, e.g.
//   snap_pixel_to_stream_order_threshold_expanding_square(...).
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <cstdint>
#include <limits>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"


namespace lsdtt
{



// -----------------------------------------------------------------------------
// lsdtt_snap_status
//
// Per-point snap outcome.
//
// success
//   A qualifying pixel was found and returned in snapped_* fields.
//
// input_out_of_bounds
//   The supplied input row/col or flat index does not refer to a raster cell
//   inside the flow grid.
//
// input_invalid_flow_cell
//   The input location lies inside the raster extent but is not a valid flow
//   cell. This status is currently used by downstream snapping, which requires
//   a valid starting cell in the routing graph.
//
// threshold_not_found_in_search_window
//   Expanding-square search completed through the requested maximum search
//   radius without finding any qualifying pixel, or no cell in the
//   largest-area or reported-area search square met the threshold.
//
// outlet_reached_before_threshold
//   Downstream tracing reached an outlet before finding a cell that met the
//   requested contributing-pixel threshold.
// -----------------------------------------------------------------------------
enum class lsdtt_snap_status
{
    success,
    input_out_of_bounds,
    input_invalid_flow_cell,
    threshold_not_found_in_search_window,
    outlet_reached_before_threshold
};

// -----------------------------------------------------------------------------
// snap_status_to_string
//
// The enumerator's own name, for reporting a failed snap to the user. Total
// over the enum: "unknown" is returned only if this switch has fallen out of
// step with the enum above, which is a bug rather than a value.
//
// Three CLIs (basin-extraction, basin-stats, divide-migration) each carried a
// byte-identical copy of this switch. The strings are part of the user-facing
// output of all three, so they must not be allowed to drift apart.
// -----------------------------------------------------------------------------
inline const char* snap_status_to_string(const lsdtt_snap_status status)
{
    switch (status)
    {
        case lsdtt_snap_status::success:
            return "success";
        case lsdtt_snap_status::input_out_of_bounds:
            return "input_out_of_bounds";
        case lsdtt_snap_status::input_invalid_flow_cell:
            return "input_invalid_flow_cell";
        case lsdtt_snap_status::threshold_not_found_in_search_window:
            return "threshold_not_found_in_search_window";
        case lsdtt_snap_status::outlet_reached_before_threshold:
            return "outlet_reached_before_threshold";
    }
    return "unknown";
}

// -----------------------------------------------------------------------------
// lsdtt_square_snap_options
//
// Options for expanding-square search.
//
// min_search_radius_pixels
//   Inclusive lower bound on search radius.
//   0 means the input pixel itself is tested first.
//
// max_search_radius_pixels
//   Inclusive upper bound on search radius.
//
// contributing_pixels_threshold
//   Minimum upslope contributing-pixel count required for an input to qualify.
// -----------------------------------------------------------------------------
struct lsdtt_square_snap_options
{
    std::uint64_t contributing_pixels_threshold = 1000;
    int           min_search_radius_pixels      = 0;
    int           max_search_radius_pixels      = 25;
};

// -----------------------------------------------------------------------------
// lsdtt_downstream_snap_options
//
// Options for downstream D8 snapping.
//
// contributing_pixels_threshold
//   Minimum upslope contributing-pixel count required for an input to qualify.
// -----------------------------------------------------------------------------
struct lsdtt_downstream_snap_options
{
    std::uint64_t contributing_pixels_threshold = 1000;
};

// -----------------------------------------------------------------------------
// lsdtt_largest_area_snap_options
//
// Options for largest-area snapping.
//
// contributing_pixels_threshold
//   Minimum upslope contributing-pixel count the chosen cell must have. Every
//   valid cell has at least 1 (itself), so 1 means "no minimum".
//
// max_search_radius_pixels
//   Inclusive half-width R of the search square; 0 searches only the input
//   pixel.
// -----------------------------------------------------------------------------
struct lsdtt_largest_area_snap_options
{
    std::uint64_t contributing_pixels_threshold = 1;
    int           max_search_radius_pixels      = 25;
};

// -----------------------------------------------------------------------------
// lsdtt_reported_area_snap_options
//
// Options for area-matched snapping.
//
// contributing_pixels_threshold
//   Only cells with at least this many contributing pixels (channel cells) are
//   candidates.
//
// max_search_radius_pixels
//   Inclusive half-width of the search square.
//
// area_ratio_tolerance
//   A candidate matches when 1/tol <= A_cell / A_reported <= tol, i.e.
//   |ln(A_cell / A_reported)| <= ln(tol). Must be > 1.
// -----------------------------------------------------------------------------
struct lsdtt_reported_area_snap_options
{
    std::uint64_t contributing_pixels_threshold = 1000;
    int           max_search_radius_pixels      = 10;
    double        area_ratio_tolerance          = 1.25;
};

// -----------------------------------------------------------------------------
// lsdtt_snap_result
//
// Result record returned for each attempted snap.
//
// input_* fields
//   Describe the supplied input location.
//
// input_is_valid_flow_cell
//   True when the input location corresponds to a valid flow cell.
//
// input_contributing_pixels
//   Contributing-pixel count at the input cell when the input is a valid flow
//   cell. Otherwise left as 0.
//
// snapped_* fields
//   Describe the chosen snapped location when status == success.
//
// search_radius_pixels
//   Radius at which the successful expanding-square search found its chosen
//   candidate. For downstream search, or for failures before a successful ring,
//   this is -1.
//
// downstream_steps
//   Number of receiver-to-receiver steps taken from the input cell during a
//   successful or failed downstream search.
//   0 means the input cell itself already met the threshold.
//   For expanding-square search this is -1.
// -----------------------------------------------------------------------------
struct lsdtt_snap_result
{
    lsdtt_snap_status status = lsdtt_snap_status::input_out_of_bounds;

    std::int64_t input_flat = -1;
    int          input_row  = -1;
    int          input_col  = -1;

    bool          input_is_valid_flow_cell   = false;
    std::uint64_t input_contributing_pixels  = 0;

    std::int64_t snapped_flat = -1;
    int          snapped_row  = -1;
    int          snapped_col  = -1;

    std::uint64_t snapped_contributing_pixels = 0;

    int           search_radius_pixels = -1;
    std::int64_t  downstream_steps     = -1;

    // Straight-line distance from the input point to the snapped cell centre
    // (m). Set by snap_point_to_contributing_pixels_threshold_nearest only.
    double snap_distance_metres = std::numeric_limits<double>::quiet_NaN();
};

// -----------------------------------------------------------------------------
// Single-pixel snapping: row/col inputs
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_snap_result
snap_pixel_to_contributing_pixels_threshold_expanding_square(
    const lsdtt_flow_base&          flow,
    int                             row,
    int                             col,
    const lsdtt_square_snap_options& options = lsdtt_square_snap_options{}
);

[[nodiscard]] lsdtt_snap_result
snap_pixel_to_contributing_pixels_threshold_downstream(
    const lsdtt_flow_base&               flow,
    int                                  row,
    int                                  col,
    const lsdtt_downstream_snap_options& options =
        lsdtt_downstream_snap_options{}
);

// -----------------------------------------------------------------------------
// snap_pixel_to_contributing_pixels_threshold_largest_area
//
// Largest-area snapping (method 3 above). Over the valid cells of the square
// of half-width max_search_radius_pixels around (row, col), take the cell with
// the largest contributing-pixel count; ties go to the smaller squared pixel
// distance from (row, col), then the smaller flat index. The input pixel need
// not be a valid flow cell. Status success when that cell's count is >=
// contributing_pixels_threshold (the same comparison as the other methods),
// threshold_not_found_in_search_window when it is below the threshold or the
// square holds no valid cell, input_out_of_bounds for a point off the grid.
// search_radius_pixels is the ring (Chebyshev distance) of the chosen cell;
// downstream_steps is -1. Throws std::invalid_argument for a zero threshold
// or a negative radius.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_snap_result
snap_pixel_to_contributing_pixels_threshold_largest_area(
    const lsdtt_flow_base&                 flow,
    int                                    row,
    int                                    col,
    const lsdtt_largest_area_snap_options& options =
        lsdtt_largest_area_snap_options{}
);

// -----------------------------------------------------------------------------
// snap_pixel_to_reported_drainage_area_expanding_square
//
// Area-matched snapping (method 4 above). reported_area_m2 is in the squared
// map units of the grid (m^2 for a projected DEM); a cell's area is
// contributing pixels * cell_size^2. Rule, over the valid cells of the square
// of radius max_search_radius_pixels around (row, col) whose contributing
// pixels are at least the threshold:
//   1. if any candidate has |ln(A_cell / A_rep)| <= ln(tol), take the nearest
//      of those (pixel distance), ties to the smaller |ln ratio|, then the
//      smaller flat index;
//   2. otherwise take the candidate with the smallest |ln ratio|, ties to the
//      smaller distance, then the smaller flat index, and still return
//      success: the caller compares the snapped area with the reported one
//      and flags the mismatch.
// search_radius_pixels is the ring (Chebyshev distance) of the chosen cell.
// Status threshold_not_found_in_search_window when there is no candidate,
// input_out_of_bounds for a point off the grid. Throws std::invalid_argument
// for a zero threshold, a negative radius, a tolerance <= 1 or a reported
// area that is not finite and > 0.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_snap_result
snap_pixel_to_reported_drainage_area_expanding_square(
    const lsdtt_flow_base&                   flow,
    int                                      row,
    int                                      col,
    double                                   reported_area_m2,
    const lsdtt_reported_area_snap_options&  options
);

// -----------------------------------------------------------------------------
// snap_point_to_contributing_pixels_threshold_nearest
//
// Nearest snapping (method 5 above). (x, y) is a point in the raster CRS. Over
// the valid cells with at least contributing_pixels_threshold contributing
// pixels and within Chebyshev distance max_search_radius_pixels of the point's
// pixel, take the cell whose centre is nearest the point (Euclidean); ties go
// to the larger contributing-pixel count, then the smaller flat index. The
// search is exact: rings are scanned outward and stop once a ring's nearest
// possible cell (ring - 0.5 cells away, since the point lies within its own
// pixel) is farther than the best found. min_search_radius_pixels is not used.
//
// Status: success; input_out_of_bounds when the point is off the grid;
// threshold_not_found_in_search_window when no qualifying cell is in range.
// snap_distance_metres is set on success; search_radius_pixels is the ring of
// the chosen cell; downstream_steps is -1. input_* describe the point's pixel.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_snap_result
snap_point_to_contributing_pixels_threshold_nearest(
    const lsdtt_flow_base&           flow,
    double                           x,
    double                           y,
    const lsdtt_square_snap_options& options = lsdtt_square_snap_options{}
);

// -----------------------------------------------------------------------------
// Single-pixel snapping: flat-index inputs
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_snap_result
snap_flat_to_contributing_pixels_threshold_expanding_square(
    const lsdtt_flow_base&          flow,
    std::int64_t                    flat,
    const lsdtt_square_snap_options& options = lsdtt_square_snap_options{}
);

[[nodiscard]] lsdtt_snap_result
snap_flat_to_contributing_pixels_threshold_downstream(
    const lsdtt_flow_base&               flow,
    std::int64_t                         flat,
    const lsdtt_downstream_snap_options& options =
        lsdtt_downstream_snap_options{}
);

// -----------------------------------------------------------------------------
// Batch snapping: row/col inputs
// -----------------------------------------------------------------------------

[[nodiscard]] std::vector<lsdtt_snap_result>
snap_pixels_to_contributing_pixels_threshold_expanding_square(
    const lsdtt_flow_base&                  flow,
    const std::vector<lsdtt_pixel_rc>&      pixels,
    const lsdtt_square_snap_options&        options = lsdtt_square_snap_options{}
);

[[nodiscard]] std::vector<lsdtt_snap_result>
snap_pixels_to_contributing_pixels_threshold_downstream(
    const lsdtt_flow_base&                  flow,
    const std::vector<lsdtt_pixel_rc>&      pixels,
    const lsdtt_downstream_snap_options&    options =
        lsdtt_downstream_snap_options{}
);

// -----------------------------------------------------------------------------
// Batch snapping: flat-index inputs
// -----------------------------------------------------------------------------

[[nodiscard]] std::vector<lsdtt_snap_result>
snap_flats_to_contributing_pixels_threshold_expanding_square(
    const lsdtt_flow_base&                  flow,
    const std::vector<std::int64_t>&        flats,
    const lsdtt_square_snap_options&        options = lsdtt_square_snap_options{}
);

[[nodiscard]] std::vector<lsdtt_snap_result>
snap_flats_to_contributing_pixels_threshold_downstream(
    const lsdtt_flow_base&                  flow,
    const std::vector<std::int64_t>&        flats,
    const lsdtt_downstream_snap_options&    options =
        lsdtt_downstream_snap_options{}
);

// -----------------------------------------------------------------------------
// snap_pixel_by_method
//
// Dispatch to the snapping routine named by `method`, so a caller that takes
// the method as a string parameter does not have to spell the if/else out.
//
//   "none"              no snapping; the input pixel is returned with
//                       status success, so a caller can treat all the
//                       methods uniformly
//   "expanding_square"  snap_pixel_to_contributing_pixels_threshold_expanding_square
//   "downstream"        snap_pixel_to_contributing_pixels_threshold_downstream
//   "largest_area"      snap_pixel_to_contributing_pixels_threshold_largest_area,
//                       taking contributing_pixels_threshold and
//                       max_search_radius_pixels from square_options
//                       (min_search_radius_pixels is not used)
//
// Throws std::invalid_argument on an unrecognised method, rather than silently
// falling through to no snapping: a typo in a parameter file would otherwise
// produce unsnapped outlets that look like a modelling choice.
//
// -----------------------------------------------------------------------------
// "reported_area" needs the reported area, so this overload throws
// std::invalid_argument for it and names the overload below.
[[nodiscard]] lsdtt_snap_result snap_pixel_by_method(
    const lsdtt_flow_base&               flow,
    int                                  row,
    int                                  col,
    const std::string&                   method,
    const lsdtt_square_snap_options&     square_options,
    const lsdtt_downstream_snap_options& downstream_options
);

// The same, with "reported_area" too:
//   "reported_area"     snap_pixel_to_reported_drainage_area_expanding_square
//                       with reported_area_m2 and area_options
// The other methods ignore reported_area_m2 and area_options.
[[nodiscard]] lsdtt_snap_result snap_pixel_by_method(
    const lsdtt_flow_base&                  flow,
    int                                     row,
    int                                     col,
    const std::string&                      method,
    const lsdtt_square_snap_options&        square_options,
    const lsdtt_downstream_snap_options&    downstream_options,
    double                                  reported_area_m2,
    const lsdtt_reported_area_snap_options& area_options
);

} // namespace lsdtt
