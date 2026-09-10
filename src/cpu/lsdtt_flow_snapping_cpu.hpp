// lsdtt_flow_snapping_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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
// Two snapping behaviours are currently supported:
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

#include <cstdint>
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
//   radius without finding any qualifying pixel.
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

} // namespace lsdtt
