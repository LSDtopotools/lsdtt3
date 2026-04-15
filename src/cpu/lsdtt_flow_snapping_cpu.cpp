// lsdtt_flow_snapping_cpu.cpp
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
// Implementation of flow-based pixel snapping helpers.
//
// These routines snap user-supplied raster pixels to cells that satisfy a
// contributing-pixel threshold based on the flow-routing s_vector:
//
//   contributing_pixels = flow.s_vector_size(flat)
//
// Two snapping behaviours are implemented:
//
//   1. expanding_square
//      Search outwards ring-by-ring in raster pixel space. At the first radius
//      that contains one or more qualifying cells, choose the candidate with:
//        a. smallest Euclidean distance in pixel units
//        b. if tied, largest contributing-pixel count
//        c. if still tied, smallest flat index
//
//   2. downstream_d8
//      Follow the D8 receiver chain until the first cell is reached whose
//      contributing-pixel count meets or exceeds the threshold.
//
// Notes
// -----
// - These routines require a flow object built with s_vector enabled.
// - Expanding-square search can begin from any in-bounds raster pixel.
// - Downstream search requires the starting pixel to be a valid flow cell.
// -----------------------------------------------------------------------------

#include "lsdtt_flow_snapping_cpu.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// Validation helpers
// -----------------------------------------------------------------------------
void validate_common_preconditions(
    const lsdtt_flow_base& flow,
    const std::string& function_name)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            function_name + ": flow object is empty"
        );
    }

    if (!flow.has_s_vector())
    {
        throw std::runtime_error(
            function_name + ": flow object must have s_vector enabled because "
            "snapping thresholds are based on contributing pixels "
            "(flow.s_vector_size)"
        );
    }
}

void validate_square_snap_options(
    const lsdtt_square_snap_options& options,
    const std::string& function_name)
{
    if (options.contributing_pixels_threshold == 0)
    {
        throw std::invalid_argument(
            function_name + ": contributing_pixels_threshold must be > 0"
        );
    }

    if (options.min_search_radius_pixels < 0)
    {
        throw std::invalid_argument(
            function_name + ": min_search_radius_pixels must be >= 0"
        );
    }

    if (options.max_search_radius_pixels < 0)
    {
        throw std::invalid_argument(
            function_name + ": max_search_radius_pixels must be >= 0"
        );
    }

    if (options.max_search_radius_pixels < options.min_search_radius_pixels)
    {
        throw std::invalid_argument(
            function_name + ": max_search_radius_pixels must be >= "
            "min_search_radius_pixels"
        );
    }
}

void validate_downstream_snap_options(
    const lsdtt_downstream_snap_options& options,
    const std::string& function_name)
{
    if (options.contributing_pixels_threshold == 0)
    {
        throw std::invalid_argument(
            function_name + ": contributing_pixels_threshold must be > 0"
        );
    }
}

// -----------------------------------------------------------------------------
// Input initialisation helpers
// -----------------------------------------------------------------------------
lsdtt_snap_result initialise_snap_result_from_row_col(
    const lsdtt_flow_base& flow,
    int row,
    int col)
{
    lsdtt_snap_result result;
    result.input_row = row;
    result.input_col = col;

    if (!flow.in_bounds(row, col))
    {
        result.status = lsdtt_snap_status::input_out_of_bounds;
        return result;
    }

    // In bounds: allow snapping routines to proceed.
    result.status    = lsdtt_snap_status::success;
    result.input_flat = flow.flat_index_unchecked(row, col);

    if (flow.is_valid_cell(result.input_flat))
    {
        result.input_is_valid_flow_cell  = true;
        result.input_contributing_pixels =
            flow.s_vector_size_unchecked(result.input_flat);
    }

    return result;
}

lsdtt_snap_result initialise_snap_result_from_flat(
    const lsdtt_flow_base& flow,
    std::int64_t flat)
{
    lsdtt_snap_result result;
    result.input_flat = flat;

    if (!flow.is_flat_in_range(flat))
    {
        result.status = lsdtt_snap_status::input_out_of_bounds;
        return result;
    }

    // In range: allow snapping routines to proceed.
    result.status    = lsdtt_snap_status::success;
    result.input_row = flow.row_from_flat_unchecked(flat);
    result.input_col = flow.col_from_flat_unchecked(flat);

    if (flow.is_valid_cell(flat))
    {
        result.input_is_valid_flow_cell  = true;
        result.input_contributing_pixels =
            flow.s_vector_size_unchecked(flat);
    }

    return result;
}


// -----------------------------------------------------------------------------
// Success helper
// -----------------------------------------------------------------------------
void write_successful_snap_to_result(
    const lsdtt_flow_base& flow,
    std::int64_t snapped_flat,
    std::uint64_t snapped_contributing_pixels,
    lsdtt_snap_result& result)
{
    result.status                     = lsdtt_snap_status::success;
    result.snapped_flat               = snapped_flat;
    result.snapped_row                = flow.row_from_flat_unchecked(snapped_flat);
    result.snapped_col                = flow.col_from_flat_unchecked(snapped_flat);
    result.snapped_contributing_pixels = snapped_contributing_pixels;
}

// -----------------------------------------------------------------------------
// Expanding-square candidate helpers
// -----------------------------------------------------------------------------
struct square_snap_candidate
{
    bool         found                = false;
    std::int64_t flat                 = -1;
    int          row                  = -1;
    int          col                  = -1;
    std::uint64_t contributing_pixels = 0;
    std::int64_t distance_squared     = std::numeric_limits<std::int64_t>::max();
};

bool candidate_is_better_than_current_best(
    const square_snap_candidate& candidate,
    const square_snap_candidate& current_best)
{
    if (!candidate.found)
    {
        return false;
    }

    if (!current_best.found)
    {
        return true;
    }

    if (candidate.distance_squared < current_best.distance_squared)
    {
        return true;
    }

    if (candidate.distance_squared > current_best.distance_squared)
    {
        return false;
    }

    if (candidate.contributing_pixels > current_best.contributing_pixels)
    {
        return true;
    }

    if (candidate.contributing_pixels < current_best.contributing_pixels)
    {
        return false;
    }

    return candidate.flat < current_best.flat;
}

void consider_square_search_candidate(
    const lsdtt_flow_base& flow,
    int input_row,
    int input_col,
    int candidate_row,
    int candidate_col,
    std::uint64_t threshold_pixels,
    square_snap_candidate& current_best)
{
    if (!flow.in_bounds(candidate_row, candidate_col))
    {
        return;
    }

    const std::int64_t candidate_flat =
        flow.flat_index_unchecked(candidate_row, candidate_col);

    if (!flow.is_valid_cell(candidate_flat))
    {
        return;
    }

    const std::uint64_t candidate_contributing_pixels =
        flow.s_vector_size_unchecked(candidate_flat);

    if (candidate_contributing_pixels < threshold_pixels)
    {
        return;
    }

    const std::int64_t dr =
        static_cast<std::int64_t>(candidate_row) -
        static_cast<std::int64_t>(input_row);

    const std::int64_t dc =
        static_cast<std::int64_t>(candidate_col) -
        static_cast<std::int64_t>(input_col);

    square_snap_candidate candidate;
    candidate.found                = true;
    candidate.flat                 = candidate_flat;
    candidate.row                  = candidate_row;
    candidate.col                  = candidate_col;
    candidate.contributing_pixels  = candidate_contributing_pixels;
    candidate.distance_squared     = dr * dr + dc * dc;

    if (candidate_is_better_than_current_best(candidate, current_best))
    {
        current_best = candidate;
    }
}

square_snap_candidate find_best_candidate_on_square_search_ring(
    const lsdtt_flow_base& flow,
    int input_row,
    int input_col,
    int radius_pixels,
    std::uint64_t threshold_pixels)
{
    square_snap_candidate best_candidate;

    if (radius_pixels == 0)
    {
        consider_square_search_candidate(
            flow,
            input_row,
            input_col,
            input_row,
            input_col,
            threshold_pixels,
            best_candidate
        );

        return best_candidate;
    }

    const int row_min = input_row - radius_pixels;
    const int row_max = input_row + radius_pixels;
    const int col_min = input_col - radius_pixels;
    const int col_max = input_col + radius_pixels;

    // Top row
    for (int col = col_min; col <= col_max; ++col)
    {
        consider_square_search_candidate(
            flow,
            input_row,
            input_col,
            row_min,
            col,
            threshold_pixels,
            best_candidate
        );
    }

    // Bottom row
    if (row_max != row_min)
    {
        for (int col = col_min; col <= col_max; ++col)
        {
            consider_square_search_candidate(
                flow,
                input_row,
                input_col,
                row_max,
                col,
                threshold_pixels,
                best_candidate
            );
        }
    }

    // Left column, excluding corners already visited
    for (int row = row_min + 1; row <= row_max - 1; ++row)
    {
        consider_square_search_candidate(
            flow,
            input_row,
            input_col,
            row,
            col_min,
            threshold_pixels,
            best_candidate
        );
    }

    // Right column, excluding corners already visited
    if (col_max != col_min)
    {
        for (int row = row_min + 1; row <= row_max - 1; ++row)
        {
            consider_square_search_candidate(
                flow,
                input_row,
                input_col,
                row,
                col_max,
                threshold_pixels,
                best_candidate
            );
        }
    }

    return best_candidate;
}

// -----------------------------------------------------------------------------
// Core snapping implementations
// These assume that common/global validation has already been done.
// -----------------------------------------------------------------------------
lsdtt_snap_result snap_preinitialised_input_to_contributing_pixels_threshold_expanding_square(
    const lsdtt_flow_base& flow,
    lsdtt_snap_result result,
    const lsdtt_square_snap_options& options)
{
    if (result.status == lsdtt_snap_status::input_out_of_bounds)
    {
        return result;
    }

    for (int radius = options.min_search_radius_pixels;
         radius <= options.max_search_radius_pixels;
         ++radius)
    {
        const square_snap_candidate best_candidate =
            find_best_candidate_on_square_search_ring(
                flow,
                result.input_row,
                result.input_col,
                radius,
                options.contributing_pixels_threshold
            );

        if (!best_candidate.found)
        {
            continue;
        }

        write_successful_snap_to_result(
            flow,
            best_candidate.flat,
            best_candidate.contributing_pixels,
            result
        );

        result.search_radius_pixels = radius;
        result.downstream_steps     = -1;
        return result;
    }

    result.status = lsdtt_snap_status::threshold_not_found_in_search_window;
    return result;
}

lsdtt_snap_result snap_preinitialised_input_to_contributing_pixels_threshold_downstream(
    const lsdtt_flow_base& flow,
    lsdtt_snap_result result,
    const lsdtt_downstream_snap_options& options)
{
    if (result.status == lsdtt_snap_status::input_out_of_bounds)
    {
        return result;
    }

    if (!result.input_is_valid_flow_cell)
    {
        result.status = lsdtt_snap_status::input_invalid_flow_cell;
        return result;
    }

    std::int64_t current_flat = result.input_flat;
    std::int64_t downstream_steps = 0;
    const std::int64_t max_steps =
        static_cast<std::int64_t>(flow.full_grid_size());

    while (true)
    {
        const std::uint64_t current_contributing_pixels =
            flow.s_vector_size_unchecked(current_flat);

        if (current_contributing_pixels >=
            options.contributing_pixels_threshold)
        {
            write_successful_snap_to_result(
                flow,
                current_flat,
                current_contributing_pixels,
                result
            );

            result.search_radius_pixels = -1;
            result.downstream_steps     = downstream_steps;
            return result;
        }

        if (flow.is_outlet_unchecked(current_flat))
        {
            result.status           =
                lsdtt_snap_status::outlet_reached_before_threshold;
            result.search_radius_pixels = -1;
            result.downstream_steps     = downstream_steps;
            return result;
        }

        const std::int64_t receiver_flat =
            flow.receiver_of_unchecked(current_flat);

        if (receiver_flat < 0)
        {
            throw std::runtime_error(
                "snap_preinitialised_input_to_contributing_pixels_threshold_downstream: "
                "encountered a negative non-outlet receiver while tracing "
                "downstream"
            );
        }

        if (!flow.is_valid_cell(receiver_flat))
        {
            throw std::runtime_error(
                "snap_preinitialised_input_to_contributing_pixels_threshold_downstream: "
                "encountered an invalid receiver cell while tracing downstream"
            );
        }

        current_flat = receiver_flat;
        ++downstream_steps;

        if (downstream_steps > max_steps)
        {
            throw std::runtime_error(
                "snap_preinitialised_input_to_contributing_pixels_threshold_downstream: "
                "exceeded downstream safety step limit; the routing graph may "
                "be cyclic or inconsistent"
            );
        }
    }
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// Single-pixel snapping: row/col inputs
// -----------------------------------------------------------------------------
lsdtt_snap_result
snap_pixel_to_contributing_pixels_threshold_expanding_square(
    const lsdtt_flow_base& flow,
    int row,
    int col,
    const lsdtt_square_snap_options& options)
{
    validate_common_preconditions(
        flow,
        "snap_pixel_to_contributing_pixels_threshold_expanding_square"
    );

    validate_square_snap_options(
        options,
        "snap_pixel_to_contributing_pixels_threshold_expanding_square"
    );

    const lsdtt_snap_result result =
        initialise_snap_result_from_row_col(flow, row, col);

    return snap_preinitialised_input_to_contributing_pixels_threshold_expanding_square(
        flow,
        result,
        options
    );
}

lsdtt_snap_result
snap_pixel_to_contributing_pixels_threshold_downstream(
    const lsdtt_flow_base& flow,
    int row,
    int col,
    const lsdtt_downstream_snap_options& options)
{
    validate_common_preconditions(
        flow,
        "snap_pixel_to_contributing_pixels_threshold_downstream"
    );

    validate_downstream_snap_options(
        options,
        "snap_pixel_to_contributing_pixels_threshold_downstream"
    );

    const lsdtt_snap_result result =
        initialise_snap_result_from_row_col(flow, row, col);

    return snap_preinitialised_input_to_contributing_pixels_threshold_downstream(
        flow,
        result,
        options
    );
}

// -----------------------------------------------------------------------------
// Single-pixel snapping: flat-index inputs
// -----------------------------------------------------------------------------
lsdtt_snap_result
snap_flat_to_contributing_pixels_threshold_expanding_square(
    const lsdtt_flow_base& flow,
    std::int64_t flat,
    const lsdtt_square_snap_options& options)
{
    validate_common_preconditions(
        flow,
        "snap_flat_to_contributing_pixels_threshold_expanding_square"
    );

    validate_square_snap_options(
        options,
        "snap_flat_to_contributing_pixels_threshold_expanding_square"
    );

    const lsdtt_snap_result result =
        initialise_snap_result_from_flat(flow, flat);

    return snap_preinitialised_input_to_contributing_pixels_threshold_expanding_square(
        flow,
        result,
        options
    );
}

lsdtt_snap_result
snap_flat_to_contributing_pixels_threshold_downstream(
    const lsdtt_flow_base& flow,
    std::int64_t flat,
    const lsdtt_downstream_snap_options& options)
{
    validate_common_preconditions(
        flow,
        "snap_flat_to_contributing_pixels_threshold_downstream"
    );

    validate_downstream_snap_options(
        options,
        "snap_flat_to_contributing_pixels_threshold_downstream"
    );

    const lsdtt_snap_result result =
        initialise_snap_result_from_flat(flow, flat);

    return snap_preinitialised_input_to_contributing_pixels_threshold_downstream(
        flow,
        result,
        options
    );
}

// -----------------------------------------------------------------------------
// Batch snapping: row/col inputs
// -----------------------------------------------------------------------------
std::vector<lsdtt_snap_result>
snap_pixels_to_contributing_pixels_threshold_expanding_square(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_pixel_rc>& pixels,
    const lsdtt_square_snap_options& options)
{
    validate_common_preconditions(
        flow,
        "snap_pixels_to_contributing_pixels_threshold_expanding_square"
    );

    validate_square_snap_options(
        options,
        "snap_pixels_to_contributing_pixels_threshold_expanding_square"
    );

    std::vector<lsdtt_snap_result> results;
    results.reserve(pixels.size());

    for (const auto& pixel : pixels)
    {
        const lsdtt_snap_result result =
            initialise_snap_result_from_row_col(
                flow,
                pixel.row,
                pixel.col
            );

        results.push_back(
            snap_preinitialised_input_to_contributing_pixels_threshold_expanding_square(
                flow,
                result,
                options
            )
        );
    }

    return results;
}

std::vector<lsdtt_snap_result>
snap_pixels_to_contributing_pixels_threshold_downstream(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_pixel_rc>& pixels,
    const lsdtt_downstream_snap_options& options)
{
    validate_common_preconditions(
        flow,
        "snap_pixels_to_contributing_pixels_threshold_downstream"
    );

    validate_downstream_snap_options(
        options,
        "snap_pixels_to_contributing_pixels_threshold_downstream"
    );

    std::vector<lsdtt_snap_result> results;
    results.reserve(pixels.size());

    for (const auto& pixel : pixels)
    {
        const lsdtt_snap_result result =
            initialise_snap_result_from_row_col(
                flow,
                pixel.row,
                pixel.col
            );

        results.push_back(
            snap_preinitialised_input_to_contributing_pixels_threshold_downstream(
                flow,
                result,
                options
            )
        );
    }

    return results;
}

// -----------------------------------------------------------------------------
// Batch snapping: flat-index inputs
// -----------------------------------------------------------------------------
std::vector<lsdtt_snap_result>
snap_flats_to_contributing_pixels_threshold_expanding_square(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& flats,
    const lsdtt_square_snap_options& options)
{
    validate_common_preconditions(
        flow,
        "snap_flats_to_contributing_pixels_threshold_expanding_square"
    );

    validate_square_snap_options(
        options,
        "snap_flats_to_contributing_pixels_threshold_expanding_square"
    );

    std::vector<lsdtt_snap_result> results;
    results.reserve(flats.size());

    for (const std::int64_t flat : flats)
    {
        const lsdtt_snap_result result =
            initialise_snap_result_from_flat(flow, flat);

        results.push_back(
            snap_preinitialised_input_to_contributing_pixels_threshold_expanding_square(
                flow,
                result,
                options
            )
        );
    }

    return results;
}

std::vector<lsdtt_snap_result>
snap_flats_to_contributing_pixels_threshold_downstream(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& flats,
    const lsdtt_downstream_snap_options& options)
{
    validate_common_preconditions(
        flow,
        "snap_flats_to_contributing_pixels_threshold_downstream"
    );

    validate_downstream_snap_options(
        options,
        "snap_flats_to_contributing_pixels_threshold_downstream"
    );

    std::vector<lsdtt_snap_result> results;
    results.reserve(flats.size());

    for (const std::int64_t flat : flats)
    {
        const lsdtt_snap_result result =
            initialise_snap_result_from_flat(flow, flat);

        results.push_back(
            snap_preinitialised_input_to_contributing_pixels_threshold_downstream(
                flow,
                result,
                options
            )
        );
    }

    return results;
}

} // namespace lsdtt
