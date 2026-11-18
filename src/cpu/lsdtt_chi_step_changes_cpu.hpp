// lsdtt_chi_step_changes_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Step changes (waterfalls, dams, steep gorge reaches) on the branches of a
// fitted chi network (lsdtt-chi-analysis, run_step_change_detection).
//
// A step is found on the longitudinal profile, not in chi: the gradient
// between adjacent channel pixels is a physical slope, independent of m/n and
// A0. On each branch path of the fit (the same nodes and D8 receivers, the
// fit's observed elevations), from the branch base (the basin outlet on the
// main stem, the junction on the parent for a tributary) to its head:
//
//   1. PAIRS. Each node and its receiver form a pair: drop dz = z_node -
//      z_receiver, length L = the D8 step (cell size, x sqrt 2 on a diagonal).
//   2. WINDOWS. Every run of consecutive pairs whose flow length is at most
//      max_window_length_metres (always at least one pair, so a pixel coarser
//      than the window still has single-pair windows) is a window, with
//      gradient sum(dz) / sum(L). A window is steep if that is at least
//      min_gradient. Every window length up to the maximum is used, so a
//      one-pixel waterfall and a step smeared over several pixels are both
//      found; flat or uphill pairs inside a window lower its gradient
//      instead of breaking it.
//   3. REACHES. Every pair in a steep window is marked; consecutive marked
//      pairs form a reach. Its end pairs whose own gradient is below
//      min_gradient are trimmed (a window holding one very steep pair marks
//      its gentle neighbours); gentle pairs inside the reach stay.
//   4. MINIMUM DROP. A reach is kept if its drop (top - bottom elevation) is
//      at least min_drop_metres, so lidar bed roughness is not reported.
//
// Detection does not depend on m/n; only the chi of a reach's ends does.
// The threshold is resolution dependent: at 30 m, 0.40 also marks steep gorge
// reaches, which is why each reach carries its drop, length, gradients and
// flags.
//
// The CSV writer uses std::ofstream; the line writer (write_chi_steps_vector)
// goes through lsdtt_data_loader.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_chi_data_maps_cpu.hpp"
#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_chi_knickpoints_cpu.hpp"
#include "lsdtt_chi_network_analysis_cpu.hpp"

namespace lsdtt
{

struct lsdtt_chi_step_options
{
    double min_gradient = 0.40;                // m/m
    double max_window_length_metres = 90.0;    // flow length
    double min_drop_metres = 5.0;
};

// A reach on one path: pairs first_pair..last_pair (pair 0 is the lowest,
// at the path's base).
struct lsdtt_chi_step_reach_on_path
{
    std::size_t first_pair = 0;
    std::size_t last_pair = 0;
    double drop_metres = 0.0;            // sum of dz over the reach
    double flow_length_metres = 0.0;     // sum of L
    double mean_gradient = 0.0;          // drop / length
    double max_window_gradient = 0.0;    // steepest window inside the reach
    std::int32_t n_pairs = 0;
    bool contains_uphill_pairs = false;  // a pair with dz <= 0
};

// Steps 2-4 on one path. dz[i] and length_metres[i] are pair i, ordered
// upstream from the base; lengths must be > 0. Throws std::invalid_argument
// on bad options (min_gradient <= 0, max_window_length_metres <= 0,
// min_drop_metres < 0) or mismatched sizes.
[[nodiscard]] std::vector<lsdtt_chi_step_reach_on_path> detect_chi_step_reaches(
    const std::vector<double>&    dz,
    const std::vector<double>&    length_metres,
    const lsdtt_chi_step_options& options
);

inline constexpr double k_chi_step_nan = std::numeric_limits<double>::quiet_NaN();

struct lsdtt_chi_step_row
{
    std::int32_t step_id = -1;
    std::int32_t basin_index = -1;
    std::int32_t basin_key = -1;
    std::int32_t source_key = -1;            // of the branch's head
    std::int32_t branch_id = -1;
    std::int64_t bottom_flat = -1, top_flat = -1;
    double bottom_x = k_chi_step_nan, bottom_y = k_chi_step_nan;
    double top_x = k_chi_step_nan, top_y = k_chi_step_nan;
    double bottom_elevation = k_chi_step_nan, top_elevation = k_chi_step_nan;
    double drop_m = k_chi_step_nan;
    double flow_length_m = k_chi_step_nan;
    double mean_gradient = k_chi_step_nan;
    double max_window_gradient = k_chi_step_nan;
    std::int32_t n_pairs = 0;
    double bottom_chi = k_chi_step_nan, top_chi = k_chi_step_nan;
    double bottom_flow_distance_m = k_chi_step_nan;     // from the chi data maps keys
    double bottom_drainage_area_m2 = k_chi_step_nan;
    double distance_to_branch_base_m = k_chi_step_nan;  // along the path
    std::uint8_t at_tributary_mouth = 0;     // tributary, bottom within the max window of the junction
    std::uint8_t at_basin_outlet = 0;        // main stem, bottom within the max window of the outlet
    std::uint8_t contains_uphill_pairs = 0;
};

// All step reaches of the fitted network. fits.basins[i] must have been
// fitted from chi_networks.basins[i]; keys from assign_chi_source_keys (may be
// empty: then source_key is -1 and flow distance / area are NaN).
[[nodiscard]] std::vector<lsdtt_chi_step_row> build_chi_step_changes(
    const lsdtt_chi_network_fit_collection_result&  fits,
    const lsdtt_chi_collection_result&              chi_networks,
    const std::vector<lsdtt_chi_basin_source_keys>& keys,
    const lsdtt_raster_meta&                        meta,
    const lsdtt_chi_step_options&                   options
);

void write_chi_steps_csv(const std::string& filename, const std::vector<lsdtt_chi_step_row>& rows);

// One polyline per step, traced node by node along the fit's D8 receivers
// from its top to its bottom (downstream, as the fit segment lines), with the
// CSV's attributes except the four end coordinates (the geometry carries
// them). fits must be the fit the steps were built from; a step whose trace
// leaves its basin is written as its two end points. fmt must be a binary
// vector format (flatgeobuf or gpkg).
void write_chi_steps_vector(
    const std::string&                              filename,
    const std::vector<lsdtt_chi_step_row>&          steps,
    const lsdtt_chi_network_fit_collection_result&  fits,
    const lsdtt_raster_meta&                        meta,
    vector_output_format                            fmt,
    bool                                            spatial_index = false
);

// PART C. The step reaches as the step-aware fit takes them
// (lsdtt_chi_network_fit_options::step_reaches): (bottom_flat, top_flat).
[[nodiscard]] std::vector<std::pair<std::int64_t, std::int64_t>> chi_step_reaches_for_fit(
    const std::vector<lsdtt_chi_step_row>& steps
);

// Fill the stage 1 knickpoints' step columns (lsdtt_chi_knickpoint_row):
//   nearest_step_id, nearest_step_chi_distance: the nearest step on the row's
//     branch, by chi (0 if the boundary lies inside it); -1 / NaN if none;
//   step_in_adjacent_segments: 1 if a step on the row's branch, or on the
//     branch of its downstream segment (the parent, at a junction), overlaps
//     the chi range of the two segments either side of the boundary.
// fits must be the fit the knickpoints were built from (stage 1).
void cross_reference_chi_knickpoints_with_steps(
    std::vector<lsdtt_chi_knickpoint_row>&          knickpoints,
    const std::vector<lsdtt_chi_step_row>&          steps,
    const lsdtt_chi_network_fit_collection_result&  fits
);

} // namespace lsdtt
