// lsdtt_divide_migration_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "../utils/lsdtt_raster_geometry.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Defined in lsdtt_chi_extraction_cpu.hpp; only referenced here.
struct lsdtt_chi_collection_result;

struct lsdtt_outlet_record
{
    int point_id = -1;
    int row = -1;
    int col = -1;
    std::int64_t input_flat = -1;
    std::int64_t snapped_flat = -1;
    bool snapped = false;

    int unit_id = -1;
    std::string unit_key = "";
    std::string cluster_value = "";
};

struct lsdtt_unit_record
{
    int unit_id = -1;
    std::string unit_key = "";
    std::vector<int> point_ids;
    std::vector<std::int64_t> snapped_outlet_flats;
};

struct lsdtt_divide_segment_record
{
    std::string pair_id = "";
    int unit_a = -1;
    int unit_b = -1;
    std::int32_t basin_junction_a = -1;
    std::int32_t basin_junction_b = -1;
    int segment_id = -1;
    std::vector<lsdtt_xy> points;
    double length_m = 0.0;
};

struct lsdtt_channel_head_metric_record
{
    int head_id = -1;
    int outlet_id = -1;
    int unit_id = -1;
    std::string pair_id_nearest = "";

    double head_x = 0.0;
    double head_y = 0.0;
    double head_elevation_m = 0.0;
    double chi = 0.0;

    double nearest_divide_x = 0.0;
    double nearest_divide_y = 0.0;
    double divide_elevation_m = 0.0;
    double distance_to_divide_m = 0.0;
    double relief_m = 0.0;
    double gradient_to_divide = 0.0;
};

struct lsdtt_pair_summary_record
{
    std::string pair_id = "";
    int unit_a = -1;
    int unit_b = -1;

    std::int64_t n_heads_a = 0;
    std::int64_t n_heads_b = 0;

    double mean_head_elevation_a = 0.0;
    double mean_head_elevation_b = 0.0;
    double delta_head_elevation = 0.0;

    double mean_gradient_a = 0.0;
    double mean_gradient_b = 0.0;
    double delta_gradient = 0.0;

    double mean_relief_a = 0.0;
    double mean_relief_b = 0.0;
    double delta_relief = 0.0;

    double mean_chi_a = 0.0;
    double mean_chi_b = 0.0;
    double delta_chi = 0.0;

    std::string uncertainty_method = "";
    double uncertainty_head_elevation = 0.0;
    double uncertainty_gradient = 0.0;
    double uncertainty_relief = 0.0;
    double uncertainty_chi = 0.0;
};

// -----------------------------------------------------------------------------
// compute_nearest_divide_metrics
//
// For every channel head, find the nearest location on the divide geometry and
// fill the nearest-divide fields of the record:
//   pair_id_nearest, nearest_divide_x, nearest_divide_y, divide_elevation_m,
//   distance_to_divide_m, relief_m, gradient_to_divide.
//
// Rules:
//   - only divide segments whose pair involves the head's own unit_id are
//     considered, so a head is never matched to a divide between two other
//     units. A head whose unit takes part in no divide keeps its NaN fields
//     and an empty pair_id_nearest,
//   - the nearest location is the closest point on a divide *segment*
//     (point-to-line-segment projection), not merely the closest vertex,
//   - divide_elevation_m is sampled from elevation (normally the filled/routed
//     DEM) at the cell containing the nearest point. Out-of-bounds or nodata
//     cells leave the elevation-derived fields as NaN,
//   - gradient_to_divide is NaN when the head sits on (or within rounding of)
//     the divide, so the division by zero never happens. Negative relief is a
//     legitimate result and is passed through unchanged.
//
// The search uses a uniform-grid index built per unit, so cost is close to
// linear in the number of heads for well-distributed divides.
//
// Returns the number of heads that were matched to a divide.
// -----------------------------------------------------------------------------
[[nodiscard]] std::size_t compute_nearest_divide_metrics(
    const std::vector<lsdtt_divide_segment_record>& divide_segments,
    const lsdtt_raster_cpu&                         elevation,
    std::vector<lsdtt_channel_head_metric_record>&  head_metrics);

// -----------------------------------------------------------------------------
// apply_common_base_level_to_chi
//
// Impose a common base level on an outlet-seeded chi collection, mirroring the
// intent of DivideTools' SetOutlet.m: every chi network is truncated below
// base_elevation_m and chi is re-zeroed at the lowest retained node of each
// remaining channel fragment. Comparisons across basins are then made from the
// same elevation datum rather than from each basin's own outlet.
//
// Elevations are read from elevation (this must be the filled/routed DEM so
// that flow paths are monotonic and the retained node set stays upslope
// closed). basin.outlet_flat is left untouched so caller-side outlet -> unit
// maps remain valid.
//
// Returns the number of chi nodes removed.
// -----------------------------------------------------------------------------
std::size_t apply_common_base_level_to_chi(
    lsdtt_chi_collection_result& chi,
    const lsdtt_raster_cpu&      elevation,
    double                       base_elevation_m);

} // namespace lsdtt
