// lsdtt_hillslope_channel_tagging_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Links every hillslope pixel to a channel pixel and a ridge pixel, and reports
// the relief between them and the pixel's relative height.
//
// THE CHANNEL LINK IS ALWAYS THE D8 PATH. A hillslope pixel is linked to the
// first channel pixel on its D8 path, so HAND (Height Above the Nearest
// Drainage; Nobre et al., Journal of Hydrology 404, 13-29, 2011) does not
// depend on how the ridge is chosen.
//
//   ridge_relief(p)    = z[ridge(p)] - z[channel(p)]
//   relative_height(p) = HAND(p) / ridge_relief(p)
//
// WHY THE RIDGE IS CHOSEN PER GROUP. D8 flow converges, so about half of the
// hillslope pixels of a typical DEM receive no flow from any neighbour, and
// most of those are part-way down a slope. If each pixel used its own farthest
// upslope pixel, each of them would be its own ridge with a relative height of
// exactly 1. Pixels therefore share a ridge within a group:
//
//   channel_pixel  all pixels draining to the same channel pixel (default)
//   entry_pixel    all pixels entering the channel through the same last
//                  hillslope pixel, which separates the slopes joining a
//                  channel pixel
//   pixel          the pixel and those draining through it (the literal
//                  per-pixel definition)
//
// THE RIDGE METHODS
//
//   farthest  the group's pixel with the longest D8 path to the channel. Ties:
//             higher elevation, then lower flat index.
//   highest   the group's highest pixel. Ties: longer D8 path, then lower flat
//             index.
//   aspect_ray  a straight ray from the pixel centre, upslope along the pixel's
//             aspect (not per group). See ASPECT RAYS below.
//
// RELATIVE HEIGHT CAN EXCEED 1 WITH farthest AND A SHARED GROUP. The farthest
// pixel need not be the highest: a pixel on a shorter, steeper branch of the
// same group can stand above it. The value is reported as it is, not clipped,
// and those pixels are counted in n_above_ridge. With highest it lies in [0, 1]
// on a surface where elevation does not rise downstream.
//
// CHANNELS STOP THE SEARCH. A ridge search never passes through a channel
// pixel, so a fragment of an ingested network upslope of a hillslope pixel is
// that fragment's own hillslope's limit.
//
// THE DIVIDE NETWORK. Each channel pixel belongs to a junction-to-junction link
// (a junction belongs to the link below it), each hillslope pixel to the link
// its D8 path reaches, and a pixel reaching no channel to its outlet. A divide
// edge is a pixel edge shared by two valid non-channel pixels with different
// labels. Nothing is filtered. Divides lie BETWEEN pixels. Edges against nodata
// or the grid boundary are not divides. No outline is traced, so a basin joined
// to itself only through diagonal contacts on a jagged coast is handled.
//
// ASPECT RAYS. The ray visits every cell it touches (Amanatides and Woo 1987),
// so it cannot slip through an 8-connected channel or divide. It stops at the
// first divide edge satisfying both:
//
//   at or above the start   the edge's higher pixel is not below the start pixel
//   a crest along the ray   the next cell beyond the crossing is not higher than
//                           that pixel (nodata, off-grid or channel: not higher)
//
// Link catchments also meet along flow-parallel lines running upslope from
// junctions, where the ground keeps rising past the divide. The crest test
// passes those. The ridge is the edge's higher pixel. If the ray reaches a
// channel, nodata or the grid boundary first, the ridge is the highest cell on
// the ray, including the start pixel. A pixel with no aspect has no ridge.
// ridge_distance is then the straight-line distance from pixel to ridge.
//
// Aspect is the DOWNSLOPE compass bearing relative to GRID north (0 = up the
// grid, 90 = increasing column). Periodic boundaries are refused.
//
// Heights are taken from the elevation passed in, which should be the routing
// surface. Undefined values are nodata, never 0, and are counted.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "../core/lsdtt_channel_attribute_input.hpp"
#include "../core/lsdtt_raster_result.hpp"
#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_chi_network_analysis_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Per-pixel record of how the ridge was found.
enum class hillslope_ridge_status : std::uint8_t
{
    not_hillslope = 0,  // nodata or a channel pixel
    no_channel    = 1,  // valid, but its D8 path reaches no channel pixel
    d8_upslope    = 2,  // ridge found on D8 paths (farthest or highest)
    crest_divide  = 3,  // aspect ray stopped at a crest on a divide
    ray_channel   = 4,  // aspect ray reached a channel: highest cell on the ray
    ray_boundary  = 5,  // aspect ray reached nodata or the grid edge: highest cell
    no_aspect     = 6   // no aspect, so no ray and no ridge
};

// Divide edges on the full grid. east[i] is the edge between cell i and the cell
// to its east; south[i] is the edge between cell i and the cell to its south.
struct hillslope_divide_edges
{
    int n_rows = 0;
    int n_cols = 0;
    std::vector<std::uint8_t> east;
    std::vector<std::uint8_t> south;
    std::int64_t n_edges = 0;
    std::int64_t n_channel_links = 0;
};

struct hillslope_relief_options
{
    std::string ridge    = "farthest";       // farthest | highest | aspect_ray
    // farthest and highest only:
    std::string grouping = "channel_pixel";  // channel_pixel | entry_pixel | pixel
    // aspect_ray only. ray_aspect is required (see compute_hillslope_ray_aspect).
    // divides is optional: null builds them from the flow and channel mask.
    const lsdtt_raster_result*    ray_aspect = nullptr;
    const hillslope_divide_edges* divides    = nullptr;
    int n_threads = 0;  // 0 = the runtime default
};

struct hillslope_relief_result
{
    // Flat indices on the full grid; -1 where undefined. Channel pixels are
    // their own channel_flat, and have no entry or ridge.
    std::vector<std::int64_t> channel_flat;
    std::vector<std::int64_t> entry_flat;
    std::vector<std::int64_t> ridge_flat;
    std::vector<std::uint8_t> ridge_status;  // hillslope_ridge_status

    // Hillslope pixels only; nodata elsewhere.
    lsdtt_raster_result hand;              // metres
    lsdtt_raster_result ridge_relief;      // metres
    lsdtt_raster_result relative_height;   // HAND / ridge_relief
    lsdtt_raster_result channel_distance;  // D8 path, pixel to channel (metres)
    // farthest, highest: D8 path, ridge to channel. aspect_ray: straight line,
    // pixel to ridge (metres).
    lsdtt_raster_result ridge_distance;

    std::int64_t n_channel      = 0;
    std::int64_t n_hillslope    = 0;  // valid non-channel pixels reaching a channel
    std::int64_t n_no_channel   = 0;  // valid non-channel pixels reaching none
    std::int64_t n_zero_relief  = 0;  // ridge_relief <= 0: relative height is nodata
    std::int64_t n_above_ridge  = 0;  // relative height > 1
    // aspect_ray only: how each hillslope pixel's ridge was found.
    std::int64_t n_crest_divide = 0;
    std::int64_t n_ray_channel  = 0;
    std::int64_t n_ray_boundary = 0;
    std::int64_t n_no_aspect    = 0;
};

// The divide edges of the channel network (see THE DIVIDE NETWORK above).
// Throws std::invalid_argument for a wrong array size or periodic boundaries,
// and std::runtime_error if the flow object has no s_vector.
[[nodiscard]] hillslope_divide_edges build_hillslope_divide_edges(
    const lsdtt_flowinfo&            flow,
    const std::vector<std::uint8_t>& elevation_mask,
    const std::vector<std::uint8_t>& channel_mask);

// Downslope aspect in compass degrees relative to grid north, for aspect rays.
//   method polyfit: the polyfit surface, length scale lambda_metres (half
//                   wavelength; <= 0 means three cells).
//   method dinf:    the D-infinity facet direction; noisier. boundary_conditions
//                   is used for the routing, and lambda_metres is ignored.
// Flat or sink cells are nodata. Throws std::invalid_argument for another method.
[[nodiscard]] lsdtt_raster_result compute_hillslope_ray_aspect(
    const lsdtt_raster_cpu& surface,
    const std::string&      method,
    double                  lambda_metres,
    const std::string&      boundary_conditions = "oooo",
    int                     n_threads = 0);

// elevation, elevation_mask and channel_mask have one value per grid cell
// (n_rows * n_cols); a non-zero channel_mask value marks a channel pixel.
//
// Throws std::invalid_argument for an unknown ridge or grouping, a wrong array
// size, a grouping given with aspect_ray, a missing or wrongly sized ray_aspect,
// or periodic boundaries with aspect_ray; std::runtime_error if the flow object
// has no s_vector or aspect_ray is used on a grid not projected in metres.
[[nodiscard]] hillslope_relief_result compute_hillslope_relief(
    const lsdtt_flowinfo&            flow,
    const std::vector<float>&        elevation,
    const std::vector<std::uint8_t>& elevation_mask,
    const std::vector<std::uint8_t>& channel_mask,
    const hillslope_relief_options&  options = {});

// =============================================================================
// CHANNEL PROPERTIES
//
// nodes / segments come from export_chi_network_fit_nodes / _segments. Within
// each (basin_index, branch_id) group a node takes the segment whose chi
// interval contains it; a node just outside every interval (endpoint rounding)
// takes the nearest segment in that branch. flow_distance and drainage_area are
// optional per-cell rasters; nullptr leaves those columns undefined, not zero.
// =============================================================================
struct lsdtt_channel_properties
{
    double       k_sn              = std::numeric_limits<double>::quiet_NaN();
    double       chi               = std::numeric_limits<double>::quiet_NaN();
    double       channel_elevation = std::numeric_limits<double>::quiet_NaN();
    double       flow_distance_m   = std::numeric_limits<double>::quiet_NaN();
    double       drainage_area_m2  = std::numeric_limits<double>::quiet_NaN();
    std::int32_t basin_index       = -1;
    std::int32_t branch_id         = -1;
    std::int32_t segment_id        = -1;
};

// A per-DEM-cell lookup of channel properties, empty where the cell is not a
// fitted channel node.
using lsdtt_channel_property_map =
    std::unordered_map<std::int64_t, lsdtt_channel_properties>;

// =============================================================================
// build_channel_property_map
//
// nodes / segments come from export_chi_network_fit_nodes / _segments.
// flow_distance and drainage_area are optional per-cell rasters; pass nullptr
// to leave those columns undefined rather than zero.
// =============================================================================
[[nodiscard]] lsdtt_channel_property_map build_channel_property_map(
    const lsdtt_chi_network_fit_node_table&    nodes,
    const lsdtt_chi_network_fit_segment_table& segments,
    const std::vector<float>*                  flow_distance,
    const std::vector<float>*                  drainage_area);

// =============================================================================
// TAGS
//
// tag(p) = channel_values[channel_flat[p]], through the D8 channel link, so a tag
// does not depend on the ridge method. Channel pixels carry their own value.
// Nothing is interpolated: a pixel whose channel pixel has no value is nodata.
// =============================================================================
struct hillslope_tag_result
{
    lsdtt_raster_result tag;
    std::int64_t n_channel_tagged   = 0;
    std::int64_t n_hillslope_tagged = 0;
    // Hillslope pixels left as nodata, by reason:
    std::int64_t n_no_value = 0;  // its channel pixel has no value
    std::int64_t n_excluded = 0;  // its channel pixel is in an excluded basin
    std::int64_t n_no_channel = 0;  // its D8 path reaches no channel
};

// All vectors have one value per grid cell. valid_mask marks valid DEM cells.
// channel_included is optional: where it is 0 the channel pixel's basin is
// excluded, so its pixels are nodata and counted in n_excluded.
// Throws std::invalid_argument for a wrong size or a channel_flat off the grid.
[[nodiscard]] hillslope_tag_result tag_hillslopes_from_channel(
    const lsdtt_raster_meta&         meta,
    const std::vector<std::int64_t>& channel_flat,
    const std::vector<std::uint8_t>& valid_mask,
    const std::vector<float>&        channel_values,
    const std::vector<std::uint8_t>& channel_value_mask,
    const std::vector<std::uint8_t>* channel_included = nullptr);

// =============================================================================
// BASINS VALID FOR CHI AND k_sn TAGS
//
// Chi needs a true baselevel, and nodata geometry cannot separate a coast from a
// clipped DEM, so the caller declares the valid basins, by exactly one of:
//
//   outlet_flats          these basins
//   max_outlet_elevation_m  every baselevel basin with its outlet at or below
//                         this height
//
// A basin with any pixel on the grid edge is always excluded. A basin below
// minimum_contributing_pixels (counted as the chi extractor does) is excluded.
// Of nested outlets that remain, only the outermost is kept, so that each
// channel pixel has one chi.
// =============================================================================
struct hillslope_tag_basin_options
{
    std::vector<std::int64_t> outlet_flats;
    double max_outlet_elevation_m = std::numeric_limits<double>::quiet_NaN();
    int    minimum_contributing_pixels = 1;
};

struct hillslope_tag_basins
{
    std::vector<std::int64_t> outlets;              // ascending flat index
    std::vector<std::uint8_t> in_selected_basin;    // one per grid cell
    std::int64_t n_candidates          = 0;
    std::int64_t n_excluded_elevation  = 0;
    std::int64_t n_excluded_grid_edge  = 0;
    std::int64_t n_excluded_too_small  = 0;
    std::int64_t n_excluded_nested     = 0;
};

// Throws std::invalid_argument unless exactly one of the two settings is given,
// for an outlet that is not a valid cell, or for a wrong array size;
// std::runtime_error if the flow object has no s_vector.
[[nodiscard]] hillslope_tag_basins select_hillslope_tag_basins(
    const lsdtt_flowinfo&              flow,
    const std::vector<float>&          elevation,
    const hillslope_tag_basin_options& options);

struct hillslope_chi_tag_options
{
    lsdtt_chi_options chi;
    float critical_divergence_m = 5.0f;
    int   minimum_contributing_pixels = 1;  // the channel source threshold
};

struct hillslope_chi_tag_result
{
    hillslope_tag_result chi;
    hillslope_tag_result k_sn;
    std::int64_t n_fitted_nodes = 0;
    std::int64_t n_segments     = 0;
};

// Fits chi on the selected basins and tags hillslopes with chi and k_sn. The
// fitted network comes from the source threshold, so channel pixels of an
// ingested network that it does not contain have no value and are counted in
// n_no_value. Selecting no basins gives all-nodata tags, not an error.
[[nodiscard]] hillslope_chi_tag_result compute_hillslope_chi_tags(
    const lsdtt_flowinfo&            flow,
    const lsdtt_raster_cpu&          routing_surface,
    const std::vector<std::int64_t>& channel_flat,
    const hillslope_tag_basins&      basins,
    const hillslope_chi_tag_options& options);

// =============================================================================
// IMPORTED CHANNEL ATTRIBUTES
//
// Joins a table (read_channel_attribute_table) to channel pixels:
//
//   exact pixel  must be a channel pixel, otherwise dropped and counted.
//   coordinates  the nearest channel pixel centre within max_snap_distance_m
//                (inclusive; ties go to the lower flat index). Further samples
//                are dropped.
//
// Several records on one pixel give their median, ignoring NaN. A line counts
// once per pixel however many of its samples land there.
//
// GAPS. A channel pixel with no record takes the nearest record on its own
// junction-to-junction link, by D8 channel distance, up to
// max_gap_fill_distance_m (0 turns filling off). Equidistant records give their median. A junction belongs to the link below
// it, so filling never crosses a junction. Filled values are never used as
// sources. Longer gaps stay nodata.
// =============================================================================
struct channel_attribute_column
{
    std::string name;
    std::vector<float>        values;  // one per grid cell
    std::vector<std::uint8_t> mask;    // 1 on channel pixels with a value
    std::int64_t n_direct           = 0;  // channel pixels with a record
    std::int64_t n_multiple_records = 0;  // of those, with more than one
    std::int64_t n_filled           = 0;  // channel pixels gap-filled
    std::int64_t n_missing          = 0;  // channel pixels left as nodata
};

struct channel_attribute_join_result
{
    std::vector<channel_attribute_column> columns;
    std::int64_t n_channel_pixels  = 0;
    std::int64_t n_records         = 0;
    std::int64_t n_records_matched = 0;  // at least one sample on a channel pixel
    std::int64_t n_records_dropped = 0;
    // Matched records with a sample moved more than one cell size.
    std::int64_t n_records_snapped_beyond_one_cell = 0;
    // Samples with an exact pixel that is not a channel pixel.
    std::int64_t n_exact_samples_off_channel = 0;
};

// Throws std::invalid_argument for a wrong array size, an inconsistent table, a
// sample off the grid, or a negative or non-finite distance limit;
// std::runtime_error if the flow object has no s_vector.
[[nodiscard]] channel_attribute_join_result join_channel_attributes(
    const lsdtt_flowinfo&            flow,
    const std::vector<std::uint8_t>& elevation_mask,
    const std::vector<std::uint8_t>& channel_mask,
    const channel_attribute_table&   table,
    double                           max_snap_distance_m,
    double                           max_gap_fill_distance_m);

// =============================================================================
// ONE ROW PER CHANNEL PIXEL
//
// How many hillslope pixels drain to each channel pixel, and the median of a
// hillslope raster (for example relative height) over those with a value.
// Rows are in ascending flat order. median is NaN where no pixel has a value.
// =============================================================================
struct hillslope_channel_pixel_summary
{
    std::vector<std::int64_t> flat;
    std::vector<std::int64_t> n_hillslope;
    std::vector<std::int64_t> n_with_value;
    std::vector<double>       median;
};

// Throws std::invalid_argument for a wrong array size or a link off the grid.
[[nodiscard]] hillslope_channel_pixel_summary summarise_hillslopes_by_channel_pixel(
    const std::vector<std::int64_t>& channel_flat,
    const lsdtt_raster_result&       values);

} // namespace lsdtt
