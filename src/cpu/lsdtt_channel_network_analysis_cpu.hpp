// lsdtt_channel_network_analysis_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#pragma once

// -----------------------------------------------------------------------------
// ANGULAR CONVENTION (project-wide)
//
// Every bearing, aspect and orientation lsdtt3 emits uses the COMPASS
// convention: 0 degrees = NORTH, increasing CLOCKWISE (90 = east, 180 = south,
// 270 = west).  This matches GDAL, ArcGIS and the geomorphology literature.
//
// It is NOT the mathematical convention (0 = +x axis, increasing
// anticlockwise), which in map view puts 0 at EAST.  To convert:
//     compass = fmod(90 - maths + 360, 360)
// Note that `atan2(y, x)` is the maths form; a compass angle is
// `atan2(east, north)` -- arguments in that order.
//
// ASPECT additionally means the DOWN-SLOPE azimuth: the direction the slope
// FACES.  With p = dz/dx (east-positive) and q = dz/dy (north-positive) that is
// `atan2(-p, -q)`.  `atan2(p, q)` is the up-slope azimuth and is 180 degrees
// wrong.
//
// AXIAL quantities (line orientations, which have no head or tail: MADk2 DMC,
// basin principal axis, skeleton dominant orientation, feature azimuths) use
// the same zero and sense but fold to [0, 180).
// -----------------------------------------------------------------------------


#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <limits>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_true_north.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt
{

enum class lsdtt_link_sampling_mode
{
    full_link,
    fixed_length,
    fixed_n_pixels
};

enum class lsdtt_link_sample_anchor
{
    upstream_end,
    downstream_end
};

struct lsdtt_link_sampling_options
{
    lsdtt_link_sampling_mode mode = lsdtt_link_sampling_mode::full_link;

    double fixed_length_metres = 0.0;
    int    fixed_n_pixels = 0;

    // Minimum number of sampled nodes required for later analysis.
    // Keep this small here; bearing code can impose stricter requirements later.
    int min_nodes = 2;
};

struct lsdtt_traced_channel_link
{
    std::int32_t upstream_junction   = k_nodata_junction;
    std::int32_t downstream_junction = k_nodata_junction;

    std::int64_t upstream_flat   = -1;
    std::int64_t downstream_flat = -1;

    std::vector<std::int64_t> flats;
    std::vector<lsdtt_xy>     xy;
    std::vector<double>       cumulative_distance_from_upstream_metres;

    double total_length_metres = 0.0;

    bool empty() const noexcept
    {
        return flats.empty();
    }

    std::size_t size() const noexcept
    {
        return flats.size();
    }
};

struct lsdtt_sampled_channel_link
{
    std::int32_t upstream_junction   = k_nodata_junction;
    std::int32_t downstream_junction = k_nodata_junction;

    std::vector<std::int64_t> flats;
    std::vector<lsdtt_xy>     xy;

    // Distance measured along the sampled polyline, in sample order.
    std::vector<double> cumulative_distance_metres;

    double sampled_length_metres = 0.0;

    bool empty() const noexcept
    {
        return flats.empty();
    }

    std::size_t size() const noexcept
    {
        return flats.size();
    }
};

// -----------------------------------------------------------------------------
// Trace the compact channel link from one junction to its immediate downstream
// receiver junction. Includes both endpoint junction pixels.
// Throws if upstream_junction is invalid, is an outlet junction, or the traced
// path leaves the channel network or fails to reach the declared receiver.
// -----------------------------------------------------------------------------
lsdtt_traced_channel_link trace_channel_link_between_junctions(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    std::int32_t                       upstream_junction
);

// -----------------------------------------------------------------------------
// Sample a traced link according to the chosen mode.
//
// full_link:
//   returns all nodes
//
// fixed_n_pixels:
//   returns a contiguous sub-chain of at most fixed_n_pixels nodes,
//   anchored at either the upstream or downstream end
//
// fixed_length:
//   returns a contiguous sub-chain whose path length is <= fixed_length_metres,
//   anchored at either the upstream or downstream end
//
// The returned nodes preserve channel order.
// -----------------------------------------------------------------------------
lsdtt_sampled_channel_link sample_traced_channel_link(
    const lsdtt_traced_channel_link&    traced_link,
    const lsdtt_link_sampling_options&  options,
    lsdtt_link_sample_anchor            anchor
);

// Convenience helper for reversing a sampled link in-place.
void reverse_sampled_channel_link(
    lsdtt_sampled_channel_link& sampled_link
);


struct lsdtt_channel_bearing_fit
{
    bool valid = false;

    std::size_t n_nodes_used = 0;
    double sampled_length_metres = std::numeric_limits<double>::quiet_NaN();

    double direction_x = std::numeric_limits<double>::quiet_NaN();
    double direction_y = std::numeric_limits<double>::quiet_NaN();

    double bearing_rad = std::numeric_limits<double>::quiet_NaN();
    double bearing_deg = std::numeric_limits<double>::quiet_NaN();

    double fit_linearity = std::numeric_limits<double>::quiet_NaN();

    bool empty() const noexcept
    {
        return !valid;
    }
};

struct lsdtt_channel_link_bearing
{
    std::int32_t upstream_junction   = k_nodata_junction;
    std::int32_t downstream_junction = k_nodata_junction;

    std::int64_t upstream_flat   = -1;
    std::int64_t downstream_flat = -1;

    std::int32_t stream_order = -1;

    double drainage_area_m2 = std::numeric_limits<double>::quiet_NaN();

    std::size_t n_nodes_total = 0;
    std::size_t n_nodes_used  = 0;

    double total_length_metres   = std::numeric_limits<double>::quiet_NaN();
    double sampled_length_metres = std::numeric_limits<double>::quiet_NaN();

    double direction_x = std::numeric_limits<double>::quiet_NaN();
    double direction_y = std::numeric_limits<double>::quiet_NaN();

    double bearing_rad = std::numeric_limits<double>::quiet_NaN();
    double bearing_deg = std::numeric_limits<double>::quiet_NaN();

    double fit_linearity = std::numeric_limits<double>::quiet_NaN();

    // Representative point used for visual bearing-vector output.
    // This is taken from the centre sampled node of the fitted link segment.
    std::int64_t representative_flat = -1;

    double representative_x = std::numeric_limits<double>::quiet_NaN();
    double representative_y = std::numeric_limits<double>::quiet_NaN();

    // The pixels the bearing was fitted on (the sampled sub-chain, in channel
    // order; size == n_nodes_used). Not written to any output; kept so that
    // per-record attributes (e.g. the lithologic unit) can be taken over
    // exactly the nodes that set the bearing.
    std::vector<std::int64_t> sampled_flats;

    bool empty() const noexcept
    {
        return upstream_junction == k_nodata_junction;
    }
};

// -----------------------------------------------------------------------------
// Fit a directed channel bearing to an already-sampled link.
//
// The sampled link must be ordered in channel order. The fitted orthogonal
// regression axis is oriented to match the chord from the first sampled point
// to the last sampled point.
// -----------------------------------------------------------------------------
lsdtt_channel_bearing_fit fit_channel_bearing_to_sampled_link(
    const lsdtt_sampled_channel_link& sampled_link
);

// -----------------------------------------------------------------------------
// Calculate one bearing record per downstream compact link.
//
// Input junctions are interpreted as the upstream ends of the links to analyse.
// Outlet junctions are skipped because they have no downstream receiver link.
//
// If upstream_junctions_to_analyse is empty, all non-outlet junctions are used.
// -----------------------------------------------------------------------------
std::vector<lsdtt_channel_link_bearing> calculate_channel_link_bearings(
    const lsdtt_flow_base&                    flow,
    const lsdtt_junction_network_base&        network,
    const std::vector<std::int32_t>&          upstream_junctions_to_analyse = {},
    const lsdtt_link_sampling_options&        options = lsdtt_link_sampling_options{},
    lsdtt_link_sample_anchor                  anchor = lsdtt_link_sample_anchor::upstream_end
);



struct lsdtt_junction_angle_record
{
    bool valid = false;

    std::int32_t junction  = k_nodata_junction;
    std::int32_t donor_1   = k_nodata_junction;
    std::int32_t donor_2   = k_nodata_junction;
    std::int32_t receiver  = k_nodata_junction;

    std::int32_t junction_order = -1;
    std::int32_t donor_1_order  = -1;
    std::int32_t donor_2_order  = -1;
    std::int32_t receiver_order = -1;

    // Drainage areas:
    // donor_*_area_m2 are taken at the penultimate pixel immediately upstream
    // of the analysed junction along each donor link.
    // junction_area_m2 is the area at the junction pixel itself.
    double donor_1_area_m2 = std::numeric_limits<double>::quiet_NaN();
    double donor_2_area_m2 = std::numeric_limits<double>::quiet_NaN();
    double junction_area_m2 = std::numeric_limits<double>::quiet_NaN();

    // Bearings
    // donor_*_bearing_* are inward-pointing, i.e. toward the junction
    // receiver_bearing_out_* points downstream away from the junction
    // receiver_bearing_in_* points back toward the junction
    double donor_1_bearing_rad = std::numeric_limits<double>::quiet_NaN();
    double donor_2_bearing_rad = std::numeric_limits<double>::quiet_NaN();
    double receiver_bearing_out_rad = std::numeric_limits<double>::quiet_NaN();
    double receiver_bearing_in_rad  = std::numeric_limits<double>::quiet_NaN();

    double donor_1_bearing_deg = std::numeric_limits<double>::quiet_NaN();
    double donor_2_bearing_deg = std::numeric_limits<double>::quiet_NaN();
    double receiver_bearing_out_deg = std::numeric_limits<double>::quiet_NaN();
    double receiver_bearing_in_deg  = std::numeric_limits<double>::quiet_NaN();

    // Angles
    double donor_opening_angle_rad = std::numeric_limits<double>::quiet_NaN();
    double donor_1_receiver_angle_rad = std::numeric_limits<double>::quiet_NaN();
    double donor_2_receiver_angle_rad = std::numeric_limits<double>::quiet_NaN();

    double donor_opening_angle_deg = std::numeric_limits<double>::quiet_NaN();
    double donor_1_receiver_angle_deg = std::numeric_limits<double>::quiet_NaN();
    double donor_2_receiver_angle_deg = std::numeric_limits<double>::quiet_NaN();

    // Bending angle (NOT Strong & Mudd's beta): the direction of the midpoint
    // of the donor opening sector, measured from the inward receiver direction,
    // bending_angle_signed = wrap(180 + (beta_1 - beta_2) / 2). A symmetric
    // junction (beta_1 = beta_2) therefore reads +/-180. abs is |signed|.
    double bending_angle_signed_rad = std::numeric_limits<double>::quiet_NaN();
    double bending_angle_abs_rad    = std::numeric_limits<double>::quiet_NaN();

    double bending_angle_signed_deg = std::numeric_limits<double>::quiet_NaN();
    double bending_angle_abs_deg    = std::numeric_limits<double>::quiet_NaN();

    // Link-fit diagnostics
    double donor_1_linearity = std::numeric_limits<double>::quiet_NaN();
    double donor_2_linearity = std::numeric_limits<double>::quiet_NaN();
    double receiver_linearity = std::numeric_limits<double>::quiet_NaN();

    std::size_t donor_1_nodes_total = 0;
    std::size_t donor_2_nodes_total = 0;
    std::size_t receiver_nodes_total = 0;

    std::size_t donor_1_nodes_used = 0;
    std::size_t donor_2_nodes_used = 0;
    std::size_t receiver_nodes_used = 0;

    double donor_1_total_length_metres = std::numeric_limits<double>::quiet_NaN();
    double donor_2_total_length_metres = std::numeric_limits<double>::quiet_NaN();
    double receiver_total_length_metres = std::numeric_limits<double>::quiet_NaN();

    double donor_1_sampled_length_metres = std::numeric_limits<double>::quiet_NaN();
    double donor_2_sampled_length_metres = std::numeric_limits<double>::quiet_NaN();
    double receiver_sampled_length_metres = std::numeric_limits<double>::quiet_NaN();

    // The pixels each arm's bearing was fitted on, in channel order: donors
    // sampled at their downstream end (ending at the junction), the receiver at
    // its upstream end (starting at the junction). size == *_nodes_used. Not
    // written to any output.
    std::vector<std::int64_t> donor_1_sampled_flats;
    std::vector<std::int64_t> donor_2_sampled_flats;
    std::vector<std::int64_t> receiver_sampled_flats;

    bool empty() const noexcept
    {
        return !valid;
    }
};

// -----------------------------------------------------------------------------
// Calculate junction angles for the selected junctions.
//
// Input junctions are analysed as confluences. Junctions are skipped if they:
//   - are invalid
//   - are outlet junctions
//   - do not have exactly two direct donor junctions
//   - fail donor/receiver tracing or bearing fitting
//
// If junctions_to_analyse is empty, all junctions in the network are considered.
// -----------------------------------------------------------------------------
// IMPORTANT TERMINOLOGY NOTE - JUNCTION-ANGLE NAMING
//
// This code uses descriptive internal names that differ from the historical
// Strong and Mudd (PNAS, 2023) naming convention.
//
// Historical convention:
//   alpha   = angle between the two donor channels
//   beta_1  = angle between donor_1 and the receiver
//   beta_2  = angle between donor_2 and the receiver
//
// In this implementation:
//   donor_opening_angle      = alpha
//   donor_1_receiver_angle   = beta_1
//   donor_2_receiver_angle   = beta_2
//
// So what is called "junction angle" in the earlier literature is here stored
// as donor_opening_angle.
//
// This code may also report:
//   bending_angle_signed
//   bending_angle_abs
//
// These are NEW derived diagnostics and are not part of the original
// alpha / beta_1 / beta_2 set, and bending_angle is NOT Strong & Mudd's beta.
// It is the direction of the midpoint of the donor opening sector, measured
// from the inward receiver direction:
//   bending_angle_signed = wrap(180 + (beta_1 - beta_2) / 2)
// so a symmetric junction (beta_1 = beta_2) reads +/-180. The columns are
// written unchanged and are kept out of the statistics.
//
// donor_1 and donor_2 are stable implementation labels only; they should not
// be interpreted as left/right tributaries unless an additional geometric
// relabelling step is performed.
// -----------------------------------------------------------------------------
std::vector<lsdtt_junction_angle_record> calculate_junction_angles(
    const lsdtt_flow_base&                    flow,
    const lsdtt_junction_network_base&        network,
    const std::vector<std::int32_t>&          junctions_to_analyse = {},
    const lsdtt_link_sampling_options&        options = lsdtt_link_sampling_options{}
);



// -----------------------------------------------------------------------------
// Write channel-link bearing records.
//
// Output format behaviour:
//   - csv         : tabular output with upstream/downstream coordinates
//   - gpkg / fgb  : one line feature per analysed channel link
//
// NOTE:
// We use point_output_format here so that csv remains available as an option.
// For non-csv formats, gpkg/fgb are converted internally to vector output.
// -----------------------------------------------------------------------------
// Optional per-record columns appended after every standard column of the
// record writers below (e.g. the lithologic unit of each record). names are the
// CSV headers, vector_names the attribute names in vector outputs (kept to 10
// characters for shapefiles), values[i] the values for records[i] (same order
// as the records, empty records included). nullptr = no extra columns, and the
// output is exactly what it was without the argument.
struct lsdtt_channel_record_extra_columns
{
    std::vector<std::string> names;
    std::vector<std::string> vector_names;
    std::vector<lsdtt_vector_field_type> types;
    std::vector<std::vector<lsdtt_vector_field_value>> values;
};

void write_channel_link_bearings(
    const std::vector<lsdtt_channel_link_bearing>& records,
    const lsdtt_flow_base&                         flow,
    const lsdtt_junction_network_base&             network,
    const std::string&                             out_path,
    point_output_format                            format,
    const std::string&                             layer_name = "channel_link_bearings",
    bool                                           spatial_index = true,
    const lsdtt_channel_record_extra_columns*      extra = nullptr
);

// -----------------------------------------------------------------------------
// Write junction-angle records.
//
// Output format behaviour:
//   - csv         : tabular output with junction latitude/longitude
//   - gpkg / fgb  : one point feature per analysed junction
// -----------------------------------------------------------------------------
void write_junction_angle_data(
    const std::vector<lsdtt_junction_angle_record>& records,
    const lsdtt_flow_base&                          flow,
    const lsdtt_junction_network_base&              network,
    const std::string&                              out_path,
    point_output_format                             format,
    const std::string&                              layer_name = "junction_angles",
    bool                                            spatial_index = true,
    const lsdtt_channel_record_extra_columns*       extra = nullptr
);


// -----------------------------------------------------------------------------
// write_channel_link_bearing_vectors
//
// Write one short visual bearing vector per analysed compact downstream link.
//
// Geometry:
//   one short line feature per analysed link
//
// Direction:
//   points downstream in the fitted link-bearing direction
//
// Origin:
//   representative sampled pixel used by the bearing fit
//
// Length:
//   min(mean_link_length / 3, link_length / 2)
//
// This is a visual QA layer rather than a full attribute export.
// -----------------------------------------------------------------------------
void write_channel_link_bearing_vectors(
    const std::vector<lsdtt_channel_link_bearing>& records,
    const lsdtt_flow_base&                         flow,
    const std::string&                             out_path,
    vector_output_format                           format,
    const std::string&                             layer_name = "channel_link_bearing_vectors",
    bool                                           spatial_index = true
);


// -----------------------------------------------------------------------------
// write_junction_angle_bearing_vectors
//
// Write three short centred orientation lines per analysed junction:
//
//   - donor_1
//   - donor_2
//   - receiver
//
// Geometry:
//   three short centred line features per junction
//
// Bearings:
//   donor_1 and donor_2 use inward-pointing donor bearings
//   receiver uses the outward/downstream receiver bearing
//
// This is a visual QA layer intended to let the user quickly inspect whether
// junction-angle bearings are being measured correctly.
// -----------------------------------------------------------------------------
void write_junction_angle_bearing_vectors(
    const std::vector<lsdtt_junction_angle_record>& records,
    const lsdtt_flow_base&                          flow,
    const lsdtt_junction_network_base&              network,
    const std::string&                              out_path,
    vector_output_format                            format,
    const std::string&                              layer_name = "junction_angle_bearing_vectors",
    bool                                            spatial_index = true,
    const lsdtt_grid_convergence*                   true_north = nullptr
);


// true_north (bearings_relative_to_true_north): the bearing_rad / bearing_deg
// attributes are converted at the junction; the arms are still drawn along
// the grid bearings, so the geometry is unchanged.

// =============================================================================
// Basin centre-of-mass bearing
//
// For a basin draining to a given outlet pixel, the centre of mass is the
// planimetric centroid of every upslope pixel. The centre-of-mass bearing is the
// geographic bearing (degrees from North, clockwise, in [0,360)) from the outlet
// to that centroid, a measure of the dominant direction in which the basin's
// area lies relative to its outlet.
// =============================================================================

// Bearing in degrees from North (clockwise) from (ox,oy) to (cx,cy) in a
// projected CRS (x = easting, y = northing). Returns NaN for a zero-length
// vector.
double bearing_from_north_degrees(double ox, double oy, double cx, double cy);

struct lsdtt_junction_com_bearing
{
    std::int32_t junction_id   = -1;
    std::int32_t order         = -1;   // Strahler order of the junction

    // Basin outlet = the pixel just before the downstream junction (or the
    // base-level pixel for the outlet junction).
    double       outlet_x      = std::numeric_limits<double>::quiet_NaN();
    double       outlet_y      = std::numeric_limits<double>::quiet_NaN();
    double       com_x         = std::numeric_limits<double>::quiet_NaN();
    double       com_y         = std::numeric_limits<double>::quiet_NaN();
    double       bearing_deg   = std::numeric_limits<double>::quiet_NaN();
    double       distance_m    = std::numeric_limits<double>::quiet_NaN();
    std::int64_t n_pixels      = 0;
};

// -----------------------------------------------------------------------------
// calculate_junction_com_bearings
//
// For every junction, defines the basin by the pixel just before the next
// downstream junction (the base-level pixel for the outlet junction), then
// computes the basin centre of mass and the outlet -> centre-of-mass bearing.
// Centroids are obtained from a single O(N) upslope accumulation of the cell
// coordinates, so the whole network is processed efficiently.
// -----------------------------------------------------------------------------
std::vector<lsdtt_junction_com_bearing> calculate_junction_com_bearings(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network
);

// Per-junction CSV (junction_id, order, outlet/com coordinates, bearing, ...).
void write_junction_com_bearings(
    const std::vector<lsdtt_junction_com_bearing>& records,
    const std::string&                             out_path
);

// =============================================================================
// bearings_relative_to_true_north (lsdtt_true_north.hpp)
//
// Copies of the records with every bearing converted from grid north to true
// north at the record's own location: a link at its representative point
// (NaN if it has none), junction arms at the junction pixel, a junction
// basin's outlet -> centre-of-mass bearing at the midpoint of the two.
// direction_x / direction_y stay grid components (they draw the vectors in
// map coordinates), and angles between arms are unchanged.
// =============================================================================
std::vector<lsdtt_channel_link_bearing> channel_link_bearings_to_true_north(
    const std::vector<lsdtt_channel_link_bearing>& records,
    const lsdtt_grid_convergence&                  true_north
);

std::vector<lsdtt_junction_angle_record> junction_angle_bearings_to_true_north(
    const std::vector<lsdtt_junction_angle_record>& records,
    const lsdtt_flow_base&                          flow,
    const lsdtt_junction_network_base&              network,
    const lsdtt_grid_convergence&                   true_north
);

std::vector<lsdtt_junction_com_bearing> junction_com_bearings_to_true_north(
    const std::vector<lsdtt_junction_com_bearing>& records,
    const lsdtt_grid_convergence&                  true_north
);


} // namespace lsdtt
