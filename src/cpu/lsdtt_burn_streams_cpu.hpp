// lsdtt_burn_streams_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_burn_streams_cpu
//
// Burns a user-supplied polyline vector into a DEM by sampling the elevation
// profile of pixels that fall beneath the vector channel and enforcing a
// monotone-decreasing (downhill) elevation constraint along the profile.
//
// Two burn modes are provided (via stream_burn_mode):
//
//   burn_in_place
//     For each polyline feature in the input vector:
//       1. For every valid DEM pixel within half a pixel diagonal of the
//          polyline, record its along-axis distance and DEM elevation.
//       2. Determine which end of the polyline is the outlet (lower elevation).
//       3. Apply a running-minimum scan from upstream to downstream so that
//          the burned elevation profile is always non-increasing downstream.
//     Output: one burned_stream_profile per feature, containing per-node
//     flow_distance (from outlet), dem_elevation, burned_elevation, and x/y.
//
//   snap_to_network
//     Not yet implemented.  Throws std::runtime_error if called.
//
// Input vector CRS must match the DEM CRS; no on-the-fly reprojection is
// performed when reading polyline vertices.
//
// Output:
//   write_burned_stream_profiles_csv  — CSV table, one row per node
//   write_burned_stream_profiles_points — vector point file (gpkg / flatgeobuf)
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// vector_network_endpoints
//
// Result of find_vector_network_source_tips.
//
// source_tips  — upstream headwater channel endpoints (degree-1 tips that are
//                not the lowest-elevation tip of the network).
// outlet_tip   — the single lowest-elevation degree-1 tip (basin outlet).
//                has_outlet is false if no valid DEM elevation could be
//                sampled at any tip (e.g. all tips fall outside the raster).
// ---------------------------------------------------------------------------
struct vector_network_endpoints
{
    std::vector<lsdtt_xy> source_tips;
    lsdtt_xy              outlet_tip;
    bool                  has_outlet = false;
};

// ---------------------------------------------------------------------------
// find_vector_network_source_tips
//
// Identifies upstream source tips from a polyline vector channel network by
// finding degree-1 endpoints — polyline vertices that are not shared by more
// than one feature at that position (within half a DEM cell of tolerance).
//
// Of all degree-1 tips, the one with the lowest DEM elevation is classified
// as the basin outlet. All others are source tips.
//
// vector_path must point to a polyline (linestring / multilinestring) layer
// in the same CRS as the DEM.
//
// Outlet identification uses the segment gradient criterion only.
// All inward vertices of the connected segment are walked and the minimum
// valid DEM elevation among them is recorded.  That minimum is compared to
// the tip elevation:
//   min_inward > tip_elev  →  all visible inward terrain is higher  →  tip
//                             is at the downstream (outlet) end.
//   min_inward < tip_elev  →  the segment descends away from the tip  →  tip
//                             is at the upstream (source) end.
// Using the minimum (rather than just the first valid vertex) prevents the
// spurious positive gradient that can arise when the first vertex after a
// source tip sits in a slightly higher DEM cell due to cell-alignment noise
// on flat terrain.  It also naturally handles segments that cross the DEM
// boundary with several nodata vertices before reaching valid terrain.
// Among multiple outlet candidates the one with the greatest upslope_count
// (s_vector_size from the flow object) is chosen; steepest gradient breaks
// ties.  If no tip satisfies the downstream-end condition has_outlet is false
// and all degree-1 tips are returned as sources — the real outlet is either
// at a junction (degree > 1) or outside the DEM.  Fallbacks to minimum
// elevation or maximum upslope count are intentionally absent: both force an
// outlet selection even when none exists, silently removing a source tip.
//
// flow — optional D8 flow object built on the same DEM.  Provides upslope
//         pixel counts used as a tiebreaker in criterion 1 and as criterion 2.
// ---------------------------------------------------------------------------
[[nodiscard]] vector_network_endpoints find_vector_network_source_tips(
    const std::string&      vector_path,
    const lsdtt_raster_cpu& dem,
    const lsdtt_flow_base*  flow = nullptr
);

// ---------------------------------------------------------------------------
// stream_burn_mode
// ---------------------------------------------------------------------------
enum class stream_burn_mode
{
    burn_in_place,
    snap_to_network
};

// Resolve from string ("burn_in_place" / "snap_to_network").
stream_burn_mode resolve_stream_burn_mode(const std::string& s);

// ---------------------------------------------------------------------------
// burned_stream_node
//
// One node in a burned channel profile.
// ---------------------------------------------------------------------------
struct burned_stream_node
{
    double flow_distance   = 0.0;   // distance from outlet along polyline (m)
    double dem_elevation   = 0.0;   // original DEM elevation (m)
    double burned_elevation = 0.0;  // enforced downhill-flow elevation (m)
    double x               = 0.0;   // projected easting (m)
    double y               = 0.0;   // projected northing (m)
};

// ---------------------------------------------------------------------------
// burned_stream_profile
//
// Profile for one polyline feature.  Nodes are ordered by flow_distance
// ascending (outlet first, headwater last).
// ---------------------------------------------------------------------------
struct burned_stream_profile
{
    std::int32_t                    feature_id = 0;
    std::vector<burned_stream_node> nodes;
};

// ---------------------------------------------------------------------------
// burn_streams_in_place
//
// For each polyline feature in `vector_path`, extracts the DEM elevation
// profile from pixels that lie beneath the vector and computes a monotone
// burned profile that always flows downhill.
//
// `vector_path` must contain polyline (linestring / multilinestring) features.
//  Other feature types are skipped with a warning.
//
// The vector CRS must match the DEM CRS; vertices are read as-is.
//
// Returns one profile per polyline feature (empty profiles omitted).
// ---------------------------------------------------------------------------
[[nodiscard]] std::vector<burned_stream_profile> burn_streams_in_place(
    const lsdtt_raster_cpu& dem,
    const std::string&      vector_path,
    int                     n_threads = 0
);

// ---------------------------------------------------------------------------
// vector_dem_junction_pair
//
// One matched pair of junctions: a first-layer junction from the supplied
// vector network and the spatially nearest junction in the DEM-extracted
// network.
// ---------------------------------------------------------------------------
struct vector_dem_junction_pair
{
    int      pair_id      = 0;
    lsdtt_xy vector_xy;          // junction location in the supplied vector
    lsdtt_xy dem_xy;             // nearest junction in the DEM network
    double   distance_m   = 0.0; // Euclidean distance between them (m)
    int      n_sources    = 0;   // number of source tips that merge here
};

// ---------------------------------------------------------------------------
// find_first_layer_vector_junctions
//
// Reads the supplied polyline vector network and returns the locations of
// "first-layer" junctions — nodes where exactly the channels coming directly
// from upstream source tips (degree-1 endpoints) converge.
//
// Orientation at each junction is determined by DEM elevation: segments whose
// far endpoint sits higher than the junction are upstream; the segment whose
// far endpoint is lower is downstream.  A junction is first-layer when every
// upstream-connected segment ends at a degree-1 source tip.
//
// endpoints is used to exclude the outlet tip from the source-tip set when
// has_outlet is true.
// ---------------------------------------------------------------------------
[[nodiscard]] std::vector<vector_dem_junction_pair>
find_first_layer_junction_pairs(
    const std::string&                     vector_path,
    const lsdtt_raster_cpu&                dem,
    const vector_network_endpoints&        endpoints,
    const lsdtt_junction_network_base&     dem_network
);

// ---------------------------------------------------------------------------
// write_junction_comparison_points
//
// Writes a single point file containing both members of every junction pair.
// Each row represents one point (vector or DEM junction).
//
// Attributes per feature:
//   pair_id        (int64)  — shared by the two points in a pair
//   is_dem_junction (int64) — 0 = vector network junction, 1 = DEM junction
//   distance_m     (double) — distance between the paired points
//   n_sources      (int64)  — source tips converging at this junction
// ---------------------------------------------------------------------------
void write_junction_comparison_points(
    const std::vector<vector_dem_junction_pair>& pairs,
    const lsdtt_raster_meta&                      meta,
    const std::string&                            filename,
    point_output_format                           fmt,
    bool                                          spatial_index = false
);

// ---------------------------------------------------------------------------
// Write helpers
// ---------------------------------------------------------------------------

// Write all profiles to a CSV file.
//
// Schema:
//   feature_id, flow_distance, x, y, dem_elevation, burned_elevation
//
// Rows are ordered per feature, then by flow_distance ascending.
void write_burned_stream_profiles_csv(
    const std::vector<burned_stream_profile>& profiles,
    const std::string&                        filename
);

// Write all profiles as vector points (gpkg or flatgeobuf).
//
// Attributes per feature:
//   feature_id (int64), flow_distance (double),
//   dem_elevation (double), burned_elevation (double)
//
// fmt must not be csv; use write_burned_stream_profiles_csv for that.
void write_burned_stream_profiles_points(
    const std::vector<burned_stream_profile>& profiles,
    const lsdtt_raster_meta&                  meta,
    const std::string&                        filename,
    point_output_format                       fmt,
    bool                                      spatial_index = false
);

} // namespace lsdtt
