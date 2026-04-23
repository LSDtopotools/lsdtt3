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
#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

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
