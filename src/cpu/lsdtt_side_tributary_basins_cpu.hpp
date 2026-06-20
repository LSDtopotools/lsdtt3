// lsdtt_side_tributary_basins_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Side tributary basin extraction from downstream traces seeded by vector points.
//
// Workflow:
//   1. Read upstream seed points from FlatGeobuf / GeoPackage / any OGR vector.
//   2. Transform point coordinates to the raster CRS and map them to DEM pixels.
//   3. Trace D8 downstream from each point, building a unique traced-channel mask.
//   4. Find non-channel pixels whose receiver lies on the traced channel.
//   5. Keep those side outlets whose D8 upslope pixel count exceeds a threshold.
//   6. Reject incomplete side basins if any upslope pixel is influenced by nodata
//      or an open/out-of-domain edge.
//   7. Write side outlet points and basin outlines.
//
// Notes:
//   - The flow object must have been built with build_s_vector = true.
//   - Basin outlines are written by constructing a labelled side-basin raster and
//     vectorising it using the binary-mask polygon tracer in
//     lsdtt_raster_binary_mask_vectorise_cpu.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt
{

struct lsdtt_side_tributary_basin_record
{
    std::int64_t outlet_flat = -1;
    int          outlet_row  = -1;
    int          outlet_col  = -1;

    std::int64_t receiver_channel_flat = -1;
    int          receiver_channel_row  = -1;
    int          receiver_channel_col  = -1;

    std::uint64_t n_pixels      = 0;
    double        drainage_area = 0.0;

    bool complete = false;
};

struct lsdtt_side_tributary_basin_result
{
    std::vector<std::uint8_t> traced_channel_mask;
    std::vector<lsdtt_side_tributary_basin_record> accepted_basins;

    // Number of valid seed (source) points that were traced. Zero means no
    // usable sources were found in the input vector.
    std::size_t n_seed_points = 0;

    // True if a stopping point was supplied and used to truncate the traces.
    bool stopping_point_applied = false;
};

struct lsdtt_side_tributary_basin_options
{
    std::string input_points_path;
    std::string input_points_layer;

    int  minimum_contributing_pixels = 1000;
    bool reject_incomplete_basins    = true;

    // If true, the final point in the input vector is treated as a stopping
    // point rather than a source. Each downstream trace is truncated when it
    // first enters a square window of half-width stopping_point_window_pixels
    // around the stopping point.
    bool include_stopping_point        = false;
    int  stopping_point_window_pixels  = 5;

    bool spatial_index = true;
};

[[nodiscard]] lsdtt_side_tributary_basin_result extract_side_tributary_basins(
    const lsdtt_flow_base& flow,
    const lsdtt_side_tributary_basin_options& opts
);

// Build a SINGLE main channel from an upstream point and extract the side
// tributary basins draining into it. The input file holds one or two points:
//   - one point  : trace D8 from it to the DEM edge (the single channel).
//   - two points : identify the upstream point by D8 reachability (the point
//     whose flow path the other joins), trace from it, and truncate at the
//     "common downstream point" where the other point's flow meets the channel.
//     If neither point's flow reaches the other (different sub-basins), a
//     warning is emitted and the longer of the two D8 paths is used.
// Side basins are then found exactly as in extract_side_tributary_basins
// (outlet = pixel just before joining; complete + minimum-size filtered).
[[nodiscard]] lsdtt_side_tributary_basin_result extract_single_channel_side_basins(
    const lsdtt_flow_base& flow,
    const lsdtt_side_tributary_basin_options& opts
);

void write_side_tributary_outlets(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_side_tributary_basin_record>& basins,
    const std::string& out_path,
    point_output_format format,
    const std::string& layer_name = "side_tributary_outlets",
    bool spatial_index = true
);

void write_side_tributary_basin_outlines(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_side_tributary_basin_record>& basins,
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name = "side_tributary_basin_outlines",
    bool spatial_index = true
);

// Write the traced main channel (the D8 flow paths from the source points down
// to the outlet or stopping point) as line features. Each traced pixel emits one
// segment to its receiver, so the result is the full set of channel links from
// which side tributaries were extracted.
void write_side_tributary_main_channel(
    const lsdtt_flow_base& flow,
    const std::vector<std::uint8_t>& traced_channel_mask,
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name = "side_tributary_main_channel",
    bool spatial_index = true
);

} // namespace lsdtt
