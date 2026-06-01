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
};

struct lsdtt_side_tributary_basin_options
{
    std::string input_points_path;
    std::string input_points_layer;

    int  minimum_contributing_pixels = 1000;
    bool reject_incomplete_basins    = true;

    bool spatial_index = true;
};

[[nodiscard]] lsdtt_side_tributary_basin_result extract_side_tributary_basins(
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

} // namespace lsdtt
