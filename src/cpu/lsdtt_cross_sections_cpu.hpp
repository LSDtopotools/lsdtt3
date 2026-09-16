// lsdtt_cross_sections_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_cross_sections_cpu
//
// CPU cross-section extraction on a channel network.
//
// For each channel pixel in a junction network:
//
//   1. Build a local along-channel window in metres
//   2. Fit an orthogonal linear regression to the channel-pixel centres
//   3. Use the fitted line as the channel tangent direction
//   4. Construct the orthogonal cross section through the channel pixel centre
//   5. Sample the DEM along that cross section
//   6. Optionally sample a vertical_excess_curvature raster at the same points
//
// Bearing estimation uses a user-specified along-channel window length in
// metres, not just the immediate upstream/downstream pixel.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct lsdtt_cross_section_sample
{
    std::int32_t sample_index = -1;

    // Signed offset from the section centre in metres.
    // Negative and positive follow the chosen normal direction.
    double offset_metres = 0.0;

    // Sample coordinate in raster CRS
    double x = 0.0;
    double y = 0.0;

    // DEM elevation sample
    float elevation = 0.0f;
    std::uint8_t valid = 0;

    // Optional sampled vertical excess curvature
    float vertical_excess_curvature = 0.0f;
    std::uint8_t vertical_excess_curvature_valid = 0;
};

struct lsdtt_cross_section
{
    std::int64_t cross_section_id = -1;

    std::int64_t centre_flat = -1;
    int          centre_row  = -1;
    int          centre_col  = -1;

    double centre_x = 0.0;
    double centre_y = 0.0;

    // Tangent unit vector, oriented downstream where possible.
    double tangent_dx = 0.0;
    double tangent_dy = 0.0;

    // Left-hand normal unit vector relative to downstream direction.
    double normal_dx = 0.0;
    double normal_dy = 0.0;

    // Compass bearings in degrees clockwise from north.
    double tangent_bearing_degrees = 0.0;
    double normal_bearing_degrees  = 0.0;

    double bearing_window_length_metres = 0.0;
    double half_length_metres           = 0.0;
    double sample_spacing_metres        = 0.0;

    lsdtt_xy start_xy;
    lsdtt_xy end_xy;

    std::vector<lsdtt_cross_section_sample> samples;
};

struct lsdtt_cross_section_options
{
    double cross_section_half_length_metres = 50.0;
    double cross_section_sample_spacing_metres = 5.0;

    // Along-channel window length used to estimate local channel bearing.
    double bearing_window_length_metres = 50.0;

    // Keep every Nth channel pixel in row-major raster order.
    int cross_section_thinning = 1;

    // If true, do not generate sections centred on junction pixels.
    bool skip_junction_pixels = false;

    // Minimum number of channel pixels needed in the fitted bearing window.
    int minimum_bearing_window_pixels = 3;
};

// -----------------------------------------------------------------------------
// extract_cross_sections
//
// dem
//   DEM raster used for elevation sampling along the section.
//
// vertical_excess_curvature_raster
//   Optional raster sampled at the same section points.
//   Pass nullptr to disable vertical excess curvature sampling.
//
// flow / network
//   Existing routed channel network infrastructure.
//
// options
//   Cross-section geometry and bearing-fit controls.
// -----------------------------------------------------------------------------
std::vector<lsdtt_cross_section> extract_cross_sections(
    const lsdtt_raster_cpu&            dem,
    const lsdtt_raster_cpu*            vertical_excess_curvature_raster,
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    const lsdtt_cross_section_options& options = lsdtt_cross_section_options{}
);

// -----------------------------------------------------------------------------
// write_cross_section_lines
//
// Writes one line per cross section.
// Intended as the simple geometry output for section locations.
// -----------------------------------------------------------------------------
void write_cross_section_lines(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_raster_meta&                meta,
    const std::string&                      out_path,
    vector_output_format                    format,
    bool                                    spatial_index = true,
    const std::string&                      layer_name = "cross_sections"
);

// -----------------------------------------------------------------------------
// write_cross_section_sample_points
//
// Writes sampled section points with elevation and, where available,
// vertical_excess_curvature.
// -----------------------------------------------------------------------------
void write_cross_section_sample_points(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_raster_meta&                meta,
    const std::string&                      out_path,
    point_output_format                     format,
    bool                                    spatial_index = true,
    const std::string&                      layer_name = "cross_section_samples"
);

// -----------------------------------------------------------------------------
// write_cross_section_samples_csv
//
// Writes sampled section points to CSV with section geometry / location,
// elevation, and optional vertical_excess_curvature values.
// -----------------------------------------------------------------------------
void write_cross_section_samples_csv(
    const std::vector<lsdtt_cross_section>& sections,
    const std::string&                      out_path
);

} // namespace lsdtt
