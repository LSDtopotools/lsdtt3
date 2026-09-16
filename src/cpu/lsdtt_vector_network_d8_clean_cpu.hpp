// lsdtt_vector_network_d8_clean_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_vector_network_d8_clean
//
// D8-guided cleaning of a clipped vector stream network so it can seed the
// junction-finding routines in lsdtt-stream-burning. Two topology problems in
// the repaired OS VectorMap water network break junction finding:
//
//   (a) Dead-ending segments. An interior degree-1 downstream tip has nowhere
//       to drain in the vector graph. We take the downstream point, walk D8
//       downhill from its DEM pixel, and add a connector to the first OTHER
//       vector segment the flow path reaches.
//
//   (b) Loops. A water polygon skeleton that splits and rejoins forms a cycle.
//       The cycle is found, its two arms are scored by mean distance to a D8
//       channel (cells whose contributing area exceeds a threshold), and the
//       arm farthest from the D8 channel is removed, leaving a single through
//       channel. Isolated rings fall back to removing the single farthest edge.
//
// The loop-tolerant builder lsdtt_vector_to_junction_network must have been run
// with break_cycles = true so loop edges are labelled (is_loop_edge) and their
// endpoints flagged (is_on_loop); no source is created on a loop.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_vector_to_junction_network.hpp"

namespace lsdtt
{

struct vector_network_d8_clean_options
{
    // Contributing-pixel threshold defining a D8 channel. A DEM cell belongs to
    // the D8 channel mask when its number of contributing pixels (include-self)
    // is at least this value. Used to score loop arms.
    std::int64_t channel_threshold_pixels = 1000;

    // Maximum D8 walk distance (metres) when reconnecting a dead-ending tip to
    // the next vector segment. If the flow path does not reach another segment
    // within this distance the dead end is left unresolved.
    double dead_end_search_max_metres = 750.0;

    // Sampling step (metres) along edge polylines for channel-distance scoring
    // and for rasterising the segment-coverage mask. If <= 0, defaults to half
    // the DEM cell size.
    double edge_sample_step_metres = -1.0;

    // Maximum radius (metres) of the bounded nearest-channel search used when
    // scoring loop arms. If <= 0, defaults to 250 m.
    double channel_search_max_metres = -1.0;

    bool resolve_loops     = true;
    bool resolve_dead_ends = true;
};

// One output segment of the cleaned network, tagged with its role.
struct vector_clean_segment
{
    std::vector<lsdtt_xy> xy;

    // "channel"             - retained network edge
    // "removed_loop_arm"    - loop edge/arm removed because it is farthest from
    //                         the D8 channel
    // "dead_end_connector"  - new D8-derived link from a dead-ending tip to the
    //                         next vector segment
    std::string role = "channel";

    std::int32_t edge_id              = -1;
    double       mean_channel_dist_m  = -1.0;   // for scored loop edges
};

struct vector_network_d8_clean_result
{
    std::vector<vector_clean_segment> segments;

    std::int64_t channel_mask_cells       = 0;
    int          n_loops_detected         = 0;
    int          n_loops_resolved         = 0;
    int          n_loop_edges_removed     = 0;
    int          n_dead_ends_detected     = 0;
    int          n_dead_ends_resolved     = 0;
    int          n_dead_ends_unresolved   = 0;
};

// Clean the vector network with D8 guidance. The network must have been built
// from `flow`'s DEM (same CRS / extent / cell size).
[[nodiscard]] vector_network_d8_clean_result clean_vector_network_with_d8(
    const lsdtt_vector_to_junction_network& network,
    const lsdtt_flow_base&                  flow,
    const vector_network_d8_clean_options&  options =
        vector_network_d8_clean_options{}
);

// Write the cleaned network to a single-layer FlatGeobuf (or gpkg).
void write_cleaned_network(
    const vector_network_d8_clean_result& result,
    const lsdtt_raster_meta&              meta,
    const std::string&                    out_path,
    vector_output_format                  format = vector_output_format::flatgeobuf,
    bool                                  spatial_index = true,
    const std::string&                    layer_name = "cleaned_network"
);

} // namespace lsdtt
