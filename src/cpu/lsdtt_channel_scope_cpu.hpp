// lsdtt_channel_scope_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_channel_scope_cpu
//
// Resolves which channel pixels are "in scope" for channel-linked analyses
// (floodplain/terrace classification, valley centrelines, flat-surface
// breaking). Three modes, selected implicitly from the request:
//
//   1. all_channels_threshold
//      source_flats empty. Channel heads come from the contributing-pixel
//      threshold; every channel pixel at or above the threshold is in scope.
//
//   2. user_sources
//      source_flats given, outlet_flats empty. Each source is (optionally)
//      snapped to the threshold channel network, then traced downstream to
//      the DEM outlet. One source = single-channel mode.
//
//   3. user_sources_bounded
//      source_flats and outlet_flats given. As mode 2, but each traced
//      channel is truncated at the outlet point matched to it. Outlets are
//      snapped onto the *traced flow paths* (nearest path node), not onto the
//      channel network, so a slightly misplaced outlet still truncates the
//      intended channel.
//
// Auto role detection
// -------------------
// When auto_detect_roles is set, all points arrive in source_flats (a single
// combined file) and the engine infers roles: a point that snaps strictly
// downstream of another supplied point is treated as an outlet (it truncates
// that upstream point's channel); a point downstream of none is an independent
// source / tributary head. This resolves to mode 2 (all independent) or mode 3
// (one or more downstream = bounded) without a separate outlet file.
//
// Per-pixel source id semantics
// -----------------------------
// Where two in-scope channels share downstream pixels, channel_source_id
// records the first source (in resolved_sources order) to claim the pixel.
//
// Node-list semantics differ by mode:
//   - mode 1: a source's node list stops where it meets a pixel already
//     claimed by an earlier source (link-like behaviour, keeps total storage
//     O(n_pixels) with thousands of heads).
//   - modes 2/3: node lists are complete per-source channel profiles and may
//     overlap downstream of confluences.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_point_input.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_flow_snapping_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

enum class channel_scope_mode
{
    all_channels_threshold,
    user_sources,
    user_sources_bounded
};

// -----------------------------------------------------------------------------
// channel_scope_request
//
// threshold_contributing_pixels
//   Channel-head threshold (include-self contributing pixel count). Used to
//   extract heads in mode 1 and as the snapping target in modes 2/3.
//
// source_flats / outlet_flats
//   Raw user point locations as flat indices (e.g. from
//   successfully_mapped_input_points_to_flats). Empty source_flats selects
//   mode 1, in which case outlet_flats must also be empty.
//
// snap_sources
//   When true (default), sources are snapped to the nearest cell meeting the
//   contributing-pixel threshold with an expanding-square search. When false,
//   sources are used as given (they must be valid flow cells).
//
// source_snap_max_radius_pixels
//   Maximum expanding-square search radius for source snapping.
//
// outlet_snap_max_distance_metres
//   Maximum allowed distance between a raw outlet point and its nearest
//   traced-path node. 0 or negative = unlimited. Outlets beyond the limit are
//   reported unmatched and do not truncate anything.
// -----------------------------------------------------------------------------
struct channel_scope_request
{
    int threshold_contributing_pixels = 1000;

    std::vector<std::int64_t> source_flats;
    std::vector<std::int64_t> outlet_flats;

    bool   snap_sources                    = true;
    int    source_snap_max_radius_pixels   = 25;
    double outlet_snap_max_distance_metres = 0.0;

    // Auto-classify source_flats into sources and outlets by downstream
    // membership (see "Auto role detection" above). Requires outlet_flats to be
    // empty. Sources are always snapped when this is set.
    bool   auto_detect_roles               = false;
};

// -----------------------------------------------------------------------------
// channel_scope_result
//
// channel_mask
//   Full-grid (n_rows * n_cols) mask, 1 = in-scope channel pixel.
//
// channel_source_id
//   Full-grid raster of the claiming source index (index into
//   resolved_sources), -1 where not an in-scope channel pixel.
//
// channel_node_lists
//   Per resolved source, the ordered (upstream -> downstream) flat indices of
//   its channel. See header comment for mode-dependent semantics.
//
// resolved_sources
//   Snapped (or threshold-derived) source flats actually used for tracing.
//   Sources whose snap failed are dropped; source_snaps retains the full
//   per-input record for diagnostics.
//
// resolved_outlets
//   Per resolved source, the flat index at which its channel was truncated,
//   or -1 when the channel runs to the DEM outlet.
//
// source_snaps
//   Per *input* source point, the snap diagnostic record (empty in mode 1 or
//   when snap_sources is false).
//
// outlet_matched_source / outlet_snap_distance_metres
//   Per *input* outlet point: the resolved-source index whose path it
//   truncates (-1 = unmatched) and the distance from the raw outlet location
//   to the chosen path node.
// -----------------------------------------------------------------------------
struct channel_scope_result
{
    channel_scope_mode mode = channel_scope_mode::all_channels_threshold;

    std::vector<std::uint8_t> channel_mask;
    std::vector<std::int32_t> channel_source_id;

    std::vector<std::vector<std::int64_t>> channel_node_lists;

    std::vector<std::int64_t> resolved_sources;
    std::vector<std::int64_t> resolved_outlets;

    std::vector<lsdtt_snap_result> source_snaps;

    std::vector<std::int32_t> outlet_matched_source;
    std::vector<double>       outlet_snap_distance_metres;

    // Auto role detection only: per *input* source point, the inferred role -
    // 0 = source, 1 = outlet, -1 = snap failed. Aligned with source_snaps.
    // Empty when auto_detect_roles was not used.
    std::vector<std::int32_t> input_point_role;
};

// -----------------------------------------------------------------------------
// resolve_channel_scope
//
// Requires a flow object built with s_vector enabled. Throws
// std::invalid_argument on an empty flow object, a missing s_vector, or
// outlet_flats supplied without source_flats.
// -----------------------------------------------------------------------------
[[nodiscard]] channel_scope_result resolve_channel_scope(
    const lsdtt_flow_base&       flow,
    const channel_scope_request& request
);

// =============================================================================
// Shared CLI-side plumbing
//
// Everything below is here rather than in a CLI because more than one CLI needs
// it and a function at file scope in a CLI is invisible to the rest of the
// project. lsdtt-flat-surface-finder, lsdtt-fluvial-surfaces and
// lsdtt-valley-centreline all resolve a channel scope from the same parameter
// names, write the same scope diagnostics, and write the same three channel
// network layers; before this they each carried their own copy.
// =============================================================================

// -----------------------------------------------------------------------------
// resolve_point_input_format
//
// Map a CLI *_points_input_format string onto the enum. Throws
// std::runtime_error naming the offending value if it matches nothing.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_point_input_format resolve_point_input_format(
    const std::string& format_name
);

// -----------------------------------------------------------------------------
// point_role_input
//
// One role's worth of point-input parameters, as every CLI spells them:
// <role>_points_fname, _input_format, _layer_name, _latitude_field_name,
// _longitude_field_name.
// -----------------------------------------------------------------------------
struct point_role_input
{
    std::string role_name = "source";   // for the progress message only
    std::string path;                   // empty = this role is not supplied
    std::string input_format = "auto_detect";
    std::string layer_name;
    std::string latitude_field_name  = "latitude";
    std::string longitude_field_name = "longitude";
};

// -----------------------------------------------------------------------------
// load_role_point_flats
//
// Load one role's points in the raster's CRS and return the flat indices of
// those that mapped onto valid cells. Returns empty for an empty path.
// Progress goes to stdout, as it did in each CLI.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::int64_t> load_role_point_flats(
    const point_role_input&  role,
    const lsdtt_raster_meta& meta
);

// -----------------------------------------------------------------------------
// channel_scope_diagnostic_request / write_channel_scope_diagnostics
//
// The scope-resolution diagnostics: which channels are in scope, which nodes
// they run through, and where each supplied point snapped to. These stay with
// every CLI that resolves a scope, because they answer "did the tool use the
// channels I meant?" - a question the dedicated channel CLIs cannot answer for
// someone else's run.
//
// output_stem is the full write path prefix; each writer appends its own
// suffix. Failures are reported to stderr and do not throw: a diagnostic that
// cannot be written must not lose the analysis it describes.
// -----------------------------------------------------------------------------
struct channel_scope_diagnostic_request
{
    bool write_scope_raster    = false;   // <stem>_channel_scope.tif
    bool write_scope_nodes_csv = false;   // <stem>_channel_scope_nodes.csv
    bool write_point_snap_csv  = false;   // <stem>_point_snaps.csv + 2 fgb layers

    raster_output_format raster_format = raster_output_format::cog;
};

void write_channel_scope_diagnostics(
    const channel_scope_result&             scope,
    const lsdtt_flow_base&                  flow,
    const lsdtt_raster_cpu&                 dem,
    lsdtt_data_loader&                      loader,
    const std::string&                      output_stem,
    const channel_scope_diagnostic_request& request
);

// -----------------------------------------------------------------------------
// channel_network_write_request / write_channel_network_layers
//
// The three channel-network vector layers. This lives here, with exactly one
// implementation, because more than one CLI offers the parameters.
//
// That a tool which classifies fluvial surfaces also writes the channel network
// it used is deliberate, not an oversight: a user should not have to chain a
// second CLI to see the network their surfaces were measured against. The
// overlap with lsdtt-channel-extraction is answered by sharing this writer, not
// by removing the parameter - the parameter appears in several CLIs, the code
// exists once. If you are about to write these three layers inline in a fourth
// place, call this instead.
//
// Failures are reported to stderr and do not throw.
// -----------------------------------------------------------------------------
struct channel_network_write_request
{
    bool write_lines           = true;    // <stem>_channel_network_lines<ext>
    bool write_channel_points  = false;   // <stem>_channel_network_points<ext>
    bool write_junction_points = false;   // <stem>_junction_network_points<ext>

    vector_output_format format = vector_output_format::flatgeobuf;

    // Optional per-pixel source id written as a line attribute. Empty = none.
    std::vector<std::int32_t> channel_source_id;
};

void write_channel_network_layers(
    const lsdtt_junction_network_base&   jnet,
    const std::string&                   output_stem,
    const channel_network_write_request& request
);

} // namespace lsdtt
