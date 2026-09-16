// lsdtt_internal_drainage_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implementation of closed-depression detection and the sink review round-trip.
// See lsdtt_internal_drainage_cpu.hpp for the contract.
// -----------------------------------------------------------------------------

#include "lsdtt_internal_drainage_cpu.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "../core/lsdtt_point_input.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "lsdtt_depression_hierarchy_cpu.hpp"
#include "lsdtt_raster_algebra.hpp"
#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_binary_mask_vectorise_cpu.hpp"

namespace lsdtt {

namespace {

// =============================================================================
// Small shared helpers
// =============================================================================

std::string to_lower_copy(const std::string& s)
{
    std::string out = s;
    for (char& c : out)
    {
        c = static_cast<char>(
            std::tolower(static_cast<unsigned char>(c)));
    }
    return out;
}

// The file extension, lower-cased, including the leading dot. Empty when the
// name has no extension after the last path separator.
std::string extension_of(const std::string& path)
{
    const std::size_t slash = path.find_last_of("/\\");
    const std::size_t start = (slash == std::string::npos) ? 0 : slash + 1;
    const std::size_t dot   = path.find_last_of('.');

    if (dot == std::string::npos || dot < start) return "";

    return to_lower_copy(path.substr(dot));
}

// Renders a double at the SHORTEST precision that reads back to exactly the
// same value. A CSV round-trip is then exact rather than merely close, which is
// what lets the CSV and vector candidate files parse back to identical records,
// while keeping the file readable for a human reviewer.
std::string exact_double_to_string(double value)
{
    char buffer[64];

    for (int precision = 6; precision < 18; ++precision)
    {
        std::snprintf(buffer, sizeof(buffer), "%.*g", precision, value);
        if (std::strtod(buffer, nullptr) == value) return std::string(buffer);
    }

    std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    return std::string(buffer);
}

// The schema of the candidate layer, in the order the fields are written. The
// CSV header and the vector schema are generated from this one list so that the
// two formats cannot drift apart.
const std::vector<std::string>& candidate_field_names()
{
    static const std::vector<std::string> names = {
        "sink_id",
        "routing_sink_id",
        "pit_row",
        "pit_col",
        "pit_flat",
        "easting",
        "northing",
        "latitude",
        "longitude",
        "pit_elevation",
        "spill_elevation",
        "depth_m",
        "area_m2",
        "volume_m3",
        "n_pixels",
        "tagged_drainage_area_m2",
        "tagged_n_pixels",
        "hierarchy_level",
        "parent_sink_id",
        "accepted"
    };
    return names;
}

const std::vector<lsdtt_vector_field_type>& candidate_field_types()
{
    static const std::vector<lsdtt_vector_field_type> types = {
        lsdtt_vector_field_type::int64_type,     // sink_id
        lsdtt_vector_field_type::int64_type,     // routing_sink_id
        lsdtt_vector_field_type::int64_type,     // pit_row
        lsdtt_vector_field_type::int64_type,     // pit_col
        lsdtt_vector_field_type::int64_type,     // pit_flat
        lsdtt_vector_field_type::double_type,    // easting
        lsdtt_vector_field_type::double_type,    // northing
        lsdtt_vector_field_type::double_type,    // latitude
        lsdtt_vector_field_type::double_type,    // longitude
        lsdtt_vector_field_type::double_type,    // pit_elevation
        lsdtt_vector_field_type::double_type,    // spill_elevation
        lsdtt_vector_field_type::double_type,    // depth_m
        lsdtt_vector_field_type::double_type,    // area_m2
        lsdtt_vector_field_type::double_type,    // volume_m3
        lsdtt_vector_field_type::int64_type,     // n_pixels
        lsdtt_vector_field_type::double_type,    // tagged_drainage_area_m2
        lsdtt_vector_field_type::int64_type,     // tagged_n_pixels
        lsdtt_vector_field_type::int64_type,     // hierarchy_level
        lsdtt_vector_field_type::int64_type,     // parent_sink_id
        lsdtt_vector_field_type::int64_type      // accepted
    };
    return types;
}

// One candidate as attribute values, in schema order. Elevations are floats in
// the struct and widen to double exactly, so they narrow back exactly too.
std::vector<lsdtt_vector_field_value>
candidate_field_values(const sink_candidate& sink)
{
    return {
        static_cast<std::int64_t>(sink.sink_id),
        static_cast<std::int64_t>(sink.routing_sink_id),
        static_cast<std::int64_t>(sink.pit_row),
        static_cast<std::int64_t>(sink.pit_col),
        sink.pit_flat,
        sink.pit_x,
        sink.pit_y,
        sink.latitude,
        sink.longitude,
        static_cast<double>(sink.pit_elevation),
        static_cast<double>(sink.spill_elevation),
        static_cast<double>(sink.depth_m),
        sink.area_m2,
        sink.volume_m3,
        static_cast<std::int64_t>(sink.n_pixels),
        sink.tagged_drainage_area_m2,
        sink.tagged_n_pixels,
        static_cast<std::int64_t>(sink.hierarchy_level),
        static_cast<std::int64_t>(sink.parent_sink_id),
        static_cast<std::int64_t>(sink.accepted ? 1 : 0)
    };
}

// Builds a transform from the raster CRS to WGS84 so that pit latitude and
// longitude go through the OGR transform rather than a hand-rolled projection.
// Returns false when the raster carries no CRS, in which case lat/lon are left
// at zero rather than invented.
bool make_raster_to_wgs84_transform(const lsdtt_raster_meta&    meta,
                                    lsdtt_coordinate_transform& transform_out)
{
    if (meta.crs_wkt.empty()) return false;

    transform_out.initialise_from_epsg(meta.crs_wkt, 4326);
    return true;
}

// The inverse: WGS84 latitude/longitude supplied by a reviewer, brought into
// the raster CRS.
bool make_wgs84_to_raster_transform(const lsdtt_raster_meta&    meta,
                                    lsdtt_coordinate_transform& transform_out)
{
    if (meta.crs_wkt.empty()) return false;

    transform_out.initialise_from_source_epsg(4326, meta.crs_wkt);
    return true;
}

// Fills pit_x / pit_y / latitude / longitude from the pit row and column. The
// pit is a CELL, so its position is the cell centre.
void set_pit_position(sink_candidate&                   sink,
                      const lsdtt_raster_meta&          meta,
                      const lsdtt_coordinate_transform* to_wgs84)
{
    sink.pit_x = cell_centre_x(meta, sink.pit_col);
    sink.pit_y = cell_centre_y(meta, sink.pit_row);
    sink.pit_flat =
        static_cast<std::int64_t>(sink.pit_row) * meta.n_cols + sink.pit_col;

    if (to_wgs84 != nullptr)
    {
        const lsdtt_lonlat ll =
            cell_center_to_lonlat(meta, sink.pit_row, sink.pit_col, *to_wgs84);
        sink.longitude = ll.longitude;
        sink.latitude  = ll.latitude;
    }
}

// =============================================================================
// Attribute lookup
//
// The candidate file may have been round-tripped through a GIS or a
// spreadsheet, so field names are matched case-insensitively and a short list
// of aliases is accepted for the fields a reviewer is most likely to rename or
// hand-create.
// =============================================================================

using attribute_map = std::map<std::string, lsdtt_point_attribute_value>;

const lsdtt_point_attribute_value*
find_attribute(const attribute_map&            attributes,
               const std::vector<std::string>& aliases)
{
    for (const std::string& alias : aliases)
    {
        for (const auto& entry : attributes)
        {
            if (to_lower_copy(entry.first) == alias &&
                entry.second.is_present())
            {
                return &entry.second;
            }
        }
    }
    return nullptr;
}

bool attribute_as_double(const attribute_map&            attributes,
                         const std::vector<std::string>& aliases,
                         double&                         value_out)
{
    const lsdtt_point_attribute_value* value = find_attribute(attributes, aliases);
    if (value == nullptr) return false;
    return value->as_double(value_out);
}

bool attribute_as_int(const attribute_map&            attributes,
                      const std::vector<std::string>& aliases,
                      int&                            value_out)
{
    const lsdtt_point_attribute_value* value = find_attribute(attributes, aliases);
    if (value == nullptr) return false;

    std::int64_t as_int = 0;
    if (!value->as_int64(as_int)) return false;

    value_out = static_cast<int>(as_int);
    return true;
}

bool attribute_as_int64(const attribute_map&            attributes,
                        const std::vector<std::string>& aliases,
                        std::int64_t&                   value_out)
{
    const lsdtt_point_attribute_value* value = find_attribute(attributes, aliases);
    return value != nullptr && value->as_int64(value_out);
}

// =============================================================================
// Depression-hierarchy detection state
//
// Everything the hierarchy method needs, computed once, so that detection and
// the depression polygons share one derivation rather than two that could
// disagree.
// =============================================================================

struct hierarchy_state
{
    lsdtt_depression_hierarchy deps;
    std::vector<dep_label_t>   label;      // leaf depression label per cell

    std::vector<std::int64_t>  pit_cell;   // per depression, deterministic
    std::vector<float>         pit_elev;
    std::vector<int>           level;      // 0 = leaf
    std::vector<char>          selected;   // survives leaves_only + the filters
};

// True when a depression node is an elementary (leaf) depression.
bool is_leaf(const depression_node& node)
{
    return node.lchild == k_no_label && node.rchild == k_no_label;
}

// Resolves the pit of every depression deterministically: the lowest cell, and
// among equally low cells the LOWEST FLAT INDEX. A flat depression floor has no
// unique lowest cell, so without an explicit tie-break the reported pit would
// depend on the order the flood happened to visit the floor.
void resolve_deterministic_pits(const lsdtt_raster_cpu& dem,
                                hierarchy_state&        state)
{
    const std::size_t n_deps = state.deps.size();

    state.pit_cell.assign(n_deps, -1);
    state.pit_elev.assign(n_deps, std::numeric_limits<float>::infinity());

    const std::vector<float>&   elevation = dem.get_data();
    const std::vector<uint8_t>& mask      = dem.get_mask();

    const std::int64_t n_cells =
        static_cast<std::int64_t>(dem.get_meta().n_rows) *
        dem.get_meta().n_cols;

    // Pass one: the lowest cell of each leaf label.
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (mask[static_cast<std::size_t>(i)] == 0u) continue;

        const dep_label_t cell_label = state.label[static_cast<std::size_t>(i)];
        if (cell_label == k_ocean_label) continue;
        if (static_cast<std::size_t>(cell_label) >= n_deps) continue;

        const float z = elevation[static_cast<std::size_t>(i)];
        float& best_z = state.pit_elev[cell_label];
        std::int64_t& best_cell = state.pit_cell[cell_label];

        if (best_cell < 0 || z < best_z || (z == best_z && i < best_cell))
        {
            best_z    = z;
            best_cell = i;
        }
    }

    // Pass two: a meta-depression's pit is the lowest pit of its children. A
    // parent always has a larger index than its children, so one forward sweep
    // is enough.
    for (std::size_t d = 1; d < n_deps; ++d)
    {
        const depression_node& node = state.deps[d];
        if (is_leaf(node)) continue;

        const dep_label_t children[2] = { node.lchild, node.rchild };

        for (const dep_label_t child : children)
        {
            if (child == k_no_label) continue;

            const std::int64_t child_cell = state.pit_cell[child];
            if (child_cell < 0) continue;

            const float child_z = state.pit_elev[child];

            if (state.pit_cell[d] < 0 ||
                child_z < state.pit_elev[d] ||
                (child_z == state.pit_elev[d] && child_cell < state.pit_cell[d]))
            {
                state.pit_elev[d] = child_z;
                state.pit_cell[d] = child_cell;
            }
        }
    }
}

// hierarchy_level: 0 for a leaf, otherwise one more than the deepest child.
void resolve_hierarchy_levels(hierarchy_state& state)
{
    state.level.assign(state.deps.size(), 0);

    for (std::size_t d = 1; d < state.deps.size(); ++d)
    {
        const depression_node& node = state.deps[d];
        if (is_leaf(node)) continue;

        int deepest = 0;
        const dep_label_t children[2] = { node.lchild, node.rchild };

        for (const dep_label_t child : children)
        {
            if (child == k_no_label) continue;
            deepest = std::max(deepest, state.level[child] + 1);
        }

        state.level[d] = deepest;
    }
}

// Builds the hierarchy, the deterministic pits, the tree levels and the
// selection of depressions that pass leaves_only and the size filters.
hierarchy_state build_hierarchy_state(const lsdtt_raster_cpu&       dem,
                                      const sink_detection_options& opts)
{
    hierarchy_state state;

    depression_hierarchy_options dh_opts;
    dh_opts.boundary_conditions = opts.boundary_conditions;

    std::vector<uint8_t> flow_direction;   // needed by fill-spill-merge, not here

    state.deps = build_depression_hierarchy(
        dem, state.label, flow_direction, dh_opts);

    resolve_deterministic_pits(dem, state);
    resolve_hierarchy_levels(state);

    const double cell_area = dem.get_meta().pixel_area();

    state.selected.assign(state.deps.size(), 0);

    // Index 0 is the ocean node and is never a sink.
    for (std::size_t d = 1; d < state.deps.size(); ++d)
    {
        const depression_node& node = state.deps[d];

        if (opts.leaves_only && !is_leaf(node)) continue;
        if (state.pit_cell[d] < 0)              continue;

        // A depression with no finite spill has no depth and cannot be
        // described; this only happens for degenerate nodes at the top of the
        // tree.
        if (!std::isfinite(node.out_elev)) continue;

        const double depth   = static_cast<double>(node.out_elev) -
                               static_cast<double>(state.pit_elev[d]);
        const double area    = static_cast<double>(node.cell_count) * cell_area;
        const double volume  = node.dep_vol * cell_area;

        if (static_cast<int>(node.cell_count) < opts.min_pixels_in_patch) continue;
        if (depth  < static_cast<double>(opts.min_depth_m))               continue;
        if (area   < opts.min_area_m2)                                    continue;
        if (volume < opts.min_volume_m3)                                  continue;

        state.selected[d] = 1;
    }

    return state;
}

// The nearest selected proper ancestor of a depression, or k_no_label.
//
// build_depression_hierarchy() sets depression_node::parent for OCEAN LINKS as
// well as for true containment, and depression_node documents ocean_linked as
// depressions that reach the ocean THROUGH this one but which are NOT
// sub-depressions of it. Following such a link would report two unrelated
// craters as nested, so a node flagged ocean_parent ends the walk: crossing an
// ocean link means there is no containing ancestor. This is the same test the
// cell walk in sink_depressions_to_polygons() makes, and the two must stay
// consistent -- do not "simplify" it back to a plain parent walk.
dep_label_t nearest_selected_ancestor(const hierarchy_state& state,
                                      std::size_t            d)
{
    dep_label_t current = static_cast<dep_label_t>(d);

    while (current != k_no_parent && current != k_ocean_label)
    {
        if (state.deps[current].ocean_parent) return k_no_label;

        const dep_label_t parent = state.deps[current].parent;

        if (parent == k_no_parent || parent == k_ocean_label) return k_no_label;
        if (state.selected[parent] != 0) return parent;

        current = parent;
    }

    return k_no_label;
}

// =============================================================================
// Fill-difference detection state
//
// The legacy path: filled minus raw, 8-connected components with a minimum
// size, one sink per component. The component labelling, area and volume all
// come from depressions_to_polygons(), so there is exactly one connected
// components pass definition in the codebase and this method cannot drift away
// from the depression polygons written by lsdtt-raster-preprocessing.
// =============================================================================

struct fill_difference_component
{
    std::int64_t deepest_cell = -1;    // legacy pit: the DEEPEST-FILL cell
    float        depth_m      = 0.0f;
    float        pit_elev     = 0.0f;
    float        spill_elev   = 0.0f;
    std::int64_t n_pixels     = 0;
    double       area_m2      = 0.0;
    double       volume_m3    = 0.0;
};

struct fill_difference_state
{
    lsdtt_raster_cpu                          labelled;
    std::map<int, fill_difference_component>  components;
    depression_polygons_result                polygons;
};

fill_difference_state
build_fill_difference_state(const lsdtt_raster_cpu&       dem,
                            const sink_detection_options& opts)
{
    fill_difference_state state;

    if (opts.min_pixels_in_patch < 1)
    {
        throw std::invalid_argument(
            "lsdtt::find_internal_drainage: min_pixels_in_patch must be >= 1");
    }

    // min_slope is zero so that the fill difference is the depression depth
    // itself, with no epsilon gradient laid across the flats to confuse the
    // deepest-fill cell.
    const lsdtt_raster_cpu filled =
        fill_pits(dem, 0.0f, opts.boundary_conditions, opts.n_threads);

    const lsdtt_raster_cpu fill_difference =
        algebra::raster_subtract(filled, dem);

    // The component definition is "any cell the fill raised", exactly as
    // depressions_to_polygons() uses it below. min_depth_m is a filter on the
    // SINK, applied further down, not on which cells make up a depression.
    const lsdtt_raster_cpu binary =
        algebra::raster_binary_true_above_but_excluding_threshold(
            fill_difference, 0.0f);

    state.labelled = label_binary_components_by_size(
        binary, opts.min_pixels_in_patch, pixel_connectivity::d8);

    state.polygons = depressions_to_polygons(
        fill_difference, opts.min_pixels_in_patch, 0.0f);

    // Per-component extremum: the deepest-fill cell, ties broken on the lowest
    // flat index so that a flat-bottomed depression gives a reproducible pit.
    const std::vector<float>&   label_data = state.labelled.get_data();
    const std::vector<uint8_t>& label_mask = state.labelled.get_mask();
    const std::vector<float>&   difference = fill_difference.get_data();
    const std::vector<float>&   elevation  = dem.get_data();
    const std::vector<float>&   filled_z   = filled.get_data();

    const std::int64_t n_cells =
        static_cast<std::int64_t>(dem.get_meta().n_rows) * dem.get_meta().n_cols;

    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t index = static_cast<std::size_t>(i);
        if (label_mask[index] == 0u) continue;

        const int component_label = static_cast<int>(label_data[index]);
        if (component_label <= 0) continue;

        fill_difference_component& component = state.components[component_label];
        ++component.n_pixels;

        const float depth = difference[index];

        if (component.deepest_cell < 0 || depth > component.depth_m)
        {
            component.deepest_cell = i;
            component.depth_m      = depth;
            component.pit_elev     = elevation[index];
            component.spill_elev   = filled_z[index];
        }
    }

    // Area and volume come from the polygon records so that they are defined in
    // exactly one place.
    for (const lsdtt_vector_record& record : state.polygons.records)
    {
        if (record.field_values.size() < 3) continue;

        const std::int64_t component_label =
            std::get<std::int64_t>(record.field_values[0]);

        auto found = state.components.find(static_cast<int>(component_label));
        if (found == state.components.end()) continue;

        found->second.area_m2   = std::get<double>(record.field_values[1]);
        found->second.volume_m3 = std::get<double>(record.field_values[2]);
    }

    return state;
}

} // anonymous namespace

// =============================================================================
// find_internal_drainage
// =============================================================================
std::vector<sink_candidate>
find_internal_drainage(const lsdtt_raster_cpu&       dem,
                       const sink_detection_options& opts)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    meta.check_projected_metres();

    lsdtt_coordinate_transform to_wgs84;
    const bool have_wgs84 = make_raster_to_wgs84_transform(meta, to_wgs84);
    const lsdtt_coordinate_transform* transform =
        have_wgs84 ? &to_wgs84 : nullptr;

    std::vector<sink_candidate> sinks;

    if (opts.method == sink_detection_method::depression_hierarchy)
    {
        const hierarchy_state state = build_hierarchy_state(dem, opts);

        // sink_id runs in ascending depression label, so a child always gets a
        // smaller id than its parent.
        std::unordered_map<dep_label_t, int> sink_id_of_depression;

        for (std::size_t d = 1; d < state.deps.size(); ++d)
        {
            if (state.selected[d] == 0) continue;

            const depression_node& node = state.deps[d];

            sink_candidate sink;
            sink.sink_id         = static_cast<int>(sinks.size());
            sink.pit_row         = static_cast<int>(state.pit_cell[d] / meta.n_cols);
            sink.pit_col         = static_cast<int>(state.pit_cell[d] % meta.n_cols);
            sink.pit_elevation   = state.pit_elev[d];
            sink.spill_elevation = node.out_elev;
            sink.depth_m         = node.out_elev - state.pit_elev[d];
            sink.n_pixels        = static_cast<int>(node.cell_count);
            sink.area_m2         = static_cast<double>(node.cell_count) *
                                   meta.pixel_area();
            sink.volume_m3       = node.dep_vol * meta.pixel_area();
            sink.hierarchy_level = state.level[d];
            sink.accepted        = true;

            set_pit_position(sink, meta, transform);

            sink_id_of_depression[static_cast<dep_label_t>(d)] = sink.sink_id;
            sinks.push_back(sink);
        }

        // Nesting. A depression whose parent was filtered out is attached to
        // its nearest surviving ancestor, so the reported tree is the tree of
        // the depressions that were actually emitted.
        std::size_t emitted = 0;
        for (std::size_t d = 1; d < state.deps.size(); ++d)
        {
            if (state.selected[d] == 0) continue;

            const dep_label_t ancestor = nearest_selected_ancestor(state, d);
            if (ancestor != k_no_label)
            {
                sinks[emitted].parent_sink_id = sink_id_of_depression[ancestor];
            }
            ++emitted;
        }
    }
    else
    {
        const fill_difference_state state = build_fill_difference_state(dem, opts);

        for (const auto& entry : state.components)
        {
            const fill_difference_component& component = entry.second;

            if (component.deepest_cell < 0) continue;
            if (component.depth_m   < opts.min_depth_m)   continue;
            if (component.area_m2   < opts.min_area_m2)   continue;
            if (component.volume_m3 < opts.min_volume_m3) continue;

            sink_candidate sink;
            sink.sink_id         = static_cast<int>(sinks.size());
            sink.pit_row         = static_cast<int>(component.deepest_cell / meta.n_cols);
            sink.pit_col         = static_cast<int>(component.deepest_cell % meta.n_cols);
            sink.pit_elevation   = component.pit_elev;
            sink.spill_elevation = component.spill_elev;
            sink.depth_m         = component.depth_m;
            sink.n_pixels        = static_cast<int>(component.n_pixels);
            sink.area_m2         = component.area_m2;
            sink.volume_m3       = component.volume_m3;

            // A fill-difference component is a single flat-topped depression:
            // it has no tree, so every sink is a leaf with no parent.
            sink.hierarchy_level = 0;
            sink.parent_sink_id  = -1;
            sink.accepted        = true;

            set_pit_position(sink, meta, transform);

            sinks.push_back(sink);
        }
    }

    return sinks;
}

// =============================================================================
// resolve_sink_io_format
// =============================================================================
sink_io_format
resolve_sink_io_format(const std::string& path, const std::string& requested)
{
    const std::string request = to_lower_copy(requested);

    if (!request.empty() && request != "same_as_extension")
    {
        if (request == "csv")                                 return sink_io_format::csv;
        if (request == "fgb" || request == "flatgeobuf")      return sink_io_format::flatgeobuf;
        if (request == "gpkg" || request == "geopackage")     return sink_io_format::gpkg;

        throw std::invalid_argument(
            "lsdtt::resolve_sink_io_format: unknown format '" + requested +
            "'. Use csv, flatgeobuf, gpkg, or same_as_extension.");
    }

    const std::string extension = extension_of(path);

    if (extension == ".csv")  return sink_io_format::csv;
    if (extension == ".fgb")  return sink_io_format::flatgeobuf;
    if (extension == ".gpkg") return sink_io_format::gpkg;

    throw std::invalid_argument(
        "lsdtt::resolve_sink_io_format: cannot resolve a sink file format from "
        "the extension of '" + path + "'. Supported extensions are .csv, .fgb "
        "and .gpkg; alternatively name the format explicitly.");
}

// =============================================================================
// write_sink_candidates
// =============================================================================
void write_sink_candidates(lsdtt_data_loader&                 loader,
                           const std::string&                 path,
                           const std::vector<sink_candidate>& sinks,
                           const lsdtt_raster_meta&           meta,
                           sink_io_format                     format,
                           const std::string&                 layer_name,
                           bool                               spatial_index)
{
    if (format == sink_io_format::csv)
    {
        std::vector<std::vector<std::string>> rows;
        rows.reserve(sinks.size());

        for (const sink_candidate& sink : sinks)
        {
            const std::vector<lsdtt_vector_field_value> values =
                candidate_field_values(sink);

            std::vector<std::string> row;
            row.reserve(values.size());

            for (const lsdtt_vector_field_value& value : values)
            {
                if (std::holds_alternative<std::int64_t>(value))
                {
                    row.push_back(std::to_string(std::get<std::int64_t>(value)));
                }
                else
                {
                    row.push_back(exact_double_to_string(std::get<double>(value)));
                }
            }

            rows.push_back(std::move(row));
        }

        loader.write_csv(path, candidate_field_names(), rows);
        return;
    }

    lsdtt_vector_schema schema;
    schema.field_names = candidate_field_names();
    schema.field_types = candidate_field_types();

    const vector_output_format vector_format =
        (format == sink_io_format::gpkg) ? vector_output_format::gpkg
                                         : vector_output_format::flatgeobuf;

    loader.begin_vector_write(
        path, layer_name, schema, meta.crs_wkt, vector_format, spatial_index);

    for (const sink_candidate& sink : sinks)
    {
        loader.write_point(sink.pit_x, sink.pit_y, candidate_field_values(sink));
    }

    loader.finish_vector_write();
}

// =============================================================================
// read_sink_candidates
// =============================================================================
std::vector<sink_candidate>
read_sink_candidates(const std::string&       path,
                     const lsdtt_raster_meta& meta,
                     sink_io_format           format,
                     const std::string&       layer_name)
{
    lsdtt_point_input_options input_options;
    input_options.vector_layer_name = layer_name;
    input_options.input_format =
        (format == sink_io_format::csv)
            ? lsdtt_point_input_format::csv_latitude_longitude_epsg4326
            : lsdtt_point_input_format::ogr_vector_layer;

    // An empty field list reads every field in the layer, which is what a
    // reviewer-edited file needs: extra columns are ignored and missing ones
    // simply come back absent.
    const lsdtt_point_attribute_read_result read_result =
        load_point_layer_with_attributes(
            path, meta, std::vector<std::string>{}, input_options);

    lsdtt_coordinate_transform to_wgs84;
    const bool have_wgs84 = make_raster_to_wgs84_transform(meta, to_wgs84);
    const lsdtt_coordinate_transform* transform =
        have_wgs84 ? &to_wgs84 : nullptr;

    // Built lazily: a candidate file positioned by latitude/longitude
    // attributes is the exception, not the rule.
    lsdtt_coordinate_transform from_wgs84;
    bool from_wgs84_ready = false;

    std::vector<sink_candidate> sinks;
    sinks.reserve(read_result.records.size());

    int highest_sink_id = -1;
    std::vector<std::size_t> sinks_needing_an_id;

    // Skipping is silent per feature but not silent overall: a reviewer who
    // loses a candidate off the edge of the DEM is told how many went.
    int n_skipped_off_raster = 0;

    for (const lsdtt_input_point_with_attributes& record : read_result.records)
    {
        const attribute_map& attributes = record.attributes;

        double x = 0.0;
        double y = 0.0;
        bool   have_position = false;

        // Tier 1: the point geometry of a vector layer. A CSV has no geometry
        // of its own -- its latitude/longitude columns are handled at tier 3 --
        // so this applies to OGR input only.
        if (record.point.source_format == lsdtt_point_input_format::ogr_vector_layer &&
            record.point.status == lsdtt_point_input_status::success)
        {
            x = record.point.x;
            y = record.point.y;
            have_position = true;
        }

        // Tier 2: easting / northing attributes, already in the raster CRS.
        if (!have_position)
        {
            double easting  = 0.0;
            double northing = 0.0;

            if (attribute_as_double(attributes, { "easting", "pit_x", "x" }, easting) &&
                attribute_as_double(attributes, { "northing", "pit_y", "y" }, northing))
            {
                x = easting;
                y = northing;
                have_position = true;
            }
        }

        // Tier 3: latitude / longitude. For a CSV these are the columns the
        // loader has already transformed into the raster CRS; for a vector
        // layer they are plain attributes and need transforming here.
        if (!have_position &&
            record.point.source_format ==
                lsdtt_point_input_format::csv_latitude_longitude_epsg4326 &&
            record.point.status == lsdtt_point_input_status::success)
        {
            x = record.point.x;
            y = record.point.y;
            have_position = true;
        }

        if (!have_position)
        {
            double latitude  = 0.0;
            double longitude = 0.0;

            if (attribute_as_double(attributes, { "latitude", "lat" }, latitude) &&
                attribute_as_double(attributes, { "longitude", "lon", "long", "lng" },
                                    longitude))
            {
                if (!from_wgs84_ready)
                {
                    from_wgs84_ready = make_wgs84_to_raster_transform(meta, from_wgs84);
                }

                if (from_wgs84_ready)
                {
                    x = longitude;
                    y = latitude;
                    from_wgs84.transform_in_place(x, y);
                    have_position = true;
                }
            }
        }

        sink_candidate sink;

        // Tier 4: row / col attributes, used as they stand.
        if (have_position)
        {
            // col_from_x / row_from_y THROW for a point off the grid, so the
            // extent is tested first: the header promises that a feature which
            // cannot be placed on the raster is skipped, and one candidate
            // dragged half a cell off the DEM in QGIS must not abort the read
            // of every other feature in the file. Tier 1 positions are already
            // safe because the loader flags them outside_raster_extent, and
            // this is the same test it makes there.
            if (!point_within_raster(meta, x, y))
            {
                ++n_skipped_off_raster;
                continue;
            }

            try
            {
                sink.pit_col = col_from_x(meta, x);
                sink.pit_row = row_from_y(meta, y);
            }
            catch (const std::exception&)
            {
                // A point exactly on the far edge rounds onto the first cell
                // outside the grid; treat it as off the raster, not an error.
                ++n_skipped_off_raster;
                continue;
            }
        }
        else
        {
            if (!attribute_as_int(attributes, { "pit_row", "row" }, sink.pit_row) ||
                !attribute_as_int(attributes, { "pit_col", "col" }, sink.pit_col))
            {
                // Nothing in this feature says where it is; skip it rather than
                // invent a location.
                ++n_skipped_off_raster;
                continue;
            }
        }

        if (sink.pit_row < 0 || sink.pit_row >= meta.n_rows ||
            sink.pit_col < 0 || sink.pit_col >= meta.n_cols)
        {
            ++n_skipped_off_raster;
            continue;
        }

        // A sink is a cell, so the resolved position is snapped to that cell's
        // centre and its latitude/longitude recomputed. This is what makes the
        // CSV and vector round-trips agree exactly whichever tier positioned
        // the feature.
        set_pit_position(sink, meta, transform);

        double value = 0.0;

        if (attribute_as_double(attributes, { "pit_elevation" }, value))
            sink.pit_elevation = static_cast<float>(value);
        if (attribute_as_double(attributes, { "spill_elevation" }, value))
            sink.spill_elevation = static_cast<float>(value);
        if (attribute_as_double(attributes, { "depth_m", "depth" }, value))
            sink.depth_m = static_cast<float>(value);
        if (attribute_as_double(attributes, { "area_m2", "area" }, value))
            sink.area_m2 = value;
        if (attribute_as_double(attributes, { "volume_m3", "volume" }, value))
            sink.volume_m3 = value;
        if (attribute_as_double(attributes, { "tagged_drainage_area_m2" }, value))
            sink.tagged_drainage_area_m2 = value;

        int integer_value = 0;

        if (attribute_as_int(attributes, { "n_pixels" }, integer_value))
            sink.n_pixels = integer_value;
        std::int64_t int64_value = 0;
        if (attribute_as_int64(attributes, { "tagged_n_pixels" }, int64_value))
            sink.tagged_n_pixels = int64_value;
        if (attribute_as_int(attributes, { "hierarchy_level" }, integer_value))
            sink.hierarchy_level = integer_value;
        if (attribute_as_int(attributes, { "parent_sink_id" }, integer_value))
            sink.parent_sink_id = integer_value;
        if (attribute_as_int(attributes, { "routing_sink_id" }, integer_value))
            sink.routing_sink_id = integer_value;

        // An absent accepted column means every feature present is accepted:
        // the reviewer deleted the ones they rejected.
        const lsdtt_point_attribute_value* accepted_value =
            find_attribute(attributes, { "accepted", "accept", "valid" });

        sink.accepted = true;
        if (accepted_value != nullptr)
        {
            bool accepted = true;
            if (accepted_value->as_bool(accepted)) sink.accepted = accepted;
        }

        if (attribute_as_int(attributes, { "sink_id", "id" }, integer_value))
        {
            sink.sink_id    = integer_value;
            highest_sink_id = std::max(highest_sink_id, integer_value);
        }
        else
        {
            sinks_needing_an_id.push_back(sinks.size());
        }

        sinks.push_back(sink);
    }

    // A hand-digitised feature carries no sink_id, so it gets a fresh one
    // continuing from the largest id already in the file. Sort by cell first:
    // indexed vector formats may reorder features, but distinct digitised points
    // must receive the same durable ids regardless of storage order.
    std::sort(sinks_needing_an_id.begin(), sinks_needing_an_id.end(),
              [&](std::size_t a, std::size_t b)
              {
                  return sinks[a].pit_flat < sinks[b].pit_flat;
              });
    for (const std::size_t index : sinks_needing_an_id)
    {
        sinks[index].sink_id = ++highest_sink_id;
    }

    if (n_skipped_off_raster > 0)
    {
        std::cerr << "Warning: read_sink_candidates: " << n_skipped_off_raster
                  << " of " << read_result.records.size()
                  << " features in '" << path
                  << "' could not be placed on the raster and were skipped\n";
    }

    return sinks;
}

// =============================================================================
// accepted_sinks
// =============================================================================
std::vector<sink_candidate>
accepted_sinks(const std::vector<sink_candidate>& sinks)
{
    std::vector<sink_candidate> accepted;

    for (const sink_candidate& sink : sinks)
    {
        if (!sink.accepted) continue;
        accepted.push_back(sink);
    }

    std::sort(accepted.begin(), accepted.end(),
              [](const sink_candidate& a, const sink_candidate& b)
              {
                  return a.sink_id < b.sink_id;
              });

    for (std::size_t i = 0; i < accepted.size(); ++i)
    {
        if (accepted[i].sink_id < 0)
        {
            throw std::invalid_argument(
                "accepted_sinks: accepted sinks must have non-negative sink_id values");
        }
        if (i > 0 && accepted[i - 1].sink_id == accepted[i].sink_id)
        {
            throw std::invalid_argument(
                "accepted_sinks: accepted sinks must have unique sink_id values");
        }

        accepted[i].routing_sink_id = static_cast<int>(i);
    }

    return accepted;
}

// =============================================================================
// sink_depressions_to_polygons
// =============================================================================
depression_polygons_result
sink_depressions_to_polygons(const lsdtt_raster_cpu&            dem,
                             const std::vector<sink_candidate>& sinks,
                             const sink_detection_options&      opts)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    meta.check_projected_metres();

    depression_polygons_result result;
    result.schema.field_names = {
        "sink_id", "pit_elevation", "spill_elevation", "depth_m",
        "area_m2", "volume_m3", "n_pixels", "hierarchy_level",
        "parent_sink_id", "tagged_drainage_area_m2", "tagged_n_pixels",
        "accepted"
    };
    result.schema.field_types = {
        lsdtt_vector_field_type::int64_type,     // sink_id
        lsdtt_vector_field_type::double_type,    // pit_elevation
        lsdtt_vector_field_type::double_type,    // spill_elevation
        lsdtt_vector_field_type::double_type,    // depth_m
        lsdtt_vector_field_type::double_type,    // area_m2
        lsdtt_vector_field_type::double_type,    // volume_m3
        lsdtt_vector_field_type::int64_type,     // n_pixels
        lsdtt_vector_field_type::int64_type,     // hierarchy_level
        lsdtt_vector_field_type::int64_type,     // parent_sink_id
        lsdtt_vector_field_type::double_type,    // tagged_drainage_area_m2
        lsdtt_vector_field_type::int64_type,     // tagged_n_pixels
        lsdtt_vector_field_type::int64_type      // accepted
    };

    if (sinks.empty()) return result;
    if (sinks.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument(
            "sink_depressions_to_polygons: too many sinks for integer component labels");
    }

    // The attributes of a polygon come from its sink, never from a second
    // measurement, so the point and polygon layers cannot disagree.
    const auto attributes_of = [](const sink_candidate& sink)
    {
        return std::vector<lsdtt_vector_field_value>{
            static_cast<std::int64_t>(sink.sink_id),
            static_cast<double>(sink.pit_elevation),
            static_cast<double>(sink.spill_elevation),
            static_cast<double>(sink.depth_m),
            sink.area_m2,
            sink.volume_m3,
            static_cast<std::int64_t>(sink.n_pixels),
            static_cast<std::int64_t>(sink.hierarchy_level),
            static_cast<std::int64_t>(sink.parent_sink_id),
            sink.tagged_drainage_area_m2,
            sink.tagged_n_pixels,
            static_cast<std::int64_t>(sink.accepted ? 1 : 0)
        };
    };

    // Footprint geometry is traced by the shared component vectoriser, which
    // wants a dense label raster whose values are 1-based. Durable sink_id
    // values may be sparse, so geometry uses local vector indices and only the
    // output attributes carry sink_id.
    const auto trace_footprints =
        [&](const std::vector<int>& sink_of_cell,
            const std::vector<const sink_candidate*>& sink_by_index)
    {
        lsdtt_raster_meta label_meta = meta;
        label_meta.no_data_value = 0.0;

        std::vector<float>   label_data(sink_of_cell.size(), 0.0f);
        std::vector<uint8_t> label_mask(sink_of_cell.size(), 0u);

        // A sink's cells need not be CONNECTED. When min_pixels_in_patch drops a
        // leaf, its cells reattach to the nearest selected ancestor, and that
        // ancestor's footprint then has an outlying patch separated from its main
        // body by the footprints of its other children.
        //
        // binary_components_to_polygon_records_with_holes() returns one record
        // per LABEL, so labelling cells by sink id would collapse those patches
        // and silently discard all but one -- terrain that is plainly on the DEM
        // would appear in no footprint at all.
        //
        // So each 8-connected COMPONENT of each sink gets its own label here, and
        // component_sink maps back afterwards. A sink with a multipart footprint
        // therefore emits one record per part, all carrying its attributes.
        std::vector<int> component_of_cell(sink_of_cell.size(), -1);
        std::vector<int> component_sink;

        {
            const int n_rows_l = meta.n_rows;
            const int n_cols_l = meta.n_cols;
            std::vector<std::size_t> stack;

            for (std::size_t start = 0; start < sink_of_cell.size(); ++start)
            {
                if (sink_of_cell[start] < 0)          continue;
                if (component_of_cell[start] >= 0)    continue;

                const int this_sink = sink_of_cell[start];
                const int component = static_cast<int>(component_sink.size());
                component_sink.push_back(this_sink);

                component_of_cell[start] = component;
                stack.push_back(start);

                while (!stack.empty())
                {
                    const std::size_t cell = stack.back();
                    stack.pop_back();

                    const int row = static_cast<int>(cell) / n_cols_l;
                    const int col = static_cast<int>(cell) % n_cols_l;

                    for (int d_row = -1; d_row <= 1; ++d_row)
                    for (int d_col = -1; d_col <= 1; ++d_col)
                    {
                        if (d_row == 0 && d_col == 0) continue;

                        const int n_row = row + d_row;
                        const int n_col = col + d_col;

                        if (n_row < 0 || n_row >= n_rows_l) continue;
                        if (n_col < 0 || n_col >= n_cols_l) continue;

                        const std::size_t neighbour =
                            static_cast<std::size_t>(n_row) * n_cols_l + n_col;

                        if (sink_of_cell[neighbour] != this_sink)  continue;
                        if (component_of_cell[neighbour] >= 0)     continue;

                        component_of_cell[neighbour] = component;
                        stack.push_back(neighbour);
                    }
                }
            }
        }

        bool any = false;

        for (std::size_t i = 0; i < sink_of_cell.size(); ++i)
        {
            if (component_of_cell[i] < 0) continue;

            label_data[i] = static_cast<float>(component_of_cell[i] + 1);
            label_mask[i] = 1u;
            any = true;
        }

        if (!any) return;

        const lsdtt_raster_cpu label_raster(
            std::move(label_data), std::move(label_mask), label_meta);

        lsdtt_vector_schema trace_schema;
        trace_schema.field_names = { "label_id", "area_m2" };
        trace_schema.field_types = {
            lsdtt_vector_field_type::int64_type,
            lsdtt_vector_field_type::double_type
        };

        const std::vector<lsdtt_vector_record> traced =
            binary_components_to_polygon_records_with_holes(
                label_raster, nullptr, trace_schema);

        for (const lsdtt_vector_record& record : traced)
        {
            if (record.field_values.empty()) continue;

            // label_id is a COMPONENT id (+1), not a sink index -- see the
            // relabelling above. Map it back to the sink that owns the part.
            const std::int64_t label_id =
                std::get<std::int64_t>(record.field_values[0]);
            const std::size_t component =
                static_cast<std::size_t>(label_id - 1);

            if (component >= component_sink.size())     continue;

            const std::size_t sink_index =
                static_cast<std::size_t>(component_sink[component]);

            if (sink_index >= sink_by_index.size())     continue;
            if (sink_by_index[sink_index] == nullptr)   continue;

            lsdtt_vector_record out;
            out.geometry     = record.geometry;
            out.field_values = attributes_of(*sink_by_index[sink_index]);
            result.records.push_back(std::move(out));
        }
    };

    std::vector<const sink_candidate*> sink_by_index(sinks.size(), nullptr);
    for (std::size_t i = 0; i < sinks.size(); ++i)
    {
        sink_by_index[i] = &sinks[i];
    }

    const std::size_t n_cells =
        static_cast<std::size_t>(meta.n_rows) * meta.n_cols;

    if (opts.method == sink_detection_method::depression_hierarchy)
    {
        hierarchy_state state = build_hierarchy_state(dem, opts);

        // Match the supplied sinks -- which may be an edited or accepted-only
        // subset -- onto depressions by pit cell, and select exactly those.
        //
        // A pit cell is NOT a unique key: a meta-depression's pit is the pit of
        // its lowest child, so a lake and the hollow at the bottom of it share
        // one cell. hierarchy_level separates them. A hand-added point carries
        // no level, so it falls back to the innermost depression at that cell,
        // which is the one a reviewer clicking on a map means.
        std::unordered_map<std::int64_t, std::vector<std::size_t>>
            depressions_at_pit;

        for (std::size_t d = 1; d < state.deps.size(); ++d)
        {
            if (state.pit_cell[d] >= 0)
                depressions_at_pit[state.pit_cell[d]].push_back(d);
        }

        std::vector<int> sink_of_depression(state.deps.size(), -1);
        std::fill(state.selected.begin(), state.selected.end(), 0);

        for (std::size_t sink_index = 0; sink_index < sinks.size(); ++sink_index)
        {
            const sink_candidate& sink = sinks[sink_index];
            const auto found = depressions_at_pit.find(sink.pit_flat);
            if (found == depressions_at_pit.end()) continue;

            std::size_t chosen   = found->second.front();
            bool        by_level = false;

            for (const std::size_t d : found->second)
            {
                if (state.level[d] == sink.hierarchy_level)
                {
                    chosen   = d;
                    by_level = true;
                    break;
                }
                if (!by_level && state.level[d] < state.level[chosen]) chosen = d;
            }

            state.selected[chosen]     = 1;
            sink_of_depression[chosen] = static_cast<int>(sink_index);
        }

        // Each cell is drawn once, into the INNERMOST selected depression that
        // contains it: first walk out of the leaf label until the depression is
        // deep enough to hold the cell, exactly as the hierarchy's own volume
        // accounting does, then walk up to the first selected depression.
        const std::vector<float>&   elevation = dem.get_data();
        const std::vector<uint8_t>& mask      = dem.get_mask();

        std::vector<int> sink_of_cell(n_cells, -1);

        for (std::size_t i = 0; i < n_cells; ++i)
        {
            if (mask[i] == 0u) continue;

            dep_label_t current = state.label[i];
            if (current == k_ocean_label) continue;

            const float z = elevation[i];

            while (current != k_no_parent && current != k_ocean_label)
            {
                if (z <= state.deps[current].out_elev) break;
                if (state.deps[current].ocean_parent)
                {
                    current = k_ocean_label;
                    break;
                }
                current = state.deps[current].parent;
            }

            while (current != k_no_parent && current != k_ocean_label)
            {
                if (state.selected[current] != 0) break;
                current = state.deps[current].parent;
            }


            if (current == k_no_parent || current == k_ocean_label) continue;

            sink_of_cell[i] = sink_of_depression[current];
        }


        // Two selected depressions overlap only when one is an ancestor of the
        // other, so grouping by the number of selected ancestors gives batches
        // in which no two footprints can touch. One tracing pass per batch is
        // then enough, however deeply the depressions nest.
        std::vector<int> batch_of_depression(state.deps.size(), -1);
        int              n_batches = 0;

        for (std::size_t d = 1; d < state.deps.size(); ++d)
        {
            if (state.selected[d] == 0) continue;

            int depth = 0;
            dep_label_t parent = state.deps[d].parent;
            while (parent != k_no_parent && parent != k_ocean_label)
            {
                if (state.selected[parent] != 0) ++depth;
                parent = state.deps[parent].parent;
            }

            batch_of_depression[d] = depth;
            n_batches = std::max(n_batches, depth + 1);
        }

        std::vector<int> batch_of_sink(sink_by_index.size(), -1);

        for (std::size_t d = 1; d < state.deps.size(); ++d)
        {
            if (state.selected[d] == 0)   continue;
            if (sink_of_depression[d] < 0) continue;

            batch_of_sink[static_cast<std::size_t>(sink_of_depression[d])] =
                batch_of_depression[d];
        }

        for (int batch = 0; batch < n_batches; ++batch)
        {
            std::vector<int> batch_cells(n_cells, -1);

            for (std::size_t i = 0; i < n_cells; ++i)
            {
                if (sink_of_cell[i] < 0) continue;

                if (batch_of_sink[static_cast<std::size_t>(sink_of_cell[i])] == batch)
                {
                    batch_cells[i] = sink_of_cell[i];
                }
            }

            trace_footprints(batch_cells, sink_by_index);
        }
    }
    else
    {
        const fill_difference_state state = build_fill_difference_state(dem, opts);

        // Fill-difference components are disjoint by construction, so every
        // footprint can be traced in one pass. Sinks are matched to components
        // by the label under the pit cell.
        const std::vector<float>&   label_data = state.labelled.get_data();
        const std::vector<uint8_t>& label_mask = state.labelled.get_mask();

        std::unordered_map<int, int> sink_of_component;

        for (std::size_t sink_index = 0; sink_index < sinks.size(); ++sink_index)
        {
            const sink_candidate& sink = sinks[sink_index];
            if (sink.pit_flat < 0 ||
                static_cast<std::size_t>(sink.pit_flat) >= n_cells) continue;

            const std::size_t index = static_cast<std::size_t>(sink.pit_flat);
            if (label_mask[index] == 0u) continue;

            const int component_label = static_cast<int>(label_data[index]);
            if (component_label <= 0) continue;

            sink_of_component[component_label] = static_cast<int>(sink_index);
        }

        std::vector<int> sink_of_cell(n_cells, -1);

        for (std::size_t i = 0; i < n_cells; ++i)
        {
            if (label_mask[i] == 0u) continue;

            const int component_label = static_cast<int>(label_data[i]);
            if (component_label <= 0) continue;

            const auto found = sink_of_component.find(component_label);
            if (found == sink_of_component.end()) continue;

            sink_of_cell[i] = found->second;
        }

        trace_footprints(sink_of_cell, sink_by_index);
    }

    std::sort(result.records.begin(), result.records.end(),
              [](const lsdtt_vector_record& a, const lsdtt_vector_record& b)
              {
                  return std::get<std::int64_t>(a.field_values[0]) <
                         std::get<std::int64_t>(b.field_values[0]);
              });

    return result;
}

// =============================================================================
// Sink-terminated routing support
// =============================================================================

namespace {

void validate_routing_sinks(const std::vector<sink_candidate>& sinks,
                            std::size_t                         n_cells,
                            int                                 n_cols,
                            const char*                         caller)
{
    constexpr std::size_t max_exact_float_label =
        static_cast<std::size_t>(1u) << 24u;
    if (!sinks.empty() && sinks.size() - 1u > max_exact_float_label)
    {
        throw std::invalid_argument(
            std::string(caller) +
            ": routing_sink_id exceeds the exact integer range of a float raster");
    }

    std::vector<std::uint8_t> routing_ids(sinks.size(), 0u);
    std::unordered_set<std::int64_t> pit_flats;

    for (const sink_candidate& sink : sinks)
    {
        if (!sink.accepted)
        {
            throw std::invalid_argument(
                std::string(caller) +
                ": input contains a rejected sink; call accepted_sinks first");
        }
        if (sink.routing_sink_id < 0 ||
            static_cast<std::size_t>(sink.routing_sink_id) >= sinks.size())
        {
            throw std::invalid_argument(
                std::string(caller) +
                ": routing_sink_id values must be contiguous from zero");
        }

        const std::size_t routing_id =
            static_cast<std::size_t>(sink.routing_sink_id);
        if (routing_ids[routing_id] != 0u)
        {
            throw std::invalid_argument(
                std::string(caller) + ": duplicate routing_sink_id " +
                std::to_string(sink.routing_sink_id));
        }
        routing_ids[routing_id] = 1u;

        if (sink.pit_flat < 0 ||
            static_cast<std::size_t>(sink.pit_flat) >= n_cells)
        {
            throw std::invalid_argument(
                std::string(caller) + ": pit_flat is outside the raster for sink_id " +
                std::to_string(sink.sink_id));
        }
        if (!pit_flats.insert(sink.pit_flat).second)
        {
            throw std::invalid_argument(
                std::string(caller) + ": duplicate pit_flat " +
                std::to_string(sink.pit_flat));
        }

        if (sink.pit_row >= 0 || sink.pit_col >= 0)
        {
            const int expected_row = static_cast<int>(sink.pit_flat / n_cols);
            const int expected_col = static_cast<int>(sink.pit_flat % n_cols);
            if (sink.pit_row != expected_row || sink.pit_col != expected_col)
            {
                throw std::invalid_argument(
                    std::string(caller) +
                    ": pit_row/pit_col do not match pit_flat for sink_id " +
                    std::to_string(sink.sink_id));
            }
        }
    }
}

struct routing_seed
{
    std::int64_t flat = -1;
    double       final_step_m = 0.0;
};

} // namespace

lsdtt_raster_cpu
punch_sink_holes(const lsdtt_raster_cpu&            dem,
                 const std::vector<sink_candidate>& sinks)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n_cells = dem.full_grid_size();

    validate_routing_sinks(sinks, n_cells, meta.n_cols, "punch_sink_holes");

    std::vector<float> data = dem.get_data();
    std::vector<std::uint8_t> mask = dem.get_mask();
    const float ndv = static_cast<float>(meta.no_data_value);

    for (const sink_candidate& sink : sinks)
    {
        const std::size_t index = static_cast<std::size_t>(sink.pit_flat);
        if (mask[index] == 0u || meta.is_nodata(data[index]))
        {
            throw std::invalid_argument(
                "punch_sink_holes: pit_flat lands on nodata or a masked cell for sink_id " +
                std::to_string(sink.sink_id));
        }

        data[index] = ndv;
        mask[index] = 0u;
    }

    return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

sink_tagging_result
tag_upslope_of_sinks(const lsdtt_flow_base&            flow,
                     const std::vector<sink_candidate>& sinks,
                     sink_mode                          mode,
                     const sink_tagging_options&        opts)
{
    if (flow.empty())
    {
        throw std::invalid_argument("tag_upslope_of_sinks: flow object is empty");
    }
    if (!flow.has_s_vector())
    {
        throw std::invalid_argument(
            "tag_upslope_of_sinks: flow object has no s_vector");
    }
    if (std::isnan(opts.max_upslope_distance_m) ||
        opts.max_upslope_distance_m < 0.0f)
    {
        throw std::invalid_argument(
            "tag_upslope_of_sinks: max_upslope_distance_m must be non-negative");
    }

    const lsdtt_raster_meta& meta = flow.meta();
    const std::size_t n_cells = flow.full_grid_size();
    if (n_cells != static_cast<std::size_t>(meta.n_pixels()))
    {
        throw std::invalid_argument(
            "tag_upslope_of_sinks: flow grid size does not match its metadata");
    }

    validate_routing_sinks(
        sinks, n_cells, meta.n_cols, "tag_upslope_of_sinks");

    std::vector<std::vector<routing_seed>> seeds(sinks.size());
    std::unordered_map<std::int64_t, int> selected_outlets;

    if (mode == sink_mode::forced_outlet)
    {
        for (const sink_candidate& sink : sinks)
        {
            if (!flow.is_valid_cell(sink.pit_flat) ||
                !flow.is_outlet(sink.pit_flat))
            {
                throw std::invalid_argument(
                    "tag_upslope_of_sinks: forced-outlet pit is not a valid outlet for sink_id " +
                    std::to_string(sink.sink_id));
            }

            seeds[static_cast<std::size_t>(sink.routing_sink_id)].push_back(
                { sink.pit_flat, 0.0 });
            selected_outlets.emplace(sink.pit_flat, sink.routing_sink_id);
        }
    }
    else
    {
        static constexpr int d_row[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
        static constexpr int d_col[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };

        for (const sink_candidate& sink : sinks)
        {
            if (flow.is_valid_cell(sink.pit_flat))
            {
                throw std::invalid_argument(
                    "tag_upslope_of_sinks: nodata-hole pit is still valid for sink_id " +
                    std::to_string(sink.sink_id));
            }

            const int pit_row = flow.row_from_flat(sink.pit_flat);
            const int pit_col = flow.col_from_flat(sink.pit_flat);
            std::vector<routing_seed>& sink_seeds =
                seeds[static_cast<std::size_t>(sink.routing_sink_id)];

            for (int direction = 0; direction < 8; ++direction)
            {
                const int row = pit_row + d_row[direction];
                const int col = pit_col + d_col[direction];
                if (!flow.in_bounds(row, col) || !flow.is_valid_cell(row, col))
                    continue;

                const std::int64_t flat = flow.flat_index_unchecked(row, col);
                if (!flow.is_outlet(flat)) continue;

                // The abstract flow API records only that this cell is an
                // outlet, not which adjacent exit caused it. An edge outlet or
                // a cell touching another void cannot be attributed uniquely
                // to this hole, so refuse it rather than importing an external
                // basin into the sink label.
                if (row == 0 || row == meta.n_rows - 1 ||
                    col == 0 || col == meta.n_cols - 1)
                {
                    throw std::invalid_argument(
                        "tag_upslope_of_sinks: nodata-hole outlet attribution is ambiguous at the raster edge");
                }
                for (int around = 0; around < 8; ++around)
                {
                    const int other_row = row + d_row[around];
                    const int other_col = col + d_col[around];
                    if (!flow.in_bounds(other_row, other_col)) continue;

                    const std::int64_t other =
                        flow.flat_index_unchecked(other_row, other_col);
                    if (other != sink.pit_flat && !flow.is_valid_cell(other))
                    {
                        throw std::invalid_argument(
                            "tag_upslope_of_sinks: nodata-hole outlet attribution is ambiguous beside another nodata cell");
                    }
                }

                const auto existing = selected_outlets.find(flat);
                if (existing != selected_outlets.end() &&
                    existing->second != sink.routing_sink_id)
                {
                    throw std::invalid_argument(
                        "tag_upslope_of_sinks: an outlet neighbour is shared by multiple nodata holes");
                }
                if (existing != selected_outlets.end()) continue;

                selected_outlets.emplace(flat, sink.routing_sink_id);
                const bool diagonal = d_row[direction] != 0 && d_col[direction] != 0;
                sink_seeds.push_back(
                    { flat, meta.cell_size * (diagonal ? std::sqrt(2.0) : 1.0) });
            }

            if (sink_seeds.empty())
            {
                throw std::invalid_argument(
                    "tag_upslope_of_sinks: nodata hole has no adjacent outlet for sink_id " +
                    std::to_string(sink.sink_id));
            }
        }
    }

    const lsdtt_raster_result distance = flow.get_distance_from_outlet_raster();
    if (distance.data.size() != n_cells || distance.mask.size() != n_cells)
    {
        throw std::runtime_error(
            "tag_upslope_of_sinks: distance raster size does not match flow grid");
    }

    lsdtt_raster_meta label_meta = meta;
    label_meta.no_data_value = -1.0;
    const float ndv = -1.0f;
    std::vector<float> tag_data(n_cells, ndv);
    std::vector<std::uint8_t> tag_mask(n_cells, 0u);
    std::vector<int> owner(n_cells, -1);

    sink_tagging_result result;
    result.area_m2.assign(sinks.size(), 0.0);
    result.n_pixels.assign(sinks.size(), 0);

    std::int64_t tagged_pixels = 0;
    std::int64_t untagged_internal_pixels = 0;

    for (std::size_t routing_id = 0; routing_id < seeds.size(); ++routing_id)
    {
        for (const routing_seed& seed : seeds[routing_id])
        {
            const lsdtt_flat_index_slice upslope = flow.upslope_nodes(seed.flat);
            for (const std::int64_t flat : upslope)
            {
                const std::size_t index = static_cast<std::size_t>(flat);
                const double distance_to_sink =
                    static_cast<double>(distance.data[index]) + seed.final_step_m;
                if (distance_to_sink >=
                    static_cast<double>(opts.max_upslope_distance_m))
                {
                    ++untagged_internal_pixels;
                    continue;
                }

                if (owner[index] >= 0 &&
                    owner[index] != static_cast<int>(routing_id))
                {
                    throw std::runtime_error(
                        "tag_upslope_of_sinks: sink upslope sets overlap");
                }
                if (owner[index] == static_cast<int>(routing_id)) continue;

                owner[index] = static_cast<int>(routing_id);
                tag_data[index] = static_cast<float>(routing_id);
                tag_mask[index] = 1u;
                ++result.n_pixels[routing_id];
                ++tagged_pixels;
            }
        }
    }

    std::int64_t non_sink_outlet_pixels = 0;
    for (std::size_t index = 0; index < n_cells; ++index)
    {
        const std::int64_t flat = static_cast<std::int64_t>(index);
        if (!flow.is_valid_cell(flat) || !flow.is_outlet(flat)) continue;
        if (selected_outlets.find(flat) != selected_outlets.end()) continue;
        non_sink_outlet_pixels +=
            static_cast<std::int64_t>(flow.s_vector_size(flat));
    }

    if (tagged_pixels + untagged_internal_pixels + non_sink_outlet_pixels !=
        flow.n_valid_cells())
    {
        throw std::runtime_error(
            "tag_upslope_of_sinks: area conservation failed");
    }

    const double pixel_area = meta.pixel_area();
    for (std::size_t i = 0; i < result.n_pixels.size(); ++i)
    {
        result.area_m2[i] =
            static_cast<double>(result.n_pixels[i]) * pixel_area;
    }
    result.non_sink_outlet_area_m2 =
        static_cast<double>(non_sink_outlet_pixels) * pixel_area;
    result.untagged_internal_area_m2 =
        static_cast<double>(untagged_internal_pixels) * pixel_area;
    result.sink_id_raster =
        lsdtt_raster_cpu(std::move(tag_data), std::move(tag_mask), label_meta);

    return result;
}


std::vector<std::int64_t> sink_channel_sources(
    const lsdtt_flow_base&     flow,
    const sink_tagging_result& tagged,
    const std::uint64_t        threshold_contributing_pixels
)
{
    const std::size_t n_cells = flow.full_grid_size();

    const std::vector<std::uint8_t>& sink_mask = tagged.sink_id_raster.get_mask();
    const std::vector<float>&        sink_id   = tagged.sink_id_raster.get_data();

    std::vector<std::uint8_t> qualifies(n_cells, 0u);
    std::vector<std::uint8_t> has_qualifying_donor(n_cells, 0u);

    for (std::size_t index = 0; index < n_cells; ++index)
    {
        const std::int64_t flat = static_cast<std::int64_t>(index);
        if (sink_mask[index] != 0u &&
            flow.s_vector_size(flat) >= threshold_contributing_pixels)
        {
            qualifies[index] = 1u;
        }
    }

    for (std::size_t index = 0; index < n_cells; ++index)
    {
        if (qualifies[index] == 0u) continue;

        const std::int64_t receiver = flow.receiver_of(static_cast<std::int64_t>(index));
        if (receiver < 0) continue;

        const std::size_t r = static_cast<std::size_t>(receiver);

        // Same sink only: see the header.
        if (qualifies[r] != 0u && sink_id[index] == sink_id[r])
        {
            has_qualifying_donor[r] = 1u;
        }
    }

    std::vector<std::int64_t> sources;
    for (std::size_t index = 0; index < n_cells; ++index)
    {
        if (qualifies[index] != 0u && has_qualifying_donor[index] == 0u)
        {
            sources.push_back(static_cast<std::int64_t>(index));
        }
    }

    return sources;
}

} // namespace lsdtt
