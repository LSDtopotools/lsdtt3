// lsdtt_river_line_fetch.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// River centrelines for an area, READ from a source (lsdtt-fetch-river-lines).
//
// Sources (river_line_source):
//   local_file         any OGR-readable line file the user has
//   gb_os_open_rivers  OS Open Rivers (GB, OGL v3). The national GeoPackage zip
//                      is downloaded once from the keyless OS Downloads API into
//                      the cache, its md5 checked against the API's, and read in
//                      place through /vsizip/ with a spatial filter.
//
// This core module resolves the area and output CRS, fetches / caches the
// source, and reads each line feature through a field map (name, id, feature
// class, direction, start / end node ids, culverts) into plain
// river_line_feature records in the output CRS. The area is a WGS84 box or a
// raster (bounds_raster_path): with a raster, the output CRS defaults to the
// raster's and the clip box is the raster's extent plus the buffer.
//
// Preparation (noding, direction vote, topology) and writing are computation
// and live in src/cpu/lsdtt_river_line_prepare_cpu, which also provides
// fetch_river_lines(): read -> prepare -> write.
//
// GDAL-free header. Never prints: notes and warnings are returned.
// -----------------------------------------------------------------------------

#pragma once

#include <map>
#include <string>
#include <vector>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_http_fetch.hpp"
#include "lsdtt_types.hpp"

namespace lsdtt
{

// One feature as read from a source, already in the output CRS (metres).
struct river_line_feature
{
    // One entry per part; a single-part feature has one.
    std::vector<std::vector<lsdtt_xy>> parts;

    std::string source_feature_id;
    std::string name;
    std::string feature_class;      // river, stream, lake_connector, canal, drain, tidal, other

    // Source node ids at the digitised start and end (Open Rivers start_node /
    // end_node). Empty when the source has none. Only used for multipart-free
    // features: for a multipart feature they are ignored.
    std::string start_node_id;
    std::string end_node_id;

    // +1: digitised downstream; -1: digitised upstream (reverse it);
    // 0: unknown (kept as digitised).
    int direction = 0;
    // attribute, digitised, osm_convention or unknown.
    std::string direction_source = "unknown";

    int culverted = -1;             // 1 yes, 0 no, -1 unknown
    int strahler_order = -1;        // -1 when the source gives none
};

// How to read one source's attributes. An empty field name means "not present".
struct lsdtt_river_line_field_map
{
    std::string id_field;
    std::string name_field;                 // empty: first of watercourse_name, name, gnis_name
    std::string feature_class_field;        // empty: every line is "river"
    std::string direction_field;
    std::map<std::string, int> direction_values;   // value (spaces as _) -> +1 / -1
    bool digitised_downstream = false;      // when there is no direction field
    std::string start_node_field;
    std::string end_node_field;
    std::string culvert_field;              // OSM "tunnel": culvert / yes -> 1
};

struct lsdtt_river_line_fetch_options
{
    std::string source = "local_file";      // local_file | gb_os_open_rivers

    // Area: a WGS84 box, or a raster whose extent and CRS are used.
    double west = 0.0, south = 0.0, east = 0.0, north = 0.0;
    std::string bounds_raster_path;

    int    output_epsg = 0;                 // 0: the raster's CRS, else UTM of the box centre
    double clip_buffer_metres = -1.0;       // < 0: 5% of the shorter side, at least 200 m

    // local_file
    std::string river_lines_path;
    std::string river_lines_layer;
    lsdtt_river_line_field_map fields;

    std::vector<std::string> keep_feature_classes = {"river", "stream", "lake_connector"};
    double node_snap_tolerance_metres = -1.0;   // < 0: 1 m
    double min_line_length_metres = 0.0;
    bool   repair_direction = true;         // direction_repair_algorithm_choice: network_outlet

    // gb_os_open_rivers
    std::string cache_dir;                  // empty: default_river_cache_dir()
    std::string os_downloads_api_url =
        "https://api.os.uk/downloads/v1/products/OpenRivers/downloads";
    lsdtt_http_request_options http;

    // outputs (used by fetch_river_lines in lsdtt_river_line_prepare_cpu)
    std::string out_lines_path;
    std::string out_nodes_path;
    vector_output_format out_format = vector_output_format::flatgeobuf;
    bool spatial_index = true;
};

struct lsdtt_river_line_class_count
{
    int n_read = 0;
    int n_kept = 0;
};

struct lsdtt_river_line_read_result
{
    std::vector<river_line_feature> features;   // kept classes, output CRS
    std::string crs_wkt;                      // of the output CRS
    int epsg = 0;                            // 0 if the output CRS has no EPSG code
    std::string crs_description;
    double area_x_min = 0, area_y_min = 0, area_x_max = 0, area_y_max = 0;   // output CRS
    double clip_x_min = 0, clip_y_min = 0, clip_x_max = 0, clip_y_max = 0;   // area + buffer
    double clip_buffer_metres = 0.0;
    double node_snap_tolerance_metres = 0.0;
    std::string read_path;                   // the dataset actually opened
    int n_features_read = 0;
    int n_non_line_features = 0;
    std::map<std::string, lsdtt_river_line_class_count> classes;   // by feature_class
    std::map<std::string, int> direction_sources;                  // kept lines' parts
    std::string attribution;
    std::vector<std::string> notes;
    std::vector<std::string> warnings;
};

// Resolve the area and output CRS, fetch / cache the source, and read the line
// features of the kept classes. Throws std::runtime_error /
// std::invalid_argument on any failure that leaves nothing to prepare.
lsdtt_river_line_read_result read_river_lines(const lsdtt_river_line_fetch_options& opts);

// ---------------------------------------------------------------------------
// Pure helpers (exposed for tests)
// ---------------------------------------------------------------------------

// The common feature class for a source value (case-insensitive):
// river, inlandriver -> river; stream, brook -> stream; lake -> lake_connector;
// canal; drain, ditch -> drain; tidalriver, tidal -> tidal; anything else -> other.
std::string river_feature_class_from_value(const std::string& value);

// "in_direction=1,in_opposite_direction=-1" -> map. Spaces in a value are
// written as underscores (a parameter value cannot contain a space); the
// attribute values are compared the same way. Throws on a malformed pair.
std::map<std::string, int> parse_direction_field_mapping(const std::string& text);

// The OS Downloads API product list: the GeoPackage entry.
struct lsdtt_os_download_entry
{
    std::string url, file_name, md5;
    long long size = 0;
};
lsdtt_os_download_entry parse_os_downloads_geopackage(const std::string& json);

// Lower-case hex MD5 of a byte string (RFC 1321).
std::string md5_hex(const std::string& bytes);

// $XDG_CACHE_HOME/lsdtt3/river_lines, else $HOME/.cache/lsdtt3/river_lines,
// else %LOCALAPPDATA%/lsdtt3/river_lines, else ./lsdtt3_river_lines_cache.
std::string default_river_cache_dir();

// The field map used for OS Open Rivers.
lsdtt_river_line_field_map os_open_rivers_field_map();

} // namespace lsdtt
