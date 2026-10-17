// lsdtt_macrostrat_fetch.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Geological map units from Macrostrat (macrostrat.org), for a bounding box.
//
// Macrostrat's map API answers per point, so an area comes from its vector
// tiles (https://tiles.macrostrat.org/carto/{z}/{x}/{y}.mvt, layer "units",
// EPSG:3857), read through GDAL's MVT driver. Each unit carries map_id,
// source_id, name, a free-text "lith", ages and its source map. The API then
// adds, per map_id, the linked lithology vocabulary ids ("liths"), resolved
// through /defs/lithologies, and the source references for attribution.
//
// The output is a polygon layer reprojected to the UTM zone of the bounding
// box centre (Web Mercator does not preserve distance, so lsdtt refuses it),
// with lith_primary: the first lithology of each unit (see
// macrostrat_primary_lithology). Data are CC BY 4.0: cite Macrostrat and each
// source map.
//
// GDAL-free header.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "lsdtt_http_fetch.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// Slippy-map tiles (EPSG:3857, XYZ, y counted from the north)
// ---------------------------------------------------------------------------
struct lsdtt_slippy_tile
{
    int z = 0, x = 0, y = 0;
};

// The tiles at zoom z covering the WGS84 box (degrees), row by row from the
// north-west. Throws std::invalid_argument for an empty or out-of-range box
// (latitudes are limited to +-85.0511, the Web Mercator limit).
[[nodiscard]] std::vector<lsdtt_slippy_tile> slippy_tiles_for_bbox(
    double west, double south, double east, double north, int z);

// A tile's bounds in EPSG:3857 metres.
void slippy_tile_bounds_3857(const lsdtt_slippy_tile& t, double& x_min, double& y_min,
                             double& x_max, double& y_max);

// ---------------------------------------------------------------------------
// Macrostrat vocabulary and parsers (pure; throw std::runtime_error on
// malformed JSON or a response without success/data)
// ---------------------------------------------------------------------------
struct lsdtt_macrostrat_lithology
{
    int id = 0;
    std::string name, type, group, lith_class;
};

struct lsdtt_macrostrat_source
{
    int id = 0;
    std::string name, authors, year, url, ref_title, ref_source, doi, scale;
};

// /defs/lithologies?all
std::map<int, lsdtt_macrostrat_lithology> parse_macrostrat_lithologies(const std::string& json);
// /geologic_units/map?map_id=a,b,...: map_id -> liths (vocabulary ids, API order)
std::map<std::int64_t, std::vector<int>> parse_macrostrat_map_unit_liths(const std::string& json);
// /defs/sources?source_id=a,b,...
std::vector<lsdtt_macrostrat_source> parse_macrostrat_sources(const std::string& json);

// The first lithology of a map unit (decided, LD6), from, in order:
//   "major": the first term of a "Major:{a,b,...}" lith text;
//   "liths": the name of the first linked vocabulary id;
//   "text":  the free-text lith (trimmed, lower case, whitespace collapsed);
//   "none":  nothing (empty value).
// The result is lower case.
struct lsdtt_macrostrat_primary
{
    std::string value;
    std::string source;
};
[[nodiscard]] lsdtt_macrostrat_primary macrostrat_primary_lithology(
    const std::string& lith_text, const std::vector<int>& lith_ids,
    const std::map<int, lsdtt_macrostrat_lithology>& vocabulary);

// WGS 84 / UTM EPSG code for a longitude and latitude (326zz north, 327zz south).
[[nodiscard]] int utm_epsg_for(double lon, double lat);

// ---------------------------------------------------------------------------
// The fetch
// ---------------------------------------------------------------------------
struct lsdtt_macrostrat_fetch_options
{
    double west = 0.0, south = 0.0, east = 0.0, north = 0.0;   // WGS84 degrees
    int    zoom = 10;
    int    max_tiles = 256;
    int    target_epsg = 0;                // 0 = UTM of the box centre
    bool   lookup_map_units = true;        // API: liths per map_id, vocabulary, sources
    int    request_delay_milliseconds = 200;
    lsdtt_http_request_options http;
    std::string tile_url_template = "https://tiles.macrostrat.org/carto/{z}/{x}/{y}.mvt";
    std::string api_base = "https://macrostrat.org/api/v2";
    std::string tile_dir;                  // non-empty: read <tile_dir>/{z}/{x}/{y}.mvt|.pbf, no HTTP
    std::string out_path;                  // vector file to write
    std::string out_driver = "FlatGeobuf"; // OGR driver name
};

struct lsdtt_macrostrat_vocabulary_row   // one per distinct lith_primary
{
    std::string lith_primary, lith_primary_source, lith_type, lith_class;
    std::size_t n_features = 0;
    double area_km2 = 0.0;
};

struct lsdtt_macrostrat_fetch_result
{
    int epsg = 0;
    std::size_t n_tiles = 0, n_tiles_missing = 0, n_features = 0, n_map_units = 0;
    std::vector<lsdtt_macrostrat_vocabulary_row> vocabulary;   // by area, largest first
    std::vector<lsdtt_macrostrat_source> sources;
    std::vector<std::string> notes, warnings;
};

// Fetches the tiles, merges their units, looks up the map units (optional),
// writes opts.out_path (layer "units") and returns the summary. Throws
// std::runtime_error on a failure that leaves no usable output (no tile read,
// too many tiles, an unwritable output); a failed API lookup is a warning and
// lith_primary then falls back to the tile's lith text.
lsdtt_macrostrat_fetch_result fetch_macrostrat_units(const lsdtt_macrostrat_fetch_options& opts);

// The attribution text for a fetch: Macrostrat and every source map (CC BY 4.0).
[[nodiscard]] std::string macrostrat_attribution_text(const lsdtt_macrostrat_fetch_result& r);

} // namespace lsdtt
