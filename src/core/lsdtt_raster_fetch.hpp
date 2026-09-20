// lsdtt_raster_fetch.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Native DEM fetching: turn a WGS84 bounding box into an analysis-ready
// projected raster, using only GDAL. No Python, no extra dependency.
//
// Testing and benchmarking notes: docs/developer_notes/testing_and_benchmarking.md
//
// Source: Copernicus GLO-30 on the AWS Open Data registry. Public,
// keyless, stable URL scheme, stored as Cloud-Optimised GeoTIFFs, so only the
// requested window travels over the wire (HTTP range requests) rather than a
// whole 1-degree tile.
//
// *** CREDENTIALS ***
// This header intentionally contains NO credential fields. The default source
// is keyless, so most users never supply one. When keyed providers arrive (docs/developer_notes/testing_and_benchmarking.md
// tiers 1-3) their secrets must come from lsdtt_secret (environment variable,
// key file, or TTY prompt) and must NEVER become:
//   - a member of lsdtt_fetch_options,
//   - a CLI parameter, or
//   - anything reachable by lsdtt_param_parser::print_parameters().
// lsdtt3 CLIs are parameter-file driven and echo every parameter at the end of
// a run, so a credential held as a parameter would be written to logs and
// committed alongside users' analyses. Keeping secrets out of this struct is
// what makes that impossible rather than merely discouraged.
//
// Pipeline:
//   1. resolve which source tiles the bounding box covers (mosaic if several)
//   2. windowed read over /vsicurl/  -> in-memory dataset
//   3. warp to a metric CRS at a fixed grid spacing, explicit nodata
//   4. optional exact centre crop to a square pixel count
//
// Reprojection happens here, at fetch time, rather than in whatever CLI later
// consumes the raster. When two implementations are compared, both must read a
// byte-identical grid or the differences are uninterpretable.
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

namespace lsdtt {

// -----------------------------------------------------------------------------
// lsdtt_dem_source
//
// Sources are declared here so callers
// and CLI help can be written against the final vocabulary, and each returns a
// clear "not implemented" result rather than silently doing something else.
// -----------------------------------------------------------------------------
// Copernicus Data Space (CDSE) is deliberately NOT offered (Simon, 2026-07-20):
// its only distinctive product is EU-DEM at 10 m, which a general public account
// does not unlock, so it would carry the only username/password credential risk
// in the tool for coverage already available keyless from cop30_aws.
enum class lsdtt_dem_source
{
    cop30_aws,            // Copernicus GLO-30, AWS Open Data. No key.
    opentopography,       // OpenTopography REST. Needs an API key.
    planetary_computer    // Keyless, but asset URLs need SAS signing.
};

// Parse a user-facing source name ("cop30_aws", ...). Returns false if unknown.
bool parse_dem_source(const std::string& name, lsdtt_dem_source& out);

// Name of a source, for messages and help text.
std::string dem_source_name(lsdtt_dem_source source);


// -----------------------------------------------------------------------------
// lsdtt_ot_dataset
//
// One entry of the OpenTopography global raster catalogue: the `demtype` values
// accepted by https://portal.opentopography.org/API/globaldem.
//
// The catalogue is carried in the binary rather than being discovered at run
// time so that `-h ot_dem_type` works offline, and so that a mistyped demtype
// fails immediately with the list of valid values instead of after a network
// round trip that returns an opaque 400.
//
// `citation` is the citation OpenTopography publishes for the dataset. It is
// printed after every successful OpenTopography download: these are DOI-minted
// distributions whose providers require attribution, and a citation nobody sees
// until they go looking for it is a citation that does not get used.
//
// Source of the fields: the globaldem OpenAPI description and the
// otCatalog API (accessed 2026-07-22).
// -----------------------------------------------------------------------------
struct lsdtt_ot_dataset
{
    std::string demtype;      // value passed as demtype=, e.g. "COP30"
    std::string name;         // human-readable dataset name
    std::string resolution;   // nominal grid spacing, e.g. "1 arc-second (~30 m)"
    std::string gdal_type;    // source band type in GTiff responses, e.g. "Int16"
    std::string coverage;     // "global", "South America", ...
    std::string description;  // one line: what it is and when to reach for it
    std::string citation;     // full citation including DOI
    double max_area_km2 = 0;  // per-request area limit imposed by the API
};

// The full catalogue, in the order it should be shown to a user.
const std::vector<lsdtt_ot_dataset>& ot_dataset_catalogue();

// Look up one entry by demtype. Returns nullptr if the demtype is not known.
// The comparison is case sensitive: the API is.
const lsdtt_ot_dataset* ot_dataset_info(const std::string& demtype);

// Comma-separated list of valid demtype values, for error messages.
std::string ot_dataset_demtype_list();


// -----------------------------------------------------------------------------
// lsdtt_fetch_options
//
// NOTE: deliberately contains no credential fields. See the header comment.
// -----------------------------------------------------------------------------
struct lsdtt_fetch_options
{
    // Bounding box in WGS84 decimal degrees.
    double west  = 0.0;
    double south = 0.0;
    double east  = 0.0;
    double north = 0.0;

    lsdtt_dem_source source = lsdtt_dem_source::cop30_aws;

    // Target CRS. 0 = derive the UTM zone from the bounding box centre.
    int target_epsg = 0;

    // Target grid spacing in projection units (metres for UTM).
    double grid_spacing = 30.0;

    // Any GDAL warp resampling string ("cubic", "bilinear", "near", ...).
    std::string resampling = "cubic";

    // Written into the output. Matches lsdtt_raster_meta's default so the rest
    // of lsdtt3 needs no special handling. An explicit nodata value also stops
    // cells outside the source footprint defaulting to 0, which would read as
    // sea level and become a large spurious depression.
    double no_data_value = -9999.0;

    // > 0 crops an exact square of this many pixels from the centre after
    // warping, so output dimensions are deterministic rather than a function of
    // GDAL's warp extent arithmetic. 0 keeps the full warped extent.
    int size_px = 0;

    // Opt-in caching. Empty = no caching; the fetch always runs and writes only
    // to out_path. When set, a previously fetched raster with the same cache key
    // is reused instead of re-downloading.
    std::string cache_dir;

    // OpenTopography dataset identifier, e.g. "COP30", "SRTMGL1", "NASADEM".
    // Ignored by keyless sources.
    std::string ot_dem_type = "COP30";

    // Planetary Computer STAC collection and asset key.
    std::string pc_collection = "cop-dem-glo-30";
    std::string pc_asset      = "data";

    // *** NOT A CREDENTIAL FIELD ***
    // This is resolved by lsdtt_secret at the call site and passed in, rather
    // than being read from a parameter. It is here only so fetch_dem() can
    // build the request URL and mask the key out of any message it produces.
    // It must never be populated from a parameter file, from argv, or from
    // anything reachable by print_parameters(). See lsdtt_secret.hpp.
    std::string api_key;

    // Seconds. Bounds how long an unreachable host can stall a caller.
    int http_timeout_seconds = 60;
    int http_connect_timeout_seconds = 15;
    int http_max_retry = 3;
};


// -----------------------------------------------------------------------------
// lsdtt_fetch_result
//
// `skipped` is separate from `ok` on purpose. A machine with no network, or a
// GDAL built without curl, is an *unavailable dependency*, not a failure -- it
// must let tests skip (CTest exit 77) rather than go red. A malformed bounding
// box, by contrast, is a genuine error: ok = false, skipped = false.
// -----------------------------------------------------------------------------
struct lsdtt_fetch_result
{
    bool ok         = false;
    bool skipped    = false;
    bool from_cache = false;

    std::string path;      // written raster (empty if not written)
    std::string message;   // human-readable outcome; safe to log

    int    epsg         = 0;
    int    n_cols       = 0;
    int    n_rows       = 0;
    double grid_spacing = 0.0;

    std::vector<std::string> tiles;   // source tiles used
};


// -----------------------------------------------------------------------------
// fetch_dem
//
// Fetch, reproject and write one DEM. Never throws for expected conditions
// (offline, no curl, bad bounding box) -- inspect the result instead.
// -----------------------------------------------------------------------------
lsdtt_fetch_result fetch_dem(const lsdtt_fetch_options& options,
                             const std::string&         out_path);


// -----------------------------------------------------------------------------
// Helpers, exposed for testing and for callers that need them directly.
// -----------------------------------------------------------------------------

// EPSG code of the UTM zone containing a point. Northern hemisphere -> 326xx,
// southern -> 327xx.
int utm_epsg_for_point(double lon, double lat);

// Copernicus GLO-30 tile names covering a bounding box. Tiles are 1 degree
// squares named by their south-west corner, e.g. "N46_00_E010_00".
std::vector<std::string> cop30_tiles_for_bbox(double west, double south,
                                              double east, double north);

// Full /vsicurl/ URL for one Copernicus GLO-30 tile.
std::string cop30_tile_url(const std::string& tile_name);

// True if GDAL can make HTTP requests at all (built with curl).
bool http_available();


// -----------------------------------------------------------------------------
// lsdtt_bounds_from_vector_result
// -----------------------------------------------------------------------------
struct lsdtt_bounds_from_vector_result
{
    bool ok = false;

    // WGS84 decimal degrees, ready to drop into lsdtt_fetch_options.
    double west  = 0.0;
    double south = 0.0;
    double east  = 0.0;
    double north = 0.0;

    std::string layer_name;      // layer actually used
    std::string geometry_kind;   // "polygon" or "points"
    int         n_features = 0;
    std::string source_crs;      // description of the input CRS, for reporting

    std::string warning;         // non-fatal advisory, safe to log
    std::string message;         // why it failed, or a summary on success
};


// -----------------------------------------------------------------------------
// bounds_from_vector_file
//
// Derive a WGS84 bounding box from a vector file (FlatGeobuf, GeoPackage, or
// anything else OGR reads). Accepts either:
//
//   - one or more polygons  -> the envelope of all of them
//   - exactly two points    -> the box they span
//
// Two points define a rectangle by its diagonal. Note that it does not matter
// which diagonal: lower-left/upper-right and lower-right/upper-left produce the
// SAME bounding box, because it is just min/max of the two x and y values. So
// the orientation is detected for free, with no ambiguity and nothing to get
// wrong.
//
// The geometry is reprojected to WGS84 when the layer has a different CRS, so a
// polygon drawn in UTM in a GIS works directly. Reprojection is applied to the
// geometry before taking the envelope, not to the envelope corners, which would
// be wrong for anything but an axis-aligned transform.
//
// layer_name may be empty to take the first (or only) layer.
// Never throws: inspect .ok and .message.
// -----------------------------------------------------------------------------
lsdtt_bounds_from_vector_result bounds_from_vector_file(
    const std::string& path,
    const std::string& layer_name = ""
);

} // namespace lsdtt
