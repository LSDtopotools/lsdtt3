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

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_oblique_frame.hpp"
#include "lsdtt_rotation_frame_file.hpp"

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
    planetary_computer,   // Keyless, but asset URLs need SAS signing.

    // Catalogue sources (lsdtt_dem_source_info): national and polar DEMs,
    // named country first. New values are appended so existing ones keep
    // their meaning.
    us_3dep,              // USGS 3DEP seamless, via Planetary Computer. No key.
    arctic_arcticdem,     // PGC ArcticDEM mosaic v4.1. No key.
    antarctic_rema,       // PGC REMA mosaic v2.0. No key.
    no_dtm,               // Kartverket NHM DTM 1 m (Norway), WCS. No key.
    nl_ahn,               // AHN 0.5 m (Netherlands), WCS. No key.
    es_mdt,               // IGN MDT05 5 m (Spain), WCS. No key.
    pl_nmt,               // GUGiK NMT 1 m (Poland), WCS. No key.
    ca_hrdem,             // NRCan HRDEM mosaic 1/2 m (Canada), STAC. No key.
    ch_swissalti3d,       // swisstopo swissALTI3D 0.5/2 m (Switzerland), STAC. No key.
    nz_lidar,             // LINZ national 1 m DEM/DSM (New Zealand), public S3. No key.
    gb_sct_lidar,         // Scottish Remote Sensing Portal lidar, public S3. No key.
    jp_gsi                // GSI elevation tiles 1/5/10 m (Japan), XYZ PNG. No key.
};

// Parse a user-facing source name ("cop30_aws", ...). Returns false if unknown.
bool parse_dem_source(const std::string& name, lsdtt_dem_source& out);

// Name of a source, for messages and help text.
std::string dem_source_name(lsdtt_dem_source source);


// -----------------------------------------------------------------------------
// lsdtt_dem_source_info
//
// One row of the source catalogue. Every source has a row, so `-h dem_source`
// can describe them offline and messages can say what a source is. The three
// original sources (cop30_aws, opentopography, planetary_computer) are marked
// legacy_path: they keep their own fetch code, and their row is used only for
// help text, the banner and provenance.
//
// The catalogue sources share one fetch path: a pre-flight check (coverage,
// resolution, surface model, size) that needs no network, discovery of the
// source files, one GDALWarp straight from them onto the target grid, and a
// provenance stamp on the output.
//
// `description` must say in words whether the source is a surface model (DSM:
// canopy and buildings included) or a terrain model (DTM: bare earth), and
// give its vertical datum: users analyse these as ground surfaces.
// -----------------------------------------------------------------------------
struct lsdtt_dem_source_info
{
    lsdtt_dem_source id = lsdtt_dem_source::cop30_aws;
    std::string name;                 // "us_3dep"
    std::string title;                // "USGS 3DEP seamless DEM"
    std::string provider;             // who distributes it
    std::string licence_attribution;  // printed after a successful fetch
    std::string coverage_text;        // "CONUS, Alaska, Hawaii, Puerto Rico"
    std::string description;          // help text: DSM or DTM in words, datum
    std::string accessed;             // when the access details were checked
    bool legacy_path = false;

    std::vector<double> native_resolutions_metres;   // nominal, e.g. {10, 30}
    double default_resolution_metres = 0.0;
    std::vector<std::string> surface_models;         // "dtm" and/or "dsm"
    std::string default_surface_model;
    std::string vertical_datum;                      // "NAVD88", "EGM2008", ...
    std::string vertical_units = "metre";            // what the source must carry

    // Coarse WGS84 boxes (west, south, east, north) for the offline coverage
    // check: a request outside all of them fails before any network request.
    std::vector<std::array<double, 4>> coverage_envelopes;

    // STAC access (catalogue sources).
    std::string stac_search_url;
    std::string stac_collection;
    std::string stac_asset;
    bool        stac_planetary_computer_signing = false;
    // Hosts that publish one collection per resolution (PGC): the collection
    // for each native resolution, used instead of stac_collection.
    std::vector<std::pair<double, std::string>> stac_collection_by_resolution;
    // The asset per surface model, where one item carries both (HRDEM:
    // "dtm" -> "dtm", "dsm" -> "dsm"); otherwise stac_asset.
    std::vector<std::pair<std::string, std::string>> stac_asset_by_model;
    // Assets are keyed by file name, one per resolution and format
    // (swisstopo): take the GeoTIFF asset whose gsd is the resolution.
    bool stac_asset_by_gsd = false;
    // Items carry no gsd (the collection fixes the resolution): no gsd filter.
    bool stac_items_carry_gsd = true;
    // Several acquisitions per tile, the year a token of the item id
    // ("swissalti3d_2024_2623-1091"): keep the latest per tile.
    bool stac_latest_per_tile = false;
    // STAC fields extension: return only these item fields. PGC items list
    // every source scene, tens of kB each; the fetch needs none of it.
    std::vector<std::string> stac_fields_include;

    // The source grid is projected (metre cells), not geographic: the
    // download estimate counts cells by area rather than by degrees.
    bool native_grid_projected = false;

    // WCS access (catalogue sources served by a Web Coverage Service; a row
    // with a wcs_url uses it instead of STAC). The native grid is requested
    // exactly, so the server does no resampling, in chunks of at most
    // wcs_max_chunk_side pixels a side that abut on native pixel edges.
    std::string wcs_url;                       // base URL, without a query
    std::string wcs_version;                   // "1.0.0" or "2.0.1"
    std::vector<std::pair<std::string, std::string>> wcs_coverage_by_model;  // "dtm" -> id
    std::string wcs_format;                    // the GetCoverage FORMAT value
    bool   wcs_multipart_ascii = false;        // response is multipart MIME around an ESRI ASCII grid
    int    wcs_native_epsg = 0;
    double wcs_origin_x = 0.0;                 // top-left corner of the grid (a pixel edge)
    double wcs_origin_y = 0.0;
    double wcs_pixel_metres = 0.0;
    std::string wcs_x_axis = "x";              // 2.0.1 SUBSET label for easting
    std::string wcs_y_axis = "y";              // 2.0.1 SUBSET label for northing
    int    wcs_max_chunk_side = 2000;
    double wcs_bytes_per_cell = 4.0;           // for the download estimate
    // The service has no nodata value and returns 0 outside the country (and
    // over sea): 0 is set as nodata, so a border is not a cliff to 0 m.
    bool   wcs_zero_is_nodata = false;
    // Grid-keyed public bucket (no catalogue): file names are computed from
    // coordinates in grid_epsg, and a file that does not exist means no data
    // there. Collections are listed oldest first; one is used per fetch
    // (source_collection_choice), never a mix.
    struct grid_collection
    {
        std::string label;         // "phase-3", recorded in the provenance
        std::string key_template;  // tokens: {model} {MODEL} {sheet} {os10} {os5} {os1}
        double      tile_width = 0.0;    // metres, easting
        double      tile_height = 0.0;   // metres, northing
        double      origin_x = 0.0;      // a tile corner (west edge)
        double      origin_y = 0.0;      // a tile corner (south edge)
    };
    std::string grid_base_url;               // https base; a file is base + key
    int         grid_epsg = 0;
    std::vector<grid_collection> grid_collections;
    // The bucket's name for a surface model, where it differs ("dtm" -> "dem").
    std::vector<std::pair<std::string, std::string>> grid_model_names;
    std::string grid_sentinel_key;           // a file that exists: tells "no data
                                             // here" from "host unreachable"
    bool requires_lerc = false;              // files are LERC-compressed

    // XYZ elevation PNG tiles (GSI encoding, EPSG:3857): one layer per
    // resolution. A tile that does not exist (404) means no data there.
    struct xyz_layer
    {
        double      resolution_metres = 0.0;   // nominal
        int         zoom = 0;
        std::string label;                     // "dem5a_png", recorded in the provenance
        std::string url_template;              // with {z} {x} {y}
    };
    std::vector<xyz_layer> xyz_layers;
};

// The whole catalogue, in the order shown to users.
const std::vector<lsdtt_dem_source_info>& dem_source_catalogue();

// The row for one source (every lsdtt_dem_source has one).
const lsdtt_dem_source_info& dem_source_info(lsdtt_dem_source source);

// Comma-separated list of every valid dem_source name.
std::string dem_source_name_list();


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

    // Native product resolution to fetch (nominal metres). 0: the source's
    // default product; for planetary_computer, 0 keeps every item the search
    // returns (the original behaviour) and a value keeps only items whose
    // STAC gsd equals it. Ignored by cop30_aws and opentopography.
    double source_resolution_metres = 0.0;

    // "source_default", "dtm" or "dsm". Asking for a product a source does
    // not have is an error.
    std::string surface_model_choice = "source_default";
    // Sources published as several survey collections (gb_sct_lidar): "latest"
    // takes the newest collection with tiles in the area; a collection label
    // ("phase-3", "national-lidar-programme") takes that one only.
    std::string source_collection_choice = "latest";

    // Size limits, checked before any download (and, where lsdtt3 does the
    // HTTP itself, while downloading). 0 disables a limit.
    std::int64_t max_output_pixels      = 250000000;
    double       max_download_megabytes = 2048.0;
    int          max_tiles              = 256;
    int          request_delay_milliseconds = 200;

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

    // Frame mode (lsdtt_oblique_frame.hpp). When use_frame is true the output
    // is warped straight from the source tiles into the frame (one resampling
    // pass; -te/-ts in frame metres, no -tap): the bounding box above is
    // replaced by the frame's WGS84 envelope (fetch_bbox_for_frame), size_px
    // must be 0, and target_epsg becomes the frame file's rotate-back target
    // instead of the output CRS (0: the analysis CRS of the frame centre).
    // frame.cell_size must equal grid_spacing.
    bool use_frame = false;
    lsdtt_oblique_frame frame;

    // Frame mode only: where to write the rotation frame file (empty: none),
    // and its provenance (rotated_dem_fname is filled from out_path).
    std::string frame_file_path;
    lsdtt_rotation_provenance frame_provenance;
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

    int    epsg         = 0;       // output CRS; 0 in frame mode
    bool   in_frame     = false;   // output is in options.frame
    int    rotate_back_epsg = 0;   // frame mode: the frame file's target
    int    n_cols       = 0;
    int    n_rows       = 0;
    double grid_spacing = 0.0;

    std::vector<std::string> tiles;   // source tiles used

    // Non-fatal advisories for the caller to print (partial coverage, mixed
    // resolutions, overview reads).
    std::vector<std::string> warnings;

    // What was written into the output's metadata (LSDTT_* items) and should
    // go into a run summary.
    std::vector<std::pair<std::string, std::string>> provenance;
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

// The WGS84 envelope of a frame's rectangle (edges densified, 50 points per
// side) padded by pad_metres on every side, for the tile request. Throws
// std::runtime_error if the transform fails.
void fetch_bbox_for_frame(const lsdtt_oblique_frame& frame, double pad_metres,
                          double& west, double& south, double& east, double& north);

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

// True if this GDAL build can decode LERC-compressed GeoTIFFs (needed by
// sources with requires_lerc, e.g. nz_lidar).
bool gdal_can_read_lerc();

// The screen warning for a source whose files are LERC-compressed: what it
// needs, whether this GDAL has it, and the alternative. Empty for others.
std::string lerc_source_warning(lsdtt_dem_source source);


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
    std::string geometry_kind;   // "polygon", "points" or "raster"
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
//   - a raster (e.g. the DEM itself, a GeoTIFF) -> the WGS84 envelope of its
//     footprint, every edge densified before reprojection; geometry_kind is
//     "raster". Tried only when the file opens as no vector layer, so every
//     vector input behaves as before.
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
