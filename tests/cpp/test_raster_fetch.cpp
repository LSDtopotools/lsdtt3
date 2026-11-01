// test_raster_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Part of the flow-routing benchmark (docs/developer_notes/testing_and_benchmarking.md): tests for native DEM fetching.
//
// Two modes, because they have very different requirements:
//
//   --offline  (default)  Pure logic plus the graceful-degradation paths. Needs
//                         no network and runs in the default CTest suite. This
//                         is where the invariants that must never regress live.
//
//   --network             Actually fetches a small window from Copernicus
//                         GLO-30 on AWS. Opt-in via LSDTT_ENABLE_NETWORK_TESTS,
//                         and exits 77 (CTest skip) when the source is
//                         unreachable, so it can never turn CI red.
//
// Usage: lsdtt-test-raster-fetch [--offline|--network] [output_dir]
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <gdal.h>
#include <cpl_conv.h>
#include <cpl_vsi.h>
#include <ogr_srs_api.h>

#include "lsdtt_raster_fetch.hpp"
#include "lsdtt_raster_fetch_detail.hpp"

namespace
{

constexpr int k_skip_exit = 77;   // CTest's skip code

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

std::string join(const std::vector<std::string>& items)
{
    std::string out;
    for (std::size_t i = 0; i < items.size(); ++i)
    {
        if (i) out += ", ";
        out += items[i];
    }
    return out;
}


// =============================================================================
// UTM zone derivation
// =============================================================================
void test_utm_zone_derivation()
{
    struct utm_case { double lon; double lat; int expected; const char* what; };

    // Zone boundaries are the interesting part: an off-by-one here silently
    // reprojects a fixture into the neighbouring zone, which distorts distances
    // rather than failing outright.
    const std::vector<utm_case> cases = {
        { 10.40,  46.50, 32632, "Italian Alps -> UTM 32N" },
        { 25.40,  61.50, 32635, "Finnish Lakeland -> UTM 35N" },
        {-3.19,   55.95, 32630, "Edinburgh -> UTM 30N" },
        {170.50, -45.87, 32759, "New Zealand -> UTM 59S" },
        {  0.00,   0.10, 32631, "prime meridian, northern -> zone 31N" },
        { -0.01,  -0.10, 32730, "just west of prime meridian, southern -> 30S" },
        {-180.00,  10.0, 32601, "western limit -> zone 1" },
        { 179.99,  10.0, 32660, "eastern limit -> zone 60" },
        {   6.00,  10.0, 32632, "exact zone boundary rounds up into zone 32" },
        {   0.00,   0.00, 32631, "equator counts as northern hemisphere" },
    };

    for (const utm_case& c : cases)
    {
        const int got = lsdtt::utm_epsg_for_point(c.lon, c.lat);
        check(std::string("utm zone: ") + c.what,
              got == c.expected,
              "expected " + std::to_string(c.expected) + ", got " + std::to_string(got));
    }

    // 180 degrees exactly must not produce a zone 61, which does not exist.
    const int at_antimeridian = lsdtt::utm_epsg_for_point(180.0, 10.0);
    check("utm zone: +180 does not overflow to zone 61",
          at_antimeridian >= 32601 && at_antimeridian <= 32660,
          "got EPSG " + std::to_string(at_antimeridian));
}


// =============================================================================
// Tile naming and coverage
// =============================================================================
void test_tile_naming()
{
    {
        const auto tiles = lsdtt::cop30_tiles_for_bbox(10.30, 46.40, 10.50, 46.60);
        check("tiles: bbox inside one tile -> 1 tile",
              tiles.size() == 1 && tiles[0] == "N46_00_E010_00",
              "got [" + join(tiles) + "]");
    }
    {
        // Crossing a boundary is normal for real study areas; a mosaic is
        // required or the fetch silently truncates.
        const auto tiles = lsdtt::cop30_tiles_for_bbox(9.80, 46.40, 10.20, 46.60);
        check("tiles: bbox crossing a longitude boundary -> 2 tiles",
              tiles.size() == 2,
              "got [" + join(tiles) + "]");
    }
    {
        const auto tiles = lsdtt::cop30_tiles_for_bbox(10.30, 45.90, 10.50, 46.10);
        check("tiles: bbox crossing a latitude boundary -> 2 tiles",
              tiles.size() == 2,
              "got [" + join(tiles) + "]");
    }
    {
        const auto tiles = lsdtt::cop30_tiles_for_bbox(9.80, 45.90, 10.20, 46.10);
        check("tiles: bbox crossing both -> 4 tiles",
              tiles.size() == 4,
              "got [" + join(tiles) + "]");
    }
    {
        // Southern/western hemispheres use S and W and must still zero-pad.
        const auto tiles = lsdtt::cop30_tiles_for_bbox(-72.40, -13.60, -72.30, -13.50);
        check("tiles: southern/western hemisphere naming",
              tiles.size() == 1 && tiles[0] == "S14_00_W073_00",
              "got [" + join(tiles) + "]");
    }
    {
        const std::string url = lsdtt::cop30_tile_url("N46_00_E010_00");
        const bool well_formed =
            url.find("/vsicurl/https://") == 0 &&
            url.find("Copernicus_DSM_COG_10_N46_00_E010_00_DEM/") != std::string::npos &&
            url.rfind("Copernicus_DSM_COG_10_N46_00_E010_00_DEM.tif") ==
                url.size() - std::string("Copernicus_DSM_COG_10_N46_00_E010_00_DEM.tif").size();
        check("tile url: well formed /vsicurl path", well_formed, url);
    }
}


// =============================================================================
// Source vocabulary
// =============================================================================
void test_source_parsing()
{
    lsdtt::lsdtt_dem_source source;

    check("source: cop30_aws parses",
          lsdtt::parse_dem_source("cop30_aws", source) &&
          source == lsdtt::lsdtt_dem_source::cop30_aws);

    check("source: unknown name rejected",
          !lsdtt::parse_dem_source("not_a_source", source));

    check("source: round-trips through its name",
          lsdtt::dem_source_name(lsdtt::lsdtt_dem_source::cop30_aws) == "cop30_aws");

    check("source: opentopography parses",
          lsdtt::parse_dem_source("opentopography", source) &&
          source == lsdtt::lsdtt_dem_source::opentopography);

    check("source: planetary_computer parses",
          lsdtt::parse_dem_source("planetary_computer", source) &&
          source == lsdtt::lsdtt_dem_source::planetary_computer);

    // CDSE was cut (docs/developer_notes/testing_and_benchmarking.md): a general public account does not unlock the
    // 10 m EU-DEM, so it would have added the only username/password credential
    // path for coverage already free from cop30_aws. It must not be silently
    // accepted and quietly served from a different provider.
    check("source: cdse is rejected, not silently substituted",
          !lsdtt::parse_dem_source("cdse", source));

    // NOTE: there is deliberately no "unimplemented enum value" check here.
    // Every value in lsdtt_dem_source is implemented, and a previous version of
    // this test used a not-yet-implemented value as its subject -- which meant
    // that once the value was implemented, this OFFLINE test silently started
    // performing a real network fetch. Assert on the parsing vocabulary
    // instead, which needs no network by construction.
}


// =============================================================================
// Keyed sources
//
// A missing key must be a clear, actionable error before any network access is
// attempted -- not a skip (which would hide a misconfiguration) and not an
// obscure HTTP failure from the provider.
// =============================================================================
void test_opentopography_requires_key()
{
    lsdtt::lsdtt_fetch_options options;
    options.west = 10.30; options.south = 46.40;
    options.east = 10.40; options.north = 46.50;
    options.source  = lsdtt::lsdtt_dem_source::opentopography;
    options.api_key = "";   // resolved by the caller; empty means none found

    const auto result = lsdtt::fetch_dem(options, "/vsimem/no_key.tif");

    check("opentopography: missing key is an error, not a skip",
          !result.ok && !result.skipped,
          "ok=" + std::to_string(result.ok) +
          " skipped=" + std::to_string(result.skipped));

    check("opentopography: missing-key message names the environment variable",
          result.message.find("OT_API_KEY") != std::string::npos,
          result.message);

    check("opentopography: missing-key message points at the keyless option",
          result.message.find("cop30_aws") != std::string::npos,
          result.message);
}


// =============================================================================
// OpenTopography dataset catalogue
//
// The catalogue is compiled in rather than fetched, so it can drift from the
// upstream /globaldem demtype list. These checks cannot detect that drift --
// only that the table is internally coherent and actually used for validation,
// which is what stops a typo reaching the API as an opaque 400.
// =============================================================================
void test_ot_dataset_catalogue()
{
    const auto& catalogue = lsdtt::ot_dataset_catalogue();

    check("ot catalogue: is populated", catalogue.size() >= 17,
          std::to_string(catalogue.size()) + " entries");

    bool all_fields_set = true;
    bool all_cited      = true;
    std::set<std::string> seen;
    bool unique = true;

    for (const auto& d : catalogue)
    {
        if (d.demtype.empty() || d.name.empty() || d.resolution.empty() ||
            d.gdal_type.empty() || d.coverage.empty() ||
            d.description.empty() || d.max_area_km2 <= 0.0)
        {
            all_fields_set = false;
        }
        // A citation without a resolvable identifier is not a citation.
        if (d.citation.find("https://doi.org/") == std::string::npos)
        {
            all_cited = false;
        }
        if (!seen.insert(d.demtype).second) unique = false;
    }

    check("ot catalogue: every entry has all descriptive fields", all_fields_set);
    check("ot catalogue: every entry carries a DOI citation", all_cited);
    check("ot catalogue: demtype values are unique", unique);

    check("ot catalogue: lookup finds a known demtype",
          lsdtt::ot_dataset_info("COP30") != nullptr);

    // The API is case sensitive, so the lookup must be too: silently accepting
    // "cop30" here would send a request that fails downstream instead.
    check("ot catalogue: lookup is case sensitive",
          lsdtt::ot_dataset_info("cop30") == nullptr);

    check("ot catalogue: lookup rejects an unknown demtype",
          lsdtt::ot_dataset_info("SRTMGL1_typo") == nullptr);

    const std::string list = lsdtt::ot_dataset_demtype_list();
    check("ot catalogue: demtype list names each entry",
          list.find("COP30") != std::string::npos &&
          list.find("GEDI_L3") != std::string::npos,
          list);

    // A bad demtype must fail locally, before the key is even used.
    lsdtt::lsdtt_fetch_options options;
    options.west = 10.30; options.south = 46.40;
    options.east = 10.40; options.north = 46.50;
    options.source      = lsdtt::lsdtt_dem_source::opentopography;
    options.api_key     = "dummy_key_not_used";
    options.ot_dem_type = "NOT_A_DATASET";

    const auto result = lsdtt::fetch_dem(options, "/vsimem/bad_demtype.tif");

    check("ot catalogue: bad demtype is an error, not a skip",
          !result.ok && !result.skipped,
          "ok=" + std::to_string(result.ok) +
          " skipped=" + std::to_string(result.skipped));

    check("ot catalogue: bad-demtype message lists the valid values",
          result.message.find("NOT_A_DATASET") != std::string::npos &&
          result.message.find("COP30")         != std::string::npos,
          result.message);
}


// =============================================================================
// Option validation
//
// A malformed request is a genuine error (ok=false, skipped=false). Only an
// unavailable dependency may set skipped -- otherwise a typo in a bounding box
// would quietly register as "skipped" and the test would pass.
// =============================================================================
void test_option_validation()
{
    struct bad_case { double w, s, e, n; const char* what; };
    const std::vector<bad_case> cases = {
        { 10.5, 46.4, 10.3, 46.6, "west >= east" },
        { 10.3, 46.6, 10.5, 46.4, "south >= north" },
        { 10.3, 46.4, 10.3, 46.6, "zero-width bbox" },
        {-190.0, 46.4, 10.5, 46.6, "longitude below -180" },
        { 10.3, 46.4, 10.5, 91.0,  "latitude above 90" },
    };

    for (const bad_case& c : cases)
    {
        lsdtt::lsdtt_fetch_options options;
        options.west = c.w; options.south = c.s;
        options.east = c.e; options.north = c.n;

        const auto result = lsdtt::fetch_dem(options, "/vsimem/bad.tif");
        check(std::string("validation: ") + c.what + " is an error",
              !result.ok && !result.skipped,
              "ok=" + std::to_string(result.ok) +
              " skipped=" + std::to_string(result.skipped) +
              " msg=" + result.message);
    }

    {
        lsdtt::lsdtt_fetch_options options;
        options.west = 10.30; options.south = 46.40;
        options.east = 10.40; options.north = 46.50;
        options.grid_spacing = 0.0;
        const auto result = lsdtt::fetch_dem(options, "/vsimem/bad.tif");
        check("validation: non-positive grid_spacing is an error",
              !result.ok && !result.skipped, result.message);
    }
    {
        lsdtt::lsdtt_fetch_options options;
        options.west = 10.30; options.south = 46.40;
        options.east = 10.40; options.north = 46.50;
        const auto result = lsdtt::fetch_dem(options, "");
        check("validation: empty out_path is an error",
              !result.ok && !result.skipped, result.message);
    }
}


// =============================================================================
// Graceful degradation
//
// The single most important behaviour in this module: an unreachable source
// must yield skipped=true, never a hard failure and never a throw. Tested
// offline by pointing at a host that cannot resolve, which fails the same way
// whether or not this machine has network.
// =============================================================================
void test_unreachable_source_skips()
{
    // Redirect the vsicurl prefix at a guaranteed-invalid TLD. .invalid is
    // reserved by RFC 2606 precisely so it can never resolve.
    CPLSetThreadLocalConfigOption("GDAL_HTTP_CONNECTTIMEOUT", "2");
    CPLSetThreadLocalConfigOption("GDAL_HTTP_TIMEOUT", "3");

    lsdtt::lsdtt_fetch_options options;
    // A bbox over open ocean far from any Copernicus land tile: the tile does
    // not exist, so the open fails without needing an unreachable host.
    options.west  = -140.00; options.south = -40.00;
    options.east  = -139.90; options.north = -39.90;
    options.http_connect_timeout_seconds = 2;
    options.http_timeout_seconds = 5;
    options.http_max_retry = 0;

    const auto result = lsdtt::fetch_dem(options, "/vsimem/unreachable.tif");

    check("degradation: unobtainable tile yields skipped, not failure",
          !result.ok && result.skipped,
          "ok=" + std::to_string(result.ok) +
          " skipped=" + std::to_string(result.skipped) +
          " msg=" + result.message);

    check("degradation: skip message names the tile(s) for diagnosis",
          result.message.find("S40_00_W140_00") != std::string::npos,
          result.message);

    CPLSetThreadLocalConfigOption("GDAL_HTTP_CONNECTTIMEOUT", nullptr);
    CPLSetThreadLocalConfigOption("GDAL_HTTP_TIMEOUT", nullptr);
}


// =============================================================================
// No credentials anywhere near this module
//
// lsdtt3 CLIs echo every parameter at the end of a run and users commit their
// parameter files, so a credential must never be reachable from options or
// results. Phase 1 has no keyed sources at all; this check exists so that
// adding one later cannot quietly route a secret through these structs.
// =============================================================================
void test_no_credential_leakage()
{
    lsdtt::lsdtt_fetch_options options;
    options.west = 10.30; options.south = 46.40;
    options.east = 10.40; options.north = 46.50;

    // Anything that looks like a secret must not survive into a loggable field.
    const auto result = lsdtt::fetch_dem(options, "");

    const std::string haystack = result.message + result.path;
    const bool clean =
        haystack.find("API_Key") == std::string::npos &&
        haystack.find("api_key") == std::string::npos &&
        haystack.find("Bearer")  == std::string::npos &&
        haystack.find("password") == std::string::npos &&
        haystack.find("token")   == std::string::npos;

    check("credentials: result carries nothing credential-shaped", clean, haystack);
}


// =============================================================================
// Bounds from a vector file
//
// The headline property: for two points, BOTH diagonals must give the same
// box. Lower-left/upper-right and lower-right/upper-left are just min/max of
// the same two x and y values, so there is no orientation to detect and no way
// for a user to hold it wrong.
// =============================================================================
void test_bounds_from_vector(const std::string& data_dir)
{
    // The extent the fixtures were built from.
    const double expect_w = 10.3373594105;
    const double expect_s = 46.4568810636;
    const double expect_e = 10.4626405895;
    const double expect_n = 46.5431189364;

    auto close = [](double a, double b, double tol) { return std::fabs(a - b) <= tol; };

    auto check_extent = [&](const std::string& label,
                            const lsdtt::lsdtt_bounds_from_vector_result& r,
                            double tol)
    {
        check(label + ": bounds match the fixture extent",
              r.ok &&
              close(r.west, expect_w, tol) && close(r.south, expect_s, tol) &&
              close(r.east, expect_e, tol) && close(r.north, expect_n, tol),
              "got [" + std::to_string(r.west) + ", " + std::to_string(r.south) +
              ", " + std::to_string(r.east) + ", " + std::to_string(r.north) +
              "] msg=" + r.message);
    };

    // --- polygon in WGS84, FlatGeobuf ---
    {
        const auto r = lsdtt::bounds_from_vector_file(data_dir + "/bounds_polygon.fgb");
        check_extent("polygon fgb", r, 1e-9);
        check("polygon fgb: reported as polygon", r.geometry_kind == "polygon",
              r.geometry_kind);
    }

    // --- polygon in UTM 32N, GeoPackage: must be reprojected to WGS84 ---
    {
        const auto r = lsdtt::bounds_from_vector_file(data_dir + "/bounds_polygon_utm.gpkg");
        // Reprojection round-trip, plus a UTM rectangle is not exactly a WGS84
        // rectangle, so the envelope is very slightly larger. Metres, not
        // degrees, is the right scale of tolerance here.
        check_extent("polygon gpkg UTM", r, 1e-4);
        check("polygon gpkg UTM: source CRS reported",
              r.source_crs.find("32632") != std::string::npos, r.source_crs);
    }

    // --- a raster (bounds_fname: <DEM>): a 10 x 10 GeoTIFF in UTM 32N whose
    // extent is the UTM envelope of bounds_polygon_utm.gpkg. Its WGS84
    // footprint is therefore larger than that polygon's; the expected edges are
    // its four corners transformed by gdaltransform (PROJ, GDAL 3.4.1):
    // SW 10.3352468907 / 46.4569057500, NW 10.3373970652 / 46.5446520523,
    // NE 10.4649546110 / 46.5430893193, SE 10.4625995289 / 46.4553477705.
    // The edges are densified, so the tolerance allows the 1e-7 degree bow
    // of a 10 km UTM edge. ---
    {
        const auto r = lsdtt::bounds_from_vector_file(data_dir + "/bounds_raster_utm.tif");
        check("raster GeoTIFF UTM: footprint matches the transformed corners",
              r.ok && close(r.west, 10.3352468907, 1e-6) && close(r.south, 46.4553477705, 1e-6) &&
              close(r.east, 10.4649546110, 1e-6) && close(r.north, 46.5446520523, 1e-6),
              "got [" + std::to_string(r.west) + ", " + std::to_string(r.south) + ", " +
              std::to_string(r.east) + ", " + std::to_string(r.north) + "] msg=" + r.message);
        check("raster GeoTIFF UTM: reported as raster with its CRS",
              r.geometry_kind == "raster" && r.source_crs.find("32632") != std::string::npos,
              r.geometry_kind + " " + r.source_crs);
    }

    // --- two points, lower-left + upper-right ---
    const auto ll_ur = lsdtt::bounds_from_vector_file(data_dir + "/bounds_points_ll_ur.gpkg");
    check_extent("points LL/UR", ll_ur, 1e-9);
    check("points LL/UR: reported as points", ll_ur.geometry_kind == "points",
          ll_ur.geometry_kind);

    // --- two points, the OTHER diagonal: lower-right + upper-left ---
    const auto lr_ul = lsdtt::bounds_from_vector_file(data_dir + "/bounds_points_lr_ul.fgb");
    check_extent("points LR/UL", lr_ul, 1e-9);

    // --- the property that makes "autodetection" unnecessary ---
    check("points: both diagonals give an identical box",
          ll_ur.ok && lr_ul.ok &&
          close(ll_ur.west,  lr_ul.west,  1e-12) &&
          close(ll_ur.south, lr_ul.south, 1e-12) &&
          close(ll_ur.east,  lr_ul.east,  1e-12) &&
          close(ll_ur.north, lr_ul.north, 1e-12),
          "LL/UR vs LR/UL differ");

    // --- three points cannot define a box unambiguously ---
    {
        const auto r = lsdtt::bounds_from_vector_file(data_dir + "/bounds_points_three.gpkg");
        check("points: three points rejected",
              !r.ok && r.message.find("exactly two") != std::string::npos,
              r.message);
    }

    // --- missing file ---
    {
        const auto r = lsdtt::bounds_from_vector_file(data_dir + "/no_such_file.fgb");
        check("bounds: missing file is a clear error",
              !r.ok && !r.message.empty(), r.message);
    }

    // --- empty path ---
    {
        const auto r = lsdtt::bounds_from_vector_file("");
        check("bounds: empty path is an error", !r.ok, r.message);
    }
}



// =============================================================================
// Frozen Copernicus path. The cache key, the Planetary Computer search body
// and the warp arguments determine what the cop30_aws / planetary_computer
// path asks GDAL to do, so pinning them pins that path on any machine. Any
// difference here is a regression of the frozen-Copernicus requirement.
// =============================================================================
lsdtt::lsdtt_fetch_options ortler_options(lsdtt::lsdtt_dem_source source)
{
    lsdtt::lsdtt_fetch_options o;
    o.west = 10.337; o.south = 46.457; o.east = 10.463; o.north = 46.543;
    o.source = source;
    o.target_epsg = 32632;
    o.grid_spacing = 30.0;
    o.size_px = 128;
    return o;
}

void test_copernicus_frozen()
{
    using lsdtt::lsdtt_dem_source;
    namespace fd = lsdtt::fetch_detail;

    const auto cop = ortler_options(lsdtt_dem_source::cop30_aws);
    const auto pc  = ortler_options(lsdtt_dem_source::planetary_computer);

    check("frozen: cop30_aws cache key",
          fd::fetch_cache_key(cop) ==
          "cop30_aws_10_337000_46_457000_10_463000_46_543000_epsg32632_res30_0000"
          "_cubic_px128_nd_9999_0000.tif", fd::fetch_cache_key(cop));
    check("frozen: planetary_computer cache key",
          fd::fetch_cache_key(pc) ==
          "planetary_computer_10_337000_46_457000_10_463000_46_543000_epsg32632"
          "_res30_0000_cubic_px128_nd_9999_0000_cop_dem_glo_30_data.tif",
          fd::fetch_cache_key(pc));
    check("frozen: planetary_computer search body",
          fd::pc_stac_search_body(pc) ==
          "{\"collections\":[\"cop-dem-glo-30\"],\"bbox\":[10.3370000000,"
          "46.4570000000,10.4630000000,46.5430000000],\"limit\":100}",
          fd::pc_stac_search_body(pc));

    const std::vector<std::string> expected_warp = {
        "-t_srs", "EPSG:32632", "-tr", "30.000000", "30.000000", "-r", "cubic",
        "-tap", "-dstnodata", "-9999.000000", "-of", "GTiff",
        "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES"};
    check("frozen: warp arguments", fd::fetch_warp_arguments(cop, 32632) == expected_warp,
          join(fd::fetch_warp_arguments(cop, 32632)));

    // Frame mode: -te/-ts in frame metres, no -tap. The WKT is PROJ's, so only
    // the structure and the numbers lsdtt3 writes are pinned.
    lsdtt::lsdtt_fetch_options framed = cop;
    framed.size_px = 0;
    framed.use_frame = true;
    framed.frame.centre.longitude = 10.4;
    framed.frame.centre.latitude  = 46.5;
    framed.frame.azimuth_degrees  = 30.0;
    framed.frame.x_min = -3000.0; framed.frame.x_max = 3000.0;
    framed.frame.y_min = -1500.0; framed.frame.y_max = 1500.0;
    framed.frame.cell_size = 30.0;
    const auto fa = fd::fetch_warp_arguments(framed, 0);
    check("frozen: frame warp arguments",
          fa.size() == 20 && fa[0] == "-t_srs" && fa[2] == "-te" &&
          fa[3] == "-3000.000000" && fa[4] == "-1500.000000" &&
          fa[5] == "3000.000000" && fa[6] == "1500.000000" &&
          fa[7] == "-ts" && fa[8] == "200" && fa[9] == "100" &&
          fa[10] == "-r" && fa[11] == "cubic" &&
          std::find(fa.begin(), fa.end(), "-tap") == fa.end(),
          join(fa));
}


// Raster content of a written DEM: what the frozen-Copernicus guard compares.
// Metadata items are deliberately not part of it (provenance is stamped on
// every output), so two files are "raster identical" when everything a user
// computes from them is identical.
struct raster_content
{
    bool ok = false;
    int cols = 0, rows = 0;
    double geotransform[6] = {0, 0, 0, 0, 0, 0};
    std::string wkt;
    int has_nodata = 0;
    double nodata = 0.0;
    int data_type = 0;
    std::vector<float> values;
};

raster_content read_raster_content(const std::string& path)
{
    raster_content c;
    GDALDatasetH ds = GDALOpen(path.c_str(), GA_ReadOnly);
    if (!ds) return c;
    c.cols = GDALGetRasterXSize(ds);
    c.rows = GDALGetRasterYSize(ds);
    GDALGetGeoTransform(ds, c.geotransform);
    const char* wkt = GDALGetProjectionRef(ds);
    c.wkt = wkt ? wkt : "";
    GDALRasterBandH band = GDALGetRasterBand(ds, 1);
    c.nodata = GDALGetRasterNoDataValue(band, &c.has_nodata);
    c.data_type = static_cast<int>(GDALGetRasterDataType(band));
    c.values.assign(static_cast<std::size_t>(c.cols) * c.rows, 0.0f);
    c.ok = GDALRasterIO(band, GF_Read, 0, 0, c.cols, c.rows, c.values.data(),
                        c.cols, c.rows, GDT_Float32, 0, 0) == CE_None;
    GDALClose(ds);
    return c;
}

bool raster_identical(const raster_content& a, const raster_content& b, std::string& why)
{
    if (!a.ok || !b.ok) { why = "unreadable"; return false; }
    if (a.cols != b.cols || a.rows != b.rows) { why = "size"; return false; }
    for (int i = 0; i < 6; ++i)
        if (a.geotransform[i] != b.geotransform[i]) { why = "geotransform"; return false; }
    if (a.wkt != b.wkt) { why = "CRS"; return false; }
    if (a.has_nodata != b.has_nodata || a.nodata != b.nodata) { why = "nodata"; return false; }
    if (a.data_type != b.data_type) { why = "data type"; return false; }
    if (std::memcmp(a.values.data(), b.values.data(), a.values.size() * sizeof(float)) != 0)
    { why = "band values"; return false; }
    return true;
}

// FNV-1a (64 bit) over the band values and the geotransform: a content hash
// that ignores metadata. Its value depends on the GDAL build (warp kernels,
// PROJ transforms), so it is pinned per GDAL release.
std::string raster_content_hash(const raster_content& c)
{
    std::uint64_t h = 1469598103934665603ull;
    auto feed = [&h](const void* data, std::size_t n) {
        const unsigned char* p = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < n; ++i) { h ^= p[i]; h *= 1099511628211ull; }
    };
    feed(c.geotransform, sizeof(c.geotransform));
    feed(c.values.data(), c.values.size() * sizeof(float));
    char buffer[17];
    std::snprintf(buffer, sizeof(buffer), "%016llx", static_cast<unsigned long long>(h));
    return buffer;
}


// =============================================================================
// The DEM source catalogue and the catalogue path (offline)
// =============================================================================
bool contains(const std::string& text, const std::string& part)
{
    return text.find(part) != std::string::npos;
}

void test_dem_source_catalogue()
{
    for (const lsdtt::lsdtt_dem_source_info& info : lsdtt::dem_source_catalogue())
    {
        lsdtt::lsdtt_dem_source parsed;
        check("catalogue: " + info.name + " parses and round-trips",
              lsdtt::parse_dem_source(info.name, parsed) && parsed == info.id &&
              lsdtt::dem_source_name(parsed) == info.name);
        check("catalogue: " + info.name + " is described",
              !info.title.empty() && !info.description.empty() &&
              !info.coverage_text.empty() && !info.licence_attribution.empty() &&
              !info.accessed.empty() && !info.coverage_envelopes.empty());
        if (info.legacy_path) continue;
        bool default_offered = false;
        for (double r : info.native_resolutions_metres)
            if (r == info.default_resolution_metres) default_offered = true;
        check("catalogue: " + info.name + " default resolution is offered", default_offered);
        check("catalogue: " + info.name + " default surface model is offered",
              std::find(info.surface_models.begin(), info.surface_models.end(),
                        info.default_surface_model) != info.surface_models.end());
        for (const std::string& m : info.surface_models)
        {
            check("catalogue: " + info.name + " says " + m + " in words",
                  m == "dsm" ? contains(info.description, "SURFACE MODEL")
                             : contains(info.description, "TERRAIN MODEL"));
        }
        check("catalogue: " + info.name + " has a datum and metres",
              !info.vertical_datum.empty() && info.vertical_units == "metre");
        check("catalogue: " + info.name + " has its STAC access",
              !info.stac_search_url.empty() && !info.stac_collection.empty() &&
              !info.stac_asset.empty());
    }
    lsdtt::lsdtt_dem_source parsed;
    check("catalogue: cdse is still not a source", !lsdtt::parse_dem_source("cdse", parsed));
    check("catalogue: the name list has us_3dep",
          contains(lsdtt::dem_source_name_list(), "us_3dep"));
    check("catalogue: cop30_aws is a DSM",
          contains(lsdtt::dem_source_info(lsdtt::lsdtt_dem_source::cop30_aws).description,
                   "SURFACE MODEL"));
}

// Pre-flight refusals need no network and are errors, never skips.
void test_catalogue_preflight()
{
    auto colorado = []() {
        lsdtt::lsdtt_fetch_options o;
        o.west = -105.30; o.south = 39.95; o.east = -105.25; o.north = 40.00;
        o.source = lsdtt::lsdtt_dem_source::us_3dep;
        o.grid_spacing = 10.0;
        return o;
    };
    auto refused = [](const lsdtt::lsdtt_fetch_result& r, const std::string& part) {
        return !r.ok && !r.skipped && contains(r.message, part);
    };

    auto o = colorado();
    o.west = 10.3; o.east = 10.4; o.south = 46.4; o.north = 46.5;
    auto r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("preflight: us_3dep outside its coverage is refused",
          refused(r, "outside us_3dep coverage") && contains(r.message, "cop30_aws"),
          r.message);

    o = colorado();
    o.source_resolution_metres = 5.0;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("preflight: an unoffered resolution lists the offered ones",
          refused(r, "offers 10, 30 m"), r.message);

    o = colorado();
    o.surface_model_choice = "dsm";
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("preflight: a missing surface model is refused",
          refused(r, "no dsm product"), r.message);

    o = colorado();
    o.source = lsdtt::lsdtt_dem_source::cop30_aws;
    o.surface_model_choice = "dtm";
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("preflight: Copernicus is a DSM only", refused(r, "DSM only"), r.message);

    o = colorado();
    o.surface_model_choice = "bare";
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("preflight: an unknown surface model is refused",
          refused(r, "source_default, dtm or dsm"), r.message);

    o = colorado();
    o.max_output_pixels = 1000;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("preflight: max_output_pixels refuses before download",
          refused(r, "max_output_pixels"), r.message);

    o = colorado();
    o.max_download_megabytes = 0.5;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("preflight: max_download_megabytes refuses before download",
          refused(r, "max_download_megabytes"), r.message);

    o = colorado();
    o.source = lsdtt::lsdtt_dem_source::cop30_aws;
    o.grid_spacing = 30.0;
    o.max_download_megabytes = 0.01;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("preflight: the size limits apply to Copernicus too",
          refused(r, "max_download_megabytes"), r.message);
}

void test_fetch_geometry_helpers()
{
    namespace fd = lsdtt::fetch_detail;
    const fd::fetch_extent e =
        fd::target_extent_for_bbox(10.337, 46.457, 10.463, 46.543, 32632, 30.0);
    auto multiple = [](double v, double g) {
        return std::abs(v / g - std::round(v / g)) < 1e-9;
    };
    check("geometry: snapped extent sits on the 30 m lattice",
          multiple(e.x_min, 30.0) && multiple(e.y_min, 30.0) &&
          multiple(e.x_max, 30.0) && multiple(e.y_max, 30.0));
    // Densified envelope of the Ortler box in UTM 32N, snapped outward.
    check("geometry: the Ortler box in UTM 32N, pinned",
          e.x_min == 602490.0 && e.x_max == 612360.0 && e.y_min == 5145690.0 &&
          e.y_max == 5155440.0,
          std::to_string(e.x_min) + " " + std::to_string(e.y_min) + " " +
          std::to_string(e.x_max) + " " + std::to_string(e.y_max));
    check("geometry: output pixel estimate is cols x rows",
          fd::estimate_output_pixels({0.0, 0.0, 300.0, 600.0}, 30.0) == 200);
    check("geometry: download estimate is cells x bytes",
          std::abs(fd::estimate_download_megabytes(1.0e6, 4.0) - 4.0) < 1e-12);
    check("geometry: one square degree at 1 arc-second is 3600^2 cells",
          std::abs(fd::geographic_cell_count(0, 0, 1, 1, 1.0 / 3600.0) - 3600.0 * 3600.0) < 1e-3);
}

void test_stac_parsing()
{
    namespace fd = lsdtt::fetch_detail;
    // Trimmed from a real 3dep-seamless search (hrefs replaced).
    const std::string json = R"({
      "features": [
        {"id": "n41w106-13", "bbox": [-106.0, 40.0, -105.0, 41.0],
         "properties": {"gsd": 10, "datetime": "2015-07-17T00:00:00Z"},
         "assets": {"data": {"href": "https://example.invalid/n41w106-13.tif"}}},
        {"id": "n41w106-1", "bbox": [-106.0, 40.0, -105.0, 41.0],
         "properties": {"gsd": 30, "datetime": "2015-07-17T00:00:00Z"},
         "assets": {"data": {"href": "https://example.invalid/n41w106-1.tif"}}},
        {"id": "n40w106-13", "bbox": [-106.0, 39.0, -105.0, 40.0],
         "properties": {"gsd": 10, "datetime": "2018-03-28T00:00:00Z"},
         "assets": {"data": {"href": "https://example.invalid/n40w106-13.tif"}}},
        {"id": "no-asset", "properties": {"gsd": 10}, "assets": {}}
      ],
      "links": [{"rel": "next", "href": "https://example.invalid/search",
                 "method": "POST", "body": {"token": "next:abc"}}]
    })";
    std::vector<fd::stac_item> items;
    std::vector<std::string> missing;
    std::string next_href, next_body, error;
    const bool ok = fd::parse_stac_items(json, "data", items, missing, next_href, next_body, error);
    check("stac: parses", ok && items.size() == 3, error);
    check("stac: an item without the asset is named, not used",
          missing.size() == 1 && missing[0] == "no-asset");
    check("stac: the next link and its body are carried",
          next_href == "https://example.invalid/search" && contains(next_body, "next:abc"),
          next_href + " " + next_body);
    check("stac: gsd and bbox read",
          items.size() == 3 && items[1].gsd == 30.0 && items[2].bbox[1] == 39.0);

    const auto all = fd::filter_stac_items_by_gsd(items, 0.0);
    check("stac: resolution 0 keeps every item in order (planetary_computer default)",
          all.size() == 3 && all[0].id == "n41w106-13" && all[1].id == "n41w106-1" &&
          all[2].id == "n40w106-13");
    auto ten = fd::filter_stac_items_by_gsd(items, 10.0);
    fd::sort_stac_items(ten);
    check("stac: gsd 10 keeps the 1/3 arc-second items, sorted by id",
          ten.size() == 2 && ten[0].id == "n40w106-13" && ten[1].id == "n41w106-13");
    check("stac: gsd 30 keeps the 1 arc-second item",
          fd::filter_stac_items_by_gsd(items, 30.0).size() == 1);
    const auto gsds = fd::distinct_gsds(items);
    check("stac: distinct resolutions are 10 and 30",
          gsds.size() == 2 && gsds[0] == 10.0 && gsds[1] == 30.0);
}

// A synthetic EPSG:4269 GeoTIFF holding a plane z = a + b*lon + c*lat, with
// overviews, written to /vsimem (nothing on disk, nothing committed).
std::string make_plane_geotiff(const std::string& path, const char* unit, bool overviews)
{
    const int n = 540;                           // 0.05 degrees at 1/3 arc-second
    const double cell = 1.0 / 10800.0;
    const double west = -105.30, north = 40.00;
    GDALDriverH gtiff = GDALGetDriverByName("GTiff");
    GDALDatasetH ds = GDALCreate(gtiff, path.c_str(), n, n, 1, GDT_Float32, nullptr);
    const double gt[6] = {west, cell, 0.0, north, 0.0, -cell};
    GDALSetGeoTransform(ds, const_cast<double*>(gt));
    OGRSpatialReferenceH srs = OSRNewSpatialReference(nullptr);
    OSRImportFromEPSG(srs, 4269);
    char* wkt = nullptr;
    OSRExportToWkt(srs, &wkt);
    GDALSetProjection(ds, wkt);
    CPLFree(wkt);
    OSRDestroySpatialReference(srs);
    std::vector<float> z(static_cast<std::size_t>(n) * n);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
        {
            const double lon = west + (c + 0.5) * cell;
            const double lat = north - (r + 0.5) * cell;
            z[static_cast<std::size_t>(r) * n + c] =
                static_cast<float>(1500.0 + 2000.0 * (lon + 105.30) + 3000.0 * (lat - 39.95));
        }
    GDALRasterBandH band = GDALGetRasterBand(ds, 1);
    GDALRasterIO(band, GF_Write, 0, 0, n, n, z.data(), n, n, GDT_Float32, 0, 0);
    GDALSetRasterNoDataValue(band, -9999.0);
    GDALSetRasterUnitType(band, unit);
    if (overviews)
    {
        int levels[2] = {2, 4};
        GDALBuildOverviews(ds, "AVERAGE", 2, levels, 0, nullptr, nullptr, nullptr);
    }
    GDALClose(ds);
    return path;
}

std::string metadata_item(const std::string& path, const char* key)
{
    GDALDatasetH ds = GDALOpen(path.c_str(), GA_ReadOnly);
    if (!ds) return "";
    const char* v = GDALGetMetadataItem(ds, key, nullptr);
    const std::string out = v ? v : "";
    GDALClose(ds);
    return out;
}

void test_catalogue_warp_offline()
{
    namespace fd = lsdtt::fetch_detail;
    const std::string src = make_plane_geotiff("/vsimem/plane_4269.tif", "m", true);

    lsdtt::lsdtt_fetch_options o;
    o.source = lsdtt::lsdtt_dem_source::us_3dep;
    o.west = -105.29; o.south = 39.96; o.east = -105.26; o.north = 39.99;
    o.target_epsg = 32613;
    o.grid_spacing = 10.0;

    const std::string out = "/vsimem/plane_out_10.tif";
    auto r = fd::warp_catalogue_files(o, {src}, out);
    check("catalogue warp: geographic source written", r.ok, r.message);
    const raster_content c = read_raster_content(out);
    check("catalogue warp: square cells at grid_spacing from an EPSG:4269 source",
          c.ok && std::abs(c.geotransform[1]) == 10.0 && std::abs(c.geotransform[5]) == 10.0 &&
          c.geotransform[2] == 0.0 && c.geotransform[4] == 0.0);
    check("catalogue warp: extent snapped to the 10 m lattice",
          std::abs(c.geotransform[0] / 10.0 - std::round(c.geotransform[0] / 10.0)) < 1e-9 &&
          std::abs(c.geotransform[3] / 10.0 - std::round(c.geotransform[3] / 10.0)) < 1e-9);
    check("catalogue warp: provenance stamped",
          metadata_item(out, "LSDTT_DEM_SOURCE") == "us_3dep" &&
          metadata_item(out, "LSDTT_SURFACE_MODEL") == "DTM" &&
          metadata_item(out, "LSDTT_VERTICAL_DATUM") == "NAVD88" &&
          metadata_item(out, "LSDTT_VERTICAL_UNITS") == "metre" &&
          metadata_item(out, "LSDTT_OVERVIEW_USED") == "false" &&
          metadata_item(out, "LSDTT_RESAMPLING") == "cubic");
    check("catalogue warp: no overview warning at the native grid", r.warnings.empty());

    // The plane is linear in the source's pixel coordinates, which cubic
    // reproduces exactly: check interior cells against the plane at their
    // longitude/latitude.
    {
        OGRSpatialReferenceH utm = OSRNewSpatialReference(nullptr);
        OGRSpatialReferenceH geo = OSRNewSpatialReference(nullptr);
        OSRImportFromEPSG(utm, 32613);
        OSRImportFromEPSG(geo, 4269);
        OSRSetAxisMappingStrategy(utm, OAMS_TRADITIONAL_GIS_ORDER);
        OSRSetAxisMappingStrategy(geo, OAMS_TRADITIONAL_GIS_ORDER);
        OGRCoordinateTransformationH ct = OCTNewCoordinateTransformation(utm, geo);
        double worst = 0.0;
        int n_checked = 0;
        for (int row = 20; row < c.rows - 20; row += 17)
            for (int col = 20; col < c.cols - 20; col += 17)
            {
                double x = c.geotransform[0] + (col + 0.5) * c.geotransform[1];
                double y = c.geotransform[3] + (row + 0.5) * c.geotransform[5];
                OCTTransform(ct, 1, &x, &y, nullptr);
                const double expected = 1500.0 + 2000.0 * (x + 105.30) + 3000.0 * (y - 39.95);
                const float got = c.values[static_cast<std::size_t>(row) * c.cols + col];
                worst = std::max(worst, std::abs(got - expected));
                ++n_checked;
            }
        OCTDestroyCoordinateTransformation(ct);
        OSRDestroySpatialReference(utm);
        OSRDestroySpatialReference(geo);
        check("catalogue warp: a plane is reproduced away from the edges (cubic)",
              n_checked > 50 && worst < 1e-2,
              "max |dz| = " + std::to_string(worst) + " over " + std::to_string(n_checked));
    }

    // Coarser than twice the native resolution: GDAL reads the overviews.
    o.grid_spacing = 30.0;
    r = fd::warp_catalogue_files(o, {src}, "/vsimem/plane_out_30.tif");
    check("catalogue warp: overview read warned and recorded",
          r.ok && r.warnings.size() == 1 && contains(r.warnings[0], "overviews") &&
          metadata_item("/vsimem/plane_out_30.tif", "LSDTT_OVERVIEW_USED") == "true",
          r.message);

    // Feet: refused, not converted.
    const std::string ft = make_plane_geotiff("/vsimem/plane_ft.tif", "ft", false);
    o.grid_spacing = 10.0;
    r = fd::warp_catalogue_files(o, {ft}, "/vsimem/plane_out_ft.tif");
    check("catalogue warp: a source in feet is refused",
          !r.ok && !r.skipped && contains(r.message, "'ft'") &&
          contains(r.message, "metres"), r.message);

    for (const char* f : {"/vsimem/plane_4269.tif", "/vsimem/plane_ft.tif",
                          "/vsimem/plane_out_10.tif", "/vsimem/plane_out_30.tif",
                          "/vsimem/plane_out_ft.tif"})
        VSIUnlink(f);
}

// =============================================================================
// Network mode
// =============================================================================
int run_network_tests(const std::string& output_dir)
{
    if (!lsdtt::http_available())
    {
        std::cout << "SKIP: this GDAL build has no HTTP support (no curl)\n";
        return k_skip_exit;
    }

    const std::string out_path = output_dir + "/fetch_alpine_128.tif";

    lsdtt::lsdtt_fetch_options options;
    options.west  = 10.337; options.south = 46.457;
    options.east  = 10.463; options.north = 46.543;
    options.target_epsg  = 32632;
    options.grid_spacing = 30.0;
    options.size_px      = 128;

    const auto result = lsdtt::fetch_dem(options, out_path);

    if (result.skipped)
    {
        std::cout << "SKIP: " << result.message << "\n";
        return k_skip_exit;
    }

    check("network: fetch succeeded", result.ok, result.message);
    if (!result.ok)
    {
        std::cout << "RESULT: FAIL\n";
        return 1;
    }

    check("network: exact centre crop honoured",
          result.n_cols == 128 && result.n_rows == 128,
          std::to_string(result.n_cols) + "x" + std::to_string(result.n_rows));

    check("network: target EPSG honoured", result.epsg == 32632,
          "got EPSG " + std::to_string(result.epsg));

    check("network: single source tile identified",
          result.tiles.size() == 1 && result.tiles[0] == "N46_00_E010_00",
          join(result.tiles));

    // Open the written raster and verify it is real topography rather than an
    // all-nodata shell, which is what a partially failed range read produces.
    GDALDatasetH ds = GDALOpen(out_path.c_str(), GA_ReadOnly);
    check("network: output opens", ds != nullptr, out_path);
    if (!ds)
    {
        std::cout << "RESULT: FAIL\n";
        return 1;
    }

    GDALRasterBandH band = GDALGetRasterBand(ds, 1);

    int has_nodata = 0;
    const double nodata = GDALGetRasterNoDataValue(band, &has_nodata);
    check("network: nodata value written", has_nodata != 0 && nodata == -9999.0,
          "has_nodata=" + std::to_string(has_nodata) +
          " value=" + std::to_string(nodata));

    std::vector<float> buffer(128 * 128, 0.0f);
    const CPLErr io_err = GDALRasterIO(band, GF_Read, 0, 0, 128, 128,
                                       buffer.data(), 128, 128, GDT_Float32, 0, 0);
    check("network: pixels readable", io_err == CE_None);

    std::int64_t n_valid = 0;
    double min_z =  1e30;
    double max_z = -1e30;
    for (float v : buffer)
    {
        if (v == static_cast<float>(nodata)) continue;
        ++n_valid;
        if (v < min_z) min_z = v;
        if (v > max_z) max_z = v;
    }

    check("network: raster is mostly valid data",
          n_valid > (128 * 128) / 2,
          "valid cells = " + std::to_string(n_valid) + " of " + std::to_string(128 * 128));

    // The Ortler window is high alpine terrain; anything outside this range
    // means the wrong window, wrong units, or a botched warp.
    check("network: elevations are plausible alpine values",
          n_valid > 0 && min_z > 500.0 && max_z < 4500.0,
          "range [" + std::to_string(min_z) + ", " + std::to_string(max_z) + "]");

    GDALClose(ds);

    // ---- the frozen-Copernicus guard ----------------------------------------
    // planetary_computer serves the same GLO-30 COGs, so with the default
    // collection it must write the same raster as cop30_aws. This holds on any
    // GDAL build, which makes it the strongest guard.
    const std::string pc_path = output_dir + "/fetch_alpine_128_pc.tif";
    lsdtt::lsdtt_fetch_options pc_options = options;
    pc_options.source = lsdtt::lsdtt_dem_source::planetary_computer;
    const auto pc_result = lsdtt::fetch_dem(pc_options, pc_path);
    if (pc_result.skipped)
    {
        std::cout << "SKIP(pc): " << pc_result.message << "\n";
    }
    else
    {
        const raster_content cop_c = read_raster_content(out_path);
        const raster_content pc_c  = read_raster_content(pc_path);
        std::string why;
        check("network: cop30_aws and planetary_computer rasters identical",
              pc_result.ok && raster_identical(cop_c, pc_c, why),
              pc_result.ok ? why : pc_result.message);
    }

    // Content hash pinned per GDAL release (the warp output depends on the
    // build). An unknown release reports the hash for adding to the table.
    const std::vector<std::pair<std::string, std::string>> pinned = {
        {"3.4.1", "9f55ea5e672c03d0"},
    };
    const std::string gdal_release = GDALVersionInfo("RELEASE_NAME");
    const std::string hash = raster_content_hash(read_raster_content(out_path));
    bool known = false;
    for (const auto& [release, expected] : pinned)
    {
        if (release != gdal_release) continue;
        known = true;
        check("network: cop30_aws content hash pinned for GDAL " + release,
              hash == expected, "got " + hash);
    }
    if (!known)
    {
        std::cout << "SKIP(pin): no content hash pinned for GDAL " << gdal_release
                  << "; this build gives " << hash << "\n";
    }

    // ---- us_3dep: a Colorado Front Range window ------------------------------
    // No content hash: the NAD83 -> WGS84 transform PROJ picks depends on the
    // grids installed, so sizes, labels and an elevation band are asserted.
    {
        lsdtt::lsdtt_fetch_options dep;
        dep.source = lsdtt::lsdtt_dem_source::us_3dep;
        dep.west = -105.30; dep.south = 39.95; dep.east = -105.25; dep.north = 40.00;
        dep.target_epsg = 32613;
        dep.grid_spacing = 10.0;
        dep.size_px = 256;
        const std::string dep_path = output_dir + "/fetch_3dep_256.tif";
        const auto d = lsdtt::fetch_dem(dep, dep_path);
        if (d.skipped)
        {
            std::cout << "SKIP(us_3dep): " << d.message << "\n";
        }
        else
        {
            check("network: us_3dep fetch succeeded", d.ok, d.message);
            bool only_10m = !d.tiles.empty();
            for (const std::string& t : d.tiles)
                if (t.size() < 3 || t.substr(t.size() - 3) != "-13") only_10m = false;
            check("network: us_3dep used only the 1/3 arc-second (10 m) items", only_10m,
                  join(d.tiles));
            const raster_content c = read_raster_content(dep_path);
            float lo = 1e30f, hi = -1e30f;
            std::int64_t valid = 0;
            for (float v : c.values)
            {
                if (v == -9999.0f) continue;
                ++valid; lo = std::min(lo, v); hi = std::max(hi, v);
            }
            check("network: us_3dep 256 x 256 square 10 m cells",
                  c.ok && c.cols == 256 && c.rows == 256 && c.geotransform[1] == 10.0 &&
                  c.geotransform[5] == -10.0);
            check("network: us_3dep all cells valid, Front Range elevations in metres",
                  valid == 256 * 256 && lo > 1500.0 && hi < 2700.0,
                  "range [" + std::to_string(lo) + ", " + std::to_string(hi) + "]");
            check("network: us_3dep provenance says DTM, NAVD88",
                  metadata_item(dep_path, "LSDTT_SURFACE_MODEL") == "DTM" &&
                  metadata_item(dep_path, "LSDTT_VERTICAL_DATUM") == "NAVD88");
        }
    }
    return 0;
}

} // namespace


int main(int argc, char** argv)
{
    GDALAllRegister();

    bool network_mode = false;
    std::string output_dir = ".";
    std::string data_dir   = "tests/data/bounds";

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if      (arg == "--network") network_mode = true;
        else if (arg == "--offline") network_mode = false;
        else if (arg.rfind("--data-dir=", 0) == 0) data_dir = arg.substr(11);
        else                          output_dir = arg;
    }

    std::cout << "lsdtt3 raster fetch tests "
                 " (docs/developer_notes/testing_and_benchmarking.md) -- "
              << (network_mode ? "network" : "offline") << " mode\n\n";

    if (network_mode)
    {
        const int rc = run_network_tests(output_dir);
        if (rc == k_skip_exit) return rc;
        std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
                  << " checks passed\n";
        std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
        return (g_n_failed == 0) ? 0 : 1;
    }

    test_utm_zone_derivation();
    std::cout << "\n";
    test_tile_naming();
    std::cout << "\n";
    test_source_parsing();
    std::cout << "\n";
    test_opentopography_requires_key();
    std::cout << "\n";
    test_ot_dataset_catalogue();
    std::cout << "\n";
    test_bounds_from_vector(data_dir);
    std::cout << "\n";
    test_option_validation();
    std::cout << "\n";
    test_unreachable_source_skips();
    std::cout << "\n";
    test_no_credential_leakage();
    std::cout << "\n";
    test_copernicus_frozen();
    std::cout << "\n";
    test_dem_source_catalogue();
    std::cout << "\n";
    test_catalogue_preflight();
    std::cout << "\n";
    test_fetch_geometry_helpers();
    std::cout << "\n";
    test_stac_parsing();
    std::cout << "\n";
    test_catalogue_warp_offline();

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return (g_n_failed == 0) ? 0 : 1;
}
