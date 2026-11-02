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
#include "lsdtt_macrostrat_fetch.hpp"

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
std::string list_number(double v)
{
    std::string text = std::to_string(v);
    while (!text.empty() && text.back() == '0') text.pop_back();
    if (!text.empty() && text.back() == '.') text.pop_back();
    return text;
}

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
        if (!info.xyz_layers.empty())
        {
            bool layers_ok = info.xyz_layers.size() == info.native_resolutions_metres.size();
            for (const auto& l : info.xyz_layers)
                if (l.zoom <= 0 || l.url_template.find("{z}/{x}/{y}") == std::string::npos)
                    layers_ok = false;
            check("catalogue: " + info.name + " has an XYZ layer per resolution", layers_ok);
            continue;
        }
        if (!info.grid_base_url.empty())
        {
            bool templates_ok = !info.grid_collections.empty();
            for (const auto& c : info.grid_collections)
                if (c.key_template.find("{model}") == std::string::npos ||
                    !(c.tile_width > 0.0 && c.tile_height > 0.0))
                    templates_ok = false;
            check("catalogue: " + info.name + " has its bucket access, keyed by model",
                  templates_ok && info.grid_epsg > 0 && !info.grid_sentinel_key.empty());
            continue;
        }
        if (!info.wcs_url.empty())
        {
            bool every_model_served = true;
            for (const std::string& m : info.surface_models)
            {
                bool found = false;
                for (const auto& entry : info.wcs_coverage_by_model)
                    if (entry.first == m && !entry.second.empty()) found = true;
                if (!found) every_model_served = false;
            }
            check("catalogue: " + info.name + " has its WCS access, a coverage per model",
                  every_model_served && info.wcs_native_epsg > 0 &&
                  info.wcs_pixel_metres > 0.0 && info.wcs_max_chunk_side > 0 &&
                  info.native_resolutions_metres ==
                      std::vector<double>{info.wcs_pixel_metres} &&
                  (info.wcs_version == "1.0.0" || info.wcs_version == "2.0.1"));
            continue;
        }
        check("catalogue: " + info.name + " has its STAC access",
              !info.stac_search_url.empty() &&
              (!info.stac_asset.empty() || info.stac_asset_by_gsd));
        bool every_resolution_searchable = true;
        for (double r : info.native_resolutions_metres)
            if (lsdtt::fetch_detail::catalogue_stac_collection(info, r).empty())
                every_resolution_searchable = false;
        check("catalogue: " + info.name + " has a STAC collection for every resolution",
              every_resolution_searchable);
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

// ArcticDEM and REMA: catalogue rows, the PGC search, the pre-flight checks.
void test_polar_sources()
{
    namespace fd = lsdtt::fetch_detail;
    using lsdtt::lsdtt_dem_source;
    const auto& arc  = lsdtt::dem_source_info(lsdtt_dem_source::arctic_arcticdem);
    const auto& rema = lsdtt::dem_source_info(lsdtt_dem_source::antarctic_rema);

    for (const auto* info : {&arc, &rema})
    {
        check("polar: " + info->name + " is a DSM on the WGS84 ellipsoid, said in words",
              info->surface_models == std::vector<std::string>{"dsm"} &&
              info->vertical_datum == "WGS84 ellipsoid" &&
              contains(info->description, "SURFACE MODEL") &&
              contains(info->description, "ELLIPSOID"));
        check("polar: " + info->name + " offers 2, 10, 32 m, default 10",
              info->native_resolutions_metres == std::vector<double>{2.0, 10.0, 32.0} &&
              info->default_resolution_metres == 10.0);
    }
    check("polar: ArcticDEM collections per resolution",
          fd::catalogue_stac_collection(arc, 2.0) == "arcticdem-mosaics-v4.1-2m" &&
          fd::catalogue_stac_collection(arc, 10.0) == "arcticdem-mosaics-v4.1-10m" &&
          fd::catalogue_stac_collection(arc, 32.0) == "arcticdem-mosaics-v4.1-32m");
    check("polar: REMA collections per resolution",
          fd::catalogue_stac_collection(rema, 2.0) == "rema-mosaics-v2.0-2m" &&
          fd::catalogue_stac_collection(rema, 10.0) == "rema-mosaics-v2.0-10m" &&
          fd::catalogue_stac_collection(rema, 32.0) == "rema-mosaics-v2.0-32m");

    lsdtt::lsdtt_fetch_options o;
    o.west = 15.40; o.south = 78.20; o.east = 15.60; o.north = 78.25;
    o.grid_spacing = 10.0;

    // us_3dep's search body is the Planetary Computer body on its collection.
    {
        lsdtt::lsdtt_fetch_options d = o;
        d.source = lsdtt_dem_source::us_3dep;
        lsdtt::lsdtt_fetch_options pc = o;
        pc.pc_collection = "3dep-seamless";
        check("polar: us_3dep search body unchanged by the fields option",
              fd::catalogue_stac_search_body(d) == fd::pc_stac_search_body(pc),
              fd::catalogue_stac_search_body(d));
    }
    o.source = lsdtt_dem_source::arctic_arcticdem;
    o.source_resolution_metres = 2.0;
    const std::string body = fd::catalogue_stac_search_body(o);
    check("polar: ArcticDEM 2 m search body",
          body == "{\"collections\":[\"arcticdem-mosaics-v4.1-2m\"],"
                  "\"bbox\":[15.4000000000,78.2000000000,15.6000000000,78.2500000000],"
                  "\"limit\":100,\"fields\":{\"include\":[\"id\",\"bbox\","
                  "\"properties.gsd\",\"properties.datetime\",\"assets.dem.href\"]}}",
          body);
    o.source_resolution_metres = 0.0;

    // A PGC search response, trimmed (fields filter applied by the server).
    const std::string json = R"({"type": "FeatureCollection",
      "features": [{"id": "34_52_10m_v4.1",
        "bbox": [12.522740954, 77.225277216, 18.44067821, 78.471965713],
        "assets": {"dem": {"href": "https://pgc-opendata-dems.s3.us-west-2.amazonaws.com/arcticdem/mosaics/v4.1/10m/34_52/34_52_10m_v4.1_dem.tif"}},
        "properties": {"gsd": 10, "datetime": "2010-03-29T00:00:00Z"}}],
      "links": [{"rel": "root", "href": "https://stac.pgc.umn.edu/api/v1/"}]})";
    std::vector<fd::stac_item> items;
    std::vector<std::string> missing;
    std::string next_href, next_body, error;
    check("polar: PGC response parsed on the dem asset",
          fd::parse_stac_items(json, "dem", items, missing, next_href, next_body, error) &&
          items.size() == 1 && items[0].gsd == 10.0 && next_href.empty() &&
          contains(items[0].href, "34_52_10m_v4.1_dem.tif"), error);

    auto refused = [](const lsdtt::lsdtt_fetch_result& r, const std::string& part) {
        return !r.ok && !r.skipped && contains(r.message, part);
    };
    lsdtt::lsdtt_fetch_options a = o;
    a.west = 10.3; a.east = 10.4; a.south = 46.4; a.north = 46.5;
    auto r = lsdtt::fetch_dem(a, "/vsimem/never_written.tif");
    check("polar: ArcticDEM over the Alps is refused",
          refused(r, "outside arctic_arcticdem coverage"), r.message);
    a = o;
    a.source = lsdtt_dem_source::antarctic_rema;
    r = lsdtt::fetch_dem(a, "/vsimem/never_written.tif");
    check("polar: REMA over Svalbard is refused",
          refused(r, "outside antarctic_rema coverage"), r.message);
    a = o;
    a.surface_model_choice = "dtm";
    r = lsdtt::fetch_dem(a, "/vsimem/never_written.tif");
    check("polar: ArcticDEM has no terrain model", refused(r, "no dtm product"), r.message);
    a = o;
    a.source_resolution_metres = 5.0;
    r = lsdtt::fetch_dem(a, "/vsimem/never_written.tif");
    check("polar: ArcticDEM resolutions listed", refused(r, "offers 2, 10, 32 m"), r.message);

    // Projected source grid: the download is counted by area. At 78 N a
    // count in degrees would say about 5 MB for this box; by area it is
    // about 1 MB (461 x 562 cells of 10 m, Float32).
    a = o;
    a.max_download_megabytes = 0.5;
    r = lsdtt::fetch_dem(a, "/vsimem/never_written.tif");
    check("polar: download estimated by area on a projected source grid",
          refused(r, "estimated at 1.0 MB"), r.message);
}

// A synthetic EPSG:3413 GeoTIFF over a WGS84 box (10 m cells, 2 km margin)
// holding the plane z = 100 + 0.01 (x - x0) + 0.02 (y - y0).
std::string make_polar_plane_geotiff(const std::string& path, double west, double south,
                                     double east, double north, double& x0, double& y0)
{
    OGRSpatialReferenceH geo = OSRNewSpatialReference(nullptr);
    OGRSpatialReferenceH ps  = OSRNewSpatialReference(nullptr);
    OSRImportFromEPSG(geo, 4326);
    OSRImportFromEPSG(ps, 3413);
    OSRSetAxisMappingStrategy(geo, OAMS_TRADITIONAL_GIS_ORDER);
    OSRSetAxisMappingStrategy(ps, OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformationH ct = OCTNewCoordinateTransformation(geo, ps);
    double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
    for (int j = 0; j <= 8; ++j)
        for (int i = 0; i <= 8; ++i)
        {
            double x = west + i * (east - west) / 8.0, y = south + j * (north - south) / 8.0;
            OCTTransform(ct, 1, &x, &y, nullptr);
            xmin = std::min(xmin, x); xmax = std::max(xmax, x);
            ymin = std::min(ymin, y); ymax = std::max(ymax, y);
        }
    OCTDestroyCoordinateTransformation(ct);
    const double cell = 10.0;
    x0 = std::floor((xmin - 2000.0) / cell) * cell;
    y0 = std::ceil((ymax + 2000.0) / cell) * cell;
    const int nx = static_cast<int>(std::ceil((xmax + 2000.0 - x0) / cell));
    const int ny = static_cast<int>(std::ceil((y0 - (ymin - 2000.0)) / cell));

    GDALDatasetH ds = GDALCreate(GDALGetDriverByName("GTiff"), path.c_str(), nx, ny, 1,
                                 GDT_Float32, nullptr);
    const double gt[6] = {x0, cell, 0.0, y0, 0.0, -cell};
    GDALSetGeoTransform(ds, const_cast<double*>(gt));
    char* wkt = nullptr;
    OSRExportToWkt(ps, &wkt);
    GDALSetProjection(ds, wkt);
    CPLFree(wkt);
    OSRDestroySpatialReference(geo);
    OSRDestroySpatialReference(ps);
    std::vector<float> z(static_cast<std::size_t>(nx) * ny);
    for (int r = 0; r < ny; ++r)
        for (int c = 0; c < nx; ++c)
            z[static_cast<std::size_t>(r) * nx + c] =
                static_cast<float>(100.0 + 0.01 * (c + 0.5) * cell - 0.02 * (r + 0.5) * cell);
    GDALRasterBandH band = GDALGetRasterBand(ds, 1);
    const CPLErr written =
        GDALRasterIO(band, GF_Write, 0, 0, nx, ny, z.data(), nx, ny, GDT_Float32, 0, 0);
    check("polar warp: synthetic EPSG:3413 plane written", written == CE_None, path);
    GDALSetRasterNoDataValue(band, -9999.0);
    GDALClose(ds);
    return path;
}

void test_polar_warp_offline()
{
    namespace fd = lsdtt::fetch_detail;
    lsdtt::lsdtt_fetch_options o;
    o.source = lsdtt::lsdtt_dem_source::arctic_arcticdem;
    o.west = 15.40; o.south = 78.20; o.east = 15.60; o.north = 78.25;
    o.grid_spacing = 10.0;
    double x0 = 0.0, y0 = 0.0;
    const std::string src = make_polar_plane_geotiff("/vsimem/polar_3413.tif", o.west,
                                                     o.south, o.east, o.north, x0, y0);
    const std::string out = "/vsimem/polar_out.tif";
    auto r = fd::warp_catalogue_files(o, {src}, out);
    check("polar warp: EPSG:3413 source written to UTM 33N", r.ok && r.epsg == 32633,
          r.message);
    check("polar warp: no warning inside the UTM latitudes", r.warnings.empty());
    const raster_content c = read_raster_content(out);
    check("polar warp: square 10 m cells",
          c.ok && c.geotransform[1] == 10.0 && c.geotransform[5] == -10.0 &&
          c.geotransform[2] == 0.0 && c.geotransform[4] == 0.0);
    check("polar warp: provenance says DSM, WGS84 ellipsoid",
          metadata_item(out, "LSDTT_DEM_SOURCE") == "arctic_arcticdem" &&
          metadata_item(out, "LSDTT_SURFACE_MODEL") == "DSM" &&
          metadata_item(out, "LSDTT_VERTICAL_DATUM") == "WGS84 ellipsoid" &&
          metadata_item(out, "LSDTT_VALID_PERCENT") == "100.0");
    {
        OGRSpatialReferenceH utm = OSRNewSpatialReference(nullptr);
        OGRSpatialReferenceH ps  = OSRNewSpatialReference(nullptr);
        OSRImportFromEPSG(utm, 32633);
        OSRImportFromEPSG(ps, 3413);
        OSRSetAxisMappingStrategy(utm, OAMS_TRADITIONAL_GIS_ORDER);
        OSRSetAxisMappingStrategy(ps, OAMS_TRADITIONAL_GIS_ORDER);
        OGRCoordinateTransformationH ct = OCTNewCoordinateTransformation(utm, ps);
        double worst = 0.0;
        int n_checked = 0;
        for (int row = 20; row < c.rows - 20; row += 17)
            for (int col = 20; col < c.cols - 20; col += 17)
            {
                double x = c.geotransform[0] + (col + 0.5) * c.geotransform[1];
                double y = c.geotransform[3] + (row + 0.5) * c.geotransform[5];
                OCTTransform(ct, 1, &x, &y, nullptr);
                const double expected = 100.0 + 0.01 * (x - x0) + 0.02 * (y - y0);
                const float got = c.values[static_cast<std::size_t>(row) * c.cols + col];
                worst = std::max(worst, std::abs(got - expected));
                ++n_checked;
            }
        OCTDestroyCoordinateTransformation(ct);
        OSRDestroySpatialReference(utm);
        OSRDestroySpatialReference(ps);
        check("polar warp: a plane in polar stereographic survives the reprojection",
              n_checked > 50 && worst < 1e-2,
              "max |dz| = " + std::to_string(worst) + " over " + std::to_string(n_checked));
    }

    // Beyond 84 N the default UTM output is warned about; target_epsg 3413
    // silences it.
    lsdtt::lsdtt_fetch_options p = o;
    p.west = 10.0; p.east = 12.0; p.south = 85.00; p.north = 85.05;
    const std::string far = make_polar_plane_geotiff("/vsimem/polar_85.tif", p.west,
                                                     p.south, p.east, p.north, x0, y0);
    r = fd::warp_catalogue_files(p, {far}, "/vsimem/polar_85_out.tif");
    check("polar warp: UTM beyond 84 N is warned about",
          r.ok && r.warnings.size() == 1 && contains(r.warnings[0], "target_epsg"),
          r.message);
    p.target_epsg = 3413;
    r = fd::warp_catalogue_files(p, {far}, "/vsimem/polar_85_out.tif");
    check("polar warp: target_epsg 3413 is not warned about",
          r.ok && r.epsg == 3413 && r.warnings.empty(), r.message);

    for (const char* f : {"/vsimem/polar_3413.tif", "/vsimem/polar_out.tif",
                          "/vsimem/polar_85.tif", "/vsimem/polar_85_out.tif"})
        VSIUnlink(f);
}

// GSI elevation tiles: the decoder, a synthetic PNG through the tile step,
// pre-flight refusals.
void test_gsi_tiles()
{
    namespace fd = lsdtt::fetch_detail;
    double h = 0.0;
    check("gsi: (128,0,0) is NoData", !fd::gsi_decode_rgb(128, 0, 0, h));
    check("gsi: (0,0,0) is 0 m", fd::gsi_decode_rgb(0, 0, 0, h) && h == 0.0);
    check("gsi: (0,185,203) is 475.63 m (a real Mt Takao pixel)",
          fd::gsi_decode_rgb(0, 185, 203, h) && std::abs(h - 475.63) < 1e-9, std::to_string(h));
    check("gsi: (255,255,255) is -0.01 m", fd::gsi_decode_rgb(255, 255, 255, h) &&
          std::abs(h + 0.01) < 1e-9, std::to_string(h));
    check("gsi: (127,255,255) is the maximum, 83886.07 m",
          fd::gsi_decode_rgb(127, 255, 255, h) && std::abs(h - 83886.07) < 1e-6);
    check("gsi: (255,254,12) is -5.00 m (x - 2^24 = -500)",
          fd::gsi_decode_rgb(255, 254, 12, h) && std::abs(h + 5.0) < 1e-9, std::to_string(h));

    // A synthetic 256 x 256 RGB PNG: heights encoded by the GSI rule, one
    // NoData pixel, one negative pixel.
    const int n = 256;
    GDALDatasetH mem = GDALCreate(GDALGetDriverByName("MEM"), "", n, n, 3, GDT_Byte, nullptr);
    std::vector<unsigned char> bands[3];
    for (auto& b : bands) b.assign(static_cast<std::size_t>(n) * n, 0);
    auto encode = [&](int idx, double height) {
        long x = std::lround(height * 100.0);
        if (x < 0) x += (1L << 24);
        bands[0][idx] = static_cast<unsigned char>((x >> 16) & 255);
        bands[1][idx] = static_cast<unsigned char>((x >> 8) & 255);
        bands[2][idx] = static_cast<unsigned char>(x & 255);
    };
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) encode(r * n + c, 100.0 + 0.5 * c + 0.25 * r);
    bands[0][0] = 128; bands[1][0] = 0; bands[2][0] = 0;   // NoData
    encode(1, -2.37);
    bool written = true;
    for (int b = 0; b < 3; ++b)
        written = written && GDALRasterIO(GDALGetRasterBand(mem, b + 1), GF_Write, 0, 0, n, n,
                                          bands[b].data(), n, n, GDT_Byte, 0, 0) == CE_None;
    GDALDatasetH png = GDALCreateCopy(GDALGetDriverByName("PNG"), "/vsimem/gsi_tile.png", mem,
                                      FALSE, nullptr, nullptr, nullptr);
    GDALClose(mem);
    if (png) GDALClose(png);
    vsi_l_offset size = 0;
    GByte* data = VSIGetMemFileBuffer("/vsimem/gsi_tile.png", &size, FALSE);
    const std::string bytes = data ? std::string(reinterpret_cast<const char*>(data),
                                                 static_cast<std::size_t>(size)) : "";
    std::string error;
    const bool ok = written && !bytes.empty() &&
                    fd::gsi_tile_to_geotiff(bytes, 15, 29057, 12909, "/vsimem/gsi_tile.tif", error);
    check("gsi: synthetic PNG tile decoded", ok, error);
    if (ok)
    {
        const raster_content c = read_raster_content("/vsimem/gsi_tile.tif");
        lsdtt::lsdtt_slippy_tile t;
        t.z = 15; t.x = 29057; t.y = 12909;
        double x0, y0, x1, y1;
        lsdtt::slippy_tile_bounds_3857(t, x0, y0, x1, y1);
        check("gsi: tile georeferenced to its Web Mercator bounds",
              c.ok && std::abs(c.geotransform[0] - x0) < 1e-6 &&
              std::abs(c.geotransform[3] - y1) < 1e-6 &&
              std::abs(c.geotransform[1] - (x1 - x0) / 256.0) < 1e-9);
        check("gsi: heights decoded, NoData and negative pixels kept",
              c.ok && c.values[0] == -9999.0f && std::abs(c.values[1] + 2.37f) < 1e-4f &&
              std::abs(c.values[n * 10 + 20] - (100.0f + 10.0f + 2.5f)) < 1e-4f);
        GDALDatasetH ds = GDALOpen("/vsimem/gsi_tile.tif", GA_ReadOnly);
        check("gsi: tile CRS is EPSG:3857", ds && contains(GDALGetProjectionRef(ds), "3857"));
        if (ds) GDALClose(ds);
    }
    VSIUnlink("/vsimem/gsi_tile.png");
    VSIUnlink("/vsimem/gsi_tile.tif");

    auto refused = [](const lsdtt::lsdtt_fetch_result& r, const std::string& part) {
        return !r.ok && !r.skipped && contains(r.message, part);
    };
    lsdtt::lsdtt_fetch_options o;
    o.source = lsdtt::lsdtt_dem_source::jp_gsi;
    o.west = 2.30; o.south = 48.80; o.east = 2.35; o.north = 48.85;
    o.grid_spacing = 10.0;
    auto r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("gsi: Paris is refused", refused(r, "outside jp_gsi coverage"), r.message);
    o.west = 139.22; o.south = 35.61; o.east = 139.26; o.north = 35.64;
    o.source_resolution_metres = 2.0;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("gsi: resolutions listed", refused(r, "jp_gsi offers 1, 5, 10 m"), r.message);
    o.source_resolution_metres = 5.0;
    o.grid_spacing = 5.0;
    o.max_tiles = 10;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("gsi: the tile count is refused before any request",
          refused(r, "dem5a_png tiles, more than max_tiles (10)"), r.message);
}

// Grid-keyed buckets: NZTopo50 sheets and OS grid references (truth points
// read from the buckets' own files), key templates, pre-flight refusals.
void test_grid_keyed_sources()
{
    namespace fd = lsdtt::fetch_detail;
    using lsdtt::lsdtt_dem_source;
    // Sheet origins (top-left) read from the LINZ files: AS21 (1492000, 6234000),
    // BQ31 (1732000, 5442000); CC11 and BE36 from their STAC bboxes.
    check("grid: NZTopo50 AS21 at its top-left corner",
          fd::nztopo50_sheet(1492000.5, 6233999.5) == "AS21", fd::nztopo50_sheet(1492000.5, 6233999.5));
    check("grid: NZTopo50 BQ31 (Wellington)", fd::nztopo50_sheet(1744000, 5424000) == "BQ31");
    check("grid: NZTopo50 CC11 (rows skip I and O)", fd::nztopo50_sheet(1264000, 4992000) == "CC11");
    check("grid: NZTopo50 BE36", fd::nztopo50_sheet(1864000, 5784000) == "BE36");
    check("grid: NZTopo50 sheet edges: west of AS21 is AS20, south is AT21",
          fd::nztopo50_sheet(1491999.5, 6233999.5) == "AS20" &&
          fd::nztopo50_sheet(1492000.5, 6197999.5) == "AT21");
    // OS references, from the SRSP files' own origins.
    check("grid: OS 10 km square NH24", fd::os_grid_reference(226500, 848500, 1) == "NH24",
          fd::os_grid_reference(226500, 848500, 1));
    check("grid: OS 1 km square NG2004", fd::os_grid_reference(120500, 804500, 2) == "NG2004",
          fd::os_grid_reference(120500, 804500, 2));
    check("grid: OS 5 km quadrant NS14NE", fd::os_grid_quadrant(217500, 647500) == "NS14NE",
          fd::os_grid_quadrant(217500, 647500));
    check("grid: OS 100 km squares HY, HU, NB", fd::os_grid_reference(320000, 1010000, 0) == "HY" &&
          fd::os_grid_reference(430000, 1130000, 0) == "HU" &&
          fd::os_grid_reference(140000, 930000, 0) == "NB");
    check("grid: outside the National Grid is empty", fd::os_grid_reference(-10, 5, 1).empty());

    const auto& sct = lsdtt::dem_source_info(lsdtt_dem_source::gb_sct_lidar);
    const auto& phase3 = sct.grid_collections[2];
    const auto keys = fd::grid_keys_for_extent(phase3, {214000, 644000, 221000, 651000}, "dtm");
    const std::vector<std::string> want = {
        "NS15SW", "NS15SE", "NS25SW", "NS14NW", "NS14NE", "NS24NW", "NS14SW", "NS14SE", "NS24SW"};
    bool keys_ok = keys.size() == want.size();
    for (std::size_t i = 0; keys_ok && i < want.size(); ++i)
        keys_ok = keys[i] == "lidar/phase-3/dtm/27700/gridded/" + want[i] + "_50CM_DTM_PHASE3.tif";
    check("grid: SRSP phase-3 keys for a 7 km box, row-major from the top left", keys_ok,
          keys.empty() ? "" : keys.front());
    check("grid: SRSP collections oldest first, national programme last",
          sct.grid_collections.front().label == "phase-1" &&
          sct.grid_collections.back().label == "national-lidar-programme");
    const auto nlp = fd::grid_keys_for_extent(sct.grid_collections.back(),
                                              {120200, 804100, 120900, 804900}, "dsm");
    check("grid: SRSP national programme 1 km key, DSM",
          nlp.size() == 1 &&
          nlp[0] == "lidar/national-lidar-programme/dsm/27700/gridded/NG2004_50cm_DSM_ScotlandNationalLiDAR.tif",
          nlp.empty() ? "" : nlp[0]);
    const auto& nz = lsdtt::dem_source_info(lsdtt_dem_source::nz_lidar);
    const auto nzk = fd::grid_keys_for_extent(nz.grid_collections[0],
                                              {1740000, 5420000, 1750000, 5430000}, "dem");
    check("grid: LINZ key for Wellington", nzk.size() == 1 &&
          nzk[0] == "new-zealand/new-zealand/dem_1m/2193/BQ31.tiff", nzk.empty() ? "" : nzk[0]);
    check("grid: LINZ dtm is the bucket's dem folder",
          nz.grid_model_names.size() == 1 && nz.grid_model_names[0].first == "dtm" &&
          nz.grid_model_names[0].second == "dem");

    auto refused = [](const lsdtt::lsdtt_fetch_result& r, const std::string& part) {
        return !r.ok && !r.skipped && contains(r.message, part);
    };
    lsdtt::lsdtt_fetch_options o;
    o.source = lsdtt_dem_source::gb_sct_lidar;
    o.west = 2.30; o.south = 48.80; o.east = 2.35; o.north = 48.85;
    o.grid_spacing = 1.0;
    auto r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("grid: Scottish lidar over Paris is refused",
          refused(r, "outside gb_sct_lidar coverage"), r.message);
    o.west = -4.93; o.south = 55.717; o.east = -4.885; o.north = 55.745;
    o.source_collection_choice = "phase-9";
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("grid: an unknown Scottish collection is refused, listing the valid ones",
          refused(r, "has no collection 'phase-9'") && contains(r.message, "phase-6") &&
          contains(r.message, "national-lidar-programme"), r.message);
    o.source = lsdtt_dem_source::us_3dep;
    o.west = -105.30; o.south = 39.95; o.east = -105.25; o.north = 40.00;
    o.source_collection_choice = "phase-3";
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("grid: source_collection_choice is refused for a source without collections",
          refused(r, "leave it at latest for us_3dep"), r.message);
    o.source_collection_choice = "latest";
    o.source = lsdtt_dem_source::nz_lidar;
    o.west = 174.75; o.south = -41.30; o.east = 174.78; o.north = -41.28;
    o.source_resolution_metres = 2.0;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("grid: LINZ offers 1 m only (checked before the codec)",
          refused(r, "nz_lidar offers 1 m"), r.message);
    const std::string lerc_warning = lsdtt::lerc_source_warning(lsdtt_dem_source::nz_lidar);
    check("grid: nz_lidar screen warning names LERC, this GDAL's status and OpenTopography",
          contains(lerc_warning, "LERC") && contains(lerc_warning, "OpenTopography") &&
          contains(lerc_warning, "API key") &&
          contains(lerc_warning, fd::gdal_has_lerc() ? "has LERC" : "does NOT have LERC"),
          lerc_warning);
    check("grid: no LERC warning for sources that do not need it",
          lsdtt::lerc_source_warning(lsdtt_dem_source::gb_sct_lidar).empty() &&
          lsdtt::lerc_source_warning(lsdtt_dem_source::cop30_aws).empty());
    if (!fd::gdal_has_lerc())
    {
        o.source_resolution_metres = 0.0;
        r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
        check("grid: without LERC in GDAL, nz_lidar is refused before any request, "
              "pointing to OpenTopography",
              refused(r, "LERC") && contains(r.message, "OpenTopography"), r.message);
    }
}

// HRDEM and swissALTI3D: asset selection, latest-per-tile, cursor paging.
void test_hrdem_swiss_sources()
{
    namespace fd = lsdtt::fetch_detail;
    using lsdtt::lsdtt_dem_source;

    // Trimmed from a real swissALTI3D search: assets keyed by file name, two
    // acquisitions of tile 2623-1091, one of 2624-1091.
    const std::string swiss = R"({"type": "FeatureCollection", "features": [
      {"id": "swissalti3d_2019_2623-1091", "bbox": [7.7353597, 45.9701488, 7.7483109, 45.9791788],
       "properties": {"datetime": "2019-01-01T00:00:00Z"},
       "assets": {
         "swissalti3d_2019_2623-1091_0.5_2056_5728.tif": {"href": "https://example.invalid/2019_2623_05.tif", "type": "image/tiff; application=geotiff; profile=cloud-optimized", "gsd": 0.5},
         "swissalti3d_2019_2623-1091_2_2056_5728.tif": {"href": "https://example.invalid/2019_2623_2.tif", "type": "image/tiff; application=geotiff; profile=cloud-optimized", "gsd": 2.0},
         "swissalti3d_2019_2623-1091_2_2056_5728.xyz.zip": {"href": "https://example.invalid/2019_2623_2.xyz.zip", "type": "application/x.ascii-xyz+zip", "gsd": 2.0}}},
      {"id": "swissalti3d_2024_2623-1091", "bbox": [7.7353597, 45.9701488, 7.7483109, 45.9791788],
       "properties": {"datetime": "2024-01-01T00:00:00Z"},
       "assets": {
         "swissalti3d_2024_2623-1091_2_2056_5728.xyz.zip": {"href": "https://example.invalid/2024_2623_2.xyz.zip", "type": "application/x.ascii-xyz+zip", "gsd": 2.0},
         "swissalti3d_2024_2623-1091_2_2056_5728.tif": {"href": "https://example.invalid/2024_2623_2.tif", "type": "image/tiff; application=geotiff; profile=cloud-optimized", "gsd": 2.0}}},
      {"id": "swissalti3d_2019_2624-1091", "bbox": [7.7482601, 45.9701126, 7.7612133, 45.9791441],
       "properties": {"datetime": "2019-01-01T00:00:00Z"},
       "assets": {
         "swissalti3d_2019_2624-1091_2_2056_5728.tif": {"href": "https://example.invalid/2019_2624_2.tif", "type": "image/tiff; application=geotiff; profile=cloud-optimized", "gsd": 2.0}}}],
      "links": [{"rel": "next", "href": "https://example.invalid/search", "method": "POST", "body": {"cursor": "cD01OTM="}}]})";
    std::vector<fd::stac_item> items;
    std::vector<std::string> missing;
    std::string next_href, next_body, error;
    bool ok = fd::parse_stac_items(swiss, "", items, missing, next_href, next_body, error, 2.0);
    bool all_tif = items.size() == 3;
    for (const auto& item : items)
        if (!contains(item.href, "_2.tif") || item.gsd != 2.0) all_tif = false;
    check("swiss: the 2 m GeoTIFF asset is picked by gsd, never the XYZ text", ok && all_tif,
          error);
    items.clear(); missing.clear();
    fd::parse_stac_items(swiss, "", items, missing, next_href, next_body, error, 0.5);
    check("swiss: 0.5 m picks the 0.5 m asset; items without it are named",
          items.size() == 1 && contains(items[0].href, "2019_2623_05.tif") &&
          missing.size() == 2);

    items.clear(); missing.clear();
    fd::parse_stac_items(swiss, "", items, missing, next_href, next_body, error, 2.0);
    const auto latest = fd::latest_stac_item_per_tile(items);
    check("swiss: latest acquisition per tile, in id order",
          latest.size() == 2 && latest[0].id == "swissalti3d_2019_2624-1091" &&
          latest[1].id == "swissalti3d_2024_2623-1091",
          latest.empty() ? "" : latest[0].id);

    // The swisstopo next link carries only a cursor: merged with the
    // original body it keeps the collection and bbox filter.
    const std::string first =
        R"({"collections":["ch.swisstopo.swissalti3d"],"bbox":[7.74,45.97,7.76,45.98],"limit":100})";
    const std::string merged = fd::merge_stac_next_body(first, next_body);
    check("stac paging: a cursor-only next body keeps the original filter",
          contains(merged, "ch.swisstopo.swissalti3d") && contains(merged, "7.74") &&
          contains(merged, "cD01OTM="), merged);
    const std::string full_next =
        R"({"collections":["3dep-seamless"],"bbox":[1,2,3,4],"limit":100,"token":"next:abc"})";
    const std::string merged_full = fd::merge_stac_next_body(
        R"({"collections":["3dep-seamless"],"bbox":[1,2,3,4],"limit":100})", full_next);
    check("stac paging: a full next body is carried unchanged in substance",
          contains(merged_full, "next:abc") && contains(merged_full, "3dep-seamless") &&
          !contains(merged_full, "cursor"), merged_full);

    // A host that ignores the filter: far-away items are dropped.
    fd::stac_item far;
    far.id = "far";
    far.bbox[0] = 8.5; far.bbox[1] = 47.3; far.bbox[2] = 8.6; far.bbox[3] = 47.4;
    items.push_back(far);
    check("stac: items outside the request bbox are dropped",
          fd::stac_items_in_bbox(items, 7.74, 45.97, 7.76, 45.98).size() == items.size() - 1);

    // HRDEM: one item, both models as assets.
    const std::string hrdem = R"({"features": [{"id": "9_2-mosaic-1m",
      "bbox": [-76.986506, 41.266347, -69.184612, 46.888783],
      "properties": {"datetime": "2025-05-05T12:00:00Z", "proj:epsg": 3979},
      "assets": {"dsm": {"href": "https://example.invalid/9_2-mosaic-1m-dsm.tif"},
                 "dtm": {"href": "https://example.invalid/9_2-mosaic-1m-dtm.tif"}}}]})";
    items.clear(); missing.clear();
    const auto& ca = lsdtt::dem_source_info(lsdtt_dem_source::ca_hrdem);
    std::string dsm_key;
    for (const auto& entry : ca.stac_asset_by_model)
        if (entry.first == "dsm") dsm_key = entry.second;
    fd::parse_stac_items(hrdem, dsm_key, items, missing, next_href, next_body, error);
    check("hrdem: the dsm asset is chosen for surface_model_choice dsm",
          items.size() == 1 && contains(items[0].href, "-dsm.tif") && items[0].gsd == 0.0);
    check("hrdem: items carry no gsd, so no gsd filter; collection per resolution",
          !ca.stac_items_carry_gsd &&
          fd::catalogue_stac_collection(ca, 1.0) == "hrdem-mosaic-1m" &&
          fd::catalogue_stac_collection(ca, 2.0) == "hrdem-mosaic-2m");

    auto refused = [](const lsdtt::lsdtt_fetch_result& r, const std::string& part) {
        return !r.ok && !r.skipped && contains(r.message, part);
    };
    lsdtt::lsdtt_fetch_options o;
    o.source = lsdtt_dem_source::ch_swissalti3d;
    o.west = 2.30; o.south = 48.80; o.east = 2.35; o.north = 48.85;
    o.grid_spacing = 2.0;
    auto r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("swiss: Paris is refused", refused(r, "outside ch_swissalti3d coverage"), r.message);
    o.west = 7.74; o.south = 45.97; o.east = 7.76; o.north = 45.98;
    o.surface_model_choice = "dsm";
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("swiss: no surface model", refused(r, "ch_swissalti3d has no dsm product"), r.message);
    o.source = lsdtt_dem_source::ca_hrdem;
    o.west = -75.72; o.south = 45.40; o.east = -75.69; o.north = 45.42;
    o.surface_model_choice = "source_default";
    o.source_resolution_metres = 5.0;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("hrdem: resolutions listed", refused(r, "ca_hrdem offers 1, 2 m"), r.message);
}

// The WCS sources: native extent, chunks, URLs, exception and multipart
// parsing, pre-flight refusals.
void test_wcs_sources()
{
    namespace fd = lsdtt::fetch_detail;
    using lsdtt::lsdtt_dem_source;
    const auto& no = lsdtt::dem_source_info(lsdtt_dem_source::no_dtm);
    const auto& nl = lsdtt::dem_source_info(lsdtt_dem_source::nl_ahn);
    const auto& es = lsdtt::dem_source_info(lsdtt_dem_source::es_mdt);
    const auto& pl = lsdtt::dem_source_info(lsdtt_dem_source::pl_nmt);
    auto on_lattice = [](double v, double origin, double p) {
        const double k = (v - origin) / p;
        return std::abs(k - std::round(k)) < 1e-6;
    };

    // ---- native extent: on the native lattice, covering the output area ----
    {
        lsdtt::lsdtt_fetch_options o;
        o.source = lsdtt_dem_source::pl_nmt;
        o.west = 19.95; o.south = 49.20; o.east = 19.98; o.north = 49.22;
        o.grid_spacing = 1.0;
        const fd::fetch_extent e = fd::wcs_native_extent(o, pl, 32634);
        check("wcs: Polish native extent sits on the grid's pixel edges",
              on_lattice(e.x_min, pl.wcs_origin_x, 1.0) &&
              on_lattice(e.x_max, pl.wcs_origin_x, 1.0) &&
              on_lattice(e.y_min, pl.wcs_origin_y, 1.0) &&
              on_lattice(e.y_max, pl.wcs_origin_y, 1.0),
              std::to_string(e.x_min) + " " + std::to_string(e.y_max));
        // Every corner of the snapped UTM output grid must fall inside it.
        const fd::fetch_extent t =
            fd::target_extent_for_bbox(o.west, o.south, o.east, o.north, 32634, 1.0);
        OGRSpatialReferenceH utm = OSRNewSpatialReference(nullptr);
        OGRSpatialReferenceH pls = OSRNewSpatialReference(nullptr);
        OSRImportFromEPSG(utm, 32634);
        OSRImportFromEPSG(pls, 2180);
        OSRSetAxisMappingStrategy(utm, OAMS_TRADITIONAL_GIS_ORDER);
        OSRSetAxisMappingStrategy(pls, OAMS_TRADITIONAL_GIS_ORDER);
        OGRCoordinateTransformationH ct = OCTNewCoordinateTransformation(utm, pls);
        double xs[4] = {t.x_min, t.x_max, t.x_min, t.x_max};
        double ys[4] = {t.y_min, t.y_min, t.y_max, t.y_max};
        OCTTransform(ct, 4, xs, ys, nullptr);
        OCTDestroyCoordinateTransformation(ct);
        OSRDestroySpatialReference(utm);
        OSRDestroySpatialReference(pls);
        bool inside = true;
        for (int k = 0; k < 4; ++k)
            if (xs[k] < e.x_min + 2.0 || xs[k] > e.x_max - 2.0 ||
                ys[k] < e.y_min + 2.0 || ys[k] > e.y_max - 2.0) inside = false;
        check("wcs: native extent covers the output grid's corners with a margin", inside);
    }

    // ---- chunks: abut exactly, cover exactly, respect max_side ----
    {
        const fd::fetch_extent e{1000.0, 2000.0, 1000.0 + 4500.0, 2000.0 + 2100.0};
        const auto chunks = fd::wcs_plan_chunks(e, 1.0, 2000);
        bool sizes_ok = chunks.size() == 6;
        double area = 0.0;
        for (const auto& c : chunks)
        {
            if (c.x_max - c.x_min > 2000.0 || c.y_max - c.y_min > 2000.0) sizes_ok = false;
            area += (c.x_max - c.x_min) * (c.y_max - c.y_min);
        }
        check("wcs: 4500 x 2100 px at max side 2000 is 3 x 2 chunks", sizes_ok,
              std::to_string(chunks.size()));
        check("wcs: chunks cover the extent exactly, no overlap",
              area == 4500.0 * 2100.0 && chunks.front().x_min == 1000.0 &&
              chunks.front().y_max == 4100.0 && chunks.back().x_max == 5500.0 &&
              chunks.back().y_min == 2000.0);
        check("wcs: neighbouring chunks share an edge",
              chunks.size() == 6 && chunks[0].x_max == chunks[1].x_min &&
              chunks[0].y_min == chunks[3].y_max);
        const auto half = fd::wcs_plan_chunks({0.0, 0.0, 2500.0, 1000.0}, 0.5, 4000);
        check("wcs: 0.5 m pixels count in pixels, not metres",
              half.size() == 2 && half[0].x_max == 2000.0);
    }

    // ---- GetCoverage URLs, pinned ----
    check("wcs: Norway URL (WCS 1.0.0, BBOX + size)",
          fd::wcs_getcoverage_url(no, "nhm_dtm_topo_25833", {140000, 6830000, 140200, 6830100}) ==
          "https://wcs.geonorge.no/skwms1/wcs.hoyde-dtm-nhm-25833?SERVICE=WCS&VERSION=1.0.0"
          "&REQUEST=GetCoverage&COVERAGE=nhm_dtm_topo_25833&CRS=EPSG:25833"
          "&BBOX=140000,6830000,140200,6830100&WIDTH=200&HEIGHT=100&FORMAT=GeoTIFF");
    check("wcs: Netherlands URL (WCS 2.0.1, x/y subsets)",
          fd::wcs_getcoverage_url(nl, "dsm_05m", {136000, 455000, 136100, 455050.5}) ==
          "https://service.pdok.nl/rws/ahn/wcs/v1_0?SERVICE=WCS&VERSION=2.0.1"
          "&REQUEST=GetCoverage&COVERAGEID=dsm_05m&FORMAT=image/tiff"
          "&SUBSET=x(136000,136100)&SUBSET=y(455000,455050.5)");
    check("wcs: Spain URL asks for the text grid",
          fd::wcs_getcoverage_url(es, "Elevacion25830_5", {420002.5, 4520002.5, 420102.5, 4520102.5}) ==
          "https://servicios.idee.es/wcs-inspire/mdt?SERVICE=WCS&VERSION=2.0.1"
          "&REQUEST=GetCoverage&COVERAGEID=Elevacion25830_5&FORMAT=application/asc"
          "&SUBSET=x(420002.5,420102.5)&SUBSET=y(4520002.5,4520102.5)");
    check("wcs: Poland URL keeps the grid's fractional origin",
          fd::wcs_getcoverage_url(pl, "DTM_PL-KRON86-NH_TIFF",
                                  {637000.343266, 486000.66941, 637100.343266, 486100.66941}) ==
          "https://mapy.geoportal.gov.pl/wss/service/PZGIK/NMT/GRID1/WCS/"
          "DigitalTerrainModelFormatTIFF?SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage"
          "&COVERAGEID=DTM_PL-KRON86-NH_TIFF&FORMAT=image/tiff"
          "&SUBSET=x(637000.343266,637100.343266)&SUBSET=y(486000.66941,486100.66941)");

    // ---- exception reports (as returned by the services, trimmed) ----
    const std::string ows = R"(<?xml version="1.0" encoding="UTF-8"?>
<ows:ExceptionReport xmlns:ows="http://www.opengis.net/ows/2.0" version="2.0.1" xml:lang="en-US">
  <ows:Exception exceptionCode="InvalidParameterValue" locator="size">
    <ows:ExceptionText>msWCSGetCoverage20(): WCS server error. Raster size out of range, width and height of resulting coverage must be no more than MAXSIZE=4000.</ows:ExceptionText>
  </ows:Exception>
</ows:ExceptionReport>)";
    const std::string wcs10 = R"(<?xml version="1.0" encoding="UTF-8" ?>
<ServiceExceptionReport version="1.2.0" xmlns="http://www.opengis.net/ogc">
  <ServiceException code="UnSupportedCombination">
Request image size exceeds allowed limit.
  </ServiceException>
</ServiceExceptionReport>)";
    check("wcs: OWS 2.0 exception text read",
          contains(fd::parse_ows_exception(ows), "MAXSIZE=4000"), fd::parse_ows_exception(ows));
    check("wcs: WCS 1.0 service exception read",
          fd::parse_ows_exception(wcs10) == "Request image size exceeds allowed limit.",
          fd::parse_ows_exception(wcs10));
    check("wcs: a GeoTIFF is not an exception", fd::parse_ows_exception("II*\0abc").empty());

    // ---- multipart (the Spanish service wraps its ASCII grid) ----
    const std::string multipart =
        "\r\n--wcs\r\nContent-Type: application/asc\r\nContent-ID: coverage/out.asc\r\n\r\n"
        "ncols 2\nnrows 1\nxllcorner 0\nyllcorner 0\ncellsize 5\n 1.5 2.5\n\r\n--wcs--\r\n";
    check("wcs: multipart payload extracted",
          fd::multipart_first_payload(multipart, "wcs") ==
          "ncols 2\nnrows 1\nxllcorner 0\nyllcorner 0\ncellsize 5\n 1.5 2.5",
          fd::multipart_first_payload(multipart, "wcs"));
    const std::string with_xml =
        "--b\r\nContent-Type: text/xml\r\n\r\n<x/>\r\n--b\r\nContent-Type: image/tiff\r\n\r\nDATA\r\n--b--";
    check("wcs: multipart skips an XML part",
          fd::multipart_first_payload(with_xml, "b") == "DATA");

    // ---- pre-flight ----
    auto refused = [](const lsdtt::lsdtt_fetch_result& r, const std::string& part) {
        return !r.ok && !r.skipped && contains(r.message, part);
    };
    lsdtt::lsdtt_fetch_options o;
    o.source = lsdtt_dem_source::no_dtm;
    o.west = 10.3; o.east = 10.4; o.south = 46.4; o.north = 46.5;
    o.grid_spacing = 1.0;
    auto r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("wcs: Norway over the Alps is refused", refused(r, "outside no_dtm coverage"),
          r.message);
    o.source = lsdtt_dem_source::es_mdt;
    o.west = -3.98; o.south = 40.78; o.east = -3.95; o.north = 40.80;
    o.surface_model_choice = "dsm";
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("wcs: Spain has no surface model", refused(r, "es_mdt has no dsm product"), r.message);
    o.source = lsdtt_dem_source::nl_ahn;
    o.west = 5.10; o.south = 52.08; o.east = 5.12; o.north = 52.09;
    o.surface_model_choice = "source_default";
    o.source_resolution_metres = 5.0;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("wcs: the AHN service has only 0.5 m", refused(r, "nl_ahn offers 0.5 m"), r.message);
    // A 10 m output from AHN still downloads 0.5 m cells: the limits count those.
    o.source_resolution_metres = 0.0;
    o.grid_spacing = 10.0;
    o.max_output_pixels = 1000000;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("wcs: max_output_pixels counts the native 0.5 m cells downloaded",
          refused(r, "native 0.5 m cells"), r.message);
    o.max_output_pixels = 250000000;
    o.max_download_megabytes = 5.0;
    r = lsdtt::fetch_dem(o, "/vsimem/never_written.tif");
    check("wcs: max_download_megabytes counts the native cells downloaded",
          refused(r, "max_download_megabytes"), r.message);
}

// Write a small GeoTIFF on a given grid and return its bytes (a stand-in
// for a GetCoverage response).
std::string geotiff_bytes(double x0, double y0, double pixel, int w, int h, float value)
{
    const std::string path = "/vsimem/wcs_response.tif";
    GDALDatasetH ds = GDALCreate(GDALGetDriverByName("GTiff"), path.c_str(), w, h, 1,
                                 GDT_Float32, nullptr);
    const double gt[6] = {x0, pixel, 0.0, y0, 0.0, -pixel};
    GDALSetGeoTransform(ds, const_cast<double*>(gt));
    std::vector<float> z(static_cast<std::size_t>(w) * h, value);
    z[0] = 0.0f;   // one "outside the country" cell
    const CPLErr written = GDALRasterIO(GDALGetRasterBand(ds, 1), GF_Write, 0, 0, w, h,
                                        z.data(), w, h, GDT_Float32, 0, 0);
    check("wcs staging: synthetic response written", written == CE_None);
    GDALClose(ds);
    vsi_l_offset size = 0;
    GByte* data = VSIGetMemFileBuffer(path.c_str(), &size, FALSE);
    const std::string bytes(reinterpret_cast<const char*>(data), static_cast<std::size_t>(size));
    VSIUnlink(path.c_str());
    return bytes;
}

void test_wcs_staging_offline()
{
    namespace fd = lsdtt::fetch_detail;
    const auto& pl = lsdtt::dem_source_info(lsdtt::lsdtt_dem_source::pl_nmt);
    const fd::fetch_extent chunk{637000.343266, 486000.66941, 637100.343266, 486050.66941};
    std::string error;

    const std::string good = geotiff_bytes(chunk.x_min, chunk.y_max, 1.0, 100, 50, 120.0f);
    bool ok = fd::wcs_stage_response(pl, good, "image/tiff", chunk, "/vsimem/wcs_ok.tif", error);
    GDALDatasetH ds = ok ? GDALOpen("/vsimem/wcs_ok.tif", GA_ReadOnly) : nullptr;
    int has_nodata = 0;
    const double nodata = ds ? GDALGetRasterNoDataValue(GDALGetRasterBand(ds, 1), &has_nodata) : -1;
    const char* wkt = ds ? GDALGetProjectionRef(ds) : "";
    check("wcs staging: an exact native chunk is staged with its CRS",
          ok && ds && contains(wkt, "2180") &&
          GDALGetRasterDataType(GDALGetRasterBand(ds, 1)) == GDT_Float32, error);
    check("wcs staging: 0 is nodata where the service returns 0 outside the country",
          has_nodata && nodata == 0.0);
    if (ds) GDALClose(ds);

    const std::string shifted = geotiff_bytes(chunk.x_min + 0.5, chunk.y_max, 1.0, 100, 50, 120.0f);
    ok = fd::wcs_stage_response(pl, shifted, "image/tiff", chunk, "/vsimem/wcs_bad.tif", error);
    check("wcs staging: a half-pixel shift (server resampling) is refused",
          !ok && contains(error, "not the requested native"), error);
    const std::string coarse = geotiff_bytes(chunk.x_min, chunk.y_max, 2.0, 50, 25, 120.0f);
    ok = fd::wcs_stage_response(pl, coarse, "image/tiff", chunk, "/vsimem/wcs_bad.tif", error);
    check("wcs staging: a coarser grid is refused", !ok && contains(error, "not the requested native"),
          error);

    // Spain: multipart ESRI ASCII, millimetre heights, no CRS in the payload.
    const auto& es = lsdtt::dem_source_info(lsdtt::lsdtt_dem_source::es_mdt);
    const fd::fetch_extent es_chunk{420002.5, 4520002.5, 420012.5, 4520012.5};
    const std::string body =
        "\r\n--wcs\r\nContent-Type: application/asc\r\nContent-ID: coverage/out.asc\r\n\r\n"
        "ncols        2\nnrows        2\nxllcorner    420002.500000000000\n"
        "yllcorner    4520002.500000000000\ncellsize     5.000000000000\n"
        " 1740.251 1739.792\n 1741.003 0.000\n\r\n--wcs--\r\n";
    ok = fd::wcs_stage_response(es, body, "multipart/related; boundary=wcs", es_chunk,
                                "/vsimem/wcs_es.tif", error);
    float v[4] = {0, 0, 0, 0};
    ds = ok ? GDALOpen("/vsimem/wcs_es.tif", GA_ReadOnly) : nullptr;
    CPLErr read = CE_Failure;
    if (ds) read = GDALRasterIO(GDALGetRasterBand(ds, 1), GF_Read, 0, 0, 2, 2, v, 2, 2,
                                GDT_Float32, 0, 0);
    check("wcs staging: Spanish text grid staged with millimetre heights and EPSG:25830",
          ok && ds && read == CE_None && contains(GDALGetProjectionRef(ds), "25830") &&
          std::abs(v[0] - 1740.251f) < 1e-3f && std::abs(v[2] - 1741.003f) < 1e-3f, error);
    if (ds) GDALClose(ds);

    for (const char* f : {"/vsimem/wcs_ok.tif", "/vsimem/wcs_bad.tif", "/vsimem/wcs_es.tif"})
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

    // ---- ArcticDEM (Svalbard) and REMA (Antarctic Peninsula), 10 m ----------
    struct polar_case
    {
        lsdtt::lsdtt_dem_source source;
        const char* label;
        double west, south, east, north;
        int epsg;
        const char* tile_suffix;
        float lo, hi;   // plausible ellipsoidal heights for the window
    };
    const polar_case cases[] = {
        {lsdtt::lsdtt_dem_source::arctic_arcticdem, "arctic_arcticdem",
         15.40, 78.20, 15.60, 78.25, 32633, "_10m_v4.1", 0.0f, 1200.0f},
        {lsdtt::lsdtt_dem_source::antarctic_rema, "antarctic_rema",
         -63.80, -64.85, -63.60, -64.75, 32720, "_10m_v2.0", 0.0f, 2500.0f}};
    for (const polar_case& pc : cases)
    {
        lsdtt::lsdtt_fetch_options po;
        po.source = pc.source;
        po.west = pc.west; po.south = pc.south; po.east = pc.east; po.north = pc.north;
        po.grid_spacing = 10.0;
        po.size_px = 256;
        const std::string path = output_dir + "/fetch_" + pc.label + "_256.tif";
        const auto p = lsdtt::fetch_dem(po, path);
        const std::string name = pc.label;
        if (p.skipped)
        {
            std::cout << "SKIP(" << name << "): " << p.message << "\n";
            continue;
        }
        check("network: " + name + " fetch succeeded", p.ok, p.message);
        bool right_tiles = !p.tiles.empty();
        for (const std::string& t : p.tiles)
            if (!contains(t, pc.tile_suffix)) right_tiles = false;
        check("network: " + name + " used the 10 m mosaic tiles, UTM output",
              right_tiles && p.epsg == pc.epsg, join(p.tiles));
        const raster_content c = read_raster_content(path);
        float lo = 1e30f, hi = -1e30f;
        std::int64_t valid = 0;
        for (float v : c.values)
        {
            if (v == -9999.0f) continue;
            ++valid; lo = std::min(lo, v); hi = std::max(hi, v);
        }
        check("network: " + name + " 256 x 256 square 10 m cells, all valid",
              c.ok && c.cols == 256 && c.rows == 256 && c.geotransform[1] == 10.0 &&
              c.geotransform[5] == -10.0 && valid == 256 * 256);
        check("network: " + name + " heights in a plausible band (metres)",
              lo > pc.lo && hi < pc.hi,
              "range [" + std::to_string(lo) + ", " + std::to_string(hi) + "]");
        check("network: " + name + " provenance says DSM, WGS84 ellipsoid",
              metadata_item(path, "LSDTT_SURFACE_MODEL") == "DSM" &&
              metadata_item(path, "LSDTT_VERTICAL_DATUM") == "WGS84 ellipsoid");
    }

    // ---- HRDEM (Ottawa) and swissALTI3D (Zermatt) ---------------------------
    struct stac_case
    {
        lsdtt::lsdtt_dem_source source;
        const char* label;
        double west, south, east, north;
        double grid;
        int epsg;
        const char* datum;
        float lo, hi;
    };
    const stac_case stac_cases[] = {
        {lsdtt::lsdtt_dem_source::ca_hrdem, "ca_hrdem", -75.72, 45.40, -75.69, 45.42, 1.0,
         32618, "CGVD2013", 30.0f, 120.0f},
        {lsdtt::lsdtt_dem_source::ch_swissalti3d, "ch_swissalti3d", 7.74, 45.97, 7.76, 45.98,
         2.0, 32632, "LN02", 2000.0f, 3000.0f}};
    for (const stac_case& sc : stac_cases)
    {
        lsdtt::lsdtt_fetch_options so;
        so.source = sc.source;
        so.west = sc.west; so.south = sc.south; so.east = sc.east; so.north = sc.north;
        so.grid_spacing = sc.grid;
        so.size_px = 256;
        const std::string path = output_dir + "/fetch_" + sc.label + "_256.tif";
        const auto f = lsdtt::fetch_dem(so, path);
        const std::string name = sc.label;
        if (f.skipped)
        {
            std::cout << "SKIP(" << name << "): " << f.message << "\n";
            continue;
        }
        check("network: " + name + " fetch succeeded", f.ok, f.message);
        const raster_content c = read_raster_content(path);
        float lo = 1e30f, hi = -1e30f;
        std::int64_t valid = 0;
        for (float v : c.values)
        {
            if (v == -9999.0f) continue;
            ++valid; lo = std::min(lo, v); hi = std::max(hi, v);
        }
        check("network: " + name + " 256 x 256 square cells, all valid, plausible heights",
              c.ok && c.cols == 256 && c.rows == 256 && c.geotransform[1] == sc.grid &&
              f.epsg == sc.epsg && valid == 256 * 256 && lo > sc.lo && hi < sc.hi,
              "valid " + std::to_string(valid) + ", range [" + std::to_string(lo) + ", " +
              std::to_string(hi) + "]");
        check("network: " + name + " provenance says DTM and its datum",
              metadata_item(path, "LSDTT_SURFACE_MODEL") == "DTM" &&
              metadata_item(path, "LSDTT_VERTICAL_DATUM") == sc.datum);
        if (sc.source == lsdtt::lsdtt_dem_source::ch_swissalti3d)
        {
            std::set<std::string> tiles;
            for (const std::string& t : f.tiles) tiles.insert(t.substr(t.rfind('_') + 1));
            check("network: swissALTI3D one acquisition per tile",
                  !f.tiles.empty() && tiles.size() == f.tiles.size(), join(f.tiles));
        }
    }

    // ---- Grid-keyed buckets: Scotland (phase 1, inland) and New Zealand -------
    {
        lsdtt::lsdtt_fetch_options go;
        go.source = lsdtt::lsdtt_dem_source::gb_sct_lidar;
        go.west = -4.894; go.south = 57.471; go.east = -4.855; go.north = 57.494;
        go.grid_spacing = 1.0;
        go.size_px = 256;
        const std::string path = output_dir + "/fetch_gb_sct_lidar_256.tif";
        const auto g = lsdtt::fetch_dem(go, path);
        if (g.skipped)
        {
            std::cout << "SKIP(gb_sct_lidar): " << g.message << "\n";
        }
        else
        {
            check("network: gb_sct_lidar fetch succeeded", g.ok, g.message);
            const raster_content c = read_raster_content(path);
            float lo = 1e30f, hi = -1e30f;
            std::int64_t valid = 0;
            for (float v : c.values)
            {
                if (v == -9999.0f) continue;
                ++valid; lo = std::min(lo, v); hi = std::max(hi, v);
            }
            check("network: gb_sct_lidar 256 x 256 at 1 m, the phase-1 tile, plausible heights",
                  c.ok && c.cols == 256 && c.rows == 256 && c.geotransform[1] == 1.0 &&
                  g.tiles.size() == 1 && g.tiles[0] == "NH24_1M_DTM_PHASE1.tif" &&
                  valid > 0 && lo > 150.0f && hi < 900.0f,
                  join(g.tiles) + " valid " + std::to_string(valid) + " [" +
                  std::to_string(lo) + ", " + std::to_string(hi) + "]");
            check("network: gb_sct_lidar provenance says DTM, ODN, phase-1",
                  metadata_item(path, "LSDTT_VERTICAL_DATUM") == "ODN" &&
                  metadata_item(path, "LSDTT_SOURCE_COLLECTIONS") == "phase-1");
        }
        // Clyde coast: phases 1, 2, 4 and 6 all have tiles; latest takes
        // phase-6 alone, a named phase that phase alone.
        for (const char* choice : {"latest", "phase-4"})
        {
            lsdtt::lsdtt_fetch_options co;
            co.source = lsdtt::lsdtt_dem_source::gb_sct_lidar;
            co.west = -4.93; co.south = 55.717; co.east = -4.885; co.north = 55.745;
            co.grid_spacing = 10.0;
            co.source_collection_choice = choice;
            const std::string cpath = output_dir + "/fetch_gb_sct_" + choice + ".tif";
            const auto cr = lsdtt::fetch_dem(co, cpath);
            if (cr.skipped) continue;
            const std::string want = std::string(choice) == "latest" ? "phase-6" : "phase-4";
            bool one_collection = !cr.tiles.empty();
            for (const std::string& t : cr.tiles)
                if (!contains(t, want == "phase-6" ? "PHASE6" : "PHASE4")) one_collection = false;
            check(std::string("network: gb_sct_lidar ") + choice + " uses " + want + " only",
                  cr.ok && one_collection &&
                  metadata_item(cpath, "LSDTT_SOURCE_COLLECTIONS") == want,
                  cr.message + " " + join(cr.tiles));
        }
    }
    if (!lsdtt::fetch_detail::gdal_has_lerc())
    {
        std::cout << "SKIP(nz_lidar): this GDAL build has no LERC codec\n";
    }
    else
    {
        lsdtt::lsdtt_fetch_options no;
        no.source = lsdtt::lsdtt_dem_source::nz_lidar;
        no.west = 174.75; no.south = -41.30; no.east = 174.78; no.north = -41.28;
        no.grid_spacing = 1.0;
        no.size_px = 256;
        const std::string path = output_dir + "/fetch_nz_lidar_256.tif";
        const auto n = lsdtt::fetch_dem(no, path);
        if (n.skipped)
        {
            std::cout << "SKIP(nz_lidar): " << n.message << "\n";
        }
        else
        {
            check("network: nz_lidar fetch succeeded, sheet BQ31", n.ok && n.tiles.size() == 1 &&
                  n.tiles[0] == "BQ31.tiff", n.message + " " + join(n.tiles));
            check("network: nz_lidar provenance says DTM, NZVD2016",
                  metadata_item(path, "LSDTT_SURFACE_MODEL") == "DTM" &&
                  metadata_item(path, "LSDTT_VERTICAL_DATUM") == "NZVD2016");
        }
    }

    // ---- GSI tiles near Tokyo (Mt Takao), 10 m and 5 m ------------------------
    for (double gsi_res : {10.0, 5.0})
    {
        lsdtt::lsdtt_fetch_options jo;
        jo.source = lsdtt::lsdtt_dem_source::jp_gsi;
        jo.west = 139.22; jo.south = 35.61; jo.east = 139.26; jo.north = 35.64;
        jo.source_resolution_metres = gsi_res;
        jo.grid_spacing = gsi_res;
        jo.size_px = 256;
        const std::string label = "jp_gsi " + list_number(gsi_res) + " m";
        const std::string path = output_dir + "/fetch_jp_gsi_" + list_number(gsi_res) + ".tif";
        const auto j = lsdtt::fetch_dem(jo, path);
        if (j.skipped)
        {
            std::cout << "SKIP(" << label << "): " << j.message << "\n";
            continue;
        }
        check("network: " + label + " fetch succeeded", j.ok, j.message);
        bool right_layer = !j.tiles.empty();
        for (const std::string& t : j.tiles)
            if (t.rfind(gsi_res == 10.0 ? "dem_png/14/" : "dem5a_png/15/", 0) != 0)
                right_layer = false;
        const raster_content c = read_raster_content(path);
        float lo = 1e30f, hi = -1e30f;
        std::int64_t valid = 0;
        for (float v : c.values)
        {
            if (v == -9999.0f) continue;
            ++valid; lo = std::min(lo, v); hi = std::max(hi, v);
        }
        check("network: " + label + " used only its layer, 256 x 256, Takao heights",
              right_layer && c.ok && c.cols == 256 && c.rows == 256 && j.epsg == 32654 &&
              valid > 0.99 * 256 * 256 && lo > 100.0f && hi < 700.0f,
              join(j.tiles) + " [" + std::to_string(lo) + ", " + std::to_string(hi) + "]");
    }

    // ---- WCS sources: one small window each ---------------------------------
    struct wcs_case
    {
        lsdtt::lsdtt_dem_source source;
        const char* label;
        double west, south, east, north;
        double grid;
        int epsg;
        const char* datum;
        float lo, hi;          // plausible heights for the window
        double min_valid;      // percent of the 256 x 256 window
    };
    const wcs_case wcs_cases[] = {
        {lsdtt::lsdtt_dem_source::no_dtm, "no_dtm", 8.30, 61.60, 8.33, 61.62, 1.0,
         32632, "NN2000", 1300.0f, 2500.0f, 99.0},
        {lsdtt::lsdtt_dem_source::nl_ahn, "nl_ahn", 5.10, 52.08, 5.12, 52.09, 0.5,
         32631, "NAP", -5.0f, 30.0f, 20.0},
        {lsdtt::lsdtt_dem_source::es_mdt, "es_mdt", -3.98, 40.78, -3.95, 40.80, 5.0,
         32630, "REDNAP (Alicante)", 1700.0f, 2500.0f, 99.0},
        {lsdtt::lsdtt_dem_source::pl_nmt, "pl_nmt", 19.97, 49.24, 19.98, 49.25, 1.0,
         32634, "PL-KRON86-NH", 800.0f, 2000.0f, 99.0}};
    for (const wcs_case& wc : wcs_cases)
    {
        lsdtt::lsdtt_fetch_options wo;
        wo.source = wc.source;
        wo.west = wc.west; wo.south = wc.south; wo.east = wc.east; wo.north = wc.north;
        wo.grid_spacing = wc.grid;
        wo.size_px = 256;
        const std::string path = output_dir + "/fetch_" + wc.label + "_256.tif";
        const auto w = lsdtt::fetch_dem(wo, path);
        const std::string name = wc.label;
        if (w.skipped)
        {
            std::cout << "SKIP(" << name << "): " << w.message << "\n";
            continue;
        }
        check("network: " + name + " fetch succeeded", w.ok, w.message);
        const raster_content c = read_raster_content(path);
        float lo = 1e30f, hi = -1e30f;
        std::int64_t valid = 0;
        for (float v : c.values)
        {
            if (v == -9999.0f) continue;
            ++valid; lo = std::min(lo, v); hi = std::max(hi, v);
        }
        check("network: " + name + " 256 x 256 square cells at the native spacing, UTM",
              c.ok && c.cols == 256 && c.rows == 256 && c.geotransform[1] == wc.grid &&
              c.geotransform[5] == -wc.grid && w.epsg == wc.epsg);
        check("network: " + name + " valid cells and heights in a plausible band (metres)",
              100.0 * valid / (256.0 * 256.0) >= wc.min_valid && lo > wc.lo && hi < wc.hi,
              "valid " + std::to_string(valid) + ", range [" + std::to_string(lo) + ", " +
              std::to_string(hi) + "]");
        check("network: " + name + " provenance says DTM and its datum",
              metadata_item(path, "LSDTT_SURFACE_MODEL") == "DTM" &&
              metadata_item(path, "LSDTT_VERTICAL_DATUM") == wc.datum &&
              !metadata_item(path, "LSDTT_DOWNLOAD_MEGABYTES").empty());
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
    std::cout << "\n";
    test_polar_sources();
    std::cout << "\n";
    test_polar_warp_offline();
    std::cout << "\n";
    test_hrdem_swiss_sources();
    std::cout << "\n";
    test_grid_keyed_sources();
    std::cout << "\n";
    test_gsi_tiles();
    std::cout << "\n";
    test_wcs_sources();
    std::cout << "\n";
    test_wcs_staging_offline();

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return (g_n_failed == 0) ? 0 : 1;
}
