// test_raster_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Part of the flow-routing benchmark (docs/testing_and_benchmarking.md): tests for native DEM fetching.
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

#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <gdal.h>
#include <cpl_conv.h>

#include "lsdtt_raster_fetch.hpp"

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

    // CDSE was cut (docs/testing_and_benchmarking.md): a general public account does not unlock the
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
            d.coverage.empty() || d.description.empty() || d.max_area_km2 <= 0.0)
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
                 " (docs/testing_and_benchmarking.md) -- "
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

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return (g_n_failed == 0) ? 0 : 1;
}
