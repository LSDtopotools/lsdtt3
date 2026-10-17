// test_macrostrat_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Tests for the Macrostrat fetch (lsdtt-fetch-lithology's library), offline.
//
//   * slippy-tile maths against hand-computed tiles;
//   * the JSON parsers, including malformed and error responses;
//   * the lith_primary rule (major -> liths -> text -> none);
//   * a fetch from the synthetic tile fixture (tests/data/macrostrat_tiles):
//     units outside the box dropped, reprojection to the local UTM zone, areas
//     against the source polygon projected independently, the vocabulary;
//   * graceful degradation: an unreachable API is a warning, not a failure;
//     too many tiles and no readable tile are failures.
//
// Usage: lsdtt-test-macrostrat-fetch <fixture_tile_dir> <work_dir>
// Exits 77 (CTest skip) when this GDAL has no MVT driver.
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_http_fetch.hpp"
#include "lsdtt_macrostrat_fetch.hpp"

namespace
{

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

template <typename F>
std::string thrown_message(F f)
{
    try { f(); }
    catch (const std::exception& e) { return e.what(); }
    return "";
}

bool contains(const std::string& s, const std::string& part)
{
    return s.find(part) != std::string::npos;
}

// Area (m2) of a WGS84 lon/lat rectangle projected to an EPSG code, densified
// so the projected edges follow the curved ones.
double projected_rectangle_area(double w, double s, double e, double n, int epsg)
{
    OGRSpatialReference ll, out;
    ll.importFromEPSG(4326);
    out.importFromEPSG(epsg);
    ll.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    out.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    std::unique_ptr<OGRCoordinateTransformation> ct(OGRCreateCoordinateTransformation(&ll, &out));
    OGRLinearRing ring;
    const int k = 50;
    for (int i = 0; i < k; ++i) ring.addPoint(w + (e - w) * i / k, s);
    for (int i = 0; i < k; ++i) ring.addPoint(e, s + (n - s) * i / k);
    for (int i = 0; i < k; ++i) ring.addPoint(e - (e - w) * i / k, n);
    for (int i = 0; i < k; ++i) ring.addPoint(w, n - (n - s) * i / k);
    ring.closeRings();
    OGRPolygon poly;
    poly.addRing(&ring);
    poly.transform(ct.get());
    return poly.get_Area();
}

void test_tiles()
{
    using lsdtt::slippy_tiles_for_bbox;
    auto t0 = slippy_tiles_for_bbox(-180.0, -85.0, 180.0, 85.0, 0);
    check("zoom 0: one tile", t0.size() == 1 && t0[0].x == 0 && t0[0].y == 0);

    // Berlin (13.40 E, 52.52 N) is tile 550/335 at zoom 10 (by hand from the
    // slippy-map formulae).
    auto tb = slippy_tiles_for_bbox(13.40, 52.52, 13.401, 52.521, 10);
    check("Berlin zoom 10 = 550/335", tb.size() == 1 && tb[0].x == 550 && tb[0].y == 335,
          tb.empty() ? "none" : std::to_string(tb[0].x) + "/" + std::to_string(tb[0].y));

    // The fixture box sits in one tile, 540/364.
    auto tf = slippy_tiles_for_bbox(10.00, 46.00, 10.08, 46.05, 10);
    check("fixture box = tile 540/364", tf.size() == 1 && tf[0].x == 540 && tf[0].y == 364);

    // Rows run north to south, columns west to east.
    auto t4 = slippy_tiles_for_bbox(-1.0, -1.0, 1.0, 1.0, 1);
    check("2x2 around (0,0) at zoom 1, NW first",
          t4.size() == 4 && t4[0].x == 0 && t4[0].y == 0 && t4[3].x == 1 && t4[3].y == 1);

    // Latitudes beyond Web Mercator's limit are clamped, not refused.
    auto tp = slippy_tiles_for_bbox(0.0, 80.0, 1.0, 90.0, 2);
    check("polar box clamped to the top row", tp.size() == 1 && tp[0].y == 0);

    double x0, y0, x1, y1;
    lsdtt::slippy_tile_bounds_3857({0, 0, 0}, x0, y0, x1, y1);
    check("zoom 0 bounds = +-20037508.34 m",
          std::abs(x0 + 20037508.342789) < 1e-3 && std::abs(y1 - 20037508.342789) < 1e-3);
    lsdtt::slippy_tile_bounds_3857({10, 540, 364}, x0, y0, x1, y1);
    check("tile 10/540/364 is 39135.76 m wide",
          std::abs((x1 - x0) - 39135.758482) < 1e-3 && std::abs((y1 - y0) - 39135.758482) < 1e-3);

    check("west >= east refused",
          contains(thrown_message([] { (void)slippy_tiles_for_bbox(2.0, 0.0, 1.0, 1.0, 5); }),
                   "west < east"));
    check("longitude beyond 180 refused",
          !thrown_message([] { (void)slippy_tiles_for_bbox(170.0, 0.0, 190.0, 1.0, 5); }).empty());
    check("zoom 23 refused",
          !thrown_message([] { (void)slippy_tiles_for_bbox(0.0, 0.0, 1.0, 1.0, 23); }).empty());

    check("UTM: Edinburgh -> 32630", lsdtt::utm_epsg_for(-3.2, 55.9) == 32630);
    check("UTM: Canterbury NZ -> 32759", lsdtt::utm_epsg_for(172.0, -43.5) == 32759);
    check("UTM: equator is north", lsdtt::utm_epsg_for(0.5, 0.0) == 32631);
    check("UTM: 180 E -> zone 60", lsdtt::utm_epsg_for(180.0, 10.0) == 32660);
}

void test_parsers()
{
    const std::string liths =
        R"({"success":{"v":2,"data":[)"
        R"({"lith_id":1,"name":"Siliciclastic","type":"siliciclastic","group":"","class":"sedimentary"},)"
        R"({"lith_id":43,"name":"granite","type":"plutonic","group":"felsic","class":"igneous"},)"
        R"({"name":"no id"}]}})";
    const auto v = lsdtt::parse_macrostrat_lithologies(liths);
    check("lithologies: rows without lith_id skipped", v.size() == 2);
    check("lithologies: names lower-cased", v.count(1) && v.at(1).name == "siliciclastic");
    check("lithologies: class kept", v.count(43) && v.at(43).lith_class == "igneous");

    const std::string units =
        R"({"success":{"data":[{"map_id":17,"liths":[43,1]},{"map_id":18,"liths":[]},)"
        R"({"map_id":5000000000,"liths":[1]}]}})";
    const auto u = lsdtt::parse_macrostrat_map_unit_liths(units);
    check("map units: liths in API order",
          u.count(17) && u.at(17) == std::vector<int>({43, 1}));
    check("map units: empty liths kept", u.count(18) && u.at(18).empty());
    check("map units: 64-bit map_id", u.count(5000000000LL) == 1);

    const std::string sources =
        R"({"success":{"data":[{"source_id":9,"name":"B","ref_year":"2001","scale":"1:24000"},)"
        R"({"source_id":3,"name":"A","authors":"Someone","isbn_doi":"10.1/x","ref_year":null}]}})";
    const auto s = lsdtt::parse_macrostrat_sources(sources);
    check("sources: sorted by id", s.size() == 2 && s[0].id == 3 && s[1].id == 9);
    check("sources: null year is empty", s.size() == 2 && s[0].year.empty() && s[1].year == "2001");

    check("malformed JSON refused",
          contains(thrown_message([] { (void)lsdtt::parse_macrostrat_lithologies("{not json"); }),
                   "not JSON"));
    check("error response refused",
          contains(thrown_message([] {
                       (void)lsdtt::parse_macrostrat_sources(R"({"error":{"message":"bad"}})");
                   }),
                   "error"));
    check("response without success/data refused",
          contains(thrown_message([] {
                       (void)lsdtt::parse_macrostrat_map_unit_liths(R"({"success":{"v":2}})");
                   }),
                   "success/data"));

    // The lith_primary rule.
    using lsdtt::macrostrat_primary_lithology;
    auto p = macrostrat_primary_lithology("Major:{ Alkali  Feldspar Granite ,diorite}, Minor:{x}",
                                          {1}, v);
    check("primary: first major term wins over liths",
          p.value == "alkali feldspar granite" && p.source == "major", p.value + "/" + p.source);
    p = macrostrat_primary_lithology("plutonic: granite", {43, 1}, v);
    check("primary: first known liths id", p.value == "granite" && p.source == "liths");
    p = macrostrat_primary_lithology("plutonic: granite", {999, 1}, v);
    check("primary: unknown liths ids skipped", p.value == "siliciclastic" && p.source == "liths");
    p = macrostrat_primary_lithology("Major:{}", {}, v);
    check("primary: empty major falls through to text",
          p.value == "major:{}" && p.source == "text", p.value + "/" + p.source);
    p = macrostrat_primary_lithology("  Red   SANDSTONE\t", {}, v);
    check("primary: text trimmed, lower case, squashed",
          p.value == "red sandstone" && p.source == "text");
    p = macrostrat_primary_lithology(" ", {}, v);
    check("primary: nothing -> none", p.value.empty() && p.source == "none");
}

struct read_back
{
    int epsg = 0;
    std::map<long long, std::string> primary, primary_source;
    std::map<long long, double> area;
    std::size_t n = 0;
};

read_back read_units(const std::string& path)
{
    read_back r;
    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpenEx(path.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!ds) return r;
    OGRLayer* l = ds->GetLayerByName("units");
    if (l)
    {
        if (const OGRSpatialReference* srs = l->GetSpatialRef())
        {
            const char* code = srs->GetAuthorityCode(nullptr);
            if (code) r.epsg = std::atoi(code);
        }
        OGRFeature* f = nullptr;
        while ((f = l->GetNextFeature()) != nullptr)
        {
            const long long id = f->GetFieldAsInteger64("map_id");
            r.primary[id] = f->GetFieldAsString("lith_primary");
            r.primary_source[id] = f->GetFieldAsString("lith_primary_source");
            if (auto* mp = dynamic_cast<const OGRMultiPolygon*>(f->GetGeometryRef()))
                r.area[id] += mp->get_Area();
            ++r.n;
            OGRFeature::DestroyFeature(f);
        }
    }
    GDALClose(ds);
    return r;
}

lsdtt::lsdtt_macrostrat_fetch_options fixture_options(const std::string& tile_dir,
                                                      const std::string& out)
{
    lsdtt::lsdtt_macrostrat_fetch_options o;
    o.west = 10.00; o.south = 46.00; o.east = 10.08; o.north = 46.05;
    o.zoom = 10;
    o.tile_dir = tile_dir;
    o.lookup_map_units = false;
    o.request_delay_milliseconds = 0;
    o.out_path = out;
    return o;
}

void test_fetch(const std::string& tile_dir, const std::string& work)
{
    const std::string out = work + "/fixture_units.fgb";
    lsdtt::lsdtt_macrostrat_fetch_result r;
    const std::string err = thrown_message([&] { r = lsdtt::fetch_macrostrat_units(fixture_options(tile_dir, out)); });
    check("fixture fetch runs", err.empty(), err);
    if (!err.empty()) return;

    check("one tile read", r.n_tiles == 1 && r.n_tiles_missing == 0);
    check("three units inside the box (limestone outside dropped)",
          r.n_features == 3 && r.n_map_units == 3,
          std::to_string(r.n_features) + " features");
    check("output EPSG is UTM 32N", r.epsg == 32632, std::to_string(r.epsg));
    check("no warnings offline", r.warnings.empty(),
          r.warnings.empty() ? "" : r.warnings.front());

    const read_back b = read_units(out);
    check("written layer: 3 features in EPSG:32632", b.n == 3 && b.epsg == 32632,
          std::to_string(b.n) + " / " + std::to_string(b.epsg));
    check("map 1: major term", b.primary.count(1) && b.primary.at(1) == "granite" &&
                                   b.primary_source.at(1) == "major");
    check("map 2: free text squashed",
          b.primary.count(2) && b.primary.at(2) == "sandstone and shale" &&
              b.primary_source.at(2) == "text");
    check("map 3 (outside the box) absent", b.primary.count(3) == 0);
    check("map 4: none", b.primary.count(4) && b.primary.at(4).empty() &&
                             b.primary_source.at(4) == "none");

    // Areas against the source rectangles projected independently. The tile
    // grid quantises vertices to 39135.76 m / 4096 = 9.6 m, so a 3 x 5.6 km
    // rectangle may differ by about 2 x 9.6 m x (3 + 5.6) km = 0.2 km2 (1 %).
    const double a1 = projected_rectangle_area(10.00, 46.00, 10.04, 46.05, 32632);
    const double a2 = projected_rectangle_area(10.04, 46.00, 10.08, 46.05, 32632);
    check("map 1 area within 1.5 % of the source polygon",
          b.area.count(1) && std::abs(b.area.at(1) / a1 - 1.0) < 0.015,
          b.area.count(1) ? std::to_string(b.area.at(1)) + " vs " + std::to_string(a1) : "");
    check("map 2 area within 1.5 % of the source polygon",
          b.area.count(2) && std::abs(b.area.at(2) / a2 - 1.0) < 0.015);

    // The vocabulary: one row per lith_primary, largest area first.
    std::map<std::string, double> vocab;
    for (const auto& v : r.vocabulary) vocab[v.lith_primary] = v.area_km2;
    check("vocabulary: granite, sandstone and shale, empty",
          r.vocabulary.size() == 3 && vocab.count("granite") && vocab.count("sandstone and shale") &&
              vocab.count(""));
    check("vocabulary sorted by area",
          r.vocabulary.size() == 3 && r.vocabulary[0].area_km2 >= r.vocabulary[1].area_km2 &&
              r.vocabulary[1].area_km2 >= r.vocabulary[2].area_km2);
    check("vocabulary area = feature areas",
          vocab.count("granite") && b.area.count(1) &&
              std::abs(vocab["granite"] - b.area.at(1) / 1e6) < 1e-6);

    // Without the API the attribution still names Macrostrat and each source id.
    const std::string attr = lsdtt::macrostrat_attribution_text(r);
    check("attribution cites Macrostrat and source 7",
          contains(attr, "doi:10.1029/2018GC007467") && contains(attr, "source_id 7") &&
              contains(attr, "CC BY 4.0"));
    check("attribution omits source 8 (only outside the box)", !contains(attr, "source_id 8"));

    // An explicit output CRS is honoured.
    auto o = fixture_options(tile_dir, work + "/fixture_units_32633.gpkg");
    o.target_epsg = 32633;
    o.out_driver = "GPKG";
    lsdtt::lsdtt_macrostrat_fetch_result r33;
    check("explicit EPSG:32633 to GeoPackage",
          thrown_message([&] { r33 = lsdtt::fetch_macrostrat_units(o); }).empty() &&
              r33.epsg == 32633 && read_units(o.out_path).epsg == 32633);

    // Failures.
    o = fixture_options(tile_dir, work + "/unused.fgb");
    o.max_tiles = 0;
    check("too many tiles refused",
          contains(thrown_message([&] { (void)lsdtt::fetch_macrostrat_units(o); }), "max_tiles"));
    std::filesystem::create_directories(work + "/empty_tiles");
    o = fixture_options(work + "/empty_tiles", work + "/unused.fgb");
    check("no readable tile refused",
          contains(thrown_message([&] { (void)lsdtt::fetch_macrostrat_units(o); }),
                   "no Macrostrat tile could be read"));
    o = fixture_options(tile_dir, work + "/unused.fgb");
    o.west = 9.90; o.east = 9.95;     // inside the tile, west of every unit
    check("box with no units refused",
          contains(thrown_message([&] { (void)lsdtt::fetch_macrostrat_units(o); }),
                   "no map units inside"));

    // An unreachable API degrades to a warning and the tiles' text.
    if (lsdtt::http_fetch_available())
    {
        o = fixture_options(tile_dir, work + "/fixture_units_noapi.fgb");
        o.lookup_map_units = true;
        o.api_base = "http://127.0.0.1:9/api/v2";
        o.http.max_retry = 0;
        o.http.timeout_seconds = 5;
        o.http.connect_timeout_seconds = 5;
        lsdtt::lsdtt_macrostrat_fetch_result rn;
        const std::string e = thrown_message([&] { rn = lsdtt::fetch_macrostrat_units(o); });
        bool warned = false;
        for (const auto& w : rn.warnings) warned = warned || contains(w, "API lookup failed");
        const read_back bn = read_units(o.out_path);
        check("unreachable API: warning, output still written", e.empty() && warned && bn.n == 3, e);
        check("unreachable API: lith_primary from the tiles' text",
              bn.primary.count(1) && bn.primary.at(1) == "granite");
    }
    else
    {
        std::cout << "SKIP  unreachable API (this GDAL has no HTTP support)\n";
    }
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "usage: lsdtt-test-macrostrat-fetch <fixture_tile_dir> <work_dir>\n";
        return 2;
    }
    const std::string tile_dir = argv[1];
    const std::string work = argv[2];
    std::filesystem::create_directories(work);
    GDALAllRegister();

    test_tiles();
    test_parsers();
    if (!GetGDALDriverManager()->GetDriverByName("MVT"))
    {
        std::cout << "SKIP  fixture fetch: this GDAL has no MVT driver\n";
        std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks << " passed\n";
        return g_n_failed == 0 ? 77 : 1;
    }
    test_fetch(tile_dir, work);

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks << " passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
