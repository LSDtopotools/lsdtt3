// test_octopus_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Tests for the OCTOPUS erosion-rate fetch (lsdtt-fetch-erosion-rates' library),
// offline; no HTTP request is made.
//
//   * the request recipe: the endpoint and layers, the latitude-first bbox,
//     the hits / page / by-id URLs as exact strings, never maxFeatures;
//   * numberMatched and ExceptionText parsing; the paging plan;
//   * the no-data sentinel predicate and the rate-field map of every
//     scheme / nuclide;
//   * offline reads of the committed WV fixtures
//     (tests/data/octopus_erosion_rates): basin and outlet counts, masked
//     values per field, the join, the CSV rows (no sentinel reaches the CSV),
//     the internal uncertainty, the polygons, the DEM box and the attribution;
//   * failures: too many basins (the canned hits response), an
//     ExceptionReport body.
//
// Usage: lsdtt-test-octopus-fetch <fixture_dir> <work_dir>
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_octopus_fetch.hpp"

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

bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

std::string lower(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string read_text(const std::string& path)
{
    std::ifstream in(path);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

void test_recipe(const std::string& dir)
{
    check("endpoint is the HTTP geoserver host",
          lsdtt::octopus_wfs_endpoint() == "http://geoserver.octopusdata.org/geoserver/wfs");
    check("layer names", lsdtt::octopus_basins_layer() == "be10-denude:crn_int_basins" &&
                             lsdtt::octopus_outlets_layer() == "be10-denude:crn_int_outlets");

    const std::string big = lsdtt::octopus_wfs_bbox_parameter(-79.80, 38.85, -79.05, 39.15);
    const std::string large = lsdtt::octopus_wfs_bbox_parameter(-80.5, 38.3, -78.5, 39.7);
    check("bbox is latitude first (S,W,N,E) with the URN", big == "38.85,-79.8,39.15,-79.05,urn:ogc:def:crs:EPSG::4326", big);
    check("bbox of the large box", large == "38.3,-80.5,39.7,-78.5,urn:ogc:def:crs:EPSG::4326", large);

    const std::string hits = lsdtt::octopus_wfs_hits_url(lsdtt::octopus_basins_layer(), large);
    check("hits URL exact", hits ==
          "http://geoserver.octopusdata.org/geoserver/wfs?service=WFS&version=2.0.0&request=GetFeature"
          "&typeNames=be10-denude:crn_int_basins&resultType=hits"
          "&bbox=38.3,-80.5,39.7,-78.5,urn:ogc:def:crs:EPSG::4326&srsName=EPSG:4326", hits);
    const std::string page = lsdtt::octopus_wfs_page_url(lsdtt::octopus_basins_layer(), large, 10, 20);
    check("page URL exact", page ==
          "http://geoserver.octopusdata.org/geoserver/wfs?service=WFS&version=2.0.0&request=GetFeature"
          "&typeNames=be10-denude:crn_int_basins&outputFormat=application/json&srsName=EPSG:4326"
          "&bbox=38.3,-80.5,39.7,-78.5,urn:ogc:def:crs:EPSG::4326&sortBy=OBSID1+A&count=10&startIndex=20",
          page);
    const std::string by_id = lsdtt::octopus_wfs_outlets_by_id_url({"S186WTS001", "S186WTS007"});
    check("outlets-by-id URL exact (CQL_FILTER, no bbox)", by_id ==
          "http://geoserver.octopusdata.org/geoserver/wfs?service=WFS&version=2.0.0&request=GetFeature"
          "&typeNames=be10-denude:crn_int_outlets&outputFormat=application/json&srsName=EPSG:4326"
          "&sortBy=OBSID1+A&count=2"
          "&CQL_FILTER=OBSID1%20IN%20%28%27S186WTS001%27%2C%27S186WTS007%27%29", by_id);
    bool clean = true;
    for (const std::string& u : {hits, page, by_id})
        clean = clean && !contains(lower(u), "maxfeatures") && u.rfind("http://geoserver.octopusdata.org", 0) == 0;
    check("no URL carries maxFeatures, and all start with http://geoserver.octopusdata.org", clean);

    const auto starts = lsdtt::octopus_wfs_page_starts(28, 10);
    check("paging plan for 28 features, page 10: startIndex 0, 10, 20",
          starts == std::vector<std::int64_t>{0, 10, 20});
    check("paging plan for 30 / 10 and 0 / 10",
          lsdtt::octopus_wfs_page_starts(30, 10) == std::vector<std::int64_t>{0, 10, 20} &&
          lsdtt::octopus_wfs_page_starts(0, 10).empty());
    check("page size 0 throws", !thrown_message([] { (void)lsdtt::octopus_wfs_page_starts(5, 0); }).empty());

    const std::string hits_xml = read_text(dir + "/wv_large_box/hits_basins.xml");
    check("numberMatched of the canned hits response is 28",
          lsdtt::parse_octopus_wfs_number_matched(hits_xml) == 28);
    check("numberMatched of a GeoJSON body",
          lsdtt::parse_octopus_wfs_number_matched(read_text(dir + "/wv_small_box/crn_int_outlets.json")) == 1);
    check("numberMatched absent -> -1", lsdtt::parse_octopus_wfs_number_matched("{}") == -1);

    const std::string ex = lsdtt::octopus_wfs_exception_text(read_text(dir + "/exception_report/crn_int_basins.json"));
    check("ExceptionText extracted", ex == "Feature type be10-denude:crn_int_basinz unknown", ex);
    check("no ExceptionText for a hits response or GeoJSON",
          lsdtt::octopus_wfs_exception_text(hits_xml).empty() &&
          lsdtt::octopus_wfs_exception_text(read_text(dir + "/wv_small_box/crn_int_outlets.json")).empty());

    bool sentinels = true;
    for (double x : {-9999.0, -9999.99, -99.99, -99.999, -9.999, -9.99999, -999.0})
        sentinels = sentinels && lsdtt::is_octopus_no_data_sentinel(x);
    bool values = true;
    for (double x : {-9998.9, 0.0, -9.9, -99.0, -9.0, 12.88, 9999.0, -9999.5, std::nan("")})
        values = values && !lsdtt::is_octopus_no_data_sentinel(x);
    check("sentinel predicate: true for the all-nines no-data values", sentinels);
    check("sentinel predicate: false for -9998.9, 0, -9.9, -99, data and NaN", values);

    struct expect { const char* scheme; const char* nuc; const char* rate; const char* unc; const char* inner; const char* mass; const char* kind; };
    const expect table[] = {
        {"st", "be10", "EBE_MMKYR", "EBE_ERR", "ERRBE_AMS", "EBE_GCMYR", "total_1sigma"},
        {"lsdn", "be10", "EBEL", "EBEL_EXT", "EBEL_INT", "", "external_1sigma"},
        {"lsdn_no_topo_shielding", "be10", "EBELNS", "EBELNS_EXT", "EBELNS_INT", "", "external_1sigma"},
        {"published", "be10", "BE10EP", "BE10EP_ERR", "", "", "published_1sigma"},
        {"st", "al26", "EAL_MMKYR", "EAL_ERR", "ERRAL_AMS", "EAL_GCMYR", "total_1sigma"},
        {"lsdn", "al26", "EALL", "EALL_EXT", "EALL_INT", "", "external_1sigma"},
        {"lsdn_no_topo_shielding", "al26", "EALLNS", "EALLNS_EXT", "EALLNS_INT", "", "external_1sigma"},
        {"published", "al26", "AL26EP", "AL26EP_ERR", "", "", "published_1sigma"}};
    bool map_ok = true;
    for (const auto& e : table)
    {
        const auto f = lsdtt::octopus_rate_fields(e.scheme, e.nuc);
        map_ok = map_ok && f.rate == e.rate && f.uncertainty == e.unc && f.internal_uncertainty == e.inner &&
                 f.internal_mass_flux_rate == e.mass && f.kind == e.kind;
    }
    check("rate-field map exact for every scheme and nuclide", map_ok);
    const std::string bad = thrown_message([] { (void)lsdtt::octopus_rate_fields("cronus", "be10"); });
    check("unknown scheme throws listing the schemes",
          contains(bad, "st, lsdn, lsdn_no_topo_shielding, published"), bad);
    check("unknown nuclide throws", contains(thrown_message([] { (void)lsdtt::octopus_rate_fields("st", "c14"); }), "be10 or al26"));
}

void test_offline(const std::string& dir, const std::string& work)
{
    lsdtt::lsdtt_octopus_fetch_options o;
    o.west = -80.5; o.east = -78.5; o.south = 38.3; o.north = 39.7;
    o.response_dir = dir + "/wv_large_box";
    const auto r = lsdtt::fetch_octopus_erosion_rates(o);

    check("large box: 28 basins, numberMatched 28", r.basins.size() == 28 && r.n_matched == 28);
    check("large box: 23 outlets, 23 basins with the outlet in the box, 23 rows",
          r.n_outlets == 23 && r.n_outlet_in_bbox == 23 && r.rows.size() == 23,
          std::to_string(r.n_outlets) + " " + std::to_string(r.n_outlet_in_bbox) + " " + std::to_string(r.rows.size()));
    check("large box: 94 fields (93 attributes + the WFS feature id), no request made",
          r.fields.size() == 94 && r.fields.front().name == "id" && r.n_requests == 0,
          std::to_string(r.fields.size()));
    std::size_t total = 0;
    for (const auto& [f, n] : r.masked_by_field) total += n;
    auto masked = [&](const char* f) { const auto it = r.masked_by_field.find(f); return it == r.masked_by_field.end() ? 0u : it->second; };
    check("masked values: 896 in 32 fields", total == 896 && r.masked_by_field.size() == 32,
          std::to_string(total) + " in " + std::to_string(r.masked_by_field.size()));
    check("masked: every Al-26 rate, the glacier and sampling-date fields (28 each)",
          masked("EAL_MMKYR") == 28 && masked("EAL_ERR") == 28 && masked("EAL_GCMYR") == 28 &&
          masked("ALTOPO") == 28 && masked("GLA_PCNT") == 28 && masked("AL26NC") == 28 &&
          masked("SMP_YEAR") == 28 && masked("EBE_MMKYR") == 0);

    const lsdtt::lsdtt_octopus_csv_row* s046 = nullptr;
    for (const auto& w : r.rows) if (w.sample_id == "S186WTS046") s046 = &w;
    check("S186WTS046 row: outlet, rate, uncertainty, area",
          s046 && close(s046->latitude, 38.88961897, 1e-9) && close(s046->longitude, -79.40200545, 1e-9) &&
          close(s046->rate, 12.88, 1e-12) && close(s046->uncertainty, 2.62, 1e-12) &&
          close(s046->area_km2, 20.14, 1e-12) && s046->outlet_in_bbox && s046->study_id == "S186" &&
          s046->ref_doi == "10.1130/b31840.1");
    double expected_internal = std::nan("");
    for (const auto& b : r.basins)
        if (lsdtt::octopus_basin_text(r, b, "OBSID1") == "S186WTS046")
            expected_internal = lsdtt::octopus_basin_number(r, b, "ERRBE_AMS") *
                                lsdtt::octopus_basin_number(r, b, "EBE_MMKYR") /
                                lsdtt::octopus_basin_number(r, b, "EBE_GCMYR");
    check("S186WTS046 internal uncertainty = ERRBE_AMS * EBE_MMKYR / EBE_GCMYR",
          s046 && std::isfinite(expected_internal) && close(s046->internal_uncertainty, expected_internal, 1e-12));
    bool s001_out = std::find(r.no_outlet_in_bbox.begin(), r.no_outlet_in_bbox.end(), "S186WTS001") != r.no_outlet_in_bbox.end();
    for (const auto& w : r.rows) s001_out = s001_out && w.sample_id != "S186WTS001";
    check("S186WTS001 (outlet outside the box) is not a row", s001_out && r.no_outlet_in_bbox.size() == 5);
    std::size_t s046_dupes = 0;
    for (const auto& w : r.rows)
        if (w.sample_id >= "S046WTS009" && w.sample_id <= "S046WTS012" && close(w.area_km2, 9.45, 1e-12)) ++s046_dupes;
    check("S046WTS009-012 are four rows, each AREA 9.45 km2", s046_dupes == 4);

    std::filesystem::create_directories(work);
    const std::string csv = work + "/off_erosion_rates.csv";
    lsdtt::write_octopus_erosion_rate_csv(csv, r);
    const std::string text = read_text(csv);
    check("CSV header exact", text.rfind(lsdtt::octopus_erosion_rate_csv_header() + "\n", 0) == 0 &&
          lsdtt::octopus_erosion_rate_csv_header() ==
          "sample_id,latitude,longitude,erosion_rate,erosion_rate_uncertainty,erosion_rate_units,"
          "erosion_rate_uncertainty_kind,erosion_rate_internal_uncertainty,reported_area_km2,"
          "rate_scheme,rate_field,nuclide,obsid2,study_id,author,pub_year,ref_doi,mean_elevation_m,outlet_in_bbox");
    std::size_t lines = 0;
    for (char c : text) if (c == '\n') ++lines;
    check("CSV: 23 data rows, mm_per_kyr, no sentinel", lines == 24 && contains(text, ",mm_per_kyr,total_1sigma,") &&
          !contains(text, "-9999") && !contains(text, "-99.99") && !contains(text, "-9.999"));

    const std::string fgb = work + "/off_octopus_basins.fgb";
    lsdtt::write_octopus_basin_polygons(fgb, "FlatGeobuf", 32617, r);
    {
        GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(fgb.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
        OGRLayer* l = ds ? ds->GetLayer(0) : nullptr;
        const OGRSpatialReference* srs = l ? l->GetSpatialRef() : nullptr;
        const char* code = srs ? srs->GetAuthorityCode(nullptr) : nullptr;
        const int n_fields = l ? l->GetLayerDefn()->GetFieldCount() : 0;
        std::size_t n_null_al = 0;
        if (l)
        {
            const int i = l->GetLayerDefn()->GetFieldIndex("EAL_MMKYR");
            OGRFeature* f = nullptr;
            l->ResetReading();
            while ((f = l->GetNextFeature()) != nullptr)
            {
                if (i >= 0 && !f->IsFieldSetAndNotNull(i)) ++n_null_al;
                OGRFeature::DestroyFeature(f);
            }
        }
        check("polygons: 28 features in EPSG:32617, 94 + 3 fields, sentinels null",
              l && l->GetFeatureCount() == 28 && code && std::string(code) == "32617" && n_fields == 97 &&
              n_null_al == 28, std::to_string(n_fields));
        if (ds) GDALClose(ds);
    }

    const auto db = lsdtt::octopus_dem_bounds(r, 0.02, 1000.0);
    std::cout << std::setprecision(12) << "DEM box: " << db.n_basins << " basins, " << db.west << " "
              << db.south << " " << db.east << " " << db.north << "\n";
    check("DEM box: 21 basins (outlet in the box, AREA <= 1000 km2), padded envelope",
          db.ok && db.n_basins == 21 && close(db.west, -79.719166666, 1e-6) && close(db.south, 38.270833333, 1e-6) &&
          close(db.east, -78.51, 1e-6) && close(db.north, 39.729166666, 1e-6));

    const std::string attr = lsdtt::octopus_attribution_text(r);
    bool dois = true;
    for (const char* d : {"10.5194/essd-10-2123-2018", "10.5194/essd-14-3695-2022", "10.5194/gchron-7-113-2025",
                          "10.71747/uow-r3gk326m.28216865.v1", "10.5194/esurf-4-655-2016", "CC BY 4.0"})
        dois = dois && contains(attr, d);
    check("attribution: every DOI and the licence", dois);
    check("attribution: one line per study",
          contains(attr, "S186, Portenga, 2019, 10.1130/b31840.1, 17") &&
          contains(attr, "S046, Duxburry, 2009, (no DOI; REFID phdthesis), 11") &&
          contains(attr, "DBVERNO): v2.3(4)"), attr);

    // lsdn: the EBEL fields.
    lsdtt::lsdtt_octopus_fetch_options ol = o;
    ol.rate_scheme = "lsdn";
    const auto rl = lsdtt::fetch_octopus_erosion_rates(ol);
    bool lsdn_ok = rl.rows.size() == 23 && rl.rate_fields.rate == "EBEL";
    for (const auto& w : rl.rows)
        if (w.sample_id == "S186WTS046") lsdn_ok = lsdn_ok && close(w.rate, 12.8, 1e-12) && close(w.uncertainty, 0.83, 1e-12);
    check("lsdn scheme: EBEL / EBEL_EXT rows (S186WTS046 12.8 +- 0.83)", lsdn_ok);

    // al26: no basin here has an Al-26 rate.
    lsdtt::lsdtt_octopus_fetch_options oa = o;
    oa.nuclide = "al26";
    const auto ra = lsdtt::fetch_octopus_erosion_rates(oa);
    check("al26: every in-box basin dropped for no usable EAL_MMKYR",
          ra.rows.empty() && ra.dropped_no_rate.size() == 23);

    // The small (wvbig) box: 5 basins, 1 outlet.
    lsdtt::lsdtt_octopus_fetch_options os;
    os.west = -79.80; os.east = -79.05; os.south = 38.85; os.north = 39.15;
    os.response_dir = dir + "/wv_small_box";
    const auto rs = lsdtt::fetch_octopus_erosion_rates(os);
    check("wvbig box: 5 basins, 1 with its outlet in the box (S186WTS046)",
          rs.basins.size() == 5 && rs.n_outlet_in_bbox == 1 && rs.rows.size() == 1 &&
          rs.rows[0].sample_id == "S186WTS046");

    // Failures.
    lsdtt::lsdtt_octopus_fetch_options om = o;
    om.max_features = 5;
    const std::string too_many = thrown_message([&] { (void)lsdtt::fetch_octopus_erosion_rates(om); });
    check("max_features 5 against the hits fixture fails",
          contains(too_many, "OCTOPUS returns 28 basins for this box, more than max_features (5)"), too_many);
    lsdtt::lsdtt_octopus_fetch_options oe = o;
    oe.response_dir = dir + "/exception_report";
    const std::string exc = thrown_message([&] { (void)lsdtt::fetch_octopus_erosion_rates(oe); });
    check("an ExceptionReport body fails with its text",
          contains(exc, "OCTOPUS WFS error: Feature type be10-denude:crn_int_basinz unknown"), exc);
    lsdtt::lsdtt_octopus_fetch_options ob = o;
    ob.west = -78.0;
    check("a box with west >= east throws", !thrown_message([&] { (void)lsdtt::fetch_octopus_erosion_rates(ob); }).empty());
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "usage: lsdtt-test-octopus-fetch <fixture_dir> <work_dir>\n";
        return 2;
    }
    GDALAllRegister();
    const std::string dir = argv[1], work = argv[2];
    std::cout << "lsdtt3 OCTOPUS fetch tests (offline)\n";
    test_recipe(dir);
    try { test_offline(dir, work); }
    catch (const std::exception& e) { check("offline fetch ran without an exception", false, e.what()); }
    std::cout << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
