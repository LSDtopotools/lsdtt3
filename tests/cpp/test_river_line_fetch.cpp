// test_river_line_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Tests for lsdtt_river_line_fetch (lsdtt-fetch-river-lines), offline.
//
//   * md5 against the RFC 1321 test suite; the direction mapping parser; the
//     OS Downloads API parser on a frozen response; the feature class map;
//   * the gb_os_open_rivers adapter on the committed Open Rivers extract
//     (tests/data/river_lines/os_open_rivers_extract), zipped into a cache
//     directory here: the network and every node type must equal a reference
//     built in this test directly from the start_node / end_node ids;
//   * local_file with the Open Rivers field names gives the same lines;
//   * a copy with every line reversed and flow_direction flipped gives the
//     same lines (orientation comes from the attribute);
//   * coordinate-only noding (no ids) gives the same node count;
//   * the class filter, the outputs read back, and a missing field refused.
//
// Usage: lsdtt-test-river-line-fetch <fixture_dir> <work_dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <cpl_vsi.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_river_line_fetch.hpp"
#include "lsdtt_river_line_prepare_cpu.hpp"

namespace fs = std::filesystem;

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

void check_int(const std::string& name, long got, long want)
{
    check(name, got == want, "got " + std::to_string(got) + ", want " + std::to_string(want));
}

template <typename F>
std::string thrown_message(F f)
{
    try { f(); }
    catch (const std::exception& e) { return e.what(); }
    return "";
}

std::string read_file(const fs::path& p)
{
    std::ifstream f(p, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

void vsi_write(const std::string& path, const std::string& bytes)
{
    VSILFILE* f = VSIFOpenL(path.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot open " + path);
    VSIFWriteL(bytes.data(), 1, bytes.size(), f);
    VSIFCloseL(f);
}

// Reference network from the fixture's node ids alone.
struct reference
{
    int n_links = 0;
    int n_nodes = 0;
    int n_loops = 0;
    std::map<std::string, int> types;
    double x0 = 1e300, y0 = 1e300, x1 = -1e300, y1 = -1e300;
};

reference build_reference(const std::string& gpkg)
{
    reference ref;
    std::unique_ptr<GDALDataset> ds(static_cast<GDALDataset*>(
        GDALOpenEx(gpkg.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr)));
    OGRLayer* l = ds->GetLayerByName("watercourse_link");
    std::map<std::string, int> din, dout;
    std::map<std::string, std::string> parent;
    std::function<std::string(const std::string&)> find = [&](const std::string& a)
    {
        if (parent[a] == a) return a;
        return parent[a] = find(parent[a]);
    };
    for (auto& f : *l)
    {
        const std::string s = f->GetFieldAsString("start_node");
        const std::string e = f->GetFieldAsString("end_node");
        const bool flip = std::string(f->GetFieldAsString("flow_direction")) == "in opposite direction";
        ++dout[flip ? e : s];
        ++din[flip ? s : e];
        if (!parent.count(s)) parent[s] = s;
        if (!parent.count(e)) parent[e] = e;
        const std::string a = find(s), b = find(e);
        if (a == b) ++ref.n_loops; else parent[a] = b;
        ++ref.n_links;
        OGREnvelope env;
        f->GetGeometryRef()->getEnvelope(&env);
        ref.x0 = std::min(ref.x0, env.MinX); ref.x1 = std::max(ref.x1, env.MaxX);
        ref.y0 = std::min(ref.y0, env.MinY); ref.y1 = std::max(ref.y1, env.MaxY);
    }
    std::set<std::string> nodes;
    for (const auto& kv : parent) nodes.insert(kv.first);
    ref.n_nodes = static_cast<int>(nodes.size());
    // The extract is one network. If its data have exactly one sink (inflow,
    // no outflow), the network drains there: that sink is the outlet even with
    // several inflows (the extract cut off its outflow).
    int n_sinks = 0;
    for (const std::string& n : nodes) if (dout[n] == 0 && din[n] >= 1) ++n_sinks;
    for (const std::string& n : nodes)
    {
        const int i = din[n], o = dout[n];
        std::string t;
        if (n_sinks == 1 && o == 0 && i >= 1) t = "outlet";
        else if (i == 0 && o == 1) t = "tip";
        else if (i == 1 && o == 0) t = "outlet";
        else if (i >= 1 && o >= 2) t = "bifurcation";
        else if (i >= 2 && o == 1) t = "confluence";
        else if (i == 1 && o == 1) t = "link";
        else t = "conflict";
        ++ref.types[t];
    }
    return ref;
}

std::map<std::string, int> type_counts(const lsdtt::river_line_network& net)
{
    std::map<std::string, int> m;
    for (const auto& n : net.nodes) ++m[n.node_type];
    return m;
}

std::string describe(const std::map<std::string, int>& m)
{
    std::string s;
    for (const auto& kv : m) s += kv.first + "=" + std::to_string(kv.second) + " ";
    return s;
}

// Lines keyed by source id, as vertex lists.
std::map<std::string, std::vector<lsdtt::lsdtt_xy>> lines_by_id(const lsdtt::river_line_network& net)
{
    std::map<std::string, std::vector<lsdtt::lsdtt_xy>> m;
    for (const auto& l : net.lines) m[l.source_feature_id] = l.vertices;
    return m;
}

bool same_lines(const lsdtt::river_line_network& a, const lsdtt::river_line_network& b, double eps,
                std::string& why)
{
    const auto la = lines_by_id(a), lb = lines_by_id(b);
    if (la.size() != lb.size()) { why = "line counts differ"; return false; }
    for (const auto& kv : la)
    {
        const auto it = lb.find(kv.first);
        if (it == lb.end()) { why = "missing id " + kv.first; return false; }
        if (it->second.size() != kv.second.size()) { why = "vertex count differs for " + kv.first; return false; }
        for (size_t i = 0; i < kv.second.size(); ++i)
        {
            if (std::hypot(kv.second[i].x - it->second[i].x, kv.second[i].y - it->second[i].y) > eps)
            {
                why = "vertex " + std::to_string(i) + " of " + kv.first + " differs";
                return false;
            }
        }
    }
    return true;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "usage: lsdtt-test-river-line-fetch <fixture_dir> <work_dir>\n";
        return 2;
    }
    GDALAllRegister();
    const fs::path fixture = fs::path(argv[1]) / "os_open_rivers_extract";
    const fs::path work = argv[2];
    fs::remove_all(work);
    fs::create_directories(work);

    // ---- pure helpers -------------------------------------------------------
    check("md5 ''", lsdtt::md5_hex("") == "d41d8cd98f00b204e9800998ecf8427e");
    check("md5 'abc'", lsdtt::md5_hex("abc") == "900150983cd24fb0d6963f7d28e17f72");
    check("md5 'message digest'", lsdtt::md5_hex("message digest") == "f96b697d7cb7938d525a2f31aaf161d0");
    check("md5 80 digits (two blocks)",
          lsdtt::md5_hex("1234567890123456789012345678901234567890"
                         "1234567890123456789012345678901234567890") == "57edf4a22be3c955ac49da2e2107b67a");

    const auto dm = lsdtt::parse_direction_field_mapping("in_direction=1, in_opposite_direction=-1");
    check("mapping: two values", dm.size() == 2 && dm.at("in_direction") == 1 &&
                                 dm.at("in_opposite_direction") == -1);
    check("mapping: a space in a key is stored as _",
          lsdtt::parse_direction_field_mapping("in direction=1").count("in_direction") == 1);
    check("mapping: empty text gives empty map", lsdtt::parse_direction_field_mapping("").empty());
    check("mapping: missing = refused",
          thrown_message([] { lsdtt::parse_direction_field_mapping("downstream"); }).find("value=1") !=
              std::string::npos);
    check("mapping: bad sign refused",
          thrown_message([] { lsdtt::parse_direction_field_mapping("a=2"); }).find("must be 1 or -1") !=
              std::string::npos);

    const std::string api_json =
        "[{\"md5\":\"12b3ee038ea7ef29b1789e040acb0f39\",\"size\":37876554,\"url\":\"https://x/shp\","
        "\"format\":\"ESRI\\u00ae Shapefile\",\"area\":\"GB\",\"fileName\":\"oprvrs_essh_gb.zip\"},"
        "{\"md5\":\"DD777B3F966FB814AE2352D9C2076CD4\",\"size\":52461242,"
        "\"url\":\"https://api.os.uk/downloads/v1/products/OpenRivers/downloads?area=GB&format=GeoPackage&redirect\","
        "\"format\":\"GeoPackage\",\"area\":\"GB\",\"fileName\":\"oprvrs_gpkg_gb.zip\"}]";
    const auto e = lsdtt::parse_os_downloads_geopackage(api_json);
    check("OS API: GeoPackage entry chosen", e.file_name == "oprvrs_gpkg_gb.zip" && e.size == 52461242);
    check("OS API: md5 lower-cased", e.md5 == "dd777b3f966fb814ae2352d9c2076cd4");
    check("OS API: no GeoPackage refused",
          thrown_message([] { lsdtt::parse_os_downloads_geopackage("[]"); }).find("no GeoPackage") !=
              std::string::npos);

    check("class: inlandRiver -> river", lsdtt::river_feature_class_from_value("inlandRiver") == "river");
    check("class: lake -> lake_connector", lsdtt::river_feature_class_from_value("lake") == "lake_connector");
    check("class: tidalRiver -> tidal", lsdtt::river_feature_class_from_value("tidalRiver") == "tidal");
    check("class: ditch -> drain", lsdtt::river_feature_class_from_value(" Ditch ") == "drain");
    check("class: unknown -> other", lsdtt::river_feature_class_from_value("waterfall") == "other");

    // ---- cache built from the fixture ---------------------------------------
    const fs::path cache = work / "cache";
    fs::create_directories(cache / "gb_os_open_rivers");
    const std::string zip = (cache / "gb_os_open_rivers" / "oprvrs_gpkg_gb.zip").string();
    vsi_write("/vsizip/" + zip + "/Data/oprvrs_gb.gpkg", read_file(fixture / "Data" / "oprvrs_gb.gpkg"));
    vsi_write("/vsizip/" + zip + "/Doc/licence.txt", read_file(fixture / "Doc" / "licence.txt"));

    const std::string gpkg = (fixture / "Data" / "oprvrs_gb.gpkg").string();
    const reference ref = build_reference(gpkg);

    // A WGS84 box around the whole extract, so nothing is clipped.
    double west, south, east, north;
    {
        OGRSpatialReference bng, ll;
        bng.importFromEPSG(27700);
        ll.importFromEPSG(4326);
        bng.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        ll.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        std::unique_ptr<OGRCoordinateTransformation> ct(OGRCreateCoordinateTransformation(&bng, &ll));
        double xs[4] = {ref.x0, ref.x1, ref.x0, ref.x1};
        double ys[4] = {ref.y0, ref.y0, ref.y1, ref.y1};
        ct->Transform(4, xs, ys);
        west = std::min(xs[0], xs[2]) - 0.02; east = std::max(xs[1], xs[3]) + 0.02;
        south = std::min(ys[0], ys[1]) - 0.02; north = std::max(ys[2], ys[3]) + 0.02;
    }

    lsdtt::lsdtt_river_line_fetch_options base;
    base.west = west; base.south = south; base.east = east; base.north = north;
    base.repair_direction = false;   // the reference is the data's own directions
    base.output_epsg = 27700;   // the source CRS: coordinates pass through exactly

    lsdtt::lsdtt_river_line_fetch_options os = base;
    os.source = "gb_os_open_rivers";
    os.cache_dir = cache.string();
    os.out_lines_path = (work / "os_river_lines.fgb").string();
    os.out_nodes_path = (work / "os_river_nodes.fgb").string();
    const auto r_os = lsdtt::fetch_river_lines(os);

    check_int("OS: features read", r_os.n_features_read, ref.n_links);
    check_int("OS: lines (already noded: no splits)", static_cast<long>(r_os.network.lines.size()), ref.n_links);
    check_int("OS: T-junction splits", r_os.network.n_t_junction_splits, 0);
    check_int("OS: nodes = distinct node ids", static_cast<long>(r_os.network.nodes.size()), ref.n_nodes);
    check("OS: node types equal the id reference", type_counts(r_os.network) == ref.types,
          describe(type_counts(r_os.network)) + "vs " + describe(ref.types));
    int loops = 0;
    for (const auto& c : r_os.network.components) loops += c.n_loops;
    check_int("OS: loops equal the id reference", loops, ref.n_loops);
    check_int("OS: EPSG", r_os.epsg, 27700);
    check("OS: attribution year from the licence file",
          r_os.attribution.find("database rights 2026") != std::string::npos, r_os.attribution);
    check("OS: every direction from the attribute",
          r_os.direction_sources.size() == 1 && r_os.direction_sources.count("attribute"));
    check("OS: lake links kept by default", r_os.classes.count("lake_connector") &&
                                            r_os.classes.at("lake_connector").n_kept > 0);
    {
        std::unique_ptr<GDALDataset> ds(static_cast<GDALDataset*>(
            GDALOpenEx(os.out_lines_path.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr)));
        OGRLayer* l = ds ? ds->GetLayer(0) : nullptr;
        check("OS: lines file read back",
              l && l->GetFeatureCount() == static_cast<GIntBig>(r_os.network.lines.size()) &&
                  l->GetLayerDefn()->GetFieldIndex("on_loop") >= 0 &&
                  l->GetLayerDefn()->GetFieldIndex("direction_source") >= 0);
        std::unique_ptr<GDALDataset> dn(static_cast<GDALDataset*>(
            GDALOpenEx(os.out_nodes_path.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr)));
        OGRLayer* ln = dn ? dn->GetLayer(0) : nullptr;
        check("OS: nodes file read back",
              ln && ln->GetFeatureCount() == static_cast<GIntBig>(r_os.network.nodes.size()));
    }

    // ---- local_file with the Open Rivers field names: the same lines --------
    lsdtt::lsdtt_river_line_fetch_options lf = base;
    lf.source = "local_file";
    lf.river_lines_path = gpkg;
    lf.river_lines_layer = "watercourse_link";
    lf.fields = lsdtt::os_open_rivers_field_map();
    const auto r_lf = lsdtt::fetch_river_lines(lf);
    std::string why;
    check("local_file = gb_os_open_rivers (lines)", same_lines(r_lf.network, r_os.network, 0.0, why), why);
    check("local_file = gb_os_open_rivers (node types)", type_counts(r_lf.network) == type_counts(r_os.network));

    // ---- reversed copy: same lines ------------------------------------------
    const std::string rev = (work / "reversed.gpkg").string();
    {
        std::unique_ptr<GDALDataset> in(static_cast<GDALDataset*>(
            GDALOpenEx(gpkg.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr)));
        OGRLayer* li = in->GetLayerByName("watercourse_link");
        GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GPKG");
        std::unique_ptr<GDALDataset> out(drv->Create(rev.c_str(), 0, 0, 0, GDT_Unknown, nullptr));
        OGRLayer* lo = out->CreateLayer("watercourse_link", li->GetSpatialRef(), wkbLineString, nullptr);
        for (int i = 0; i < li->GetLayerDefn()->GetFieldCount(); ++i)
            lo->CreateField(li->GetLayerDefn()->GetFieldDefn(i));
        for (auto& f : *li)
        {
            std::unique_ptr<OGRFeature, void (*)(OGRFeature*)> g(
                OGRFeature::CreateFeature(lo->GetLayerDefn()), OGRFeature::DestroyFeature);
            g->SetFrom(f.get());
            auto* ls = g->GetGeometryRef()->toLineString();
            ls->reversePoints();
            const std::string s = f->GetFieldAsString("start_node");
            const std::string en = f->GetFieldAsString("end_node");
            g->SetField("start_node", en.c_str());
            g->SetField("end_node", s.c_str());
            const bool opp = std::string(f->GetFieldAsString("flow_direction")) == "in opposite direction";
            g->SetField("flow_direction", opp ? "in direction" : "in opposite direction");
            if (lo->CreateFeature(g.get()) != OGRERR_NONE) throw std::runtime_error("could not write the reversed copy");
        }
    }
    lsdtt::lsdtt_river_line_fetch_options lr = lf;
    lr.river_lines_path = rev;
    const auto r_rev = lsdtt::fetch_river_lines(lr);
    check("reversed + flipped attribute = original (lines)", same_lines(r_rev.network, r_os.network, 0.0, why), why);
    check_int("reversed: every part reversed", r_rev.network.n_parts_reversed, ref.n_links);

    // ---- coordinates only: same node count ----------------------------------
    lsdtt::lsdtt_river_line_fetch_options lc = lf;
    lc.fields.start_node_field.clear();
    lc.fields.end_node_field.clear();
    const auto r_lc = lsdtt::fetch_river_lines(lc);
    check_int("no node ids: nodes from coordinates = id nodes",
              static_cast<long>(r_lc.network.nodes.size()), ref.n_nodes);

    // ---- class filter ---------------------------------------------------------
    lsdtt::lsdtt_river_line_fetch_options lk = lf;
    lk.keep_feature_classes = {"river"};
    const auto r_lk = lsdtt::fetch_river_lines(lk);
    const int n_lake = r_os.classes.at("lake_connector").n_read;
    check_int("filter: river only drops the lake links",
              static_cast<long>(r_lk.network.lines.size()), ref.n_links - n_lake);

    // ---- negative ---------------------------------------------------------------
    lsdtt::lsdtt_river_line_fetch_options bad = lf;
    bad.fields.feature_class_field = "nope";
    check("missing field refused, naming it",
          thrown_message([&] { lsdtt::fetch_river_lines(bad); }).find("'nope' is not a field") !=
              std::string::npos);
    lsdtt::lsdtt_river_line_fetch_options bad2 = base;
    bad2.source = "hydrorivers";
    check("unknown source refused",
          thrown_message([&] { lsdtt::fetch_river_lines(bad2); }).find("is not one of") != std::string::npos);
    lsdtt::lsdtt_river_line_fetch_options bad3 = os;
    bad3.output_epsg = 4326;
    check("geographic output CRS refused",
          thrown_message([&] { lsdtt::fetch_river_lines(bad3); }).find("geographic") != std::string::npos);

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
