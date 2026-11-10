// test_multipart_geometry.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Multipart geometry on read.
//
// Until 2026-08-26 `lsdtt_data_loader` had NO branch for any multipart type --
// not MultiPoint, MultiLineString, MultiPolygon or GeometryCollection. An
// unrecognised geometry returned `std::monostate` with no throw and no warning,
// and every consumer then did `get_if<lsdtt_line_geom>(...)`, got nullptr, and
// skipped the feature. So a multipart input produced an empty output and a ZERO
// EXIT CODE.
//
// Demonstrated before the fix, on shipped CLIs:
//   lsdtt-linear-stats            LineString -> 1 feature;  MultiLineString -> empty layer
//   lsdtt-raster-preprocessing    Point x3   -> 3 features; MultiPoint      -> 0-byte file
//
// Both exited 0. THE FAILURE MODE IS SILENCE, so every check here asserts a
// positive -- a part count, a coordinate, an attribute -- and never a clean
// return.
//
//   M1  MultiLineString explodes into one record per part, coordinates intact
//   M2  MultiPoint likewise
//   M3  MultiPolygon likewise, rings intact
//   M4  ATTRIBUTES ARE COPIED ONTO EVERY PART, so an exploded feature does not
//       lose the fields that identified it
//   5   `source_feature_index` recovers the original grouping, and is NOT a
//       schema field -- a synthetic column would make every read-then-write
//       round trip emit something the source did not have
//   M6  single-part geometry is completely unaffected
//   M7  GEOPACKAGE. NOTE the promotion is CONDITIONAL, not universal: a GPKG
//       layer whose content is all single-part is declared Point/LineString and
//       always read correctly. It is declared Multi -- and then every feature in
//       it reads back as Multi, single-part ones included -- only when the
//       content is mixed or already multipart
//   M8  SHAPEFILE, which preserves the distinction and so yields a genuinely
//       MIXED layer. Public datasets arrive in it constantly; lsdtt3 never
//       writes it, so this is an ingest-only concern -- and the pre-fix failure
//       was PARTIAL (single-part features survived, multi-part vanished), which
//       is worse than total loss because the output looks plausible
//
// Usage: lsdtt-test-multipart-geometry <out_dir>
// -----------------------------------------------------------------------------

#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"

// This test builds its vector fixtures with the OGR API directly, so it needs
// the GDAL headers in its own right. It used to get them transitively through
// lsdtt_data_loader.hpp; that header no longer includes them, which is the
// point of the change -- a file that uses GDAL now has to say so.
#include <gdal.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <ogr_spatialref.h>

namespace
{

int g_checks   = 0;
int g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

std::string g_dir;

std::string path_for(const std::string& name)
{
    return g_dir + "/" + name;
}

// Write one layer through OGR directly, so the test does not depend on the
// writer under test to build its own fixtures.
void write_wkt_layer(const std::string& file,
                     const std::string& layer,
                     OGRwkbGeometryType type,
                     const std::vector<std::string>& wkts,
                     const std::vector<std::int64_t>& ids,
                     const char* driver = "FlatGeobuf")
{
    GDALAllRegister();
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName(driver);
    if (!drv) { check(false, "FlatGeobuf driver available"); return; }

    GDALDataset* ds = drv->Create(file.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) { check(false, "created " + file); return; }

    OGRSpatialReference srs;
    srs.importFromEPSG(32630);
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRLayer* lyr = ds->CreateLayer(layer.c_str(), &srs, type, nullptr);
    OGRFieldDefn fd("id", OFTInteger64);
    lyr->CreateField(&fd);

    for (std::size_t i = 0; i < wkts.size(); ++i)
    {
        OGRFeature* f = OGRFeature::CreateFeature(lyr->GetLayerDefn());
        f->SetField("id", static_cast<GIntBig>(ids[i]));
        OGRGeometry* g = nullptr;
        OGRGeometryFactory::createFromWkt(wkts[i].c_str(), &srs, &g);
        f->SetGeometry(g);
        // Checked, not ignored: a fixture that silently fails to write a
        // feature would make every assertion below vacuous.
        if (lyr->CreateFeature(f) != OGRERR_NONE)
            check(false, "wrote feature to " + file);
        OGRGeometryFactory::destroyGeometry(g);
        OGRFeature::DestroyFeature(f);
    }
    GDALClose(ds);
}

lsdtt::lsdtt_vector_layer read_back(const std::string& file,
                                    const std::string& layer)
{
    lsdtt::lsdtt_data_loader loader;
    loader.open_vector(file);
    lsdtt::lsdtt_vector_layer out = loader.read_vector_layer(layer);
    loader.close_vector();
    return out;
}

}  // namespace


int main(int argc, char** argv)
{
    using namespace lsdtt;

    g_dir = (argc > 1) ? argv[1] : ".";

    std::cout << "=====================================================\n"
              << " Multipart geometry on read\n"
              << "=====================================================\n";

    // -------------------------------------------------------------------------
    std::cout << "\nM1: MultiLineString explodes into parts\n";
    {
        const std::string f = path_for("mls.fgb");
        write_wkt_layer(f, "mls", wkbMultiLineString,
                        {"MULTILINESTRING((0 0,10 10),(100 100,110 110,120 130))"},
                        {7});
        const auto layer = read_back(f, "mls");

        check(layer.records.size() == 2,
              "one feature -> 2 records (got " +
              std::to_string(layer.records.size()) + ")");

        if (layer.records.size() == 2)
        {
            const auto* a = std::get_if<lsdtt_line_geom>(&layer.records[0].geometry);
            const auto* b = std::get_if<lsdtt_line_geom>(&layer.records[1].geometry);
            check(a && b, "both records carry line geometry");
            if (a && b)
            {
                check(a->points.size() == 2 && b->points.size() == 3,
                      "part vertex counts are 2 and 3");
                // The parts must NOT have been joined: that would invent a
                // segment from (10,10) to (100,100), which is the artefact
                // multipart exists to avoid.
                check(a->points.back().x == 10.0 && b->points.front().x == 100.0,
                      "the parts are SEPARATE -- no joining segment invented");
            }
        }
    }

    // -------------------------------------------------------------------------
    std::cout << "\nM2: MultiPoint explodes into parts\n";
    {
        const std::string f = path_for("mp.fgb");
        write_wkt_layer(f, "mp", wkbMultiPoint,
                        {"MULTIPOINT((5 5),(25 25),(45 45))"}, {3});
        const auto layer = read_back(f, "mp");

        check(layer.records.size() == 3,
              "one feature -> 3 records (got " +
              std::to_string(layer.records.size()) + ")");
        int n_pt = 0;
        for (const auto& r : layer.records)
            if (std::holds_alternative<lsdtt_point_geom>(r.geometry)) ++n_pt;
        check(n_pt == 3, "all three carry point geometry");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nM3: MultiPolygon explodes into parts\n";
    {
        const std::string f = path_for("mpoly.fgb");
        write_wkt_layer(f, "mpoly", wkbMultiPolygon,
                        {"MULTIPOLYGON(((0 0,10 0,10 10,0 10,0 0)),"
                         "((50 50,60 50,60 60,50 60,50 50)))"}, {1});
        const auto layer = read_back(f, "mpoly");

        check(layer.records.size() == 2,
              "one feature -> 2 records (got " +
              std::to_string(layer.records.size()) + ")");
        const auto* p = layer.records.empty() ? nullptr
            : std::get_if<lsdtt_polygon_geom>(&layer.records[0].geometry);
        check(p != nullptr, "record carries polygon geometry");
        if (p) check(p->exterior.size() >= 4, "the exterior ring survived");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nM4: attributes are copied onto every part\n";
    {
        const auto layer = read_back(path_for("mp.fgb"), "mp");
        bool all_have_id = !layer.records.empty();
        for (const auto& r : layer.records)
        {
            if (r.field_values.empty()) { all_have_id = false; break; }
            const auto* v = std::get_if<std::int64_t>(&r.field_values[0]);
            if (!v || *v != 3) all_have_id = false;
        }
        check(all_have_id,
              "every exploded part kept id = 3 -- an exploded feature must not "
              "lose the fields that identified it");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nM5: source_feature_index, and it is NOT a schema field\n";
    {
        const auto layer = read_back(path_for("mp.fgb"), "mp");
        check(layer.source_feature_index.size() == layer.records.size(),
              "one index per record");
        bool all_zero = !layer.source_feature_index.empty();
        for (auto v : layer.source_feature_index) if (v != 0) all_zero = false;
        check(all_zero, "all three parts trace back to source feature 0");

        // A synthetic schema column would make every read-then-write round trip
        // emit something the source file did not have. unit_internal_drainage_-
        // detect caught exactly that when this was first written as a field.
        bool present = false;
        for (const auto& n : layer.schema.field_names)
            if (n == "source_feature_index") present = true;
        check(!present,
              "source_feature_index is NOT in the schema -- the schema must keep "
              "describing the file");
        check(layer.schema.field_names.size() == 1,
              "the schema still has exactly the file's own 1 field");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nM6: single-part geometry is unaffected\n";
    {
        const std::string f = path_for("ls.fgb");
        write_wkt_layer(f, "ls", wkbLineString,
                        {"LINESTRING(0 0,10 10,20 30)"}, {42});
        const auto layer = read_back(f, "ls");

        check(layer.records.size() == 1, "one feature -> one record");
        const auto* a = layer.records.empty() ? nullptr
            : std::get_if<lsdtt_line_geom>(&layer.records[0].geometry);
        check(a && a->points.size() == 3, "all three vertices present");
        check(layer.source_feature_index.size() == 1 &&
              layer.source_feature_index[0] == 0,
              "and it still records its source feature");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nM7: GEOPACKAGE -- every feature promoted to MultiLineString\n";
    {
        // Where an unguarded std::get used to turn a silent skip into a crash:
        // one file, one type single-part and one multipart.
        // GEOPACKAGE, mixed content.
        //
        // The promotion is CONDITIONAL and it is worth being precise, because
        // the loose version of this claim is wrong. Measured:
        //
        //   all single-part  -> layer declared Point / LineString  -> always read
        //   already multipart-> layer declared Multi*              -> dropped pre-fix
        //   MIXED            -> layer declared Multi*, and EVERY feature reads
        //                       back as Multi, single-part ones included
        //                       -> the whole layer dropped pre-fix
        //
        // So a GPKG only loses everything when its content is mixed or
        // multipart. A plain single-part GPKG was never affected, which is why
        // ordinary point layers from QGIS have always imported without trouble.
        const std::string f = path_for("mixed.gpkg");
        write_wkt_layer(f, "mixed", wkbUnknown,
                        {"LINESTRING(0 0,5 5)",
                         "MULTILINESTRING((10 10,20 20),(30 30,40 40))"},
                        {1, 2}, "GPKG");
        const auto layer = read_back(f, "mixed");

        check(layer.records.size() == 3,
              "2 features (1 single + 1 two-part) -> 3 records (got " +
              std::to_string(layer.records.size()) + ")");
        int n_line = 0;
        for (const auto& r : layer.records)
            if (std::holds_alternative<lsdtt_line_geom>(r.geometry)) ++n_line;
        check(n_line == 3, "every record carries usable line geometry");

        if (layer.source_feature_index.size() == 3)
            check(layer.source_feature_index[0] == 0 &&
                  layer.source_feature_index[1] == 1 &&
                  layer.source_feature_index[2] == 1,
                  "the grouping is 0, 1, 1");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nM8: SHAPEFILE -- a genuinely mixed layer\n";
    {
        // Shapefile has no LineString/MultiLineString distinction in the format:
        // a PolyLine with one part reads back as LineString, one with several as
        // MultiLineString, IN THE SAME FILE. It is the ONLY one of the four
        // formats tested that preserves the distinction -- GPKG, FlatGeobuf and
        // GeoJSON all promote singles to Multi in a mixed layer -- so it is the
        // one that produces a genuinely mixed layer. Public datasets arrive in
        // it constantly.
        //
        // lsdtt3 never WRITES shapefile -- this is purely an ingest concern.
        //
        // The pre-fix failure here was PARTIAL: the single-part feature survived
        // and the multi-part one vanished, so the output looked plausible. That
        // is harder to notice than the total loss GPKG produced.
        const std::string f = path_for("mixed.shp");
        write_wkt_layer(f, "mixed", wkbLineString,
                        {"LINESTRING(0 0,5 5)",
                         "MULTILINESTRING((10 10,20 20),(30 30,40 40))"},
                        {1, 2}, "ESRI Shapefile");
        const auto layer = read_back(f, "mixed");

        check(layer.records.size() == 3,
              "2 shapefile features (1 single + 1 two-part) -> 3 records (got " +
              std::to_string(layer.records.size()) + ")");
        int n_line = 0;
        for (const auto& r : layer.records)
            if (std::holds_alternative<lsdtt_line_geom>(r.geometry)) ++n_line;
        check(n_line == 3,
              "every record carries usable line geometry -- pre-fix only the "
              "single-part one did, and the output still looked plausible");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nM9: features with no geometry are kept AND reported\n";
    {
        // A feature with a null geometry is kept as a record (so the feature
        // count and attributes survive) and must be reported, like the
        // multipart and unsupported cases: before, it was counted (n_empty)
        // but the count was never printed, so it vanished from every consumer
        // without a word.
        const std::string f = path_for("empty_geom.gpkg");
        write_wkt_layer(f, "empty_geom", wkbLineString,
                        {"LINESTRING(0 0,5 5)", ""},   // "" -> null geometry
                        {1, 2}, "GPKG");

        std::ostringstream captured;
        std::streambuf* old = std::cout.rdbuf(captured.rdbuf());
        const auto layer = read_back(f, "empty_geom");
        std::cout.rdbuf(old);
        const std::string out = captured.str();

        check(layer.records.size() == 2,
              "both features kept as records (got " +
              std::to_string(layer.records.size()) + ")");
        check(out.find("WARNING - 1 feature(s) in layer 'empty_geom' have no geometry")
                  != std::string::npos,
              "the empty feature is reported (loader printed: [" + out + "])");
        check(out.find("geometry type this reader does not support") == std::string::npos,
              "and is not mis-reported as an unsupported geometry type");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
