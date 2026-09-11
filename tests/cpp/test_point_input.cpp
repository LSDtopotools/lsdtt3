// test_point_input.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Unit tests for src/core/lsdtt_point_input
// (docs/plan_docs/internal_drainage_plan.md, sections 6.3 and 6.4).
//
// Phase 2a factored the point-input module into one CSV parser and one OGR read
// loop, added load_point_layer_with_attributes(), and rewrote
// load_drainage_injection_points() as a thin wrapper over it.  That function has
// two live callers (src/cpu/lsdtt_flowinfo_cpu.cpp, src/cpu/lsdtt_flowgraph_cpu.cpp)
// and the module had no CTest coverage at all, so the refactor went in
// unguarded.  This driver is the guard.
//
// What is under test
// ------------------
//   load_drainage_injection_points()
//     the exact text of the missing-field and open-failure messages, which are
//     reconstructed deliberately by the wrapper and are therefore the contract;
//     silent skipping of non-point / unmappable / non-finite / negative /
//     absent records; point_id continuing to count skipped features
//
//   load_point_layer_with_attributes()
//     cross-format equivalence (.fgb / .gpkg / .csv), fields_missing yielding
//     kind::absent rather than throwing, empty field_names reading every field,
//     CSV latitude/longitude being geometry and not attributes, order
//     preservation, per-record status, dataset-level throws
//
//   lsdtt_point_attribute_value accessors
//     every accepted as_bool spelling, and failure-by-return-value rather than
//     by exception on junk input
//
//   the deliberate divergences of plan section 6.4
//     6.4(1) a non-numeric string-typed drainage_area is now SKIPPED
//     6.4(4) attribute kind is not identical across formats
//
// All fixtures are written at run time into the directory given as argv[1]
// (default ./point_input_test_output), so nothing is committed to tests/data.
//
// Usage: lsdtt-test-point-input [output_dir]
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

#include "lsdtt_point_input.hpp"
#include "lsdtt_raster_meta.hpp"

namespace
{

// =============================================================================
// Tiny check harness (same shape as tests/cpp/test_internal_drainage.cpp)
// =============================================================================
int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

// Runs f and requires that it threw std::runtime_error whose message is
// EXACTLY expected_message.  The text is asserted, not just the throw: the
// wrapper in load_drainage_injection_points() goes out of its way to rebuild
// these strings, so the strings are the contract.
template <typename F>
void check_throws_exact(const std::string& name,
                        const F&           f,
                        const std::string& expected_message)
{
    bool        threw   = false;
    bool        matched = false;
    std::string what;

    try
    {
        f();
    }
    catch (const std::runtime_error& e)
    {
        threw   = true;
        what    = e.what();
        matched = (what == expected_message);
    }
    catch (const std::exception& e)
    {
        what = std::string("wrong exception type: ") + e.what();
    }

    check(name, threw && matched,
          threw ? ("got: \"" + what + "\"") : ("did not throw; " + what));
}

bool close_enough(double a, double b, double tol = 1.0e-9)
{
    return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b));
}

// =============================================================================
// Test raster and coordinate helpers
//
// The raster is a small projected (UTM 33N) grid.  The vector and CSV fixtures
// are written in EPSG:4326, so every read path exercises a real coordinate
// transformation into the raster CRS rather than an identity.
// =============================================================================
lsdtt::lsdtt_raster_meta make_test_raster_meta()
{
    lsdtt::lsdtt_raster_meta meta;

    meta.n_rows    = 100;
    meta.n_cols    = 100;
    meta.x_min     = 500000.0;
    meta.y_min     = 4000000.0;
    meta.cell_size = 10.0;
    meta.epsg_code = 32633;

    OGRSpatialReference srs;
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (srs.importFromEPSG(32633) != OGRERR_NONE)
    {
        throw std::runtime_error("test_point_input: could not build EPSG:32633");
    }

    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    meta.crs_wkt = (wkt != nullptr) ? std::string(wkt) : std::string();
    CPLFree(wkt);

    return meta;
}

struct lon_lat
{
    double lon = 0.0;
    double lat = 0.0;
};

// Converts a raster-CRS coordinate into EPSG:4326 lon/lat, so the fixtures can
// be specified in the pixel coordinates whose row/col the test asserts.
lon_lat raster_xy_to_lon_lat(const lsdtt::lsdtt_raster_meta& meta,
                             double x,
                             double y)
{
    OGRSpatialReference src_srs;
    OGRSpatialReference dst_srs;

    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (src_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE ||
        dst_srs.importFromEPSG(4326) != OGRERR_NONE)
    {
        throw std::runtime_error("test_point_input: could not build CRS pair");
    }

    OGRCoordinateTransformation* tx =
        OGRCreateCoordinateTransformation(&src_srs, &dst_srs);

    if (tx == nullptr)
    {
        throw std::runtime_error("test_point_input: could not build transform");
    }

    double xx = x;
    double yy = y;
    const int ok = tx->Transform(1, &xx, &yy);
    OCTDestroyCoordinateTransformation(tx);

    if (ok == 0)
    {
        throw std::runtime_error("test_point_input: transform to 4326 failed");
    }

    return lon_lat{xx, yy};
}

// Pixel centres of the grid above.  Chosen so the expected row/col/flat can be
// written down by hand and a transform round trip cannot move them.
//   A -> row 4,  col 5,  flat 405
//   B -> row 12, col 12, flat 1212
//   C -> row 20, col 20, flat 2020
//   D -> row 30, col 30, flat 3030
//   E -> row 40, col 40, flat 4040
const double k_ax = 500055.0, k_ay = 4000955.0;
const double k_bx = 500125.0, k_by = 4000875.0;
const double k_cx = 500205.0, k_cy = 4000795.0;
const double k_dx = 500305.0, k_dy = 4000695.0;
const double k_ex = 500405.0, k_ey = 4000595.0;

// Comfortably outside the 1 km x 1 km grid.
const double k_outx = 600000.0, k_outy = 4100000.0;

const std::int64_t k_a_flat = 405;
const std::int64_t k_b_flat = 1212;
const std::int64_t k_c_flat = 2020;
const std::int64_t k_e_flat = 4040;

// =============================================================================
// Fixture writers
// =============================================================================
GDALDataset* create_vector_dataset(const std::string& driver_name,
                                   const std::string& path)
{
    std::filesystem::remove_all(path);

    GDALDriver* driver =
        GetGDALDriverManager()->GetDriverByName(driver_name.c_str());

    if (driver == nullptr)
    {
        throw std::runtime_error(
            "test_point_input: GDAL driver not available: " + driver_name);
    }

    GDALDataset* ds =
        driver->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);

    if (ds == nullptr)
    {
        throw std::runtime_error(
            "test_point_input: could not create dataset: " + path);
    }

    return ds;
}

void set_point_geometry(OGRFeature*                     feature,
                        const lsdtt::lsdtt_raster_meta& meta,
                        double                          x,
                        double                          y)
{
    const lon_lat ll = raster_xy_to_lon_lat(meta, x, y);
    OGRPoint point(ll.lon, ll.lat);
    feature->SetGeometry(&point);
}

// The mixed injection layer.  Feature order matters: only features 1 and 7 are
// usable, so point_id 7 on the second returned record is what proves that
// point_id counts skipped features.
void write_injection_layer(const lsdtt::lsdtt_raster_meta& meta,
                           const std::string&              path)
{
    GDALDataset* ds = create_vector_dataset("GPKG", path);

    OGRSpatialReference srs;
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    srs.importFromEPSG(4326);

    // wkbUnknown so a non-point feature can live in the same layer.
    OGRLayer* layer = ds->CreateLayer("sinks", &srs, wkbUnknown, nullptr);
    if (layer == nullptr)
    {
        GDALClose(ds);
        throw std::runtime_error("test_point_input: could not create layer");
    }

    OGRFieldDefn area_defn("drainage_area", OFTReal);
    layer->CreateField(&area_defn);

    OGRFieldDefn label_defn("label", OFTString);
    layer->CreateField(&label_defn);

    struct injection_fixture_feature
    {
        double      x;
        double      y;
        bool        is_line;
        bool        set_area;
        double      area;
        const char* label;
    };

    const injection_fixture_feature features[] = {
        // 1: usable
        {k_ax,   k_ay,   false, true,  1000.0,                            "valid_a"},
        // 2: a line, not a point -> skipped
        {k_bx,   k_by,   true,  true,  2000.0,                            "line"},
        // 3: outside the raster -> skipped
        {k_outx, k_outy, false, true,  3000.0,                           "outside"},
        // 4: negative area -> skipped
        {k_bx,   k_by,   false, true,  -5.0,                             "negative"},
        // 5: non-finite area -> skipped (GPKG may also store this as NULL,
        //    which is skipped for the same reason: no usable double)
        {k_cx,   k_cy,   false, true,  std::numeric_limits<double>::quiet_NaN(),
                                                                          "nan"},
        // 6: field left unset -> absent -> skipped
        {k_dx,   k_dy,   false, false, 0.0,                              "null"},
        // 7: usable
        {k_ex,   k_ey,   false, true,  4000.0,                            "valid_e"}
    };

    for (const injection_fixture_feature& spec : features)
    {
        OGRFeature* feature = OGRFeature::CreateFeature(layer->GetLayerDefn());

        if (spec.set_area)
        {
            feature->SetField("drainage_area", spec.area);
        }
        feature->SetField("label", spec.label);

        if (spec.is_line)
        {
            const lon_lat p0 = raster_xy_to_lon_lat(meta, spec.x, spec.y);
            const lon_lat p1 = raster_xy_to_lon_lat(meta, spec.x + 20.0,
                                                    spec.y + 20.0);
            OGRLineString line;
            line.addPoint(p0.lon, p0.lat);
            line.addPoint(p1.lon, p1.lat);
            feature->SetGeometry(&line);
        }
        else
        {
            set_point_geometry(feature, meta, spec.x, spec.y);
        }

        const OGRErr create_error = layer->CreateFeature(feature);
        OGRFeature::DestroyFeature(feature);

        if (create_error != OGRERR_NONE)
        {
            GDALClose(ds);
            throw std::runtime_error(
                "test_point_input: could not write fixture feature");
        }
    }

    GDALClose(ds);
}

// A layer whose drainage_area is a STRING field.  Feature 1 holds text that is
// not a number; feature 2 holds numeric text.  This is the plan 6.4(1) case.
void write_text_area_layer(const lsdtt::lsdtt_raster_meta& meta,
                           const std::string&              path)
{
    GDALDataset* ds = create_vector_dataset("GPKG", path);

    OGRSpatialReference srs;
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    srs.importFromEPSG(4326);

    OGRLayer* layer = ds->CreateLayer("sinks", &srs, wkbPoint, nullptr);
    if (layer == nullptr)
    {
        GDALClose(ds);
        throw std::runtime_error("test_point_input: could not create layer");
    }

    OGRFieldDefn area_defn("drainage_area", OFTString);
    layer->CreateField(&area_defn);

    struct text_fixture_feature
    {
        double      x;
        double      y;
        const char* area_text;
    };

    const text_fixture_feature features[] = {
        {k_ax, k_ay, "not_a_number"},
        {k_bx, k_by, "2500.5"}
    };

    for (const text_fixture_feature& spec : features)
    {
        OGRFeature* feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
        feature->SetField("drainage_area", spec.area_text);
        set_point_geometry(feature, meta, spec.x, spec.y);
        const OGRErr create_error = layer->CreateFeature(feature);
        OGRFeature::DestroyFeature(feature);

        if (create_error != OGRERR_NONE)
        {
            GDALClose(ds);
            throw std::runtime_error(
                "test_point_input: could not write fixture feature");
        }
    }

    GDALClose(ds);
}

// Three points carrying one integer, one real and one text attribute, written
// through OGR.  The same points and values are written as CSV by
// write_attribute_csv(), and the two must read back equivalently.
void write_attribute_layer(const lsdtt::lsdtt_raster_meta& meta,
                           const std::string&              driver_name,
                           const std::string&              path)
{
    GDALDataset* ds = create_vector_dataset(driver_name, path);

    OGRSpatialReference srs;
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    srs.importFromEPSG(4326);

    // SPATIAL_INDEX=NO for FlatGeobuf: the driver's default packed Hilbert
    // R-tree REORDERS features on write, which would make a cross-format order
    // comparison meaningless.  The reader takes features in layer order, so
    // this keeps layer order equal to the order the fixture was written in.
    const char* const layer_options[] = {"SPATIAL_INDEX=NO", nullptr};

    OGRLayer* layer = ds->CreateLayer(
        "sink_points", &srs, wkbPoint,
        const_cast<char**>(layer_options));
    if (layer == nullptr)
    {
        GDALClose(ds);
        throw std::runtime_error("test_point_input: could not create layer");
    }

    OGRFieldDefn id_defn("sink_id", OFTInteger64);
    layer->CreateField(&id_defn);

    OGRFieldDefn area_defn("area_m2", OFTReal);
    layer->CreateField(&area_defn);

    OGRFieldDefn label_defn("label", OFTString);
    layer->CreateField(&label_defn);

    struct attribute_fixture_feature
    {
        double       x;
        double       y;
        std::int64_t sink_id;
        double       area_m2;
        const char*  label;
    };

    const attribute_fixture_feature features[] = {
        {k_ax, k_ay, 7, 1234.5, "alpha"},
        // area_m2 = 0 is the plan 6.4(4) case: whole numbers come back
        // kind::integer from CSV and kind::real from OGR.
        {k_bx, k_by, 8,    0.0, "beta"},
        {k_cx, k_cy, 9,  99.25, "gamma"}
    };

    for (const attribute_fixture_feature& spec : features)
    {
        OGRFeature* feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
        feature->SetField("sink_id", static_cast<GIntBig>(spec.sink_id));
        feature->SetField("area_m2", spec.area_m2);
        feature->SetField("label", spec.label);
        set_point_geometry(feature, meta, spec.x, spec.y);
        const OGRErr create_error = layer->CreateFeature(feature);
        OGRFeature::DestroyFeature(feature);

        if (create_error != OGRERR_NONE)
        {
            GDALClose(ds);
            throw std::runtime_error(
                "test_point_input: could not write fixture feature");
        }
    }

    GDALClose(ds);
}

void write_attribute_csv(const lsdtt::lsdtt_raster_meta& meta,
                         const std::string&              path)
{
    const lon_lat a = raster_xy_to_lon_lat(meta, k_ax, k_ay);
    const lon_lat b = raster_xy_to_lon_lat(meta, k_bx, k_by);
    const lon_lat c = raster_xy_to_lon_lat(meta, k_cx, k_cy);

    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("test_point_input: could not write " + path);
    }

    out << std::setprecision(15);
    out << "latitude,longitude,sink_id,area_m2,label\n";
    out << a.lat << "," << a.lon << ",007,1234.5,alpha\n";
    out << b.lat << "," << b.lon << ",8,0,beta\n";
    out << c.lat << "," << c.lon << ",9,99.25,gamma\n";

    out.close();
}

// =============================================================================
// A. The load_drainage_injection_points() contract
// =============================================================================

// A missing drainage_area_field is a DATASET-level error, and the message is
// the contract: it names the field, the resolved layer and the dataset so a
// user can fix their file.  Asserted exactly, because phase 2a rebuilds this
// string by hand inside the wrapper.
void test_injection_missing_field_message(const lsdtt::lsdtt_raster_meta& meta,
                                          const std::string& injection_path)
{
    check_throws_exact(
        "injection: missing field throws with the exact contract message",
        [&]() {
            (void)lsdtt::load_drainage_injection_points(
                injection_path, meta, "no_such_field");
        },
        "load_drainage_injection_points: field 'no_such_field' not found in "
        "layer 'sinks' of dataset: " + injection_path);
}

// An unopenable dataset must surface this function's OWN prefix, not the
// generic lsdtt_point_input: one it now delegates to.
void test_injection_open_failure_message(const lsdtt::lsdtt_raster_meta& meta,
                                         const std::string& missing_path)
{
    check_throws_exact(
        "injection: open failure keeps the load_drainage_injection_points prefix",
        [&]() {
            (void)lsdtt::load_drainage_injection_points(missing_path, meta);
        },
        "load_drainage_injection_points: could not open vector dataset: " +
        missing_path);
}

// Non-point, unmappable, non-finite, negative and absent records are dropped
// with no diagnostic; point_id keeps counting the features that were dropped.
void test_injection_skips_and_point_ids(const lsdtt::lsdtt_raster_meta& meta,
                                        const std::string& injection_path)
{
    const std::vector<lsdtt::lsdtt_drainage_injection_record> records =
        lsdtt::load_drainage_injection_points(injection_path, meta,
                                              "drainage_area");

    check("injection: only the two usable features are returned",
          records.size() == 2,
          "n_records = " + std::to_string(records.size()));

    if (records.size() != 2)
    {
        return;
    }

    check("injection: point_id counts skipped features (1 and 7, not 1 and 2)",
          records[0].point_id == 1 && records[1].point_id == 7,
          "point_ids = " + std::to_string(records[0].point_id) + ", " +
          std::to_string(records[1].point_id));

    check("injection: flats are the mapped pixels of the usable features",
          records[0].flat == k_a_flat && records[1].flat == k_e_flat,
          "flats = " + std::to_string(records[0].flat) + ", " +
          std::to_string(records[1].flat));

    check("injection: drainage areas are read from the named field",
          close_enough(records[0].drainage_area_m2, 1000.0) &&
          close_enough(records[1].drainage_area_m2, 4000.0),
          "areas = " + std::to_string(records[0].drainage_area_m2) + ", " +
          std::to_string(records[1].drainage_area_m2));
}

// -----------------------------------------------------------------------------
// D. The pinned divergence, plan section 6.4(1)
//
// DELIBERATE BEHAVIOUR CHANGE, DO NOT "FIX" THIS BACK.  The pre-phase-2a
// implementation read the drainage area with OGRFeature::GetFieldAsDouble(),
// so a string field holding non-numeric text silently became 0.0 and produced
// a record.  The typed accessor now refuses to convert it and the feature is
// skipped instead.  Injecting 0.0 was a no-op for both callers, so only the
// record count changes -- but it is a change, and this test pins it.
// -----------------------------------------------------------------------------
void test_injection_non_numeric_string_area_is_skipped(
    const lsdtt::lsdtt_raster_meta& meta,
    const std::string&              text_area_path)
{
    const std::vector<lsdtt::lsdtt_drainage_injection_record> records =
        lsdtt::load_drainage_injection_points(text_area_path, meta,
                                              "drainage_area");

    check("injection 6.4(1): non-numeric text in a string area field is skipped, "
          "not injected as 0.0",
          records.size() == 1,
          "n_records = " + std::to_string(records.size()));

    if (records.size() != 1)
    {
        return;
    }

    check("injection: numeric text in a string area field still parses",
          records[0].point_id == 2 &&
          records[0].flat == k_b_flat &&
          close_enough(records[0].drainage_area_m2, 2500.5),
          "point_id = " + std::to_string(records[0].point_id) +
          ", flat = " + std::to_string(records[0].flat) +
          ", area = " + std::to_string(records[0].drainage_area_m2));
}

// =============================================================================
// B. load_point_layer_with_attributes()
// =============================================================================

// Pulls the three attributes through the accessors.  Values are NEVER compared
// via kind: plan 6.4(4) says a whole number is kind::integer from CSV and
// kind::real from OGR by design, so comparing type across formats would be
// asserting the wrong thing.
struct attribute_triplet
{
    std::int64_t sink_id = 0;
    double       area_m2 = 0.0;
    std::string  label   = "";
    bool         ok      = false;
};

attribute_triplet read_triplet(
    const lsdtt::lsdtt_input_point_with_attributes& record)
{
    attribute_triplet out;

    const auto id_it    = record.attributes.find("sink_id");
    const auto area_it  = record.attributes.find("area_m2");
    const auto label_it = record.attributes.find("label");

    if (id_it    == record.attributes.end() ||
        area_it  == record.attributes.end() ||
        label_it == record.attributes.end())
    {
        return out;
    }

    out.ok = id_it->second.as_int64(out.sink_id) &&
             area_it->second.as_double(out.area_m2) &&
             label_it->second.as_string(out.label);

    return out;
}

bool triplets_match(const attribute_triplet& a, const attribute_triplet& b)
{
    return a.ok && b.ok &&
           a.sink_id == b.sink_id &&
           close_enough(a.area_m2, b.area_m2) &&
           a.label == b.label;
}

// Reading the same three points and attributes from .fgb, .gpkg and .csv must
// give equivalent records: same order, same mapped pixels, same attribute
// names, same attribute VALUES through the accessors.
void test_cross_format_equivalence(const lsdtt::lsdtt_raster_meta& meta,
                                   const std::string&              fgb_path,
                                   const std::string&              gpkg_path,
                                   const std::string&              csv_path)
{
    const std::vector<std::string> all_fields;

    const lsdtt::lsdtt_point_attribute_read_result from_fgb =
        lsdtt::load_point_layer_with_attributes(fgb_path, meta, all_fields);
    const lsdtt::lsdtt_point_attribute_read_result from_gpkg =
        lsdtt::load_point_layer_with_attributes(gpkg_path, meta, all_fields);
    const lsdtt::lsdtt_point_attribute_read_result from_csv =
        lsdtt::load_point_layer_with_attributes(csv_path, meta, all_fields);

    check("attributes: all three formats return three records",
          from_fgb.records.size() == 3 &&
          from_gpkg.records.size() == 3 &&
          from_csv.records.size() == 3,
          "fgb " + std::to_string(from_fgb.records.size()) +
          ", gpkg " + std::to_string(from_gpkg.records.size()) +
          ", csv " + std::to_string(from_csv.records.size()));

    if (from_fgb.records.size() != 3 ||
        from_gpkg.records.size() != 3 ||
        from_csv.records.size() != 3)
    {
        return;
    }

    // Empty field_names reads every field in the layer.  For the CSV this also
    // proves plan 6.4(3): latitude and longitude are geometry, not attributes,
    // so the CSV attribute set is the same size as the OGR one.
    bool schema_ok = true;
    for (std::size_t i = 0; i < 3; ++i)
    {
        schema_ok = schema_ok &&
                    from_fgb.records[i].attributes.size() == 3 &&
                    from_gpkg.records[i].attributes.size() == 3 &&
                    from_csv.records[i].attributes.size() == 3;
    }

    check("attributes: empty field_names reads every field in the layer",
          schema_ok,
          "attribute-map sizes fgb/gpkg/csv = " +
          std::to_string(from_fgb.records[0].attributes.size()) + "/" +
          std::to_string(from_gpkg.records[0].attributes.size()) + "/" +
          std::to_string(from_csv.records[0].attributes.size()));

    check("attributes 6.4(3): CSV latitude/longitude are geometry, so the CSV "
          "attribute names match the OGR ones",
          from_csv.records[0].attributes.count("sink_id") == 1 &&
          from_csv.records[0].attributes.count("area_m2") == 1 &&
          from_csv.records[0].attributes.count("label") == 1 &&
          from_csv.records[0].attributes.count("latitude") == 0 &&
          from_csv.records[0].attributes.count("longitude") == 0);

    const std::int64_t expected_flats[3] = {k_a_flat, k_b_flat, k_c_flat};

    bool geometry_ok = true;
    for (std::size_t i = 0; i < 3; ++i)
    {
        geometry_ok = geometry_ok &&
            from_fgb.records[i].point.status ==
                lsdtt::lsdtt_point_input_status::success &&
            from_gpkg.records[i].point.status ==
                lsdtt::lsdtt_point_input_status::success &&
            from_csv.records[i].point.status ==
                lsdtt::lsdtt_point_input_status::success &&
            from_fgb.records[i].point.flat  == expected_flats[i] &&
            from_gpkg.records[i].point.flat == expected_flats[i] &&
            from_csv.records[i].point.flat  == expected_flats[i] &&
            from_fgb.records[i].point.point_id == static_cast<int>(i) + 1 &&
            from_gpkg.records[i].point.point_id == static_cast<int>(i) + 1 &&
            from_csv.records[i].point.point_id == static_cast<int>(i) + 1;
    }

    check("attributes: all three formats map to the same pixels, in input order",
          geometry_ok);

    bool values_ok = true;
    for (std::size_t i = 0; i < 3; ++i)
    {
        const attribute_triplet fgb  = read_triplet(from_fgb.records[i]);
        const attribute_triplet gpkg = read_triplet(from_gpkg.records[i]);
        const attribute_triplet csv  = read_triplet(from_csv.records[i]);

        values_ok = values_ok &&
                    triplets_match(fgb, gpkg) &&
                    triplets_match(fgb, csv);
    }

    const attribute_triplet first = read_triplet(from_fgb.records[0]);

    check("attributes: values agree across .fgb / .gpkg / .csv via the accessors",
          values_ok,
          "first record: sink_id = " + std::to_string(first.sink_id) +
          ", area_m2 = " + std::to_string(first.area_m2) +
          ", label = " + first.label);

    check("attributes: the written values are the values read back",
          first.sink_id == 7 && close_enough(first.area_m2, 1234.5) &&
          first.label == "alpha");

    std::string csv_source_text;
    check("attributes: CSV numeric source spelling is retained for identifiers",
          from_csv.records[0].attributes.at("sink_id").as_source_string(
              csv_source_text) && csv_source_text == "007",
          "source text = " + csv_source_text);

    // Plan 6.4(4), pinned so nobody starts comparing kind across formats: the
    // whole number 0 is typed by inspection from CSV (integer) and by schema
    // from OGR (real).  as_double bridges them; as_int64 accepts both because
    // the real is exactly integral.
    const lsdtt::lsdtt_point_attribute_value& csv_zero =
        from_csv.records[1].attributes.at("area_m2");
    const lsdtt::lsdtt_point_attribute_value& ogr_zero =
        from_fgb.records[1].attributes.at("area_m2");

    double csv_zero_value = -1.0;
    double ogr_zero_value = -1.0;

    check("attributes 6.4(4): a whole number is kind::integer from CSV and "
          "kind::real from OGR, and as_double bridges them",
          csv_zero.type == lsdtt::lsdtt_point_attribute_value::kind::integer &&
          ogr_zero.type == lsdtt::lsdtt_point_attribute_value::kind::real &&
          csv_zero.as_double(csv_zero_value) &&
          ogr_zero.as_double(ogr_zero_value) &&
          csv_zero_value == 0.0 && ogr_zero_value == 0.0);
}

// A requested field that is not in the schema must be reported in
// fields_missing and appear as kind::absent in every record.  Phase 2b relies
// on this for its optional "accepted" column, so it must NOT throw.
void test_missing_field_is_absent_not_a_throw(
    const lsdtt::lsdtt_raster_meta& meta,
    const std::string&              fgb_path,
    const std::string&              csv_path)
{
    const std::vector<std::string> requested = {"sink_id", "accepted"};

    bool threw = false;
    lsdtt::lsdtt_point_attribute_read_result from_fgb;
    lsdtt::lsdtt_point_attribute_read_result from_csv;

    try
    {
        from_fgb = lsdtt::load_point_layer_with_attributes(fgb_path, meta,
                                                           requested);
        from_csv = lsdtt::load_point_layer_with_attributes(csv_path, meta,
                                                           requested);
    }
    catch (const std::exception& e)
    {
        threw = true;
        check("fields_missing: an absent field must not throw", false, e.what());
    }

    if (threw)
    {
        return;
    }

    check("fields_missing: an absent field does not throw", true);

    check("fields_missing: the absent field is reported, the present one is not",
          from_fgb.fields_missing.size() == 1 &&
          from_fgb.fields_missing[0] == "accepted" &&
          from_fgb.fields_found.size() == 1 &&
          from_fgb.fields_found[0] == "sink_id" &&
          from_csv.fields_missing.size() == 1 &&
          from_csv.fields_missing[0] == "accepted" &&
          from_csv.fields_found.size() == 1 &&
          from_csv.fields_found[0] == "sink_id");

    bool absent_everywhere = (from_fgb.records.size() == 3 &&
                              from_csv.records.size() == 3);

    for (const lsdtt::lsdtt_input_point_with_attributes& record :
         from_fgb.records)
    {
        const auto found = record.attributes.find("accepted");
        absent_everywhere = absent_everywhere &&
            found != record.attributes.end() &&
            found->second.type ==
                lsdtt::lsdtt_point_attribute_value::kind::absent &&
            !found->second.is_present();
    }

    for (const lsdtt::lsdtt_input_point_with_attributes& record :
         from_csv.records)
    {
        const auto found = record.attributes.find("accepted");
        absent_everywhere = absent_everywhere &&
            found != record.attributes.end() &&
            found->second.type ==
                lsdtt::lsdtt_point_attribute_value::kind::absent;
    }

    check("fields_missing: the absent field is kind::absent in every record",
          absent_everywhere);

    // The requested field that IS present must still carry its value.
    std::int64_t sink_id = 0;
    const bool sibling_ok =
        from_fgb.records.size() == 3 &&
        from_fgb.records[0].attributes.at("sink_id").as_int64(sink_id) &&
        sink_id == 7;

    check("fields_missing: a present sibling field is unaffected",
          sibling_ok,
          "sink_id = " + std::to_string(sink_id));
}

// Order is preserved with one record per feature, per-point problems land in
// record.point.status, and only dataset-level problems throw.
void test_order_status_and_dataset_throws(const lsdtt::lsdtt_raster_meta& meta,
                                          const std::string& injection_path,
                                          const std::string& output_dir)
{
    const lsdtt::lsdtt_point_attribute_read_result result =
        lsdtt::load_point_layer_with_attributes(
            injection_path, meta, std::vector<std::string>{});

    check("order: one record per feature, including the unusable ones",
          result.records.size() == 7,
          "n_records = " + std::to_string(result.records.size()));

    if (result.records.size() != 7)
    {
        return;
    }

    const lsdtt::lsdtt_point_input_status expected[7] = {
        lsdtt::lsdtt_point_input_status::success,
        lsdtt::lsdtt_point_input_status::unsupported_geometry_type,
        lsdtt::lsdtt_point_input_status::outside_raster_extent,
        lsdtt::lsdtt_point_input_status::success,
        lsdtt::lsdtt_point_input_status::success,
        lsdtt::lsdtt_point_input_status::success,
        lsdtt::lsdtt_point_input_status::success
    };

    bool status_ok = true;
    bool order_ok  = true;
    for (std::size_t i = 0; i < 7; ++i)
    {
        status_ok = status_ok && result.records[i].point.status == expected[i];
        order_ok  = order_ok &&
                    result.records[i].point.point_id == static_cast<int>(i) + 1;
    }

    check("status: per-point problems are recorded, not thrown", status_ok);
    check("order: point_id is 1-based and follows input order", order_ok);

    check("layer: the resolved layer name is reported",
          result.resolved_layer_name == "sinks",
          "resolved_layer_name = " + result.resolved_layer_name);

    // The unusable features still carry their attributes, which is what lets a
    // reviewer round-trip a file without losing rows.
    std::string line_label;
    const bool line_label_ok =
        result.records[1].attributes.at("label").as_string(line_label) &&
        line_label == "line";

    check("attributes: an unusable feature still carries its attributes",
          line_label_ok,
          "label = " + line_label);

    const std::string missing_vector = output_dir + "/does_not_exist.gpkg";
    check_throws_exact(
        "dataset: an unopenable vector dataset throws",
        [&]() {
            (void)lsdtt::load_point_layer_with_attributes(
                missing_vector, meta, std::vector<std::string>{});
        },
        "lsdtt_point_input: could not open vector dataset: " + missing_vector);

    // A CSV without the geometry columns is a dataset-level failure too.
    const std::string bad_csv = output_dir + "/no_lat_lon.csv";
    {
        std::ofstream out(bad_csv);
        out << "sink_id,area_m2\n1,10.0\n";
    }

    // The message names BOTH accepted positionings and the offending file, so a
    // user whose CSV is in the raster's CRS is told what to do rather than being
    // pushed into converting a file that was already usable.
    check_throws_exact(
        "dataset: a CSV with neither latitude/longitude nor easting/northing throws",
        [&]() {
            (void)lsdtt::load_point_layer_with_attributes(
                bad_csv, meta, std::vector<std::string>{});
        },
        "lsdtt_point_input: CSV input must contain columns named 'latitude' "
        "and 'longitude', or 'easting' and 'northing' (or 'x' and 'y') if the "
        "points are already in the raster's coordinate system: " + bad_csv);
}

// =============================================================================
// A CSV positioned in the RASTER'S OWN CRS
//
// Projected point files are at least as common as geographic ones: anything
// exported from a GIS beside the DEM carries easting and northing. Requiring
// latitude/longitude made users convert a file that was already in the right
// coordinate system -- and the OGR path had accepted easting/northing all along,
// so CSV was the odd one out.
//
// The three checks below are the contract: projected columns position a point
// identically to the lat/long form, the aliases x/y work, and latitude/longitude
// still wins when a file carries both.
// =============================================================================
void test_csv_in_raster_crs(const lsdtt::lsdtt_raster_meta& meta,
                            const std::string&              dir)
{
    const lon_lat a = raster_xy_to_lon_lat(meta, k_ax, k_ay);

    // Same point, written three ways.
    const std::string ll_csv   = dir + "/crs_latlon.csv";
    const std::string en_csv   = dir + "/crs_easting.csv";
    const std::string xy_csv   = dir + "/crs_xy.csv";
    const std::string both_csv = dir + "/crs_both.csv";

    {
        std::ofstream out(ll_csv);
        out << std::setprecision(15) << "latitude,longitude,tag\n"
            << a.lat << "," << a.lon << ",one\n";
    }
    {
        std::ofstream out(en_csv);
        out << std::setprecision(15) << "easting,northing,tag\n"
            << k_ax << "," << k_ay << ",one\n";
    }
    {
        std::ofstream out(xy_csv);
        out << std::setprecision(15) << "x,y,tag\n"
            << k_ax << "," << k_ay << ",one\n";
    }
    {
        // Both present: latitude/longitude must win. The easting/northing here
        // deliberately point somewhere else, so a wrong choice is visible.
        std::ofstream out(both_csv);
        out << std::setprecision(15) << "latitude,longitude,easting,northing,tag\n"
            << a.lat << "," << a.lon << "," << k_bx << "," << k_by << ",one\n";
    }

    const auto load = [&](const std::string& path) {
        return lsdtt::load_point_layer_with_attributes(
            path, meta, std::vector<std::string>{});
    };

    const auto ll   = load(ll_csv);
    const auto en   = load(en_csv);
    const auto xy   = load(xy_csv);
    const auto both = load(both_csv);

    const bool have = ll.records.size() == 1 && en.records.size() == 1 &&
                      xy.records.size() == 1 && both.records.size() == 1;

    check("raster-CRS CSV: every variant reads one record", have,
          "sizes " + std::to_string(ll.records.size()) + "/" +
              std::to_string(en.records.size()) + "/" +
              std::to_string(xy.records.size()) + "/" +
              std::to_string(both.records.size()));
    if (!have) return;

    const auto& p_ll   = ll.records[0].point;
    const auto& p_en   = en.records[0].point;
    const auto& p_xy   = xy.records[0].point;
    const auto& p_both = both.records[0].point;

    check("raster-CRS CSV: easting/northing lands on the same pixel as lat/long",
          p_en.row == p_ll.row && p_en.col == p_ll.col,
          "easting/northing (" + std::to_string(p_en.row) + "," +
              std::to_string(p_en.col) + ") vs lat/long (" +
              std::to_string(p_ll.row) + "," + std::to_string(p_ll.col) + ")");

    check("raster-CRS CSV: x/y are accepted as aliases",
          p_xy.row == p_ll.row && p_xy.col == p_ll.col,
          "x/y (" + std::to_string(p_xy.row) + "," + std::to_string(p_xy.col) + ")");

    check("raster-CRS CSV: latitude/longitude wins when a file carries both",
          p_both.row == p_ll.row && p_both.col == p_ll.col,
          "both-columns file landed at (" + std::to_string(p_both.row) + "," +
              std::to_string(p_both.col) + "), lat/long is (" +
              std::to_string(p_ll.row) + "," + std::to_string(p_ll.col) + ")");
}

// =============================================================================
// C. Typed accessors
// =============================================================================
lsdtt::lsdtt_point_attribute_value make_text_value(const std::string& text)
{
    lsdtt::lsdtt_point_attribute_value value;
    value.type = lsdtt::lsdtt_point_attribute_value::kind::text;
    value.s    = text;
    return value;
}

lsdtt::lsdtt_point_attribute_value make_integer_value(std::int64_t i)
{
    lsdtt::lsdtt_point_attribute_value value;
    value.type = lsdtt::lsdtt_point_attribute_value::kind::integer;
    value.i    = i;
    return value;
}

lsdtt::lsdtt_point_attribute_value make_real_value(double d)
{
    lsdtt::lsdtt_point_attribute_value value;
    value.type = lsdtt::lsdtt_point_attribute_value::kind::real;
    value.d    = d;
    return value;
}

// Every spelling a reviewer might type into an "accepted" column.
void test_as_bool_spellings()
{
    struct bool_case
    {
        const char* text;
        bool        expected;
    };

    const bool_case cases[] = {
        {"1", true},     {"0", false},
        {"true", true},  {"false", false},
        {"TRUE", true},  {"False", false},
        {"yes", true},   {"no", false},
        {"YES", true},   {"No", false},
        {"t", true},     {"f", false},
        {"T", true},     {"F", false},
        {"y", true},     {"n", false},
        {"Y", true},     {"N", false},
        {"  TrUe  ", true}, {"  n  ", false}
    };

    bool all_ok = true;
    std::string failures;

    for (const bool_case& c : cases)
    {
        bool value = !c.expected;
        const bool ok = make_text_value(c.text).as_bool(value);

        if (!ok || value != c.expected)
        {
            all_ok = false;
            failures += std::string(" '") + c.text + "'";
        }
    }

    check("as_bool: every accepted spelling, case insensitively", all_ok,
          all_ok ? "" : "failed on:" + failures);

    // Any non-zero number is true, from any representation.
    bool from_integer   = false;
    bool from_zero      = true;
    bool from_real      = false;
    bool from_real_zero = true;
    bool from_numeric_text = false;

    const bool numeric_ok =
        make_integer_value(5).as_bool(from_integer) &&
        make_integer_value(0).as_bool(from_zero) &&
        make_real_value(-2.5).as_bool(from_real) &&
        make_real_value(0.0).as_bool(from_real_zero) &&
        make_text_value("3.5").as_bool(from_numeric_text);

    check("as_bool: any non-zero number is true, zero is false",
          numeric_ok && from_integer && !from_zero && from_real &&
          !from_real_zero && from_numeric_text);
}

// This data is human-edited, so the accessors must report failure through the
// return value.  A throw here would take down a whole reviewer round trip.
void test_accessors_report_failure_without_throwing()
{
    bool threw = false;

    bool         b = false;
    std::int64_t i = 0;
    double       d = 0.0;
    std::string  s = "";

    bool junk_rejected     = false;
    bool partial_rejected  = false;
    bool absent_rejected   = false;
    bool fractional_reject = false;
    bool conversions_ok    = false;

    try
    {
        const lsdtt::lsdtt_point_attribute_value junk = make_text_value("maybe");
        junk_rejected = !junk.as_bool(b) && !junk.as_int64(i) &&
                        !junk.as_double(d);

        const lsdtt::lsdtt_point_attribute_value partial =
            make_text_value("12abc");
        partial_rejected = !partial.as_int64(i) && !partial.as_double(d);

        const lsdtt::lsdtt_point_attribute_value absent;
        absent_rejected = !absent.is_present() && !absent.as_bool(b) &&
                          !absent.as_int64(i) && !absent.as_double(d) &&
                          !absent.as_string(s);

        // A real is only an int64 when it is exactly integral.
        fractional_reject = !make_real_value(2.5).as_int64(i);

        std::int64_t whole_real = 0;
        std::int64_t from_text  = 0;
        double       widened    = 0.0;
        std::string  rendered   = "";

        conversions_ok =
            make_real_value(4.0).as_int64(whole_real) && whole_real == 4 &&
            make_text_value(" 42 ").as_int64(from_text) && from_text == 42 &&
            make_integer_value(42).as_double(widened) && widened == 42.0 &&
            make_integer_value(42).as_string(rendered) && rendered == "42" &&
            !make_real_value(std::numeric_limits<double>::quiet_NaN())
                 .as_double(widened);
    }
    catch (const std::exception& e)
    {
        threw = true;
        check("accessors: junk input must not throw", false, e.what());
    }

    if (threw)
    {
        return;
    }

    check("accessors: junk text is rejected by return value, not by throwing",
          junk_rejected);
    check("accessors: a partially numeric token is rejected", partial_rejected);
    check("accessors: an absent value fails every accessor", absent_rejected);
    check("accessors: a non-integral real is not an int64", fractional_reject);
    check("accessors: the sane conversions still work, non-finite reals do not",
          conversions_ok);
}

} // namespace


int main(int argc, char** argv)
{
    const std::string output_dir =
        (argc > 1) ? std::string(argv[1]) : std::string("point_input_test_output");

    std::cout << "lsdtt3 point-input unit tests "
                 "(internal-drainage plan sections 6.3 and 6.4)\n\n";

    try
    {
        GDALAllRegister();
        std::filesystem::create_directories(output_dir);

        const lsdtt::lsdtt_raster_meta meta = make_test_raster_meta();

        const std::string injection_path = output_dir + "/injection_points.gpkg";
        const std::string text_area_path = output_dir + "/injection_text.gpkg";
        const std::string fgb_path       = output_dir + "/sink_points.fgb";
        const std::string gpkg_path      = output_dir + "/sink_points.gpkg";
        const std::string csv_path       = output_dir + "/sink_points.csv";

        write_injection_layer(meta, injection_path);
        write_text_area_layer(meta, text_area_path);
        write_attribute_layer(meta, "FlatGeobuf", fgb_path);
        write_attribute_layer(meta, "GPKG", gpkg_path);
        write_attribute_csv(meta, csv_path);

        std::cout << "--- load_drainage_injection_points ---\n";
        test_injection_missing_field_message(meta, injection_path);
        test_injection_open_failure_message(
            meta, output_dir + "/not_a_dataset.gpkg");
        test_injection_skips_and_point_ids(meta, injection_path);
        test_injection_non_numeric_string_area_is_skipped(meta, text_area_path);

        std::cout << "\n--- load_point_layer_with_attributes ---\n";
        test_cross_format_equivalence(meta, fgb_path, gpkg_path, csv_path);
        test_missing_field_is_absent_not_a_throw(meta, fgb_path, csv_path);
        test_order_status_and_dataset_throws(meta, injection_path, output_dir);
        test_csv_in_raster_crs(meta, output_dir);

        std::cout << "\n--- typed attribute accessors ---\n";
        test_as_bool_spellings();
        test_accessors_report_failure_without_throwing();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_checks;
        ++g_n_failed;
    }

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";

    return (g_n_failed == 0) ? 0 : 1;
}
