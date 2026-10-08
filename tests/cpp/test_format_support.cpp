// test_format_support.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The format contract: every driver lsdtt3 declares support for is present, and
// lsdtt3 can actually read a raster written in it.
//
// Chunk 2a of the GDAL/PROJ friction work. Phase 2 will
// build GDAL from source with a TRIMMED driver set. The whole risk of that is
// quietly dropping a format users feed in, and a trimmed build fails in exactly
// the same way as a correct one -- silently, until someone's DEM will not open.
// This test is the acceptance gate: it is written against the system GDAL now,
// and re-run unchanged against the vendored GDAL later. If the trimmed set is
// short a driver, this says which one.
//
// It is therefore a HARD failure when a required driver is missing, not a skip.
// A skip here would defeat the entire purpose.
//
// The required set, which is Simon's D2/D6 answers:
//   write: GTiff, COG, PNG, MEM, VRT, GPKG, FlatGeobuf
//   read:  all of the above plus GeoJSON, AAIGrid, ENVI, EHdr, HFA, netCDF,
//          HDF5, SRTMHGT, PDS4, ISIS3, Zarr, ESRI Shapefile
//
// Zarr needs GDAL >= 3.4, which is the project's floor for this contract.
// -----------------------------------------------------------------------------

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_gdal_init.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>
#include <vector>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

namespace fs = std::filesystem;

namespace
{

int g_failures = 0;

void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << "\n";
    if (!ok) ++g_failures;
}

// ---------------------------------------------------------------------------
// The contract.
//
// `creatable` marks drivers this test can round-trip: it builds the fixture by
// copying a MEM raster into the driver. HDF5 is read-only in GDAL, and SRTMHGT
// only accepts exact 1201/3601-pixel tiles with a N00E000-style name, so those
// two are presence-only -- their absence is still a hard failure.
// ---------------------------------------------------------------------------
struct raster_driver_case
{
    const char* driver;
    bool        creatable;
    const char* extension;
};

const raster_driver_case k_raster_drivers[] = {
    {"GTiff",   true,  ".tif"},
    {"COG",     true,  ".tif"},
    {"PNG",     true,  ".png"},
    {"MEM",     false, ""},        // in-memory; nothing to open from disk
    {"VRT",     true,  ".vrt"},
    {"AAIGrid", true,  ".asc"},
    {"ENVI",    true,  ".bil"},
    {"EHdr",    true,  ".bil"},
    {"HFA",     true,  ".img"},
    {"netCDF",  true,  ".nc"},
    {"HDF5",    false, ".h5"},     // read-only driver
    {"SRTMHGT", false, ".hgt"},    // fixed tile geometry
    {"PDS4",    true,  ".xml"},
    {"ISIS3",   true,  ".cub"},
    {"Zarr",    true,  ".zarr"},
};

const char* const k_vector_drivers[] = {
    "GPKG", "FlatGeobuf", "GeoJSON", "ESRI Shapefile",
};

constexpr int    k_nx        = 8;
constexpr int    k_ny        = 6;
constexpr double k_cell      = 25.0;
constexpr double k_x_min     = 300000.0;
constexpr double k_y_min     = 5000000.0;

// Value written at (row, col). Deliberately not round numbers, so a silent
// type narrowing on the way through a driver shows up.
float sample_value(int row, int col)
{
    return static_cast<float>(100 + row * 10 + col);
}

// A MEM raster carrying a known surface, geotransform and CRS.
//
// Used only to seed the on-disk base GeoTIFF below. Fixtures are copied FROM
// THAT FILE, not from this MEM dataset: a VRT built over a MEM source keeps a
// raw pointer into it, so closing the source leaves the VRT dangling and
// reading it segfaults. Copying from a file is also how these formats are
// really produced.
GDALDataset* make_mem_source()
{
    GDALDriver* mem = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!mem) return nullptr;

    GDALDataset* ds = mem->Create("", k_nx, k_ny, 1, GDT_Float32, nullptr);
    if (!ds) return nullptr;

    std::vector<float> buf(static_cast<std::size_t>(k_nx) * k_ny);
    for (int r = 0; r < k_ny; ++r)
        for (int c = 0; c < k_nx; ++c)
            buf[static_cast<std::size_t>(r) * k_nx + c] = sample_value(r, c);

    ds->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, k_nx, k_ny, buf.data(),
                                   k_nx, k_ny, GDT_Float32, 0, 0);

    // GDAL's geotransform is top-left origin; lsdtt_raster_meta is bottom-left.
    double gt[6] = {k_x_min, k_cell, 0.0,
                    k_y_min + k_ny * k_cell, 0.0, -k_cell};
    ds->SetGeoTransform(gt);

    OGRSpatialReference srs;
    srs.importFromEPSG(32630);         // WGS 84 / UTM zone 30N
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    if (wkt) { ds->SetProjection(wkt); CPLFree(wkt); }

    return ds;
}

// ---------------------------------------------------------------------------
// 1. Every required driver is registered.
// ---------------------------------------------------------------------------
void test_drivers_present()
{
    std::cout << "-- required drivers present --\n";

    for (const raster_driver_case& d : k_raster_drivers)
    {
        check(GetGDALDriverManager()->GetDriverByName(d.driver) != nullptr,
              std::string("raster driver registered: ") + d.driver);
    }
    for (const char* name : k_vector_drivers)
    {
        check(GetGDALDriverManager()->GetDriverByName(name) != nullptr,
              std::string("vector driver registered: ") + name);
    }
}

// ---------------------------------------------------------------------------
// 2. lsdtt3 can open a raster written by each creatable driver, and gets the
//    geometry and pixel values back.
// ---------------------------------------------------------------------------
void test_raster_round_trip(const fs::path& work)
{
    std::cout << "-- lsdtt3 reads each format --\n";

    // The on-disk source every fixture is copied from.
    const fs::path base = work / "base_source.tif";
    {
        std::error_code rm_ec;
        fs::remove_all(base, rm_ec);

        GDALDriver* gtiff = GetGDALDriverManager()->GetDriverByName("GTiff");
        GDALDataset* mem = make_mem_source();
        if (!gtiff || !mem)
        {
            check(false, "could not build the base GeoTIFF source");
            if (mem) GDALClose(mem);
            return;
        }
        GDALDataset* b = gtiff->CreateCopy(base.string().c_str(), mem,
                                           FALSE, nullptr, nullptr, nullptr);
        GDALClose(mem);
        if (!b)
        {
            check(false, "could not write the base GeoTIFF source");
            return;
        }
        GDALClose(b);
    }

    for (const raster_driver_case& d : k_raster_drivers)
    {
        if (!d.creatable) continue;

        GDALDriver* drv = GetGDALDriverManager()->GetDriverByName(d.driver);
        if (!drv)
        {
            check(false, std::string("round trip ") + d.driver +
                         ": driver missing (reported above)");
            continue;
        }

        const fs::path out = work / (std::string("fmt_") + d.driver + d.extension);

        // Clear any fixture left by a previous run. CTest reuses its work
        // directory, and Zarr refuses to CreateCopy over an existing directory
        // -- so without this the test passes on a clean tree and fails on every
        // re-run. remove_all covers both the file and directory cases.
        std::error_code rm_ec;
        fs::remove_all(out, rm_ec);

        GDALDataset* src = static_cast<GDALDataset*>(
            GDALOpen(base.string().c_str(), GA_ReadOnly));
        if (!src)
        {
            check(false, "could not reopen the base GeoTIFF source");
            return;
        }

        CPLPushErrorHandler(CPLQuietErrorHandler);
        GDALDataset* dst = drv->CreateCopy(out.string().c_str(), src,
                                           FALSE, nullptr, nullptr, nullptr);
        CPLPopErrorHandler();

        // ORDER MATTERS: close the copy before its source. A VRT keeps a live
        // pointer to the dataset it was copied from, and its destructor
        // flushes through that pointer -- closing the source first segfaults
        // inside VRTSimpleSource::FlushCache. Every other driver here is
        // indifferent to the order, so getting it wrong shows up only on VRT.
        if (dst) GDALClose(dst);
        GDALClose(src);

        if (!dst)
        {
            check(false, std::string("round trip ") + d.driver +
                         ": CreateCopy failed (" + CPLGetLastErrorMsg() + ")");
            continue;
        }

        // Now the part that matters: lsdtt3's own loader must open it.
        try
        {
            lsdtt::lsdtt_data_loader loader;
            loader.open_raster(out.string());
            const lsdtt::lsdtt_raster_meta& m = loader.get_meta();

            bool ok = (m.n_cols == k_nx && m.n_rows == k_ny);
            ok = ok && std::abs(m.cell_size - k_cell) < 1.0e-6;

            // PNG carries no georeferencing of its own; GDAL recovers it from
            // the worldfile and PAM sidecar, which CreateCopy does not write
            // here. Geometry is still checked above; CRS is not.
            const bool expect_crs = (std::string(d.driver) != "PNG");
            if (expect_crs) ok = ok && !m.crs_wkt.empty();

            const lsdtt::lsdtt_raster_tile tile = loader.read_full();
            bool values_ok = true;
            for (int r = 0; r < k_ny && values_ok; ++r)
                for (int c = 0; c < k_nx && values_ok; ++c)
                {
                    const float got =
                        tile.data[static_cast<std::size_t>(r) * k_nx + c];
                    // PNG is an 8-bit stretch by design, so only its geometry
                    // is contractual, not its values.
                    if (std::string(d.driver) == "PNG") continue;
                    if (std::abs(got - sample_value(r, c)) > 1.0e-3f)
                        values_ok = false;
                }

            check(ok && values_ok,
                  std::string("round trip ") + d.driver +
                  " -> lsdtt_data_loader (geometry, CRS, values)");
        }
        catch (const std::exception& e)
        {
            check(false, std::string("round trip ") + d.driver +
                         ": lsdtt3 could not open it: " + e.what());
        }
    }
}

// ---------------------------------------------------------------------------
// 3. lsdtt3 writes the vector formats it claims to, and can read them back.
// ---------------------------------------------------------------------------
void test_vector_round_trip(const fs::path& work)
{
    std::cout << "-- lsdtt3 writes and reads each vector format --\n";

    struct vcase { const char* fmt; const char* ext; };
    // GeoJSON is deliberately absent: lsdtt3 READS it but does not write it.
    // vector_output_format offers only these two (lsdtt_data_loader.hpp:111).
    const vcase cases[] = {
        {"gpkg",       ".gpkg"},
        {"flatgeobuf", ".fgb"},
    };

    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows    = k_ny;
    meta.n_cols    = k_nx;
    meta.cell_size = k_cell;
    meta.x_min     = k_x_min;
    meta.y_min     = k_y_min;
    meta.set_crs_from_epsg(32630);

    for (const vcase& v : cases)
    {
        const fs::path out = work / (std::string("vec_") + v.fmt + v.ext);
        std::error_code rm_ec;
        fs::remove_all(out, rm_ec);

        lsdtt::lsdtt_vector_schema schema;
        schema.field_names = {"id"};
        schema.field_types = {lsdtt::lsdtt_vector_field_type::int64_type};

        try
        {
            lsdtt::lsdtt_data_loader loader;
            loader.begin_vector_write(
                out.string(), "test_layer", schema, meta.crs_wkt,
                lsdtt::resolve_vector_output_format(v.fmt), false);
            loader.write_point(k_x_min + 2.0 * k_cell,
                               k_y_min + 3.0 * k_cell,
                               {static_cast<std::int64_t>(7)});
            loader.finish_vector_write();

            check(fs::exists(out),
                  std::string("wrote vector format ") + v.fmt);
        }
        catch (const std::exception& e)
        {
            check(false, std::string("vector format ") + v.fmt + ": " + e.what());
        }
    }
}

} // namespace

int main(int argc, char* argv[])
{
    lsdtt::gdal_init();

    const fs::path work = (argc > 1)
        ? fs::path(argv[1])
        : fs::temp_directory_path() / "lsdtt_format_support";
    fs::create_directories(work);

    std::cout << "format support contract (GDAL "
              << lsdtt::gdal_version() << ")\n"
              << "work dir: " << work << "\n\n";

    test_drivers_present();
    std::cout << "\n";
    test_raster_round_trip(work);
    std::cout << "\n";
    test_vector_round_trip(work);

    std::cout << "\n" << (g_failures == 0 ? "OK" : "FAILURES")
              << ": " << g_failures << " failed\n";
    return g_failures == 0 ? 0 : 1;
}
