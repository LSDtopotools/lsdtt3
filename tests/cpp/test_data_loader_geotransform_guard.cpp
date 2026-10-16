// test_data_loader_geotransform_guard.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// lsdtt_data_loader refuses geotransforms its meta cannot represent.
//
// lsdtt_raster_meta keeps only the origin and one cell size, so a rotated,
// non-square or south-up raster would be misread by every CLI without a word.
// The loader refuses them at open_raster:
//
//   G1  north-up, square, rectangular (3x5), nodata border: opens, and the
//       meta equals the geotransform.
//   G2  gt[2] = 0.5 -> "rotated".       G3  gt[4] = 0.5 -> "rotated".
//   G4  dx 10, dy 12 -> "not square".   G5  dx 10, dy 10.000000000001 -> opens.
//   G6  gt[5] = +10 -> "not north-up".
//   G7  UTM (projected) with dx 10, dy 12 -> "not square".
//   G8  EPSG:4326 with 1.5" x 1" pixels (Copernicus above 50 degrees) -> opens,
//       so lsdtt-raster-preprocessing can reproject it; the analysis CLIs
//       refuse geographic rasters themselves.
//
// Usage: lsdtt-test-data-loader-geotransform-guard <scratch_dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_gdal_init.hpp"

namespace
{

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

// A 3-row x 5-col float GeoTIFF with the given geotransform and a nodata
// border (only the middle row's three inner pixels are data).
std::string write_tif(const std::string& dir, const std::string& name, const double gt[6],
                      int epsg = 0)
{
    const std::string path = dir + "/" + name;
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset* ds = drv->Create(path.c_str(), 5, 3, 1, GDT_Float32, nullptr);
    double g[6] = {gt[0], gt[1], gt[2], gt[3], gt[4], gt[5]};
    ds->SetGeoTransform(g);
    if (epsg > 0)
    {
        OGRSpatialReference srs;
        srs.importFromEPSG(epsg);
        ds->SetSpatialRef(&srs);
    }
    GDALRasterBand* b = ds->GetRasterBand(1);
    b->SetNoDataValue(-9999.0);
    std::vector<float> v(15, -9999.0f);
    v[6] = 1.0f;
    v[7] = 2.0f;
    v[8] = 3.0f;
    if (b->RasterIO(GF_Write, 0, 0, 5, 3, v.data(), 5, 3, GDT_Float32, 0, 0) != CE_None)
        std::cerr << "write failed: " << path << "\n";
    GDALClose(ds);
    return path;
}

// "" if it opened; the exception text otherwise.
std::string open_error(const std::string& path, lsdtt::lsdtt_raster_meta* meta = nullptr)
{
    try
    {
        lsdtt::lsdtt_data_loader loader;
        loader.open_raster(path);
        if (meta) *meta = loader.get_meta();
        return "";
    }
    catch (const std::exception& e)
    {
        return e.what();
    }
}

bool contains(const std::string& s, const std::string& needle)
{
    return s.find(needle) != std::string::npos;
}

} // namespace

int main(int argc, char* argv[])
{
    const std::string dir = (argc > 1) ? argv[1] : ".";
    std::filesystem::create_directories(dir);
    lsdtt::gdal_init();
    std::cout << "== lsdtt_data_loader geotransform guard ==\n";

    {
        const double gt[6] = {500000.0, 10.0, 0.0, 8700030.0, 0.0, -10.0};
        lsdtt::lsdtt_raster_meta m;
        const std::string err = open_error(write_tif(dir, "gt_ok.tif", gt), &m);
        check("G1 north-up, square, 3x5, nodata border: opens", err.empty(), err);
        check("G1 meta equals the geotransform",
              m.n_rows == 3 && m.n_cols == 5 && m.cell_size == 10.0 && m.x_min == 500000.0 &&
                  m.y_min == 8700000.0);
    }
    {
        const double gt[6] = {500000.0, 10.0, 0.5, 8700030.0, 0.0, -10.0};
        const std::string err = open_error(write_tif(dir, "gt_rot2.tif", gt));
        check("G2 gt[2] = 0.5: refused as rotated", contains(err, "rotated"), err);
    }
    {
        const double gt[6] = {500000.0, 10.0, 0.0, 8700030.0, 0.5, -10.0};
        const std::string err = open_error(write_tif(dir, "gt_rot4.tif", gt));
        check("G3 gt[4] = 0.5: refused as rotated", contains(err, "rotated"), err);
    }
    {
        const double gt[6] = {500000.0, 10.0, 0.0, 8700036.0, 0.0, -12.0};
        const std::string err = open_error(write_tif(dir, "gt_nonsquare.tif", gt));
        check("G4 dx 10, dy 12: refused as not square", contains(err, "not square"), err);
    }
    {
        const double gt[6] = {500000.0, 10.0, 0.0, 8700030.0, 0.0, -10.000000000001};
        const std::string err = open_error(write_tif(dir, "gt_rounding.tif", gt));
        check("G5 dy 10.000000000001 (print rounding): opens", err.empty(), err);
    }
    {
        const double gt[6] = {500000.0, 10.0, 0.0, 8700000.0, 0.0, 10.0};
        const std::string err = open_error(write_tif(dir, "gt_southup.tif", gt));
        check("G6 gt[5] = +10: refused as not north-up", contains(err, "not north-up"), err);
    }

    {
        const double gt[6] = {500000.0, 10.0, 0.0, 8700036.0, 0.0, -12.0};
        const std::string err = open_error(write_tif(dir, "gt_utm_nonsquare.tif", gt, 32758));
        check("G7 projected (UTM 58S), dx 10, dy 12: refused as not square",
              contains(err, "not square"), err);
    }
    {
        const double gt[6] = {-5.5, 1.5 / 3600.0, 0.0, 56.5, 0.0, -1.0 / 3600.0};
        lsdtt::lsdtt_raster_meta m;
        const std::string err = open_error(write_tif(dir, "gt_geographic.tif", gt, 4326), &m);
        check("G8 EPSG:4326 with 1.5\" x 1\" pixels: opens (for reprojection)",
              err.empty(), err);
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
