// test_dem_difference_alignment.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// The grid handling behind DEMs of difference in lsdtt-raster-calculator:
//
//   * non-square geographic pixels (Copernicus 30 m above 50 degrees) keep their
//     geometry through the loader, the geometry helpers and a warp;
//   * interpolating an Int16 DEM can keep its fractional values;
//   * a warp can refuse cells whose kernel touched nodata or left the source;
//   * vertical datum parsing and refusals;
//   * Nuth and Kaab: the initial vertical guess does not bias the result, and a
//     stable-terrain mask excludes real change.
//
// Usage:
//   lsdtt-test-dem-difference-alignment
//   lsdtt-test-dem-difference-alignment --write-cli-fixture <dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal_priv.h>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/core/lsdtt_raster_warp.hpp"
#include "../../src/core/lsdtt_vertical_datum.hpp"
#include "../../src/cpu/lsdtt_raster_coregistration_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/utils/lsdtt_raster_geometry.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;
int g_checks = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_checks;
    if (!ok) ++g_failures;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
}

const double k_arcsec = 1.0 / 3600.0;

lsdtt_raster_meta geographic_meta(int nr, int nc, double x0, double y0, double dx, double dy)
{
    lsdtt_raster_meta m;
    m.n_rows = nr;
    m.n_cols = nc;
    m.x_min = x0;
    m.y_min = y0;
    m.cell_size = dx;
    m.cell_size_y = std::abs(dx - dy) > 1e-12 ? dy : 0.0;
    m.no_data_value = -9999.0;
    m.set_crs_from_epsg(4326);
    return m;
}

// A plane in metres of height per degree: any error in a pixel's position shows
// up as a height error.
std::vector<float> plane(const lsdtt_raster_meta& m)
{
    std::vector<float> z(static_cast<std::size_t>(m.n_rows) * m.n_cols);
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c)
        {
            const lsdtt_xy p = cell_center_from_row_col(m, r, c);
            z[static_cast<std::size_t>(r) * m.n_cols + c] =
                static_cast<float>(1000.0 + 4000.0 * (p.x - 104.0) + 9000.0 * (p.y - 31.0));
        }
    return z;
}

void write_int16(const std::string& path, const lsdtt_raster_meta& m, const std::vector<std::int16_t>& v)
{
    gdal_init();
    GDALDataset* d = GetGDALDriverManager()->GetDriverByName("GTiff")
                         ->Create(path.c_str(), m.n_cols, m.n_rows, 1, GDT_Int16, nullptr);
    if (d == nullptr) throw std::runtime_error("cannot create " + path);
    double gt[6] = {m.x_min, m.cell_size, 0.0, m.y_max(), 0.0, -m.pixel_height()};
    d->SetGeoTransform(gt);
    d->SetProjection(m.crs_wkt.c_str());
    d->GetRasterBand(1)->SetNoDataValue(-32768);
    if (d->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, m.n_cols, m.n_rows,
                                      const_cast<std::int16_t*>(v.data()), m.n_cols, m.n_rows,
                                      GDT_Int16, 0, 0) != CE_None)
        throw std::runtime_error("cannot write " + path);
    GDALClose(d);
}

void test_non_square(const std::string& dir)
{
    std::cout << "\n--- non-square geographic pixels ---\n";
    const lsdtt_raster_meta m = geographic_meta(40, 30, 104.0, 31.0, 1.5 * k_arcsec, k_arcsec);
    const std::string path = dir + "/non_square.tif";
    lsdtt_data_loader writer;
    writer.write_full(plane(m), m, path, raster_output_format::gtiff);

    lsdtt_data_loader loader;
    loader.open_raster(path);
    const lsdtt_raster_meta r = loader.get_meta();
    check("loader keeps the pixel height", std::abs(r.pixel_height() - k_arcsec) < 1e-12);
    check("top edge is right", std::abs(r.y_max() - (31.0 + 40.0 * k_arcsec)) < 1e-9);
    check("first row centre is right",
          std::abs(cell_centre_y(r, 0) - (31.0 + 39.5 * k_arcsec)) < 1e-9);

    lsdtt_raster_meta square = r;
    square.cell_size_y = 0.0;
    check("a different pixel height is a different grid", !same_grid(r, square));

    lsdtt_warp_options options;
    options.resampling_method = "bilinear";
    options.force_float32 = true;
    options.invalidate_incomplete_kernels = true;
    const lsdtt_raster_meta target = geographic_meta(40, 30, 104.0 + 2.0 * k_arcsec, 31.0 + 2.5 * k_arcsec, k_arcsec, k_arcsec);
    const lsdtt_raster_tile tile = warp_raster_to_grid_in_memory(path, target, options);
    const std::vector<float> truth = plane(target);
    double worst = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < truth.size(); ++i)
        if (tile.mask[i]) { worst = std::max(worst, std::abs(static_cast<double>(tile.data[i] - truth[i]))); ++n; }
    check("a warp from non-square pixels lands on the true plane", n > 300 && worst < 0.01,
          std::to_string(n) + " cells, worst " + std::to_string(worst));
}

void test_int16(const std::string& dir)
{
    std::cout << "\n--- Int16 sources ---\n";
    const lsdtt_raster_meta m = geographic_meta(20, 20, 104.0, 31.0, k_arcsec, k_arcsec);
    std::vector<std::int16_t> stripes(400);
    for (int r = 0; r < 20; ++r)
        for (int c = 0; c < 20; ++c) stripes[static_cast<std::size_t>(r) * 20 + c] = static_cast<std::int16_t>(c % 2);
    const std::string path = dir + "/int16.tif";
    write_int16(path, m, stripes);

    lsdtt_raster_meta half = m;
    half.x_min += 0.5 * k_arcsec;
    half.n_cols = 19;
    lsdtt_warp_options options;
    options.resampling_method = "bilinear";
    const lsdtt_raster_tile rounded = warp_raster_to_grid_in_memory(path, half, options);
    options.force_float32 = true;
    const lsdtt_raster_tile kept = warp_raster_to_grid_in_memory(path, half, options);
    const std::size_t i = 10 * 19 + 9;
    check("without force_float32 the value is rounded (the old behaviour)",
          rounded.data[i] == std::round(rounded.data[i]));
    check("with force_float32 a half-pixel value is 0.5", std::abs(kept.data[i] - 0.5f) < 1e-5,
          std::to_string(kept.data[i]));
}

void test_kernel_validity(const std::string& dir)
{
    std::cout << "\n--- kernels touching missing data ---\n";
    const lsdtt_raster_meta m = geographic_meta(30, 30, 104.0, 31.0, k_arcsec, k_arcsec);
    std::vector<float> z = plane(m);
    z[15 * 30 + 15] = -9999.0f;
    const std::string path = dir + "/hole.tif";
    lsdtt_data_loader writer;
    writer.write_full(z, m, path, raster_output_format::gtiff);

    lsdtt_raster_meta half = m;
    half.x_min += 0.5 * k_arcsec;
    half.y_min += 0.5 * k_arcsec;
    half.n_rows = 29;
    half.n_cols = 29;
    lsdtt_warp_options options;
    options.resampling_method = "bilinear";
    options.force_float32 = true;
    const lsdtt_raster_tile permissive = warp_raster_to_grid_in_memory(path, half, options);
    options.invalidate_incomplete_kernels = true;
    const lsdtt_raster_tile strict = warp_raster_to_grid_in_memory(path, half, options);

    // Row 15 and col 15 of the source: the four target cells around the hole.
    int n_lost = 0, n_gdal_kept = 0;
    for (int r = 14; r <= 15; ++r)
        for (int c = 14; c <= 15; ++c)
        {
            if (!strict.mask[static_cast<std::size_t>(r) * 29 + c]) ++n_lost;
            if (permissive.mask[static_cast<std::size_t>(r) * 29 + c]) ++n_gdal_kept;
        }
    check("GDAL alone keeps cells whose kernel touched the hole", n_gdal_kept >= 3,
          std::to_string(n_gdal_kept));
    check("the strict warp drops all four cells beside a hole", n_lost == 4, std::to_string(n_lost));
    check("the next cell is kept", strict.mask[12 * 29 + 12] != 0u);

    lsdtt_raster_meta wide = m;
    wide.x_min -= 5.5 * k_arcsec;
    const lsdtt_raster_tile edge = warp_raster_to_grid_in_memory(path, wide, options);
    bool edge_ok = true;
    for (int r = 2; r < 28; ++r)
        edge_ok = edge_ok && edge.mask[static_cast<std::size_t>(r) * 30 + 5] == 0u &&
                  edge.mask[static_cast<std::size_t>(r) * 30 + 6] != 0u;
    check("a kernel running off the source is dropped", edge_ok);
}

void test_datum()
{
    std::cout << "\n--- vertical datum ---\n";
    check("parses egm96", parse_vertical_datum("egm96") == lsdtt_vertical_datum::egm96);
    check("parses egm2008", parse_vertical_datum("egm2008") == lsdtt_vertical_datum::egm2008);
    bool bad = false, unspecified = false;
    try { (void)parse_vertical_datum("wgs84"); } catch (const std::invalid_argument&) { bad = true; }
    const lsdtt_raster_meta m = geographic_meta(4, 4, 104.0, 31.0, k_arcsec, k_arcsec);
    try { (void)vertical_datum_shift(m, lsdtt_vertical_datum::unspecified, lsdtt_vertical_datum::egm96); }
    catch (const std::invalid_argument&) { unspecified = true; }
    check("an unknown datum is refused", bad);
    check("an unspecified datum is refused", unspecified);
    const std::vector<float> none = vertical_datum_shift(m, lsdtt_vertical_datum::egm96, lsdtt_vertical_datum::egm96);
    check("the same datum needs no grid and shifts by 0", none.size() == 16 && none[5] == 0.0f);
}

// Rolling terrain on a projected grid, sampled with an offset.
lsdtt_raster_cpu rolling(int n, double dx, double dy, double dz)
{
    lsdtt_raster_meta m;
    m.n_rows = n;
    m.n_cols = n;
    m.x_min = 500000.0;
    m.y_min = 4000000.0;
    m.cell_size = 10.0;
    m.no_data_value = -9999.0;
    m.set_crs_from_epsg(32633);
    std::vector<float> z(static_cast<std::size_t>(n) * n);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
        {
            const lsdtt_xy p = cell_center_from_row_col(m, r, c);
            const double x = p.x - 500000.0 + dx, y = p.y - 4000000.0 + dy;
            z[static_cast<std::size_t>(r) * n + c] = static_cast<float>(
                500.0 + 60.0 * std::sin(x / 170.0) * std::cos(y / 130.0) + 35.0 * std::sin((x + y) / 90.0) + dz);
        }
    return lsdtt_raster_cpu(std::move(z), std::vector<std::uint8_t>(static_cast<std::size_t>(n) * n, 1u), m);
}

void test_nuth_kaab()
{
    std::cout << "\n--- Nuth and Kaab ---\n";
    const lsdtt_raster_cpu reference = rolling(160, 0.0, 0.0, 0.0);
    const lsdtt_raster_cpu moving = rolling(160, -4.0, 3.0, -2.5);

    dem_coreg_options options;
    options.method = dem_coreg_method::nuth_kaab;
    options.nuth_kaab.nuth_kaab_coreg_max_iterations = 12;
    const dem_coreg_result plain = dem_coregister_to_reference(reference, moving, options);
    check("recovers the horizontal shift", std::abs(plain.dx_m - 4.0) < 0.1 && std::abs(plain.dy_m + 3.0) < 0.1,
          std::to_string(plain.dx_m) + ", " + std::to_string(plain.dy_m));
    check("recovers the vertical shift", std::abs(plain.dz_m - 2.5) < 0.05, std::to_string(plain.dz_m));

    options.nuth_kaab.nuth_kaab_coreg_initial_dz_m = 40.0;
    const dem_coreg_result guessed = dem_coregister_to_reference(reference, moving, options);
    check("an initial vertical guess does not bias the result",
          std::abs(guessed.dz_m - plain.dz_m) < 1e-6, std::to_string(guessed.dz_m));
    options.nuth_kaab.nuth_kaab_coreg_initial_dz_m = 0.0;

    // Real change over 60 % of the grid defeats the median without a mask.
    std::vector<float> changed = moving.get_data();
    std::vector<std::uint8_t> stable(changed.size(), 1u);
    for (int r = 0; r < 160; ++r)
        for (int c = 0; c < 96; ++c)
        {
            changed[static_cast<std::size_t>(r) * 160 + c] -= 25.0f;
            stable[static_cast<std::size_t>(r) * 160 + c] = 0u;
        }
    const lsdtt_raster_cpu moved(std::move(changed), std::vector<std::uint8_t>(moving.get_mask()), moving.get_meta());
    const dem_coreg_result unmasked = dem_coregister_to_reference(reference, moved, options);
    options.nuth_kaab.stable_terrain_mask = &stable;
    const dem_coreg_result masked = dem_coregister_to_reference(reference, moved, options);
    check("without a mask, change over most of the grid biases the vertical shift",
          std::abs(unmasked.dz_m - 2.5) > 5.0, std::to_string(unmasked.dz_m));
    check("with a stable mask the vertical shift is right", std::abs(masked.dz_m - 2.5) < 0.1,
          std::to_string(masked.dz_m));
    check("with a stable mask the horizontal shift is right",
          std::abs(masked.dx_m - 4.0) < 0.15 && std::abs(masked.dy_m + 3.0) < 0.15,
          std::to_string(masked.dx_m) + ", " + std::to_string(masked.dy_m));

    std::vector<std::uint8_t> wrong(3, 1u);
    options.nuth_kaab.stable_terrain_mask = &wrong;
    bool refused = false;
    try { (void)dem_coregister_to_reference(reference, moved, options); }
    catch (const std::invalid_argument&) { refused = true; }
    check("a stable mask of the wrong size is refused", refused);
}

// Planes in degrees on the grids of the real products: 1 arc-second cells, on
// and off whole arc-seconds, and 1.5 x 1 arc-second above 50 degrees.
void write_cli_fixture(const std::string& dir)
{
    std::filesystem::create_directories(dir);
    lsdtt_data_loader w;
    const lsdtt_raster_meta first = geographic_meta(120, 120, 104.0, 31.0, k_arcsec, k_arcsec);
    w.write_full(plane(first), first, dir + "/plane_first.tif", raster_output_format::gtiff);
    const lsdtt_raster_meta half = geographic_meta(120, 120, 104.0 + 0.5 * k_arcsec, 31.0 + 0.5 * k_arcsec, k_arcsec, k_arcsec);
    std::vector<float> second = plane(half);
    w.write_full(second, half, dir + "/plane_half_pixel.tif", raster_output_format::gtiff);
    for (float& v : second) v = static_cast<float>(v / 0.3048);
    w.write_full(second, half, dir + "/plane_half_pixel_feet.tif", raster_output_format::gtiff);
    const lsdtt_raster_meta wide = geographic_meta(100, 90, 104.0, 31.0, 1.5 * k_arcsec, k_arcsec);
    w.write_full(plane(wide), wide, dir + "/plane_non_square.tif", raster_output_format::gtiff);
    const lsdtt_raster_meta away = geographic_meta(20, 20, 110.0, 40.0, k_arcsec, k_arcsec);
    w.write_full(plane(away), away, dir + "/plane_elsewhere.tif", raster_output_format::gtiff);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc == 3 && std::string(argv[1]) == "--write-cli-fixture")
        {
            write_cli_fixture(argv[2]);
            return 0;
        }
        const std::filesystem::path dir = std::filesystem::temp_directory_path() /
            ("lsdtt_test_dem_difference_" + std::to_string(std::random_device{}()));
        std::filesystem::create_directories(dir);
        test_non_square(dir.string());
        test_int16(dir.string());
        test_kernel_validity(dir.string());
        test_datum();
        test_nuth_kaab();
        std::error_code ignored;
        std::filesystem::remove_all(dir, ignored);
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL  unexpected exception: " << e.what() << "\n";
        ++g_failures;
    }
    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
