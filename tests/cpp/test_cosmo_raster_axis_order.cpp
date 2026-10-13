// test_cosmo_raster_axis_order.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// The production raster must find the same latitude for a pixel whatever the
// axis order the DEM's CRS declares.
//
// calculate_production_raster and calculate_pressure_raster convert pixel
// centres to WGS84 latitude/longitude. They used to put only the EPSG:4326
// destination in traditional GIS (x = easting) axis order. A projected CRS
// whose WKT declares NORTHING first (EPSG:2193 NZTM, EPSG:3006 SWEREF99 TM)
// then had its pixel (x, y) read as (northing, easting): an NZTM point at
// (1750000, 5900000), which is near Wellington at (174.69 E, 37.03 S), came out
// at (117.45 W, 51.20 S), so the Stone scaling was evaluated at the wrong
// latitude.
//
// Each case puts a 3x3 DEM of constant elevation, centred on the same point,
// once in the northing-first CRS and once in the UTM zone around it (easting
// first), and requires the centre pixel's production to agree. Stone scaling
// with the standard atmosphere depends on the location only through latitude,
// so a swapped-axis latitude shows up as a difference of several percent.
//
// Usage: lsdtt-test-cosmo-raster-axis-order
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/core/lsdtt_vector_utils.hpp"
#include "../../src/cpu/lsdtt_cosmo_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

// 3x3 DEM at elevation_m, in CRS epsg, whose centre pixel centre is (x, y).
lsdtt_raster_cpu make_dem(int epsg, double x, double y, double cell_size, float elevation_m)
{
    lsdtt_raster_meta meta;
    meta.n_rows = 3;
    meta.n_cols = 3;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.x_min = x - 1.5 * cell_size;
    meta.y_min = y - 1.5 * cell_size;
    meta.set_crs_from_epsg(epsg);
    return lsdtt_raster_cpu(std::vector<float>(9, elevation_m), std::vector<std::uint8_t>(9, 1), meta);
}

float centre_production(const lsdtt_raster_cpu& dem)
{
    return calculate_production_raster(dem).get_data()[4];
}

void run_case(const std::string& label, int northing_first_epsg, double x, double y, int utm_epsg)
{
    const lsdtt_raster_cpu dem_nf = make_dem(northing_first_epsg, x, y, 30.0, 500.0f);

    // The same point in the UTM zone. lsdtt_coordinate_transform puts both CRSs
    // in traditional GIS order, so this conversion is independent of the code
    // under test.
    lsdtt_raster_meta utm_meta;
    utm_meta.set_crs_from_epsg(utm_epsg);
    const lsdtt_coordinate_transform to_utm(dem_nf.get_meta().crs_wkt, utm_meta.crs_wkt);
    const lsdtt_xy p = to_utm.transform_xy(x, y);
    const lsdtt_raster_cpu dem_utm = make_dem(utm_epsg, p.x, p.y, 30.0, 500.0f);

    const lsdtt_lonlat ll = cell_center_to_lonlat(dem_nf.get_meta(), 1, 1);
    const float prod_nf = centre_production(dem_nf);
    const float prod_utm = centre_production(dem_utm);
    const double rel = std::abs(static_cast<double>(prod_nf) - prod_utm) / std::abs(static_cast<double>(prod_utm));
    check(label + ": centre production in EPSG:" + std::to_string(northing_first_epsg) +
              " matches EPSG:" + std::to_string(utm_epsg) + " (point at lon " + std::to_string(ll.longitude) +
              ", lat " + std::to_string(ll.latitude) + ")",
          rel <= 1e-6,
          "EPSG:" + std::to_string(northing_first_epsg) + " " + std::to_string(prod_nf) + ", EPSG:" +
              std::to_string(utm_epsg) + " " + std::to_string(prod_utm) + ", relative difference " +
              std::to_string(rel));
}

} // namespace

int main()
{
    std::cout << "== cosmo raster: WGS84 conversion independent of the CRS axis order ==\n";
    run_case("NZTM (northing first) vs UTM 60S", 2193, 1750000.0, 5900000.0, 32760);
    run_case("SWEREF99 TM (northing first) vs UTM 33N", 3006, 674000.0, 6580000.0, 32633);
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
