// test_crs_semantic_matching.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Rasters in one CRS whose WKT text differs (an lsdtt3-written DEM and a mask
// made by gdalwarp in the same oblique Mercator frame) must be accepted
// together. Checks crs_equal_semantically and the core entry points that use
// it for a user-supplied raster: resample_raster and the hierarchy's
// set_active.
//
//   C1  identical strings, two empty strings: equal.
//   C2  one empty string: not equal.
//   C3  omerc as WKT1 and as WKT2_2019 (different text): equal.
//   C4  omerc with a different azimuth: not equal.
//   C5  UTM 12N from EPSG and from a PROJ string (no authority): equal.
//   C6  resample_raster accepts WKT1/WKT2 omerc with no EPSG code; refuses a
//       different azimuth.
//   C7  hierarchy set_active accepts a WKT2 replacement for a WKT1 level;
//       refuses a different azimuth.
// -----------------------------------------------------------------------------

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cpl_conv.h>
#include <ogr_spatialref.h>

#include "lsdtt_raster_hierarchy_cpu.hpp"
#include "lsdtt_raster_resampling_cpu.hpp"
#include "lsdtt_raster_warp.hpp"

namespace {

int checks = 0;
int failures = 0;

void check(const std::string& name, bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
}

// An oblique Mercator frame along a feature at azimuth alpha (degrees).
std::string omerc_proj(double alpha)
{
    return "+proj=omerc +lat_0=38.4 +lonc=-110.0 +alpha=" + std::to_string(alpha) +
           " +gamma=90 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs";
}

std::string wkt_from_proj(const std::string& proj, const char* format)
{
    OGRSpatialReference srs;
    if (srs.importFromProj4(proj.c_str()) != OGRERR_NONE)
        throw std::runtime_error("could not parse " + proj);
    const char* options[] = {format, nullptr};
    char* raw = nullptr;
    if (srs.exportToWkt(&raw, options) != OGRERR_NONE || raw == nullptr)
        throw std::runtime_error("could not export " + proj);
    std::string wkt(raw);
    CPLFree(raw);
    return wkt;
}

std::string wkt_from_epsg(int epsg)
{
    OGRSpatialReference srs;
    if (srs.importFromEPSG(epsg) != OGRERR_NONE)
        throw std::runtime_error("could not import EPSG:" + std::to_string(epsg));
    char* raw = nullptr;
    srs.exportToWkt(&raw);
    std::string wkt(raw);
    CPLFree(raw);
    return wkt;
}

lsdtt::lsdtt_raster_cpu raster(const std::string& wkt, float value)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 4;
    meta.n_cols = 4;
    meta.cell_size = 30.0;
    meta.x_min = -60.0;
    meta.y_min = -60.0;
    meta.no_data_value = -9999.0;
    meta.crs_wkt = wkt;
    meta.epsg_code = 0;
    return lsdtt::lsdtt_raster_cpu(std::vector<float>(16, value),
                                   std::vector<std::uint8_t>(16, 1), meta);
}

template <typename F>
bool throws(F&& f)
{
    try
    {
        f();
    }
    catch (const std::invalid_argument&)
    {
        return true;
    }
    return false;
}

} // namespace

int main()
{
    const std::string wkt1 = wkt_from_proj(omerc_proj(30.0), "FORMAT=WKT1");
    const std::string wkt2 = wkt_from_proj(omerc_proj(30.0), "FORMAT=WKT2_2019");
    const std::string other = wkt_from_proj(omerc_proj(31.0), "FORMAT=WKT2_2019");

    check("C1 identical strings are equal", lsdtt::crs_equal_semantically(wkt1, wkt1));
    check("C1 two empty strings are equal", lsdtt::crs_equal_semantically("", ""));
    check("C2 one empty string is not equal", !lsdtt::crs_equal_semantically(wkt1, ""));
    check("C3 WKT1 and WKT2 texts differ", wkt1 != wkt2);
    check("C3 omerc WKT1 equals WKT2", lsdtt::crs_equal_semantically(wkt1, wkt2));
    check("C4 different azimuth is not equal", !lsdtt::crs_equal_semantically(wkt1, other));
    check("C5 UTM 12N from EPSG equals PROJ string",
          lsdtt::crs_equal_semantically(
              wkt_from_epsg(32612),
              wkt_from_proj("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs +type=crs",
                            "FORMAT=WKT2_2019")));

    {
        const lsdtt::lsdtt_raster_cpu source = raster(wkt1, 5.0f);
        lsdtt::lsdtt_raster_cpu destination = raster(wkt2, -9999.0f);
        check("C6 resample accepts WKT1 source and WKT2 destination",
              !throws([&] {
                  lsdtt::resample_raster(source, destination,
                                         lsdtt::raster_resampling_method::mean);
              }));
        check("C6 resampled values are the source values",
              destination.get_data() == std::vector<float>(16, 5.0f));
        lsdtt::lsdtt_raster_cpu rotated = raster(other, -9999.0f);
        check("C6 resample refuses a different azimuth", throws([&] {
                  lsdtt::resample_raster(source, rotated,
                                         lsdtt::raster_resampling_method::mean);
              }));
    }

    {
        lsdtt::lsdtt_raster_hierarchy_cpu hierarchy(raster(wkt1, 1.0f));
        check("C7 set_active accepts a WKT2 replacement",
              !throws([&] { hierarchy.set_active(raster(wkt2, 2.0f)); }));
        check("C7 set_active refuses a different azimuth",
              throws([&] { hierarchy.set_active(raster(other, 3.0f)); }));
    }

    std::cout << "RESULT: " << (failures == 0 ? "PASS" : "FAIL") << " (" << checks - failures
              << "/" << checks << ")\n";
    return failures == 0 ? 0 : 1;
}
