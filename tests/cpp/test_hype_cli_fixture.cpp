// test_hype_cli_fixture.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// fixture_hype_cli: writes the input set for the column_model: hype tests of
// lsdtt-distributed-hydrology. A CTest
// FIXTURES_SETUP; nothing here is committed, every file is generated.
//
// Usage: lsdtt-test-hype-cli-fixture <output_dir>
//
// THE DEM (hype_dem.tif) is 60 x 60 at 10 m, on the default synthetic
// georeference (UTM 58S), so it is projected in metres as the CLI requires.
// Elevation rises 0.08 m/m to the north and 0.25 m/m away from a Y-shaped
// network of three straight valley lines:
//
//     trunk      (305, 0)   -> (305, 300)    column 30, draining off the
//                                            southern edge at row 59
//     west arm   (105, 590) -> (305, 300)
//     east arm   (505, 590) -> (305, 300)
//
// (local metres from the south-west corner; 305 m is the centre of column 30,
// so the trunk lies on cell centres and has no elevation ties across it). The
// two arms are first-order channels that meet in a second-order trunk, so the
// Strahler outlet modes have something to select, and a gauge on an arm plus
// one on the trunk make a two-sub-basin chain for the user-outlet mode.
//
// THE CLASSES. Land use (hype_landuse.tif) is 1 west of column 30 and 2 from
// column 30 east; soil (hype_soil.tif) is 1 north of row 30 and 2 from row 30
// south. Four classes cover the four pairs. Parameter values are those of the
// oracle case tests/expected/hype/default_all/lsdtt_input, which are known to
// run, plus a fourth class for the (1, 2) pair.
//
// THE FORCING is 120 days from 2001-01-01: temperature rising from about -6 to
// +12 deg C with a 9-day wobble, so there is snowfall, a snowpack and a melt.
// Rainfall carries a date column, which the CLI must ignore.
//
// FAILURE INPUTS. Three deliberately broken tables for the failure tests:
// fsceff != 0, an unknown land-use parameter column, and a class table with no
// class for the (1, 2) pair.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <ogr_spatialref.h>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_raster_meta.hpp"

namespace
{

constexpr int    k_n_rows = 60;
constexpr int    k_n_cols = 60;
constexpr double k_dx     = 10.0;
constexpr int    k_n_days = 120;

lsdtt::lsdtt_raster_meta fixture_meta()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = k_n_rows;
    meta.n_cols        = k_n_cols;
    meta.cell_size     = k_dx;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    return meta;
}

// Distance from (px, py) to the segment (ax, ay) -> (bx, by).
double segment_distance(double px, double py, double ax, double ay, double bx, double by)
{
    const double vx = bx - ax, vy = by - ay;
    const double t  = std::clamp(((px - ax) * vx + (py - ay) * vy) / (vx * vx + vy * vy), 0.0, 1.0);
    const double dx = px - (ax + t * vx), dy = py - (ay + t * vy);
    return std::sqrt(dx * dx + dy * dy);
}

std::vector<float> make_dem()
{
    std::vector<float> z(static_cast<std::size_t>(k_n_rows) * k_n_cols);
    for (int r = 0; r < k_n_rows; ++r)
    {
        for (int c = 0; c < k_n_cols; ++c)
        {
            const double x = (c + 0.5) * k_dx;
            const double y = (k_n_rows - r - 0.5) * k_dx;
            const double d = std::min({segment_distance(x, y, 305.0, 0.0, 305.0, 300.0),
                                       segment_distance(x, y, 105.0, 590.0, 305.0, 300.0),
                                       segment_distance(x, y, 505.0, 590.0, 305.0, 300.0)});
            z[static_cast<std::size_t>(r) * k_n_cols + c] = static_cast<float>(100.0 + 0.08 * y + 0.25 * d);
        }
    }
    return z;
}

void write_raster(const std::vector<float>& data, const lsdtt::lsdtt_raster_meta& meta, const std::string& path)
{
    lsdtt::lsdtt_data_loader writer;
    writer.write_full(data, meta, path, lsdtt::raster_output_format::gtiff);
    std::cout << "Wrote " << path << "\n";
}

void write_text(const std::string& path, const std::string& text)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("could not write " + path);
    f << text;
    std::cout << "Wrote " << path << "\n";
}

// Two gauges, as row/col cell centres: one on the west arm, one on the trunk
// above the southern edge. The CLI snaps them downstream onto the channel.
void write_outlets(const lsdtt::lsdtt_raster_meta& meta, const std::string& path)
{
    const int rows[2] = {15, 50};
    const int cols[2] = {20, 30};

    OGRSpatialReference utm;
    OGRSpatialReference wgs84;
    if (utm.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE || wgs84.importFromEPSG(4326) != OGRERR_NONE)
        throw std::runtime_error("could not build the fixture CRS transform");
    utm.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformation* ct = OGRCreateCoordinateTransformation(&utm, &wgs84);
    if (ct == nullptr) throw std::runtime_error("could not create the fixture CRS transform");

    std::ofstream csv(path);
    if (!csv) throw std::runtime_error("could not write " + path);
    csv << "latitude,longitude,point_id\n" << std::setprecision(12);
    for (int i = 0; i < 2; ++i)
    {
        double x = meta.x_min + (cols[i] + 0.5) * meta.cell_size;
        double y = meta.y_min + (meta.n_rows - rows[i] - 0.5) * meta.cell_size;
        if (!ct->Transform(1, &x, &y)) throw std::runtime_error("fixture outlet transform failed");
        csv << y << ',' << x << ',' << (i + 1) << '\n';
    }
    OCTDestroyCoordinateTransformation(ct);
    std::cout << "Wrote " << path << "\n";
}

void write_forcing(const std::string& dir)
{
    constexpr double k_pi = 3.14159265358979323846;
    std::ostringstream rain, temp, pet;
    rain << "date,rainfall_mm\n";
    temp << "temperature_c\n";
    pet << "pet_mm\n";
    rain << std::setprecision(10);
    temp << std::setprecision(10);
    pet << std::setprecision(10);
    // 2001-01-01 onward; 2001 is not a leap year.
    const int month_days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int month = 0, day = 1;
    for (int t = 0; t < k_n_days; ++t)
    {
        const double p = (t % 5 == 0) ? 14.0 : ((t % 2 == 0) ? 2.5 : 0.0);
        const double tc = -6.0 + 18.0 * t / (k_n_days - 1.0) + 3.0 * std::sin(2.0 * k_pi * t / 9.0);
        const double e = std::max(0.0, 0.15 * tc);
        char date[32];
        std::snprintf(date, sizeof(date), "2001-%02d-%02d", month + 1, day);
        rain << date << ',' << p << '\n';
        temp << tc << '\n';
        pet << e << '\n';
        if (++day > month_days[month]) { day = 1; ++month; }
    }
    write_text(dir + "/hype_rainfall.csv", rain.str());
    write_text(dir + "/hype_temperature.csv", temp.str());
    write_text(dir + "/hype_pet.csv", pet.str());
}

const char* k_general =
    "name,value\n"
    "cevpam,0.3\n"
    "cevpph,20\n"
    "damp,0.3\n"
    "epotdist,4\n"
    "fsceff,0\n"
    "lp,0.75\n"
    "rivvel,0.7\n"
    "rrcs3,0.002\n"
    "tcalt,0.6\n"
    "ttpd,0.5\n"
    "ttpi,1.5\n"
    "whcsnow,0.08\n";

const char* k_landuse =
    "code,cevp,cmlt,pcluse,srrcs,ttmp\n"
    "1,0.2,3.5,0.05,0.15,0\n"
    "2,0.25,2.5,0.15,0.1,0.5\n";

const char* k_soil =
    "code,macrate,mactrinf,mactrsm,mperc1,mperc2,rrcs1,rrcs2,srrate,wcep,wcfc,wcwp\n"
    "1,0.15,12,0.6,20,8,0.25,0.02,0.1,0.14,0.22,0.08\n"
    "2,0.05,6,0.7,8,3,0.15,0.01,0.2,0.08,0.25,0.15\n";

const char* k_classes_header =
    "class_id,landuse,soil,n_layers,layer1_bottom_m,layer2_bottom_m,layer3_bottom_m,stream_depth_m,tile_depth_m\n";
const char* k_classes_rows_123 =
    "1,1,1,3,0.1,0.5,1.5,1,0\n"
    "2,2,1,2,0.15,0.6,,0.5,0\n"
    "3,2,2,1,0.3,,,0.25,0\n";
const char* k_class_row_4 =
    "4,1,2,2,0.2,0.8,,0.5,0\n";

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-hype-cli-fixture <output_dir>\n";
        return 2;
    }
    const std::string dir = argv[1];
    try
    {
        std::filesystem::create_directories(dir);
        const lsdtt::lsdtt_raster_meta meta = fixture_meta();

        write_raster(make_dem(), meta, dir + "/hype_dem.tif");

        std::vector<float> landuse(static_cast<std::size_t>(k_n_rows) * k_n_cols);
        std::vector<float> soil(landuse.size());
        for (int r = 0; r < k_n_rows; ++r)
            for (int c = 0; c < k_n_cols; ++c)
            {
                landuse[static_cast<std::size_t>(r) * k_n_cols + c] = (c < 30) ? 1.0f : 2.0f;
                soil[static_cast<std::size_t>(r) * k_n_cols + c]    = (r < 30) ? 1.0f : 2.0f;
            }
        write_raster(landuse, meta, dir + "/hype_landuse.tif");
        write_raster(soil, meta, dir + "/hype_soil.tif");

        write_outlets(meta, dir + "/hype_outlets.csv");
        write_forcing(dir);

        write_text(dir + "/hype_general_parameters.csv", k_general);
        write_text(dir + "/hype_landuse_parameters.csv", k_landuse);
        write_text(dir + "/hype_soil_parameters.csv", k_soil);
        write_text(dir + "/hype_classes.csv",
                   std::string(k_classes_header) + k_classes_rows_123 + k_class_row_4);

        // Broken tables for the failure tests.
        std::string fsceff = k_general;
        fsceff.replace(fsceff.find("fsceff,0"), 8, "fsceff,0.5");
        write_text(dir + "/hype_general_parameters_fsceff.csv", fsceff);
        write_text(dir + "/hype_landuse_parameters_unknown.csv",
                   "code,cevp,cmlt,pcluse,srrcs,ttmp,cmelt\n"
                   "1,0.2,3.5,0.05,0.15,0,1\n"
                   "2,0.25,2.5,0.15,0.1,0.5,1\n");
        write_text(dir + "/hype_classes_missing_pair.csv",
                   std::string(k_classes_header) + k_classes_rows_123);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FAIL: writing the HYPE CLI fixture: " << e.what() << "\n";
        std::cout << "RESULT: FAIL\n";
        return 1;
    }
    std::cout << "RESULT: PASS\n";
    return 0;
}
