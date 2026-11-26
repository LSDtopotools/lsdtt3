// test_proj_grid_check.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// PROJ grid availability check (src/core/lsdtt_proj_grid_check):
//   - WGS 84 -> UTM needs no grid: nothing missing, no warning;
//   - British National Grid -> UTM over Glen Coe: the most accurate operation
//     is the OSTN15 one, whether or not the grid is installed; if it is not,
//     it is the missing grid and the fallback is a different, less accurate
//     operation; if it is, the operation used is the best one;
//   - a warning is printed once per CRS pair, only when a grid is missing;
//   - the default grid list and the banner text.
// Which grids are installed depends on the machine, so the checks hold in
// both cases.
// -----------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <string>

#include <ogr_spatialref.h>

#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/core/lsdtt_proj_grid_check.hpp"

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

std::string capture_warning(const char* src, const char* dst, double x0, double y0, double x1, double y1)
{
    OGRSpatialReference a, b;
    a.SetFromUserInput(src);
    b.SetFromUserInput(dst);
    std::ostringstream buf;
    std::streambuf* old = std::cerr.rdbuf(buf.rdbuf());
    warn_if_proj_grids_missing(a, b, "test step", x0, y0, x1, y1);
    std::cerr.rdbuf(old);
    return buf.str();
}

} // namespace

int main()
{
    gdal_init();
    const std::string status = proj_grid_check_status();
    std::cout << "grid check " << status << "; " << proj_grids_description() << "\n";
    check("default grid list has OSTN15 and EGM96",
          default_proj_grids().size() == 4 && default_proj_grids()[0] == "uk_os_OSTN15_NTv2_OSGBtoETRS.tif" &&
              default_proj_grids()[1] == "us_nga_egm96_15.tif");
    check("banner text names the default count", proj_grids_description().find("default") != std::string::npos ||
                                                    proj_grids_description().find("no bundled") != std::string::npos);
    if (status != "on")
    {
        check("check is off: nothing evaluated",
              !check_proj_transformation_grids("EPSG:27700", "EPSG:32630").evaluated);
        std::cout << g_checks << " checks, " << g_failures << " failed\nRESULT: "
                  << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
        return g_failures == 0 ? 0 : 1;
    }

    const double glencoe[4] = {-5.2, 56.6, -4.9, 56.75};
    const proj_grid_check_result w = check_proj_transformation_grids("EPSG:4326", "EPSG:32630", glencoe);
    check("WGS 84 -> UTM: evaluated, no grid needed", w.evaluated && w.missing_grids.empty(), w.best_operation);

    const proj_grid_check_result b = check_proj_transformation_grids("EPSG:27700", "EPSG:32630", glencoe);
    check("BNG -> UTM: the best operation is the OSTN15 one",
          b.evaluated && b.best_operation.find("OSGB36 to WGS 84 (9)") != std::string::npos, b.best_operation);
    if (b.missing_grids.empty())
        check("OSTN15 installed: the operation used is the best one", b.used_operation == b.best_operation);
    else
        check("OSTN15 missing: it is named, and a less accurate fallback is used",
              b.missing_grids.size() == 1 && b.missing_grids[0] == "uk_os_OSTN15_NTv2_OSGBtoETRS.tif" &&
                  b.used_operation != b.best_operation && b.used_accuracy_m > b.best_accuracy_m,
              b.used_operation);

    check("no warning for WGS 84 -> UTM",
          capture_warning("EPSG:4326", "EPSG:32630", -5.2, 56.6, -4.9, 56.75).empty());
    const std::string first = capture_warning("EPSG:27700", "EPSG:32630", 220000, 750000, 230000, 760000);
    const std::string second = capture_warning("EPSG:27700", "EPSG:32630", 220000, 750000, 230000, 760000);
    if (b.missing_grids.empty())
        check("no warning when OSTN15 is installed", first.empty());
    else
        check("warning names the grid and the script",
              first.find("uk_os_OSTN15_NTv2_OSGBtoETRS.tif") != std::string::npos &&
                  first.find("fetch_proj_grids.py") != std::string::npos && first.find("test step") == 9,
              first.substr(0, 80));
    check("warned at most once per CRS pair", second.empty());

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
