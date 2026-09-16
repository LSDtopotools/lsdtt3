// test_gdal_init.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Tests for lsdtt::gdal_init() -- the PROJ/GDAL data pinning in
// src/core/lsdtt_gdal_init.cpp.
//
// set_bundled_data_root() only takes effect before the first gdal_init(), and
// gdal_init() runs once per process, so each scenario needs its own process.
// This driver therefore takes a mode on argv and CTest runs it several times.
//
//   test_gdal_init idempotent
//       gdal_init() is safe to call repeatedly and registers drivers.
//
//   test_gdal_init pinned <staging_dir>
//       Stages a proj.db into <staging_dir>/share/proj, points the bundled-data
//       root at <staging_dir>, sets PROJ_DATA to an EMPTY directory, and checks
//       that a datum transformation still gives the right answer. This is the
//       case the whole module exists for: the user's environment is wrong and
//       must be ignored.
//
//   test_gdal_init fallthrough
//       No bundled proj.db. Asserts we do NOT claim to have pinned anything and
//       that the inherited environment still works -- a conda-forge install must
//       be unaffected.
//
// THE NUMERIC CHECK
//
// The transform is WGS84 (EPSG:4326) -> OSGB36 / British National Grid
// (EPSG:27700), chosen deliberately over a UTM zone. A UTM projection is pure
// mathematics on the same datum and keeps working even when proj.db cannot be
// read, so it would not detect the failure this module guards against.
// EPSG:27700 needs the OSGB36 datum shift out of the database; without it PROJ
// silently applies a null shift and the answer is wrong by ~100 m while still
// looking like a plausible coordinate. That is exactly the "plausible and
// wrong" failure mode described in lsdtt_gdal_init.hpp.
// -----------------------------------------------------------------------------

#include "../../src/core/lsdtt_gdal_init.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

#include <gdal.h>
#include <ogr_spatialref.h>

namespace fs = std::filesystem;

namespace
{

int g_failures = 0;

void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << "\n";
    if (!ok) ++g_failures;
}

// A point in central Edinburgh, and its EPSG:27700 coordinates as PROJ
// computes them with a working proj.db. MEASURED, not derived: produced by
// this test on PROJ 8.2.1 / GDAL 3.4.1 and pinned as an oracle. The assertion
// is that the answer does not MOVE -- these digits are not independently
// authoritative, and the sanity check on them is that central Edinburgh is
// around National Grid NT 25 73.
constexpr double k_test_lon      = -3.1883;      // degrees east
constexpr double k_test_lat      = 55.9533;      // degrees north
constexpr double k_expect_east   = 325897.0;     // EPSG:27700 easting, metres
constexpr double k_expect_north  = 674001.0;     // EPSG:27700 northing, metres

// Tolerance is 50 m. Far tighter than the ~100 m error a missing OSGB36 datum
// shift produces, which is what this test is for; and far looser than the ~3 m
// spread between PROJ's OSTN15 grid and the 7-parameter Helmert fallback, so a
// PROJ upgrade that changes which one is picked does not fail the build.
constexpr double k_tol_m = 50.0;

// Transform the test point to EPSG:27700. Returns false if PROJ refuses.
bool transform_test_point(double& east, double& north)
{
    OGRSpatialReference src, dst;
    if (src.importFromEPSG(4326) != OGRERR_NONE) return false;
    if (dst.importFromEPSG(27700) != OGRERR_NONE) return false;

    // Force lon,lat argument order regardless of the CRS's declared axis order.
    src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&src, &dst);
    if (!ct) return false;

    east  = k_test_lon;
    north = k_test_lat;
    const bool ok = ct->Transform(1, &east, &north) != FALSE;
    OGRCoordinateTransformation::DestroyCT(ct);
    return ok;
}

void check_transform(const std::string& label)
{
    double east = 0.0, north = 0.0;
    if (!transform_test_point(east, north))
    {
        check(false, label + ": transform to EPSG:27700 failed outright");
        return;
    }

    std::cout << "        easting=" << east << " northing=" << north << "\n";

    const double de = std::fabs(east  - k_expect_east);
    const double dn = std::fabs(north - k_expect_north);
    check(de < k_tol_m && dn < k_tol_m,
          label + ": EPSG:27700 coordinates within " +
          std::to_string(static_cast<int>(k_tol_m)) + " m of the pinned value");
}

// Point PROJ_DATA and PROJ_LIB at a directory with no proj.db in it.
void poison_proj_env(const fs::path& empty_dir)
{
    fs::create_directories(empty_dir);
#if defined(_WIN32)
    _putenv_s("PROJ_DATA", empty_dir.string().c_str());
    _putenv_s("PROJ_LIB",  empty_dir.string().c_str());
#else
    ::setenv("PROJ_DATA", empty_dir.string().c_str(), 1);
    ::setenv("PROJ_LIB",  empty_dir.string().c_str(), 1);
#endif
}

int run_idempotent()
{
    lsdtt::gdal_init();
    const int n1 = GDALGetDriverCount();
    lsdtt::gdal_init();
    lsdtt::gdal_init();
    const int n2 = GDALGetDriverCount();

    check(n1 > 0, "gdal_init() registers drivers");
    check(n1 == n2, "repeated gdal_init() does not re-register");
    check(!lsdtt::gdal_version().empty(), "gdal_version() is non-empty");

    // set_bundled_data_root must refuse once initialisation has happened,
    // rather than silently pretending to take effect.
    check(!lsdtt::set_bundled_data_root("/nonexistent"),
          "set_bundled_data_root() returns false after init");
    return 0;
}

int run_pinned(const fs::path& staging)
{
    // Find a real proj.db to stage. The build host has one; if it cannot be
    // found there is nothing to test, so say so rather than passing silently.
    const char* env_src = std::getenv("LSDTT_TEST_PROJ_DB");
    if (!env_src || !*env_src || !fs::exists(env_src))
    {
        std::cerr << "SKIP: LSDTT_TEST_PROJ_DB not set to an existing proj.db\n";
        return 77;  // CTest SKIP_RETURN_CODE
    }

    const fs::path share_proj = staging / "share" / "proj";
    fs::create_directories(share_proj);

    // Copy the whole directory: proj.db alone is enough for EPSG:27700, but a
    // real bundle carries the transformation grids beside it.
    const fs::path src_dir = fs::path(env_src).parent_path();
    for (const auto& e : fs::directory_iterator(src_dir))
    {
        if (e.is_regular_file())
            fs::copy_file(e.path(), share_proj / e.path().filename(),
                          fs::copy_options::overwrite_existing);
    }

    poison_proj_env(staging / "empty_proj_dir");

    check(lsdtt::set_bundled_data_root(staging.string()),
          "set_bundled_data_root() accepted before init");

    const std::string desc = lsdtt::proj_data_description();
    std::cout << "        proj data: " << desc << "\n";
    check(desc.rfind("bundled:", 0) == 0,
          "bundled proj.db wins over a poisoned PROJ_DATA");

    check_transform("pinned");
    return g_failures == 0 ? 0 : 1;
}

int run_fallthrough()
{
    // No set_bundled_data_root(), and the test binary does not sit next to a
    // share/proj, so nothing should be pinned and the inherited environment
    // must still work.
    const std::string desc = lsdtt::proj_data_description();
    std::cout << "        proj data: " << desc << "\n";
    check(desc.rfind("bundled:", 0) != 0,
          "nothing claimed as bundled when no bundled proj.db exists");

    check_transform("fallthrough");
    return g_failures == 0 ? 0 : 1;
}

} // namespace

int main(int argc, char* argv[])
{
    const std::string mode = (argc > 1) ? argv[1] : "";

    if (mode == "idempotent")
    {
        run_idempotent();
    }
    else if (mode == "pinned")
    {
        if (argc < 3)
        {
            std::cerr << "usage: test_gdal_init pinned <staging_dir>\n";
            return 2;
        }
        const int rc = run_pinned(argv[2]);
        if (rc == 77) return 77;
    }
    else if (mode == "fallthrough")
    {
        run_fallthrough();
    }
    else
    {
        std::cerr << "usage: test_gdal_init "
                     "{idempotent|pinned <dir>|fallthrough}\n";
        return 2;
    }

    std::cout << (g_failures == 0 ? "OK" : "FAILURES") << ": " << mode
              << " (" << g_failures << " failed)\n";
    return g_failures == 0 ? 0 : 1;
}
