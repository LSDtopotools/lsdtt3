// lsdtt_usgs_elevation.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// See lsdtt_usgs_elevation.hpp for why this exists and why it is US-only.

#include "lsdtt_usgs_elevation.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <sstream>

#include <cpl_conv.h>
#include <cpl_error.h>
#include <cpl_http.h>
#include <cpl_string.h>
#include <ogr_spatialref.h>

namespace lsdtt
{

namespace
{

// The international foot. us_survey_feet differs by 2 ppm, which is far below
// anything this comparison can resolve, so the verdict names `feet` for both --
// the point is the factor of 3.28, not which foot.
constexpr double kMetresPerFoot = 0.3048;

// A ratio must land within this fraction of a hypothesis to be accepted. 12% is
// wide enough to absorb 3DEP-vs-lidar disagreement on steep ground and horizontal
// registration error, and far narrower than the 228% gap between the two
// hypotheses it has to separate.
constexpr double kRatioTolerance = 0.12;

// Silence GDAL's own error printing for the duration of a scope: a failed
// diagnostic query must not look like a failed run.
struct quiet_gdal
{
    quiet_gdal()  { CPLPushErrorHandler(CPLQuietErrorHandler); }
    ~quiet_gdal() { CPLPopErrorHandler(); }
};

// Pull "value":<number> out of an EPQS JSON reply without a JSON dependency.
// The service returns a flat object and has done so across versions; anything
// unexpected simply fails to parse and the point is dropped.
bool parse_epqs_value(const char* body, double& out)
{
    if (body == nullptr) return false;
    const std::string s(body);

    const std::string key = "\"value\"";
    std::size_t p = s.find(key);
    if (p == std::string::npos) return false;
    p = s.find(':', p + key.size());
    if (p == std::string::npos) return false;
    ++p;

    // The value may be bare (12.34) or quoted ("12.34").
    while (p < s.size() && (s[p] == ' ' || s[p] == '"')) ++p;

    try
    {
        std::size_t used = 0;
        const double v = std::stod(s.substr(p, 64), &used);
        if (used == 0) return false;
        // EPQS reports -1000000 for a point outside coverage.
        if (v < -100000.0) return false;
        out = v;
        return true;
    }
    catch (const std::exception&) { return false; }
}

} // namespace


bool in_usgs_coverage(double lon_deg, double lat_deg)
{
    const auto box = [&](double w, double s, double e, double n)
    { return lon_deg >= w && lon_deg <= e && lat_deg >= s && lat_deg <= n; };

    return box(-125.0, 24.0, -66.5, 49.5)     // CONUS
        || box(-170.0, 51.0, -129.0, 72.0)    // Alaska
        || box(-161.0, 18.5, -154.5, 22.5)    // Hawaii
        || box(-68.0, 17.8, -65.2, 18.6);     // Puerto Rico / USVI
}


usgs_unit_check check_vertical_units_against_usgs(const lsdtt_raster_cpu& dem,
                                                 int n_points,
                                                 int timeout_s)
{
    usgs_unit_check out;

    const lsdtt_raster_meta& meta = dem.meta();

    if (meta.crs_wkt.empty())
    {
        out.message = "the DEM carries no CRS, so its location is unknown and "
                      "3DEP cannot be queried.";
        return out;
    }

    // ---- projected -> WGS84 lon/lat --------------------------------------
    OGRSpatialReference src, wgs84;
    {
        quiet_gdal quiet;
        if (src.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE ||
            wgs84.importFromEPSG(4326) != OGRERR_NONE)
        {
            out.message = "the DEM's CRS could not be interpreted, so 3DEP "
                          "cannot be queried.";
            return out;
        }
    }
    src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* tf = nullptr;
    {
        quiet_gdal quiet;
        tf = OGRCreateCoordinateTransformation(&src, &wgs84);
    }
    if (tf == nullptr)
    {
        out.message = "the DEM's CRS could not be transformed to WGS84, so "
                      "3DEP cannot be queried.";
        return out;
    }

    // Coverage is decided on the CENTRE, once, before any request goes out.
    {
        double cx = meta.x_min + 0.5 * meta.n_cols * meta.cell_size;
        double cy = meta.y_min + 0.5 * meta.n_rows * meta.cell_size;
        if (!tf->Transform(1, &cx, &cy) || !in_usgs_coverage(cx, cy))
        {
            OGRCoordinateTransformation::DestroyCT(tf);
            out.message = "this DEM is outside USGS 3DEP coverage, so the "
                          "vertical unit cannot be checked automatically. "
                          "Note that feet-on-a-metre-grid is a US delivery "
                          "convention, so this is unlikely to be the problem "
                          "here.";
            return out;
        }
    }

    if (!CPLHTTPEnabled())
    {
        OGRCoordinateTransformation::DestroyCT(tf);
        out.message = "GDAL was built without HTTP support (no curl), so 3DEP "
                      "cannot be queried.";
        return out;
    }

    // ---- pick valid pixels ------------------------------------------------
    const std::vector<float>&   d = dem.get_data();
    const std::vector<uint8_t>& m = dem.get_mask();

    std::vector<std::size_t> valid;
    valid.reserve(d.size() / 8);
    for (std::size_t i = 0; i < d.size(); ++i)
        if (m[i]) valid.push_back(i);

    if (valid.empty())
    {
        OGRCoordinateTransformation::DestroyCT(tf);
        out.message = "the DEM has no valid pixels to compare.";
        return out;
    }

    // Fixed seed: a diagnostic that gives a different answer on a rerun is not
    // a diagnostic.
    std::mt19937 rng(12345);
    std::shuffle(valid.begin(), valid.end(), rng);
    const int n_try = std::min<int>(n_points, static_cast<int>(valid.size()));
    out.n_requested = n_try;

    // ---- query ------------------------------------------------------------
    std::vector<double> ratios, err_m, err_ft;

    for (int k = 0; k < n_try; ++k)
    {
        const std::size_t idx = valid[static_cast<std::size_t>(k)];
        const int r = static_cast<int>(idx / meta.n_cols);
        const int c = static_cast<int>(idx % meta.n_cols);

        // Pixel centre. Row 0 is the TOP row, so northing counts down from y_max.
        double x = meta.x_min + (c + 0.5) * meta.cell_size;
        double y = meta.y_max() - (r + 0.5) * meta.cell_size;
        if (!tf->Transform(1, &x, &y)) continue;

        std::ostringstream url;
        url.precision(8);
        url << std::fixed
            << "https://epqs.nationalmap.gov/v1/json?x=" << x
            << "&y=" << y << "&units=Meters&wkid=4326";

        char** opts = nullptr;
        opts = CSLSetNameValue(opts, "TIMEOUT", std::to_string(timeout_s).c_str());
        opts = CSLSetNameValue(opts, "MAX_RETRY", "1");

        CPLHTTPResult* resp = nullptr;
        {
            quiet_gdal quiet;
            resp = CPLHTTPFetch(url.str().c_str(), opts);
        }
        CSLDestroy(opts);

        if (resp == nullptr) continue;

        double usgs_m = 0.0;
        const bool ok = (resp->nStatus == 0) && (resp->nDataLen > 0) &&
                        parse_epqs_value(
                            reinterpret_cast<const char*>(resp->pabyData), usgs_m);
        CPLHTTPDestroyResult(resp);

        if (!ok || std::fabs(usgs_m) < 1e-6) continue;

        const double dem_v = static_cast<double>(d[idx]);
        ratios.push_back(dem_v / usgs_m);
        err_m .push_back(std::fabs(dem_v - usgs_m));
        err_ft.push_back(std::fabs(dem_v * kMetresPerFoot - usgs_m));
    }

    OGRCoordinateTransformation::DestroyCT(tf);

    out.n_returned = static_cast<int>(ratios.size());

    // Fewer than three replies is not a verdict. Say so rather than deciding on
    // one or two points -- an unreachable service must not become a conclusion,
    // and three is the smallest median that can outvote a single bad pixel.
    if (out.n_returned < 3)
    {
        out.verdict = vertical_unit_verdict::no_response;
        out.message = "USGS 3DEP returned too few points (" +
                      std::to_string(out.n_returned) + " of " +
                      std::to_string(out.n_requested) +
                      ") to decide the vertical unit. This is usually the "
                      "service being unreachable rather than anything about "
                      "the DEM.";
        return out;
    }

    out.median_ratio           = compute_median_upper_of_copy(ratios, 0.0);
    out.median_abs_err_as_metres = compute_median_upper_of_copy(err_m, 0.0);
    out.median_abs_err_as_feet   = compute_median_upper_of_copy(err_ft, 0.0);

    const double feet_ratio = 1.0 / kMetresPerFoot;   // 3.28084
    std::ostringstream msg;
    msg.precision(3);
    msg << std::fixed;

    if (std::fabs(out.median_ratio - 1.0) <= kRatioTolerance)
    {
        out.verdict = vertical_unit_verdict::metres;
        msg << "USGS 3DEP agrees with this DEM as it stands ("
            << out.n_returned << " points, median ratio "
            << out.median_ratio << ", median error "
            << out.median_abs_err_as_metres << " m).\n"
            << "  THE ELEVATIONS ARE ALREADY IN METRES. Leave vertical_units "
               "at metres.\n"
            << "  The steepness flagged above is therefore real terrain, not a "
               "unit error.";
    }
    else if (std::fabs(out.median_ratio - feet_ratio) <=
             kRatioTolerance * feet_ratio)
    {
        out.verdict = vertical_unit_verdict::feet;
        msg << "USGS 3DEP says THIS DEM'S ELEVATIONS ARE IN FEET ("
            << out.n_returned << " points, median ratio "
            << out.median_ratio << " vs 3.281 expected for feet).\n"
            << "    median |DEM - 3DEP|           = "
            << out.median_abs_err_as_metres << " m\n"
            << "    median |DEM x 0.3048 - 3DEP|  = "
            << out.median_abs_err_as_feet << " m   <- agrees\n"
            << "  SET  vertical_units: feet  AND RERUN. Until you do, every "
               "slope, curvature\n"
            << "  and relief in this run is too large by a factor of 3.28.";
    }
    else
    {
        out.verdict = vertical_unit_verdict::inconsistent;
        msg << "USGS 3DEP does not match this DEM under either hypothesis ("
            << out.n_returned << " points, median ratio "
            << out.median_ratio << "; 1.0 would mean metres, 3.281 feet).\n"
            << "  Do NOT set vertical_units from this. Something else differs "
               "-- a vertical\n"
            << "  datum offset, a horizontal registration error, or a DEM that "
               "is not ground\n"
            << "  surface (first-return canopy). Investigate before converting "
               "anything.";
    }

    out.message = msg.str();
    return out;
}

} // namespace lsdtt
