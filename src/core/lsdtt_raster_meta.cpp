// lsdtt_raster_meta.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// The CRS methods of lsdtt_raster_meta.
//
// These live here rather than inline in the header because they need
// <ogr_spatialref.h>, and lsdtt_raster_meta.hpp is included by every compute
// object in the tree. Keeping the OGR include here is what stops it reaching
// several hundred translation units that only carry the struct around.

#include "lsdtt_raster_meta.hpp"

#include "lsdtt_gdal_init.hpp"

#include <cpl_conv.h>        // CPLFree
#include <ogr_spatialref.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace lsdtt {

// Sets crs_wkt and epsg_code from an EPSG code (e.g. 32601 for UTM zone 1
// North on WGS84). Used to give synthetically generated rasters a sensible
// default projected CRS. Throws std::runtime_error if the EPSG code cannot
// be resolved.
void lsdtt_raster_meta::set_crs_from_epsg(int epsg)
{
    gdal_init();
    OGRSpatialReference srs;
    if (srs.importFromEPSG(epsg) != OGRERR_NONE)
    {
        throw std::runtime_error(
            std::string("lsdtt_raster_meta::set_crs_from_epsg: failed to "
                        "import EPSG:") + std::to_string(epsg));
    }
    char* wkt_raw = nullptr;
    if (srs.exportToWkt(&wkt_raw) != OGRERR_NONE || !wkt_raw)
    {
        CPLFree(wkt_raw);
        throw std::runtime_error(
            "lsdtt_raster_meta::set_crs_from_epsg: failed to export WKT");
    }
    crs_wkt   = wkt_raw;
    epsg_code = epsg;
    CPLFree(wkt_raw);
}

// Places a synthetic grid at an arbitrary geographic anchor. The output
// grid is centred on (centre_lat_deg, centre_lon_deg) in the selected
// projected CRS. If epsg == 0, a WGS84 UTM EPSG is derived from
// (centre_lon_deg, centre_lat_deg).
//
// Throws std::runtime_error when the anchor is non-finite or when
// epsg == 0 and |centre_lat_deg| > 84 (UTM undefined at polar latitudes).
// n_rows, n_cols and cell_size must already be set.
void lsdtt_raster_meta::set_synthetic_georeference(double centre_lat_deg,
                                double centre_lon_deg,
                                int    epsg)
{
    if (!std::isfinite(centre_lat_deg) || !std::isfinite(centre_lon_deg))
    {
        throw std::runtime_error(
            "lsdtt_raster_meta::set_synthetic_georeference: anchor latitude "
            "and longitude must be finite");
    }

    int target_epsg = epsg;
    if (target_epsg == 0)
    {
        if (std::abs(centre_lat_deg) > 84.0)
        {
            throw std::runtime_error(
                "lsdtt_raster_meta::set_synthetic_georeference: "
                "epsg=0 derives UTM, but UTM is undefined for |latitude| > 84. "
                "Set projection_epsg explicitly (e.g. 3413 north, 3976 south).");
        }

        // Keep this logic aligned with utm_epsg_for_point in
        // lsdtt_raster_fetch.cpp.
        double normalised = std::fmod(centre_lon_deg + 180.0, 360.0);
        if (normalised < 0.0) normalised += 360.0;

        int zone = static_cast<int>(normalised / 6.0) + 1;
        if (zone > 60) zone = 60;
        if (zone < 1)  zone = 1;
        target_epsg = (centre_lat_deg >= 0.0 ? 32600 : 32700) + zone;
    }

    set_crs_from_epsg(target_epsg);

    OGRSpatialReference wgs84;
    OGRSpatialReference target;
    if (wgs84.importFromEPSG(4326) != OGRERR_NONE ||
        target.importFromEPSG(target_epsg) != OGRERR_NONE)
    {
        return;
    }
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    target.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&wgs84, &target);
    if (!ct) return;

    double e = centre_lon_deg;
    double n = centre_lat_deg;
    const bool ok = ct->Transform(1, &e, &n) != 0;
    OGRCoordinateTransformation::DestroyCT(ct);
    if (!ok) return;

    x_min = e - 0.5 * static_cast<double>(n_cols) * cell_size;
    y_min = n - 0.5 * static_cast<double>(n_rows) * cell_size;
}

// Gives a synthetically generated raster (no input DEM) a sensible default
// georeference: the grid is centred on Teanu Island, Solomon Islands
// (11.633 S, 166.967 E) - the real-world location to which Zwift's fictional
// "Watopia" is geocoded for Strava. This sits comfortably inside UTM zone
// 58S (EPSG:32758), well away from the antimeridian and the poles, so
// generated DEMs carry realistic projected coordinates and round-trip
// cleanly to and from lat/lon. n_rows, n_cols and cell_size must already be
// set. Falls back to no CRS if the transform fails.
void lsdtt_raster_meta::set_default_synthetic_georeference()
{
    constexpr double watopia_lon = 166.967;  // Teanu Island
    constexpr double watopia_lat = -11.633;
    constexpr int    epsg        = 32758;    // WGS 84 / UTM zone 58S

    try
    {
        set_synthetic_georeference(watopia_lat, watopia_lon, epsg);
    }
    catch (const std::exception&)
    {
        crs_wkt = "";
        return;
    }
}

// Returns true if the CRS is geographic (latitude/longitude).
// Returns false if projected or if no CRS is defined.
//
bool lsdtt_raster_meta::is_geographic() const
{
    gdal_init();
    if (crs_wkt.empty()) return false;
    OGRSpatialReference srs;
    if (srs.importFromWkt(crs_wkt.c_str()) != OGRERR_NONE) return false;
    return srs.IsGeographic() != 0;
}

// Returns the name of the linear unit of the CRS (e.g. "metre", "foot",
// "US survey foot"). Returns empty string if no CRS or if geographic.
std::string lsdtt_raster_meta::linear_unit_name() const
{
    gdal_init();
    if (crs_wkt.empty()) return "";
    OGRSpatialReference srs;
    if (srs.importFromWkt(crs_wkt.c_str()) != OGRERR_NONE) return "";
    if (srs.IsGeographic()) return "";
    const char* unit_name = nullptr;
    srs.GetLinearUnits(&unit_name);
    return unit_name ? std::string(unit_name) : "";
}

// Returns the conversion factor from the CRS linear unit to metres.
// Returns 1.0 if the unit is metres, ~0.3048 for feet, etc.
// Returns 0.0 if the CRS is geographic or unknown.
double lsdtt_raster_meta::linear_unit_to_metres() const
{
    gdal_init();
    if (crs_wkt.empty()) return 0.0;
    OGRSpatialReference srs;
    if (srs.importFromWkt(crs_wkt.c_str()) != OGRERR_NONE) return 0.0;
    if (srs.IsGeographic()) return 0.0;
    double to_metres = srs.GetLinearUnits(nullptr);
    return to_metres;
}

// Returns true if the linear unit is metres (within floating point tolerance).
bool lsdtt_raster_meta::is_in_metres() const
{
    double factor = linear_unit_to_metres();
    return std::abs(factor - 1.0) < 1e-6;
}

// Throws std::runtime_error with a clear message if the CRS is geographic
// or if the linear unit is not metres.
// Call this at the start of any analysis function that requires projected
// data in metres.
void lsdtt_raster_meta::check_projected_metres() const
{
    if (is_geographic())
    {
        throw std::runtime_error(
            "You are using a geographic coordinate system but this program "
            "only works on projected data. Please convert your raster to a "
            "projected coordinate system like UTM."
        );
    }
    if (!crs_wkt.empty() && !is_in_metres())
    {
        throw std::runtime_error(
            "Your units are not in metres. Are you using feet? Please "
            "convert your raster to a unit that is not only used by lunatics."
        );
    }
}

// -----------------------------------------------------------------------------
// max_linear_scale_error
// -----------------------------------------------------------------------------
namespace
{

// Vincenty (1975) inverse: ellipsoidal geodesic length (metres) between two
// points given in degrees. a = semi-major axis, f = flattening. Converges for
// the short (10-cell) baselines used here; never near-antipodal.
double vincenty_distance(double lat1, double lon1, double lat2, double lon2, double a, double f)
{
    const double deg = M_PI / 180.0;
    const double b = a * (1.0 - f);
    const double L = (lon2 - lon1) * deg;
    const double U1 = std::atan((1.0 - f) * std::tan(lat1 * deg));
    const double U2 = std::atan((1.0 - f) * std::tan(lat2 * deg));
    const double sinU1 = std::sin(U1), cosU1 = std::cos(U1);
    const double sinU2 = std::sin(U2), cosU2 = std::cos(U2);
    double lambda = L, sin_sigma = 0, cos_sigma = 0, sigma = 0, cos2_alpha = 0, cos_2sm = 0;
    for (int it = 0; it < 200; ++it)
    {
        const double sin_l = std::sin(lambda), cos_l = std::cos(lambda);
        sin_sigma = std::sqrt((cosU2 * sin_l) * (cosU2 * sin_l) +
                              (cosU1 * sinU2 - sinU1 * cosU2 * cos_l) *
                                  (cosU1 * sinU2 - sinU1 * cosU2 * cos_l));
        if (sin_sigma == 0.0) return 0.0;
        cos_sigma = sinU1 * sinU2 + cosU1 * cosU2 * cos_l;
        sigma = std::atan2(sin_sigma, cos_sigma);
        const double sin_alpha = cosU1 * cosU2 * sin_l / sin_sigma;
        cos2_alpha = 1.0 - sin_alpha * sin_alpha;
        cos_2sm = (cos2_alpha != 0.0) ? cos_sigma - 2.0 * sinU1 * sinU2 / cos2_alpha : 0.0;
        const double C = f / 16.0 * cos2_alpha * (4.0 + f * (4.0 - 3.0 * cos2_alpha));
        const double prev = lambda;
        lambda = L + (1.0 - C) * f * sin_alpha *
                         (sigma + C * sin_sigma *
                                      (cos_2sm + C * cos_sigma * (-1.0 + 2.0 * cos_2sm * cos_2sm)));
        if (std::fabs(lambda - prev) < 1e-13) break;
    }
    const double u2 = cos2_alpha * (a * a - b * b) / (b * b);
    const double A = 1.0 + u2 / 16384.0 * (4096.0 + u2 * (-768.0 + u2 * (320.0 - 175.0 * u2)));
    const double B = u2 / 1024.0 * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)));
    const double d_sigma =
        B * sin_sigma *
        (cos_2sm + B / 4.0 *
                       (cos_sigma * (-1.0 + 2.0 * cos_2sm * cos_2sm) -
                        B / 6.0 * cos_2sm * (-3.0 + 4.0 * sin_sigma * sin_sigma) *
                            (-3.0 + 4.0 * cos_2sm * cos_2sm)));
    return b * A * (sigma - d_sigma);
}

} // namespace

double lsdtt_raster_meta::max_linear_scale_error() const
{
    gdal_init();
    if (crs_wkt.empty()) return -1.0;
    OGRSpatialReference proj;
    if (proj.importFromWkt(crs_wkt.c_str()) != OGRERR_NONE)
        throw std::runtime_error("max_linear_scale_error: cannot parse the CRS");
    if (proj.IsGeographic()) return -1.0;
    proj.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    std::unique_ptr<OGRSpatialReference> geog(proj.CloneGeogCS());
    if (!geog) throw std::runtime_error("max_linear_scale_error: CRS has no geographic CRS");
    geog->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRErr err = OGRERR_NONE;
    const double a = geog->GetSemiMajor(&err);
    const double inv_f = geog->GetInvFlattening(&err);
    const double f = inv_f > 0.0 ? 1.0 / inv_f : 0.0;
    std::unique_ptr<OGRCoordinateTransformation> ct(
        OGRCreateCoordinateTransformation(&proj, geog.get()));
    if (!ct) throw std::runtime_error("max_linear_scale_error: cannot transform to geographic");

    // Distances on the grid are in the CRS linear unit; compare in metres.
    const double unit = proj.GetLinearUnits();
    const double step = 10.0 * cell_size;
    const double xs[3] = {x_min, x_min + 0.5 * n_cols * cell_size, x_min + n_cols * cell_size};
    const double ys[3] = {y_min, y_min + 0.5 * n_rows * cell_size, y_min + n_rows * cell_size};
    const int pts[5][2] = {{1, 1}, {0, 0}, {2, 0}, {0, 2}, {2, 2}};

    double worst = 0.0;
    for (const auto& p : pts)
    {
        // Step inwards from the corners so the baseline stays on the grid.
        const double x0 = xs[p[0]], y0 = ys[p[1]];
        const double dx = (p[0] == 2) ? -step : step;
        const double dy = (p[1] == 2) ? -step : step;
        double px[3] = {x0, x0 + dx, x0};
        double py[3] = {y0, y0, y0 + dy};
        if (!ct->Transform(3, px, py))
            throw std::runtime_error("max_linear_scale_error: coordinate transform failed");
        // After the transform px = longitude, py = latitude (traditional order).
        const double k = std::fabs(dx) * unit / vincenty_distance(py[0], px[0], py[1], px[1], a, f);
        const double h = std::fabs(dy) * unit / vincenty_distance(py[0], px[0], py[2], px[2], a, f);
        worst = std::max({worst, std::fabs(k - 1.0), std::fabs(h - 1.0)});
    }
    return worst;
}

std::string lsdtt_raster_meta::crs_name() const
{
    gdal_init();
    if (crs_wkt.empty()) return "none";
    OGRSpatialReference srs;
    if (srs.importFromWkt(crs_wkt.c_str()) != OGRERR_NONE) return "unparsed";
    const char* n = srs.GetName();
    return n ? std::string(n) : std::string("unnamed");
}

} // namespace lsdtt
