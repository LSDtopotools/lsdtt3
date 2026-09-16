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

} // namespace lsdtt
