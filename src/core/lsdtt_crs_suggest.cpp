// lsdtt_crs_suggest.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// See lsdtt_crs_suggest.hpp for full documentation.

#include "lsdtt_crs_suggest.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include <ogr_spatialref.h>
#include <ogr_geometry.h>

namespace lsdtt {

// =============================================================================
// Internal helpers
// =============================================================================

namespace {

// -----------------------------------------------------------------------------
// lon_to_utm_zone
//
// Returns the standard UTM zone number (1–60) for a given longitude.
// Does not handle the Norway/Svalbard UTM exceptions — those areas will use
// the polar stereographic path anyway.
// -----------------------------------------------------------------------------
static int lon_to_utm_zone(double lon)
{
    // Normalise to [-180, 180)
    while (lon >= 180.0)  lon -= 360.0;
    while (lon < -180.0)  lon += 360.0;
    return static_cast<int>(std::floor((lon + 180.0) / 6.0)) % 60 + 1;
}

// -----------------------------------------------------------------------------
// utm_epsg
//
// Returns the EPSG code for a UTM zone + hemisphere.
//   North: 32601–32660
//   South: 32701–32760
// -----------------------------------------------------------------------------
static int utm_epsg(int zone, bool north)
{
    return (north ? 32600 : 32700) + zone;
}

// -----------------------------------------------------------------------------
// make_utm_wkt
//
// Build WKT for a UTM CRS via OGR, given zone and hemisphere.
// Throws std::runtime_error on failure.
// -----------------------------------------------------------------------------
static std::string make_utm_wkt(int zone, bool north)
{
    OGRSpatialReference srs;
    srs.SetWellKnownGeogCS("WGS84");
    srs.SetUTM(zone, north ? TRUE : FALSE);

    char* wkt_raw = nullptr;
    if (srs.exportToWkt(&wkt_raw) != OGRERR_NONE || !wkt_raw)
    {
        CPLFree(wkt_raw);
        throw std::runtime_error("lsdtt_suggest_crs: failed to build UTM WKT");
    }
    std::string result(wkt_raw);
    CPLFree(wkt_raw);
    return result;
}

// -----------------------------------------------------------------------------
// make_polar_wkt
//
// Return WKT for NSIDC polar stereographic CRS by importing from EPSG.
// EPSG:3413 — North  (true scale at 70°N, lon of origin 45°W)
// EPSG:3976 — South  (true scale at 71°S, lon of origin 0°)
// -----------------------------------------------------------------------------
static std::string make_polar_wkt(bool north)
{
    int epsg = north ? 3413 : 3976;
    OGRSpatialReference srs;
    if (srs.importFromEPSG(epsg) != OGRERR_NONE)
    {
        throw std::runtime_error(
            std::string("lsdtt_suggest_crs: failed to import polar EPSG:")
            + std::to_string(epsg)
        );
    }
    char* wkt_raw = nullptr;
    if (srs.exportToWkt(&wkt_raw) != OGRERR_NONE || !wkt_raw)
    {
        CPLFree(wkt_raw);
        throw std::runtime_error("lsdtt_suggest_crs: failed to export polar WKT");
    }
    std::string result(wkt_raw);
    CPLFree(wkt_raw);
    return result;
}

// -----------------------------------------------------------------------------
// make_laea_wkt
//
// Build a custom LAEA WKT centred on (centre_lon, centre_lat), based on
// WGS84, with metres as the linear unit.
// -----------------------------------------------------------------------------
static std::string make_laea_wkt(double centre_lon, double centre_lat)
{
    OGRSpatialReference srs;
    srs.SetWellKnownGeogCS("WGS84");
    srs.SetLAEA(
        centre_lat,   // latitude of origin
        centre_lon,   // central meridian
        0.0,          // false easting
        0.0           // false northing
    );
    srs.SetLinearUnits("metre", 1.0);

    char* wkt_raw = nullptr;
    if (srs.exportToWkt(&wkt_raw) != OGRERR_NONE || !wkt_raw)
    {
        CPLFree(wkt_raw);
        throw std::runtime_error("lsdtt_suggest_crs: failed to build LAEA WKT");
    }
    std::string result(wkt_raw);
    CPLFree(wkt_raw);
    return result;
}

// -----------------------------------------------------------------------------
// unproject_corners_to_geographic
//
// Given the meta for a raster (which may be in any CRS, including geographic),
// fills in the six geographic extent fields on the suggestion struct.
//
// Strategy:
//   • If already geographic, corners are trivially computed from x_min/y_min
//     and x_max()/y_max().
//   • If projected (including Web Mercator), each corner is transformed to
//     WGS84 lon/lat using OGR coordinate transformation.
//
// Populated fields: min_lon, max_lon, min_lat, max_lat, centre_lon, centre_lat
// -----------------------------------------------------------------------------
static void unproject_corners_to_geographic(
    const lsdtt_raster_meta& meta,
    lsdtt_crs_suggestion&    sug)
{
    // Four corners in raster pixel space → projected coordinates
    const double x0 = meta.x_min;
    const double y0 = meta.y_min;
    const double x1 = meta.x_max();
    const double y1 = meta.y_max();

    std::array<double, 4> xs = { x0, x1, x0, x1 };
    std::array<double, 4> ys = { y0, y0, y1, y1 };

    OGRSpatialReference src_srs;
    if (src_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_suggest_crs: cannot parse input CRS WKT");
    }

    // If already geographic, corners are already lon/lat
    if (src_srs.IsGeographic())
    {
        // x_min/x_max = lon, y_min/y_max = lat  (standard for EPSG:4326)
        sug.min_lon = x0;
        sug.max_lon = x1;
        sug.min_lat = y0;
        sug.max_lat = y1;
    }
    else
    {
        // Transform projected corners to WGS84
        OGRSpatialReference wgs84;
        wgs84.SetWellKnownGeogCS("WGS84");
        wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        OGRCoordinateTransformation* ct =
            OGRCreateCoordinateTransformation(&src_srs, &wgs84);
        if (!ct)
        {
            throw std::runtime_error(
                "lsdtt_suggest_crs: cannot create coordinate transformation "
                "to WGS84");
        }

        for (int i = 0; i < 4; ++i)
        {
            double lon = xs[i];
            double lat = ys[i];
            if (!ct->Transform(1, &lon, &lat))
            {
                OGRCoordinateTransformation::DestroyCT(ct);
                throw std::runtime_error(
                    "lsdtt_suggest_crs: corner coordinate transformation failed");
            }
            xs[i] = lon;
            ys[i] = lat;
        }
        OGRCoordinateTransformation::DestroyCT(ct);

        sug.min_lon = *std::min_element(xs.begin(), xs.end());
        sug.max_lon = *std::max_element(xs.begin(), xs.end());
        sug.min_lat = *std::min_element(ys.begin(), ys.end());
        sug.max_lat = *std::max_element(ys.begin(), ys.end());
    }

    sug.centre_lon = 0.5 * (sug.min_lon + sug.max_lon);
    sug.centre_lat = 0.5 * (sug.min_lat + sug.max_lat);
}

// -----------------------------------------------------------------------------
// count_utm_zones
//
// Populates utm_zone_min, utm_zone_max, utm_zone_count on the suggestion.
// Handles anti-meridian crossing by checking both a coarse sweep and the
// corner zones, but does not attempt to handle split-zone edge cases
// perfectly — the warning is conservative.
// -----------------------------------------------------------------------------
static void count_utm_zones(lsdtt_crs_suggestion& sug)
{
    std::set<int> zones;

    // Sample lon at 6° intervals across the raster extent
    double lon = sug.min_lon;
    while (lon <= sug.max_lon + 1e-9)
    {
        zones.insert(lon_to_utm_zone(lon));
        lon += 6.0;
    }
    // Ensure the exact max is included
    zones.insert(lon_to_utm_zone(sug.max_lon));

    sug.utm_zone_min   = *zones.begin();
    sug.utm_zone_max   = *zones.rbegin();
    sug.utm_zone_count = static_cast<int>(zones.size());
}

} // anonymous namespace

// =============================================================================
// Public API
// =============================================================================

// -----------------------------------------------------------------------------
// lsdtt_suggest_crs
// -----------------------------------------------------------------------------
lsdtt_crs_suggestion lsdtt_suggest_crs(
    const lsdtt_raster_meta& meta,
    bool offer_laea)
{
    if (meta.crs_wkt.empty())
    {
        throw std::runtime_error(
            "lsdtt_suggest_crs: input raster has no CRS defined. "
            "Assign a CRS before calling this function.");
    }

    lsdtt_crs_suggestion sug;

    // ------------------------------------------------------------------
    // 1. Parse input CRS
    // ------------------------------------------------------------------
    OGRSpatialReference src_srs;
    if (src_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_suggest_crs: cannot parse input CRS WKT");
    }

    // ------------------------------------------------------------------
    // 2. Check whether the raster is already in a suitable projected CRS
    //    (not geographic and not Web Mercator, which severely distorts
    //    areas and distances at high latitudes).
    //
    //    Web Mercator = EPSG:3857.  We detect it by checking the EPSG code
    //    embedded in the WKT.
    // ------------------------------------------------------------------
    bool is_web_mercator = false;
    {
        const char* auth_name = src_srs.GetAuthorityName(nullptr);
        const char* auth_code = src_srs.GetAuthorityCode(nullptr);
        if (auth_name && auth_code)
        {
            if (std::string(auth_name) == "EPSG" &&
                std::string(auth_code) == "3857")
            {
                is_web_mercator = true;
            }
        }
    }

    if (!src_srs.IsGeographic() && !is_web_mercator)
    {
        // Already a sensible projected CRS — nothing to recommend
        sug.already_projected = true;
        sug.type              = CrsRecommendationType::already_projected;

        // Still compute geographic extents for diagnostics
        unproject_corners_to_geographic(meta, sug);

        sug.recommended_name = "Already projected — no reprojection needed";
        return sug;
    }

    // ------------------------------------------------------------------
    // 3. Unproject corners to geographic coordinates
    // ------------------------------------------------------------------
    unproject_corners_to_geographic(meta, sug);

    // ------------------------------------------------------------------
    // 4. Polar test:  |centre_lat| > 75°
    // ------------------------------------------------------------------
    constexpr double k_polar_threshold = 75.0;

    if (std::abs(sug.centre_lat) > k_polar_threshold)
    {
        bool north = sug.centre_lat > 0.0;

        sug.type              = north
                                    ? CrsRecommendationType::polar_north
                                    : CrsRecommendationType::polar_south;
        sug.recommended_epsg  = north ? 3413 : 3976;
        sug.recommended_wkt   = make_polar_wkt(north);
        sug.recommended_name  = north
            ? "NSIDC North Polar Stereographic (EPSG:3413)"
            : "NSIDC South Polar Stereographic (EPSG:3976)";

        return sug;
    }

    // ------------------------------------------------------------------
    // 5. UTM selection
    //    Zone from centre longitude; hemisphere from centre latitude.
    // ------------------------------------------------------------------
    int  zone  = lon_to_utm_zone(sug.centre_lon);
    bool north = sug.centre_lat >= 0.0;
    int  epsg  = utm_epsg(zone, north);

    sug.type              = CrsRecommendationType::utm;
    sug.recommended_epsg  = epsg;
    sug.recommended_wkt   = make_utm_wkt(zone, north);

    std::ostringstream name;
    name << "UTM Zone " << zone << (north ? "N" : "S")
         << " (EPSG:" << epsg << ")";
    sug.recommended_name = name.str();

    // ------------------------------------------------------------------
    // 6. Multi-zone warning
    // ------------------------------------------------------------------
    count_utm_zones(sug);

    if (sug.utm_zone_count >= 3)
    {
        sug.warn_multi_zone = true;

        std::ostringstream warn;
        warn << "WARNING: This raster spans " << sug.utm_zone_count
             << " UTM zones (zones " << sug.utm_zone_min
             << "–" << sug.utm_zone_max << "). "
             << "Projecting into UTM Zone " << zone
             << " will introduce increasing distortion toward the raster edges. "
             << "Consider using a Lambert Azimuthal Equal Area projection "
             << "centred on the raster instead (see laea_wkt).";
        sug.warning_message = warn.str();

        if (offer_laea)
        {
            sug.laea_wkt = make_laea_wkt(sug.centre_lon, sug.centre_lat);

            std::ostringstream laea_name;
            laea_name << std::fixed << std::setprecision(2)
                      << "Custom LAEA centred at ("
                      << sug.centre_lon << "°, "
                      << sug.centre_lat << "°)";
            sug.laea_name = laea_name.str();
        }
    }

    return sug;
}

// -----------------------------------------------------------------------------
// lsdtt_print_crs_suggestion
// -----------------------------------------------------------------------------
void lsdtt_print_crs_suggestion(const lsdtt_crs_suggestion& sug)
{
    std::cout << "\n=== lsdtt_suggest_crs ===\n";

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Centre       : " << sug.centre_lat << "°N  "
              << sug.centre_lon << "°E\n";
    std::cout << "  Extent lat   : " << sug.min_lat << " → " << sug.max_lat << "\n";
    std::cout << "  Extent lon   : " << sug.min_lon << " → " << sug.max_lon << "\n";

    if (sug.already_projected)
    {
        std::cout << "  Status       : already in a projected CRS — "
                     "no reprojection needed\n";
        std::cout << "=========================\n\n";
        return;
    }

    std::cout << "  Recommended  : " << sug.recommended_name << "\n";

    if (sug.type == CrsRecommendationType::utm)
    {
        std::cout << "  UTM zones    : " << sug.utm_zone_min
                  << " → " << sug.utm_zone_max
                  << " (" << sug.utm_zone_count << " zones spanned)\n";
    }

    if (sug.warn_multi_zone)
    {
        std::cout << "\n  [!] " << sug.warning_message << "\n";
        if (!sug.laea_wkt.empty())
        {
            std::cout << "  LAEA alt     : " << sug.laea_name << "\n";
        }
    }

    std::cout << "=========================\n\n";
}

} // namespace lsdtt
