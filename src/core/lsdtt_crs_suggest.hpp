// lsdtt_crs_suggest.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_crs_suggest
//
// Inspects a raster (via its lsdtt_raster_meta) and recommends a projected
// coordinate system suitable for analysis.
// For the time being this works on rasters on Earth. You will need different
// logic for Mars or Titan or any other non-Earth coordinate system.
//
// Decision logic
// 
//   1. Raster must **either** be in a geographic CRS (EPSG:4326 or Web Mercator are the
//      common cases) **or** the user must call `force_reproject` 
//
//   2. The four corners of the raster are unprojected to geographic
//      coordinates. The geographic centre is taken as the mean lat/lon.
//
//   3. Polar test  |centre_lat| > 75 deg
//      -> North:  NSIDC North Polar Stereographic  (EPSG:3413)
//      -> South:  NSIDC South Polar Stereographic  (EPSG:3976)
//
//   4. UTM test  otherwise
//      -> UTM zone chosen from centre longitude; hemisphere from centre lat.
//      -> If the raster spans >= 3 UTM zones a warning is emitted and the
//        caller may request a custom Lambert Azimuthal Equal Area (LAEA)
//        centred on the raster centre instead.
//
// Usage
// 
//   lsdtt_crs_suggestion sug = lsdtt_suggest_crs(meta);
//
//   if (sug.already_projected)
//       // nothing to do
//
//   if (sug.warn_multi_zone)
//       // print sug.warning_message, offer LAEA option
//
//   // Reproject using sug.recommended_wkt  (or sug.laea_wkt if preferred)
//
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include "lsdtt_raster_meta.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// CrsRecommendationType - what kind of CRS was selected
//  in future we may include a "force_projection" so you can change from
//  on projection to another and update the gird spacing
// -----------------------------------------------------------------------------
enum class CrsRecommendationType
{
    already_projected,          // Input was already in a projected CRS
    polar_north,                // NSIDC North Polar Stereographic (EPSG:3413)
    polar_south,                // NSIDC South Polar Stereographic (EPSG:3976)
    utm,                        // UTM zone centred on the raster
    laea_custom                 // Custom Lambert Azimuthal Equal Area
};

// -----------------------------------------------------------------------------
// lsdtt_crs_suggestion
//
// Returned by lsdtt_suggest_crs(). Contains the recommended WKT, diagnostic
// information, and an optional LAEA alternative for wide-extent rasters.
// -----------------------------------------------------------------------------
struct lsdtt_crs_suggestion
{
    // -------------------------------------------------------------------------
    // Primary recommendation
    // -------------------------------------------------------------------------

    CrsRecommendationType type = CrsRecommendationType::utm;

    // WKT string for the recommended CRS. Pass this to your reprojection
    // pipeline. Empty only when already_projected == true.
    std::string recommended_wkt;

    // Human-readable name for logging / UI display.
    // e.g. "UTM Zone 30N (EPSG:32630)"
    std::string recommended_name;

    // EPSG code for the recommended CRS (0 if custom / not applicable).
    int recommended_epsg = 0;

    // -------------------------------------------------------------------------
    // Geographic diagnostics
    //
    // All values are in decimal degrees.
    // Populated even if the raster is already projected, using whatever
    // geographic coordinates can be inferred.
    // -------------------------------------------------------------------------

    double centre_lon  = 0.0;   // geographic centre longitude
    double centre_lat  = 0.0;   // geographic centre latitude

    double min_lon     = 0.0;   // westernmost extent
    double max_lon     = 0.0;   // easternmost extent
    double min_lat     = 0.0;   // southernmost extent
    double max_lat     = 0.0;   // northernmost extent

    int    utm_zone_min = 0;    // westernmost UTM zone intersected
    int    utm_zone_max = 0;    // easternmost UTM zone intersected
    int    utm_zone_count = 0;  // number of distinct UTM zones spanned

    // -------------------------------------------------------------------------
    // Status flags
    // -------------------------------------------------------------------------

    // True if the input was already in a non-geographic projected CRS.
    // When true, recommended_wkt is empty and no reprojection is needed.
    bool already_projected = false;

    // True if the raster spans 3 or more UTM zones. Only set when
    // type == utm. When true, consider using laea_wkt instead.
    bool warn_multi_zone = false;

    // Human-readable warning / suggestion text for multi-zone rasters.
    std::string warning_message;

    // -------------------------------------------------------------------------
    // LAEA alternative (multi-zone rasters only)
    //
    // A custom Lambert Azimuthal Equal Area projection centred on
    // (centre_lon, centre_lat). Populated when warn_multi_zone == true.
    // Use this when you need area-preserving properties across a wide extent.
    // -------------------------------------------------------------------------

    std::string laea_wkt;   // WKT for the custom LAEA CRS
    std::string laea_name;  // e.g. "Custom LAEA centred at (lon, lat)"
};

// -----------------------------------------------------------------------------
// lsdtt_suggest_crs
//
// Inspect the raster described by `meta` and return a CRS recommendation.
//
// Parameters
//   meta              - raster metadata (must contain a valid crs_wkt)
//   offer_laea        - if true and the raster spans >= 3 UTM zones, also
//                       build the laea_wkt alternative (default true)
//   force_recommend_projected - if true, continue through the recommendation
//                               logic even when the input raster is already in
//                               a projected CRS. This is used when a caller
//                               wants to force reprojection into the preferred
//                               analysis CRS.
//
// Throws std::runtime_error if meta.crs_wkt is empty or cannot be parsed.
// -----------------------------------------------------------------------------
lsdtt_crs_suggestion lsdtt_suggest_crs(
    const lsdtt_raster_meta& meta,
    bool offer_laea = true,
    bool force_recommend_projected = false
);


// -----------------------------------------------------------------------------
// lsdtt_print_crs_suggestion
//
// Pretty-print a lsdtt_crs_suggestion to stdout. Useful for CLI tools.
// -----------------------------------------------------------------------------
void lsdtt_print_crs_suggestion(const lsdtt_crs_suggestion& sug);

} // namespace lsdtt
