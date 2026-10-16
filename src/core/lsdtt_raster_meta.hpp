// lsdtt_raster_meta.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_meta
//
// Plain struct carrying all geospatial metadata for a raster dataset.
// Populated exclusively by lsdtt_data_loader and passed read-only to all
// compute objects (lsdtt_raster_cpu, lsdtt_raster_gpu).
//
// Compute objects must never modify this struct. Only lsdtt_data_loader
// creates or updates lsdtt_raster_meta instances.
//
// The coordinate convention used throughout is:
//   x_min  - easting  of the left   edge of the leftmost  column
//   y_min  - northing of the bottom edge of the bottommost row
//   cell_size - pixel width and height (assumed square)
//
// This matches the convention used in LSDTopoTools1 and is consistent with
// the bottom-left origin used in most GIS software when displaying rasters.
// Note that GDAL's geotransform uses a top-left origin internally; the
// conversion is handled in lsdtt_data_loader.
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

#include "lsdtt_pixel_type.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// lsdtt_raster_meta
// -----------------------------------------------------------------------------
struct lsdtt_raster_meta
{
    // -------------------------------------------------------------------------
    // Spatial dimensions
    // -------------------------------------------------------------------------

    int n_rows = 0;   // number of rows    (y direction)
    int n_cols = 0;   // number of columns (x direction)

    // -------------------------------------------------------------------------
    // Georeferencing
    //
    // Origin is the bottom-left corner of the bottom-left pixel.
    // Cell size is assumed square (equal x and y resolution).
    // -------------------------------------------------------------------------

    double x_min      = 0.0;   // easting  of left   edge (metres or degrees)
    double y_min      = 0.0;   // northing of bottom edge (metres or degrees)
    double cell_size  = 1.0;   // pixel size in projection units

    // Derived extents - computed from the above, provided for convenience
    double x_max() const { return x_min + n_cols * cell_size; }
    double y_max() const { return y_min + n_rows * cell_size; }

    // -------------------------------------------------------------------------
    // NoData
    // -------------------------------------------------------------------------

    double no_data_value = -9999.0;

    // -------------------------------------------------------------------------
    // Coordinate reference system
    //
    // Stored as a full OGC WKT string. Empty string means no CRS defined.
    // The loader also retains the EPSG code where available, but WKT is the
    // authoritative representation used for writing outputs.
    // -------------------------------------------------------------------------

    std::string crs_wkt;       // full WKT projection string
    int         epsg_code = 0; // EPSG code if known, 0 if unknown

    // -------------------------------------------------------------------------
    // Data type
    //
    // Sample type of the source band. Compute objects work internally in
    // float32 regardless of the source type; this field is retained so that
    // outputs can be written in the same type as the input if desired.
    //
    // lsdtt_pixel_type rather than GDALDataType so that this header -- which
    // every compute object includes -- needs no GDAL header at all.
    // -------------------------------------------------------------------------

    lsdtt_pixel_type data_type = lsdtt_pixel_type::float32;

    // -------------------------------------------------------------------------
    // Utility
    // -------------------------------------------------------------------------

    // Returns true if the metadata describes a valid (non-empty) raster.
    bool is_valid() const
    {
        return n_rows > 0 && n_cols > 0 && cell_size > 0.0;
    }

    // Returns true if the value v should be treated as nodata.
    // Handles both sentinel values (e.g. -9999) and NaN, since some DEMs
    // use NaN as their nodata representation. NaN == NaN is always false
    // in IEEE 754, so it must be tested explicitly with std::isnan.
    bool is_nodata(double v) const
    {
        if (std::isnan(no_data_value))
        {
            return std::isnan(v);
        }
        return v == no_data_value;
    }

    // Returns the total number of pixels.
    long long n_pixels() const
    {
        return static_cast<long long>(n_rows) * n_cols;
    }

    // Returns the ground area of one pixel in projection units squared.
    double pixel_area() const
    {
        return cell_size * cell_size;
    }

    // -------------------------------------------------------------------------
    // CRS validation
    //
    // These methods parse the WKT string to determine whether the coordinate
    // system is geographic (lat/lon) or projected, and what the linear units
    // are. They are DECLARED here and defined in lsdtt_raster_meta.cpp, so that
    // the OGR headers they need stay out of every translation unit that merely
    // carries an lsdtt_raster_meta around.
    //
    // They are used by analysis functions to enforce the requirement that all
    // lsdtt3 analyses operate on projected data in metres.
    // -------------------------------------------------------------------------

    // Sets crs_wkt and epsg_code from an EPSG code (e.g. 32601 for UTM zone 1
    // North on WGS84). Used to give synthetically generated rasters a sensible
    // default projected CRS. Throws std::runtime_error if the EPSG code cannot
    // be resolved.
    void set_crs_from_epsg(int epsg);

    // Places a synthetic grid at an arbitrary geographic anchor. The output
    // grid is centred on (centre_lat_deg, centre_lon_deg) in the selected
    // projected CRS. If epsg == 0, a WGS84 UTM EPSG is derived from
    // (centre_lon_deg, centre_lat_deg).
    //
    // Throws std::runtime_error when the anchor is non-finite or when
    // epsg == 0 and |centre_lat_deg| > 84 (UTM undefined at polar latitudes).
    // n_rows, n_cols and cell_size must already be set.
    void set_synthetic_georeference(double centre_lat_deg,
                                    double centre_lon_deg,
                                    int    epsg = 0);

    // Gives a synthetically generated raster (no input DEM) a sensible default
    // georeference: the grid is centred on Teanu Island, Solomon Islands
    // (11.633 S, 166.967 E) - the real-world location to which Zwift's fictional
    // "Watopia" is geocoded for Strava. This sits comfortably inside UTM zone
    // 58S (EPSG:32758), well away from the antimeridian and the poles, so
    // generated DEMs carry realistic projected coordinates and round-trip
    // cleanly to and from lat/lon. n_rows, n_cols and cell_size must already be
    // set. Falls back to no CRS if the transform fails.
    void set_default_synthetic_georeference();

    // Returns true if the CRS is geographic (latitude/longitude).
    // Returns false if projected or if no CRS is defined.
    bool is_geographic() const;

    // Returns the name of the linear unit of the CRS (e.g. "metre", "foot",
    // "US survey foot"). Returns empty string if no CRS or if geographic.
    std::string linear_unit_name() const;

    // Returns the conversion factor from the CRS linear unit to metres.
    // Returns 1.0 if the unit is metres, ~0.3048 for feet, etc.
    // Returns 0.0 if the CRS is geographic or unknown.
    double linear_unit_to_metres() const;

    // Returns true if the linear unit is metres (within floating point tolerance).
    bool is_in_metres() const;

    // Throws std::runtime_error with a clear message if the CRS is geographic
    // or if the linear unit is not metres.
    // Call this at the start of any analysis function that requires projected
    // data in metres.
    void check_projected_metres() const;

    // How far the projection is from preserving distance on this grid: the
    // largest |h - 1| or |k - 1| over the grid's centre and four corners, where
    // h and k are the point scale factors along the grid's y and x axes
    // (projected distance of a 10-cell step divided by its ellipsoidal
    // geodesic length, Vincenty inverse on the CRS's own ellipsoid). 0.0004 for
    // UTM on its central meridian; about 0.49 for Web Mercator at 48 degrees.
    // Returns -1 if there is no CRS or it is geographic. Throws
    // std::runtime_error if the CRS cannot be transformed to its geographic CRS.
    [[nodiscard]] double max_linear_scale_error() const;

    // The CRS name from its WKT (e.g. "WGS 84 / UTM zone 32N"), or "none" when
    // there is no CRS, or "unparsed" when the WKT cannot be read.
    [[nodiscard]] std::string crs_name() const;
};

} // namespace lsdtt
