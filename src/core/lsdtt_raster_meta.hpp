// lsdtt_raster_meta.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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
#include <gdal.h>
#include <ogr_spatialref.h>

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
    // GDAL data type of the source band. Compute objects work internally in
    // float32 regardless of the source type; this field is retained so that
    // outputs can be written in the same type as the input if desired.
    // -------------------------------------------------------------------------

    GDALDataType data_type = GDT_Float32;

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
    // These methods parse the WKT string using OGR to determine whether the
    // coordinate system is geographic (lat/lon) or projected, and what the
    // linear units are.
    //
    // They are used by analysis functions to enforce the requirement that all
    // lsdtt3 analyses operate on projected data in metres.
    // -------------------------------------------------------------------------

    // Sets crs_wkt and epsg_code from an EPSG code (e.g. 32601 for UTM zone 1
    // North on WGS84). Used to give synthetically generated rasters a sensible
    // default projected CRS. Throws std::runtime_error if the EPSG code cannot
    // be resolved. Requires OGR - include <ogr_spatialref.h> before using.
    void set_crs_from_epsg(int epsg)
    {
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

    // Gives a synthetically generated raster (no input DEM) a sensible default
    // georeference: the grid is centred on Teanu Island, Solomon Islands
    // (11.633 S, 166.967 E) - the real-world location to which Zwift's fictional
    // "Watopia" is geocoded for Strava. This sits comfortably inside UTM zone
    // 58S (EPSG:32758), well away from the antimeridian and the poles, so
    // generated DEMs carry realistic projected coordinates and round-trip
    // cleanly to and from lat/lon. n_rows, n_cols and cell_size must already be
    // set. Requires OGR. Falls back to no CRS if the transform fails.
    void set_default_synthetic_georeference()
    {
        constexpr double watopia_lon = 166.967;  // Teanu Island
        constexpr double watopia_lat = -11.633;
        constexpr int    epsg        = 32758;    // WGS 84 / UTM zone 58S

        try
        {
            set_crs_from_epsg(epsg);
        }
        catch (const std::exception&)
        {
            crs_wkt = "";
            return;
        }

        OGRSpatialReference wgs84;
        OGRSpatialReference utm;
        if (wgs84.importFromEPSG(4326) != OGRERR_NONE ||
            utm.importFromEPSG(epsg)   != OGRERR_NONE)
        {
            return;
        }
        wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        utm.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        OGRCoordinateTransformation* ct =
            OGRCreateCoordinateTransformation(&wgs84, &utm);
        if (!ct) return;

        double e = watopia_lon;
        double n = watopia_lat;
        const bool ok = ct->Transform(1, &e, &n) != 0;
        OGRCoordinateTransformation::DestroyCT(ct);
        if (!ok) return;

        // Centre the grid on the anchor point.
        x_min = e - 0.5 * static_cast<double>(n_cols) * cell_size;
        y_min = n - 0.5 * static_cast<double>(n_rows) * cell_size;
    }

    // Returns true if the CRS is geographic (latitude/longitude).
    // Returns false if projected or if no CRS is defined.
    // Requires OGR - include <ogr_spatialref.h> before using.
    bool is_geographic() const
    {
        if (crs_wkt.empty()) return false;
        OGRSpatialReference srs;
        if (srs.importFromWkt(crs_wkt.c_str()) != OGRERR_NONE) return false;
        return srs.IsGeographic() != 0;
    }

    // Returns the name of the linear unit of the CRS (e.g. "metre", "foot",
    // "US survey foot"). Returns empty string if no CRS or if geographic.
    std::string linear_unit_name() const
    {
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
    double linear_unit_to_metres() const
    {
        if (crs_wkt.empty()) return 0.0;
        OGRSpatialReference srs;
        if (srs.importFromWkt(crs_wkt.c_str()) != OGRERR_NONE) return 0.0;
        if (srs.IsGeographic()) return 0.0;
        double to_metres = srs.GetLinearUnits(nullptr);
        return to_metres;
    }

    // Returns true if the linear unit is metres (within floating point tolerance).
    bool is_in_metres() const
    {
        double factor = linear_unit_to_metres();
        return std::abs(factor - 1.0) < 1e-6;
    }

    // Throws std::runtime_error with a clear message if the CRS is geographic
    // or if the linear unit is not metres.
    // Call this at the start of any analysis function that requires projected
    // data in metres.
    void check_projected_metres() const
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
};

} // namespace lsdtt
