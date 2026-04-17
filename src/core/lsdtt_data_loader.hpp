// lsdtt_data_loader.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_data_loader
//
// The sole GDAL/OGR interface in lsdtt3. Responsible for all geospatial
// file I/O — reading and writing rasters, vector layers, and tabular data.
//
// Compute objects (lsdtt_raster_cpu, lsdtt_raster_gpu) never touch GDAL directly.
// They receive raw data arrays and lsdtt_raster_meta from this class, and return
// result arrays to it for writing.
//
// Supported input formats: any format supported by the installed GDAL build.
// Tested formats: GeoTIFF, Cloud Optimised GeoTIFF (COG), VRT, HDF5, NetCDF, ENVI bil.
//
// Supported output formats:
//   Raster  — Cloud Optimised GeoTIFF (.tif)
//   Vector  — GeoPackage (.gpkg)
//   Tabular — CSV (.csv, geopandas compatible)
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <variant>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_raster_meta.hpp"
#include "lsdtt_raster_result.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// lsdtt_tile_window
//
// Describes a single tile within a raster, including its halo (overlap region).
// The halo is read from the file but should not be written to the output — it
// exists solely to eliminate edge artefacts in compute operations.
//
//   x_off, y_off   — pixel offset of the tile origin (without halo) in the
//                    full raster
//   tile_width,
//   tile_height    — dimensions of the valid (non-halo) tile region in pixels
//   halo           — number of pixels of overlap on each edge
//
// The actual array read from disk is:
//   width  = tile_width  + 2 * halo  (clamped at raster edges)
//   height = tile_height + 2 * halo
// -----------------------------------------------------------------------------
struct lsdtt_tile_window {
    int x_off      = 0;
    int y_off      = 0;
    int tile_width  = 0;
    int tile_height = 0;
    int halo        = 0;

    // Actual read dimensions (halo included, clamped at edges)
    int read_x_off    = 0;
    int read_y_off    = 0;
    int read_width    = 0;
    int read_height   = 0;
};

// -----------------------------------------------------------------------------
// lsdtt_raster_tile
//
// A single tile of raster data returned by lsdtt_data_loader::read_tile().
// Owns its data and mask arrays.
//
//   data   — float32 elevation (or other) values, row-major, size = read_width
//            * read_height. NoData pixels are set to meta.no_data_value.
//   mask   — uint8 validity mask, same dimensions. 1 = valid, 0 = nodata.
//   window — the tile geometry used to read this tile
//   meta   — full raster metadata (shared for all tiles from the same file)
// -----------------------------------------------------------------------------
struct lsdtt_raster_tile {
    std::vector<float>    data;
    std::vector<uint8_t>  mask;
    lsdtt_tile_window       window;
    lsdtt_raster_meta       meta;
};


// -----------------------------------------------------------------------------
// vector_output_format
// point_output_format
// -----------------------------------------------------------------------------
enum class vector_output_format
{
    flatgeobuf, 
    gpkg        
};

enum class point_output_format
{ 
    csv,
    flatgeobuf,
    gpkg
};

// This takes a format string and returns the vector output format
vector_output_format resolve_vector_output_format(const std::string& fmt_str);

// This turns a point output format (not csv!) into a vector output format
vector_output_format convert_point_to_vector_output_format(const point_output_format& format);

// This takes a format string and returns the vector output format
point_output_format resolve_point_output_format(const std::string& fmt_str);

// These overloaded functions get the output extension for a 
// vector output or point output format 
inline std::string output_extension(const vector_output_format fmt)
{
    switch (fmt)
    {
        case vector_output_format::gpkg:       return ".gpkg";
        case vector_output_format::flatgeobuf: return ".fgb";
        default:
            throw std::invalid_argument("Unknown vector_output_format");
    }
}

inline std::string output_extension(const point_output_format fmt)
{
    switch (fmt)
    {
        case point_output_format::csv:        return ".csv";
        case point_output_format::gpkg:       return ".gpkg";
        case point_output_format::flatgeobuf: return ".fgb";
        default:
            throw std::invalid_argument("Unknown point_output_format");
    }
}
// =============================================================================
// Geometry value types
//
// lsdtt_xy is defined in lsdtt_raster_geometry.hpp already.
// =============================================================================
struct lsdtt_point_geom
{
    double x = 0.0;
    double y = 0.0;
};
 
struct lsdtt_line_geom
{
    std::vector<lsdtt_xy> points;   // >= 2 elements required when writing
};
 
struct lsdtt_polygon_geom
{
    std::vector<lsdtt_xy>              exterior;  // outer ring
    std::vector<std::vector<lsdtt_xy>> holes;     // inner rings (may be empty)
};
 
// std::monostate = null / no geometry
using lsdtt_geom = std::variant<
    std::monostate,
    lsdtt_point_geom,
    lsdtt_line_geom,
    lsdtt_polygon_geom
>;
 
// =============================================================================
// Vector attribute types
// =============================================================================
enum class lsdtt_vector_field_type
{
    int64_type,
    double_type,
    string_type
};
 
// std::monostate = SQL NULL
using lsdtt_vector_field_value = std::variant<
    std::monostate,
    std::int64_t,
    double,
    std::string
>;
 
// =============================================================================
// lsdtt_vector_schema
//
// Field names and types for a layer, stored once — not per feature.
// =============================================================================
struct lsdtt_vector_schema
{
    std::vector<std::string>             field_names;
    std::vector<lsdtt_vector_field_type> field_types;
 
    std::size_t size()  const noexcept { return field_names.size(); }
    bool        empty() const noexcept { return field_names.empty(); }
 
    // Returns the index of a named field, or -1 if not found.
    int index_of(const std::string& name) const noexcept
    {
        for (std::size_t i = 0; i < field_names.size(); ++i)
            if (field_names[i] == name) return static_cast<int>(i);
        return -1;
    }
};
 
// =============================================================================
// lsdtt_vector_record
//
// One feature: geometry + positional attribute values matching a schema.
// =============================================================================
struct lsdtt_vector_record
{
    lsdtt_geom                            geometry;
    std::vector<lsdtt_vector_field_value> field_values;
};
 
// =============================================================================
// lsdtt_vector_layer
//
// A complete in-memory layer: schema stored once, then records.
// Used as the return type for read_vector_layer() and as the buffer
// for spatial-index writes.
// =============================================================================
struct lsdtt_vector_layer
{
    lsdtt_vector_schema              schema;
    std::vector<lsdtt_vector_record> records;
 
    std::size_t size()  const noexcept { return records.size(); }
    bool        empty() const noexcept { return records.empty(); }
 
    // Typed geometry accessors — throw std::bad_variant_access if wrong type.
    static lsdtt_point_geom   as_point(const lsdtt_vector_record& r)
    {
        return std::get<lsdtt_point_geom>(r.geometry);
    }
    static lsdtt_line_geom    as_line(const lsdtt_vector_record& r)
    {
        return std::get<lsdtt_line_geom>(r.geometry);
    }
    static lsdtt_polygon_geom as_polygon(const lsdtt_vector_record& r)
    {
        return std::get<lsdtt_polygon_geom>(r.geometry);
    }
};

// -----------------------------------------------------------------------------
// raster_output_format
//
// Controls how output rasters are written.
//
//   cog            — Cloud Optimised GeoTIFF with DEFLATE compression and
//                    overview pyramid. Best for large rasters and cloud
//                    storage. Automatically falls back to gtiff for small
//                    rasters (below k_cog_min_pixels threshold).
//
//   gtiff          — Plain compressed GeoTIFF (DEFLATE, tiled). No overview
//                    pyramid. Faster to write, slightly larger than COG for
//                    large rasters. Suitable for intermediate results or
//                    small datasets.
//
//   same_as_input  — Uses the same GDAL driver as the input raster. Lets
//                    the output format match whatever format the user provided
//                    (e.g. ENVI bil, NetCDF, plain GeoTIFF). The file
//                    extension in out_path must match the driver.
//
// Threshold for automatic COG → GeoTIFF fallback:
//   COG overview pyramid construction has a fixed overhead that outweighs
//   the benefit for small rasters. Below k_cog_min_pixels pixels the loader
//   automatically writes a plain GeoTIFF even when cog is requested.
//   Benchmark crossover is approximately 512×512 = 262,144 pixels.
// -----------------------------------------------------------------------------
enum class raster_output_format
{
    cog,            // Cloud Optimised GeoTIFF (default; falls back to gtiff below threshold)
    gtiff,          // Plain compressed GeoTIFF
    same_as_input   // Match the input raster driver
};

// This resolves the raster output format from a string 
raster_output_format resolve_raster_output_format(const std::string& fmt_str);


// Minimum pixel count for COG overhead to be worthwhile.
// Below this threshold, write_full() uses plain GeoTIFF even when cog is requested.
static constexpr long long k_cog_min_pixels = 512LL * 512LL;  // 262,144 pixels

// -----------------------------------------------------------------------------
// lsdtt_data_loader
//
// Main I/O class. One instance per input raster or vector dataset.
//
// Typical raster usage:
//
//   lsdtt_data_loader loader;
//   loader.open_raster("dem.tif");
//
//   auto tiles = loader.make_tile_windows(512, 16);
//   for (auto& win : tiles)
//   {
//       auto tile = loader.read_tile(win);
//       // pass tile.data, tile.mask, tile.meta to a compute engine
//       // receive result_data back from compute engine
//       loader.write_tile(result_data, win, "slope.tif");
//   }
//   loader.finalise_cog("slope.tif");
//
// -----------------------------------------------------------------------------
class lsdtt_data_loader {
public:

    // -------------------------------------------------------------------------
    // Construction / destruction
    // -------------------------------------------------------------------------

    lsdtt_data_loader();
    ~lsdtt_data_loader();

    // Non-copyable — owns GDAL handles
    lsdtt_data_loader(const lsdtt_data_loader&)            = delete;
    lsdtt_data_loader& operator=(const lsdtt_data_loader&) = delete;

    // Movable
    lsdtt_data_loader(lsdtt_data_loader&&)            noexcept;
    lsdtt_data_loader& operator=(lsdtt_data_loader&&) noexcept;

    // -------------------------------------------------------------------------
    // Raster input
    // -------------------------------------------------------------------------

    // Open any GDAL-supported raster. Reads metadata; does not load pixel data.
    // Throws std::runtime_error if the file cannot be opened.
    void open_raster(const std::string& filepath);

    // Returns true if a raster is currently open.
    bool is_open() const;

    // Close the open raster and release all GDAL resources.
    void close_raster();

    // Return the metadata for the currently open raster.
    const lsdtt_raster_meta& get_meta() const;

    // -------------------------------------------------------------------------
    // Tile management
    // -------------------------------------------------------------------------

    // Partition the full raster into a grid of tiles.
    //
    //   tile_size  — target tile width and height in pixels (default 512)
    //   halo       — overlap in pixels added to each edge (default 16)
    //
    // Returns a vector of lsdtt_tile_window, one per tile, in row-major order.
    // Tiles at the raster edges are automatically clamped to the raster extent.
    std::vector<lsdtt_tile_window> make_tile_windows(
        int tile_size = 512,
        int halo      = 16
    ) const;

    // Read one tile from the open raster into an lsdtt_raster_tile.
    // The returned data array is float32, row-major, size = read_width *
    // read_height. NoData pixels appear at meta.no_data_value in data[] and
    // as 0 in mask[].
    lsdtt_raster_tile read_tile(const lsdtt_tile_window& window) const;

    // Read the entire raster into a single lsdtt_raster_tile (no tiling).
    // Use only when the full raster is known to fit in RAM.
    lsdtt_raster_tile read_full() const;

    // -------------------------------------------------------------------------
    // Raster output
    // -------------------------------------------------------------------------

    // Initialise an output raster for tiled writing.
    //
    //   out_path   — output file path
    //   meta       — metadata to embed (typically from get_meta())
    //   format     — output format (default cog; see raster_output_format)
    //   n_bands    — number of bands to create (default 1)
    //   out_type   — GDAL data type for the output band (default GDT_Float32).
    //                Pass meta.data_type only when you explicitly want to match
    //                the source type. Analysis results are float32 regardless
    //                of the source type; only use a different type when the
    //                data is genuinely of that type (e.g. GDT_UInt32 for
    //                integer labels produced by label_binary_components_by_size).
    //
    // For cog and gtiff, out_path should end in .tif.
    // For same_as_input, out_path extension must match the input driver.
    void create_output_raster(
        const std::string&       out_path,
        const lsdtt_raster_meta& meta,
        raster_output_format     format   = raster_output_format::cog,
        int                      n_bands  = 1,
        GDALDataType             out_type = GDT_Float32
    );

    // Write a result tile to a previously created output raster.
    // data must be float32, row-major, size = window.tile_width *
    // window.tile_height (the valid region only — do not include the halo).
    // GDAL converts the float32 values to the band's data type at write time.
    //
    //   data       — result values for the valid tile region
    //   window     — the tile window (used for pixel offsets)
    //   out_path   — must match a path previously passed to create_output_raster
    //   band       — 1-based band index (default 1)
    void write_tile(
        const std::vector<float>& data,
        const lsdtt_tile_window&  window,
        const std::string&        out_path,
        int                       band = 1
    );

    // Write a full raster in one call.
    // Convenience wrapper around create_output_raster + write_tile + finalise.
    //
    // For cog format: automatically falls back to gtiff for rasters smaller
    // than k_cog_min_pixels pixels (COG overhead not worthwhile below this).
    //
    //   out_type   — GDAL data type for the output band (default GDT_Float32).
    //                Use meta.data_type when the raster carries an explicit
    //                output type (e.g. the UInt32 label raster from
    //                label_binary_components_by_size). For all standard analysis
    //                outputs leave this as the default.
    void write_full(
        const std::vector<float>& data,
        const lsdtt_raster_meta&  meta,
        const std::string&        out_path,
        raster_output_format      format   = raster_output_format::cog,
        GDALDataType              out_type = GDT_Float32
    );

    // Convenience overload for raster results returned by core-layer engines.
    // Uses meta.data_type as the output type so that typed results (e.g. the
    // UInt32 label raster) are written correctly without an extra argument.
    void write_full(
        const lsdtt_raster_result& result,
        const std::string&         out_path,
        raster_output_format       format = raster_output_format::cog
    );

    // Finalise a tiled output file.
    // For cog format: adds overviews and re-orders blocks for COG compliance.
    // For gtiff and same_as_input: flushes and closes the dataset.
    // Call after all write_tile() calls for a given output file are complete.
    void finalise_output(const std::string& out_path);

    // -------------------------------------------------------------------------
    // Vector input (OGR)
    // -------------------------------------------------------------------------

    // Open a vector dataset (GeoPackage, Shapefile, etc.).
    // Throws std::runtime_error if the file cannot be opened.
    void open_vector(const std::string& filepath);

    // Return the names of all layers in the open vector dataset.
    std::vector<std::string> get_layer_names() const;

    // Read a vector layer. This is loaded in as an lsdtt_vector_layer,
    // which contains the schema (attribute table names and field types) and
    // the indiviudal features (including geometry and attribute values).
    [[nodiscard]] lsdtt_vector_layer read_vector_layer(
        const std::string& layer_name
    ) const;
 
    // Close the open vector dataset.
    void close_vector();

    // -----------------------------------------------------------------------
    // Vector output — session interface
    //
    // 1. Call begin_vector_write() to open the file and declare the schema.
    // 2. Call write_point / write_line / write_polygon for each feature.
    // 3. Call finish_vector_write() to flush and close.
    //
    // spatial_index = false  →  streaming, no index (low memory)
    // spatial_index = true   →  buffer all features, then write with GDAL
    //                           spatial index (higher memory, but will the 
    //                           index speeds up spatial queries when post-
    //                           processing in a GIS system.
    //
    // LSDTopoTools2 (the legacy codebase) printed data to csv and geojson
    //  geojson has been retired as it was only implemented to be GIS readable
    //  but has very large file sizes. csv output is still an option here
    //  for backward compatibility with lsdviztools. Although we will favour binary file 
    //  formats in the future (because of storage)
    //
    // The two binary formats available are GeoPackage (.gpkg) and FlatGeobuf. 
    //  The formats are passed on to the vector write through the format parameter.
    //  GeoPackage are backed by a full sqlite database and are a richer data format
    //  but they are slower to load and sometimes have write issues (for example, 
    //  this will fail if you are running in WSL2 and trying to write
    //  to /mnt/c)
    //  For this reason the FlatGeobuf is preferred. The FlatGeobuf can only hold one
    //  layer, you cannot combine point and polylines in the same file, for example. 
    // -------------------------------------------------------------------------
    void begin_vector_write(
        const std::string&          out_path,
        const std::string&          layer_name,
        const lsdtt_vector_schema&  schema,
        const std::string&          crs_wkt,
        vector_output_format  format,
        bool                        spatial_index = false
    );
 
    void write_point(
        double x, double y,
        const std::vector<lsdtt_vector_field_value>& field_values
    );
 
    void write_line(
        const std::vector<lsdtt_xy>&                 points,
        const std::vector<lsdtt_vector_field_value>& field_values
    );
 
    void write_polygon(
        const std::vector<lsdtt_xy>&                 exterior,
        const std::vector<lsdtt_vector_field_value>& field_values
    );
 
    void write_polygon(
        const std::vector<lsdtt_xy>&                 exterior,
        const std::vector<std::vector<lsdtt_xy>>&    holes,
        const std::vector<lsdtt_vector_field_value>& field_values
    );
 
    // Flush and close the active write session. Must be called after all
    // write_* calls. For the buffered path this triggers the GDAL write that
    // includes spatial index.
    void finish_vector_write();

    // -------------------------------------------------------------------------
    // Tabular output (CSV)
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Vector clipping
    // -------------------------------------------------------------------------

    // Clip a vector file to the bounding box of a raster.
    //
    // The raster bbox is transformed into the vector's native CRS before the
    // spatial filter is applied, so the two datasets need not share a CRS.
    // All layers in the source file are processed and written to the output.
    //
    //   input_path   — path to any OGR-readable vector file
    //   output_path  — path for the output file (extension must match out_fmt)
    //   raster_meta  — metadata for the raster whose bbox is used as the clip
    //                  rectangle; supplies both the extent and the CRS WKT
    //   out_fmt      — output format (flatgeobuf or gpkg)
    static void clip_vector_to_raster_bbox(
        const std::string&       input_path,
        const std::string&       output_path,
        const lsdtt_raster_meta& raster_meta,
        vector_output_format     out_fmt
    );

    // Write a table to a CSV file. The first row is a header.
    //
    //   out_path   — output .csv file path
    //   headers    — column names
    //   rows       — data rows; each inner vector must match headers.size()
    void write_csv(
        const std::string&                        out_path,
        const std::vector<std::string>&           headers,
        const std::vector<std::vector<std::string>>& rows
    );

    // -------------------------------------------------------------------------
    // Utility
    // -------------------------------------------------------------------------

    // Print a summary of the open raster to stdout (dimensions, CRS, nodata).
    void print_raster_info() const;

    // Return the GDAL version string.
    static std::string gdal_version();

private:

    // RAII wrapper for GDALDataset* using a custom deleter
    struct gdal_dataset_deleter
    {
        void operator()(GDALDataset* ds) const
        {
            if (ds) GDALClose(ds);
        }
    };
    using gdal_dataset_ptr = std::unique_ptr<GDALDataset, gdal_dataset_deleter>;

    // Open raster dataset
    gdal_dataset_ptr  raster_dataset_;

    // Open vector dataset
    gdal_dataset_ptr  vector_dataset_;

    // Output raster datasets keyed by file path
    // (multiple outputs can be open simultaneously during tiled writing)
    std::unordered_map<std::string, gdal_dataset_ptr>       output_datasets_;

    // Output format per dataset — needed by finalise_output()
    std::unordered_map<std::string, raster_output_format>   output_formats_;

    // Driver name of the open input raster — used for same_as_input output
    std::string input_driver_name_;

    // Metadata for the currently open raster
    lsdtt_raster_meta   meta_;

    // Whether a raster is currently open
    bool              raster_open_ = false;

    // Whether a vector dataset is currently open
    bool              vector_open_ = false;

    // The state of a vector write sessoin
    struct vector_write_session;
    std::unique_ptr<vector_write_session> write_session_;

    // Internal helpers
    void              read_meta_from_dataset();
    lsdtt_tile_window   clamp_window_to_raster(const lsdtt_tile_window& win) const;
    void              add_cog_overviews(GDALDataset* ds);

    // Write a single lsdtt_vector_record to an already-open OGRLayer.
    // Used by both the streaming path (one at a time) and the buffered
    // path (bulk write in finish_vector_write).
    void write_record_to_layer(
        const lsdtt_vector_record& record,
        void*                      ogr_layer  // OGRLayer* — opaque to avoid
    ) const;                                  // GDAL headers leaking into .hpp
};

} // namespace lsdtt
