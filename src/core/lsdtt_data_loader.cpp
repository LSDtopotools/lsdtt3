// lsdtt_data_loader.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_data_loader.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <type_traits>
#include <variant>



#include <gdal.h>
#include <gdal_priv.h>
#include <gdal_utils.h>
#include <cpl_conv.h>
#include <ogrsf_frmts.h>

namespace lsdtt {

// -----------------------------------------------------------------------------
// GDAL initialisation helper - called once via static local
// -----------------------------------------------------------------------------
static void ensure_gdal_init()
{
    static bool initialised = false;
    if (!initialised)
    {
        GDALAllRegister();
        initialised = true;
    }
}

// some helper functions
namespace
{

OGRFieldType to_ogr_field_type(lsdtt::lsdtt_vector_field_type type)
{
    switch (type)
    {
        case lsdtt::lsdtt_vector_field_type::int64_type:
            return OFTInteger64;

        case lsdtt::lsdtt_vector_field_type::double_type:
            return OFTReal;

        case lsdtt::lsdtt_vector_field_type::string_type:
        default:
            return OFTString;
    }
}

lsdtt::lsdtt_vector_field_type from_ogr_field_type(OGRFieldType type)
{
    switch (type)
    {
        case OFTInteger:
        case OFTInteger64:
            return lsdtt::lsdtt_vector_field_type::int64_type;

        case OFTReal:
            return lsdtt::lsdtt_vector_field_type::double_type;

        case OFTString:
        case OFTWideString:
        default:
            return lsdtt::lsdtt_vector_field_type::string_type;
    }
}

lsdtt::lsdtt_vector_field_value get_typed_field_value(
    const OGRFeature*                 feature,
    int                               field_index,
    lsdtt::lsdtt_vector_field_type    field_type
)
{
    if (!feature->IsFieldSetAndNotNull(field_index))
    {
        return std::monostate{};
    }

    switch (field_type)
    {
        case lsdtt::lsdtt_vector_field_type::int64_type:
            return static_cast<std::int64_t>(
                feature->GetFieldAsInteger64(field_index)
            );

        case lsdtt::lsdtt_vector_field_type::double_type:
            return feature->GetFieldAsDouble(field_index);

        case lsdtt::lsdtt_vector_field_type::string_type:
        default:
            return std::string(feature->GetFieldAsString(field_index));
    }
}

void set_typed_field_value(
    OGRFeature*                               feature,
    int                                       field_index,
    const lsdtt::lsdtt_vector_field_value&    value
)
{
    std::visit(
        [&](const auto& v)
        {
            using value_type = std::decay_t<decltype(v)>;

            if constexpr (std::is_same_v<value_type, std::monostate>)
            {
                feature->SetFieldNull(field_index);
            }
            else if constexpr (std::is_same_v<value_type, std::int64_t>)
            {
                feature->SetField(
                    field_index,
                    static_cast<GIntBig>(v)
                );
            }
            else if constexpr (std::is_same_v<value_type, double>)
            {
                feature->SetField(field_index, v);
            }
            else if constexpr (std::is_same_v<value_type, std::string>)
            {
                feature->SetField(field_index, v.c_str());
            }
        },
        value
    );
}

// Build an OGRGeometry* from an lsdtt_geom. Caller takes ownership via
// OGRFeature::SetGeometryDirectly().
OGRGeometry* build_ogr_geometry(const lsdtt_geom& geom)
{
    return std::visit(
        [](const auto& g) -> OGRGeometry*
        {
            using T = std::decay_t<decltype(g)>;
 
            if constexpr (std::is_same_v<T, std::monostate>)
            {
                return nullptr;
            }
            else if constexpr (std::is_same_v<T, lsdtt_point_geom>)
            {
                return new OGRPoint(g.x, g.y);
            }
            else if constexpr (std::is_same_v<T, lsdtt_line_geom>)
            {
                OGRLineString* ls = new OGRLineString();
                ls->setNumPoints(static_cast<int>(g.points.size()));
                for (int i = 0; i < static_cast<int>(g.points.size()); ++i)
                    ls->setPoint(i, g.points[i].x, g.points[i].y);
                return ls;
            }
            else if constexpr (std::is_same_v<T, lsdtt_polygon_geom>)
            {
                OGRPolygon*    poly = new OGRPolygon();
                OGRLinearRing* ring = new OGRLinearRing();
                for (const auto& pt : g.exterior) ring->addPoint(pt.x, pt.y);
                ring->closeRings();
                poly->addRingDirectly(ring);
                for (const auto& hole : g.holes)
                {
                    OGRLinearRing* hr = new OGRLinearRing();
                    for (const auto& pt : hole) hr->addPoint(pt.x, pt.y);
                    hr->closeRings();
                    poly->addRingDirectly(hr);
                }
                return poly;
            }
        },
        geom
    );
}
 
// Read a SINGLE-PART OGRGeometry into an lsdtt_geom. Returns monostate for
// anything else -- callers should go through read_ogr_geometry_parts(), which
// handles multipart and reports what it could not read.
lsdtt_geom read_ogr_geometry_single(OGRGeometry* ogr_geom)
{
    if (!ogr_geom) return std::monostate{};
 
    const OGRwkbGeometryType gt = wkbFlatten(ogr_geom->getGeometryType());
 
    if (gt == wkbPoint)
    {
        OGRPoint* pt = static_cast<OGRPoint*>(ogr_geom);
        return lsdtt_point_geom{ pt->getX(), pt->getY() };
    }
 
    if (gt == wkbLineString)
    {
        OGRLineString* ls = static_cast<OGRLineString*>(ogr_geom);
        lsdtt_line_geom lg;
        const int n = ls->getNumPoints();
        lg.points.reserve(n);
        for (int i = 0; i < n; ++i)
            lg.points.push_back({ ls->getX(i), ls->getY(i) });
        return lg;
    }
 
    if (gt == wkbPolygon)
    {
        OGRPolygon* poly = static_cast<OGRPolygon*>(ogr_geom);
        lsdtt_polygon_geom pg;
 
        auto read_ring = [](const OGRLinearRing* ring)
        {
            std::vector<lsdtt_xy> pts;
            const int n = ring->getNumPoints();
            pts.reserve(n);
            for (int i = 0; i < n; ++i)
                pts.push_back({ ring->getX(i), ring->getY(i) });
            return pts;
        };
 
        if (const OGRLinearRing* ext = poly->getExteriorRing())
            pg.exterior = read_ring(ext);
        for (int h = 0; h < poly->getNumInteriorRings(); ++h)
            pg.holes.push_back(read_ring(poly->getInteriorRing(h)));
 
        return pg;
    }
 
    return std::monostate{};  // unsupported geometry type
}

// -----------------------------------------------------------------------------
// read_ogr_geometry_parts
//
// Explode a geometry into its single parts.
//
// WHY EXPLODE RATHER THAN ADD A MULTIPART VARIANT. `lsdtt_line_geom` holds one
// vector of points, so a multipart line has to become either several records or
// a new variant type. Flattening the parts into ONE line is not an option: it
// invents a segment joining part N's end to part N+1's start, which is exactly
// the artefact that makes multipart necessary in the first place (a channel
// crossing a periodic boundary, an island, a road interrupted by a bridge).
//
// Exploding gives every existing consumer the right answer with no change: they
// all iterate records and handle single-part geometry. The cost is that the
// record count no longer matches the input feature count, which is why the
// caller stamps `source_feature_index` on each part.
//
// Until 2026-08-26 there was no multipart branch AT ALL, in any of the four
// types, and an unrecognised geometry returned monostate silently. Every line
// consumer then did `get_if<lsdtt_line_geom>(...)`, got nullptr, and skipped the
// feature -- so a MultiLineString input produced an empty output and a zero exit
// code. Demonstrated on lsdtt-linear-stats (1 feature vs empty layer) and on
// lsdtt-raster-preprocessing point sampling (3 features vs a 0-byte file).
// -----------------------------------------------------------------------------
std::vector<lsdtt_geom> read_ogr_geometry_parts(OGRGeometry* ogr_geom,
                                                std::string&  unsupported_name)
{
    std::vector<lsdtt_geom> out;
    if (!ogr_geom) return out;

    const OGRwkbGeometryType gt = wkbFlatten(ogr_geom->getGeometryType());

    if (gt == wkbMultiPoint || gt == wkbMultiLineString ||
        gt == wkbMultiPolygon || gt == wkbGeometryCollection)
    {
        OGRGeometryCollection* gc =
            static_cast<OGRGeometryCollection*>(ogr_geom);

        for (int i = 0; i < gc->getNumGeometries(); ++i)
        {
            // Recurse: a GeometryCollection may nest, and a nested collection
            // is not a part we can hand to a consumer.
            std::vector<lsdtt_geom> inner =
                read_ogr_geometry_parts(gc->getGeometryRef(i), unsupported_name);
            out.insert(out.end(), inner.begin(), inner.end());
        }
        return out;
    }

    lsdtt_geom single = read_ogr_geometry_single(ogr_geom);
    if (std::holds_alternative<std::monostate>(single))
    {
        // Report the type rather than dropping it in silence. A reader that
        // discards data without saying so is the defect this function exists to
        // fix; leaving the behaviour for the NEXT unsupported type repeats it.
        if (unsupported_name.empty())
        {
            const char* name = OGRGeometryTypeToName(ogr_geom->getGeometryType());
            unsupported_name = name ? name : "unknown";
        }
        return out;
    }

    out.push_back(std::move(single));
    return out;
}

} // unnamed namespace



// =============================================================================
// Resolving output formats 
// =============================================================================
vector_output_format resolve_vector_output_format(
    const std::string& fmt_str)
{
    if (fmt_str == "gpkg")
    {
        return vector_output_format::gpkg;
    }
    if (fmt_str == "flatgeobuf" || fmt_str == "fgb")
    {
        return vector_output_format::flatgeobuf;
    }
    std::cerr << "Warning: unknown vector_output_format "
              << fmt_str << "' - defaulting to flatgeobuf.\n";

    return vector_output_format::flatgeobuf;
}

point_output_format resolve_point_output_format(
    const std::string& fmt_str)
{
    if (fmt_str == "csv")
    {
        return point_output_format::csv;
    }
    if (fmt_str == "gpkg")
    {
        return point_output_format::gpkg;
    }
    if (fmt_str == "flatgeobuf" || fmt_str == "fgb")
    {
        return point_output_format::flatgeobuf;
    }
    std::cerr << "Warning: unknown vector_output_format '"
              << fmt_str << "' - defaulting to flatgeobuf.\n";

    return point_output_format::flatgeobuf;
}

vector_output_format convert_point_to_vector_output_format(
    const point_output_format& format)
{
   
    if (format == point_output_format::gpkg)
    {
        return vector_output_format::gpkg;
    }
    if (format == point_output_format::flatgeobuf)
    {
        return vector_output_format::flatgeobuf;
    }
    
    std::cerr << "Warning: cannot convert point_output_format "
                 "to vector_output_format - defaulting to flatgeobuf.\n";
    return vector_output_format::flatgeobuf;
}

raster_output_format resolve_raster_output_format(
    const std::string& fmt_str)
{
    if (fmt_str == "cog")          return raster_output_format::cog;
    if (fmt_str == "gtiff")        return raster_output_format::gtiff;
    if (fmt_str == "same_as_input") return raster_output_format::same_as_input;
    if (fmt_str == "png")          return raster_output_format::png;

    std::cerr << "Warning: unknown raster_output_format '" << fmt_str
              << "' - defaulting to cog.\n";
    return raster_output_format::cog;
}

// =============================================================================
// vector_write_session
//
// Holds all state for an open vector write session. Lives inside lsdtt_data_loader
// as a unique_ptr so GDAL types don't leak into the header.
// =============================================================================
struct lsdtt_data_loader::vector_write_session
{
    // Configuration (set at begin_vector_write)
    std::string                out_path;
    std::string                layer_name;
    lsdtt_vector_schema        schema;
    std::string                crs_wkt;
    vector_output_format format;
    bool                       spatial_index;
 
    // Streaming path: dataset and layer are open throughout the session.
    struct gdal_dataset_deleter_local
    {
        void operator()(GDALDataset* ds) const { if (ds) GDALClose(ds); }
    };
    using ds_ptr = std::unique_ptr<GDALDataset, gdal_dataset_deleter_local>;
 
    ds_ptr     dataset;
    OGRLayer*  layer = nullptr;
 
    // Buffered path: accumulate records here; write all in finish().
    lsdtt_vector_layer buffer;
 
    // -----------------------------------------------------------------------
    // Open the GDAL dataset and create the layer.
    // -----------------------------------------------------------------------
    void open_dataset_and_layer(OGRwkbGeometryType geom_type)
    {
        const char* driver_name =
            (format == vector_output_format::flatgeobuf)
                ? "FlatGeobuf"
                : "GPKG";
 
        GDALDriver* driver =
            GetGDALDriverManager()->GetDriverByName(driver_name);
        if (!driver)
            throw std::runtime_error(
                std::string("lsdtt_data_loader: driver not available: ") +
                driver_name
            );
 
        GDALDataset* ds = nullptr;
 
        if (format == vector_output_format::flatgeobuf)
        {
            // FlatGeobuf is a single-layer single-file format: always recreate.
            VSIUnlink(out_path.c_str());
            ds = driver->Create(out_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        }
        else
        {
            // GeoPackage: open existing file for update, or create new.
            VSIStatBufL stat;
            if (VSIStatL(out_path.c_str(), &stat) == 0)
            {
                ds = static_cast<GDALDataset*>(
                    GDALOpenEx(
                        out_path.c_str(),
                        GDAL_OF_VECTOR | GDAL_OF_UPDATE,
                        nullptr, nullptr, nullptr
                    )
                );
            }
            if (!ds)
                ds = driver->Create(out_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        }
 
        if (!ds)
            throw std::runtime_error(
                "lsdtt_data_loader: could not create/open file: " + out_path
            );
 
        dataset.reset(ds);
 
        // Replace an existing layer of the same name (GeoPackage only;
        // FlatGeobuf was already unlinked above).
        for (int i = 0; i < ds->GetLayerCount(); ++i)
        {
            OGRLayer* l = ds->GetLayer(i);
            if (l && layer_name == l->GetName()) { ds->DeleteLayer(i); break; }
        }
 
        OGRSpatialReference  srs;
        OGRSpatialReference* srs_ptr = nullptr;
        if (!crs_wkt.empty())
        {
            srs.importFromWkt(crs_wkt.c_str());
            srs_ptr = &srs;
        }
 
        // Pass SPATIAL_INDEX=NO when not requested.
        const char* opts_no_index[]   = { "SPATIAL_INDEX=NO", nullptr };
        const char* opts_with_index[] = { nullptr };
        char** layer_opts = const_cast<char**>(
            spatial_index ? opts_with_index : opts_no_index
        );
 
        layer = ds->CreateLayer(
            layer_name.c_str(), srs_ptr, geom_type, layer_opts
        );
        if (!layer)
            throw std::runtime_error(
                "lsdtt_data_loader: could not create layer: " + layer_name
            );
 
        // Create attribute fields from schema.
        for (std::size_t f = 0; f < schema.size(); ++f)
        {
            OGRFieldDefn field(
                schema.field_names[f].c_str(),
                to_ogr_field_type(schema.field_types[f])
            );
            if (layer->CreateField(&field) != OGRERR_NONE)
                throw std::runtime_error(
                    "lsdtt_data_loader: failed to create field: " +
                    schema.field_names[f]
                );
        }
 
        // For GeoPackage streaming writes, open a transaction now so that
        // all feature inserts are batched into a single SQLite commit.
        // (For the buffered path the transaction is opened in finish_vector_write
        // immediately before the bulk insert loop.)
        if (format == vector_output_format::gpkg && !spatial_index)
        {
            if (layer->StartTransaction() != OGRERR_NONE)
                throw std::runtime_error(
                    "lsdtt_data_loader: failed to start transaction on layer: " +
                    layer_name
                );
        }
    }
};

// -----------------------------------------------------------------------------
// Construction / destruction
// -----------------------------------------------------------------------------

lsdtt_data_loader::lsdtt_data_loader()
{
    ensure_gdal_init();
}

lsdtt_data_loader::~lsdtt_data_loader() = default;

lsdtt_data_loader::lsdtt_data_loader(lsdtt_data_loader&&) noexcept = default;
lsdtt_data_loader& lsdtt_data_loader::operator=(lsdtt_data_loader&&) noexcept = default;

// -----------------------------------------------------------------------------
// Raster input
// -----------------------------------------------------------------------------

void lsdtt_data_loader::open_raster(const std::string& filepath)
{
    if (raster_open_)
    {
        close_raster();
    }

    GDALDataset* raw = static_cast<GDALDataset*>(
        GDALOpen(filepath.c_str(), GA_ReadOnly)
    );

    if (!raw)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: could not open raster: " + filepath +
            "\nGDAL error: " + std::string(CPLGetLastErrorMsg())
        );
    }

    raster_dataset_.reset(raw);
    read_meta_from_dataset();

    // Record the driver name for same_as_input output
    input_driver_name_ = raw->GetDriver()
                         ? std::string(raw->GetDriver()->GetDescription())
                         : "GTiff";

    raster_open_ = true;
}

bool lsdtt_data_loader::is_open() const
{
    return raster_open_;
}

void lsdtt_data_loader::close_raster()
{
    raster_dataset_.reset();
    raster_open_ = false;
    meta_ = lsdtt_raster_meta{};
}

const lsdtt_raster_meta& lsdtt_data_loader::get_meta() const
{
    if (!raster_open_)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: get_meta() called with no raster open"
        );
    }
    return meta_;
}

// -----------------------------------------------------------------------------
// Tile management
// -----------------------------------------------------------------------------

std::vector<lsdtt_tile_window> lsdtt_data_loader::make_tile_windows(
    int tile_size,
    int halo
) const
{
    if (!raster_open_)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: make_tile_windows() called with no raster open"
        );
    }

    std::vector<lsdtt_tile_window> windows;

    int n_tiles_x = static_cast<int>(
        std::ceil(static_cast<double>(meta_.n_cols) / tile_size)
    );
    int n_tiles_y = static_cast<int>(
        std::ceil(static_cast<double>(meta_.n_rows) / tile_size)
    );

    windows.reserve(n_tiles_x * n_tiles_y);

    for (int ty = 0; ty < n_tiles_y; ++ty)
    {
        for (int tx = 0; tx < n_tiles_x; ++tx)
        {
            lsdtt_tile_window win;
            win.x_off       = tx * tile_size;
            win.y_off       = ty * tile_size;
            win.tile_width  = std::min(tile_size, meta_.n_cols - win.x_off);
            win.tile_height = std::min(tile_size, meta_.n_rows - win.y_off);
            win.halo        = halo;

            windows.push_back(clamp_window_to_raster(win));
        }
    }

    return windows;
}

lsdtt_raster_tile lsdtt_data_loader::read_tile(const lsdtt_tile_window& window) const
{
    if (!raster_open_)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: read_tile() called with no raster open"
        );
    }

    GDALRasterBand* band = raster_dataset_->GetRasterBand(1);
    if (!band)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: could not get raster band 1"
        );
    }

    int buf_size = window.read_width * window.read_height;

    lsdtt_raster_tile tile;
    tile.window = window;
    tile.meta   = meta_;
    tile.data.resize(buf_size);
    tile.mask.resize(buf_size, 1);  // default all valid

    CPLErr err = band->RasterIO(
        GF_Read,
        window.read_x_off,
        window.read_y_off,
        window.read_width,
        window.read_height,
        tile.data.data(),
        window.read_width,
        window.read_height,
        GDT_Float32,
        0, 0
    );

    if (err != CE_None)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: RasterIO read failed: " +
            std::string(CPLGetLastErrorMsg())
        );
    }

    // Build nodata mask.
    // Treat non-finite values (NaN / Inf) as nodata as well: some DEMs
    // (e.g. ENVI .bil exports) encode their no-data border as IEEE NaN while
    // the header still declares a sentinel such as -9999. If these are not
    // caught here they survive threshold masking and poison flow
    // preprocessing - a NaN rim around the raster reads as a dam with no
    // outlet, so the priority-flood fill pancakes the whole interior.
    // Normalise detected nodata pixels to the sentinel so the
    // mask==0 => data holds no_data_value convention is preserved downstream.
    float nodata = static_cast<float>(meta_.no_data_value);
    for (int i = 0; i < buf_size; ++i)
    {
        float v = tile.data[i];
        if (v == nodata || !std::isfinite(v))
        {
            tile.mask[i] = 0;
            tile.data[i] = nodata;
        }
    }

    return tile;
}

lsdtt_raster_tile lsdtt_data_loader::read_full() const
{
    if (!raster_open_)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: read_full() called with no raster open"
        );
    }

    lsdtt_tile_window win;
    win.x_off       = 0;
    win.y_off       = 0;
    win.tile_width  = meta_.n_cols;
    win.tile_height = meta_.n_rows;
    win.halo        = 0;
    win.read_x_off  = 0;
    win.read_y_off  = 0;
    win.read_width  = meta_.n_cols;
    win.read_height = meta_.n_rows;

    return read_tile(win);
}

// -----------------------------------------------------------------------------
// Raster output
// -----------------------------------------------------------------------------

void lsdtt_data_loader::create_output_raster(
    const std::string&       out_path,
    const lsdtt_raster_meta& meta,
    raster_output_format     format,
    int                      n_bands,
    GDALDataType             out_type
)
{
    // Resolve effective format - COG falls back to GTiff for small rasters
    raster_output_format effective_format = format;
    if (format == raster_output_format::cog &&
        meta.n_pixels() < k_cog_min_pixels)
    {
        std::cout << "lsdtt_data_loader: raster has " << meta.n_pixels()
                  << " pixels (below COG threshold of " << k_cog_min_pixels
                  << ") - writing plain GeoTIFF instead of COG.\n";
        effective_format = raster_output_format::gtiff;
    }

    // Choose driver
    const char* driver_name = nullptr;
    if (effective_format == raster_output_format::same_as_input)
    {
        driver_name = input_driver_name_.empty()
                      ? "GTiff"
                      : input_driver_name_.c_str();
    }
    else
    {
        driver_name = "GTiff";  // both cog and gtiff use the GTiff driver at creation
    }

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(driver_name);
    if (!driver)
    {
        // Fallback to GTiff if requested driver is not available
        std::cerr << "lsdtt_data_loader: driver '" << driver_name
                  << "' not available, falling back to GTiff.\n";
        driver = GetGDALDriverManager()->GetDriverByName("GTiff");
        effective_format = raster_output_format::gtiff;
    }
    if (!driver)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: GTiff driver not available"
        );
    }

    // Creation options - tiled DEFLATE for GTiff/COG, minimal for other formats
    char** options = nullptr;
    if (effective_format == raster_output_format::cog ||
        effective_format == raster_output_format::gtiff)
    {
        options = CSLSetNameValue(options, "TILED",      "YES");
        options = CSLSetNameValue(options, "BLOCKXSIZE", "512");
        options = CSLSetNameValue(options, "BLOCKYSIZE", "512");
        options = CSLSetNameValue(options, "COMPRESS",   "DEFLATE");
        options = CSLSetNameValue(options, "PREDICTOR",  "2");
    }

    GDALDataset* raw = driver->Create(
        out_path.c_str(),
        meta.n_cols,
        meta.n_rows,
        n_bands,
        out_type,
        options
    );
    CSLDestroy(options);

    if (!raw)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: could not create output raster: " + out_path
        );
    }

    // Write geotransform
    double gt[6] = {
        meta.x_min,
        meta.cell_size,
        0.0,
        meta.y_min + meta.n_rows * meta.cell_size,  // y_max (top-left northing)
        0.0,
        -meta.cell_size
    };
    raw->SetGeoTransform(gt);

    // Write CRS
    if (!meta.crs_wkt.empty())
    {
        OGRSpatialReference srs;
        srs.importFromWkt(meta.crs_wkt.c_str());
        char* wkt = nullptr;
        srs.exportToWkt(&wkt);
        raw->SetProjection(wkt);
        CPLFree(wkt);
    }

    // Set nodata on all bands
    for (int b = 1; b <= n_bands; ++b)
    {
        raw->GetRasterBand(b)->SetNoDataValue(meta.no_data_value);
    }

    output_datasets_[out_path].reset(raw);
    output_formats_[out_path] = effective_format;
}

void lsdtt_data_loader::write_tile(
    const std::vector<float>& data,
    const lsdtt_tile_window&  window,
    const std::string&        out_path,
    int                       band
)
{
    auto it = output_datasets_.find(out_path);
    if (it == output_datasets_.end())
    {
        throw std::runtime_error(
            "lsdtt_data_loader: write_tile() called for unknown output path: " +
            out_path + " - call create_output_raster() first"
        );
    }

    int expected = window.tile_width * window.tile_height;
    if (static_cast<int>(data.size()) != expected)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: write_tile() data size mismatch - "
            "expected " + std::to_string(expected) +
            ", got "    + std::to_string(data.size())
        );
    }

    GDALRasterBand* out_band = it->second->GetRasterBand(band);
    if (!out_band)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: write_tile() invalid band index: " +
            std::to_string(band)
        );
    }

    CPLErr err = out_band->RasterIO(
        GF_Write,
        window.x_off,
        window.y_off,
        window.tile_width,
        window.tile_height,
        const_cast<float*>(data.data()),
        window.tile_width,
        window.tile_height,
        GDT_Float32,
        0, 0
    );

    if (err != CE_None)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: RasterIO write failed: " +
            std::string(CPLGetLastErrorMsg())
        );
    }
}

void lsdtt_data_loader::write_png(
    const std::vector<float>& data,
    const lsdtt_raster_meta&  meta,
    const std::string&        out_path)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const float ndv = static_cast<float>(meta.no_data_value);

    // -------------------------------------------------------------------------
    // Find the valid-data range (exclude nodata and non-finite values).
    // -------------------------------------------------------------------------
    float dmin =  std::numeric_limits<float>::max();
    float dmax = -std::numeric_limits<float>::max();
    for (const float v : data)
    {
        if (!std::isfinite(v) || v == ndv) continue;
        if (v < dmin) dmin = v;
        if (v > dmax) dmax = v;
    }

    if (!std::isfinite(dmin) || dmax <= dmin)
        dmax = dmin + 1.0f;   // degenerate: all nodata or uniform

    // Maps dmin -> 1, dmax -> 255;  0 is reserved for nodata.
    const float scale = 254.0f / (dmax - dmin);

    // -------------------------------------------------------------------------
    // Build uint8 image array.
    // -------------------------------------------------------------------------
    const auto n_pix = static_cast<std::size_t>(nrows) * ncols;
    std::vector<uint8_t> img(n_pix, 0u);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        const float v = data[i];
        if (!std::isfinite(v) || v == ndv) continue;
        const float norm = 1.0f + (v - dmin) * scale;
        img[i] = static_cast<uint8_t>(std::max(1.0f, std::min(255.0f, norm)));
    }

    // -------------------------------------------------------------------------
    // Write to an in-memory GDAL dataset with full georeferencing.
    // -------------------------------------------------------------------------
    GDALDriver* mem_drv = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!mem_drv)
        throw std::runtime_error("write_png: GDAL MEM driver not available");

    GDALDataset* mem_ds = mem_drv->Create("", ncols, nrows, 1, GDT_Byte, nullptr);
    if (!mem_ds)
        throw std::runtime_error("write_png: failed to create MEM dataset");

    // Same geotransform convention as create_output_raster (top-left origin).
    double gt[6] = {
        meta.x_min,
        meta.cell_size,
        0.0,
        meta.y_max(),        // = y_min + n_rows * cell_size
        0.0,
        -meta.cell_size
    };
    mem_ds->SetGeoTransform(gt);

    if (!meta.crs_wkt.empty())
        mem_ds->SetProjection(meta.crs_wkt.c_str());

    GDALRasterBand* band = mem_ds->GetRasterBand(1);
    band->SetNoDataValue(0.0);

    // Store the stretch range as band metadata so the mapping can be reversed.
    band->SetMetadataItem("LSDTT_DATA_MIN", std::to_string(dmin).c_str());
    band->SetMetadataItem("LSDTT_DATA_MAX", std::to_string(dmax).c_str());

    const CPLErr werr = band->RasterIO(
        GF_Write, 0, 0, ncols, nrows,
        img.data(), ncols, nrows, GDT_Byte, 0, 0);

    if (werr != CE_None)
    {
        GDALClose(mem_ds);
        throw std::runtime_error("write_png: RasterIO to MEM dataset failed");
    }

    // -------------------------------------------------------------------------
    // Copy to PNG.  WORLDFILE=YES writes a .pgw for the geotransform.
    // Closing the PNG dataset triggers the PAM mechanism which writes a
    // .png.aux.xml sidecar carrying the full CRS WKT.
    // -------------------------------------------------------------------------
    GDALDriver* png_drv = GetGDALDriverManager()->GetDriverByName("PNG");
    if (!png_drv)
    {
        GDALClose(mem_ds);
        throw std::runtime_error(
            "write_png: GDAL PNG driver not available - "
            "rebuild GDAL with PNG/libpng support.");
    }

    char** opts = nullptr;
    opts = CSLSetNameValue(opts, "WORLDFILE", "YES");

    GDALDataset* png_ds = png_drv->CreateCopy(
        out_path.c_str(), mem_ds, FALSE, opts, nullptr, nullptr);

    CSLDestroy(opts);
    GDALClose(mem_ds);

    if (!png_ds)
        throw std::runtime_error("write_png: CreateCopy failed for: " + out_path);

    GDALClose(png_ds);   // flushes PNG and triggers .aux.xml write
}

void lsdtt_data_loader::write_full(
    const std::vector<float>& data,
    const lsdtt_raster_meta&  meta,
    const std::string&        out_path,
    raster_output_format      format,
    GDALDataType              out_type
)
{
    if (format == raster_output_format::png)
    {
        write_png(data, meta, out_path);
        return;
    }

    create_output_raster(out_path, meta, format, 1, out_type);

    lsdtt_tile_window win;
    win.x_off       = 0;
    win.y_off       = 0;
    win.tile_width  = meta.n_cols;
    win.tile_height = meta.n_rows;
    win.halo        = 0;
    win.read_x_off  = 0;
    win.read_y_off  = 0;
    win.read_width  = meta.n_cols;
    win.read_height = meta.n_rows;

    write_tile(data, win, out_path);
    finalise_output(out_path);
}

void lsdtt_data_loader::write_full(
    const lsdtt_raster_result& result,
    const std::string&         out_path,
    raster_output_format       format
)
{
    write_full(result.data, result.meta, out_path, format, result.meta.data_type);
}

void lsdtt_data_loader::finalise_output(const std::string& out_path)
{
    auto ds_it  = output_datasets_.find(out_path);
    auto fmt_it = output_formats_.find(out_path);

    if (ds_it == output_datasets_.end())
    {
        throw std::runtime_error(
            "lsdtt_data_loader: finalise_output() called for unknown path: " + out_path
        );
    }

    const raster_output_format fmt = (fmt_it != output_formats_.end())
                                     ? fmt_it->second
                                     : raster_output_format::gtiff;

    if (fmt == raster_output_format::cog)
    {
        // Add overviews then re-copy as proper COG
        add_cog_overviews(ds_it->second.get());

        // Flush and close so GDAL can write the footer
        ds_it->second.reset();
        output_datasets_.erase(ds_it);
        if (fmt_it != output_formats_.end()) output_formats_.erase(fmt_it);

        GDALDriver* cog_driver = GetGDALDriverManager()->GetDriverByName("COG");
        if (cog_driver)
        {
            std::string tmp_path = out_path + ".tmp_cog.tif";

            GDALDataset* src = static_cast<GDALDataset*>(
                GDALOpen(out_path.c_str(), GA_ReadOnly)
            );
            if (!src) return;  // non-fatal - tiled GTiff already written

            char** options = nullptr;
            options = CSLSetNameValue(options, "COMPRESS",  "DEFLATE");
            options = CSLSetNameValue(options, "PREDICTOR", "2");

            GDALDataset* cog = cog_driver->CreateCopy(
                tmp_path.c_str(), src, FALSE, options, nullptr, nullptr
            );
            CSLDestroy(options);
            GDALClose(src);

            if (cog)
            {
                GDALClose(cog);
                VSIUnlink(out_path.c_str());
                VSIRename(tmp_path.c_str(), out_path.c_str());
            }
        }
        // If COG driver unavailable: tiled DEFLATE GTiff is COG-compatible
    }
    else
    {
        // gtiff or same_as_input - just flush and close
        ds_it->second.reset();
        output_datasets_.erase(ds_it);
        if (fmt_it != output_formats_.end()) output_formats_.erase(fmt_it);
    }
}

// -----------------------------------------------------------------------------
// Vector input
// -----------------------------------------------------------------------------
 
void lsdtt_data_loader::open_vector(const std::string& filepath)
{
    if (vector_open_)
    {
        close_vector();
    }
 
    GDALDataset* raw = static_cast<GDALDataset*>(
        GDALOpenEx(filepath.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr)
    );
 
    if (!raw)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: could not open vector dataset: " + filepath +
            "\nGDAL error: " + std::string(CPLGetLastErrorMsg())
        );
    }
 
    vector_dataset_.reset(raw);
    vector_open_ = true;
}
 
std::vector<std::string> lsdtt_data_loader::get_layer_names() const
{
    if (!vector_open_)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: get_layer_names() called with no vector dataset open"
        );
    }
 
    std::vector<std::string> names;
    int n = vector_dataset_->GetLayerCount();
    names.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        OGRLayer* layer = vector_dataset_->GetLayer(i);
        if (layer) names.emplace_back(layer->GetName());
    }
    return names;
}
 
lsdtt_vector_layer lsdtt_data_loader::read_vector_layer(
    const std::string& layer_name
) const
{
    if (!vector_open_)
        throw std::runtime_error(
            "lsdtt_data_loader: read_vector_layer() called with no vector dataset open"
        );
 
    OGRLayer* layer = vector_dataset_->GetLayerByName(layer_name.c_str());
    if (!layer)
        throw std::runtime_error(
            "lsdtt_data_loader: layer not found: " + layer_name
        );
 
    lsdtt_vector_layer result;
 
    // Build schema once from the layer definition.
    OGRFeatureDefn* defn     = layer->GetLayerDefn();
    const int       n_fields = defn->GetFieldCount();
    result.schema.field_names.reserve(n_fields);
    result.schema.field_types.reserve(n_fields);
    for (int f = 0; f < n_fields; ++f)
    {
        OGRFieldDefn* fd = defn->GetFieldDefn(f);
        result.schema.field_names.emplace_back(fd->GetNameRef());
        result.schema.field_types.push_back(from_ogr_field_type(fd->GetType()));
    }
 
    // Read features.
    layer->ResetReading();
    OGRFeature* raw_feat = nullptr;
    std::int64_t feature_index  = 0;
    std::int64_t n_empty        = 0;
    std::int64_t n_exploded     = 0;
    std::int64_t n_unsupported  = 0;
    std::string  unsupported_name;

    while ((raw_feat = layer->GetNextFeature()) != nullptr)
    {
        std::vector<lsdtt_geom> parts =
            read_ogr_geometry_parts(raw_feat->GetGeometryRef(), unsupported_name);

        std::vector<lsdtt_vector_field_value> fields;
        fields.reserve(n_fields);
        for (int f = 0; f < n_fields; ++f)
            fields.push_back(
                get_typed_field_value(raw_feat, f, result.schema.field_types[f])
            );

        if (parts.empty())
        {
            // Nothing readable. The record is still emitted with an empty
            // geometry so the feature count and the field values survive, which
            // is the historical behaviour; what is NEW is that it is counted and
            // reported below rather than vanishing without trace.
            if (raw_feat->GetGeometryRef() == nullptr) ++n_empty;
            else                                       ++n_unsupported;

            lsdtt_vector_record record;
            record.field_values = std::move(fields);
            result.records.push_back(std::move(record));
            result.source_feature_index.push_back(feature_index);
        }
        else
        {
            if (parts.size() > 1) ++n_exploded;
            for (lsdtt_geom& g : parts)
            {
                lsdtt_vector_record record;
                record.geometry     = std::move(g);
                record.field_values = fields;      // copied onto every part
                result.records.push_back(std::move(record));
                result.source_feature_index.push_back(feature_index);
            }
        }

        ++feature_index;
        OGRFeature::DestroyFeature(raw_feat);
    }

    if (n_exploded > 0)
        std::cout << "lsdtt_data_loader: " << n_exploded
                  << " multipart feature(s) in layer '" << layer_name
                  << "' were split into single parts ("
                  << result.records.size() << " records from "
                  << feature_index << " features). Attributes are copied onto "
                     "each part; layer.source_feature_index recovers the "
                     "grouping.\n";

    if (n_unsupported > 0)
        std::cout << "lsdtt_data_loader: WARNING - " << n_unsupported
                  << " feature(s) in layer '" << layer_name
                  << "' have a geometry type this reader does not support"
                  << (unsupported_name.empty() ? std::string()
                                               : " (" + unsupported_name + ")")
                  << " and carry NO geometry. They are counted in the record "
                     "total but every consumer will skip them.\n";

    return result;
}
 

std::string lsdtt_data_loader::get_vector_layer_crs_wkt(
    const std::string& layer_name
) const
{
    if (!vector_open_)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: get_vector_layer_crs_wkt() called with no vector dataset open"
        );
    }

    OGRLayer* layer = nullptr;

    if (!layer_name.empty())
    {
        layer = vector_dataset_->GetLayerByName(layer_name.c_str());
        if (!layer)
        {
            throw std::runtime_error(
                "lsdtt_data_loader: layer not found: " + layer_name
            );
        }
    }
    else
    {
        if (vector_dataset_->GetLayerCount() == 0)
        {
            throw std::runtime_error(
                "lsdtt_data_loader: vector dataset contains no layers"
            );
        }

        layer = vector_dataset_->GetLayer(0);
        if (!layer)
        {
            throw std::runtime_error(
                "lsdtt_data_loader: failed to access first vector layer"
            );
        }
    }

    const OGRSpatialReference* srs = layer->GetSpatialRef();
    if (!srs)
    {
        return "";
    }

    OGRSpatialReference srs_copy(*srs);
    srs_copy.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    char* wkt = nullptr;
    if (srs_copy.exportToWkt(&wkt) != OGRERR_NONE || !wkt)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: failed to export vector layer CRS to WKT"
        );
    }

    std::string out_wkt(wkt);
    CPLFree(wkt);
    return out_wkt;
}




void lsdtt_data_loader::close_vector()
{
    vector_dataset_.reset();
    vector_open_ = false;
}
 
// -----------------------------------------------------------------------------
// Vector output
// -----------------------------------------------------------------------------
 
void lsdtt_data_loader::begin_vector_write(
    const std::string&          out_path,
    const std::string&          layer_name,
    const lsdtt_vector_schema&  schema,
    const std::string&          crs_wkt,
    vector_output_format  format,
    bool                        spatial_index
)
{
    if (write_session_)
        throw std::runtime_error(
            "lsdtt_data_loader: begin_vector_write() called while a write "
            "session is already active - call finish_vector_write() first"
        );
 
    if (format == vector_output_format::gpkg)
    {
        std::cerr << "Warning: writing to GeoPackage. "
                  << "SQLite/GeoPackage I/O can be problematic on some "
                  << "filesystems (e.g. WSL2 /mnt/c). "
                  << "If writing fails, use FlatGeobuf instead.\n";
    }
 
    write_session_ = std::make_unique<vector_write_session>();
    write_session_->out_path      = out_path;
    write_session_->layer_name    = layer_name;
    write_session_->schema        = schema;
    write_session_->crs_wkt       = crs_wkt;
    write_session_->format        = format;
    write_session_->spatial_index = spatial_index;
 
    if (!spatial_index)
    {
        // Streaming path: open the file and layer immediately.
        // We don't know the geometry type yet so use wkbUnknown - this costs
        // a per-feature type tag in FlatGeobuf but avoids a two-pass design.
        write_session_->open_dataset_and_layer(wkbUnknown);
    }
    // Buffered path: defer file creation until finish_vector_write(), when
    // we can infer the geometry type from the first buffered record.
}
 
void lsdtt_data_loader::write_record_to_layer(
    const lsdtt_vector_record& record,
    void*                      ogr_layer_ptr
) const
{
    OGRLayer* layer = static_cast<OGRLayer*>(ogr_layer_ptr);
 
    OGRFeature* feat = OGRFeature::CreateFeature(layer->GetLayerDefn());
 
    OGRGeometry* geom = build_ogr_geometry(record.geometry);
    if (geom) feat->SetGeometryDirectly(geom);
 
    for (int f = 0; f < static_cast<int>(record.field_values.size()); ++f)
        set_typed_field_value(feat, f, record.field_values[static_cast<std::size_t>(f)]);
 
    if (layer->CreateFeature(feat) != OGRERR_NONE)
    {
        OGRFeature::DestroyFeature(feat);
        throw std::runtime_error("lsdtt_data_loader: failed to write feature");
    }
    OGRFeature::DestroyFeature(feat);
}
 
void lsdtt_data_loader::write_point(
    double x, double y,
    const std::vector<lsdtt_vector_field_value>& field_values
)
{
    if (!write_session_)
        throw std::runtime_error(
            "lsdtt_data_loader: write_point() called outside a write session"
        );
 
    lsdtt_vector_record record;
    record.geometry     = lsdtt_point_geom{ x, y };
    record.field_values = field_values;
 
    if (!write_session_->spatial_index)
        write_record_to_layer(record, write_session_->layer);
    else
        write_session_->buffer.records.push_back(std::move(record));
}
 
void lsdtt_data_loader::write_line(
    const std::vector<lsdtt_xy>&                 points,
    const std::vector<lsdtt_vector_field_value>& field_values
)
{
    if (!write_session_)
        throw std::runtime_error(
            "lsdtt_data_loader: write_line() called outside a write session"
        );
 
    lsdtt_vector_record record;
    record.geometry     = lsdtt_line_geom{ points };
    record.field_values = field_values;
 
    if (!write_session_->spatial_index)
        write_record_to_layer(record, write_session_->layer);
    else
        write_session_->buffer.records.push_back(std::move(record));
}
 
void lsdtt_data_loader::write_polygon(
    const std::vector<lsdtt_xy>&                 exterior,
    const std::vector<lsdtt_vector_field_value>& field_values
)
{
    write_polygon(exterior, {}, field_values);
}
 
void lsdtt_data_loader::write_polygon(
    const std::vector<lsdtt_xy>&                 exterior,
    const std::vector<std::vector<lsdtt_xy>>&    holes,
    const std::vector<lsdtt_vector_field_value>& field_values
)
{
    if (!write_session_)
        throw std::runtime_error(
            "lsdtt_data_loader: write_polygon() called outside a write session"
        );
 
    lsdtt_vector_record record;
    record.geometry     = lsdtt_polygon_geom{ exterior, holes };
    record.field_values = field_values;
 
    if (!write_session_->spatial_index)
        write_record_to_layer(record, write_session_->layer);
    else
        write_session_->buffer.records.push_back(std::move(record));
}
 
void lsdtt_data_loader::finish_vector_write()
{
    if (!write_session_)
        throw std::runtime_error(
            "lsdtt_data_loader: finish_vector_write() called with no active session"
        );
 
    if (write_session_->spatial_index)
    {
        // Buffered path: infer geometry type from first record, open the
        // file, then write all buffered records in one pass. GDAL builds
        // its Hilbert R-tree spatial index over the complete dataset.
 
        OGRwkbGeometryType geom_type = wkbUnknown;
        if (!write_session_->buffer.empty())
        {
            const lsdtt_geom& first_geom =
                write_session_->buffer.records[0].geometry;
 
            if (std::holds_alternative<lsdtt_point_geom>(first_geom))
                geom_type = wkbPoint;
            else if (std::holds_alternative<lsdtt_line_geom>(first_geom))
                geom_type = wkbLineString;
            else if (std::holds_alternative<lsdtt_polygon_geom>(first_geom))
                geom_type = wkbPolygon;
        }
 
        write_session_->open_dataset_and_layer(geom_type);
 
        std::cout << "lsdtt_data_loader: writing "
                  << write_session_->buffer.size()
                  << " buffered features with spatial index...\n";
 
        if (write_session_->layer->StartTransaction() != OGRERR_NONE)
            throw std::runtime_error(
                "lsdtt_data_loader: failed to start transaction on layer: " +
                write_session_->layer_name
            );
 
        for (const auto& record : write_session_->buffer.records)
            write_record_to_layer(record, write_session_->layer);
 
        if (write_session_->layer->CommitTransaction() != OGRERR_NONE)
            throw std::runtime_error(
                "lsdtt_data_loader: failed to commit transaction on layer: " +
                write_session_->layer_name
            );
 
        // Free buffer memory before closing the dataset.
        write_session_->buffer.records.clear();
        write_session_->buffer.records.shrink_to_fit();
    }
    else if (write_session_->format == vector_output_format::gpkg)
    {
        // Streaming path, GeoPackage: commit the transaction opened in
        // open_dataset_and_layer().
        if (write_session_->layer->CommitTransaction() != OGRERR_NONE)
            throw std::runtime_error(
                "lsdtt_data_loader: failed to commit transaction on layer: " +
                write_session_->layer_name
            );
    }
 
    // Close dataset (flushes to disk for both paths).
    write_session_->dataset.reset();
    write_session_->layer = nullptr;
 
    std::cout << "lsdtt_data_loader: finished writing to "
              << write_session_->out_path << "\n";
 
    write_session_.reset();
}

// -----------------------------------------------------------------------------
// Tabular output
// -----------------------------------------------------------------------------

void lsdtt_data_loader::write_csv(
    const std::string&                           out_path,
    const std::vector<std::string>&              headers,
    const std::vector<std::vector<std::string>>& rows
)
{
    std::ofstream ofs(out_path);
    if (!ofs.is_open())
    {
        throw std::runtime_error(
            "lsdtt_data_loader: could not open CSV for writing: " + out_path
        );
    }

    auto write_row = [&](const std::vector<std::string>& fields)
    {
        for (std::size_t i = 0; i < fields.size(); ++i)
        {
            // Quote fields containing commas, quotes, or newlines
            const std::string& f = fields[i];
            bool needs_quote = f.find(',')  != std::string::npos ||
                               f.find('"')  != std::string::npos ||
                               f.find('\n') != std::string::npos;
            if (needs_quote)
            {
                ofs << '"';
                for (char c : f)
                {
                    if (c == '"') ofs << '"';  // escape embedded quotes
                    ofs << c;
                }
                ofs << '"';
            }
            else
            {
                ofs << f;
            }
            if (i + 1 < fields.size()) ofs << ',';
        }
        ofs << '\n';
    };

    write_row(headers);
    for (const auto& row : rows)
    {
        if (row.size() != headers.size())
        {
            throw std::runtime_error(
                "lsdtt_data_loader: write_csv() row has " +
                std::to_string(row.size()) + " fields but header has " +
                std::to_string(headers.size())
            );
        }
        write_row(row);
    }
}

// -----------------------------------------------------------------------------
// Utility
// -----------------------------------------------------------------------------

void lsdtt_data_loader::print_raster_info() const
{
    if (!raster_open_)
    {
        std::cout << "lsdtt_data_loader: no raster currently open\n";
        return;
    }
    std::cout << "Raster info:\n"
              << "  Dimensions : " << meta_.n_cols << " x " << meta_.n_rows << "\n"
              << "  Cell size  : " << meta_.cell_size << "\n"
              << "  Origin     : (" << meta_.x_min << ", " << meta_.y_min << ")\n"
              << "  NoData     : " << meta_.no_data_value << "\n"
              << "  Data type  : " << GDALGetDataTypeName(meta_.data_type) << "\n"
              << "  CRS (WKT)  : " << meta_.crs_wkt.substr(0, 80)
              << (meta_.crs_wkt.size() > 80 ? "..." : "") << "\n";
}

std::string lsdtt_data_loader::gdal_version()
{
    return std::string(GDALVersionInfo("RELEASE_NAME"));
}

// -----------------------------------------------------------------------------
// Private helpers
// -----------------------------------------------------------------------------

void lsdtt_data_loader::read_meta_from_dataset()
{
    meta_.n_cols = raster_dataset_->GetRasterXSize();
    meta_.n_rows = raster_dataset_->GetRasterYSize();

    double gt[6];
    if (raster_dataset_->GetGeoTransform(gt) == CE_None)
    {
        meta_.x_min     = gt[0];
        meta_.cell_size = gt[1];
        // gt[3] is the top-left y (northing of first row)
        // y_min is the bottom of the raster
        meta_.y_min     = gt[3] + meta_.n_rows * gt[5];
    }
    else
    {
        meta_.x_min     = 0.0;
        meta_.cell_size = 1.0;
        meta_.y_min     = 0.0;
    }

    // CRS
    const char* wkt = raster_dataset_->GetProjectionRef();
    meta_.crs_wkt = wkt ? std::string(wkt) : "";

    // NoData and data type from band 1
    GDALRasterBand* band = raster_dataset_->GetRasterBand(1);
    if (band)
    {
        int    has_nodata = 0;
        double nd         = band->GetNoDataValue(&has_nodata);
        meta_.no_data_value = has_nodata ? nd : -9999.0;
        meta_.data_type     = band->GetRasterDataType();
    }
    else
    {
        meta_.no_data_value = -9999.0;
        meta_.data_type     = GDT_Float32;
    }
}

lsdtt_tile_window lsdtt_data_loader::clamp_window_to_raster(
    const lsdtt_tile_window& win
) const
{
    lsdtt_tile_window clamped = win;

    // Read origin - pull back by halo, clamp to raster edge
    clamped.read_x_off = std::max(0, win.x_off - win.halo);
    clamped.read_y_off = std::max(0, win.y_off - win.halo);

    // Read extent - push forward by halo, clamp to raster edge
    int read_x_end = std::min(meta_.n_cols, win.x_off + win.tile_width  + win.halo);
    int read_y_end = std::min(meta_.n_rows, win.y_off + win.tile_height + win.halo);

    clamped.read_width  = read_x_end - clamped.read_x_off;
    clamped.read_height = read_y_end - clamped.read_y_off;

    return clamped;
}

void lsdtt_data_loader::add_cog_overviews(GDALDataset* ds)
{
    // Build power-of-2 overview levels down to ~256px
    int min_dim = std::min(ds->GetRasterXSize(), ds->GetRasterYSize());
    std::vector<int> levels;
    int level = 2;
    while (min_dim / level > 256)
    {
        levels.push_back(level);
        level *= 2;
    }
    if (levels.empty()) return;

    ds->BuildOverviews(
        "NEAREST",
        static_cast<int>(levels.size()),
        levels.data(),
        0, nullptr,
        GDALDummyProgress, nullptr
    );
}

// =============================================================================
// clip_vector_to_raster_bbox
// =============================================================================
void lsdtt_data_loader::clip_vector_to_raster_bbox(
    const std::string&       input_path,
    const std::string&       output_path,
    const lsdtt_raster_meta& raster_meta,
    vector_output_format     out_fmt
)
{
    ensure_gdal_init();

    const char* driver_name =
        (out_fmt == vector_output_format::gpkg) ? "GPKG" : "FlatGeobuf";

    struct ds_closer
    {
        void operator()(GDALDataset* d) const { if (d) GDALClose(d); }
    };

    // Open input vector
    GDALDataset* raw_src = static_cast<GDALDataset*>(
        GDALOpenEx(input_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr)
    );
    if (!raw_src)
        throw std::runtime_error(
            "clip_vector_to_raster_bbox: cannot open input vector: " +
            input_path + "\nGDAL error: " + std::string(CPLGetLastErrorMsg())
        );
    std::unique_ptr<GDALDataset, ds_closer> src_guard(raw_src);

    if (raw_src->GetLayerCount() == 0)
        throw std::runtime_error(
            "clip_vector_to_raster_bbox: no layers in: " + input_path
        );

    // Create output dataset
    GDALDriver* out_driver =
        GetGDALDriverManager()->GetDriverByName(driver_name);
    if (!out_driver)
        throw std::runtime_error(
            std::string("clip_vector_to_raster_bbox: GDAL driver not available: ")
            + driver_name
        );

    VSIStatBufL stat_buf;
    if (VSIStatL(output_path.c_str(), &stat_buf) == 0)
        out_driver->Delete(output_path.c_str());

    GDALDataset* raw_dst = out_driver->Create(
        output_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr
    );
    if (!raw_dst)
        throw std::runtime_error(
            "clip_vector_to_raster_bbox: cannot create output: " + output_path
        );
    std::unique_ptr<GDALDataset, ds_closer> dst_guard(raw_dst);

    const int n_layers = raw_src->GetLayerCount();
    for (int li = 0; li < n_layers; ++li)
    {
        OGRLayer* src_layer = raw_src->GetLayer(li);
        if (!src_layer) continue;

        // Compute the filter bbox in the layer's CRS
        double filter_minx = raster_meta.x_min;
        double filter_miny = raster_meta.y_min;
        double filter_maxx = raster_meta.x_max();
        double filter_maxy = raster_meta.y_max();

        const OGRSpatialReference* layer_srs = src_layer->GetSpatialRef();
        if (layer_srs && !raster_meta.crs_wkt.empty())
        {
            OGRSpatialReference raster_srs;
            raster_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

            if (raster_srs.importFromWkt(raster_meta.crs_wkt.c_str()) == OGRERR_NONE)
            {
                OGRSpatialReference dst_srs = *layer_srs;
                dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

                OGRCoordinateTransformation* ct =
                    OGRCreateCoordinateTransformation(&raster_srs, &dst_srs);

                if (ct)
                {
                    double xs[4] = {raster_meta.x_min,  raster_meta.x_max(),
                                    raster_meta.x_min,  raster_meta.x_max()};
                    double ys[4] = {raster_meta.y_min,  raster_meta.y_min,
                                    raster_meta.y_max(), raster_meta.y_max()};

                    if (ct->Transform(4, xs, ys))
                    {
                        filter_minx = *std::min_element(xs, xs + 4);
                        filter_maxx = *std::max_element(xs, xs + 4);
                        filter_miny = *std::min_element(ys, ys + 4);
                        filter_maxy = *std::max_element(ys, ys + 4);
                    }
                    OCTDestroyCoordinateTransformation(ct);
                }
            }
        }

        src_layer->SetSpatialFilterRect(filter_minx, filter_miny,
                                        filter_maxx, filter_maxy);

        // Create output layer matching source schema and CRS
        OGRLayer* dst_layer = raw_dst->CreateLayer(
            src_layer->GetName(), const_cast<OGRSpatialReference*>(layer_srs), src_layer->GetGeomType(), nullptr
        );
        if (!dst_layer)
            throw std::runtime_error(
                "clip_vector_to_raster_bbox: failed to create output layer"
            );

        OGRFeatureDefn* defn = src_layer->GetLayerDefn();
        for (int fi = 0; fi < defn->GetFieldCount(); ++fi)
        {
            if (dst_layer->CreateField(defn->GetFieldDefn(fi)) != OGRERR_NONE)
                throw std::runtime_error(
                    "clip_vector_to_raster_bbox: failed to copy field definition"
                );
        }

        // Stream filtered features to output
        src_layer->ResetReading();
        OGRFeature* feat = nullptr;
        while ((feat = src_layer->GetNextFeature()) != nullptr)
        {
            OGRFeature* out_feat =
                OGRFeature::CreateFeature(dst_layer->GetLayerDefn());
            out_feat->SetFrom(feat);
            if (dst_layer->CreateFeature(out_feat) != OGRERR_NONE)
            {
                OGRFeature::DestroyFeature(out_feat);
                OGRFeature::DestroyFeature(feat);
                throw std::runtime_error(
                    "clip_vector_to_raster_bbox: failed to write feature"
                );
            }
            OGRFeature::DestroyFeature(out_feat);
            OGRFeature::DestroyFeature(feat);
        }

        src_layer->SetSpatialFilter(nullptr);
    }
}

// count_features_in_first_layer
// =============================================================================
int lsdtt_data_loader::count_features_in_first_layer(
    const std::string& vector_path)
{
    GDALAllRegister();

    GDALDataset* raw = static_cast<GDALDataset*>(
        GDALOpenEx(vector_path.c_str(),
                   GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr)
    );

    if (!raw)
    {
        throw std::runtime_error(
            "count_features_in_first_layer: cannot open: " + vector_path +
            "\nGDAL error: " + std::string(CPLGetLastErrorMsg())
        );
    }

    struct ds_closer { void operator()(GDALDataset* d) const { if (d) GDALClose(d); } };
    std::unique_ptr<GDALDataset, ds_closer> guard(raw);

    OGRLayer* layer = raw->GetLayer(0);
    if (!layer)
    {
        throw std::runtime_error(
            "count_features_in_first_layer: no layers in: " + vector_path
        );
    }

    return static_cast<int>(layer->GetFeatureCount(/*bForce=*/TRUE));
}

// =============================================================================
// resolve_writable_path
// =============================================================================
std::string lsdtt_data_loader::resolve_writable_path(
    const std::string& base,
    const std::string& ext,
    bool*              was_renamed)
{
    if (was_renamed) *was_renamed = false;

    auto exists = [](const std::string& p) -> bool
    {
        std::ifstream f(p);
        return f.good();
    };

    const std::string path = base + ext;
    if (!exists(path))
    {
        return path;   // nothing there: free to create
    }

    if (std::remove(path.c_str()) == 0)
    {
        return path;   // removed (not locked): caller will recreate in place
    }

    // Could not remove - treat as locked; find the next free fallback name.
    for (int n = 2; n < 1000; ++n)
    {
        const std::string cand = base + "_" + std::to_string(n) + ext;
        if (!exists(cand))
        {
            if (was_renamed) *was_renamed = true;
            return cand;
        }
    }
    // Exhausted fallbacks: return the original (write will likely fail loudly).
    return path;
}

} // namespace lsdtt
