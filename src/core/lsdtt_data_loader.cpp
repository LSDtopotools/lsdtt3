// lsdtt_data_loader.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_data_loader.hpp"

#include <algorithm>
#include <cmath>
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
// GDAL initialisation helper — called once via static local
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

void validate_vector_feature_schema(
    const std::vector<lsdtt::lsdtt_vector_feature>& features,
    const std::string&                              where
)
{
    if (features.empty())
    {
        return;
    }

    const auto& names0 = features[0].field_names;
    const auto& types0 = features[0].field_types;
    const auto& vals0  = features[0].field_values;

    if (names0.size() != types0.size() || names0.size() != vals0.size())
    {
        throw std::runtime_error(
            where + ": first feature has inconsistent field schema lengths"
        );
    }

    for (std::size_t i = 1; i < features.size(); ++i)
    {
        const auto& feat = features[i];

        if (feat.field_names.size() != feat.field_types.size() ||
            feat.field_names.size() != feat.field_values.size())
        {
            throw std::runtime_error(
                where + ": feature has inconsistent field schema lengths"
            );
        }

        if (feat.field_names != names0)
        {
            throw std::runtime_error(
                where + ": all features in a layer must share identical field names"
            );
        }

        if (feat.field_types != types0)
        {
            throw std::runtime_error(
                where + ": all features in a layer must share identical field types"
            );
        }
    }
}

} // unnamed namespace


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

    // Build nodata mask
    float nodata = static_cast<float>(meta_.no_data_value);
    for (int i = 0; i < buf_size; ++i)
    {
        if (tile.data[i] == nodata)
        {
            tile.mask[i] = 0;
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
    int                      n_bands
)
{
    // Resolve effective format — COG falls back to GTiff for small rasters
    raster_output_format effective_format = format;
    if (format == raster_output_format::cog &&
        meta.n_pixels() < k_cog_min_pixels)
    {
        std::cout << "lsdtt_data_loader: raster has " << meta.n_pixels()
                  << " pixels (below COG threshold of " << k_cog_min_pixels
                  << ") — writing plain GeoTIFF instead of COG.\n";
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

    // Creation options — tiled DEFLATE for GTiff/COG, minimal for other formats
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
        GDT_Float32,
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
            out_path + " — call create_output_raster() first"
        );
    }

    int expected = window.tile_width * window.tile_height;
    if (static_cast<int>(data.size()) != expected)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: write_tile() data size mismatch — "
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

void lsdtt_data_loader::write_full(
    const std::vector<float>& data,
    const lsdtt_raster_meta&  meta,
    const std::string&        out_path,
    raster_output_format      format
)
{
    create_output_raster(out_path, meta, format);

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
            if (!src) return;  // non-fatal — tiled GTiff already written

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
        // gtiff or same_as_input — just flush and close
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

std::vector<lsdtt_vector_feature> lsdtt_data_loader::read_layer(
    const std::string& layer_name
) const
{
    if (!vector_open_)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: read_layer() called with no vector dataset open"
        );
    }

    OGRLayer* layer = vector_dataset_->GetLayerByName(layer_name.c_str());
    if (!layer)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: layer not found: " + layer_name
        );
    }

    std::vector<lsdtt_vector_feature> features;
    layer->ResetReading();

    OGRFeatureDefn* defn     = layer->GetLayerDefn();
    int             n_fields = defn->GetFieldCount();

    OGRFeature* raw_feat = nullptr;
    while ((raw_feat = layer->GetNextFeature()) != nullptr)
    {
        lsdtt_vector_feature feat;

        // Geometry as WKT
        OGRGeometry* geom = raw_feat->GetGeometryRef();
        if (geom)
        {
            char* wkt = nullptr;
            geom->exportToWkt(&wkt);
            feat.geometry_wkt = wkt ? wkt : "";
            CPLFree(wkt);
        }

        // Typed fields
        feat.field_names.reserve(n_fields);
        feat.field_types.reserve(n_fields);
        feat.field_values.reserve(n_fields);

        for (int f = 0; f < n_fields; ++f)
        {
            OGRFieldDefn* field_defn = defn->GetFieldDefn(f);

            feat.field_names.emplace_back(field_defn->GetNameRef());

            const lsdtt_vector_field_type field_type =
                from_ogr_field_type(field_defn->GetType());

            feat.field_types.push_back(field_type);
            feat.field_values.push_back(
                get_typed_field_value(raw_feat, f, field_type)
            );
        }

        features.push_back(std::move(feat));
        OGRFeature::DestroyFeature(raw_feat);
    }

    return features;
}


void lsdtt_data_loader::close_vector()
{
    vector_dataset_.reset();
    vector_open_ = false;
}

// -----------------------------------------------------------------------------
// Vector output
// -----------------------------------------------------------------------------
void lsdtt_data_loader::write_gpkg(
    const std::string&                       out_path,
    const std::string&                       layer_name,
    const std::vector<lsdtt_vector_feature>& features,
    const std::string&                       crs_wkt
)
{
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GPKG");
    if (!driver)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: GeoPackage (GPKG) driver not available"
        );
    }

    validate_vector_feature_schema(features, "lsdtt_data_loader::write_gpkg");

    // Open existing or create new GPKG
    GDALDataset* ds = nullptr;
    VSIStatBufL stat;
    if (VSIStatL(out_path.c_str(), &stat) == 0)
    {
        ds = static_cast<GDALDataset*>(
            GDALOpenEx(out_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE,
                       nullptr, nullptr, nullptr)
        );
    }
    if (!ds)
    {
        ds = driver->Create(out_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    }
    if (!ds)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: could not create or open GeoPackage: " + out_path
        );
    }

    gdal_dataset_ptr ds_ptr(ds);

    OGRSpatialReference srs;
    if (!crs_wkt.empty())
    {
        srs.importFromWkt(crs_wkt.c_str());
    }

    int existing = -1;
    for (int i = 0; i < ds->GetLayerCount(); ++i)
    {
        OGRLayer* l = ds->GetLayer(i);
        if (l && layer_name == l->GetName())
        {
            existing = i;
            break;
        }
    }
    if (existing >= 0)
    {
        ds->DeleteLayer(existing);
    }

    OGRLayer* layer = ds->CreateLayer(
        layer_name.c_str(),
        crs_wkt.empty() ? nullptr : &srs,
        wkbUnknown,
        nullptr
    );
    if (!layer)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: could not create layer: " + layer_name
        );
    }

    // Create typed fields from first feature
    if (!features.empty())
    {
        for (std::size_t f = 0; f < features[0].field_names.size(); ++f)
        {
            OGRFieldDefn field(
                features[0].field_names[f].c_str(),
                to_ogr_field_type(features[0].field_types[f])
            );

            if (layer->CreateField(&field) != OGRERR_NONE)
            {
                throw std::runtime_error(
                    "lsdtt_data_loader: failed to create field: " +
                    features[0].field_names[f]
                );
            }
        }
    }

    for (const auto& feat : features)
    {
        OGRFeature* ogr_feat = OGRFeature::CreateFeature(layer->GetLayerDefn());

        if (!feat.geometry_wkt.empty())
        {
            OGRGeometry* geom = nullptr;
            OGRGeometryFactory::createFromWkt(
                feat.geometry_wkt.c_str(), nullptr, &geom
            );
            if (geom)
            {
                ogr_feat->SetGeometryDirectly(geom);
            }
        }

        for (int f = 0; f < static_cast<int>(feat.field_values.size()); ++f)
        {
            set_typed_field_value(
                ogr_feat,
                f,
                feat.field_values[static_cast<std::size_t>(f)]
            );
        }

        if (layer->CreateFeature(ogr_feat) != OGRERR_NONE)
        {
            OGRFeature::DestroyFeature(ogr_feat);
            throw std::runtime_error(
                "lsdtt_data_loader: failed to write feature to layer: " + layer_name
            );
        }

        OGRFeature::DestroyFeature(ogr_feat);
    }
}




void lsdtt_data_loader::write_flatgeobuf(
    const std::string&                       out_path,
    const std::string&                       layer_name,
    const std::vector<lsdtt_vector_feature>& features,
    const std::string&                       crs_wkt
)
{
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
    if (!driver)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: FlatGeobuf driver not available"
        );
    }

    validate_vector_feature_schema(features, "lsdtt_data_loader::write_flatgeobuf");

    VSIUnlink(out_path.c_str());

    GDALDataset* ds = driver->Create(out_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: could not create FlatGeobuf file: " + out_path
        );
    }

    gdal_dataset_ptr ds_ptr(ds);

    OGRSpatialReference srs;
    if (!crs_wkt.empty())
    {
        srs.importFromWkt(crs_wkt.c_str());
    }

    OGRLayer* layer = ds->CreateLayer(
        layer_name.c_str(),
        crs_wkt.empty() ? nullptr : &srs,
        wkbUnknown,
        nullptr
    );

    if (!layer)
    {
        throw std::runtime_error(
            "lsdtt_data_loader: could not create FlatGeobuf layer: " + layer_name
        );
    }

    if (!features.empty())
    {
        for (std::size_t f = 0; f < features[0].field_names.size(); ++f)
        {
            OGRFieldDefn field(
                features[0].field_names[f].c_str(),
                to_ogr_field_type(features[0].field_types[f])
            );

            if (layer->CreateField(&field) != OGRERR_NONE)
            {
                throw std::runtime_error(
                    "lsdtt_data_loader: failed creating field '" +
                    features[0].field_names[f] + "' in FlatGeobuf layer"
                );
            }
        }
    }

    for (const auto& feat : features)
    {
        OGRFeature* ogr_feat = OGRFeature::CreateFeature(layer->GetLayerDefn());

        if (!feat.geometry_wkt.empty())
        {
            OGRGeometry* geom = nullptr;
            OGRGeometryFactory::createFromWkt(
                feat.geometry_wkt.c_str(), nullptr, &geom
            );
            if (geom)
            {
                ogr_feat->SetGeometryDirectly(geom);
            }
        }

        for (int f = 0; f < static_cast<int>(feat.field_values.size()); ++f)
        {
            set_typed_field_value(
                ogr_feat,
                f,
                feat.field_values[static_cast<std::size_t>(f)]
            );
        }

        if (layer->CreateFeature(ogr_feat) != OGRERR_NONE)
        {
            OGRFeature::DestroyFeature(ogr_feat);
            throw std::runtime_error(
                "lsdtt_data_loader: failed to write feature to FlatGeobuf"
            );
        }

        OGRFeature::DestroyFeature(ogr_feat);
    }
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

    // Read origin — pull back by halo, clamp to raster edge
    clamped.read_x_off = std::max(0, win.x_off - win.halo);
    clamped.read_y_off = std::max(0, win.y_off - win.halo);

    // Read extent — push forward by halo, clamp to raster edge
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

} // namespace lsdtt
