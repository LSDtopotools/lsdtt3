// lsdtt_raster_contours_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_contours_cpu.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <cpl_conv.h>
#include <cpl_vsi.h>
#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// ensure_gdal_init
// -----------------------------------------------------------------------------
void ensure_gdal_init()
{
    static bool initialised = false;
    if (!initialised)
    {
        GDALAllRegister();
        initialised = true;
    }
}

// -----------------------------------------------------------------------------
// GDAL dataset deleter
// -----------------------------------------------------------------------------
struct gdal_dataset_deleter
{
    void operator()(GDALDataset* ds) const
    {
        if (ds) GDALClose(ds);
    }
};

using gdal_dataset_ptr = std::unique_ptr<GDALDataset, gdal_dataset_deleter>;

// -----------------------------------------------------------------------------
// driver_name_from_format
// -----------------------------------------------------------------------------
const char* driver_name_from_format(vector_output_format format)
{
    switch (format)
    {
        case vector_output_format::gpkg:
            return "GPKG";

        case vector_output_format::flatgeobuf:
            return "FlatGeobuf";
    }

    throw std::invalid_argument(
        "lsdtt::write_raster_contours: unsupported vector_output_format"
    );
}

// -----------------------------------------------------------------------------
// validate_request
// -----------------------------------------------------------------------------
void validate_request(
    const lsdtt_raster_cpu&      raster,
    const std::string&           out_path,
    const raster_contour_request& request)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: raster metadata is invalid"
        );
    }

    if (out_path.empty())
    {
        throw std::invalid_argument(
            "lsdtt::write_raster_contours: out_path must not be empty"
        );
    }

    if (request.layer_name.empty())
    {
        throw std::invalid_argument(
            "lsdtt::write_raster_contours: layer_name must not be empty"
        );
    }

    if (request.elevation_field_name.empty())
    {
        throw std::invalid_argument(
            "lsdtt::write_raster_contours: elevation_field_name must not be empty"
        );
    }

    if (request.level_mode == contour_level_mode::fixed_levels)
    {
        if (request.fixed_levels.empty())
        {
            throw std::invalid_argument(
                "lsdtt::write_raster_contours: fixed_levels must not be empty"
            );
        }
    }
    else
    {
        if (!(request.contour_interval > 0.0))
        {
            throw std::invalid_argument(
                "lsdtt::write_raster_contours: contour_interval must be > 0"
            );
        }
    }

    const std::string expected_ext = output_extension(request.output_format);
    if (out_path.size() < expected_ext.size() ||
        out_path.substr(out_path.size() - expected_ext.size()) != expected_ext)
    {
        throw std::invalid_argument(
            "lsdtt::write_raster_contours: output extension does not match output_format"
        );
    }
}

// -----------------------------------------------------------------------------
// sorted_unique_levels
// -----------------------------------------------------------------------------
std::vector<double> sorted_unique_levels(
    const std::vector<double>& levels)
{
    std::vector<double> out = levels;
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

// -----------------------------------------------------------------------------
// make_mem_raster_from_lsdtt_raster
//
// Creates a GDAL MEM raster dataset from an lsdtt_raster_cpu.
// Invalid mask pixels are written as nodata if use_nodata = true.
// -----------------------------------------------------------------------------
gdal_dataset_ptr make_mem_raster_from_lsdtt_raster(
    const lsdtt_raster_cpu& raster,
    bool                    use_nodata)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    GDALDriver* mem_driver =
        GetGDALDriverManager()->GetDriverByName("MEM");
    if (!mem_driver)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: MEM driver not available"
        );
    }

    GDALDataset* raw_ds =
        mem_driver->Create("", ncols, nrows, 1, GDT_Float32, nullptr);

    if (!raw_ds)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: could not create MEM raster dataset"
        );
    }

    gdal_dataset_ptr ds(raw_ds);

    // Geotransform
    double gt[6] =
    {
        meta.x_min,
        meta.cell_size,
        0.0,
        meta.y_min + static_cast<double>(meta.n_rows) * meta.cell_size,
        0.0,
        -meta.cell_size
    };

    if (ds->SetGeoTransform(gt) != CE_None)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: failed to set MEM raster geotransform"
        );
    }

    // Projection
    if (!meta.crs_wkt.empty())
    {
        if (ds->SetProjection(meta.crs_wkt.c_str()) != CE_None)
        {
            throw std::runtime_error(
                "lsdtt::write_raster_contours: failed to set MEM raster projection"
            );
        }
    }

    GDALRasterBand* band = ds->GetRasterBand(1);
    if (!band)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: could not get MEM raster band"
        );
    }

    const std::vector<float>&   in_data = raster.get_data();
    const std::vector<uint8_t>& in_mask = raster.get_mask();

    if (in_data.size() != in_mask.size())
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: raster data/mask size mismatch"
        );
    }

    std::vector<float> work = in_data;

    if (use_nodata)
    {
        const float ndv = static_cast<float>(meta.no_data_value);

        for (std::size_t i = 0; i < work.size(); ++i)
        {
            if (in_mask[i] == 0)
            {
                work[i] = ndv;
            }
        }

        if (band->SetNoDataValue(meta.no_data_value) != CE_None)
        {
            throw std::runtime_error(
                "lsdtt::write_raster_contours: failed to set MEM raster nodata"
            );
        }
    }

    CPLErr err = band->RasterIO(
        GF_Write,
        0,
        0,
        ncols,
        nrows,
        work.data(),
        ncols,
        nrows,
        GDT_Float32,
        0,
        0
    );

    if (err != CE_None)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: failed to write MEM raster band"
        );
    }

    return ds;
}

// -----------------------------------------------------------------------------
// create_output_vector_dataset
// -----------------------------------------------------------------------------
gdal_dataset_ptr create_output_vector_dataset(
    const std::string&   out_path,
    vector_output_format format,
    bool                 overwrite)
{
    const char* driver_name = driver_name_from_format(format);

    GDALDriver* driver =
        GetGDALDriverManager()->GetDriverByName(driver_name);

    if (!driver)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: output vector driver not available"
        );
    }

    VSIStatBufL stat_buf;
    const bool exists = VSIStatL(out_path.c_str(), &stat_buf) == 0;

    if (exists)
    {
        if (!overwrite)
        {
            throw std::runtime_error(
                "lsdtt::write_raster_contours: output file exists and overwrite=false"
            );
        }

        // Prefer driver delete where available.
        if (driver->Delete(out_path.c_str()) != CE_None)
        {
            // Fallback to VSI unlink for single-file outputs.
            VSIUnlink(out_path.c_str());
        }
    }

    GDALDataset* raw_ds =
        driver->Create(out_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);

    if (!raw_ds)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: could not create output vector dataset"
        );
    }

    return gdal_dataset_ptr(raw_ds);
}

// -----------------------------------------------------------------------------
// create_output_layer
// -----------------------------------------------------------------------------
OGRLayer* create_output_layer(
    GDALDataset*                 ds,
    const lsdtt_raster_meta&     meta,
    const raster_contour_request& request,
    int&                         id_field_index,
    int&                         elev_field_index)
{
    OGRSpatialReference srs;
    OGRSpatialReference* srs_ptr = nullptr;

    if (!meta.crs_wkt.empty())
    {
        if (srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE)
        {
            throw std::runtime_error(
                "lsdtt::write_raster_contours: failed to parse CRS WKT for output layer"
            );
        }
        srs_ptr = &srs;
    }

    OGRLayer* layer = ds->CreateLayer(
        request.layer_name.c_str(),
        srs_ptr,
        wkbLineString,
        nullptr
    );

    if (!layer)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: could not create output layer"
        );
    }

    id_field_index = -1;
    elev_field_index = -1;

    if (!request.id_field_name.empty())
    {
        OGRFieldDefn id_field(
            request.id_field_name.c_str(),
            OFTInteger64
        );

        if (layer->CreateField(&id_field) != OGRERR_NONE)
        {
            throw std::runtime_error(
                "lsdtt::write_raster_contours: failed to create id field"
            );
        }

        id_field_index =
            layer->GetLayerDefn()->GetFieldIndex(request.id_field_name.c_str());
    }

    OGRFieldDefn elev_field(
        request.elevation_field_name.c_str(),
        OFTReal
    );

    if (layer->CreateField(&elev_field) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: failed to create elevation field"
        );
    }

    elev_field_index =
        layer->GetLayerDefn()->GetFieldIndex(request.elevation_field_name.c_str());

    if (elev_field_index < 0)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: invalid elevation field index"
        );
    }

    return layer;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// write_raster_contours
// -----------------------------------------------------------------------------
raster_contour_result write_raster_contours(
    const lsdtt_raster_cpu&      raster,
    const std::string&           out_path,
    const raster_contour_request& request)
{
    ensure_gdal_init();
    validate_request(raster, out_path, request);

    gdal_dataset_ptr mem_ds =
        make_mem_raster_from_lsdtt_raster(raster, request.use_nodata);

    GDALRasterBand* band = mem_ds->GetRasterBand(1);
    if (!band)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: missing MEM raster band"
        );
    }

    gdal_dataset_ptr out_ds =
        create_output_vector_dataset(
            out_path,
            request.output_format,
            request.overwrite
        );

    int id_field_index   = -1;
    int elev_field_index = -1;

    OGRLayer* layer =
        create_output_layer(
            out_ds.get(),
            raster.get_meta(),
            request,
            id_field_index,
            elev_field_index
        );

    std::vector<double> fixed_levels;
    int n_fixed_levels = 0;
    double contour_interval = 0.0;
    double contour_base     = 0.0;

    if (request.level_mode == contour_level_mode::fixed_levels)
    {
        fixed_levels = sorted_unique_levels(request.fixed_levels);
        n_fixed_levels = static_cast<int>(fixed_levels.size());
    }
    else
    {
        contour_interval = request.contour_interval;
        contour_base     = request.contour_base;
    }

    int has_nodata = FALSE;
    const double nodata = band->GetNoDataValue(&has_nodata);

    GDALProgressFunc progress_func =
        request.verbose ? GDALTermProgress : GDALDummyProgress;

    CPLErr err = GDALContourGenerate(
        GDALRasterBand::ToHandle(band),
        contour_interval,
        contour_base,
        n_fixed_levels,
        n_fixed_levels > 0 ? fixed_levels.data() : nullptr,
        request.use_nodata && has_nodata,
        nodata,
        OGRLayer::ToHandle(layer),
        id_field_index,
        elev_field_index,
        progress_func,
        nullptr
    );

    if (err != CE_None)
    {
        throw std::runtime_error(
            "lsdtt::write_raster_contours: GDALContourGenerate failed"
        );
    }

    raster_contour_result result;
    result.out_path      = out_path;
    result.layer_name    = request.layer_name;
    result.output_format = request.output_format;
    result.n_features    = static_cast<std::int64_t>(layer->GetFeatureCount());

    return result;
}

// -----------------------------------------------------------------------------
// write_raster_contours_fixed_levels
// -----------------------------------------------------------------------------
raster_contour_result write_raster_contours_fixed_levels(
    const lsdtt_raster_cpu&    raster,
    const std::vector<double>& fixed_levels,
    const std::string&         out_path,
    vector_output_format       output_format,
    const std::string&         layer_name)
{
    raster_contour_request request;
    request.level_mode    = contour_level_mode::fixed_levels;
    request.fixed_levels  = fixed_levels;
    request.output_format = output_format;
    request.layer_name    = layer_name;

    return write_raster_contours(raster, out_path, request);
}

// -----------------------------------------------------------------------------
// write_raster_contours_interval
// -----------------------------------------------------------------------------
raster_contour_result write_raster_contours_interval(
    const lsdtt_raster_cpu& raster,
    double                  contour_interval,
    double                  contour_base,
    const std::string&      out_path,
    vector_output_format    output_format,
    const std::string&      layer_name)
{
    raster_contour_request request;
    request.level_mode       = contour_level_mode::interval;
    request.contour_interval = contour_interval;
    request.contour_base     = contour_base;
    request.output_format    = output_format;
    request.layer_name       = layer_name;

    return write_raster_contours(raster, out_path, request);
}

} // namespace lsdtt
