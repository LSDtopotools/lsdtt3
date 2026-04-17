// lsdtt_vector_rasterise_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_vector_rasterise_cpu.hpp"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal.h>
#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt {
namespace {

// Returns the GDAL geotransform for a raster described by `meta`.
// GDAL convention: GT[0] = x of left edge, GT[3] = y of top edge,
// GT[1] = +cell_size, GT[5] = -cell_size.
void fill_geotransform(const lsdtt_raster_meta& meta, double gt[6])
{
    gt[0] = meta.x_min;
    gt[1] = meta.cell_size;
    gt[2] = 0.0;
    gt[3] = meta.y_min + static_cast<double>(meta.n_rows) * meta.cell_size;
    gt[4] = 0.0;
    gt[5] = -meta.cell_size;
}

// Creates a single-band GDT_Byte in-memory raster of size (ncols x nrows)
// with the given geotransform.  Initialised to all-zero.
// Caller owns the returned dataset (must GDALClose).
GDALDataset* make_mem_byte_raster(int nrows, int ncols, const double gt[6])
{
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!drv)
        throw std::runtime_error("lsdtt_vector_rasterise: MEM driver unavailable");

    GDALDataset* ds = drv->Create("", ncols, nrows, 1, GDT_Byte, nullptr);
    if (!ds)
        throw std::runtime_error("lsdtt_vector_rasterise: failed to create MEM raster");

    ds->SetGeoTransform(const_cast<double*>(gt));

    GDALRasterBand* band = ds->GetRasterBand(1);
    if (band->Fill(0.0) != CE_None)
    {
        GDALClose(ds);
        throw std::runtime_error("lsdtt_vector_rasterise: failed to initialise MEM band");
    }
    return ds;
}

// Opens `vector_path` (read-only), returns its first layer.
// `ds_out` receives ownership of the GDALDataset — caller must GDALClose it.
// Throws on failure.
OGRLayer* open_first_layer(const std::string& vector_path, GDALDataset*& ds_out)
{
    GDALAllRegister();
    ds_out = static_cast<GDALDataset*>(
        GDALOpenEx(vector_path.c_str(),
                   GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr)
    );
    if (!ds_out)
        throw std::runtime_error(
            "lsdtt_vector_rasterise: cannot open vector file: " + vector_path);

    OGRLayer* layer = ds_out->GetLayer(0);
    if (!layer)
    {
        GDALClose(ds_out);
        ds_out = nullptr;
        throw std::runtime_error(
            "lsdtt_vector_rasterise: no layers in: " + vector_path);
    }
    return layer;
}

// Returns true if `t` is (or could contain) a polygon geometry.
bool is_polygon_type(OGRwkbGeometryType t)
{
    t = wkbFlatten(t);
    return t == wkbPolygon ||
           t == wkbMultiPolygon ||
           t == wkbGeometryCollection ||
           t == wkbUnknown;
}

// Returns true if `t` is (or could contain) a line geometry.
bool is_line_type(OGRwkbGeometryType t)
{
    t = wkbFlatten(t);
    return t == wkbLineString ||
           t == wkbMultiLineString ||
           t == wkbGeometryCollection ||
           t == wkbUnknown;
}

// Collects all geometry pointers from `layer` (one per feature), optionally
// transforming with `ct`.  Returns cloned OGRGeometry* objects — caller
// owns them.  Features whose geometry is missing or of the wrong wkb class
// are skipped.  `accept_type` is called to filter individual geometries.
std::vector<OGRGeometry*> collect_geometries(
    OGRLayer*                      layer,
    OGRCoordinateTransformation*   ct,
    bool (*accept_type)(OGRwkbGeometryType)
)
{
    std::vector<OGRGeometry*> out;
    layer->ResetReading();
    OGRFeature* feat;
    while ((feat = layer->GetNextFeature()) != nullptr)
    {
        OGRGeometry* geom = feat->GetGeometryRef();
        if (geom && accept_type(geom->getGeometryType()))
        {
            OGRGeometry* clone = geom->clone();
            if (ct)
            {
                if (clone->transform(ct) != OGRERR_NONE)
                {
                    delete clone;
                    OGRFeature::DestroyFeature(feat);
                    continue;
                }
            }
            out.push_back(clone);
        }
        OGRFeature::DestroyFeature(feat);
    }
    return out;
}

// Frees a vector of cloned OGRGeometry*.
void free_geometries(std::vector<OGRGeometry*>& geoms)
{
    for (auto* g : geoms) delete g;
    geoms.clear();
}

// Converts OGRGeometry* to OGRGeometryH (void*) for the C rasterise API.
OGRGeometryH to_handle(OGRGeometry* g)
{
    return reinterpret_cast<OGRGeometryH>(g);
}

// Calls GDALRasterizeGeometries with burn value 1 into band 1 of `mem_ds`.
// `all_touched` sets the ALL_TOUCHED option.
void burn_geometries(
    GDALDataset*                    mem_ds,
    const std::vector<OGRGeometry*>& geoms,
    bool                            all_touched)
{
    if (geoms.empty()) return;

    std::vector<OGRGeometryH> handles;
    handles.reserve(geoms.size());
    for (auto* g : geoms) handles.push_back(to_handle(g));

    int    band_no   = 1;
    double burn_val  = 1.0;

    char** opts = nullptr;
    if (all_touched)
        opts = CSLSetNameValue(nullptr, "ALL_TOUCHED", "TRUE");

    GDALRasterizeGeometries(
        mem_ds,
        1, &band_no,
        static_cast<int>(handles.size()), handles.data(),
        nullptr, nullptr,
        &burn_val,
        opts,
        nullptr, nullptr
    );

    CSLDestroy(opts);
}

// Reads the single byte band of `mem_ds` into a vector<uint8_t> sized
// nrows * ncols.
std::vector<uint8_t> read_byte_band(GDALDataset* mem_ds, int nrows, int ncols)
{
    std::vector<uint8_t> buf(static_cast<std::size_t>(nrows) * ncols, 0u);
    if (mem_ds->GetRasterBand(1)->RasterIO(
            GF_Read, 0, 0, ncols, nrows,
            buf.data(), ncols, nrows, GDT_Byte,
            0, 0) != CE_None)
    {
        throw std::runtime_error(
            "lsdtt_vector_rasterise: failed to read rasterised band");
    }
    return buf;
}

// Builds the output lsdtt_raster_cpu binary mask.
// `raw[i]` > 0  →  output data[i] = 1.0f (if reference mask is set).
lsdtt_raster_cpu build_output(
    const lsdtt_raster_cpu&     reference_raster,
    const std::vector<uint8_t>& raw)
{
    const lsdtt_raster_meta&     meta = reference_raster.get_meta();
    const std::vector<uint8_t>&  ref_mask = reference_raster.get_mask();
    const std::size_t npix = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;

    std::vector<float>   data(npix, 0.0f);
    std::vector<uint8_t> mask(npix, 0u);

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!ref_mask[i]) continue;
        mask[i] = 1u;
        data[i] = (raw[i] > 0u) ? 1.0f : 0.0f;
    }

    return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// Sets up OGRCoordinateTransformation from `layer_srs` to `raster_srs` if they
// differ.  Returns nullptr when no transform is needed.
// Throws on failure to create the transform.
std::unique_ptr<OGRCoordinateTransformation> make_ct_if_needed(
    OGRSpatialReference* layer_srs,
    const std::string&   raster_wkt)
{
    if (!layer_srs || raster_wkt.empty()) return nullptr;

    OGRSpatialReference raster_srs;
    if (raster_srs.importFromWkt(raster_wkt.c_str()) != OGRERR_NONE) return nullptr;

    if (layer_srs->IsSame(&raster_srs)) return nullptr;

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(layer_srs, &raster_srs);
    if (!ct)
        throw std::runtime_error(
            "lsdtt_vector_rasterise: CRS transform failed — "
            "check that vector and raster have compatible projections");

    return std::unique_ptr<OGRCoordinateTransformation>(ct);
}

} // anonymous namespace

// =============================================================================
// rasterise_polygons
// =============================================================================
lsdtt_raster_cpu rasterise_polygons(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path,
    polygon_rasterise_rule  rule)
{
    const lsdtt_raster_meta& meta = reference_raster.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error("rasterise_polygons: reference raster metadata is invalid");

    // Open vector file
    GDALDataset* vec_ds = nullptr;
    OGRLayer*    layer  = open_first_layer(vector_path, vec_ds);

    // Validate geometry type
    const OGRwkbGeometryType layer_type = layer->GetGeomType();
    if (!is_polygon_type(layer_type))
    {
        GDALClose(vec_ds);
        throw std::runtime_error(
            "rasterise_polygons: first layer of '" + vector_path +
            "' is not a polygon layer (type=" +
            std::string(OGRGeometryTypeToName(layer_type)) + ")");
    }

    // CRS transform
    auto ct = make_ct_if_needed(layer->GetSpatialRef(), meta.crs_wkt);

    // Collect geometries
    auto geoms = collect_geometries(layer, ct.get(), is_polygon_type);
    GDALClose(vec_ds);

    if (geoms.empty())
    {
        // No polygon features — return all-zero mask
        const std::size_t npix = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
        std::vector<float>   data(npix, 0.0f);
        std::vector<uint8_t> mask = reference_raster.get_mask();
        return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
    }

    // -------------------------------------------------------------------------
    // For full_coverage: rasterise at 4× and require all 16 sub-pixels covered.
    // For any_coverage:  rasterise with ALL_TOUCHED=TRUE.
    // For majority_coverage: rasterise with centre-point rule (default).
    // -------------------------------------------------------------------------
    std::vector<uint8_t> raw;

    if (rule == polygon_rasterise_rule::full_coverage)
    {
        constexpr int os = 4;
        const int burn_nrows = meta.n_rows * os;
        const int burn_ncols = meta.n_cols * os;

        double gt[6];
        fill_geotransform(meta, gt);
        gt[1] /= os;
        gt[5] /= os;

        GDALDataset* mem_ds = make_mem_byte_raster(burn_nrows, burn_ncols, gt);
        burn_geometries(mem_ds, geoms, /*all_touched=*/false);
        const auto hi_res = read_byte_band(mem_ds, burn_nrows, burn_ncols);
        GDALClose(mem_ds);

        constexpr int sub = os * os;
        raw.resize(static_cast<std::size_t>(meta.n_rows) * meta.n_cols, 0u);
        for (int r = 0; r < meta.n_rows; ++r)
        {
            for (int c = 0; c < meta.n_cols; ++c)
            {
                int count = 0;
                for (int sr = 0; sr < os; ++sr)
                    for (int sc = 0; sc < os; ++sc)
                        count += hi_res[static_cast<std::size_t>(r * os + sr) * burn_ncols +
                                        (c * os + sc)];
                raw[static_cast<std::size_t>(r) * meta.n_cols + c] =
                    (count == sub) ? 1u : 0u;
            }
        }
    }
    else
    {
        double gt[6];
        fill_geotransform(meta, gt);

        GDALDataset* mem_ds = make_mem_byte_raster(meta.n_rows, meta.n_cols, gt);
        burn_geometries(mem_ds, geoms,
                        rule == polygon_rasterise_rule::any_coverage);
        raw = read_byte_band(mem_ds, meta.n_rows, meta.n_cols);
        GDALClose(mem_ds);
    }

    free_geometries(geoms);
    return build_output(reference_raster, raw);
}

// =============================================================================
// rasterise_polylines
// =============================================================================
lsdtt_raster_cpu rasterise_polylines(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path)
{
    const lsdtt_raster_meta& meta = reference_raster.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error("rasterise_polylines: reference raster metadata is invalid");

    // Open vector file
    GDALDataset* vec_ds = nullptr;
    OGRLayer*    layer  = open_first_layer(vector_path, vec_ds);

    // Validate geometry type
    const OGRwkbGeometryType layer_type = layer->GetGeomType();
    if (!is_line_type(layer_type))
    {
        GDALClose(vec_ds);
        throw std::runtime_error(
            "rasterise_polylines: first layer of '" + vector_path +
            "' is not a line layer (type=" +
            std::string(OGRGeometryTypeToName(layer_type)) + ")");
    }

    // CRS transform
    auto ct = make_ct_if_needed(layer->GetSpatialRef(), meta.crs_wkt);

    // Collect geometries
    auto geoms = collect_geometries(layer, ct.get(), is_line_type);
    GDALClose(vec_ds);

    if (geoms.empty())
    {
        const std::size_t npix = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
        std::vector<float>   data(npix, 0.0f);
        std::vector<uint8_t> mask = reference_raster.get_mask();
        return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
    }

    double gt[6];
    fill_geotransform(meta, gt);

    GDALDataset* mem_ds = make_mem_byte_raster(meta.n_rows, meta.n_cols, gt);
    burn_geometries(mem_ds, geoms, /*all_touched=*/true);
    const auto raw = read_byte_band(mem_ds, meta.n_rows, meta.n_cols);
    GDALClose(mem_ds);

    free_geometries(geoms);
    return build_output(reference_raster, raw);
}

} // namespace lsdtt
