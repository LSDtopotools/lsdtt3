// lsdtt_vector_rasterise_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_vector_rasterise_cpu.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "../core/lsdtt_gdal_init.hpp"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <gdal.h>
#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt {
namespace {

// Thrown when the vector and raster CRSs are both valid but GDAL cannot build
// a transform between them.
const char* const k_crs_transform_failure =
    "lsdtt_vector_rasterise: CRS transform failed - check that vector and raster have compatible projections";


struct macrostrat_units_field_info
{
    const char*            field_name = "";
    macrostrat_field_kind  kind       = macrostrat_field_kind::integer;
};

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
// with the given geotransform. Initialised to all-zero.
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
// `ds_out` receives ownership of the GDALDataset - caller must GDALClose it.
// Throws on failure.
OGRLayer* open_first_layer(const std::string& vector_path, GDALDataset*& ds_out)
{
    lsdtt::gdal_init();

    ds_out = static_cast<GDALDataset*>(
        GDALOpenEx(vector_path.c_str(),
                   GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr)
    );

    if (!ds_out)
    {
        throw std::runtime_error(
            "lsdtt_vector_rasterise: cannot open vector file: " + vector_path);
    }

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

// Opens `vector_path` (read-only), returns a named layer.
// `ds_out` receives ownership of the GDALDataset - caller must GDALClose it.
// Throws on failure.
OGRLayer* open_named_layer(
    const std::string& vector_path,
    const std::string& layer_name,
    GDALDataset*&      ds_out)
{
    lsdtt::gdal_init();

    ds_out = static_cast<GDALDataset*>(
        GDALOpenEx(vector_path.c_str(),
                   GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr)
    );

    if (!ds_out)
    {
        throw std::runtime_error(
            "lsdtt_vector_rasterise: cannot open vector file: " + vector_path);
    }

    OGRLayer* layer = ds_out->GetLayerByName(layer_name.c_str());
    if (!layer)
    {
        GDALClose(ds_out);
        ds_out = nullptr;
        throw std::runtime_error(
            "lsdtt_vector_rasterise: layer '" + layer_name +
            "' not found in: " + vector_path);
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

// Collects geometry pointers from `layer`, optionally transforming with `ct`.
// Returns cloned OGRGeometry* objects - caller owns them. Features whose
// geometry is missing or of the wrong wkb class are skipped.
// max_features == 0 means collect all; otherwise stop after max_features.
std::vector<OGRGeometry*> collect_geometries(
    OGRLayer*                      layer,
    OGRCoordinateTransformation*   ct,
    bool (*accept_type)(OGRwkbGeometryType),
    int                            max_features = 0)
{
    std::vector<OGRGeometry*> out;

    layer->ResetReading();
    OGRFeature* feat;
    while ((feat = layer->GetNextFeature()) != nullptr)
    {
        if (max_features > 0 && static_cast<int>(out.size()) >= max_features)
        {
            OGRFeature::DestroyFeature(feat);
            break;
        }

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
    GDALDataset*                     mem_ds,
    const std::vector<OGRGeometry*>& geoms,
    bool                             all_touched)
{
    if (geoms.empty()) return;

    std::vector<OGRGeometryH> handles;
    handles.reserve(geoms.size());
    for (auto* g : geoms) handles.push_back(to_handle(g));

    int    band_no  = 1;
    double burn_val = 1.0;

    char** opts = nullptr;
    if (all_touched)
        opts = CSLSetNameValue(nullptr, "ALL_TOUCHED", "TRUE");

    const CPLErr err = GDALRasterizeGeometries(
        mem_ds,
        1, &band_no,
        static_cast<int>(handles.size()), handles.data(),
        nullptr, nullptr,
        &burn_val,
        opts,
        nullptr, nullptr
    );

    CSLDestroy(opts);

    if (err != CE_None)
    {
        throw std::runtime_error(
            "lsdtt_vector_rasterise: GDALRasterizeGeometries failed");
    }
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
// `raw[i]` > 0  -> output data[i] = 1.0f (if reference mask is set).
lsdtt_raster_cpu build_output(
    const lsdtt_raster_cpu&     reference_raster,
    const std::vector<uint8_t>& raw)
{
    const lsdtt_raster_meta&    meta     = reference_raster.get_meta();
    const std::vector<uint8_t>& ref_mask = reference_raster.get_mask();
    const std::size_t           npix =
        static_cast<std::size_t>(meta.n_rows) * meta.n_cols;

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

macrostrat_units_field_info get_macrostrat_units_field_info(
    macrostrat_units_field field)
{
    switch (field)
    {
        case macrostrat_units_field::map_id:
            return {"map_id", macrostrat_field_kind::integer};

        case macrostrat_units_field::legeln_id:
            return {"legeln_id", macrostrat_field_kind::integer};

        case macrostrat_units_field::best_age_top:
            return {"best_age_top", macrostrat_field_kind::real};

        case macrostrat_units_field::best_age_bottom:
            return {"best_age_bottom", macrostrat_field_kind::real};

        case macrostrat_units_field::lith:
            return {"lith", macrostrat_field_kind::string};

        case macrostrat_units_field::t_int_id:
            return {"t_int_id", macrostrat_field_kind::integer};

        case macrostrat_units_field::t_int:
            return {"t_int", macrostrat_field_kind::string};

        case macrostrat_units_field::b_int_id:
            return {"b_int_id", macrostrat_field_kind::integer};

        case macrostrat_units_field::b_int:
            return {"b_int", macrostrat_field_kind::string};
    }

    throw std::runtime_error("rasterise_macrostrat_geology: unknown field");
}

int find_field_index_or_throw(OGRLayer* layer, const char* field_name)
{
    OGRFeatureDefn* defn = layer->GetLayerDefn();
    if (!defn)
    {
        throw std::runtime_error(
            "rasterise_macrostrat_geology: layer has no valid schema");
    }

    const int idx = defn->GetFieldIndex(field_name);
    if (idx < 0)
    {
        throw std::runtime_error(
            "rasterise_macrostrat_geology: field '" +
            std::string(field_name) + "' not found in Macrostrat units layer");
    }

    return idx;
}

float float_from_int64_checked(std::int64_t v, const std::string& field_name)
{
    static constexpr std::int64_t k_float_exact_int_limit = 16777216LL;

    if (v < -k_float_exact_int_limit || v > k_float_exact_int_limit)
    {
        throw std::runtime_error(
            "rasterise_macrostrat_geology: integer field '" + field_name +
            "' contains a value outside the exact float32 range");
    }

    return static_cast<float>(v);
}

bool get_numeric_field_value(
    const OGRFeature*         feat,
    int                       field_index,
    macrostrat_field_kind     kind,
    const std::string&        field_name,
    float&                    out_value)
{
    if (!feat->IsFieldSetAndNotNull(field_index))
        return false;

    if (kind == macrostrat_field_kind::integer)
    {
        const std::int64_t v =
            static_cast<std::int64_t>(feat->GetFieldAsInteger64(field_index));
        out_value = float_from_int64_checked(v, field_name);
        return true;
    }

    if (kind == macrostrat_field_kind::real)
    {
        out_value = static_cast<float>(feat->GetFieldAsDouble(field_index));
        return true;
    }

    return false;
}

bool get_string_field_value(
    const OGRFeature*  feat,
    int                field_index,
    std::string&       out_value)
{
    if (!feat->IsFieldSetAndNotNull(field_index))
        return false;

    out_value = feat->GetFieldAsString(field_index);
    return true;
}

bool raster_window_from_envelope(
    const lsdtt_raster_meta& meta,
    const OGREnvelope&       env,
    int&                     row0,
    int&                     row1,
    int&                     col0,
    int&                     col1)
{
    const double raster_minx = meta.x_min;
    const double raster_maxx = meta.x_max();
    const double raster_miny = meta.y_min;
    const double raster_maxy = meta.y_max();
    const double cs          = meta.cell_size;

    const double minx = std::max(env.MinX, raster_minx);
    const double maxx = std::min(env.MaxX, raster_maxx);
    const double miny = std::max(env.MinY, raster_miny);
    const double maxy = std::min(env.MaxY, raster_maxy);

    if (!(minx < maxx && miny < maxy))
        return false;

    col0 = std::max(
        0,
        std::min(meta.n_cols - 1,
                 static_cast<int>(std::floor((minx - raster_minx) / cs)))
    );

    col1 = std::max(
        0,
        std::min(meta.n_cols - 1,
                 static_cast<int>(std::ceil((maxx - raster_minx) / cs)) - 1)
    );

    row0 = std::max(
        0,
        std::min(meta.n_rows - 1,
                 static_cast<int>(std::floor((raster_maxy - maxy) / cs)))
    );

    row1 = std::max(
        0,
        std::min(meta.n_rows - 1,
                 static_cast<int>(std::ceil((raster_maxy - miny) / cs)) - 1)
    );

    return row0 <= row1 && col0 <= col1;
}

void fill_window_geotransform(
    const lsdtt_raster_meta& meta,
    int                      row0,
    int                      col0,
    int                      oversample,
    double                   gt[6])
{
    const double raster_top = meta.y_min + static_cast<double>(meta.n_rows) * meta.cell_size;

    gt[0] = meta.x_min + static_cast<double>(col0) * meta.cell_size;
    gt[1] = meta.cell_size / static_cast<double>(oversample);
    gt[2] = 0.0;
    gt[3] = raster_top - static_cast<double>(row0) * meta.cell_size;
    gt[4] = 0.0;
    gt[5] = -meta.cell_size / static_cast<double>(oversample);
}

std::vector<macrostrat_lookup_entry> build_lookup_table(
    OGRLayer*          layer,
    int                field_index)
{
    std::set<std::string> unique_values;

    layer->ResetReading();
    OGRFeature* feat = nullptr;
    while ((feat = layer->GetNextFeature()) != nullptr)
    {
        OGRGeometry* geom = feat->GetGeometryRef();
        if (!geom || !is_polygon_type(geom->getGeometryType()))
        {
            OGRFeature::DestroyFeature(feat);
            continue;
        }

        std::string v;
        if (get_string_field_value(feat, field_index, v))
        {
            unique_values.insert(v);
        }

        OGRFeature::DestroyFeature(feat);
    }

    std::vector<macrostrat_lookup_entry> table;
    table.reserve(unique_values.size());

    std::uint32_t code = 1;
    for (const std::string& v : unique_values)
    {
        table.push_back(macrostrat_lookup_entry{code, v});
        ++code;
    }

    return table;
}

std::unordered_map<std::string, std::uint32_t> build_lookup_map(
    const std::vector<macrostrat_lookup_entry>& table)
{
    std::unordered_map<std::string, std::uint32_t> out;
    out.reserve(table.size());

    for (const auto& entry : table)
    {
        out.emplace(entry.value, entry.code);
    }

    return out;
}

void csv_write_escaped(std::ofstream& ofs, const std::string& s)
{
    bool needs_quotes = false;
    for (char c : s)
    {
        if (c == ',' || c == '"' || c == '\n' || c == '\r')
        {
            needs_quotes = true;
            break;
        }
    }

    if (!needs_quotes)
    {
        ofs << s;
        return;
    }

    ofs << '"';
    for (char c : s)
    {
        if (c == '"') ofs << '"';
        ofs << c;
    }
    ofs << '"';
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
    {
        throw std::runtime_error(
            "rasterise_polygons: reference raster metadata is invalid");
    }

    GDALDataset* vec_ds = nullptr;
    OGRLayer*    layer  = open_first_layer(vector_path, vec_ds);

    const OGRwkbGeometryType layer_type = layer->GetGeomType();
    if (!is_polygon_type(layer_type))
    {
        GDALClose(vec_ds);
        throw std::runtime_error(
            "rasterise_polygons: first layer of '" + vector_path +
            "' is not a polygon layer (type=" +
            std::string(OGRGeometryTypeToName(layer_type)) + ")");
    }

    auto ct = make_layer_to_raster_transform(layer->GetSpatialRef(), meta.crs_wkt, k_crs_transform_failure);

    auto geoms = collect_geometries(layer, ct.get(), is_polygon_type);
    GDALClose(vec_ds);

    if (geoms.empty())
    {
        const std::size_t npix =
            static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
        std::vector<float>   data(npix, 0.0f);
        std::vector<uint8_t> mask = reference_raster.get_mask();
        return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
    }

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
                {
                    for (int sc = 0; sc < os; ++sc)
                    {
                        count += hi_res[
                            static_cast<std::size_t>(r * os + sr) * burn_ncols +
                            (c * os + sc)
                        ];
                    }
                }

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
// rasterise_polygons (max_features overload)
// =============================================================================
lsdtt_raster_cpu rasterise_polygons(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path,
    polygon_rasterise_rule  rule,
    int                     max_features)
{
    const lsdtt_raster_meta& meta = reference_raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "rasterise_polygons: reference raster metadata is invalid");
    }

    GDALDataset* vec_ds = nullptr;
    OGRLayer*    layer  = open_first_layer(vector_path, vec_ds);

    const OGRwkbGeometryType layer_type = layer->GetGeomType();
    if (!is_polygon_type(layer_type))
    {
        GDALClose(vec_ds);
        throw std::runtime_error(
            "rasterise_polygons: first layer of '" + vector_path +
            "' is not a polygon layer (type=" +
            std::string(OGRGeometryTypeToName(layer_type)) + ")");
    }

    auto ct = make_layer_to_raster_transform(layer->GetSpatialRef(), meta.crs_wkt, k_crs_transform_failure);

    auto geoms = collect_geometries(layer, ct.get(), is_polygon_type,
                                    max_features);
    GDALClose(vec_ds);

    if (geoms.empty())
    {
        const std::size_t npix =
            static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
        std::vector<float>   data(npix, 0.0f);
        std::vector<uint8_t> mask = reference_raster.get_mask();
        return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
    }

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
                        count += hi_res[
                            static_cast<std::size_t>(r * os + sr) * burn_ncols +
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
// get_first_polygon_bbox_in_raster_crs
// =============================================================================
polygon_bbox get_first_polygon_bbox_in_raster_crs(
    const std::string& vector_path,
    const std::string& raster_crs_wkt)
{
    GDALDataset* vec_ds = nullptr;
    OGRLayer*    layer  = open_first_layer(vector_path, vec_ds);

    struct ds_guard_t {
        GDALDataset* ds;
        ~ds_guard_t() { if (ds) GDALClose(ds); }
    } ds_guard{vec_ds};

    auto ct = make_layer_to_raster_transform(layer->GetSpatialRef(), raster_crs_wkt, k_crs_transform_failure);

    layer->ResetReading();
    OGRFeature* feat = layer->GetNextFeature();

    if (!feat)
    {
        throw std::runtime_error(
            "get_first_polygon_bbox_in_raster_crs: no features in: " +
            vector_path);
    }

    struct feat_guard_t {
        OGRFeature* f;
        ~feat_guard_t() { if (f) OGRFeature::DestroyFeature(f); }
    } feat_guard{feat};

    OGRGeometry* geom = feat->GetGeometryRef();
    if (!geom)
    {
        throw std::runtime_error(
            "get_first_polygon_bbox_in_raster_crs: first feature has no "
            "geometry in: " + vector_path);
    }

    // Work on a clone if we need to transform.
    OGRGeometry* work = ct ? geom->clone() : geom;
    struct work_guard_t {
        OGRGeometry* g;
        bool owned;
        ~work_guard_t() { if (owned && g) delete g; }
    } work_guard{work, ct != nullptr};

    if (ct)
    {
        if (work->transform(ct.get()) != OGRERR_NONE)
        {
            throw std::runtime_error(
                "get_first_polygon_bbox_in_raster_crs: CRS transform failed "
                "for: " + vector_path);
        }
    }

    OGREnvelope env;
    work->getEnvelope(&env);

    return polygon_bbox{env.MinX, env.MinY, env.MaxX, env.MaxY};
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
    {
        throw std::runtime_error(
            "rasterise_polylines: reference raster metadata is invalid");
    }

    GDALDataset* vec_ds = nullptr;
    OGRLayer*    layer  = open_first_layer(vector_path, vec_ds);

    const OGRwkbGeometryType layer_type = layer->GetGeomType();
    if (!is_line_type(layer_type))
    {
        GDALClose(vec_ds);
        throw std::runtime_error(
            "rasterise_polylines: first layer of '" + vector_path +
            "' is not a line layer (type=" +
            std::string(OGRGeometryTypeToName(layer_type)) + ")");
    }

    auto ct = make_layer_to_raster_transform(layer->GetSpatialRef(), meta.crs_wkt, k_crs_transform_failure);

    auto geoms = collect_geometries(layer, ct.get(), is_line_type);
    GDALClose(vec_ds);

    if (geoms.empty())
    {
        const std::size_t npix =
            static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
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

// =============================================================================
// rasterise_macrostrat_geology
// =============================================================================
macrostrat_geology_rasterise_result rasterise_macrostrat_geology(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      mvt_path,
    macrostrat_units_field  field,
    int                     oversample)
{
    const lsdtt_raster_meta& meta = reference_raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "rasterise_macrostrat_geology: reference raster metadata is invalid");
    }

    if (oversample <= 0)
    {
        throw std::runtime_error(
            "rasterise_macrostrat_geology: oversample must be > 0");
    }

    GDALDataset* vec_ds = nullptr;
    OGRLayer*    layer  = open_named_layer(mvt_path, "units", vec_ds);

    const OGRwkbGeometryType layer_type = layer->GetGeomType();
    if (!is_polygon_type(layer_type))
    {
        GDALClose(vec_ds);
        throw std::runtime_error(
            "rasterise_macrostrat_geology: layer 'units' in '" + mvt_path +
            "' is not a polygon layer (type=" +
            std::string(OGRGeometryTypeToName(layer_type)) + ")");
    }

    set_layer_spatial_filter_from_raster_bbox(layer, meta);

    const macrostrat_units_field_info info = get_macrostrat_units_field_info(field);
    const int field_index = find_field_index_or_throw(layer, info.field_name);

    std::vector<macrostrat_lookup_entry> lookup_table;
    std::unordered_map<std::string, std::uint32_t> lookup_map;

    if (info.kind == macrostrat_field_kind::string)
    {
        lookup_table = build_lookup_table(layer, field_index);
        lookup_map   = build_lookup_map(lookup_table);
    }

    auto ct = make_layer_to_raster_transform(layer->GetSpatialRef(), meta.crs_wkt, k_crs_transform_failure);

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) * meta.n_cols;

    std::vector<float>   data(npix, 0.0f);
    std::vector<uint8_t> mask(npix, 0u);
    std::vector<int>     best_count(npix, -1);

    const std::vector<uint8_t>& ref_mask = reference_raster.get_mask();
    for (std::size_t i = 0; i < npix; ++i)
    {
        mask[i] = ref_mask[i] ? 1u : 0u;
    }

    layer->ResetReading();
    OGRFeature* feat = nullptr;
    while ((feat = layer->GetNextFeature()) != nullptr)
    {
        OGRGeometry* geom = feat->GetGeometryRef();
        if (!geom || !is_polygon_type(geom->getGeometryType()))
        {
            OGRFeature::DestroyFeature(feat);
            continue;
        }

        float raster_value = 0.0f;

        if (info.kind == macrostrat_field_kind::string)
        {
            std::string s;
            if (!get_string_field_value(feat, field_index, s))
            {
                OGRFeature::DestroyFeature(feat);
                continue;
            }

            auto it = lookup_map.find(s);
            if (it == lookup_map.end())
            {
                OGRFeature::DestroyFeature(feat);
                continue;
            }

            raster_value = static_cast<float>(it->second);
        }
        else
        {
            if (!get_numeric_field_value(
                    feat, field_index, info.kind, info.field_name, raster_value))
            {
                OGRFeature::DestroyFeature(feat);
                continue;
            }
        }

        OGRGeometry* clone = geom->clone();
        if (ct)
        {
            if (clone->transform(ct.get()) != OGRERR_NONE)
            {
                delete clone;
                OGRFeature::DestroyFeature(feat);
                continue;
            }
        }

        OGREnvelope env;
        clone->getEnvelope(&env);

        int row0, row1, col0, col1;
        if (!raster_window_from_envelope(meta, env, row0, row1, col0, col1))
        {
            delete clone;
            OGRFeature::DestroyFeature(feat);
            continue;
        }

        const int coarse_rows = row1 - row0 + 1;
        const int coarse_cols = col1 - col0 + 1;
        const int burn_rows   = coarse_rows * oversample;
        const int burn_cols   = coarse_cols * oversample;

        double gt[6];
        fill_window_geotransform(meta, row0, col0, oversample, gt);

        GDALDataset* mem_ds = make_mem_byte_raster(burn_rows, burn_cols, gt);

        {
            std::vector<OGRGeometry*> one_geom;
            one_geom.push_back(clone);
            burn_geometries(mem_ds, one_geom, /*all_touched=*/false);
        }

        const std::vector<uint8_t> hi_res = read_byte_band(mem_ds, burn_rows, burn_cols);
        GDALClose(mem_ds);
        delete clone;

        for (int wr = 0; wr < coarse_rows; ++wr)
        {
            const int global_r = row0 + wr;

            for (int wc = 0; wc < coarse_cols; ++wc)
            {
                const int global_c = col0 + wc;
                const std::size_t idx =
                    static_cast<std::size_t>(global_r) * meta.n_cols + global_c;

                if (!ref_mask[idx])
                    continue;

                int count = 0;
                for (int sr = 0; sr < oversample; ++sr)
                {
                    const int hi_r = wr * oversample + sr;
                    const std::size_t row_base =
                        static_cast<std::size_t>(hi_r) * burn_cols;

                    for (int sc = 0; sc < oversample; ++sc)
                    {
                        const int hi_c = wc * oversample + sc;
                        count += static_cast<int>(hi_res[row_base + hi_c]);
                    }
                }

                if (count <= 0)
                    continue;

                if (count > best_count[idx])
                {
                    best_count[idx] = count;
                    data[idx]       = raster_value;
                }
            }
        }

        OGRFeature::DestroyFeature(feat);
    }

    layer->SetSpatialFilter(nullptr);
    GDALClose(vec_ds);

    macrostrat_geology_rasterise_result result;
    result.raster       = lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
    result.field        = field;
    result.kind         = info.kind;
    result.lookup_table = std::move(lookup_table);

    return result;
}

// =============================================================================
// write_macrostrat_lookup_csv
// =============================================================================
void write_macrostrat_lookup_csv(
    const std::vector<macrostrat_lookup_entry>& lookup_table,
    const std::string&                          out_csv_path)
{
    std::ofstream ofs(out_csv_path);
    if (!ofs.is_open())
    {
        throw std::runtime_error(
            "write_macrostrat_lookup_csv: could not open CSV for writing: " +
            out_csv_path);
    }

    ofs << "code,value\n";
    for (const auto& entry : lookup_table)
    {
        ofs << entry.code << ',';
        csv_write_escaped(ofs, entry.value);
        ofs << '\n';
    }
}

// =============================================================================
// macrostrat_rasterise_faults
// =============================================================================
lsdtt_raster_cpu macrostrat_rasterise_faults(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      mvt_path)
{
    (void)reference_raster;
    (void)mvt_path;

    throw std::runtime_error(
        "macrostrat_rasterise_faults: not yet implemented");
}

} // namespace lsdtt
