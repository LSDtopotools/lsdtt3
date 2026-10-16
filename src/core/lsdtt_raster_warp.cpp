// lsdtt_raster_warp.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_raster_warp.hpp for full documentation.

#include "lsdtt_raster_warp.hpp"
#include "lsdtt_pixel_type_gdal.hpp"
#include "lsdtt_gdal_init.hpp"

#include "../utils/lsdtt_raster_geometry.hpp"

#include <array>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <gdal.h>
#include <gdal_priv.h>
#include <gdal_utils.h>
#include <ogr_spatialref.h>
#include <cpl_string.h>

namespace lsdtt {

// =============================================================================
// Internal helpers
// =============================================================================

namespace
{

// -----------------------------------------------------------------------------
// warp_request
// -----------------------------------------------------------------------------
struct warp_request
{
    bool          set_target_srs  = false;
    std::string   target_wkt;

    bool          set_extent      = false;
    lsdtt_extent  target_extent;

    bool          set_dimensions  = false;
    int           target_cols     = 0;
    int           target_rows     = 0;

    std::string          resampling_method = "cubic";
    raster_output_format output_format     = raster_output_format::cog;

    double src_nodata = -9999.0;
    double dst_nodata = -9999.0;

    int n_threads   = 0;
    int warp_mem_mb = 512;
};

// -----------------------------------------------------------------------------
// normalise_resampling_method
// -----------------------------------------------------------------------------
std::string normalise_resampling_method(const std::string& s)
{
    if (s == "near"         || s == "nearest")     return "near";
    if (s == "bilinear")                           return "bilinear";
    if (s == "cubic")                              return "cubic";
    if (s == "cubicspline")                        return "cubicspline";
    if (s == "lanczos")                            return "lanczos";
    if (s == "average")                            return "average";
    if (s == "mode")                               return "mode";
    if (s == "max")                                return "max";
    if (s == "min")                                return "min";
    if (s == "med")                                return "med";
    if (s == "q1")                                 return "q1";
    if (s == "q3")                                 return "q3";
    if (s == "sum")                                return "sum";

    std::cerr << "lsdtt_raster_warp: unknown resampling_method '" << s
              << "' - falling back to cubic.\n";
    return "cubic";
}

// -----------------------------------------------------------------------------
// single_line_wkt
// -----------------------------------------------------------------------------
std::string single_line_wkt(const std::string& wkt)
{
    std::string safe = wkt;
    for (char& c : safe)
    {
        if (c == '\n' || c == '\r')
        {
            c = ' ';
        }
    }
    return safe;
}


// -----------------------------------------------------------------------------
// transformed_extent_in_target_crs
//
// Transform the four source-raster corners into the target CRS and return the
// axis-aligned bounding extent of those transformed corners.
// -----------------------------------------------------------------------------
lsdtt_extent transformed_extent_in_target_crs(
    const lsdtt_raster_meta& source_meta,
    const std::string&       target_wkt)
{
    if (target_wkt.empty())
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: target CRS WKT is empty"
        );
    }

    if (crs_equal_semantically(source_meta.crs_wkt, target_wkt))
    {
        return raster_extent(source_meta);
    }

    OGRSpatialReference src_srs;
    OGRSpatialReference tgt_srs;

    if (src_srs.importFromWkt(source_meta.crs_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: cannot parse source CRS WKT"
        );
    }
    if (tgt_srs.importFromWkt(target_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: cannot parse target CRS WKT"
        );
    }

    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    tgt_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&src_srs, &tgt_srs);

    if (!ct)
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: cannot create source->target CRS transform"
        );
    }

    std::array<double, 4> xs = {
        source_meta.x_min,
        source_meta.x_max(),
        source_meta.x_min,
        source_meta.x_max()
    };

    std::array<double, 4> ys = {
        source_meta.y_min,
        source_meta.y_min,
        source_meta.y_max(),
        source_meta.y_max()
    };

    for (int i = 0; i < 4; ++i)
    {
        if (!ct->Transform(1, &xs[i], &ys[i]))
        {
            OGRCoordinateTransformation::DestroyCT(ct);
            throw std::runtime_error(
                "lsdtt_raster_warp: corner transform failed while checking overlap"
            );
        }
    }

    OGRCoordinateTransformation::DestroyCT(ct);

    lsdtt_extent out;
    out.x_min = *std::min_element(xs.begin(), xs.end());
    out.x_max = *std::max_element(xs.begin(), xs.end());
    out.y_min = *std::min_element(ys.begin(), ys.end());
    out.y_max = *std::max_element(ys.begin(), ys.end());
    return out;
}

// -----------------------------------------------------------------------------
// build_gdalwarp_options
// -----------------------------------------------------------------------------
char** build_gdalwarp_options(
    const warp_request& req,
    const std::string&  source_driver_name,
    // Band type of the SOURCE. GDALWarp preserves it unless -ot is given, and
    // -ot is not given here, so this is also the type of the output -- which is
    // what the compression predictor has to match. Passed in rather than read
    // from the file again: the caller already holds the open dataset.
    GDALDataType        source_type)
{
    char** opts = nullptr;

    if (req.set_target_srs)
    {
        opts = CSLAddString(opts, "-t_srs");
        const std::string safe_wkt = single_line_wkt(req.target_wkt);
        opts = CSLAddString(opts, safe_wkt.c_str());
    }

    if (req.set_extent)
    {
        // Full double precision: the default 6 significant digits would turn a
        // UTM northing of 5000202 into 5.0002e+06 (5000200) and shift the grid.
        std::ostringstream x0, y0, x1, y1;
        x0.precision(17);
        y0.precision(17);
        x1.precision(17);
        y1.precision(17);
        x0 << req.target_extent.x_min;
        y0 << req.target_extent.y_min;
        x1 << req.target_extent.x_max;
        y1 << req.target_extent.y_max;

        opts = CSLAddString(opts, "-te");
        opts = CSLAddString(opts, x0.str().c_str());
        opts = CSLAddString(opts, y0.str().c_str());
        opts = CSLAddString(opts, x1.str().c_str());
        opts = CSLAddString(opts, y1.str().c_str());
    }

    if (req.set_dimensions)
    {
        std::ostringstream nx, ny;
        nx << req.target_cols;
        ny << req.target_rows;

        opts = CSLAddString(opts, "-ts");
        opts = CSLAddString(opts, nx.str().c_str());
        opts = CSLAddString(opts, ny.str().c_str());
    }

    opts = CSLAddString(opts, "-r");
    const std::string resampling = normalise_resampling_method(req.resampling_method);
    opts = CSLAddString(opts, resampling.c_str());

    {
        std::ostringstream src_nd, dst_nd;
        src_nd << req.src_nodata;
        dst_nd << req.dst_nodata;

        opts = CSLAddString(opts, "-srcnodata");
        opts = CSLAddString(opts, src_nd.str().c_str());

        opts = CSLAddString(opts, "-dstnodata");
        opts = CSLAddString(opts, dst_nd.str().c_str());
    }

    // Output format
    std::string driver_name = "GTiff";
    if (req.output_format == raster_output_format::cog)
    {
        driver_name = "COG";
    }
    else if (req.output_format == raster_output_format::same_as_input)
    {
        driver_name = source_driver_name.empty() ? "GTiff" : source_driver_name;
    }

    opts = CSLAddString(opts, "-of");
    opts = CSLAddString(opts, driver_name.c_str());

    // See tiff_predictor_for() in lsdtt_pixel_type.hpp: the predictor depends
    // on the band type, and 2 on float data is spec-invalid and reads as NaN in
    // strict TIFF readers even though GDAL round-trips it happily.
    const std::string predictor_co =
        "PREDICTOR=" + std::to_string(
            tiff_predictor_for(pixel_type_from_gdal(source_type)));

    // Creation options
    if (driver_name == "COG")
    {
        opts = CSLAddString(opts, "-co");
        opts = CSLAddString(opts, "COMPRESS=DEFLATE");
        opts = CSLAddString(opts, "-co");
        opts = CSLAddString(opts, predictor_co.c_str());
    }
    else if (driver_name == "GTiff")
    {
        opts = CSLAddString(opts, "-co");
        opts = CSLAddString(opts, "TILED=YES");
        opts = CSLAddString(opts, "-co");
        opts = CSLAddString(opts, "BLOCKXSIZE=512");
        opts = CSLAddString(opts, "-co");
        opts = CSLAddString(opts, "BLOCKYSIZE=512");
        opts = CSLAddString(opts, "-co");
        opts = CSLAddString(opts, "COMPRESS=DEFLATE");
        opts = CSLAddString(opts, "-co");
        opts = CSLAddString(opts, predictor_co.c_str());
    }

    // Multi-threading and memory
    opts = CSLAddString(opts, "-multi");

    {
        std::ostringstream wm;
        wm << req.warp_mem_mb;
        opts = CSLAddString(opts, "-wm");
        opts = CSLAddString(opts, wm.str().c_str());
    }

    {
        const std::string thread_val =
            (req.n_threads > 0) ? std::to_string(req.n_threads) : "ALL_CPUS";
        opts = CSLAddString(opts, "-wo");
        const std::string wo = "NUM_THREADS=" + thread_val;
        opts = CSLAddString(opts, wo.c_str());
    }

    return opts;
}

// -----------------------------------------------------------------------------
// run_gdalwarp
// -----------------------------------------------------------------------------
void run_gdalwarp(
    const std::string& in_path,
    const std::string& out_path,
    const warp_request& req)
{
    lsdtt::gdal_init();

    GDALDatasetH src_ds = GDALOpen(in_path.c_str(), GA_ReadOnly);
    if (!src_ds)
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: cannot open source raster: " + in_path
        );
    }

    std::string source_driver_name = "GTiff";
    if (GDALGetDatasetDriver(src_ds))
    {
        const char* desc = GDALGetDriverShortName(GDALGetDatasetDriver(src_ds));
        if (desc) source_driver_name = desc;
    }

    const GDALDataType source_type =
        (GDALGetRasterCount(src_ds) > 0)
            ? GDALGetRasterDataType(GDALGetRasterBand(src_ds, 1))
            : GDT_Float32;

    char** warp_opts =
        build_gdalwarp_options(req, source_driver_name, source_type);

    int usage_error = CE_None;
    GDALWarpAppOptions* app_opts = GDALWarpAppOptionsNew(warp_opts, nullptr);
    CSLDestroy(warp_opts);

    if (!app_opts)
    {
        GDALClose(src_ds);
        throw std::runtime_error(
            "lsdtt_raster_warp: failed to create GDALWarp options"
        );
    }

    GDALDatasetH src_arr[1] = { src_ds };

    GDALDatasetH dst_ds = GDALWarp(
        out_path.c_str(),
        nullptr,
        1,
        src_arr,
        app_opts,
        &usage_error
    );

    GDALWarpAppOptionsFree(app_opts);
    GDALClose(src_ds);

    if (!dst_ds || usage_error != CE_None)
    {
        if (dst_ds) GDALClose(dst_ds);
        throw std::runtime_error(
            "lsdtt_raster_warp: GDALWarp failed writing " + out_path
        );
    }

    GDALClose(dst_ds);
}

// -----------------------------------------------------------------------------
// load_meta
// -----------------------------------------------------------------------------
lsdtt_raster_meta load_meta(const std::string& path)
{
    lsdtt_data_loader loader;
    loader.open_raster(path);
    return loader.get_meta();
}

// -----------------------------------------------------------------------------
// summarise_result
// -----------------------------------------------------------------------------
lsdtt_warp_result summarise_result(
    const std::string&       out_path,
    const lsdtt_raster_meta& source_meta,
    const lsdtt_raster_meta& reference_meta,
    bool                     matched_reference_grid_requested)
{
    const lsdtt_raster_meta out_meta = load_meta(out_path);

    lsdtt_warp_result result;
    result.wrote_output = true;
    result.out_path     = out_path;

    result.changed_crs =
        !crs_equal_semantically(source_meta.crs_wkt, out_meta.crs_wkt);

    result.changed_grid =
        !same_grid(source_meta, out_meta);

    const lsdtt_extent src_extent = raster_extent(source_meta);
    const lsdtt_extent out_extent = raster_extent(out_meta);

    result.changed_extent =
        !nearly_equal(src_extent.x_min, out_extent.x_min) ||
        !nearly_equal(src_extent.y_min, out_extent.y_min) ||
        !nearly_equal(src_extent.x_max, out_extent.x_max) ||
        !nearly_equal(src_extent.y_max, out_extent.y_max);

    result.matched_reference_grid =
        matched_reference_grid_requested &&
        same_grid(out_meta, reference_meta) &&
        crs_equal_semantically(out_meta.crs_wkt, reference_meta.crs_wkt);

    return result;
}

// -----------------------------------------------------------------------------
// grid_warp_request
//
// The request that places a source exactly on a target grid (CRS, extent,
// rows, cols), shared by warp_raster_to_reference_grid and
// warp_raster_to_grid_in_memory. Checks overlap first when asked.
// -----------------------------------------------------------------------------
warp_request grid_warp_request(
    const lsdtt_raster_meta&  source_meta,
    const lsdtt_raster_meta&  target_grid,
    const lsdtt_warp_options& opts,
    const char*               target_name)
{
    if (opts.error_if_no_overlap)
    {
        const lsdtt_extent source_in_target =
            transformed_extent_in_target_crs(source_meta, target_grid.crs_wkt);
        if (!extents_intersect(source_in_target, raster_extent(target_grid)))
        {
            throw std::runtime_error(
                std::string("lsdtt_raster_warp: source and ") + target_name +
                " rasters do not overlap in the " + target_name + " CRS"
            );
        }
    }

    warp_request req;
    req.set_target_srs    = true;
    req.target_wkt        = target_grid.crs_wkt;
    req.set_extent        = true;
    req.target_extent     = raster_extent(target_grid);
    req.set_dimensions    = true;
    req.target_cols       = target_grid.n_cols;
    req.target_rows       = target_grid.n_rows;
    req.resampling_method = opts.resampling_method;
    req.output_format     = opts.output_format;
    req.src_nodata        = source_meta.no_data_value;
    req.dst_nodata        = source_meta.no_data_value;
    req.n_threads         = opts.n_threads;
    req.warp_mem_mb       = opts.warp_mem_mb;
    return req;
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// crs_equal_semantically
// -----------------------------------------------------------------------------
bool crs_equal_semantically(
    const std::string& a_wkt,
    const std::string& b_wkt)
{
    if (a_wkt.empty() && b_wkt.empty())
    {
        return true;
    }
    if (a_wkt.empty() || b_wkt.empty())
    {
        return false;
    }

    OGRSpatialReference a_srs;
    OGRSpatialReference b_srs;

    if (a_srs.importFromWkt(a_wkt.c_str()) != OGRERR_NONE) return false;
    if (b_srs.importFromWkt(b_wkt.c_str()) != OGRERR_NONE) return false;

    a_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    b_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    return a_srs.IsSame(&b_srs) != 0;
}

// =============================================================================
// Public API
// =============================================================================

// -----------------------------------------------------------------------------
// clip_raster_to_reference_extent
// -----------------------------------------------------------------------------
lsdtt_warp_result clip_raster_to_reference_extent(
    const std::string&        in_path,
    const std::string&        reference_path,
    const std::string&        out_path,
    const lsdtt_warp_options& opts)
{
    const lsdtt_raster_meta source_meta    = load_meta(in_path);
    const lsdtt_raster_meta reference_meta = load_meta(reference_path);

    if (!source_meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: source raster metadata are invalid"
        );
    }
    if (!reference_meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: reference raster metadata are invalid"
        );
    }

    const bool same_crs =
        crs_equal_semantically(source_meta.crs_wkt, reference_meta.crs_wkt);

    if (!same_crs && opts.require_same_crs_for_clip)
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: source and reference rasters do not share the "
            "same CRS and require_same_crs_for_clip is true"
        );
    }

    // Overlap check performed in the clipping CRS:
    //   - source CRS if same_crs
    //   - reference CRS if source will be warped there
    if (opts.error_if_no_overlap)
    {
        lsdtt_extent source_extent_for_check;
        lsdtt_extent reference_extent_for_check;

        if (same_crs)
        {
            source_extent_for_check    = raster_extent(source_meta);
            reference_extent_for_check = raster_extent(reference_meta);
        }
        else
        {
            source_extent_for_check =
                transformed_extent_in_target_crs(source_meta, reference_meta.crs_wkt);
            reference_extent_for_check = raster_extent(reference_meta);
        }

        if (!extents_intersect(source_extent_for_check, reference_extent_for_check))
        {
            throw std::runtime_error(
                "lsdtt_raster_warp: source and reference rasters do not overlap"
            );
        }
    }

    warp_request req;
    req.set_target_srs   = !same_crs;
    req.target_wkt       = reference_meta.crs_wkt;
    req.set_extent       = true;
    req.target_extent    = raster_extent(reference_meta);
    req.set_dimensions   = false;
    req.resampling_method = opts.resampling_method;
    req.output_format     = opts.output_format;
    req.src_nodata        = source_meta.no_data_value;
    req.dst_nodata        = source_meta.no_data_value;
    req.n_threads         = opts.n_threads;
    req.warp_mem_mb       = opts.warp_mem_mb;

    run_gdalwarp(in_path, out_path, req);

    return summarise_result(
        out_path,
        source_meta,
        reference_meta,
        false
    );
}

// -----------------------------------------------------------------------------
// warp_raster_to_reference_grid
// -----------------------------------------------------------------------------
lsdtt_warp_result warp_raster_to_reference_grid(
    const std::string&        in_path,
    const std::string&        reference_path,
    const std::string&        out_path,
    const lsdtt_warp_options& opts)
{
    const lsdtt_raster_meta source_meta    = load_meta(in_path);
    const lsdtt_raster_meta reference_meta = load_meta(reference_path);

    if (!source_meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: source raster metadata are invalid"
        );
    }
    if (!reference_meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt_raster_warp: reference raster metadata are invalid"
        );
    }

    const warp_request req = grid_warp_request(source_meta, reference_meta, opts, "reference");

    run_gdalwarp(in_path, out_path, req);

    return summarise_result(
        out_path,
        source_meta,
        reference_meta,
        true
    );
}

// -----------------------------------------------------------------------------
// warp_raster_to_grid_in_memory
// -----------------------------------------------------------------------------
lsdtt_raster_tile warp_raster_to_grid_in_memory(
    const std::string&        in_path,
    const lsdtt_raster_meta&  target_grid,
    const lsdtt_warp_options& opts)
{
    const lsdtt_raster_meta source_meta = load_meta(in_path);
    if (!source_meta.is_valid())
        throw std::runtime_error("lsdtt_raster_warp: source raster metadata are invalid");
    if (!target_grid.is_valid())
        throw std::runtime_error("lsdtt_raster_warp: target grid description is invalid");
    if (source_meta.crs_wkt.empty() || target_grid.crs_wkt.empty())
        throw std::runtime_error(
            "lsdtt_raster_warp: warping onto a grid needs a CRS on both the source (" + in_path +
            ") and the target grid");

    lsdtt_warp_options o = opts;
    o.output_format = raster_output_format::gtiff;   // no overviews, nothing averaged
    warp_request req = grid_warp_request(source_meta, target_grid, o, "target");

    static int counter = 0;
    std::ostringstream name;
    name << "/vsimem/lsdtt_warp_to_grid_" << reinterpret_cast<std::uintptr_t>(&req) << "_"
         << ++counter << ".tif";
    const std::string mem_path = name.str();

    lsdtt_raster_tile tile;
    try
    {
        run_gdalwarp(in_path, mem_path, req);
        lsdtt_data_loader loader;
        loader.open_raster(mem_path);
        tile = loader.read_full();
    }
    catch (...)
    {
        VSIUnlink(mem_path.c_str());
        throw;
    }
    VSIUnlink(mem_path.c_str());
    return tile;
}

} // namespace lsdtt
