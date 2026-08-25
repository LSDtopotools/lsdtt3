// lsdtt_raster_reproject.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// See lsdtt_raster_reproject.hpp for full documentation.

#include "lsdtt_raster_reproject.hpp"
#include "lsdtt_crs_suggest.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <gdal.h>
#include <gdal_version.h>
#include <gdal_priv.h>
#include <gdal_utils.h>      // GDALWarpAppOptions
#include <gdalwarper.h>      // GDALResampleAlg, GRA_* constants
#include <ogr_spatialref.h>
#include <cpl_string.h>      // CSLAddString / CSLDestroy

namespace lsdtt {

// =============================================================================
// Internal helpers
// =============================================================================

namespace {

// -----------------------------------------------------------------------------
// gdal_resample_string_to_algo
//
// Map a user-facing string to a GDALResampleAlg enum value.
// Falls back to GRA_CubicSpline with a warning if unrecognised.
// -----------------------------------------------------------------------------
static GDALResampleAlg resample_algo(const std::string& s)
{
    if (s == "near"         || s == "nearest")  return GRA_NearestNeighbour;
    if (s == "bilinear")                         return GRA_Bilinear;
    if (s == "cubic")                            return GRA_Cubic;
    if (s == "cubicspline")                      return GRA_CubicSpline;
    if (s == "lanczos")                          return GRA_Lanczos;
    if (s == "average")                          return GRA_Average;
    if (s == "mode")                             return GRA_Mode;
    if (s == "max")                              return GRA_Max;
    if (s == "min")                              return GRA_Min;
    if (s == "med")                              return GRA_Med;
    if (s == "q1")                               return GRA_Q1;
    if (s == "q3")                               return GRA_Q3;
    if (s == "sum")
    {
#if GDAL_VERSION_NUM >= 3010000
        return GRA_Sum;
#else
        throw std::invalid_argument(
            "lsdtt_reproject: 'sum' resampling requires GDAL 3.1 or newer");
#endif
    }

    std::cerr << "lsdtt_reproject: unknown resampling_method '"
              << s << "' - falling back to cubic.\n";
    return GRA_Cubic;
}

// -----------------------------------------------------------------------------
// estimate_projected_pixel_size
//
// Given the source meta (assumed geographic or Web Mercator) and a target WKT,
// estimate what the pixel size will be in the target projection in metres.
//
// Strategy: transform the raster centre + one pixel offset to the target CRS
// and measure the resulting distance.
// -----------------------------------------------------------------------------
static double estimate_projected_pixel_size(
    const lsdtt_raster_meta& meta,
    const lsdtt_crs_suggestion& sug,
    const std::string& target_wkt)
{
    // Parse source and target SRS
    OGRSpatialReference src_srs;
    if (src_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_reproject: cannot parse source CRS for pixel size estimation");
    }

    OGRSpatialReference tgt_srs;
    if (tgt_srs.importFromWkt(target_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_reproject: cannot parse target CRS for pixel size estimation");
    }

    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    tgt_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&src_srs, &tgt_srs);
    if (!ct)
    {
        throw std::runtime_error(
            "lsdtt_reproject: cannot create coordinate transformation "
            "for pixel size estimation");
    }

    // Centre of raster in source coordinates
    double cx = 0.5 * (meta.x_min + meta.x_max());
    double cy = 0.5 * (meta.y_min + meta.y_max());

    // One pixel offset in source coordinates
    double ox = cx + meta.cell_size;
    double oy = cy;

    if (!ct->Transform(1, &cx, &cy) || !ct->Transform(1, &ox, &oy))
    {
        OGRCoordinateTransformation::DestroyCT(ct);
        throw std::runtime_error(
            "lsdtt_reproject: pixel size estimation transform failed");
    }
    OGRCoordinateTransformation::DestroyCT(ct);

    double dx = ox - cx;
    double dy = oy - cy;
    return std::sqrt(dx * dx + dy * dy);
}

// -----------------------------------------------------------------------------
// projection_suffix
// -----------------------------------------------------------------------------
static std::string projection_suffix(CrsRecommendationType type)
{
    switch (type)
    {
        case CrsRecommendationType::polar_north:
        case CrsRecommendationType::polar_south:
            return "_polar";
        case CrsRecommendationType::laea_custom:
            return "_lambert";
        case CrsRecommendationType::utm:
        default:
            return "_utm";
    }
}

// -----------------------------------------------------------------------------
// run_gdalwarp
//
// Wrap GDALWarp() to reproject in_path -> out_path with the given WKT and
// pixel size. Writes a Cloud Optimised GeoTIFF.
// -----------------------------------------------------------------------------
static void run_gdalwarp(
    const std::string& in_path,
    const std::string& out_path,
    const std::string& target_wkt,
    double             pixel_size,
    double             no_data_value,
    GDALResampleAlg    algo,
    int                n_threads  = 0,   // 0 = ALL_CPUS
    int                warp_mem_mb = 512) // working memory per thread in MB
{
    GDALAllRegister();

    // Build WKT as a single-line string safe to pass as a command option
    // (replace embedded newlines just in case)
    std::string safe_wkt = target_wkt;
    for (char& c : safe_wkt)
        if (c == '\n' || c == '\r') c = ' ';

    // Build argv-style options for GDALWarpAppOptions
    // We use the string-array form so we don't need to spawn a subprocess.
    char** warp_opts = nullptr;

    // Target SRS
    warp_opts = CSLAddString(warp_opts, "-t_srs");
    warp_opts = CSLAddString(warp_opts, safe_wkt.c_str());

    // Pixel size
    {
        std::ostringstream ps;
        ps << pixel_size;
        warp_opts = CSLAddString(warp_opts, "-tr");
        warp_opts = CSLAddString(warp_opts, ps.str().c_str());
        warp_opts = CSLAddString(warp_opts, ps.str().c_str());
    }

    // Resampling algorithm
    const char* algo_str = "cubic";
    switch (algo)
    {
        case GRA_NearestNeighbour: algo_str = "near";        break;
        case GRA_Bilinear:         algo_str = "bilinear";    break;
        case GRA_Cubic:            algo_str = "cubic";       break;
        case GRA_CubicSpline:      algo_str = "cubicspline"; break;
        case GRA_Lanczos:          algo_str = "lanczos";     break;
        case GRA_Average:          algo_str = "average";     break;
        case GRA_Mode:             algo_str = "mode";        break;
        default:                   algo_str = "cubic";       break;
    }
    warp_opts = CSLAddString(warp_opts, "-r");
    warp_opts = CSLAddString(warp_opts, algo_str);

    // NoData
    {
        std::ostringstream nd;
        nd << no_data_value;
        warp_opts = CSLAddString(warp_opts, "-dstnodata");
        warp_opts = CSLAddString(warp_opts, nd.str().c_str());
        warp_opts = CSLAddString(warp_opts, "-srcnodata");
        warp_opts = CSLAddString(warp_opts, nd.str().c_str());
    }

    // Output format - COG
    warp_opts = CSLAddString(warp_opts, "-of");
    warp_opts = CSLAddString(warp_opts, "COG");

    // COG creation options
    warp_opts = CSLAddString(warp_opts, "-co");
    warp_opts = CSLAddString(warp_opts, "COMPRESS=DEFLATE");
    warp_opts = CSLAddString(warp_opts, "-co");
    warp_opts = CSLAddString(warp_opts, "RESAMPLING=AVERAGE");

    // -multi: two-thread I/O pipeline (reader + writer run concurrently).
    // This is separate from the warp kernel thread pool below.
    warp_opts = CSLAddString(warp_opts, "-multi");

    // Warp working memory.  Larger values mean fewer tile read passes for
    // large rasters.  512 MB is a safe default; expose via n_reproject_mem_mb
    // if users need to tune it on memory-constrained systems.
    {
        std::ostringstream wm;
        wm << warp_mem_mb;
        warp_opts = CSLAddString(warp_opts, "-wm");
        warp_opts = CSLAddString(warp_opts, wm.str().c_str());
    }

    // Warp kernel thread pool.
    // -wo NUM_THREADS controls how many threads the actual resampling kernel
    // uses.  Without this flag the kernel is single-threaded even with -multi.
    // ALL_CPUS tells GDAL to query the CPU count at runtime.
    // If n_threads > 0 the user has specified an explicit count.
    {
        std::string thread_val =
            (n_threads > 0)
                ? std::to_string(n_threads)
                : "ALL_CPUS";
        warp_opts = CSLAddString(warp_opts, "-wo");
        std::string wo_str = "NUM_THREADS=" + thread_val;
        warp_opts = CSLAddString(warp_opts, wo_str.c_str());
    }

    int usage_error = CE_None;
    GDALWarpAppOptions* app_opts =
        GDALWarpAppOptionsNew(warp_opts, nullptr);
    CSLDestroy(warp_opts);

    if (!app_opts)
    {
        throw std::runtime_error(
            "lsdtt_reproject: GDALWarpAppOptionsNew failed");
    }

    // Open source
    GDALDatasetH src_ds = GDALOpen(in_path.c_str(), GA_ReadOnly);
    if (!src_ds)
    {
        GDALWarpAppOptionsFree(app_opts);
        throw std::runtime_error(
            "lsdtt_reproject: cannot open source raster: " + in_path);
    }

    GDALDatasetH src_arr[1] = { src_ds };

    GDALDatasetH dst_ds = GDALWarp(
        out_path.c_str(),   // output path
        nullptr,            // no pre-existing destination dataset
        1,                  // number of source datasets
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
            "lsdtt_reproject: GDALWarp failed writing " + out_path);
    }

    GDALClose(dst_ds);
}

} // anonymous namespace

// =============================================================================
// Public API
// =============================================================================

// -----------------------------------------------------------------------------
// snap_grid_spacing
// -----------------------------------------------------------------------------
double snap_grid_spacing(double pixel_m)
{
    if (pixel_m >= 20.0)
    {
        // Round to nearest 5 m
        return std::round(pixel_m / 5.0) * 5.0;
    }
    else if (pixel_m >= 5.0)
    {
        // Round to nearest 2 m
        return std::round(pixel_m / 2.0) * 2.0;
    }
    else
    {
        // Round to nearest 1 m
        double snapped = std::round(pixel_m);
        // Guard: rounding can produce 0 for sub-metre inputs.
        // In that case preserve the raw pixel size (do not snap).
        return snapped > 0.0 ? snapped : pixel_m;
    }
}

// -----------------------------------------------------------------------------
// lsdtt_reproject_if_needed
// -----------------------------------------------------------------------------
lsdtt_reproject_result lsdtt_reproject_if_needed(
    const std::string&             in_path,
    const std::string&             write_path,
    const std::string&             write_stem,
    const lsdtt_raster_meta&       meta,
    const lsdtt_reproject_options& opts)
{
    lsdtt_reproject_result result;

    // ------------------------------------------------------------------
    // 1. Get a CRS recommendation
    // ------------------------------------------------------------------
    result.suggestion = lsdtt_suggest_crs(
        meta,
        opts.allow_lambert_projection,
        opts.force_reproject
    );
    lsdtt_print_crs_suggestion(result.suggestion);

    // ------------------------------------------------------------------
    // 2. Already projected  check if we need to force the reprojection
    // ------------------------------------------------------------------
    if (result.suggestion.already_projected && !opts.force_reproject)
    {
        std::cout << "lsdtt_reproject: raster is already in a projected CRS "
                     "with metre units - skipping reprojection.\n";
        result.reprojected = false;
        return result;
    }

    if (result.suggestion.already_projected && opts.force_reproject)
    {
        std::cout << "lsdtt_reproject: force_reproject is true - warping "
                     "projected input to the recommended analysis CRS.\n";
    }


    // ------------------------------------------------------------------
    // 3. enforce_projection == false - bail out
    // ------------------------------------------------------------------
    if (!opts.enforce_projection)
    {
        throw std::runtime_error(
            "lsdtt_reproject: raster is in a geographic or Web Mercator CRS "
            "and enforce_projection is false. "
            "Please reproject your raster before running lsdtt3 analyses.");
    }

    // ------------------------------------------------------------------
    // 4. Choose target WKT
    //    If allow_lambert_projection is true and the raster spans >= 3
    //    UTM zones, prefer the LAEA WKT.
    // ------------------------------------------------------------------
    std::string target_wkt = result.suggestion.recommended_wkt;
    CrsRecommendationType chosen_type = result.suggestion.type;

    if (opts.allow_lambert_projection
        && result.suggestion.warn_multi_zone
        && !result.suggestion.laea_wkt.empty())
    {
        std::cout << "lsdtt_reproject: raster spans "
                  << result.suggestion.utm_zone_count
                  << " UTM zones - using custom Lambert Azimuthal Equal Area "
                     "projection instead of UTM.\n";
        target_wkt   = result.suggestion.laea_wkt;
        chosen_type  = CrsRecommendationType::laea_custom;
    }
    else if (result.suggestion.warn_multi_zone)
    {
        std::cout << "lsdtt_reproject: WARNING - raster spans "
                  << result.suggestion.utm_zone_count
                  << " UTM zones. Consider setting allow_lambert_projection: true "
                     "for a less distorted result.\n";
    }

    result.target_wkt = target_wkt;

    // ------------------------------------------------------------------
    // 5. Determine grid spacing
    // ------------------------------------------------------------------
    double grid_m = -1.0;

    if (opts.force_resampling)
    {
        if (opts.forced_grid_spacing <= 0.0)
        {
            throw std::runtime_error(
                "lsdtt_reproject: force_resampling is true but "
                "forced_grid_spacing is not positive");
        }

        grid_m = opts.forced_grid_spacing;

        std::cout << "lsdtt_reproject: force_resampling is true - using "
                  << "forced_grid_spacing = " << grid_m << " m\n";
    }
    else if (opts.grid_spacing > 0.0)
    {
        grid_m = opts.grid_spacing;

        std::cout << "lsdtt_reproject: using user-supplied grid spacing = "
                  << grid_m << " m\n";
    }
    else
    {
        double raw_pixel_m = estimate_projected_pixel_size(
            meta, result.suggestion, target_wkt);

        grid_m = snap_grid_spacing(raw_pixel_m);

        std::cout << "lsdtt_reproject: estimated projected pixel size = "
                  << raw_pixel_m << " m  -> snapped grid spacing = "
                  << grid_m << " m\n";
    }

    if (grid_m <= 0.0)
    {
        std::cerr << "lsdtt_reproject: WARNING - computed grid spacing is "
                  << grid_m << " m; clamping to 0.01 m.\n";
        grid_m = 0.01;
    }

    result.grid_spacing_used = grid_m;

    // ------------------------------------------------------------------
    // 6. Build output path
    // ------------------------------------------------------------------
    const std::string suffix   = projection_suffix(chosen_type);
    const std::string out_file = write_path + write_stem + suffix + ".tif";

    std::cout << "lsdtt_reproject: writing reprojected raster to "
              << out_file << "\n";

    // ------------------------------------------------------------------
    // 7. Run gdalwarp
    // ------------------------------------------------------------------
    GDALResampleAlg algo = resample_algo(opts.resampling_method);

    run_gdalwarp(
        in_path,
        out_file,
        target_wkt,
        grid_m,
        meta.no_data_value,
        algo,
        opts.n_threads,
        opts.warp_mem_mb
    );

    std::cout << "lsdtt_reproject: reprojection complete.\n";

    result.reprojected = true;
    result.out_path    = out_file;
    return result;
}

} // namespace lsdtt
