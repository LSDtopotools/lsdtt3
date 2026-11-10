// lsdtt_raster_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_raster_fetch.hpp for full documentation.

#include "lsdtt_raster_fetch.hpp"
#include "lsdtt_raster_fetch_detail.hpp"
#include "lsdtt_gdal_init.hpp"
#include "lsdtt_secret.hpp"
#include "lsdtt_macrostrat_fetch.hpp"   // slippy_tiles_for_bbox (XYZ tiles)

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cpl_conv.h>        // CPLSetThreadLocalConfigOption, CPLSPrintf
#include <cpl_error.h>       // CPLPushErrorHandler / CPLQuietErrorHandler
#include <cpl_http.h>        // CPLHTTPEnabled
#include <cpl_string.h>      // CSLAddString / CSLDestroy
#include <cpl_vsi.h>         // VSIStatL, VSIUnlink
#include <gdal.h>
#include <gdal_utils.h>      // GDALTranslate / GDALWarp
#include <ogr_spatialref.h>
#include "lsdtt_http_fetch.hpp"
#include <cpl_json.h>        // CPLJSONDocument (STAC + SAS responses)
#include <cpl_multiproc.h>   // CPLSleep (WCS request delay)
#include <ogr_api.h>         // OGR_L_* / OGR_G_* vector access
#include <ogr_srs_api.h>     // OSR* / OCT* coordinate transforms

namespace lsdtt {

namespace {

constexpr const char* k_cop30_url_prefix =
    "/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_";

// -----------------------------------------------------------------------------
// gdal_dataset_guard / gdal_string_list_guard
//
// Minimal RAII so early returns cannot leak GDAL handles. The fetch path has
// many failure exits (offline, bad bbox, warp failure), which is exactly the
// situation where manual cleanup goes wrong.
// -----------------------------------------------------------------------------
class gdal_dataset_guard
{
public:
    explicit gdal_dataset_guard(GDALDatasetH ds = nullptr) : ds_(ds) {}
    ~gdal_dataset_guard() { reset(); }

    gdal_dataset_guard(const gdal_dataset_guard&)            = delete;
    gdal_dataset_guard& operator=(const gdal_dataset_guard&) = delete;

    void reset(GDALDatasetH ds = nullptr)
    {
        if (ds_ && ds_ != ds) GDALClose(ds_);
        ds_ = ds;
    }

    GDALDatasetH get() const noexcept { return ds_; }
    explicit operator bool() const noexcept { return ds_ != nullptr; }

private:
    GDALDatasetH ds_;
};


class gdal_string_list_guard
{
public:
    gdal_string_list_guard() = default;
    ~gdal_string_list_guard() { CSLDestroy(list_); }

    gdal_string_list_guard(const gdal_string_list_guard&)            = delete;
    gdal_string_list_guard& operator=(const gdal_string_list_guard&) = delete;

    void add(const char* value) { list_ = CSLAddString(list_, value); }
    void add(const std::string& value) { add(value.c_str()); }

    char** get() const noexcept { return list_; }

private:
    char** list_ = nullptr;
};


// -----------------------------------------------------------------------------
// quiet_gdal_errors
//
// GDAL writes failures to stderr by default. A probe that is *expected* to fail
// when offline should not spray CURL errors over test output, so silence the
// handler for the probe's lifetime only.
//
// This is also where credential masking will hook in once keyed providers land
// (docs/developer_notes/testing_and_benchmarking.md rule 4): GDAL error messages include the failing URL, so a key in
// a query string would otherwise leak to stderr on any failure.
// -----------------------------------------------------------------------------
class quiet_gdal_errors
{
public:
    quiet_gdal_errors() { CPLPushErrorHandler(CPLQuietErrorHandler); }
    ~quiet_gdal_errors() { CPLPopErrorHandler(); }

    quiet_gdal_errors(const quiet_gdal_errors&)            = delete;
    quiet_gdal_errors& operator=(const quiet_gdal_errors&) = delete;
};


// -----------------------------------------------------------------------------
// http_config_guard
//
// Sets /vsicurl/ tuning for the duration of a fetch and restores the previous
// values. Thread-local so a fetch on one thread cannot disturb another.
//
// GDAL_DISABLE_READDIR_ON_OPEN suppresses a bucket listing on every open, which
// is a large latency win against S3. The timeouts turn an unreachable host into
// a prompt bounded failure instead of a hung caller.
// -----------------------------------------------------------------------------
class http_config_guard
{
public:
    http_config_guard(int timeout_s, int connect_timeout_s, int max_retry)
    {
        set("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR");
        // Deliberately NOT setting CPL_VSIL_CURL_ALLOWED_EXTENSIONS: Planetary
        // Computer asset hrefs carry no file extension at all
        // (".../Copernicus_DSM_COG_10_N46_00_E010_00_DEM"), so an extension
        // allow-list silently blocks every one of them. READDIR_ON_OPEN above
        // already provides the bulk of the latency win.
        set("GDAL_HTTP_TIMEOUT",        std::to_string(timeout_s).c_str());
        set("GDAL_HTTP_CONNECTTIMEOUT", std::to_string(connect_timeout_s).c_str());
        set("GDAL_HTTP_MAX_RETRY",      std::to_string(max_retry).c_str());
        set("GDAL_HTTP_RETRY_DELAY",    "2");
    }

    ~http_config_guard()
    {
        for (const auto& kv : saved_)
        {
            CPLSetThreadLocalConfigOption(
                kv.first.c_str(), kv.second.empty() ? nullptr : kv.second.c_str());
        }
    }

    http_config_guard(const http_config_guard&)            = delete;
    http_config_guard& operator=(const http_config_guard&) = delete;

private:
    void set(const char* key, const char* value)
    {
        const char* previous = CPLGetConfigOption(key, nullptr);
        saved_.emplace_back(key, previous ? previous : "");
        CPLSetThreadLocalConfigOption(key, value);
    }

    std::vector<std::pair<std::string, std::string>> saved_;
};


std::string format_double(double value, int precision)
{
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.*f", precision, value);
    return std::string(buffer);
}


bool file_exists(const std::string& path)
{
    VSIStatBufL stat_buffer;
    return VSIStatL(path.c_str(), &stat_buffer) == 0;
}


std::string cache_key_for(const lsdtt_fetch_options& o)
{
    return fetch_detail::fetch_cache_key(o);
}


// -----------------------------------------------------------------------------
// opentopography_url
//
// Builds the OpenTopography global DEM request. The API key travels in the
// query string, which is what the service requires -- so EVERY message derived
// from this URL must be passed through mask_secret() before it is shown.
// -----------------------------------------------------------------------------
std::string opentopography_url(const lsdtt_fetch_options& o)
{
    std::ostringstream url;
    url << "https://portal.opentopography.org/API/globaldem"
        << "?demtype="      << (o.ot_dem_type.empty() ? "COP30" : o.ot_dem_type)
        << "&south="        << format_double(o.south, 10)
        << "&north="        << format_double(o.north, 10)
        << "&west="         << format_double(o.west,  10)
        << "&east="         << format_double(o.east,  10)
        << "&outputFormat=GTiff"
        << "&API_Key="      << o.api_key;
    return url.str();
}


// -----------------------------------------------------------------------------
// download_to_vsimem
//
// GET a URL into an in-memory file. Returns false on any failure, with a
// message that has already had the secret masked out.
//
// CPLHTTPFetch is used rather than /vsicurl/ because the OpenTopography
// endpoint is a dynamic request, not a Cloud-Optimised GeoTIFF: there is no
// range-request advantage, and a failed request returns an error document that
// must be inspected rather than handed to GDAL as if it were a raster.
// -----------------------------------------------------------------------------
// CPLHTTPDestroyResult frees pabyData; ours points into a capped_download,
// so detach it first.
void destroy_detached(CPLHTTPResult* response)
{
    response->pabyData = nullptr;
    response->nDataLen = 0;
    CPLHTTPDestroyResult(response);
}

// The body of a download, collected by a write callback so the byte limit
// (max_download_megabytes) can stop the transfer as soon as it is passed.
struct capped_download
{
    std::vector<unsigned char> data;
    double max_bytes = 0.0;   // 0: no limit
    bool   over_limit = false;
};

size_t capped_download_write(void* buffer, size_t size, size_t count, void* arg)
{
    auto* d = static_cast<capped_download*>(arg);
    const size_t n = size * count;
    if (d->max_bytes > 0.0 &&
        static_cast<double>(d->data.size() + n) > d->max_bytes)
    {
        d->over_limit = true;
        return 0;   // tells curl to abort the transfer
    }
    const unsigned char* p = static_cast<const unsigned char*>(buffer);
    d->data.insert(d->data.end(), p, p + n);
    return n;
}

bool download_to_vsimem(const std::string& url,
                        const std::string& secret,
                        const std::string& vsimem_path,
                        int timeout_seconds,
                        int max_retry,
                        double max_download_megabytes,
                        std::string& error)
{
    gdal_string_list_guard http_opts;
    http_opts.add("TIMEOUT="   + std::to_string(timeout_seconds));
    http_opts.add("MAX_RETRY=" + std::to_string(max_retry));
    http_opts.add("RETRY_DELAY=2");

    capped_download body;
    body.max_bytes = max_download_megabytes * 1.0e6;

    CPLHTTPResult* response = nullptr;
    {
        quiet_gdal_errors quiet;   // GDAL would otherwise print the keyed URL
        response = CPLHTTPFetchEx(url.c_str(), http_opts.get(), nullptr, nullptr,
                                  capped_download_write, &body);
    }

    if (body.over_limit)
    {
        if (response != nullptr) destroy_detached(response);
        error = "the download passed max_download_megabytes (" +
                format_double(max_download_megabytes, 1) + " MB) and was stopped; "
                "raise max_download_megabytes, fetch a smaller area, or use a "
                "coarser product";
        return false;
    }
    if (response != nullptr)
    {
        response->pabyData = body.data.empty() ? nullptr : body.data.data();
        response->nDataLen = static_cast<int>(body.data.size());
    }

    if (response == nullptr)
    {
        error = "no response from the OpenTopography API";
        return false;
    }

    if (response->nStatus != 0)
    {
        error = "HTTP request failed: " +
                mask_secret(response->pszErrBuf ? response->pszErrBuf
                                                : "unknown transport error",
                            secret);
        destroy_detached(response);
        return false;
    }

    if (response->nDataLen <= 0 || response->pabyData == nullptr)
    {
        error = "the OpenTopography API returned an empty response";
        destroy_detached(response);
        return false;
    }

    // On a bad key, an out-of-coverage box, or a quota breach the API returns
    // 200 with a short text/plain or text/html body. Handing that to GDAL would
    // produce an unhelpful "not a supported file format"; report the actual
    // message instead, masked.
    const std::string content_type =
        response->pszContentType ? response->pszContentType : "";

    if (content_type.find("text") != std::string::npos ||
        content_type.find("json") != std::string::npos ||
        content_type.find("xml")  != std::string::npos)
    {
        std::string snippet(reinterpret_cast<const char*>(response->pabyData),
                            static_cast<std::size_t>(response->nDataLen));
        if (snippet.size() > 400) snippet = snippet.substr(0, 400) + "...";
        error = "the OpenTopography API returned an error: " +
                mask_secret(snippet, secret);
        destroy_detached(response);
        return false;
    }

    VSILFILE* handle = VSIFileFromMemBuffer(
        vsimem_path.c_str(), response->pabyData,
        static_cast<vsi_l_offset>(response->nDataLen),
        FALSE);   // FALSE: /vsimem does not take ownership of the buffer

    if (handle == nullptr)
    {
        error = "could not stage the downloaded raster in memory";
        destroy_detached(response);
        return false;
    }
    VSIFCloseL(handle);

    // The buffer must outlive the /vsimem file, so copy before freeing the
    // HTTP result. Reading it back is cheap next to the download itself.
    gdal_dataset_guard staged;
    {
        quiet_gdal_errors quiet;
        staged.reset(GDALOpen(vsimem_path.c_str(), GA_ReadOnly));
    }
    if (!staged)
    {
        VSIUnlink(vsimem_path.c_str());
        destroy_detached(response);
        error = "the downloaded data is not a readable raster";
        return false;
    }

    const std::string owned_path = vsimem_path + ".owned.tif";
    gdal_dataset_guard owned;
    {
        quiet_gdal_errors quiet;
        GDALDriverH gtiff = GDALGetDriverByName("GTiff");
        if (gtiff != nullptr)
        {
            owned.reset(GDALCreateCopy(gtiff, owned_path.c_str(), staged.get(),
                                       FALSE, nullptr, nullptr, nullptr));
        }
    }
    staged.reset();
    VSIUnlink(vsimem_path.c_str());
    destroy_detached(response);

    if (!owned)
    {
        error = "could not copy the downloaded raster into memory";
        return false;
    }
    owned.reset();

    // Hand the caller the self-owned copy.
    VSIRename(owned_path.c_str(), vsimem_path.c_str());
    return true;
}


lsdtt_fetch_result make_error(const std::string& message);
lsdtt_fetch_result make_skip(const std::string& message);
lsdtt_fetch_result warp_and_write_sources(const std::vector<GDALDatasetH>& sources,
                                          const lsdtt_fetch_options& options,
                                          lsdtt_fetch_result result,
                                          const std::string& out_path,
                                          const fetch_detail::fetch_extent* explicit_extent);


// -----------------------------------------------------------------------------
// http_get_body / http_post_body
//
// Small wrappers returning the response body as a string. Errors are reported
// with the URL query masked, because Planetary Computer signed hrefs carry a
// SAS token in the query string.
// -----------------------------------------------------------------------------
bool http_request_body(const std::string& url,
                       const std::string& post_body,   // empty = GET
                       int timeout_seconds,
                       int max_retry,
                       std::string& body,
                       std::string& error)
{
    lsdtt_http_request_options o;
    o.timeout_seconds = timeout_seconds;
    o.max_retry = max_retry;
    return lsdtt::http_request_body(url, post_body, o, body, error);
}


// -----------------------------------------------------------------------------
// pc_sign_href
//
// Planetary Computer asset URLs are not publicly readable -- an unsigned GET
// returns HTTP 409. The signing endpoint needs no account and no key, and
// returns a time-limited SAS token appended as a query string.
//
// The returned URL is therefore a short-lived credential: never log it
// unmasked.
// -----------------------------------------------------------------------------
bool pc_sign_href(const std::string& href, int timeout_seconds, int max_retry,
                  std::string& signed_href, std::string& error)
{
    const std::string sign_url =
        "https://planetarycomputer.microsoft.com/api/sas/v1/sign?href=" + href;

    std::string body;
    if (!http_request_body(sign_url, "", timeout_seconds, max_retry, body, error))
    {
        return false;
    }

    CPLJSONDocument document;
    if (!document.LoadMemory(body))
    {
        error = "could not parse the signing response as JSON";
        return false;
    }

    signed_href = document.GetRoot().GetString("href", "");
    if (signed_href.empty())
    {
        error = "the signing response contained no 'href' field";
        return false;
    }
    return true;
}


// -----------------------------------------------------------------------------
// pc_collection_token
//
// One SAS token for a whole Planetary Computer collection (anonymous, no key),
// appended to every asset href. One request instead of one per file, which
// the per-href signing endpoint sometimes refuses when called in quick
// succession. The token is a short-lived credential: never log it.
// -----------------------------------------------------------------------------
bool pc_collection_token(const std::string& collection, int timeout_seconds, int max_retry,
                         std::string& token, std::string& error)
{
    std::string body;
    if (!http_request_body("https://planetarycomputer.microsoft.com/api/sas/v1/token/" +
                           collection, "", timeout_seconds, max_retry, body, error))
    {
        return false;
    }
    CPLJSONDocument document;
    if (!document.LoadMemory(body))
    {
        error = "could not parse the token response as JSON";
        return false;
    }
    token = document.GetRoot().GetString("token", "");
    if (token.empty())
    {
        error = "the token response contained no 'token' field";
        return false;
    }
    return true;
}


// -----------------------------------------------------------------------------
// pc_search_assets
//
// STAC search for items intersecting the bounding box, returning the signed
// href of the chosen asset for each item.
// -----------------------------------------------------------------------------
bool pc_search_assets(const lsdtt_fetch_options& o,
                      std::vector<std::string>& signed_hrefs,
                      std::vector<std::string>& item_ids,
                      std::vector<std::string>& warnings,
                      bool& is_error,
                      std::string& error)
{
    is_error = false;
    const std::string collection =
        o.pc_collection.empty() ? "cop-dem-glo-30" : o.pc_collection;
    const std::string asset_key =
        o.pc_asset.empty() ? "data" : o.pc_asset;

    std::string body;
    if (!http_request_body(
            "https://planetarycomputer.microsoft.com/api/stac/v1/search",
            fetch_detail::pc_stac_search_body(o), o.http_timeout_seconds, o.http_max_retry,
            body, error))
    {
        return false;
    }

    CPLJSONDocument document;
    if (!document.LoadMemory(body))
    {
        error = "could not parse the STAC search response as JSON";
        return false;
    }

    CPLJSONArray features = document.GetRoot().GetArray("features");
    if (!features.IsValid() || features.Size() == 0)
    {
        error = "the STAC search returned no items for this bounding box in "
                "collection '" + collection + "'";
        return false;
    }

    // gsd of every item, for the source_resolution_metres filter and the
    // mixed-resolution warning. With the filter at 0 every item is kept, in
    // the order returned: the original behaviour.
    const double wanted_gsd = o.source_resolution_metres;
    std::vector<double> seen_gsd;
    for (int i = 0; i < features.Size(); ++i)
    {
        CPLJSONObject feature = features[i];
        const std::string id = feature.GetString("id", "");
        const double gsd = feature.GetDouble("properties/gsd", 0.0);
        if (std::find(seen_gsd.begin(), seen_gsd.end(), gsd) == seen_gsd.end())
            seen_gsd.push_back(gsd);
        if (wanted_gsd > 0.0 && std::abs(gsd - wanted_gsd) > 1e-6 * wanted_gsd)
            continue;

        CPLJSONObject asset = feature.GetObj("assets/" + asset_key);
        if (!asset.IsValid())
        {
            error = "item '" + id + "' has no asset named '" + asset_key + "'";
            return false;
        }

        const std::string href = asset.GetString("href", "");
        if (href.empty())
        {
            error = "asset '" + asset_key + "' of item '" + id + "' has no href";
            return false;
        }

        std::string signed_href;
        if (!pc_sign_href(href, o.http_timeout_seconds, o.http_max_retry,
                          signed_href, error))
        {
            return false;
        }

        signed_hrefs.push_back(signed_href);
        item_ids.push_back(id.empty() ? "(unnamed item)" : id);
    }

    std::sort(seen_gsd.begin(), seen_gsd.end());
    std::string seen_list;
    for (double g : seen_gsd)
        seen_list += (seen_list.empty() ? "" : ", ") + format_double(g, 1);
    if (wanted_gsd > 0.0 && signed_hrefs.empty())
    {
        is_error = true;
        error = "no item of collection '" + collection + "' over this area has gsd " +
                format_double(wanted_gsd, 1) + " m (gsd values found: " + seen_list +
                "); change source_resolution_metres";
        return false;
    }
    if (!(wanted_gsd > 0.0) && seen_gsd.size() > 1)
    {
        warnings.push_back("collection '" + collection + "' returned items of "
                           "several resolutions (gsd " + seen_list + " m) and all "
                           "of them were mosaicked; set source_resolution_metres "
                           "to use one");
    }
    return true;
}

// -----------------------------------------------------------------------------
// frame_source_kernel_scale
//
// The resampling-kernel scale of a warp into a frame: 1 when the frame cell is
// no coarser than the source's ground pixel (the geometric mean of its two
// sides, in metres at the frame centre's latitude for a geographic source),
// otherwise source / frame cell.
// -----------------------------------------------------------------------------
double frame_source_kernel_scale(GDALDatasetH ds, const lsdtt_oblique_frame& f)
{
    double gt[6] = {0.0, 1.0, 0.0, 0.0, 0.0, -1.0};
    if (!ds || GDALGetGeoTransform(ds, gt) != CE_None || !(f.cell_size > 0.0)) return 1.0;
    double sx = std::abs(gt[1]);
    double sy = std::abs(gt[5]);
    OGRSpatialReference srs;
    const char* wkt = GDALGetProjectionRef(ds);
    if (wkt && *wkt && srs.importFromWkt(wkt) == OGRERR_NONE && srs.IsGeographic())
    {
        constexpr double metres_per_degree = 111320.0;
        sx *= metres_per_degree * std::cos(f.centre.latitude * 3.14159265358979323846 / 180.0);
        sy *= metres_per_degree;
    }
    const double pixel = std::sqrt(sx * sy);
    if (!(pixel > 0.0)) return 1.0;
    return std::min(1.0, pixel / f.cell_size);
}

// -----------------------------------------------------------------------------
// warp_and_write
//
// Steps 3 and 4 of the pipeline, shared by every source: reproject to the
// target CRS at a fixed grid spacing, then optionally crop an exact square
// from the centre. Sources differ only in how they obtain the input dataset.
// -----------------------------------------------------------------------------
lsdtt_fetch_result warp_and_write_sources(const std::vector<GDALDatasetH>& sources,
                                          const lsdtt_fetch_options& options,
                                          lsdtt_fetch_result result,
                                          const std::string& out_path,
                                          const fetch_detail::fetch_extent* explicit_extent)
{
    // ---------------------------------------------------------------------
    // Step 3: warp to a metric CRS at a fixed grid spacing
    // ---------------------------------------------------------------------
    const std::string target_label = options.use_frame
                                   ? std::string("the oblique frame")
                                   : "EPSG:" + std::to_string(result.epsg);
    gdal_string_list_guard warp_args;
    for (const std::string& arg :
         fetch_detail::fetch_warp_arguments(options, result.epsg, explicit_extent))
    {
        warp_args.add(arg);
    }
    if (options.use_frame && !sources.empty())
    {
        // A rotated grid: fix the resampling kernel scale and use the exact
        // transform, or GDAL widens the kernel per warp chunk and smooths the
        // DEM by a chunk-dependent amount (see warp_request in
        // lsdtt_raster_warp.cpp). Scale 1 unless the frame cell is coarser
        // than the source's ground pixel; with several mosaic tiles, the first
        // tile's pixel size is used.
        warp_args.add("-et");     warp_args.add("0");
        const double scale = frame_source_kernel_scale(sources.front(), options.frame);
        warp_args.add("-wo");     warp_args.add("XSCALE=" + format_double(scale, 17));
        warp_args.add("-wo");     warp_args.add("YSCALE=" + format_double(scale, 17));
    }

    GDALWarpAppOptions* warp_opts = GDALWarpAppOptionsNew(warp_args.get(), nullptr);
    if (!warp_opts)
    {
        return make_error("lsdtt::fetch_dem: GDALWarpAppOptionsNew failed");
    }

    // When no centre crop is requested the warp writes straight to out_path;
    // otherwise it lands in memory and step 4 does the final write.
    const bool need_crop = (options.size_px > 0);
    const std::string warp_target =
        need_crop ? std::string("/vsimem/lsdtt_fetch_warp.tif") : out_path;

    gdal_dataset_guard warped;
    {
        int usage_error = CE_None;
        std::vector<GDALDatasetH> inputs = sources;
        warped.reset(GDALWarp(warp_target.c_str(), nullptr,
                              static_cast<int>(inputs.size()), inputs.data(),
                              warp_opts, &usage_error));
        GDALWarpAppOptionsFree(warp_opts);
        if (!warped || usage_error != CE_None)
        {
            return make_error("lsdtt::fetch_dem: warp to " + target_label + " failed");
        }
    }

    // ---------------------------------------------------------------------
    // Step 4: optional exact centre crop
    // ---------------------------------------------------------------------
    if (need_crop)
    {
        const int warped_cols = GDALGetRasterXSize(warped.get());
        const int warped_rows = GDALGetRasterYSize(warped.get());

        if (warped_cols < options.size_px || warped_rows < options.size_px)
        {
            warped.reset();
            VSIUnlink("/vsimem/lsdtt_fetch_warp.tif");
            return make_error(
                "lsdtt::fetch_dem: warped raster is " +
                std::to_string(warped_cols) + "x" + std::to_string(warped_rows) +
                ", smaller than the requested " + std::to_string(options.size_px) +
                "x" + std::to_string(options.size_px) + " centre crop -- widen "
                "the bounding box or reduce size_px");
        }

        const int x_off = (warped_cols - options.size_px) / 2;
        const int y_off = (warped_rows - options.size_px) / 2;

        gdal_string_list_guard crop_args;
        crop_args.add("-srcwin");
        crop_args.add(std::to_string(x_off));
        crop_args.add(std::to_string(y_off));
        crop_args.add(std::to_string(options.size_px));
        crop_args.add(std::to_string(options.size_px));
        crop_args.add("-a_nodata"); crop_args.add(format_double(options.no_data_value, 6));
        crop_args.add("-co");       crop_args.add("COMPRESS=DEFLATE");
        crop_args.add("-co");       crop_args.add("TILED=YES");

        GDALTranslateOptions* crop_opts =
            GDALTranslateOptionsNew(crop_args.get(), nullptr);
        if (!crop_opts)
        {
            warped.reset();
            VSIUnlink("/vsimem/lsdtt_fetch_warp.tif");
            return make_error("lsdtt::fetch_dem: GDALTranslateOptionsNew failed");
        }

        int usage_error = CE_None;
        gdal_dataset_guard cropped(
            GDALTranslate(out_path.c_str(), warped.get(), crop_opts, &usage_error));
        GDALTranslateOptionsFree(crop_opts);
        warped.reset();
        VSIUnlink("/vsimem/lsdtt_fetch_warp.tif");

        if (!cropped || usage_error != CE_None)
        {
            return make_error("lsdtt::fetch_dem: centre crop failed writing " +
                              out_path);
        }

        result.n_cols = GDALGetRasterXSize(cropped.get());
        result.n_rows = GDALGetRasterYSize(cropped.get());
    }
    else
    {
        result.n_cols = GDALGetRasterXSize(warped.get());
        result.n_rows = GDALGetRasterYSize(warped.get());
        warped.reset();
    }

    // Caching is the caller's concern: warp_and_write only produces the raster.
    result.ok      = true;
    result.path    = out_path;
    result.message = "fetched " + std::to_string(result.n_cols) + "x" +
                     std::to_string(result.n_rows) + " raster to " + target_label;
    return result;
}


// The single-source form every original source uses.
lsdtt_fetch_result warp_and_write(gdal_dataset_guard& source_ds,
                                  const lsdtt_fetch_options& options,
                                  lsdtt_fetch_result result,
                                  const std::string& out_path)
{
    return warp_and_write_sources({source_ds.get()}, options, result, out_path, nullptr);
}


lsdtt_fetch_result make_error(const std::string& message)
{
    lsdtt_fetch_result result;
    result.ok      = false;
    result.skipped = false;
    result.message = message;
    return result;
}


lsdtt_fetch_result make_skip(const std::string& message)
{
    lsdtt_fetch_result result;
    result.ok      = false;
    result.skipped = true;
    result.message = message;
    return result;
}


bool validate_options(const lsdtt_fetch_options& o, std::string& error)
{
    if (!(o.west < o.east))
    {
        error = "west must be less than east (got west=" +
                format_double(o.west, 6) + ", east=" + format_double(o.east, 6) + ")";
        return false;
    }
    if (!(o.south < o.north))
    {
        error = "south must be less than north (got south=" +
                format_double(o.south, 6) + ", north=" + format_double(o.north, 6) + ")";
        return false;
    }
    if (o.west < -180.0 || o.east > 180.0 || o.south < -90.0 || o.north > 90.0)
    {
        error = "bounding box must lie within WGS84 limits "
                "(-180..180 longitude, -90..90 latitude)";
        return false;
    }
    if (!(o.grid_spacing > 0.0))
    {
        error = "grid_spacing must be > 0";
        return false;
    }
    if (o.size_px < 0)
    {
        error = "size_px must be >= 0";
        return false;
    }
    if (o.resampling.empty())
    {
        error = "resampling must not be empty";
        return false;
    }
    if (o.use_frame)
    {
        if (o.size_px != 0)
        {
            error = "size_px cannot be used with a frame (the frame sets the extent)";
            return false;
        }
        if (!(o.frame.cell_size > 0.0) ||
            std::abs(o.frame.cell_size - o.grid_spacing) > 1e-9 * o.grid_spacing)
        {
            error = "the frame's cell size (" + format_double(o.frame.cell_size, 6) +
                    ") must equal grid_spacing (" + format_double(o.grid_spacing, 6) + ")";
            return false;
        }
    }
    return true;
}

} // anonymous namespace

// =============================================================================
// fetch_detail: pure internals, pinned by the tests
// =============================================================================
namespace fetch_detail {

// Cache key: everything that changes the pixels. Omitting any of these would
// let a cached raster be reused for a request it does not answer.
std::string fetch_cache_key(const lsdtt_fetch_options& o)
{
    std::ostringstream key;
    key << dem_source_name(o.source)
        << '_' << format_double(o.west,  6)
        << '_' << format_double(o.south, 6)
        << '_' << format_double(o.east,  6)
        << '_' << format_double(o.north, 6)
        << "_epsg" << o.target_epsg
        << "_res"  << format_double(o.grid_spacing, 4)
        << '_' << o.resampling
        << "_px" << o.size_px
        << "_nd" << format_double(o.no_data_value, 4);

    // Which dataset was asked for changes the pixels just as much as the bbox
    // does, so it has to be in the key: otherwise a COP30 fetch and a later
    // SRTMGL1 fetch over the same area would collide and the second would
    // silently return the first. Added per-source so keys for the single-dataset
    // cop30_aws source are unaffected.
    if (o.source == lsdtt_dem_source::opentopography)
    {
        key << '_' << (o.ot_dem_type.empty() ? "COP30" : o.ot_dem_type);
    }
    else if (o.source == lsdtt_dem_source::planetary_computer)
    {
        key << '_' << o.pc_collection << '_' << o.pc_asset;
        // Only when filtering, so the default keys are unchanged.
        if (o.source_resolution_metres > 0.0)
            key << "_gsd" << format_double(o.source_resolution_metres, 3);
    }
    else if (!dem_source_info(o.source).legacy_path)
    {
        const lsdtt_dem_source_info& info = dem_source_info(o.source);
        const double res = o.source_resolution_metres > 0.0
                         ? o.source_resolution_metres : info.default_resolution_metres;
        const std::string model = (o.surface_model_choice.empty() ||
                                   o.surface_model_choice == "source_default")
                                ? info.default_surface_model : o.surface_model_choice;
        key << "_src" << format_double(res, 3) << '_' << model;
        if (!o.source_collection_choice.empty() && o.source_collection_choice != "latest")
            key << '_' << o.source_collection_choice;
    }

    // Frame mode: the frame defines the grid. Added only in frame mode so
    // plain fetch keys are unaffected.
    if (o.use_frame)
    {
        key << "_frame_" << format_double(o.frame.centre.longitude, 9)
            << '_' << format_double(o.frame.centre.latitude, 9)
            << '_' << format_double(o.frame.azimuth_degrees, 9)
            << '_' << format_double(o.frame.x_min, 3) << '_' << format_double(o.frame.x_max, 3)
            << '_' << format_double(o.frame.y_min, 3) << '_' << format_double(o.frame.y_max, 3);
    }

    // Keep it filesystem-safe.
    std::string out = key.str();
    for (char& c : out)
    {
        if (c == '.' || c == '-' || c == ' ') c = '_';
    }
    return out + ".tif";
}



std::string pc_stac_search_body(const lsdtt_fetch_options& o)
{
    const std::string collection =
        o.pc_collection.empty() ? "cop-dem-glo-30" : o.pc_collection;
    std::ostringstream request;
    request << "{\"collections\":[\"" << collection << "\"],"
            << "\"bbox\":[" << format_double(o.west,  10) << ','
                            << format_double(o.south, 10) << ','
                            << format_double(o.east,  10) << ','
                            << format_double(o.north, 10) << "],"
            << "\"limit\":100}";
    return request.str();
}


namespace {

// Envelope in `to_epsg` of a box in `from_epsg` with every edge densified.
fetch_extent densified_envelope(const fetch_extent& box, int from_epsg, int to_epsg)
{
    OGRSpatialReferenceH from = OSRNewSpatialReference(nullptr);
    OGRSpatialReferenceH to   = OSRNewSpatialReference(nullptr);
    const OGRErr e1 = OSRImportFromEPSG(from, from_epsg);
    const OGRErr e2 = OSRImportFromEPSG(to, to_epsg);
    OSRSetAxisMappingStrategy(from, OAMS_TRADITIONAL_GIS_ORDER);
    OSRSetAxisMappingStrategy(to, OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformationH ct = nullptr;
    if (e1 == OGRERR_NONE && e2 == OGRERR_NONE)
    {
        quiet_gdal_errors quiet;
        ct = OCTNewCoordinateTransformation(from, to);
    }
    OSRDestroySpatialReference(from);
    OSRDestroySpatialReference(to);
    if (ct == nullptr)
        throw std::runtime_error("cannot transform EPSG:" + std::to_string(from_epsg) +
                                 " to EPSG:" + std::to_string(to_epsg));

    const int per_side = 64;
    std::vector<double> xs, ys;
    for (int i = 0; i <= per_side; ++i)
    {
        const double t = static_cast<double>(i) / per_side;
        const double x = box.x_min + t * (box.x_max - box.x_min);
        const double y = box.y_min + t * (box.y_max - box.y_min);
        xs.insert(xs.end(), {x, x, box.x_min, box.x_max});
        ys.insert(ys.end(), {box.y_min, box.y_max, y, y});
    }
    std::vector<int> ok(xs.size(), 0);
    OCTTransformEx(ct, static_cast<int>(xs.size()), xs.data(), ys.data(), nullptr,
                   ok.data());
    OCTDestroyCoordinateTransformation(ct);
    fetch_extent out{1e300, 1e300, -1e300, -1e300};
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        if (!ok[i]) continue;
        out.x_min = std::min(out.x_min, xs[i]); out.x_max = std::max(out.x_max, xs[i]);
        out.y_min = std::min(out.y_min, ys[i]); out.y_max = std::max(out.y_max, ys[i]);
    }
    if (!(out.x_max > out.x_min && out.y_max > out.y_min))
        throw std::runtime_error("the area does not transform to EPSG:" +
                                 std::to_string(to_epsg));
    return out;
}

std::string url_number(double v)
{
    std::string text = format_double(v, 6);
    while (!text.empty() && text.back() == '0') text.pop_back();
    if (!text.empty() && text.back() == '.') text.pop_back();
    return text;
}

std::string xml_text_of(const std::string& xml, const std::string& tag)
{
    // <tag ...>text</tag>, with or without a namespace prefix.
    std::size_t pos = 0;
    while ((pos = xml.find(tag, pos)) != std::string::npos)
    {
        const bool opening = pos > 0 && (xml[pos - 1] == '<' || xml[pos - 1] == ':');
        const std::size_t after = pos + tag.size();
        if (opening && after < xml.size() && (xml[after] == '>' || xml[after] == ' '))
        {
            std::size_t lt = pos;
            while (lt > 0 && xml[lt] != '<') --lt;
            if (xml[lt] == '<' && xml[lt + 1] != '/')
            {
                const std::size_t start = xml.find('>', after);
                const std::size_t end = xml.find("</", start);
                if (start != std::string::npos && end != std::string::npos)
                {
                    std::string text = xml.substr(start + 1, end - start - 1);
                    const std::size_t a = text.find_first_not_of(" \t\r\n");
                    const std::size_t b = text.find_last_not_of(" \t\r\n");
                    return a == std::string::npos ? std::string() : text.substr(a, b - a + 1);
                }
            }
        }
        pos = after;
    }
    return "";
}

} // namespace


namespace {

// The output area (the snapped target extent, or the bbox in frame mode) as
// an envelope in native_epsg.
fetch_extent output_area_in(const lsdtt_fetch_options& o, int native_epsg, int target_epsg)
{
    if (o.use_frame || target_epsg <= 0)
        return densified_envelope({o.west, o.south, o.east, o.north}, 4326, native_epsg);
    const fetch_extent target =
        target_extent_for_bbox(o.west, o.south, o.east, o.north, target_epsg, o.grid_spacing);
    return densified_envelope(target, target_epsg, native_epsg);
}

} // namespace


fetch_extent wcs_native_extent(const lsdtt_fetch_options&   o,
                               const lsdtt_dem_source_info& info,
                               int                          target_epsg)
{
    const fetch_extent area = output_area_in(o, info.wcs_native_epsg, target_epsg);
    // Pad by the resampling kernel: two output cells, at least two native
    // pixels, plus one.
    const double p = info.wcs_pixel_metres;
    const double pad = p * (std::ceil(2.0 * std::max(1.0, o.grid_spacing / p)) + 1.0);
    fetch_extent out;
    out.x_min = info.wcs_origin_x + std::floor((area.x_min - pad - info.wcs_origin_x) / p) * p;
    out.x_max = info.wcs_origin_x + std::ceil((area.x_max + pad - info.wcs_origin_x) / p) * p;
    out.y_max = info.wcs_origin_y - std::floor((info.wcs_origin_y - area.y_max - pad) / p) * p;
    out.y_min = info.wcs_origin_y - std::ceil((info.wcs_origin_y - area.y_min + pad) / p) * p;
    return out;
}


std::string nztopo50_sheet(double x, double y)
{
    static const std::string letters = "ABCDEFGHJKLMNPQRSTUVWXYZ";   // no I, O
    const double col = std::floor((x - 1012000.0) / 24000.0) + 1.0;
    const double row = std::floor((6234000.0 - y) / 36000.0) + 16.0;   // AS = 0*24 + 16
    if (col < 1.0 || col > 99.0 || row < 0.0 || row >= 24.0 * 24.0) return "";
    const int r = static_cast<int>(row);
    const int c = static_cast<int>(col);
    std::string out;
    out += letters[static_cast<std::size_t>(r / 24)];
    out += letters[static_cast<std::size_t>(r % 24)];
    out += (c < 10 ? "0" : "") + std::to_string(c);
    return out;
}


std::string os_grid_reference(double easting, double northing, int digits_per_axis)
{
    if (easting < 0.0 || northing < 0.0 || easting >= 700000.0 || northing >= 1300000.0 ||
        digits_per_axis < 0 || digits_per_axis > 5)
        return "";
    static const std::string letters = "ABCDEFGHJKLMNOPQRSTUVWXYZ";   // no I
    const int e100k = static_cast<int>(easting / 100000.0);
    const int n100k = static_cast<int>(northing / 100000.0);
    const int l1 = (19 - n100k) - (19 - n100k) % 5 + (e100k + 10) / 5;
    const int l2 = (19 - n100k) * 5 % 25 + e100k % 5;
    std::string out;
    out += letters[static_cast<std::size_t>(l1)];
    out += letters[static_cast<std::size_t>(l2)];
    if (digits_per_axis > 0)
    {
        const double unit = std::pow(10.0, 5 - digits_per_axis);
        const long long e = static_cast<long long>(std::floor(std::fmod(easting, 100000.0) / unit));
        const long long n = static_cast<long long>(std::floor(std::fmod(northing, 100000.0) / unit));
        char buffer[16];
        std::snprintf(buffer, sizeof(buffer), "%0*lld%0*lld", digits_per_axis, e,
                      digits_per_axis, n);
        out += buffer;
    }
    return out;
}


std::string os_grid_quadrant(double easting, double northing)
{
    const std::string square = os_grid_reference(easting, northing, 1);
    if (square.empty()) return "";
    const bool north = std::fmod(northing, 10000.0) >= 5000.0;
    const bool east  = std::fmod(easting, 10000.0) >= 5000.0;
    return square + (north ? "N" : "S") + (east ? "E" : "W");
}


std::vector<std::string> grid_keys_for_extent(
    const lsdtt_dem_source_info::grid_collection& c, const fetch_extent& extent,
    const std::string& model)
{
    std::vector<std::string> keys;
    if (!(c.tile_width > 0.0 && c.tile_height > 0.0)) return keys;
    const long long i0 = static_cast<long long>(std::floor((extent.x_min - c.origin_x) / c.tile_width));
    const long long i1 = static_cast<long long>(std::ceil((extent.x_max - c.origin_x) / c.tile_width));
    const long long j0 = static_cast<long long>(std::floor((extent.y_min - c.origin_y) / c.tile_height));
    const long long j1 = static_cast<long long>(std::ceil((extent.y_max - c.origin_y) / c.tile_height));
    std::string upper = model;
    for (char& ch : upper) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    for (long long j = j1 - 1; j >= j0; --j)
    {
        for (long long i = i0; i < i1; ++i)
        {
            const double xc = c.origin_x + (i + 0.5) * c.tile_width;
            const double yc = c.origin_y + (j + 0.5) * c.tile_height;
            std::string key = c.key_template;
            bool ok = true;
            auto fill = [&](const std::string& token, const std::string& value) {
                const std::size_t at = key.find(token);
                if (at == std::string::npos) return;
                if (value.empty()) ok = false;
                key.replace(at, token.size(), value);
            };
            fill("{model}", model);
            fill("{MODEL}", upper);
            fill("{sheet}", nztopo50_sheet(xc, yc));
            fill("{os10}", os_grid_reference(xc, yc, 1));
            fill("{os5}", os_grid_quadrant(xc, yc));
            fill("{os1}", os_grid_reference(xc, yc, 2));
            if (ok) keys.push_back(key);
        }
    }
    return keys;
}


fetch_extent grid_search_extent(const lsdtt_fetch_options&   o,
                                const lsdtt_dem_source_info& info,
                                int                          target_epsg)
{
    fetch_extent area = output_area_in(o, info.grid_epsg, target_epsg);
    const double pad = 2.0 * o.grid_spacing;
    area.x_min -= pad; area.y_min -= pad; area.x_max += pad; area.y_max += pad;
    return area;
}


bool gsi_decode_rgb(int r, int g, int b, double& height_metres)
{
    const long x = 65536L * r + 256L * g + b;
    if (x == (1L << 23)) return false;
    height_metres = 0.01 * static_cast<double>(x < (1L << 23) ? x : x - (1L << 24));
    return true;
}


bool gsi_tile_to_geotiff(const std::string& png_bytes, int z, int x, int y,
                         const std::string& out_path, std::string& error)
{
    const std::string png_path = out_path + ".png";
    {
        VSILFILE* f = VSIFOpenL(png_path.c_str(), "wb");
        if (f == nullptr ||
            VSIFWriteL(png_bytes.data(), 1, png_bytes.size(), f) != png_bytes.size())
        {
            if (f) VSIFCloseL(f);
            error = "could not stage the tile in memory";
            return false;
        }
        VSIFCloseL(f);
    }
    GDALDatasetH png = nullptr;
    {
        quiet_gdal_errors quiet;
        png = GDALOpen(png_path.c_str(), GA_ReadOnly);
    }
    if (png == nullptr || GDALGetRasterCount(png) < 3)
    {
        if (png) GDALClose(png);
        VSIUnlink(png_path.c_str());
        error = "the tile is not an RGB PNG";
        return false;
    }
    const int w = GDALGetRasterXSize(png);
    const int h = GDALGetRasterYSize(png);
    std::vector<unsigned char> rgb(static_cast<std::size_t>(w) * h * 3);
    int bands[3] = {1, 2, 3};
    const CPLErr read = GDALDatasetRasterIO(png, GF_Read, 0, 0, w, h, rgb.data(), w, h,
                                            GDT_Byte, 3, bands, 3, 3 * w, 1);
    GDALClose(png);
    VSIUnlink(png_path.c_str());
    if (read != CE_None)
    {
        error = "could not read the tile's pixels";
        return false;
    }
    std::vector<float> z_values(static_cast<std::size_t>(w) * h);
    for (std::size_t i = 0; i < z_values.size(); ++i)
    {
        double height = 0.0;
        z_values[i] = gsi_decode_rgb(rgb[3 * i], rgb[3 * i + 1], rgb[3 * i + 2], height)
                    ? static_cast<float>(height) : -9999.0f;
    }

    lsdtt_slippy_tile tile;
    tile.z = z; tile.x = x; tile.y = y;
    double x_min = 0, y_min = 0, x_max = 0, y_max = 0;
    slippy_tile_bounds_3857(tile, x_min, y_min, x_max, y_max);
    GDALDriverH gtiff = GDALGetDriverByName("GTiff");
    GDALDatasetH ds = gtiff ? GDALCreate(gtiff, out_path.c_str(), w, h, 1, GDT_Float32, nullptr)
                            : nullptr;
    if (ds == nullptr)
    {
        error = "could not create the decoded tile";
        return false;
    }
    double gt[6] = {x_min, (x_max - x_min) / w, 0.0, y_max, 0.0, -(y_max - y_min) / h};
    GDALSetGeoTransform(ds, gt);
    OGRSpatialReferenceH srs = OSRNewSpatialReference(nullptr);
    OSRImportFromEPSG(srs, 3857);
    GDALSetSpatialRef(ds, srs);
    OSRDestroySpatialReference(srs);
    GDALRasterBandH band = GDALGetRasterBand(ds, 1);
    GDALSetRasterNoDataValue(band, -9999.0);
    const CPLErr written = GDALRasterIO(band, GF_Write, 0, 0, w, h, z_values.data(), w, h,
                                        GDT_Float32, 0, 0);
    GDALClose(ds);
    if (written != CE_None)
    {
        VSIUnlink(out_path.c_str());
        error = "could not write the decoded tile";
        return false;
    }
    return true;
}


fetch_extent output_area_wgs84(const lsdtt_fetch_options& o, int target_epsg)
{
    fetch_extent area = output_area_in(o, 4326, target_epsg);
    const double pad = 2.0 * o.grid_spacing / 111320.0;
    area.x_min -= pad; area.y_min -= pad; area.x_max += pad; area.y_max += pad;
    return area;
}


bool gdal_has_lerc()
{
    GDALDriverH gtiff = GDALGetDriverByName("GTiff");
    const char* options = gtiff ? GDALGetMetadataItem(gtiff, GDAL_DMD_CREATIONOPTIONLIST, nullptr)
                                : nullptr;
    return options != nullptr && std::strstr(options, "LERC") != nullptr;
}


std::vector<fetch_extent> wcs_plan_chunks(const fetch_extent& extent, double pixel,
                                          int max_side)
{
    const long long cols = std::llround((extent.x_max - extent.x_min) / pixel);
    const long long rows = std::llround((extent.y_max - extent.y_min) / pixel);
    std::vector<fetch_extent> chunks;
    if (cols <= 0 || rows <= 0 || max_side <= 0) return chunks;
    for (long long r0 = 0; r0 < rows; r0 += max_side)
    {
        const long long r1 = std::min(rows, r0 + max_side);
        for (long long c0 = 0; c0 < cols; c0 += max_side)
        {
            const long long c1 = std::min(cols, c0 + max_side);
            chunks.push_back({extent.x_min + c0 * pixel, extent.y_max - r1 * pixel,
                              extent.x_min + c1 * pixel, extent.y_max - r0 * pixel});
        }
    }
    return chunks;
}


std::string wcs_getcoverage_url(const lsdtt_dem_source_info& info,
                                const std::string&           coverage_id,
                                const fetch_extent&          c)
{
    std::ostringstream url;
    url << info.wcs_url << "?SERVICE=WCS&VERSION=" << info.wcs_version
        << "&REQUEST=GetCoverage";
    if (info.wcs_version == "1.0.0")
    {
        const long long w = std::llround((c.x_max - c.x_min) / info.wcs_pixel_metres);
        const long long h = std::llround((c.y_max - c.y_min) / info.wcs_pixel_metres);
        url << "&COVERAGE=" << coverage_id << "&CRS=EPSG:" << info.wcs_native_epsg
            << "&BBOX=" << url_number(c.x_min) << ',' << url_number(c.y_min) << ','
            << url_number(c.x_max) << ',' << url_number(c.y_max)
            << "&WIDTH=" << w << "&HEIGHT=" << h << "&FORMAT=" << info.wcs_format;
    }
    else
    {
        url << "&COVERAGEID=" << coverage_id << "&FORMAT=" << info.wcs_format
            << "&SUBSET=" << info.wcs_x_axis << '(' << url_number(c.x_min) << ','
            << url_number(c.x_max) << ')'
            << "&SUBSET=" << info.wcs_y_axis << '(' << url_number(c.y_min) << ','
            << url_number(c.y_max) << ')';
    }
    return url.str();
}


std::string parse_ows_exception(const std::string& xml)
{
    if (xml.find("ExceptionReport") == std::string::npos) return "";
    std::string text = xml_text_of(xml, "ExceptionText");
    if (text.empty()) text = xml_text_of(xml, "ServiceException");
    return text.empty() ? std::string("the service returned an exception report") : text;
}


std::string multipart_first_payload(const std::string& body, const std::string& boundary)
{
    const std::string marker = "--" + boundary;
    std::size_t pos = body.find(marker);
    while (pos != std::string::npos)
    {
        std::size_t headers = pos + marker.size();
        if (body.compare(headers, 2, "--") == 0) break;   // closing boundary
        std::size_t sep = body.find("\r\n\r\n", headers);
        std::size_t sep_len = 4;
        const std::size_t lf = body.find("\n\n", headers);
        if (lf != std::string::npos && (sep == std::string::npos || lf < sep))
        {
            sep = lf;
            sep_len = 2;
        }
        if (sep == std::string::npos) break;
        const std::size_t next = body.find(marker, sep);
        const std::string head = body.substr(headers, sep - headers);
        std::size_t end = next == std::string::npos ? body.size() : next;
        while (end > sep + sep_len && (body[end - 1] == '\n' || body[end - 1] == '\r')) --end;
        std::string lower = head;
        for (char& ch : lower)
            ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        if (lower.find("xml") == std::string::npos)
            return body.substr(sep + sep_len, end - sep - sep_len);
        pos = next;
    }
    return "";
}


std::string catalogue_stac_collection(const lsdtt_dem_source_info& info,
                                      double                       resolution_metres)
{
    for (const auto& entry : info.stac_collection_by_resolution)
    {
        if (std::abs(entry.first - resolution_metres) <= 1e-6 * entry.first)
            return entry.second;
    }
    return info.stac_collection;
}


std::string catalogue_stac_search_body(const lsdtt_fetch_options& o)
{
    const lsdtt_dem_source_info& info = dem_source_info(o.source);
    const double res = o.source_resolution_metres > 0.0 ? o.source_resolution_metres
                                                        : info.default_resolution_metres;
    lsdtt_fetch_options search = o;
    search.pc_collection = catalogue_stac_collection(info, res);
    std::string body = pc_stac_search_body(search);
    if (!info.stac_fields_include.empty())
    {
        std::string fields;
        for (const std::string& f : info.stac_fields_include)
            fields += (fields.empty() ? "\"" : ",\"") + f + "\"";
        body.insert(body.size() - 1, ",\"fields\":{\"include\":[" + fields + "]}");
    }
    return body;
}


std::vector<std::string> fetch_warp_arguments(const lsdtt_fetch_options& options,
                                              int                        epsg,
                                              const fetch_extent*        explicit_extent)
{
    std::vector<std::string> args;
    const std::string spacing = format_double(options.grid_spacing, 6);
    if (options.use_frame)
    {
        // Straight into the frame: its own extent and cell count, so the grid
        // is exactly the frame's (no -tap, which would snap to the frame
        // origin's multiples but not to the frame extent).
        const lsdtt_oblique_frame& f = options.frame;
        args.insert(args.end(), {"-t_srs", f.to_wkt(), "-te",
                                 format_double(f.x_min, 6), format_double(f.y_min, 6),
                                 format_double(f.x_max, 6), format_double(f.y_max, 6),
                                 "-ts", std::to_string(f.n_cols()), std::to_string(f.n_rows()),
                                 "-r", options.resampling});
    }
    else if (explicit_extent != nullptr)
    {
        // The catalogue sources: the extent is already snapped to the grid
        // (target_extent_for_bbox), which is what -tap would do, but -tap
        // cannot be combined with -te.
        args.insert(args.end(), {"-t_srs", "EPSG:" + std::to_string(epsg),
                                 "-te",
                                 format_double(explicit_extent->x_min, 6),
                                 format_double(explicit_extent->y_min, 6),
                                 format_double(explicit_extent->x_max, 6),
                                 format_double(explicit_extent->y_max, 6),
                                 "-tr", spacing, spacing,
                                 "-r", options.resampling});
    }
    else
    {
        args.insert(args.end(), {"-t_srs", "EPSG:" + std::to_string(epsg),
                                 "-tr", spacing, spacing,
                                 "-r", options.resampling,
                                 "-tap"});
    }
    args.insert(args.end(), {"-dstnodata", format_double(options.no_data_value, 6),
                             "-of", "GTiff",
                             "-co", "COMPRESS=DEFLATE",
                             "-co", "TILED=YES"});
    return args;
}


fetch_extent target_extent_for_bbox(double west, double south, double east, double north,
                                    int epsg, double grid_spacing)
{
    OGRSpatialReferenceH wgs84  = OSRNewSpatialReference(nullptr);
    OGRSpatialReferenceH target = OSRNewSpatialReference(nullptr);
    OSRImportFromEPSG(wgs84, 4326);
    OSRSetAxisMappingStrategy(wgs84, OAMS_TRADITIONAL_GIS_ORDER);
    const OGRErr err = OSRImportFromEPSG(target, epsg);
    OSRSetAxisMappingStrategy(target, OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformationH ct = nullptr;
    if (err == OGRERR_NONE)
    {
        quiet_gdal_errors quiet;
        ct = OCTNewCoordinateTransformation(wgs84, target);
    }
    if (ct == nullptr)
    {
        OSRDestroySpatialReference(wgs84);
        OSRDestroySpatialReference(target);
        throw std::runtime_error("cannot transform WGS84 to EPSG:" + std::to_string(epsg));
    }

    // Densify every edge: a projected box is not the image of the corners.
    const int per_side = 64;
    std::vector<double> xs, ys;
    for (int side = 0; side < 4; ++side)
    {
        for (int i = 0; i <= per_side; ++i)
        {
            const double t = static_cast<double>(i) / per_side;
            switch (side)
            {
                case 0: xs.push_back(west + t * (east - west)); ys.push_back(south); break;
                case 1: xs.push_back(east); ys.push_back(south + t * (north - south)); break;
                case 2: xs.push_back(east - t * (east - west)); ys.push_back(north); break;
                default: xs.push_back(west); ys.push_back(north - t * (north - south)); break;
            }
        }
    }
    std::vector<int> ok(xs.size(), 0);
    OCTTransformEx(ct, static_cast<int>(xs.size()), xs.data(), ys.data(), nullptr, ok.data());
    OCTDestroyCoordinateTransformation(ct);
    OSRDestroySpatialReference(wgs84);
    OSRDestroySpatialReference(target);

    fetch_extent e;
    e.x_min = e.y_min = 1e300;
    e.x_max = e.y_max = -1e300;
    bool any = false;
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        if (!ok[i]) continue;
        any = true;
        e.x_min = std::min(e.x_min, xs[i]); e.x_max = std::max(e.x_max, xs[i]);
        e.y_min = std::min(e.y_min, ys[i]); e.y_max = std::max(e.y_max, ys[i]);
    }
    if (!any) throw std::runtime_error("no box edge point transforms to EPSG:" +
                                       std::to_string(epsg));

    // Snap outward to the grid: the -tap lattice.
    e.x_min = std::floor(e.x_min / grid_spacing) * grid_spacing;
    e.y_min = std::floor(e.y_min / grid_spacing) * grid_spacing;
    e.x_max = std::ceil (e.x_max / grid_spacing) * grid_spacing;
    e.y_max = std::ceil (e.y_max / grid_spacing) * grid_spacing;
    return e;
}


bool parse_stac_items(const std::string&        json,
                      const std::string&        asset_key,
                      std::vector<stac_item>&   items,
                      std::vector<std::string>& missing_asset,
                      std::string&              next_href,
                      std::string&              next_body,
                      std::string&              error,
                      double                    asset_gsd)
{
    next_href.clear();
    next_body.clear();
    CPLJSONDocument document;
    if (!document.LoadMemory(json))
    {
        error = "could not parse the STAC search response as JSON";
        return false;
    }
    const CPLJSONObject root = document.GetRoot();
    const CPLJSONArray features = root.GetArray("features");
    for (int i = 0; features.IsValid() && i < features.Size(); ++i)
    {
        const CPLJSONObject feature = features[i];
        stac_item item;
        item.id = feature.GetString("id", "");
        item.gsd = feature.GetDouble("properties/gsd", 0.0);
        item.datetime = feature.GetString("properties/datetime", "");
        const CPLJSONArray bbox = feature.GetArray("bbox");
        if (bbox.IsValid() && bbox.Size() >= 4)
        {
            for (int k = 0; k < 4; ++k) item.bbox[k] = bbox[k].ToDouble();
        }
        CPLJSONObject asset = feature.GetObj("assets/" + asset_key);
        if (asset_key.empty())
        {
            asset = CPLJSONObject();
            asset.Deinit();
            for (const CPLJSONObject& candidate : feature.GetObj("assets").GetChildren())
            {
                double g = candidate.GetDouble("gsd", 0.0);
                if (g <= 0.0) g = candidate.GetDouble("eo:gsd", 0.0);
                const std::string type = candidate.GetString("type", "");
                if (std::abs(g - asset_gsd) <= 1e-6 * asset_gsd &&
                    type.find("geotiff") != std::string::npos)
                {
                    asset = candidate;
                    item.gsd = g;
                    break;
                }
            }
        }
        item.href = asset.IsValid() ? asset.GetString("href", "") : "";
        if (item.href.empty())
        {
            missing_asset.push_back(item.id.empty() ? "(unnamed item)" : item.id);
            continue;
        }
        if (item.gsd <= 0.0 && asset.IsValid())
            item.gsd = asset.GetDouble("gsd", 0.0);
        items.push_back(item);
    }
    const CPLJSONArray links = root.GetArray("links");
    for (int i = 0; links.IsValid() && i < links.Size(); ++i)
    {
        const CPLJSONObject link = links[i];
        if (link.GetString("rel", "") != "next") continue;
        next_href = link.GetString("href", "");
        const CPLJSONObject body = link.GetObj("body");
        if (body.IsValid()) next_body = body.Format(CPLJSONObject::PrettyFormat::Plain);
    }
    return true;
}


std::string merge_stac_next_body(const std::string& original, const std::string& next)
{
    CPLJSONDocument base, over;
    if (next.empty() || !over.LoadMemory(next)) return next;
    if (!base.LoadMemory(original)) return next;
    CPLJSONObject root = base.GetRoot();
    for (const CPLJSONObject& member : over.GetRoot().GetChildren())
    {
        root.Delete(member.GetName());
        root.Add(member.GetName(), member);
    }
    return root.Format(CPLJSONObject::PrettyFormat::Plain);
}


std::vector<stac_item> latest_stac_item_per_tile(const std::vector<stac_item>& items)
{
    auto tile_of = [](const std::string& id) {
        std::string out;
        std::size_t start = 0;
        while (start <= id.size())
        {
            const std::size_t end = std::min(id.find('_', start), id.size());
            const std::string token = id.substr(start, end - start);
            const bool year = token.size() == 4 &&
                              std::all_of(token.begin(), token.end(),
                                          [](char c) { return c >= '0' && c <= '9'; });
            if (!year) out += (out.empty() ? "" : "_") + token;
            start = end + 1;
        }
        return out;
    };
    std::vector<std::pair<std::string, stac_item>> best;
    for (const stac_item& item : items)
    {
        const std::string tile = tile_of(item.id);
        auto it = std::find_if(best.begin(), best.end(),
                               [&](const auto& b) { return b.first == tile; });
        if (it == best.end())
            best.emplace_back(tile, item);
        else if (std::make_pair(item.datetime, item.id) >
                 std::make_pair(it->second.datetime, it->second.id))
            it->second = item;
    }
    std::vector<stac_item> out;
    for (const auto& b : best) out.push_back(b.second);
    sort_stac_items(out);
    return out;
}


std::vector<stac_item> stac_items_in_bbox(const std::vector<stac_item>& items, double west,
                                          double south, double east, double north)
{
    std::vector<stac_item> out;
    for (const stac_item& item : items)
    {
        const bool no_bbox = item.bbox[0] == 0.0 && item.bbox[1] == 0.0 &&
                             item.bbox[2] == 0.0 && item.bbox[3] == 0.0;
        if (no_bbox || (item.bbox[0] < east && item.bbox[2] > west &&
                        item.bbox[1] < north && item.bbox[3] > south))
            out.push_back(item);
    }
    return out;
}


std::vector<stac_item> filter_stac_items_by_gsd(const std::vector<stac_item>& items,
                                                double resolution_metres)
{
    if (!(resolution_metres > 0.0)) return items;
    std::vector<stac_item> kept;
    for (const stac_item& item : items)
    {
        // gsd is a nominal product resolution: compare with a tolerance so
        // 10 and 10.0000001 are the same product.
        if (std::abs(item.gsd - resolution_metres) <= 1e-6 * resolution_metres)
            kept.push_back(item);
    }
    return kept;
}


std::vector<double> distinct_gsds(const std::vector<stac_item>& items)
{
    std::vector<double> out;
    for (const stac_item& item : items)
    {
        if (std::find(out.begin(), out.end(), item.gsd) == out.end())
            out.push_back(item.gsd);
    }
    std::sort(out.begin(), out.end());
    return out;
}


void sort_stac_items(std::vector<stac_item>& items)
{
    std::stable_sort(items.begin(), items.end(),
                     [](const stac_item& a, const stac_item& b) { return a.id < b.id; });
}


std::int64_t estimate_output_pixels(const fetch_extent& extent, double grid_spacing)
{
    const double cols = std::ceil((extent.x_max - extent.x_min) / grid_spacing - 1e-9);
    const double rows = std::ceil((extent.y_max - extent.y_min) / grid_spacing - 1e-9);
    return static_cast<std::int64_t>(std::max(0.0, cols) * std::max(0.0, rows));
}


double estimate_download_megabytes(double n_cells, double bytes_per_cell)
{
    return n_cells * bytes_per_cell / 1.0e6;
}


double geographic_cell_count(double west, double south, double east, double north,
                             double cell_degrees)
{
    return ((east - west) / cell_degrees) * ((north - south) / cell_degrees);
}


bool stamp_provenance(const std::string&                                       path,
                      const std::vector<std::pair<std::string, std::string>>& items,
                      std::string&                                             error)
{
    gdal_dataset_guard ds;
    {
        quiet_gdal_errors quiet;
        ds.reset(GDALOpen(path.c_str(), GA_Update));
    }
    if (!ds)
    {
        error = "could not open " + path + " to record its provenance";
        return false;
    }
    for (const auto& kv : items)
    {
        GDALSetMetadataItem(ds.get(), kv.first.c_str(), kv.second.c_str(), nullptr);
    }
    ds.reset();
    return true;
}

} // namespace fetch_detail


// =============================================================================
// Public helpers
// =============================================================================

bool parse_dem_source(const std::string& name, lsdtt_dem_source& out)
{
    if (name == "cop30_aws")          { out = lsdtt_dem_source::cop30_aws;          return true; }
    if (name == "opentopography")     { out = lsdtt_dem_source::opentopography;     return true; }
    if (name == "planetary_computer") { out = lsdtt_dem_source::planetary_computer; return true; }
    for (const lsdtt_dem_source_info& info : dem_source_catalogue())
    {
        if (!info.legacy_path && info.name == name) { out = info.id; return true; }
    }
    return false;
}


std::string dem_source_name(lsdtt_dem_source source)
{
    switch (source)
    {
        case lsdtt_dem_source::cop30_aws:          return "cop30_aws";
        case lsdtt_dem_source::opentopography:     return "opentopography";
        case lsdtt_dem_source::planetary_computer: return "planetary_computer";
        case lsdtt_dem_source::us_3dep:            return "us_3dep";
        case lsdtt_dem_source::arctic_arcticdem:   return "arctic_arcticdem";
        case lsdtt_dem_source::antarctic_rema:     return "antarctic_rema";
        case lsdtt_dem_source::no_dtm:             return "no_dtm";
        case lsdtt_dem_source::nl_ahn:             return "nl_ahn";
        case lsdtt_dem_source::es_mdt:             return "es_mdt";
        case lsdtt_dem_source::pl_nmt:             return "pl_nmt";
        case lsdtt_dem_source::ca_hrdem:           return "ca_hrdem";
        case lsdtt_dem_source::ch_swissalti3d:     return "ch_swissalti3d";
        case lsdtt_dem_source::nz_lidar:           return "nz_lidar";
        case lsdtt_dem_source::gb_sct_lidar:       return "gb_sct_lidar";
        case lsdtt_dem_source::jp_gsi:             return "jp_gsi";
    }
    return "unknown";
}



// -----------------------------------------------------------------------------
// The DEM source catalogue. Access details checked live on the date in
// `accessed`; national endpoints change, so the network tests re-check them.
// -----------------------------------------------------------------------------
const std::vector<lsdtt_dem_source_info>& dem_source_catalogue()
{
    static const std::vector<lsdtt_dem_source_info> catalogue = [] {
        std::vector<lsdtt_dem_source_info> rows;

        lsdtt_dem_source_info cop;
        cop.id = lsdtt_dem_source::cop30_aws;
        cop.name = "cop30_aws";
        cop.title = "Copernicus GLO-30 (AWS Open Data)";
        cop.provider = "ESA / Copernicus, hosted on the AWS Open Data registry";
        cop.licence_attribution =
            "Copernicus DEM GLO-30: produced using Copernicus WorldDEM-30 "
            "(c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH "
            "2014-2018, provided under COPERNICUS by the European Union and ESA.";
        cop.coverage_text = "global";
        cop.description =
            "Copernicus GLO-30 -- a DIGITAL SURFACE MODEL (canopy and buildings "
            "included) from TanDEM-X radar, nominal 30 m (1 arc-second). Heights "
            "in metres above EGM2008 (geoid). Global and keyless: start here.";
        cop.accessed = "2026-07-20";
        cop.legacy_path = true;
        cop.native_resolutions_metres = {30.0};
        cop.default_resolution_metres = 30.0;
        cop.surface_models = {"dsm"};
        cop.default_surface_model = "dsm";
        cop.vertical_datum = "EGM2008";
        cop.coverage_envelopes = {{-180.0, -90.0, 180.0, 90.0}};
        rows.push_back(cop);

        lsdtt_dem_source_info ot;
        ot.id = lsdtt_dem_source::opentopography;
        ot.name = "opentopography";
        ot.title = "OpenTopography global datasets (needs an API key)";
        ot.provider = "OpenTopography";
        ot.licence_attribution = "See the citation printed for the ot_dem_type used.";
        ot.coverage_text = "depends on ot_dem_type";
        ot.description =
            "OpenTopography global DEMs (ot_dem_type picks the dataset; run -h "
            "ot_dem_type). Surface or terrain model and datum depend on the "
            "dataset. Needs a free API key.";
        ot.accessed = "2026-07-22";
        ot.legacy_path = true;
        ot.coverage_envelopes = {{-180.0, -90.0, 180.0, 90.0}};
        rows.push_back(ot);

        lsdtt_dem_source_info pc;
        pc.id = lsdtt_dem_source::planetary_computer;
        pc.name = "planetary_computer";
        pc.title = "Microsoft Planetary Computer STAC collection";
        pc.provider = "Microsoft Planetary Computer";
        pc.licence_attribution = "See the licence of the pc_collection used.";
        pc.coverage_text = "depends on pc_collection";
        pc.description =
            "Any Planetary Computer raster collection (pc_collection; default "
            "cop-dem-glo-30, the same data as cop30_aws). Keyless: asset URLs "
            "are signed automatically. Surface or terrain model and datum depend "
            "on the collection; source_resolution_metres picks items by gsd.";
        pc.accessed = "2026-07-22";
        pc.legacy_path = true;
        pc.coverage_envelopes = {{-180.0, -90.0, 180.0, 90.0}};
        rows.push_back(pc);

        lsdtt_dem_source_info dep;
        dep.id = lsdtt_dem_source::us_3dep;
        dep.name = "us_3dep";
        dep.title = "USGS 3DEP seamless DEM (via Planetary Computer)";
        dep.provider = "U.S. Geological Survey, 3D Elevation Program; hosted by "
                       "Microsoft Planetary Computer";
        dep.licence_attribution =
            "U.S. Geological Survey, 3D Elevation Program (3DEP) seamless DEM, "
            "public domain; accessed through Microsoft Planetary Computer "
            "(collection 3dep-seamless).";
        dep.coverage_text = "conterminous US, Alaska, Hawaii, Puerto Rico and "
                            "US territories (patchy outside CONUS)";
        dep.description =
            "USGS 3DEP seamless DEM -- a DIGITAL TERRAIN MODEL (bare earth, "
            "mostly from lidar). Nominal 10 m (1/3 arc-second, the default) or "
            "30 m (1 arc-second) on a geographic NAD83 grid, so it is always "
            "reprojected to square cells. Heights in metres above NAVD88 "
            "(geoid).";
        dep.accessed = "2026-11-01";
        dep.native_resolutions_metres = {10.0, 30.0};
        dep.default_resolution_metres = 10.0;
        dep.surface_models = {"dtm"};
        dep.default_surface_model = "dtm";
        dep.vertical_datum = "NAVD88";
        dep.coverage_envelopes = {
            {-125.0,  24.0,  -66.0,  50.0},   // conterminous US
            {-180.0,  51.0, -129.0,  72.0},   // Alaska
            { 172.0,  51.0,  180.0,  54.0},   // western Aleutians
            {-161.0,  18.5, -154.0,  22.5},   // Hawaii
            { -68.0,  17.5,  -64.5,  18.6},   // Puerto Rico, US Virgin Islands
            { 144.5,  13.2,  146.2,  20.6},   // Guam, Northern Marianas
            {-171.0, -14.6, -168.0, -11.0}};  // American Samoa
        dep.stac_search_url = "https://planetarycomputer.microsoft.com/api/stac/v1/search";
        dep.stac_collection = "3dep-seamless";
        dep.stac_asset = "data";
        dep.stac_planetary_computer_signing = true;
        rows.push_back(dep);

        // ArcticDEM and REMA: PGC's mosaics, searched through the PGC STAC API
        // (one collection per resolution), read straight from the public
        // pgc-opendata-dems bucket. No key, no signing.
        const std::vector<std::string> pgc_fields = {
            "id", "bbox", "properties.gsd", "properties.datetime", "assets.dem.href"};

        lsdtt_dem_source_info arc;
        arc.id = lsdtt_dem_source::arctic_arcticdem;
        arc.name = "arctic_arcticdem";
        arc.title = "ArcticDEM mosaic v4.1 (Polar Geospatial Center)";
        arc.provider = "Polar Geospatial Center, University of Minnesota; hosted on "
                       "the AWS Open Data registry";
        arc.licence_attribution =
            "ArcticDEM mosaic v4.1, Polar Geospatial Center, University of Minnesota "
            "(Porter et al., 2023, doi:10.7910/DVN/3VDC4W), CC-BY-4.0. DEMs provided "
            "by the Polar Geospatial Center under NSF-OPP awards 1043681, 1559691, "
            "1542736, 1810976 and 2129685.";
        arc.coverage_text = "land north of 60 N, plus Greenland, Alaska and Kamchatka";
        arc.description =
            "ArcticDEM mosaic v4.1 -- a DIGITAL SURFACE MODEL (canopy and buildings "
            "included) from stereo satellite imagery, composited from scenes of "
            "different dates. Nominal 10 m (the default), 2 m or 32 m on a polar "
            "stereographic grid (EPSG:3413). Heights in metres above the WGS84 "
            "ELLIPSOID, not a geoid: they differ from sea-level heights by tens "
            "of metres.";
        arc.accessed = "2026-11-01";
        arc.native_resolutions_metres = {2.0, 10.0, 32.0};
        arc.default_resolution_metres = 10.0;
        arc.surface_models = {"dsm"};
        arc.default_surface_model = "dsm";
        arc.vertical_datum = "WGS84 ellipsoid";
        arc.coverage_envelopes = {{-180.0, 50.0, 180.0, 90.0}};
        arc.stac_search_url = "https://stac.pgc.umn.edu/api/v1/search";
        arc.stac_asset = "dem";
        arc.stac_collection_by_resolution = {
            {2.0,  "arcticdem-mosaics-v4.1-2m"},
            {10.0, "arcticdem-mosaics-v4.1-10m"},
            {32.0, "arcticdem-mosaics-v4.1-32m"}};
        arc.stac_fields_include = pgc_fields;
        arc.native_grid_projected = true;
        rows.push_back(arc);

        lsdtt_dem_source_info rema;
        rema.id = lsdtt_dem_source::antarctic_rema;
        rema.name = "antarctic_rema";
        rema.title = "REMA mosaic v2.0 (Polar Geospatial Center)";
        rema.provider = "Polar Geospatial Center, University of Minnesota; hosted on "
                        "the AWS Open Data registry";
        rema.licence_attribution =
            "Reference Elevation Model of Antarctica (REMA) mosaic v2.0, Polar "
            "Geospatial Center, University of Minnesota (Howat et al., 2022, "
            "doi:10.7910/DVN/EBW8UC), CC-BY-4.0. DEMs provided by the Polar "
            "Geospatial Center under NSF-OPP awards 1043681, 1542736, 1543501, "
            "1559691, 1810976 and 2129685.";
        rema.coverage_text = "Antarctica and its islands, south of 60 S";
        rema.description =
            "REMA mosaic v2.0 -- a DIGITAL SURFACE MODEL from stereo satellite "
            "imagery, composited from scenes of different dates. Nominal 10 m "
            "(the default), 2 m or 32 m on a polar stereographic grid "
            "(EPSG:3031). Heights in metres above the WGS84 ELLIPSOID, not a "
            "geoid: they differ from sea-level heights by tens of metres.";
        rema.accessed = "2026-11-01";
        rema.native_resolutions_metres = {2.0, 10.0, 32.0};
        rema.default_resolution_metres = 10.0;
        rema.surface_models = {"dsm"};
        rema.default_surface_model = "dsm";
        rema.vertical_datum = "WGS84 ellipsoid";
        rema.coverage_envelopes = {{-180.0, -90.0, 180.0, -60.0}};
        rema.stac_search_url = "https://stac.pgc.umn.edu/api/v1/search";
        rema.stac_asset = "dem";
        rema.stac_collection_by_resolution = {
            {2.0,  "rema-mosaics-v2.0-2m"},
            {10.0, "rema-mosaics-v2.0-10m"},
            {32.0, "rema-mosaics-v2.0-32m"}};
        rema.stac_fields_include = pgc_fields;
        rema.native_grid_projected = true;
        rows.push_back(rema);

        // National lidar DEMs served by WCS. Each row pins the native grid
        // from the service's DescribeCoverage (origin = top-left pixel edge),
        // so every request is an exact window of stored pixels. Where the
        // service has no nodata value, sea and areas outside the country come
        // back as 0; kept as data, a border becomes a cliff to 0 m that cubic
        // resampling overshoots (-229 m measured in the Tatra), so 0 is nodata.
        lsdtt_dem_source_info no;
        no.id = lsdtt_dem_source::no_dtm;
        no.name = "no_dtm";
        no.title = "Norway national height model DTM 1 m (Kartverket)";
        no.provider = "Kartverket (Norwegian Mapping Authority), via the Geonorge WCS";
        no.licence_attribution =
            "Nasjonal hoydemodell (NHM) DTM 1 m, (c) Kartverket, from the Geonorge WCS "
            "(free, no access constraints declared by the service). Credit Kartverket.";
        no.coverage_text = "mainland Norway (not Svalbard)";
        no.description =
            "Kartverket national height model -- a DIGITAL TERRAIN MODEL (bare earth, "
            "from airborne lidar) at 1 m on UTM 33N (EPSG:25833), the whole country "
            "in one grid. Heights in metres above NN2000. Sea and areas outside "
            "Norway are NoData (the service returns 0 there, so cells of exactly "
            "0 m are NoData).";
        no.accessed = "2026-11-01";
        no.native_resolutions_metres = {1.0};
        no.default_resolution_metres = 1.0;
        no.surface_models = {"dtm"};
        no.default_surface_model = "dtm";
        no.vertical_datum = "NN2000";
        no.coverage_envelopes = {{4.0, 57.9, 31.2, 71.25}};
        no.native_grid_projected = true;
        no.wcs_url = "https://wcs.geonorge.no/skwms1/wcs.hoyde-dtm-nhm-25833";
        no.wcs_version = "1.0.0";
        no.wcs_coverage_by_model = {{"dtm", "nhm_dtm_topo_25833"}};
        no.wcs_format = "GeoTIFF";
        no.wcs_native_epsg = 25833;
        no.wcs_origin_x = -100275.0;
        no.wcs_origin_y = 8000275.0;
        no.wcs_pixel_metres = 1.0;
        no.wcs_max_chunk_side = 2000;   // 4000 x 4000 timed out (HTTP 504)
        no.wcs_zero_is_nodata = true;
        rows.push_back(no);

        lsdtt_dem_source_info nl;
        nl.id = lsdtt_dem_source::nl_ahn;
        nl.name = "nl_ahn";
        nl.title = "AHN, Actueel Hoogtebestand Nederland 0.5 m (PDOK)";
        nl.provider = "Rijkswaterstaat and the water boards (AHN), via the PDOK WCS";
        nl.licence_attribution =
            "Actueel Hoogtebestand Nederland (AHN), from the PDOK WCS, CC0 1.0 "
            "(public domain).";
        nl.coverage_text = "the Netherlands";
        nl.description =
            "AHN at 0.5 m on the Dutch RD grid (EPSG:28992): a DIGITAL TERRAIN MODEL "
            "(dtm, the default: bare earth from airborne lidar, with gaps where "
            "buildings and water were removed) or a DIGITAL SURFACE MODEL (dsm). "
            "Heights in metres above NAP. Only the 0.5 m product is served, so "
            "large areas download a lot of data.";
        nl.accessed = "2026-11-01";
        nl.native_resolutions_metres = {0.5};
        nl.default_resolution_metres = 0.5;
        nl.surface_models = {"dtm", "dsm"};
        nl.default_surface_model = "dtm";
        nl.vertical_datum = "NAP";
        nl.coverage_envelopes = {{3.2, 50.7, 7.3, 53.6}};
        nl.native_grid_projected = true;
        nl.wcs_url = "https://service.pdok.nl/rws/ahn/wcs/v1_0";
        nl.wcs_version = "2.0.1";
        nl.wcs_coverage_by_model = {{"dtm", "dtm_05m"}, {"dsm", "dsm_05m"}};
        nl.wcs_format = "image/tiff";
        nl.wcs_native_epsg = 28992;
        nl.wcs_origin_x = 10000.0;
        nl.wcs_origin_y = 618750.0;
        nl.wcs_pixel_metres = 0.5;
        nl.wcs_max_chunk_side = 4000;   // the service's MAXSIZE
        rows.push_back(nl);

        lsdtt_dem_source_info es;
        es.id = lsdtt_dem_source::es_mdt;
        es.name = "es_mdt";
        es.title = "Spain MDT05 5 m (Instituto Geografico Nacional)";
        es.provider = "Instituto Geografico Nacional (IGN), via the IDEE WCS";
        es.licence_attribution =
            "MDT05, (c) Instituto Geografico Nacional de Espana, from the IDEE WCS, "
            "CC BY 4.0 scne.es.";
        es.coverage_text = "mainland Spain and the Balearic Islands (not the Canaries)";
        es.description =
            "IGN MDT05 -- a DIGITAL TERRAIN MODEL (bare earth, from airborne lidar) "
            "at 5 m on UTM 30N (EPSG:25830) for the whole peninsula. Heights in "
            "metres above mean sea level at Alicante (REDNAP). Fetched as text to "
            "keep millimetre heights (the service's GeoTIFF is whole metres). Sea "
            "and areas outside Spain are NoData (the service returns 0 there, so "
            "cells of exactly 0 m are NoData).";
        es.accessed = "2026-11-01";
        es.native_resolutions_metres = {5.0};
        es.default_resolution_metres = 5.0;
        es.surface_models = {"dtm"};
        es.default_surface_model = "dtm";
        es.vertical_datum = "REDNAP (Alicante)";
        es.coverage_envelopes = {{-9.4, 35.9, 4.4, 43.9}};
        es.native_grid_projected = true;
        es.wcs_url = "https://servicios.idee.es/wcs-inspire/mdt";
        es.wcs_version = "2.0.1";
        es.wcs_coverage_by_model = {{"dtm", "Elevacion25830_5"}};
        es.wcs_format = "application/asc";
        es.wcs_multipart_ascii = true;
        es.wcs_native_epsg = 25830;
        es.wcs_origin_x = -19452.5;
        es.wcs_origin_y = 4865682.5;
        es.wcs_pixel_metres = 5.0;
        es.wcs_max_chunk_side = 2000;   // MAXSIZE is 4096; text is ~9 bytes a cell
        es.wcs_bytes_per_cell = 9.0;
        es.wcs_zero_is_nodata = true;
        rows.push_back(es);

        lsdtt_dem_source_info pl;
        pl.id = lsdtt_dem_source::pl_nmt;
        pl.name = "pl_nmt";
        pl.title = "Poland NMT 1 m grid (GUGiK)";
        pl.provider = "Glowny Urzad Geodezji i Kartografii (GUGiK), via the "
                      "geoportal.gov.pl WCS";
        pl.licence_attribution =
            "Numeryczny Model Terenu (NMT) 1 m, Glowny Urzad Geodezji i Kartografii "
            "(GUGiK), from the geoportal.gov.pl WCS (no fees or access constraints "
            "declared by the service). Credit GUGiK.";
        pl.coverage_text = "Poland";
        pl.description =
            "GUGiK NMT -- a DIGITAL TERRAIN MODEL (bare earth, from airborne lidar) "
            "at 1 m on the PL-1992 grid (EPSG:2180). Heights in metres in the "
            "Kronstadt 1986 system (PL-KRON86-NH). Areas outside Poland are "
            "NoData (the service returns 0 there, so cells of exactly 0 m are "
            "NoData). The service is slow (about 1 MB/s).";
        pl.accessed = "2026-11-01";
        pl.native_resolutions_metres = {1.0};
        pl.default_resolution_metres = 1.0;
        pl.surface_models = {"dtm"};
        pl.default_surface_model = "dtm";
        pl.vertical_datum = "PL-KRON86-NH";
        pl.coverage_envelopes = {{14.1, 49.0, 24.2, 54.9}};
        pl.native_grid_projected = true;
        pl.wcs_url = "https://mapy.geoportal.gov.pl/wss/service/PZGIK/NMT/GRID1/WCS/"
                     "DigitalTerrainModelFormatTIFF";
        pl.wcs_version = "2.0.1";
        pl.wcs_coverage_by_model = {{"dtm", "DTM_PL-KRON86-NH_TIFF"}};
        pl.wcs_format = "image/tiff";
        pl.wcs_native_epsg = 2180;
        pl.wcs_origin_x = 160828.343266;
        pl.wcs_origin_y = 796521.669410;
        pl.wcs_pixel_metres = 1.0;
        pl.wcs_max_chunk_side = 2000;
        pl.wcs_zero_is_nodata = true;
        rows.push_back(pl);

        lsdtt_dem_source_info ca;
        ca.id = lsdtt_dem_source::ca_hrdem;
        ca.name = "ca_hrdem";
        ca.title = "Canada HRDEM mosaic 1 m / 2 m (NRCan)";
        ca.provider = "Natural Resources Canada, via the CCMEO datacube STAC";
        ca.licence_attribution =
            "High Resolution Digital Elevation Model (HRDEM) mosaic, Natural Resources "
            "Canada. Contains information licensed under the Open Government Licence "
            "- Canada.";
        ca.coverage_text = "Canada where airborne lidar exists (patchy; gaps are NoData)";
        ca.description =
            "NRCan HRDEM mosaic: a DIGITAL TERRAIN MODEL (dtm, the default: bare earth "
            "from airborne lidar) or a DIGITAL SURFACE MODEL (dsm), at 1 m (the "
            "default) or 2 m on Canada Atlas Lambert (EPSG:3979). Heights in metres "
            "above CGVD2013. Coverage follows lidar campaigns, so expect NoData gaps.";
        ca.accessed = "2026-11-01";
        ca.native_resolutions_metres = {1.0, 2.0};
        ca.default_resolution_metres = 1.0;
        ca.surface_models = {"dtm", "dsm"};
        ca.default_surface_model = "dtm";
        ca.vertical_datum = "CGVD2013";
        ca.coverage_envelopes = {{-141.1, 41.6, -52.5, 83.2}};
        ca.native_grid_projected = true;
        ca.stac_search_url = "https://datacube.services.geo.ca/stac/api/search";
        ca.stac_collection_by_resolution = {{1.0, "hrdem-mosaic-1m"},
                                            {2.0, "hrdem-mosaic-2m"}};
        ca.stac_asset = "dtm";
        ca.stac_asset_by_model = {{"dtm", "dtm"}, {"dsm", "dsm"}};
        ca.stac_items_carry_gsd = false;
        rows.push_back(ca);

        lsdtt_dem_source_info ch;
        ch.id = lsdtt_dem_source::ch_swissalti3d;
        ch.name = "ch_swissalti3d";
        ch.title = "swissALTI3D 0.5 m / 2 m (swisstopo)";
        ch.provider = "Federal Office of Topography swisstopo, via the data.geo.admin.ch STAC";
        ch.licence_attribution =
            "swissALTI3D, Federal Office of Topography swisstopo (free geodata; the "
            "source reference is required).";
        ch.coverage_text = "Switzerland and Liechtenstein";
        ch.description =
            "swisstopo swissALTI3D -- a DIGITAL TERRAIN MODEL (bare earth, from lidar "
            "and stereo imagery) at 2 m (the default) or 0.5 m on the Swiss LV95 grid "
            "(EPSG:2056), in 1 km tiles; the latest acquisition of each tile is used. "
            "Heights in metres above LN02.";
        ch.accessed = "2026-11-01";
        ch.native_resolutions_metres = {0.5, 2.0};
        ch.default_resolution_metres = 2.0;
        ch.surface_models = {"dtm"};
        ch.default_surface_model = "dtm";
        ch.vertical_datum = "LN02";
        ch.coverage_envelopes = {{5.9, 45.8, 10.5, 47.85}};
        ch.native_grid_projected = true;
        ch.stac_search_url = "https://data.geo.admin.ch/api/stac/v1/search";
        ch.stac_collection = "ch.swisstopo.swissalti3d";
        ch.stac_asset_by_gsd = true;
        ch.stac_latest_per_tile = true;
        rows.push_back(ch);

        // Grid-keyed public buckets: file names computed from coordinates.
        lsdtt_dem_source_info nz;
        nz.id = lsdtt_dem_source::nz_lidar;
        nz.name = "nz_lidar";
        nz.title = "New Zealand national lidar 1 m DEM / DSM (LINZ)";
        nz.provider = "Toitu Te Whenua Land Information New Zealand, on the AWS Open "
                      "Data registry (nz-elevation)";
        nz.licence_attribution =
            "New Zealand LiDAR 1m DEM / DSM, Toitu Te Whenua Land Information New "
            "Zealand, CC BY 4.0.";
        nz.coverage_text = "New Zealand where lidar has been flown (nearly all of it)";
        nz.description =
            "LINZ national elevation: a DIGITAL TERRAIN MODEL (dtm, the default: bare "
            "earth from airborne lidar) or a DIGITAL SURFACE MODEL (dsm) at 1 m on "
            "NZTM2000 (EPSG:2193), in NZTopo50 sheets. Heights in metres above "
            "NZVD2016. The files are LERC-compressed: needs a GDAL built with LERC "
            "(otherwise the LINZ lidar DEMs are available from OpenTopography with "
            "an API key).";
        nz.accessed = "2026-11-01";
        nz.native_resolutions_metres = {1.0};
        nz.default_resolution_metres = 1.0;
        nz.surface_models = {"dtm", "dsm"};
        nz.default_surface_model = "dtm";
        nz.vertical_datum = "NZVD2016";
        nz.coverage_envelopes = {{166.0, -47.4, 178.7, -34.3}};
        nz.native_grid_projected = true;
        nz.grid_base_url = "https://nz-elevation.s3.ap-southeast-2.amazonaws.com/";
        nz.grid_epsg = 2193;
        nz.grid_model_names = {{"dtm", "dem"}};   // the bucket's folders say dem/dsm
        nz.grid_collections = {{"new-zealand", "new-zealand/new-zealand/{model}_1m/2193/{sheet}.tiff",
                                24000.0, 36000.0, 1012000.0, 6234000.0}};
        nz.grid_sentinel_key = "new-zealand/new-zealand/dem_1m/2193/collection.json";
        nz.requires_lerc = true;
        rows.push_back(nz);

        lsdtt_dem_source_info sct;
        sct.id = lsdtt_dem_source::gb_sct_lidar;
        sct.name = "gb_sct_lidar";
        sct.title = "Lidar for Scotland (Scottish Remote Sensing Portal)";
        sct.provider = "Scottish Government and partners, Scottish Remote Sensing Portal "
                       "(srsp-open-data on AWS)";
        sct.licence_attribution =
            "Lidar for Scotland, Scottish Remote Sensing Portal, Open Government "
            "Licence. Attribution differs by phase, e.g. Phase 1: Crown copyright "
            "Scottish Government, SEPA and Scottish Water (2012); Phase 3: Crown "
            "copyright Scottish Government, SEPA and Fugro (2020). See each phase's "
            "record on spatialdata.gov.scot.";
        sct.coverage_text = "parts of Scotland (survey phases 1-6 and the national "
                            "lidar programme; patchy)";
        sct.description =
            "Lidar for Scotland: a DIGITAL TERRAIN MODEL (dtm, the default) or a "
            "DIGITAL SURFACE MODEL (dsm) on the British National Grid (EPSG:27700), "
            "from one survey collection at a time: phases 1-2 (1 m), phases 3-6 and "
            "the national lidar programme (0.5 m). By default the newest collection "
            "with tiles in the area; source_collection_choice names another. "
            "Heights in metres above ODN (Newlyn).";
        sct.accessed = "2026-11-01";
        sct.native_resolutions_metres = {0.5};
        sct.default_resolution_metres = 0.5;
        sct.surface_models = {"dtm", "dsm"};
        sct.default_surface_model = "dtm";
        sct.vertical_datum = "ODN";
        sct.coverage_envelopes = {{-8.7, 54.6, -0.7, 60.9}};
        sct.native_grid_projected = true;
        sct.grid_base_url = "https://srsp-open-data.s3.eu-west-2.amazonaws.com/";
        sct.grid_epsg = 27700;
        // Oldest first (phase numbers follow release order; the national
        // programme is the newest): "latest" takes the last one with tiles.
        sct.grid_collections = {
            {"phase-1", "lidar/phase-1/{model}/27700/gridded/{os10}_1M_{MODEL}_PHASE1.tif",
             10000.0, 10000.0, 0.0, 0.0},
            {"phase-2", "lidar/phase-2/{model}/27700/gridded/{os10}_1M_{MODEL}_PHASE2.tif",
             10000.0, 10000.0, 0.0, 0.0},
            {"phase-3", "lidar/phase-3/{model}/27700/gridded/{os5}_50CM_{MODEL}_PHASE3.tif",
             5000.0, 5000.0, 0.0, 0.0},
            {"phase-4", "lidar/phase-4/{model}/27700/gridded/{os5}_50CM_{MODEL}_PHASE4.tif",
             5000.0, 5000.0, 0.0, 0.0},
            {"phase-5", "lidar/phase-5/{model}/27700/gridded/{os5}_50CM_{MODEL}_PHASE5.tif",
             5000.0, 5000.0, 0.0, 0.0},
            {"phase-6", "lidar/phase-6/{model}/27700/gridded/{os5}_50CM_{MODEL}_PHASE6.tif",
             5000.0, 5000.0, 0.0, 0.0},
            {"national-lidar-programme",
             "lidar/national-lidar-programme/{model}/27700/gridded/{os1}_50cm_{MODEL}_ScotlandNationalLiDAR.tif",
             1000.0, 1000.0, 0.0, 0.0}};
        sct.grid_sentinel_key = "lidar/phase-1/dtm/27700/gridded/NH24_1M_DTM_PHASE1.tif";
        rows.push_back(sct);

        lsdtt_dem_source_info jp;
        jp.id = lsdtt_dem_source::jp_gsi;
        jp.name = "jp_gsi";
        jp.title = "Japan GSI elevation tiles (DEM10B / DEM5A / DEM1A)";
        jp.provider = "Geospatial Information Authority of Japan (GSI), elevation tiles";
        jp.licence_attribution =
            "Source: Geospatial Information Authority of Japan, GSI tiles (elevation "
            "tiles), https://maps.gsi.go.jp/development/ichiran.html";
        jp.coverage_text = "Japan (10 m everywhere; 5 m and 1 m where lidar exists)";
        jp.description =
            "GSI elevation tiles -- a DIGITAL TERRAIN MODEL from the national base "
            "DEMs: 10 m (DEM10B, the default, national), 5 m (DEM5A, airborne lidar, "
            "patchy) or 1 m (DEM1A, patchy). Served as Web Mercator PNG tiles that GSI "
            "has already interpolated from the base DEM, so lsdtt3's warp is a second "
            "resampling. Heights in metres above Tokyo Bay mean sea level (T.P.).";
        jp.accessed = "2026-11-01";
        jp.native_resolutions_metres = {1.0, 5.0, 10.0};
        jp.default_resolution_metres = 10.0;
        jp.surface_models = {"dtm"};
        jp.default_surface_model = "dtm";
        jp.vertical_datum = "TP (Tokyo Bay mean sea level)";
        jp.coverage_envelopes = {{122.9, 24.0, 154.0, 45.6}};
        jp.native_grid_projected = true;
        jp.xyz_layers = {
            {1.0, 17, "dem1a_png", "https://cyberjapandata.gsi.go.jp/xyz/dem1a_png/{z}/{x}/{y}.png"},
            {5.0, 15, "dem5a_png", "https://cyberjapandata.gsi.go.jp/xyz/dem5a_png/{z}/{x}/{y}.png"},
            {10.0, 14, "dem_png", "https://cyberjapandata.gsi.go.jp/xyz/dem_png/{z}/{x}/{y}.png"}};
        rows.push_back(jp);

        return rows;
    }();
    return catalogue;
}


const lsdtt_dem_source_info& dem_source_info(lsdtt_dem_source source)
{
    for (const lsdtt_dem_source_info& info : dem_source_catalogue())
    {
        if (info.id == source) return info;
    }
    return dem_source_catalogue().front();   // unreachable: every id has a row
}


std::string dem_source_name_list()
{
    std::string out;
    for (const lsdtt_dem_source_info& info : dem_source_catalogue())
    {
        out += (out.empty() ? "" : ", ") + info.name;
    }
    return out;
}

// -----------------------------------------------------------------------------
// The OpenTopography global raster catalogue.
//
// Per-request area limits are those stated by the globaldem endpoint: 500e6 km2
// for GEDI L3; 125e6 km2 for the 500 m bathymetry products; 4.05e6 km2 for the
// 90 m products; 450e3 km2 for everything else.
//
// Ordered roughly by how often a geomorphologist reaches for them, not
// alphabetically: the Copernicus and SRTM families first, regional and
// specialist products last.
// -----------------------------------------------------------------------------
const std::vector<lsdtt_ot_dataset>& ot_dataset_catalogue()
{
    static const std::vector<lsdtt_ot_dataset> catalogue = {
        {"COP30", "Copernicus GLO-30", "1 arc-second (~30 m)", "Float32", "global",
         "TanDEM-X derived DSM with edited water bodies and coastlines. The "
         "best general-purpose global 30 m surface; the same data the keyless "
         "cop30_aws source serves.",
         "European Space Agency (2024). Copernicus Global Digital Elevation "
         "Model. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9028PQB",
         450000.0},

        {"COP90", "Copernicus GLO-90", "3 arc-second (~90 m)", "Float32", "global",
         "The 90 m Copernicus DSM. Use it when the area is too large for the "
         "30 m request limit.",
         "European Space Agency (2024). Copernicus Global Digital Elevation "
         "Model. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9028PQB",
         4050000.0},

        {"SRTMGL1", "SRTM GL1", "1 arc-second (~30 m)", "Int16", "60 N to 56 S",
         "NASA Shuttle Radar Topography Mission, 1 arc-second, EGM96 heights. "
         "The long-standing comparison baseline in the literature.",
         "NASA Shuttle Radar Topography Mission (SRTM) (2013). Shuttle Radar "
         "Topography Mission (SRTM) Global. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9445JDF",
         450000.0},

        {"SRTMGL1_E", "SRTM GL1 Ellipsoidal", "1 arc-second (~30 m)", "Float32", "60 N to 56 S",
         "SRTM GL1 with heights referenced to the WGS84 ellipsoid rather than "
         "EGM96. Use it when co-registering with GNSS or lidar point clouds.",
         "NASA Shuttle Radar Topography Mission (SRTM) (2013). Shuttle Radar "
         "Topography Mission (SRTM) Global. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9445JDF",
         450000.0},

        {"SRTMGL3", "SRTM GL3", "3 arc-second (~90 m)", "Int16", "60 N to 56 S",
         "SRTM at 3 arc-seconds. Larger request limit than GL1.",
         "NASA Shuttle Radar Topography Mission (SRTM) (2013). Shuttle Radar "
         "Topography Mission (SRTM) Global. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9445JDF",
         4050000.0},

        {"NASADEM", "NASADEM", "1 arc-second (~30 m)", "Int16", "60 N to 56 S",
         "SRTM reprocessed with improved geolocation and void filling from ASTER "
         "GDEM and ICESat. Generally preferable to raw SRTMGL1.",
         "NASA JPL (2021). NASADEM Merged DEM Global 1 arc second V001. "
         "Distributed by OpenTopography. https://doi.org/10.5069/G93T9FD9",
         450000.0},

        {"AW3D30", "ALOS World 3D", "1 arc-second (~30 m)", "Int16", "global",
         "JAXA optical-stereo DSM (PRISM). An independent global 30 m surface, "
         "useful as a cross-check on the radar products.",
         "Japan Aerospace Exploration Agency (2021). ALOS World 3D 30 meter "
         "DEM. V3.2, Jan 2021. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G94M92HB",
         450000.0},

        {"AW3D30_E", "ALOS World 3D Ellipsoidal", "1 arc-second (~30 m)", "Float32", "global",
         "AW3D30 with ellipsoidal rather than geoid heights.",
         "Japan Aerospace Exploration Agency (2021). ALOS World 3D 30 meter "
         "DEM. V3.2, Jan 2021. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G94M92HB",
         450000.0},

        {"GEDTM30", "Global Ensemble DTM", "1 arc-second (~30 m)", "Float32", "global",
         "Machine-learning ensemble bare-earth DTM. A terrain model rather than "
         "a surface model, so vegetation and buildings are removed.",
         "Ho, Y.-F., & Hengl, T. (2025). Global ensemble digital terrain model "
         "30m (GEDTM30) (Version v20250619). Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9BV7DT1",
         450000.0},

        {"EU_DTM", "Continental Europe DTM", "30 m", "Float32", "Europe",
         "Bare-earth DTM for continental Europe, built from national lidar and "
         "photogrammetric surveys.",
         "Hengl, T., Leal Parente, L., Krizan, J., and Bonannella, C. (2022). "
         "Continental Europe Digital Terrain Model. Distributed by "
         "OpenTopography. https://doi.org/10.5069/G99021ZF",
         450000.0},

        {"ANADEM", "ANADEM", "~30 m", "Float32", "South America",
         "Bare-earth DTM for South America, with vegetation bias removed from "
         "Copernicus GLO-30.",
         "Agencia Nacional de Aguas e Saneamento Basico (2025). ANADEM: A "
         "Digital Terrain Model for South America. Distributed by "
         "OpenTopography. https://doi.org/10.5069/G9736P4G",
         450000.0},

        {"CA_MRDEM_DTM", "Canada MRDEM (terrain)", "30 m", "Float32", "Canada",
         "Natural Resources Canada medium-resolution bare-earth model.",
         "Natural Resources Canada (2025). Medium Resolution Digital Elevation "
         "Model of Canada. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9N8780J",
         450000.0},

        {"CA_MRDEM_DSM", "Canada MRDEM (surface)", "30 m", "Float32", "Canada",
         "The surface-model companion to CA_MRDEM_DTM.",
         "Natural Resources Canada (2025). Medium Resolution Digital Elevation "
         "Model of Canada. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9N8780J",
         450000.0},

        {"SRTM15Plus", "SRTM15+ V2.1", "15 arc-second (~500 m)", "Float32", "global, land and sea",
         "Merged topography and bathymetry at 15 arc-seconds. Use it when the "
         "domain has to continue offshore.",
         "Tozer, B., Sandwell, D. T., Smith, W. H. F., Olson, C., Beale, J. R., "
         "& Wessel, P. (2019). Global bathymetry and topography at 15 arc sec: "
         "SRTM15+. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G92R3PT9",
         125000000.0},

        {"GEBCOIceTopo", "GEBCO 2023 (ice surface)", "15 arc-second (~500 m)",
         "Int16", "global, land and sea",
         "GEBCO global bathymetry and topography, with the ice-sheet surface as "
         "the elevation over Greenland and Antarctica.",
         "GEBCO Compilation Group (2023). GEBCO 2023 Grid. Distributed by "
         "OpenTopography. https://doi.org/10.5069/G9D21VTT",
         125000000.0},

        {"GEBCOSubIceTopo", "GEBCO 2023 (bed under ice)", "15 arc-second (~500 m)",
         "Int16", "global, land and sea",
         "As GEBCOIceTopo, but the ice sheets are replaced by the underlying "
         "bed topography.",
         "GEBCO Compilation Group (2023). GEBCO 2023 Grid. Distributed by "
         "OpenTopography. https://doi.org/10.5069/G9D21VTT",
         125000000.0},

        {"GEDI_L3", "GEDI L3 gridded metrics", "1000 m nominal", "Float32", "52 N to 52 S",
         "Gridded GEDI lidar land-surface metrics. A canopy and elevation "
         "metrics product, not a DEM for flow routing.",
         "Dubayah, R. O., Luthcke, S. B., Sabaka, T. J., Nicholas, J. B., "
         "Preaux, S., and Hofton, M. A. (2021). GEDI L3 Gridded Land Surface "
         "Metrics, Version 2. ORNL DAAC, Oak Ridge, Tennessee, USA. Distributed "
         "by OpenTopography. https://doi.org/10.3334/ORNLDAAC/1952",
         500000000.0}
    };
    return catalogue;
}


const lsdtt_ot_dataset* ot_dataset_info(const std::string& demtype)
{
    for (const lsdtt_ot_dataset& entry : ot_dataset_catalogue())
    {
        if (entry.demtype == demtype) return &entry;
    }
    return nullptr;
}


std::string ot_dataset_demtype_list()
{
    std::string out;
    for (const lsdtt_ot_dataset& entry : ot_dataset_catalogue())
    {
        if (!out.empty()) out += ", ";
        out += entry.demtype;
    }
    return out;
}


int utm_epsg_for_point(double lon, double lat)
{
    // Normalise so that lon = 180 lands in zone 60 rather than a non-existent 61.
    double normalised = std::fmod(lon + 180.0, 360.0);
    if (normalised < 0.0) normalised += 360.0;

    int zone = static_cast<int>(normalised / 6.0) + 1;
    if (zone > 60) zone = 60;
    if (zone < 1)  zone = 1;

    return (lat >= 0.0 ? 32600 : 32700) + zone;
}


std::vector<std::string> cop30_tiles_for_bbox(double west, double south,
                                              double east, double north)
{
    std::vector<std::string> tiles;

    const int lat_begin = static_cast<int>(std::floor(south));
    const int lat_end   = static_cast<int>(std::floor(north));
    const int lon_begin = static_cast<int>(std::floor(west));
    const int lon_end   = static_cast<int>(std::floor(east));

    for (int lat = lat_begin; lat <= lat_end; ++lat)
    {
        for (int lon = lon_begin; lon <= lon_end; ++lon)
        {
            char buffer[32];
            std::snprintf(buffer, sizeof(buffer), "%c%02d_00_%c%03d_00",
                          lat >= 0 ? 'N' : 'S', std::abs(lat),
                          lon >= 0 ? 'E' : 'W', std::abs(lon));
            tiles.emplace_back(buffer);
        }
    }
    return tiles;
}


std::string cop30_tile_url(const std::string& tile_name)
{
    return std::string(k_cop30_url_prefix) + tile_name + "_DEM/" +
           "Copernicus_DSM_COG_10_" + tile_name + "_DEM.tif";
}


namespace
{

// The WGS84 envelope of a raster's footprint, or false when `path` is not a
// raster with a geotransform. Every edge is densified (64 steps) and
// transformed, so a projected footprint is covered fully.
bool bounds_from_raster_footprint(const std::string& path, lsdtt_bounds_from_vector_result& result)
{
    GDALDatasetH ds = nullptr;
    {
        quiet_gdal_errors quiet;
        ds = GDALOpenEx(path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY, nullptr, nullptr, nullptr);
    }
    if (ds == nullptr) return false;
    double gt[6];
    if (GDALGetGeoTransform(ds, gt) != CE_None)
    {
        GDALClose(ds);
        return false;
    }
    const int nx = GDALGetRasterXSize(ds), ny = GDALGetRasterYSize(ds);

    OGRSpatialReferenceH wgs84 = OSRNewSpatialReference(nullptr);
    OSRImportFromEPSG(wgs84, 4326);
    OSRSetAxisMappingStrategy(wgs84, OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformationH transform = nullptr;
    bool transform_failed = false;
    const char* wkt = GDALGetProjectionRef(ds);
    if (wkt != nullptr && wkt[0] != '\0')
    {
        OGRSpatialReferenceH src = OSRNewSpatialReference(wkt);
        OSRSetAxisMappingStrategy(src, OAMS_TRADITIONAL_GIS_ORDER);
        const char* auth = OSRGetAuthorityCode(src, nullptr);
        result.source_crs = auth ? (std::string("EPSG:") + auth) : "(named CRS)";
        if (OSRIsSame(src, wgs84) == FALSE)
        {
            quiet_gdal_errors quiet;
            transform = OCTNewCoordinateTransformation(src, wgs84);
            transform_failed = transform == nullptr;
        }
        OSRDestroySpatialReference(src);
    }
    else
    {
        result.source_crs = "(none declared)";
        result.warning = "'" + path + "' declares no CRS; assuming the coordinates are WGS84 "
                         "degrees.";
    }

    constexpr int k_steps = 64;
    std::vector<double> xs, ys;
    for (int k = 0; k <= k_steps; ++k)
    {
        const double t = static_cast<double>(k) / k_steps;
        const double px[4] = {t * nx, t * nx, 0.0, static_cast<double>(nx)};
        const double py[4] = {0.0, static_cast<double>(ny), t * ny, t * ny};
        for (int e = 0; e < 4; ++e)
        {
            xs.push_back(gt[0] + px[e] * gt[1] + py[e] * gt[2]);
            ys.push_back(gt[3] + px[e] * gt[4] + py[e] * gt[5]);
        }
    }
    bool ok = !transform_failed;
    if (transform != nullptr)
    {
        std::vector<int> success(xs.size(), 0);
        {
            quiet_gdal_errors quiet;
            OCTTransformEx(transform, static_cast<int>(xs.size()), xs.data(), ys.data(), nullptr,
                           success.data());
        }
        ok = std::all_of(success.begin(), success.end(), [](int v) { return v != 0; });
        OCTDestroyCoordinateTransformation(transform);
    }
    OSRDestroySpatialReference(wgs84);
    GDALClose(ds);

    result.layer_name = "";
    result.geometry_kind = "raster";
    result.n_features = 1;
    if (!ok)
    {
        result.message = "could not transform the footprint of raster '" + path + "' (" +
                         result.source_crs + ") to WGS84";
        return true;
    }
    result.west = *std::min_element(xs.begin(), xs.end());
    result.east = *std::max_element(xs.begin(), xs.end());
    result.south = *std::min_element(ys.begin(), ys.end());
    result.north = *std::max_element(ys.begin(), ys.end());
    if (!(result.west < result.east) || !(result.south < result.north))
    {
        result.message = "raster '" + path + "' has a zero-area footprint";
        return true;
    }
    result.ok = true;
    result.message = "bounds from the footprint of raster '" + path + "'";
    return true;
}

} // namespace

lsdtt_bounds_from_vector_result bounds_from_vector_file(
    const std::string& path,
    const std::string& layer_name)
{
    lsdtt_bounds_from_vector_result result;

    if (path.empty())
    {
        result.message = "no bounds file given";
        return result;
    }

    // A CLI may resolve its bounds before anything else has registered the
    // GDAL drivers (lsdtt-fetch-lithology and lsdtt-fetch-erosion-rates do).
    gdal_init();

    GDALDatasetH dataset = nullptr;
    {
        quiet_gdal_errors quiet;
        dataset = GDALOpenEx(path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                             nullptr, nullptr, nullptr);
    }
    if (dataset == nullptr)
    {
        // Not a vector: a raster's footprint (e.g. the DEM itself) also
        // defines a box.
        if (bounds_from_raster_footprint(path, result)) return result;
        result.message = "could not open '" + path +
                         "' as a vector dataset (expected FlatGeobuf, "
                         "GeoPackage, or another OGR-readable format)";
        return result;
    }

    // ---- select the layer ----
    OGRLayerH layer = nullptr;
    const int n_layers = GDALDatasetGetLayerCount(dataset);

    if (n_layers <= 0)
    {
        GDALClose(dataset);
        if (bounds_from_raster_footprint(path, result)) return result;
        result.message = "'" + path + "' contains no layers";
        return result;
    }

    if (!layer_name.empty())
    {
        layer = GDALDatasetGetLayerByName(dataset, layer_name.c_str());
        if (layer == nullptr)
        {
            GDALClose(dataset);
            result.message = "layer '" + layer_name + "' not found in '" + path + "'";
            return result;
        }
    }
    else
    {
        layer = GDALDatasetGetLayer(dataset, 0);
        if (n_layers > 1)
        {
            result.warning = "'" + path + "' has " + std::to_string(n_layers) +
                             " layers; using the first. Set bounds_layer to "
                             "choose a different one.";
        }
    }
    if (layer == nullptr)
    {
        GDALClose(dataset);
        result.message = "could not access a layer in '" + path + "'";
        return result;
    }

    const char* raw_layer_name = OGR_L_GetName(layer);
    result.layer_name = raw_layer_name ? raw_layer_name : "";

    // ---- build a transform to WGS84 ----
    // Geometry is transformed before the envelope is taken. Transforming the
    // envelope corners instead would under-cover the true extent for any
    // rotated or strongly curved projection.
    OGRSpatialReferenceH source_srs = OGR_L_GetSpatialRef(layer);
    OGRCoordinateTransformationH transform = nullptr;
    OGRSpatialReferenceH wgs84 = OSRNewSpatialReference(nullptr);
    OSRImportFromEPSG(wgs84, 4326);
    // Without this, EPSG:4326 yields lat/lon rather than lon/lat and every
    // bounding box comes out transposed.
    OSRSetAxisMappingStrategy(wgs84, OAMS_TRADITIONAL_GIS_ORDER);

    if (source_srs != nullptr)
    {
        const char* auth = OSRGetAuthorityCode(source_srs, nullptr);
        result.source_crs = auth ? (std::string("EPSG:") + auth) : "(named CRS)";

        if (OSRIsSame(source_srs, wgs84) == FALSE)
        {
            OGRSpatialReferenceH source_copy = OSRClone(source_srs);
            OSRSetAxisMappingStrategy(source_copy, OAMS_TRADITIONAL_GIS_ORDER);
            {
                quiet_gdal_errors quiet;
                transform = OCTNewCoordinateTransformation(source_copy, wgs84);
            }
            OSRDestroySpatialReference(source_copy);

            if (transform == nullptr)
            {
                OSRDestroySpatialReference(wgs84);
                GDALClose(dataset);
                result.message = "could not build a transform from the layer CRS (" +
                                 result.source_crs + ") to WGS84";
                return result;
            }
        }
    }
    else
    {
        // A file with no CRS is ambiguous. Assuming WGS84 is the least
        // surprising choice, but say so loudly: silently misinterpreting UTM
        // metres as degrees would produce a nonsensical bounding box.
        result.source_crs = "(none declared)";
        result.warning += (result.warning.empty() ? "" : " ");
        result.warning += "'" + path + "' declares no CRS; assuming the "
                          "coordinates are WGS84 degrees.";
    }

    // ---- accumulate the envelope ----
    double min_x =  1e300, min_y =  1e300;
    double max_x = -1e300, max_y = -1e300;

    int n_polygons = 0;
    int n_points   = 0;
    int n_other    = 0;

    OGR_L_ResetReading(layer);
    OGRFeatureH feature = nullptr;
    while ((feature = OGR_L_GetNextFeature(layer)) != nullptr)
    {
        OGRGeometryH geometry = OGR_F_GetGeometryRef(feature);
        if (geometry == nullptr)
        {
            OGR_F_Destroy(feature);
            continue;
        }

        const OGRwkbGeometryType type = wkbFlatten(OGR_G_GetGeometryType(geometry));
        if (type == wkbPolygon || type == wkbMultiPolygon)      ++n_polygons;
        else if (type == wkbPoint || type == wkbMultiPoint)     ++n_points;
        else                                                    ++n_other;

        OGRGeometryH working = OGR_G_Clone(geometry);
        if (transform != nullptr)
        {
            quiet_gdal_errors quiet;
            if (OGR_G_Transform(working, transform) != OGRERR_NONE)
            {
                OGR_G_DestroyGeometry(working);
                OGR_F_Destroy(feature);
                continue;
            }
        }

        OGREnvelope envelope;
        OGR_G_GetEnvelope(working, &envelope);
        min_x = std::min(min_x, envelope.MinX);
        min_y = std::min(min_y, envelope.MinY);
        max_x = std::max(max_x, envelope.MaxX);
        max_y = std::max(max_y, envelope.MaxY);

        OGR_G_DestroyGeometry(working);
        OGR_F_Destroy(feature);
        ++result.n_features;
    }

    if (transform != nullptr) OCTDestroyCoordinateTransformation(transform);
    OSRDestroySpatialReference(wgs84);
    GDALClose(dataset);

    // ---- validate what we found ----
    if (result.n_features == 0)
    {
        result.message = "layer '" + result.layer_name + "' in '" + path +
                         "' contains no usable geometry";
        return result;
    }

    if (n_other > 0)
    {
        result.message = "layer '" + result.layer_name + "' contains " +
                         std::to_string(n_other) + " geometries that are "
                         "neither polygons nor points. Provide either polygons, "
                         "or exactly two points marking opposite corners.";
        return result;
    }

    if (n_polygons > 0 && n_points > 0)
    {
        result.message = "layer '" + result.layer_name + "' mixes polygons and "
                         "points. Provide one kind or the other, so the "
                         "intended extent is unambiguous.";
        return result;
    }

    if (n_points > 0)
    {
        // Two points define the box by a diagonal. Which diagonal does not
        // matter: lower-left/upper-right and lower-right/upper-left give the
        // same min/max, so orientation needs no detection. One point cannot
        // define a box at all.
        if (n_points != 2)
        {
            result.message = "layer '" + result.layer_name + "' contains " +
                             std::to_string(n_points) + " points; exactly two "
                             "are needed, marking opposite corners of the area "
                             "(either diagonal works).";
            return result;
        }
        result.geometry_kind = "points";
    }
    else
    {
        result.geometry_kind = "polygon";
    }

    if (!(min_x < max_x) || !(min_y < max_y))
    {
        result.message = "the geometry in '" + path + "' has zero width or "
                         "height after reprojection, so it does not define an "
                         "area to fetch";
        return result;
    }

    result.ok    = true;
    result.west  = min_x;
    result.south = min_y;
    result.east  = max_x;
    result.north = max_y;
    result.message = "bounds from " + std::to_string(result.n_features) + " " +
                     result.geometry_kind + " feature(s) in layer '" +
                     result.layer_name + "'";
    return result;
}


bool gdal_can_read_lerc()
{
    return fetch_detail::gdal_has_lerc();
}


std::string lerc_source_warning(lsdtt_dem_source source)
{
    const lsdtt_dem_source_info& info = dem_source_info(source);
    if (!info.requires_lerc) return "";
    const std::string gdal = GDALVersionInfo("RELEASE_NAME");
    return info.name + " files are LERC-compressed, which needs a GDAL built with LERC. "
           "This GDAL (" + gdal + ") " +
           (gdal_can_read_lerc() ? std::string("has LERC, so the fetch can go ahead.")
                                 : std::string("does NOT have LERC, so the fetch will be "
                                               "refused; use an lsdtt3 built against a "
                                               "GDAL with LERC (e.g. conda-forge).")) +
           " The LINZ lidar DEMs are also available from OpenTopography "
           "(opentopography.org) with an API key.";
}


bool http_available()
{
    // CPLHTTPEnabled() reports whether GDAL was built with curl at all. It says
    // nothing about whether the network is reachable -- that is probed
    // separately by trying to open the source.
    return CPLHTTPEnabled() != 0;
}


// =============================================================================
// fetch_dem_impl: the fetch itself; fetch_dem adds frame mode around it
// =============================================================================

namespace {

// =============================================================================
// Pre-flight checks, provenance and the catalogue path
// =============================================================================

double resolved_source_resolution(const lsdtt_fetch_options& o, const lsdtt_dem_source_info& info)
{
    return o.source_resolution_metres > 0.0 ? o.source_resolution_metres
                                            : info.default_resolution_metres;
}

std::string resolved_surface_model(const lsdtt_fetch_options& o, const lsdtt_dem_source_info& info)
{
    return (o.surface_model_choice.empty() || o.surface_model_choice == "source_default")
         ? info.default_surface_model : o.surface_model_choice;
}

std::string list_of(const std::vector<double>& values)
{
    std::string out;
    for (double v : values)
    {
        std::string text = format_double(v, 3);
        while (!text.empty() && text.back() == '0') text.pop_back();
        if (!text.empty() && text.back() == '.') text.pop_back();
        out += (out.empty() ? "" : ", ") + text;
    }
    return out;
}

// Source, resolution, surface model and coverage: everything that can be
// refused without the network. Errors, never skips: none of them is the
// network's fault.
bool preflight_source(const lsdtt_fetch_options& o, std::string& error)
{
    const lsdtt_dem_source_info& info = dem_source_info(o.source);
    const std::string model_choice = o.surface_model_choice.empty()
                                   ? std::string("source_default") : o.surface_model_choice;
    if (model_choice != "source_default" && model_choice != "dtm" && model_choice != "dsm")
    {
        error = "surface_model_choice must be source_default, dtm or dsm (got '" +
                model_choice + "')";
        return false;
    }
    const std::string collection_choice = o.source_collection_choice.empty()
                                        ? std::string("latest") : o.source_collection_choice;
    if (collection_choice != "latest")
    {
        if (info.grid_collections.size() < 2)
        {
            error = "source_collection_choice applies to sources published as several "
                    "survey collections (gb_sct_lidar); leave it at latest for " + info.name;
            return false;
        }
        std::string labels;
        bool known = false;
        for (const auto& c : info.grid_collections)
        {
            labels += (labels.empty() ? "" : ", ") + c.label;
            if (c.label == collection_choice) known = true;
        }
        if (!known)
        {
            error = info.name + " has no collection '" + collection_choice +
                    "' (it has: latest, " + labels + ")";
            return false;
        }
    }
    if (info.legacy_path)
    {
        if (model_choice != "source_default")
        {
            error = (o.source == lsdtt_dem_source::cop30_aws)
                  ? "Copernicus GLO-30 is a DSM only; surface_model_choice must be "
                    "source_default for cop30_aws"
                  : "surface_model_choice applies to catalogue sources only; leave it "
                    "at source_default for " + info.name;
            return false;
        }
        return true;
    }

    const double res = resolved_source_resolution(o, info);
    bool offered = false;
    for (double r : info.native_resolutions_metres)
        if (std::abs(r - res) <= 1e-6 * r) offered = true;
    if (!offered)
    {
        error = info.name + " offers " + list_of(info.native_resolutions_metres) +
                " m (nominal); source_resolution_metres " + format_double(res, 3) +
                " is not one of them";
        return false;
    }

    const std::string model = resolved_surface_model(o, info);
    if (std::find(info.surface_models.begin(), info.surface_models.end(), model) ==
        info.surface_models.end())
    {
        std::string offered_models;
        for (const std::string& m : info.surface_models)
            offered_models += (offered_models.empty() ? "" : ", ") + m;
        error = info.name + " has no " + model + " product (it offers: " +
                offered_models + ")";
        return false;
    }

    bool inside = false;
    for (const auto& box : info.coverage_envelopes)
    {
        if (o.west < box[2] && o.east > box[0] && o.south < box[3] && o.north > box[1])
            inside = true;
    }
    if (!inside)
    {
        error = "the area (" + format_double(o.west, 4) + ", " + format_double(o.south, 4) +
                ", " + format_double(o.east, 4) + ", " + format_double(o.north, 4) +
                ") lies outside " + info.name + " coverage (" + info.coverage_text +
                "). No data was fetched; use cop30_aws for global 30 m coverage";
        return false;
    }
    if (!info.xyz_layers.empty() && GDALGetDriverByName("PNG") == nullptr)
    {
        error = info.name + " tiles are PNG and this GDAL build has no PNG driver";
        return false;
    }
    if (info.requires_lerc && !fetch_detail::gdal_has_lerc())
    {
        error = info.name + " files are LERC-compressed and this GDAL build (" +
                std::string(GDALVersionInfo("RELEASE_NAME")) + ") cannot decode LERC; "
                "use an lsdtt3 built against a GDAL with LERC support (e.g. conda-forge), "
                "or get the data from OpenTopography (opentopography.org) with an API key";
        return false;
    }
    return true;
}

// The source-cell size GDAL will actually read: the native resolution, or
// the coarsest power-of-two overview not coarser than the output grid.
double effective_read_resolution(double native, double grid_spacing)
{
    double r = native;
    while (r * 2.0 <= grid_spacing + 1e-9) r *= 2.0;
    return r;
}

bool check_download_estimate(const lsdtt_fetch_options& o, double megabytes,
                             std::string& error)
{
    if (o.max_download_megabytes > 0.0 && megabytes > o.max_download_megabytes)
    {
        error = "the download is estimated at " + format_double(megabytes, 1) +
                " MB, more than max_download_megabytes (" +
                format_double(o.max_download_megabytes, 1) + " MB); fetch a smaller "
                "area, use a coarser grid_spacing or source_resolution_metres, or "
                "raise max_download_megabytes";
        return false;
    }
    return true;
}

// Output pixels and download size, before any request. The download estimate
// assumes Float32 cells and ignores compression, so it errs high.
bool preflight_size(const lsdtt_fetch_options& o, int epsg, std::string& error)
{
    std::int64_t out_pixels = 0;
    if (o.use_frame)
    {
        out_pixels = static_cast<std::int64_t>(o.frame.n_cols()) * o.frame.n_rows();
    }
    else
    {
        try
        {
            out_pixels = fetch_detail::estimate_output_pixels(
                fetch_detail::target_extent_for_bbox(o.west, o.south, o.east, o.north,
                                                     epsg, o.grid_spacing),
                o.grid_spacing);
        }
        catch (const std::exception& e)
        {
            error = e.what();
            return false;
        }
    }
    if (o.max_output_pixels > 0 && out_pixels > o.max_output_pixels)
    {
        error = "the output grid would be " + std::to_string(out_pixels) +
                " pixels, more than max_output_pixels (" +
                std::to_string(o.max_output_pixels) + "); use a coarser grid_spacing, "
                "a smaller area, or raise max_output_pixels";
        return false;
    }

    const lsdtt_dem_source_info& info = dem_source_info(o.source);
    double cell_degrees = 1.0 / 3600.0;   // GLO-30 and the 1 arc-second products
    if (!info.wcs_url.empty())
    {
        // The native cells requested, exactly as the chunk planner will ask.
        double cells = 0.0;
        try
        {
            const fetch_detail::fetch_extent native =
                fetch_detail::wcs_native_extent(o, info, o.use_frame ? 0 : epsg);
            cells = ((native.x_max - native.x_min) / info.wcs_pixel_metres) *
                    ((native.y_max - native.y_min) / info.wcs_pixel_metres);
        }
        catch (const std::exception& e)
        {
            error = e.what();
            return false;
        }
        if (o.max_output_pixels > 0 && cells > static_cast<double>(o.max_output_pixels))
        {
            error = "the " + info.name + " download would be " +
                    format_double(cells, 0) + " native " +
                    list_of({info.wcs_pixel_metres}) + " m cells, more than "
                    "max_output_pixels (" + std::to_string(o.max_output_pixels) +
                    "); fetch a smaller area or raise max_output_pixels";
            return false;
        }
        return check_download_estimate(
            o, fetch_detail::estimate_download_megabytes(cells, info.wcs_bytes_per_cell),
            error);
    }
    if (!info.legacy_path && info.native_grid_projected)
    {
        // Metre cells: count them over the output extent's area.
        const double read_res = effective_read_resolution(
            resolved_source_resolution(o, info), o.grid_spacing);
        const double megabytes = fetch_detail::estimate_download_megabytes(
            static_cast<double>(out_pixels) * (o.grid_spacing / read_res) *
                (o.grid_spacing / read_res),
            4.0);
        return check_download_estimate(o, megabytes, error);
    }
    if (!info.legacy_path)
    {
        const double read_res = effective_read_resolution(
            resolved_source_resolution(o, info), o.grid_spacing);
        cell_degrees = read_res / 111320.0;
    }
    else if (o.source == lsdtt_dem_source::planetary_computer &&
             o.source_resolution_metres > 0.0)
    {
        cell_degrees = o.source_resolution_metres / 111320.0;
    }
    const double megabytes = fetch_detail::estimate_download_megabytes(
        fetch_detail::geographic_cell_count(o.west, o.south, o.east, o.north, cell_degrees),
        4.0);
    return check_download_estimate(o, megabytes, error);
}

// The LSDTT_* items for a legacy-path output.
std::vector<std::pair<std::string, std::string>> legacy_provenance(
    const lsdtt_fetch_options& o, const lsdtt_fetch_result& result)
{
    const lsdtt_dem_source_info& info = dem_source_info(o.source);
    std::string surface = "not recorded (depends on the dataset)";
    std::string datum   = surface;
    std::string resolution = "not filtered";
    std::string dataset;
    if (o.source == lsdtt_dem_source::cop30_aws ||
        (o.source == lsdtt_dem_source::planetary_computer &&
         (o.pc_collection.empty() || o.pc_collection == "cop-dem-glo-30")))
    {
        surface = "DSM";
        datum = "EGM2008";
        resolution = "30";
    }
    if (o.source == lsdtt_dem_source::planetary_computer)
    {
        dataset = o.pc_collection.empty() ? "cop-dem-glo-30" : o.pc_collection;
        if (o.source_resolution_metres > 0.0)
            resolution = format_double(o.source_resolution_metres, 3);
    }
    else if (o.source == lsdtt_dem_source::opentopography)
    {
        dataset = o.ot_dem_type.empty() ? "COP30" : o.ot_dem_type;
        const lsdtt_ot_dataset* d = ot_dataset_info(dataset);
        resolution = d ? d->resolution : "not recorded";
    }
    std::vector<std::pair<std::string, std::string>> items = {
        {"LSDTT_DEM_SOURCE", info.name},
        {"LSDTT_SURFACE_MODEL", surface},
        {"LSDTT_VERTICAL_DATUM", datum},
        {"LSDTT_VERTICAL_UNITS", "metre"},
        {"LSDTT_SOURCE_RESOLUTION_METRES", resolution},
        {"LSDTT_SOURCE_ITEM_COUNT", std::to_string(result.tiles.size())},
        {"LSDTT_RESAMPLING", o.resampling},
        {"LSDTT_OVERVIEW_USED", "false"}};
    if (!dataset.empty()) items.push_back({"LSDTT_SOURCE_DATASET", dataset});
    return items;
}

// Stamp a written output and record what was stamped. A failure to stamp
// fails the fetch: the provenance is part of the product.
lsdtt_fetch_result stamp_result(lsdtt_fetch_result result,
                                const std::vector<std::pair<std::string, std::string>>& items)
{
    if (!result.ok) return result;
    std::string error;
    if (!fetch_detail::stamp_provenance(result.path, items, error))
        return make_error("lsdtt::fetch_dem: " + error);
    result.provenance = items;
    return result;
}

// Percentage of output cells that are not nodata, read row by row so a large
// output is never held whole.
double valid_percent(const std::string& path)
{
    gdal_dataset_guard ds(GDALOpen(path.c_str(), GA_ReadOnly));
    if (!ds) return 0.0;
    GDALRasterBandH band = GDALGetRasterBand(ds.get(), 1);
    int has_nodata = 0;
    const double nodata = GDALGetRasterNoDataValue(band, &has_nodata);
    const int cols = GDALGetRasterXSize(ds.get());
    const int rows = GDALGetRasterYSize(ds.get());
    std::vector<float> line(static_cast<std::size_t>(cols));
    std::int64_t valid = 0;
    for (int r = 0; r < rows; ++r)
    {
        if (GDALRasterIO(band, GF_Read, 0, r, cols, 1, line.data(), cols, 1,
                         GDT_Float32, 0, 0) != CE_None) return 0.0;
        for (float v : line)
            if (!has_nodata || v != static_cast<float>(nodata)) ++valid;
    }
    const double total = static_cast<double>(cols) * rows;
    return total > 0.0 ? 100.0 * static_cast<double>(valid) / total : 0.0;
}

// A declared vertical unit that is not metres is refused: elevations must be
// metres, and a feet DEM is otherwise indistinguishable from steep terrain.
bool vertical_units_are_metres(GDALDatasetH ds, std::string& found)
{
    auto is_metre = [](std::string u) {
        for (char& c : u) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return u.empty() || u == "m" || u == "metre" || u == "meter" ||
               u == "metres" || u == "meters";
    };
    const char* band_unit = GDALGetRasterUnitType(GDALGetRasterBand(ds, 1));
    found = band_unit ? band_unit : "";
    if (!is_metre(found)) return false;

    OGRSpatialReferenceH srs = GDALGetSpatialRef(ds);
    if (srs != nullptr && (OSRIsCompound(srs) || OSRIsVertical(srs)))
    {
        char* name = nullptr;
        const double factor = OSRGetTargetLinearUnits(srs, "VERT_CS", &name);
        if (name != nullptr && factor > 0.0 && std::abs(factor - 1.0) > 1e-9)
        {
            found = name;
            return false;
        }
    }
    return true;
}

// Fraction (percent) of a 64 x 64 lattice over the box covered by any item
// bbox: the partial-coverage check for sources whose footprints are known.
double footprint_coverage_percent(const lsdtt_fetch_options& o,
                                  const std::vector<fetch_detail::stac_item>& items)
{
    const int n = 64;
    int covered = 0;
    for (int j = 0; j < n; ++j)
    {
        const double lat = o.south + (j + 0.5) * (o.north - o.south) / n;
        for (int i = 0; i < n; ++i)
        {
            const double lon = o.west + (i + 0.5) * (o.east - o.west) / n;
            for (const auto& item : items)
            {
                if (lon >= item.bbox[0] && lon <= item.bbox[2] &&
                    lat >= item.bbox[1] && lat <= item.bbox[3]) { ++covered; break; }
            }
        }
    }
    return 100.0 * covered / (n * n);
}

lsdtt_fetch_result warp_catalogue_files(const lsdtt_fetch_options&      options,
                                        const lsdtt_dem_source_info&    info,
                                        double                          res,
                                        const std::string&              model,
                                        const std::vector<std::string>& paths,
                                        const std::string&              out_path,
                                        lsdtt_fetch_result              result,
                                        double                          coverage);

// One WCS GetCoverage request: the body, capped at max_bytes (0 = no cap).
// An OGC exception report, whatever the HTTP status, becomes the error.
bool wcs_download(const std::string& url, int timeout_seconds, int max_retry,
                  double max_bytes, std::string& body, std::string& content_type,
                  bool& over_limit, std::string& error)
{
    gdal_string_list_guard http_opts;
    http_opts.add("TIMEOUT="   + std::to_string(timeout_seconds));
    http_opts.add("MAX_RETRY=" + std::to_string(max_retry));
    http_opts.add("RETRY_DELAY=2");
    capped_download download;
    download.max_bytes = max_bytes;
    CPLHTTPResult* response = nullptr;
    {
        quiet_gdal_errors quiet;
        response = CPLHTTPFetchEx(url.c_str(), http_opts.get(), nullptr, nullptr,
                                  capped_download_write, &download);
    }
    over_limit = download.over_limit;
    body.assign(download.data.begin(), download.data.end());
    content_type = (response && response->pszContentType) ? response->pszContentType : "";
    const int status = response ? response->nStatus : -1;
    const std::string transport = (response && response->pszErrBuf) ? response->pszErrBuf : "";
    if (response != nullptr) destroy_detached(response);
    if (over_limit) return false;

    const std::string exception = fetch_detail::parse_ows_exception(
        body.substr(0, std::min<std::size_t>(body.size(), 8192)));
    if (!exception.empty())
    {
        error = "the service refused the request: " + exception;
        return false;
    }
    if (response == nullptr || status != 0)
    {
        error = "HTTP request failed: " + (transport.empty() ? std::string("no response")
                                                             : transport);
        return false;
    }
    if (body.empty())
    {
        error = "the service returned an empty response";
        return false;
    }
    return true;
}

// WCS sources: plan native-grid chunks, download each, then the common warp.
lsdtt_fetch_result fetch_wcs_source(const lsdtt_fetch_options&   options,
                                    const lsdtt_dem_source_info& info,
                                    double                       res,
                                    const std::string&           model,
                                    const std::string&           out_path,
                                    lsdtt_fetch_result           result)
{
    std::string coverage_id;
    for (const auto& entry : info.wcs_coverage_by_model)
        if (entry.first == model) coverage_id = entry.second;

    std::vector<fetch_detail::fetch_extent> chunks;
    try
    {
        chunks = fetch_detail::wcs_plan_chunks(
            fetch_detail::wcs_native_extent(options, info,
                                            options.use_frame ? 0 : result.epsg),
            info.wcs_pixel_metres, info.wcs_max_chunk_side);
    }
    catch (const std::exception& e)
    {
        return make_error(std::string("lsdtt::fetch_dem: ") + e.what());
    }
    if (options.max_tiles > 0 && static_cast<int>(chunks.size()) > options.max_tiles)
    {
        return make_error("lsdtt::fetch_dem: this area needs " +
                          std::to_string(chunks.size()) + " " + info.name +
                          " requests, more than max_tiles (" +
                          std::to_string(options.max_tiles) + "); raise max_tiles or "
                          "fetch a smaller area");
    }

    // Fraction of the area inside the coverage envelopes (coarse).
    std::vector<fetch_detail::stac_item> envelopes;
    for (const auto& box : info.coverage_envelopes)
    {
        fetch_detail::stac_item item;
        for (int k = 0; k < 4; ++k) item.bbox[k] = box[k];
        envelopes.push_back(item);
    }
    const double coverage = footprint_coverage_percent(options, envelopes);

    std::vector<std::string> paths;
    auto remove_chunks = [&paths]() {
        for (const std::string& path : paths) VSIUnlink(path.c_str());
    };
    double bytes_so_far = 0.0;
    const double max_bytes = options.max_download_megabytes * 1.0e6;
    for (std::size_t i = 0; i < chunks.size(); ++i)
    {
        if (i > 0 && options.request_delay_milliseconds > 0)
            CPLSleep(options.request_delay_milliseconds / 1000.0);
        const std::string url = fetch_detail::wcs_getcoverage_url(info, coverage_id, chunks[i]);
        std::string body, content_type, error;
        bool over_limit = false;
        const double remaining = max_bytes > 0.0 ? std::max(1.0, max_bytes - bytes_so_far) : 0.0;
        const bool ok = wcs_download(url, options.http_timeout_seconds, options.http_max_retry,
                                     remaining, body, content_type, over_limit, error);
        const std::string which = info.name + " request " + std::to_string(i + 1) + " of " +
                                  std::to_string(chunks.size());
        if (over_limit)
        {
            remove_chunks();
            return make_error("lsdtt::fetch_dem: " + which + " passed max_download_megabytes (" +
                              format_double(options.max_download_megabytes, 1) +
                              " MB) and the fetch was stopped; raise max_download_megabytes "
                              "or fetch a smaller area");
        }
        if (!ok)
        {
            remove_chunks();
            // A refusal is the service's answer, not the network's fault.
            if (error.rfind("the service refused", 0) == 0)
                return make_error("lsdtt::fetch_dem: " + which + ": " + error);
            return make_skip("lsdtt::fetch_dem: " + which + ": " + error);
        }
        bytes_so_far += static_cast<double>(body.size());
        const std::string chunk_path =
            "/vsimem/lsdtt_wcs_" + info.name + "_" + std::to_string(i) + ".tif";
        if (!fetch_detail::wcs_stage_response(info, body, content_type, chunks[i], chunk_path,
                                                error))
        {
            remove_chunks();
            return make_error("lsdtt::fetch_dem: " + which + ": " + error);
        }
        paths.push_back(chunk_path);
        result.tiles.push_back(coverage_id + "[" + std::to_string(i + 1) + "/" +
                               std::to_string(chunks.size()) + "]");
    }
    if (coverage < 99.95)
    {
        result.warnings.push_back(info.name + " covers about " + format_double(coverage, 1) +
                                  " % of the requested area (" + info.coverage_text + ")");
    }
    result.provenance.push_back({"LSDTT_DOWNLOAD_MEGABYTES",
                                 format_double(bytes_so_far / 1.0e6, 1)});
    lsdtt_fetch_result written = warp_catalogue_files(options, info, res, model, paths,
                                                      out_path, result, coverage);
    remove_chunks();
    return written;
}

// XYZ elevation tiles (GSI): fetch each tile of the layer, decode it to a
// Float32 GeoTIFF on its Web Mercator bounds (a 404 means no data there),
// then the common warp. No fallback to a coarser layer.
lsdtt_fetch_result fetch_xyz_source(const lsdtt_fetch_options&   options,
                                    const lsdtt_dem_source_info& info,
                                    double                       res,
                                    const std::string&           model,
                                    const std::string&           out_path,
                                    lsdtt_fetch_result           result)
{
    const lsdtt_dem_source_info::xyz_layer* layer = nullptr;
    for (const auto& l : info.xyz_layers)
        if (std::abs(l.resolution_metres - res) <= 1e-6 * res) layer = &l;
    if (layer == nullptr)
        return make_error("lsdtt::fetch_dem: " + info.name + " has no layer at " +
                          list_of({res}) + " m");

    std::vector<lsdtt_slippy_tile> tiles;
    try
    {
        const fetch_detail::fetch_extent box =
            fetch_detail::output_area_wgs84(options, options.use_frame ? 0 : result.epsg);
        tiles = slippy_tiles_for_bbox(box.x_min, box.y_min, box.x_max, box.y_max, layer->zoom);
    }
    catch (const std::exception& e)
    {
        return make_error(std::string("lsdtt::fetch_dem: ") + e.what());
    }
    if (options.max_tiles > 0 && static_cast<int>(tiles.size()) > options.max_tiles)
    {
        return make_error("lsdtt::fetch_dem: this area needs " + std::to_string(tiles.size()) +
                          " " + layer->label + " tiles, more than max_tiles (" +
                          std::to_string(options.max_tiles) + "); raise max_tiles, use a "
                          "coarser source_resolution_metres, or fetch a smaller area");
    }

    std::vector<std::string> paths;
    auto remove_tiles = [&paths]() {
        for (const std::string& path : paths) VSIUnlink(path.c_str());
    };
    double bytes_so_far = 0.0;
    const double max_bytes = options.max_download_megabytes * 1.0e6;
    std::size_t n_missing = 0;
    for (std::size_t i = 0; i < tiles.size(); ++i)
    {
        if (i > 0 && options.request_delay_milliseconds > 0)
            CPLSleep(options.request_delay_milliseconds / 1000.0);
        const lsdtt_slippy_tile& t = tiles[i];
        std::string url = layer->url_template;
        for (const auto& token : {std::make_pair(std::string("{z}"), t.z),
                                  std::make_pair(std::string("{x}"), t.x),
                                  std::make_pair(std::string("{y}"), t.y)})
        {
            const std::size_t at = url.find(token.first);
            if (at != std::string::npos) url.replace(at, token.first.size(), std::to_string(token.second));
        }
        std::string body, content_type, error;
        bool over_limit = false;
        const double remaining = max_bytes > 0.0 ? std::max(1.0, max_bytes - bytes_so_far) : 0.0;
        const bool ok = wcs_download(url, options.http_timeout_seconds, options.http_max_retry,
                                     remaining, body, content_type, over_limit, error);
        const std::string tile_name = layer->label + "/" + std::to_string(t.z) + "/" +
                                      std::to_string(t.x) + "/" + std::to_string(t.y);
        if (over_limit)
        {
            remove_tiles();
            return make_error("lsdtt::fetch_dem: " + info.name + " passed "
                              "max_download_megabytes (" +
                              format_double(options.max_download_megabytes, 1) +
                              " MB) and the fetch was stopped");
        }
        // A missing tile is an HTTP 404 with S3's NoSuchKey XML; GDAL may hand
        // the body over as a successful transfer.
        const bool not_image = content_type.find("image/") == std::string::npos;
        if ((ok && not_image && body.find("NoSuchKey") != std::string::npos) ||
            (!ok && error.find("404") != std::string::npos))
        {
            ++n_missing;   // no data in this tile
            continue;
        }
        if (ok && not_image)
        {
            remove_tiles();
            return make_error("lsdtt::fetch_dem: " + info.name + " tile " + tile_name +
                              ": the service returned '" + content_type + "', not a PNG: " +
                              body.substr(0, 200));
        }
        if (!ok)
        {
            remove_tiles();
            return make_skip("lsdtt::fetch_dem: " + info.name + " tile " + tile_name + ": " +
                             error);
        }
        bytes_so_far += static_cast<double>(body.size());
        const std::string tile_path =
            "/vsimem/lsdtt_xyz_" + info.name + "_" + std::to_string(i) + ".tif";
        if (!fetch_detail::gsi_tile_to_geotiff(body, t.z, t.x, t.y, tile_path, error))
        {
            remove_tiles();
            return make_error("lsdtt::fetch_dem: " + info.name + " tile " + tile_name + ": " +
                              error);
        }
        paths.push_back(tile_path);
        result.tiles.push_back(tile_name);
    }
    if (paths.empty())
    {
        return make_error("lsdtt::fetch_dem: no " + layer->label + " tiles exist in this area (" +
                          std::to_string(tiles.size()) + " checked)" +
                          (res < 10.0 ? std::string("; the 10 m layer covers all of Japan "
                                                    "(source_resolution_metres: 10)")
                                      : std::string()));
    }
    const double coverage = 100.0 * static_cast<double>(paths.size()) /
                            static_cast<double>(tiles.size());
    if (n_missing > 0)
    {
        result.warnings.push_back(std::to_string(n_missing) + " of " +
                                  std::to_string(tiles.size()) + " " + layer->label +
                                  " tiles have no data (no fallback to a coarser layer)");
    }
    result.provenance.push_back({"LSDTT_SOURCE_COLLECTIONS", layer->label});
    result.provenance.push_back({"LSDTT_DOWNLOAD_MEGABYTES", format_double(bytes_so_far / 1.0e6, 1)});
    lsdtt_fetch_result written = warp_catalogue_files(options, info, res, model, paths,
                                                      out_path, result, coverage);
    remove_tiles();
    return written;
}

// Grid-keyed buckets: compute the candidate file names, keep those that
// exist (a missing file means no data there), pick one collection, then the
// common warp.
lsdtt_fetch_result fetch_grid_source(const lsdtt_fetch_options&   options,
                                     const lsdtt_dem_source_info& info,
                                     double                       res,
                                     const std::string&           model,
                                     const std::string&           out_path,
                                     lsdtt_fetch_result           result)
{
    std::string bucket_model = model;
    for (const auto& entry : info.grid_model_names)
        if (entry.first == model) bucket_model = entry.second;

    const std::string choice = options.source_collection_choice.empty()
                             ? std::string("latest") : options.source_collection_choice;
    std::vector<std::pair<std::string, std::string>> candidates;   // (collection, key)
    try
    {
        const fetch_detail::fetch_extent area = fetch_detail::grid_search_extent(
            options, info, options.use_frame ? 0 : result.epsg);
        for (const auto& collection : info.grid_collections)
            if (choice == "latest" || collection.label == choice)
            for (const std::string& key :
                 fetch_detail::grid_keys_for_extent(collection, area, bucket_model))
                candidates.emplace_back(collection.label, key);
    }
    catch (const std::exception& e)
    {
        return make_error(std::string("lsdtt::fetch_dem: ") + e.what());
    }
    if (options.max_tiles > 0 && static_cast<int>(candidates.size()) > options.max_tiles)
    {
        return make_error("lsdtt::fetch_dem: this area needs " +
                          std::to_string(candidates.size()) + " " + info.name +
                          " tiles checked, more than max_tiles (" +
                          std::to_string(options.max_tiles) + "); raise max_tiles or "
                          "fetch a smaller area");
    }

    // Which collections have tiles here: the named one, or all of them (for
    // "latest", the newest with any tile here is used).
    std::vector<std::pair<std::string, std::vector<std::string>>> found;   // oldest first
    for (const auto& collection : info.grid_collections)
    {
        if (choice != "latest" && collection.label != choice) continue;
        std::vector<std::string> existing;
        for (const auto& candidate : candidates)
        {
            if (candidate.first != collection.label) continue;
            const std::string path = "/vsicurl/" + info.grid_base_url + candidate.second;
            VSIStatBufL stat;
            bool exists = false;
            {
                quiet_gdal_errors quiet;
                exists = VSIStatExL(path.c_str(), &stat, VSI_STAT_EXISTS_FLAG) == 0;
            }
            if (exists) existing.push_back(candidate.second);
        }
        if (!existing.empty()) found.emplace_back(collection.label, existing);
    }
    std::vector<std::string> paths;
    std::string used;
    if (!found.empty())
    {
        used = found.back().first;
        for (const std::string& key : found.back().second)
        {
            paths.push_back("/vsicurl/" + info.grid_base_url + key);
            result.tiles.push_back(key.substr(key.rfind('/') + 1));
        }
    }
    if (paths.empty())
    {
        VSIStatBufL stat;
        bool reachable = false;
        {
            quiet_gdal_errors quiet;
            reachable = VSIStatExL(("/vsicurl/" + info.grid_base_url + info.grid_sentinel_key).c_str(),
                                   &stat, VSI_STAT_EXISTS_FLAG) == 0;
        }
        if (!reachable)
            return make_skip("lsdtt::fetch_dem: the " + info.name + " bucket is unreachable");
        return make_error("lsdtt::fetch_dem: no " + info.name + " " + model + " tiles " +
                          (choice == "latest" ? std::string() : "of " + choice + " ") +
                          "exist in this area (" + std::to_string(candidates.size()) +
                          " checked; coverage: " + info.coverage_text + ")");
    }

    std::vector<fetch_detail::stac_item> envelopes;
    for (const auto& box : info.coverage_envelopes)
    {
        fetch_detail::stac_item item;
        for (int k = 0; k < 4; ++k) item.bbox[k] = box[k];
        envelopes.push_back(item);
    }
    const double coverage = footprint_coverage_percent(options, envelopes);
    if (info.grid_collections.size() > 1)
    {
        result.provenance.push_back({"LSDTT_SOURCE_COLLECTIONS", used});
        std::string others;
        for (std::size_t k = 0; k + 1 < found.size(); ++k)
            others += (others.empty() ? "" : ", ") + found[k].first;
        if (!others.empty())
            result.warnings.push_back(info.name + " used " + used + " (the newest collection "
                                      "with tiles here); " + others + " also cover parts of "
                                      "this area: set source_collection_choice to use one");
    }
    return warp_catalogue_files(options, info, res, model, paths, out_path, result, coverage);
}

// The catalogue path: discover the source files, warp straight from them onto
// the snapped target grid, check the result, stamp its provenance.
lsdtt_fetch_result fetch_catalogue_source(const lsdtt_fetch_options& options,
                                          const std::string&         out_path,
                                          lsdtt_fetch_result         result)
{
    const lsdtt_dem_source_info& info = dem_source_info(options.source);
    const double res = resolved_source_resolution(options, info);
    const std::string model = resolved_surface_model(options, info);
    if (!info.wcs_url.empty())
        return fetch_wcs_source(options, info, res, model, out_path, result);
    if (!info.grid_base_url.empty())
        return fetch_grid_source(options, info, res, model, out_path, result);
    if (!info.xyz_layers.empty())
        return fetch_xyz_source(options, info, res, model, out_path, result);

    // ---- discovery: STAC search, following rel=next ----
    std::vector<fetch_detail::stac_item> items;
    std::vector<std::string> missing_asset;
    {
        std::string asset_key = info.stac_asset_by_gsd ? std::string() : info.stac_asset;
        for (const auto& entry : info.stac_asset_by_model)
            if (entry.first == model) asset_key = entry.second;
        std::string url  = info.stac_search_url;
        const std::string first_body = fetch_detail::catalogue_stac_search_body(options);
        std::string body = first_body;
        for (;;)    // one STAC page per pass, until there is no "next" link
        {
            std::string response, error, next_href, next_body;
            if (!http_request_body(url, body, options.http_timeout_seconds,
                                   options.http_max_retry, response, error))
            {
                return make_skip("lsdtt::fetch_dem: " + info.name + " search failed: " + error);
            }
            if (!fetch_detail::parse_stac_items(response, asset_key, items, missing_asset,
                                                next_href, next_body, error, res))
            {
                return make_skip("lsdtt::fetch_dem: " + info.name + ": " + error);
            }
            if (next_href.empty()) break;
            if (options.max_tiles > 0 &&
                static_cast<int>(items.size()) > options.max_tiles * 2)
                break;   // already far past the limit; the count check below refuses
            url = next_href;
            body = fetch_detail::merge_stac_next_body(first_body, next_body);
        }
    }

    // A host that ignored the filter must not pull in far-away tiles.
    items = fetch_detail::stac_items_in_bbox(items, options.west, options.south,
                                             options.east, options.north);
    const std::vector<double> gsds = fetch_detail::distinct_gsds(items);
    if (info.stac_items_carry_gsd) items = fetch_detail::filter_stac_items_by_gsd(items, res);
    if (info.stac_latest_per_tile) items = fetch_detail::latest_stac_item_per_tile(items);
    if (items.empty())
    {
        return make_error("lsdtt::fetch_dem: no " + info.name + " data at " +
                          list_of({res}) + " m intersects this area (resolutions "
                          "found: " + (gsds.empty() ? std::string("none") : list_of(gsds)) +
                          " m; coverage: " + info.coverage_text + ")");
    }
    if (options.max_tiles > 0 && static_cast<int>(items.size()) > options.max_tiles)
    {
        return make_error("lsdtt::fetch_dem: this area needs " +
                          std::to_string(items.size()) + " " + info.name +
                          " source files, more than max_tiles (" +
                          std::to_string(options.max_tiles) + "); raise max_tiles or "
                          "fetch a smaller area");
    }
    fetch_detail::sort_stac_items(items);

    const double coverage = footprint_coverage_percent(options, items);
    if (coverage < 99.95)
    {
        result.warnings.push_back(info.name + " data covers " + format_double(coverage, 1) +
                                  " % of the requested area; the rest is NoData");
    }

    // ---- sign and open the sources ----
    std::string token;
    if (info.stac_planetary_computer_signing)
    {
        std::string error;
        if (!pc_collection_token(info.stac_collection, options.http_timeout_seconds,
                                 options.http_max_retry, token, error))
        {
            return make_skip("lsdtt::fetch_dem: " + info.name + ": " + error);
        }
    }
    std::vector<std::string> paths;
    for (const auto& item : items)
    {
        std::string href = item.href;
        if (!token.empty())
            href += (href.find('?') == std::string::npos ? "?" : "&") + token;
        paths.push_back("/vsicurl/" + href);
        result.tiles.push_back(item.id);
    }
    return warp_catalogue_files(options, info, res, model, paths, out_path, result,
                                coverage);
}

// After discovery: open, check units, warp onto the snapped grid, check the
// output, stamp its provenance. Shared with the offline tests through
// fetch_detail::warp_catalogue_files.
lsdtt_fetch_result warp_catalogue_files(const lsdtt_fetch_options&      options,
                                        const lsdtt_dem_source_info&    info,
                                        double                          res,
                                        const std::string&              model,
                                        const std::vector<std::string>& paths,
                                        const std::string&              out_path,
                                        lsdtt_fetch_result              result,
                                        double                          coverage)
{
    fetch_detail::fetch_extent extent;
    if (!options.use_frame)
    {
        try
        {
            extent = fetch_detail::target_extent_for_bbox(options.west, options.south,
                                                          options.east, options.north,
                                                          result.epsg, options.grid_spacing);
        }
        catch (const std::exception& e)
        {
            return make_error(std::string("lsdtt::fetch_dem: ") + e.what());
        }
    }
    // UTM, the default output CRS, is not defined beyond 84 N and 80 S.
    const double centre_lat = 0.5 * (options.south + options.north);
    if (!options.use_frame && options.target_epsg <= 0 &&
        (centre_lat > 84.0 || centre_lat < -80.0))
    {
        result.warnings.push_back(
            "the area is centred beyond the UTM latitude limits (84 N / 80 S), so "
            "the default UTM output is strongly distorted; set target_epsg to 3413 "
            "(Arctic) or 3031 (Antarctic) for a polar stereographic grid");
    }

    const std::vector<std::string>& item_names = result.tiles;

    std::vector<GDALDatasetH> datasets;
    auto close_all = [&datasets]() {
        for (GDALDatasetH d : datasets) GDALClose(d);
        datasets.clear();
    };
    for (std::size_t i = 0; i < paths.size(); ++i)
    {
        GDALDatasetH d = nullptr;
        {
            quiet_gdal_errors quiet;
            d = GDALOpen(paths[i].c_str(), GA_ReadOnly);
        }
        const std::string name = i < item_names.size() ? item_names[i] : paths[i];
        if (d == nullptr)
        {
            close_all();
            return make_skip("lsdtt::fetch_dem: could not open " + info.name +
                             " source file " + name + " -- the host may be unreachable");
        }
        std::string unit;
        if (!vertical_units_are_metres(d, unit))
        {
            GDALClose(d);
            close_all();
            return make_error("lsdtt::fetch_dem: " + info.name + " source file " +
                              name + " declares vertical units '" + unit +
                              "'; elevations must be in metres");
        }
        datasets.push_back(d);
    }

    // ---- overviews: GDALWarp reads them when the grid is coarse enough ----
    const int n_overviews = GDALGetOverviewCount(GDALGetRasterBand(datasets.front(), 1));
    const bool overview_used = n_overviews > 0 && options.grid_spacing >= 2.0 * res - 1e-9;
    if (overview_used)
    {
        result.warnings.push_back(
            "grid_spacing " + list_of({options.grid_spacing}) + " m is coarser than the " +
            list_of({res}) + " m source; GDAL will read the source's overviews (built "
            "by the provider, kernel not controlled), and " + options.resampling +
            " is applied only to the final step");
    }

    lsdtt_fetch_result written = warp_and_write_sources(
        datasets, options, result, out_path, options.use_frame ? nullptr : &extent);
    close_all();
    if (!written.ok) return written;

    const double valid = valid_percent(written.path);
    if (!(valid > 0.0))
    {
        VSIUnlink(written.path.c_str());
        return make_error("lsdtt::fetch_dem: the " + info.name + " output has no valid "
                          "cells; the area may lie in a gap of its coverage (" +
                          info.coverage_text + ")");
    }

    if (valid < 99.95)
    {
        written.warnings.push_back(format_double(100.0 - valid, 1) + " % of the output "
                                   "cells are NoData (outside " + info.name +
                                   " coverage, or gaps in the source)");
    }

    std::string surface = model;
    for (char& c : surface) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    std::vector<std::pair<std::string, std::string>> items_out = {
        {"LSDTT_DEM_SOURCE", info.name},
        {"LSDTT_SURFACE_MODEL", surface},
        {"LSDTT_VERTICAL_DATUM", info.vertical_datum},
        {"LSDTT_VERTICAL_UNITS", "metre"},
        {"LSDTT_SOURCE_RESOLUTION_METRES", list_of({res})},
        {"LSDTT_SOURCE_ITEM_COUNT", std::to_string(paths.size())},
        {"LSDTT_RESAMPLING", options.resampling},
        {"LSDTT_OVERVIEW_USED", overview_used ? "true" : "false"},
        {"LSDTT_COVERAGE_PERCENT", format_double(coverage, 1)},
        {"LSDTT_VALID_PERCENT", format_double(valid, 1)}};
    // Items the discovery step recorded (e.g. the WCS download size).
    for (const auto& item : written.provenance) items_out.push_back(item);
    return stamp_result(written, items_out);
}

lsdtt_fetch_result fetch_dem_impl(const lsdtt_fetch_options& options,
                                  const std::string&         out_path)
{
    std::string error;
    if (!validate_options(options, error))
    {
        return make_error("lsdtt::fetch_dem: " + error);
    }

    if (dem_source_name(options.source) == "unknown")
    {
        // Explicit rather than silently substituting a different source.
        return make_error(
            "lsdtt::fetch_dem: source is not supported (available: " +
            dem_source_name_list() + ")");
    }
    if (!preflight_source(options, error))
    {
        return make_error("lsdtt::fetch_dem: " + error);
    }

    if (options.source == lsdtt_dem_source::opentopography &&
        options.api_key.empty())
    {
        return make_error(
            "lsdtt::fetch_dem: the opentopography source needs an API key. "
            "Set OT_API_KEY, or set OT_API_KEY_FILE to a file containing it. "
            "(cop30_aws needs no key and gives the same global 30 m coverage.)");
    }

    // Caught here rather than at the API, which answers a bad demtype with a
    // bare 400 that says nothing about what would have been valid.
    if (options.source == lsdtt_dem_source::opentopography &&
        !options.ot_dem_type.empty() &&
        ot_dataset_info(options.ot_dem_type) == nullptr)
    {
        return make_error(
            "lsdtt::fetch_dem: unknown ot_dem_type '" + options.ot_dem_type +
            "'. Available datasets: " + ot_dataset_demtype_list() +
            " (run lsdtt-fetch-raster -h ot_dem_type for descriptions)");
    }

    if (out_path.empty())
    {
        return make_error("lsdtt::fetch_dem: out_path must not be empty");
    }

    lsdtt_fetch_result result;
    if (options.source == lsdtt_dem_source::cop30_aws)
    {
        result.tiles = cop30_tiles_for_bbox(options.west, options.south,
                                            options.east, options.north);
    }
    else if (options.source == lsdtt_dem_source::opentopography)
    {
        // OpenTopography serves an arbitrary window; there are no tiles to name.
        result.tiles.push_back(options.ot_dem_type.empty() ? "COP30"
                                                           : options.ot_dem_type);
    }
    // Planetary Computer tiles are discovered by the STAC search below.
    if (options.use_frame)
    {
        result.epsg = 0;
        result.in_frame = true;
        try
        {
            result.rotate_back_epsg = (options.target_epsg > 0)
                                    ? options.target_epsg
                                    : analysis_epsg_for_frame_centre(options.frame);
        }
        catch (const std::exception& e)
        {
            return make_error(std::string("lsdtt::fetch_dem: ") + e.what());
        }
    }
    else
    {
        result.epsg  = (options.target_epsg > 0)
                     ? options.target_epsg
                     : utm_epsg_for_point(0.5 * (options.west + options.east),
                                          0.5 * (options.south + options.north));
    }
    result.grid_spacing = options.grid_spacing;

    if (!preflight_size(options, options.use_frame ? 0 : result.epsg, error))
    {
        return make_error("lsdtt::fetch_dem: " + error);
    }

    if (!http_available())
    {
        return make_skip(
            "lsdtt::fetch_dem: this GDAL build has no HTTP support (built "
            "without curl), so remote DEMs cannot be fetched");
    }


    // ---------------------------------------------------------------------
    // Cache lookup (opt-in)
    // ---------------------------------------------------------------------
    std::string cached_path;
    if (!options.cache_dir.empty())
    {
        cached_path = options.cache_dir + "/" + cache_key_for(options);
        if (file_exists(cached_path))
        {
            gdal_dataset_guard cached;
            {
                quiet_gdal_errors quiet;
                cached.reset(GDALOpen(cached_path.c_str(), GA_ReadOnly));
            }
            if (cached)
            {
                result.n_cols = GDALGetRasterXSize(cached.get());
                result.n_rows = GDALGetRasterYSize(cached.get());
                // The cached raster carries the provenance stamped when it
                // was fetched.
                // CSLConstList, not char**: GDAL >= 3.13 returns a const list
                // here (older GDAL's char** converts to it implicitly).
                CSLConstList md = GDALGetMetadata(cached.get(), nullptr);
                for (int i = 0; md != nullptr && md[i] != nullptr; ++i)
                {
                    const std::string entry = md[i];
                    const std::size_t eq = entry.find('=');
                    if (entry.rfind("LSDTT_", 0) == 0 && eq != std::string::npos)
                        result.provenance.emplace_back(entry.substr(0, eq), entry.substr(eq + 1));
                }
                cached.reset();

                if (cached_path != out_path)
                {
                    if (CPLCopyFile(out_path.c_str(), cached_path.c_str()) != 0)
                    {
                        return make_error(
                            "lsdtt::fetch_dem: could not copy cached raster to " +
                            out_path);
                    }
                }
                result.ok         = true;
                result.from_cache = true;
                result.path       = out_path;
                result.message    = "reused cached raster " + cached_path;
                return result;
            }
            // A corrupt cache entry is not fatal: fall through and re-fetch.
        }
    }

    http_config_guard http_guard(options.http_timeout_seconds,
                                 options.http_connect_timeout_seconds,
                                 options.http_max_retry);

    // ---------------------------------------------------------------------
    // Catalogue sources: their own discovery and warp (no -projwin cut).
    // ---------------------------------------------------------------------
    if (!dem_source_info(options.source).legacy_path)
    {
        const lsdtt_fetch_result catalogue_result =
            fetch_catalogue_source(options, out_path, result);
        if (catalogue_result.ok && !cached_path.empty() && cached_path != out_path)
        {
            quiet_gdal_errors quiet;
            CPLCopyFile(cached_path.c_str(), out_path.c_str());
        }
        return catalogue_result;
    }

    // ---------------------------------------------------------------------
    // Step 1: open the source (mosaicking when the bbox spans several tiles)
    //
    // A bounding box crossing a 1-degree boundary is normal for real study
    // areas, so mosaic rather than failing.
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // OpenTopography: the service performs the windowing, so the download is
    // already the requested extent and step 2 (projwin) is skipped.
    // ---------------------------------------------------------------------
    if (options.source == lsdtt_dem_source::opentopography)
    {
        const std::string url = opentopography_url(options);
        const std::string staged_path = "/vsimem/lsdtt_fetch_ot.tif";

        std::string download_error;
        if (!download_to_vsimem(url, options.api_key, staged_path,
                                options.http_timeout_seconds,
                                options.http_max_retry,
                                options.max_download_megabytes,
                                download_error))
        {
            // download_error is already masked. Never interpolate `url` here.
            // Passing the byte limit is the user's limit, not the network.
            if (download_error.find("max_download_megabytes") != std::string::npos)
                return make_error("lsdtt::fetch_dem: " + download_error);
            return make_skip("lsdtt::fetch_dem: " + download_error);
        }

        gdal_dataset_guard downloaded;
        {
            quiet_gdal_errors quiet;
            downloaded.reset(GDALOpen(staged_path.c_str(), GA_ReadOnly));
        }
        if (!downloaded)
        {
            VSIUnlink(staged_path.c_str());
            return make_error("lsdtt::fetch_dem: downloaded raster is unreadable");
        }

        lsdtt_fetch_result warped_result =
            warp_and_write(downloaded, options, result, out_path);

        downloaded.reset();
        VSIUnlink(staged_path.c_str());
        warped_result = stamp_result(warped_result, legacy_provenance(options, warped_result));

        if (warped_result.ok && !cached_path.empty() && cached_path != out_path)
        {
            quiet_gdal_errors quiet;
            CPLCopyFile(cached_path.c_str(), out_path.c_str());
        }
        return warped_result;
    }

    // ---------------------------------------------------------------------
    // Planetary Computer: discover items via STAC, sign each asset, then hand
    // the signed URLs to the same windowed-read path as the AWS source.
    // ---------------------------------------------------------------------
    std::vector<std::string> pc_urls;
    if (options.source == lsdtt_dem_source::planetary_computer)
    {
        std::vector<std::string> signed_hrefs;
        std::string search_error;
        bool search_is_error = false;

        if (!pc_search_assets(options, signed_hrefs, result.tiles, result.warnings,
                              search_is_error, search_error))
        {
            if (search_is_error)
                return make_error("lsdtt::fetch_dem: " + search_error);
            // Any failure here is "data not obtainable right now" from the
            // caller's point of view, so skip rather than fail: an offline
            // machine and an empty search look the same to a test.
            return make_skip("lsdtt::fetch_dem: " + search_error);
        }

        pc_urls.reserve(signed_hrefs.size());
        for (const std::string& href : signed_hrefs)
        {
            pc_urls.push_back("/vsicurl/" + href);
        }
    }

    gdal_dataset_guard source;
    {
        quiet_gdal_errors quiet;   // an offline probe must not spray stderr

        if (!pc_urls.empty())
        {
            if (pc_urls.size() == 1)
            {
                source.reset(GDALOpen(pc_urls[0].c_str(), GA_ReadOnly));
            }
            else
            {
                std::vector<const char*> url_ptrs;
                url_ptrs.reserve(pc_urls.size());
                for (const std::string& url : pc_urls) url_ptrs.push_back(url.c_str());

                source.reset(GDALBuildVRT(
                    "/vsimem/lsdtt_fetch_mosaic.vrt",
                    static_cast<int>(url_ptrs.size()),
                    nullptr, url_ptrs.data(), nullptr, nullptr));
            }
        }
        else if (result.tiles.size() == 1)
        {
            source.reset(GDALOpen(cop30_tile_url(result.tiles[0]).c_str(),
                                  GA_ReadOnly));
        }
        else
        {
            std::vector<std::string> urls;
            std::vector<const char*> url_ptrs;
            urls.reserve(result.tiles.size());
            url_ptrs.reserve(result.tiles.size());
            for (const std::string& tile : result.tiles)
            {
                urls.push_back(cop30_tile_url(tile));
            }
            for (const std::string& url : urls) url_ptrs.push_back(url.c_str());

            source.reset(GDALBuildVRT(
                "/vsimem/lsdtt_fetch_mosaic.vrt",
                static_cast<int>(url_ptrs.size()),
                nullptr,
                url_ptrs.data(),
                nullptr,
                nullptr));
        }
    }

    if (!source)
    {
        // Cannot distinguish "offline" from "tile genuinely absent" without
        // more probing, and both mean the same thing to a caller: the data is
        // not obtainable right now. Treat as skip so tests degrade rather than
        // fail; the message names the tiles so a real 404 is diagnosable.
        std::ostringstream msg;
        msg << "lsdtt::fetch_dem: could not open source tile(s) [";
        for (std::size_t i = 0; i < result.tiles.size(); ++i)
        {
            if (i) msg << ", ";
            msg << result.tiles[i];
        }
        msg << "] -- the host may be unreachable, or the bounding box may lie "
               "outside Copernicus GLO-30 coverage";
        return make_skip(msg.str());
    }

    // ---------------------------------------------------------------------
    // Step 2: windowed read over /vsicurl (HTTP range requests)
    // ---------------------------------------------------------------------
    gdal_string_list_guard cut_args;
    cut_args.add("-projwin");
    cut_args.add(format_double(options.west,  10));
    cut_args.add(format_double(options.north, 10));
    cut_args.add(format_double(options.east,  10));
    cut_args.add(format_double(options.south, 10));

    GDALTranslateOptions* cut_opts =
        GDALTranslateOptionsNew(cut_args.get(), nullptr);
    if (!cut_opts)
    {
        return make_error("lsdtt::fetch_dem: GDALTranslateOptionsNew failed");
    }

    gdal_dataset_guard cut;
    {
        int usage_error = CE_None;
        cut.reset(GDALTranslate("/vsimem/lsdtt_fetch_cut.tif", source.get(),
                                cut_opts, &usage_error));
        GDALTranslateOptionsFree(cut_opts);
        if (!cut || usage_error != CE_None)
        {
            return make_skip(
                "lsdtt::fetch_dem: windowed read failed -- the connection may "
                "have dropped mid-transfer");
        }
    }
    source.reset();

    // ---------------------------------------------------------------------
    // Steps 3 and 4: warp, optional centre crop, write
    // ---------------------------------------------------------------------
    lsdtt_fetch_result final_result =
        warp_and_write(cut, options, result, out_path);

    cut.reset();
    VSIUnlink("/vsimem/lsdtt_fetch_cut.tif");
    final_result = stamp_result(final_result, legacy_provenance(options, final_result));

    // ---------------------------------------------------------------------
    // Populate the cache
    // ---------------------------------------------------------------------
    if (final_result.ok && !cached_path.empty() && cached_path != out_path)
    {
        // A cache write failure must not fail the fetch: the caller already
        // has the raster they asked for.
        quiet_gdal_errors quiet;
        CPLCopyFile(cached_path.c_str(), out_path.c_str());
    }

    return final_result;
}

} // namespace

bool fetch_detail::wcs_stage_response(const lsdtt_dem_source_info& info,
                                      const std::string&           body,
                                      const std::string&           content_type,
                                      const fetch_extent&          chunk,
                                      const std::string&           out_path,
                                      std::string&                 error)
{
    std::string payload = body;
    std::string staged = out_path + ".download";
    if (info.wcs_multipart_ascii)
    {
        std::string boundary;
        const std::size_t b = content_type.find("boundary=");
        if (b != std::string::npos)
        {
            boundary = content_type.substr(b + 9);
            const std::size_t end = boundary.find(';');
            if (end != std::string::npos) boundary = boundary.substr(0, end);
            if (!boundary.empty() && boundary.front() == '"')
                boundary = boundary.substr(1, boundary.size() - 2);
        }
        payload = boundary.empty() ? body
                                   : fetch_detail::multipart_first_payload(body, boundary);
        if (payload.empty())
        {
            error = "the multipart response held no grid";
            return false;
        }
        staged += ".asc";
    }
    else
    {
        staged += ".tif";
    }
    {
        VSILFILE* f = VSIFOpenL(staged.c_str(), "wb");
        if (f == nullptr ||
            VSIFWriteL(payload.data(), 1, payload.size(), f) != payload.size())
        {
            if (f) VSIFCloseL(f);
            error = "could not stage the response in memory";
            return false;
        }
        VSIFCloseL(f);
    }
    GDALDatasetH src = nullptr;
    {
        quiet_gdal_errors quiet;
        src = GDALOpen(staged.c_str(), GA_ReadOnly);
    }
    if (src == nullptr)
    {
        VSIUnlink(staged.c_str());
        error = "the response is not a readable raster (content type '" + content_type + "')";
        return false;
    }

    double gt[6] = {0, 0, 0, 0, 0, 0};
    GDALGetGeoTransform(src, gt);
    const double p = info.wcs_pixel_metres;
    const double tol = 1e-4 * p;   // origins are given to 1e-6 m
    const int w = static_cast<int>(std::llround((chunk.x_max - chunk.x_min) / p));
    const int h = static_cast<int>(std::llround((chunk.y_max - chunk.y_min) / p));
    if (std::abs(gt[1] - p) > 1e-9 * p || std::abs(gt[5] + p) > 1e-9 * p ||
        gt[2] != 0.0 || gt[4] != 0.0 || std::abs(gt[0] - chunk.x_min) > tol ||
        std::abs(gt[3] - chunk.y_max) > tol ||
        GDALGetRasterXSize(src) != w || GDALGetRasterYSize(src) != h)
    {
        error = "the service returned a " + std::to_string(GDALGetRasterXSize(src)) + " x " +
                std::to_string(GDALGetRasterYSize(src)) + " grid at " +
                list_of({gt[1]}) + " m from (" + list_of({gt[0]}) + ", " +
                list_of({gt[3]}) + "), not the requested native " + std::to_string(w) +
                " x " + std::to_string(h) + " grid at " + list_of({p}) + " m from (" +
                list_of({chunk.x_min}) + ", " + list_of({chunk.y_max}) + ")";
        GDALClose(src);
        VSIUnlink(staged.c_str());
        return false;
    }

    gdal_string_list_guard args;
    for (const char* a : {"-of", "GTiff", "-ot", "Float32"}) args.add(a);
    args.add("-a_srs");
    args.add("EPSG:" + std::to_string(info.wcs_native_epsg));
    if (info.wcs_zero_is_nodata)
    {
        args.add("-a_nodata");
        args.add("0");
    }
    GDALTranslateOptions* opts = GDALTranslateOptionsNew(args.get(), nullptr);
    GDALDatasetH copy = nullptr;
    {
        quiet_gdal_errors quiet;
        copy = GDALTranslate(out_path.c_str(), src, opts, nullptr);
    }
    GDALTranslateOptionsFree(opts);
    GDALClose(src);
    VSIUnlink(staged.c_str());
    VSIUnlink((staged + ".aux.xml").c_str());
    if (copy == nullptr)
    {
        error = "could not convert the response to GeoTIFF";
        return false;
    }
    GDALClose(copy);
    return true;
}


lsdtt_fetch_result fetch_detail::warp_catalogue_files(const lsdtt_fetch_options&      options,
                                                      const std::vector<std::string>& paths,
                                                      const std::string&              out_path)
{
    std::string error;
    if (!validate_options(options, error)) return make_error("lsdtt::fetch_dem: " + error);
    const lsdtt_dem_source_info& info = dem_source_info(options.source);
    lsdtt_fetch_result result;
    result.epsg = options.target_epsg > 0
                ? options.target_epsg
                : utm_epsg_for_point(0.5 * (options.west + options.east),
                                     0.5 * (options.south + options.north));
    result.grid_spacing = options.grid_spacing;
    return ::lsdtt::warp_catalogue_files(options, info,
                                         resolved_source_resolution(options, info),
                                         resolved_surface_model(options, info),
                                         paths, out_path, result, 100.0);
}

// -----------------------------------------------------------------------------
// fetch_bbox_for_frame
// -----------------------------------------------------------------------------
void fetch_bbox_for_frame(const lsdtt_oblique_frame& frame, double pad_metres,
                          double& west, double& south, double& east, double& north)
{
    lsdtt_oblique_frame_projector projector(frame);
    const int per_side = 50;
    west = 1e300; east = -1e300; south = 1e300; north = -1e300;
    for (int side = 0; side < 4; ++side)
    {
        for (int i = 0; i < per_side; ++i)
        {
            const double t = static_cast<double>(i) / per_side;
            double x = 0.0, y = 0.0;
            switch (side)
            {
            case 0: x = frame.x_min + t * (frame.x_max - frame.x_min); y = frame.y_min; break;
            case 1: x = frame.x_max; y = frame.y_min + t * (frame.y_max - frame.y_min); break;
            case 2: x = frame.x_max - t * (frame.x_max - frame.x_min); y = frame.y_max; break;
            default: x = frame.x_min; y = frame.y_max - t * (frame.y_max - frame.y_min); break;
            }
            const lsdtt_lonlat p = projector.to_geographic(x, y);
            west = std::min(west, p.longitude); east = std::max(east, p.longitude);
            south = std::min(south, p.latitude); north = std::max(north, p.latitude);
        }
    }
    if (east - west > 180.0)
        throw std::runtime_error(
            "fetch_bbox_for_frame: the frame crosses the antimeridian, which fetching does not "
            "support");
    const double metres_per_degree = 111320.0;
    const double max_abs_lat = std::min(89.0, std::max(std::abs(south), std::abs(north)));
    const double dlat = pad_metres / metres_per_degree;
    const double dlon =
        pad_metres / (metres_per_degree * std::cos(max_abs_lat * 3.14159265358979323846 / 180.0));
    west = std::max(-180.0, west - dlon);
    east = std::min(180.0, east + dlon);
    south = std::max(-90.0, south - dlat);
    north = std::min(90.0, north + dlat);
}

// =============================================================================
// fetch_dem
// =============================================================================
lsdtt_fetch_result fetch_dem(const lsdtt_fetch_options& options,
                             const std::string&         out_path)
{
    if (!options.use_frame) return fetch_dem_impl(options, out_path);

    // Frame mode: the tile request is the frame's envelope plus five cells
    // (room for the resampling kernel and any datum offset of the base CRS).
    lsdtt_fetch_options o = options;
    try
    {
        fetch_bbox_for_frame(o.frame, 5.0 * o.grid_spacing, o.west, o.south, o.east, o.north);
    }
    catch (const std::exception& e)
    {
        return make_error(std::string("lsdtt::fetch_dem: ") + e.what());
    }

    lsdtt_fetch_result result = fetch_dem_impl(o, out_path);
    if (!result.ok || o.frame_file_path.empty()) return result;

    lsdtt_rotation_target target;
    target.epsg = result.rotate_back_epsg;
    target.cell_size = o.grid_spacing;
    target.origin_x = 0.0;   // the -tap lattice of a plain fetch
    target.origin_y = 0.0;
    lsdtt_rotation_provenance provenance = o.frame_provenance;
    // The file name only: the frame file sits next to the DEM, so the pair
    // can be moved together.
    if (provenance.rotated_dem_fname.empty())
        provenance.rotated_dem_fname = std::filesystem::path(out_path).filename().string();
    try
    {
        write_rotation_frame(o.frame_file_path, o.frame, target, provenance);
    }
    catch (const std::exception& e)
    {
        lsdtt_fetch_result failed = make_error(
            std::string("lsdtt::fetch_dem: the DEM was written but the frame file was not: ") +
            e.what());
        failed.path = out_path;
        return failed;
    }
    result.message += "; frame file " + o.frame_file_path;
    return result;
}

} // namespace lsdtt
