// lsdtt_raster_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// See lsdtt_raster_fetch.hpp for full documentation.

#include "lsdtt_raster_fetch.hpp"
#include "lsdtt_secret.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include <cpl_conv.h>        // CPLSetThreadLocalConfigOption, CPLSPrintf
#include <cpl_error.h>       // CPLPushErrorHandler / CPLQuietErrorHandler
#include <cpl_http.h>        // CPLHTTPEnabled
#include <cpl_string.h>      // CSLAddString / CSLDestroy
#include <cpl_vsi.h>         // VSIStatL, VSIUnlink
#include <gdal.h>
#include <gdal_utils.h>      // GDALTranslate / GDALWarp
#include <cpl_json.h>        // CPLJSONDocument (STAC + SAS responses)
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
// (docs/testing_and_benchmarking.md rule 4): GDAL error messages include the failing URL, so a key in
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


// Cache key: everything that changes the pixels. Omitting any of these would
// let a cached raster be reused for a request it does not answer.
std::string cache_key_for(const lsdtt_fetch_options& o)
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
    }

    // Keep it filesystem-safe.
    std::string out = key.str();
    for (char& c : out)
    {
        if (c == '.' || c == '-' || c == ' ') c = '_';
    }
    return out + ".tif";
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
bool download_to_vsimem(const std::string& url,
                        const std::string& secret,
                        const std::string& vsimem_path,
                        int timeout_seconds,
                        int max_retry,
                        std::string& error)
{
    gdal_string_list_guard http_opts;
    http_opts.add("TIMEOUT="   + std::to_string(timeout_seconds));
    http_opts.add("MAX_RETRY=" + std::to_string(max_retry));
    http_opts.add("RETRY_DELAY=2");

    CPLHTTPResult* response = nullptr;
    {
        quiet_gdal_errors quiet;   // GDAL would otherwise print the keyed URL
        response = CPLHTTPFetch(url.c_str(), http_opts.get());
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
        CPLHTTPDestroyResult(response);
        return false;
    }

    if (response->nDataLen <= 0 || response->pabyData == nullptr)
    {
        error = "the OpenTopography API returned an empty response";
        CPLHTTPDestroyResult(response);
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
        std::string body(reinterpret_cast<const char*>(response->pabyData),
                         static_cast<std::size_t>(response->nDataLen));
        if (body.size() > 400) body = body.substr(0, 400) + "...";
        error = "the OpenTopography API returned an error: " +
                mask_secret(body, secret);
        CPLHTTPDestroyResult(response);
        return false;
    }

    VSILFILE* handle = VSIFileFromMemBuffer(
        vsimem_path.c_str(), response->pabyData,
        static_cast<vsi_l_offset>(response->nDataLen),
        FALSE);   // FALSE: /vsimem does not take ownership of the buffer

    if (handle == nullptr)
    {
        error = "could not stage the downloaded raster in memory";
        CPLHTTPDestroyResult(response);
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
        CPLHTTPDestroyResult(response);
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
    CPLHTTPDestroyResult(response);

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
    gdal_string_list_guard opts;
    opts.add("TIMEOUT="     + std::to_string(timeout_seconds));
    opts.add("MAX_RETRY="   + std::to_string(max_retry));
    opts.add("RETRY_DELAY=2");
    if (!post_body.empty())
    {
        opts.add("POSTFIELDS=" + post_body);
        opts.add("HEADERS=Content-Type: application/json");
    }

    CPLHTTPResult* response = nullptr;
    {
        quiet_gdal_errors quiet;
        response = CPLHTTPFetch(url.c_str(), opts.get());
    }

    if (response == nullptr)
    {
        error = "no response from " + mask_url_query(url);
        return false;
    }
    if (response->nStatus != 0)
    {
        error = "request to " + mask_url_query(url) + " failed: " +
                (response->pszErrBuf ? response->pszErrBuf : "transport error");
        CPLHTTPDestroyResult(response);
        return false;
    }
    if (response->nDataLen <= 0 || response->pabyData == nullptr)
    {
        error = "empty response from " + mask_url_query(url);
        CPLHTTPDestroyResult(response);
        return false;
    }

    body.assign(reinterpret_cast<const char*>(response->pabyData),
                static_cast<std::size_t>(response->nDataLen));
    CPLHTTPDestroyResult(response);
    return true;
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
// pc_search_assets
//
// STAC search for items intersecting the bounding box, returning the signed
// href of the chosen asset for each item.
// -----------------------------------------------------------------------------
bool pc_search_assets(const lsdtt_fetch_options& o,
                      std::vector<std::string>& signed_hrefs,
                      std::vector<std::string>& item_ids,
                      std::string& error)
{
    const std::string collection =
        o.pc_collection.empty() ? "cop-dem-glo-30" : o.pc_collection;
    const std::string asset_key =
        o.pc_asset.empty() ? "data" : o.pc_asset;

    std::ostringstream request;
    request << "{\"collections\":[\"" << collection << "\"],"
            << "\"bbox\":[" << format_double(o.west,  10) << ','
                            << format_double(o.south, 10) << ','
                            << format_double(o.east,  10) << ','
                            << format_double(o.north, 10) << "],"
            << "\"limit\":100}";

    std::string body;
    if (!http_request_body(
            "https://planetarycomputer.microsoft.com/api/stac/v1/search",
            request.str(), o.http_timeout_seconds, o.http_max_retry,
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

    for (int i = 0; i < features.Size(); ++i)
    {
        CPLJSONObject feature = features[i];
        const std::string id = feature.GetString("id", "");

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

    return true;
}

// -----------------------------------------------------------------------------
// warp_and_write
//
// Steps 3 and 4 of the pipeline, shared by every source: reproject to the
// target CRS at a fixed grid spacing, then optionally crop an exact square
// from the centre. Sources differ only in how they obtain the input dataset.
// -----------------------------------------------------------------------------
lsdtt_fetch_result warp_and_write(gdal_dataset_guard& source_ds,
                                  const lsdtt_fetch_options& options,
                                  lsdtt_fetch_result result,
                                  const std::string& out_path)
{
    // ---------------------------------------------------------------------
    // Step 3: warp to a metric CRS at a fixed grid spacing
    // ---------------------------------------------------------------------
    const std::string spacing = format_double(options.grid_spacing, 6);

    gdal_string_list_guard warp_args;
    warp_args.add("-t_srs");     warp_args.add("EPSG:" + std::to_string(result.epsg));
    warp_args.add("-tr");        warp_args.add(spacing); warp_args.add(spacing);
    warp_args.add("-r");         warp_args.add(options.resampling);
    warp_args.add("-tap");
    warp_args.add("-dstnodata"); warp_args.add(format_double(options.no_data_value, 6));
    warp_args.add("-of");        warp_args.add("GTiff");
    warp_args.add("-co");        warp_args.add("COMPRESS=DEFLATE");
    warp_args.add("-co");        warp_args.add("TILED=YES");

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
        GDALDatasetH sources[1] = { source_ds.get() };
        warped.reset(GDALWarp(warp_target.c_str(), nullptr, 1, sources,
                              warp_opts, &usage_error));
        GDALWarpAppOptionsFree(warp_opts);
        if (!warped || usage_error != CE_None)
        {
            return make_error("lsdtt::fetch_dem: warp to EPSG:" +
                              std::to_string(result.epsg) + " failed");
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
                     std::to_string(result.n_rows) + " raster to EPSG:" +
                     std::to_string(result.epsg);
    return result;
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
    return true;
}

} // anonymous namespace


// =============================================================================
// Public helpers
// =============================================================================

bool parse_dem_source(const std::string& name, lsdtt_dem_source& out)
{
    if (name == "cop30_aws")          { out = lsdtt_dem_source::cop30_aws;          return true; }
    if (name == "opentopography")     { out = lsdtt_dem_source::opentopography;     return true; }
    if (name == "planetary_computer") { out = lsdtt_dem_source::planetary_computer; return true; }
    return false;
}


std::string dem_source_name(lsdtt_dem_source source)
{
    switch (source)
    {
        case lsdtt_dem_source::cop30_aws:          return "cop30_aws";
        case lsdtt_dem_source::opentopography:     return "opentopography";
        case lsdtt_dem_source::planetary_computer: return "planetary_computer";
    }
    return "unknown";
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
        {"COP30", "Copernicus GLO-30", "30 m", "global",
         "TanDEM-X derived DSM with edited water bodies and coastlines. The "
         "best general-purpose global 30 m surface; the same data the keyless "
         "cop30_aws source serves.",
         "European Space Agency (2024). Copernicus Global Digital Elevation "
         "Model. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9028PQB",
         450000.0},

        {"COP90", "Copernicus GLO-90", "90 m", "global",
         "The 90 m Copernicus DSM. Use it when the area is too large for the "
         "30 m request limit.",
         "European Space Agency (2024). Copernicus Global Digital Elevation "
         "Model. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9028PQB",
         4050000.0},

        {"SRTMGL1", "SRTM GL1", "30 m", "60 N to 56 S",
         "NASA Shuttle Radar Topography Mission, 1 arc-second, EGM96 heights. "
         "The long-standing comparison baseline in the literature.",
         "NASA Shuttle Radar Topography Mission (SRTM) (2013). Shuttle Radar "
         "Topography Mission (SRTM) Global. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9445JDF",
         450000.0},

        {"SRTMGL1_E", "SRTM GL1 Ellipsoidal", "30 m", "60 N to 56 S",
         "SRTM GL1 with heights referenced to the WGS84 ellipsoid rather than "
         "EGM96. Use it when co-registering with GNSS or lidar point clouds.",
         "NASA Shuttle Radar Topography Mission (SRTM) (2013). Shuttle Radar "
         "Topography Mission (SRTM) Global. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9445JDF",
         450000.0},

        {"SRTMGL3", "SRTM GL3", "90 m", "60 N to 56 S",
         "SRTM at 3 arc-seconds. Larger request limit than GL1.",
         "NASA Shuttle Radar Topography Mission (SRTM) (2013). Shuttle Radar "
         "Topography Mission (SRTM) Global. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9445JDF",
         4050000.0},

        {"NASADEM", "NASADEM", "30 m", "60 N to 56 S",
         "SRTM reprocessed with improved geolocation and void filling from ASTER "
         "GDEM and ICESat. Generally preferable to raw SRTMGL1.",
         "NASA JPL (2021). NASADEM Merged DEM Global 1 arc second V001. "
         "Distributed by OpenTopography. https://doi.org/10.5069/G93T9FD9",
         450000.0},

        {"AW3D30", "ALOS World 3D", "30 m", "global",
         "JAXA optical-stereo DSM (PRISM). An independent global 30 m surface, "
         "useful as a cross-check on the radar products.",
         "Japan Aerospace Exploration Agency (2021). ALOS World 3D 30 meter "
         "DEM. V3.2, Jan 2021. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G94M92HB",
         450000.0},

        {"AW3D30_E", "ALOS World 3D Ellipsoidal", "30 m", "global",
         "AW3D30 with ellipsoidal rather than geoid heights.",
         "Japan Aerospace Exploration Agency (2021). ALOS World 3D 30 meter "
         "DEM. V3.2, Jan 2021. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G94M92HB",
         450000.0},

        {"GEDTM30", "Global Ensemble DTM", "30 m", "global",
         "Machine-learning ensemble bare-earth DTM. A terrain model rather than "
         "a surface model, so vegetation and buildings are removed.",
         "Ho, Y.-F., & Hengl, T. (2025). Global ensemble digital terrain model "
         "30m (GEDTM30) (Version v20250619). Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9BV7DT1",
         450000.0},

        {"EU_DTM", "Continental Europe DTM", "30 m", "Europe",
         "Bare-earth DTM for continental Europe, built from national lidar and "
         "photogrammetric surveys.",
         "Hengl, T., Leal Parente, L., Krizan, J., and Bonannella, C. (2022). "
         "Continental Europe Digital Terrain Model. Distributed by "
         "OpenTopography. https://doi.org/10.5069/G99021ZF",
         450000.0},

        {"ANADEM", "ANADEM", "30 m", "South America",
         "Bare-earth DTM for South America, with vegetation bias removed from "
         "Copernicus GLO-30.",
         "Agencia Nacional de Aguas e Saneamento Basico (2025). ANADEM: A "
         "Digital Terrain Model for South America. Distributed by "
         "OpenTopography. https://doi.org/10.5069/G9736P4G",
         450000.0},

        {"CA_MRDEM_DTM", "Canada MRDEM (terrain)", "30 m", "Canada",
         "Natural Resources Canada medium-resolution bare-earth model.",
         "Natural Resources Canada (2025). Medium Resolution Digital Elevation "
         "Model of Canada. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9N8780J",
         450000.0},

        {"CA_MRDEM_DSM", "Canada MRDEM (surface)", "30 m", "Canada",
         "The surface-model companion to CA_MRDEM_DTM.",
         "Natural Resources Canada (2025). Medium Resolution Digital Elevation "
         "Model of Canada. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G9N8780J",
         450000.0},

        {"SRTM15Plus", "SRTM15+ V2.1", "500 m", "global, land and sea",
         "Merged topography and bathymetry at 15 arc-seconds. Use it when the "
         "domain has to continue offshore.",
         "Tozer, B., Sandwell, D. T., Smith, W. H. F., Olson, C., Beale, J. R., "
         "& Wessel, P. (2019). Global bathymetry and topography at 15 arc sec: "
         "SRTM15+. Distributed by OpenTopography. "
         "https://doi.org/10.5069/G92R3PT9",
         125000000.0},

        {"GEBCOIceTopo", "GEBCO 2023 (ice surface)", "500 m",
         "global, land and sea",
         "GEBCO global bathymetry and topography, with the ice-sheet surface as "
         "the elevation over Greenland and Antarctica.",
         "GEBCO Compilation Group (2023). GEBCO 2023 Grid. Distributed by "
         "OpenTopography. https://doi.org/10.5069/G9D21VTT",
         125000000.0},

        {"GEBCOSubIceTopo", "GEBCO 2023 (bed under ice)", "500 m",
         "global, land and sea",
         "As GEBCOIceTopo, but the ice sheets are replaced by the underlying "
         "bed topography.",
         "GEBCO Compilation Group (2023). GEBCO 2023 Grid. Distributed by "
         "OpenTopography. https://doi.org/10.5069/G9D21VTT",
         125000000.0},

        {"GEDI_L3", "GEDI L3 gridded metrics", "1000 m", "52 N to 52 S",
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

    GDALDatasetH dataset = nullptr;
    {
        quiet_gdal_errors quiet;
        dataset = GDALOpenEx(path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                             nullptr, nullptr, nullptr);
    }
    if (dataset == nullptr)
    {
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


bool http_available()
{
    // CPLHTTPEnabled() reports whether GDAL was built with curl at all. It says
    // nothing about whether the network is reachable -- that is probed
    // separately by trying to open the source.
    return CPLHTTPEnabled() != 0;
}


// =============================================================================
// fetch_dem
// =============================================================================

lsdtt_fetch_result fetch_dem(const lsdtt_fetch_options& options,
                             const std::string&         out_path)
{
    std::string error;
    if (!validate_options(options, error))
    {
        return make_error("lsdtt::fetch_dem: " + error);
    }

    if (options.source != lsdtt_dem_source::cop30_aws &&
        options.source != lsdtt_dem_source::opentopography &&
        options.source != lsdtt_dem_source::planetary_computer)
    {
        // Explicit rather than silently substituting a different source.
        return make_error(
            "lsdtt::fetch_dem: source '" + dem_source_name(options.source) +
            "' is not supported (available: cop30_aws, opentopography, "
            "planetary_computer; see docs/testing_and_benchmarking.md)");
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

    if (!http_available())
    {
        return make_skip(
            "lsdtt::fetch_dem: this GDAL build has no HTTP support (built "
            "without curl), so remote DEMs cannot be fetched");
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
    result.epsg  = (options.target_epsg > 0)
                 ? options.target_epsg
                 : utm_epsg_for_point(0.5 * (options.west + options.east),
                                      0.5 * (options.south + options.north));
    result.grid_spacing = options.grid_spacing;

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
                                download_error))
        {
            // download_error is already masked. Never interpolate `url` here.
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

        const lsdtt_fetch_result warped_result =
            warp_and_write(downloaded, options, result, out_path);

        downloaded.reset();
        VSIUnlink(staged_path.c_str());

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

        if (!pc_search_assets(options, signed_hrefs, result.tiles, search_error))
        {
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

} // namespace lsdtt
