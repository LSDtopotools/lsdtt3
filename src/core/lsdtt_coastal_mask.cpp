// lsdtt_coastal_mask.cpp

#include "lsdtt_coastal_mask.hpp"
#include "lsdtt_gdal_init.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <cpl_conv.h>
#include <cpl_http.h>
#include <cpl_string.h>
#include <cpl_vsi.h>
#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

namespace lsdtt {
namespace {

std::string to_lower_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string choose_layer_name(GDALDataset* dataset,
                              const std::string& requested_layer)
{
    if (!requested_layer.empty()) return requested_layer;
    if (dataset == nullptr) return "";
    OGRLayer* layer = dataset->GetLayer(0);
    if (layer == nullptr) return "";
    const char* name = layer->GetName();
    return name == nullptr ? "" : std::string(name);
}

bool dataset_has_polygonal_layer(GDALDataset* dataset,
                                 const std::string& requested_layer,
                                 std::string& resolved_layer,
                                 std::string& error)
{
    if (dataset == nullptr)
    {
        error = "internal error: null vector dataset";
        return false;
    }

    OGRLayer* layer = nullptr;
    if (!requested_layer.empty())
    {
        layer = dataset->GetLayerByName(requested_layer.c_str());
        if (layer == nullptr)
        {
            error = "coastal mask layer not found: '" + requested_layer + "'";
            return false;
        }
    }
    else
    {
        layer = dataset->GetLayer(0);
        if (layer == nullptr)
        {
            error = "coastal mask dataset contains no layers";
            return false;
        }
    }

    resolved_layer = choose_layer_name(dataset, requested_layer);
    layer->ResetReading();

    OGRFeature* feature = nullptr;
    while ((feature = layer->GetNextFeature()) != nullptr)
    {
        OGRGeometry* geometry = feature->GetGeometryRef();
        if (geometry != nullptr)
        {
            const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
            OGRFeature::DestroyFeature(feature);
            if (type == wkbPolygon || type == wkbMultiPolygon)
            {
                return true;
            }
            error = "coastal mask layer must contain polygon geometry";
            return false;
        }
        OGRFeature::DestroyFeature(feature);
    }

    error = "coastal mask layer contains no features";
    return false;
}

bool dataset_is_readable_polygon_mask(const std::string& path,
                                      const std::string& requested_layer,
                                      std::string& resolved_layer,
                                      std::string& error)
{
    lsdtt::gdal_init();

    GDALDataset* raw = static_cast<GDALDataset*>(
        GDALOpenEx(path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr));
    if (raw == nullptr)
    {
        error = "could not open vector coastal mask dataset: " + path;
        return false;
    }

    std::unique_ptr<GDALDataset, void(*)(GDALDataset*)> dataset(raw, [](GDALDataset* ds)
    {
        if (ds) GDALClose(ds);
    });

    return dataset_has_polygonal_layer(dataset.get(), requested_layer,
                                       resolved_layer, error);
}

std::vector<std::string> default_dataset_candidates(
    const std::string& dataset_id,
    const std::string& read_path)
{
    std::vector<std::string> out;
    if (dataset_id != "natural_earth_110m_land") return out;

    const std::vector<std::string> names = {
        "ne_110m_land.fgb",
        "ne_110m_land.gpkg",
        "ne_110m_land.shp",
        "ne_110m_land.zip"
    };

    const std::filesystem::path root(read_path);
    for (const std::string& n : names)
    {
        out.push_back((root / n).string());
    }
    for (const std::string& n : names)
    {
        out.push_back((root / "coastal_masks" / n).string());
    }
    for (const std::string& n : names)
    {
        out.push_back((root / "masks" / n).string());
    }
    return out;
}

bool download_file_https(const std::string& url,
                         const std::string& out_path,
                         std::string& error)
{
    char** http_opts = nullptr;
    http_opts = CSLSetNameValue(http_opts, "TIMEOUT", "60");
    std::unique_ptr<CPLHTTPResult, decltype(&CPLHTTPDestroyResult)> response(
        CPLHTTPFetch(url.c_str(), http_opts),
        &CPLHTTPDestroyResult);
    CSLDestroy(http_opts);

    if (!response)
    {
        error = "HTTP fetch returned no response";
        return false;
    }
    if (response->nStatus != 0)
    {
        error = response->pszErrBuf ? std::string(response->pszErrBuf)
                                    : "HTTP fetch failed";
        return false;
    }
    if (response->pabyData == nullptr || response->nDataLen <= 0)
    {
        error = "HTTP fetch returned empty payload";
        return false;
    }

    std::error_code ec;
    std::filesystem::create_directories(
        std::filesystem::path(out_path).parent_path(), ec);
    if (ec)
    {
        error = "could not create output directory: " + ec.message();
        return false;
    }

    VSILFILE* f = VSIFOpenL(out_path.c_str(), "wb");
    if (f == nullptr)
    {
        error = "could not open output file for writing";
        return false;
    }

    const size_t written = VSIFWriteL(response->pabyData,
                                      1,
                                      static_cast<size_t>(response->nDataLen),
                                      f);
    const int close_status = VSIFCloseL(f);

    if (written != static_cast<size_t>(response->nDataLen) || close_status != 0)
    {
        error = "could not write full downloaded payload to disk";
        return false;
    }

    return true;
}

constexpr double k_pi = 3.14159265358979323846;

double km_to_deg_lat(double km)
{
    return km / 111.32;
}

double km_to_deg_lon(double km, double latitude_deg)
{
    double cos_lat = std::cos(latitude_deg * (k_pi / 180.0));
    if (std::fabs(cos_lat) < 1.0e-6) cos_lat = (cos_lat < 0.0 ? -1.0 : 1.0) * 1.0e-6;
    return km / (111.32 * std::fabs(cos_lat));
}

std::unique_ptr<OGRPolygon> make_bbox_polygon(const lsdtt_wgs84_bbox& bbox)
{
    auto poly = std::make_unique<OGRPolygon>();
    OGRLinearRing ring;
    ring.addPoint(bbox.west, bbox.south);
    ring.addPoint(bbox.west, bbox.north);
    ring.addPoint(bbox.east, bbox.north);
    ring.addPoint(bbox.east, bbox.south);
    ring.addPoint(bbox.west, bbox.south);
    poly->addRing(&ring);
    return poly;
}

void expand_envelope(const OGRGeometry* geom, OGREnvelope& env, bool& has_env)
{
    if (geom == nullptr || geom->IsEmpty()) return;
    OGREnvelope current;
    geom->getEnvelope(&current);
    if (!has_env)
    {
        env = current;
        has_env = true;
        return;
    }
    env.MinX = std::min(env.MinX, current.MinX);
    env.MaxX = std::max(env.MaxX, current.MaxX);
    env.MinY = std::min(env.MinY, current.MinY);
    env.MaxY = std::max(env.MaxY, current.MaxY);
}

} // unnamed namespace

lsdtt_coastal_mask_resolution resolve_coastal_land_mask(
    const lsdtt_coastal_mask_request& request)
{
    lsdtt_coastal_mask_resolution out;

    const std::string mode = to_lower_copy(request.mode);
    out.requested = (mode == "clip_bbox_edges");
    out.dataset_id = request.dataset;

    if (mode == "off")
    {
        out.available = false;
        out.message = "coastal mask mode is off";
        return out;
    }

    if (mode != "clip_bbox_edges")
    {
        out.available = false;
        out.message = "unsupported coastal mask mode: '" + request.mode + "'";
        return out;
    }

    if (!request.mask_fname.empty())
    {
        std::filesystem::path p(request.mask_fname);
        if (!p.is_absolute())
        {
            const std::filesystem::path base(
                request.mask_path.empty() ? request.read_path : request.mask_path);
            p = base / p;
        }

        out.path = p.string();
        std::string error;
        if (!dataset_is_readable_polygon_mask(out.path, request.mask_layer,
                                              out.layer_name, error))
        {
            out.available = false;
            out.message = error;
            return out;
        }

        out.available = true;
        out.message = "using user-provided coastal mask";
        return out;
    }

    const std::string dataset_id = to_lower_copy(request.dataset);
    if (dataset_id != "natural_earth_110m_land")
    {
        out.available = false;
        out.message = "unsupported coastal_mask_dataset '" + request.dataset +
                      "' (supported: natural_earth_110m_land)";
        return out;
    }

    const std::vector<std::string> candidates =
        default_dataset_candidates(dataset_id, request.read_path);
    for (const std::string& candidate : candidates)
    {
        std::string layer_name;
        std::string error;
        if (dataset_is_readable_polygon_mask(candidate, request.mask_layer,
                                             layer_name, error))
        {
            out.available = true;
            out.path = candidate;
            out.layer_name = layer_name;
            out.message = "using default coastal mask dataset";
            return out;
        }
    }

    if (request.allow_dataset_download)
    {
        const std::filesystem::path local_zip =
            std::filesystem::path(request.read_path) /
            "coastal_masks" /
            "ne_110m_land.zip";

        std::string download_error;
        const std::string url =
            "https://naturalearth.s3.amazonaws.com/110m_physical/ne_110m_land.zip";
        if (download_file_https(url, local_zip.string(), download_error))
        {
            std::string layer_name;
            std::string validate_error;
            if (dataset_is_readable_polygon_mask(local_zip.string(),
                                                 request.mask_layer,
                                                 layer_name,
                                                 validate_error))
            {
                out.available = true;
                out.path = local_zip.string();
                out.layer_name = layer_name;
                out.message = "downloaded and using default coastal mask dataset";
                out.warning =
                    "downloaded Natural Earth ne_110m_land.zip into read_path/coastal_masks/";
                return out;
            }

            out.warning =
                "downloaded Natural Earth land-mask archive, but it was not readable: " +
                validate_error;
        }
        else
        {
            out.warning =
                "automatic download of Natural Earth land-mask dataset failed: " +
                download_error;
        }
    }

    out.available = false;
    if (out.warning.empty())
    {
        out.warning =
            "default coastal mask dataset 'natural_earth_110m_land' not found in "
            "read_path search locations (read_path/, read_path/coastal_masks/, "
            "read_path/masks/); outlet mode will proceed without pre-download "
            "coastal clipping";
    }
    out.message =
        "set coastal_mask_fname to a polygon mask file (e.g. ne_110m_land.*), "
        "or place ne_110m_land in read_path";
    if (request.allow_dataset_download)
    {
        out.message += " (auto-download attempted)";
    }
    return out;
}

lsdtt_coastal_bbox_clip_result clip_bbox_edges_to_land(
    const std::string& mask_path,
    const std::string& mask_layer,
    const lsdtt_wgs84_bbox& input_bbox,
    double outlet_lon,
    double outlet_lat,
    double edge_buffer_km)
{
    lsdtt_coastal_bbox_clip_result out;
    out.bbox = input_bbox;

    if (mask_path.empty())
    {
        out.message = "coastal mask path is empty";
        return out;
    }
    if (!(input_bbox.west < input_bbox.east && input_bbox.south < input_bbox.north))
    {
        out.message = "input bbox is invalid";
        return out;
    }

    lsdtt::gdal_init();
    GDALDataset* raw = static_cast<GDALDataset*>(
        GDALOpenEx(mask_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr));
    if (raw == nullptr)
    {
        out.message = "could not open coastal mask dataset: " + mask_path;
        return out;
    }

    std::unique_ptr<GDALDataset, void(*)(GDALDataset*)> dataset(raw, [](GDALDataset* ds)
    {
        if (ds) GDALClose(ds);
    });

    OGRLayer* layer = nullptr;
    if (!mask_layer.empty())
    {
        layer = dataset->GetLayerByName(mask_layer.c_str());
        if (layer == nullptr)
        {
            out.message = "coastal mask layer not found: '" + mask_layer + "'";
            return out;
        }
    }
    else
    {
        layer = dataset->GetLayer(0);
        if (layer == nullptr)
        {
            out.message = "coastal mask dataset contains no layers";
            return out;
        }
    }

    OGRSpatialReference dst_srs;
    dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (dst_srs.importFromEPSG(4326) != OGRERR_NONE)
    {
        out.message = "could not initialise EPSG:4326 for coastal clipping";
        return out;
    }

    OGRCoordinateTransformation* tx = nullptr;
    const OGRSpatialReference* layer_srs = layer->GetSpatialRef();
    OGRSpatialReference src_srs;
    if (layer_srs != nullptr)
    {
        src_srs = *layer_srs;
        src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        tx = OGRCreateCoordinateTransformation(&src_srs, &dst_srs);
        if (tx == nullptr)
        {
            out.message = "could not create coastal-mask transform to EPSG:4326";
            return out;
        }
    }
    else
    {
        out.warning = "coastal mask layer has no CRS; assuming EPSG:4326";
    }

    const std::unique_ptr<OGRPolygon> bbox_poly = make_bbox_polygon(input_bbox);

    OGREnvelope land_env;
    bool has_land_env = false;

    layer->ResetReading();
    OGRFeature* feature = nullptr;
    while ((feature = layer->GetNextFeature()) != nullptr)
    {
        OGRGeometry* geom = feature->GetGeometryRef();
        if (geom == nullptr)
        {
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        std::unique_ptr<OGRGeometry> geom_copy(geom->clone());
        if (!geom_copy)
        {
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        if (tx != nullptr)
        {
            if (geom_copy->transform(tx) != OGRERR_NONE)
            {
                OGRFeature::DestroyFeature(feature);
                continue;
            }
        }

        std::unique_ptr<OGRGeometry> intersection(geom_copy->Intersection(bbox_poly.get()));
        expand_envelope(intersection.get(), land_env, has_land_env);

        OGRFeature::DestroyFeature(feature);
    }

    if (tx != nullptr) OCTDestroyCoordinateTransformation(tx);

    if (!has_land_env)
    {
        out.ok = true;
        out.clipped = false;
        out.message = "no land intersects bbox; keeping radius bounds";
        return out;
    }

    const double buffer_lat = km_to_deg_lat(edge_buffer_km);
    const double buffer_lon = km_to_deg_lon(edge_buffer_km, outlet_lat);

    lsdtt_wgs84_bbox clipped = input_bbox;
    clipped.west  = std::max(input_bbox.west,  land_env.MinX - buffer_lon);
    clipped.east  = std::min(input_bbox.east,  land_env.MaxX + buffer_lon);
    clipped.south = std::max(input_bbox.south, land_env.MinY - buffer_lat);
    clipped.north = std::min(input_bbox.north, land_env.MaxY + buffer_lat);

    const double eps = 1.0e-6;
    clipped.west  = std::min(clipped.west, outlet_lon - eps);
    clipped.east  = std::max(clipped.east, outlet_lon + eps);
    clipped.south = std::min(clipped.south, outlet_lat - eps);
    clipped.north = std::max(clipped.north, outlet_lat + eps);

    clipped.west = std::max(-180.0, clipped.west);
    clipped.east = std::min(180.0, clipped.east);
    clipped.south = std::max(-90.0, clipped.south);
    clipped.north = std::min(90.0, clipped.north);

    if (!(clipped.west < clipped.east && clipped.south < clipped.north))
    {
        out.ok = true;
        out.clipped = false;
        out.message = "land-mask clipping produced invalid bbox; keeping radius bounds";
        return out;
    }

    out.ok = true;
    out.bbox = clipped;
    out.clipped =
        std::fabs(clipped.west - input_bbox.west) > 1.0e-12 ||
        std::fabs(clipped.south - input_bbox.south) > 1.0e-12 ||
        std::fabs(clipped.east - input_bbox.east) > 1.0e-12 ||
        std::fabs(clipped.north - input_bbox.north) > 1.0e-12;
    out.message = out.clipped ? "bbox clipped to land-intersecting envelope"
                              : "bbox unchanged by land mask";
    return out;
}

} // namespace lsdtt
