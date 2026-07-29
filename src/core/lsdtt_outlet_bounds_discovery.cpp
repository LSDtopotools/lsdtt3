// lsdtt_outlet_bounds_discovery.cpp

#include "lsdtt_outlet_bounds_discovery.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

namespace lsdtt {
namespace {

constexpr double k_pi = 3.14159265358979323846;

std::string trim(const std::string& s)
{
    std::size_t begin = 0;
    while (begin < s.size() &&
           std::isspace(static_cast<unsigned char>(s[begin])))
    {
        ++begin;
    }
    std::size_t end = s.size();
    while (end > begin &&
           std::isspace(static_cast<unsigned char>(s[end - 1])))
    {
        --end;
    }
    return s.substr(begin, end - begin);
}

std::string to_lower_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string strip_utf8_bom(const std::string& s)
{
    if (s.size() >= 3 &&
        static_cast<unsigned char>(s[0]) == 0xEF &&
        static_cast<unsigned char>(s[1]) == 0xBB &&
        static_cast<unsigned char>(s[2]) == 0xBF)
    {
        return s.substr(3);
    }
    return s;
}

std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ','))
    {
        out.push_back(trim(token));
    }
    return out;
}

std::string lowercase_extension(const std::string& path)
{
    const std::size_t dot = path.find_last_of('.');
    if (dot == std::string::npos) return "";
    return to_lower_copy(path.substr(dot));
}

void ensure_gdal_init()
{
    static bool initialised = false;
    if (!initialised)
    {
        GDALAllRegister();
        initialised = true;
    }
}

void bbox_from_radius_km(double lon, double lat, double radius_km,
                         double& west, double& south,
                         double& east, double& north)
{
    const double radius_deg_lat = radius_km / 111.32;
    const double lat_rad = lat * (k_pi / 180.0);
    double cos_lat = std::cos(lat_rad);
    if (std::fabs(cos_lat) < 1.0e-6) cos_lat = (cos_lat < 0.0 ? -1.0 : 1.0) * 1.0e-6;
    const double radius_deg_lon = radius_km / (111.32 * std::fabs(cos_lat));

    west = lon - radius_deg_lon;
    east = lon + radius_deg_lon;
    south = lat - radius_deg_lat;
    north = lat + radius_deg_lat;

    if (west < -180.0) west = -180.0;
    if (east > 180.0) east = 180.0;
    if (south < -90.0) south = -90.0;
    if (north > 90.0) north = 90.0;
}

} // unnamed namespace

std::vector<lsdtt_outlet_point_wgs84>
load_outlet_points_wgs84(const std::string& input_path,
                         const std::string& layer_name,
                         const std::string& latitude_field_name,
                         const std::string& longitude_field_name)
{
    if (input_path.empty())
    {
        throw std::runtime_error("no outlet points file given");
    }

    ensure_gdal_init();
    std::vector<lsdtt_outlet_point_wgs84> points;
    const std::string ext = lowercase_extension(input_path);

    if (ext == ".csv")
    {
        std::ifstream in(input_path);
        if (!in)
        {
            throw std::runtime_error("could not open outlet points CSV: " + input_path);
        }

        std::string header;
        if (!std::getline(in, header))
        {
            throw std::runtime_error("outlet points CSV is empty: " + input_path);
        }
        header = strip_utf8_bom(header);
        const std::vector<std::string> cols = split_csv_line(header);

        const std::string lat_key = to_lower_copy(trim(latitude_field_name));
        const std::string lon_key = to_lower_copy(trim(longitude_field_name));
        int lat_idx = -1;
        int lon_idx = -1;
        for (std::size_t i = 0; i < cols.size(); ++i)
        {
            const std::string key = to_lower_copy(cols[i]);
            if (key == lat_key) lat_idx = static_cast<int>(i);
            if (key == lon_key) lon_idx = static_cast<int>(i);
        }
        if (lat_idx < 0 || lon_idx < 0)
        {
            throw std::runtime_error(
                "outlet points CSV must contain columns named '" +
                latitude_field_name + "' and '" + longitude_field_name + "'");
        }

        std::string line;
        int next_id = 1;
        while (std::getline(in, line))
        {
            if (trim(line).empty()) continue;
            const std::vector<std::string> tok = split_csv_line(line);
            if (static_cast<int>(tok.size()) <= std::max(lat_idx, lon_idx))
            {
                continue;
            }

            try
            {
                const double lat = std::stod(tok[static_cast<std::size_t>(lat_idx)]);
                const double lon = std::stod(tok[static_cast<std::size_t>(lon_idx)]);
                lsdtt_outlet_point_wgs84 p;
                p.point_id = next_id++;
                p.latitude = lat;
                p.longitude = lon;
                points.push_back(p);
            }
            catch (const std::exception&)
            {
                continue;
            }
        }
        return points;
    }

    GDALDataset* raw_dataset = static_cast<GDALDataset*>(
        GDALOpenEx(input_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr));
    if (!raw_dataset)
    {
        throw std::runtime_error("could not open outlet points vector: " + input_path);
    }

    std::unique_ptr<GDALDataset, void(*)(GDALDataset*)> dataset(raw_dataset, [](GDALDataset* ds)
    {
        if (ds) GDALClose(ds);
    });

    OGRLayer* layer = nullptr;
    if (!layer_name.empty())
    {
        layer = dataset->GetLayerByName(layer_name.c_str());
        if (!layer)
        {
            throw std::runtime_error("outlet points layer not found: " + layer_name);
        }
    }
    else
    {
        layer = dataset->GetLayer(0);
        if (!layer)
        {
            throw std::runtime_error("could not access first outlet-points layer");
        }
    }

    OGRSpatialReference src_srs;
    OGRSpatialReference dst_srs;
    dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (dst_srs.importFromEPSG(4326) != OGRERR_NONE)
    {
        throw std::runtime_error("could not initialise EPSG:4326 for outlet-point transform");
    }

    OGRCoordinateTransformation* tx = nullptr;
    OGRSpatialReference* layer_srs = layer->GetSpatialRef();
    if (layer_srs)
    {
        src_srs = *layer_srs;
        src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        tx = OGRCreateCoordinateTransformation(&src_srs, &dst_srs);
        if (!tx)
        {
            throw std::runtime_error("could not create transform from outlet layer CRS to EPSG:4326");
        }
    }

    int next_id = 1;
    layer->ResetReading();
    OGRFeature* feat = nullptr;
    while ((feat = layer->GetNextFeature()) != nullptr)
    {
        OGRGeometry* geom = feat->GetGeometryRef();
        if (!geom)
        {
            OGRFeature::DestroyFeature(feat);
            continue;
        }

        const OGRwkbGeometryType geom_type = wkbFlatten(geom->getGeometryType());
        if (geom_type != wkbPoint)
        {
            OGRFeature::DestroyFeature(feat);
            continue;
        }

        OGRPoint* p = geom->toPoint();
        double x = p->getX();
        double y = p->getY();
        if (tx)
        {
            if (!tx->Transform(1, &x, &y))
            {
                OGRFeature::DestroyFeature(feat);
                continue;
            }
        }

        lsdtt_outlet_point_wgs84 out;
        out.point_id = next_id++;
        out.longitude = x;
        out.latitude = y;
        points.push_back(out);
        OGRFeature::DestroyFeature(feat);
    }

    if (tx)
    {
        OCTDestroyCoordinateTransformation(tx);
    }

    return points;
}

lsdtt_outlet_bounds_discovery_result
discover_outlet_bounds(const std::vector<lsdtt_outlet_point_wgs84>& points,
                       const lsdtt_outlet_bounds_discovery_options& options)
{
    lsdtt_outlet_bounds_discovery_result result;

    if (points.empty())
    {
        result.ok = false;
        result.message = "no valid outlet points found";
        return result;
    }

    double merged_w = std::numeric_limits<double>::infinity();
    double merged_s = std::numeric_limits<double>::infinity();
    double merged_e = -std::numeric_limits<double>::infinity();
    double merged_n = -std::numeric_limits<double>::infinity();
    int n_ok = 0;

    for (const lsdtt_outlet_point_wgs84& p : points)
    {
        lsdtt_outlet_bounds_record rec;
        rec.point_id = p.point_id;
        rec.longitude = p.longitude;
        rec.latitude = p.latitude;

        if (!(p.longitude >= -180.0 && p.longitude <= 180.0 &&
              p.latitude >= -90.0 && p.latitude <= 90.0))
        {
            rec.status = "invalid_coordinate";
            rec.message = "outlet coordinate is outside WGS84 range";
            result.records.push_back(rec);
            continue;
        }

        double radius = options.initial_radius_km;
        rec.iterations = 1;
        for (int it = 0; it < options.max_iters; ++it)
        {
            bbox_from_radius_km(p.longitude, p.latitude, radius,
                                rec.west, rec.south, rec.east, rec.north);

            // E1 placeholder acceptance: full basin-completeness checks arrive
            // in E2+ when DEM delineation is wired into this loop.
            rec.status = "ok";
            rec.message = "radius-based first guess (E1)";
            rec.iterations = it + 1;
            break;
        }

        if (rec.status == "ok")
        {
            ++n_ok;
            merged_w = std::min(merged_w, rec.west);
            merged_s = std::min(merged_s, rec.south);
            merged_e = std::max(merged_e, rec.east);
            merged_n = std::max(merged_n, rec.north);
        }
        result.records.push_back(rec);
    }

    if (n_ok == 0)
    {
        result.ok = false;
        result.message = "no outlet points produced usable first-guess bounds";
        return result;
    }

    result.ok = true;
    result.merged_west = merged_w;
    result.merged_south = merged_s;
    result.merged_east = merged_e;
    result.merged_north = merged_n;
    result.message = "discovered first-guess bounds for " + std::to_string(n_ok) +
                     " outlet point(s)";
    return result;
}

bool write_outlet_bounds_status_csv(
    const std::string& out_csv,
    const lsdtt_outlet_bounds_discovery_result& result,
    std::string& error)
{
    std::ofstream out(out_csv);
    if (!out)
    {
        error = "could not open status CSV for writing: " + out_csv;
        return false;
    }

    out << "point_id,longitude,latitude,status,iterations,west,south,east,north,message\n";
    for (const lsdtt_outlet_bounds_record& rec : result.records)
    {
        out << rec.point_id << ','
            << rec.longitude << ','
            << rec.latitude << ','
            << rec.status << ','
            << rec.iterations << ','
            << rec.west << ','
            << rec.south << ','
            << rec.east << ','
            << rec.north << ','
            << rec.message << '\n';
    }

    return true;
}

} // namespace lsdtt
