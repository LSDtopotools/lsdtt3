// lsdtt_river_line_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_river_line_fetch.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_river_line_fetch.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

#include <cpl_conv.h>
#include <cpl_json.h>
#include <cpl_vsi.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_gdal_init.hpp"
#include "lsdtt_macrostrat_fetch.hpp"   // utm_epsg_for
#include "lsdtt_param_parser.hpp"      // lsdtt3_package_version
#include "lsdtt_raster_fetch.hpp"      // bounds_from_vector_file

namespace lsdtt
{

namespace
{

std::string lower(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Attribute values are compared with spaces written as underscores, because a
// parameter value cannot contain a space ("in direction" -> "in_direction").
std::string underscore_spaces(std::string s)
{
    for (char& c : s) if (c == ' ') c = '_';
    return s;
}

std::string trim(const std::string& s)
{
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

struct dataset_closer
{
    void operator()(GDALDataset* d) const { if (d) GDALClose(d); }
};
using dataset_ptr = std::unique_ptr<GDALDataset, dataset_closer>;

struct ct_deleter
{
    void operator()(OGRCoordinateTransformation* c) const
    {
        if (c) OGRCoordinateTransformation::DestroyCT(c);
    }
};
using ct_ptr = std::unique_ptr<OGRCoordinateTransformation, ct_deleter>;

ct_ptr make_ct(const OGRSpatialReference& from, const OGRSpatialReference& to)
{
    ct_ptr ct(OGRCreateCoordinateTransformation(&from, &to));
    if (!ct) throw std::runtime_error("could not build a coordinate transformation");
    return ct;
}

// Envelope, in "to", of the box [x0,x1]x[y0,y1] in "from", densified so the
// projected edges follow the curved ones.
void transform_box(const OGRSpatialReference& from, const OGRSpatialReference& to,
                   double x0, double y0, double x1, double y1,
                   double& ox0, double& oy0, double& ox1, double& oy1)
{
    ct_ptr ct = make_ct(from, to);
    const int k = 50;
    std::vector<double> xs, ys;
    for (int i = 0; i <= k; ++i)
    {
        const double f = static_cast<double>(i) / k;
        xs.push_back(x0 + (x1 - x0) * f); ys.push_back(y0);
        xs.push_back(x0 + (x1 - x0) * f); ys.push_back(y1);
        xs.push_back(x0); ys.push_back(y0 + (y1 - y0) * f);
        xs.push_back(x1); ys.push_back(y0 + (y1 - y0) * f);
    }
    std::vector<int> ok(xs.size(), 0);
    ct->Transform(static_cast<int>(xs.size()), xs.data(), ys.data(), nullptr, ok.data());
    ox0 = oy0 = std::numeric_limits<double>::infinity();
    ox1 = oy1 = -std::numeric_limits<double>::infinity();
    int n_ok = 0;
    for (size_t i = 0; i < xs.size(); ++i)
    {
        if (!ok[i] || !std::isfinite(xs[i]) || !std::isfinite(ys[i])) continue;
        ++n_ok;
        ox0 = std::min(ox0, xs[i]); ox1 = std::max(ox1, xs[i]);
        oy0 = std::min(oy0, ys[i]); oy1 = std::max(oy1, ys[i]);
    }
    if (n_ok == 0) throw std::runtime_error("the area could not be transformed between CRSs");
}

std::string srs_description(const OGRSpatialReference& srs, int& epsg)
{
    OGRSpatialReference s(srs);
    epsg = 0;
    const char* auth = s.GetAuthorityName(nullptr);
    const char* code = s.GetAuthorityCode(nullptr);
    if ((!auth || !code) && s.AutoIdentifyEPSG() == OGRERR_NONE)
    {
        auth = s.GetAuthorityName(nullptr);
        code = s.GetAuthorityCode(nullptr);
    }
    if (auth && code && std::string(auth) == "EPSG") epsg = std::atoi(code);
    const char* name = s.GetName();
    std::string d = name ? name : "unnamed CRS";
    if (epsg > 0) d += " (EPSG:" + std::to_string(epsg) + ")";
    return d;
}

// --- MD5 (RFC 1321) ----------------------------------------------------------
std::uint32_t rotl(std::uint32_t x, int c) { return (x << c) | (x >> (32 - c)); }

} // namespace

std::string md5_hex(const std::string& bytes)
{
    static const std::array<int, 64> s = {
        7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
        5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20,
        4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
        6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
    std::array<std::uint32_t, 64> K{};
    for (int i = 0; i < 64; ++i)
        K[static_cast<size_t>(i)] = static_cast<std::uint32_t>(
            std::floor(std::fabs(std::sin(static_cast<double>(i + 1))) * 4294967296.0));

    std::uint32_t a0 = 0x67452301, b0 = 0xefcdab89, c0 = 0x98badcfe, d0 = 0x10325476;
    std::string msg = bytes;
    const std::uint64_t bit_len = static_cast<std::uint64_t>(bytes.size()) * 8ULL;
    msg.push_back(static_cast<char>(0x80));
    while (msg.size() % 64 != 56) msg.push_back('\0');
    for (int i = 0; i < 8; ++i) msg.push_back(static_cast<char>((bit_len >> (8 * i)) & 0xff));

    for (size_t off = 0; off < msg.size(); off += 64)
    {
        std::uint32_t M[16];
        for (int i = 0; i < 16; ++i)
        {
            const auto* p = reinterpret_cast<const unsigned char*>(msg.data() + off + 4 * i);
            M[i] = static_cast<std::uint32_t>(p[0]) | (static_cast<std::uint32_t>(p[1]) << 8) |
                   (static_cast<std::uint32_t>(p[2]) << 16) | (static_cast<std::uint32_t>(p[3]) << 24);
        }
        std::uint32_t A = a0, B = b0, C = c0, D = d0;
        for (int i = 0; i < 64; ++i)
        {
            std::uint32_t F;
            int g;
            if (i < 16)      { F = (B & C) | (~B & D); g = i; }
            else if (i < 32) { F = (D & B) | (~D & C); g = (5 * i + 1) % 16; }
            else if (i < 48) { F = B ^ C ^ D;          g = (3 * i + 5) % 16; }
            else             { F = C ^ (B | ~D);       g = (7 * i) % 16; }
            F = F + A + K[static_cast<size_t>(i)] + M[g];
            A = D; D = C; C = B;
            B = B + rotl(F, s[static_cast<size_t>(i)]);
        }
        a0 += A; b0 += B; c0 += C; d0 += D;
    }
    std::ostringstream os;
    for (std::uint32_t v : {a0, b0, c0, d0})
        for (int i = 0; i < 4; ++i)
            os << std::hex << std::setw(2) << std::setfill('0') << ((v >> (8 * i)) & 0xff);
    return os.str();
}

std::string river_feature_class_from_value(const std::string& value)
{
    const std::string v = lower(trim(value));
    if (v == "river" || v == "inlandriver") return "river";
    if (v == "stream" || v == "brook") return "stream";
    if (v == "lake" || v == "lake_connector" || v == "artificialpath") return "lake_connector";
    if (v == "canal") return "canal";
    if (v == "drain" || v == "ditch") return "drain";
    if (v == "tidalriver" || v == "tidal") return "tidal";
    return "other";
}

std::map<std::string, int> parse_direction_field_mapping(const std::string& text)
{
    std::map<std::string, int> m;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        if (trim(item).empty()) continue;
        const auto eq = item.rfind('=');
        if (eq == std::string::npos)
            throw std::invalid_argument("direction_field_mapping: '" + trim(item) +
                                        "' is not value=1 or value=-1");
        const std::string key = underscore_spaces(trim(item.substr(0, eq)));
        const std::string val = trim(item.substr(eq + 1));
        if (val == "1" || val == "+1") m[key] = 1;
        else if (val == "-1") m[key] = -1;
        else
            throw std::invalid_argument("direction_field_mapping: the sign for '" + key +
                                        "' must be 1 or -1, not '" + val + "'");
    }
    return m;
}

lsdtt_os_download_entry parse_os_downloads_geopackage(const std::string& json)
{
    CPLJSONDocument doc;
    if (!doc.LoadMemory(json)) throw std::runtime_error("OS Downloads API: response is not JSON");
    const CPLJSONArray arr = doc.GetRoot().ToArray();
    if (!arr.IsValid()) throw std::runtime_error("OS Downloads API: expected a JSON array");
    for (int i = 0; i < arr.Size(); ++i)
    {
        const CPLJSONObject o = arr[i];
        if (o.GetString("format") != "GeoPackage") continue;
        lsdtt_os_download_entry e;
        e.url = o.GetString("url");
        e.file_name = o.GetString("fileName");
        e.md5 = lower(o.GetString("md5"));
        e.size = o.GetLong("size");
        if (e.url.empty() || e.file_name.empty())
            throw std::runtime_error("OS Downloads API: GeoPackage entry without url or fileName");
        return e;
    }
    throw std::runtime_error("OS Downloads API: no GeoPackage download listed");
}

std::string default_river_cache_dir()
{
    namespace fs = std::filesystem;
    if (const char* x = std::getenv("XDG_CACHE_HOME"); x && *x)
        return (fs::path(x) / "lsdtt3" / "river_lines").string();
    if (const char* h = std::getenv("HOME"); h && *h)
        return (fs::path(h) / ".cache" / "lsdtt3" / "river_lines").string();
    if (const char* l = std::getenv("LOCALAPPDATA"); l && *l)
        return (fs::path(l) / "lsdtt3" / "river_lines").string();
    return "lsdtt3_river_lines_cache";
}

lsdtt_river_line_field_map os_open_rivers_field_map()
{
    lsdtt_river_line_field_map f;
    f.id_field = "id";
    f.name_field = "watercourse_name";
    f.feature_class_field = "form";
    f.direction_field = "flow_direction";
    f.direction_values = {{"in_direction", 1}, {"in_opposite_direction", -1}};
    f.start_node_field = "start_node";
    f.end_node_field = "end_node";
    return f;
}

namespace
{

// Download (once) and return the local path of the Open Rivers zip.
std::string ensure_os_open_rivers_cached(const lsdtt_river_line_fetch_options& o,
                                         lsdtt_river_line_read_result& r)
{
    namespace fs = std::filesystem;
    const fs::path dir = fs::path(o.cache_dir.empty() ? default_river_cache_dir() : o.cache_dir) /
                         "gb_os_open_rivers";
    // Any oprvrs_gpkg*.zip already there is used as it is.
    if (fs::is_directory(dir))
    {
        for (const auto& e : fs::directory_iterator(dir))
        {
            const std::string n = e.path().filename().string();
            if (e.is_regular_file() && n.rfind("oprvrs_gpkg", 0) == 0 && e.path().extension() == ".zip")
            {
                r.notes.push_back("OS Open Rivers from the cache: " + e.path().string());
                return e.path().string();
            }
        }
    }
    if (!http_fetch_available())
        throw std::runtime_error("OS Open Rivers is not in the cache (" + dir.string() +
                                 ") and this GDAL cannot make HTTP requests");
    std::string body, error;
    if (!http_request_body(o.os_downloads_api_url, "", o.http, body, error))
        throw std::runtime_error("OS Downloads API: " + error);
    const lsdtt_os_download_entry e = parse_os_downloads_geopackage(body);
    r.notes.push_back("Downloading OS Open Rivers (" + e.file_name + ", " +
                      std::to_string(e.size / 1000000) + " MB) into " + dir.string() +
                      "; later runs read the cache");
    std::string data;
    lsdtt_http_request_options big = o.http;
    big.timeout_seconds = std::max(big.timeout_seconds, 600);
    if (!http_request_body(e.url, "", big, data, error))
        throw std::runtime_error("OS Open Rivers download: " + error);
    if (e.size > 0 && static_cast<long long>(data.size()) != e.size)
        throw std::runtime_error("OS Open Rivers download: got " + std::to_string(data.size()) +
                                 " bytes, the API lists " + std::to_string(e.size));
    if (!e.md5.empty())
    {
        const std::string got = md5_hex(data);
        if (got != e.md5)
            throw std::runtime_error("OS Open Rivers download: md5 " + got + " does not match the "
                                     "API's " + e.md5);
    }
    fs::create_directories(dir);
    const fs::path out = dir / e.file_name;
    const fs::path tmp = dir / (e.file_name + ".part");
    {
        std::ofstream f(tmp, std::ios::binary);
        f.write(data.data(), static_cast<std::streamsize>(data.size()));
        if (!f) throw std::runtime_error("could not write " + tmp.string());
    }
    fs::rename(tmp, out);
    r.notes.push_back("md5 verified: " + e.md5);
    return out.string();
}

// The first .gpkg inside a zip (VSI path), e.g. /vsizip/x.zip/Data/oprvrs_gb.gpkg.
std::string first_gpkg_in_zip(const std::string& zip)
{
    const std::string root = "/vsizip/" + zip;
    char** list = VSIReadDirRecursive(root.c_str());
    std::string found;
    for (char** p = list; p && *p; ++p)
    {
        const std::string n = *p;
        if (n.size() > 5 && lower(n.substr(n.size() - 5)) == ".gpkg") { found = root + "/" + n; break; }
    }
    CSLDestroy(list);
    if (found.empty()) throw std::runtime_error("no .gpkg inside " + zip);
    return found;
}

std::string read_vsi_text(const std::string& path)
{
    VSILFILE* f = VSIFOpenL(path.c_str(), "rb");
    if (!f) return "";
    std::string s;
    char buf[4096];
    size_t n;
    while ((n = VSIFReadL(buf, 1, sizeof(buf), f)) > 0) s.append(buf, n);
    VSIFCloseL(f);
    return s;
}

std::string os_attribution(const std::string& zip)
{
    std::string year;
    const std::string lic = read_vsi_text("/vsizip/" + zip + "/Doc/licence.txt");
    std::smatch m;
    if (std::regex_search(lic, m, std::regex("database rights ([0-9]{4})"))) year = m[1];
    std::ostringstream os;
    os << "OS Open Rivers, Ordnance Survey.\n"
       << "Contains OS data (c) Crown copyright and database rights "
       << (year.empty() ? "[year: see Doc/licence.txt in the download]" : year) << ".\n"
       << "Licensed under the Open Government Licence v3.0 (OS OpenData licence,\n"
       << "http://os.uk/opendata/licence).\n"
       << "Source file: " << zip << "\n";
    return os.str();
}

int field_index(OGRFeatureDefn* defn, const std::string& name, const std::string& what)
{
    if (name.empty()) return -1;
    const int i = defn->GetFieldIndex(name.c_str());
    if (i < 0)
    {
        std::string avail;
        for (int k = 0; k < defn->GetFieldCount(); ++k)
            avail += (k ? ", " : "") + std::string(defn->GetFieldDefn(k)->GetNameRef());
        throw std::invalid_argument(what + " '" + name + "' is not a field of the layer (fields: " +
                                    avail + ")");
    }
    return i;
}

std::string field_text(OGRFeature* f, int i)
{
    if (i < 0 || !f->IsFieldSetAndNotNull(i)) return "";
    return f->GetFieldAsString(i);
}

void append_parts(const OGRGeometry* g, std::vector<std::vector<lsdtt_xy>>& parts)
{
    const OGRwkbGeometryType t = wkbFlatten(g->getGeometryType());
    if (t == wkbLineString)
    {
        const auto* ls = g->toLineString();
        std::vector<lsdtt_xy> v;
        v.reserve(static_cast<size_t>(ls->getNumPoints()));
        for (int i = 0; i < ls->getNumPoints(); ++i) v.push_back({ls->getX(i), ls->getY(i)});
        parts.push_back(std::move(v));
    }
    else if (t == wkbMultiLineString || t == wkbGeometryCollection)
    {
        const auto* gc = g->toGeometryCollection();
        for (int i = 0; i < gc->getNumGeometries(); ++i) append_parts(gc->getGeometryRef(i), parts);
    }
}

} // namespace

lsdtt_river_line_read_result read_river_lines(const lsdtt_river_line_fetch_options& o)
{
    gdal_init();
    lsdtt_river_line_read_result r;

    if (o.source != "local_file" && o.source != "gb_os_open_rivers")
        throw std::invalid_argument("river_line_source '" + o.source +
                                    "' is not one of: local_file, gb_os_open_rivers");

    // ---- area and output CRS -------------------------------------------------
    OGRSpatialReference wgs84;
    wgs84.importFromEPSG(4326);
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRSpatialReference out_srs;
    if (!o.bounds_raster_path.empty())
    {
        dataset_ptr rs(static_cast<GDALDataset*>(GDALOpen(o.bounds_raster_path.c_str(), GA_ReadOnly)));
        if (!rs) throw std::runtime_error("could not open raster " + o.bounds_raster_path);
        double gt[6];
        if (rs->GetGeoTransform(gt) != CE_None || gt[2] != 0.0 || gt[4] != 0.0)
            throw std::runtime_error("raster " + o.bounds_raster_path +
                                     " has no north-up geotransform");
        const OGRSpatialReference* rsrs = rs->GetSpatialRef();
        if (!rsrs) throw std::runtime_error("raster " + o.bounds_raster_path + " has no CRS");
        OGRSpatialReference raster_srs(*rsrs);
        raster_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        const double x0 = gt[0], x1 = gt[0] + gt[1] * rs->GetRasterXSize();
        const double y1 = gt[3], y0 = gt[3] + gt[5] * rs->GetRasterYSize();
        if (o.output_epsg > 0)
        {
            if (out_srs.importFromEPSG(o.output_epsg) != OGRERR_NONE)
                throw std::invalid_argument("unknown EPSG:" + std::to_string(o.output_epsg));
            out_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
            transform_box(raster_srs, out_srs, std::min(x0, x1), std::min(y0, y1),
                          std::max(x0, x1), std::max(y0, y1),
                          r.area_x_min, r.area_y_min, r.area_x_max, r.area_y_max);
        }
        else
        {
            out_srs = raster_srs;
            r.area_x_min = std::min(x0, x1); r.area_x_max = std::max(x0, x1);
            r.area_y_min = std::min(y0, y1); r.area_y_max = std::max(y0, y1);
        }
    }
    else
    {
        if (!(o.west < o.east) || !(o.south < o.north) || o.west < -180 || o.east > 180 ||
            o.south < -90 || o.north > 90)
            throw std::invalid_argument("the bounding box must have west < east and south < north, "
                                        "in WGS84 degrees");
        const int epsg = o.output_epsg > 0
                             ? o.output_epsg
                             : utm_epsg_for(0.5 * (o.west + o.east), 0.5 * (o.south + o.north));
        if (out_srs.importFromEPSG(epsg) != OGRERR_NONE)
            throw std::invalid_argument("unknown EPSG:" + std::to_string(epsg));
        out_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        transform_box(wgs84, out_srs, o.west, o.south, o.east, o.north,
                      r.area_x_min, r.area_y_min, r.area_x_max, r.area_y_max);
    }
    if (out_srs.IsGeographic())
        throw std::invalid_argument("the output CRS is geographic; river lines need a projected CRS "
                                    "in metres (set output_epsg)");
    if (std::fabs(out_srs.GetLinearUnits() - 1.0) > 1e-6)
        throw std::invalid_argument("the output CRS is not in metres");
    r.crs_description = srs_description(out_srs, r.epsg);
    {
        char* wkt = nullptr;
        out_srs.exportToWkt(&wkt);
        r.crs_wkt = wkt ? wkt : "";
        CPLFree(wkt);
    }

    const double shorter = std::min(r.area_x_max - r.area_x_min, r.area_y_max - r.area_y_min);
    r.clip_buffer_metres = o.clip_buffer_metres >= 0.0 ? o.clip_buffer_metres
                                                       : std::max(200.0, 0.05 * shorter);
    r.clip_x_min = r.area_x_min - r.clip_buffer_metres;
    r.clip_y_min = r.area_y_min - r.clip_buffer_metres;
    r.clip_x_max = r.area_x_max + r.clip_buffer_metres;
    r.clip_y_max = r.area_y_max + r.clip_buffer_metres;
    r.node_snap_tolerance_metres = o.node_snap_tolerance_metres >= 0.0 ? o.node_snap_tolerance_metres : 1.0;

    // ---- the source ------------------------------------------------------------
    lsdtt_river_line_field_map fields = o.fields;
    std::string layer_name = o.river_lines_layer;
    if (o.source == "local_file")
    {
        if (o.river_lines_path.empty())
            throw std::invalid_argument("river_line_source local_file needs river_lines_fname");
        r.read_path = o.river_lines_path;
        r.attribution = "River lines from a local file: " + o.river_lines_path +
                        "\nThe licence and citation are those of that file's source.\n";
    }
    else
    {
        const std::string zip = ensure_os_open_rivers_cached(o, r);
        r.read_path = first_gpkg_in_zip(zip);
        layer_name = "watercourse_link";
        fields = os_open_rivers_field_map();
        r.attribution = os_attribution(zip);
    }

    // GDAL 3.4 warns that GeoPackage 1.4 "may only be partially supported";
    // reading works, so that warning is not shown.
    CPLPushErrorHandler(CPLQuietErrorHandler);
    dataset_ptr src(static_cast<GDALDataset*>(
        GDALOpenEx(r.read_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, nullptr, nullptr, nullptr)));
    CPLPopErrorHandler();
    if (!src) throw std::runtime_error("could not open vector " + r.read_path);
    OGRLayer* layer = layer_name.empty() ? src->GetLayer(0) : src->GetLayerByName(layer_name.c_str());
    if (!layer)
        throw std::runtime_error("no layer " + (layer_name.empty() ? std::string("0") : layer_name) +
                                 " in " + r.read_path);

    ct_ptr to_out;
    const OGRSpatialReference* lsrs = layer->GetSpatialRef();
    if (lsrs)
    {
        OGRSpatialReference layer_srs(*lsrs);
        layer_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        double fx0, fy0, fx1, fy1;
        transform_box(out_srs, layer_srs, r.clip_x_min, r.clip_y_min, r.clip_x_max, r.clip_y_max,
                      fx0, fy0, fx1, fy1);
        layer->SetSpatialFilterRect(fx0, fy0, fx1, fy1);
        if (!layer_srs.IsSame(&out_srs)) to_out = make_ct(layer_srs, out_srs);
    }
    else
    {
        r.warnings.push_back("the river line layer has no CRS; it is assumed to be in the output CRS");
        layer->SetSpatialFilterRect(r.clip_x_min, r.clip_y_min, r.clip_x_max, r.clip_y_max);
    }

    OGRFeatureDefn* defn = layer->GetLayerDefn();
    const int i_id = field_index(defn, fields.id_field, "id field");
    int i_name = field_index(defn, fields.name_field, "name field");
    if (fields.name_field.empty())
    {
        for (const char* n : {"watercourse_name", "name", "gnis_name"})
        {
            i_name = defn->GetFieldIndex(n);
            if (i_name >= 0) break;
        }
    }
    const int i_class = field_index(defn, fields.feature_class_field, "feature_class_field_name");
    const int i_dir = field_index(defn, fields.direction_field, "direction_field_name");
    const int i_start = field_index(defn, fields.start_node_field, "start_node_field_name");
    const int i_end = field_index(defn, fields.end_node_field, "end_node_field_name");
    const int i_culvert = field_index(defn, fields.culvert_field, "culvert field");
    if (i_dir >= 0 && fields.direction_values.empty())
        throw std::invalid_argument("direction_field_name is set but direction_field_mapping is empty");

    const std::set<std::string> keep(o.keep_feature_classes.begin(), o.keep_feature_classes.end());
    std::vector<river_line_feature> features;
    layer->ResetReading();
    for (auto& feat : *layer)
    {
        ++r.n_features_read;
        const OGRGeometry* g = feat->GetGeometryRef();
        if (!g) { ++r.n_non_line_features; continue; }
        const OGRwkbGeometryType t = wkbFlatten(g->getGeometryType());
        if (t != wkbLineString && t != wkbMultiLineString) { ++r.n_non_line_features; continue; }

        river_line_feature rf;
        rf.feature_class = i_class >= 0 ? river_feature_class_from_value(field_text(feat.get(), i_class))
                                        : std::string("river");
        auto& cc = r.classes[rf.feature_class];
        ++cc.n_read;
        if (!keep.count(rf.feature_class)) continue;
        ++cc.n_kept;

        std::unique_ptr<OGRGeometry> gc(g->clone());
        if (to_out && gc->transform(to_out.get()) != OGRERR_NONE)
        {
            r.warnings.push_back("a feature could not be reprojected and was skipped");
            continue;
        }
        append_parts(gc.get(), rf.parts);
        rf.source_feature_id = i_id >= 0 ? field_text(feat.get(), i_id) : std::to_string(feat->GetFID());
        rf.name = field_text(feat.get(), i_name);
        rf.start_node_id = field_text(feat.get(), i_start);
        rf.end_node_id = field_text(feat.get(), i_end);
        if (i_dir >= 0)
        {
            const std::string v = underscore_spaces(trim(field_text(feat.get(), i_dir)));
            const auto it = fields.direction_values.find(v);
            if (it != fields.direction_values.end())
            {
                rf.direction = it->second;
                rf.direction_source = "attribute";
            }
        }
        else if (fields.digitised_downstream)
        {
            rf.direction = 1;
            rf.direction_source = "digitised";
        }
        if (i_culvert >= 0)
        {
            const std::string v = lower(trim(field_text(feat.get(), i_culvert)));
            rf.culverted = (v == "culvert" || v == "yes" || v == "1" || v == "true") ? 1 : 0;
        }
        r.direction_sources[rf.direction_source] += static_cast<int>(rf.parts.size());
        features.push_back(std::move(rf));
    }
    if (features.empty())
        throw std::runtime_error("no river lines of the kept classes in the area (read " +
                                 std::to_string(r.n_features_read) + " features from " + r.read_path + ")");

    r.features = std::move(features);

    // ---- optional roads ---------------------------------------------------------
    if (!o.road_vector_path.empty())
    {
        CPLPushErrorHandler(CPLQuietErrorHandler);
        dataset_ptr rds(static_cast<GDALDataset*>(GDALOpenEx(
            o.road_vector_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, nullptr, nullptr, nullptr)));
        CPLPopErrorHandler();
        if (!rds) throw std::runtime_error("could not open road vector " + o.road_vector_path);
        OGRLayer* rl = o.road_layer.empty() ? rds->GetLayer(0) : rds->GetLayerByName(o.road_layer.c_str());
        if (!rl) throw std::runtime_error("no road layer " + o.road_layer + " in " + o.road_vector_path);
        ct_ptr road_to_out;
        if (const OGRSpatialReference* rsrs = rl->GetSpatialRef())
        {
            OGRSpatialReference road_srs(*rsrs);
            road_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
            double fx0, fy0, fx1, fy1;
            transform_box(out_srs, road_srs, r.clip_x_min, r.clip_y_min, r.clip_x_max, r.clip_y_max,
                          fx0, fy0, fx1, fy1);
            rl->SetSpatialFilterRect(fx0, fy0, fx1, fy1);
            if (!road_srs.IsSame(&out_srs)) road_to_out = make_ct(road_srs, out_srs);
        }
        else
        {
            r.warnings.push_back("the road layer has no CRS; it is assumed to be in the output CRS");
            rl->SetSpatialFilterRect(r.clip_x_min, r.clip_y_min, r.clip_x_max, r.clip_y_max);
        }
        const int i_road_name = rl->GetLayerDefn()->GetFieldIndex(o.road_name_field.c_str());
        r.roads.schema.field_names = {o.road_name_field};
        r.roads.schema.field_types = {lsdtt_vector_field_type::string_type};
        int n_roads = 0;
        rl->ResetReading();
        for (auto& feat : *rl)
        {
            const OGRGeometry* g = feat->GetGeometryRef();
            if (!g) continue;
            std::unique_ptr<OGRGeometry> gc(g->clone());
            if (road_to_out && gc->transform(road_to_out.get()) != OGRERR_NONE) continue;
            std::vector<std::vector<lsdtt_xy>> parts;
            append_parts(gc.get(), parts);
            const std::string name = field_text(feat.get(), i_road_name);
            for (auto& p : parts)
            {
                if (p.size() < 2) continue;
                lsdtt_vector_record rec;
                rec.geometry = lsdtt_line_geom{std::move(p)};
                rec.field_values = {name};
                r.roads.records.push_back(std::move(rec));
                ++n_roads;
            }
        }
        r.notes.push_back("Roads: " + std::to_string(n_roads) + " line(s) in the clip box from " +
                          o.road_vector_path);
    }
    return r;
}

// ---------------------------------------------------------------------------
// build_river_line_fetch_options (moved from lsdtt-fetch-river-lines)
// ---------------------------------------------------------------------------
namespace
{
bool parse_degrees_text(const std::string& s, double& out)
{
    if (s.empty()) return false;
    std::size_t pos = 0;
    try { out = std::stod(s, &pos); }
    catch (const std::exception&) { return false; }
    return pos == s.size() && std::isfinite(out);
}
} // namespace

lsdtt_river_line_fetch_options build_river_line_fetch_options(
    const lsdtt_river_line_parameters& p,
    const std::function<std::string(const std::string&)>& resolve_path,
    std::vector<std::string>& warnings,
    std::string* area_note)
{
    auto S = [&](const char* k, const char* d = "") -> std::string
    { auto it = p.s.find(k); return it == p.s.end() ? std::string(d) : it->second; };
    auto F = [&](const char* k, double d) { auto it = p.f.find(k); return it == p.f.end() ? d : it->second; };
    auto I = [&](const char* k, int d) { auto it = p.i.find(k); return it == p.i.end() ? d : it->second; };
    auto B = [&](const char* k, bool d) { auto it = p.b.find(k); return it == p.b.end() ? d : it->second; };
    auto file = [&](const std::string& f) { return f.rfind("/vsi", 0) == 0 ? f : resolve_path(f); };

    lsdtt_river_line_fetch_options o;
    o.source = S("river_line_source");
    if (o.source.empty())
        throw std::invalid_argument("set river_line_source: local_file or gb_os_open_rivers");

    // Area
    const bool edges_given = !S("west").empty() || !S("east").empty() || !S("south").empty() || !S("north").empty();
    if (!S("bounds_fname").empty())
    {
        if (edges_given)
            throw std::invalid_argument("give either bounds_fname or west/south/east/north, not both");
        std::string path = S("bounds_fname");
        if (!std::filesystem::path(path).is_absolute())
            path = S("bounds_path").empty() ? resolve_path(path)
                                            : (std::filesystem::path(S("bounds_path")) / path).string();
        const lsdtt_bounds_from_vector_result b = bounds_from_vector_file(path, S("bounds_layer"));
        if (!b.ok) throw std::invalid_argument("bounds_fname: " + b.message);
        if (!b.warning.empty()) warnings.push_back(b.warning);
        o.west = b.west; o.south = b.south; o.east = b.east; o.north = b.north;
        if (b.geometry_kind == "raster") o.bounds_raster_path = path;
        if (area_note) *area_note = "Area from " + path + " (" + b.geometry_kind + ", " + b.source_crs + ")";
    }
    else if (!parse_degrees_text(S("west"), o.west) || !parse_degrees_text(S("east"), o.east) ||
             !parse_degrees_text(S("south"), o.south) || !parse_degrees_text(S("north"), o.north))
    {
        throw std::invalid_argument("give the area: west, south, east and north (WGS84 degrees), or bounds_fname");
    }

    o.output_epsg = I("output_epsg", 0);
    if (o.output_epsg == 3857 || o.output_epsg == 3785 || o.output_epsg == 900913)
        throw std::invalid_argument("output_epsg: Web Mercator distorts distances and is refused by lsdtt; "
                                    "use 0 or another metric CRS");
    o.clip_buffer_metres = F("clip_buffer_metres", -1.0);
    o.node_snap_tolerance_metres = F("node_snap_tolerance_metres", -1.0);
    o.min_line_length_metres = F("min_line_length_metres", 0.0);
    const std::string gb = S("gap_bridge_algorithm_choice", "none");
    if (gb == "facing_tips") o.gap_bridge.enabled = true;
    else if (gb != "none") throw std::invalid_argument("gap_bridge_algorithm_choice must be none or facing_tips");
    o.gap_bridge.max_gap_length_metres = F("max_gap_length_metres", 25.0);
    o.gap_bridge.max_nonroad_gap_metres = F("max_nonroad_gap_metres", 30.0);
    o.gap_bridge.max_facing_angle_degrees = F("max_facing_angle_degrees", 45.0);
    o.gap_bridge.max_heading_misalignment_degrees = F("max_heading_misalignment_degrees", 50.0);
    o.gap_bridge.road_max_facing_angle_degrees = F("road_max_facing_angle_degrees", 60.0);
    o.gap_bridge.tangent_fit_length_metres = F("tangent_fit_length_metres", 15.0);
    o.gap_bridge.road_intersection_tolerance_metres = F("road_intersection_tolerance_metres", 10.0);
    o.gap_bridge.require_mutual_best = B("require_mutual_best", true);
    if (!S("road_vector_fname").empty())
    {
        if (!o.gap_bridge.enabled)
            warnings.push_back("road_vector_fname is only used with gap_bridge_algorithm_choice: facing_tips; ignored");
        else
        {
            o.road_vector_path = file(S("road_vector_fname"));
            o.road_layer = S("road_layer");
        }
    }
    o.road_name_field = S("road_name_field", "roadName");
    const std::string dr = S("direction_repair_algorithm_choice", "network_outlet");
    if (dr == "network_outlet") o.repair_direction = true;
    else if (dr == "none") o.repair_direction = false;
    else throw std::invalid_argument("direction_repair_algorithm_choice must be network_outlet or none");

    o.keep_feature_classes.clear();
    {
        static const std::set<std::string> known = {"river", "stream", "lake_connector", "canal",
                                                    "drain", "tidal", "other"};
        std::stringstream ss(S("river_feature_classes", "river,stream,lake_connector"));
        std::string c;
        while (std::getline(ss, c, ','))
        {
            c.erase(0, c.find_first_not_of(" \t"));
            c.erase(c.find_last_not_of(" \t") + 1);
            if (c.empty()) continue;
            if (!known.count(c))
                throw std::invalid_argument("river_feature_classes: '" + c + "' is not one of "
                                            "river, stream, lake_connector, canal, drain, tidal, other");
            o.keep_feature_classes.push_back(c);
        }
        if (o.keep_feature_classes.empty()) throw std::invalid_argument("river_feature_classes is empty");
    }

    if (o.source == "local_file")
    {
        if (S("river_lines_fname").empty())
            throw std::invalid_argument("river_line_source local_file needs river_lines_fname");
        o.river_lines_path = file(S("river_lines_fname"));
        o.river_lines_layer = S("river_lines_layer");
        o.fields.feature_class_field = S("feature_class_field_name");
        o.fields.direction_field = S("direction_field_name");
        o.fields.direction_values = parse_direction_field_mapping(S("direction_field_mapping"));
        o.fields.digitised_downstream = B("digitised_downstream", false);
        o.fields.start_node_field = S("start_node_field_name");
        o.fields.end_node_field = S("end_node_field_name");
        o.fields.culvert_field = S("culvert_field_name");
        if (!o.fields.direction_field.empty() && o.fields.digitised_downstream)
            throw std::invalid_argument("give direction_field_name or digitised_downstream, not both");
        if (o.fields.start_node_field.empty() != o.fields.end_node_field.empty())
            throw std::invalid_argument("give both start_node_field_name and end_node_field_name, or neither");
    }
    else
    {
        for (const char* k : {"river_lines_fname", "feature_class_field_name", "direction_field_name",
                              "direction_field_mapping", "start_node_field_name", "end_node_field_name",
                              "culvert_field_name"})
            if (!S(k).empty())
                warnings.push_back(std::string(k) + " is only used with river_line_source: local_file; ignored for " +
                                   o.source);
    }
    if (!S("river_cache_dir").empty()) o.cache_dir = resolve_path(S("river_cache_dir"));

    o.http.timeout_seconds = I("http_timeout_seconds", 60);
    o.http.connect_timeout_seconds = I("http_connect_timeout_seconds", 15);
    o.http.max_retry = I("http_max_retry", 3);
    o.http.user_agent = "lsdtt3 lsdtt-fetch-river-lines/" + lsdtt3_package_version();
    return o;
}

} // namespace lsdtt
