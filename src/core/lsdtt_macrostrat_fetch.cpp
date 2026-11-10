// lsdtt_macrostrat_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_macrostrat_fetch.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include <cpl_conv.h>
#include <cpl_json.h>
#include <cpl_vsi.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_gdal_init.hpp"

namespace lsdtt
{

namespace
{

constexpr double k_pi = 3.14159265358979323846;
constexpr double k_max_lat = 85.0511287798066;
constexpr double k_r = 6378137.0;   // Web Mercator sphere

std::string lower(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string squash(const std::string& s)
{
    std::string out;
    bool space = false;
    for (char c : s)
    {
        if (std::isspace(static_cast<unsigned char>(c))) { space = !out.empty(); continue; }
        if (space) { out += ' '; space = false; }
        out += c;
    }
    return out;
}

CPLJSONArray success_data(const std::string& json, const char* what)
{
    CPLJSONDocument doc;
    if (!doc.LoadMemory(json))
        throw std::runtime_error(std::string("Macrostrat ") + what + ": the response is not JSON");
    const CPLJSONObject root = doc.GetRoot();
    const CPLJSONObject err = root.GetObj("error");
    if (err.IsValid())
        throw std::runtime_error(std::string("Macrostrat ") + what + ": error: " +
                                 err.ToString());
    const CPLJSONArray data = root.GetArray("success/data");
    if (!data.IsValid())
        throw std::runtime_error(std::string("Macrostrat ") + what +
                                 ": the response has no success/data");
    return data;
}

std::string str(const CPLJSONObject& o, const char* key)
{
    const CPLJSONObject v = o.GetObj(key);
    if (!v.IsValid() || v.GetType() == CPLJSONObject::Type::Null) return "";
    if (v.GetType() == CPLJSONObject::Type::String) return v.ToString();
    return v.Format(CPLJSONObject::PrettyFormat::Plain);
}

std::string join(const std::vector<std::string>& v, const char* sep)
{
    std::string out;
    for (const auto& s : v) out += (out.empty() ? "" : sep) + s;
    return out;
}

} // namespace

// ---------------------------------------------------------------------------
// Tiles
// ---------------------------------------------------------------------------
std::vector<lsdtt_slippy_tile> slippy_tiles_for_bbox(double west, double south, double east,
                                                     double north, int z)
{
    if (z < 0 || z > 22) throw std::invalid_argument("tile zoom must be in 0..22");
    if (!(west < east) || !(south < north) || west < -180.0 || east > 180.0 || south < -90.0 ||
        north > 90.0)
        throw std::invalid_argument("bounding box must have west < east and south < north, "
                                    "within -180..180 and -90..90");
    south = std::max(south, -k_max_lat);
    north = std::min(north, k_max_lat);
    const double n = std::ldexp(1.0, z);
    auto tx = [&](double lon) {
        return std::min(static_cast<int>(n) - 1,
                        std::max(0, static_cast<int>(std::floor((lon + 180.0) / 360.0 * n))));
    };
    auto ty = [&](double lat) {
        const double r = lat * k_pi / 180.0;
        const double y = (1.0 - std::asinh(std::tan(r)) / k_pi) / 2.0 * n;
        return std::min(static_cast<int>(n) - 1, std::max(0, static_cast<int>(std::floor(y))));
    };
    std::vector<lsdtt_slippy_tile> out;
    for (int y = ty(north); y <= ty(south); ++y)
        for (int x = tx(west); x <= tx(east); ++x) out.push_back({z, x, y});
    return out;
}

void slippy_tile_bounds_3857(const lsdtt_slippy_tile& t, double& x_min, double& y_min,
                             double& x_max, double& y_max)
{
    const double size = 2.0 * k_pi * k_r / std::ldexp(1.0, t.z);
    const double origin = k_pi * k_r;
    x_min = -origin + t.x * size;
    x_max = x_min + size;
    y_max = origin - t.y * size;
    y_min = y_max - size;
}

int utm_epsg_for(double lon, double lat)
{
    const int zone = std::min(60, std::max(1, static_cast<int>(std::floor((lon + 180.0) / 6.0)) + 1));
    return (lat >= 0.0 ? 32600 : 32700) + zone;
}

// ---------------------------------------------------------------------------
// Parsers
// ---------------------------------------------------------------------------
std::map<int, lsdtt_macrostrat_lithology> parse_macrostrat_lithologies(const std::string& json)
{
    std::map<int, lsdtt_macrostrat_lithology> out;
    const CPLJSONArray data = success_data(json, "lithology definitions");
    for (int i = 0; i < data.Size(); ++i)
    {
        const CPLJSONObject r = data[i];
        lsdtt_macrostrat_lithology l;
        l.id = r.GetInteger("lith_id", -1);
        if (l.id < 0) continue;
        l.name = lower(str(r, "name"));
        l.type = str(r, "type");
        l.group = str(r, "group");
        l.lith_class = str(r, "class");
        out[l.id] = l;
    }
    return out;
}

std::map<std::int64_t, std::vector<int>> parse_macrostrat_map_unit_liths(const std::string& json)
{
    std::map<std::int64_t, std::vector<int>> out;
    const CPLJSONArray data = success_data(json, "map units");
    for (int i = 0; i < data.Size(); ++i)
    {
        const CPLJSONObject r = data[i];
        const std::int64_t id = r.GetLong("map_id", -1);
        if (id < 0) continue;
        std::vector<int> liths;
        const CPLJSONArray a = r.GetArray("liths");
        if (a.IsValid())
            for (int k = 0; k < a.Size(); ++k) liths.push_back(a[k].ToInteger());
        out[id] = liths;
    }
    return out;
}

std::vector<lsdtt_macrostrat_source> parse_macrostrat_sources(const std::string& json)
{
    std::vector<lsdtt_macrostrat_source> out;
    const CPLJSONArray data = success_data(json, "sources");
    for (int i = 0; i < data.Size(); ++i)
    {
        const CPLJSONObject r = data[i];
        lsdtt_macrostrat_source s;
        s.id = r.GetInteger("source_id", -1);
        s.name = str(r, "name");
        s.authors = str(r, "authors");
        s.year = str(r, "ref_year");
        s.url = str(r, "url");
        s.ref_title = str(r, "ref_title");
        s.ref_source = str(r, "ref_source");
        s.doi = str(r, "isbn_doi");
        s.scale = str(r, "scale");
        out.push_back(s);
    }
    std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
    return out;
}

lsdtt_macrostrat_primary macrostrat_primary_lithology(
    const std::string& lith_text, const std::vector<int>& lith_ids,
    const std::map<int, lsdtt_macrostrat_lithology>& vocabulary)
{
    const std::string low = lower(lith_text);
    const auto major = low.find("major:");
    if (major != std::string::npos)
    {
        const auto open = low.find('{', major);
        const auto close = low.find('}', open == std::string::npos ? major : open);
        if (open != std::string::npos && close != std::string::npos && close > open + 1)
        {
            std::string inside = low.substr(open + 1, close - open - 1);
            const std::string first = squash(inside.substr(0, inside.find(',')));
            if (!first.empty()) return {first, "major"};
        }
    }
    for (int id : lith_ids)
    {
        const auto it = vocabulary.find(id);
        if (it != vocabulary.end() && !it->second.name.empty()) return {it->second.name, "liths"};
    }
    const std::string text = squash(low);
    if (!text.empty()) return {text, "text"};
    return {"", "none"};
}

// ---------------------------------------------------------------------------
// The fetch
// ---------------------------------------------------------------------------
namespace
{

struct unit_feature
{
    std::unique_ptr<OGRGeometry> geom;   // EPSG:3857
    std::int64_t map_id = -1;
    int source_id = -1;
    std::string name, lith, ref_name;
    double age_top = std::nan(""), age_bottom = std::nan("");
};

std::string tile_url(const std::string& tmpl, const lsdtt_slippy_tile& t)
{
    std::string u = tmpl;
    auto rep = [&](const std::string& k, int v) {
        for (auto p = u.find(k); p != std::string::npos; p = u.find(k)) u.replace(p, k.size(), std::to_string(v));
    };
    rep("{z}", t.z);
    rep("{x}", t.x);
    rep("{y}", t.y);
    return u;
}

void pause(int ms)
{
    if (ms > 0) CPLSleep(ms / 1000.0);
}

bool get_with_delay(const std::string& url, const lsdtt_macrostrat_fetch_options& o,
                    std::string& body, std::string& error)
{
    const bool ok = http_request_body(url, "", o.http, body, error);
    pause(o.request_delay_milliseconds);
    return ok;
}

} // namespace

lsdtt_macrostrat_fetch_result fetch_macrostrat_units(const lsdtt_macrostrat_fetch_options& o)
{
    gdal_init();
    lsdtt_macrostrat_fetch_result res;
    const std::vector<lsdtt_slippy_tile> tiles =
        slippy_tiles_for_bbox(o.west, o.south, o.east, o.north, o.zoom);
    if (static_cast<int>(tiles.size()) > o.max_tiles)
        throw std::runtime_error(
            "the bounding box needs " + std::to_string(tiles.size()) + " tiles at zoom " +
            std::to_string(o.zoom) + ", more than max_tiles = " + std::to_string(o.max_tiles) +
            ". Use a smaller box, a lower macrostrat_tile_zoom, or raise max_tiles");
    if (o.tile_dir.empty() && !http_fetch_available())
        throw std::runtime_error("this GDAL has no HTTP support (built without curl)");
    if (!GetGDALDriverManager()->GetDriverByName("MVT"))
        throw std::runtime_error("this GDAL has no MVT (Mapbox Vector Tile) driver, which "
                                 "reading Macrostrat tiles needs");
    res.n_tiles = tiles.size();

    // The box in EPSG:3857, to drop features that lie wholly outside it.
    OGRSpatialReference wgs84, merc;
    wgs84.importFromEPSG(4326);
    merc.importFromEPSG(3857);
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    merc.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    std::unique_ptr<OGRCoordinateTransformation> to_merc(
        OGRCreateCoordinateTransformation(&wgs84, &merc));
    double bx[2] = {o.west, o.east}, by[2] = {std::max(o.south, -k_max_lat), std::min(o.north, k_max_lat)};
    to_merc->Transform(2, bx, by);
    OGREnvelope box;
    box.MinX = bx[0]; box.MaxX = bx[1]; box.MinY = by[0]; box.MaxY = by[1];

    std::vector<unit_feature> features;
    for (const auto& t : tiles)
    {
        std::string path;
        std::string mem;
        if (!o.tile_dir.empty())
        {
            for (const char* ext : {".mvt", ".pbf"})
            {
                const std::string p = o.tile_dir + "/" + std::to_string(t.z) + "/" +
                                      std::to_string(t.x) + "/" + std::to_string(t.y) + ext;
                if (std::filesystem::exists(p)) { path = p; break; }
            }
            if (path.empty()) { ++res.n_tiles_missing; continue; }
        }
        else
        {
            std::string body, error;
            if (!get_with_delay(tile_url(o.tile_url_template, t), o, body, error))
            {
                ++res.n_tiles_missing;
                res.warnings.push_back("tile " + std::to_string(t.z) + "/" + std::to_string(t.x) +
                                       "/" + std::to_string(t.y) + ": " + error);
                continue;
            }
            // The /vsimem path carries z/x/y, which is how the MVT driver
            // georeferences a single tile.
            mem = "/vsimem/lsdtt_macrostrat/" + std::to_string(t.z) + "/" + std::to_string(t.x) +
                  "/" + std::to_string(t.y) + ".mvt";
            VSILFILE* w = VSIFOpenL(mem.c_str(), "wb");
            if (!w) throw std::runtime_error("could not create " + mem);
            VSIFWriteL(body.data(), 1, body.size(), w);
            VSIFCloseL(w);
            path = mem;
        }
        const char* drivers[] = {"MVT", nullptr};
        GDALDataset* ds = static_cast<GDALDataset*>(
            GDALOpenEx(path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, drivers, nullptr, nullptr));
        if (!ds)
        {
            ++res.n_tiles_missing;
            res.warnings.push_back("tile " + path + " could not be read as MVT");
            if (!mem.empty()) VSIUnlink(mem.c_str());
            continue;
        }
        OGRLayer* layer = ds->GetLayerByName("units");
        if (layer)
        {
            OGRFeature* feat = nullptr;
            layer->ResetReading();
            while ((feat = layer->GetNextFeature()) != nullptr)
            {
                std::unique_ptr<OGRFeature, void (*)(OGRFeature*)> fg(feat, OGRFeature::DestroyFeature);
                OGRGeometry* g = feat->GetGeometryRef();
                if (!g) continue;
                OGREnvelope env;
                g->getEnvelope(&env);
                if (!env.Intersects(box)) continue;
                unit_feature u;
                u.geom.reset(g->clone());
                const OGRFeatureDefn* d = feat->GetDefnRef();
                auto idx = [&](const char* n) { return d->GetFieldIndex(n); };
                if (idx("map_id") >= 0) u.map_id = feat->GetFieldAsInteger64(idx("map_id"));
                if (idx("source_id") >= 0) u.source_id = feat->GetFieldAsInteger(idx("source_id"));
                if (idx("name") >= 0) u.name = feat->GetFieldAsString(idx("name"));
                if (idx("lith") >= 0) u.lith = feat->GetFieldAsString(idx("lith"));
                if (idx("ref_name") >= 0) u.ref_name = feat->GetFieldAsString(idx("ref_name"));
                if (idx("best_age_top") >= 0 && feat->IsFieldSetAndNotNull(idx("best_age_top")))
                    u.age_top = feat->GetFieldAsDouble(idx("best_age_top"));
                if (idx("best_age_bottom") >= 0 && feat->IsFieldSetAndNotNull(idx("best_age_bottom")))
                    u.age_bottom = feat->GetFieldAsDouble(idx("best_age_bottom"));
                features.push_back(std::move(u));
            }
        }
        GDALClose(ds);
        if (!mem.empty()) VSIUnlink(mem.c_str());
    }
    if (res.n_tiles_missing == res.n_tiles)
        throw std::runtime_error("no Macrostrat tile could be read (" +
                                 std::to_string(res.n_tiles) + " tiles tried)" +
                                 (res.warnings.empty() ? "" : ": " + res.warnings.front()));
    if (features.empty())
        throw std::runtime_error("the Macrostrat tiles hold no map units inside the bounding box");
    res.n_features = features.size();

    // Map-unit lookups: liths per map_id, the vocabulary, and the sources.
    std::set<std::int64_t> map_ids;
    std::set<int> source_ids;
    for (const auto& f : features)
    {
        if (f.map_id >= 0) map_ids.insert(f.map_id);
        if (f.source_id >= 0) source_ids.insert(f.source_id);
    }
    res.n_map_units = map_ids.size();
    std::map<int, lsdtt_macrostrat_lithology> vocabulary;
    std::map<std::int64_t, std::vector<int>> liths;
    if (o.lookup_map_units)
    {
        std::string body, error;
        try
        {
            if (!get_with_delay(o.api_base + "/defs/lithologies?all", o, body, error))
                throw std::runtime_error(error);
            vocabulary = parse_macrostrat_lithologies(body);
            std::vector<std::int64_t> ids(map_ids.begin(), map_ids.end());
            for (std::size_t i = 0; i < ids.size(); i += 100)
            {
                std::string list;
                for (std::size_t k = i; k < std::min(ids.size(), i + 100); ++k)
                    list += (list.empty() ? "" : ",") + std::to_string(ids[k]);
                if (!get_with_delay(o.api_base + "/geologic_units/map?map_id=" + list, o, body, error))
                    throw std::runtime_error(error);
                for (auto& [id, l] : parse_macrostrat_map_unit_liths(body)) liths[id] = l;
            }
            std::string list;
            for (int s : source_ids) list += (list.empty() ? "" : ",") + std::to_string(s);
            if (!list.empty())
            {
                if (!get_with_delay(o.api_base + "/defs/sources?source_id=" + list, o, body, error))
                    throw std::runtime_error(error);
                res.sources = parse_macrostrat_sources(body);
            }
        }
        catch (const std::exception& e)
        {
            res.warnings.push_back(std::string("Macrostrat API lookup failed (") + e.what() +
                                   "); lith_primary falls back to the tiles' lith text and the "
                                   "attribution lists source ids only");
        }
    }
    if (res.sources.empty())
        for (int s : source_ids)
        {
            lsdtt_macrostrat_source src;
            src.id = s;
            res.sources.push_back(src);
        }

    // Reproject to UTM (or the requested CRS) and write.
    res.epsg = o.target_epsg > 0
                   ? o.target_epsg
                   : utm_epsg_for(0.5 * (o.west + o.east), 0.5 * (o.south + o.north));
    if (o.target_epsg <= 0 && o.east - o.west > 6.0)
        res.warnings.push_back("the bounding box is wider than one UTM zone (6 degrees); "
                               "distances far from EPSG:" + std::to_string(res.epsg) +
                               "'s central meridian are distorted");
    OGRSpatialReference out_srs;
    if (out_srs.importFromEPSG(res.epsg) != OGRERR_NONE)
        throw std::runtime_error("unknown EPSG:" + std::to_string(res.epsg));
    out_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    std::unique_ptr<OGRCoordinateTransformation> to_out(
        OGRCreateCoordinateTransformation(&merc, &out_srs));

    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName(o.out_driver.c_str());
    if (!drv) throw std::runtime_error("no OGR driver " + o.out_driver);
    std::filesystem::remove(o.out_path);
    GDALDataset* out = drv->Create(o.out_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!out) throw std::runtime_error("could not create " + o.out_path);
    std::unique_ptr<GDALDataset, void (*)(GDALDataset*)> out_guard(out, [](GDALDataset* d) { GDALClose(d); });
    OGRLayer* layer = out->CreateLayer("units", &out_srs, wkbMultiPolygon, nullptr);
    if (!layer) throw std::runtime_error("could not create layer 'units' in " + o.out_path);
    auto field = [&](const char* n, OGRFieldType t) {
        OGRFieldDefn fd(n, t);
        if (layer->CreateField(&fd) != OGRERR_NONE)
            throw std::runtime_error(std::string("could not create field '") + n +
                                     "' in " + o.out_path);
    };
    field("map_id", OFTInteger64);
    field("source_id", OFTInteger);
    field("name", OFTString);
    field("lith", OFTString);
    field("lith_primary", OFTString);
    field("lith_primary_source", OFTString);
    field("lith_ids", OFTString);
    field("lith_names", OFTString);
    field("lith_type", OFTString);
    field("lith_class", OFTString);
    field("best_age_top", OFTReal);
    field("best_age_bottom", OFTReal);
    field("ref_name", OFTString);

    std::map<std::string, std::string> name_to_type, name_to_class;
    for (const auto& [id, l] : vocabulary)
    {
        name_to_type[l.name] = l.type;
        name_to_class[l.name] = l.lith_class;
    }
    std::map<std::string, lsdtt_macrostrat_vocabulary_row> vocab_rows;
    std::size_t n_written = 0;
    for (auto& f : features)
    {
        const auto li = liths.find(f.map_id);
        const std::vector<int> ids = li == liths.end() ? std::vector<int>() : li->second;
        const lsdtt_macrostrat_primary p = macrostrat_primary_lithology(f.lith, ids, vocabulary);
        std::vector<std::string> id_s, names;
        for (int id : ids)
        {
            id_s.push_back(std::to_string(id));
            const auto v = vocabulary.find(id);
            if (v != vocabulary.end()) names.push_back(v->second.name);
        }
        std::unique_ptr<OGRGeometry> g(OGRGeometryFactory::forceToMultiPolygon(f.geom.release()));
        if (!g || g->transform(to_out.get()) != OGRERR_NONE) continue;
        OGRFeature* of = OGRFeature::CreateFeature(layer->GetLayerDefn());
        std::unique_ptr<OGRFeature, void (*)(OGRFeature*)> og(of, OGRFeature::DestroyFeature);
        of->SetField("map_id", static_cast<GIntBig>(f.map_id));
        of->SetField("source_id", f.source_id);
        of->SetField("name", f.name.c_str());
        of->SetField("lith", f.lith.c_str());
        of->SetField("lith_primary", p.value.c_str());
        of->SetField("lith_primary_source", p.source.c_str());
        of->SetField("lith_ids", join(id_s, ",").c_str());
        of->SetField("lith_names", join(names, ",").c_str());
        const std::string type = name_to_type.count(p.value) ? name_to_type[p.value] : "";
        const std::string cls = name_to_class.count(p.value) ? name_to_class[p.value] : "";
        of->SetField("lith_type", type.c_str());
        of->SetField("lith_class", cls.c_str());
        if (!std::isnan(f.age_top)) of->SetField("best_age_top", f.age_top);
        if (!std::isnan(f.age_bottom)) of->SetField("best_age_bottom", f.age_bottom);
        of->SetField("ref_name", f.ref_name.c_str());
        double area = 0.0;
        if (auto* mp = dynamic_cast<OGRMultiPolygon*>(g.get())) area = mp->get_Area();
        of->SetGeometryDirectly(g.release());
        if (layer->CreateFeature(of) != OGRERR_NONE)
            throw std::runtime_error("could not write a feature to " + o.out_path);
        ++n_written;
        auto& row = vocab_rows[p.value];
        row.lith_primary = p.value;
        row.lith_primary_source = p.source;
        row.lith_type = type;
        row.lith_class = cls;
        ++row.n_features;
        row.area_km2 += area / 1.0e6;
    }
    res.n_features = n_written;
    for (auto& [k, row] : vocab_rows) res.vocabulary.push_back(row);
    std::sort(res.vocabulary.begin(), res.vocabulary.end(),
              [](const auto& a, const auto& b) { return a.area_km2 > b.area_km2; });
    res.notes.push_back("units overlap where tiles' buffers overlap; polygons of one unit "
                        "are kept as separate features (a rasteriser adds them up)");
    return res;
}

std::string macrostrat_attribution_text(const lsdtt_macrostrat_fetch_result& r)
{
    std::ostringstream s;
    s << "Geological map data from Macrostrat (https://macrostrat.org), licensed CC BY 4.0.\n"
      << "Cite Macrostrat: Peters, S. E., Husson, J. M., & Czaplewski, J. (2018). Macrostrat: A\n"
      << "platform for geological data integration and deep-time Earth crust research.\n"
      << "Geochemistry, Geophysics, Geosystems, 19, 1393-1409. doi:10.1029/2018GC007467\n\n"
      << "and each source map used:\n";
    for (const auto& src : r.sources)
    {
        s << "\n- source_id " << src.id;
        if (!src.name.empty()) s << ": " << src.name;
        s << "\n";
        if (!src.authors.empty() || !src.year.empty())
            s << "  " << src.authors << (src.year.empty() ? "" : " (" + src.year + ")") << "\n";
        if (!src.ref_title.empty()) s << "  " << src.ref_title << "\n";
        if (!src.ref_source.empty()) s << "  " << src.ref_source << "\n";
        if (!src.doi.empty()) s << "  " << src.doi << "\n";
        if (!src.url.empty()) s << "  " << src.url << "\n";
    }
    return s.str();
}

} // namespace lsdtt
