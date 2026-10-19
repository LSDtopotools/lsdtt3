// lsdtt_octopus_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_octopus_fetch.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

#include <cpl_conv.h>
#include <cpl_vsi.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_gdal_init.hpp"
#include "lsdtt_macrostrat_fetch.hpp"   // utm_epsg_for
#include "lsdtt_secret.hpp"             // mask_url_query

namespace lsdtt
{

namespace
{

constexpr const char* k_endpoint = "http://geoserver.octopusdata.org/geoserver/wfs";
constexpr const char* k_basins_layer = "be10-denude:crn_int_basins";
constexpr const char* k_outlets_layer = "be10-denude:crn_int_outlets";
constexpr double k_nan = std::numeric_limits<double>::quiet_NaN();

// Shortest decimal form that reads back to the same double.
std::string shortest(double x)
{
    for (int p = 1; p <= 17; ++p)
    {
        std::ostringstream s;
        s << std::setprecision(p) << x;
        if (std::stod(s.str()) == x) return s.str();
    }
    std::ostringstream s;
    s << std::setprecision(17) << x;
    return s.str();
}

std::string num(double v, int precision = 10)
{
    if (!std::isfinite(v)) return "";
    std::ostringstream s;
    s << std::setprecision(precision) << v;
    return s.str();
}

std::string csv_text(const std::string& s)
{
    if (s.find_first_of(",\"\n\r") == std::string::npos) return s;
    std::string q = "\"";
    for (char c : s) { if (c == '"') q += "\"\""; else q += c; }
    return q + "\"";
}

std::string trim(const std::string& s)
{
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

std::string xml_unescape(std::string s)
{
    const std::pair<const char*, const char*> ents[] = {
        {"&lt;", "<"}, {"&gt;", ">"}, {"&quot;", "\""}, {"&apos;", "'"}, {"&amp;", "&"}};
    for (const auto& [from, to] : ents)
        for (auto p = s.find(from); p != std::string::npos; p = s.find(from, p + std::strlen(to)))
            s.replace(p, std::strlen(from), to);
    return s;
}

std::string first_bytes(const std::string& body)
{
    std::string head = body.substr(0, 200);
    for (char& c : head) if (c == '\n' || c == '\r') c = ' ';
    return head;
}

std::string today_utc()
{
    const std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[16];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d", &tm);
    return buf;
}

std::string url_encode_value(const std::string& v)
{
    std::ostringstream o;
    for (unsigned char c : v)
    {
        if (std::isalnum(c) || c == '-' || c == '_' || c == '.' || c == '~')
            o << c;
        else
            o << '%' << std::uppercase << std::hex << std::setw(2) << std::setfill('0')
              << static_cast<int>(c) << std::nouppercase << std::dec;
    }
    return o.str();
}

// ---------------------------------------------------------------------------
// One GeoJSON body -> features
// ---------------------------------------------------------------------------
struct parsed_feature
{
    std::map<std::string, lsdtt_octopus_value> values;
    std::vector<unsigned char> wkb;
    double x = k_nan, y = k_nan;   // first point, for outlets
};

struct parsed_page
{
    std::vector<std::string> field_order;
    std::vector<parsed_feature> features;
};

void check_body(const std::string& body, const std::string& what)
{
    const std::string ex = octopus_wfs_exception_text(body);
    if (!ex.empty()) throw std::runtime_error("OCTOPUS WFS error: " + ex);
    const std::string t = trim(body.substr(0, 64));
    if (t.empty() || (t[0] != '{' && t[0] != '<'))
        throw std::runtime_error("OCTOPUS WFS: the " + what +
                                 " response is neither JSON nor XML; it starts: " + first_bytes(body));
    if (t[0] == '<')
        throw std::runtime_error("OCTOPUS WFS: the " + what +
                                 " response is XML, not GeoJSON; it starts: " + first_bytes(body));
}

parsed_page parse_geojson_page(const std::string& body, const std::string& mem_name)
{
    parsed_page page;
    const std::string mem = "/vsimem/lsdtt_octopus/" + mem_name + ".json";
    VSILFILE* w = VSIFOpenL(mem.c_str(), "wb");
    if (!w) throw std::runtime_error("could not create " + mem);
    VSIFWriteL(body.data(), 1, body.size(), w);
    VSIFCloseL(w);
    const char* drivers[] = {"GeoJSON", nullptr};
    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpenEx(mem.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, drivers, nullptr, nullptr));
    if (!ds)
    {
        VSIUnlink(mem.c_str());
        throw std::runtime_error("OCTOPUS WFS: a response could not be read as GeoJSON; it starts: " +
                                 first_bytes(body));
    }
    std::unique_ptr<GDALDataset, void (*)(GDALDataset*)> guard(ds, [](GDALDataset* d) { GDALClose(d); });
    if (ds->GetLayerCount() > 0)
    {
        OGRLayer* layer = ds->GetLayer(0);
        OGRFeatureDefn* defn = layer->GetLayerDefn();
        for (int i = 0; i < defn->GetFieldCount(); ++i)
            page.field_order.push_back(defn->GetFieldDefn(i)->GetNameRef());
        layer->ResetReading();
        OGRFeature* f = nullptr;
        while ((f = layer->GetNextFeature()) != nullptr)
        {
            std::unique_ptr<OGRFeature, void (*)(OGRFeature*)> fg(f, OGRFeature::DestroyFeature);
            parsed_feature pf;
            for (int i = 0; i < defn->GetFieldCount(); ++i)
            {
                const OGRFieldDefn* fd = defn->GetFieldDefn(i);
                lsdtt_octopus_value v;
                if (f->IsFieldSetAndNotNull(i))
                {
                    switch (fd->GetType())
                    {
                        case OFTInteger:
                        case OFTInteger64:
                            v.type = lsdtt_octopus_value::kind::integer;
                            v.i = f->GetFieldAsInteger64(i);
                            v.d = static_cast<double>(v.i);
                            break;
                        case OFTReal:
                            v.type = lsdtt_octopus_value::kind::real;
                            v.d = f->GetFieldAsDouble(i);
                            break;
                        default:
                            v.type = lsdtt_octopus_value::kind::text;
                            v.s = f->GetFieldAsString(i);
                            break;
                    }
                }
                pf.values[fd->GetNameRef()] = v;
            }
            if (const OGRGeometry* g = f->GetGeometryRef())
            {
                pf.wkb.resize(static_cast<std::size_t>(g->WkbSize()));
                g->exportToWkb(wkbNDR, pf.wkb.data());
                if (wkbFlatten(g->getGeometryType()) == wkbPoint)
                {
                    const auto* p = g->toPoint();
                    pf.x = p->getX();
                    pf.y = p->getY();
                }
            }
            page.features.push_back(std::move(pf));
        }
    }
    guard.reset();
    VSIUnlink(mem.c_str());
    return page;
}

std::string read_file(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::invalid_argument("octopus_response_dir: cannot read " + path);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

struct http_session
{
    const lsdtt_octopus_fetch_options& o;
    lsdtt_octopus_fetch_result& r;
    std::string get(const std::string& url)
    {
        std::string body, error;
        const bool ok = http_request_body(url, "", o.http, body, error);
        ++r.n_requests;
        if (o.request_delay_milliseconds > 0) CPLSleep(o.request_delay_milliseconds / 1000.0);
        if (!ok)
            throw std::runtime_error("OCTOPUS WFS unreachable (" + mask_url_query(url) + "): " + error);
        r.n_bytes += body.size();
        return body;
    }
};

std::int64_t hits(http_session& http, const std::string& layer, const std::string& bbox)
{
    const std::string body = http.get(octopus_wfs_hits_url(layer, bbox));
    const std::string ex = octopus_wfs_exception_text(body);
    if (!ex.empty()) throw std::runtime_error("OCTOPUS WFS error: " + ex);
    const std::int64_t n = parse_octopus_wfs_number_matched(body);
    if (n < 0)
        throw std::runtime_error("OCTOPUS WFS: the hits response has no numberMatched; it starts: " +
                                 first_bytes(body));
    return n;
}

std::vector<parsed_feature> fetch_all_pages(http_session& http, const std::string& layer,
                                            const std::string& bbox, std::int64_t n_matched,
                                            int page_size, std::vector<std::string>& field_order,
                                            const std::string& tag)
{
    std::vector<parsed_feature> out;
    std::size_t k = 0;
    for (std::int64_t start : octopus_wfs_page_starts(n_matched, page_size))
    {
        const std::string body = http.get(octopus_wfs_page_url(layer, bbox, page_size, start));
        check_body(body, tag);
        parsed_page p = parse_geojson_page(body, tag + "_" + std::to_string(k++));
        ++http.r.n_pages;
        for (const auto& f : p.field_order)
            if (std::find(field_order.begin(), field_order.end(), f) == field_order.end())
                field_order.push_back(f);
        const std::size_t got = p.features.size();
        for (auto& f : p.features) out.push_back(std::move(f));
        if (static_cast<int>(got) < page_size) break;
    }
    if (static_cast<std::int64_t>(out.size()) != n_matched)
        http.r.warnings.push_back("OCTOPUS layer " + layer + ": numberMatched was " +
                                  std::to_string(n_matched) + " but " + std::to_string(out.size()) +
                                  " features were read (the layer changed during paging?)");
    return out;
}

std::string value_text(const lsdtt_octopus_value& v)
{
    switch (v.type)
    {
        case lsdtt_octopus_value::kind::integer: return std::to_string(v.i);
        case lsdtt_octopus_value::kind::real: return shortest(v.d);
        case lsdtt_octopus_value::kind::text: return v.s;
        default: return "";
    }
}

int field_index(const lsdtt_octopus_fetch_result& r, const std::string& name)
{
    for (std::size_t i = 0; i < r.fields.size(); ++i)
        if (r.fields[i].name == name) return static_cast<int>(i);
    return -1;
}

} // namespace

// ---------------------------------------------------------------------------
// Constants and pure recipe
// ---------------------------------------------------------------------------
std::string octopus_wfs_endpoint() { return k_endpoint; }
std::string octopus_basins_layer() { return k_basins_layer; }
std::string octopus_outlets_layer() { return k_outlets_layer; }

std::string octopus_wfs_bbox_parameter(double west, double south, double east, double north)
{
    return shortest(south) + "," + shortest(west) + "," + shortest(north) + "," + shortest(east) +
           ",urn:ogc:def:crs:EPSG::4326";
}

std::string octopus_wfs_hits_url(const std::string& layer, const std::string& bbox)
{
    return std::string(k_endpoint) + "?service=WFS&version=2.0.0&request=GetFeature&typeNames=" +
           layer + "&resultType=hits&bbox=" + bbox + "&srsName=EPSG:4326";
}

std::string octopus_wfs_page_url(const std::string& layer, const std::string& bbox, int count,
                                 std::int64_t start_index)
{
    return std::string(k_endpoint) + "?service=WFS&version=2.0.0&request=GetFeature&typeNames=" +
           layer + "&outputFormat=application/json&srsName=EPSG:4326&bbox=" + bbox +
           "&sortBy=OBSID1+A&count=" + std::to_string(count) +
           "&startIndex=" + std::to_string(start_index);
}

std::string octopus_wfs_outlets_by_id_url(const std::vector<std::string>& obsid1)
{
    std::string list;
    for (const auto& id : obsid1)
    {
        std::string q;
        for (char c : id) { if (c == '\'') q += "''"; else q += c; }
        list += (list.empty() ? "" : ",") + ("'" + q + "'");
    }
    return std::string(k_endpoint) + "?service=WFS&version=2.0.0&request=GetFeature&typeNames=" +
           k_outlets_layer + "&outputFormat=application/json&srsName=EPSG:4326&sortBy=OBSID1+A" +
           "&count=" + std::to_string(obsid1.size()) +
           "&CQL_FILTER=" + url_encode_value("OBSID1 IN (" + list + ")");
}

std::vector<std::int64_t> octopus_wfs_page_starts(std::int64_t n_matched, int page_size)
{
    if (page_size <= 0) throw std::invalid_argument("octopus_page_size must be > 0");
    std::vector<std::int64_t> out;
    for (std::int64_t s = 0; s < n_matched; s += page_size) out.push_back(s);
    return out;
}

std::int64_t parse_octopus_wfs_number_matched(const std::string& body)
{
    static const std::regex xml_re("numberMatched\\s*=\\s*\"([0-9]+)\"");
    static const std::regex json_re("\"numberMatched\"\\s*:\\s*([0-9]+)");
    std::smatch m;
    if (std::regex_search(body, m, xml_re)) return std::stoll(m[1].str());
    if (std::regex_search(body, m, json_re)) return std::stoll(m[1].str());
    return -1;
}

std::string octopus_wfs_exception_text(const std::string& body)
{
    const std::string head = body.substr(0, 4096);
    if (head.find("ExceptionReport") == std::string::npos) return "";
    static const std::regex text_re(
        "<(?:[A-Za-z0-9_]+:)?ExceptionText[^>]*>([\\s\\S]*?)</(?:[A-Za-z0-9_]+:)?ExceptionText>");
    std::string out;
    for (auto it = std::sregex_iterator(body.begin(), body.end(), text_re); it != std::sregex_iterator();
         ++it)
    {
        const std::string t = trim(xml_unescape((*it)[1].str()));
        if (!t.empty()) out += (out.empty() ? "" : "; ") + t;
    }
    return out.empty() ? "an ExceptionReport without ExceptionText" : out;
}

bool is_octopus_no_data_sentinel(double x)
{
    if (!std::isfinite(x) || !(x < 0.0)) return false;
    const double a = -x;
    for (int k = 1; k <= 6; ++k)          // nines before the decimal point
        for (int d = 0; d <= 6; ++d)      // nines after it
        {
            if (k + d < 3) continue;
            const double s = std::pow(10.0, k) - std::pow(10.0, -d);
            if (std::fabs(a - s) <= 1e-6) return true;
        }
    return false;
}

lsdtt_octopus_rate_fields octopus_rate_fields(const std::string& scheme, const std::string& nuclide)
{
    if (nuclide != "be10" && nuclide != "al26")
        throw std::invalid_argument("octopus_nuclide must be be10 or al26 (got '" + nuclide + "')");
    const bool be = nuclide == "be10";
    lsdtt_octopus_rate_fields f;
    if (scheme == "st")
    {
        f.rate = be ? "EBE_MMKYR" : "EAL_MMKYR";
        f.uncertainty = be ? "EBE_ERR" : "EAL_ERR";
        f.internal_uncertainty = be ? "ERRBE_AMS" : "ERRAL_AMS";
        f.internal_mass_flux_rate = be ? "EBE_GCMYR" : "EAL_GCMYR";
        f.kind = "total_1sigma";
    }
    else if (scheme == "lsdn")
    {
        f.rate = be ? "EBEL" : "EALL";
        f.uncertainty = be ? "EBEL_EXT" : "EALL_EXT";
        f.internal_uncertainty = be ? "EBEL_INT" : "EALL_INT";
        f.kind = "external_1sigma";
    }
    else if (scheme == "lsdn_no_topo_shielding")
    {
        f.rate = be ? "EBELNS" : "EALLNS";
        f.uncertainty = be ? "EBELNS_EXT" : "EALLNS_EXT";
        f.internal_uncertainty = be ? "EBELNS_INT" : "EALLNS_INT";
        f.kind = "external_1sigma";
    }
    else if (scheme == "published")
    {
        f.rate = be ? "BE10EP" : "AL26EP";
        f.uncertainty = be ? "BE10EP_ERR" : "AL26EP_ERR";
        f.kind = "published_1sigma";
    }
    else
        throw std::invalid_argument("octopus_rate_scheme must be one of st, lsdn, "
                                    "lsdn_no_topo_shielding, published (got '" + scheme + "')");
    return f;
}

// ---------------------------------------------------------------------------
// The fetch
// ---------------------------------------------------------------------------
lsdtt_octopus_fetch_result fetch_octopus_erosion_rates(const lsdtt_octopus_fetch_options& o)
{
    gdal_init();
    if (!(o.west < o.east) || !(o.south < o.north) || o.west < -180.0 || o.east > 180.0 ||
        o.south < -90.0 || o.north > 90.0)
        throw std::invalid_argument("bounding box must have west < east and south < north, within "
                                    "-180..180 and -90..90");
    if (o.page_size <= 0) throw std::invalid_argument("octopus_page_size must be > 0");
    lsdtt_octopus_fetch_result r;
    r.rate_fields = octopus_rate_fields(o.rate_scheme, o.nuclide);
    r.rate_scheme = o.rate_scheme;
    r.nuclide = o.nuclide;
    r.west = o.west; r.south = o.south; r.east = o.east; r.north = o.north;
    r.access_date = today_utc();
    if (!GetGDALDriverManager()->GetDriverByName("GeoJSON"))
        throw std::runtime_error("this GDAL has no GeoJSON driver, which reading OCTOPUS needs");

    const std::string bbox = octopus_wfs_bbox_parameter(o.west, o.south, o.east, o.north);
    std::vector<std::string> field_order;
    std::vector<parsed_feature> basins, outlets;
    http_session http{o, r};

    auto max_check = [&](std::int64_t n) {
        if (n > o.max_features)
            throw std::runtime_error("OCTOPUS returns " + std::to_string(n) +
                                     " basins for this box, more than max_features (" +
                                     std::to_string(o.max_features) +
                                     "). Use a smaller box or raise max_features");
    };

    if (!o.response_dir.empty())
    {
        r.source = o.response_dir;
        const std::string hits_path = o.response_dir + "/hits_basins.xml";
        if (std::filesystem::exists(hits_path))
        {
            const std::string h = read_file(hits_path);
            const std::string ex = octopus_wfs_exception_text(h);
            if (!ex.empty()) throw std::runtime_error("OCTOPUS WFS error: " + ex);
            const std::int64_t n = parse_octopus_wfs_number_matched(h);
            if (n >= 0) max_check(n);
        }
        const std::string b = read_file(o.response_dir + "/crn_int_basins.json");
        r.n_bytes += b.size();
        check_body(b, "basins");
        parsed_page pb = parse_geojson_page(b, "basins_offline");
        field_order = pb.field_order;
        basins = std::move(pb.features);
        const std::int64_t nm = parse_octopus_wfs_number_matched(b);
        r.n_matched = nm >= 0 ? nm : static_cast<std::int64_t>(basins.size());
        max_check(r.n_matched);
        const std::string op = o.response_dir + "/crn_int_outlets.json";
        if (std::filesystem::exists(op))
        {
            const std::string ob = read_file(op);
            r.n_bytes += ob.size();
            check_body(ob, "outlets");
            outlets = parse_geojson_page(ob, "outlets_offline").features;
        }
        r.n_pages = 2;
    }
    else
    {
        r.source = k_endpoint;
        if (!http_fetch_available())
            throw std::runtime_error("this GDAL has no HTTP support (built without curl)");
        r.n_matched = hits(http, k_basins_layer, bbox);
        max_check(r.n_matched);
        if (r.n_matched > 0)
        {
            basins = fetch_all_pages(http, k_basins_layer, bbox, r.n_matched, o.page_size,
                                     field_order, "basins");
            const std::int64_t n_out = hits(http, k_outlets_layer, bbox);
            std::vector<std::string> outlet_fields;
            if (n_out > 0)
                outlets = fetch_all_pages(http, k_outlets_layer, bbox, n_out, o.page_size,
                                          outlet_fields, "outlets");
        }
    }

    // The schema: the union of the pages' fields, typed by the widest
    // non-null value (text > real > integer).
    for (const auto& name : field_order)
    {
        bool any_text = false, any_real = false, any_int = false;
        for (const auto& f : basins)
        {
            const auto it = f.values.find(name);
            if (it == f.values.end()) continue;
            any_text |= it->second.type == lsdtt_octopus_value::kind::text;
            any_real |= it->second.type == lsdtt_octopus_value::kind::real;
            any_int |= it->second.type == lsdtt_octopus_value::kind::integer;
        }
        lsdtt_octopus_field fd;
        fd.name = name;
        fd.type = any_text ? lsdtt_octopus_value::kind::text
                  : any_real ? lsdtt_octopus_value::kind::real
                  : any_int  ? lsdtt_octopus_value::kind::integer
                             : lsdtt_octopus_value::kind::text;
        r.fields.push_back(fd);
    }

    // Outlets by OBSID1 (the first of repeated ids).
    std::map<std::string, std::pair<double, double>> outlet_at;
    for (const auto& f : outlets)
    {
        const auto it = f.values.find("OBSID1");
        if (it == f.values.end() || it->second.type == lsdtt_octopus_value::kind::null) continue;
        if (!std::isfinite(f.x) || !std::isfinite(f.y)) continue;
        outlet_at.emplace(value_text(it->second), std::make_pair(f.x, f.y));
    }
    r.n_outlets = outlet_at.size();

    // Basins: mask the sentinels, join the outlets.
    std::vector<std::string> missing_outlets;
    for (auto& f : basins)
    {
        lsdtt_octopus_basin b;
        for (const auto& fd : r.fields)
        {
            const auto it = f.values.find(fd.name);
            lsdtt_octopus_value v = it == f.values.end() ? lsdtt_octopus_value{} : it->second;
            if ((v.type == lsdtt_octopus_value::kind::integer ||
                 v.type == lsdtt_octopus_value::kind::real) &&
                is_octopus_no_data_sentinel(v.d))
            {
                v = lsdtt_octopus_value{};
                ++r.masked_by_field[fd.name];
            }
            b.values.push_back(v);
        }
        b.wkb = std::move(f.wkb);
        const std::string id = octopus_basin_text(r, b, "OBSID1");
        const auto oi = outlet_at.find(id);
        if (oi != outlet_at.end())
        {
            b.has_outlet = true;
            b.outlet_longitude = oi->second.first;
            b.outlet_latitude = oi->second.second;
            b.outlet_in_bbox = b.outlet_longitude >= o.west && b.outlet_longitude <= o.east &&
                               b.outlet_latitude >= o.south && b.outlet_latitude <= o.north;
        }
        else
            missing_outlets.push_back(id);
        r.basins.push_back(std::move(b));
    }

    // Outlets outside the box (opt-in): one id query per 50 basins.
    if (o.include_outlets_outside_bbox && o.response_dir.empty() && !missing_outlets.empty())
    {
        for (std::size_t i = 0; i < missing_outlets.size(); i += 50)
        {
            const std::vector<std::string> ids(
                missing_outlets.begin() + static_cast<std::ptrdiff_t>(i),
                missing_outlets.begin() + static_cast<std::ptrdiff_t>(std::min(missing_outlets.size(), i + 50)));
            const std::string body = http.get(octopus_wfs_outlets_by_id_url(ids));
            check_body(body, "outlets by id");
            ++r.n_pages;
            for (const auto& f : parse_geojson_page(body, "outlets_by_id_" + std::to_string(i)).features)
            {
                const auto it = f.values.find("OBSID1");
                if (it == f.values.end() || !std::isfinite(f.x)) continue;
                outlet_at.emplace(value_text(it->second), std::make_pair(f.x, f.y));
            }
        }
        for (auto& b : r.basins)
        {
            if (b.has_outlet) continue;
            const auto oi = outlet_at.find(octopus_basin_text(r, b, "OBSID1"));
            if (oi == outlet_at.end()) continue;
            b.has_outlet = true;
            b.outlet_longitude = oi->second.first;
            b.outlet_latitude = oi->second.second;
            b.outlet_in_bbox = false;
        }
    }

    // Calibration rows and the studies.
    std::set<std::string> versions;
    std::map<std::string, lsdtt_octopus_study> studies;
    for (const auto& b : r.basins)
    {
        if (b.outlet_in_bbox) ++r.n_outlet_in_bbox;
        const std::string id = octopus_basin_text(r, b, "OBSID1");
        const std::string v = octopus_basin_text(r, b, "DBVERNO");
        if (!v.empty()) versions.insert(v);
        const std::string sid = octopus_basin_text(r, b, "STUDYID");
        auto& st = studies[sid];
        st.study_id = sid;
        st.author = octopus_basin_text(r, b, "AUTHOR");
        st.pub_year = octopus_basin_text(r, b, "PUBYEAR");
        st.ref_doi = octopus_basin_text(r, b, "REFDOI");
        st.ref_id = octopus_basin_text(r, b, "REFID");
        ++st.n_basins;

        if (!b.has_outlet || (!b.outlet_in_bbox && !o.include_outlets_outside_bbox))
        {
            r.no_outlet_in_bbox.push_back(id);
            continue;
        }
        lsdtt_octopus_csv_row row;
        row.rate = octopus_basin_number(r, b, r.rate_fields.rate);
        if (!std::isfinite(row.rate) || !(row.rate > 0.0))
        {
            r.dropped_no_rate.push_back(id);
            continue;
        }
        row.sample_id = id;
        row.obsid2 = octopus_basin_text(r, b, "OBSID2");
        row.study_id = sid;
        row.author = st.author;
        row.pub_year = st.pub_year;
        row.ref_doi = st.ref_doi;
        row.latitude = b.outlet_latitude;
        row.longitude = b.outlet_longitude;
        row.outlet_in_bbox = b.outlet_in_bbox;
        const double u = octopus_basin_number(r, b, r.rate_fields.uncertainty);
        row.uncertainty = std::isfinite(u) && u > 0.0 ? u : k_nan;
        if (!r.rate_fields.internal_uncertainty.empty())
        {
            const double iu = octopus_basin_number(r, b, r.rate_fields.internal_uncertainty);
            if (!r.rate_fields.internal_mass_flux_rate.empty())
            {
                const double g = octopus_basin_number(r, b, r.rate_fields.internal_mass_flux_rate);
                row.internal_uncertainty = std::isfinite(iu) && std::isfinite(g) && g > 0.0
                                               ? iu * row.rate / g : k_nan;
            }
            else
                row.internal_uncertainty = iu;
        }
        row.area_km2 = octopus_basin_number(r, b, "AREA");
        row.mean_elevation_m = octopus_basin_number(r, b, "ELEV_AVE");
        r.rows.push_back(row);
    }
    for (const auto& v : versions) r.db_version += (r.db_version.empty() ? "" : ";") + v;
    for (auto& [k, s] : studies) r.studies.push_back(s);
    return r;
}

std::string octopus_basin_text(const lsdtt_octopus_fetch_result& r, const lsdtt_octopus_basin& b,
                               const std::string& field)
{
    const int i = field_index(r, field);
    if (i < 0 || static_cast<std::size_t>(i) >= b.values.size()) return "";
    return value_text(b.values[static_cast<std::size_t>(i)]);
}

double octopus_basin_number(const lsdtt_octopus_fetch_result& r, const lsdtt_octopus_basin& b,
                            const std::string& field)
{
    const int i = field_index(r, field);
    if (i < 0 || static_cast<std::size_t>(i) >= b.values.size()) return k_nan;
    const lsdtt_octopus_value& v = b.values[static_cast<std::size_t>(i)];
    if (v.type == lsdtt_octopus_value::kind::integer || v.type == lsdtt_octopus_value::kind::real)
        return v.d;
    if (v.type == lsdtt_octopus_value::kind::text)
    {
        try
        {
            std::size_t pos = 0;
            const double d = std::stod(v.s, &pos);
            if (pos == v.s.size() && !is_octopus_no_data_sentinel(d)) return d;
        }
        catch (const std::exception&) {}
    }
    return k_nan;
}

// ---------------------------------------------------------------------------
// Outputs
// ---------------------------------------------------------------------------
std::string octopus_erosion_rate_csv_header()
{
    return "sample_id,latitude,longitude,erosion_rate,erosion_rate_uncertainty,erosion_rate_units,"
           "erosion_rate_uncertainty_kind,erosion_rate_internal_uncertainty,reported_area_km2,"
           "rate_scheme,rate_field,nuclide,obsid2,study_id,author,pub_year,ref_doi,"
           "mean_elevation_m,outlet_in_bbox";
}

void write_octopus_erosion_rate_csv(const std::string& path, const lsdtt_octopus_fetch_result& r)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << octopus_erosion_rate_csv_header() << "\n";
    for (const auto& w : r.rows)
        out << csv_text(w.sample_id) << "," << num(w.latitude, 12) << "," << num(w.longitude, 12)
            << "," << num(w.rate) << "," << num(w.uncertainty) << ",mm_per_kyr,"
            << r.rate_fields.kind << "," << num(w.internal_uncertainty, 6) << ","
            << num(w.area_km2) << "," << r.rate_scheme << "," << r.rate_fields.rate << ","
            << r.nuclide << "," << csv_text(w.obsid2) << "," << csv_text(w.study_id) << ","
            << csv_text(w.author) << "," << csv_text(w.pub_year) << "," << csv_text(w.ref_doi)
            << "," << num(w.mean_elevation_m) << "," << (w.outlet_in_bbox ? 1 : 0) << "\n";
}

void write_octopus_basin_polygons(const std::string& path, const std::string& driver, int epsg,
                                  const lsdtt_octopus_fetch_result& r)
{
    gdal_init();
    OGRSpatialReference wgs84, out_srs;
    wgs84.importFromEPSG(4326);
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (out_srs.importFromEPSG(epsg) != OGRERR_NONE)
        throw std::runtime_error("unknown EPSG:" + std::to_string(epsg));
    out_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    std::unique_ptr<OGRCoordinateTransformation> to_out(
        OGRCreateCoordinateTransformation(&wgs84, &out_srs));
    if (!to_out) throw std::runtime_error("no transform from EPSG:4326 to EPSG:" + std::to_string(epsg));

    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
    if (!drv) throw std::runtime_error("no OGR driver " + driver);
    std::filesystem::remove(path);
    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) throw std::runtime_error("could not create " + path);
    std::unique_ptr<GDALDataset, void (*)(GDALDataset*)> guard(ds, [](GDALDataset* d) { GDALClose(d); });
    OGRLayer* layer = ds->CreateLayer("octopus_basins", &out_srs, wkbMultiPolygon, nullptr);
    if (!layer) throw std::runtime_error("could not create layer 'octopus_basins' in " + path);
    for (const auto& fd : r.fields)
    {
        const OGRFieldType t = fd.type == lsdtt_octopus_value::kind::integer ? OFTInteger64
                               : fd.type == lsdtt_octopus_value::kind::real  ? OFTReal
                                                                             : OFTString;
        OGRFieldDefn d(fd.name.c_str(), t);
        layer->CreateField(&d);
    }
    {
        OGRFieldDefn a("outlet_in_bbox", OFTInteger);
        layer->CreateField(&a);
        OGRFieldDefn b("outlet_longitude", OFTReal);
        layer->CreateField(&b);
        OGRFieldDefn c("outlet_latitude", OFTReal);
        layer->CreateField(&c);
    }
    for (const auto& b : r.basins)
    {
        OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        std::unique_ptr<OGRFeature, void (*)(OGRFeature*)> fg(f, OGRFeature::DestroyFeature);
        for (std::size_t i = 0; i < r.fields.size() && i < b.values.size(); ++i)
        {
            const lsdtt_octopus_value& v = b.values[i];
            const char* name = r.fields[i].name.c_str();
            if (v.type == lsdtt_octopus_value::kind::null) { f->SetFieldNull(f->GetFieldIndex(name)); continue; }
            switch (r.fields[i].type)
            {
                case lsdtt_octopus_value::kind::integer: f->SetField(name, static_cast<GIntBig>(v.i)); break;
                case lsdtt_octopus_value::kind::real: f->SetField(name, v.d); break;
                default: f->SetField(name, value_text(v).c_str()); break;
            }
        }
        f->SetField("outlet_in_bbox", b.outlet_in_bbox ? 1 : 0);
        if (b.has_outlet)
        {
            f->SetField("outlet_longitude", b.outlet_longitude);
            f->SetField("outlet_latitude", b.outlet_latitude);
        }
        if (!b.wkb.empty())
        {
            OGRGeometry* g = nullptr;
            if (OGRGeometryFactory::createFromWkb(b.wkb.data(), nullptr, &g, b.wkb.size()) == OGRERR_NONE && g)
            {
                g->assignSpatialReference(&wgs84);
                std::unique_ptr<OGRGeometry> mp(OGRGeometryFactory::forceToMultiPolygon(g));
                if (mp && mp->transform(to_out.get()) == OGRERR_NONE) f->SetGeometryDirectly(mp.release());
            }
        }
        if (layer->CreateFeature(f) != OGRERR_NONE)
            throw std::runtime_error("could not write a feature to " + path);
    }
}

lsdtt_octopus_dem_bounds octopus_dem_bounds(const lsdtt_octopus_fetch_result& r, double pad_degrees,
                                            double max_area_km2)
{
    lsdtt_octopus_dem_bounds out;
    OGREnvelope env;
    for (const auto& b : r.basins)
    {
        if (!b.outlet_in_bbox || b.wkb.empty()) continue;
        const double area = octopus_basin_number(r, b, "AREA");
        if (std::isfinite(area) && area > max_area_km2) continue;
        OGRGeometry* g = nullptr;
        if (OGRGeometryFactory::createFromWkb(b.wkb.data(), nullptr, &g, b.wkb.size()) != OGRERR_NONE || !g)
            continue;
        std::unique_ptr<OGRGeometry> guard(g);
        OGREnvelope e;
        g->getEnvelope(&e);
        env.Merge(e);
        ++out.n_basins;
    }
    if (out.n_basins == 0) return out;
    out.ok = true;
    out.west = std::max(-180.0, env.MinX - pad_degrees);
    out.east = std::min(180.0, env.MaxX + pad_degrees);
    out.south = std::max(-90.0, env.MinY - pad_degrees);
    out.north = std::min(90.0, env.MaxY + pad_degrees);
    return out;
}

void write_octopus_dem_bounds(const std::string& path, const lsdtt_octopus_dem_bounds& bnd)
{
    gdal_init();
    if (!bnd.ok) throw std::invalid_argument("no DEM bounds to write");
    OGRSpatialReference wgs84;
    wgs84.importFromEPSG(4326);
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
    if (!drv) throw std::runtime_error("no OGR driver FlatGeobuf");
    std::filesystem::remove(path);
    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) throw std::runtime_error("could not create " + path);
    std::unique_ptr<GDALDataset, void (*)(GDALDataset*)> guard(ds, [](GDALDataset* d) { GDALClose(d); });
    OGRLayer* layer = ds->CreateLayer("dem_bounds", &wgs84, wkbPolygon, nullptr);
    if (!layer) throw std::runtime_error("could not create layer 'dem_bounds' in " + path);
    OGRFieldDefn n("n_basins", OFTInteger);
    layer->CreateField(&n);
    OGRLinearRing ring;
    ring.addPoint(bnd.west, bnd.south);
    ring.addPoint(bnd.east, bnd.south);
    ring.addPoint(bnd.east, bnd.north);
    ring.addPoint(bnd.west, bnd.north);
    ring.addPoint(bnd.west, bnd.south);
    OGRPolygon poly;
    poly.addRing(&ring);
    OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    std::unique_ptr<OGRFeature, void (*)(OGRFeature*)> fg(f, OGRFeature::DestroyFeature);
    f->SetField("n_basins", static_cast<int>(bnd.n_basins));
    f->SetGeometry(&poly);
    if (layer->CreateFeature(f) != OGRERR_NONE) throw std::runtime_error("could not write " + path);
}

std::string octopus_attribution_text(const lsdtt_octopus_fetch_result& r)
{
    std::ostringstream s;
    s << "Basin-averaged cosmogenic erosion rates from OCTOPUS (https://octopusdata.org),\n"
      << "provided by the University of Wollongong and licensed CC BY 4.0.\n\n"
      << "Source: " << r.source << "\n"
      << "Layers: " << k_basins_layer << ", " << k_outlets_layer << "\n"
      << "Accessed: " << r.access_date << " (UTC)\n"
      << "Database version (DBVERNO): " << (r.db_version.empty() ? "(none returned)" : r.db_version)
      << "\n"
      << "Bounding box (WGS84): west " << shortest(r.west) << ", south " << shortest(r.south)
      << ", east " << shortest(r.east) << ", north " << shortest(r.north) << "\n"
      << "Rate scheme: " << r.rate_scheme << " (" << r.nuclide << "; rate " << r.rate_fields.rate
      << ", uncertainty " << r.rate_fields.uncertainty << ", mm/kyr)\n"
      << "Transfer: the OCTOPUS WFS answers over HTTP only; the data are public and no credential\n"
      << "is sent, but their integrity in transit is not protected.\n\n"
      << "Cite OCTOPUS:\n"
      << "  Codilean, A. T., et al. (2018). Earth System Science Data, 10, 2123-2139.\n"
      << "  doi:10.5194/essd-10-2123-2018\n"
      << "  Codilean, A. T., et al. (2022). Earth System Science Data, 14, 3695-3713.\n"
      << "  doi:10.5194/essd-14-3695-2022\n"
      << "  Codilean, A. T. and Munack, H. (2025). Geochronology, 7, 113-122.\n"
      << "  doi:10.5194/gchron-7-113-2025\n"
      << "  Dataset (CRN Global v2.3): doi:10.71747/uow-r3gk326m.28216865.v1\n";
    if (r.rate_scheme == "st")
        s << "\nThe st rates (" << r.rate_fields.rate << ") were recalculated by OCTOPUS with CAIRN:\n"
          << "  Mudd, S. M., et al. (2016). Earth Surface Dynamics, 4, 655-674.\n"
          << "  doi:10.5194/esurf-4-655-2016\n";
    else if (r.rate_scheme == "lsdn" || r.rate_scheme == "lsdn_no_topo_shielding")
        s << "\nThe LSDn rates (" << r.rate_fields.rate << ") are a point evaluation (Balco v3 / "
          << "RIVERSAND style) at the basin's\nCENTR_LAT and ATM_PRESS, as published by OCTOPUS.\n";
    else
        s << "\nThe published rates (" << r.rate_fields.rate << ") are the source papers' own values, "
          << "with their own calculators\nand standards.\n";
    s << "\nOriginal studies (STUDYID, AUTHOR, PUBYEAR, REFDOI, n samples):\n";
    for (const auto& st : r.studies)
        s << "  " << st.study_id << ", " << st.author << ", " << st.pub_year << ", "
          << (st.ref_doi.empty() ? "(no DOI; REFID " + st.ref_id + ")" : st.ref_doi) << ", "
          << st.n_basins << "\n";
    s << "\nCite the original studies above and OCTOPUS when using these rates.\n";
    return s.str();
}

} // namespace lsdtt
