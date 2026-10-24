// lsdtt_rotation_frame_file.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_rotation_frame_file.hpp for the format.

#include "lsdtt_rotation_frame_file.hpp"

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "lsdtt_crs_suggest.hpp"
#include "lsdtt_data_loader.hpp"
#include "lsdtt_param_parser.hpp"
#include "lsdtt_raster_warp.hpp"

namespace lsdtt {

namespace {

// -----------------------------------------------------------------------------
// Writing
// -----------------------------------------------------------------------------

// 17 significant digits, always with a decimal point before any exponent.
std::string format_number(double v)
{
    std::ostringstream os;
    os.imbue(std::locale::classic());
    os << std::setprecision(17) << v;
    std::string s = os.str();
    if (s.find_first_of(".ni") != std::string::npos) return s;   // has a point, or inf/nan
    const std::size_t e = s.find('e');
    if (e == std::string::npos) return s + ".0";
    return s.substr(0, e) + ".0" + s.substr(e);
}

// One-line text: newlines and the indentation after them are dropped (WKT
// whitespace between tokens is not significant).
std::string one_line(const std::string& text)
{
    std::string out;
    out.reserve(text.size());
    bool skip_indent = false;
    for (char c : text)
    {
        if (c == '\n' || c == '\r') { skip_indent = true; continue; }
        if (skip_indent && (c == ' ' || c == '\t')) continue;
        skip_indent = false;
        out += c;
    }
    return out;
}

std::string quoted(const std::string& text)
{
    std::string out = "'";
    for (char c : one_line(text))
    {
        if (c == '\'') out += "''";
        else out += c;
    }
    return out + "'";
}

std::string definition_token(oblique_frame_definition d)
{
    switch (d)
    {
    case oblique_frame_definition::two_points: return "two_points";
    case oblique_frame_definition::line:       return "line";
    case oblique_frame_definition::polygon:    return "polygon";
    case oblique_frame_definition::crs:        return "crs";
    }
    return "line";
}

std::string now_utc()
{
    const std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

// -----------------------------------------------------------------------------
// Reading
// -----------------------------------------------------------------------------

std::string trim(const std::string& s)
{
    const std::size_t a = s.find_first_not_of(" \t");
    if (a == std::string::npos) return "";
    const std::size_t b = s.find_last_not_of(" \t");
    return s.substr(a, b - a + 1);
}

// True if what follows a quoted value is blank or a comment.
bool blank_or_comment(const std::string& rest)
{
    const std::string t = trim(rest);
    return t.empty() || t[0] == '#';
}

// Parse the text after "key:". Returns false with err set on a malformed value.
bool parse_value(const std::string& raw, std::string& value, std::string& err)
{
    const std::string v = trim(raw);
    value.clear();
    if (v.empty()) return true;
    if (v[0] == '\'')
    {
        std::size_t i = 1;
        for (; i < v.size(); ++i)
        {
            if (v[i] != '\'') { value += v[i]; continue; }
            if (i + 1 < v.size() && v[i + 1] == '\'') { value += '\''; ++i; continue; }
            break;
        }
        if (i >= v.size()) { err = "unterminated single-quoted value"; return false; }
        if (!blank_or_comment(v.substr(i + 1))) { err = "text after the closing quote"; return false; }
        return true;
    }
    if (v[0] == '"')
    {
        std::size_t i = 1;
        for (; i < v.size(); ++i)
        {
            if (v[i] == '"') break;
            if (v[i] == '\\')
            {
                if (i + 1 >= v.size()) break;
                const char e = v[++i];
                if (e == '"' || e == '\\' || e == '/') value += e;
                else { err = std::string("unsupported escape \\") + e; return false; }
                continue;
            }
            value += v[i];
        }
        if (i >= v.size()) { err = "unterminated double-quoted value"; return false; }
        if (!blank_or_comment(v.substr(i + 1))) { err = "text after the closing quote"; return false; }
        return true;
    }
    // Plain: a YAML comment starts at " #".
    const std::size_t hash = v.find(" #");
    value = trim(hash == std::string::npos ? v : v.substr(0, hash));
    return true;
}

class key_reader
{
public:
    explicit key_reader(const std::map<std::string, std::string>& kv) : kv_(kv) {}

    const std::string* find(const std::string& key)
    {
        used_.insert(key);
        const auto it = kv_.find(key);
        return it == kv_.end() ? nullptr : &it->second;
    }

    std::string text(const std::string& key)
    {
        const std::string* v = find(key);
        return v ? *v : std::string();
    }

    std::string required_text(const std::string& key)
    {
        const std::string* v = find(key);
        if (!v || v->empty()) missing_.push_back(key);
        return v ? *v : std::string();
    }

    double number(const std::string& key)
    {
        const std::string* v = find(key);
        if (!v || v->empty()) { missing_.push_back(key); return 0.0; }
        std::istringstream is(*v);
        is.imbue(std::locale::classic());
        double d = 0.0;
        is >> d;
        if (is.fail() || !(is >> std::ws).eof() || !std::isfinite(d))
            bad_.push_back(key + " ('" + *v + "' is not a finite number)");
        return d;
    }

    int integer(const std::string& key)
    {
        const std::string* v = find(key);
        if (!v || v->empty()) { missing_.push_back(key); return 0; }
        std::istringstream is(*v);
        is.imbue(std::locale::classic());
        long long i = 0;
        is >> i;
        if (is.fail() || !(is >> std::ws).eof() || i < 0 || i > 2147483647LL)
            bad_.push_back(key + " ('" + *v + "' is not a non-negative integer)");
        return static_cast<int>(i);
    }

    bool boolean(const std::string& key)
    {
        const std::string* v = find(key);
        if (!v || v->empty()) { missing_.push_back(key); return false; }
        if (*v == "true") return true;
        if (*v == "false") return false;
        bad_.push_back(key + " ('" + *v + "' is not true or false)");
        return false;
    }

    std::vector<std::string> unknown() const
    {
        std::vector<std::string> out;
        for (const auto& kv : kv_)
            if (!used_.count(kv.first)) out.push_back(kv.first);
        return out;
    }

    std::vector<std::string> missing_;
    std::vector<std::string> bad_;

private:
    const std::map<std::string, std::string>& kv_;
    std::set<std::string> used_;
};

std::string join(const std::vector<std::string>& items)
{
    std::string out;
    for (const auto& s : items) out += (out.empty() ? "" : ", ") + s;
    return out;
}

double snap_out(double half_extent, double cell_size)
{
    return std::ceil(half_extent / cell_size - 1e-9) * cell_size;
}

} // namespace

// -----------------------------------------------------------------------------
// write_rotation_frame
// -----------------------------------------------------------------------------
void write_rotation_frame(const std::string& path, const lsdtt_oblique_frame& f,
                          const lsdtt_rotation_target& target,
                          const lsdtt_rotation_provenance& provenance)
{
    if (target.epsg < 0 || (target.epsg == 0 && target.crs_wkt.empty()))
        throw std::invalid_argument(
            "write_rotation_frame: the rotate-back target needs an EPSG code or a CRS WKT");
    if (!std::isfinite(target.cell_size) || target.cell_size <= 0.0)
        throw std::invalid_argument("write_rotation_frame: target cell_size must be > 0");

    const std::string frame_wkt = f.to_wkt();
    const std::string proj = f.to_proj();

    std::ostringstream os;
    os << "# lsdtt3 rotation frame file. key: value lines; a flat YAML mapping.\n"
       << "# Keep each value on one line. See the lsdtt3 docs on rotated frames.\n"
       << "# frame_proj: " << one_line(proj) << "\n"
       << "\n# Provenance\n"
       << "frame_file_version: " << rotation_frame_file_version << "\n"
       << "lsdtt3_version: "
       << quoted(provenance.lsdtt3_version.empty() ? lsdtt3_package_version()
                                                   : provenance.lsdtt3_version) << "\n"
       << "created_utc: "
       << quoted(provenance.created_utc.empty() ? now_utc() : provenance.created_utc) << "\n"
       << "created_by: " << quoted(provenance.created_by) << "\n"
       << "rotated_dem_fname: " << quoted(provenance.rotated_dem_fname) << "\n"
       << "\n# Frame\n"
       << "centre_lon: " << format_number(f.centre.longitude) << "\n"
       << "centre_lat: " << format_number(f.centre.latitude) << "\n"
       << "azimuth_deg: " << format_number(f.azimuth_degrees) << "\n"
       << "feature_axis: x\n"
       << "cell_size: " << format_number(f.cell_size) << "\n"
       << "frame_x_min: " << format_number(f.x_min) << "\n"
       << "frame_x_max: " << format_number(f.x_max) << "\n"
       << "frame_y_min: " << format_number(f.y_min) << "\n"
       << "frame_y_max: " << format_number(f.y_max) << "\n"
       << "base_crs_wkt: " << quoted(f.base_geographic_wkt) << "\n"
       << "frame_crs_wkt: " << quoted(frame_wkt) << "\n"
       << "\n# Definition\n"
       << "frame_source: " << definition_token(f.definition) << "\n"
       << "feature_fname: " << quoted(provenance.feature_fname) << "\n"
       << "feature_layer: " << quoted(provenance.feature_layer) << "\n"
       << "end1_lon: " << format_number(f.end1.longitude) << "\n"
       << "end1_lat: " << format_number(f.end1.latitude) << "\n"
       << "end2_lon: " << format_number(f.end2.longitude) << "\n"
       << "end2_lat: " << format_number(f.end2.latitude) << "\n"
       << "axis_length_m: " << format_number(f.axis_length_metres) << "\n"
       << "width_m: " << format_number(f.width_metres) << "\n"
       << "end_buffer_m: " << format_number(f.end_buffer_metres) << "\n"
       << "buffer_all_sides: " << (f.buffer_all_sides ? "true" : "false") << "\n"
       << "\n# Rotate-back target\n"
       << "target_crs_epsg: " << target.epsg << "\n"
       << "target_crs_wkt: " << quoted(target.epsg == 0 ? target.crs_wkt : std::string()) << "\n"
       << "target_cell_size: " << format_number(target.cell_size) << "\n"
       << "target_origin_x: " << format_number(target.origin_x) << "\n"
       << "target_origin_y: " << format_number(target.origin_y) << "\n"
       << "source_dem_fname: " << quoted(target.source_dem_fname) << "\n";

    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("write_rotation_frame: cannot open " + path);
    out << os.str();
    out.close();
    if (!out) throw std::runtime_error("write_rotation_frame: failed writing " + path);
}

// -----------------------------------------------------------------------------
// read_rotation_frame
// -----------------------------------------------------------------------------
lsdtt_rotation_frame_file read_rotation_frame(const std::string& path)
{
    lsdtt_rotation_frame_file r;
    const std::string where = "rotation frame file " + path + ": ";

    std::ifstream in(path, std::ios::binary);
    if (!in)
    {
        r.message = where + "cannot be opened";
        return r;
    }

    std::map<std::string, std::string> kv;
    std::string line;
    int line_no = 0;
    while (std::getline(in, line))
    {
        ++line_no;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line_no == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;
        const std::size_t colon = t.find(':');
        if (colon == std::string::npos || colon == 0)
        {
            r.message = where + "line " + std::to_string(line_no) + " is not 'key: value'";
            return r;
        }
        const std::string key = trim(t.substr(0, colon));
        std::string value, err;
        if (!parse_value(t.substr(colon + 1), value, err))
        {
            r.message = where + "line " + std::to_string(line_no) + " (" + key + "): " + err;
            return r;
        }
        if (!kv.emplace(key, value).second)
        {
            r.message = where + "key '" + key + "' appears more than once";
            return r;
        }
    }

    key_reader k(kv);

    // Version first: a newer file may use keys this build does not know.
    if (!k.find("frame_file_version"))
    {
        r.message = where + "no frame_file_version; this is not a rotation frame file";
        return r;
    }
    r.version = k.integer("frame_file_version");
    if (!k.bad_.empty())
    {
        r.message = where + "bad value: " + join(k.bad_);
        return r;
    }
    if (r.version < 1 || r.version > rotation_frame_file_version)
    {
        r.message = where + "frame_file_version " + std::to_string(r.version) +
                    " is not supported (this build reads 1 to " +
                    std::to_string(rotation_frame_file_version) + "); use a newer lsdtt3";
        return r;
    }

    lsdtt_rotation_provenance& p = r.provenance;
    p.lsdtt3_version = k.text("lsdtt3_version");
    p.created_utc = k.text("created_utc");
    p.created_by = k.text("created_by");
    p.rotated_dem_fname = k.text("rotated_dem_fname");
    p.feature_fname = k.text("feature_fname");
    p.feature_layer = k.text("feature_layer");

    lsdtt_oblique_frame& f = r.frame;
    f.centre.longitude = k.number("centre_lon");
    f.centre.latitude = k.number("centre_lat");
    f.azimuth_degrees = k.number("azimuth_deg");
    const std::string axis = k.required_text("feature_axis");
    f.cell_size = k.number("cell_size");
    f.x_min = k.number("frame_x_min");
    f.x_max = k.number("frame_x_max");
    f.y_min = k.number("frame_y_min");
    f.y_max = k.number("frame_y_max");
    f.base_geographic_wkt = k.required_text("base_crs_wkt");
    const std::string frame_wkt = k.required_text("frame_crs_wkt");
    const std::string source = k.required_text("frame_source");
    f.end1.longitude = k.number("end1_lon");
    f.end1.latitude = k.number("end1_lat");
    f.end2.longitude = k.number("end2_lon");
    f.end2.latitude = k.number("end2_lat");
    f.axis_length_metres = k.number("axis_length_m");
    f.width_metres = k.number("width_m");
    f.end_buffer_metres = k.number("end_buffer_m");
    f.buffer_all_sides = k.boolean("buffer_all_sides");

    lsdtt_rotation_target& t = r.target;
    t.epsg = k.integer("target_crs_epsg");
    t.crs_wkt = k.text("target_crs_wkt");
    t.cell_size = k.number("target_cell_size");
    t.origin_x = k.number("target_origin_x");
    t.origin_y = k.number("target_origin_y");
    t.source_dem_fname = k.text("source_dem_fname");

    if (!k.missing_.empty())
    {
        r.message = where + "missing: " + join(k.missing_);
        return r;
    }
    if (!k.bad_.empty())
    {
        r.message = where + "bad value: " + join(k.bad_);
        return r;
    }
    for (const auto& key : k.unknown())
        r.warnings.push_back(where + "unknown key '" + key + "' ignored");

    if (axis != "x")
    {
        r.message = where + "feature_axis is '" + axis + "'; only x is supported";
        return r;
    }
    if (source == "two_points") f.definition = oblique_frame_definition::two_points;
    else if (source == "line") f.definition = oblique_frame_definition::line;
    else if (source == "polygon") f.definition = oblique_frame_definition::polygon;
    else if (source == "crs") f.definition = oblique_frame_definition::crs;
    else
    {
        r.message = where + "frame_source '" + source +
                    "' is not one of two_points, line, polygon, crs";
        return r;
    }
    if (t.epsg == 0 && t.crs_wkt.empty())
    {
        r.message = where + "target_crs_epsg is 0 and target_crs_wkt is empty; no rotate-back CRS";
        return r;
    }
    if (t.epsg != 0 && !t.crs_wkt.empty())
        r.warnings.push_back(where + "target_crs_epsg is set, so target_crs_wkt is ignored");
    if (f.cell_size <= 0.0 || t.cell_size <= 0.0)
    {
        r.message = where + "cell_size and target_cell_size must be > 0";
        return r;
    }
    if (!(f.x_max > f.x_min && f.y_max > f.y_min))
    {
        r.message = where + "the frame extent is empty";
        return r;
    }

    // The stored frame made the raster, so it is what is used; these checks
    // only confirm it is self-consistent.
    try
    {
        if (fold_axis_azimuth(f.azimuth_degrees) != f.azimuth_degrees)
            throw std::runtime_error("azimuth_deg is not folded into (-90, 90]");
        if (!crs_equal_semantically(frame_wkt, f.to_wkt()))
            throw std::runtime_error(
                "frame_crs_wkt does not match the frame built from centre_lon, centre_lat, "
                "azimuth_deg and base_crs_wkt");
        check_oblique_frame_axis(f);
    }
    catch (const std::exception& e)
    {
        r.message = where + e.what();
        return r;
    }

    if (f.definition != oblique_frame_definition::crs)
    {
        try
        {
            const lsdtt_oblique_frame g =
                frame_from_points(f.end1, f.end2, f.width_metres, f.end_buffer_metres,
                                  f.cell_size, f.base_geographic_wkt);
            const double y_max = snap_out(
                0.5 * f.width_metres + (f.buffer_all_sides ? f.end_buffer_metres : 0.0),
                f.cell_size);
            const double tol_deg = 1e-9;
            const double tol_m = 1e-6;
            if (std::abs(g.centre.longitude - f.centre.longitude) > tol_deg ||
                std::abs(g.centre.latitude - f.centre.latitude) > tol_deg ||
                std::abs(g.azimuth_degrees - f.azimuth_degrees) > tol_deg ||
                std::abs(g.axis_length_metres - f.axis_length_metres) > tol_m ||
                std::abs(g.x_min - f.x_min) > tol_m || std::abs(g.x_max - f.x_max) > tol_m ||
                std::abs(-y_max - f.y_min) > tol_m || std::abs(y_max - f.y_max) > tol_m)
            {
                std::ostringstream msg;
                msg << std::setprecision(12) << where
                    << "the frame rebuilt from the end points, width, buffer and cell size "
                       "differs from the stored one (centre " << g.centre.longitude << ", "
                    << g.centre.latitude << " vs " << f.centre.longitude << ", "
                    << f.centre.latitude << "; azimuth " << g.azimuth_degrees << " vs "
                    << f.azimuth_degrees << "; x_max " << g.x_max << " vs " << f.x_max
                    << "; y_max " << y_max << " vs " << f.y_max
                    << "); the stored frame is used (a hand edit or a different PROJ?)";
                r.warnings.push_back(msg.str());
            }
        }
        catch (const std::exception& e)
        {
            r.warnings.push_back(where + "the frame could not be rebuilt from its definition (" +
                                 e.what() + "); the stored frame is used");
        }
    }

    std::ostringstream msg;
    msg << std::setprecision(10) << "rotation frame: " << source << ", centre ("
        << f.centre.longitude << ", " << f.centre.latitude << "), azimuth "
        << f.azimuth_degrees << " deg, " << f.n_cols() << " x " << f.n_rows() << " cells of "
        << f.cell_size << " m; rotate back to "
        << (t.epsg != 0 ? "EPSG:" + std::to_string(t.epsg) : std::string("the WKT target"))
        << " at " << t.cell_size << " m";
    r.message = msg.str();
    r.ok = true;
    return r;
}

// -----------------------------------------------------------------------------
// analysis_epsg_for_frame_centre
// -----------------------------------------------------------------------------
int analysis_epsg_for_frame_centre(const lsdtt_oblique_frame& frame)
{
    // A one-cell geographic raster at the centre, so lsdtt_suggest_crs applies
    // its usual rule.
    lsdtt_raster_meta centre;
    centre.crs_wkt = frame.base_geographic_wkt.empty() ? base_geographic_wkt_of("")
                                                       : frame.base_geographic_wkt;
    centre.cell_size = 1e-6;
    centre.n_cols = 1;
    centre.n_rows = 1;
    centre.x_min = frame.centre.longitude - 0.5e-6;
    centre.y_min = frame.centre.latitude - 0.5e-6;
    const lsdtt_crs_suggestion sug = lsdtt_suggest_crs(centre, false, true);
    if (sug.recommended_epsg <= 0)
        throw std::runtime_error("no analysis CRS for the frame centre");
    return sug.recommended_epsg;
}

// -----------------------------------------------------------------------------
// rotation_frame_from_raster
// -----------------------------------------------------------------------------
lsdtt_rotation_frame_file rotation_frame_from_raster(const std::string& raster_path)
{
    lsdtt_rotation_frame_file r;
    const std::string where = "frame from raster " + raster_path + ": ";
    try
    {
        lsdtt_data_loader loader;
        loader.open_raster(raster_path);
        const lsdtt_raster_meta meta = loader.get_meta();
        r.frame = frame_from_crs_wkt(meta.crs_wkt, meta.cell_size, meta.x_min, meta.x_max(),
                                     meta.y_min, meta.y_max());
        r.warnings = r.frame.warnings;

        r.target.epsg = analysis_epsg_for_frame_centre(r.frame);
        r.target.cell_size = meta.cell_size;
        r.target.origin_x = 0.0;
        r.target.origin_y = 0.0;
        r.provenance.rotated_dem_fname = raster_path;
        r.version = rotation_frame_file_version;
        r.warnings.push_back(
            where + "no frame file; rotating back to EPSG:" + std::to_string(r.target.epsg) +
            " on the lattice anchored at (0, 0). If the source DEM's origin was not on that "
            "lattice, its exact pixels cannot be recovered without the frame file");
    }
    catch (const std::exception& e)
    {
        r.message = where + e.what();
        return r;
    }
    std::ostringstream msg;
    msg << std::setprecision(10) << "rotation frame rebuilt from " << raster_path
        << ": centre (" << r.frame.centre.longitude << ", " << r.frame.centre.latitude
        << "), azimuth " << r.frame.azimuth_degrees << " deg; rotate back to EPSG:"
        << r.target.epsg;
    r.message = msg.str();
    r.ok = true;
    return r;
}

} // namespace lsdtt
