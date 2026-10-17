// lsdtt_lithology_map_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_lithology_map_cpu.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

#include "../core/lsdtt_csv_table.hpp"

namespace lsdtt
{

namespace
{

constexpr int k_max_code = 65534;   // 65535 is reserved as nodata in uint16 code rasters

std::string lower(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string at(const std::string& path, int line)
{
    return path + " line " + std::to_string(line);
}

bool parse_int(const std::string& s, long& out)
{
    if (s.empty()) return false;
    std::size_t pos = 0;
    try { out = std::stol(s, &pos); }
    catch (const std::exception&) { return false; }
    return pos == s.size();
}

bool parse_double(const std::string& s, double& out)
{
    if (s.empty()) return false;
    std::size_t pos = 0;
    try { out = std::stod(s, &pos); }
    catch (const std::exception&) { return false; }
    return pos == s.size();
}

int checked_code(const std::string& s, const std::string& where)
{
    long v = 0;
    if (!parse_int(s, v) || v < 0 || v > k_max_code)
        throw std::invalid_argument(where + ": lithology code '" + s +
                                    "' must be an integer in 0.." + std::to_string(k_max_code));
    return static_cast<int>(v);
}

std::string num(double v)
{
    std::ostringstream s;
    s << std::setprecision(7) << v;
    return s.str();
}

} // namespace

lsdtt_lithology_property_table load_lithology_property_table(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::invalid_argument("cannot open lithology property table: " + path);

    lsdtt_lithology_property_table t;
    t.path = path;
    std::map<std::string, std::string> header;
    std::string line;
    int lineno = 0;
    std::vector<std::string> columns;
    int header_line = 0;
    while (std::getline(in, line))
    {
        ++lineno;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string trimmed = csv_trim(line);
        if (trimmed.empty()) continue;
        if (trimmed[0] == '#')
        {
            const std::string body = csv_trim(trimmed.substr(1));
            const auto colon = body.find(':');
            if (colon != std::string::npos)
                header[lower(csv_trim(body.substr(0, colon)))] = csv_trim(body.substr(colon + 1));
            continue;
        }
        columns = csv_split_row_trimmed(trimmed);
        header_line = lineno;
        break;
    }
    if (columns.empty()) throw std::invalid_argument(path + ": no column header row");

    for (const char* key : {"erodibility_convention", "m_exponent", "n_exponent"})
        if (!header.count(key))
            throw std::invalid_argument(path + ": required header line '# " + std::string(key) +
                                        ": ...' is missing (K is only meaningful with its area "
                                        "convention, m and n)");
    try
    {
        t.calibration.area_convention =
            parse_erodibility_area_convention(header["erodibility_convention"]);
    }
    catch (const std::exception& e)
    {
        throw std::invalid_argument(path + ": erodibility_convention: " + e.what());
    }
    double v = 0.0;
    if (!parse_double(header["m_exponent"], v))
        throw std::invalid_argument(path + ": m_exponent '" + header["m_exponent"] + "' is not a number");
    t.calibration.m_exponent = static_cast<float>(v);
    if (!parse_double(header["n_exponent"], v))
        throw std::invalid_argument(path + ": n_exponent '" + header["n_exponent"] + "' is not a number");
    t.calibration.n_exponent = static_cast<float>(v);
    if (header.count("calibration_cell_size_metres"))
    {
        if (!parse_double(header["calibration_cell_size_metres"], v) || !(v > 0.0))
            throw std::invalid_argument(path + ": calibration_cell_size_metres must be a number > 0");
        t.calibration.cell_size = v;
    }
    if (header.count("provenance")) t.provenance = header["provenance"];

    auto col = [&](const std::string& name, bool required) {
        for (std::size_t i = 0; i < columns.size(); ++i)
            if (lower(columns[i]) == name) return static_cast<int>(i);
        if (required)
            throw std::invalid_argument(at(path, header_line) + ": required column '" + name +
                                        "' is missing (need code, name, erodibility)");
        return -1;
    };
    const int c_code = col("code", true), c_name = col("name", true), c_k = col("erodibility", true);
    const int c_ref = col("reference", false);
    std::vector<std::string> ignored, hillslope;
    for (std::size_t i = 0; i < columns.size(); ++i)
    {
        const int ii = static_cast<int>(i);
        if (ii == c_code || ii == c_name || ii == c_k || ii == c_ref) continue;
        const std::string n = lower(columns[i]);
        if (n == "critical_hillslope_slope" || n == "hillslope_diffusivity" || n == "sc" ||
            n == "diffusivity")
            hillslope.push_back(columns[i]);
        else
            ignored.push_back(columns[i]);
    }
    for (const auto& h : hillslope)
        t.warnings.push_back(path + ": column '" + h + "' is not applied in this version; "
                             "hillslope parameters stay uniform");
    if (!ignored.empty())
    {
        std::string list;
        for (const auto& c : ignored) list += (list.empty() ? "" : ", ") + c;
        t.notes.push_back(path + ": ignored columns: " + list);
    }

    std::map<std::string, int> name_line;
    while (std::getline(in, line))
    {
        ++lineno;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string trimmed = csv_trim(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;
        const std::vector<std::string> f = csv_split_row_trimmed(trimmed);
        const auto field = [&](int c) { return c < static_cast<int>(f.size()) ? f[c] : std::string(); };
        lsdtt_lithology_unit u;
        u.line = lineno;
        u.code = checked_code(field(c_code), at(path, lineno));
        u.name = field(c_name);
        if (u.name.empty())
            throw std::invalid_argument(at(path, lineno) + ": code " + std::to_string(u.code) +
                                        " has no name");
        double k = 0.0;
        if (!parse_double(field(c_k), k))
            throw std::invalid_argument(at(path, lineno) + ", column erodibility: '" + field(c_k) +
                                        "' is not a number");
        u.erodibility = static_cast<float>(k);
        if (!(u.erodibility > 0.0f) || !std::isfinite(u.erodibility))
            throw std::invalid_argument(at(path, lineno) + ": K must be > 0; code " +
                                        std::to_string(u.code) + " (" + u.name + ") has K = " +
                                        field(c_k));
        if (c_ref >= 0) u.reference = field(c_ref);
        const auto dup = t.units.find(u.code);
        if (dup != t.units.end())
            throw std::invalid_argument(path + ": duplicate code " + std::to_string(u.code) +
                                        " on lines " + std::to_string(dup->second.line) + " and " +
                                        std::to_string(lineno));
        const auto dn = name_line.find(lower(u.name));
        if (dn != name_line.end())
            throw std::invalid_argument(path + ": duplicate name '" + u.name + "' on lines " +
                                        std::to_string(dn->second) + " and " + std::to_string(lineno));
        name_line[lower(u.name)] = lineno;
        t.units[u.code] = u;
    }
    if (t.units.empty()) throw std::invalid_argument(path + ": the table has no lithology rows");
    return t;
}

std::map<std::string, int> load_lithology_crosswalk(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::invalid_argument("cannot open lithology crosswalk: " + path);
    std::string line;
    int lineno = 0, c_value = -1, c_code = -1;
    std::map<std::string, int> out;
    std::map<std::string, int> seen;
    while (std::getline(in, line))
    {
        ++lineno;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string trimmed = csv_trim(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;
        const std::vector<std::string> f = csv_split_row_trimmed(trimmed);
        if (c_value < 0)
        {
            for (std::size_t i = 0; i < f.size(); ++i)
            {
                if (lower(f[i]) == "value") c_value = static_cast<int>(i);
                if (lower(f[i]) == "code") c_code = static_cast<int>(i);
            }
            if (c_value < 0 || c_code < 0)
                throw std::invalid_argument(at(path, lineno) + ": the crosswalk needs columns "
                                            "'value' and 'code'");
            continue;
        }
        const std::string value = c_value < static_cast<int>(f.size()) ? f[c_value] : "";
        const std::string code = c_code < static_cast<int>(f.size()) ? f[c_code] : "";
        if (seen.count(value))
            throw std::invalid_argument(path + ": duplicate value '" + value + "' on lines " +
                                        std::to_string(seen[value]) + " and " + std::to_string(lineno));
        seen[value] = lineno;
        out[value] = checked_code(code, at(path, lineno));
    }
    if (c_value < 0) throw std::invalid_argument(path + ": the crosswalk is empty");
    return out;
}

std::vector<int> lithology_codes_for_values(
    const std::vector<std::string>&       values,
    const std::vector<std::size_t>&       value_pixel_counts,
    const std::map<std::string, int>*     crosswalk,
    const lsdtt_lithology_property_table& table)
{
    std::map<std::string, int> by_name;
    for (const auto& [code, u] : table.units) by_name[lower(u.name)] = code;

    std::vector<int> codes(values.size(), -1);
    std::vector<std::size_t> missing;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        const std::string& v = values[i];
        if (crosswalk)
        {
            const auto it = crosswalk->find(v);
            if (it != crosswalk->end()) codes[i] = it->second;
        }
        else
        {
            long c = 0;
            if (parse_int(v, c) && c >= 0 && c <= k_max_code) codes[i] = static_cast<int>(c);
            else
            {
                const auto it = by_name.find(lower(csv_trim(v)));
                if (it != by_name.end()) codes[i] = it->second;
            }
        }
        if (codes[i] < 0) missing.push_back(i);
    }
    if (!missing.empty())
    {
        std::ostringstream msg;
        msg << missing.size() << " lithology value" << (missing.size() > 1 ? "s" : "")
            << " could not be given a code "
            << (crosswalk ? "(not in the crosswalk)"
                          : "(neither an integer code nor a name in the property table)")
            << ":";
        for (std::size_t k = 0; k < missing.size() && k < 25; ++k)
        {
            const std::size_t i = missing[k];
            msg << "\n  '" << values[i] << "'";
            if (i < value_pixel_counts.size()) msg << " (" << value_pixel_counts[i] << " pixels)";
        }
        if (missing.size() > 25) msg << "\n  ... and " << missing.size() - 25 << " more";
        msg << "\nAdd them to the property table by name, or map them with "
               "lithology_crosswalk_fname (CSV value,code).";
        throw std::invalid_argument(msg.str());
    }
    return codes;
}

lsdtt_lithology_erodibility_result make_lithology_erodibility(
    const lsdtt_raster_cpu&                    grid,
    const std::vector<std::int32_t>&           value_index,
    const std::vector<std::string>&            values,
    const std::vector<int>&                    codes_for_values,
    const lsdtt_lithology_property_table&      table,
    const lsdtt_lithology_erodibility_options& opts)
{
    const lsdtt_raster_meta& meta = grid.get_meta();
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    if (value_index.size() != n)
        throw std::invalid_argument("make_lithology_erodibility: value_index has " +
                                    std::to_string(value_index.size()) + " entries, the grid " +
                                    std::to_string(n));
    if (codes_for_values.size() != values.size())
        throw std::invalid_argument("make_lithology_erodibility: one code per value is needed");
    if (opts.unmapped_code >= 0 && !table.units.count(opts.unmapped_code))
        throw std::invalid_argument("lithology_unmapped_code " + std::to_string(opts.unmapped_code) +
                                    " is not in the property table");

    const std::vector<std::uint8_t>& domain = grid.get_mask();
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> codes(n, nd), k(n, nd);
    std::vector<std::uint8_t> mask(n, 0);
    std::map<int, lsdtt_lithology_legend_row> legend;
    lsdtt_lithology_erodibility_result out;
    std::size_t first_unmapped = n;

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!domain[i]) continue;
        int code = -1;
        if (value_index[i] >= 0)
            code = codes_for_values[static_cast<std::size_t>(value_index[i])];
        else
        {
            ++out.n_unmapped;
            if (first_unmapped == n) first_unmapped = i;
            if (opts.unmapped_code >= 0) code = opts.unmapped_code;
            else if (!opts.unmapped_inside_domain) continue;   // footprint: stays nodata
            else continue;                                     // refused below
        }
        const auto u = table.units.find(code);
        if (u == table.units.end())
        {
            std::string where = values.empty() || value_index[i] < 0
                                    ? std::string("")
                                    : " (map value '" + values[static_cast<std::size_t>(value_index[i])] + "')";
            throw std::invalid_argument("lithology code " + std::to_string(code) + where +
                                        " is not in the property table " + table.path);
        }
        codes[i] = static_cast<float>(code);
        k[i] = u->second.erodibility;
        mask[i] = 1;
        auto& row = legend[code];
        row.code = code;
        row.name = u->second.name;
        row.erodibility = u->second.erodibility;
        ++row.n_pixels;
    }

    if (out.n_unmapped > 0 && opts.unmapped_code < 0 && opts.unmapped_inside_domain)
    {
        const std::size_t c = static_cast<std::size_t>(meta.n_cols);
        throw std::invalid_argument(
            std::to_string(out.n_unmapped) + " pixels inside the model domain have no lithology; "
            "first at (" + std::to_string(first_unmapped / c) + ", " +
            std::to_string(first_unmapped % c) + "). Name a table unit for them with "
            "lithology_unmapped_code, or take the grid from the map so its footprint is the domain.");
    }
    if (out.n_unmapped > 0)
        out.notes.push_back(std::to_string(out.n_unmapped) + " pixels had no lithology: " +
                            (opts.unmapped_code >= 0
                                 ? "given code " + std::to_string(opts.unmapped_code)
                                 : std::string("left as nodata (outside the map's footprint)")));

    for (std::size_t v = 0; v < values.size(); ++v)
    {
        auto it = legend.find(codes_for_values[v]);
        if (it != legend.end()) it->second.map_values.push_back(values[v]);
    }
    for (auto& [code, row] : legend) out.legend.push_back(row);
    out.code_raster = lsdtt_raster_cpu(std::move(codes), std::vector<std::uint8_t>(mask), meta);
    out.erodibility_raster = lsdtt_raster_cpu(std::move(k), std::move(mask), meta);
    return out;
}

std::string lithology_erodibility_lem_parameters(const lsdtt_lithology_property_table& table,
                                                 const std::string& erodibility_fname)
{
    const auto& c = table.calibration;
    std::ostringstream s;
    s << "erodibility_fname: " << erodibility_fname << "\n"
      << "erodibility_raster_m_exponent: " << num(c.m_exponent) << "\n"
      << "erodibility_raster_n_exponent: " << num(c.n_exponent) << "\n"
      << "erodibility_raster_area_convention: "
      << (c.area_convention == erodibility_area_convention::map_area ? "map_area" : "specific_area")
      << "\n";
    if (c.cell_size > 0.0)
        s << "erodibility_raster_calibration_cell_size_metres: " << num(c.cell_size) << "\n";
    return s.str();
}

} // namespace lsdtt
