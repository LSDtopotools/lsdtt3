// lsdtt_erodibility_calibration_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_erodibility_calibration_cpu.hpp"

#include <cctype>
#include <cmath>
#include <filesystem>
#include <map>
#include <stdexcept>

namespace lsdtt
{

namespace
{

std::string lower_trimmed(const std::string& s)
{
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    std::string out = s.substr(a, b - a);
    for (char& c : out) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

// The attribute of a record whose key equals `column` ignoring case; nullptr
// when there is none.
const lsdtt_point_attribute_value* find_attribute(const lsdtt_input_point_with_attributes& rec,
                                                  const std::string& column)
{
    if (column.empty()) return nullptr;
    const auto exact = rec.attributes.find(column);
    if (exact != rec.attributes.end()) return &exact->second;
    const std::string want = lower_trimmed(column);
    for (const auto& [k, v] : rec.attributes)
        if (lower_trimmed(k) == want) return &v;
    return nullptr;
}

bool has_column(const std::vector<std::string>& columns, const std::string& column)
{
    if (column.empty()) return false;
    const std::string want = lower_trimmed(column);
    for (const auto& c : columns)
        if (lower_trimmed(c) == want) return true;
    return false;
}

std::string join_columns(const std::vector<std::string>& columns)
{
    std::string out;
    for (const auto& c : columns) out += (out.empty() ? "" : ", ") + c;
    return out.empty() ? "(none)" : out;
}

double numeric_or_nan(const lsdtt_point_attribute_value* v)
{
    double d = std::numeric_limits<double>::quiet_NaN();
    if (v == nullptr || !v->is_present()) return d;
    if (!v->as_double(d)) return std::numeric_limits<double>::quiet_NaN();
    return d;
}

} // namespace

std::string accepted_erosion_rate_units_list()
{
    return "mm_per_kyr, m_per_myr, cm_per_kyr, mm_per_yr, m_per_yr";
}

std::string canonical_erosion_rate_units(const std::string& units)
{
    const std::string u = lower_trimmed(units);
    if (u == "mm_per_kyr" || u == "m_per_myr") return "mm_per_kyr";
    if (u == "cm_per_kyr") return "cm_per_kyr";
    if (u == "mm_per_yr") return "mm_per_yr";
    if (u == "m_per_yr") return "m_per_yr";
    const bool mass_flux = u.find("g_per") != std::string::npos ||
                           u.find("g/cm") != std::string::npos ||
                           u.find("g cm") != std::string::npos ||
                           u.find("gcm") != std::string::npos;
    if (mass_flux)
        throw std::invalid_argument("erosion_rate_units '" + units +
                                    "' is a mass flux; convert it to a rate first by dividing by "
                                    "the rock density (e.g. 2.65 g cm^-3), then give one of " +
                                    accepted_erosion_rate_units_list());
    throw std::invalid_argument("erosion_rate_units must be one of " +
                                accepted_erosion_rate_units_list() + " (got '" + units + "')");
}

double erosion_rate_units_to_m_per_yr(const std::string& units)
{
    const std::string c = canonical_erosion_rate_units(units);
    if (c == "mm_per_kyr") return 1e-6;
    if (c == "cm_per_kyr") return 1e-5;
    if (c == "mm_per_yr") return 1e-3;
    return 1.0;
}

lsdtt_erosion_rate_read_result read_erosion_rate_samples(const lsdtt_raster_meta& meta,
                                                         const lsdtt_erosion_rate_input_options& o)
{
    lsdtt_erosion_rate_read_result res;
    res.units = canonical_erosion_rate_units(o.units);
    res.units_to_m_per_yr = erosion_rate_units_to_m_per_yr(o.units);
    if (o.path.empty() || !std::filesystem::exists(o.path))
        throw std::invalid_argument("erosion_rate_fname does not exist: " + o.path);

    lsdtt_point_input_options po;
    po.vector_layer_name = o.layer_name;
    po.csv_latitude_field_name = o.latitude_column;
    po.csv_longitude_field_name = o.longitude_column;
    po.keep_points_outside_raster = true;
    const lsdtt_point_attribute_read_result read = load_point_layer_with_attributes(o.path, meta, {}, po);
    res.columns = read.fields_found;

    if (!has_column(res.columns, o.rate_column))
        throw std::invalid_argument("erosion_rate_column '" + o.rate_column + "' is not a column of " +
                                    o.path + "; the columns are: " + join_columns(res.columns));
    if (!has_column(res.columns, o.uncertainty_column))
        throw std::invalid_argument("erosion_rate_uncertainty_column '" + o.uncertainty_column +
                                    "' is not a column of " + o.path + "; the columns are: " +
                                    join_columns(res.columns));
    res.has_sample_id = has_column(res.columns, o.sample_id_column);
    res.has_reported_area = has_column(res.columns, o.reported_area_column);
    res.has_units_column = has_column(res.columns, "erosion_rate_units");

    int row = 0;
    for (const auto& rec : read.records)
    {
        ++row;
        lsdtt_erosion_rate_sample s;
        s.input_row = row;
        s.point = rec.point;

        std::string id;
        const lsdtt_point_attribute_value* idv = find_attribute(rec, o.sample_id_column);
        if (idv && idv->is_present()) (void)idv->as_source_string(id);
        s.sample_id = lower_trimmed(id).empty() ? "row_" + std::to_string(row) : id;

        if (res.has_units_column)
        {
            std::string u;
            const lsdtt_point_attribute_value* uv = find_attribute(rec, "erosion_rate_units");
            if (uv && uv->is_present()) (void)uv->as_string(u);
            std::string canon;
            try { canon = canonical_erosion_rate_units(u); }
            catch (const std::invalid_argument&) { canon.clear(); }
            if (canon != res.units)
                throw std::invalid_argument(
                    "row " + std::to_string(row) + " (sample " + s.sample_id +
                    ") of " + o.path + ": its erosion_rate_units column says '" + u +
                    "', which disagrees with erosion_rate_units: " + o.units);
        }

        s.rate_input = numeric_or_nan(find_attribute(rec, o.rate_column));
        s.uncertainty_input = numeric_or_nan(find_attribute(rec, o.uncertainty_column));
        if (std::isfinite(s.rate_input)) s.E_m_per_yr = s.rate_input * res.units_to_m_per_yr;
        if (std::isfinite(s.uncertainty_input))
            s.sigma_E_m_per_yr = s.uncertainty_input * res.units_to_m_per_yr;
        if (res.has_reported_area)
        {
            const double a = numeric_or_nan(find_attribute(rec, o.reported_area_column));
            if (std::isfinite(a) && a > 0.0) s.reported_area_m2 = a * 1.0e6;
        }

        if (!std::isfinite(s.rate_input) || !(s.rate_input > 0.0))
        {
            s.status = "bad_rate";
            s.status_detail = "erosion rate missing, not a number or <= 0";
        }
        else if (!std::isfinite(s.uncertainty_input) || !(s.uncertainty_input > 0.0))
        {
            s.status = "missing_uncertainty";
            s.status_detail = "uncertainty missing, not a number or <= 0";
        }
        res.samples.push_back(std::move(s));
    }
    return res;
}

} // namespace lsdtt
