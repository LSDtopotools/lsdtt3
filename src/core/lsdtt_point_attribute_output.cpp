// lsdtt_point_attribute_output.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_point_attribute_output.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_point_attribute_output.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

#include "lsdtt_vector_utils.hpp"

namespace lsdtt
{
namespace
{
std::string lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

bool name_taken(const lsdtt_vector_schema& schema, const std::string& name)
{
    const std::string key = lower(name);
    for (const std::string& n : schema.field_names)
        if (lower(n) == key) return true;
    return false;
}
} // namespace

lsdtt_vector_schema point_attribute_schema(const lsdtt_point_attribute_read_result& read)
{
    lsdtt_vector_schema schema;
    for (const std::string& name : read.fields_found)
    {
        bool any_text = false, any_real = false, any_int = false, int_spelling_changes = false;
        for (const lsdtt_input_point_with_attributes& rec : read.records)
        {
            const auto it = rec.attributes.find(name);
            if (it == rec.attributes.end()) continue;
            const lsdtt_point_attribute_value& v = it->second;
            switch (v.type)
            {
            case lsdtt_point_attribute_value::kind::text: any_text = true; break;
            case lsdtt_point_attribute_value::kind::real: any_real = true; break;
            case lsdtt_point_attribute_value::kind::integer:
                any_int = true;
                if (v.has_source_text && v.source_text != std::to_string(v.i)) int_spelling_changes = true;
                break;
            case lsdtt_point_attribute_value::kind::absent: break;
            }
        }
        lsdtt_vector_field_type t = lsdtt_vector_field_type::string_type;
        if (any_text || (any_int && int_spelling_changes && !any_real))
            t = lsdtt_vector_field_type::string_type;
        else if (any_real)
            t = lsdtt_vector_field_type::double_type;
        else if (any_int)
            t = lsdtt_vector_field_type::int64_type;
        schema.field_names.push_back(name);
        schema.field_types.push_back(t);
    }
    return schema;
}

std::vector<lsdtt_vector_field_value> point_attribute_values(
    const lsdtt_input_point_with_attributes& record,
    const lsdtt_vector_schema&               schema,
    bool                                     as_text)
{
    std::vector<lsdtt_vector_field_value> out;
    out.reserve(schema.size());
    for (std::size_t k = 0; k < schema.size(); ++k)
    {
        const auto it = record.attributes.find(schema.field_names[k]);
        if (it == record.attributes.end() || !it->second.is_present())
        {
            out.emplace_back(std::monostate{});
            continue;
        }
        const lsdtt_point_attribute_value& v = it->second;
        if (as_text)
        {
            std::string s;
            if (v.as_source_string(s)) out.emplace_back(std::move(s));
            else out.emplace_back(std::monostate{});
            continue;
        }
        switch (schema.field_types[k])
        {
        case lsdtt_vector_field_type::int64_type:
        {
            std::int64_t i = 0;
            if (v.as_int64(i)) out.emplace_back(i);
            else out.emplace_back(std::monostate{});
            break;
        }
        case lsdtt_vector_field_type::double_type:
        {
            double d = 0.0;
            if (v.as_double(d)) out.emplace_back(d);
            else out.emplace_back(std::monostate{});
            break;
        }
        case lsdtt_vector_field_type::string_type:
        {
            std::string s;
            if (v.as_source_string(s)) out.emplace_back(std::move(s));
            else out.emplace_back(std::monostate{});
            break;
        }
        }
    }
    return out;
}

std::vector<std::pair<std::string, std::string>> append_point_fields(
    lsdtt_vector_schema&                                                 schema,
    const std::vector<std::pair<std::string, lsdtt_vector_field_type>>& fields)
{
    std::vector<std::pair<std::string, std::string>> renamed;
    for (const auto& f : fields)
    {
        std::string name = f.first;
        if (name_taken(schema, name))
        {
            name += "_lsdtt";
            while (name_taken(schema, name)) name += "_lsdtt";
            renamed.emplace_back(f.first, name);
        }
        schema.field_names.push_back(name);
        schema.field_types.push_back(f.second);
    }
    return renamed;
}

std::size_t write_points_as_layer(
    lsdtt_data_loader&                                        loader,
    const std::string&                                        crs_wkt,
    const std::string&                                        out_path,
    const std::string&                                        layer_name,
    const lsdtt_vector_schema&                                schema,
    point_output_format                                       format,
    bool                                                      spatial_index,
    const std::vector<lsdtt_xy>&                              xy,
    const std::vector<std::vector<lsdtt_vector_field_value>>& values)
{
    if (xy.size() != values.size())
        throw std::invalid_argument("write_points_as_layer: xy and values differ in length");
    for (const auto& row : values)
        if (row.size() != schema.size())
            throw std::invalid_argument("write_points_as_layer: a values row does not match the schema");

    if (format == point_output_format::csv)
    {
        if (crs_wkt.empty())
            throw std::runtime_error("write_points_as_layer: CSV output requires a valid CRS");
        const lsdtt_coordinate_transform to_wgs84 = lsdtt_coordinate_transform::from_epsg(crs_wkt, 4326);
        std::vector<std::vector<std::string>> rows;
        rows.reserve(xy.size());
        for (std::size_t i = 0; i < xy.size(); ++i)
        {
            const lsdtt_lonlat ll = xy_to_lonlat(xy[i].x, xy[i].y, to_wgs84);
            std::vector<std::string> cells;
            cells.reserve(schema.size() + 2);
            cells.push_back(csv_number(ll.latitude));
            cells.push_back(csv_number(ll.longitude));
            for (const lsdtt_vector_field_value& v : values[i]) cells.push_back(vector_field_value_to_string(v));
            rows.push_back(std::move(cells));
        }
        std::vector<std::string> headers = {"latitude", "longitude"};
        headers.insert(headers.end(), schema.field_names.begin(), schema.field_names.end());
        loader.write_csv(out_path, headers, rows);
        return xy.size();
    }

    loader.begin_vector_write(out_path, layer_name, schema, crs_wkt,
                              convert_point_to_vector_output_format(format), spatial_index);
    for (std::size_t i = 0; i < xy.size(); ++i) loader.write_point(xy[i].x, xy[i].y, values[i]);
    loader.finish_vector_write();
    return xy.size();
}

} // namespace lsdtt
