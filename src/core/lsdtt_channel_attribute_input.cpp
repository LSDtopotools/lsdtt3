// lsdtt_channel_attribute_input.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_channel_attribute_input.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <variant>

#include "../utils/lsdtt_raster_geometry.hpp"
#include "lsdtt_csv_table.hpp"
#include "lsdtt_data_loader.hpp"
#include "lsdtt_vector_utils.hpp"

namespace lsdtt
{

namespace
{

const double k_nan = std::numeric_limits<double>::quiet_NaN();

std::string lower_ascii(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

bool has_csv_suffix(const std::string& path)
{
    return path.size() >= 4 && lower_ascii(path.substr(path.size() - 4)) == ".csv";
}

// The header spelled as requested, ignoring case; empty if absent.
std::string find_header(const std::vector<std::string>& headers, const std::string& wanted)
{
    for (const auto& h : headers)
        if (lower_ascii(h) == wanted) return h;
    return "";
}

// A blank, NaN or infinite token is NaN; any other non-number is an error.
double parse_value(const std::string& token, const std::string& what, std::size_t line)
{
    if (token.empty()) return k_nan;
    char* end = nullptr;
    const double v = std::strtod(token.c_str(), &end);
    if (end == token.c_str() || *end != '\0')
        throw std::invalid_argument("read_channel_attribute_table: " + what + " on line " +
                                    std::to_string(line) + " is not numeric ('" + token + "')");
    return std::isfinite(v) ? v : k_nan;
}

void check_exact_pixel(const lsdtt_raster_meta& meta, double row, double col, double flat)
{
    const double nr = meta.n_rows, nc = meta.n_cols;
    if (!(row >= 0 && row < nr && col >= 0 && col < nc) || row != std::floor(row) ||
        col != std::floor(col) || flat != row * nc + col)
        throw std::invalid_argument(
            "read_channel_attribute_table: row, col and flat disagree with this " +
            std::to_string(meta.n_rows) + " x " + std::to_string(meta.n_cols) +
            " grid. The table was written on a different grid.");
}

void check_on_pixel(const lsdtt_raster_meta& meta, int row, int col, double x, double y)
{
    const lsdtt_xy centre = cell_center_from_row_col(meta, row, col);
    const double limit = 0.75 * meta.cell_size;
    if (std::abs(x - centre.x) > limit || std::abs(y - centre.y) > limit)
        throw std::invalid_argument(
            "read_channel_attribute_table: a record's coordinates are not on its row and "
            "col. The table was written on a different grid.");
}

void add_line_samples(channel_attribute_table& t, std::int64_t record,
                      const std::vector<lsdtt_xy>& vertices, double spacing)
{
    const auto add = [&](double x, double y) {
        t.sample_record.push_back(record);
        t.sample_flat.push_back(-1);
        t.sample_x.push_back(x);
        t.sample_y.push_back(y);
    };
    if (!vertices.empty()) add(vertices.front().x, vertices.front().y);
    for (std::size_t k = 1; k < vertices.size(); ++k)
    {
        const lsdtt_xy& a = vertices[k - 1];
        const lsdtt_xy& b = vertices[k];
        const double length = std::hypot(b.x - a.x, b.y - a.y);
        const int steps = std::max(1, static_cast<int>(std::ceil(length / spacing)));
        for (int s = 1; s <= steps; ++s)
        {
            const double w = static_cast<double>(s) / steps;
            add(a.x + w * (b.x - a.x), a.y + w * (b.y - a.y));
        }
    }
}

void read_csv(const std::string& path, const lsdtt_raster_meta& meta, channel_attribute_table& t)
{
    const lsdtt_csv_table csv = load_csv_table(path);
    for (const auto& name : t.columns)
        if (!csv_has_column(csv, name))
            throw std::invalid_argument("read_channel_attribute_table: column '" + name +
                                        "' is not in " + path);

    const std::string h_row = find_header(csv.headers, "row");
    const std::string h_col = find_header(csv.headers, "col");
    const std::string h_flat = find_header(csv.headers, "flat");
    const std::string h_lat = find_header(csv.headers, "latitude");
    const std::string h_lon = find_header(csv.headers, "longitude");
    const bool exact = !h_row.empty() && !h_col.empty() && !h_flat.empty();
    const bool coords = !h_lat.empty() && !h_lon.empty();
    if (!exact && !coords)
        throw std::invalid_argument("read_channel_attribute_table: " + path +
                                    " needs row, col and flat, or latitude and longitude");

    lsdtt_coordinate_transform transform;
    if (coords)
    {
        if (meta.crs_wkt.empty())
            throw std::invalid_argument(
                "read_channel_attribute_table: the raster has no CRS for latitude and longitude");
        transform.initialise_from_source_epsg(4326, meta.crs_wkt);
    }

    t.n_records = static_cast<std::int64_t>(csv.rows.size());
    for (auto& column : t.values) column.reserve(csv.rows.size());

    std::size_t line = 1;
    for (const auto& row : csv.rows)
    {
        ++line;
        const std::int64_t record = static_cast<std::int64_t>(line - 2);
        for (std::size_t c = 0; c < t.columns.size(); ++c)
            t.values[c].push_back(
                parse_value(csv_row_string(row, t.columns[c]), "column " + t.columns[c], line));

        double x = k_nan, y = k_nan;
        if (coords)
        {
            const double lat = parse_value(csv_row_string(row, h_lat), "latitude", line);
            const double lon = parse_value(csv_row_string(row, h_lon), "longitude", line);
            if (std::isfinite(lat) && std::isfinite(lon))
            {
                const lsdtt_xy p = transform.transform_xy(lon, lat);
                x = p.x;
                y = p.y;
            }
        }

        std::int64_t flat = -1;
        if (exact)
        {
            const double r = parse_value(csv_row_string(row, h_row), "row", line);
            const double c = parse_value(csv_row_string(row, h_col), "col", line);
            const double f = parse_value(csv_row_string(row, h_flat), "flat", line);
            check_exact_pixel(meta, r, c, f);
            if (std::isfinite(x)) check_on_pixel(meta, static_cast<int>(r), static_cast<int>(c), x, y);
            flat = static_cast<std::int64_t>(f);
        }
        else if (!std::isfinite(x))
        {
            continue;  // no location: the record stays unmatched and is counted
        }
        t.sample_record.push_back(record);
        t.sample_flat.push_back(flat);
        t.sample_x.push_back(x);
        t.sample_y.push_back(y);
    }
}

double field_as_number(const lsdtt_vector_field_value& v)
{
    if (const auto* i = std::get_if<std::int64_t>(&v)) return static_cast<double>(*i);
    if (const auto* d = std::get_if<double>(&v)) return std::isfinite(*d) ? *d : k_nan;
    return k_nan;
}

void read_vector(const std::string& path, const lsdtt_raster_meta& meta, channel_attribute_table& t)
{
    lsdtt_data_loader loader;
    loader.open_vector(path);
    const auto names = loader.get_layer_names();
    if (names.empty())
        throw std::invalid_argument("read_channel_attribute_table: " + path + " has no layers");
    const lsdtt_vector_layer layer = loader.read_vector_layer(names.front());
    const std::string layer_wkt = loader.get_vector_layer_crs_wkt(names.front());

    std::vector<int> index;
    for (const auto& name : t.columns)
    {
        const int k = layer.schema.index_of(name);
        if (k < 0)
            throw std::invalid_argument("read_channel_attribute_table: column '" + name +
                                        "' is not in " + path);
        if (layer.schema.field_types[static_cast<std::size_t>(k)] ==
            lsdtt_vector_field_type::string_type)
            throw std::invalid_argument("read_channel_attribute_table: column '" + name +
                                        "' is text, not numeric");
        index.push_back(k);
    }

    const int i_row = layer.schema.index_of("row");
    const int i_col = layer.schema.index_of("col");
    const int i_flat = layer.schema.index_of("flat");
    const bool exact_columns = i_row >= 0 && i_col >= 0 && i_flat >= 0;

    bool needs_crs = false;
    for (const auto& rec : layer.records)
    {
        if (std::holds_alternative<lsdtt_polygon_geom>(rec.geometry))
            throw std::invalid_argument(
                "read_channel_attribute_table: polygons cannot be joined to channel pixels; "
                "use points or lines");
        if (!std::holds_alternative<std::monostate>(rec.geometry)) needs_crs = true;
    }
    if (needs_crs && layer_wkt.empty())
        throw std::invalid_argument("read_channel_attribute_table: " + path + " has no CRS");
    if (needs_crs && meta.crs_wkt.empty())
        throw std::invalid_argument("read_channel_attribute_table: the raster has no CRS");

    lsdtt_coordinate_transform transform;
    if (needs_crs) transform.initialise(layer_wkt, meta.crs_wkt);
    const auto to_raster = [&](double x, double y) { return transform.transform_xy(x, y); };
    const auto to_raster_all = [&](const std::vector<lsdtt_xy>& in) {
        std::vector<lsdtt_xy> out;
        out.reserve(in.size());
        for (const auto& p : in) out.push_back(to_raster(p.x, p.y));
        return out;
    };

    // Parts of a multipart feature are one record.
    std::int64_t previous = -1, record = -1;
    const double spacing = 0.5 * meta.cell_size;
    for (std::size_t k = 0; k < layer.records.size(); ++k)
    {
        const lsdtt_vector_record& rec = layer.records[k];
        const std::int64_t feature = layer.source_feature_index.empty()
                                         ? static_cast<std::int64_t>(k)
                                         : layer.source_feature_index[k];
        if (k == 0 || feature != previous)
        {
            ++record;
            previous = feature;
            for (std::size_t c = 0; c < index.size(); ++c)
                t.values[c].push_back(
                    field_as_number(rec.field_values[static_cast<std::size_t>(index[c])]));
        }

        if (const auto* line = std::get_if<lsdtt_line_geom>(&rec.geometry))
        {
            add_line_samples(t, record, to_raster_all(line->points), spacing);
            continue;
        }
        if (const auto* multi = std::get_if<lsdtt_multiline_geom>(&rec.geometry))
        {
            for (const auto& part : multi->parts)
                add_line_samples(t, record, to_raster_all(part), spacing);
            continue;
        }

        double x = k_nan, y = k_nan;
        if (const auto* point = std::get_if<lsdtt_point_geom>(&rec.geometry))
        {
            const lsdtt_xy p = to_raster(point->x, point->y);
            x = p.x;
            y = p.y;
        }

        std::int64_t flat = -1;
        if (exact_columns)
        {
            const double r = field_as_number(rec.field_values[static_cast<std::size_t>(i_row)]);
            const double c = field_as_number(rec.field_values[static_cast<std::size_t>(i_col)]);
            const double f = field_as_number(rec.field_values[static_cast<std::size_t>(i_flat)]);
            check_exact_pixel(meta, r, c, f);
            if (std::isfinite(x)) check_on_pixel(meta, static_cast<int>(r), static_cast<int>(c), x, y);
            flat = static_cast<std::int64_t>(f);
        }
        else if (!std::isfinite(x))
        {
            continue;
        }
        t.sample_record.push_back(record);
        t.sample_flat.push_back(flat);
        t.sample_x.push_back(x);
        t.sample_y.push_back(y);
    }
    t.n_records = record + 1;
}

} // namespace

channel_attribute_table read_channel_attribute_table(
    const std::string&              path,
    const lsdtt_raster_meta&        meta,
    const std::vector<std::string>& columns,
    double                          nodata_value)
{
    if (columns.empty())
        throw std::invalid_argument("read_channel_attribute_table: no columns requested");
    std::vector<std::string> sorted = columns;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
        throw std::invalid_argument("read_channel_attribute_table: a column is repeated");

    channel_attribute_table t;
    t.columns = columns;
    t.values.resize(columns.size());
    if (has_csv_suffix(path)) read_csv(path, meta, t);
    else                      read_vector(path, meta, t);

    const float nodata = static_cast<float>(nodata_value);
    for (auto& column : t.values)
        for (double& v : column)
            if (std::isfinite(v) && static_cast<float>(v) == nodata) v = k_nan;
    return t;
}

} // namespace lsdtt
