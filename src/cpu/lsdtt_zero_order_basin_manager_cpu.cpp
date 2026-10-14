// lsdtt_zero_order_basin_manager_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_zero_order_basin_manager_cpu.hpp"

#include "lsdtt_basin_outline_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <variant>

namespace lsdtt
{

lsdtt_zero_order_basin_manager_cpu::lsdtt_zero_order_basin_manager_cpu(
    const lsdtt_flow_base& flow)
    : flow_(&flow)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "lsdtt_zero_order_basin_manager_cpu: input flow object is empty");
    }
}

void lsdtt_zero_order_basin_manager_cpu::select_all_zero_order_basins(
    const std::vector<std::int64_t>& source_flats)
{
    source_flats_ = source_flats;
    std::sort(source_flats_.begin(), source_flats_.end());
    source_flats_.erase(
        std::unique(source_flats_.begin(), source_flats_.end()),
        source_flats_.end());

    basins_.clear();
    basins_.reserve(source_flats_.size());
    for (const std::int64_t source_flat : source_flats_)
    {
        if (!flow_->is_valid_cell(source_flat)) continue;
        basins_.emplace_back(source_flat, *flow_);
    }
}

const lsdtt_basin_cpu&
lsdtt_zero_order_basin_manager_cpu::basin_at(std::size_t index) const
{
    if (index >= basins_.size())
    {
        throw std::out_of_range(
            "lsdtt_zero_order_basin_manager_cpu::basin_at: index out of range");
    }
    return basins_[index];
}

void lsdtt_zero_order_basin_manager_cpu::write_basins_to_vector_lines(
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name,
    bool spatial_index) const
{
    if (basins_.empty())
    {
        throw std::runtime_error(
            "lsdtt_zero_order_basin_manager_cpu::write_basins_to_vector_lines: "
            "no zero-order basins have been selected");
    }

    lsdtt_vector_schema schema;
    schema.field_names = {"zero_order_id", "source_flat"};
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type};

    lsdtt_data_loader loader;
    loader.begin_vector_write(
        out_path, layer_name, schema, flow_->meta().crs_wkt, format,
        spatial_index);

    const double max_dx = std::abs(flow_->meta().cell_size) *
        flow_->meta().n_cols * 0.5;
    const double max_dy = std::abs(flow_->meta().cell_size) *
        flow_->meta().n_rows * 0.5;

    for (std::size_t i = 0; i < basins_.size(); ++i)
    {
        const lsdtt_basin_cpu& basin = basins_[i];
        const lsdtt_basin_outline_cpu outline = basin.get_basin_outline();
        const lsdtt_vector_record record =
            outline.get_outline_as_line_record(static_cast<std::int32_t>(i));

        if (std::holds_alternative<std::monostate>(record.geometry)) continue;
        const lsdtt_line_geom& line = std::get<lsdtt_line_geom>(record.geometry);
        if (line.points.size() < 2) continue;

        std::vector<lsdtt_vector_field_value> fields = {
            static_cast<std::int64_t>(i),
            basin.get_basin_outlet_flat_index()};

        if (outline.is_periodic_ew() || outline.is_periodic_ns())
        {
            std::vector<std::vector<lsdtt_xy>> paths;
            std::vector<lsdtt_xy> current;
            current.push_back(line.points.front());
            for (std::size_t p = 1; p < line.points.size(); ++p)
            {
                const double dx = std::abs(line.points[p].x - line.points[p - 1].x);
                const double dy = std::abs(line.points[p].y - line.points[p - 1].y);
                if (dx > max_dx || dy > max_dy)
                {
                    paths.push_back(std::move(current));
                    current.clear();
                }
                current.push_back(line.points[p]);
            }
            if (!current.empty()) paths.push_back(std::move(current));
            for (const auto& path : paths)
                if (path.size() >= 2) loader.write_line(path, fields);
        }
        else
        {
            loader.write_line(line.points, fields);
        }
    }

    loader.finish_vector_write();
}

} // namespace lsdtt
