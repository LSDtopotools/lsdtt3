// lsdtt_raster_hierarchy_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Boris Gailleton
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_raster_hierarchy_cpu.hpp"

#include "../core/lsdtt_raster_warp.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace lsdtt {
namespace {

std::size_t maximum_dyadic_levels(const lsdtt_raster_meta& meta)
{
    std::size_t count = 1;
    int rows = meta.n_rows;
    int cols = meta.n_cols;
    while (rows > 1 && cols > 1 && rows % 2 == 0 && cols % 2 == 0)
    {
        rows /= 2;
        cols /= 2;
        ++count;
    }
    return count;
}

} // namespace

raster_hierarchy_layout_cpu::raster_hierarchy_layout_cpu(
    const lsdtt_raster_meta& finest, std::size_t requested_levels)
{
    if (!finest.is_valid())
        throw std::invalid_argument("raster hierarchy needs valid finest metadata");

    const std::size_t available = maximum_dyadic_levels(finest);
    const std::size_t count = requested_levels == 0 ? available : requested_levels;
    if (count == 0 || count > available)
        throw std::invalid_argument(
            "requested hierarchy levels exceed exact factor-two levels available "
            "for this square-cell raster");

    levels_.reserve(count);
    auto meta = finest;
    for (std::size_t level = 0; level < count; ++level)
    {
        levels_.push_back(meta);
        if (level + 1 < count)
        {
            meta.n_rows /= 2;
            meta.n_cols /= 2;
            meta.cell_size *= 2.0;
        }
    }
}

const lsdtt_raster_meta& raster_hierarchy_layout_cpu::level_meta(std::size_t level) const
{
    if (level >= levels_.size()) throw std::out_of_range("hierarchy level is out of range");
    return levels_[level];
}

std::size_t raster_hierarchy_layout_cpu::add_level_for_factor(std::size_t factor)
{
    if (factor == 0) throw std::invalid_argument("hierarchy factor must be positive");
    const auto& finest = levels_.front();
    if (finest.n_rows % static_cast<int>(factor) != 0 ||
        finest.n_cols % static_cast<int>(factor) != 0)
        throw std::invalid_argument("hierarchy factor must divide both finest dimensions");

    const int rows = finest.n_rows / static_cast<int>(factor);
    const int cols = finest.n_cols / static_cast<int>(factor);
    for (std::size_t i = 0; i < levels_.size(); ++i)
        if (levels_[i].n_rows == rows && levels_[i].n_cols == cols) return i;

    auto meta = finest;
    meta.n_rows = rows;
    meta.n_cols = cols;
    meta.cell_size *= static_cast<double>(factor);
    levels_.push_back(std::move(meta));
    return levels_.size() - 1;
}

lsdtt_raster_hierarchy_cpu::lsdtt_raster_hierarchy_cpu(
    const lsdtt_raster_cpu& finest, raster_hierarchy_storage storage,
    bool keep_finest, std::size_t requested_levels)
    : layout_(finest.get_meta(), requested_levels), storage_(storage), keep_finest_(keep_finest),
      levels_(layout_.n_levels())
{
    levels_[0] = finest;
    if (storage_ == raster_hierarchy_storage::preallocated)
        allocate_to(layout_.n_levels() - 1);
}

bool lsdtt_raster_hierarchy_cpu::has_level(std::size_t level) const noexcept
{
    return level < levels_.size() && levels_[level].has_value();
}

void lsdtt_raster_hierarchy_cpu::validate_level(std::size_t level) const
{
    if (level >= levels_.size()) throw std::out_of_range("hierarchy level is out of range");
}

lsdtt_raster_cpu lsdtt_raster_hierarchy_cpu::make_empty(std::size_t level) const
{
    const auto& meta = layout_.level_meta(level);
    const std::size_t size = static_cast<std::size_t>(meta.n_pixels());
    return lsdtt_raster_cpu(
        std::vector<float>(size, static_cast<float>(meta.no_data_value)),
        std::vector<std::uint8_t>(size, 0), meta);
}

lsdtt_raster_cpu& lsdtt_raster_hierarchy_cpu::ensure_level(std::size_t level)
{
    validate_level(level);
    if (!levels_[level]) levels_[level] = make_empty(level);
    return *levels_[level];
}

void lsdtt_raster_hierarchy_cpu::allocate_to(std::size_t level)
{
    validate_level(level);
    for (std::size_t i = 0; i <= level; ++i) ensure_level(i);
}

std::size_t lsdtt_raster_hierarchy_cpu::add_level_for_factor(std::size_t factor)
{
    const std::size_t level = layout_.add_level_for_factor(factor);
    if (level == levels_.size()) levels_.emplace_back();
    if (storage_ == raster_hierarchy_storage::preallocated) ensure_level(level);
    return level;
}

lsdtt_raster_cpu& lsdtt_raster_hierarchy_cpu::level(std::size_t level)
{
    validate_level(level);
    if (!levels_[level])
        throw std::runtime_error("hierarchy level is not allocated; call allocate_to or move_to");
    return *levels_[level];
}

const lsdtt_raster_cpu& lsdtt_raster_hierarchy_cpu::level(std::size_t level) const
{
    validate_level(level);
    if (!levels_[level])
        throw std::runtime_error("hierarchy level is not allocated; call allocate_to or move_to");
    return *levels_[level];
}

lsdtt_raster_cpu& lsdtt_raster_hierarchy_cpu::move_to(
    std::size_t target, raster_resampling_method method)
{
    validate_level(target);
    if (target == active_level_) return level(target);

    const std::size_t source = active_level_;
    lsdtt_raster_cpu destination = make_empty(target);
    resample_raster(level(source), destination, method);
    levels_[target] = std::move(destination);
    active_level_ = target;

    if (storage_ == raster_hierarchy_storage::active_only)
    {
        if (!(keep_finest_ && source == 0)) levels_[source].reset();
        for (std::size_t i = 0; i < levels_.size(); ++i)
            if (i != target && !(keep_finest_ && i == 0)) levels_[i].reset();
    }
    return *levels_[target];
}

lsdtt_raster_cpu& lsdtt_raster_hierarchy_cpu::restrict_to(
    std::size_t target, raster_resampling_method method)
{
    if (layout_.level_meta(target).cell_size <= layout_.level_meta(active_level_).cell_size)
        throw std::invalid_argument("restriction target must be coarser than the active level");
    return move_to(target, method);
}

lsdtt_raster_cpu& lsdtt_raster_hierarchy_cpu::prolong_to(
    std::size_t target, raster_resampling_method method)
{
    if (layout_.level_meta(target).cell_size >= layout_.level_meta(active_level_).cell_size)
        throw std::invalid_argument("prolongation target must be finer than the active level");
    return move_to(target, method);
}

void lsdtt_raster_hierarchy_cpu::set_active(const lsdtt_raster_cpu& raster)
{
    const auto& expected = layout_.level_meta(active_level_);
    const auto& actual = raster.get_meta();
    if (actual.n_rows != expected.n_rows || actual.n_cols != expected.n_cols ||
        actual.x_min != expected.x_min || actual.y_min != expected.y_min ||
        actual.cell_size != expected.cell_size ||
        !crs_equal_semantically(actual.crs_wkt, expected.crs_wkt))
        throw std::invalid_argument("replacement raster does not match active hierarchy geometry");
    levels_[active_level_] = raster;
}

void lsdtt_raster_hierarchy_cpu::clear_except(std::size_t level)
{
    validate_level(level);
    if (!levels_[level])
        throw std::runtime_error("cannot retain a hierarchy level that is not allocated");
    for (std::size_t i = 0; i < levels_.size(); ++i)
        if (i != level) levels_[i].reset();
    active_level_ = level;
}

void lsdtt_raster_hierarchy_cpu::clear()
{
    for (auto& entry : levels_) entry.reset();
}

} // namespace lsdtt
