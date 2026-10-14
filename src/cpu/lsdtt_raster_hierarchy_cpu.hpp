// lsdtt_raster_hierarchy_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Boris Gailleton
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

// In-memory resolution hierarchy for existing LSDTT float rasters.

#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_resampling_cpu.hpp"

namespace lsdtt {

enum class raster_hierarchy_storage
{
    active_only,
    retain_allocated,
    preallocated
};

// Geometry only. It deliberately has no value type or raster ownership.
class raster_hierarchy_layout_cpu
{
public:
    explicit raster_hierarchy_layout_cpu(const lsdtt_raster_meta& finest,
                                         std::size_t requested_levels = 0);

    [[nodiscard]] std::size_t n_levels() const noexcept { return levels_.size(); }
    [[nodiscard]] const lsdtt_raster_meta& level_meta(std::size_t level) const;
    // Register an exact same-extent level at finest resolution / factor.
    // Repeated factors return the existing level index.
    std::size_t add_level_for_factor(std::size_t factor);

private:
    std::vector<lsdtt_raster_meta> levels_;
};

// Future typed rasters can reuse this level storage and the layout above.
// It is intentionally independent of lsdtt_raster_cpu, which remains float.
template <typename T>
struct raster_hierarchy_typed_level_cpu
{
    std::vector<T> values;
    std::vector<std::uint8_t> mask;
    lsdtt_raster_meta metadata;
};

// A hierarchy owns one float field. Related fields (DEM, rainfall, friction)
// intentionally use separate hierarchies so their transfer policies stay explicit.
class lsdtt_raster_hierarchy_cpu
{
public:
    explicit lsdtt_raster_hierarchy_cpu(
        const lsdtt_raster_cpu& finest,
        raster_hierarchy_storage storage = raster_hierarchy_storage::retain_allocated,
        bool keep_finest = true,
        std::size_t requested_levels = 0);

    [[nodiscard]] std::size_t n_levels() const noexcept { return layout_.n_levels(); }
    [[nodiscard]] bool has_level(std::size_t level) const noexcept;
    [[nodiscard]] std::size_t active_level() const noexcept { return active_level_; }
    [[nodiscard]] const lsdtt_raster_meta& level_meta(std::size_t level) const
    {
        return layout_.level_meta(level);
    }
    std::size_t add_level_for_factor(std::size_t factor);

    // Allocate blank level buffers through `level`; no values are transferred.
    void allocate_to(std::size_t level);

    // Existing level access. Throws when the level is not retained/allocated.
    lsdtt_raster_cpu& level(std::size_t level);
    const lsdtt_raster_cpu& level(std::size_t level) const;

    // Restriction/prolongation move from the active level and make `target`
    // active. Restriction requires target > active; prolongation target < active.
    lsdtt_raster_cpu& restrict_to(
        std::size_t target,
        raster_resampling_method method = raster_resampling_method::conservative_mean);
    lsdtt_raster_cpu& prolong_to(
        std::size_t target,
        raster_resampling_method method = raster_resampling_method::linear);

    // Direction-neutral jump; the caller supplies the transfer method.
    lsdtt_raster_cpu& move_to(std::size_t target,
                              raster_resampling_method method);

    // Replace the active values after an external solver has advanced them.
    void set_active(const lsdtt_raster_cpu& raster);

    void clear_except(std::size_t level);
    void clear();

private:
    [[nodiscard]] lsdtt_raster_cpu make_empty(std::size_t level) const;
    lsdtt_raster_cpu& ensure_level(std::size_t level);
    void validate_level(std::size_t level) const;

    raster_hierarchy_layout_cpu layout_;
    raster_hierarchy_storage storage_;
    bool keep_finest_;
    std::vector<std::optional<lsdtt_raster_cpu>> levels_;
    std::size_t active_level_ = 0;
};

} // namespace lsdtt
