#pragma once

#include "lsdtt_basin_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"

#include <cstdint>
#include <vector>

namespace lsdtt
{

// Manages basins whose outlets are channel-head source pixels. These are
// deliberately separate from junction basins: a zero-order basin is the
// complete upslope area of a source pixel, not the basin of a downstream
// Strahler junction.
class lsdtt_zero_order_basin_manager_cpu
{
public:
    lsdtt_zero_order_basin_manager_cpu(const lsdtt_flow_base& flow);

    void select_all_zero_order_basins(
        const std::vector<std::int64_t>& source_flats
    );

    [[nodiscard]] std::size_t size() const noexcept { return basins_.size(); }
    [[nodiscard]] bool empty() const noexcept { return basins_.empty(); }

    [[nodiscard]] const lsdtt_basin_cpu& basin_at(std::size_t index) const;

    void write_basins_to_vector_lines(
        const std::string& out_path,
        vector_output_format format,
        const std::string& layer_name = "zero_order_basin_outlines",
        bool spatial_index = true
    ) const;

private:
    const lsdtt_flow_base* flow_;
    std::vector<std::int64_t> source_flats_;
    std::vector<lsdtt_basin_cpu> basins_;
};

} // namespace lsdtt
