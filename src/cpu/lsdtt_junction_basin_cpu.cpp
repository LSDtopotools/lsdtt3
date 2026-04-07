// lsdtt3 — lsdtt_junction_basin_cpu.cpp
// CPU implementation of a junction-network-aware drainage basin.
//
// Copyright Simon M. Mudd, GPL-3.0

#include "lsdtt_junction_basin_cpu.hpp"

namespace lsdtt
{

lsdtt_junction_basin_cpu::lsdtt_junction_basin_cpu(
    std::int64_t                     outlet_flat_index,
    const lsdtt_flowinfo_cpu&         flowinfo,
    const lsdtt_junction_network_cpu& junction_network)
    : lsdtt_basin_cpu(outlet_flat_index, flowinfo)
    , junction_number_(-1)
    , stream_order_(-1)
{
    // TODO: walk upstream from outlet_flat_index through flowinfo until the
    // first junction pixel is found in junction_network. Record its junction
    // index as junction_number_ and look up its stream order from the network.
    // If the outlet is itself a junction, junction_number_ is that junction.
}


std::int32_t lsdtt_junction_basin_cpu::get_junction_number() const noexcept
{
    return junction_number_;
}

std::int32_t lsdtt_junction_basin_cpu::get_stream_order() const noexcept
{
    return stream_order_;
}

} // namespace lsdtt