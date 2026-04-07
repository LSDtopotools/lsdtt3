// lsdtt3 — lsdtt_junction_basin_cpu.cpp
// CPU implementation of a junction-network-aware drainage basin.
//
// Copyright Simon M. Mudd, GPL-3.0

#include "lsdtt_junction_basin_cpu.hpp"

namespace lsdtt
{

lsdtt_junction_basin_cpu::lsdtt_junction_basin_cpu(
   std::int64_t                       basin_outlet_flat_index,
   const lsdtt_flow_base&             flow,
   const lsdtt_junction_network_base& junction_network)
   : lsdtt_basin_cpu(basin_outlet_flat_index, flow)
   , junction_number_(-1)
   , stream_order_(-1)
   , junction_network_(&junction_network)
{
   junction_number_ = resolve_junction_number(junction_network, basin_outlet_flat_index);
   stream_order_    = resolve_stream_order(junction_network, junction_number_);
}

std::int32_t lsdtt_junction_basin_cpu::get_junction_number() const noexcept
{
   return junction_number_;
}

std::int32_t lsdtt_junction_basin_cpu::get_stream_order() const noexcept
{
   return stream_order_;
}

const lsdtt_junction_network_base&
lsdtt_junction_basin_cpu::get_junction_network() const noexcept
{
   return *junction_network_;
}

std::int32_t lsdtt_junction_basin_cpu::resolve_junction_number(
   const lsdtt_junction_network_base& junction_network,
   std::int64_t basin_outlet_flat_index)
{
   // This lookup deliberately belongs to the junction-network interface,
   // not to the basin. The network owns the channel/junction topology.
   return junction_network.junction_upstream_of_flat(basin_outlet_flat_index);
}

std::int32_t lsdtt_junction_basin_cpu::resolve_stream_order(
   const lsdtt_junction_network_base& junction_network,
   std::int32_t junction_number)
{
   if (junction_number < 0)
   {
       return -1;
   }

   return junction_network.stream_order_of_junction(junction_number);
}

} // namespace lsdtt
