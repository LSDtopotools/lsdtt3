#pragma once

// lsdtt3 — lsdtt_junction_basin_cpu.hpp
// A basin linked to the junction network.
// Inherits lsdtt_basin_cpu and adds junction-network context:
// the junction number and associated stream order.
//
// An lsdtt_junction_basin_cpu can be selected by any outlet flat index,
// not only a junction pixel. On construction the first junction
// encountered in the upstream direction is recorded and its stream
// order is stored. If the outlet is itself a junction the two coincide.
//
// lsdtt_cosmo_basin_cpu is a sibling child of lsdtt_basin_cpu (not of
// this class) because cosmogenic samples are not tied to junctions.
//
// Copyright Simon M. Mudd, GPL-3.0

#include "lsdtt_basin_cpu.hpp"
#include "lsdtt_junction_network_cpu.hpp"

#include <cstdint>
#include <string>

namespace lsdtt
{

// -----------------------------------------------------------------------
// lsdtt_junction_basin_cpu
// -----------------------------------------------------------------------

class lsdtt_junction_basin_cpu : public lsdtt_basin_cpu
{
public:

    // Constructs a junction basin from a validated outlet flat index.
    // Both flowinfo and junction_network must outlive this object.
    lsdtt_junction_basin_cpu(std::int64_t                      outlet_flat_index,
                             const lsdtt_flowinfo_cpu&          flowinfo,
                             const lsdtt_junction_network_cpu&  junction_network);

    // -----------------------------------------------------------------------
    // Junction-specific accessors
    // -----------------------------------------------------------------------

    // The junction number of the first junction encountered upstream of
    // (or at) the outlet. Used by the basin manager for network-aware
    // filtering and for linking basins to the junction graph.
    [[nodiscard]] std::int32_t get_junction_number() const noexcept;

    // The Strahler stream order associated with junction_number_.
    [[nodiscard]] std::int32_t get_stream_order() const noexcept;

private:

    std::int32_t junction_number_;   // index into junction network
    std::int32_t stream_order_;      // Strahler order at junction_number_
};

} // namespace lsdtt