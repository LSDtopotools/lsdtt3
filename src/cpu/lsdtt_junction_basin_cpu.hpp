#pragma once

// lsdtt3 - lsdtt_basin_cpu.hpp
// CPU implementation of a single drainage basin.
// Holds a zero-copy slice into the flow-routing s_vector
// (lsdtt_flat_index_slice) - no pixel data is copied on construction.
//
// Member functions are limited to operations that require access to
// private data (the slice and outlet index). Statistical analysis on
// extracted value vectors belongs in free functions in src/cpu/,
// following the design rule: "write free functions unless you need
// private data members".
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_basin_cpu.hpp"
#include "../core/lsdtt_junction_network_base.hpp"

#include <cstdint>

namespace lsdtt
{

// -----------------------------------------------------------------------
// lsdtt_junction_basin_cpu
// -----------------------------------------------------------------------
class lsdtt_junction_basin_cpu : public lsdtt_basin_cpu
{
public:

    // Constructs a junction-aware basin from an outlet flat index,
    // a parent flow-routing object, and a junction network defined
    // on the same flow domain.
    //
    // Both flow and junction_network must outlive this basin.
    lsdtt_junction_basin_cpu(
        std::int64_t                        basin_outlet_flat_index,
        const lsdtt_flow_base&              flow,
        const lsdtt_junction_network_base&  junction_network
    );

    // Returns the linked junction number for this basin.
    // Returns -1 if no junction could be associated.
    [[nodiscard]] std::int32_t get_junction_number() const noexcept;

    // Returns the stream order associated with the linked junction.
    // Returns -1 if no valid junction is associated.
    [[nodiscard]] std::int32_t get_stream_order() const noexcept;

    // Returns the non-owning linked junction network.
    [[nodiscard]] const lsdtt_junction_network_base&
        get_junction_network() const noexcept;

private:

    [[nodiscard]] static std::int32_t resolve_junction_number(
        const lsdtt_junction_network_base& junction_network,
        std::int64_t basin_outlet_flat_index);

    [[nodiscard]] static std::int32_t resolve_stream_order(
        const lsdtt_junction_network_base& junction_network,
        std::int32_t junction_number);

private:

    std::int32_t junction_number_;
    std::int32_t stream_order_;
    const lsdtt_junction_network_base* junction_network_;   // non-owning; outlives basin
};

} // namespace lsdtt
