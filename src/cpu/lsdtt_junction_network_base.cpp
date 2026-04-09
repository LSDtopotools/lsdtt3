// lsdtt_junction_network_base.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Factory implementation for make_junction_network_base().
//
// Declared in src/core/lsdtt_junction_network_base.hpp.
// Lives in src/cpu/ so that core does not depend on cpu — the same pattern
// used by make_flow_base() in src/cpu/lsdtt_flow_base.cpp.
// -----------------------------------------------------------------------------

#include "../core/lsdtt_junction_network_base.hpp"
#include "lsdtt_junction_network_cpu.hpp"

namespace lsdtt
{

std::unique_ptr<lsdtt_junction_network_base> make_junction_network_base(
    const lsdtt_flow_base&                    flow,
    const std::vector<std::int64_t>&          source_flats,
    const lsdtt_junction_network_options&     options
)
{
    auto network = std::make_unique<lsdtt_junction_network_cpu>();
    network->create(flow, source_flats, options);
    return network;
}

} // namespace lsdtt
