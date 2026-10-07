// lsdtt_junction_graph_orders.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Traversal orders on a compact junction graph: a forest of int32 junction ids
// 0..n-1 in which every junction has at most one receiver (-1 at an outlet).
//
// junction_graph_topological_order
//   Kahn's algorithm, donors before receivers. Junctions with no donors are
//   queued in ascending id, then each receiver is appended when its last donor
//   has been taken. remaining_donors holds each junction's donor count (taken
//   by value and consumed).
//
// junction_graph_preorder
//   Depth-first preorder from each outlet in the order given, visiting each
//   junction's donors in the order donors_of(j) returns them. Produces the
//   s_vector plus, for each junction, the start and size of its contiguous
//   subtree interval (upstream junctions occupy s_vector[start, start + size)).
//   start is -1 and size 0 for a junction that is never reached.
//
// Both report a failure as a status rather than throwing, so each caller keeps
// its own error message. Shared by lsdtt_junction_network_cpu and
// lsdtt_vector_to_junction_network, which had the same loops
// (tests/cpp/test_junction_graph_orders.cpp keeps them as the reference). The
// D8 pixel graph in lsdtt_flowinfo and the multi-receiver graph in
// lsdtt_flowgraph_cpu have their own, deliberately different, versions.
//
// Header-only: integer graph traversal, no floating arithmetic.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace lsdtt
{

enum class junction_graph_order_status
{
    ok,
    cyclic_or_incomplete,   // some junction was never released (a cycle, or a broken graph)
    donor_count_underflow,  // a receiver was released more times than it has donors
    duplicate_visit,        // preorder reached a junction twice
    incomplete_cover        // preorder did not reach every junction
};

// receiver_of(j) -> std::int32_t, -1 at an outlet.
template <class ReceiverOf>
junction_graph_order_status junction_graph_topological_order(
    std::vector<std::int32_t>  remaining_donors,
    ReceiverOf                 receiver_of,
    std::vector<std::int32_t>& order)
{
    const std::size_t n = remaining_donors.size();
    order.clear();
    order.reserve(n);

    std::vector<std::int32_t> ready;
    ready.reserve(n);

    for (std::int32_t j = 0; j < static_cast<std::int32_t>(n); ++j)
    {
        if (remaining_donors[static_cast<std::size_t>(j)] == 0)
        {
            ready.push_back(j);
        }
    }

    std::size_t head = 0;
    while (head < ready.size())
    {
        const std::int32_t j = ready[head++];
        order.push_back(j);

        const std::int32_t rec = receiver_of(j);
        if (rec >= 0)
        {
            if (remaining_donors[static_cast<std::size_t>(rec)] <= 0)
            {
                return junction_graph_order_status::donor_count_underflow;
            }

            --remaining_donors[static_cast<std::size_t>(rec)];

            if (remaining_donors[static_cast<std::size_t>(rec)] == 0)
            {
                ready.push_back(rec);
            }
        }
    }

    if (order.size() != n)
    {
        return junction_graph_order_status::cyclic_or_incomplete;
    }

    return junction_graph_order_status::ok;
}

// donors_of(j) -> a container of std::int32_t in visiting order.
template <class DonorsOf>
junction_graph_order_status junction_graph_preorder(
    std::size_t                      n_junctions,
    const std::vector<std::int32_t>& outlets,
    DonorsOf                         donors_of,
    std::vector<std::int32_t>&       s_vector,
    std::vector<std::int32_t>&       start,
    std::vector<std::int32_t>&       size)
{
    s_vector.clear();
    s_vector.reserve(n_junctions);
    start.assign(n_junctions, -1);
    size.assign(n_junctions, 0);

    struct dfs_frame
    {
        std::int32_t junction   = -1;
        std::size_t  next_donor = 0;
    };

    std::vector<dfs_frame> stack;
    stack.reserve(256);

    auto push_junction = [&](std::int32_t junction)
    {
        start[static_cast<std::size_t>(junction)] =
            static_cast<std::int32_t>(s_vector.size());
        size[static_cast<std::size_t>(junction)] = 1;
        s_vector.push_back(junction);
        stack.push_back(dfs_frame{junction, 0});
    };

    for (const std::int32_t outlet : outlets)
    {
        if (start[static_cast<std::size_t>(outlet)] >= 0)
        {
            continue;
        }

        stack.clear();
        push_junction(outlet);

        while (!stack.empty())
        {
            // Copied, not referenced: push_junction can reallocate the stack.
            const std::int32_t junction = stack.back().junction;
            const auto& donor_list = donors_of(junction);

            if (stack.back().next_donor < donor_list.size())
            {
                const std::int32_t child = donor_list[stack.back().next_donor++];

                if (start[static_cast<std::size_t>(child)] >= 0)
                {
                    return junction_graph_order_status::duplicate_visit;
                }

                push_junction(child);
                continue;
            }

            const std::int32_t completed_size = size[static_cast<std::size_t>(junction)];

            stack.pop_back();

            if (!stack.empty())
            {
                size[static_cast<std::size_t>(stack.back().junction)] += completed_size;
            }
        }
    }

    if (s_vector.size() != n_junctions)
    {
        return junction_graph_order_status::incomplete_cover;
    }

    return junction_graph_order_status::ok;
}

} // namespace lsdtt
