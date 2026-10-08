// test_junction_graph_orders.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for src/utils/lsdtt_junction_graph_orders.hpp
// (consolidation C05).
//
// lsdtt_junction_network_cpu and lsdtt_vector_to_junction_network each had a
// Kahn topological order and a depth-first s_vector builder over their junction
// graph. The four removed member functions are kept below with their bodies
// verbatim; only the member variables they used are passed in as parameters.
// Each is compared with the shared routine called exactly as the migrated
// member now calls it (written below), on:
//
//   - random forests, n = 0..300, with several outlets, outlets listed in
//     ascending and in shuffled order, donors ascending (as both classes sort
//     them) and in arbitrary order;
//   - corrupt graphs: a cycle, a donor count one too low (underflow), a donor
//     listed twice (duplicate visit), and a missing outlet (incomplete cover).
//
// Compared: the order, s_vector, start and size arrays exactly, and the
// exception message exactly. One deliberate difference: the vector network had
// no donor-underflow check. With a donor list shorter than the receivers imply
// it either reported "cyclic or incomplete" or returned an order without error,
// and that order was not topological whenever the unlisted donor was taken after
// its receiver (115 of 160 such cases in this test). It now throws
// the underflow message the raster network already had (the plan's option to
// add the check to both). Only an internally inconsistent graph can reach it.
//
// Usage: lsdtt-test-junction-graph-orders
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_types.hpp"
#include "../../src/utils/lsdtt_junction_graph_orders.hpp"

using namespace lsdtt;

namespace
{

// ============================================================================
// Removed copies (references); member variables passed as parameters
// ============================================================================

// src/cpu/lsdtt_junction_network_cpu.cpp:1605-1668 (as of 2026-10-07)
void jn_build_topological_order(const std::vector<std::int64_t>& junction_flat_,
                                const std::vector<std::int32_t>& donor_count_junction_,
                                const std::vector<std::int32_t>& receiver_junction_,
                                std::vector<std::int32_t>& topological_order_)
{
    topological_order_.clear();
    topological_order_.reserve(junction_flat_.size());

    const std::int32_t n_junctions =
        static_cast<std::int32_t>(junction_flat_.size());

    if (n_junctions == 0)
    {
        return;
    }

    std::vector<std::int32_t> remaining_donors(
        donor_count_junction_.begin(),
        donor_count_junction_.end()
    );

    std::vector<std::int32_t> ready_queue;
    ready_queue.reserve(static_cast<std::size_t>(n_junctions));

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        if (remaining_donors[static_cast<std::size_t>(j)] == 0)
        {
            ready_queue.push_back(j);
        }
    }

    std::size_t head = 0;
    while (head < ready_queue.size())
    {
        const std::int32_t j = ready_queue[head++];
        topological_order_.push_back(j);

        const std::int32_t rec = receiver_junction_[static_cast<std::size_t>(j)];
        if (rec >= 0)
        {
            if (remaining_donors[static_cast<std::size_t>(rec)] <= 0)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: donor counter underflow while building topological order"
                );
            }

            --remaining_donors[static_cast<std::size_t>(rec)];

            if (remaining_donors[static_cast<std::size_t>(rec)] == 0)
            {
                ready_queue.push_back(rec);
            }
        }
    }

    if (topological_order_.size() != junction_flat_.size())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: compact junction graph is cyclic or incomplete"
        );
    }
}

// src/cpu/lsdtt_junction_network_cpu.cpp:1674-1784 (as of 2026-10-07)
void jn_build_s_vector(const std::vector<std::int64_t>& junction_flat_,
                       const std::vector<std::int32_t>& receiver_junction_,
                       const std::vector<std::int32_t>& outlets_,
                       std::vector<std::int32_t>& s_vector_,
                       std::vector<std::int32_t>& s_vector_start_,
                       std::vector<std::int32_t>& s_vector_size_)
{
    const std::int32_t n_junctions =
        static_cast<std::int32_t>(junction_flat_.size());

    s_vector_.clear();
    s_vector_start_.assign(static_cast<std::size_t>(n_junctions), k_nodata_junction);
    s_vector_size_.assign(static_cast<std::size_t>(n_junctions), 0);

    if (n_junctions == 0)
    {
        return;
    }

    std::vector<std::vector<std::int32_t>> donors(
        static_cast<std::size_t>(n_junctions)
    );

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        const std::int32_t rec = receiver_junction_[static_cast<std::size_t>(j)];
        if (rec >= 0)
        {
            donors[static_cast<std::size_t>(rec)].push_back(j);
        }
    }

    // Sort donors by ascending junction id.
    // Because junction ids are assigned in row-major order of their flat pixel,
    // this is also ascending row-major order of junction location.
    for (auto& donor_list : donors)
    {
        std::sort(donor_list.begin(), donor_list.end());
    }

    struct dfs_frame
    {
        std::int32_t junction   = k_nodata_junction;
        std::size_t  next_donor = 0;
    };

    std::vector<dfs_frame> stack;
    stack.reserve(256);

    auto push_junction = [&](std::int32_t junction)
    {
        s_vector_start_[static_cast<std::size_t>(junction)] =
            static_cast<std::int32_t>(s_vector_.size());

        s_vector_size_[static_cast<std::size_t>(junction)] = 1;
        s_vector_.push_back(junction);

        stack.push_back(dfs_frame{junction, 0});
    };

    for (const std::int32_t outlet : outlets_)
    {
        if (s_vector_start_[static_cast<std::size_t>(outlet)] >= 0)
        {
            continue;
        }

        stack.clear();
        push_junction(outlet);

        while (!stack.empty())
        {
            dfs_frame& frame = stack.back();
            const std::int32_t junction = frame.junction;
            const auto& donor_list = donors[static_cast<std::size_t>(junction)];

            if (frame.next_donor < donor_list.size())
            {
                const std::int32_t child =
                    donor_list[frame.next_donor++];

                if (s_vector_start_[static_cast<std::size_t>(child)] >= 0)
                {
                    throw std::runtime_error(
                        "lsdtt_junction_network_cpu: duplicate visit while building s_vector"
                    );
                }

                push_junction(child);
                continue;
            }

            const std::int32_t completed_junction = junction;
            const std::int32_t completed_size =
                s_vector_size_[static_cast<std::size_t>(completed_junction)];

            stack.pop_back();

            if (!stack.empty())
            {
                s_vector_size_[static_cast<std::size_t>(stack.back().junction)] +=
                    completed_size;
            }
        }
    }

    if (s_vector_.size() != static_cast<std::size_t>(n_junctions))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: s_vector does not cover all junctions"
        );
    }
}

// The vector network's per-junction record, reduced to the fields used.
struct vjunction
{
    std::int32_t receiver_junction = k_nodata_junction;
    std::vector<std::int32_t> donor_junctions;
    std::int32_t s_vector_start = -1;
    std::int32_t s_vector_size  = 0;
};

// src/cpu/lsdtt_vector_to_junction_network.cpp:1792-1845 (as of 2026-10-07)
void vn_build_topological_order(const std::vector<vjunction>& junctions_,
                                std::vector<std::int32_t>& topological_order_)
{
    topological_order_.clear();
    topological_order_.reserve(junctions_.size());

    std::vector<std::int32_t> remaining(
        junctions_.size(),
        0
    );

    for (std::size_t i = 0; i < junctions_.size(); ++i)
    {
        remaining[i] = static_cast<std::int32_t>(
            junctions_[i].donor_junctions.size()
        );
    }

    std::vector<std::int32_t> ready;
    ready.reserve(junctions_.size());

    for (std::int32_t j = 0; j < static_cast<std::int32_t>(junctions_.size()); ++j)
    {
        if (remaining[static_cast<std::size_t>(j)] == 0)
        {
            ready.push_back(j);
        }
    }

    std::size_t head = 0;
    while (head < ready.size())
    {
        const std::int32_t j = ready[head++];
        topological_order_.push_back(j);

        const std::int32_t rec = junctions_[static_cast<std::size_t>(j)].receiver_junction;
        if (rec >= 0)
        {
            --remaining[static_cast<std::size_t>(rec)];
            if (remaining[static_cast<std::size_t>(rec)] == 0)
            {
                ready.push_back(rec);
            }
        }
    }

    if (topological_order_.size() != junctions_.size())
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: vector junction graph is cyclic or incomplete"
        );
    }
}

// src/cpu/lsdtt_vector_to_junction_network.cpp:1936-2013 (as of 2026-10-07)
void vn_build_s_vector(std::vector<vjunction>& junctions_,
                       const std::vector<std::int32_t>& outlets_,
                       std::vector<std::int32_t>& s_vector_)
{
    s_vector_.clear();
    s_vector_.reserve(junctions_.size());

    struct dfs_frame
    {
        std::int32_t junction = k_nodata_junction;
        std::size_t next_donor = 0;
    };

    std::vector<dfs_frame> stack;
    stack.reserve(128);

    auto push_junction =
        [&](std::int32_t junction)
        {
            junctions_[static_cast<std::size_t>(junction)].s_vector_start =
                static_cast<std::int32_t>(s_vector_.size());

            junctions_[static_cast<std::size_t>(junction)].s_vector_size = 1;
            s_vector_.push_back(junction);
            stack.push_back({junction, 0});
        };

    for (const std::int32_t outlet : outlets_)
    {
        if (junctions_[static_cast<std::size_t>(outlet)].s_vector_start >= 0)
        {
            continue;
        }

        stack.clear();
        push_junction(outlet);

        while (!stack.empty())
        {
            dfs_frame& frame = stack.back();
            auto& donors =
                junctions_[static_cast<std::size_t>(frame.junction)].donor_junctions;

            if (frame.next_donor < donors.size())
            {
                const std::int32_t child = donors[frame.next_donor++];
                if (junctions_[static_cast<std::size_t>(child)].s_vector_start >= 0)
                {
                    throw std::runtime_error(
                        "lsdtt_vector_to_junction_network: duplicate visit while building s_vector"
                    );
                }

                push_junction(child);
                continue;
            }

            const std::int32_t done = frame.junction;
            const std::int32_t done_size =
                junctions_[static_cast<std::size_t>(done)].s_vector_size;

            stack.pop_back();

            if (!stack.empty())
            {
                junctions_[static_cast<std::size_t>(stack.back().junction)].s_vector_size +=
                    done_size;
            }
        }
    }

    if (s_vector_.size() != junctions_.size())
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: s_vector does not cover all junctions"
        );
    }
}

// ============================================================================
// The migrated members, written exactly as they now appear in the sources
// ============================================================================

void jn_new_topological_order(const std::vector<std::int32_t>& donor_count_junction_,
                              const std::vector<std::int32_t>& receiver_junction_,
                              std::vector<std::int32_t>& topological_order_)
{
    const junction_graph_order_status status = junction_graph_topological_order(
        donor_count_junction_,
        [&](std::int32_t j) { return receiver_junction_[static_cast<std::size_t>(j)]; },
        topological_order_);

    if (status == junction_graph_order_status::donor_count_underflow)
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: donor counter underflow while building topological order");
    if (status != junction_graph_order_status::ok)
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: compact junction graph is cyclic or incomplete");
}

void jn_new_s_vector(std::size_t n_junctions,
                     const std::vector<std::int32_t>& receiver_junction_,
                     const std::vector<std::int32_t>& outlets_,
                     std::vector<std::int32_t>& s_vector_,
                     std::vector<std::int32_t>& s_vector_start_,
                     std::vector<std::int32_t>& s_vector_size_)
{
    std::vector<std::vector<std::int32_t>> donors(n_junctions);
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(n_junctions); ++j)
    {
        const std::int32_t rec = receiver_junction_[static_cast<std::size_t>(j)];
        if (rec >= 0) donors[static_cast<std::size_t>(rec)].push_back(j);
    }
    for (auto& donor_list : donors) std::sort(donor_list.begin(), donor_list.end());

    const junction_graph_order_status status = junction_graph_preorder(
        n_junctions, outlets_,
        [&](std::int32_t j) -> const std::vector<std::int32_t>& {
            return donors[static_cast<std::size_t>(j)]; },
        s_vector_, s_vector_start_, s_vector_size_);

    if (status == junction_graph_order_status::duplicate_visit)
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: duplicate visit while building s_vector");
    if (status != junction_graph_order_status::ok)
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: s_vector does not cover all junctions");
}

void vn_new_topological_order(const std::vector<vjunction>& junctions_,
                              std::vector<std::int32_t>& topological_order_)
{
    std::vector<std::int32_t> donor_counts(junctions_.size());
    for (std::size_t i = 0; i < junctions_.size(); ++i)
        donor_counts[i] = static_cast<std::int32_t>(junctions_[i].donor_junctions.size());

    const junction_graph_order_status status = junction_graph_topological_order(
        std::move(donor_counts),
        [&](std::int32_t j) { return junctions_[static_cast<std::size_t>(j)].receiver_junction; },
        topological_order_);

    if (status == junction_graph_order_status::donor_count_underflow)
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: donor counter underflow while building topological order");
    if (status != junction_graph_order_status::ok)
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: vector junction graph is cyclic or incomplete");
}

void vn_new_s_vector(std::vector<vjunction>& junctions_,
                     const std::vector<std::int32_t>& outlets_,
                     std::vector<std::int32_t>& s_vector_)
{
    std::vector<std::int32_t> start, size;
    const junction_graph_order_status status = junction_graph_preorder(
        junctions_.size(), outlets_,
        [&](std::int32_t j) -> const std::vector<std::int32_t>& {
            return junctions_[static_cast<std::size_t>(j)].donor_junctions; },
        s_vector_, start, size);

    for (std::size_t i = 0; i < junctions_.size(); ++i)
    {
        junctions_[i].s_vector_start = start[i];
        junctions_[i].s_vector_size  = size[i];
    }

    if (status == junction_graph_order_status::duplicate_visit)
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: duplicate visit while building s_vector");
    if (status != junction_graph_order_status::ok)
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: s_vector does not cover all junctions");
}

// ============================================================================
// Harness
// ============================================================================

int g_failures = 0;

void check(const std::string& what, bool ok, long long n_cases, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << " (" << n_cases << " cases)";
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

// True when every junction comes before its receiver in order.
bool is_topological(const std::vector<std::int32_t>& order, const std::vector<std::int32_t>& receiver)
{
    std::vector<std::size_t> pos(receiver.size());
    for (std::size_t k = 0; k < order.size(); ++k) pos[static_cast<std::size_t>(order[k])] = k;
    for (std::size_t j = 0; j < receiver.size(); ++j)
        if (receiver[j] >= 0 && pos[j] > pos[static_cast<std::size_t>(receiver[j])]) return false;
    return true;
}

template <class F>
std::string outcome(F f)
{
    try { f(); return ""; }
    catch (const std::exception& e) { return e.what(); }
}

struct forest
{
    std::vector<std::int32_t> receiver;       // -1 at an outlet
    std::vector<std::int32_t> donor_count;
    std::vector<std::vector<std::int32_t>> donors;  // in vector-network order
    std::vector<std::int32_t> outlets;
};

// corruption: 0 none, 1 cycle, 2 donor count too low, 3 donor listed twice,
// 4 outlet missing from the outlet list
forest make_forest(std::size_t n, std::mt19937& rng, bool shuffle_outlets, bool sort_donors,
                   int corruption)
{
    forest f;
    f.receiver.assign(n, -1);
    // Receivers point to a lower id with probability ~0.85, so ids are a valid
    // forest; a random permutation then scrambles the id order.
    std::vector<std::int32_t> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (std::size_t k = 1; k < n; ++k)
        if (u(rng) < 0.85)
        {
            std::uniform_int_distribution<std::size_t> parent(0, k - 1);
            f.receiver[static_cast<std::size_t>(perm[k])] = perm[parent(rng)];
        }
    if (corruption == 1 && n >= 3)
    {
        // close a loop between two junctions that are not outlets
        f.receiver[static_cast<std::size_t>(perm[1])] = perm[2];
        f.receiver[static_cast<std::size_t>(perm[2])] = perm[1];
    }
    f.donors.assign(n, {});
    for (std::size_t j = 0; j < n; ++j)
        if (f.receiver[j] >= 0) f.donors[static_cast<std::size_t>(f.receiver[j])].push_back(static_cast<std::int32_t>(j));
    for (auto& d : f.donors)
        if (sort_donors) std::sort(d.begin(), d.end());
        else std::shuffle(d.begin(), d.end(), rng);
    if (corruption == 3)
        for (auto& d : f.donors)
            if (!d.empty()) { d.push_back(d.front()); break; }
    f.donor_count.assign(n, 0);
    for (std::size_t j = 0; j < n; ++j) f.donor_count[j] = static_cast<std::int32_t>(f.donors[j].size());
    if (corruption == 2)
        for (std::size_t j = 0; j < n; ++j)
            if (f.donor_count[j] > 0) { --f.donor_count[j]; f.donors[j].pop_back(); break; }
    for (std::size_t j = 0; j < n; ++j)
        if (f.receiver[j] < 0) f.outlets.push_back(static_cast<std::int32_t>(j));
    if (shuffle_outlets) std::shuffle(f.outlets.begin(), f.outlets.end(), rng);
    if (corruption == 4 && f.outlets.size() >= 2) f.outlets.pop_back();
    return f;
}

} // namespace

int main()
{
    std::cout << "== C05 junction-graph orders vs the removed member functions ==\n";
    std::mt19937 rng(2026);

    for (int corruption = 0; corruption <= 4; ++corruption)
    {
        long long n_cases = 0, bad_jn_k = 0, bad_jn_s = 0, bad_vn_k = 0, bad_vn_s = 0, n_threw = 0;
        long long n_vn_underflow_now_reported = 0, n_vn_old_order_invalid = 0;
        std::string first;
        for (int rep = 0; rep < 160; ++rep)
        {
            const std::size_t n = (corruption == 0) ? static_cast<std::size_t>(rep * 2 % 301)
                                                    : 3 + static_cast<std::size_t>(rep % 120);
            const bool shuffle_outlets = (rep % 2 == 1);
            const bool sort_donors = (rep % 3 != 2);
            const forest f = make_forest(n, rng, shuffle_outlets, sort_donors, corruption);
            ++n_cases;

            // ---- raster junction network: it builds its own (sorted) donors ----
            const std::vector<std::int64_t> flat(n, 0);
            std::vector<std::int32_t> o_old, o_new;
            const std::string e_old = outcome([&] { jn_build_topological_order(flat, f.donor_count, f.receiver, o_old); });
            const std::string e_new = outcome([&] { jn_new_topological_order(f.donor_count, f.receiver, o_new); });
            if (e_old != e_new || (e_old.empty() && o_old != o_new))
            {
                ++bad_jn_k;
                if (first.empty()) first = "jn kahn n=" + std::to_string(n) + " old{" + e_old + "} new{" + e_new + "}";
            }
            if (!e_old.empty()) ++n_threw;

            std::vector<std::int32_t> s1, st1, sz1, s2, st2, sz2;
            const std::string d_old = outcome([&] { jn_build_s_vector(flat, f.receiver, f.outlets, s1, st1, sz1); });
            const std::string d_new = outcome([&] { jn_new_s_vector(n, f.receiver, f.outlets, s2, st2, sz2); });
            if (d_old != d_new || (d_old.empty() && (s1 != s2 || st1 != st2 || sz1 != sz2)))
            {
                ++bad_jn_s;
                if (first.empty()) first = "jn dfs n=" + std::to_string(n) + " old{" + d_old + "} new{" + d_new + "}";
            }

            // ---- vector network: donors as stored ----
            std::vector<vjunction> va(n), vb(n);
            for (std::size_t j = 0; j < n; ++j)
            {
                va[j].receiver_junction = vb[j].receiver_junction = f.receiver[j];
                va[j].donor_junctions = vb[j].donor_junctions = f.donors[j];
            }
            std::vector<std::int32_t> vo_old, vo_new;
            const std::string ve_old = outcome([&] { vn_build_topological_order(va, vo_old); });
            const std::string ve_new = outcome([&] { vn_new_topological_order(vb, vo_new); });
            // The one intended difference: the new underflow check. Without it
            // the old code either reported "cyclic" or, worse, returned an
            // order that is not topological (a receiver released before its
            // last donor). Require that the old result really was one of those.
            const bool underflow_case =
                ve_new == "lsdtt_vector_to_junction_network: donor counter underflow while building topological order" &&
                (ve_old == "lsdtt_vector_to_junction_network: vector junction graph is cyclic or incomplete" ||
                 ve_old.empty());
            if (underflow_case)
            {
                ++n_vn_underflow_now_reported;
                if (ve_old.empty() && !is_topological(vo_old, f.receiver)) ++n_vn_old_order_invalid;
            }
            else if (ve_old != ve_new || (ve_old.empty() && vo_old != vo_new))
            {
                ++bad_vn_k;
                if (first.empty()) first = "vn kahn n=" + std::to_string(n) + " old{" + ve_old + "} new{" + ve_new + "}";
            }

            std::vector<std::int32_t> vs_old, vs_new;
            const std::string vd_old = outcome([&] { vn_build_s_vector(va, f.outlets, vs_old); });
            const std::string vd_new = outcome([&] { vn_new_s_vector(vb, f.outlets, vs_new); });
            bool same_records = true;
            if (vd_old.empty())
                for (std::size_t j = 0; j < n; ++j)
                    same_records = same_records && va[j].s_vector_start == vb[j].s_vector_start &&
                                   va[j].s_vector_size == vb[j].s_vector_size;
            if (vd_old != vd_new || (vd_old.empty() && (vs_old != vs_new || !same_records)))
            {
                ++bad_vn_s;
                if (first.empty()) first = "vn dfs n=" + std::to_string(n) + " old{" + vd_old + "} new{" + vd_new + "}";
            }
        }
        const char* names[] = {"valid forests", "cycle", "donor count too low", "donor listed twice",
                               "outlet missing"};
        const std::string tag = std::string(" [") + names[corruption] + "]";
        check("raster network topological order" + tag + "; " + std::to_string(n_threw) + " threw",
              bad_jn_k == 0, n_cases, first);
        check("raster network s_vector, start, size" + tag, bad_jn_s == 0, n_cases, first);
        check("vector network topological order" + tag +
                  (n_vn_underflow_now_reported ? "; " + std::to_string(n_vn_underflow_now_reported) +
                                                     " now report underflow; the old code returned a non-topological order in " +
                                                     std::to_string(n_vn_old_order_invalid) + " of them"
                                               : std::string()),
              bad_vn_k == 0, n_cases, first);
        check("vector network s_vector and per-junction start, size" + tag, bad_vn_s == 0, n_cases, first);
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
