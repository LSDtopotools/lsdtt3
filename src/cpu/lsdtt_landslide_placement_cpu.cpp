// lsdtt_landslide_placement_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// See lsdtt_landslide_placement_cpu.hpp for the design.

#include "lsdtt_landslide_placement_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace lsdtt
{

namespace
{

// xorshift64*, so a test can reproduce a placement exactly without dragging in
// <random> and its implementation-defined engines.
inline double next_uniform(std::uint64_t& state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    const std::uint64_t v = state * 2685821657736338717ULL;
    // 53 significant bits, in [0,1)
    return static_cast<double>(v >> 11) * (1.0 / 9007199254740992.0);
}

}  // namespace


LandslidePlacer::LandslidePlacer(const lsdtt_flow_base& flow,
                                 const lsdtt_raster_meta& meta,
                                 const LandslidePlacementConfig& config,
                                 const lsdtt_raster_cpu* mask)
  : flow_(flow), meta_(meta), cfg_(config)
{
    const std::int64_t n_pixels = meta_.n_pixels();
    anneal_blocked_.assign(static_cast<std::size_t>(n_pixels), 0);
    event_blocked_.assign(static_cast<std::size_t>(n_pixels), 0);

    const bool have_mask = (mask != nullptr);
    const float mask_ndv = have_mask
        ? static_cast<float>(mask->meta().no_data_value)
        : 0.0f;

    hosts_.reserve(static_cast<std::size_t>(n_pixels) / 4);

    for (std::int64_t flat = 0; flat < n_pixels; ++flat)
    {
        if (!flow_.is_valid_cell(flat)) continue;

        // Hillslope only: a landslide may not start in a channel. Because
        // contributing area decreases upslope, this single test also
        // guarantees every pixel in the upslope set is a hillslope pixel.
        const std::int64_t n_up =
            static_cast<std::int64_t>(flow_.s_vector_size_unchecked(flat));
        if (n_up <= 0) continue;
        if (n_up >= cfg_.channel_threshold_pixels) continue;

        if (have_mask)
        {
            // The mask must hold over the WHOLE landslide, not just its seed.
            // A landslide occupies the host's upslope set, and testing only the
            // host lets a scar seeded just inside the mask spill across the
            // boundary -- measured at 40 of 1941 excavated pixels on a
            // step-change tag mask, 35 of them on ground UPSTREAM of the
            // knickpoint, which is the exact population the mask existed to
            // exclude. A mask that leaks is not a mask.
            //
            // Rejecting the host (rather than trimming the scar) is deliberate:
            // trimming would silently shrink landslides near the boundary and
            // so distort the imposed area distribution, which is the one thing
            // this module guarantees. The cost is that hosts within a
            // landslide's reach of the edge are unusable, which is the honest
            // consequence of demanding containment.
            const auto& mmask = mask->get_mask_array();
            const auto& mdata = mask->get_data_array();
            auto inside = [&](std::int64_t f) -> bool
            {
                const int rr = static_cast<int>(f / meta_.n_cols);
                const int cc = static_cast<int>(f % meta_.n_cols);
                if (!mmask(rr, cc)) return false;
                const float mv = mdata(rr, cc);
                return !(mv == mask_ndv || mv == 0.0f);
            };

            if (!inside(flat)) continue;
            bool contained = true;
            for (std::int64_t up : flow_.upslope_nodes_unchecked(flat))
            {
                if (!inside(up)) { contained = false; break; }
            }
            if (!contained) continue;
        }

        hosts_.push_back(Host{flat, n_up});
        max_hostable_ = std::max(max_hostable_, n_up);
    }

    // Descending by upslope count, so hosts able to take an n-pixel landslide
    // form a prefix and can be located by binary search.
    std::sort(hosts_.begin(), hosts_.end(),
              [](const Host& a, const Host& b) { return a.n_up > b.n_up; });

    diag_.max_hostable_pixels = max_hostable_;
}


void LandslidePlacer::apply_annealing(const lsdtt_raster_cpu& time_of_last_failure,
                                      double t_now_yr,
                                      double anneal_time_yr)
{
    // Rebuild from scratch: annealing state is a pure function of the raster
    // and the current time, so a stale entry from an earlier call must not
    // survive.
    std::fill(anneal_blocked_.begin(), anneal_blocked_.end(), 0);

    const auto& data = time_of_last_failure.get_data_array();
    const auto& msk  = time_of_last_failure.get_mask_array();
    const float ndv  = static_cast<float>(time_of_last_failure.meta().no_data_value);

    for (int r = 0; r < meta_.n_rows; ++r)
    {
        for (int c = 0; c < meta_.n_cols; ++c)
        {
            const float t_fail = data(r, c);
            // Never failed: nodata, or the caller's negative sentinel.
            if (!msk(r, c) || t_fail == ndv) continue;

            const double elapsed = t_now_yr - static_cast<double>(t_fail);
            if (elapsed < anneal_time_yr)
            {
                const std::int64_t flat =
                    static_cast<std::int64_t>(r) * meta_.n_cols + c;
                anneal_blocked_[static_cast<std::size_t>(flat)] = 1;
            }
        }
    }
}


void LandslidePlacer::begin_event()
{
    // Clear only what THIS event excavated. Annealing state is left alone: it
    // is owned by apply_annealing, which the caller invokes with the event time.
    std::fill(event_blocked_.begin(), event_blocked_.end(), 0);

    diag_ = LandslidePlacementDiagnostics{};
    diag_.max_hostable_pixels = max_hostable_;
    fill_ratios_.clear();
}


std::size_t LandslidePlacer::n_eligible_hosts(int n_target_pixels) const
{
    // hosts_ is sorted descending by n_up, so the eligible prefix ends at the
    // first host with n_up < n_target.
    const auto it = std::lower_bound(
        hosts_.begin(), hosts_.end(), static_cast<std::int64_t>(n_target_pixels),
        [](const Host& h, std::int64_t target) { return h.n_up > target - 1; });
    // The comparator above finds the first element for which !(n_up >= target).
    return static_cast<std::size_t>(it - hosts_.begin());
}


bool LandslidePlacer::try_host(std::int64_t init_flat, int n_target,
                               LandslideScar& scar)
{
    const lsdtt_flat_index_slice slice = flow_.upslope_nodes(init_flat);
    const std::size_t n_up = slice.size();
    if (n_up < static_cast<std::size_t>(n_target)) return false;

    // Flow distance from the initiation pixel, in one O(n_up) pass.
    //
    // The s_vector is a preorder DFS of the drainage tree rooted here, so a
    // node's receiver always appears earlier in the slice than the node itself.
    // s_vector_start(v) is v's index in the GLOBAL s_vector, so v's position
    // within this slice is s_vector_start(v) - s_vector_start(init) -- an O(1)
    // lookup, which is what makes the exact flow distance affordable.
    const std::int64_t start_init = flow_.s_vector_start_unchecked(init_flat);
    const double cell = meta_.cell_size;
    const double diag = cell * grid::k_sqrt2_d;

    // The initiation pixel itself must survive, or the scar is not anchored
    // where we said it was.
    if (cfg_.forbid_overlap && !is_available(init_flat))
        return false;

    scratch_dist_.assign(n_up, 0.0);
    scratch_rec_lp_.assign(n_up, -1);

    for (std::size_t k = 1; k < n_up; ++k)   // k = 0 is init itself, distance 0
    {
        const std::int64_t v = slice.data()[k];
        const std::int64_t rec = flow_.receiver_of_unchecked(v);
        if (rec < 0) { scratch_dist_[k] = std::numeric_limits<double>::max(); continue; }

        const std::int64_t lp = flow_.s_vector_start_unchecked(rec) - start_init;
        // Guard: a receiver outside this slice would mean the tree structure is
        // not what the contract promises. Fall back to a large distance so such
        // a node is trimmed first rather than corrupting the scar.
        if (lp < 0 || static_cast<std::size_t>(lp) >= n_up)
        {
            scratch_dist_[k] = std::numeric_limits<double>::max();
            continue;
        }
        scratch_rec_lp_[k] = static_cast<std::int64_t>(lp);

        const int r_v = static_cast<int>(v / meta_.n_cols);
        const int c_v = static_cast<int>(v % meta_.n_cols);
        const int r_r = static_cast<int>(rec / meta_.n_cols);
        const int c_r = static_cast<int>(rec % meta_.n_cols);
        const bool is_diag = (r_v != r_r) && (c_v != c_r);

        scratch_dist_[k] = scratch_dist_[static_cast<std::size_t>(lp)]
                         + (is_diag ? diag : cell);
    }

    // Order every node by flow distance, nearest first.
    scratch_order_.resize(n_up);
    std::iota(scratch_order_.begin(), scratch_order_.end(), 0u);
    std::sort(scratch_order_.begin(), scratch_order_.end(),
              [this](std::uint32_t a, std::uint32_t b) {
                  return scratch_dist_[a] < scratch_dist_[b];
              });

    // Accept in that order, but ONLY if the node's receiver has already been
    // accepted. This is what keeps the scar connected.
    //
    // Distance alone is not enough. dist(v) = dist(receiver(v)) + length, so on
    // a clean landscape the n nearest nodes are automatically closed under
    // "receiver" and the scar is a connected subtree for free. But blocked
    // pixels break that: if a node's receiver was excavated by an earlier
    // landslide in this event and is therefore unavailable, taking the node
    // anyway would leave a scar detached from its own initiation pixel.
    // Requiring the receiver first skips the whole subtree beyond a blockage,
    // which is the correct behaviour -- material cannot detach across ground
    // that is not part of the failure.
    scratch_accepted_.assign(n_up, 0);
    scratch_accepted_[0] = 1;                      // the initiation pixel
    std::size_t n_accepted = 1;

    scar.pixels.clear();
    scar.pixels.reserve(static_cast<std::size_t>(n_target));
    scar.pixels.push_back(init_flat);

    for (std::uint32_t k : scratch_order_)
    {
        if (n_accepted >= static_cast<std::size_t>(n_target)) break;
        if (k == 0) continue;                      // already taken

        const std::int64_t lp = scratch_rec_lp_[k];
        if (lp < 0) continue;                      // no usable receiver
        if (!scratch_accepted_[static_cast<std::size_t>(lp)]) continue;

        const std::int64_t v = slice.data()[k];
        if (cfg_.forbid_overlap && !is_available(v)) continue;

        scratch_accepted_[k] = 1;
        scar.pixels.push_back(v);
        ++n_accepted;
    }

    if (n_accepted < static_cast<std::size_t>(n_target))
    {
        scar.pixels.clear();
        return false;
    }

    scar.initiation_flat = init_flat;
    scar.n_pixels        = n_target;

    scar.host_fill_ratio = static_cast<double>(n_target) / static_cast<double>(n_up);

    if (cfg_.forbid_overlap)
        for (std::int64_t f : scar.pixels)
            event_blocked_[static_cast<std::size_t>(f)] = 1;

    return true;
}


bool LandslidePlacer::place_scar_seeded(int n_target_pixels,
                                        std::uint64_t& rng_state,
                                        LandslideScar& scar)
{
    return place_scar(n_target_pixels,
                      [&rng_state]() { return next_uniform(rng_state); },
                      scar);
}


void LandslidePlacer::finalise_diagnostics()
{
    diag_.median_host_fill_ratio = compute_median_of_copy(fill_ratios_, 0.0);
}

}  // namespace lsdtt
