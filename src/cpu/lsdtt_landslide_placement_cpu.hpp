// lsdtt_landslide_placement_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Landslide PLACEMENT: turning a target pixel count into a scar on a landscape.
//
// The size distribution lives in lsdtt_landslide_size_cpu.hpp and knows nothing
// about topography. This file is the other half: given "I want a landslide of
// n pixels", it finds somewhere that can host one and returns the pixels.
//
// ---------------------------------------------------------------------------
// THE ORDER MATTERS: SAMPLE THE AREA FIRST, THEN CHOOSE A HOST.
//
// A pixel can host an n-pixel landslide only if it has at least n pixels
// upslope of it. lsdtt_flowinfo gives that count in O(1) via
// s_vector_size_unchecked, so eligibility is tested BEFORE any work is done and
// a scar always reaches its target area exactly. The obvious alternative --
// pick a location, then grow a scar outwards until it hits the target -- can
// run out of hillslope part way and leave a truncated scar, which silently
// biases the realised size distribution low at the large end. That failure mode
// is designed out here rather than detected later.
// ---------------------------------------------------------------------------
//
// Three properties follow from building the scar out of an upslope set:
//
//   1. The hillslope constraint is automatic. Contributing area decreases
//      monotonically upslope, so if the initiation pixel is below the channel
//      threshold then every pixel above it is too. There is no per-pixel
//      channel test anywhere in this file -- the invariant is structural, and
//      the test suite guards it rather than the code enforcing it.
//
//   2. Connectivity is guaranteed. An upslope set is a connected subtree.
//
//   3. Cost is O(n_up) per scar, against O(n * 8) neighbour scans with a
//      priority comparison for an incremental flood fill.
//
// TRIMMING: the upslope set is usually larger than the target, so the FARTHEST
// pixels are removed until the count matches. "Farthest" is flow distance from
// the initiation pixel, computed exactly in one O(n) pass -- the s_vector is a
// preorder DFS, so a node's receiver is always processed before the node
// itself, and s_vector_start gives each node's position within the slice in
// O(1). The result is a compact scar hugging the initiation point rather than a
// thin strip up a single tributary (which is what simply truncating the DFS
// order would give).
//
// A CONSEQUENCE TO BE AWARE OF: hosting an n-pixel landslide needs n pixels
// upslope, and contributing area is itself power-law distributed, so the
// eligible set shrinks rapidly as n grows. The largest landslides can only
// initiate just above channel heads and will occupy most of a zero-order basin.
// That is defensible -- large landslides really do evacuate whole zero-order
// basins -- but it is an imposed spatial correlation, not a chosen one, so
// host_fill_ratio is reported for every scar.
// -----------------------------------------------------------------------------

#ifndef LSDTT_LANDSLIDE_PLACEMENT_CPU_HPP
#define LSDTT_LANDSLIDE_PLACEMENT_CPU_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// One placed landslide.
struct LandslideScar
{
    std::int64_t initiation_flat = -1;
    int          n_pixels        = 0;
    std::vector<std::int64_t> pixels;   // includes the initiation pixel

    // n_pixels / upslope count of the host. 1.0 means the scar filled its host
    // basin completely; small values mean it was trimmed hard. The median of
    // this across an event is the diagnostic for channel-head clustering.
    double host_fill_ratio = 0.0;
};


struct LandslidePlacementConfig
{
    // A landslide may not be larger than the hillslope it comes from, so
    // initiation is restricted to pixels below this. Expressed in PIXELS to
    // match lsdtt3's threshold_contributing_pixels convention.
    std::int64_t channel_threshold_pixels = 1000;

    // Attempts to find a usable host before giving up on one landslide. Needed
    // because previously placed scars remove pixels from later hosts: a pixel
    // can pass the O(1) upslope-count test and still have too few AVAILABLE
    // pixels once overlaps are excluded.
    int max_host_attempts = 64;

    // When true, a scar may not reuse any pixel already excavated in this event
    // or still annealing. Turning it off allows overlapping scars, which double
    // counts eroded volume -- provided only for experiments.
    bool forbid_overlap = true;
};


struct LandslidePlacementDiagnostics
{
    int    n_placed              = 0;
    int    n_failed              = 0;   // no host found within max_host_attempts
    long   n_host_attempts       = 0;   // total draws, including rejected ones
    double median_host_fill_ratio = 0.0;
    // Largest n_target that could be hosted anywhere on this landscape, in
    // pixels. Reported so a caller can see when the distribution is asking for
    // landslides the landscape cannot hold.
    std::int64_t max_hostable_pixels = 0;
};


// -----------------------------------------------------------------------------
// LandslidePlacer
//
// Holds the eligible-host list for one landscape, plus the per-event record of
// which pixels are unavailable. Build once, call place_scar() repeatedly.
//
// Takes the ABSTRACT lsdtt_flow_base, not the concrete lsdtt_flowinfo, so that
// a CLI can build routing through make_flow_base() and pass it here while
// including only ../core/lsdtt_flow_base.hpp. lsdtt_flowinfo_cpu.hpp must never
// reach the CLI layer. Every query used here -- s_vector_size_unchecked,
// s_vector_start_unchecked, upslope_nodes, receiver_of_unchecked,
// is_valid_cell -- is on the base interface, so nothing is given up; the same
// call sites work unchanged against a GPU backend when one exists.
//
// The flow object must outlive the placer. The DEM is not needed here:
// placement selects pixels, it does not move elevation.
// -----------------------------------------------------------------------------
class LandslidePlacer
{
  public:
    // `mask` is optional (nullptr = whole DEM). A pixel is eligible only where
    // the mask is valid and non-zero.
    LandslidePlacer(const lsdtt_flow_base& flow,
                    const lsdtt_raster_meta& meta,
                    const LandslidePlacementConfig& config,
                    const lsdtt_raster_cpu* mask = nullptr);

    // Mark pixels as unavailable because they are still annealing. Call before
    // an event. `time_of_last_failure` uses the same no_data convention as any
    // raster; cells that have never failed must be nodata or a large negative
    // sentinel.
    void apply_annealing(const lsdtt_raster_cpu& time_of_last_failure,
                         double t_now_yr,
                         double anneal_time_yr);

    // Clear the per-event record of excavated pixels (but not the annealing
    // state, which comes from apply_annealing).
    void begin_event();

    // Place one landslide of exactly n_target pixels.
    //
    // `draw_uniform` must return values in [0,1); it is called once per host
    // attempt. Returns false if no host was found within max_host_attempts, in
    // which case `scar` is left untouched -- a caller must treat that as a
    // reportable shortfall, never as a silent zero.
    template <typename UniformFn>
    bool place_scar(int n_target_pixels, UniformFn&& draw_uniform,
                    LandslideScar& scar);

    // Non-template convenience overload taking a seeded generator.
    bool place_scar_seeded(int n_target_pixels, std::uint64_t& rng_state,
                           LandslideScar& scar);

    const LandslidePlacementDiagnostics& diagnostics() const noexcept { return diag_; }
    void finalise_diagnostics();

    // Largest landslide, in pixels, that any pixel on this landscape could host
    // given the channel threshold. Computed at construction.
    std::int64_t max_hostable_pixels() const noexcept { return max_hostable_; }

    // Number of pixels eligible to host a landslide of n pixels, ignoring
    // availability. Exposed for tests and for reporting.
    std::size_t n_eligible_hosts(int n_target_pixels) const;

    bool is_available(std::int64_t flat) const
    {
        const std::size_t i = static_cast<std::size_t>(flat);
        return !anneal_blocked_[i] && !event_blocked_[i];
    }

  private:
    // Candidate hosts, sorted by upslope pixel count DESCENDING, so the hosts
    // able to take an n-pixel slide are always a prefix of this vector and can
    // be found by binary search.
    struct Host { std::int64_t flat; std::int64_t n_up; };

    // Attempt one specific host. Returns false if it cannot supply n_target
    // available pixels.
    bool try_host(std::int64_t init_flat, int n_target, LandslideScar& scar);

    const lsdtt_flow_base&   flow_;
    lsdtt_raster_meta        meta_;
    LandslidePlacementConfig cfg_;

    std::vector<Host>          hosts_;
    // Two SEPARATE blocking records, because they have different lifetimes.
    // anneal_blocked_ comes from the time-of-last-failure raster and is rebuilt
    // whenever apply_annealing is called; event_blocked_ records what THIS
    // event has already excavated and is cleared by begin_event. Merging them
    // into one array leaks excavation state from one event into the next, which
    // starves later events of hosts.
    std::vector<char>          anneal_blocked_;
    std::vector<char>          event_blocked_;
    std::int64_t               max_hostable_ = 0;
    LandslidePlacementDiagnostics diag_{};
    std::vector<double>        fill_ratios_;  // for the median

    // Scratch, reused across calls to avoid reallocating per landslide.
    std::vector<double>        scratch_dist_;
    std::vector<std::int64_t>  scratch_rec_lp_;   // receiver's position in slice
    std::vector<std::uint32_t> scratch_order_;
    std::vector<char>          scratch_accepted_;
};


// -----------------------------------------------------------------------------
// Template implementation
// -----------------------------------------------------------------------------
template <typename UniformFn>
bool LandslidePlacer::place_scar(int n_target_pixels, UniformFn&& draw_uniform,
                                 LandslideScar& scar)
{
    if (n_target_pixels < 1) return false;

    const std::size_t n_hosts = n_eligible_hosts(n_target_pixels);
    if (n_hosts == 0)
    {
        ++diag_.n_failed;
        return false;
    }

    for (int attempt = 0; attempt < cfg_.max_host_attempts; ++attempt)
    {
        ++diag_.n_host_attempts;
        const double u = draw_uniform();
        std::size_t idx = static_cast<std::size_t>(u * static_cast<double>(n_hosts));
        if (idx >= n_hosts) idx = n_hosts - 1;

        if (try_host(hosts_[idx].flat, n_target_pixels, scar))
        {
            ++diag_.n_placed;
            fill_ratios_.push_back(scar.host_fill_ratio);
            return true;
        }
    }

    ++diag_.n_failed;
    return false;
}

}  // namespace lsdtt

#endif  // LSDTT_LANDSLIDE_PLACEMENT_CPU_HPP
