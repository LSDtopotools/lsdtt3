// lsdtt_valley_centreline_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_valley_centreline_cpu
//
// Implementation of broad per-channel valley identification and the
// trough-routed valley centreline. See header for the algorithm outline.
// -----------------------------------------------------------------------------

#include "lsdtt_valley_centreline_cpu.hpp"

#include "lsdtt_raster_binary_mask_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_proximity_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"

namespace lsdtt
{
namespace
{

constexpr std::uint8_t k_off = 0u;
constexpr std::uint8_t k_on  = 1u;

// 8-neighbour offsets.
constexpr int k_dr[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };
constexpr int k_dc[8] = { -1,  0,  1, -1, 1, -1, 0, 1 };

// Dilate a binary mask by n iterations of 8-neighbour growth, staying inside
// the valid domain.
void dilate_mask(
    std::vector<std::uint8_t>&       mask,
    const std::vector<std::uint8_t>& domain,
    int iterations, int n_rows, int n_cols)
{
    for (int it = 0; it < iterations; ++it)
    {
        std::vector<std::uint8_t> grown = mask;
        for (int r = 0; r < n_rows; ++r)
        {
            for (int c = 0; c < n_cols; ++c)
            {
                const std::size_t i =
                    static_cast<std::size_t>(r) * n_cols + c;
                if (mask[i] || !domain[i]) continue;
                for (int k = 0; k < 8; ++k)
                {
                    const int rr = r + k_dr[k], cc = c + k_dc[k];
                    if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                        continue;
                    if (mask[static_cast<std::size_t>(rr) * n_cols + cc])
                    {
                        grown[i] = k_on;
                        break;
                    }
                }
            }
        }
        mask.swap(grown);
    }
}


} // namespace


// -----------------------------------------------------------------------------
// fit_pinned_monotonic_baseline (public - see header)
// -----------------------------------------------------------------------------
std::vector<float> fit_pinned_monotonic_baseline(
    const std::vector<double>& along_distance,
    const std::vector<float>&  pin_elevation)
{
    const std::size_t n = along_distance.size();
    std::vector<float> baseline(
        n, std::numeric_limits<float>::quiet_NaN());
    if (n == 0 || pin_elevation.size() != n) return baseline;

    // Collect pins in node order.
    std::vector<std::size_t> pin_idx;
    for (std::size_t i = 0; i < n; ++i)
        if (std::isfinite(pin_elevation[i])) pin_idx.push_back(i);
    if (pin_idx.empty()) return baseline;

    // Monotonic envelope: each pin is lowered to the minimum of all pins at
    // or upstream of it (forward prefix-min). The result is non-increasing
    // downstream and <= every original pin, so the baseline sits AT the
    // channel elevation where the channel descends normally and BELOW it
    // wherever the projected channel elevation locally rises downstream (a
    // baseline that rose above the channel would give the channel a negative
    // height-above-channel, which must never happen).
    std::vector<float> pin_adj(pin_idx.size());
    float running = std::numeric_limits<float>::infinity();
    for (std::size_t p = 0; p < pin_idx.size(); ++p)
    {
        running = std::min(running, pin_elevation[pin_idx[p]]);
        pin_adj[p] = running;
    }

    // Flat before the first pin, linear between pins, flat after the last.
    std::size_t p = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (i <= pin_idx.front())      { baseline[i] = pin_adj.front(); continue; }
        if (i >= pin_idx.back())       { baseline[i] = pin_adj.back();  continue; }
        while (p + 1 < pin_idx.size() && pin_idx[p + 1] < i) ++p;
        const std::size_t ia = pin_idx[p], ib = pin_idx[p + 1];
        const double da = along_distance[ia], db = along_distance[ib];
        const double f = (db > da)
            ? (along_distance[i] - da) / (db - da)
            : 0.0;
        baseline[i] = static_cast<float>(
            pin_adj[p] + f * (pin_adj[p + 1] - pin_adj[p]));
    }
    return baseline;
}


valley_centreline_result extract_valley_centrelines(
    const lsdtt_raster_cpu&          dem,
    const channel_scope_result&      scope,
    const valley_centreline_request& request)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>&   elev = dem.get_data();
    const std::vector<uint8_t>& dmask = dem.get_mask();
    const int n_rows = meta.n_rows, n_cols = meta.n_cols;
    const std::size_t n = elev.size();
    const double cell = meta.cell_size;

    if (scope.channel_node_lists.empty())
        throw std::invalid_argument(
            "extract_valley_centrelines: channel scope has no channels.");
    if (scope.channel_mask.size() != n)
        throw std::invalid_argument(
            "extract_valley_centrelines: scope / DEM size mismatch.");

    // (The medial-axis least-cost centreline is confined to the valley mask,
    // which excludes nodata, so reservoir/water-surface voids are rounded
    // rather than routed into - no nodata plugging is needed.)

    valley_centreline_result result;
    std::vector<float>   vm_data(n, 0.0f);
    std::vector<uint8_t> vm_mask = dmask;


    for (std::size_t s = 0; s < scope.channel_node_lists.size(); ++s)
    {
        const std::vector<std::int64_t>& nodes = scope.channel_node_lists[s];
        valley_channel_centreline rec;
        rec.source_id = static_cast<std::int32_t>(s);
        if (nodes.size() < 2)
        {
            result.channels.push_back(rec);
            continue;
        }

        // ---------------------------------------------------------------------
        // 1. Corridor and Euclidean HAND relative to THIS channel only.
        // ---------------------------------------------------------------------
        std::vector<float> chan_data(n, 0.0f);
        for (const std::int64_t f : nodes)
            chan_data[static_cast<std::size_t>(f)] = 1.0f;
        lsdtt_raster_cpu chan_binary(std::move(chan_data),
                                     std::vector<uint8_t>(dmask), meta);

        lsdtt_raster_cpu dist =
            euclidean_distance_transform(chan_binary, request.n_threads);
        lsdtt_raster_cpu hand = difference_to_nearest_true_mask_value(
            chan_binary, dem, request.n_threads);

        const std::vector<float>&   dd = dist.get_data();
        const std::vector<uint8_t>& dm = dist.get_mask();
        const std::vector<float>&   hd = hand.get_data();
        const std::vector<uint8_t>& hm = hand.get_mask();

        // The break levels, and each corridor pixel's window index, are read
        // again by the baseline refinement and the reach-id raster further
        // down. They are materialised as plain data rather than left as a
        // lambda over step-2 locals, which would dangle once that block ends.
        // Both stay empty on the supplied-mask path, where there are no
        // windows and no breaks.
        std::vector<double> window_break;
        std::vector<int>    window_index;

        // The corridor test is shared by both paths below and by the baseline
        // refinement further down, so it is declared before the branch.
        auto in_corridor = [&](std::size_t i) -> bool
        {
            return dmask[i] && dm[i] && hm[i] &&
                   std::isfinite(dd[i]) && std::isfinite(hd[i]) &&
                   static_cast<double>(dd[i]) <= request.corridor_width_metres;
        };

        // Signed along-channel coordinate, for the trough's upstream test below.
        //
        // WHAT IT IS. Each corridor pixel is assigned the along-channel distance
        // of its nearest channel node, PLUS its own offset from that node
        // projected onto the local channel direction. The projection is what
        // makes it signed: a pixel beyond the channel head projects backwards
        // and lands negative, and the size of the negative number says how far
        // beyond the head it sits.
        //
        // WHY THE TROUGH NEEDS IT. Terrain upstream of the channel head is
        // still valley ground and stays in the valley mask, but it must not be
        // carved into the trough. A broad pocket behind the head forms a deep
        // distance-to-wall bowl whose lowest spill is the upstream tile edge,
        // so the least-cost walk escapes the wrong way and the centreline runs
        // off the top of the DEM instead of down the valley.
        //
        // This used to be computed inside the corridor HAND break search, which
        // is why it went inert when identification moved to
        // lsdtt-fluvial-surfaces at phase 5d. It never actually depended on
        // that search - it needs only the channel nodes and the corridor - so it
        // is computed here, next to the one thing that reads it.
        std::vector<float>   along_tilt(n, 0.0f);
        std::vector<uint8_t> along_mask(n, 0u);
        {
            // Cumulative along-channel distance at each node, from the head.
            std::vector<float> node_along(nodes.size(), 0.0f);
            for (std::size_t k = 1; k < nodes.size(); ++k)
            {
                const std::int64_t a2 = nodes[k - 1], b2 = nodes[k];
                const int dr2 = static_cast<int>(a2 / n_cols) -
                                static_cast<int>(b2 / n_cols);
                const int dc2 = static_cast<int>(a2 % n_cols) -
                                static_cast<int>(b2 % n_cols);
                node_along[k] = node_along[k - 1] + static_cast<float>(
                    cell * std::sqrt(
                        static_cast<double>(dr2 * dr2 + dc2 * dc2)));
            }

            // Inflate the nearest channel NODE INDEX across the corridor.
            std::vector<float>   sv(n, 0.0f);
            std::vector<uint8_t> sm(n, 0u);
            std::vector<float>   pd(n, 1.0f);
            std::vector<uint8_t> pm(n, 0u);
            for (std::size_t k = 0; k < nodes.size(); ++k)
            {
                const std::size_t f = static_cast<std::size_t>(nodes[k]);
                sv[f] = static_cast<float>(k);
                sm[f] = 1u;
            }
            for (std::size_t i = 0; i < n; ++i)
                pm[i] = (in_corridor(i) || sm[i]) ? 1u : 0u;

            const lsdtt_raster_cpu inflated = inflate_values_to_nearest_source(
                lsdtt_raster_cpu(std::move(sv), std::move(sm), meta),
                lsdtt_raster_cpu(std::move(pd), std::move(pm), meta),
                request.n_threads);
            const std::vector<float>&   idx_d = inflated.get_data();
            const std::vector<uint8_t>& idx_m = inflated.get_mask();

            // Local channel direction over +/- span nodes, so a single kinked
            // node cannot flip the sign of everything beside it.
            const int span = 10;
            const int last = static_cast<int>(nodes.size()) - 1;
            for (std::size_t i = 0; i < n; ++i)
            {
                if (!idx_m[i]) continue;
                const int k = std::clamp(
                    static_cast<int>(std::llround(idx_d[i])), 0, last);
                along_mask[i] = 1u;

                const int ka = std::max(k - span, 0);
                const int kb = std::min(k + span, last);
                const std::int64_t fa = nodes[static_cast<std::size_t>(ka)];
                const std::int64_t fb = nodes[static_cast<std::size_t>(kb)];
                const std::int64_t fk = nodes[static_cast<std::size_t>(k)];

                double dx = static_cast<double>((fb % n_cols) - (fa % n_cols));
                double dy = static_cast<double>((fb / n_cols) - (fa / n_cols));
                const double len = std::hypot(dx, dy);
                double proj = 0.0;
                if (len > 0.0)
                {
                    dx /= len; dy /= len;
                    const double px = static_cast<double>(
                        static_cast<std::int64_t>(i % n_cols) - (fk % n_cols));
                    const double py = static_cast<double>(
                        static_cast<std::int64_t>(i / n_cols) - (fk / n_cols));
                    proj = (px * dx + py * dy) * cell;
                }
                along_tilt[i] = node_along[static_cast<std::size_t>(k)] +
                                static_cast<float>(proj);
            }
        }

        // ---------------------------------------------------------------------
        // 2/3. The valley for this channel.
        //
        // Either taken from a supplied mask, or derived here by searching the
        // corridor HAND population for a break. The derived path is the whole
        // of what phase 5 of the fluvial-surfaces plan moves out of this
        // engine; it is kept for one release and deleted with the valley_*
        // break parameters that drive it.
        // ---------------------------------------------------------------------
        std::vector<std::uint8_t> valley(n, k_off);

        if (request.supplied_valley_mask != nullptr)
        {
            // The mask says which ground is valley; the corridor says which
            // channel owns it. Intersecting them is what keeps a shared valley
            // floor from being claimed in full by every channel running through
            // it, exactly as the derived path does.
            const std::vector<float>&   sv  = 
                request.supplied_valley_mask->get_data();
            const std::vector<uint8_t>& svm = 
                request.supplied_valley_mask->get_mask();
            for (std::size_t i = 0; i < n; ++i)
            {
                if (!in_corridor(i)) continue;
                if (svm[i] && sv[i] > 0.5f) valley[i] = k_on;
            }

            // No break was found, so none is reported. NaN rather than 0.0:
            // a zero would read as a real break at the channel's own height.
            rec.valley_break_derived     = false;
            rec.valley_hand_break_metres =
                std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            // The corridor HAND break search used to live here: about 360
            // lines that found a break per along-channel window and thresholded
            // the corridor against it. It moved to lsdtt-fluvial-surfaces at
            // phase 5d, where the same job is done by banding a reference field
            // and taking band 0 - a method that also works on a fan, which a
            // corridor HAND break does not.
            throw std::invalid_argument(
                "extract_valley_centrelines: no valley mask supplied. This tool "
                "no longer identifies the valley; run lsdtt-fluvial-surfaces "
                "first and pass its <stem>_valley_mask.tif.");
        }

        for (const std::int64_t f : nodes)
            if (dmask[static_cast<std::size_t>(f)])
                valley[static_cast<std::size_t>(f)] = k_on;

        // Generic constrained morphology, seeded by this channel: see
        // lsdtt_raster_binary_mask_cpu.
        filter_components_by_seed_reachability_in_place(
            valley, dd, dm, nodes,
            request.valley_dilation_pixels,
            request.restrict_to_channel_connected,
            request.valley_max_detour_ratio,
            cell, n_rows, n_cols);

        rec.n_valley_pixels = 0;
        for (std::size_t i = 0; i < n; ++i)
            if (valley[i]) ++rec.n_valley_pixels;

        if (rec.n_valley_pixels <
            static_cast<std::int64_t>(request.min_centreline_nodes))
        {
            result.channels.push_back(rec);
            continue;
        }

        // (The combined valley raster is claimed at the end of the channel
        // loop, after the optional baseline refinement of the mask.)

        // ---------------------------------------------------------------------
        // 4. Trough surface: hole = valley dilated by the buffer; walls keep
        //    the DEM. Hole pixels get (min bank elevation - k * distance to
        //    wall), iteratively breached/filled.
        // ---------------------------------------------------------------------
        std::vector<std::uint8_t> hole = valley;
        // Exclude terrain upstream of the channel head (negative signed
        // along-channel coordinate) from the trough. A broad valley pocket
        // behind the head otherwise forms a deep distance-to-wall bowl whose
        // spill is the upstream tile edge, and the routed path escapes the
        // wrong way. The valley MASK keeps those pixels; only the trough
        // ignores them.
        std::int64_t n_upstream_excluded = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (hole[i] && along_mask[i] &&
                along_tilt[i] < -2.0f * static_cast<float>(cell))
            {
                hole[i] = k_off;
                ++n_upstream_excluded;
            }
        }
        rec.n_upstream_of_head_excluded = n_upstream_excluded;
        if (request.nodata_buffer_pixels > 0)
            dilate_mask(hole, dmask, request.nodata_buffer_pixels,
                        n_rows, n_cols);

        // Distance to the valley wall (sources = valid non-hole pixels).
        std::vector<float> wall_src(n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
            if (dmask[i] && !hole[i]) wall_src[i] = 1.0f;
        lsdtt_raster_cpu wall_dist = euclidean_distance_transform(
            lsdtt_raster_cpu(std::move(wall_src),
                             std::vector<uint8_t>(dmask), meta),
            request.n_threads);
        const std::vector<float>& wd = wall_dist.get_data();

        // ---------------------------------------------------------------------
        // 5. Centreline = least-cost path down the medial axis.
        // ---------------------------------------------------------------------
        // The legacy trough (min_bank - k*dist_to_wall - tilt*along), filled
        // and D8-routed, CANNOT trace the centreline of a WIDE valley: filling
        // flattens the broad-topped distance-to-wall basin into a lake whose
        // priority-flood drainage runs straight to the pour point with
        // right-angle artifacts, and the path is invariant to k_scale, tilt
        // and loop count (verified). Instead walk a Dijkstra LEAST-COST PATH
        // from the channel head to the outlet whose per-step cost penalises
        // proximity to a wall, (Wmax - dist_to_wall): the cheapest route hugs
        // the locus of MAXIMUM distance-to-wall - the medial axis / geometric
        // centreline (the point equidistant from the two nearest walls) - at
        // any valley width. It is confined to the valley (hole) mask, so it
        // rounds reservoir voids and never leaves the valley; the head/outlet
        // anchors fix the downstream direction, so no synthetic tilt is
        // needed and the path cannot drain out the wrong tile edge.

        // Snap the head and terminus onto the nearest hole pixel.
        auto snap_to_hole = [&](std::int64_t f) -> std::int64_t
        {
            if (f >= 0 && f < static_cast<std::int64_t>(n) &&
                hole[static_cast<std::size_t>(f)]) return f;
            const int fr = f >= 0 ? static_cast<int>(f / n_cols) : 0;
            const int fc = f >= 0 ? static_cast<int>(f % n_cols) : 0;
            std::int64_t best = -1;
            double best_d2 = std::numeric_limits<double>::max();
            for (std::size_t i = 0; i < n; ++i)
            {
                if (!hole[i]) continue;
                const int rr = static_cast<int>(i / n_cols);
                const int cc = static_cast<int>(i % n_cols);
                const double d2 = static_cast<double>(rr - fr) * (rr - fr) +
                                  static_cast<double>(cc - fc) * (cc - fc);
                if (d2 < best_d2)
                {
                    best_d2 = d2;
                    best = static_cast<std::int64_t>(i);
                }
            }
            return best;
        };

        const bool outlet_bounded =
            s < scope.resolved_outlets.size() &&
            scope.resolved_outlets[s] >= 0;
        const std::int64_t start_f = snap_to_hole(nodes.front());
        const std::int64_t target_f = snap_to_hole(
            outlet_bounded ? scope.resolved_outlets[s] : nodes.back());

        std::vector<std::int64_t> path;
        if (start_f >= 0 && target_f >= 0 && start_f != target_f)
        {
            float wmax = 0.0f;
            for (std::size_t i = 0; i < n; ++i)
                if (hole[i]) wmax = std::max(wmax, wd[i]);
            // Small positive length penalty (m): keeps every step cost > 0 so
            // the path does not take absurd detours to gain a millimetre of
            // centrality, and breaks ties toward the shorter central route.
            const double eps = 0.5;

            std::vector<double> cost(
                n, std::numeric_limits<double>::infinity());
            std::vector<std::int64_t> prev(n, -1);
            using qnode = std::pair<double, std::int64_t>;
            std::priority_queue<qnode, std::vector<qnode>,
                                std::greater<qnode>> pq;
            cost[static_cast<std::size_t>(start_f)] = 0.0;
            pq.emplace(0.0, start_f);
            const double straight = cell;
            const double diag = cell * std::sqrt(2.0);
            while (!pq.empty())
            {
                const auto [d, flat] = pq.top();
                pq.pop();
                const std::size_t fi = static_cast<std::size_t>(flat);
                if (d > cost[fi]) continue;
                if (flat == target_f) break;
                const int r = static_cast<int>(flat / n_cols);
                const int c = static_cast<int>(flat % n_cols);
                for (int k = 0; k < 8; ++k)
                {
                    const int rr = r + k_dr[k], cc = c + k_dc[k];
                    if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                        continue;
                    const std::size_t nf =
                        static_cast<std::size_t>(rr) * n_cols + cc;
                    if (!hole[nf]) continue;
                    const double step =
                        (k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight;
                    const double w = step *
                        (static_cast<double>(wmax - wd[nf]) + eps);
                    const double nd = d + w;
                    if (nd < cost[nf])
                    {
                        cost[nf] = nd;
                        prev[nf] = flat;
                        pq.emplace(nd, static_cast<std::int64_t>(nf));
                    }
                }
            }
            if (std::isfinite(cost[static_cast<std::size_t>(target_f)]))
            {
                std::int64_t cur = target_f;
                while (cur != -1)
                {
                    path.push_back(cur);
                    if (cur == start_f) break;
                    cur = prev[static_cast<std::size_t>(cur)];
                }
                std::reverse(path.begin(), path.end());
            }
        }

        rec.routing_succeeded =
            static_cast<int>(path.size()) >= request.min_centreline_nodes;
        if (!path.empty())
        {
            // Synthetic monotonic thalweg for the trough-elevation field: the
            // running minimum of the DEM downstream along the centreline (the
            // real profile is the pinned baseline computed in step 6).
            rec.centreline_trough_elevation.reserve(path.size());
            float run = std::numeric_limits<float>::infinity();
            for (const std::int64_t f : path)
            {
                const std::size_t fi = static_cast<std::size_t>(f);
                if (dmask[fi] && elev[fi] < run) run = elev[fi];
                rec.centreline_trough_elevation.push_back(
                    std::isfinite(run) ? run : elev[fi]);
            }
            rec.centreline_flats = std::move(path);
        }

        // ---------------------------------------------------------------------
        // 6. Pinned baseline along the centreline. The centreline, resampled
        //    on the pre-carve topography, rides up and down risers and onto
        //    higher surfaces; the baseline pins to the original channel
        //    elevations wherever the channel crosses (or touches) the
        //    centreline and interpolates linearly between pins, monotonic
        //    (non-increasing) downstream and never above the channel (see
        //    fit_pinned_monotonic_baseline).
        // ---------------------------------------------------------------------
        if (!rec.centreline_flats.empty())
        {
            const std::size_t np = rec.centreline_flats.size();
            std::vector<double> along_path(np, 0.0);
            std::vector<float>  pins(
                np, std::numeric_limits<float>::quiet_NaN());
            rec.centreline_is_pin.assign(np, 0);
            for (std::size_t p = 0; p < np; ++p)
            {
                const std::int64_t f = rec.centreline_flats[p];
                if (p > 0)
                {
                    const std::int64_t g = rec.centreline_flats[p - 1];
                    const int dr2 = static_cast<int>(f / n_cols) -
                                    static_cast<int>(g / n_cols);
                    const int dc2 = static_cast<int>(f % n_cols) -
                                    static_cast<int>(g % n_cols);
                    along_path[p] = along_path[p - 1] + cell *
                        std::sqrt(static_cast<double>(
                            dr2 * dr2 + dc2 * dc2));
                }
            }

            // Pins: every channel node projects onto its nearest centreline
            // node and pins the baseline to the channel elevation there. A
            // crossing pins exactly; a laterally offset channel reach (wide
            // valley, channel meandering away from the geometric centre)
            // still pins the baseline at its along-valley position - this is
            // what strips the meander noise. Where several channel nodes
            // project to one centreline node, the LOWEST elevation wins: the
            // baseline is a thalweg reference and must sit at or below the
            // channel (a higher pin would raise it above the channel and make
            // the channel's own height-above-channel negative).
            for (const std::int64_t cf : nodes)
            {
                const double cr = static_cast<double>(cf / n_cols);
                const double cc2 = static_cast<double>(cf % n_cols);
                std::size_t best = 0;
                double best_d2 = std::numeric_limits<double>::max();
                for (std::size_t p = 0; p < np; ++p)
                {
                    const std::int64_t f = rec.centreline_flats[p];
                    const double dr2 =
                        static_cast<double>(f / n_cols) - cr;
                    const double dc2 =
                        static_cast<double>(f % n_cols) - cc2;
                    const double d2 = dr2 * dr2 + dc2 * dc2;
                    if (d2 < best_d2) { best_d2 = d2; best = p; }
                }
                const float ce = elev[static_cast<std::size_t>(cf)];
                if (!std::isfinite(pins[best]) || ce < pins[best])
                    pins[best] = ce;
                rec.centreline_is_pin[best] = 1;
            }
            rec.centreline_baseline_elevation =
                fit_pinned_monotonic_baseline(along_path, pins);
        }

        // ---------------------------------------------------------------------
        // 7. Baseline refinement of the valley mask (see request comment):
        //    valley = (elevation - baseline of the nearest centreline node)
        //    <= the reach's break. Kills the stale-reference wedges that the
        //    nearest-channel-pixel test grows onto downstream-sloping
        //    terraces wherever the channel veers away.
        // ---------------------------------------------------------------------
        // The refinement re-thresholds height-above-baseline. It used to reuse
        // the per-window break this engine derived; identification has moved
        // out, so the threshold is now an explicit request field and the
        // refinement simply does not run without one. Refusing rather than
        // inventing a number is the same rule the ladder applies to an
        // unresolvable level count.
        //
        // NO FEED-BACK: the refined mask is written for the user and is never
        // returned to lsdtt-fluvial-surfaces.
        if (request.refine_valley_with_baseline &&
            request.refine_max_height_metres > 0.0 &&
            rec.routing_succeeded &&
            !rec.centreline_baseline_elevation.empty())
        {
            std::vector<float>   sv(n, 0.0f);
            std::vector<uint8_t> sm(n, 0u);
            std::vector<float>   pd3(n, 1.0f);
            std::vector<uint8_t> pm3(n, 0u);
            for (std::size_t p = 0; p < rec.centreline_flats.size(); ++p)
            {
                const float b = rec.centreline_baseline_elevation[p];
                if (!std::isfinite(b)) continue;
                const std::size_t f =
                    static_cast<std::size_t>(rec.centreline_flats[p]);
                sv[f] = b;
                sm[f] = 1u;
            }
            for (std::size_t i = 0; i < n; ++i)
                pm3[i] = (in_corridor(i) || sm[i]) ? 1u : 0u;
            try
            {
                lsdtt_raster_cpu base_infl = inflate_values_to_nearest_source(
                    lsdtt_raster_cpu(std::move(sv), std::move(sm), meta),
                    lsdtt_raster_cpu(std::move(pd3), std::move(pm3), meta),
                    request.n_threads);
                const std::vector<float>&   bd = base_infl.get_data();
                const std::vector<uint8_t>& bm = base_infl.get_mask();

                std::vector<std::uint8_t> refined(n, k_off);
                for (std::size_t i = 0; i < n; ++i)
                {
                    if (!in_corridor(i) || !bm[i]) continue;
                    const double rel =
                        static_cast<double>(elev[i]) -
                        static_cast<double>(bd[i]);
                    if (rel <= request.refine_max_height_metres)
                        refined[i] = k_on;
                }
                // The channel itself is always valley (and anchors the
                // patch-adjacency filter).
                for (const std::int64_t f : nodes)
                    if (dmask[static_cast<std::size_t>(f)])
                        refined[static_cast<std::size_t>(f)] = k_on;

                filter_components_by_seed_reachability_in_place(
                                     refined, dd, dm, nodes,
                                     request.valley_dilation_pixels,
                                     request.restrict_to_channel_connected,
                                     request.valley_max_detour_ratio,
                                     cell, n_rows, n_cols);
                valley.swap(refined);
                rec.n_valley_pixels = 0;
                for (std::size_t i = 0; i < n; ++i)
                    if (valley[i]) ++rec.n_valley_pixels;
            }
            catch (const std::exception&)
            {
                // Keep the first-pass mask on any failure.
            }
        }

        // Claim into the combined valley raster (first claim wins).
        for (std::size_t i = 0; i < n; ++i)
            if (valley[i] && vm_data[i] == 0.0f)
                vm_data[i] = static_cast<float>(s + 1);

        if (request.store_per_channel_valley_masks)
        {
            std::vector<float> md(n, 0.0f);
            for (std::size_t i = 0; i < n; ++i)
                if (valley[i]) md[i] = 1.0f;
            rec.valley_mask = lsdtt_raster_cpu(
                std::move(md), std::vector<uint8_t>(dmask), meta);
        }

        result.channels.push_back(std::move(rec));
    }

    result.valley_mask =
        lsdtt_raster_cpu(std::move(vm_data), std::move(vm_mask), meta);
    return result;
}

} // namespace lsdtt
