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

// Valley-mask cleanup (see request comment in the header). This is a
// morphological opening constrained by distance to the channel:
// 1. erode by n_dilate iterations (the channel and its immediate neighbours
//    are protected, so the trunk valley cannot be disconnected) — this
//    severs the corridors that creep along tributary floodplains and connect
//    far-off flats (e.g. a neighbouring channel's floodplain) into the
//    valley;
// 2. label 8-connected patches and drop those with no pixel on or 8-adjacent
//    to the channel;
// 3. dilate 2 * n_dilate times (regrow the eroded rim + fill interior
//    holes), retaining a dilated pixel only when its distance to the channel
//    does not exceed the farthest pixel of the patch it grows from, plus a
//    margin covering the erosion depth — so patches regrow and holes fill,
//    but nothing extends outward away from the channel.
void clean_valley_patches(
    std::vector<std::uint8_t>&       valley,
    const std::vector<float>&        dist_to_channel,
    const std::vector<std::uint8_t>& dist_mask,
    const std::vector<std::int64_t>& channel_nodes,
    int n_dilate, bool remove_non_adjacent, double max_detour_ratio,
    double cell, int n_rows, int n_cols)
{
    const std::size_t n = valley.size();
    const int n_erode = n_dilate;

    // 0. Erode (channel + neighbours protected) to sever connecting
    //    corridors narrower than ~2 * n_erode + 1 pixels.
    std::vector<std::uint8_t> protect(n, k_off);
    for (const std::int64_t cf : channel_nodes)
    {
        if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
        protect[static_cast<std::size_t>(cf)] = k_on;
        const int r = static_cast<int>(cf / n_cols);
        const int c = static_cast<int>(cf % n_cols);
        for (int k = 0; k < 8; ++k)
        {
            const int rr = r + k_dr[k], cc = c + k_dc[k];
            if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) continue;
            protect[static_cast<std::size_t>(rr) * n_cols + cc] = k_on;
        }
    }
    for (int it = 0; it < n_erode; ++it)
    {
        std::vector<std::uint8_t> eroded = valley;
        bool changed = false;
        for (int r = 0; r < n_rows; ++r)
        {
            for (int c = 0; c < n_cols; ++c)
            {
                const std::size_t i =
                    static_cast<std::size_t>(r) * n_cols + c;
                if (!valley[i] || protect[i]) continue;
                for (int k = 0; k < 8; ++k)
                {
                    const int rr = r + k_dr[k], cc = c + k_dc[k];
                    if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols ||
                        !valley[static_cast<std::size_t>(rr) * n_cols + cc])
                    {
                        eroded[i] = k_off;
                        changed = true;
                        break;
                    }
                }
            }
        }
        valley.swap(eroded);
        if (!changed) break;
    }
    n_dilate = 2 * n_dilate;   // regrow the rim + fill holes

    // 1. Label patches.
    std::vector<std::int32_t> lbl(n, -1);
    std::vector<double> patch_max_dist;
    std::deque<std::int64_t> queue;
    for (std::size_t seed = 0; seed < n; ++seed)
    {
        if (!valley[seed] || lbl[seed] >= 0) continue;
        const std::int32_t id =
            static_cast<std::int32_t>(patch_max_dist.size());
        patch_max_dist.push_back(0.0);
        lbl[seed] = id;
        queue.push_back(static_cast<std::int64_t>(seed));
        while (!queue.empty())
        {
            const std::int64_t flat = queue.front();
            queue.pop_front();
            const std::size_t fi = static_cast<std::size_t>(flat);
            if (dist_mask[fi] && std::isfinite(dist_to_channel[fi]))
                patch_max_dist[id] = std::max(
                    patch_max_dist[id],
                    static_cast<double>(dist_to_channel[fi]));
            const int r = static_cast<int>(flat / n_cols);
            const int c = static_cast<int>(flat % n_cols);
            for (int k = 0; k < 8; ++k)
            {
                const int rr = r + k_dr[k], cc = c + k_dc[k];
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                    continue;
                const std::size_t nf =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (valley[nf] && lbl[nf] < 0)
                {
                    lbl[nf] = id;
                    queue.push_back(static_cast<std::int64_t>(nf));
                }
            }
        }
    }

    // 2. Remove patches not adjacent to the channel (before regrowth, so a
    //    severed far-off patch cannot re-bridge during dilation).
    if (remove_non_adjacent && !patch_max_dist.empty())
    {
        std::vector<char> adjacent(patch_max_dist.size(), 0);
        for (const std::int64_t cf : channel_nodes)
        {
            const int r = static_cast<int>(cf / n_cols);
            const int c = static_cast<int>(cf % n_cols);
            for (int k = -1; k < 8; ++k)
            {
                const int rr = (k < 0) ? r : r + k_dr[k];
                const int cc = (k < 0) ? c : c + k_dc[k];
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                    continue;
                const std::size_t nf =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (valley[nf] && lbl[nf] >= 0) adjacent[lbl[nf]] = 1;
            }
        }
        for (std::size_t i = 0; i < n; ++i)
            if (valley[i] && (lbl[i] < 0 || !adjacent[lbl[i]]))
            {
                valley[i] = k_off;
                lbl[i] = -1;
            }
    }

    // 2b. Detour-ratio pruning: geodesic distance from the channel THROUGH
    //     the valley mask (Dijkstra, 8-connected) must not exceed
    //     max_detour_ratio times the Euclidean distance to the channel.
    //     Genuine floodplain is directly reachable (ratio ~1); a
    //     neighbouring channel's floodplain hangs on a long corridor detour
    //     and prunes away. A ten-cell floor keeps near-channel pixels behind
    //     small obstructions.
    if (max_detour_ratio > 0.0)
    {
        constexpr float k_unreached = std::numeric_limits<float>::infinity();
        std::vector<float> geo(n, k_unreached);
        using qnode = std::pair<float, std::int64_t>;
        std::priority_queue<qnode, std::vector<qnode>,
                            std::greater<qnode>> pq;
        for (const std::int64_t cf : channel_nodes)
        {
            if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
            const std::size_t fi = static_cast<std::size_t>(cf);
            if (!valley[fi]) continue;
            geo[fi] = 0.0f;
            pq.emplace(0.0f, cf);
        }
        const float straight = static_cast<float>(cell);
        const float diag = static_cast<float>(cell * std::sqrt(2.0));
        while (!pq.empty())
        {
            const auto [d, flat] = pq.top();
            pq.pop();
            const std::size_t fi = static_cast<std::size_t>(flat);
            if (d > geo[fi]) continue;
            const int r = static_cast<int>(flat / n_cols);
            const int c = static_cast<int>(flat % n_cols);
            for (int k = 0; k < 8; ++k)
            {
                const int rr = r + k_dr[k], cc = c + k_dc[k];
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                    continue;
                const std::size_t nf =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (!valley[nf]) continue;
                const float nd =
                    d + ((k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight);
                if (nd < geo[nf])
                {
                    geo[nf] = nd;
                    pq.emplace(nd, static_cast<std::int64_t>(nf));
                }
            }
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!valley[i]) continue;
            if (!std::isfinite(geo[i])) { valley[i] = k_off; continue; }
            const double eu =
                (dist_mask[i] && std::isfinite(dist_to_channel[i]))
                    ? static_cast<double>(dist_to_channel[i]) : 0.0;
            const double allowed =
                std::max(max_detour_ratio * eu, 10.0 * cell);
            if (static_cast<double>(geo[i]) > allowed) valley[i] = k_off;
        }
    }

    // 2c. Line-of-sight pruning: a genuine valley pixel sees its nearest
    //     channel pixel across low ground — the straight line between them
    //     stays inside the valley mask (small gaps tolerated: roads,
    //     embankments). A neighbouring valley's flat that is merely
    //     CONNECTED to this valley (e.g. around a broad glacial col) has
    //     intervening high ground on the direct line and prunes away.
    {
        // Nearest channel pixel per grid cell: multi-source Dijkstra over
        // the full valid grid (grid distance ~ Euclidean), carrying the
        // seed's flat index.
        std::vector<float>        gdist(n,
            std::numeric_limits<float>::infinity());
        std::vector<std::int64_t> gsrc(n, -1);
        using qnode = std::pair<float, std::int64_t>;
        std::priority_queue<qnode, std::vector<qnode>,
                            std::greater<qnode>> pq;
        for (const std::int64_t cf : channel_nodes)
        {
            if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
            gdist[static_cast<std::size_t>(cf)] = 0.0f;
            gsrc[static_cast<std::size_t>(cf)]  = cf;
            pq.emplace(0.0f, cf);
        }
        const float straight = static_cast<float>(cell);
        const float diag = static_cast<float>(cell * std::sqrt(2.0));
        while (!pq.empty())
        {
            const auto [d, flat] = pq.top();
            pq.pop();
            const std::size_t fi = static_cast<std::size_t>(flat);
            if (d > gdist[fi]) continue;
            const int r = static_cast<int>(flat / n_cols);
            const int c = static_cast<int>(flat % n_cols);
            for (int k = 0; k < 8; ++k)
            {
                const int rr = r + k_dr[k], cc = c + k_dc[k];
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                    continue;
                const std::size_t nf =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (!dist_mask[nf]) continue;
                const float nd =
                    d + ((k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight);
                if (nd < gdist[nf])
                {
                    gdist[nf] = nd;
                    gsrc[nf]  = gsrc[fi];
                    pq.emplace(nd, static_cast<std::int64_t>(nf));
                }
            }
        }

        for (std::size_t i = 0; i < n; ++i)
        {
            if (!valley[i] || gsrc[i] < 0) continue;
            const std::int64_t tgt = gsrc[i];
            const int r0 = static_cast<int>(i / n_cols);
            const int c0 = static_cast<int>(i % n_cols);
            const int r1 = static_cast<int>(tgt / n_cols);
            const int c1 = static_cast<int>(tgt % n_cols);
            const int steps = std::max(std::abs(r1 - r0),
                                       std::abs(c1 - c0));
            if (steps <= 2) continue;
            int miss = 0;
            for (int t = 1; t < steps; ++t)
            {
                const int rr = r0 + (r1 - r0) * t / steps;
                const int cc = c0 + (c1 - c0) * t / steps;
                if (!valley[static_cast<std::size_t>(rr) * n_cols + cc])
                    ++miss;
            }
            // Tolerate small gaps: two cells or 10% of the line.
            if (miss > std::max(2, steps / 10)) valley[i] = k_off;
        }
    }

    // 3. Constrained dilation (frontier-based, label-propagating), with a
    //    margin on the cap covering the erosion depth so the rim can regrow.
    const double cap_margin = (n_erode + 0.5) * cell;
    for (int it = 0; it < n_dilate; ++it)
    {
        std::vector<std::pair<std::size_t, std::int32_t>> grown;
        for (int r = 0; r < n_rows; ++r)
        {
            for (int c = 0; c < n_cols; ++c)
            {
                const std::size_t i =
                    static_cast<std::size_t>(r) * n_cols + c;
                if (valley[i]) continue;
                if (!dist_mask[i] || !std::isfinite(dist_to_channel[i]))
                    continue;
                // Grow from the neighbouring patch with the largest reach.
                std::int32_t best = -1;
                double best_max = -1.0;
                for (int k = 0; k < 8; ++k)
                {
                    const int rr = r + k_dr[k], cc = c + k_dc[k];
                    if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                        continue;
                    const std::size_t nf =
                        static_cast<std::size_t>(rr) * n_cols + cc;
                    if (valley[nf] && lbl[nf] >= 0 &&
                        patch_max_dist[lbl[nf]] > best_max)
                    {
                        best_max = patch_max_dist[lbl[nf]];
                        best = lbl[nf];
                    }
                }
                if (best >= 0 &&
                    static_cast<double>(dist_to_channel[i]) <=
                        best_max + cap_margin)
                    grown.emplace_back(i, best);
            }
        }
        if (grown.empty()) break;
        for (const auto& [i, id] : grown)
        {
            valley[i] = k_on;
            lbl[i] = id;
        }
    }
}

// Valley break from a corridor HANC histogram (log10 bins over [lo_log,
// hi_log]). Corridor populations differ from flat-mask populations: the
// valley-floor mode sits AT the low edge (values clamp at the log floor) and
// hillslopes form a continuum, so the generic antimode finder fails. Rule
// used here: identify the dominant wall/upland peak (the largest-HANC peak
// whose density reaches wall_peak_fraction of the histogram maximum) and
// return the LAST local minimum below it whose density dips under
// valley_fraction of that wall peak — i.e. the valley is everything below
// the rise into the valley walls. Returns <= 0 when no such structure exists.
double find_valley_break(
    const std::vector<std::int64_t>& hist,
    double lo_log, double hi_log,
    double valley_fraction, double wall_peak_fraction)
{
    const int nb = static_cast<int>(hist.size());
    if (nb < 8 || hi_log - lo_log < 1e-9) return -1.0;
    const double bw = (hi_log - lo_log) / nb;

    std::vector<double> sm(hist.begin(), hist.end()), tmp(nb, 0.0);
    for (int pass = 0; pass < 2; ++pass)
    {
        for (int b = 0; b < nb; ++b)
        {
            double acc = 0.0; int c = 0;
            for (int k = -2; k <= 2; ++k)
            {
                const int j = b + k;
                if (j >= 0 && j < nb) { acc += sm[j]; ++c; }
            }
            tmp[b] = acc / std::max(c, 1);
        }
        sm.swap(tmp);
    }

    const double gmax = *std::max_element(sm.begin(), sm.end());
    if (gmax <= 0.0) return -1.0;

    // Wall peak: largest-HANC local maximum (edge-aware) with density above
    // wall_peak_fraction * gmax.
    int wall = -1;
    for (int b = nb - 1; b >= 1; --b)
    {
        const bool is_peak =
            (b == nb - 1) ? sm[b] > sm[b - 1]
                          : (sm[b] >= sm[b - 1] && sm[b] > sm[b + 1]);
        if (is_peak && sm[b] >= wall_peak_fraction * gmax) { wall = b; break; }
    }
    if (wall <= 0) return -1.0;

    // Last local minimum below the wall peak with a genuine dip.
    int best = -1;
    for (int b = 1; b < wall; ++b)
    {
        const bool is_min = sm[b] <= sm[b - 1] && sm[b] < sm[b + 1];
        if (is_min && sm[b] <= valley_fraction * sm[wall]) best = b;
    }
    if (best < 0) return -1.0;
    return std::pow(10.0, lo_log + (best + 0.5) * bw);
}

} // namespace


// -----------------------------------------------------------------------------
// fit_pinned_monotone_baseline (public — see header)
// -----------------------------------------------------------------------------
std::vector<float> fit_pinned_monotone_baseline(
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

    // Minimal monotone envelope: each pin is raised to the maximum of all
    // pins at or downstream of it (reverse suffix-max). The result is
    // non-increasing downstream and >= every original pin.
    std::vector<float> pin_adj(pin_idx.size());
    float running = -std::numeric_limits<float>::infinity();
    for (std::size_t p = pin_idx.size(); p-- > 0;)
    {
        running = std::max(running, pin_elevation[pin_idx[p]]);
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

    // Nodata voids (reservoirs, water-surface gaps) can swallow the routed
    // centreline: a path reaching a void shore terminates (nodata = outlet).
    // Label the nodata components once; per channel, any component the traced
    // channel approaches closely gets plugged with a high dome in the trough
    // surface (its tile-border pixels stay nodata so real outlets survive),
    // and the routed path goes around it instead of dying at its shore.
    std::vector<std::int32_t> nd_lbl(n, -1);
    std::int32_t n_nd_components = 0;
    float void_plug_elev = 0.0f;
    {
        std::deque<std::int64_t> queue;
        for (std::size_t seed = 0; seed < n; ++seed)
        {
            if (dmask[seed] || nd_lbl[seed] >= 0) continue;
            const std::int32_t id = n_nd_components++;
            nd_lbl[seed] = id;
            queue.push_back(static_cast<std::int64_t>(seed));
            while (!queue.empty())
            {
                const std::int64_t flat = queue.front();
                queue.pop_front();
                const int r = static_cast<int>(flat / n_cols);
                const int c = static_cast<int>(flat % n_cols);
                const int rr4[4] = { r - 1, r + 1, r, r };
                const int cc4[4] = { c, c, c - 1, c + 1 };
                for (int k = 0; k < 4; ++k)
                {
                    if (rr4[k] < 0 || rr4[k] >= n_rows ||
                        cc4[k] < 0 || cc4[k] >= n_cols) continue;
                    const std::size_t nf =
                        static_cast<std::size_t>(rr4[k]) * n_cols + cc4[k];
                    if (!dmask[nf] && nd_lbl[nf] < 0)
                    {
                        nd_lbl[nf] = id;
                        queue.push_back(static_cast<std::int64_t>(nf));
                    }
                }
            }
        }
        float zmax = -std::numeric_limits<float>::infinity();
        for (std::size_t i = 0; i < n; ++i)
            if (dmask[i]) zmax = std::max(zmax, elev[i]);
        void_plug_elev = zmax + 100.0f;
    }

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
        // 1. Corridor and Euclidean HANC relative to THIS channel only.
        // ---------------------------------------------------------------------
        std::vector<float> chan_data(n, 0.0f);
        for (const std::int64_t f : nodes)
            chan_data[static_cast<std::size_t>(f)] = 1.0f;
        lsdtt_raster_cpu chan_binary(std::move(chan_data),
                                     std::vector<uint8_t>(dmask), meta);

        lsdtt_raster_cpu dist =
            euclidean_distance_transform(chan_binary, request.n_threads);
        lsdtt_raster_cpu hanc = difference_to_nearest_true_mask_value(
            chan_binary, dem, request.n_threads);

        const std::vector<float>&   dd = dist.get_data();
        const std::vector<uint8_t>& dm = dist.get_mask();
        const std::vector<float>&   hd = hanc.get_data();
        const std::vector<uint8_t>& hm = hanc.get_mask();

        // ---------------------------------------------------------------------
        // 2. Break search in the corridor HANC population.
        // ---------------------------------------------------------------------
        const int    nbins   = std::max(request.hanc_hist_bins, 16);
        const double floor_m = std::max(request.min_hanc_metres, 1e-3);
        const double lo_log  = std::log10(floor_m);
        double hi_log = lo_log;

        auto in_corridor = [&](std::size_t i) -> bool
        {
            return dmask[i] && dm[i] && hm[i] &&
                   std::isfinite(dd[i]) && std::isfinite(hd[i]) &&
                   static_cast<double>(dd[i]) <= request.corridor_width_metres;
        };

        for (std::size_t i = 0; i < n; ++i)
            if (in_corridor(i))
                hi_log = std::max(hi_log,
                    std::log10(std::max<double>(hd[i], floor_m)));

        // Along-channel position of every corridor pixel: inflate each
        // channel node's cumulative along-channel distance across the
        // corridor (nearest-node attribution), so breaks can adapt from
        // confined upland reaches to broad valley bottoms.
        std::vector<float> node_along(nodes.size(), 0.0f);
        for (std::size_t k = 1; k < nodes.size(); ++k)
        {
            const std::int64_t a = nodes[k - 1], b = nodes[k];
            const int dr2 = static_cast<int>(a / n_cols) -
                            static_cast<int>(b / n_cols);
            const int dc2 = static_cast<int>(a % n_cols) -
                            static_cast<int>(b % n_cols);
            node_along[k] = node_along[k - 1] + static_cast<float>(
                cell * std::sqrt(static_cast<double>(dr2 * dr2 + dc2 * dc2)));
        }
        const double channel_length =
            static_cast<double>(node_along.back());

        // Inflate the nearest channel NODE INDEX across the corridor, then
        // derive (a) the unsigned along-channel distance for the reach
        // windows and (b) a SIGNED along coordinate for the trough tilt: the
        // node's along-distance plus the projection of the pixel offset onto
        // the local downstream direction. Pixels lateral/behind the channel
        // head get a negative projection, so the tilt raises them and the
        // routed path cannot escape upstream out of the head's pocket.
        std::vector<float> along_data(n, 0.0f);   // unsigned, for windows
        std::vector<float> along_tilt(n, 0.0f);   // signed, for the trough
        std::vector<uint8_t> along_mask(n, 0u);
        {
            std::vector<float>   sv(n, 0.0f);
            std::vector<uint8_t> sm(n, 0u);
            std::vector<float>   pd2(n, 1.0f);
            std::vector<uint8_t> pm2(n, 0u);
            for (std::size_t k = 0; k < nodes.size(); ++k)
            {
                const std::size_t f = static_cast<std::size_t>(nodes[k]);
                sv[f] = static_cast<float>(k);
                sm[f] = 1u;
            }
            for (std::size_t i = 0; i < n; ++i)
                pm2[i] = (in_corridor(i) || sm[i]) ? 1u : 0u;
            lsdtt_raster_cpu inflated = inflate_values_to_nearest_source(
                lsdtt_raster_cpu(std::move(sv), std::move(sm), meta),
                lsdtt_raster_cpu(std::move(pd2), std::move(pm2), meta),
                request.n_threads);
            const std::vector<float>&   idx_d = inflated.get_data();
            const std::vector<uint8_t>& idx_m = inflated.get_mask();

            const int span = 10;   // node window for the local direction
            const int last = static_cast<int>(nodes.size()) - 1;
            for (std::size_t i = 0; i < n; ++i)
            {
                if (!idx_m[i]) continue;
                int k = std::clamp(static_cast<int>(
                    std::llround(idx_d[i])), 0, last);
                along_mask[i] = 1u;
                along_data[i] = node_along[static_cast<std::size_t>(k)];

                const int ka = std::max(k - span, 0);
                const int kb = std::min(k + span, last);
                const std::int64_t fa = nodes[static_cast<std::size_t>(ka)];
                const std::int64_t fb = nodes[static_cast<std::size_t>(kb)];
                const std::int64_t fk = nodes[static_cast<std::size_t>(k)];
                double dx = static_cast<double>(
                    (fb % n_cols) - (fa % n_cols));
                double dy = static_cast<double>(
                    (fb / n_cols) - (fa / n_cols));
                const double len = std::hypot(dx, dy);
                double proj = 0.0;
                if (len > 0.0)
                {
                    dx /= len; dy /= len;
                    const double px2 = static_cast<double>(
                        static_cast<std::int64_t>(i % n_cols) -
                        (fk % n_cols));
                    const double py2 = static_cast<double>(
                        static_cast<std::int64_t>(i / n_cols) -
                        (fk / n_cols));
                    proj = (px2 * dx + py2 * dy) * cell;
                }
                along_tilt[i] = along_data[i] + static_cast<float>(proj);
            }
        }

        double wstep = request.break_reach_window_metres;
        bool reach_adaptive = wstep >= 0.0;
        if (wstep == 0.0)
            wstep = std::max(2000.0, channel_length / 8.0);
        const int n_windows = reach_adaptive
            ? std::max(1, static_cast<int>(channel_length / wstep) + 1)
            : 1;
        auto window_of = [&](std::size_t i) -> int
        {
            if (!reach_adaptive || !along_mask[i]) return 0;
            return std::clamp(
                static_cast<int>(along_data[i] / wstep), 0, n_windows - 1);
        };

        // Whole-corridor + per-window histograms.
        double break_hanc = request.fallback_valley_hanc_metres;
        rec.used_fallback_break = true;
        std::vector<double> window_break(n_windows, -1.0);
        if (hi_log > lo_log + 1e-6)
        {
            hi_log += 1e-6;
            std::vector<std::int64_t> ghist(nbins, 0);
            std::vector<std::int64_t> whist(
                static_cast<std::size_t>(n_windows) * nbins, 0);
            std::vector<std::int64_t> wcount(n_windows, 0);
            for (std::size_t i = 0; i < n; ++i)
            {
                if (!in_corridor(i)) continue;
                const double l =
                    std::log10(std::max<double>(hd[i], floor_m));
                const int b = std::clamp(static_cast<int>(
                    (l - lo_log) / (hi_log - lo_log) * nbins), 0, nbins - 1);
                ++ghist[b];
                const int w = window_of(i);
                ++whist[static_cast<std::size_t>(w) * nbins + b];
                ++wcount[w];
            }

            const double gbrk = find_valley_break(
                ghist, lo_log, hi_log,
                request.break_valley_fraction, request.wall_peak_fraction);
            if (gbrk > 0.0)
            {
                break_hanc =
                    std::max(gbrk, request.min_valley_hanc_metres);
                rec.used_fallback_break = false;
            }

            std::vector<char> own(n_windows, 0);
            for (int w = 0; w < n_windows; ++w)
            {
                // Pool each window with its neighbours (sliding window).
                std::vector<std::int64_t> h(nbins, 0);
                std::int64_t cnt = 0;
                for (int k = std::max(w - 1, 0);
                     k <= std::min(w + 1, n_windows - 1); ++k)
                {
                    cnt += wcount[k];
                    const std::int64_t* src =
                        &whist[static_cast<std::size_t>(k) * nbins];
                    for (int b = 0; b < nbins; ++b) h[b] += src[b];
                }
                if (cnt >= request.break_min_window_pixels)
                {
                    const double brk = find_valley_break(
                        h, lo_log, hi_log,
                        request.break_valley_fraction,
                        request.wall_peak_fraction);
                    if (brk > 0.0)
                    {
                        window_break[w] =
                            std::max(brk, request.min_valley_hanc_metres);
                        own[w] = 1;
                    }
                }
            }
            // Windows without their own break inherit the nearest window's
            // break, then the whole-corridor break.
            for (int w = 0; w < n_windows; ++w)
            {
                if (own[w]) continue;
                int best = -1;
                for (int k = 0; k < n_windows; ++k)
                    if (own[k] && (best < 0 ||
                                   std::abs(k - w) < std::abs(best - w)))
                        best = k;
                window_break[w] =
                    (best >= 0) ? window_break[best] : break_hanc;
            }
            // Median-of-three across windows: a single window whose own
            // break detection locked onto the wrong population mode (e.g. a
            // terrace level) otherwise admits a huge block of false valley;
            // genuine confined -> broad transitions survive the median.
            if (n_windows >= 3)
            {
                std::vector<double> smoothed(window_break);
                for (int w = 1; w + 1 < n_windows; ++w)
                {
                    double a = window_break[w - 1];
                    double b = window_break[w];
                    double c = window_break[w + 1];
                    smoothed[w] = std::max(std::min(a, b),
                                           std::min(std::max(a, b), c));
                }
                window_break.swap(smoothed);
            }
            for (int w = 0; w < n_windows; ++w)
            {
                if (request.max_valley_hanc_metres > 0.0)
                    window_break[w] = std::min(
                        window_break[w], request.max_valley_hanc_metres);
                rec.window_centre_m.push_back((w + 0.5) * wstep);
                rec.window_break_m.push_back(window_break[w]);
                rec.window_break_own.push_back(own[w]);
            }

            rec.hanc_hist   = std::move(ghist);
            rec.hist_lo_log = lo_log;
            rec.hist_hi_log = hi_log;
        }
        else
        {
            window_break.assign(n_windows, break_hanc);
        }
        if (request.max_valley_hanc_metres > 0.0)
            break_hanc = std::min(break_hanc, request.max_valley_hanc_metres);
        rec.valley_hanc_break_metres = break_hanc;

        // ---------------------------------------------------------------------
        // 3. Valley mask: below the reach's break, connected to the channel,
        //    holes filled.
        // ---------------------------------------------------------------------
        std::vector<std::uint8_t> valley(n, k_off);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!in_corridor(i)) continue;
            const double brk_i = window_break[window_of(i)];
            if (static_cast<double>(hd[i]) <= brk_i)
                valley[i] = k_on;
        }
        for (const std::int64_t f : nodes)
            if (dmask[static_cast<std::size_t>(f)])
                valley[static_cast<std::size_t>(f)] = k_on;

        clean_valley_patches(valley, dd, dm, nodes,
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
        for (std::size_t i = 0; i < n; ++i)
            if (hole[i] && along_mask[i] &&
                along_tilt[i] < -2.0f * static_cast<float>(cell))
                hole[i] = k_off;
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

        // Nearest-wall elevation, inflated across the hole.
        std::vector<float>   sv(n, 0.0f);
        std::vector<uint8_t> sm(n, 0u);
        std::vector<float>   pd(n, 1.0f);
        std::vector<uint8_t> pm(n, 0u);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!dmask[i]) continue;
            pm[i] = 1u;
            if (!hole[i]) { sv[i] = elev[i]; sm[i] = 1u; }
        }
        lsdtt_raster_cpu wall_elev = inflate_values_to_nearest_source(
            lsdtt_raster_cpu(std::move(sv), std::move(sm), meta),
            lsdtt_raster_cpu(std::move(pd), std::move(pm), meta),
            request.n_threads);
        const std::vector<float>&   we = wall_elev.get_data();
        const std::vector<uint8_t>& wem = wall_elev.get_mask();

        // Circular local-minimum filter -> "minimum bank elevation".
        const int r_px = std::max(1, static_cast<int>(std::lround(
            request.min_bank_elevation_window_radius_metres / cell)));
        std::vector<float> min_bank(n, 0.0f);
        for (int r = 0; r < n_rows; ++r)
        {
            for (int c = 0; c < n_cols; ++c)
            {
                const std::size_t i =
                    static_cast<std::size_t>(r) * n_cols + c;
                if (!hole[i] || !wem[i]) continue;
                float mn = we[i];
                for (int dr2 = -r_px; dr2 <= r_px; ++dr2)
                {
                    const int rr = r + dr2;
                    if (rr < 0 || rr >= n_rows) continue;
                    const int half = static_cast<int>(std::sqrt(
                        static_cast<double>(r_px) * r_px -
                        static_cast<double>(dr2) * dr2));
                    for (int dc2 = -half; dc2 <= half; ++dc2)
                    {
                        const int cc = c + dc2;
                        if (cc < 0 || cc >= n_cols) continue;
                        const std::size_t j =
                            static_cast<std::size_t>(rr) * n_cols + cc;
                        if (wem[j]) mn = std::min(mn, we[j]);
                    }
                }
                min_bank[i] = mn;
            }
        }

        // ---------------------------------------------------------------------
        // 5. Route on the trough; retry with a deeper trough on failure.
        // ---------------------------------------------------------------------
        // REACTIVE void plugging. Proactive plugging of nodata regions is
        // unsafe: it can seal the trough's only natural exit (the sea / a
        // collar contact beyond the outlet) and the whole filled bowl then
        // drains backwards out the upstream end. Plug nothing up front; when
        // a routing attempt terminates at a void that is not the channel's
        // terminus, plug THAT component and retry.
        std::vector<char> plug_component(
            static_cast<std::size_t>(std::max(n_nd_components, 1)), 0);

        double k_scale = request.trough_scaling_factor;
        double tilt    = request.trough_along_channel_slope;
        for (int attempt = 0; attempt <= request.max_trough_retries; ++attempt)
        {
            std::vector<float>   td(elev);
            std::vector<uint8_t> tm(dmask);
            for (std::size_t i = 0; i < n; ++i)
            {
                if (hole[i] && wem[i])
                    td[i] = min_bank[i] -
                            static_cast<float>(k_scale) * wd[i] -
                            static_cast<float>(tilt *
                                (along_mask[i] ? along_tilt[i] : 0.0f));
                else if (!dmask[i] && nd_lbl[i] >= 0 &&
                         plug_component[nd_lbl[i]])
                {
                    // Plug reservoirs / water-surface voids the channel runs
                    // into so the routed path goes around, not into, them —
                    // but only nodata DEEPER than ~300 m from the tile
                    // border: the outer collar ribbon stays nodata so real
                    // outlets at the tile edge survive.
                    const int r = static_cast<int>(i / n_cols);
                    const int c = static_cast<int>(i % n_cols);
                    const double border_dist = cell * std::min(
                        std::min(r, n_rows - 1 - r),
                        std::min(c, n_cols - 1 - c));
                    if (border_dist > 300.0)
                    {
                        td[i] = void_plug_elev;
                        tm[i] = 1u;
                    }
                }
            }
            lsdtt_raster_cpu trough(std::move(td), std::move(tm), meta);

            lsdtt_raster_cpu routed = trough;
            for (int loop = 0; loop < request.centreline_loops; ++loop)
            {
                std::vector<float>& rd = routed.get_data();
                for (std::size_t i = 0; i < n; ++i)
                    if (hole[i] && wem[i])
                        rd[i] -= static_cast<float>(k_scale) * wd[i];
                // NOTE: no breaching. breach_depressions carves from a pit
                // to the NEAREST low edge, which on a through-valley (open at
                // both tile ends) can trench the trough out the UPSTREAM end
                // and defeat the downstream tilt. Filling alone raises pits
                // to their spill level, and with the tilt the spill is always
                // downstream.
                routed = fill_pits(routed,
                                   static_cast<float>(request.fill_min_slope),
                                   request.n_threads);
            }

            lsdtt_flowinfo_options fopts;
            fopts.build_topological_order = false;
            fopts.build_s_vector          = false;
            fopts.n_threads               = request.n_threads;
            std::unique_ptr<lsdtt_flow_base> tflow;
            try
            {
                tflow = make_flow_base(routed.get_data_array(),
                                       routed.get_mask_array(),
                                       routed.get_meta(), fopts);
            }
            catch (const std::exception&)
            {
                tflow.reset();
            }

            std::vector<std::int64_t> path;
            if (tflow && !tflow->empty())
            {
                std::int64_t flat = nodes.front();   // channel head
                std::size_t guard = 0;
                while (tflow->is_valid_cell(flat) && guard++ < n)
                {
                    path.push_back(flat);
                    if (tflow->is_outlet_unchecked(flat)) break;
                    const std::int64_t rcv =
                        tflow->receiver_of_unchecked(flat);
                    if (rcv == flat || rcv < 0) break;
                    flat = rcv;
                }
            }

            // Success requires a long-enough path AND, for an outlet-bounded
            // channel, that the path actually approaches its outlet. A path
            // that drained out the wrong side of the valley (e.g. the nearest
            // tile edge) is a routing failure and triggers a retry with a
            // deeper, more strongly tilted trough.
            bool ok = static_cast<int>(path.size()) >=
                      request.min_centreline_nodes;
            std::size_t truncate_at = path.empty() ? 0 : path.size() - 1;
            const bool outlet_bounded =
                s < scope.resolved_outlets.size() &&
                scope.resolved_outlets[s] >= 0;
            // Generous tolerance: half the corridor width.
            const double tol = 0.5 * request.corridor_width_metres / cell;
            if (ok && outlet_bounded)
            {
                const std::int64_t of = scope.resolved_outlets[s];
                const int orow = static_cast<int>(of / n_cols);
                const int ocol = static_cast<int>(of % n_cols);
                double best_d2 = std::numeric_limits<double>::max();
                for (std::size_t p = 0; p < path.size(); ++p)
                {
                    const int pr = static_cast<int>(path[p] / n_cols);
                    const int pc = static_cast<int>(path[p] % n_cols);
                    const double d2 =
                        static_cast<double>(pr - orow) * (pr - orow) +
                        static_cast<double>(pc - ocol) * (pc - ocol);
                    if (d2 < best_d2) { best_d2 = d2; truncate_at = p; }
                }
                if (best_d2 > tol * tol) ok = false;
            }
            else if (ok)
            {
                // Unbounded channel: the routed path must finish near the
                // traced channel's own end. A path that terminated mid-way
                // (e.g. at a reservoir void) is a failure — the reactive
                // plugging below then closes that void and retries.
                const std::int64_t tf = nodes.back();
                const int trow = static_cast<int>(tf / n_cols);
                const int tcol = static_cast<int>(tf % n_cols);
                const int pr = static_cast<int>(path.back() / n_cols);
                const int pc = static_cast<int>(path.back() % n_cols);
                const double d2 =
                    static_cast<double>(pr - trow) * (pr - trow) +
                    static_cast<double>(pc - tcol) * (pc - tcol);
                if (d2 > tol * tol) ok = false;
            }

            if (ok)
            {
                path.resize(truncate_at + 1);
                const std::vector<float>& rd2 = routed.get_data();
                rec.centreline_trough_elevation.reserve(path.size());
                for (const std::int64_t f : path)
                    rec.centreline_trough_elevation.push_back(
                        rd2[static_cast<std::size_t>(f)]);
                rec.centreline_flats = std::move(path);
                rec.routing_succeeded = true;
                break;
            }

            rec.trough_retries = attempt + 1;

            // Reactive void plugging: if this attempt died at a nodata shore
            // (reservoir, water-surface gap), plug that component and retry
            // with the same trough parameters. Otherwise deepen and tilt the
            // trough harder.
            bool plugged_new = false;
            if (!path.empty())
            {
                const std::int64_t ef = path.back();
                const int er = static_cast<int>(ef / n_cols);
                const int ec = static_cast<int>(ef % n_cols);
                for (int dr2 = -3; dr2 <= 3; ++dr2)
                {
                    const int rr = er + dr2;
                    if (rr < 0 || rr >= n_rows) continue;
                    for (int dc2 = -3; dc2 <= 3; ++dc2)
                    {
                        const int cc = ec + dc2;
                        if (cc < 0 || cc >= n_cols) continue;
                        const std::size_t nf =
                            static_cast<std::size_t>(rr) * n_cols + cc;
                        if (nd_lbl[nf] >= 0 && !plug_component[nd_lbl[nf]])
                        {
                            plug_component[nd_lbl[nf]] = 1;
                            plugged_new = true;
                        }
                    }
                }
            }
            if (!plugged_new)
            {
                k_scale *= 2.0;
                tilt    *= 2.0;
            }

            // Keep the final attempt's path as a flagged diagnostic so a
            // failure can be inspected on the map.
            if (attempt == request.max_trough_retries && !path.empty())
            {
                const std::vector<float>& rd2 = routed.get_data();
                rec.centreline_trough_elevation.reserve(path.size());
                for (const std::int64_t f : path)
                    rec.centreline_trough_elevation.push_back(
                        rd2[static_cast<std::size_t>(f)]);
                rec.centreline_flats = std::move(path);
            }
        }

        // ---------------------------------------------------------------------
        // 6. Pinned monotone baseline along the centreline. The centreline,
        //    resampled on the pre-carve topography, rides up and down risers
        //    and onto higher surfaces; the baseline pins to the original
        //    channel elevations wherever the channel crosses (or touches)
        //    the centreline and interpolates linearly between pins, monotone
        //    downstream (see fit_pinned_monotone_baseline).
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
            // still pins the baseline at its along-valley position — this is
            // what strips the meander noise. Where several channel nodes
            // project to one centreline node, the HIGHEST elevation wins, so
            // the baseline never dips below the channel at a crossing.
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
                if (!std::isfinite(pins[best]) || ce > pins[best])
                    pins[best] = ce;
                rec.centreline_is_pin[best] = 1;
            }
            rec.centreline_baseline_elevation =
                fit_pinned_monotone_baseline(along_path, pins);
        }

        // ---------------------------------------------------------------------
        // 7. Baseline refinement of the valley mask (see request comment):
        //    valley = (elevation - baseline of the nearest centreline node)
        //    <= the reach's break. Kills the stale-reference wedges that the
        //    nearest-channel-pixel test grows onto downstream-sloping
        //    terraces wherever the channel veers away.
        // ---------------------------------------------------------------------
        if (request.refine_valley_with_baseline &&
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
                    if (rel <= window_break[window_of(i)])
                        refined[i] = k_on;
                }
                // The channel itself is always valley (and anchors the
                // patch-adjacency filter).
                for (const std::int64_t f : nodes)
                    if (dmask[static_cast<std::size_t>(f)])
                        refined[static_cast<std::size_t>(f)] = k_on;

                clean_valley_patches(refined, dd, dm, nodes,
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

        result.channels.push_back(std::move(rec));
    }

    result.valley_mask =
        lsdtt_raster_cpu(std::move(vm_data), std::move(vm_mask), meta);
    return result;
}

} // namespace lsdtt
