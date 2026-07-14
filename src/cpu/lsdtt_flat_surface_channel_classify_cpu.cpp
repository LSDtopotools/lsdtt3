// lsdtt_flat_surface_channel_classify_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_flat_surface_channel_classify_cpu
//
// Implementation of floodplain / terrace classification. See header for the
// classification rules, the leak backstop, and the reach-relative HANC banding.
// -----------------------------------------------------------------------------

#include "lsdtt_flat_surface_channel_classify_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "lsdtt_proximity_cpu.hpp"

namespace lsdtt
{
namespace
{

struct label_accumulator
{
    std::vector<float>                  heights;    // finite HANC values
    std::vector<float>                  donor_dist; // donor along-channel distances
    double                              min_dist = std::numeric_limits<double>::max();
    std::int64_t                        n_pixels = 0;
    std::unordered_map<std::int32_t, std::int64_t> source_votes;
};

struct hanc_leak_split
{
    bool   is_leak = false;
    double split   = 0.0;   // antimode HANC (m) separating floodplain / terrace
};

// Detect a label whose HANC straddles the floodplain and terrace bands and, if
// so, return the antimode HANC at which to split it. See header (leak backstop).
//
// A leak requires: substantial mass at/below the floodplain band AND at/above
// the terrace band; a low-HANC mode and a high-HANC mode separated by at least
// min_mode_gap; and a genuine valley (antimode) between them (the valley count
// must be well below the smaller of the two modal peaks).
hanc_leak_split detect_hanc_leak(
    const std::vector<float>& heights,
    double fp_max, double terr_min,
    double min_side_fraction, double min_mode_gap, int nbins)
{
    hanc_leak_split out;
    const std::size_t n = heights.size();
    if (n < 40) return out;
    if (nbins < 8) nbins = 8;

    std::vector<float> s(heights);
    std::sort(s.begin(), s.end());
    auto quantile = [&](double f) -> double {
        std::size_t idx = static_cast<std::size_t>(
            f * static_cast<double>(s.size() - 1) + 0.5);
        if (idx >= s.size()) idx = s.size() - 1;
        return static_cast<double>(s[idx]);
    };
    const double lo = quantile(0.02);
    const double hi = quantile(0.98);
    if (hi - lo < 1e-6) return out;
    const double bw = (hi - lo) / nbins;

    // Mass fractions on each side of the class bands (over all finite HANC).
    double below = 0.0, above = 0.0;
    for (const float v : s)
    {
        if (v <= fp_max)   below += 1.0;
        if (v >= terr_min) above += 1.0;
    }
    const double tot = static_cast<double>(s.size());
    if (below / tot < min_side_fraction || above / tot < min_side_fraction)
        return out;

    // Smoothed histogram over [lo, hi].
    std::vector<double> cnt(nbins, 0.0);
    for (const float v : s)
    {
        if (v < lo || v > hi) continue;
        int b = static_cast<int>((v - lo) / bw);
        if (b >= nbins) b = nbins - 1;
        if (b < 0) b = 0;
        cnt[b] += 1.0;
    }
    std::vector<double> sm(nbins, 0.0);
    for (int b = 0; b < nbins; ++b)
    {
        double acc = 0.0; int c = 0;
        for (int k = -1; k <= 1; ++k)
        {
            const int j = b + k;
            if (j >= 0 && j < nbins) { acc += cnt[j]; ++c; }
        }
        sm[b] = acc / std::max(c, 1);
    }
    auto center = [&](int b) { return lo + (b + 0.5) * bw; };

    // Low mode: highest bin whose centre sits in the floodplain band; high mode:
    // highest bin in the terrace band. Bands can touch (fp_max == terr_min).
    int mA = -1, mB = -1; double pA = -1.0, pB = -1.0;
    for (int b = 0; b < nbins; ++b)
    {
        const double c = center(b);
        if (c <= std::max(fp_max, lo + bw)   && sm[b] > pA) { pA = sm[b]; mA = b; }
        if (c >= std::min(terr_min, hi - bw) && sm[b] > pB) { pB = sm[b]; mB = b; }
    }
    if (mA < 0 || mB < 0 || mB <= mA) return out;

    // Antimode: minimum-count bin between the two modes.
    int vb = mA; double vmin = sm[mA];
    for (int b = mA; b <= mB; ++b)
        if (sm[b] < vmin) { vmin = sm[b]; vb = b; }

    const double gap = (min_mode_gap > 0.0)
        ? min_mode_gap : std::max(2.0, terr_min - fp_max);
    if (center(mB) - center(mA) < gap) return out;      // modes too close
    if (vmin > 0.6 * std::min(pA, pB)) return out;       // valley not a real dip
    out.is_leak = true;
    out.split   = center(vb);
    return out;
}

} // namespace

// -----------------------------------------------------------------------------
// find_population_breaks (public — see header; shared with the valley engine)
// -----------------------------------------------------------------------------
std::vector<double> find_population_breaks(
    const std::vector<std::int64_t>& hist,
    double lo_log, double hi_log, double valley_fraction)
{
    std::vector<double> breaks;
    const int nb = static_cast<int>(hist.size());
    if (nb < 8 || hi_log - lo_log < 1e-9) return breaks;
    const double bw = (hi_log - lo_log) / nb;

    // Two moving-average passes (radius 2) approximate a Gaussian smooth.
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

    const double peak_floor =
        0.08 * *std::max_element(sm.begin(), sm.end());
    if (peak_floor <= 0.0) return breaks;

    // Local maxima above the noise floor.
    std::vector<int> maxima;
    for (int b = 1; b + 1 < nb; ++b)
        if (sm[b] >= sm[b - 1] && sm[b] > sm[b + 1] && sm[b] > peak_floor)
            maxima.push_back(b);
    if (maxima.size() < 2) return breaks;

    // Antimode between each consecutive pair of accepted maxima.
    for (std::size_t m = 0; m + 1 < maxima.size(); ++m)
    {
        const int a = maxima[m], b2 = maxima[m + 1];
        int vb = a; double vmin = sm[a];
        for (int b = a; b <= b2; ++b)
            if (sm[b] < vmin) { vmin = sm[b]; vb = b; }
        if (vmin <= valley_fraction * std::min(sm[a], sm[b2]))
            breaks.push_back(std::pow(10.0, lo_log + (vb + 0.5) * bw));
    }
    return breaks;
}

namespace
{

// Robust spread (p98 - p02) of a value population.
double robust_spread(const std::vector<float>& values)
{
    if (values.size() < 2) return 0.0;
    std::vector<float> s(values);
    auto q = [&](double f) {
        std::size_t idx = static_cast<std::size_t>(
            f * static_cast<double>(s.size() - 1) + 0.5);
        if (idx >= s.size()) idx = s.size() - 1;
        std::nth_element(s.begin(), s.begin() + idx, s.end());
        return static_cast<double>(s[idx]);
    };
    return q(0.98) - q(0.02);
}

double median_of(std::vector<float>& v)
{
    if (v.empty()) return -1.0;
    const std::size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    return static_cast<double>(v[mid]);
}

// Per-source along-channel window layout: windows of width `wstep` starting at
// dmin; break detection for window w pools hists of w-1, w, w+1 (sliding
// window of ~1.5x the requested length centred on w).
struct source_band_bins
{
    double dmin  = 0.0;
    double wstep = 1.0;
    int    n_windows = 0;
    // n_windows x nbins histograms, flattened.
    std::vector<std::int64_t> hists;
    std::vector<std::int64_t> counts;   // per window
};

} // namespace


flat_surface_classify_result classify_flat_surfaces_by_channel(
    const lsdtt_raster_cpu&              labelled_flat_mask,
    const lsdtt_raster_cpu&              height_above_channel,
    const std::vector<std::uint8_t>&     channel_mask,
    const std::vector<std::int32_t>&     channel_source_id,
    const flat_surface_classify_request& request,
    const std::vector<std::int64_t>*     donor_flat,
    const lsdtt_raster_cpu*              channel_flow_distance)
{
    const lsdtt_raster_meta& meta = labelled_flat_mask.get_meta();
    const std::vector<float>&   lbl  = labelled_flat_mask.get_data();
    const std::vector<uint8_t>& lmask = labelled_flat_mask.get_mask();
    const std::vector<float>&   hanc = height_above_channel.get_data();
    const std::vector<uint8_t>& hmask = height_above_channel.get_mask();

    const std::size_t n = lbl.size();
    if (hanc.size() != n || channel_mask.size() != n ||
        channel_source_id.size() != n)
    {
        throw std::invalid_argument(
            "classify_flat_surfaces_by_channel: raster / mask size mismatch.");
    }

    bool banding = request.enable_hanc_banding &&
                   donor_flat != nullptr && channel_flow_distance != nullptr;
    if (banding &&
        (donor_flat->size() != n ||
         channel_flow_distance->get_data().size() != n))
    {
        throw std::invalid_argument(
            "classify_flat_surfaces_by_channel: donor_flat / "
            "channel_flow_distance size mismatch.");
    }

    const double cell = meta.cell_size;
    double adj_radius = request.channel_adjacency_radius_metres;
    if (adj_radius <= 0.0) adj_radius = 1.5 * cell;   // 8-neighbour touch

    // -------------------------------------------------------------------------
    // Euclidean distance to the nearest in-scope channel pixel (adjacency), and
    // the nearest in-scope channel's source id inflated over the domain.
    // -------------------------------------------------------------------------
    std::vector<float> dist_data;
    bool have_dist = false;
    {
        std::vector<float>   cd(n, 0.0f);
        std::vector<uint8_t> cm(n, 1);   // full domain valid for the transform
        for (std::size_t i = 0; i < n; ++i)
            cd[i] = channel_mask[i] ? 1.0f : 0.0f;
        try
        {
            lsdtt_raster_cpu chan_binary(std::move(cd), std::move(cm), meta);
            lsdtt_raster_cpu dist =
                euclidean_distance_transform(chan_binary, request.n_threads);
            dist_data = dist.get_data();
            have_dist = true;
        }
        catch (const std::exception&)
        {
            have_dist = false;   // fall back to no adjacency below
        }
    }

    std::vector<float>   src_data;
    std::vector<uint8_t> src_mask;
    bool have_src = false;
    {
        std::vector<float>   sv(n, 0.0f);
        std::vector<uint8_t> sm(n, 0);
        std::vector<float>   pd(n, 1.0f);
        std::vector<uint8_t> pm(n, 0);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (channel_mask[i] && channel_source_id[i] >= 0)
            {
                sv[i] = static_cast<float>(channel_source_id[i]);
                sm[i] = 1;
            }
            // Fill domain must include the channel source pixels themselves plus
            // the surrounding valid footprint, so the source can propagate out to
            // set-back labels. HANC is defined across the drained DEM, so its
            // mask (union with flats and channels) is the natural domain.
            pm[i] = (lmask[i] || channel_mask[i] || hmask[i]) ? 1 : 0;
        }
        try
        {
            lsdtt_raster_cpu source_vals(std::move(sv), std::move(sm), meta);
            lsdtt_raster_cpu priority(std::move(pd), std::move(pm), meta);
            lsdtt_raster_cpu inflated = inflate_values_to_nearest_source(
                source_vals, priority, request.n_threads);
            src_data = inflated.get_data();
            src_mask = inflated.get_mask();
            have_src = true;
        }
        catch (const std::exception&)
        {
            have_src = false;
        }
    }

    // -------------------------------------------------------------------------
    // Reach-relative banding pre-pass (see header). Produces, per flat pixel:
    // band[i] = HANC band index from the population breaks of the pixel's
    // along-channel window (-1 = undefined), and per parent label whether it
    // gets split into per-band sub-labels.
    // -------------------------------------------------------------------------
    flat_surface_classify_result result;

    std::vector<std::int16_t> band;                 // per pixel, -1 undefined
    std::unordered_map<std::int64_t, char> split_label;   // labels to sub-label

    const std::vector<float>* fd_data = nullptr;
    const std::vector<uint8_t>* fd_mask = nullptr;
    if (banding)
    {
        fd_data = &channel_flow_distance->get_data();
        fd_mask = &channel_flow_distance->get_mask();
    }

    // Donor source id + along-channel distance for a pixel; false when the
    // pixel has no usable donor.
    auto donor_reach = [&](std::size_t i, std::int32_t& s, double& d) -> bool
    {
        const std::int64_t df = (*donor_flat)[i];
        if (df < 0) return false;
        const std::size_t dfi = static_cast<std::size_t>(df);
        if (!(*fd_mask)[dfi] || !std::isfinite((*fd_data)[dfi])) return false;
        s = channel_source_id[dfi];
        d = static_cast<double>((*fd_data)[dfi]);
        return true;
    };

    if (banding)
    {
        const int    nbins = std::max(request.band_hist_bins, 16);
        const double floor_m = std::max(request.band_min_hanc_metres, 1e-3);
        const double lo_log = std::log10(floor_m);

        // Population extent: global log range + per-source distance ranges.
        double hi_log = lo_log;
        std::unordered_map<std::int32_t, std::pair<double, double>> drange;
        std::int64_t n_pop = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!lmask[i] || lbl[i] <= 0.5f) continue;
            if (!hmask[i] || !std::isfinite(hanc[i])) continue;
            std::int32_t s; double d;
            if (!donor_reach(i, s, d)) continue;
            hi_log = std::max(hi_log,
                              std::log10(std::max<double>(hanc[i], floor_m)));
            auto it = drange.find(s);
            if (it == drange.end()) drange[s] = {d, d};
            else
            {
                it->second.first  = std::min(it->second.first, d);
                it->second.second = std::max(it->second.second, d);
            }
            ++n_pop;
        }

        if (n_pop >= 100 && hi_log > lo_log + 1e-6)
        {
            hi_log += 1e-6;
            auto log_bin = [&](float v) -> int
            {
                const double l = std::log10(std::max<double>(v, floor_m));
                int b = static_cast<int>((l - lo_log) / (hi_log - lo_log) * nbins);
                return std::clamp(b, 0, nbins - 1);
            };

            // Global histogram + per-source windowed histograms.
            std::vector<std::int64_t> ghist(nbins, 0);
            std::unordered_map<std::int32_t, source_band_bins> bins;
            for (auto& [s, r] : drange)
            {
                source_band_bins sb;
                sb.dmin  = r.first;
                const double span = std::max(r.second - r.first, 1.0);
                double w = request.band_reach_window_metres;
                if (w <= 0.0) w = std::max(2000.0, span / 8.0);
                sb.wstep = std::max(w * 0.5, cell);
                sb.n_windows = static_cast<int>(span / sb.wstep) + 1;
                sb.hists.assign(
                    static_cast<std::size_t>(sb.n_windows) * nbins, 0);
                sb.counts.assign(sb.n_windows, 0);
                bins.emplace(s, std::move(sb));
            }

            auto window_of = [](const source_band_bins& sb, double d) -> int
            {
                const int w = static_cast<int>((d - sb.dmin) / sb.wstep);
                return std::clamp(w, 0, sb.n_windows - 1);
            };

            for (std::size_t i = 0; i < n; ++i)
            {
                if (!lmask[i] || lbl[i] <= 0.5f) continue;
                if (!hmask[i] || !std::isfinite(hanc[i])) continue;
                std::int32_t s; double d;
                if (!donor_reach(i, s, d)) continue;
                const int b = log_bin(hanc[i]);
                ++ghist[b];
                source_band_bins& sb = bins.at(s);
                const int w = window_of(sb, d);
                ++sb.hists[static_cast<std::size_t>(w) * nbins + b];
                ++sb.counts[w];
            }

            result.global_breaks_m = find_population_breaks(
                ghist, lo_log, hi_log, request.band_max_valley_fraction);

            // Per-(source, window) breaks; sparse or break-less windows fall
            // back to the nearest window with breaks on the same source, then
            // to the global set.
            std::unordered_map<std::int32_t,
                               std::vector<std::vector<double>>> window_breaks;
            for (auto& [s, sb] : bins)
            {
                std::vector<std::vector<double>>& wb = window_breaks[s];
                wb.resize(sb.n_windows);
                std::vector<char> has_own(sb.n_windows, 0);
                const std::size_t first_rec = result.band_windows.size();
                for (int w = 0; w < sb.n_windows; ++w)
                {
                    // Pool this window with its neighbours (sliding window).
                    std::vector<std::int64_t> h(nbins, 0);
                    std::int64_t cnt = 0;
                    for (int k = std::max(w - 1, 0);
                         k <= std::min(w + 1, sb.n_windows - 1); ++k)
                    {
                        cnt += sb.counts[k];
                        const std::int64_t* src =
                            &sb.hists[static_cast<std::size_t>(k) * nbins];
                        for (int b = 0; b < nbins; ++b) h[b] += src[b];
                    }

                    hanc_band_window rec;
                    rec.source_id       = s;
                    rec.window_centre_m = sb.dmin + (w + 0.5) * sb.wstep;
                    rec.n_pixels        = cnt;
                    bool own = false;
                    if (cnt >= request.band_min_window_pixels)
                    {
                        std::vector<double> brk = find_population_breaks(
                            h, lo_log, hi_log, request.band_max_valley_fraction);
                        if (!brk.empty()) { wb[w] = std::move(brk); own = true; }
                    }
                    has_own[w] = own ? 1 : 0;
                    rec.used_global = !own;
                    result.band_windows.push_back(std::move(rec));
                }

                // Fallback for windows without their own breaks: inherit the
                // nearest window with breaks on this source (an upland-swamped
                // or sparse reach still gets a sensible local break level),
                // else the global set.
                for (int w = 0; w < sb.n_windows; ++w)
                {
                    if (has_own[w]) continue;
                    int best = -1;
                    for (int k = 0; k < sb.n_windows; ++k)
                        if (has_own[k] &&
                            (best < 0 ||
                             std::abs(k - w) < std::abs(best - w)))
                            best = k;
                    wb[w] = (best >= 0) ? wb[best] : result.global_breaks_m;
                }
                for (int w = 0; w < sb.n_windows; ++w)
                    result.band_windows[first_rec + w].breaks_m = wb[w];
            }

            // Per-pixel band index.
            band.assign(n, -1);
            for (std::size_t i = 0; i < n; ++i)
            {
                if (!lmask[i] || lbl[i] <= 0.5f) continue;
                if (!hmask[i] || !std::isfinite(hanc[i])) continue;
                std::int32_t s; double d;
                if (!donor_reach(i, s, d)) continue;
                const source_band_bins& sb = bins.at(s);
                const std::vector<double>& brk =
                    window_breaks.at(s)[window_of(sb, d)];
                int b = 0;
                for (const double lvl : brk)
                    if (static_cast<double>(hanc[i]) > lvl) ++b;
                band[i] = static_cast<std::int16_t>(std::min(b, 200));
            }
        }
        else
        {
            banding = false;   // population too small / degenerate
        }
    }

    // -------------------------------------------------------------------------
    // Decide which parent labels get band-split: robust HANC spread above the
    // trigger. (A split label whose pixels share one band re-forms unchanged.)
    // -------------------------------------------------------------------------
    if (banding)
    {
        std::unordered_map<std::int64_t, std::vector<float>> label_heights;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!lmask[i]) continue;
            const std::int64_t label =
                static_cast<std::int64_t>(std::llround(lbl[i]));
            if (label <= 0) continue;
            if (hmask[i] && std::isfinite(hanc[i]))
                label_heights[label].push_back(hanc[i]);
        }
        for (auto& [label, h] : label_heights)
        {
            if (robust_spread(h) >
                std::max(request.band_split_min_label_spread_metres, 0.0))
                split_label[label] = 1;
        }
    }

    // -------------------------------------------------------------------------
    // Accumulate per classification unit: the parent label, or (label, band)
    // for split labels. Key = label * 512 + (band + 2); band -1 (undefined)
    // maps to slot 1, non-split labels to slot 0.
    // -------------------------------------------------------------------------
    auto unit_key = [&](std::size_t i, std::int64_t label) -> std::int64_t
    {
        if (!banding || split_label.find(label) == split_label.end())
            return label * 512;
        return label * 512 + 2 + band[i];
    };

    std::map<std::int64_t, label_accumulator> acc;

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!lmask[i]) continue;
        const std::int64_t label = static_cast<std::int64_t>(std::llround(lbl[i]));
        if (label <= 0) continue;

        label_accumulator& a = acc[unit_key(i, label)];
        ++a.n_pixels;

        if (hmask[i] && std::isfinite(hanc[i]))
            a.heights.push_back(hanc[i]);

        if (banding)
        {
            std::int32_t s; double d;
            if (hmask[i] && std::isfinite(hanc[i]) && donor_reach(i, s, d))
                a.donor_dist.push_back(static_cast<float>(d));
        }

        if (have_dist && std::isfinite(dist_data[i]))
            a.min_dist = std::min(a.min_dist, static_cast<double>(dist_data[i]));

        if (have_src && src_mask[i])
        {
            const std::int32_t s =
                static_cast<std::int32_t>(std::llround(src_data[i]));
            if (s >= 0) ++a.source_votes[s];
        }
    }

    // -------------------------------------------------------------------------
    // Classify each unit and build records.
    // -------------------------------------------------------------------------
    result.labels.reserve(acc.size());

    // Per-unit classification decision, consumed by the raster build below.
    // For a leak-split unit the class is decided per pixel from split + the
    // low/high eligibility flags rather than a single class.
    struct label_decision
    {
        surface_class cls     = surface_class::none;  // non-leak units
        bool          is_leak = false;
        double        split   = 0.0;
        bool          fp_ok   = false;   // low-HANC side -> floodplain
        bool          terr_ok = false;   // high-HANC side -> terrace
    };
    std::unordered_map<std::int64_t, label_decision> decision;

    std::unordered_map<std::int64_t, char> band_split_counted;

    for (auto& [key, a] : acc)
    {
        const std::int64_t parent = key / 512;
        const int slot = static_cast<int>(key % 512);

        flat_surface_label_class rec;
        rec.label    = parent;
        rec.n_pixels = a.n_pixels;
        if (slot >= 2)
        {
            rec.hanc_band  = slot - 2;
            rec.band_split = true;
            if (band_split_counted.emplace(parent, 1).second)
                ++result.n_band_split_labels;
        }
        rec.min_distance_to_channel_metres =
            (a.min_dist == std::numeric_limits<double>::max()) ? -1.0 : a.min_dist;
        rec.adjacent_to_channel =
            have_dist && rec.min_distance_to_channel_metres >= 0.0 &&
            rec.min_distance_to_channel_metres <= adj_radius;

        if (!a.heights.empty())
        {
            rec.has_height = true;
            std::vector<float>& h = a.heights;
            const std::size_t mid = h.size() / 2;
            std::nth_element(h.begin(), h.begin() + mid, h.end());
            rec.median_height = h[mid];
            if (h.size() % 2 == 0)
            {
                const float upper = *std::max_element(h.begin(), h.begin() + mid);
                rec.median_height = 0.5 * (rec.median_height + upper);
            }
            double sum = 0.0;
            rec.min_height = h[0];
            rec.max_height = h[0];
            for (const float v : h)
            {
                sum += v;
                rec.min_height = std::min<double>(rec.min_height, v);
                rec.max_height = std::max<double>(rec.max_height, v);
            }
            rec.mean_height = sum / static_cast<double>(h.size());
            rec.hanc_spread_metres = robust_spread(h);
        }

        rec.median_channel_distance_metres = median_of(a.donor_dist);

        std::int32_t best_src = -1;
        std::int64_t best_votes = 0;
        for (const auto& [s, v] : a.source_votes)
            if (v > best_votes) { best_votes = v; best_src = s; }
        rec.source_id = best_src;

        // Terrace eligibility (linked source + inside the valley corridor +
        // sits at one coherent level above its channel).
        const bool corridor_ok =
            request.terrace_max_distance_to_channel_metres <= 0.0 ||
            (rec.min_distance_to_channel_metres >= 0.0 &&
             rec.min_distance_to_channel_metres <=
                 request.terrace_max_distance_to_channel_metres);
        const bool spread_ok =
            request.terrace_max_hanc_spread_metres <= 0.0 ||
            rec.hanc_spread_metres <= request.terrace_max_hanc_spread_metres;
        const bool terr_ok = rec.source_id >= 0 && corridor_ok && spread_ok;

        label_decision d;

        // Leak backstop: does this unit straddle both bands with a bimodal HANC?
        hanc_leak_split leak;
        if (request.enable_leak_split && rec.has_height)
            leak = detect_hanc_leak(
                a.heights,
                request.floodplain_max_height_metres,
                request.terrace_min_height_metres,
                request.leak_split_min_side_fraction,
                request.leak_split_min_mode_gap_metres,
                request.leak_split_hist_bins);

        if (leak.is_leak)
        {
            rec.split_leak         = true;
            rec.split_hanc_metres  = leak.split;
            d.is_leak           = true;
            d.split             = leak.split;
            d.fp_ok             = rec.adjacent_to_channel;  // low side floodplain
            d.terr_ok           = terr_ok;                  // high side terrace
            rec.cls             = surface_class::floodplain;// nominal record class
            ++result.n_leak_split_labels;
            if (d.fp_ok)   ++result.n_floodplain_labels;
            if (d.terr_ok) ++result.n_terrace_labels;
        }
        else if (rec.adjacent_to_channel && rec.has_height &&
                 rec.median_height <= request.floodplain_max_height_metres)
        {
            rec.cls = surface_class::floodplain;
            d.cls   = surface_class::floodplain;
            ++result.n_floodplain_labels;
        }
        else if (rec.has_height && terr_ok &&
                 rec.median_height >= request.terrace_min_height_metres)
        {
            rec.cls = surface_class::terrace;
            d.cls   = surface_class::terrace;
            ++result.n_terrace_labels;
        }
        else
        {
            rec.cls = surface_class::none;
            d.cls   = surface_class::none;
            rec.spread_guard_failed =
                rec.has_height && !spread_ok && rec.source_id >= 0 &&
                corridor_ok &&
                rec.median_height >= request.terrace_min_height_metres;
        }

        decision[key] = d;
        result.labels.push_back(rec);
    }

    // -------------------------------------------------------------------------
    // Build the class raster, floodplain mask, and (when banding) band raster.
    // -------------------------------------------------------------------------
    std::vector<float>   cls_data(n, 0.0f);
    std::vector<uint8_t> cls_mask = lmask;
    std::vector<float>   fp_data(n, 0.0f);
    std::vector<uint8_t> fp_mask = lmask;

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!lmask[i]) continue;
        const std::int64_t label = static_cast<std::int64_t>(std::llround(lbl[i]));
        if (label <= 0) continue;
        auto it = decision.find(unit_key(i, label));
        if (it == decision.end()) continue;
        const label_decision& d = it->second;

        surface_class c;
        if (d.is_leak)
        {
            // Decide per pixel from its HANC relative to the antimode split.
            if (hmask[i] && std::isfinite(hanc[i]))
                c = (hanc[i] <= d.split)
                        ? (d.fp_ok   ? surface_class::floodplain : surface_class::none)
                        : (d.terr_ok ? surface_class::terrace    : surface_class::none);
            else
                c = surface_class::none;
        }
        else
        {
            c = d.cls;
        }

        // Per-pixel valley-corridor clip: an eligible label can still stretch
        // beyond the corridor (it only needs one pixel inside to qualify).
        if (c == surface_class::terrace &&
            request.terrace_corridor_per_pixel &&
            request.terrace_max_distance_to_channel_metres > 0.0 &&
            have_dist &&
            (!std::isfinite(dist_data[i]) ||
             static_cast<double>(dist_data[i]) >
                 request.terrace_max_distance_to_channel_metres))
        {
            c = surface_class::none;
        }

        cls_data[i] = static_cast<float>(static_cast<std::uint8_t>(c));
        if (c == surface_class::floodplain) fp_data[i] = 1.0f;
    }

    result.surface_class_raster =
        lsdtt_raster_cpu(std::move(cls_data), std::move(cls_mask), meta);
    result.floodplain_mask =
        lsdtt_raster_cpu(std::move(fp_data), std::move(fp_mask), meta);

    if (banding)
    {
        // The DEM's nodata can be 0 (a legitimate band index), so the band
        // raster gets its own nodata of -1.
        lsdtt_raster_meta band_meta = meta;
        band_meta.no_data_value = -1.0;
        std::vector<float>   bd(n, -1.0f);
        std::vector<uint8_t> bm(n, 0u);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (band[i] < 0) continue;
            bd[i] = static_cast<float>(band[i]);
            bm[i] = 1u;
        }
        result.hanc_band_raster =
            lsdtt_raster_cpu(std::move(bd), std::move(bm), band_meta);
    }

    return result;
}

} // namespace lsdtt
