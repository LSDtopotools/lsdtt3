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
// classification rules and I/O conventions.
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
    std::vector<float>                  heights;   // finite HAND values
    double                              min_dist = std::numeric_limits<double>::max();
    std::int64_t                        n_pixels = 0;
    std::unordered_map<std::int32_t, std::int64_t> source_votes;
};

struct hand_leak_split
{
    bool   is_leak = false;
    double split   = 0.0;   // antimode HAND (m) separating floodplain / terrace
};

// Detect a label whose HAND straddles the floodplain and terrace bands and, if
// so, return the antimode HAND at which to split it. See header (leak backstop).
//
// A leak requires: substantial mass at/below the floodplain band AND at/above
// the terrace band; a low-HAND mode and a high-HAND mode separated by at least
// min_mode_gap; and a genuine valley (antimode) between them (the valley count
// must be well below the smaller of the two modal peaks).
hand_leak_split detect_hand_leak(
    const std::vector<float>& heights,
    double fp_max, double terr_min,
    double min_side_fraction, double min_mode_gap, int nbins)
{
    hand_leak_split out;
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

    // Mass fractions on each side of the class bands (over all finite HAND).
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


flat_surface_classify_result classify_flat_surfaces_by_channel(
    const lsdtt_raster_cpu&              labelled_flat_mask,
    const lsdtt_raster_cpu&              height_above_channel,
    const std::vector<std::uint8_t>&     channel_mask,
    const std::vector<std::int32_t>&     channel_source_id,
    const flat_surface_classify_request& request)
{
    const lsdtt_raster_meta& meta = labelled_flat_mask.get_meta();
    const std::vector<float>&   lbl  = labelled_flat_mask.get_data();
    const std::vector<uint8_t>& lmask = labelled_flat_mask.get_mask();
    const std::vector<float>&   hand = height_above_channel.get_data();
    const std::vector<uint8_t>& hmask = height_above_channel.get_mask();

    const std::size_t n = lbl.size();
    if (hand.size() != n || channel_mask.size() != n ||
        channel_source_id.size() != n)
    {
        throw std::invalid_argument(
            "classify_flat_surfaces_by_channel: raster / mask size mismatch.");
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
            // set-back labels. HAND is defined across the drained DEM, so its
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
    // Accumulate per-label statistics.
    // -------------------------------------------------------------------------
    std::map<std::int64_t, label_accumulator> acc;

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!lmask[i]) continue;
        const std::int64_t label = static_cast<std::int64_t>(std::llround(lbl[i]));
        if (label <= 0) continue;

        label_accumulator& a = acc[label];
        ++a.n_pixels;

        if (hmask[i] && std::isfinite(hand[i]))
            a.heights.push_back(hand[i]);

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
    // Classify each label and build records.
    // -------------------------------------------------------------------------
    flat_surface_classify_result result;
    result.labels.reserve(acc.size());

    // Per-label classification decision, consumed by the raster build below.
    // For a leak-split label the class is decided per pixel from split + the
    // low/high eligibility flags rather than a single label class.
    struct label_decision
    {
        surface_class cls     = surface_class::none;  // non-leak labels
        bool          is_leak = false;
        double        split   = 0.0;
        bool          fp_ok   = false;   // low-HAND side -> floodplain
        bool          terr_ok = false;   // high-HAND side -> terrace
    };
    std::unordered_map<std::int64_t, label_decision> decision;

    for (auto& [label, a] : acc)
    {
        flat_surface_label_class rec;
        rec.label    = label;
        rec.n_pixels = a.n_pixels;
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
        }

        std::int32_t best_src = -1;
        std::int64_t best_votes = 0;
        for (const auto& [s, v] : a.source_votes)
            if (v > best_votes) { best_votes = v; best_src = s; }
        rec.source_id = best_src;

        // Terrace eligibility (linked source + inside the valley corridor).
        const bool corridor_ok =
            request.terrace_max_distance_to_channel_metres <= 0.0 ||
            (rec.min_distance_to_channel_metres >= 0.0 &&
             rec.min_distance_to_channel_metres <=
                 request.terrace_max_distance_to_channel_metres);
        const bool terr_ok = rec.source_id >= 0 && corridor_ok;

        label_decision d;

        // Leak backstop: does this label straddle both bands with a bimodal HAND?
        hand_leak_split leak;
        if (request.enable_leak_split && rec.has_height)
            leak = detect_hand_leak(
                a.heights,
                request.floodplain_max_height_metres,
                request.terrace_min_height_metres,
                request.leak_split_min_side_fraction,
                request.leak_split_min_mode_gap_metres,
                request.leak_split_hist_bins);

        if (leak.is_leak)
        {
            rec.split_leak         = true;
            rec.split_hand_metres  = leak.split;
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
        }

        decision[label] = d;
        result.labels.push_back(rec);
    }

    // -------------------------------------------------------------------------
    // Build the class raster and floodplain mask.
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
        auto it = decision.find(label);
        if (it == decision.end()) continue;
        const label_decision& d = it->second;

        surface_class c;
        if (d.is_leak)
        {
            // Decide per pixel from its HAND relative to the antimode split.
            if (hmask[i] && std::isfinite(hand[i]))
                c = (hand[i] <= d.split)
                        ? (d.fp_ok   ? surface_class::floodplain : surface_class::none)
                        : (d.terr_ok ? surface_class::terrace    : surface_class::none);
            else
                c = surface_class::none;
        }
        else
        {
            c = d.cls;
        }

        cls_data[i] = static_cast<float>(static_cast<std::uint8_t>(c));
        if (c == surface_class::floodplain) fp_data[i] = 1.0f;
    }

    result.surface_class_raster =
        lsdtt_raster_cpu(std::move(cls_data), std::move(cls_mask), meta);
    result.floodplain_mask =
        lsdtt_raster_cpu(std::move(fp_data), std::move(fp_mask), meta);

    return result;
}

} // namespace lsdtt
