// lsdtt_floodplain_thresholds_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_floodplain_thresholds_cpu.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace lsdtt
{

namespace
{

bool defined(const lsdtt_raster_cpu& r, std::size_t i)
{
    const float v = r.get_data()[i];
    return r.get_mask()[i] != 0 && std::isfinite(v) &&
           !r.get_meta().is_nodata(v);
}

// Labels a 1/0 candidate raster (min size 1, so nothing is dropped here) and
// returns the label per pixel, 0 = background.
std::vector<std::int64_t> label_candidates(const lsdtt_raster_cpu& candidates,
                                           pixel_connectivity      connectivity)
{
    const lsdtt_raster_cpu lab =
        label_binary_components_by_size(candidates, 1, connectivity);
    const std::vector<float>&   d = lab.get_data();
    const std::vector<uint8_t>& m = lab.get_mask();
    std::vector<std::int64_t> out(d.size(), 0);
    for (std::size_t i = 0; i < d.size(); ++i)
        if (m[i] && d[i] > 0.5f) out[i] = static_cast<std::int64_t>(d[i]);
    return out;
}

} // namespace

floodplain_threshold_result classify_floodplain_terraces_by_thresholds(
    const lsdtt_raster_cpu&             admitted,
    const lsdtt_raster_cpu&             hand,
    const lsdtt_raster_cpu&             gradient,
    const std::vector<std::uint8_t>&    channel_mask,
    const lsdtt_raster_cpu*             terrace_field,
    const floodplain_threshold_request& request)
{
    static const std::string fn =
        "classify_floodplain_terraces_by_thresholds: ";
    const lsdtt_raster_meta& meta = admitted.get_meta();
    const std::size_t n = admitted.get_data().size();
    if (hand.get_data().size() != n || gradient.get_data().size() != n ||
        channel_mask.size() != n)
        throw std::invalid_argument(fn + "rasters and channel mask differ in size");
    if (request.find_terraces)
    {
        if (terrace_field == nullptr || terrace_field->get_data().size() != n)
            throw std::invalid_argument(fn + "find_terraces needs a co-registered "
                                             "terrace field");
        if (!(request.min_terrace_height_above_floodplain_metres >= 0.0))
            throw std::invalid_argument(
                fn + "min_terrace_height_above_floodplain_metres must be >= 0: a "
                     "terrace starts above the floodplain bound, never inside it");
        const double start = request.max_floodplain_hand_metres +
                             request.min_terrace_height_above_floodplain_metres;
        if (request.max_terrace_height_metres > 0.0 &&
            request.max_terrace_height_metres <= start)
            throw std::invalid_argument(
                fn + "max_terrace_height_metres (" +
                std::to_string(request.max_terrace_height_metres) +
                ") is at or below where terraces start (" +
                std::to_string(start) + ")");
    }

    const std::vector<float>&   ad = admitted.get_data();
    const std::vector<uint8_t>& am = admitted.get_mask();
    const bool use_gradient = request.max_floodplain_gradient > 0.0;
    const float hand_bound  = static_cast<float>(request.max_floodplain_hand_metres);
    const float grad_bound  = static_cast<float>(request.max_floodplain_gradient);

    auto passes_gradient = [&](std::size_t i)
    {
        if (!use_gradient) return true;
        return defined(gradient, i) && gradient.get_data()[i] < grad_bound;
    };

    floodplain_threshold_result out;

    // ---- floodplain candidates (plus channel seeds) -------------------------
    std::vector<float> cand(n, 0.0f);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!am[i]) continue;
        const bool is_channel = channel_mask[i] != 0;
        if (request.seed_channel_pixels && is_channel)
        {
            cand[i] = 1.0f;
            ++out.n_floodplain_candidates;
            continue;
        }
        if (!(ad[i] > 0.5f)) continue;
        if (!defined(hand, i) || !(hand.get_data()[i] < hand_bound)) continue;
        if (!passes_gradient(i)) continue;
        cand[i] = 1.0f;
        ++out.n_floodplain_candidates;
    }

    std::vector<float> fp(n, 0.0f);
    {
        const lsdtt_raster_cpu cand_r(std::vector<float>(cand),
                                      std::vector<uint8_t>(am), meta);
        const std::vector<std::int64_t> lab =
            label_candidates(cand_r, request.connectivity);

        std::unordered_map<std::int64_t, std::int64_t> size;
        std::unordered_map<std::int64_t, bool>         touches;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (lab[i] == 0) continue;
            ++size[lab[i]];
            if (channel_mask[i]) touches[lab[i]] = true;
        }
        std::unordered_map<std::int64_t, bool> keep;
        for (const auto& [l, sz] : size)
        {
            if (request.min_component_pixels > 1 &&
                sz < request.min_component_pixels)
            {
                ++out.n_floodplain_dropped_small;
                keep[l] = false;
            }
            else if (request.require_channel_contact && !touches.count(l))
            {
                ++out.n_floodplain_dropped_no_channel;
                keep[l] = false;
            }
            else
            {
                ++out.n_floodplain_components_kept;
                keep[l] = true;
            }
        }
        for (std::size_t i = 0; i < n; ++i)
            if (lab[i] != 0 && keep[lab[i]]) fp[i] = 1.0f;
    }
    lsdtt_raster_cpu fp_r(std::move(fp), std::vector<uint8_t>(am), meta);
    if (request.fill_holes) fp_r = binary_fill_holes(fp_r);
    for (const float v : fp_r.get_data())
        if (v > 0.5f) ++out.n_floodplain_pixels;

    // ---- terraces -------------------------------------------------------------
    std::vector<float> tr(n, 0.0f);
    if (request.find_terraces)
    {
        const float start = static_cast<float>(
            request.max_floodplain_hand_metres +
            request.min_terrace_height_above_floodplain_metres);
        const bool  capped  = request.max_terrace_height_metres > 0.0;
        const float ceiling = static_cast<float>(request.max_terrace_height_metres);
        const std::vector<float>& fpd = fp_r.get_data();

        std::vector<float> tcand(n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!am[i] || !(ad[i] > 0.5f) || fpd[i] > 0.5f) continue;
            if (!defined(*terrace_field, i)) continue;
            const float t = terrace_field->get_data()[i];
            if (!(t > start)) continue;
            if (capped && t > ceiling) continue;
            if (!passes_gradient(i)) continue;
            tcand[i] = 1.0f;
        }
        const lsdtt_raster_cpu tcand_r(std::move(tcand),
                                       std::vector<uint8_t>(am), meta);
        const int min_px = std::max(1, request.min_component_pixels);
        const lsdtt_raster_cpu kept =
            filter_binary_components_by_size(tcand_r, min_px, request.connectivity);
        tr = kept.get_data();
        const std::vector<std::int64_t> lab =
            label_candidates(kept, request.connectivity);
        std::int64_t max_label = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (tr[i] > 0.5f) ++out.n_terrace_pixels;
            max_label = std::max(max_label, lab[i]);
        }
        out.n_terrace_components = max_label;
    }

    std::vector<float> bi(n, -1.0f);
    const std::vector<float>& fpd = fp_r.get_data();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (fpd[i] > 0.5f)     bi[i] = 0.0f;
        else if (tr[i] > 0.5f) bi[i] = 1.0f;
    }

    out.floodplain_mask = std::move(fp_r);
    out.terrace_mask    = lsdtt_raster_cpu(std::move(tr),
                                           std::vector<uint8_t>(am), meta);
    out.band_index      = lsdtt_raster_cpu(std::move(bi),
                                           std::vector<uint8_t>(am), meta);
    return out;
}

} // namespace lsdtt
