// lsdtt_valley_terrace_profile_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// See lsdtt_valley_terrace_profile_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_valley_terrace_profile_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>

#include "../utils/lsdtt_point_bucket_index.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include <memory>

namespace lsdtt
{

namespace
{

// Linear-interpolated quantile of a sorted vector (q in [0, 1]).
float quantile_sorted(const std::vector<float>& v, double q)
{
    if (v.empty()) return std::numeric_limits<float>::quiet_NaN();
    const double pos = q * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    const double t = pos - static_cast<double>(lo);
    return static_cast<float>(v[lo] * (1.0 - t) + v[hi] * t);
}

five_numbers summary(std::vector<float> v)
{
    v.erase(std::remove_if(v.begin(), v.end(),
                           [](float x) { return !std::isfinite(x); }), v.end());
    std::sort(v.begin(), v.end());
    return {quantile_sorted(v, 0.0), quantile_sorted(v, 0.25), quantile_sorted(v, 0.5),
            quantile_sorted(v, 0.75), quantile_sorted(v, 1.0)};
}

} // namespace

terrace_profile_result compute_terrace_profile(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_raster_cpu&                   labels,
    const std::vector<valley_rer_centreline>& centrelines,
    const valley_rer_result&                  rer,
    const terrace_profile_request&            request)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int n_cols = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * n_cols;
    if (labels.get_data().size() != n || rer.owner_centreline.size() != n ||
        rer.owner_node.size() != n)
        throw std::invalid_argument(
            "compute_terrace_profile: labels / relative elevation / DEM size mismatch.");
    if (!(request.bin_metres > 0.0))
        throw std::invalid_argument("compute_terrace_profile: bin_metres must be > 0.");

    const std::vector<float>&        elev = dem.get_data();
    const std::vector<std::uint8_t>& em   = dem.get_mask();
    const std::vector<float>&        ld   = labels.get_data();
    const std::vector<std::uint8_t>& lm   = labels.get_mask();

    // Node positions and along-centreline distances.
    std::vector<std::vector<lsdtt_xy>> nxy(centrelines.size());
    std::vector<std::vector<double>>   along(centrelines.size());
    for (std::size_t c = 0; c < centrelines.size(); ++c)
    {
        const auto& f = centrelines[c].flats;
        nxy[c].resize(f.size());
        along[c].assign(f.size(), 0.0);
        for (std::size_t k = 0; k < f.size(); ++k)
        {
            nxy[c][k] = cell_center_from_row_col(meta, static_cast<int>(f[k] / n_cols),
                                                 static_cast<int>(f[k] % n_cols));
            if (k > 0)
                along[c][k] = along[c][k - 1] + std::hypot(nxy[c][k].x - nxy[c][k - 1].x,
                                                           nxy[c][k].y - nxy[c][k - 1].y);
        }
    }

    // Nearest-node index per centreline, built on first use (only centrelines
    // a patch is reassigned to need one).
    std::vector<std::unique_ptr<point_bucket_index>> node_index(centrelines.size());
    auto index_of = [&](std::size_t c) -> const point_bucket_index& {
        if (!node_index[c])
        {
            std::vector<double> xs(nxy[c].size()), ys(nxy[c].size());
            for (std::size_t k = 0; k < nxy[c].size(); ++k) { xs[k] = nxy[c][k].x; ys[k] = nxy[c][k].y; }
            node_index[c] = std::make_unique<point_bucket_index>(xs, ys, 4.0 * meta.cell_size);
        }
        return *node_index[c]; };

    // Pixels of each patch, and owner counts per patch.
    std::map<std::int64_t, std::vector<std::size_t>> pix;
    std::unordered_map<std::int64_t, std::unordered_map<std::int32_t, std::int64_t>> votes;
    terrace_profile_result res;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!lm[i] || !em[i]) continue;
        const std::int64_t lab = std::llround(ld[i]);
        if (lab <= 0) continue;
        // A pixel with no owner (beyond the distance cap, or unreachable) is
        // left out of the profile, not tied to some node.
        if (rer.owner_centreline[i] < 0) { ++res.n_patch_pixels_unowned; continue; }
        pix[lab].push_back(i);
        ++votes[lab][rer.owner_centreline[i]];
    }

    for (const auto& [lab, cells] : pix)
    {
        const auto vit = votes.find(lab);
        if (vit == votes.end() || vit->second.empty()) continue;
        // The centreline owning most of the patch (ties: lower index).
        std::int32_t cl = -1;
        std::int64_t best = -1;
        for (const auto& [c, cnt] : vit->second)
            if (cnt > best || (cnt == best && c < cl)) { best = cnt; cl = c; }
        const valley_rer_centreline& C = centrelines[static_cast<std::size_t>(cl)];
        if (C.flats.empty()) continue;

        terrace_profile_patch pp;
        pp.patch_label = lab;
        pp.centreline_id = C.id;
        std::map<std::int64_t, std::vector<float>> bin_z, bin_r;
        std::vector<float> all_z, all_r;
        double amin = std::numeric_limits<double>::infinity(), amax = -amin;
        for (const std::size_t i : cells)
        {
            std::int64_t k;
            if (rer.owner_centreline[i] == cl) k = rer.owner_node[i];
            else
            {
                // Owned by another centreline: nearest node of this patch's
                // centreline.
                const lsdtt_xy p = cell_center_from_row_col(
                    meta, static_cast<int>(i / n_cols), static_cast<int>(i % n_cols));
                k = static_cast<std::int64_t>(
                    index_of(static_cast<std::size_t>(cl)).nearest(p.x, p.y));
                ++pp.n_pixels_reassigned;
            }
            const double a = along[static_cast<std::size_t>(cl)][static_cast<std::size_t>(k)];
            const float b = C.baseline[static_cast<std::size_t>(k)];
            const float z = elev[i];
            const float r = std::isfinite(b) ? z - b : std::numeric_limits<float>::quiet_NaN();
            const std::int64_t bin = static_cast<std::int64_t>(std::floor(a / request.bin_metres));
            bin_z[bin].push_back(z);
            bin_r[bin].push_back(r);
            all_z.push_back(z);
            all_r.push_back(r);
            amin = std::min(amin, a);
            amax = std::max(amax, a);
            ++pp.n_pixels;
        }
        pp.along_min_metres = amin;
        pp.along_max_metres = amax;
        pp.elevation_median = summary(all_z)[2];
        pp.relative_elevation_median = summary(all_r)[2];
        res.patches.push_back(pp);
        for (const auto& [bin, zs] : bin_z)
        {
            terrace_profile_bin tb;
            tb.centreline_id = C.id;
            tb.patch_label = lab;
            tb.bin_centre_metres = (static_cast<double>(bin) + 0.5) * request.bin_metres;
            tb.n_pixels = static_cast<std::int64_t>(zs.size());
            tb.elevation = summary(zs);
            tb.relative_elevation = summary(bin_r[bin]);
            res.bins.push_back(tb);
        }
    }
    std::stable_sort(res.bins.begin(), res.bins.end(),
                     [](const terrace_profile_bin& a, const terrace_profile_bin& b) {
                         return a.centreline_id < b.centreline_id; });
    return res;
}

} // namespace lsdtt
