// lsdtt_swath_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_swath_cpu.hpp"
#include "lsdtt_proximity_cpu.hpp"      // compute_swath_distances

#include "../utils/lsdtt_statistics.hpp"  // compute_percentile, compute_mean
#include "../utils/lsdtt_thread_utils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lsdtt {

// =============================================================================
// compute_swath_profile
// =============================================================================
swath_profile_result compute_swath_profile(
    const lsdtt_raster_cpu&      reference_raster,
    const lsdtt_raster_cpu&      values_raster,
    const std::vector<lsdtt_xy>& axis_points,
    double                       half_width_metres,
    double                       bin_width_metres,
    int                          n_threads)
{
    // -------------------------------------------------------------------------
    // Validate
    // -------------------------------------------------------------------------
    const lsdtt_raster_meta& meta = reference_raster.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::compute_swath_profile: invalid raster metadata");
    meta.check_projected_metres();

    if (axis_points.size() < 2)
        throw std::runtime_error(
            "lsdtt::compute_swath_profile: at least 2 axis points required");

    if (bin_width_metres <= 0.0)
        throw std::invalid_argument(
            "lsdtt::compute_swath_profile: bin_width_metres must be > 0");

    // -------------------------------------------------------------------------
    // Compute per-pixel swath distances
    // (along_axis, perp_distance, signed_perp over the full valid extent)
    // -------------------------------------------------------------------------
    swath_distances dists = compute_swath_distances(
        reference_raster, axis_points, n_threads);

    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    const std::size_t npix = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);

    const std::vector<float>&   val_data = values_raster.get_data();
    const std::vector<uint8_t>& val_mask = values_raster.get_mask();

    const std::vector<float>& along_data  = dists.along_axis.get_data();
    const std::vector<float>& perp_data   = dists.perp_distance.get_data();
    const std::vector<float>& signed_data = dists.signed_perp.get_data();
    const std::vector<uint8_t>& dist_mask = dists.along_axis.get_mask();

    const float ndv = static_cast<float>(meta.no_data_value);

    // -------------------------------------------------------------------------
    // Apply half-width clipping.
    // Build output distance rasters that carry nodata for pixels outside the
    // half-width constraint.  Also collect (along_axis, value) pairs for
    // binning.
    // -------------------------------------------------------------------------
    const bool do_clip = (half_width_metres > 0.0);
    const float hw_f   = static_cast<float>(half_width_metres);

    std::vector<float>   v_along(npix, ndv);
    std::vector<float>   v_perp(npix, ndv);
    std::vector<float>   v_signed(npix, ndv);
    std::vector<uint8_t> omsk(npix, uint8_t{0});

    // Collect (along_axis_metres, value) for binning.
    // Reserve conservatively.
    std::vector<std::pair<float, float>> samples;
    samples.reserve(npix / 4);

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (dist_mask[i] == 0) continue;  // nodata in distance rasters
        if (val_mask[i]  == 0) continue;  // nodata in values raster

        const float perp = perp_data[i];
        if (do_clip && perp > hw_f) continue;  // outside half-width

        const float v = val_data[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;

        omsk[i]     = uint8_t{1};
        v_along[i]  = along_data[i];
        v_perp[i]   = perp;
        v_signed[i] = signed_data[i];

        samples.emplace_back(along_data[i], v);
    }

    // -------------------------------------------------------------------------
    // Determine total baseline arc-length and bin layout.
    // -------------------------------------------------------------------------
    double total_arc = 0.0;
    for (std::size_t s = 0; s + 1 < axis_points.size(); ++s)
    {
        const double dx = axis_points[s + 1].x - axis_points[s].x;
        const double dy = axis_points[s + 1].y - axis_points[s].y;
        total_arc += std::sqrt(dx * dx + dy * dy);
    }

    const int n_bins = static_cast<int>(std::ceil(total_arc / bin_width_metres));

    // -------------------------------------------------------------------------
    // Bucket samples into bins.
    // -------------------------------------------------------------------------
    std::vector<std::vector<float>> bucket_values(static_cast<std::size_t>(n_bins));
    for (const auto& [along, val] : samples)
    {
        int bin_idx = static_cast<int>(std::floor(static_cast<double>(along) / bin_width_metres));
        // Guard: pixel projected exactly at total_arc lands in bin n_bins-1
        if (bin_idx < 0)         bin_idx = 0;
        if (bin_idx >= n_bins)   bin_idx = n_bins - 1;
        bucket_values[static_cast<std::size_t>(bin_idx)].push_back(val);
    }

    // -------------------------------------------------------------------------
    // Compute per-bin statistics.
    // -------------------------------------------------------------------------
    constexpr std::array<double, 21> q_levels = {
         0.0,  5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
        55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0
    };
    const float nan_f = std::numeric_limits<float>::quiet_NaN();

    std::vector<swath_bin> bins(static_cast<std::size_t>(n_bins));
    for (int b = 0; b < n_bins; ++b)
    {
        swath_bin& sb = bins[static_cast<std::size_t>(b)];
        sb.along_axis_centre = (static_cast<double>(b) + 0.5) * bin_width_metres;

        std::vector<float>& bv = bucket_values[static_cast<std::size_t>(b)];
        sb.n_pixels = static_cast<int>(bv.size());

        if (bv.empty())
        {
            sb.mean = nan_f;
            sb.quantiles.fill(nan_f);
            continue;
        }

        // Sort once; compute_percentile uses sorted data via the overload that
        // accepts a pre-sorted flag, or we sort here and call directly.
        std::sort(bv.begin(), bv.end());

        double sum = 0.0;
        for (float x : bv) sum += static_cast<double>(x);
        sb.mean = static_cast<float>(sum / static_cast<double>(bv.size()));

        for (std::size_t q = 0; q < 21; ++q)
        {
            sb.quantiles[q] = static_cast<float>(
                compute_percentile<float>(bv, q_levels[q])
            );
        }
    }

    // -------------------------------------------------------------------------
    // Assemble result.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> omsk2 = omsk;
    std::vector<uint8_t> omsk3 = omsk;

    swath_profile_result result;
    result.along_axis    = lsdtt_raster_cpu(std::move(v_along),  std::move(omsk),  meta);
    result.perp_distance = lsdtt_raster_cpu(std::move(v_perp),   std::move(omsk2), meta);
    result.signed_perp   = lsdtt_raster_cpu(std::move(v_signed), std::move(omsk3), meta);
    result.bins          = std::move(bins);

    return result;
}

} // namespace lsdtt
