// lsdtt_swath_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

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
    // Baseline ends.
    // compute_swath_distances measures each pixel to the nearest point on the
    // baseline, end points included, so a pixel beyond either end (in the
    // half-disc behind the start or past the end) has its foot clamped to that
    // end point: along = 0 or along = total arc. Binned, those pixels piled
    // into the first and last bins (a 30 m wide straight section in 2 m bins
    // had ~350 extra pixels in each end bin against 60 in the others). They
    // are not alongside the baseline, so they are left out of the swath.
    // A pixel is beyond an end when its foot is that end point AND it lies
    // ahead of / behind the end along the direction of the end segment.
    // -------------------------------------------------------------------------
    double total_arc = 0.0;
    for (std::size_t s = 0; s + 1 < axis_points.size(); ++s)
    {
        const double dx = axis_points[s + 1].x - axis_points[s].x;
        const double dy = axis_points[s + 1].y - axis_points[s].y;
        total_arc += std::sqrt(dx * dx + dy * dy);
    }

    // Direction of the first and last segments of non-zero length
    // (repeated vertices give zero-length segments).
    double start_dx = 0.0, start_dy = 0.0, end_dx = 0.0, end_dy = 0.0;
    for (std::size_t s = 0; s + 1 < axis_points.size(); ++s)
    {
        start_dx = axis_points[s + 1].x - axis_points[s].x;
        start_dy = axis_points[s + 1].y - axis_points[s].y;
        if (start_dx != 0.0 || start_dy != 0.0) break;
    }
    for (std::size_t s = axis_points.size() - 1; s > 0; --s)
    {
        end_dx = axis_points[s].x - axis_points[s - 1].x;
        end_dy = axis_points[s].y - axis_points[s - 1].y;
        if (end_dx != 0.0 || end_dy != 0.0) break;
    }
    const lsdtt_xy& start_pt = axis_points.front();
    const lsdtt_xy& end_pt   = axis_points.back();
    // along_axis is stored as float: allow for its rounding at the far end.
    const double end_tol = std::max(1e-3, 1e-6 * total_arc);

    auto beyond_an_end = [&](std::size_t i, float along) -> bool
    {
        const double a = static_cast<double>(along);
        const bool at_start = (a <= end_tol);
        const bool at_end   = (a >= total_arc - end_tol);
        if (!at_start && !at_end) return false;
        const std::size_t r = i / static_cast<std::size_t>(nc);
        const std::size_t c = i % static_cast<std::size_t>(nc);
        const double px = meta.x_min + (static_cast<double>(c) + 0.5) * meta.cell_size;
        const double py = meta.y_max() - (static_cast<double>(r) + 0.5) * meta.cell_size;
        if (at_start && ((px - start_pt.x) * start_dx + (py - start_pt.y) * start_dy) < 0.0)
            return true;
        if (at_end && ((px - end_pt.x) * end_dx + (py - end_pt.y) * end_dy) > 0.0)
            return true;
        return false;
    };

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
        if (beyond_an_end(i, along_data[i])) continue;  // not alongside the baseline

        const float v = val_data[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;

        omsk[i]     = uint8_t{1};
        v_along[i]  = along_data[i];
        v_perp[i]   = perp;
        v_signed[i] = signed_data[i];

        samples.emplace_back(along_data[i], v);
    }

    // -------------------------------------------------------------------------
    // Bin layout (total_arc was computed above, with the baseline ends).
    // -------------------------------------------------------------------------
    // The tolerance stops rounding in the baseline length from adding an empty
    // bin: a 960 m line built from coordinates measures 960.0000000001 m, and
    // ceil(480.00000000005) = 481 would add a bin past its end. Pixels at the
    // very end still land in the last bin (bin_idx is clamped below).
    const int n_bins = std::max(1, static_cast<int>(
        std::ceil(total_arc / bin_width_metres - 1e-6)));

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
