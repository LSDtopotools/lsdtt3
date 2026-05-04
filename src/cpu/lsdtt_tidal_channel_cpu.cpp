// lsdtt_tidal_channel_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_tidal_channel_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>


#ifdef _OPENMP
#include <omp.h>
#endif

#include "lsdtt_proximity_cpu.hpp"
#include "lsdtt_raster_algebra.hpp"
#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_skeleton_graph_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

inline bool is_valid_scalar_(
    const lsdtt_raster_meta&    meta,
    const std::vector<float>&   data,
    const std::vector<uint8_t>& mask,
    std::size_t                 idx)
{
    return mask[idx] != 0 &&
           !meta.is_nodata(data[idx]) &&
           std::isfinite(data[idx]);
}

static void validate_request_(const tidal_channel_request& request)
{
    if (request.median_window_radii_metres.empty())
    {
        throw std::invalid_argument(
            "lsdtt::find_tidal_channels: median_window_radii_metres must not be empty"
        );
    }

    for (double r : request.median_window_radii_metres)
    {
        if (!(r > 0.0))
        {
            throw std::invalid_argument(
                "lsdtt::find_tidal_channels: all median_window_radii_metres must be > 0"
            );
        }
    }

    if (request.residual_thresholds_metres.empty())
    {
        throw std::invalid_argument(
            "lsdtt::find_tidal_channels: residual_thresholds_metres must not be empty"
        );
    }

    if (!(request.residual_thresholds_metres.size() == 1 ||
          request.residual_thresholds_metres.size() == request.median_window_radii_metres.size()))
    {
        throw std::invalid_argument(
            "lsdtt::find_tidal_channels: residual_thresholds_metres must have size 1 or match median_window_radii_metres"
        );
    }

    if (request.apply_binary_closing &&
        !(request.binary_closing_radius_metres > 0.0))
    {
        throw std::invalid_argument(
            "lsdtt::find_tidal_channels: binary_closing_radius_metres must be > 0 when apply_binary_closing is true"
        );
    }

    if (request.min_channel_component_pixels < 1)
    {
        throw std::invalid_argument(
            "lsdtt::find_tidal_channels: min_channel_component_pixels must be >= 1"
        );
    }

    if (request.extract_skeleton)
    {
        if (request.skeleton_spur_max_pixels < 0)
        {
            throw std::invalid_argument(
                "lsdtt::find_tidal_channels: skeleton_spur_max_pixels must be >= 0"
            );
        }

        if (request.skeleton_terminal_loop_max_pixels < 0)
        {
            throw std::invalid_argument(
                "lsdtt::find_tidal_channels: skeleton_terminal_loop_max_pixels must be >= 0"
            );
        }

        if (request.skeleton_min_component_pixels < 1)
        {
            throw std::invalid_argument(
                "lsdtt::find_tidal_channels: skeleton_min_component_pixels must be >= 1"
            );
        }
    }

    if (request.connect_skeleton_segments && !request.extract_skeleton)
    {
        throw std::invalid_argument(
            "lsdtt::find_tidal_channels: connect_skeleton_segments requires extract_skeleton = true"
        );
    }

    if (request.connect_skeleton_segments)
    {
        if (!(request.connect_tangent_fit_length_pixels > 0.0))
        {
            throw std::invalid_argument(
                "lsdtt::find_tidal_channels: connect_tangent_fit_length_pixels must be > 0"
            );
        }

        if (request.connect_tangent_min_support_pixels < 2)
        {
            throw std::invalid_argument(
                "lsdtt::find_tidal_channels: connect_tangent_min_support_pixels must be >= 2"
            );
        }

        if (!(request.connect_max_gap_pixels > 0.0))
        {
            throw std::invalid_argument(
                "lsdtt::find_tidal_channels: connect_max_gap_pixels must be > 0"
            );
        }

        if (!(request.connect_route_corridor_half_width_pixels > 0.0))
        {
            throw std::invalid_argument(
                "lsdtt::find_tidal_channels: connect_route_corridor_half_width_pixels must be > 0"
            );
        }
    }
}

static std::vector<double> resolve_thresholds_(
    const tidal_channel_request& request)
{
    if (request.residual_thresholds_metres.size() == 1)
    {
        return std::vector<double>(
            request.median_window_radii_metres.size(),
            request.residual_thresholds_metres.front()
        );
    }

    return request.residual_thresholds_metres;
}

// -----------------------------------------------------------------------------
// Combine threshold masks across scales.
//
// any_scale_pass:
//   output valid where any scale is valid
//   output foreground where any valid scale is foreground
//
// all_scales_pass:
//   output valid only where all scales are valid
//   output foreground only where all scales are foreground
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu combine_threshold_masks_(
    const std::vector<lsdtt_raster_cpu>& masks,
    tidal_channel_scale_mode             mode,
    int                                  n_threads)
{
    if (masks.empty())
    {
        throw std::runtime_error(
            "lsdtt_tidal_channel_cpu: no threshold masks to combine"
        );
    }

    const lsdtt_raster_meta& meta = masks.front().get_meta();
    const std::size_t npix = masks.front().get_data().size();
    const float ndv = static_cast<float>(meta.no_data_value);

    for (std::size_t i = 1; i < masks.size(); ++i)
    {
        if (!same_grid(meta, masks[i].get_meta()))
        {
            throw std::invalid_argument(
                "lsdtt_tidal_channel_cpu: threshold masks are not co-registered"
            );
        }
    }

    std::vector<float> out_data(npix, ndv);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for if(npix > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(npix); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mode == tidal_channel_scale_mode::any_scale_pass)
        {
            bool any_valid = false;
            bool any_fg = false;

            for (const auto& m : masks)
            {
                const auto& md = m.get_data();
                const auto& mm = m.get_mask();

                if (mm[idx] == 0) continue;

                any_valid = true;
                if (md[idx] > 0.5f)
                {
                    any_fg = true;
                }
            }

            if (any_valid)
            {
                out_mask[idx] = uint8_t{1};
                out_data[idx] = any_fg ? 1.0f : 0.0f;
            }
        }
        else
        {
            bool all_valid = true;
            bool all_fg = true;

            for (const auto& m : masks)
            {
                const auto& md = m.get_data();
                const auto& mm = m.get_mask();

                if (mm[idx] == 0)
                {
                    all_valid = false;
                    break;
                }

                if (!(md[idx] > 0.5f))
                {
                    all_fg = false;
                }
            }

            if (all_valid)
            {
                out_mask[idx] = uint8_t{1};
                out_data[idx] = all_fg ? 1.0f : 0.0f;
            }
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data),
        std::move(out_mask),
        meta
    );
}

// -----------------------------------------------------------------------------
// Build a non-negative support raster for least-cost bridge routing.
//
// We take the per-scale residual rasters, clamp them to non-negative values,
// normalise each to [0,1] by its valid maximum, then keep the union maximum
// across scales.  This gives strongest support where any scale strongly
// indicates a channel low.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu build_connection_support_raster_(
    const std::vector<lsdtt_raster_cpu>& residual_rasters)
{
    if (residual_rasters.empty())
    {
        throw std::runtime_error(
            "lsdtt_tidal_channel_cpu: cannot build connection-support raster from an empty residual list"
        );
    }

    lsdtt_raster_cpu support =
        algebra::raster_normalise_nonnegative(residual_rasters.front());

    for (std::size_t i = 1; i < residual_rasters.size(); ++i)
    {
        const lsdtt_raster_cpu next_support =
            algebra::raster_normalise_nonnegative(residual_rasters[i]);

        support = algebra::raster_max_union(support, next_support);
    }

    return support;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// find_tidal_channels
// -----------------------------------------------------------------------------
tidal_channel_result find_tidal_channels(
    const lsdtt_raster_cpu&      dem,
    const tidal_channel_request& request)
{
    validate_request_(request);

    const lsdtt_raster_meta& meta = dem.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::find_tidal_channels: DEM metadata are invalid"
        );
    }

    meta.check_projected_metres();

    if (dem.get_data().empty())
    {
        throw std::runtime_error(
            "lsdtt::find_tidal_channels: DEM raster is empty"
        );
    }

    tidal_channel_result result;
    result.scales_metres     = request.median_window_radii_metres;
    result.thresholds_metres = resolve_thresholds_(request);

    const std::size_t n_scales = request.median_window_radii_metres.size();
    result.median_rasters.reserve(n_scales);
    result.residual_rasters.reserve(n_scales);
    result.threshold_masks.reserve(n_scales);

    // -------------------------------------------------------------------------
    // Per-scale median residuals and threshold masks
    // -------------------------------------------------------------------------
    for (std::size_t i = 0; i < n_scales; ++i)
    {
        const double radius = request.median_window_radii_metres[i];
        const double threshold = result.thresholds_metres[i];

        std::cout << "  Scale " << (i + 1) << "/" << n_scales
                  << " (radius = " << radius << " m)..." << std::endl;

        lsdtt_raster_cpu median_raster =
            moving_window_median(
                dem,
                radius,
                window_shape::square
            );

        lsdtt_raster_cpu residual_raster =
            algebra::raster_subtract(
                median_raster,
                dem
            );

        lsdtt_raster_cpu threshold_mask =
            algebra::raster_binary_true_above_threshold(
                residual_raster,
                static_cast<float>(threshold)
            );

        result.median_rasters.push_back(std::move(median_raster));
        result.residual_rasters.push_back(std::move(residual_raster));
        result.threshold_masks.push_back(std::move(threshold_mask));
    }

    // -------------------------------------------------------------------------
    // Multi-scale combination
    // -------------------------------------------------------------------------
    std::cout << "  Combining scale masks..." << std::endl;
    result.channel_mask =
        combine_threshold_masks_(
            result.threshold_masks,
            request.scale_mode,
            request.n_threads
        );

    // -------------------------------------------------------------------------
    // Optional binary cleanup
    // -------------------------------------------------------------------------
    if (request.apply_binary_closing || request.fill_holes ||
        request.min_channel_component_pixels > 1)
    {
        std::cout << "  Cleaning channel mask..." << std::endl;
    }

    if (request.apply_binary_closing)
    {
        result.channel_mask =
            binary_morphological_closing(
                result.channel_mask,
                request.binary_closing_radius_metres
            );
    }

    if (request.fill_holes)
    {
        result.channel_mask =
            binary_fill_holes(result.channel_mask);
    }

    if (request.min_channel_component_pixels > 1)
    {
        result.channel_mask =
            filter_binary_components_by_size(
                result.channel_mask,
                request.min_channel_component_pixels,
                pixel_connectivity::d8
            );
    }

    // -------------------------------------------------------------------------
    // Optional UPL
    // -------------------------------------------------------------------------
    if (request.compute_upl)
    {
        const lsdtt_raster_meta& cm_meta = result.channel_mask.get_meta();
        const std::vector<float>& cm_data = result.channel_mask.get_data();
        const std::vector<uint8_t>& cm_mask = result.channel_mask.get_mask();

        std::size_t n_total   = cm_data.size();
        std::size_t n_valid   = 0;
        std::size_t n_invalid = 0;
        std::size_t n_fg      = 0;
        std::size_t n_bg      = 0;

        for (std::size_t i = 0; i < cm_data.size(); ++i)
        {
            if (cm_mask[i] == 0 ||
                cm_meta.is_nodata(cm_data[i]) ||
                !std::isfinite(cm_data[i]))
            {
                ++n_invalid;
                continue;
            }

            ++n_valid;

            if (cm_data[i] > 0.5f)
            {
                ++n_fg;
            }
            else
            {
                ++n_bg;
            }
        }

        std::cout
            << "  UPL diagnostics:\n"
            << "    total pixels              : " << n_total << "\n"
            << "    valid mask pixels         : " << n_valid << "\n"
            << "    invalid/nodata pixels     : " << n_invalid << "\n"
            << "    foreground channel pixels : " << n_fg << "\n"
            << "    background valid pixels   : " << n_bg << "\n";

        if (n_fg == 0)
        {
            throw std::runtime_error(
                "lsdtt::find_tidal_channels: channel mask contains no foreground pixels; cannot compute UPL"
            );
        }

        // ---------------------------------------------------------------------
        // Build an all-valid binary raster for EDT.
        //
        // Some EDT implementations do not cope safely with nodata / masked holes.
        // Here we convert every pixel to valid:
        //   foreground channel -> 1
        //   everything else    -> 0
        // then restore the original valid-domain mask after EDT.
        // ---------------------------------------------------------------------
        std::vector<float> edt_data(n_total, 0.0f);
        std::vector<uint8_t> edt_mask(n_total, uint8_t{1});

        for (std::size_t i = 0; i < n_total; ++i)
        {
            const bool is_valid =
                cm_mask[i] != 0 &&
                !cm_meta.is_nodata(cm_data[i]) &&
                std::isfinite(cm_data[i]);

            edt_data[i] = (is_valid && cm_data[i] > 0.5f) ? 1.0f : 0.0f;
        }

        lsdtt_raster_cpu edt_input(
            std::move(edt_data),
            std::move(edt_mask),
            cm_meta
        );

        std::cout << "  Calling euclidean_distance_transform on sanitised all-valid mask...\n";

        lsdtt_raster_cpu edt_upl =
            euclidean_distance_transform(
                edt_input,
                request.n_threads
            );

        std::cout << "  EDT returned; restoring original valid-domain mask...\n";

        std::vector<float> upl_data = edt_upl.get_data();
        std::vector<uint8_t> upl_mask(n_total, uint8_t{0});

        for (std::size_t i = 0; i < n_total; ++i)
        {
            const bool is_valid =
                cm_mask[i] != 0 &&
                !cm_meta.is_nodata(cm_data[i]) &&
                std::isfinite(cm_data[i]);

            if (is_valid)
            {
                upl_mask[i] = uint8_t{1};
            }
            else
            {
                upl_data[i] = static_cast<float>(cm_meta.no_data_value);
                upl_mask[i] = uint8_t{0};
            }
        }

        result.upl = lsdtt_raster_cpu(
            std::move(upl_data),
            std::move(upl_mask),
            cm_meta
        );

        std::cout << "  UPL complete.\n";
    }



    // -------------------------------------------------------------------------
    // Optional skeletonisation
    // -------------------------------------------------------------------------
    if (request.extract_skeleton)
    {
        std::cout << "  Extracting skeleton with thin medial ridge routine..." << std::endl;
        result.channel_skeleton =
            thin_distance_guided(
                result.channel_mask,
                0.5   // preserve_distance_pixels
            );



        if (request.skeleton_terminal_loop_max_pixels > 0 ||
            request.skeleton_spur_max_pixels > 0 ||
            request.skeleton_min_component_pixels > 1)
        {
            std::cout << "  Pruning skeleton..." << std::endl;
        }

        if (request.skeleton_terminal_loop_max_pixels > 0)
        {
            result.channel_skeleton =
                prune_terminal_loops(
                    result.channel_skeleton,
                    request.skeleton_terminal_loop_max_pixels
                );
        }

        if (request.skeleton_spur_max_pixels > 0)
        {
            result.channel_skeleton =
                prune_skeleton_spurs(
                    result.channel_skeleton,
                    request.skeleton_spur_max_pixels
                );
        }

        if (request.skeleton_min_component_pixels > 1)
        {
            result.channel_skeleton =
                filter_binary_components_by_size(
                    result.channel_skeleton,
                    request.skeleton_min_component_pixels,
                    pixel_connectivity::d8
                );
        }

        // ---------------------------------------------------------------------
        // Optional least-cost reconnection of nearby skeleton segments
        // ---------------------------------------------------------------------
        if (request.connect_skeleton_segments &&
            !result.channel_skeleton.get_data().empty())
        {
            std::cout << "  Reconnecting skeleton segments..." << std::endl;

            const lsdtt_raster_cpu support_raster =
                build_connection_support_raster_(result.residual_rasters);

            const std::vector<skeleton_segment> segments =
                extract_skeleton_segments(result.channel_skeleton);

            const std::vector<skeleton_endpoint_tangent> tangents =
                estimate_skeleton_endpoint_tangents(
                    segments,
                    request.connect_tangent_fit_length_pixels,
                    request.connect_tangent_min_support_pixels
                );

            endpoint_pair_request pair_request;
            pair_request.max_gap_pixels                   = request.connect_max_gap_pixels;
            pair_request.max_facing_angle_degrees         = request.connect_max_facing_angle_degrees;
            pair_request.max_heading_misalignment_degrees = request.connect_max_heading_misalignment_degrees;
            pair_request.require_strong_segments          = false;
            pair_request.reject_same_component            = true;

            const std::vector<skeleton_endpoint_pair_candidate> candidates =
                find_skeleton_endpoint_pair_candidates(
                    segments,
                    tangents,
                    pair_request
                );

            endpoint_pair_selection_request selection_request;
            selection_request.max_score           = request.connect_max_candidate_score;
            selection_request.require_mutual_best = false;
            selection_request.enforce_one_to_one  = true;

            const std::vector<skeleton_endpoint_pair_candidate> accepted =
                select_skeleton_endpoint_pairs(
                    candidates,
                    tangents.size(),
                    selection_request
                );

            skeleton_bridge_route_request route_request;
            route_request.corridor_half_width_pixels = request.connect_route_corridor_half_width_pixels;
            route_request.base_step_cost             = 1.0;
            route_request.support_penalty_weight     = 4.0;
            route_request.foreground_cost_scale      = 0.25;

            const skeleton_bridge_route_result routed =
                route_skeleton_endpoint_bridges_least_cost(
                    result.channel_skeleton,
                    support_raster,
                    accepted,
                    route_request
                );

            result.skeleton_bridge_mask = routed.bridge_mask;
            result.connected_skeleton   = routed.connected_skeleton;
        }
    }

    return result;
}

} // namespace lsdtt
