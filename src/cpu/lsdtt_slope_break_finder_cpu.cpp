// lsdtt_slope_break_finder_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_slope_break_finder_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#ifdef LSDTT_HAS_OPENCV
#include "lsdtt_raster_opencv_cpu.hpp"
#endif

#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_skeleton_graph_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_algebra.hpp"
#include "../core/lsdtt_data_loader.hpp"



namespace lsdtt {

// -----------------------------------------------------------------------------
// keep_adjacent_to
//
// Labels D8-connected components of `candidate` foreground pixels.
// Keeps an entire component if:
//   (a) ANY pixel in it is D8-adjacent to a foreground pixel in `prev`, OR
//   (b) the component has >= min_isolated_pixels pixels.
// Criterion (b) retains large isolated curvilinear features that are separated
// from the seeded network by a small gap.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu keep_adjacent_to(
    const lsdtt_raster_cpu& candidate,
    const lsdtt_raster_cpu& prev,
    int                     min_isolated_pixels)
{
    const auto& cd = candidate.get_data();
    const auto& cm = candidate.get_mask();
    const auto& pd = prev.get_data();
    const auto& pm = prev.get_mask();
    const lsdtt_raster_meta& meta = candidate.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const int npix  = nrows * ncols;
    const float ndv = static_cast<float>(meta.no_data_value);

    // --- Step 1: label D8-connected components of candidate foreground pixels ---
    // label 0 = background / unlabelled; labels start at 1
    std::vector<int> labels(npix, 0);
    int next_label = 1;

    // BFS flood fill
    std::vector<int> queue;
    queue.reserve(1024);

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const int idx = r * ncols + c;
            if (cm[idx] == 0 || cd[idx] < 0.5f || labels[idx] != 0) continue;

            const int lbl = next_label++;
            labels[idx] = lbl;
            queue.clear();
            queue.push_back(idx);

            for (std::size_t qi = 0; qi < queue.size(); ++qi)
            {
                const int cur  = queue[qi];
                const int cr   = cur / ncols;
                const int cc   = cur % ncols;
                for (int dr = -1; dr <= 1; ++dr)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;
                        const int nr = cr + dr;
                        const int nc = cc + dc;
                        if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                        const int ni = nr * ncols + nc;
                        if (cm[ni] == 0 || cd[ni] < 0.5f || labels[ni] != 0) continue;
                        labels[ni] = lbl;
                        queue.push_back(ni);
                    }
                }
            }
        }
    }

    // --- Step 2: count component sizes and mark which labels touch prev ---
    const std::size_t nlabels = static_cast<std::size_t>(next_label);
    std::vector<bool> touched(nlabels, false);
    std::vector<int>  comp_size(nlabels, 0);

    for (int i = 0; i < npix; ++i)
    {
        const int lbl = labels[i];
        if (lbl != 0) ++comp_size[static_cast<std::size_t>(lbl)];
    }

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const int idx = r * ncols + c;
            if (labels[idx] == 0) continue;

            for (int dr = -1; dr <= 1; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    const int nr = r + dr;
                    const int nc = c + dc;
                    if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                    const int ni = nr * ncols + nc;
                    if (pm[ni] != 0 && pd[ni] > 0.5f)
                        touched[static_cast<std::size_t>(labels[idx])] = true;
                }
            }
        }
    }

    // --- Step 3: build output — keep pixels whose component passes either criterion ---
    std::vector<float>   out_data(npix, ndv);
    std::vector<uint8_t> out_mask(npix, 0);

    for (int i = 0; i < npix; ++i)
    {
        if (cm[i] == 0) continue;
        out_mask[i] = 1;
        const int lbl = labels[i];
        const bool keep = lbl != 0 &&
            (touched[static_cast<std::size_t>(lbl)] ||
             comp_size[static_cast<std::size_t>(lbl)] >= min_isolated_pixels);
        out_data[i] = keep ? 1.0f : 0.0f;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// binary_or
//
// Returns a binary raster: 1 where either a or b is 1 (union semantics).
// Nodata follows intersection of masks.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu binary_or(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b)
{
    const auto& ad = a.get_data();
    const auto& am = a.get_mask();
    const auto& bd = b.get_data();
    const auto& bm = b.get_mask();
    const lsdtt_raster_meta& meta = a.get_meta();
    const float ndv = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(ad.size(), ndv);
    std::vector<uint8_t> out_mask(am.size(), 0);

    for (std::size_t i = 0; i < ad.size(); ++i)
    {
        if (am[i] == 0 || bm[i] == 0) continue;
        out_mask[i] = 1;
        out_data[i] = (ad[i] > 0.5f || bd[i] > 0.5f) ? 1.0f : 0.0f;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// mask_to_flat
//
// Restricts a curvature raster to pixels where the flat mask is foreground.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu mask_to_flat(
    const lsdtt_raster_cpu& curv,
    const lsdtt_raster_cpu& flat_mask)
{
    const auto& cd = curv.get_data();
    const auto& cm = curv.get_mask();
    const auto& fd = flat_mask.get_data();
    const auto& fm = flat_mask.get_mask();
    const lsdtt_raster_meta& meta = curv.get_meta();
    const float ndv = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(cd.size(), ndv);
    std::vector<uint8_t> out_mask(cm.size(), 0);

    for (std::size_t i = 0; i < cd.size(); ++i)
    {
        if (cm[i] == 0) continue;
        if (fm[i] == 0 || fd[i] < 0.5f) continue;
        out_data[i] = cd[i];
        out_mask[i] = 1;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// count_foreground
// -----------------------------------------------------------------------------
static int count_foreground(const lsdtt_raster_cpu& r)
{
    const auto& d = r.get_data();
    const auto& m = r.get_mask();
    int n = 0;
    for (std::size_t i = 0; i < d.size(); ++i)
        if (m[i] != 0 && d[i] > 0.5f) ++n;
    return n;
}

// =============================================================================
// find_slope_breaks
// =============================================================================
slope_break_result find_slope_breaks(
    const lsdtt_raster_cpu&    dem,
    const lsdtt_raster_cpu&    mask,
    const slope_break_request& request)
{
    if (request.canny_thresholds.empty())
        throw std::invalid_argument("find_slope_breaks: canny_thresholds must not be empty");

    for (std::size_t i = 1; i < request.canny_thresholds.size(); ++i)
    {
        if (request.canny_thresholds[i] >= request.canny_thresholds[i - 1])
            throw std::invalid_argument(
                "find_slope_breaks: canny_thresholds must be strictly descending");
    }

#ifndef LSDTT_HAS_OPENCV
    throw std::runtime_error(
        "find_slope_breaks: OpenCV is required for Canny edge detection but was "
        "not compiled in (LSDTT_HAS_OPENCV not defined)");
#else

    slope_break_result result;

    // ------------------------------------------------------------------
    // 1. Polyfit: compute vertical_excess_curvature and profile_curvature
    // ------------------------------------------------------------------
    std::cout << "\n--- Slope break finder: polyfit (lambda="
              << request.polyfit_lambda_metres << " m) ---\n";

    poly_fit_coeffs coeffs = direct_polyfit(
        dem,
        request.polyfit_lambda_metres,
        request.n_threads);

    polyfit_metric_request pmr;
    pmr.metrics   = { polyfit_metric::vertical_excess_curvature,
                      polyfit_metric::profile_curvature };
    pmr.n_threads = request.n_threads;

    polyfit_metric_result pmres = compute_polyfit_metrics(coeffs, pmr);

    result.vertical_excess_curvature =
        mask_to_flat(pmres.at(polyfit_metric::vertical_excess_curvature), mask);
    result.profile_curvature =
        mask_to_flat(pmres.at(polyfit_metric::profile_curvature), mask);

    std::cout << "  vertical_excess_curvature: "
              << count_foreground(result.vertical_excess_curvature) << " valid flat pixels\n";
    std::cout << "  profile_curvature:         "
              << count_foreground(result.profile_curvature) << " valid flat pixels\n";

    // ------------------------------------------------------------------
    // 2. Progressive Canny annealing
    // ------------------------------------------------------------------
    std::cout << "\n--- Progressive Canny annealing ("
              << request.canny_thresholds.size() << " steps) ---\n";

    // Accumulated growing masks (empty = no previous seed yet)
    lsdtt_raster_cpu acc_vert;
    lsdtt_raster_cpu acc_prof;
    bool first_step = true;

    for (std::size_t s = 0; s < request.canny_thresholds.size(); ++s)
    {
        const double thr = request.canny_thresholds[s];
        std::cout << "\n  Step " << (s + 1) << "/" << request.canny_thresholds.size()
                  << "  Canny threshold=" << thr << "\n";

        slope_break_step step;
        step.canny_threshold = thr;

        // Canny on vertical_excess_curvature
        step.canny_vert = opencv_cpu::detect_edges_canny(
            result.vertical_excess_curvature,
            thr,
            request.canny_ratio,
            request.canny_aperture_size,
            request.canny_gaussian_sigma,
            request.n_threads);

        // Canny on profile_curvature
        step.canny_prof = opencv_cpu::detect_edges_canny(
            result.profile_curvature,
            thr,
            request.canny_ratio,
            request.canny_aperture_size,
            request.canny_gaussian_sigma,
            request.n_threads);

        std::cout << "    canny_vert edges: " << count_foreground(step.canny_vert) << "\n";
        std::cout << "    canny_prof edges: " << count_foreground(step.canny_prof) << "\n";

        if (first_step)
        {
            // Seed: all Canny edges at the highest threshold
            acc_vert  = step.canny_vert;
            acc_prof  = step.canny_prof;
            first_step = false;
        }
        else
        {
            // Grow: keep components adjacent to previous mask or large enough
            const lsdtt_raster_cpu grown_vert =
                keep_adjacent_to(step.canny_vert, acc_vert, request.min_isolated_pixels);
            const lsdtt_raster_cpu grown_prof =
                keep_adjacent_to(step.canny_prof, acc_prof, request.min_isolated_pixels);

            // Union of previous mask and newly grown edges
            acc_vert = binary_or(acc_vert, grown_vert);
            acc_prof = binary_or(acc_prof, grown_prof);
        }

        step.break_mask_vert = acc_vert;
        step.break_mask_prof = acc_prof;

        std::cout << "    break_mask_vert pixels: " << count_foreground(acc_vert) << "\n";
        std::cout << "    break_mask_prof pixels: " << count_foreground(acc_prof) << "\n";

        result.steps.push_back(std::move(step));
    }

    // ------------------------------------------------------------------
    // 3. Combine: intersection of the two final growing masks
    // ------------------------------------------------------------------
    {
        const auto& vd = acc_vert.get_data();
        const auto& vm = acc_vert.get_mask();
        const auto& pd = acc_prof.get_data();
        const auto& pm = acc_prof.get_mask();
        const lsdtt_raster_meta& meta = acc_vert.get_meta();
        const float ndv = static_cast<float>(meta.no_data_value);
        const int npix = static_cast<int>(vd.size());

        std::vector<float>   out_data(npix, ndv);
        std::vector<uint8_t> out_mask(npix, 0);

        for (int i = 0; i < npix; ++i)
        {
            if (vm[i] == 0 || pm[i] == 0) continue;
            out_mask[i] = 1;
            out_data[i] = (vd[i] > 0.5f && pd[i] > 0.5f) ? 1.0f : 0.0f;
        }

        result.combined_break_mask =
            lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);

        std::cout << "\n--- Combined break mask (intersection) ---\n";
        std::cout << "    combined_break_mask pixels: "
                  << count_foreground(result.combined_break_mask) << "\n";
    }

    // ------------------------------------------------------------------
    // 4. Close the combined mask to seal narrow gaps
    // ------------------------------------------------------------------
    if (request.closing_radius_metres > 0.0)
    {
        std::cout << "\n--- Closing combined mask (radius="
                  << request.closing_radius_metres << " m) ---\n";
        result.combined_closed =
            binary_closing_by_reconstruction(
                result.combined_break_mask,
                request.closing_radius_metres);
        std::cout << "    combined_closed pixels: "
                  << count_foreground(result.combined_closed) << "\n";
    }

    // ------------------------------------------------------------------
    // 5. Build a belt mask from the rail mask
    //
    // Seed from combined_closed when present, otherwise from combined_break_mask.
    // Then:
    //   - optional isotropic closing to merge paired rails into a band
    //   - optional isotropic dilation to thicken slightly
    //   - optional component filtering
    // ------------------------------------------------------------------
    {
        std::cout << "\n--- Belt construction from break-mask rails ---\n";

        lsdtt_raster_cpu belt_seed =
            (!result.combined_closed.get_data().empty())
            ? result.combined_closed
            : result.combined_break_mask;

        lsdtt_raster_cpu belt = belt_seed;

        if (request.belt_closing_radius_metres > 0.0)
        {
            std::cout << "    belt closing radius: "
                      << request.belt_closing_radius_metres << " m\n";
            belt = binary_morphological_closing(
                belt,
                request.belt_closing_radius_metres);
        }

        if (request.belt_dilation_radius_metres > 0.0)
        {
            std::cout << "    belt dilation radius: "
                      << request.belt_dilation_radius_metres << " m\n";
            belt = binary_morphological_dilation(
                belt,
                request.belt_dilation_radius_metres);
        }

        if (request.belt_min_component_pixels > 1)
        {
            std::cout << "    belt min component size: "
                      << request.belt_min_component_pixels << " px\n";
            belt = filter_binary_components_by_size(
                belt,
                request.belt_min_component_pixels);
        }

        result.break_belt = std::move(belt);

        std::cout << "    break_belt pixels: "
                  << count_foreground(result.break_belt) << "\n";
    }

    // ------------------------------------------------------------------
    // 6. Extract a 1-pixel centreline from the belt
    // ------------------------------------------------------------------
    if (request.extract_belt_centreline)
    {
        std::cout << "\n--- Centreline extraction from break belt ---\n";

        lsdtt_raster_cpu centreline = thin_zhang_suen(result.break_belt);
        std::cout << "    raw centreline pixels: "
                  << count_foreground(centreline) << "\n";

        if (request.centreline_max_spur_pixels > 0)
        {
            std::cout << "    pruning spurs <= "
                      << request.centreline_max_spur_pixels
                      << " px\n";
            centreline = prune_skeleton_spurs(
                centreline,
                request.centreline_max_spur_pixels);
        }

        if (request.centreline_max_terminal_loop_pixels > 0)
        {
            std::cout << "    pruning spurs <= "
                      << request.centreline_max_terminal_loop_pixels
                      << " px\n";
            centreline = prune_terminal_loops(
                centreline,
                request.centreline_max_terminal_loop_pixels);
        }


        if (request.centreline_min_component_pixels > 1)
        {
            std::cout << "    centreline min component size: "
                      << request.centreline_min_component_pixels
                      << " px\n";
            centreline = filter_binary_components_by_size(
                centreline,
                request.centreline_min_component_pixels);
        }

        result.break_centreline = std::move(centreline);

        std::cout << "    final centreline pixels: "
                  << count_foreground(result.break_centreline) << "\n";
    }

    // ------------------------------------------------------------------
    // 7. Strength-gated connection of centreline segments
    // ------------------------------------------------------------------
    if (request.connect_centreline_segments &&
        !result.break_centreline.get_data().empty())
    {
        std::cout << "\n--- Strength-gated centreline reconnection ---\n";

        // Build a non-negative connection-support raster from the two
        // curvature rasters.  Use the minimum of the two normalised magnitudes
        // so that support is strong only where both metrics are strong.
        lsdtt_raster_cpu vert_strength =
            algebra::raster_normalise_nonnegative(
                algebra::raster_abs(result.vertical_excess_curvature));

        lsdtt_raster_cpu prof_strength =
            algebra::raster_normalise_nonnegative(
                algebra::raster_abs(result.profile_curvature));

        result.break_connection_strength =
            algebra::raster_min(vert_strength, prof_strength);

        std::vector<skeleton_segment> segments =
            extract_skeleton_segments(result.break_centreline);

        std::vector<skeleton_segment_strength> segment_strengths =
            compute_skeleton_segment_strength(
                segments,
                result.break_connection_strength);

        std::vector<skeleton_endpoint_tangent> tangents =
            estimate_skeleton_endpoint_tangents(
                segments,
                request.connect_tangent_fit_length_pixels,
                request.connect_tangent_min_support_pixels);

        std::cout << "    segments: " << segments.size() << "\n";
        std::cout << "    tangents: " << tangents.size() << "\n";

        endpoint_pair_request epr;
        epr.max_gap_pixels                   = request.connect_max_gap_pixels;
        epr.max_facing_angle_degrees         = request.connect_max_facing_angle_degrees;
        epr.max_heading_misalignment_degrees = request.connect_max_heading_misalignment_degrees;
        epr.require_strong_segments          = true;
        epr.min_segment_mean_strength        = request.connect_min_segment_mean_strength;
        epr.reject_same_component            = true;

        std::vector<skeleton_endpoint_pair_candidate> candidates =
            find_skeleton_endpoint_pair_candidates(
                segments,
                tangents,
                epr,
                segment_strengths);

        std::cout << "    strength-gated candidates: "
                  << candidates.size() << "\n";

        endpoint_pair_selection_request epsr;
        epsr.max_score           = request.connect_max_candidate_score;
        epsr.require_mutual_best = false;
        epsr.enforce_one_to_one  = true;

        std::vector<skeleton_endpoint_pair_candidate> accepted =
            select_skeleton_endpoint_pairs(
                candidates,
                tangents.size(),
                epsr);

        std::cout << "    accepted pairs: "
                  << accepted.size() << "\n";

        skeleton_bridge_route_request brr;
        brr.corridor_half_width_pixels = request.connect_route_corridor_half_width_pixels;
        brr.base_step_cost             = 1.0;
        brr.support_penalty_weight     = 4.0;
        brr.foreground_cost_scale      = 0.25;

        skeleton_bridge_route_result routed =
            route_skeleton_endpoint_bridges_least_cost(
                result.break_centreline,
                result.break_connection_strength,
                accepted,
                brr);

        result.break_bridge_mask = routed.bridge_mask;
        result.break_centreline_connected = routed.connected_skeleton;

        std::cout << "    routed bridges: "
                  << routed.n_routed << "\n";
    }


    return result;

#endif // LSDTT_HAS_OPENCV
}


// =============================================================================
// find_vec_skeleton
// =============================================================================
vec_skeleton_result find_vec_skeleton(
    const lsdtt_raster_cpu&     dem,
    const lsdtt_raster_cpu&     mask,
    const vec_skeleton_request& request)
{
#ifndef LSDTT_HAS_OPENCV
    throw std::runtime_error(
        "find_vec_skeleton: OpenCV is required for adaptive thresholding "
        "but was not compiled in (LSDTT_HAS_OPENCV not defined)");
#else
    vec_skeleton_result result;

    // ------------------------------------------------------------------
    // 1. Direct polyfit → VEC masked to flat pixels
    // ------------------------------------------------------------------
    std::cout << "\n--- VEC skeleton: polyfit (lambda="
              << request.polyfit_lambda_metres << " m) ---\n";

    poly_fit_coeffs coeffs = direct_polyfit(
        dem, request.polyfit_lambda_metres, request.n_threads);

    polyfit_metric_request pmr;
    pmr.metrics   = { polyfit_metric::vertical_excess_curvature };
    pmr.n_threads = request.n_threads;

    polyfit_metric_result pmres = compute_polyfit_metrics(coeffs, pmr);

    result.vertical_excess_curvature = mask_to_flat(
        pmres.at(polyfit_metric::vertical_excess_curvature), mask);

    std::cout << "  VEC valid flat pixels: "
              << count_foreground(result.vertical_excess_curvature) << "\n";

    // ------------------------------------------------------------------
    // 2. Adaptive threshold → binary mask
    // ------------------------------------------------------------------
    std::cout << "\n--- Adaptive thresholding ---\n";
    std::cout << "  method     : " << request.adaptive_method << "\n";
    std::cout << "  block size : " << request.adaptive_block_size << " px\n";
    std::cout << "  C          : " << request.adaptive_C << "\n";
    std::cout << "  pre-blur σ : " << request.adaptive_gaussian_sigma << "\n";

    if (request.adaptive_method == "mean")
    {
        result.adaptive_binary = opencv_cpu::threshold_adaptive_mean(
            result.vertical_excess_curvature,
            request.adaptive_block_size,
            request.adaptive_C,
            request.adaptive_gaussian_sigma,
            request.n_threads);
    }
    else
    {
        result.adaptive_binary = opencv_cpu::threshold_adaptive_gaussian(
            result.vertical_excess_curvature,
            request.adaptive_block_size,
            request.adaptive_C,
            request.adaptive_gaussian_sigma,
            request.n_threads);
    }

    std::cout << "  adaptive_binary foreground pixels: "
              << count_foreground(result.adaptive_binary) << "\n";

    // ------------------------------------------------------------------
    // 3. Optional closing to seal small gaps before skeletonisation
    // ------------------------------------------------------------------
    if (request.closing_radius_metres > 0.0)
    {
        std::cout << "  Closing binary mask (radius="
                  << request.closing_radius_metres << " m)...\n";
        result.adaptive_binary = binary_morphological_closing(
            result.adaptive_binary, request.closing_radius_metres);
        std::cout << "  after closing: "
                  << count_foreground(result.adaptive_binary) << "\n";
    }

    // ------------------------------------------------------------------
    // 3b. Connected-component filter on the binary mask
    // ------------------------------------------------------------------
    if (request.binary_min_component_pixels > 1)
    {
        std::cout << "  Filtering binary components (min="
                  << request.binary_min_component_pixels << " px)...\n";
        result.adaptive_binary = filter_binary_components_by_size(
            result.adaptive_binary, request.binary_min_component_pixels);
        std::cout << "  after binary component filter: "
                  << count_foreground(result.adaptive_binary) << "\n";
    }

    // ------------------------------------------------------------------
    // 4. Zhang–Suen thinning → spur pruning → component filtering
    // ------------------------------------------------------------------
    std::cout << "\n--- Skeletonisation ---\n";

    result.skeleton = thin_zhang_suen(result.adaptive_binary);
    std::cout << "  raw skeleton pixels: "
              << count_foreground(result.skeleton) << "\n";

    // Fill tiny false holes / loop interiors in the 1-pixel centreline, then
    // re-skeletonise.
    result.skeleton = fill_dense_false_pixels_3x3(
        result.skeleton,
        5,   // min_true_neighbours
        1    // n_iterations
    );
    std::cout << "  after 3x3 dense-false fill: "
              << count_foreground(result.skeleton) << "\n";

    result.skeleton = thin_zhang_suen(result.skeleton);
    std::cout << "  after re-thinning: "
              << count_foreground(result.skeleton) << "\n";

    if (request.terminal_loop_max_pixels > 0)
    {
        result.skeleton = prune_terminal_loops(
            result.skeleton, request.terminal_loop_max_pixels);
        std::cout << "  after terminal loop pruning: "
                  << count_foreground(result.skeleton) << "\n";
    }

    if (request.spur_max_pixels > 0)
    {
        result.skeleton = prune_skeleton_spurs(
            result.skeleton, request.spur_max_pixels);
        std::cout << "  after spur pruning: "
                  << count_foreground(result.skeleton) << "\n";
    }



    if (request.min_component_pixels > 1)
    {
        result.skeleton = filter_binary_components_by_size(
            result.skeleton, request.min_component_pixels);
        std::cout << "  after component filter: "
                  << count_foreground(result.skeleton) << "\n";
    }

    // ------------------------------------------------------------------
    // 5. Elevation-based contour connection of skeleton endpoints
    // ------------------------------------------------------------------
    if (request.elevation_connect.search_radius_metres > 0.0 &&
        !result.skeleton.get_data().empty())
    {
        elevation_connect_result ecres =
            connect_skeleton_by_elevation(result.skeleton, dem,
                                          request.elevation_connect);
        result.break_bridge_mask  = std::move(ecres.bridge_mask);
        result.skeleton_connected = std::move(ecres.skeleton_connected);
    }

    return result;

#endif // LSDTT_HAS_OPENCV
}


// =============================================================================
// extract_curvature_skeletons
// =============================================================================
curvature_skeleton_result extract_curvature_skeletons(
    const lsdtt_raster_cpu&           vec_raster,
    const lsdtt_raster_cpu&           profile_raster,
    const curvature_skeleton_request& request,
    lsdtt_data_loader*                loader)
{
    if (vec_raster.get_data().empty())
        throw std::runtime_error(
            "extract_curvature_skeletons: vertical_excess_curvature raster is empty");
    if (profile_raster.get_data().empty())
        throw std::runtime_error(
            "extract_curvature_skeletons: profile_curvature raster is empty");

    curvature_skeleton_result result;
    const lsdtt_raster_meta& meta = vec_raster.get_meta();
    const float ndv  = static_cast<float>(meta.no_data_value);
    const int npix   = meta.n_rows * meta.n_cols;

    // ------------------------------------------------------------------
    // 1. VEC skeleton
    // ------------------------------------------------------------------
    std::cout << "\n--- VEC skeleton ---\n";

    // Always print the CDF divergence threshold for threshold guidance
    double vec_thr = request.vec_threshold;
    try
    {
        const double cdf_thr = algebra::raster_cdf_divergence_threshold(
            vec_raster,
            /*verbose=*/true,
            request.vec_write_cdf && loader != nullptr,
            loader,
            request.vec_cdf_path);

        if (request.vec_use_cdf_threshold)
        {
            std::cout << "  Using CDF divergence threshold: " << cdf_thr << "\n";
            vec_thr = cdf_thr;
        }
        else
        {
            std::cout << "  CDF divergence threshold (reference): " << cdf_thr
                      << "  user threshold: " << vec_thr << "\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "  Warning: CDF threshold computation failed: "
                  << e.what() << " — using vec_threshold=" << vec_thr << "\n";
    }

    // Threshold: 1 where VEC >= vec_thr
    {
        const auto& vd = vec_raster.get_data();
        const auto& vm = vec_raster.get_mask();
        std::vector<float>   bd(npix, ndv);
        std::vector<uint8_t> bm(npix, 0);

        for (int i = 0; i < npix; ++i)
        {
            if (vm[i] == 0) continue;
            bm[i] = 1;
            bd[i] = (vd[i] >= static_cast<float>(vec_thr)) ? 1.0f : 0.0f;
        }

        result.vec_binary = lsdtt_raster_cpu(std::move(bd), std::move(bm), meta);
        std::cout << "  vec_binary foreground pixels: "
                  << count_foreground(result.vec_binary) << "\n";
    }

    // Skeletonise
    result.vec_skeleton = thin_zhang_suen(result.vec_binary);
    std::cout << "  vec_skeleton raw pixels: "
              << count_foreground(result.vec_skeleton) << "\n";

    if (request.spur_max_pixels > 0)
    {
        result.vec_skeleton = prune_skeleton_spurs(
            result.vec_skeleton, request.spur_max_pixels);
        std::cout << "  vec_skeleton after spur pruning: "
                  << count_foreground(result.vec_skeleton) << "\n";
    }

    if (request.terminal_loop_max_pixels > 0)
    {
        result.vec_skeleton = prune_terminal_loops(
            result.vec_skeleton, request.terminal_loop_max_pixels);
        std::cout << "  after terminal loop pruning: "
                  << count_foreground(result.vec_skeleton) << "\n";
    }

    if (request.min_component_pixels > 1)
    {
        result.vec_skeleton = filter_binary_components_by_size(
            result.vec_skeleton, request.min_component_pixels);
        std::cout << "  vec_skeleton after component filter: "
                  << count_foreground(result.vec_skeleton) << "\n";
    }

    // ------------------------------------------------------------------
    // 2. Profile curvature relief skeleton
    // ------------------------------------------------------------------
    std::cout << "\n--- Profile curvature relief skeleton ---\n";
    std::cout << "  Relief window: " << request.profile_relief_window_metres << " m\n";
    std::cout << "  Relief threshold: " << request.profile_relief_threshold << "\n";

    // Local relief of profile_curvature captures sign-change magnitude at
    // slope breaks (positive and negative curvature adjacent to each other
    // means a large within-window range).
    result.profile_relief = moving_window_local_relief(
        profile_raster,
        request.profile_relief_window_metres,
        window_shape::circular);

    std::cout << "  profile_relief computed: "
              << count_foreground(result.profile_relief) << " valid pixels\n";

    // Threshold
    {
        const auto& rd = result.profile_relief.get_data();
        const auto& rm = result.profile_relief.get_mask();
        const int rn   = static_cast<int>(rd.size());
        std::vector<float>   bd(rn, ndv);
        std::vector<uint8_t> bm(rn, 0);

        const float thr_f = static_cast<float>(request.profile_relief_threshold);
        for (int i = 0; i < rn; ++i)
        {
            if (rm[i] == 0) continue;
            bm[i] = 1;
            bd[i] = (rd[i] >= thr_f) ? 1.0f : 0.0f;
        }

        result.profile_binary = lsdtt_raster_cpu(
            std::move(bd), std::move(bm), result.profile_relief.get_meta());
        std::cout << "  profile_binary foreground pixels: "
                  << count_foreground(result.profile_binary) << "\n";
    }

    // Skeletonise
    result.profile_skeleton = thin_zhang_suen(result.profile_binary);
    std::cout << "  profile_skeleton raw pixels: "
              << count_foreground(result.profile_skeleton) << "\n";

    if (request.spur_max_pixels > 0)
    {
        result.profile_skeleton = prune_skeleton_spurs(
            result.profile_skeleton, request.spur_max_pixels);
        std::cout << "  profile_skeleton after spur pruning: "
                  << count_foreground(result.profile_skeleton) << "\n";
    }

    if (request.terminal_loop_max_pixels > 0)
    {
        result.profile_skeleton = prune_terminal_loops(
            result.profile_skeleton, request.terminal_loop_max_pixels);
        std::cout << "  after terminal loop pruning: "
                  << count_foreground(result.profile_skeleton) << "\n";
    }


    if (request.min_component_pixels > 1)
    {
        result.profile_skeleton = filter_binary_components_by_size(
            result.profile_skeleton, request.min_component_pixels);
        std::cout << "  profile_skeleton after component filter: "
                  << count_foreground(result.profile_skeleton) << "\n";
    }

    return result;
}


// =============================================================================
// connect_skeleton_by_elevation
// =============================================================================

namespace {

// BFS connected-component labelling on a binary raster.
// Returns label image (0 = background, 1..n = components).
std::vector<int> label_components(const std::vector<float>& data,
                                   int nrows, int ncols)
{
    std::vector<int> labels(nrows * ncols, 0);
    int next_label = 1;
    for (int r = 0; r < nrows; ++r)
    for (int c = 0; c < ncols; ++c)
    {
        int idx = r * ncols + c;
        if (data[idx] <= 0.5f || labels[idx] != 0) continue;
        std::queue<std::pair<int,int>> q;
        q.push({r, c});
        labels[idx] = next_label;
        while (!q.empty())
        {
            auto [cr, cc] = q.front(); q.pop();
            for (int dr = -1; dr <= 1; ++dr)
            for (int dc = -1; dc <= 1; ++dc)
            {
                if (!dr && !dc) continue;
                int nr = cr + dr, nc = cc + dc;
                if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                int ni = nr * ncols + nc;
                if (data[ni] <= 0.5f || labels[ni] != 0) continue;
                labels[ni] = next_label;
                q.push({nr, nc});
            }
        }
        ++next_label;
    }
    return labels;
}

} // anonymous namespace

elevation_connect_result connect_skeleton_by_elevation(
    const lsdtt_raster_cpu&          skeleton,
    const lsdtt_raster_cpu&          dem,
    const elevation_connect_request& request)
{
    const lsdtt_raster_meta& meta = skeleton.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    if (dem.get_meta().n_rows != nrows || dem.get_meta().n_cols != ncols)
        throw std::invalid_argument(
            "connect_skeleton_by_elevation: skeleton and DEM dimensions must match");

    const double cell = meta.cell_size;
    const int search_px = static_cast<int>(std::ceil(
        request.search_radius_metres / cell));
    const float dem_nd = static_cast<float>(dem.get_meta().no_data_value);

    const std::vector<float>& skel_data = skeleton.get_data();
    const std::vector<float>& dem_data  = dem.get_data();

    // ------------------------------------------------------------------
    // 1. Label connected components and compute per-component pixel counts
    // ------------------------------------------------------------------
    std::vector<int> comp_labels;
    std::unordered_map<int,int> comp_size;

    if (!request.connect_same_component)
    {
        comp_labels = label_components(skel_data, nrows, ncols);
        for (int lbl : comp_labels)
            if (lbl > 0) ++comp_size[lbl];
    }

    // ------------------------------------------------------------------
    // 2. Extract endpoints: skeleton pixels with exactly one 8-neighbour
    // ------------------------------------------------------------------
    // 2. Extract endpoints using crossing-number classification.
    //    classify_skeleton_nodes uses CN=1 (not raw degree=1) so corner
    //    endpoints — e.g. south + southeast neighbours — are not missed.
    // ------------------------------------------------------------------
    struct Endpoint { int r, c; double elev; int label; int seg_size; };
    std::vector<Endpoint> eps;

    {
        const skeleton_node_result snr = classify_skeleton_nodes(skeleton);
        const std::vector<skeleton_node>& ep_nodes = snr.endpoints;
        for (const skeleton_node& nd : ep_nodes)
        {
            float elev = dem_data[nd.row * ncols + nd.col];
            if (elev == dem_nd) continue;
            int lbl = request.connect_same_component ? 0
                                                     : comp_labels[nd.row * ncols + nd.col];
            int sz  = request.connect_same_component ? 0 : comp_size[lbl];
            eps.push_back({nd.row, nd.col, static_cast<double>(elev), lbl, sz});
        }
    }

    std::cout << "\n--- Elevation-based skeleton connection ---\n";
    std::cout << "  endpoints found : " << eps.size() << "\n";
    std::cout << "  search radius   : " << request.search_radius_metres << " m ("
              << search_px << " px)\n";
    std::cout << "  elev tolerance  : " << request.elevation_tolerance_metres << " m\n";

    // ------------------------------------------------------------------
    // 3. Build candidate pairs
    //    Sort: largest segment in the pair first (descending), then by
    //    elevation difference (ascending).  This ensures the longest —
    //    i.e., clearest — terrace-edge traces get first pick of matches.
    // ------------------------------------------------------------------
    struct Candidate { int i, j; double elev_diff; int max_seg_size; };
    std::vector<Candidate> candidates;
    const int n = static_cast<int>(eps.size());
    const double search_px2 = static_cast<double>(search_px) * search_px;

    for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
    {
        if (!request.connect_same_component && eps[i].label == eps[j].label)
            continue;
        double dr = eps[i].r - eps[j].r;
        double dc = eps[i].c - eps[j].c;
        if (dr*dr + dc*dc > search_px2) continue;
        double de = std::abs(eps[i].elev - eps[j].elev);
        if (de > request.elevation_tolerance_metres) continue;
        int msz = std::max(eps[i].seg_size, eps[j].seg_size);
        candidates.push_back({i, j, de, msz});
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b)
              {
                  if (a.max_seg_size != b.max_seg_size)
                      return a.max_seg_size > b.max_seg_size; // longer first
                  return a.elev_diff < b.elev_diff;           // then closer elevation
              });

    std::cout << "  candidates      : " << candidates.size() << "\n";

    // ------------------------------------------------------------------
    // 4. Greedy matching (each endpoint at most one bridge)
    // ------------------------------------------------------------------
    struct Bridge { int i, j; double z_target; };
    std::vector<bool>   used(n, false);
    std::vector<Bridge> bridges;
    for (const auto& cand : candidates)
    {
        if (used[cand.i] || used[cand.j]) continue;
        used[cand.i] = used[cand.j] = true;
        bridges.push_back({cand.i, cand.j,
                           (eps[cand.i].elev + eps[cand.j].elev) * 0.5});
    }

    std::cout << "  bridges accepted: " << bridges.size() << "\n";

    // ------------------------------------------------------------------
    // 5. Route each bridge along the target contour via A* least-cost path.
    //
    //    Per-bridge support raster:
    //      support(r,c) = max(0, 1 - |DEM(r,c) - z_target| / tolerance)
    //    → 1.0 exactly on the target contour, 0.0 at ±tolerance from it.
    //    The A* router minimises  base_cost + support_penalty * (1 - support),
    //    so it naturally threads along the target elevation.
    // ------------------------------------------------------------------
    elevation_connect_result result;

    // Initialise the cumulative bridge mask (all zeros, same geometry)
    {
        std::vector<float>   bd(nrows * ncols, 0.0f);
        std::vector<uint8_t> bm = skeleton.get_mask();
        result.bridge_mask = lsdtt_raster_cpu(std::move(bd), std::move(bm), meta);
    }

    // working_skeleton accumulates routed bridges so successive pairs see
    // each other's pixels as already-foreground (lower traversal cost)
    lsdtt_raster_cpu working_skeleton = skeleton;

    const double tol = request.elevation_tolerance_metres;
    const float  dem_nd_f = static_cast<float>(dem.get_meta().no_data_value);

    skeleton_bridge_route_request route_req;
    route_req.corridor_half_width_pixels = request.corridor_half_width_pixels;
    route_req.base_step_cost             = request.base_step_cost;
    route_req.support_penalty_weight     = request.support_penalty_weight;
    route_req.foreground_cost_scale      = request.foreground_cost_scale;

    // Pre-allocate ONE support raster (zeros everywhere).  Per bridge we fill
    // only the corridor bounding box, route, then reset those pixels.  This
    // avoids allocating + zeroing a full nrows×ncols array 646 times.
    lsdtt_raster_cpu support_raster(
        std::vector<float>(nrows * ncols, 0.0f),
        std::vector<uint8_t>(skeleton.get_mask()),
        meta);
    std::vector<float>& sup_data = support_raster.get_data();

    const int corridor_pad =
        static_cast<int>(std::ceil(request.corridor_half_width_pixels)) + 1;

    const int n_bridges = static_cast<int>(bridges.size());
    int n_routed = 0, n_failed = 0;
    for (int bi = 0; bi < n_bridges; ++bi)
    {
        const Bridge& br = bridges[bi];
        if (bi % 10 == 0)
        {
            std::cout << "  routing bridge " << (bi + 1) << " / " << n_bridges
                      << "  (routed=" << n_routed
                      << " failed=" << n_failed << ")\r" << std::flush;
        }

        const double z_target = br.z_target;

        // Bounding box of the A* corridor (same extent the router will search)
        const int r_lo = std::max(0,        std::min(eps[br.i].r, eps[br.j].r) - corridor_pad);
        const int r_hi = std::min(nrows - 1, std::max(eps[br.i].r, eps[br.j].r) + corridor_pad);
        const int c_lo = std::max(0,        std::min(eps[br.i].c, eps[br.j].c) - corridor_pad);
        const int c_hi = std::min(ncols - 1, std::max(eps[br.i].c, eps[br.j].c) + corridor_pad);

        // Fill support within the corridor bounding box only
        for (int r = r_lo; r <= r_hi; ++r)
        for (int c = c_lo; c <= c_hi; ++c)
        {
            int k = r * ncols + c;
            if (dem_data[k] == dem_nd_f) { sup_data[k] = 0.0f; continue; }
            double dev = std::abs(static_cast<double>(dem_data[k]) - z_target) / tol;
            sup_data[k] = static_cast<float>(std::max(0.0, 1.0 - dev));
        }

        skeleton_endpoint_pair_candidate pair{};
        pair.row_a = eps[br.i].r;
        pair.col_a = eps[br.i].c;
        pair.row_b = eps[br.j].r;
        pair.col_b = eps[br.j].c;

        skeleton_bridge_route_result rr =
            route_skeleton_endpoint_bridges_least_cost(
                working_skeleton, support_raster, {pair}, route_req);

        // Reset the corridor bounding box to zero before the next bridge
        for (int r = r_lo; r <= r_hi; ++r)
        for (int c = c_lo; c <= c_hi; ++c)
            sup_data[r * ncols + c] = 0.0f;

        if (rr.n_routed > 0)
        {
            ++n_routed;
            const auto& new_b = rr.bridge_mask.get_data();
            auto& bdata = result.bridge_mask.get_data();
            auto& wdata = working_skeleton.get_data();
            for (int k = 0; k < nrows * ncols; ++k)
                if (new_b[k] > 0.5f) { bdata[k] = 1.0f; wdata[k] = 1.0f; }
        }
        else
        {
            ++n_failed;
        }
    }

    std::cout << "\n  routed OK : " << n_routed << "  failed: " << n_failed << "\n";

    result.skeleton_connected = working_skeleton;
    return result;
}

} // namespace lsdtt
