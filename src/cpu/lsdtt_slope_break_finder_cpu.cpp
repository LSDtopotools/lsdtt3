// lsdtt_slope_break_finder_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_slope_break_finder_cpu.hpp"

#include <iostream>
#include <stdexcept>
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

    if (request.min_component_pixels > 1)
    {
        result.profile_skeleton = filter_binary_components_by_size(
            result.profile_skeleton, request.min_component_pixels);
        std::cout << "  profile_skeleton after component filter: "
                  << count_foreground(result.profile_skeleton) << "\n";
    }

    return result;
}

} // namespace lsdtt
