// lsdtt_slope_break_finder_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_slope_break_finder_cpu.hpp"
#include "lsdtt_raster_cell_predicates_cpu.hpp"
#include "../utils/lsdtt_connected_components.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "lsdtt_raster_image_ops_cpu.hpp"

#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_skeleton_graph_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_algebra.hpp"
#include "../core/lsdtt_data_loader.hpp"



namespace lsdtt {

// -----------------------------------------------------------------------------
// mask_to_flat
//
// Restricts a curvature raster to pixels where the flat mask is foreground.
// -----------------------------------------------------------------------------
[[maybe_unused]] static lsdtt_raster_cpu mask_to_flat(
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

// -----------------------------------------------------------------------------
// parse_skeleton_edge_metric
// -----------------------------------------------------------------------------
skeleton_edge_metric parse_skeleton_edge_metric(const std::string& name)
{
    if (name == "vertical_excess_curvature")  return skeleton_edge_metric::vertical_excess_curvature;
    if (name == "curvedness")                 return skeleton_edge_metric::curvedness;
    if (name == "maximum_curvature")          return skeleton_edge_metric::maximum_curvature;
    if (name == "multidirectional_hillshade") return skeleton_edge_metric::multidirectional_hillshade;
    throw std::invalid_argument(
        "Unknown edge metric '" + name + "': use vertical_excess_curvature, "
        "curvedness, maximum_curvature or multidirectional_hillshade.");
}

// -----------------------------------------------------------------------------
// threshold_local_mean_std
//
// Foreground (1) where the field exceeds mean + k * max(sd, sd_floor) over the
// valid cells of the block x block square centred on it; summed-area tables of
// the value, its square and the valid count, so the cost does not depend on
// the block size. Computed in double: the field's own units, no rescaling.
//
// sd_floor = (field max - field min) / 255, one grey level of the 8-bit
// methods. Without it the threshold has no scale: on a smooth tread the local
// sd is rounding noise (1e-9 on a synthetic plane) and half the tread passes.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu threshold_local_mean_std(const lsdtt_raster_cpu& field,
                                                 int block, double k)
{
    if (block < 3 || block % 2 == 0)
        throw std::invalid_argument("threshold_local_mean_std: block must be odd and >= 3");
    const lsdtt_raster_meta& meta = field.get_meta();
    const int R = meta.n_rows, C = meta.n_cols, h = block / 2;
    const std::vector<float>&   d = field.get_data();
    const std::vector<uint8_t>& m = field.get_mask();
    const std::size_t W = static_cast<std::size_t>(C) + 1;
    std::vector<double> s1((static_cast<std::size_t>(R) + 1) * W, 0.0), s2(s1.size(), 0.0);
    std::vector<std::int64_t> sn(s1.size(), 0);
    double vmin = std::numeric_limits<double>::infinity(), vmax = -vmin;
    for (std::size_t i = 0; i < d.size(); ++i)
        if (m[i] && std::isfinite(d[i])) { vmin = std::min(vmin, double(d[i])); vmax = std::max(vmax, double(d[i])); }
    const double sd_floor = (vmax > vmin) ? (vmax - vmin) / 255.0 : 0.0;
    for (int r = 0; r < R; ++r)
    {
        double a1 = 0.0, a2 = 0.0;
        std::int64_t an = 0;
        for (int c = 0; c < C; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * C + c;
            if (m[i] && std::isfinite(d[i])) { a1 += d[i]; a2 += double(d[i]) * d[i]; ++an; }
            const std::size_t o = (static_cast<std::size_t>(r) + 1) * W + c + 1;
            s1[o] = s1[o - W] + a1;
            s2[o] = s2[o - W] + a2;
            sn[o] = sn[o - W] + an;
        }
    }
    std::vector<float> out(d.size(), 0.0f);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int r = 0; r < R; ++r)
    {
        const std::size_t r0 = static_cast<std::size_t>(std::max(0, r - h));
        const std::size_t r1 = static_cast<std::size_t>(std::min(R, r + h + 1));
        for (int c = 0; c < C; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * C + c;
            if (!m[i] || !std::isfinite(d[i])) continue;
            const std::size_t c0 = static_cast<std::size_t>(std::max(0, c - h));
            const std::size_t c1 = static_cast<std::size_t>(std::min(C, c + h + 1));
            auto box = [&](const auto& t) {
                return t[r1 * W + c1] - t[r0 * W + c1] - t[r1 * W + c0] + t[r0 * W + c0]; };
            const double n = static_cast<double>(box(sn));
            const double mean = box(s1) / n;
            const double var  = std::max(0.0, box(s2) / n - mean * mean);
            if (d[i] > mean + k * std::max(std::sqrt(var), sd_floor)) out[i] = 1.0f;
        }
    }
    return lsdtt_raster_cpu(std::move(out), std::vector<uint8_t>(m), meta);
}

// =============================================================================
// find_vec_skeleton
// =============================================================================
vec_skeleton_result find_vec_skeleton(
    const lsdtt_raster_cpu&     dem,
    const lsdtt_raster_cpu&     mask,
    const vec_skeleton_request& request)
{
    vec_skeleton_result result;

    // ------------------------------------------------------------------
    // 1. The edge field, masked to flat pixels
    // ------------------------------------------------------------------
    lsdtt_raster_cpu field;
    if (request.edge_metric == skeleton_edge_metric::multidirectional_hillshade)
    {
        std::cout << "\n--- Edge field: multidirectional hillshade ("
                  << request.hillshade_n_directions << " directions) ---\n";
        lsdtt_raster_cpu hs = compute_multidirectional_hillshade(
            dem, request.hillshade_n_directions);
        std::vector<float> v = hs.get_data();
        for (float& x : v) x = -x;   // dark (steep) is high
        field = lsdtt_raster_cpu(std::move(v), std::vector<uint8_t>(hs.get_mask()),
                                 hs.get_meta());
    }
    else
    {
        polyfit_metric metric = polyfit_metric::vertical_excess_curvature;
        const char* label = "vertical excess curvature";
        if (request.edge_metric == skeleton_edge_metric::curvedness)
        { metric = polyfit_metric::curvedness; label = "curvedness"; }
        else if (request.edge_metric == skeleton_edge_metric::maximum_curvature)
        { metric = polyfit_metric::maximum_curvature; label = "|maximum curvature|"; }
        std::cout << "\n--- Edge field: " << label << ", polyfit (lambda="
                  << request.polyfit_lambda_metres << " m) ---\n";

        poly_fit_coeffs coeffs = direct_polyfit(
            dem, request.polyfit_lambda_metres, request.n_threads);
        polyfit_metric_request pmr;
        pmr.metrics   = { metric };
        pmr.n_threads = request.n_threads;
        polyfit_metric_result pmres = compute_polyfit_metrics(coeffs, pmr);
        field = pmres.at(metric);
        if (request.edge_metric == skeleton_edge_metric::maximum_curvature)
        {
            std::vector<float> v = field.get_data();
            for (float& x : v) x = std::fabs(x);
            field = lsdtt_raster_cpu(std::move(v), std::vector<uint8_t>(field.get_mask()),
                                     field.get_meta());
        }
    }
    result.edge_field = mask_to_flat(field, mask);

    std::cout << "  edge field valid flat pixels: "
              << count_foreground(result.edge_field) << "\n";

    // ------------------------------------------------------------------
    // 2. Adaptive threshold -> binary mask
    // ------------------------------------------------------------------
    std::cout << "\n--- Adaptive thresholding ---\n";
    std::cout << "  method     : " << request.adaptive_method << "\n";
    std::cout << "  block size : " << request.adaptive_block_size << " px\n";

    if (request.adaptive_method == "mean_std")
    {
        std::cout << "  std factor : " << request.adaptive_std_factor << "\n";
        result.adaptive_binary = threshold_local_mean_std(
            result.edge_field, request.adaptive_block_size, request.adaptive_std_factor);
    }
    else if (request.adaptive_method == "mean" || request.adaptive_method == "gaussian")
    {
        std::cout << "  C          : " << request.adaptive_C << " grey levels\n";
        std::cout << "  pre-blur sigma : " << request.adaptive_gaussian_sigma << "\n";
        result.adaptive_binary = (request.adaptive_method == "mean")
            ? image_ops_cpu::threshold_adaptive_mean(
                  result.edge_field, request.adaptive_block_size, request.adaptive_C,
                  request.adaptive_gaussian_sigma, request.n_threads)
            : image_ops_cpu::threshold_adaptive_gaussian(
                  result.edge_field, request.adaptive_block_size, request.adaptive_C,
                  request.adaptive_gaussian_sigma, request.n_threads);
    }
    else
    {
        throw std::invalid_argument("find_vec_skeleton: adaptive_method '" +
                                    request.adaptive_method +
                                    "' is not mean_std, mean or gaussian.");
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
    // 4. Zhang-Suen thinning -> spur pruning -> component filtering
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
                  << e.what() << " - using vec_threshold=" << vec_thr << "\n";
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
std::vector<int> label_components(const lsdtt_raster_cpu& skeleton)
{
    // Foreground = a valid cell (mask set, not nodata, finite) above 0.5
    // (decision D3, 2026-10-08). Before, any cell not <= 0.5 counted, so a
    // masked-out cell holding the raster's nodata value counted as skeleton
    // whenever that value was above 0.5 (prune_terminal_loops and
    // prune_skeleton_spurs fill invalid cells with it), and NaN always did.
    const lsdtt_raster_meta&    meta = skeleton.get_meta();
    const std::vector<float>&   data = skeleton.get_data();
    const std::vector<uint8_t>& mask = skeleton.get_mask();
    std::vector<std::int32_t> label;
    std::vector<std::int64_t> size;
    label_connected_components(meta.n_rows, meta.n_cols, 8,
        [&](std::size_t i) { return raster_cell_is_foreground(meta, data, mask, i); },
        label, size);
    return std::vector<int>(label.begin(), label.end());
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

    const std::vector<float>& dem_data  = dem.get_data();

    // ------------------------------------------------------------------
    // 1. Label connected components and compute per-component pixel counts
    // ------------------------------------------------------------------
    std::vector<int> comp_labels;
    std::unordered_map<int,int> comp_size;

    if (!request.connect_same_component)
    {
        comp_labels = label_components(skeleton);
        for (int lbl : comp_labels)
            if (lbl > 0) ++comp_size[lbl];
    }

    // ------------------------------------------------------------------
    // 2. Extract endpoints: skeleton pixels with exactly one 8-neighbour
    // ------------------------------------------------------------------
    // 2. Extract endpoints using crossing-number classification.
    //    classify_skeleton_nodes uses CN=1 (not raw degree=1) so corner
    //    endpoints - e.g. south + southeast neighbours - are not missed.
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
    //    elevation difference (ascending).  This ensures the longest -
    //    i.e., clearest - terrace-edge traces get first pick of matches.
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
    //    -> 1.0 exactly on the target contour, 0.0 at +/-tolerance from it.
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
    // avoids allocating + zeroing a full nrowsxncols array 646 times.
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
