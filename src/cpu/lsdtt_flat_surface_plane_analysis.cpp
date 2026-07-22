// lsdtt_flat_surface_plane_analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_flat_surface_plane_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lsdtt {
namespace {

static constexpr double k_rad2deg = 57.2957795130823208768;
static constexpr double k_nan     = std::numeric_limits<double>::quiet_NaN();

// =============================================================================
// label_accum  - per-label sufficient statistics accumulated across passes
// =============================================================================
struct label_accum
{
    int    n     = 0;
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    double mean_x = 0.0, mean_y = 0.0, mean_z = 0.0;

    // Centered second-order sums (filled in pass 2)
    double Sxx = 0.0, Syy = 0.0, Szz = 0.0;
    double Sxy = 0.0, Sxz = 0.0, Syz = 0.0;
};

// =============================================================================
// solve_plane_2x2
//
// Solves the centred normal equations for a least-squares plane
//   z' = a*x' + b*y'
// where primes denote coordinates centred on the label mean.
//
// The system is:
//   [ Sxx  Sxy ] [ a ]   [ Sxz ]
//   [ Sxy  Syy ] [ b ] = [ Syz ]
//
// Returns true and writes a, b on success.
// Returns false if the system is degenerate (|det| < eps).
// =============================================================================
static bool solve_plane_2x2(
    double Sxx, double Sxy, double Syy,
    double Sxz, double Syz,
    double& a_out, double& b_out)
{
    const double det = Sxx * Syy - Sxy * Sxy;
    constexpr double eps = 1e-12;
    if (std::fabs(det) < eps) return false;

    a_out = (Sxz * Syy - Syz * Sxy) / det;
    b_out = (Syz * Sxx - Sxz * Sxy) / det;
    return true;
}

// =============================================================================
// detect_bimodal
//
// Peak-valley analysis on a histogram of per-label residuals.
// Finds the two tallest peaks, measures the valley between them.
//
// Returns (is_bimodal, confidence) where
//   confidence = 1 - valley_min / min(peak1_count, peak2_count)
//   is_bimodal  = confidence > bimodal_threshold
// =============================================================================
static std::pair<bool, double> detect_bimodal(
    const std::vector<int>& hist,
    double                  bimodal_threshold)
{
    const int nb = static_cast<int>(hist.size());
    if (nb < 4) return {false, 0.0};

    std::vector<int> peaks;
    peaks.reserve(8);

    if (hist[0] > hist[1] && hist[0] > 0)
        peaks.push_back(0);

    for (int i = 1; i < nb - 1; ++i)
    {
        if (hist[i] >= hist[i-1] && hist[i] >= hist[i+1] && hist[i] > 0)
            peaks.push_back(i);
    }

    if (hist[nb-1] > hist[nb-2] && hist[nb-1] > 0)
        peaks.push_back(nb - 1);

    if (static_cast<int>(peaks.size()) < 2) return {false, 0.0};

    std::partial_sort(peaks.begin(),
                      peaks.begin() + 2,
                      peaks.end(),
                      [&](int a, int b){ return hist[a] > hist[b]; });

    int pk_lo = std::min(peaks[0], peaks[1]);
    int pk_hi = std::max(peaks[0], peaks[1]);

    if (pk_lo == pk_hi) return {false, 0.0};

    int valley_min = hist[pk_lo];
    for (int i = pk_lo + 1; i < pk_hi; ++i)
        valley_min = std::min(valley_min, hist[i]);

    const int lower_peak = std::min(hist[pk_lo], hist[pk_hi]);
    if (lower_peak == 0) return {false, 0.0};

    const double confidence = 1.0 - static_cast<double>(valley_min)
                                   / static_cast<double>(lower_peak);

    return {confidence > bimodal_threshold, confidence};
}

// =============================================================================
// Cross-surface residual helpers
//
// For plane i:
//   z = a_i*x + b_i*y + c_i
//
// evaluated on label j, define the mean residual:
//
//   mean_rji = mean_zj - a_i*mean_xj - b_i*mean_yj - c_i
//
// and the SSE:
//
//   SSE_{j<-i} =
//       Szz_j
//     + a_i^2 Sxx_j
//     + b_i^2 Syy_j
//     + 2 a_i b_i Sxy_j
//     - 2 a_i Sxz_j
//     - 2 b_i Syz_j
//     + n_j mean_rji^2
//
// where all S terms are centred on label j's own mean coordinates.
// This gives O(1) cross-surface RMSE per pair without rescanning pixels.
// =============================================================================
static double cross_mean_residual(
    const label_plane_fit& plane,
    const label_accum&     surface)
{
    return surface.mean_z
         - plane.a * surface.mean_x
         - plane.b * surface.mean_y
         - plane.c;
}

static double cross_sse(
    const label_plane_fit& plane,
    const label_accum&     surface)
{
    if (!plane.valid || surface.n <= 0) return k_nan;

    const double mean_r = cross_mean_residual(plane, surface);

    const double sse =
          surface.Szz
        + plane.a * plane.a * surface.Sxx
        + plane.b * plane.b * surface.Syy
        + 2.0 * plane.a * plane.b * surface.Sxy
        - 2.0 * plane.a * surface.Sxz
        - 2.0 * plane.b * surface.Syz
        + static_cast<double>(surface.n) * mean_r * mean_r;

    return std::max(0.0, sse);
}

static double cross_rmse(
    const label_plane_fit& plane,
    const label_accum&     surface)
{
    const double sse = cross_sse(plane, surface);
    if (!std::isfinite(sse) || surface.n <= 0) return k_nan;
    return std::sqrt(sse / static_cast<double>(surface.n));
}

// =============================================================================
// union-find helpers
// =============================================================================
static int uf_find(std::vector<int>& parent, int x)
{
    while (parent[x] != x)
    {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static void uf_unite(std::vector<int>& parent, std::vector<int>& rank,
                     int x, int y)
{
    x = uf_find(parent, x);
    y = uf_find(parent, y);
    if (x == y) return;
    if (rank[x] < rank[y]) std::swap(x, y);
    parent[y] = x;
    if (rank[x] == rank[y]) ++rank[x];
}

} // anonymous namespace


// =============================================================================
// analyse_flat_surface_planes
// =============================================================================
plane_analysis_result analyse_flat_surface_planes(
    const lsdtt_raster_cpu&       dem,
    const lsdtt_raster_cpu&       labelled_mask,
    const plane_analysis_request& request)
{
    // -------------------------------------------------------------------------
    // Validate
    // -------------------------------------------------------------------------
    const lsdtt_raster_meta& meta = dem.get_meta();

    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::analyse_flat_surface_planes: DEM metadata is invalid"
        );

    const lsdtt_raster_meta& lmeta = labelled_mask.get_meta();
    if (lmeta.n_rows != meta.n_rows || lmeta.n_cols != meta.n_cols)
        throw std::runtime_error(
            "lsdtt::analyse_flat_surface_planes: DEM and labelled_mask "
            "have different dimensions"
        );

    if (request.min_pixels_for_analysis < 3)
        throw std::invalid_argument(
            "lsdtt::analyse_flat_surface_planes: "
            "min_pixels_for_analysis must be >= 3"
        );
    if (request.merge_tolerance_rmse_factor <= 0.0)
        throw std::invalid_argument(
            "lsdtt::analyse_flat_surface_planes: "
            "merge_tolerance_rmse_factor must be > 0"
        );
    if (request.split_histogram_n_bins < 4)
        throw std::invalid_argument(
            "lsdtt::analyse_flat_surface_planes: "
            "split_histogram_n_bins must be >= 4"
        );

    const int    nrows     = meta.n_rows;
    const int    ncols     = meta.n_cols;
    const double cell_size = meta.cell_size;
    const std::size_t npix = static_cast<std::size_t>(nrows) *
                             static_cast<std::size_t>(ncols);

    const std::vector<float>&   dem_data   = dem.get_data();
    const std::vector<uint8_t>& dem_mask   = dem.get_mask();
    const std::vector<float>&   label_data = labelled_mask.get_data();
    const std::vector<uint8_t>& label_mask = labelled_mask.get_mask();

    // -------------------------------------------------------------------------
    // Find the maximum label so we know how many accumulators to allocate.
    // -------------------------------------------------------------------------
    int max_label = 0;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (label_mask[i] == 0) continue;
        const int lbl = static_cast<int>(std::round(label_data[i]));
        if (lbl > max_label) max_label = lbl;
    }

    if (max_label == 0)
    {
        plane_analysis_result empty;
        return empty;
    }

    // -------------------------------------------------------------------------
    // Pass 1 - accumulate n, sum_x, sum_y, sum_z per label.
    // -------------------------------------------------------------------------
    std::vector<label_accum> accums(static_cast<std::size_t>(max_label + 1));

    for (int row = 0; row < nrows; ++row)
    {
        const double y_rel = (static_cast<double>(nrows - row) - 0.5) * cell_size;

        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(ncols) +
                                    static_cast<std::size_t>(col);

            if (dem_mask[idx] == 0 || label_mask[idx] == 0) continue;

            const int lbl = static_cast<int>(std::round(label_data[idx]));
            if (lbl <= 0 || lbl > max_label) continue;

            const float z = dem_data[idx];
            if (!std::isfinite(z) || meta.is_nodata(z)) continue;

            const double x_rel = (static_cast<double>(col) + 0.5) * cell_size;

            label_accum& acc = accums[static_cast<std::size_t>(lbl)];
            acc.n++;
            acc.sum_x += x_rel;
            acc.sum_y += y_rel;
            acc.sum_z += static_cast<double>(z);
        }
    }

    for (int lbl = 1; lbl <= max_label; ++lbl)
    {
        label_accum& acc = accums[static_cast<std::size_t>(lbl)];
        if (acc.n == 0) continue;
        const double inv_n = 1.0 / static_cast<double>(acc.n);
        acc.mean_x = acc.sum_x * inv_n;
        acc.mean_y = acc.sum_y * inv_n;
        acc.mean_z = acc.sum_z * inv_n;
    }

    // -------------------------------------------------------------------------
    // Pass 2 - accumulate centred second-order sums per label.
    // -------------------------------------------------------------------------
    for (int row = 0; row < nrows; ++row)
    {
        const double y_rel = (static_cast<double>(nrows - row) - 0.5) * cell_size;

        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(ncols) +
                                    static_cast<std::size_t>(col);

            if (dem_mask[idx] == 0 || label_mask[idx] == 0) continue;

            const int lbl = static_cast<int>(std::round(label_data[idx]));
            if (lbl <= 0 || lbl > max_label) continue;

            const float z = dem_data[idx];
            if (!std::isfinite(z) || meta.is_nodata(z)) continue;

            const double x_rel = (static_cast<double>(col) + 0.5) * cell_size;

            label_accum& acc = accums[static_cast<std::size_t>(lbl)];
            const double xp = x_rel - acc.mean_x;
            const double yp = y_rel - acc.mean_y;
            const double zp = static_cast<double>(z) - acc.mean_z;

            acc.Sxx += xp * xp;
            acc.Syy += yp * yp;
            acc.Szz += zp * zp;
            acc.Sxy += xp * yp;
            acc.Sxz += xp * zp;
            acc.Syz += yp * zp;
        }
    }

    // -------------------------------------------------------------------------
    // Fit planes - solve 2x2 Cramer per label.
    // -------------------------------------------------------------------------
    std::vector<label_plane_fit> plane_fits(
        static_cast<std::size_t>(max_label + 1));

    for (int lbl = 1; lbl <= max_label; ++lbl)
    {
        label_plane_fit& pf    = plane_fits[static_cast<std::size_t>(lbl)];
        const label_accum& acc = accums[static_cast<std::size_t>(lbl)];

        pf.label    = lbl;
        pf.n_pixels = acc.n;
        pf.mean_z   = acc.mean_z;

        if (acc.n < request.min_pixels_for_analysis)
        {
            pf.valid = false;
            continue;
        }

        double a = 0.0, b = 0.0;
        if (!solve_plane_2x2(acc.Sxx, acc.Sxy, acc.Syy,
                             acc.Sxz, acc.Syz, a, b))
        {
            pf.valid = false;
            continue;
        }

        const double c    = acc.mean_z - a * acc.mean_x - b * acc.mean_y;
        const double sse  = acc.Szz - a * acc.Sxz - b * acc.Syz;
        const double rmse = std::sqrt(std::max(0.0, sse) /
                                      static_cast<double>(acc.n));

        pf.a               = a;
        pf.b               = b;
        pf.c               = c;
        pf.rmse            = rmse;
        pf.slope_magnitude = std::sqrt(a * a + b * b);
        // COMPASS down-slope azimuth: 0 = north, clockwise, [0, 360).
        // a = dz/dx (east-positive), b = dz/dy (north-positive), so the
        // direction the plane FACES is atan2(-a, -b).  Was atan2(b, a) until
        // 2026-07-22, which was the maths convention (0 = east, anticlockwise)
        // AND the up-slope direction - wrong twice over, and it shares the
        // column name `aspect_degrees` with lsdtt-surface-metrics.
        {
            double asp = std::atan2(-a, -b) * k_rad2deg;
            if (asp < 0.0) asp += 360.0;
            pf.aspect_degrees = asp;
        }
        pf.valid           = true;
    }

    // -------------------------------------------------------------------------
    // Pass 3 - per-label residual histograms for bimodal split detection.
    // -------------------------------------------------------------------------
    const int n_bins = request.split_histogram_n_bins;

    std::vector<int> histograms(
        static_cast<std::size_t>(max_label + 1) *
        static_cast<std::size_t>(n_bins), 0);

    for (int row = 0; row < nrows; ++row)
    {
        const double y_rel = (static_cast<double>(nrows - row) - 0.5) * cell_size;

        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(ncols) +
                                    static_cast<std::size_t>(col);

            if (dem_mask[idx] == 0 || label_mask[idx] == 0) continue;

            const int lbl = static_cast<int>(std::round(label_data[idx]));
            if (lbl <= 0 || lbl > max_label) continue;

            const label_plane_fit& pf = plane_fits[static_cast<std::size_t>(lbl)];
            if (!pf.valid || pf.rmse < 1e-12) continue;

            const float z = dem_data[idx];
            if (!std::isfinite(z) || meta.is_nodata(z)) continue;

            const double x_rel = (static_cast<double>(col) + 0.5) * cell_size;
            const double residual = static_cast<double>(z) -
                                    (pf.a * x_rel + pf.b * y_rel + pf.c);

            const double half_range = request.split_histogram_sigma * pf.rmse;
            const double bin_width  = 2.0 * half_range / static_cast<double>(n_bins);
            const int    bin        = static_cast<int>(
                std::floor((residual + half_range) / bin_width));

            if (bin < 0 || bin >= n_bins) continue;

            ++histograms[static_cast<std::size_t>(lbl) *
                         static_cast<std::size_t>(n_bins) +
                         static_cast<std::size_t>(bin)];
        }
    }

    for (int lbl = 1; lbl <= max_label; ++lbl)
    {
        label_plane_fit& pf = plane_fits[static_cast<std::size_t>(lbl)];
        if (!pf.valid) continue;

        const std::size_t hist_offset =
            static_cast<std::size_t>(lbl) * static_cast<std::size_t>(n_bins);
        const std::vector<int> hist(
            histograms.begin() + static_cast<std::ptrdiff_t>(hist_offset),
            histograms.begin() + static_cast<std::ptrdiff_t>(
                hist_offset + static_cast<std::size_t>(n_bins))
        );

        const auto [is_bimodal, confidence] =
            detect_bimodal(hist, request.split_bimodal_threshold);

        pf.split_candidate  = is_bimodal;
        pf.split_confidence = confidence;
    }

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Per-pixel plane residual raster.
    // -------------------------------------------------------------------------
    const float nodata_f = static_cast<float>(meta.no_data_value);

    std::vector<float>   res_data(npix, nodata_f);
    std::vector<uint8_t> res_mask(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        const double y_rel = (static_cast<double>(nrows - row) - 0.5) * cell_size;

        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(ncols) +
                                    static_cast<std::size_t>(col);

            if (dem_mask[idx] == 0 || label_mask[idx] == 0) continue;

            const int lbl = static_cast<int>(std::round(label_data[idx]));
            if (lbl <= 0 || lbl > max_label) continue;

            const label_plane_fit& pf = plane_fits[static_cast<std::size_t>(lbl)];
            if (!pf.valid) continue;

            const float z = dem_data[idx];
            if (!std::isfinite(z) || meta.is_nodata(z)) continue;

            const double x_rel = (static_cast<double>(col) + 0.5) * cell_size;
            const double z_hat = pf.a * x_rel + pf.b * y_rel + pf.c;

            res_data[idx] = static_cast<float>(static_cast<double>(z) - z_hat);
            res_mask[idx] = uint8_t{1};
        }
    }

    // -------------------------------------------------------------------------
    // Dense cross-surface RMSE matrix
    //
    // plane_on_label_rmse[i][j] = RMSE of plane_i evaluated on label_j pixels.
    // -------------------------------------------------------------------------
    std::vector<std::vector<double>> plane_on_label_rmse(
        static_cast<std::size_t>(max_label + 1),
        std::vector<double>(static_cast<std::size_t>(max_label + 1), k_nan)
    );

    for (int i = 1; i <= max_label; ++i)
    {
        const label_plane_fit& pfi = plane_fits[static_cast<std::size_t>(i)];
        if (!pfi.valid) continue;

        for (int j = 1; j <= max_label; ++j)
        {
            const label_accum& aj = accums[static_cast<std::size_t>(j)];
            if (aj.n <= 0) continue;

            plane_on_label_rmse[static_cast<std::size_t>(i)]
                               [static_cast<std::size_t>(j)] =
                cross_rmse(pfi, aj);
        }
    }

    // -------------------------------------------------------------------------
    // Pairwise merge check - O(1) per pair using sufficient statistics.
    // -------------------------------------------------------------------------
    std::vector<label_merge_candidate> merge_candidates;
    merge_candidates.reserve(
        static_cast<std::size_t>(max_label) *
        static_cast<std::size_t>(max_label) / 2
    );

    for (int i = 1; i <= max_label; ++i)
    {
        const label_plane_fit& pfi = plane_fits[static_cast<std::size_t>(i)];
        if (!pfi.valid) continue;

        const label_accum& ai = accums[static_cast<std::size_t>(i)];

        for (int j = i + 1; j <= max_label; ++j)
        {
            const label_plane_fit& pfj = plane_fits[static_cast<std::size_t>(j)];
            if (!pfj.valid) continue;

            const label_accum& aj = accums[static_cast<std::size_t>(j)];

            const double mean_rji = cross_mean_residual(pfi, aj);
            const double mean_rij = cross_mean_residual(pfj, ai);

            const double rmse_ji =
                plane_on_label_rmse[static_cast<std::size_t>(i)]
                                   [static_cast<std::size_t>(j)];
            const double rmse_ij =
                plane_on_label_rmse[static_cast<std::size_t>(j)]
                                   [static_cast<std::size_t>(i)];

            double norm_ji = (pfi.rmse > 1e-12)
                             ? std::fabs(mean_rji) / pfi.rmse
                             : std::fabs(mean_rji);
            double norm_ij = (pfj.rmse > 1e-12)
                             ? std::fabs(mean_rij) / pfj.rmse
                             : std::fabs(mean_rij);

            const double norm_offset = std::max(norm_ji, norm_ij);
            const bool   same_level  =
                norm_offset < request.merge_tolerance_rmse_factor;

            merge_candidates.push_back({
                i, j,
                mean_rji, mean_rij,
                rmse_ji, rmse_ij,
                norm_offset,
                same_level
            });
        }
    }

    std::sort(merge_candidates.begin(), merge_candidates.end(),
              [](const label_merge_candidate& a,
                 const label_merge_candidate& b)
              { return a.normalised_offset < b.normalised_offset; });

    // -------------------------------------------------------------------------
    // Union-find connected components on same_level pairs -> merge groups.
    // -------------------------------------------------------------------------
    std::vector<int> uf_parent(static_cast<std::size_t>(max_label + 1));
    std::vector<int> uf_rank  (static_cast<std::size_t>(max_label + 1), 0);
    std::iota(uf_parent.begin(), uf_parent.end(), 0);

    for (const auto& mc : merge_candidates)
    {
        if (mc.same_level)
            uf_unite(uf_parent, uf_rank, mc.label_i, mc.label_j);
    }

    std::vector<int> root_to_group(static_cast<std::size_t>(max_label + 1), -1);
    std::vector<int> merge_group  (static_cast<std::size_t>(max_label + 1), -1);
    int n_groups = 0;

    for (int lbl = 1; lbl <= max_label; ++lbl)
    {
        if (!plane_fits[static_cast<std::size_t>(lbl)].valid) continue;
        const int root = uf_find(uf_parent, lbl);
        if (root_to_group[static_cast<std::size_t>(root)] < 0)
            root_to_group[static_cast<std::size_t>(root)] = n_groups++;
        merge_group[static_cast<std::size_t>(lbl)] =
            root_to_group[static_cast<std::size_t>(root)];
    }

    // -------------------------------------------------------------------------
    // Assemble result
    // -------------------------------------------------------------------------
    plane_analysis_result result;
    result.plane_fits          = std::move(plane_fits);
    result.merge_candidates    = std::move(merge_candidates);
    result.merge_group         = std::move(merge_group);
    result.n_groups            = n_groups;
    result.plane_on_label_rmse = std::move(plane_on_label_rmse);
    result.plane_residual_raster =
        lsdtt_raster_cpu(std::move(res_data), std::move(res_mask), meta);

    return result;
}



// =============================================================================
// merge_labelled_flat_mask_from_plane_analysis
// =============================================================================
lsdtt_raster_cpu merge_labelled_flat_mask_from_plane_analysis(
    const lsdtt_raster_cpu&    labelled_mask,
    const plane_analysis_result& analysis_result)
{
    const lsdtt_raster_meta& meta = labelled_mask.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::merge_labelled_flat_mask_from_plane_analysis: "
            "labelled_mask metadata is invalid"
        );
    }

    const std::vector<float>&   in_data = labelled_mask.get_data();
    const std::vector<uint8_t>& in_mask = labelled_mask.get_mask();
    const std::size_t npix = in_data.size();

    if (in_mask.size() != npix)
    {
        throw std::runtime_error(
            "lsdtt::merge_labelled_flat_mask_from_plane_analysis: "
            "data/mask size mismatch"
        );
    }

    // Find the largest label present in the raster.
    int max_label = 0;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        const int lbl = static_cast<int>(std::round(in_data[i]));
        if (lbl > max_label) max_label = lbl;
    }

    // Default mapping: identity.
    std::vector<int> label_map(static_cast<std::size_t>(max_label + 1), 0);
    for (int lbl = 0; lbl <= max_label; ++lbl)
        label_map[static_cast<std::size_t>(lbl)] = lbl;

    // Determine, for each merge group, the minimum original label in that group.
    int max_group = -1;
    for (std::size_t lbl = 1; lbl < analysis_result.merge_group.size(); ++lbl)
    {
        if (analysis_result.merge_group[lbl] > max_group)
            max_group = analysis_result.merge_group[lbl];
    }

    if (max_group >= 0)
    {
        std::vector<int> group_min_label(
            static_cast<std::size_t>(max_group + 1),
            std::numeric_limits<int>::max()
        );

        for (std::size_t lbl = 1; lbl < analysis_result.merge_group.size(); ++lbl)
        {
            const int grp = analysis_result.merge_group[lbl];
            if (grp >= 0)
            {
                if (static_cast<int>(lbl) < group_min_label[static_cast<std::size_t>(grp)])
                    group_min_label[static_cast<std::size_t>(grp)] = static_cast<int>(lbl);
            }
        }

        // Build label -> merged-label mapping.
        for (int lbl = 1; lbl <= max_label; ++lbl)
        {
            if (static_cast<std::size_t>(lbl) < analysis_result.merge_group.size())
            {
                const int grp = analysis_result.merge_group[static_cast<std::size_t>(lbl)];
                if (grp >= 0)
                {
                    label_map[static_cast<std::size_t>(lbl)] =
                        group_min_label[static_cast<std::size_t>(grp)];
                }
            }
        }
    }

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask = in_mask;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;

        const int in_label = static_cast<int>(std::round(in_data[i]));
        if (in_label <= 0)
        {
            out_data[i] = 0.0f;
            continue;
        }

        int out_label = in_label;
        if (in_label <= max_label)
            out_label = label_map[static_cast<std::size_t>(in_label)];

        out_data[i] = static_cast<float>(out_label);
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}



} // namespace lsdtt