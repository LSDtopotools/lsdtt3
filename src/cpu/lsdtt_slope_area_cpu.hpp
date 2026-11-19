// lsdtt_slope_area_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Slope-area routines for lsdtt-slope-area
// (the method: docs/design_and_algorithms/slope_area_and_low_concavity.md).
//
// d8_flow_path_gradient: the gradient at EVERY valid pixel (hillslopes too)
// along its D8 flow path, over a window CENTRED on the pixel of total flow
// length window_metres:
//   - downstream half: the D8 receivers, to window/2 or an outlet;
//   - upstream half: at each step the D8 donor of LARGEST drainage area (the
//     area raster passed in, so a D-infinity area when that is the run's
//     choice; ties: the smaller flat index), to window/2 or a divide;
//   - the slope is the least-squares fit of z against flow distance
//     (fit_along_path_slope, via channel_slope_at, so lsdtt3 has one
//     definition): at least 3 valid cells, spanning at least window/2 less
//     one diagonal step (each side's walk stops before passing window/2, so
//     on coarse cells window/2 itself is out of reach; plan section 10, F1),
//     slope > 0. A one-sided window (outlet, divide) is used if it alone
//     spans that;
//   - exclude (optional, nonzero = excluded, e.g. cells the conditioning
//     carved or filled): those cells are walked through but not fitted, and
//     their own gradient is nodata.
// Pixels without a usable fit are nodata. The window must be at least four
// cell sizes (channel_slope_at needs a half-length of two). Parallel over
// pixels; the result does not depend on n_threads.
//
// GDAL-free, except write_low_concavity_segments (lsdtt_data_loader).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_histogram.hpp"

namespace lsdtt
{

// slope_area_log_bins: lsdtt-slope-area's binning of a slope raster by a
// drainage-area raster (same grid). The pixels used are valid in both, with
// area > 0 and both finite, and not excluded (exclude: nonzero = excluded,
// e.g. carved / filled pixels). They are binned in n_bins bins of equal width
// in ln(area) with the d05, d25, d50, d75 and d95 slope quantiles
// (compute_log_binned_quantiles, float form, so the CLI's lambda CSVs stay
// byte-identical). Throws std::runtime_error("No valid slope-area pixels were
// found") if no pixel is usable or all areas are equal, and
// std::invalid_argument if the grids differ or n_bins < 1.
inline const std::vector<double> k_slope_area_quantile_levels{0.05, 0.25, 0.50, 0.75, 0.95};

[[nodiscard]] lsdtt_log_binned_quantiles slope_area_log_bins(
    const lsdtt_raster_cpu&          drainage_area,
    const lsdtt_raster_cpu&          slope,
    int                              n_bins,
    const std::vector<std::uint8_t>* exclude = nullptr);

struct lsdtt_d8_flow_path_gradient_options
{
    double                           window_metres = 150.0;
    const std::vector<std::uint8_t>* exclude = nullptr;   // one byte per cell; nonzero = excluded
    int                              n_threads = 0;       // 0 = OpenMP default
};

// routing_dem: the surface the flow was routed on (elevations of the fit).
// drainage_area: same grid; ranks the upstream donors.
// Throws std::invalid_argument if the window is under four cell sizes or a
// raster is not on the flow's grid.
[[nodiscard]] lsdtt_raster_cpu d8_flow_path_gradient(
    const lsdtt_flow_base&                     flow,
    const lsdtt_raster_cpu&                    routing_dem,
    const lsdtt_raster_cpu&                    drainage_area,
    const lsdtt_d8_flow_path_gradient_options& options);

// CONDITIONED pixels (plan 3A, section 10 F2): the cells the conditioning
// changed. FILLED (k_pixel_filled = 1): the routing surface is above the
// input DEM (a pit fill, or a single-cell pit the least-cost breach raised);
// CARVED (k_pixel_carved = 2): below it (a breach channel). Both are excluded
// from lsdtt-slope-area's measurements (Simon, 2026-11-18): a carved reach
// has the breach's own slope, not the landscape's. Strict comparisons; cells
// invalid in either raster are 0. One byte per cell. Throws
// std::invalid_argument if the grids differ.
inline constexpr std::uint8_t k_pixel_untouched = 0;
inline constexpr std::uint8_t k_pixel_filled = 1;
inline constexpr std::uint8_t k_pixel_carved = 2;

[[nodiscard]] std::vector<std::uint8_t> conditioned_pixel_mask(
    const lsdtt_raster_cpu& input_dem,
    const lsdtt_raster_cpu& routing_dem);

// -----------------------------------------------------------------------------
// LOW CONCAVITY (plan section 4): reaches whose gradient changes little with
// drainage area -- threshold slopes, debris-flow channels, alluvial fans.
//
// fit_log_slope_area: least squares log10 S = a + b log10 A over the pairs
// with A > 0, S > 0 (both finite); theta = -b, so S ~ A^-theta.
//   valid: at least min_pixels pairs and the log10 A values not all equal;
//   r2: 1 - SS_res / SS_tot of log10 S (NaN when log10 S is constant, as for
//       a constant gradient, where theta = 0 exactly);
//   log10_area_range: max - min log10 A of the pairs used;
//   theta_se: the standard error of theta (of the fitted slope b),
//       sqrt(SS_res / (n - 2) / sum (x - mean x)^2) with x = log10 A; 0 for a
//       constant gradient. Reported, not used by the flag (Simon 2026-11-18).
struct lsdtt_log_slope_area_fit
{
    double       theta = std::numeric_limits<double>::quiet_NaN();
    double       theta_se = std::numeric_limits<double>::quiet_NaN();
    double       r2 = std::numeric_limits<double>::quiet_NaN();
    double       log10_area_range = std::numeric_limits<double>::quiet_NaN();
    std::int32_t n_pixels = 0;
    bool         valid = false;
};

[[nodiscard]] lsdtt_log_slope_area_fit fit_log_slope_area(
    const std::vector<double>& area,
    const std::vector<double>& slope,
    int                        min_pixels);

// find_low_concavity: at each CHANNEL pixel (D8 contributing pixels >=
// threshold_contributing_pixels; the flow needs its s-vector), take the pixel
// and the next ones DOWNSTREAM along the D8 receivers -- window_pixels in all,
// fewer at an outlet -- and fit log S against log A over those with a valid
// gradient that are not excluded (carved / filled). Downstream because the
// path upstream branches. The window runs through tributary junctions;
// n_junctions counts the window pixels after the first that have two or more
// channel donors (where area jumps). An excluded channel pixel gets no fit.
// flag: 1 low concavity (valid and |theta| < max_abs_concavity, so a weakly
// convex reach counts too), 0 not, -1 no valid fit.
struct lsdtt_low_concavity_options
{
    std::uint64_t                    threshold_contributing_pixels = 1000;
    int                              window_pixels = 10;
    int                              min_pixels = 5;
    double                           max_abs_concavity = 0.15;
    const std::vector<std::uint8_t>* exclude = nullptr;   // nonzero = excluded
    int                              n_threads = 0;
};

// One entry per channel pixel, in ascending flat order.
struct lsdtt_low_concavity_result
{
    std::vector<std::int64_t> flat;
    std::vector<double>       area_m2;
    std::vector<double>       slope;               // the pixel's own gradient (NaN if none)
    std::vector<double>       theta;
    std::vector<double>       theta_se;            // standard error of theta
    std::vector<double>       r2;
    std::vector<double>       log10_area_range;
    std::vector<std::int32_t> n_pixels;            // pairs in the fit
    std::vector<std::int32_t> n_junctions;
    std::vector<std::int8_t>  flag;                // 1, 0, -1
};

// gradient: the D8 flow-path gradient raster (d8_flow_path_gradient);
// drainage_area: the run's area raster (m^2). Throws std::invalid_argument on
// bad options (window_pixels < 2, min_pixels < 3 or > window_pixels,
// max_abs_concavity < 0) or rasters off the grid, or a flow without its
// s-vector.
[[nodiscard]] lsdtt_low_concavity_result find_low_concavity(
    const lsdtt_flow_base&             flow,
    const lsdtt_raster_cpu&            drainage_area,
    const lsdtt_raster_cpu&            gradient,
    const lsdtt_low_concavity_options& options);

// -----------------------------------------------------------------------------
// LOW-CONCAVITY SEGMENTS (plan section 4, step 4).
//
// The channel pixels of `result` are cut into LINKS: a link starts at a
// source (no channel donor) or just below a junction (a pixel with two or
// more channel donors starts its own link) and runs down the D8 receivers to
// the pixel above the next junction, or an outlet. Every channel pixel is in
// exactly one link. Each link is then split where the flag changes, so each
// segment is a run of consecutive pixels with one flag (1, 0 or -1). Before
// splitting, a run of flag 1 shorter than min_run_pixels is relabelled 0
// (merged with its neighbours, not dropped). Pixels with include[f] == 0
// (e.g. outside a basin clip) are left out and break runs.
//
// Per segment: its pixels in downstream order; flag; Strahler order of its
// pixels (the largest; computed over the channel network, from the smallest
// contributing area up); median theta (valid fits only, NaN if none) and
// median slope; min / max area; length_m (D8 steps from the first pixel to
// the receiver of the last, so a one-pixel segment has a length);
// upstream_flat / downstream_flat (first / last pixel).
struct lsdtt_low_concavity_segment
{
    std::int32_t              segment_id = -1;
    std::int32_t              link_id = -1;
    std::int8_t               flag = -1;
    std::int32_t              stream_order = 0;
    std::vector<std::int64_t> flats;           // downstream order
    std::int64_t              end_flat = -1;   // receiver of the last pixel (-1 at an outlet)
    double                    length_m = 0.0;
    double                    median_theta = std::numeric_limits<double>::quiet_NaN();
    double                    median_theta_se = std::numeric_limits<double>::quiet_NaN();
    double                    median_slope = std::numeric_limits<double>::quiet_NaN();
    double                    min_area_m2 = std::numeric_limits<double>::quiet_NaN();
    double                    max_area_m2 = std::numeric_limits<double>::quiet_NaN();
};

[[nodiscard]] std::vector<lsdtt_low_concavity_segment> build_low_concavity_segments(
    const lsdtt_flow_base&             flow,
    const lsdtt_low_concavity_result&  result,
    int                                min_run_pixels = 1,
    const std::vector<std::uint8_t>*   include = nullptr);

// One polyline per segment through the pixel centres (and the receiver of the
// last pixel), attributes: segment_id, link_id, low_concavity (the flag),
// n_pixels, length_m, median_theta, median_theta_se, median_slope, min_area_m2, max_area_m2,
// stream_order, upstream_flat, downstream_flat. fmt: flatgeobuf or gpkg.
// (This writer goes through lsdtt_data_loader; the rest of the module is
// GDAL-free.)
void write_low_concavity_segments(
    const std::string&                              filename,
    const std::vector<lsdtt_low_concavity_segment>& segments,
    const lsdtt_raster_meta&                        meta,
    vector_output_format                            fmt,
    bool                                            spatial_index = false);

} // namespace lsdtt
