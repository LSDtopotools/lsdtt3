// lsdtt_raster_linear_features_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Linear-feature primitives — CPU free functions.
//
// Purpose
// -------
// Generic image-processing operators for detecting and enhancing narrow,
// elongated ("curvilinear") structures in an arbitrary continuous scalar
// raster — for example a flow-divergence response, a Frangi vesselness map, a
// curvature magnitude, or any other per-pixel score where the real features of
// interest are continuous threads while the background is isolated speckle.
//
// These operators score (or keep) a pixel by how much support it receives from
// its neighbours ALONG a consistent line, not by its own pointwise shape. That
// is what separates a continuous road/channel/drain from same-brightness
// hillslope speckle.
//
// This module is intentionally independent of any specific geomorphic workflow.
// Feature-scoring pipelines that *compose* these primitives (the road / channel
// / drain detectors) live in lsdtt_raster_feature_finder_cpu.hpp.
//
// What belongs here
// -----------------
//   - grayscale path opening / closing (length-based continuity)
//   - structure tensor of a scalar field: coherence + orientation       (later)
//   - coherence-enhancing anisotropic diffusion                         (later)
//   - oriented non-maximum suppression + hysteresis linking             (later)
//
// What does NOT belong here
// -------------------------
//   - binary 0/1 morphology .............. lsdtt_raster_binary_mask_cpu.hpp
//   - DEM surface-geometry metrics ....... lsdtt_raster_surface_metrics_cpu.hpp
//   - FFT / periodicity / spacing ........ lsdtt_raster_spectral_cpu.hpp
//   - DEM feature scoring pipelines ...... lsdtt_raster_feature_finder_cpu.hpp
//
// Input convention
// ----------------
//   - valid pixels are those with mask != 0 and a finite, non-nodata value
//   - invalid / nodata pixels are excluded from every operation and written
//     back as nodata
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// grayscale_path_opening
//
// Complete grayscale path opening (Appleton & Talbot 2007). Keeps, at each
// pixel, the strongest length-constrained path of high values that passes
// through it; isolated speckle — which lies on no long path — is suppressed to
// the local background level. The result is anti-extensive: output <= input
// everywhere on valid pixels.
//
// Definition
// ----------
// For a path family P (a set of directed near-straight paths on the grid), the
// path opening is
//
//     g(p) = max over paths pi in P with |pi| >= L and p in pi
//                of  ( min over q in pi of f(q) )
//
// i.e. the brightest level at which p still lies on a connected path of at
// least L pixels. A pixel on a long bright thread keeps a high value; a bright
// speckle with no long bright path collapses to background.
//
// Path family (completeness)
// --------------------------
// Two orientation families are used and the results combined pixelwise by max
// (a union of openings):
//   - "horizontalish": each step advances one column and may deviate by -1, 0,
//     or +1 row  -> covers every line within +/-45 deg of horizontal
//   - "verticalish":   each step advances one row    and may deviate by -1, 0,
//     or +1 column -> covers every line within +/-45 deg of vertical
// Together they cover all orientations (the "complete" path opening), and the
// +/-1 transverse wobble lets gently curving features (real roads, meandering
// channels) survive without being broken at every bend.
//
// Implementation
// --------------
// Grayscale is handled by threshold decomposition: the value range over valid
// pixels is split into n_quantization_levels thresholds; at each threshold the
// binary path opening is evaluated via a longest-path-in-DAG computation
// (length through p = longest path ending at p + longest path starting at p
// minus 1); each pixel is assigned the highest threshold at which it still lies
// on a path of length >= L. Levels are evaluated in parallel across n_threads.
//
// Parameters
// ----------
//   response                 — input scalar raster (any continuous field)
//   min_path_length_pixels   — L, minimum path length in pixels. Set near the
//                              shortest real feature you want to keep; larger L
//                              removes more speckle but also drops short
//                              fragments. Must be >= 2.
//   n_quantization_levels    — number of grayscale threshold levels (default
//                              256). Higher = finer value preservation at
//                              higher cost. Must be >= 2.
//   n_threads                — OpenMP thread count; 0 = all available.
//
// Returns
// -------
//   A raster the same shape as the input. Valid pixels hold the path-opened
//   value (<= input); invalid pixels are nodata.
//
// Throws
// ------
//   std::runtime_error    if metadata is invalid.
//   std::invalid_argument if min_path_length_pixels < 2 or
//                         n_quantization_levels < 2.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu grayscale_path_opening(
    const lsdtt_raster_cpu& response,
    int                     min_path_length_pixels,
    int                     n_quantization_levels = 256,
    int                     n_threads             = 0
);

// -----------------------------------------------------------------------------
// grayscale_path_closing
//
// Dual of grayscale_path_opening: a path closing of f, computed as
//   -grayscale_path_opening(-f)
// It fills short dark gaps that lie across an otherwise continuous bright
// thread (e.g. a one-pixel break in a road) while leaving large dark regions
// unchanged. The result is extensive: output >= input on valid pixels.
//
// Parameters and exceptions are identical to grayscale_path_opening.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu grayscale_path_closing(
    const lsdtt_raster_cpu& response,
    int                     min_path_length_pixels,
    int                     n_quantization_levels = 256,
    int                     n_threads             = 0
);

// -----------------------------------------------------------------------------
// structure_tensor_result
//
// Per-pixel outputs of the 2D structure-tensor analysis of a scalar field.
//
//   coherence       — anisotropy / line-likeness in [0, 1]:
//                       (lam1 - lam2) / (lam1 + lam2).
//                       0 = isotropic (speckle, flat, or a corner);
//                       1 = a perfectly oriented edge/line. This is the robust,
//                       neighbourhood-pooled replacement for the pointwise
//                       Hessian elongation: speckle scores low because its
//                       gradients point every way and cancel under averaging,
//                       while a continuous thread scores high.
//
//   orientation_deg — azimuth ALONG which the linear feature runs (its tangent),
//                       in degrees clockwise from north, folded into [0, 180).
//                       Meaningful only where coherence is high. This is the
//                       orientation field consumed by oriented non-maximum
//                       suppression and the regular-spacing detector.
//
//   energy          — lam1 + lam2, the local gradient energy (tensor trace).
//                       Use to distinguish a strong oriented edge from a weak
//                       one: high coherence + low energy is a faint structure.
// -----------------------------------------------------------------------------
struct structure_tensor_result
{
    lsdtt_raster_cpu coherence;
    lsdtt_raster_cpu orientation_deg;
    lsdtt_raster_cpu energy;
};

// -----------------------------------------------------------------------------
// compute_structure_tensor_coherence
//
// Builds the 2D structure (second-moment) tensor of a scalar field and returns
// its coherence, orientation, and energy. This scores a pixel by how
// consistently its NEIGHBOURHOOD is oriented, not by its own pointwise shape —
// which is what separates a continuous road/channel from same-brightness
// hillslope speckle.
//
// Algorithm
// ---------
//   1. Optionally pre-smooth the field at derivative_radius_metres (a circular
//      moving-window mean) to set the gradient noise scale.
//   2. Central-difference gradient (gx, gy) in metres. A pixel is invalid if any
//      pixel used in its difference is nodata.
//   3. Form the per-pixel outer product J = [[gx*gx, gx*gy], [gx*gy, gy*gy]].
//   4. Average each component over a circular window of
//      integration_radius_metres (THE key scale — this pooling is what makes the
//      tensor robust). Use ~ the feature half-width plus a few pixels.
//   5. Eigen-decompose the smoothed 2x2 tensor analytically:
//        lam1,2 = tr/2 +/- sqrt( ((Jxx-Jyy)/2)^2 + Jxy^2 ),  lam1 >= lam2 >= 0.
//      coherence = (lam1 - lam2)/(lam1 + lam2); the feature tangent is
//      perpendicular to the dominant gradient direction.
//
// Parameters
// ----------
//   field                     — input scalar raster (any continuous response).
//   derivative_radius_metres  — pre-smoothing radius for the gradient; 0 = none.
//   integration_radius_metres — tensor-averaging radius; must be > 0.
//   n_threads                 — OpenMP thread count; 0 = all available.
//
// Throws
// ------
//   std::runtime_error    if metadata is invalid.
//   std::invalid_argument if integration_radius_metres <= 0 or
//                         derivative_radius_metres < 0.
// -----------------------------------------------------------------------------
[[nodiscard]] structure_tensor_result compute_structure_tensor_coherence(
    const lsdtt_raster_cpu& field,
    double                  derivative_radius_metres,
    double                  integration_radius_metres,
    int                     n_threads = 0
);

// -----------------------------------------------------------------------------
// coherence_enhancing_diffusion
//
// Weickert (1999) coherence-enhancing anisotropic diffusion. Iteratively smooths
// a scalar field ALONG its dominant local orientation and barely at all ACROSS
// it. The effect is to bridge broken / patchy linear structures (e.g. a road or
// channel that is high but interrupted in the divergence field) and to erase
// speckle perpendicular to them — a denoising pre-filter that makes the
// downstream path opening / thresholding far cleaner.
//
// Method
// ------
// Each iteration:
//   1. Compute the structure tensor of the current field (via
//      compute_structure_tensor_coherence) -> coherence kappa and the local
//      feature orientation.
//   2. Build the diffusion tensor D sharing those eigenvectors, with
//      eigenvalues
//        across the feature : lam_perp  = min_diffusivity            (~0, no smoothing across)
//        along  the feature : lam_par   = min_diffusivity
//                                       + (1 - min_diffusivity) * exp(-contrast / kappa^2)
//      so smoothing along the feature switches on only where the field is
//      coherent (kappa near 1) and stays off in isotropic speckle (kappa near 0).
//   3. Explicit update  u <- u + time_step * div(D grad u)  with no-flux
//      (Neumann) conditions at the grid edge and at any nodata boundary.
//
// Note: kappa here is the NORMALISED coherence in [0,1] (= the coherence raster
// from compute_structure_tensor_coherence), so `contrast` is scale-invariant —
// raise it to be more selective (smooth only the most coherent threads), lower
// it to smooth more of the field.
//
// Parameters
// ----------
//   field                     — input scalar raster.
//   n_iterations              — number of diffusion steps (>= 1). More =
//                               longer-range bridging and smoothing.
//   time_step                 — explicit step tau; must be in (0, 0.25] for
//                               stability (0.15 is a safe default).
//   derivative_radius_metres  — gradient pre-smoothing for the structure tensor;
//                               0 = none.
//   integration_radius_metres — structure-tensor integration radius; must be > 0.
//   contrast                  — coherence contrast C (> 0). Selectivity of the
//                               along-feature smoothing.
//   min_diffusivity           — floor eigenvalue (across-feature diffusivity) in
//                               (0, 1); small, e.g. 0.001.
//   n_threads                 — OpenMP thread count; 0 = all available.
//
// Throws
// ------
//   std::runtime_error    if metadata is invalid.
//   std::invalid_argument on out-of-range parameters.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu coherence_enhancing_diffusion(
    const lsdtt_raster_cpu& field,
    int                     n_iterations,
    double                  time_step,
    double                  derivative_radius_metres,
    double                  integration_radius_metres,
    double                  contrast,
    double                  min_diffusivity = 0.001,
    int                     n_threads       = 0
);

// -----------------------------------------------------------------------------
// oriented_ridge_result
//
//   thinned — the response after oriented non-maximum suppression: unchanged at
//             ridge crests, 0 where suppressed. A ~1-pixel-wide centreline map.
//   mask    — binary result of hysteresis linking the thinned response: 1 on
//             retained centreline pixels, 0 elsewhere.
// -----------------------------------------------------------------------------
struct oriented_ridge_result
{
    lsdtt_raster_cpu thinned;
    lsdtt_raster_cpu mask;
};

// -----------------------------------------------------------------------------
// oriented_nms_hysteresis
//
// Oriented non-maximum suppression followed by hysteresis linking — the
// Canny-style edge logic, but applied to a LINE response using an externally
// supplied orientation field (e.g. the structure-tensor feature orientation).
// Turns a thick, fragmented ridge response into thin, connected centrelines.
//
// Non-maximum suppression
// -----------------------
// At each pixel the response is compared with the response a unit step to each
// side ALONG the feature normal (perpendicular to orientation_deg), sampled by
// bilinear interpolation. The pixel is kept only if it is >= both neighbours;
// otherwise it is set to 0. This thins the ridge to its crest. Off-grid / nodata
// neighbours count as -infinity (the centre is kept), so coverage is not eroded
// at boundaries.
//
// Hysteresis
// ----------
// Seed pixels with thinned response >= high_threshold; then flood-fill (8-
// connected) through pixels with thinned response >= low_threshold. Only ridge
// chains that contain at least one strong seed are kept — isolated weak maxima
// are dropped, while genuine threads survive their faint stretches.
//
// Parameters
// ----------
//   response        — the ridge response to thin (e.g. a diffused divergence,
//                     a vesselness, or a brightness map). Higher = more ridge.
//   orientation_deg — feature-tangent azimuth in [0,180) deg clockwise from
//                     north, as produced by compute_structure_tensor_coherence.
//                     Must match `response` in size and CRS.
//   high_threshold  — strong-seed threshold on the thinned response.
//   low_threshold   — weak-link threshold; must be <= high_threshold.
//   n_threads       — OpenMP thread count; 0 = all available (NMS only; the
//                     hysteresis flood-fill is serial).
//
// Throws
// ------
//   std::runtime_error    if metadata is invalid or the two rasters differ in
//                         dimensions.
//   std::invalid_argument if low_threshold > high_threshold.
// -----------------------------------------------------------------------------
[[nodiscard]] oriented_ridge_result oriented_nms_hysteresis(
    const lsdtt_raster_cpu& response,
    const lsdtt_raster_cpu& orientation_deg,
    double                  high_threshold,
    double                  low_threshold,
    int                     n_threads = 0
);

// -----------------------------------------------------------------------------
// trim_valid_border
//
// Sets to nodata every pixel within border_pixels of the boundary of the valid
// domain — i.e. it erodes the validity mask by a square structuring element of
// radius border_pixels. The "boundary" includes BOTH the grid edge AND the
// edge of any irregular nodata footprint (out-of-grid neighbours are treated as
// invalid), so a raster with a rounded / clipped data area is trimmed inward
// from its real data edge, not just the rectangular frame.
//
// Why
// ---
// Polyfit-, flow-, and FFT-derived response rasters carry strong artefacts in a
// margin near their data boundary. Because those artefacts are often laterally
// continuous, they survive curvilinear filters such as grayscale_path_opening
// and masquerade as features. Trimming the margin before detection removes them.
//
// Implementation
// --------------
// Separable square erosion of the boolean validity mask: a horizontal pass
// (a pixel survives only if every column in [c-N, c+N] within the grid is
// valid) followed by a vertical pass over the result. Out-of-grid neighbours
// count as invalid, which trims the grid edge for free.
//
// Parameters
// ----------
//   raster        — input raster (any field)
//   border_pixels — N, margin width in pixels to remove. Must be >= 1.
//
// Returns
// -------
//   A copy of the input with the margin pixels set to nodata (data and mask).
//   Interior pixel values are unchanged.
//
// Throws
// ------
//   std::runtime_error    if metadata is invalid.
//   std::invalid_argument if border_pixels < 1.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu trim_valid_border(
    const lsdtt_raster_cpu& raster,
    int                     border_pixels
);

} // namespace lsdtt
