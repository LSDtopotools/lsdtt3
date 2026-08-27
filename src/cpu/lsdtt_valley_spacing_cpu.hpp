// =============================================================================
// lsdtt_valley_spacing_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Characteristic valley spacing of a ridge-and-valley landscape.
//
// Three measures are returned together, deliberately, because each fails in a
// different way and the disagreement between them is itself the diagnostic:
//
//   1. SPECTRAL, Perron et al. (2008, JGR 113, F04016). Detrend along the flow
//      direction, take the power spectrum across it, and where several peaks
//      exist report the POWER-WEIGHTED MEAN of the peak frequencies rather than
//      the single strongest. The weighting is the point: valley spacing occupies
//      a band rather than a single value, so an argmax hops between domain
//      harmonics while the weighted mean moves smoothly.
//
//   2. CONVEX HALF-WIDTH. The mean run length of convex (hillslope) cells along
//      the measurement axis. Thousands of samples per landscape, so it is immune
//      to the harmonic quantisation that makes a spectral peak useless on a
//      domain holding only a few valleys. On a 200x200 test grid the spectral
//      peak could only return L/2, L/3 or L/4 while this measure resolved a clean
//      L_h ~ D^0.64 scaling over a 25x range in diffusivity.
//
//   3. SPECTRAL EXPONENT alpha, Bonetti et al. (2020, PNAS 117:1375) and
//      Hooshyar et al. (2021, GRL). The small-scale power-law slope of the same
//      spectrum, predicted to be alpha = 2m - 3.
//
// Alongside them the CHANNELIZATION INDEX of Bonetti et al. (2020):
//
//      CI = K l^(m+n) / (D^n U^(1-n))
//
// which for n = 1 is independent of uplift rate. CI_c ~ 32-37 marks the first
// channel instability, so a landscape below that should show no regular valleys
// at all and one far above it should show a cascade of them. Reporting CI beside
// the measured spacing is what makes a result comparable to the literature.
//
// SPACING IS A BAND, NOT A NUMBER. Repeat runs and neighbouring parameter values
// scatter within it. A single landscape gives one draw from that band; do not
// read a trend from two points.
// =============================================================================

#ifndef LSDTT_VALLEY_SPACING_CPU_HPP
#define LSDTT_VALLEY_SPACING_CPU_HPP

#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_valley_spacing_options
// -----------------------------------------------------------------------------
struct lsdtt_valley_spacing_options
{
    // Which axis the spacing is measured ALONG. Valleys run in the flow
    // direction, so for a landscape draining north-south the spacing is measured
    // east-west, i.e. along rows.
    enum class measurement_axis { along_rows, along_columns };
    measurement_axis axis = measurement_axis::along_rows;

    // Cells to drop at each end of the FLOW direction. Base level pins elevation
    // there, so those rows carry the boundary's imprint rather than the
    // landscape's own scale.
    int margin_cells = 20;

    // Parameters for the channelization index and the alpha prediction. Leave K,
    // D or U at zero to skip CI; leave m at zero to skip the alpha prediction.
    double K = 0.0;
    double D = 0.0;
    double U = 0.0;
    double m = 0.45;
    double n = 1.0;

    // Domain length used in CI. Zero means "the extent along the flow direction".
    double domain_length_m = 0.0;

    // Window for the alpha fit, as multiples of the MEASURED valley wavelength.
    // Fitting at a fixed fraction of Nyquist instead samples the hillslope
    // roll-off at high CI and returns a slope that means nothing, so the window
    // is set relative to the scale actually present.
    double alpha_fit_min_wavelength_frac = 0.08;
    double alpha_fit_max_wavelength_frac = 0.60;

    int n_threads = 0;
};


// -----------------------------------------------------------------------------
// lsdtt_valley_spacing_result
// -----------------------------------------------------------------------------
struct lsdtt_valley_spacing_result
{
    // --- the three measures ---
    double perron_wavelength_m   = 0.0;  // power-weighted mean of spectral peaks
    double peak_wavelength_m     = 0.0;  // single argmax, reported for contrast
    double convex_halfwidth_m    = 0.0;  // mean convex run
    double convex_wavelength_m   = 0.0;  // 2 * half-width
    double spectral_exponent     = 0.0;  // alpha
    double spectral_exponent_predicted = 0.0;   // 2m - 3

    // --- context ---
    double channelization_index  = 0.0;  // CI; 0 when K/D/U were not supplied
    double domain_length_m       = 0.0;
    double valleys_in_domain     = 0.0;  // domain_length / convex_wavelength

    int    n_spectral_peaks      = 0;
    long   n_convex_runs         = 0;

    // Is the domain big enough for the SPECTRAL measure to mean anything? With
    // only a few valleys the peak can only land on domain harmonics.
    bool   spectrum_under_sampled = false;

    // The averaged spectrum, so a chart can show where the peaks were found
    // rather than asking the reader to trust a scalar.
    std::vector<double> frequency_cycles_per_m;
    std::vector<double> power;

    std::string warning;   // empty unless something needs saying
};


// -----------------------------------------------------------------------------
// compute_valley_spacing
//
// `dem` should be a filled, equilibrated landscape. Nodata cells are excluded
// from the convex-run measure; a row containing any nodata is excluded from the
// spectrum, because a gap there is a discontinuity that leaks across all
// frequencies.
// -----------------------------------------------------------------------------
lsdtt_valley_spacing_result compute_valley_spacing(
    const lsdtt_raster_cpu&             dem,
    const lsdtt_valley_spacing_options& opts);


// -----------------------------------------------------------------------------
// Critical channelization index of Bonetti et al. (2020): CI_c ~ 32-37.
// -----------------------------------------------------------------------------
constexpr double k_channelization_index_critical_low  = 32.0;
constexpr double k_channelization_index_critical_high = 37.0;

// The Bonetti channelization index, exposed on its own because it is a
// parameter-only quantity: it can be evaluated before a run, not just after one.
[[nodiscard]] double channelization_index(
    double K, double D, double U, double domain_length_m, double m, double n);

} // namespace lsdtt

#endif // LSDTT_VALLEY_SPACING_CPU_HPP
