// lsdtt_raster_spectral_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Spectral raster analysis — FFTW-based CPU free functions.
//
// Current contents:
//   wavelet_quadratic_fit
//   various spectral filters
//
// The fitted quadratic surface is:
//
//   z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
//
// Coordinate convention:
//   x  — easting offset from the pixel centre, positive east
//   y  — northing offset from the pixel centre, positive north
//
// Notes:
//   - Gaussian kernel only.
//   - Uses FFTW double API (fftw_*) internally. Not threaded FFTW.
//   - If lsdtt3 is built without FFTW, wavelet_quadratic_fit() throws at
//     runtime with a clear message.
//   - Output is conservative: pixels whose support overlaps nodata or the
//     raster edge are set to nodata.
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"

#include <vector>

namespace lsdtt
{

// -----------------------------------------------------------------------------
// poly_fit_coeffs
//
// Six coefficient rasters from the local quadratic fit.
// All share the same metadata as the input raster.
// -----------------------------------------------------------------------------
struct poly_fit_coeffs
{
    lsdtt_raster_cpu a;  // x^2 coefficient
    lsdtt_raster_cpu b;  // y^2 coefficient
    lsdtt_raster_cpu c;  // x*y coefficient
    lsdtt_raster_cpu d;  // x coefficient
    lsdtt_raster_cpu e;  // y coefficient
    lsdtt_raster_cpu f;  // constant term
};

// -----------------------------------------------------------------------------
// wavelet_quadratic_fit
//
// Parameters:
//   raster    — input DEM; must be projected in metres
//   lambda    — target wavelength in metres; must be >= cell_size
//   n_threads — OpenMP threads for post-FFT coefficient solve (0 = all)
//
// Returns:
//   poly_fit_coeffs containing coefficient rasters a..f.
//
// Throws:
//   std::runtime_error on invalid metadata, invalid lambda, missing FFTW,
//   FFTW allocation/plan failures, or if no valid pixels remain after
//   conservative support masking.
// -----------------------------------------------------------------------------
poly_fit_coeffs wavelet_quadratic_fit(
    const lsdtt_raster_cpu& raster,
    double lambda    = 100.0,
    int    n_threads = 0
);

// =============================================================================
// 2D FFT-based spectral filters
//
// All filters operate on a full DEM raster through a shared pipeline:
//   1. Detrend  — fit and subtract an OLS planar trend (all valid pixels)
//   2. Window   — 2D elliptical taper to suppress edge leakage (optional)
//   3. Pad      — zero-pad to next power-of-2 dimensions
//   4. Forward FFT (real-to-complex, r2c)
//   5. Apply filter kernel in frequency domain
//   6. Inverse FFT (complex-to-real, c2r), normalise by Ly * Lx
//   7. Restore planar trend
//   8. Apply original nodata mask
//
// Frequency convention:
//   Radial frequency  f = sqrt(fx^2 + fy^2)   [cycles / metre]
//   Wavelength        lambda = 1 / f           [metres]
//   Nyquist           f_ny = 1 / (2 * cell_size)
//
// Filter shape — Gaussian rolloff, -3 dB at the specified cutoff:
//
//   Lowpass   H(f)  = exp( -ln(2) * (f / f_c)^2 )
//   Highpass  H(f)  = 1 - exp( -ln(2) * (f / f_c)^2 )
//   Bandpass  H(f)  = H_highpass(f, f_low) * H_lowpass(f, f_high)
//
// where f_c = 1 / lambda_cutoff_metres.
//
// Edge nodata handling
// --------------------
// DEMs often have nodata around the edges — either triangular slices from
// reprojection out of a geographic CRS, or irregular boundaries from
// basin-clipping.  The detrend step sets these pixels to the trend plane level
// and the windowing taper then drives them to zero before zero-padding.  For
// mild edge fringing this is adequate.  For heavily clipped DEMs the transition
// from data to zero can still introduce low-frequency artefacts.
//
// Set clip_to_valid_rectangle = true to first crop the DEM to the largest
// interior rectangle that contains no nodata (via clip_to_maximum_rectangle),
// eliminating edge nodata before filtering.  The returned raster will be
// smaller than the input; the caller is responsible for mosaicking or
// interpreting the cropped extent.
//
// For spectral analysis (Perron et al. 2008) see lsdtt_spectral_analysis_result
// which is in this same file (to be added alongside full_spectral_analysis).
//
// If lsdtt3 is built without FFTW, all filter functions throw at runtime with
// a clear message.
// =============================================================================

// -----------------------------------------------------------------------------
// spectral_window_type
//
// none    — unity window; maximum frequency resolution but prone to edge leakage
// hann    — 2D elliptical raised-cosine window; recommended for most analyses
// hamming — 2D elliptical Hamming window; slightly different sidelobe profile
// -----------------------------------------------------------------------------
enum class spectral_window_type { none, hann, hamming };

// -----------------------------------------------------------------------------
// lsdtt_spectral_binned_result
//
// Log-binned radial spectrum product used by full_spectral_analysis().
// beta_by_bin and r_sq_by_bin are the rollover-fit traces evaluated from
// repeated suffix regressions on the binned spectrum.
//
// Units:
//   frequency      — cycles/m  (= 1/m)
//   wavelength     — metres    (= 1/frequency)
//   psd            — m^2 / (cycles/m)
//   background_psd — m^2 / (cycles/m)
//   normalised_psd — dimensionless  (psd / background_psd)
//   beta_by_bin    — dimensionless spectral slope from rollover suffix regression
//   r_sq_by_bin    — R^2 of that regression (dimensionless, 0–1)
// -----------------------------------------------------------------------------
struct lsdtt_spectral_binned_result
{
    std::vector<double> frequency;       // cycles/m
    std::vector<double> wavelength;      // m
    std::vector<double> psd;             // m^2/(cycles/m)
    std::vector<double> background_psd;  // m^2/(cycles/m)
    std::vector<double> normalised_psd;  // dimensionless
    std::vector<double> beta_by_bin;     // spectral slope
    std::vector<double> r_sq_by_bin;     // R^2
};

// -----------------------------------------------------------------------------
// lsdtt_spectral_analysis_options
//
// Options controlling the Perron-style full spectral analysis workflow.
// -----------------------------------------------------------------------------
struct lsdtt_spectral_analysis_options
{
    double               log_bin_width           = 0.1;
    int                  n_background_iterations = 100;
    spectral_window_type window                  = spectral_window_type::hann;
    bool                 clip_to_valid_rectangle = false;
    int                  n_threads               = 0;
    // Base seed for Monte Carlo background iterations.
    // Iteration i uses seed = background_seed + i.
    // Negative (default) uses the built-in sequence (base 1000), which is
    // always reproducible. Set an explicit non-negative value to pin a
    // specific realisation and document it alongside the results.
    int                  background_seed         = -1;
};

// -----------------------------------------------------------------------------
// lsdtt_spectral_analysis_result
//
// Full spectral-analysis product.
//
// 2D arrays (psd_2d, background_psd_2d):
//   - Stored in fftshift order: zero frequency (DC) at row padded_rows/2,
//     col padded_cols/2 (integer division).
//   - Row index r has y-frequency (r - padded_rows/2) * dfy  [cycles/m].
//   - Col index c has x-frequency (c - padded_cols/2) * dfx  [cycles/m].
//   - When written as a GeoTIFF via lsdtt_data_loader::write_full(), set the
//     origin so that the DC pixel centre maps to coordinate (0, 0) in
//     cycles/m. See lsdtt-spectral-analysis CLI for the derivation.
//
// Radial vectors (radial_frequency, radial_psd, etc.):
//   - One entry per non-redundant (row, col) pair in the left half of the
//     shifted spectrum. Radial frequency is in cycles/m.
//
// Units summary:
//   dfx, dfy, nyquist_frequency  — cycles/m
//   radial_frequency             — cycles/m
//   radial_wavelength            — m
//   radial_psd, background_psd   — m^2/(cycles/m)
//   normalised_psd               — dimensionless
//   rollover_frequency           — cycles/m
//   rollover_beta, sub_rollover_beta — dimensionless spectral slope
// -----------------------------------------------------------------------------
struct lsdtt_spectral_analysis_result
{
    int padded_rows = 0;
    int padded_cols = 0;

    double dfx = 0.0;               // x-frequency spacing, cycles/m
    double dfy = 0.0;               // y-frequency spacing, cycles/m
    double nyquist_frequency = 0.0; // cycles/m

    double window_sum_squared = 0.0;
    double detrended_variance = 0.0;

    lsdtt_array2d_f psd_2d;
    lsdtt_array2d_f background_psd_2d;

    std::vector<double> radial_frequency;
    std::vector<double> radial_wavelength;
    std::vector<double> radial_psd;

    std::vector<double> background_psd;
    std::vector<double> normalised_psd;

    std::vector<double> ci95;
    std::vector<double> norm_ci95;
    std::vector<double> norm_ci99;

    double rollover_frequency = 0.0;
    double rollover_beta = 0.0;
    double sub_rollover_beta = 0.0;

    lsdtt_spectral_binned_result binned;
};

// -----------------------------------------------------------------------------
// full_spectral_analysis
//
// Performs the legacy-style full spectral analysis workflow:
//
//   1. compute observed periodogram
//   2. estimate rollover frequency / beta from log-binned radial PSD
//   3. generate power-law synthetic reference DEMs
//   4. average their spectra into a background spectrum
//   5. normalise the observed spectrum by the background spectrum
//
// The dem must be in a projected coordinate system with metres as the linear
// unit (enforced via lsdtt_raster_meta::check_projected_metres()).
//
// All frequencies in the result are in cycles/m (= 1/m).
// All wavelengths are in metres.
// A frequency f is significant at the 95% confidence level when
//   normalised_psd > 2.9955  (= 0.5 * chi^2_0.95 with 2 DOF).
//
// Returns a full result struct and performs no file writing.
// -----------------------------------------------------------------------------
lsdtt_spectral_analysis_result full_spectral_analysis(
    const lsdtt_raster_cpu&              dem,
    const lsdtt_spectral_analysis_options& opts = {}
);



// -----------------------------------------------------------------------------
// fftw2D_lowpass_filter
//
// Returns the long-wavelength component of the DEM.
// Wavelengths >> lambda_cutoff_metres pass through; shorter wavelengths are
// suppressed with a Gaussian rolloff.  -3 dB point at lambda_cutoff_metres.
//
// Typical use: recover the gently-sloping piedmont / alluvial-fan surface
// from a DEM that also contains incised channels and short-wavelength relief.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fftw2D_lowpass_filter(
    const lsdtt_raster_cpu& dem,
    double               lambda_cutoff_metres,
    spectral_window_type window                = spectral_window_type::hann,
    bool                 clip_to_valid_rectangle = false,
    int                  n_threads              = 0
);

// -----------------------------------------------------------------------------
// fftw2D_highpass_filter
//
// Returns the short-wavelength residual (DEM minus lowpass component).
// -3 dB point at lambda_cutoff_metres.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fftw2D_highpass_filter(
    const lsdtt_raster_cpu& dem,
    double               lambda_cutoff_metres,
    spectral_window_type window                = spectral_window_type::hann,
    bool                 clip_to_valid_rectangle = false,
    int                  n_threads              = 0
);

// -----------------------------------------------------------------------------
// fftw2D_bandpass_filter
//
// Passes wavelengths between lambda_short_metres and lambda_long_metres.
// lambda_short_metres must be < lambda_long_metres.
// Gaussian rolloff at both edges.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fftw2D_bandpass_filter(
    const lsdtt_raster_cpu& dem,
    double               lambda_short_metres,
    double               lambda_long_metres,
    spectral_window_type window                = spectral_window_type::hann,
    bool                 clip_to_valid_rectangle = false,
    int                  n_threads              = 0
);

// -----------------------------------------------------------------------------
// fftw2D_wiener_filter
//
// Adaptive noise-removal filter following Milodowski / Perron et al.:
//   - Fits a power law S(f) = c * f^m to the spectrum in the mid-frequency
//     band (signal_wavelength_long_m to signal_wavelength_short_m).
//   - Estimates white noise amplitude N from the near-Nyquist spectrum.
//   - Applies Wiener weighting: phi(f) = S(f) / (S(f) + N).
//
// At low frequencies (signal dominates) phi -> 1; at high frequencies (noise
// dominates) phi -> 0.  The result is a DEM with noise attenuated while
// preserving the dominant topographic signal.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fftw2D_wiener_filter(
    const lsdtt_raster_cpu& dem,
    double               signal_wavelength_long_m  = 1000.0,
    double               signal_wavelength_short_m = 100.0,
    spectral_window_type window                    = spectral_window_type::hann,
    bool                 clip_to_valid_rectangle   = false,
    int                  n_threads                 = 0
);

// -----------------------------------------------------------------------------
// fftw2D_masked_surface_recovery
//
// Recovers a smooth long-wavelength background surface from only the pixels
// where mask > 0.5, by:
//   1. OLS-detrend the DEM using all valid pixels
//   2. Build a masked field: detrended elevation where mask > 0.5, 0 elsewhere
//   3. Apply 2D window taper and zero-pad to next power-of-2 dimensions
//   4. Forward FFT → lowpass filter at lambda_cutoff_metres → inverse FFT
//   5. Restore the planar trend
//
// The output is valid at every pixel where the input DEM is valid, including
// pixels not in the mask.  The lowpass filter interpolates a smooth surface
// across the masked-out gap.
//
// Intended use:
//   Recover the smooth piedmont/alluvial-fan surface from a DEM where mountain
//   pixels have been excluded from the mask.  The gap in the masked field is
//   filled by the spectral interpolation rather than the actual terrain, so the
//   recovered surface continues smoothly under the mountain region.
//
// Note:
//   Single-pass implementation.  An iterative Gerchberg-Papoulis refinement
//   can be added later as a separate function.
//
// Throws:
//   std::runtime_error if FFTW is not available, lambda_cutoff_metres <= 0,
//   or the mask contains no valid pixels.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fftw2D_masked_surface_recovery(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask,
    double               lambda_cutoff_metres,
    spectral_window_type window    = spectral_window_type::hann,
    int                  n_threads = 0
);


// -----------------------------------------------------------------------------
// directional_periodicity_result
//
// Per-pixel outputs of the local directional-periodicity (regular-spacing)
// detector. Each pixel inherits the result of the most confident overlapping
// analysis window covering it.
//
//   strength       — peak prominence of the dominant spectral peak within the
//                    searched spacing band: peak power / mean band power.
//                    ~1 = no periodicity (peak no stronger than background);
//                    large = a sharp, regular spacing. The score to threshold
//                    when deciding "is this a regularly-spaced field".
//   spacing_metres — dominant repeat spacing (1 / peak frequency), in metres.
//                    THE classifier: ditches, forest tracks, and wave dams each
//                    occupy a distinctive spacing band.
//   orientation_deg— azimuth ALONG which the parallel features run (perpendicular
//                    to the spacing direction), in [0,180) deg clockwise from
//                    north.
// -----------------------------------------------------------------------------
struct directional_periodicity_result
{
    lsdtt_raster_cpu strength;
    lsdtt_raster_cpu spacing_metres;
    lsdtt_raster_cpu orientation_deg;
};

// -----------------------------------------------------------------------------
// detect_directional_periodicity
//
// Finds families of regularly-spaced parallel features (ditches, forest tracks,
// "wave dams", cultivation ridges) by their characteristic spacing. A set of
// equally spaced parallel lines produces a sharp off-DC peak in the 2D power
// spectrum at frequency 1/spacing, oriented perpendicular to the lines. The
// detector slides an analysis window over the raster, and in each window:
//
//   1. subtract the window mean and apply a Hann taper (reduce spectral leakage);
//   2. take the 2D FFT and form the power spectrum;
//   3. within the frequency annulus corresponding to spacings in
//      [min_spacing_metres, max_spacing_metres], find the strongest peak;
//   4. record peak prominence (peak / mean band power), spacing (1/peak freq),
//      and the perpendicular line orientation.
//
// Windows overlap (stride = step_pixels); each pixel keeps the result of the
// window in which it scored the highest prominence.
//
// Parameters
// ----------
//   field               — input raster (e.g. a filled DEM such as sd_fill.tif,
//                         or a curvature/divergence response).
//   window_size_pixels  — analysis window side length. Must comfortably span
//                         several periods of the largest spacing of interest
//                         (a few * max_spacing_metres / cell_size). >= 8.
//   step_pixels         — window stride; smaller = smoother, slower. Typically
//                         window_size_pixels / 2. >= 1.
//   min_spacing_metres  — shortest spacing to search for (> 2*cell_size).
//   max_spacing_metres  — longest spacing to search for (< window extent).
//   min_valid_fraction  — minimum fraction of valid pixels for a window to be
//                         analysed (default 0.6); sparser windows are skipped.
//   n_threads           — OpenMP thread count; 0 = all available.
//
// Throws
// ------
//   std::runtime_error    if FFTW is not available, metadata is invalid, or the
//                         spacing band is empty for the given window.
//   std::invalid_argument on out-of-range parameters.
// -----------------------------------------------------------------------------
directional_periodicity_result detect_directional_periodicity(
    const lsdtt_raster_cpu& field,
    int                     window_size_pixels,
    int                     step_pixels,
    double                  min_spacing_metres,
    double                  max_spacing_metres,
    double                  min_valid_fraction = 0.6,
    int                     n_threads          = 0
);

} // namespace lsdtt
