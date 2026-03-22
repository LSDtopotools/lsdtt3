// lsdtt_raster_synthetic_surfaces.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Synthetic raster generators for testing analysis routines.
//
// Current contents:
//   - Periodic wave surfaces built from seeded random cosine modes
//   - Exact global quadratic surfaces with analytical local coefficients
//
// Coefficient convention:
//   The returned analytical coefficients are the exact local quadratic
//   coefficients at each pixel centre in the form:
//
//     z(dx,dy) ≈ a*dx^2 + b*dy^2 + c*dx*dy + d*dx + e*dy + f
//
//   where dx and dy are local coordinate offsets from the pixel centre,
//   with dx positive east and dy positive north.
//
// Therefore:
//   a = 0.5 * z_xx
//   b = 0.5 * z_yy
//   c = z_xy
//   d = z_x
//   e = z_y
//   f = z
//
// Notes:
//   - The periodic surface is ideal for testing repeatability, sign conventions,
//     and convergence behaviour.
//   - The global quadratic surface is the strict exact oracle for
//     wavelet_quadratic_fit() away from invalid edges because the fitted model
//     matches the true surface exactly.
// -----------------------------------------------------------------------------

#pragma once

#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"

namespace lsdtt
{


// -----------------------------------------------------------------------------
// Stochastic synthetic raster generators
//
// These return only a DEM raster because there is no simple closed-form
// analytical coefficient field for the stochastic cases.
//
// Determinism:
//   All stochastic generators below use a hash-based seeded pseudo-random
//   mapping keyed by pixel and/or algorithm stage, so their output is
//   deterministic for a given seed and independent of thread scheduling.
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// generate_uniform_random_noise_surface
//
// Pixelwise independent uniform random noise in [min_value, max_value].
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_uniform_random_noise_surface(
    const lsdtt_raster_meta& meta,
    int    seed       = 1,
    double min_value  = -1.0,
    double max_value  =  1.0,
    int    n_threads  = 0
);

lsdtt_raster_cpu generate_uniform_random_noise_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed       = 1,
    double min_value  = -1.0,
    double max_value  =  1.0,
    int    n_threads  = 0
);

// -----------------------------------------------------------------------------
// generate_red_noise_surface
//
// Low-frequency enhanced noise obtained by repeated periodic binomial smoothing
// of deterministic white noise, followed by mean/std normalisation.
//
// target_mean / target_stddev define the output distribution after smoothing.
// blur_passes controls the degree of correlation.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_red_noise_surface(
    const lsdtt_raster_meta& meta,
    int    seed          = 1,
    double target_mean   = 0.0,
    double target_stddev = 1.0,
    int    blur_passes   = 6,
    int    n_threads     = 0
);

lsdtt_raster_cpu generate_red_noise_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed          = 1,
    double target_mean   = 0.0,
    double target_stddev = 1.0,
    int    blur_passes   = 6,
    int    n_threads     = 0
);

// -----------------------------------------------------------------------------
// generate_blue_noise_surface
//
// High-frequency enhanced noise formed by subtracting a low-pass filtered
// version of deterministic white noise from the white noise itself, followed
// by mean/std normalisation.
//
// Note:
//   This is a blue-ish spectral tilt generator, not a strict Poisson-disk
//   blue-noise sampler.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_blue_noise_surface(
    const lsdtt_raster_meta& meta,
    int    seed          = 1,
    double target_mean   = 0.0,
    double target_stddev = 1.0,
    int    blur_passes   = 6,
    int    n_threads     = 0
);

lsdtt_raster_cpu generate_blue_noise_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed          = 1,
    double target_mean   = 0.0,
    double target_stddev = 1.0,
    int    blur_passes   = 6,
    int    n_threads     = 0
);

// -----------------------------------------------------------------------------
// generate_diamond_square_surface
//
// Fractal terrain generated by the diamond-square algorithm on the smallest
// enclosing (2^k + 1) square grid, then cropped to the requested raster size.
//
// Parameters:
//   roughness      — multiplicative decay of perturbation amplitude per level
//   initial_scale  — initial perturbation amplitude
//   target_mean    — output mean after normalisation
//   target_stddev  — output standard deviation after normalisation
//
// Determinism:
//   Random perturbations are generated from a seed + (step,row,col,phase)
//   hash, so the result is independent of thread scheduling.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_diamond_square_surface(
    const lsdtt_raster_meta& meta,
    int    seed          = 1,
    double roughness     = 0.55,
    double initial_scale = 1.0,
    double target_mean   = 0.0,
    double target_stddev = 1.0,
    int    n_threads     = 0
);

lsdtt_raster_cpu generate_diamond_square_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed          = 1,
    double roughness     = 0.55,
    double initial_scale = 1.0,
    double target_mean   = 0.0,
    double target_stddev = 1.0,
    int    n_threads     = 0
);





// -----------------------------------------------------------------------------
// lsdtt_periodic_wave_component
//
// One periodic component:
//
//   z_n(x,y) = amplitude * cos( 2*pi*(mx*x/Lx + my*y/Ly) + phase )
//
// where:
//   mx, my  — integer mode numbers
//   Lx, Ly  — full raster lengths in x and y
// -----------------------------------------------------------------------------
struct lsdtt_periodic_wave_component
{
    int    mode_x    = 0;
    int    mode_y    = 0;
    double amplitude = 0.0;
    double phase     = 0.0;
};

// -----------------------------------------------------------------------------
// lsdtt_periodic_wave_surface
//
// dem              — synthetic DEM
// analytic_taylor  — exact local Taylor coefficients at each pixel centre
// components       — actual random wave components used
// offset           — constant offset added to the DEM
// -----------------------------------------------------------------------------
struct lsdtt_periodic_wave_surface
{
    lsdtt_raster_cpu                           dem;
    poly_fit_coeffs                            analytic_taylor;
    std::vector<lsdtt_periodic_wave_component> components;
    double                                     offset = 0.0;
};

// -----------------------------------------------------------------------------
// lsdtt_global_quadratic_coefficients
//
// Global coefficients of the synthetic surface expressed in centred raster
// coordinates:
//
//   z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
//
// where x and y are measured from the raster centre, with x positive east
// and y positive north.
// -----------------------------------------------------------------------------
struct lsdtt_global_quadratic_coefficients
{
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double d = 0.0;
    double e = 0.0;
    double f = 0.0;
};

// -----------------------------------------------------------------------------
// lsdtt_global_quadratic_surface
//
// dem              — synthetic DEM
// analytic_taylor  — exact local quadratic coefficients at each pixel centre
// global_coeffs    — coefficients of the underlying global quadratic surface
// -----------------------------------------------------------------------------
struct lsdtt_global_quadratic_surface
{
    lsdtt_raster_cpu                  dem;
    poly_fit_coeffs                   analytic_taylor;
    lsdtt_global_quadratic_coefficients global_coeffs;
};

// -----------------------------------------------------------------------------
// generate_periodic_wave_test_surface
//
// Generates a periodic smooth surface over the raster domain using a seeded
// pseudorandom set of cosine modes with random amplitudes and phases.
//
// Notes:
//   - Same seed => identical surface
//   - The returned analytic_taylor rasters are exact local Taylor
//     coefficients at each pixel centre
//   - n_threads controls only the raster fill loops; random generation is
//     serial and deterministic
// -----------------------------------------------------------------------------
lsdtt_periodic_wave_surface generate_periodic_wave_test_surface(
    const lsdtt_raster_meta& meta,
    int    n_components  = 8,
    int    seed          = 1,
    double amplitude_min = 1.0,
    double amplitude_max = 10.0,
    int    mode_min      = 1,
    int    mode_max      = 6,
    bool   random_sign   = true,
    int    n_threads     = 0
);

// Convenience overload
lsdtt_periodic_wave_surface generate_periodic_wave_test_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    n_components  = 8,
    int    seed          = 1,
    double amplitude_min = 1.0,
    double amplitude_max = 10.0,
    int    mode_min      = 1,
    int    mode_max      = 6,
    bool   random_sign   = true,
    int    n_threads     = 0
);

// -----------------------------------------------------------------------------
// generate_global_quadratic_test_surface
//
// Generates a globally quadratic surface:
//
//   z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
//
// in centred raster coordinates (x,y) measured from the raster centre.
//
// This is the exact oracle for wavelet_quadratic_fit() away from edges because
// the fitted model matches the true surface class exactly.
// -----------------------------------------------------------------------------
lsdtt_global_quadratic_surface generate_global_quadratic_test_surface(
    const lsdtt_raster_meta&               meta,
    const lsdtt_global_quadratic_coefficients& coeffs,
    int                                    n_threads = 0
);

// Convenience overload
lsdtt_global_quadratic_surface generate_global_quadratic_test_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    const lsdtt_global_quadratic_coefficients& coeffs,
    int    n_threads = 0
);

// -----------------------------------------------------------------------------
// generate_random_quadratic_test_surface
//
// Generates a globally quadratic surface using seeded random coefficients.
// Same seed => identical surface.
//
// Parameter ranges:
//   quadratic_abs_max — max absolute value of a, b, c
//   linear_abs_max    — max absolute value of d, e
//   constant_abs_max  — max absolute value of f
// -----------------------------------------------------------------------------
lsdtt_global_quadratic_surface generate_random_quadratic_test_surface(
    const lsdtt_raster_meta& meta,
    int    seed              = 1,
    double quadratic_abs_max = 1e-3,
    double linear_abs_max    = 1.0,
    double constant_abs_max  = 100.0,
    bool   random_sign       = true,
    int    n_threads         = 0
);

// Convenience overload
lsdtt_global_quadratic_surface generate_random_quadratic_test_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed              = 1,
    double quadratic_abs_max = 1e-3,
    double linear_abs_max    = 1.0,
    double constant_abs_max  = 100.0,
    bool   random_sign       = true,
    int    n_threads         = 0
);

} // namespace lsdtt
