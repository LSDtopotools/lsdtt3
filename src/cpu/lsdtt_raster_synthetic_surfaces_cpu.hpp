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
//     z(dx,dy) ~= a*dx^2 + b*dy^2 + c*dx*dy + d*dx + e*dy + f
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

#include <cstdint>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"

namespace lsdtt
{


// -----------------------------------------------------------------------------
// superimpose_ns_parabola
//
// Adds a N-S parabolic tilt to an existing raster in-place equivalent fashion
// (returns a new raster).
//
// The parabola peaks at the centre row with value `peak_elevation` and is zero
// at row 0 and row nrows-1 (the north and south boundaries).  Adding this to a
// random surface ensures a general slope toward both boundaries so that flow
// drains outward rather than pooling in internal basins.
//
//   z_parabola(row) = peak * (1 - ((row - c) / c)^2)
//   where c = (nrows - 1) / 2  (centre row, fractional for even nrows)
// -----------------------------------------------------------------------------
lsdtt_raster_cpu superimpose_ns_parabola(
    const lsdtt_raster_cpu& dem,
    float peak_elevation
);

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
//   roughness      - multiplicative decay of perturbation amplitude per level
//   initial_scale  - initial perturbation amplitude
//   target_mean    - output mean after normalisation
//   target_stddev  - output standard deviation after normalisation
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
// generate_power_law_fractal_surface_fft
//
// FFT-based synthetic fractal DEM generator for spectral background modelling.
//
// Intended use:
//   - generate reference DEMs with an approximate power-law spectrum
//   - used by full spectral analysis background-spectrum construction
//
// Parameters:
//   beta          - target PSD slope exponent (> 0)
//   target_mean   - mean of the output raster after normalisation
//   target_stddev - standard deviation of the output raster after normalisation
//
// Notes:
//   - Same seed => identical surface
//   - Implementation will use FFTW when added
//   - For now this is a declared public API used by spectral analysis
// -----------------------------------------------------------------------------
lsdtt_raster_cpu generate_power_law_fractal_surface_fft(
    const lsdtt_raster_meta& meta,
    int    seed          = 1,
    double beta          = 2.0,
    double target_mean   = 0.0,
    double target_stddev = 1.0,
    int    n_threads     = 0
);

lsdtt_raster_cpu generate_power_law_fractal_surface_fft(
    int    n_rows,
    int    n_cols,
    double cell_size,
    int    seed          = 1,
    double beta          = 2.0,
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
//   mx, my  - integer mode numbers
//   Lx, Ly  - full raster lengths in x and y
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
// dem              - synthetic DEM
// analytic_taylor  - exact local Taylor coefficients at each pixel centre
// components       - actual random wave components used
// offset           - constant offset added to the DEM
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
// dem              - synthetic DEM
// analytic_taylor  - exact local quadratic coefficients at each pixel centre
// global_coeffs    - coefficients of the underlying global quadratic surface
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
//   quadratic_abs_max - max absolute value of a, b, c
//   linear_abs_max    - max absolute value of d, e
//   constant_abs_max  - max absolute value of f
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

// -----------------------------------------------------------------------------
// lsdtt_range_edge
//
// Which DEM edge a half-space or arc mountain range is associated with.
// -----------------------------------------------------------------------------
enum class lsdtt_range_edge { north, south, east, west };

// -----------------------------------------------------------------------------
// lsdtt_mountain_range_surface
//
// dem        - synthetic elevation raster
// front_mask - ground-truth zone labels: 1.0 = mountain, 0.0 = piedmont
//
// All four mountain-range generators return this type.
// The front is the boundary of the mask and is what mountain-front detection
// algorithms should recover.
// -----------------------------------------------------------------------------
struct lsdtt_mountain_range_surface
{
    lsdtt_raster_cpu dem;
    lsdtt_raster_cpu front_mask;
};

// -----------------------------------------------------------------------------
// generate_linear_range_surface
//
// A Gaussian ridge running the full length of the DEM interior, with symmetrical
// piedmont zones on both flanks.  The ridge position and half-width are expressed
// as fractions of the cross-axis extent so the geometry is resolution-independent.
//
// ridge_runs_north_south
//   true  - ridge runs parallel to the y-axis (N-S); cross-axis is E-W (x)
//   false - ridge runs parallel to the x-axis (E-W); cross-axis is N-S (y)
//
// ridge_position_fraction
//   Position of the ridge crest along the cross-axis (0 = west/south edge,
//   1 = east/north edge).  Default 0.5 places the ridge at the centre.
//
// ridge_half_width_fraction
//   Half-width of the mountain zone as a fraction of the cross-axis extent.
//   Pixels within this distance of the ridge are labelled mountain.
//   Must be in (0, 0.5].
// -----------------------------------------------------------------------------
lsdtt_mountain_range_surface generate_linear_range_surface(
    const lsdtt_raster_meta& meta,
    double ridge_elevation            = 2000.0,
    double piedmont_elevation         = 0.0,
    double ridge_position_fraction    = 0.5,
    double ridge_half_width_fraction  = 0.15,
    double mountain_noise_amplitude   = 200.0,
    double piedmont_noise_amplitude   = 20.0,
    bool   ridge_runs_north_south     = true,
    double noise_roughness            = 0.55,
    int    noise_seed                 = 1,
    int    n_threads                  = 0
);

lsdtt_mountain_range_surface generate_linear_range_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    double ridge_elevation            = 2000.0,
    double piedmont_elevation         = 0.0,
    double ridge_position_fraction    = 0.5,
    double ridge_half_width_fraction  = 0.15,
    double mountain_noise_amplitude   = 200.0,
    double piedmont_noise_amplitude   = 20.0,
    bool   ridge_runs_north_south     = true,
    double noise_roughness            = 0.55,
    int    noise_seed                 = 1,
    int    n_threads                  = 0
);

// -----------------------------------------------------------------------------
// generate_halfspace_range_surface
//
// The ridge occupies one DEM edge entirely - analogous to the south face of the
// Himalayas draining to the Ganges plains: the high terrain runs off the edge
// and only one piedmont side is represented.
//
// edge
//   Which DEM edge the mountains are against.
//
// front_fraction
//   Distance of the mountain front from that edge, expressed as a fraction of
//   the perpendicular DEM extent (e.g. 0.3 places the front 30 % of the way
//   across the DEM).  The Gaussian decay places the exp(-0.5) contour exactly
//   at this depth.  Must be in (0, 1).
// -----------------------------------------------------------------------------
lsdtt_mountain_range_surface generate_halfspace_range_surface(
    const lsdtt_raster_meta& meta,
    lsdtt_range_edge edge                 = lsdtt_range_edge::north,
    double ridge_elevation                = 2000.0,
    double piedmont_elevation             = 0.0,
    double front_fraction                 = 0.3,
    double mountain_noise_amplitude       = 200.0,
    double piedmont_noise_amplitude       = 20.0,
    double noise_roughness                = 0.55,
    int    noise_seed                     = 1,
    int    n_threads                      = 0
);

lsdtt_mountain_range_surface generate_halfspace_range_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    lsdtt_range_edge edge                 = lsdtt_range_edge::north,
    double ridge_elevation                = 2000.0,
    double piedmont_elevation             = 0.0,
    double front_fraction                 = 0.3,
    double mountain_noise_amplitude       = 200.0,
    double piedmont_noise_amplitude       = 20.0,
    double noise_roughness                = 0.55,
    int    noise_seed                     = 1,
    int    n_threads                      = 0
);

// -----------------------------------------------------------------------------
// generate_circular_range_surface
//
// An elliptical mountain dome in the DEM interior, with a continuous piedmont
// surrounding it.  Useful for testing detection of circular or basin-and-range
// type geometries.
//
// centre_{x,y}_fraction
//   Position of the dome centre as fractions of Lx and Ly (0 = west/south edge,
//   1 = east/north edge).
//
// radius_{x,y}_fraction
//   Semi-axis lengths of the ellipse as fractions of Lx and Ly respectively.
//   Pixels within the unit ellipse (normalised radius < 1) are labelled mountain.
// -----------------------------------------------------------------------------
lsdtt_mountain_range_surface generate_circular_range_surface(
    const lsdtt_raster_meta& meta,
    double ridge_elevation            = 2000.0,
    double piedmont_elevation         = 0.0,
    double centre_x_fraction         = 0.5,
    double centre_y_fraction         = 0.5,
    double radius_x_fraction         = 0.25,
    double radius_y_fraction         = 0.25,
    double mountain_noise_amplitude   = 200.0,
    double piedmont_noise_amplitude   = 20.0,
    double noise_roughness            = 0.55,
    int    noise_seed                 = 1,
    int    n_threads                  = 0
);

lsdtt_mountain_range_surface generate_circular_range_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    double ridge_elevation            = 2000.0,
    double piedmont_elevation         = 0.0,
    double centre_x_fraction         = 0.5,
    double centre_y_fraction         = 0.5,
    double radius_x_fraction         = 0.25,
    double radius_y_fraction         = 0.25,
    double mountain_noise_amplitude   = 200.0,
    double piedmont_noise_amplitude   = 20.0,
    double noise_roughness            = 0.55,
    int    noise_seed                 = 1,
    int    n_threads                  = 0
);

// -----------------------------------------------------------------------------
// generate_arc_range_surface
//
// A mountain range occupying one DEM edge whose front curves parabolically into
// the piedmont - deepest at the arc centre, shallowing toward the sides.  This
// approximates a concave range front as seen from the foreland (e.g. an eroded
// dome draining outward, or a range segment with a salient).
//
// edge
//   Which DEM edge the mountains are against.
//
// front_fraction
//   Maximum depth of the arc front from the edge (at the arc centre), expressed
//   as a fraction of the perpendicular DEM extent.
//
// arc_centre_fraction
//   Position of the arc centre along the edge (0..1).  Default 0.5.
//
// arc_half_width_fraction
//   Half-width of the arc as a fraction of the parallel DEM extent.  Controls
//   how far across the DEM the range extends.  Setting this to >= 1.0 makes the
//   mountains span the full DEM width with a gradual taper at the ends.
//   Must be > 0.
//
// The ground-truth front is the parabola
//   depth(u) = front_fraction * L_perp * max(0, 1 - u^2)
// where u = (parallel_coord - arc_centre) / arc_half_width.
// -----------------------------------------------------------------------------
lsdtt_mountain_range_surface generate_arc_range_surface(
    const lsdtt_raster_meta& meta,
    lsdtt_range_edge edge                   = lsdtt_range_edge::north,
    double ridge_elevation                  = 2000.0,
    double piedmont_elevation               = 0.0,
    double front_fraction                   = 0.3,
    double arc_centre_fraction              = 0.5,
    double arc_half_width_fraction          = 1.0,
    double mountain_noise_amplitude         = 200.0,
    double piedmont_noise_amplitude         = 20.0,
    double noise_roughness                  = 0.55,
    int    noise_seed                       = 1,
    int    n_threads                        = 0
);

lsdtt_mountain_range_surface generate_arc_range_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    lsdtt_range_edge edge                   = lsdtt_range_edge::north,
    double ridge_elevation                  = 2000.0,
    double piedmont_elevation               = 0.0,
    double front_fraction                   = 0.3,
    double arc_centre_fraction              = 0.5,
    double arc_half_width_fraction          = 1.0,
    double mountain_noise_amplitude         = 200.0,
    double piedmont_noise_amplitude         = 20.0,
    double noise_roughness                  = 0.55,
    int    noise_seed                       = 1,
    int    n_threads                        = 0
);

// -----------------------------------------------------------------------------
// Hillslope-strip test-domain generator
//
// Produces the 1-D hillslope-strip domains used to validate the transient
// Roering hillslope solver (see docs/testing_and_benchmarking.md).  The strip varies along a single cross-axis; the orthogonal axis
// is uniform, so every row (or column) is an identical 1-D hillslope profile.
//
// Channels are geometry-defined Dirichlet sinks held at baselevel_elevation
// either at both cross-axis ends with the divide at the centre (the
// symmetric test-1 parabola case) or at a single end with a no-flux divide at
// the opposite end.  The returned channel_mask is 1 on channel cells and 0
// elsewhere, exactly the is_channel convention the solver consumes.
// -----------------------------------------------------------------------------

// Where the Dirichlet channels sit along the cross-axis.
//   both_ends       - channels at both ends, no-flux divide at the centre (x=0).
//   single_end_low  - channel at the low-index end, no-flux divide at the high end.
//   single_end_high - channel at the high-index end, no-flux divide at the low end.
enum class lsdtt_strip_channel_config
{
    both_ends,
    single_end_low,
    single_end_high
};

// Initial (non-channel) elevation profile.
//   flat              - every hillslope cell at baselevel_elevation.
//   analytic_parabola - the steady-state profile (linear parabola, or the
//                       nonlinear Roering 2007 profile when nonlinear_profile
//                       is true); a discrete fixed point of the steady solver
//                       when half_width is auto-derived (half_width <= 0).
//   raised_plateau    - hillslope cells at baselevel_elevation + init_amplitude
//                       (a raised block to relax: test-2 transient relaxation).
//   perturbed         - flat baselevel plus seeded uniform noise of amplitude
//                       init_amplitude (a noisy start: test-2 transient).
//   linear_ramp       - a triangular ridge linear from a crest of init_amplitude
//                       above baselevel at the divide down to baselevel at the
//                       channel(s).
enum class lsdtt_strip_init_profile
{
    flat,
    analytic_parabola,
    raised_plateau,
    perturbed,
    linear_ramp
};

// -----------------------------------------------------------------------------
// lsdtt_hillslope_strip_surface
//
// dem          - synthetic elevation raster (the hillslope strip)
// channel_mask - 1.0 on geometry-defined channel cells, 0.0 elsewhere
//
// Analogous to lsdtt_mountain_range_surface{dem, front_mask}.
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_strip_surface
{
    lsdtt_raster_cpu dem;
    lsdtt_raster_cpu channel_mask;
};

// -----------------------------------------------------------------------------
// generate_hillslope_strip_surface
//
// cross_axis_is_x
//   true  - elevation varies along columns (E-W cross-axis x); channels are
//           full N-S lines at the E/W end column(s).
//   false - elevation varies along rows (N-S cross-axis); channels are full
//           E-W lines at the N/S end row(s).
//
// baselevel_elevation
//   Elevation of the channel cells (the Dirichlet value, m).
//
// half_width
//   The hillslope half-width L (m), i.e. the divide-to-channel distance used by
//   the analytic profiles.  If <= 0 it is auto-derived from the grid geometry so
//   the analytic profile is a discrete fixed point (recommended).  If > 0 it
//   overrides the L used in the profile formulas while channels remain at the
//   extreme cross-axis cells; a value that does not match the grid extent will
//   make analytic_parabola only approximately a fixed point.
//
// uplift_U, diffusivity_D, critical_slope_Sc
//   Physical scalars used only to build the analytic init profile
//   (analytic_parabola).  Sc is used only when nonlinear_profile is true.
//
// nonlinear_profile
//   false - analytic_parabola uses the linear steady parabola
//           z(x) = z_ch + (U/2D)(L^2 - x^2).
//   true  - analytic_parabola uses the nonlinear Roering (2007) steady profile
//           obtained by integrating D*S/(1-(S/Sc)^2) = U*x.
//
// init_amplitude
//   Height scale for raised_plateau (plateau height above baselevel), perturbed
//   (noise amplitude), and linear_ramp (crest height above baselevel).
//
// perturb_seed
//   Seed for the perturbed init profile (deterministic).
// -----------------------------------------------------------------------------
lsdtt_hillslope_strip_surface generate_hillslope_strip_surface(
    const lsdtt_raster_meta&   meta,
    lsdtt_strip_channel_config channel_config      = lsdtt_strip_channel_config::both_ends,
    lsdtt_strip_init_profile   init_profile        = lsdtt_strip_init_profile::flat,
    bool                       cross_axis_is_x     = true,
    double                     baselevel_elevation = 0.0,
    double                     half_width          = 0.0,
    double                     uplift_U            = 1e-3,
    double                     diffusivity_D       = 1e-2,
    double                     critical_slope_Sc   = 0.8,
    bool                       nonlinear_profile   = false,
    double                     init_amplitude      = 50.0,
    int                        perturb_seed        = 1,
    int                        n_threads           = 0
);

// Convenience overload
lsdtt_hillslope_strip_surface generate_hillslope_strip_surface(
    int                        n_rows,
    int                        n_cols,
    double                     cell_size,
    lsdtt_strip_channel_config channel_config      = lsdtt_strip_channel_config::both_ends,
    lsdtt_strip_init_profile   init_profile        = lsdtt_strip_init_profile::flat,
    bool                       cross_axis_is_x     = true,
    double                     baselevel_elevation = 0.0,
    double                     half_width          = 0.0,
    double                     uplift_U            = 1e-3,
    double                     diffusivity_D       = 1e-2,
    double                     critical_slope_Sc   = 0.8,
    bool                       nonlinear_profile   = false,
    double                     init_amplitude      = 50.0,
    int                        perturb_seed        = 1,
    int                        n_threads           = 0
);

// -----------------------------------------------------------------------------
// lsdtt_crater_spec
//
// One closed depression, specified analytically.
//
// centre_x, centre_y
//   Centre of the crater in the RASTER CRS (easting / northing), not in pixels
//   and not relative to the grid corner. For a raster built by the convenience
//   overload below the origin is the default synthetic georeference, so the
//   natural way to place a crater is meta.x_min + fraction * n_cols * cell_size.
//   Placing the centre exactly on a PIXEL CENTRE is strongly recommended when
//   the analytic pit is used as an exact oracle -- see lsdtt_crater_field_surface.
//
// radius_m
//   Rim radius. The rim is the circle at which the crater surface rejoins the
//   background plane, so the depression is exactly this wide at its top.
//
// depth_m
//   Rim elevation minus pit elevation, measured against the background plane
//   AT THE CRATER CENTRE. With a tilted background the usable (spill minus pit)
//   relief is smaller -- see spill_elevation below.
//
// floor_fraction
//   0    - the crater is a cone converging on a single point, so the pit is a
//          unique cell.
//   > 0  - the crater has an exactly flat circular floor of radius
//          floor_fraction * radius_m, all of it at pit_elevation. Every floor
//          cell is an equally valid pit, which is the intended input for a
//          deterministic tie-break test. Must be < 1.
// -----------------------------------------------------------------------------
struct lsdtt_crater_spec
{
    double centre_x       = 0.0;   // raster CRS
    double centre_y       = 0.0;
    double radius_m       = 0.0;   // rim radius
    double depth_m        = 0.0;   // rim elevation minus pit elevation
    double floor_fraction = 0.0;   // 0 = cone to a point; >0 = flat floor of
                                   // this fraction of radius
};

// -----------------------------------------------------------------------------
// lsdtt_crater_field_surface
//
// A synthetic DEM containing closed depressions whose geometry is known in
// closed form, plus the ground truth for each. This is the oracle for internal
// drainage work: the pit locations, depths, spill elevations, areas and volumes
// are all derived FROM THE CRATER SPEC, never measured off the raster, so a
// measurement taken off the raster can be checked against them.
//
// dem
//   The elevation surface. Every cell is valid (no nodata).
//
// sink_mask
//   1.0 at each analytic pit cell, 0.0 everywhere else. Suitable as the sink
//   set for fill_pits_with_sinks() once converted to flat indices -- pit_flats
//   below already is that set.
//
// crater_id
//   The index of the crater (0-based, as a float) at every cell strictly inside
//   that crater's rim; nodata elsewhere. Craters cannot overlap, so this is
//   unambiguous.
//
// pit_rc, pit_flats
//   The cell CONTAINING the analytic apex, one entry per crater, in spec order.
//   Because the containing cell of a point on a square grid is also the cell
//   whose centre is nearest that point, this is the raster minimum of the
//   crater whenever the centre sits at or near a pixel centre. If the centre is
//   placed close to a cell corner AND the background is tilted, the true raster
//   minimum can be an adjacent cell, so put centres on pixel centres for exact
//   tests.
//
// pit_elevation
//   background_elevation + background_slope * (centre_y - meta.y_min)
//     - depth_m
//   i.e. the background plane evaluated at the crater centre, less the depth.
//   Exact. The DEM value at the pit CELL is >= this, exceeding it by at most
//   the crater's radial gradient times half a cell diagonal, i.e. about
//   depth_m / (radius_m * (1 - floor_fraction)) * 0.71 * cell_size. For a flat
//   floor wider than one cell the pit cell value is exact.
//
// spill_elevation
//   The lowest point of the rim, which with a background tilted by
//   background_slope is the down-tilt side of the rim circle:
//     pit_elevation + depth_m - background_slope * radius_m
//   With background_slope == 0 this is simply pit_elevation + depth_m.
//   Water standing at this level is on the point of escaping the depression.
//
// area_m2
//   Planform area of the CLOSED DEPRESSION, i.e. the region where the surface
//   lies below spill_elevation. With background_slope == 0 this is exactly
//   pi * radius_m^2. With a tilt the depression is slightly smaller than the
//   rim circle on the down-tilt side and the boundary is no longer a circle;
//   the value returned is then the exact area of that region, obtained by
//   integrating the closed-form boundary radius over azimuth with a fixed
//   deterministic quadrature (relative error ~1e-12, and independent of the
//   raster grid).
//
// volume_m3
//   Fill volume of the closed depression, i.e. the volume of water it holds
//   when filled to spill_elevation. With background_slope == 0 and
//   floor_fraction == 0 this is exactly the cone volume
//   pi * radius_m^2 * depth_m / 3. Otherwise it is the same closed-form radial
//   integral, integrated over azimuth as for area_m2.
//
// TOLERANCE a caller should expect against a raster measurement: the raster
// resolves the depression boundary only to the nearest cell, so a pixel-counted
// area differs from area_m2 by roughly the boundary length times the cell size,
// i.e. of order 2 * pi * radius_m * cell_size, and a summed fill volume differs
// from volume_m3 by roughly that area times the local depth. One cell of slop
// on every linear dimension is the right expectation; do not compare to better
// than that.
// -----------------------------------------------------------------------------
struct lsdtt_crater_field_surface
{
    lsdtt_raster_cpu dem;
    lsdtt_raster_cpu sink_mask;   // 1 at each analytic pit cell, else 0
    lsdtt_raster_cpu crater_id;   // per-pixel crater id inside each rim, else nodata

    std::vector<lsdtt_pixel_rc> pit_rc;         // analytic pit, one per crater
    std::vector<std::int64_t>   pit_flats;      // same, as flat indices
    std::vector<double>         pit_elevation;  // analytic
    std::vector<double>         spill_elevation;
    std::vector<double>         area_m2;
    std::vector<double>         volume_m3;
};

// -----------------------------------------------------------------------------
// generate_crater_field_surface
//
// A tilted background plane into which a set of non-overlapping closed
// depressions has been cut. Deterministic: there is no RNG anywhere in this
// generator, so the same specs always give the same surface bit for bit.
//
// The background plane rises to the north:
//
//   z_bg(x, y) = background_elevation + background_slope * (y - meta.y_min)
//
// so with a positive background_slope the surface drains south and every cell
// outside a crater has a downhill path to the southern DEM edge.
//
// Inside a crater of radius R with a floor radius F = floor_fraction * R, let d
// be the distance from the crater centre and
//
//   w(d) = 0                for d <= F
//        = (d - F)/(R - F)  for F < d < R
//        = 1                for d >= R
//
// then
//
//   z(x, y) = z_bg(cx, cy) + (z_bg(x, y) - z_bg(cx, cy)) * w(d)
//                          - depth_m * (1 - w(d))
//
// which is continuous everywhere, equals the background plane at and beyond the
// rim, and equals pit_elevation over the whole floor. The tilt is faded out
// towards the floor by the same weight that fades the crater in, which is what
// makes a floor_fraction > 0 floor EXACTLY flat rather than merely nearly so.
//
// Constraints, all of which throw std::invalid_argument rather than silently
// producing a surface that is not what it claims to be:
//   * radius_m > 0, depth_m > 0, floor_fraction in [0, 1)
//   * background_slope >= 0
//   * depth_m > background_slope * radius_m * (2 - floor_fraction), so the
//     surface rises strictly outwards along EVERY ray and the depression really
//     is closed. Without this the tilt can outrun the crater wall on the
//     down-tilt side and the "crater" is just a notch in a hillslope.
//   * no two rims may intersect: the centres of any two craters must be more
//     than radius_i + radius_j apart. Overlapping craters would merge into one
//     depression with a pit and a spill neither spec describes.
//   * every rim must lie inside the raster, with one cell of margin, or the
//     depression is not closed.
//
// The generated DEM carries meta.set_default_synthetic_georeference() when the
// raw-dimension overload is used, as every other synthetic generator does.
// -----------------------------------------------------------------------------
lsdtt_crater_field_surface generate_crater_field_surface(
    const lsdtt_raster_meta&              meta,
    const std::vector<lsdtt_crater_spec>& craters,
    double background_slope     = 0.0,     // gentle regional tilt, m/m
    double background_elevation = 1000.0
);

// Convenience overload
lsdtt_crater_field_surface generate_crater_field_surface(
    int    n_rows,
    int    n_cols,
    double cell_size,
    const std::vector<lsdtt_crater_spec>& craters,
    double background_slope     = 0.0,
    double background_elevation = 1000.0
);

} // namespace lsdtt
