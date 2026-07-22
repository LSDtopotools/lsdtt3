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

#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
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

} // namespace lsdtt
