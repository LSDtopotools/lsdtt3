# CLI reference: lsdtt-synthetic-surfaces

Generates one-shot analytical, stochastic, and geomorphic fixture surfaces. Use
this tool when the generated raster or its ground-truth companion layer is the
product, rather than an optimised drainage network or an evolved landscape.

## Usage

```
lsdtt-synthetic-surfaces  <path> <parameter_file>
lsdtt-synthetic-surfaces  <parameter_file>        # path defaults to ./
lsdtt-synthetic-surfaces  -h                      # print help and exit
lsdtt-synthetic-surfaces  -v                      # print version and exit
lsdtt-synthetic-surfaces  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

Keys use `snake_case`. All keys and their defaults are listed below.

---

## How generation works

Surface families are selected by independent `write_*_surface` switches. More
than one family can be generated in a run on the same dimensions, cell size,
and random seed. Generated rasters use lsdtt3's default synthetic projected
georeference.

Available families are:

- uniform, red, and blue noise;
- diamond-square fractal terrain;
- periodic cosine waves and global quadratic analytical surfaces;
- terrace staircases with tread IDs and a trunk-channel mask;
- alluvial fans with inset-lobe IDs;
- hillslope strips with channel masks and controlled initial profiles.

The two analytical families can optionally be passed through
`wavelet_quadratic_fit`. The tool reports coefficient errors and can write true,
fitted, and residual coefficient rasters. These diagnostics warn rather than
acting as assertion-based unit tests.

## Principal outputs

| Family | Principal output suffixes |
|---|---|
| Uniform/red/blue noise | `_uniform_random`, `_red_noise`, `_blue_noise` |
| Diamond-square | `_diamond_square` |
| Periodic wave | `_periodic_wave`, `_periodic_wave_components.csv` |
| Global quadratic | `_global_quadratic`, `_global_quadratic_coefficients.csv` |
| Terrace staircase | `_terrace_staircase`, `_terrace_staircase_tread_id`, `_terrace_staircase_channel_mask` |
| Alluvial fan | `_alluvial_fan`, `_alluvial_fan_lobe_id` |
| Hillslope strip | `_hillslope_strip`, `_hillslope_strip_channel_mask` |
| Wavelet validation | coefficient rasters and `_wavelet_test_summary.csv` |

Raster extensions follow `write_output_format`. See
[Choosing a synthetic landscape tool](../design_and_algorithms/synthetic_cli_cluster.md)
for the distinction from `lsdtt-synthetic-networks` and `lsdtt-lem`.

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### Domain and output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_rows` | `int` | `128` | `param` | - | Number of rows in the synthetic raster. | - |
| `n_cols` | `int` | `128` | `param` | - | Number of columns in the synthetic raster. | - |
| `cell_size_metres` | `float` | `1.0` | `param` | - | Cell size in metres. | Generated rasters use the default synthetic projected georeference. wavelet_quadratic_fit() is called with the supplied cell_size_metres. |
| `random_seed` | `int` | `33` | `param` | - | Random seed for all stochastic generators. | - |
| `n_threads` | `int` | `0` | `io` | - | Number of CPU threads (0 = all available). | - |
| `write_output_format` | `string` | `cog` | `io` | - | Output format: cog, gtiff, same_as_input. | cog falls back to gtiff for rasters below ~262k pixels. |

### Uniform random surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_uniform_random_surface` | `bool` | `false` | `switch` | - | Write a uniform random noise surface. | - |
| `uniform_min_value` | `float` | `-1.0` | `param` | `write_uniform_random_surface` | Minimum value of the uniform random surface. | - |
| `uniform_max_value` | `float` | `1.0` | `param` | `write_uniform_random_surface` | Maximum value of the uniform random surface. | - |

### Red noise surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_red_noise_surface` | `bool` | `false` | `switch` | - | Write a red (low-frequency) noise surface. | - |
| `red_noise_mean` | `float` | `0.0` | `param` | `write_red_noise_surface` | Mean value of the red noise surface. | - |
| `red_noise_stddev` | `float` | `1.0` | `param` | `write_red_noise_surface` | Standard deviation of the red noise surface. | - |
| `red_noise_blur_passes` | `int` | `6` | `param` | `write_red_noise_surface` | Number of smoothing passes used to make red noise. | - |

### Blue noise surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_blue_noise_surface` | `bool` | `false` | `switch` | - | Write a blue (high-frequency) noise surface. | - |
| `blue_noise_mean` | `float` | `0.0` | `param` | `write_blue_noise_surface` | Mean value of the blue noise surface. | - |
| `blue_noise_stddev` | `float` | `1.0` | `param` | `write_blue_noise_surface` | Standard deviation of the blue noise surface. | - |
| `blue_noise_blur_passes` | `int` | `6` | `param` | `write_blue_noise_surface` | Number of smoothing passes used to make blue noise. | - |

### Diamond-square surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_diamond_square_surface` | `bool` | `false` | `switch` | - | Write a fractal diamond-square terrain surface. | - |
| `diamond_square_roughness` | `float` | `0.55` | `param` | `write_diamond_square_surface` | Roughness decay per level for diamond-square. | 0.5 gives Brownian noise; lower = smoother; higher = rougher. |
| `diamond_square_initial_scale` | `float` | `1.0` | `param` | `write_diamond_square_surface` | Initial perturbation amplitude for diamond-square. | - |
| `diamond_square_mean` | `float` | `0.0` | `param` | `write_diamond_square_surface` | Target mean of the diamond-square surface. | - |
| `diamond_square_stddev` | `float` | `1.0` | `param` | `write_diamond_square_surface` | Target standard deviation of the diamond-square surface. | - |
| `diamond_square_parabola_relief` | `float` | `0.0` | `param` | `write_diamond_square_surface` | Peak relief (m) of an N-S parabola superimposed on the diamond-square surface. 0 disables it. | A raw diamond-square surface has no regional gradient, so it is full of closed depressions; filling them creates flats, and routing across a flat gives feathered parallel drainage that is an artefact of the fill. The parabola is zero at the north and south boundaries and peaks at the centre row, giving a regional slope that drains to both boundaries so most depressions never form. Pair with diamond_square_noise_relief. CALIBRATION: measured on a 400x400 diamond-square of 537 m relief, cells lying in depressions went 16.2% (no parabola) -&gt; 12.0% (300 m) -&gt; 2.7% (1075 m), and the MEAN FILL DEPTH, which is what actually governs feathering, went 10.3 m -&gt; 3.2 m -&gt; 1.3 m. But 1075 m (2x the fractal relief) over-corrects: the parabola then dominates and the channels run near-parallel down the regional slope, losing the dendritic structure the fractal was there to provide. About 0.5x the fractal relief was the useful setting - enough regional gradient to organise the drainage, not enough to impose its direction. |
| `diamond_square_noise_relief` | `float` | `0.0` | `param` | `write_diamond_square_surface` | Amplitude (m) of white noise added to the diamond-square surface. 0 disables it. | Roughens whatever depressions survive the parabola, so that when they are filled the result still has a gradient to route on rather than a flat. Applied everywhere, not only inside depressions. SIZE IT AGAINST THE MEAN FILL DEPTH, not against the relief: the noise has to be able to break up a filled pool, so it must be comparable to or larger than the depth the fill adds. In the calibration above, 5 m of noise cleared the feathering when the parabola had brought the mean fill to 3.2 m, while 2 m did not when the mean fill was 10.3 m. Small enough not to disturb the fractal structure otherwise. |

### Periodic wave surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_periodic_wave_surface` | `bool` | `false` | `switch` | - | Write a periodic cosine-mode test surface. | - |
| `periodic_n_components` | `int` | `8` | `param` | `write_periodic_wave_surface` | Number of cosine modes in the periodic test surface. | - |
| `periodic_mode_min` | `int` | `1` | `param` | `write_periodic_wave_surface` | Minimum mode number for the periodic surface. | - |
| `periodic_mode_max` | `int` | `6` | `param` | `write_periodic_wave_surface` | Maximum mode number for the periodic surface. | - |
| `periodic_amplitude_min` | `float` | `1.0` | `param` | `write_periodic_wave_surface` | Minimum mode amplitude. | - |
| `periodic_amplitude_max` | `float` | `10.0` | `param` | `write_periodic_wave_surface` | Maximum mode amplitude. | - |
| `periodic_random_sign` | `bool` | `true` | `param` | `write_periodic_wave_surface` | Allow random sign on mode amplitudes. | - |

### Global quadratic surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_global_quadratic_surface` | `bool` | `false` | `switch` | - | Write a globally quadratic test surface. | - |
| `quadratic_abs_max` | `float` | `1e-3` | `param` | `write_global_quadratic_surface` | Max absolute value for quadratic coefficients a, b, c. | - |
| `quadratic_linear_abs_max` | `float` | `1.0` | `param` | `write_global_quadratic_surface` | Max absolute value for linear coefficients d, e. | - |
| `quadratic_constant_abs_max` | `float` | `100.0` | `param` | `write_global_quadratic_surface` | Max absolute value for constant coefficient f. | - |
| `quadratic_random_sign` | `bool` | `true` | `param` | `write_global_quadratic_surface` | Allow random sign on quadratic coefficients. | - |

### Wavelet validation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `test_wavelet_on_periodic_wave` | `bool` | `false` | `switch` | - | Run wavelet_quadratic_fit() on the periodic surface and report errors. | The periodic surface is an approximate (not exact) oracle. |
| `test_wavelet_on_global_quadratic` | `bool` | `false` | `switch` | - | Run wavelet_quadratic_fit() on the global quadratic surface and report errors. | The global quadratic is the exact oracle for wavelet_quadratic_fit(). |
| `write_true_coefficient_rasters` | `bool` | `false` | `switch` | - | Write analytic coefficient rasters alongside the DEM. | - |
| `write_wavelet_fitted_coefficient_rasters` | `bool` | `false` | `switch` | - | Write wavelet-fitted coefficient rasters. | - |
| `write_wavelet_residual_rasters` | `bool` | `false` | `switch` | - | Write residual (fitted - analytic) rasters for each coefficient. | - |
| `polyfit_lambda_metres` | `float` | `100` | `param` | - | Wavelength lambda in metres for wavelet_quadratic_fit(). | Should be &gt;= cell_size_metres. Larger values give smoother fits. |

### Terrace staircase surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_terrace_staircase_surface` | `bool` | `false` | `switch` | - | Write a symmetrical flight of fluvial terraces either side of a trunk channel, for height-stratification tests. | Terrace j sits j riser heights above the local channel and occupies the band [(j-1)*tread_width, j*tread_width) from the channel on both flanks. |
| `write_terrace_staircase_tread_id_raster` | `bool` | `true` | `switch` | `write_terrace_staircase_surface` | Write the ground-truth tread id raster (1..n_terraces on treads, 0 on the channel, on risers, and inside a riser gap). | - |
| `write_terrace_staircase_channel_mask` | `bool` | `true` | `switch` | `write_terrace_staircase_surface` | Write the trunk channel mask (1 = channel). | - |
| `terrace_n_terraces` | `int` | `4` | `param` | `write_terrace_staircase_surface` | Number of terrace treads on each flank. | - |
| `terrace_tread_width_metres` | `float` | `300` | `param` | `write_terrace_staircase_surface` | Cross-valley width of each tread in metres. | - |
| `terrace_riser_height_metres` | `float` | `8` | `param` | `write_terrace_staircase_surface` | Height of each riser in metres. | Set this BELOW the flat threshold to make every riser pass the flat test, which reproduces wholesale under-segmentation with no riser gaps at all. |
| `terrace_valley_gradient` | `float` | `0.002` | `param` | `write_terrace_staircase_surface` | Downstream fall of the channel in m/m. | - |
| `terrace_riser_height_downstream_taper` | `float` | `0.0` | `param` | `write_terrace_staircase_surface` | Fraction of the riser height lost from the upstream to the downstream end; 0 gives a flight of constant height. | 0.667 with an 12 m riser runs 12 m upstream to 4 m downstream, so a reach-relative ladder can be compared with a global one. |
| `terrace_tread_cross_gradient` | `float` | `0.002` | `param` | `write_terrace_staircase_surface` | Gradient (m/m) at which each tread falls toward the channel. | NOT cosmetic: at 0 the treads drain downstream off the grid instead of to the channel, and height above channel is then undefined over most of the DEM. Must be less than riser_height_m / tread_width_m. |
| `terrace_n_riser_gaps` | `int` | `0` | `param` | `write_terrace_staircase_surface` | Number of road- or spill-like cuts that break riser continuity. | At each gap the stepped cross-section is replaced by a continuous ramp, so the gap reads as flat and connects treads that are metres apart in height. |
| `terrace_riser_gap_width_metres` | `float` | `120` | `param` | `write_terrace_staircase_surface` | Along-valley length of each riser gap in metres. | Must comfortably exceed the roughness window of whatever detector is being tested, or a window centred in the gap still reaches the stepped rows either side and never reads as flat. At a 40 m window a 30 m gap returns 4.5 m of local relief; a 120 m gap returns about 1 m. |
| `terrace_noise_amplitude_metres` | `float` | `0` | `param` | `write_terrace_staircase_surface` | Symmetric per-pixel white noise amplitude, metres. | 0 gives an exact surface whose tread elevations can be asserted to the metre. |
| `terrace_channel_elevation` | `float` | `100` | `param` | `write_terrace_staircase_surface` | Channel elevation at the upstream edge, metres. | - |
| `terrace_valley_wall_gradient` | `float` | `0.5` | `param` | `write_terrace_staircase_surface` | Gradient of the valley wall beyond the outermost tread, m/m. | - |

### Alluvial fan surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_alluvial_fan_surface` | `bool` | `false` | `switch` | - | Write a radial alluvial fan with younger lobes inset below the older surface. | A fan climbs continuously in height above its channel from toe to apex, so a height-above-channel ladder cuts it into arcs; the lobes sit a constant depth below the older surface, so a residual to a fitted fan separates them cleanly. |
| `write_alluvial_fan_lobe_id_raster` | `bool` | `true` | `switch` | `write_alluvial_fan_surface` | Write the ground-truth lobe id raster (0 = older surface, 1..n_inset_lobes = inset lobes). | - |
| `fan_apex_x_fraction` | `float` | `0.5` | `param` | `write_alluvial_fan_surface` | Apex easting as a fraction of the grid width. | - |
| `fan_apex_y_fraction` | `float` | `0.95` | `param` | `write_alluvial_fan_surface` | Apex northing as a fraction of the grid height, increasing southward. | The default puts the apex near the south edge so the fan spreads north across the grid. |
| `fan_apex_elevation` | `float` | `500` | `param` | `write_alluvial_fan_surface` | Elevation at the apex, metres. | - |
| `fan_gradient` | `float` | `0.03` | `param` | `write_alluvial_fan_surface` | Radial fall at the apex in m/m. | - |
| `fan_concavity` | `float` | `0.0` | `param` | `write_alluvial_fan_surface` | Reduction of the radial gradient with distance; 0 gives a straight cone. | Must be less than fan_gradient. A positive value gives the concave-up profile a real fan has, which a degree-2 polynomial fits and a plane does not. |
| `fan_n_inset_lobes` | `int` | `2` | `param` | `write_alluvial_fan_surface` | Number of inset lobes radiating from the apex. | - |
| `fan_lobe_inset_metres` | `float` | `4` | `param` | `write_alluvial_fan_surface` | Depth of each lobe below the older surface, metres. | - |
| `fan_lobe_half_angle_degrees` | `float` | `25` | `param` | `write_alluvial_fan_surface` | Half-angle of each lobe wedge in degrees. | - |
| `fan_noise_amplitude_metres` | `float` | `0` | `param` | `write_alluvial_fan_surface` | Symmetric per-pixel white noise amplitude, metres. | - |

### Hillslope strip surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_hillslope_strip_surface` | `bool` | `false` | `switch` | - | Write a 1-D hillslope-strip DEM (and channel mask) for the transient hillslope solver validation tests. | - |
| `write_hillslope_strip_channel_mask` | `bool` | `true` | `switch` | `write_hillslope_strip_surface` | Write the geometry-defined channel mask (1 = channel/Dirichlet sink, 0 = hillslope) alongside the strip DEM. | - |
| `strip_channel_config` | `string` | `both_ends` | `param` | `write_hillslope_strip_surface` | Channel placement: both_ends (divide at centre), single_end_low, or single_end_high (no-flux divide at the opposite end). | - |
| `strip_init_profile` | `string` | `flat` | `param` | `write_hillslope_strip_surface` | Initial profile: flat, analytic_parabola, raised_plateau, perturbed, or linear_ramp. | analytic_parabola is the steady fixed point; raised_plateau / perturbed seed transient-relaxation tests. |
| `strip_cross_axis_is_x` | `bool` | `true` | `param` | `write_hillslope_strip_surface` | true: elevation varies E-W (channels are N-S end columns); false: varies N-S (channels are E-W end rows). | - |
| `strip_nonlinear_profile` | `bool` | `false` | `param` | `write_hillslope_strip_surface` | For analytic_parabola: true uses the nonlinear Roering (2007) steady profile; false uses the linear parabola. | - |
| `strip_baselevel_elevation` | `float` | `0.0` | `param` | `write_hillslope_strip_surface` | Channel (Dirichlet) baselevel elevation in metres. | - |
| `strip_half_width_metres` | `float` | `0.0` | `param` | `write_hillslope_strip_surface` | Hillslope half-width L (divide-to-channel distance, m) for the analytic profile. &lt;= 0 auto-derives L from the grid (recommended; makes analytic_parabola a discrete fixed point). | - |
| `strip_uplift` | `float` | `1e-3` | `param` | `write_hillslope_strip_surface` | Uplift rate U (m/yr) used to build analytic_parabola. | - |
| `strip_diffusivity` | `float` | `1e-2` | `param` | `write_hillslope_strip_surface` | Diffusivity D (m^2/yr) used to build analytic_parabola. | - |
| `strip_critical_slope` | `float` | `0.8` | `param` | `write_hillslope_strip_surface` | Critical slope Sc (m/m); used only for the nonlinear profile. | - |
| `strip_init_amplitude` | `float` | `50.0` | `param` | `write_hillslope_strip_surface` | Height scale (m): plateau height (raised_plateau), noise amplitude (perturbed), or crest height (linear_ramp). | - |

<!-- END GENERATED CLI PARAMETERS -->
