# CLI reference: lsdtt-spectral-analysis

Spectral analysis of a raster after Perron, Kirchner and Dietrich (2008): finds characteristic spatial scales, such as ridge-valley spacing, dune or vegetation-pattern spacing, and their orientations. It computes the 2D power spectrum, compares it with the spectrum of synthetic fractal surfaces, and reports the frequencies and bearings that stand out. It is analysis only: it does not filter or write a filtered raster. Requires FFTW3. The method, its equations and its tests are described in [Spectral analysis](../design_and_algorithms/spectral_analysis.md).

## Usage

```
lsdtt-spectral-analysis  <path> <parameter_file>
lsdtt-spectral-analysis  <parameter_file>        # path defaults to ./
lsdtt-spectral-analysis  -h                      # print help and exit
lsdtt-spectral-analysis  -v                      # print version and exit
lsdtt-spectral-analysis  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

Keys use `snake_case`. The full parameter list is in [Parameters](#parameters),
generated from the CLI's own help.

---

## Requirements

Requires FFTW3 library. If FFTW3 is not available at compile time, this CLI is not built. Verify availability: `lsdtt-spectral-analysis -h` should execute without error.

---

## Pipeline overview

Steps run in this order:

| Step | Description |
|---|---|
| 1. Prepare | Load the raster, remove a best-fit plane, apply a Hann window, zero-pad to a power of two in each direction |
| 2. Periodogram | 2D FFT and periodogram; collapse to a radial spectrum |
| 3. Rollover | Log-bin the radial spectrum, find the rollover and the power-law exponent beta |
| 4. Background | Average the periodograms of `n_background_iterations` synthetic fractal surfaces with that exponent |
| 5. Normalise and test | Observed / background; flag binned frequencies above 2.9955 and 2D cells above `-ln(1 - confidence)` |
| 6. Output | Write the 2D spectrum rasters and the CSVs below |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Spectral analysis

`n_threads` defaults to `0`, which uses all available hardware threads. Use
`1` for serial execution, or set an explicit limit no greater than the allocated
physical cores on a shared machine.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | OpenMP threads; `0` uses all available hardware threads. |
| `spectral_log_bin_width` | float | `0.05` | Log10 bin width for the radial spectrum. |
| `n_background_iterations` | int | `100` | Number of synthetic fractal surfaces averaged into the background. More iterations bring the 2D test's false-alarm rate closer to nominal: at 9.21 it is about 1.5 in 10,000 cells at 100, 5 in 10,000 at 20 (the default before 2026-10-21). |
| `background_seed` | int | `-1` | Background seed; `-1` selects the reproducible built-in sequence (base 1000). For a given seed every output is bit-identical, whatever `n_threads` is. |
| `significant_frequencies_2d_confidence` | float | `0.9999` | Confidence level for `_significant_frequencies_2d.csv`; the normalised-PSD threshold is `-ln(1 - confidence)` (9.2103 at 0.9999). |

---

## Valley spacing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_valley_spacing` | bool | `false` | Measure characteristic valley spacing. Former `compute_valley_spacing` is an unchanged-value alias. |
| `valley_spacing_axis` | string | `rows` | Axis along which spacing is measured: `rows` or `columns`. |
| `valley_spacing_margin` | int | `20` | Flow-direction boundary cells omitted at each end. |
| `valley_spacing_K` | float | `0` | Erodibility for the channelization index; `0` skips it. |
| `valley_spacing_D` | float | `0` | Hillslope diffusivity for the channelization index. |
| `valley_spacing_U` | float | `0` | Uplift rate for the channelization index. |
| `valley_spacing_area_exponent` | float | `0.45` | Dimensionless drainage-area exponent m. Former `valley_spacing_m` is an unchanged-value alias. |
| `valley_spacing_n` | float | `1.0` | Dimensionless slope exponent n. |

---

## Outputs

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_2d_spectrum_rasters` | bool | `true` | Write observed and background 2D PSD rasters. Former `write_2d_spectrum` is an unchanged-value alias. |
| `write_radial_spectrum` | bool | `false` | Write the log-binned radial spectrum CSV. |
| `write_significant_frequencies` | bool | `true` | Write significant frequencies and rollover summaries. |
| `write_significant_frequencies_2d` | bool | `true` | Write the individual significant 2D spectrum cells, with bearings. |

### Significant 2D frequencies with bearings

`_significant_frequencies_2d.csv` lists every cell of the 2D spectrum whose
normalised PSD (observed / background) exceeds `-ln(1 - confidence)`, the
χ²(2 DOF) quantile divided by 2 (Perron et al. 2008, eq. 9). It replaces the
legacy LSDTopoTools `_significant_freq.csv`, at the same default 99.99% level.

**Where 9.21 comes from.** For a surface with no characteristic scale, each
periodogram cell divided by its expected value (the normalised PSD) follows an
exponential distribution with mean 1: it is a sum of two squared Gaussians (the
real and imaginary parts of the Fourier coefficient), i.e. χ² with 2 degrees of
freedom, divided by 2. The chance that such a cell exceeds `t` is `exp(-t)`.
Setting that chance to 1 in 10,000 (confidence 0.9999) gives
`t = -ln(0.0001) = 9.2103`. The legacy code wrote the same value as
`18.42185 / 2`, half the χ²₂ value at 99.99 %. The binned test's 2.9955 is the
same construction at 95 % (`5.991 / 2`). See
[Spectral analysis, section 2.11](../design_and_algorithms/spectral_analysis.md)
for the derivation and for the effect of a finite number of background
iterations.

| Column | Units | Meaning |
|---|---|---|
| `bearing_rad`, `bearing_deg` | rad, ° | Direction of the wave vector, clockwise from north, in [0, 2π) / [0, 360). The crests run perpendicular to it. |
| `frequency_cycles_per_m` | cycles/m | Radial frequency of the cell. |
| `wavelength_m` | m | 1 / frequency. |
| `normalised_psd` | - | Observed PSD / background PSD. |

Five `#` comment lines at the top record the confidence, threshold,
effective background seed, number of background iterations and
`n_significant_frequencies`. When nothing passes the threshold, a further `#`
line under the column header says so in words:

```
# n_significant_frequencies,0
bearing_rad,bearing_deg,frequency_cycles_per_m,wavelength_m,normalised_psd
# No significant frequencies were found in this raster: no cell of the 2D spectrum exceeded the normalised PSD threshold of 9.21034 (confidence 0.9999). This is a valid result, not an error.
```

`_significant_frequencies.csv` does the same, with its threshold (2.9955)
and `n_significant_frequencies` added to its rollover lines. Read either file
with `pandas.read_csv(path, comment="#")`.

- The spectrum of a real surface is point-symmetric, so a significant
  frequency normally appears twice, at θ and θ + 180°. Use `bearing_deg % 180`
  for an axial (orientation) analysis, or keep one of each pair.
- The exception is the Nyquist edges of the padded grid, i.e. wavelengths near
  two cells. +Nyquist and −Nyquist are the same frequency, so the twin of a cell
  on the north Nyquist row is on that row too, at 360° − θ, and the twin of a
  cell on the east Nyquist column is on that column, at 180° − θ. Three cells
  (bearings 0°, 270° and 315° at the Nyquist corners) are their own twins and
  appear once.
- The bearing uses the physical frequency components, so it is correct when
  the padded grid is not square. The legacy code used index offsets, which
  agree only for square grids.
- The 2D spectrum has thousands of cells, so a few pass the threshold by
  chance alone. 1 in 10,000 at 0.9999 is the limit for infinitely many
  background iterations; with the default 100 the rate is about 1.5 in 10,000,
  i.e. about 2.4 entries (one or two pairs) on a 128 × 128 padded grid. Repeat
  a run with another `background_seed` and keep detections that persist.

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format: `cog`, `gtiff`, `same_as_input`, `png`. `png` writes an 8-bit scaled image with a worldfile and PAM sidecar for quick visualisation only, not for analysis reuse. Output files take the extension of the chosen format. |

---

## Notes

- The DFT treats the raster as periodic. To limit the leakage from the jump
  between opposite edges, the tile is detrended and multiplied by a Hann window
  that falls to zero at the edges, then zero-padded, before the FFT.
- PSD values are variance per DFT cell, in the input units squared (m² for an
  elevation raster); they are not divided by the frequency-cell area.
- The tile must hold several wavelengths of a pattern for its peak to stand out
  from the lowest-frequency bins.
- The 2D spectrum rasters are in frequency space, with the zero frequency at
  (0, 0) and axes in cycles/m. Each axis is padded to its own power of two, so
  when the padded row and column counts differ (e.g. a 200 × 500 input pads to
  256 × 512) the north and east frequency spacings differ and the rasters are
  written with non-square pixels. This is deliberate: both axes are exact. They
  are the only lsdtt3 outputs with non-square pixels.

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

**Status: stable.**

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### Threading

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all hardware threads (default). 1 = serial. In shared environments, set an explicit limit no greater than the physical-core allocation. |

### Output format

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | Allowed values: cog, gtiff, same_as_input, png. cog: Cloud Optimised GeoTIFF (default). gtiff: plain compressed GeoTIFF. same_as_input: match the input raster format. png: 8-bit scaled image with a worldfile and PAM sidecar for georeferencing - for quick visualisation only, not for analysis reuse. Output files take the extension of the chosen format. |

### Spectral options

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `spectral_log_bin_width` | `float` | `0.05` | `param` | - | Log10 bin width for the radial spectrum. | Controls the width of each frequency bin in log10 space. Smaller values give finer frequency resolution. Default is 0.05. |
| `n_background_iterations` | `int` | `100` | `param` | - | Number of synthetic fractal surfaces averaged into the background spectrum. | Each iteration makes a synthetic surface (white noise filtered to a power-law spectrum with the fitted exponent beta, scaled to the tile's standard deviation) and takes its periodogram with the same window and padding; the background is their mean. Its scatter in each cell is about 1/sqrt(N), and the 2D test's threshold assumes the background is exact: at 9.21 the per-cell false-alarm rate is (1 + 9.21/N)^-N, about 1.5e-4 at the default N = 100 (nominal 1e-4) and 5e-4 at N = 20. Default is 100. |
| `background_seed` | `int` | `-1` | `param` | - | Base seed for the Monte Carlo background spectrum. | Iteration i uses seed = background_seed + i. Default (-1) uses the built-in sequence (base 1000), which is always reproducible. Set an explicit non-negative integer to pin a specific realisation - record the value alongside results to reproduce them exactly. |
| `significant_frequencies_2d_confidence` | `float` | `0.9999` | `param` | - | Confidence level for a 2D spectrum cell to be written as significant. | A cell is significant when its normalised PSD exceeds -ln(1 - confidence), the chi-squared (2 DOF) quantile divided by 2 (Perron et al. 2008, eq. 9). 0.9999 gives 9.2103, as in the legacy LSDTopoTools output; 0.95 gives 2.9957. Must lie strictly between 0 and 1. The spectrum has thousands of cells, so at 0.95 hundreds pass by chance alone. |

### Output toggles

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_2d_spectrum_rasters` | `bool` | `true` | `switch` | - | Write the 2D observed and background PSD rasters. | Writes _psd_2d.&lt;ext&gt; and _background_psd_2d.&lt;ext&gt; in frequency space; write_output_format selects the extension. |
| `write_radial_spectrum` | `bool` | `false` | `switch` | - | Write the log-binned radial spectrum as a CSV. | Writes _radial_spectrum.csv with columns: frequency, wavelength, psd, background_psd, normalised_psd, beta_by_bin, r_sq_by_bin. |
| `write_significant_frequencies` | `bool` | `true` | `switch` | - | Write a CSV of frequencies with significant spectral power. | Writes _significant_frequencies.csv. A frequency is significant if its normalised PSD exceeds the chi-squared 95% CI threshold (2.9955). Lines starting with # give the rollover summary, the threshold and n_significant_frequencies; when that is 0 a line under the column header says no significant frequencies were found. |
| `write_significant_frequencies_2d` | `bool` | `true` | `switch` | - | Write a CSV of individual 2D spectrum cells with significant power, with bearings. | Writes _significant_frequencies_2d.csv with columns bearing_rad, bearing_deg, frequency_cycles_per_m, wavelength_m, normalised_psd. The bearing is the direction of the wave vector, clockwise from grid north (from true north with bearings_relative_to_true_north, converted with the convergence at the grid centre, since the spectrum belongs to the whole raster); the crests run perpendicular to it. The spectrum of a real surface is point-symmetric, so each significant frequency appears twice, at bearings theta and theta + 180, except cells on the Nyquist edges of the padded grid (wavelengths near two cells), whose twin is at 360 - theta (first row) or 180 - theta (first column), or which appear once. Equivalent to the legacy _significant_freq.csv. The threshold is set by significant_frequencies_2d_confidence. Lines starting with # record the confidence, threshold, seed, iterations and n_significant_frequencies; when that is 0 a line under the column header says no significant frequencies were found, so an empty table is never mistaken for a failed run. Read with pandas.read_csv(path, comment='#'). |

### Valley spacing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_valley_spacing` | `bool` | `false` | `switch` | - | Measure characteristic valley spacing (Perron 2008 / Bonetti 2020). | Reports three independent measures because each fails differently: the POWER-WEIGHTED mean of spectral peaks (Perron et al. 2008 -- weighted, not the argmax, because spacing occupies a band); the CONVEX HALF-WIDTH, the mean run of convex cells across the valleys, which uses thousands of samples and so survives a domain too small for a spectral peak; and the small-scale spectral exponent alpha, predicted to be 2m-3 (Bonetti et al. 2020, Hooshyar et al. 2021). Writes _valley_spacing.csv and, with the spectrum, _valley_spacing_spectrum.csv. |
| `valley_spacing_axis` | `string` | `rows` | `param` | - | Axis the spacing is measured ALONG: rows or columns. | Valleys run in the flow direction, so a landscape draining north-south has its spacing measured east-west, i.e. along rows. |
| `valley_spacing_margin` | `int` | `20` | `param` | - | Cells dropped at each end of the FLOW direction before measuring. | Base level pins elevation there, so those lines carry the boundary's imprint rather than the landscape's own scale. |
| `valley_spacing_K` | `float` | `0` | `param` | - | Erodibility K, for the Bonetti channelization index. 0 = skip CI. | CI = K l^(m+n) / (D^n U^(1-n)); for n = 1 it is independent of U. CI_c ~ 32-37 marks first channelization, so reporting CI beside the measured spacing is what makes a result comparable to the literature. |
| `valley_spacing_D` | `float` | `0` | `param` | - | Hillslope diffusivity D, for the channelization index. | See valley_spacing_K. |
| `valley_spacing_U` | `float` | `0` | `param` | - | Uplift rate U, for the channelization index (unused when n = 1). | See valley_spacing_K. |
| `valley_spacing_area_exponent` | `float` | `0.45` | `param` | - | Drainage-area exponent m, for CI and for the alpha = 2m-3 prediction. | See valley_spacing_K. |
| `valley_spacing_n` | `float` | `1.0` | `param` | - | Slope exponent n, for the channelization index. | See valley_spacing_K. |

### Orientations

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `bearings_relative_to_true_north` | `bool` | `false` | `param` | - | Write bearings, aspects and orientations clockwise from true north. | Off (default): they are clockwise from grid north (grid up), so in an oblique frame with the feature along +x, 90 degrees means along the feature. On: each value is converted with the grid convergence at its own location (each cell of a raster, each row of a table), which also removes the up to about 3 degree difference between UTM grid north and true north. Angles between two directions (junction angles, deflections, contrasts) are unchanged, and hillshade illumination always stays relative to grid up. Needs a raster with a CRS. |

<!-- END GENERATED CLI PARAMETERS -->
