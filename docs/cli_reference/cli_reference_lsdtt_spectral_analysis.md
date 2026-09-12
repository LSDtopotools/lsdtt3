# CLI reference: lsdtt-spectral-analysis

Spectral surface analysis using Fast Fourier Transform and frequency-domain filtering. Computes power spectra, applies wavelet-based smoothing, and performs spectral filters on DEMs. Requires FFTW3.

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

Keys use `snake_case`. The table below documents the current Phase 5 interface;
full reference regeneration is scheduled for Phase 6.

---

## Requirements

Requires FFTW3 library. If FFTW3 is not available at compile time, this CLI is not built. Verify availability: `lsdtt-spectral-analysis -h` should execute without error.

---

## Pipeline overview

Steps run in this order:

| Step | Description |
|---|---|
| 1. DEM preparation | Load DEM, verify projection |
| 2. FFT | Compute 2D Fast Fourier Transform |
| 3. Frequency-domain analysis | Compute power spectrum, apply filters |
| 4. IFFT | Transform back to spatial domain |
| 5. Output | Write filtered rasters and spectral diagnostics |

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
| `n_background_iterations` | int | `20` | Synthetic-background Monte Carlo iterations. |
| `background_seed` | int | `-1` | Background seed; `-1` selects the reproducible built-in sequence. |

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

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format: `cog`, `gtiff`, `same_as_input`, `png`. `png` writes an 8-bit scaled image with a worldfile and PAM sidecar for quick visualisation only, not for analysis reuse. Output files take the extension of the chosen format. |

---

## Notes

- FFT assumes periodic boundary conditions on the DEM; edge artifacts may appear.
- Power spectrum is in log-log space; useful for identifying characteristic wavelengths in topography.

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
| `n_background_iterations` | `int` | `20` | `param` | - | Number of iterations for synthetic background spectrum estimation. | Each iteration generates a random phase spectrum with the same amplitude envelope. More iterations give a smoother background. Default is 20. |
| `background_seed` | `int` | `-1` | `param` | - | Base seed for the Monte Carlo background spectrum. | Iteration i uses seed = background_seed + i. Default (-1) uses the built-in sequence (base 1000), which is always reproducible. Set an explicit non-negative integer to pin a specific realisation - record the value alongside results to reproduce them exactly. |

### Output toggles

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_2d_spectrum_rasters` | `bool` | `true` | `switch` | - | Write the 2D observed and background PSD rasters. | Writes _psd_2d.&lt;ext&gt; and _background_psd_2d.&lt;ext&gt; in frequency space; write_output_format selects the extension. |
| `write_radial_spectrum` | `bool` | `false` | `switch` | - | Write the log-binned radial spectrum as a CSV. | Writes _radial_spectrum.csv with columns: frequency, wavelength, psd, background_psd, normalised_psd, beta_by_bin, r_sq_by_bin. |
| `write_significant_frequencies` | `bool` | `true` | `switch` | - | Write a CSV of frequencies with significant spectral power. | Writes _significant_frequencies.csv. A frequency is significant if its normalised PSD exceeds the chi-squared 95% CI threshold (2.9955). Also includes rollover summary lines. |

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

<!-- END GENERATED CLI PARAMETERS -->
