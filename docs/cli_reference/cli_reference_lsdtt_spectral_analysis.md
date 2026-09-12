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
