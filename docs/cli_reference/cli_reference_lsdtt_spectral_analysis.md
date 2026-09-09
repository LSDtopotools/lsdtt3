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

Keys use `snake_case`. All keys and their defaults are listed below.

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

## Spectral filtering

| Parameter | Type | Default | Description |
|---|---|---|---|
| `apply_lowpass_filter` | bool | `false` | Apply low-pass (smoothing) filter in frequency domain. |
| `lowpass_wavelength_metres` | float | `100.0` | Wavelength (metres) of low-pass cutoff. Removes shorter-wavelength features. |
| `apply_highpass_filter` | bool | `false` | Apply high-pass (edge enhancement) filter. |
| `highpass_wavelength_metres` | float | `50.0` | Wavelength (metres) of high-pass cutoff. Removes longer-wavelength features. |
| `apply_bandpass_filter` | bool | `false` | Apply band-pass filter. |
| `bandpass_low_wavelength_metres` | float | `50.0` | Lower wavelength cutoff (metres). |
| `bandpass_high_wavelength_metres` | float | `500.0` | Upper wavelength cutoff (metres). |

---

## Wavelet analysis

| Parameter | Type | Default | Description |
|---|---|---|---|
| `compute_wavelet_transform` | bool | `false` | Compute continuous wavelet transform (CWT). |
| `wavelet_type` | string | `morlet` | Wavelet basis: `morlet`, `mexican_hat`, `gaussian`. |
| `wavelet_scales` | vector of floats | `[10, 20, 50, 100, 200]` | Scales (metres) to compute CWT at. |

---

## Power spectrum

| Parameter | Type | Default | Description |
|---|---|---|---|
| `compute_power_spectrum` | bool | `false` | Compute and write radial power spectrum. |
| `power_spectrum_radial_bins` | int | `256` | Number of radial frequency bins for power spectrum. |

---

## Output: Filtered rasters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_lowpass_filtered_dem` | bool | `false` | Write low-pass filtered DEM. |
| `write_highpass_filtered_dem` | bool | `false` | Write high-pass filtered DEM (edges). |
| `write_bandpass_filtered_dem` | bool | `false` | Write band-pass filtered DEM. |

Output files:
- `<write_fname>_lowpass_filtered.tif`
- `<write_fname>_highpass_filtered.tif`
- `<write_fname>_bandpass_filtered.tif`

---

## Output: Spectral diagnostics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_power_spectrum_csv` | bool | `false` | Write radial power spectrum CSV (wavenumber vs power). |
| `write_wavelet_transform_csv` | bool | `false` | Write wavelet transform coefficients CSV. |

Output files:
- `<write_fname>_power_spectrum.csv`
- `<write_fname>_wavelet_coefficients.csv`

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format: `cog`, `gtiff`, `same_as_input`, `png`. `png` writes an 8-bit scaled image with a worldfile and PAM sidecar for quick visualisation only, not for analysis reuse. Output files take the extension of the chosen format. |

---

## Notes

- FFT assumes periodic boundary conditions on the DEM; edge artifacts may appear.
- Power spectrum is in log-log space; useful for identifying characteristic wavelengths in topography.
- Wavelet transform is computationally expensive; consider reducing number of scales or using on subsets of large DEMs.
