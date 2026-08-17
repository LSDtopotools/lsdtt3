# CLI reference: lsdtt-raster-synthetic-tests

Generates synthetic raster test surfaces for algorithm validation and benchmarking. Creates idealized geometries including planes, cones, ridges, mountain ranges, and fractal surfaces with optional noise.

## Usage

```
lsdtt-raster-synthetic-tests  <path> <parameter_file>
lsdtt-raster-synthetic-tests  <parameter_file>        # path defaults to ./
lsdtt-raster-synthetic-tests  -h                      # print help and exit
lsdtt-raster-synthetic-tests  -v                      # print version and exit
lsdtt-raster-synthetic-tests  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

Keys use `snake_case`. All keys and their defaults are listed below.

---

## Pipeline overview

Steps run in this order:

| Step | Description |
|---|---|
| 1. Grid creation | Initialise raster grid with specified dimensions and CRS |
| 2. Surface generation | Create base synthetic surface using selected geometry |
| 3. Noise addition *(optional)* | Add white noise, red noise, or fractal perturbations |
| 4. Scaling | Apply amplitude and vertical offset |
| 5. Output | Write synthetic DEM with metadata |

---

## Output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_path` | string | current directory | Directory for output files. |
| `write_fname` | string | `synthetic_dem` | Output filename stem. |

---

## Grid specification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_rows` | int | `512` | Number of rows in output raster. |
| `n_cols` | int | `512` | Number of columns in output raster. |
| `cell_size_metres` | float | `10.0` | Pixel size (metres). |
| `set_default_synthetic_georeference` | bool | `true` | Use default georeference (Teanu Island, UTM 58S). If `false`, raster has no CRS. |

---

## Surface type

| Parameter | Type | Default | Description |
|---|---|---|---|
| `surface_type` | string | `plane` | Type of surface to generate: `plane`, `cone`, `ridge`, `range`, `diamond_square`, `perlin_noise`. |

---

## Surface type: Plane

Simple flat or tilted plane.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `plane_slope_radians` | float | `0.0` | Slope angle (radians). Set to `0` for flat plane. |
| `plane_aspect_radians` | float | `0.0` | Aspect direction (radians, 0=North, π/2=East). |
| `plane_base_elevation_metres` | float | `1000.0` | Elevation at grid centre (metres). |

---

## Surface type: Cone

Conical mountain.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `cone_base_elevation_metres` | float | `1000.0` | Base elevation (metres). |
| `cone_peak_elevation_metres` | float | `2000.0` | Peak elevation (metres). |
| `cone_centre_x_fraction` | float | `0.5` | Peak location, x-fraction of grid (0–1). |
| `cone_centre_y_fraction` | float | `0.5` | Peak location, y-fraction of grid (0–1). |

---

## Surface type: Ridge

Linear mountain ridge.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ridge_base_elevation_metres` | float | `1000.0` | Base elevation (metres). |
| `ridge_peak_elevation_metres` | float | `2000.0` | Peak elevation (metres). |
| `ridge_orientation_radians` | float | `0.0` | Ridge strike direction (radians). |

---

## Surface type: Mountain range

Multiple peaks with valley floors.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `range_base_elevation_metres` | float | `1000.0` | Valley floor elevation (metres). |
| `range_peak_elevation_metres` | float | `2000.0` | Peak elevation (metres). |
| `range_n_peaks` | int | `5` | Number of peaks. |
| `range_peak_spacing_fraction` | float | `0.2` | Spacing between peaks (fraction of grid diagonal). |

---

## Surface type: Diamond-Square (fractal)

Fractal surface generated via diamond-square algorithm.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `diamond_square_base_elevation_metres` | float | `1000.0` | Base elevation (metres). |
| `diamond_square_roughness` | float | `0.7` | Roughness parameter (0–1). Higher = rougher surface. |
| `diamond_square_max_displacement_metres` | float | `100.0` | Maximum elevation variation (metres). |

---

## Surface type: Perlin noise

Procedural noise-based surface.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `perlin_base_elevation_metres` | float | `1000.0` | Base elevation (metres). |
| `perlin_scale_metres` | float | `500.0` | Feature scale (metres). Larger = smoother features. |
| `perlin_amplitude_metres` | float | `100.0` | Amplitude of elevation variation (metres). |
| `perlin_octaves` | int | `4` | Number of noise octaves (more = more detail, slower). |

---

## Noise addition

| Parameter | Type | Default | Description |
|---|---|---|---|
| `add_white_noise` | bool | `false` | Add white noise. |
| `white_noise_amplitude_metres` | float | `10.0` | Amplitude of white noise (metres). |
| `add_red_noise` | bool | `false` | Add red (correlated) noise. |
| `red_noise_amplitude_metres` | float | `20.0` | Amplitude of red noise (metres). |
| `red_noise_correlation_length_pixels` | int | `50` | Correlation length (pixels) for red noise. |

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format: `cog` (preferred), `gtiff`, `png`. |

---

## Output files

- `<write_fname>_dem.tif` — main synthetic DEM
- `<write_fname>_dem.png` *(optional)* — PNG hillshade preview (if `write_output_format: png`)
