# CLI reference: lsdtt-synthetic-networks

Generates synthetic Optimal Channel Networks (OCNs) with tunable branching statistics by simulated-annealing edge flips that minimise the energy `E = sum_i A_i^gamma`. The cooling schedule sweeps the network from random / Shreve-like (high temperature) to strongly aggregated (temperature → 0). Elevations are then painted on the resulting spanning tree as a chi-based fluvial steady state, using the same core solver as `lsdtt-lem`.

This is a **generator**: no input raster is required. The only file it can read is an optional polygon used to define the domain outline.

Design references: Carraro et al. (2020) OCNet; Rodriguez-Iturbe & Rinaldo (1997).

## Usage

```
lsdtt-synthetic-networks <path> <parameter_file>
lsdtt-synthetic-networks <parameter_file>    # path defaults to ./
lsdtt-synthetic-networks -h                  # print help and exit
lsdtt-synthetic-networks -v                  # print version and exit
lsdtt-synthetic-networks -c                  # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

| Stage | Description |
|---|---|
| 1 | Validate base level and georeference parameters |
| 2 | Build the domain mask (rectangle, island, or imported polygon) |
| 3 | Choose outlets (`outlet_mode`, `auto` resolved from the domain) |
| 4 | Anneal the spanning tree, minimising `E = sum(A^gamma)` |
| 5 | Synthesise elevations as a chi-based fluvial steady state |
| 6 | Write the DEM and drainage-area rasters |
| 7 | Route the DEM through D8 and report round-trip graph fidelity (`verify_routing`) |
| 8 | Write `_ocn_summary.csv` and the README |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory searched for `polygon_file`. No input raster is read. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | parameter file stem | Output filename stem. |
| `write_output_format` | string | `cog` | Output raster format: `cog`, `gtiff`, or `same_as_input`. |

---

## Lattice geometry

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_rows` | int | `128` | Number of lattice rows. Every pixel is a node of the spanning tree. |
| `n_cols` | int | `128` | Number of lattice columns. |
| `cell_size` | float | `1.0` | Pixel size in map units (metres). |
| `outlet_col` | int | `-1` | Column of the single outlet on the south (bottom) edge. `-1` places it at the centre column. Only used when the resolved outlet mode is `single`. |

Note that with `domain_type: polygon` and `polygon_sizing: true_area`, `n_rows` and `n_cols` are **computed from the polygon extent and `cell_size`**, not taken from the parameter file.

---

## Domain shape

| Parameter | Type | Default | Description |
|---|---|---|---|
| `domain_type` | string | `rectangle` | `rectangle`: the full `n_rows` × `n_cols` grid. `efd_island`: a smooth random island from an Elliptical Fourier descriptor. `noise_island`: a radially-tapered red-noise island (fractal coast). `polygon`: a user polygon vector scaled to the grid (requires `polygon_file`). An unrecognised value warns and falls back to `rectangle`. |
| `outlet_mode` | string | `auto` | Where drainage leaves the domain. `auto`: islands use `all_coast`, the rectangle uses `single`. `single`: one outlet, so the whole domain drains through one network (most sensitive to boundary shape). `all_coast`: every coastal land pixel is a base level. |
| `island_fill_fraction` | float | `0.7` | Approximate fraction of the grid occupied by land. Span (`efd_island`) or area fraction (`noise_island`). |
| `efd_n_harmonics` | int | `6` | Number of Elliptical Fourier harmonics for `efd_island`. 1 = ellipse; more harmonics add lobes and bays. Try 4–12. |
| `efd_roughness` | float | `1.5` | Harmonic amplitude decay exponent for `efd_island`. Higher = smoother coast; lower = rougher. |
| `noise_blur_passes` | int | `6` | Smoothing passes for the `noise_island` coastline. More passes = smoother; fewer = rougher / fractal. |
| `polygon_file` | string | `NULL` | Vector polygon file defining the domain when `domain_type: polygon`. The first polygon's exterior ring is used as the coastline. Resolved against `read_path` if not an absolute path. |
| `polygon_sizing` | string | `scaled` | How an imported polygon is fitted to a grid. `scaled`: scale the polygon to fill the `n_rows` × `n_cols` grid (shape preserved, size set by the grid). `true_area`: derive the grid from the polygon's real extent and `cell_size`, preserving real area and the polygon's CRS. |

---

## Energy and branching

| Parameter | Type | Default | Description |
|---|---|---|---|
| `energy_exponent_gamma` | float | `0.5` | Energy exponent γ in `E = sum(A^gamma)` (OCNet `expEnergy`). **Must be in the open interval (0,1)**; a value outside it is a fatal error from the generator. Classical value 0.5. Smaller (more concave) γ rewards aggregation more strongly — larger Tokunaga *a*, *c* and Horton ratios; γ → 1 removes the aggregation incentive. |

---

## Annealing schedule

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_iter` | int | `0` | Number of annealing iterations (proposed flips). `0` selects the OCNet default `40 * n_rows * n_cols`. |
| `cooling_rate` | float | `1.0` | Cooling rate, OCNet convention: **larger = faster cooling**. After the no-cooling phase, `T_k = T0 * exp(-cooling_rate * k / N)` with `N = n_rows * n_cols`. Faster cooling freezes a more random network; slower cooling reaches lower-energy, more aggregated states. |
| `initial_no_cooling_fraction` | float | `0.0` | Fraction of iterations held at the initial temperature before cooling begins. |
| `initial_temperature` | float | `0.0` | Manual initial temperature. `<= 0` means auto, in which case `T0 = temperature_scale * E0` with `E0` the initial-tree energy. |
| `temperature_scale` | float | `1.0` | Multiplier on the auto initial temperature. Continuous knob for sweeping branching statistics: raise toward random, lower toward 0 for the optimised low-temperature limit (greedy descent at 0). `temperature_scale = 1` reproduces OCNet's native `T0 = E0`. |
| `random_seed` | int | `1` | Seed for the annealing RNG **and** for the random island generators. Fixing it makes a run reproducible. |

> **`n_iter` matters more than `temperature_scale` for routing fidelity.** Truncating the anneal costs more graph fidelity than the whole temperature range does at fixed `n_iter`, and a short, hot run may never cool at all. If you are sweeping `temperature_scale`, hold `n_iter` at the default or scale it with the grid. See `verify_routing` below and the measured table in `docs/testing_and_benchmarking.md` §8b.1.

---

## Elevation synthesis

Elevations are a chi-based fluvial steady state on the annealed tree, identical to `lsdtt-lem`'s steady-state solver: `chi[i] = chi[receiver] + dx_eff * A^(-m/n)` and `z_i = baselevel_elevation + f * chi_i`. **Fluvial relief only** — there is no hillslope component and no diffusive smoothing, so ridges are chi discontinuities between adjacent trees.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `steady_state_mode` | string | `relief` | `relief`: total relief is set to `target_relief` and K is back-calculated. `K`: erodibility is given directly and relief follows from `uplift_rate` / `erodibility`. |
| `target_relief` | float | `1000` | Target fluvial relief, `max(z) - baselevel`, in metres. Used when `steady_state_mode: relief`. |
| `uplift_rate` | float | `1e-4` | Uplift rate U (m/yr) for the steady-state stream-power solution. In relief mode only the resulting steepness matters. |
| `erodibility` | float | `1e-3` | Erodibility K. Used when `steady_state_mode: K`; ignored in relief mode. |
| `m_exponent` | float | `0.45` | Drainage-area exponent m in the stream-power / chi law. Concavity is m/n. Set the same m and n as any LEM run you compare against. |
| `n_exponent` | float | `1.0` | Slope exponent n. Concavity is m/n. |
| `baselevel_elevation` | float | `0.0` | Elevation (m) assigned to the outlet(s). Chi is integrated upstream from this value. |

### Base-level validation

`baselevel_elevation` is checked before anything is generated, and a failure is fatal:

| Condition | Behaviour |
|---|---|
| not finite | `FATAL: baselevel_elevation must be finite.` |
| `<= -9000` | `FATAL: baselevel_elevation must be > -9000.` |

The second check exists because the raster nodata value is `-9999.0`: a base level at or below that would silently turn the outlets into nodata rather than producing a visibly wrong DEM.

For cosmogenic work note that the base-level elevation is what sets the atmospheric pressure over the generated landscape, and therefore the production rate — it is a physically meaningful parameter, not just an offset.

---

## Georeference

By default the generated grid is placed at the historical default synthetic anchor. These three parameters place it somewhere specific instead, which matters whenever the latitude is physically meaningful (again, cosmogenic production).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `centre_latitude` | float | `-9999` | Latitude of the synthetic grid centre, decimal degrees. `-9999` keeps the default synthetic anchor. |
| `centre_longitude` | float | `-9999` | Longitude of the synthetic grid centre, decimal degrees. `-9999` keeps the default synthetic anchor. |
| `projection_epsg` | int | `0` | Projected CRS for the output. `0` derives the WGS84 UTM zone from `centre_latitude` / `centre_longitude`. Non-zero uses the supplied EPSG code directly. |

Validation:

| Condition | Behaviour |
|---|---|
| only one of `centre_latitude` / `centre_longitude` set | `FATAL: centre_latitude and centre_longitude must be set together.` |
| either is set but not finite | `FATAL: centre_latitude and centre_longitude must be finite.` |
| `domain_type: polygon` **and** `polygon_sizing: true_area` **and** any of `centre_latitude`, `centre_longitude`, `projection_epsg` set | **Warning**, and the anchor parameters are ignored. The grid takes the polygon's own CRS and extent. |

The `true_area` case is a genuine conflict rather than an oversight: `true_area` sizing derives the grid, the cell origin and the CRS from the polygon file, so there is nothing left for the anchor to set. The run continues with the polygon's georeference and says so on stderr.

---

## Outputs

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_dem_raster` | bool | `true` | Write the synthetic elevation raster, `<write_fname>_ocn_dem.tif`. |
| `write_drainage_area_raster` | bool | `true` | Write the OCN drainage-area raster, `<write_fname>_ocn_drainage_cells.tif`. Units are **pixel counts**, not m². |
| `verify_routing` | bool | `true` | Route the synthetic DEM through D8 and report the round-trip graph fidelity. |

### Output files

| File | Description |
|---|---|
| `<write_fname>_ocn_dem.tif` | Synthetic elevation (float32). |
| `<write_fname>_ocn_drainage_cells.tif` | OCN drainage area in pixel counts. |
| `<write_fname>_ocn_summary.csv` | One row: `n_rows, n_cols, cell_size, gamma, n_iter, cooling_rate, temperature_scale, initial_no_cooling_fraction, seed, initial_energy, final_energy, routing_fidelity`. |
| `<write_fname>_README.*` | Standard parameter README pair. |

**No flow-direction, receiver or network-topology product is written.** Downstream tools must re-route the DEM — which is exactly why `verify_routing` exists.

### Routing fidelity: read the number, do not assume it

`verify_routing` computes the fraction of interior, non-outlet land pixels whose D8 receiver on the synthesised DEM equals their OCN receiver, prints it as `D8 round-trip fidelity (interior)` and records it as `routing_fidelity` in the summary CSV. A failure inside the verification step warns and leaves `routing_fidelity` at `-1`; it does not stop the run.

This number is a **gate, not a diagnostic**. Fidelity is high for relaxed (low-temperature) networks and falls for random (high-temperature) ones, but it does not reach 1 at any temperature or anneal length: it plateaus near 0.935 even at more than four times the default `n_iter`, so several per cent of interior pixels always re-route. If your workflow re-routes the DEM — and every lsdtt3 tool downstream of this one does — then the network you analyse is not exactly the network you generated. The measured fidelity-versus-temperature-versus-`n_iter` table, and what to do about it, are in `docs/testing_and_benchmarking.md` §8b.1.

---

## Example parameter files

### A relaxed rectangular OCN

```
write_path:  /data/output/
write_fname: ocn_cold

n_rows:            128
n_cols:            128
cell_size:         50

domain_type:       rectangle
outlet_mode:       single

energy_exponent_gamma: 0.5
temperature_scale:     0.0
random_seed:           42

steady_state_mode: relief
target_relief:     1000
m_exponent:        0.45
n_exponent:        1.0
baselevel_elevation: 0

verify_routing:      true
write_output_format: gtiff
```

### An island draining to its whole coastline, placed at a real latitude

```
write_path:  /data/output/
write_fname: ocn_island

n_rows:               256
n_cols:               256
cell_size:            30

domain_type:          noise_island
outlet_mode:          all_coast
island_fill_fraction: 0.6
noise_blur_passes:    6

energy_exponent_gamma: 0.5
temperature_scale:     0.1
random_seed:           7

centre_latitude:   -41.5
centre_longitude:  174.0
projection_epsg:   0

steady_state_mode:   relief
target_relief:       1500
baselevel_elevation: 20

write_output_format: cog
```

### A domain imported from a polygon at its true area

`n_rows`, `n_cols` and the CRS all come from the polygon here; do not also set the anchor parameters, or they will be ignored with a warning.

```
read_path:   /data/vectors/
write_path:  /data/output/
write_fname: ocn_catchment

domain_type:    polygon
polygon_file:   catchment_outline.gpkg
polygon_sizing: true_area
cell_size:      30

energy_exponent_gamma: 0.5
temperature_scale:     0.0
random_seed:           1

steady_state_mode: relief
target_relief:     800

write_output_format: gtiff
```
