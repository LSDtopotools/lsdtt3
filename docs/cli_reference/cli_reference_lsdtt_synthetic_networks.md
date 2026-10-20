# CLI reference: lsdtt-synthetic-networks

Generates synthetic Optimal Channel Networks (OCNs) with tunable branching statistics by simulated-annealing edge flips that minimise the energy `E = sum_i A_i^gamma * L_i` (length-weighted by default; see `energy_length_weighted`). The cooling schedule sweeps the network from random / Shreve-like (high temperature) to strongly aggregated (temperature → 0). Elevations are then painted on the resulting spanning tree as a chi-based fluvial steady state, using the same core solver as `lsdtt-lem`.

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
| 4 | Anneal the spanning tree, minimising `E = sum(A^gamma * L)` |
| 5 | Synthesise elevations as a chi-based fluvial steady state |
| 6 | Write the DEM and drainage-area rasters |
| 7 | Route the DEM through D8 and report round-trip graph fidelity (`verify_routing`) |
| 8 | Write `_ocn_summary.csv` and the README |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory searched for `polygon_fname`. No input raster is read. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | parameter file stem | Output filename stem. |
| `write_output_format` | string | `cog` | Output raster format: `cog`, `gtiff`, or `same_as_input`. |

---

## Lattice geometry

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_rows` | int | `128` | Number of lattice rows. Every pixel is a node of the spanning tree. |
| `n_cols` | int | `128` | Number of lattice columns. |
| `cell_size_metres` | float | `1.0` | Pixel size in metres. Former `cell_size` is an unchanged-value alias. |
| `outlet_col` | int | `-1` | Column of the single outlet on the south (bottom) edge. `-1` places it at the centre column. Only used when the resolved outlet mode is `single`. |

Note that with `domain_type: polygon` and `polygon_sizing: true_area`, `n_rows` and `n_cols` are **computed from the polygon extent and `cell_size_metres`**, not taken from the parameter file.

---

## Domain shape

| Parameter | Type | Default | Description |
|---|---|---|---|
| `domain_type` | string | `rectangle` | `rectangle`: the full `n_rows` × `n_cols` grid. `efd_island`: a smooth random island from an Elliptical Fourier descriptor. `noise_island`: a radially-tapered red-noise island (fractal coast). `polygon`: a user polygon vector scaled to the grid (requires `polygon_fname`). An unrecognised value warns and falls back to `rectangle`. |
| `outlet_mode` | string | `auto` | Where drainage leaves the domain. `auto`: islands use `all_coast`, the rectangle uses `single`. `single`: one outlet, so the whole domain drains through one network (most sensitive to boundary shape). `all_coast`: every coastal land pixel is a base level. |
| `island_fill_fraction` | float | `0.7` | Approximate fraction of the grid occupied by land. Span (`efd_island`) or area fraction (`noise_island`). |
| `efd_n_harmonics` | int | `6` | Number of Elliptical Fourier harmonics for `efd_island`. 1 = ellipse; more harmonics add lobes and bays. Try 4–12. |
| `efd_roughness` | float | `1.5` | Harmonic amplitude decay exponent for `efd_island`. Higher = smoother coast; lower = rougher. |
| `coastline_blur_passes` | int | `6` | Smoothing passes for the `noise_island` coastline. Former `noise_blur_passes` is an unchanged-value alias. |
| `polygon_fname` | string | `NULL` | Vector polygon file defining the domain when `domain_type: polygon`. The first polygon's exterior ring is used as the coastline. Resolved against `read_path` if not an absolute path. |
| `polygon_sizing` | string | `scaled` | How an imported polygon is fitted to a grid. `scaled`: scale the polygon to fill the `n_rows` × `n_cols` grid (shape preserved, size set by the grid). `true_area`: derive the grid from the polygon's real extent and `cell_size_metres`, preserving real area and the polygon's CRS. |

---

## Energy and branching

| Parameter | Type | Default | Description |
|---|---|---|---|
| `energy_length_weighted` | bool | `true` | Minimise `E = sum(A^gamma * L)` rather than `sum(A^gamma)`, where `L` is the link length (`cell_size_metres`, or `cell_size_metres*sqrt(2)` for a diagonal). **This is the default and you should normally leave it on.** Link energy expenditure is `E_i ~ L_i Q_i S_i` (Strong & Mudd 2022 PNAS eq. 4; Goren & Shelef 2024 *Earth Surf. Dynam.* 12:1347–1369 eq. 5). On a 4-connected lattice `L` factors out, which is why the literature writes `sum(A^gamma)` — but on this 8-connected lattice a diagonal link is 1.414x longer for the same cost, so lateral movement is free and the trunk channels run diagonally regardless of where base level is. Setting this also selects a length-aware (Dijkstra) seed tree in place of the Chebyshev breadth-first one. Set `false` to reproduce OCNet (Carraro et al. 2020) bit-for-bit — defensible if you want the flow *graph* rather than a geometric network. |
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
| `random_seed` | int | `33` | Seed for the annealing RNG **and** for the random island generators. Fixing it makes a run reproducible. |

> **`n_iter` matters more than `temperature_scale` for routing fidelity.** Truncating the anneal costs more graph fidelity than the whole temperature range does at fixed `n_iter`, and a short, hot run may never cool at all. If you are sweeping `temperature_scale`, hold `n_iter` at the default or scale it with the grid. See `verify_routing` below and the measured table in `docs/developer_notes/testing_and_benchmarking.md` §8b.1.

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
| `write_drainage_area_raster` | bool | `true` | Write the OCN drainage-area raster. Enabled because area is a primary OCN result; analysis CLIs treat the same output as an optional diagnostic. Units are **pixel counts**, not m². |
| `verify_routing` | bool | `true` | Route the synthetic DEM through D8 and report the round-trip graph fidelity. |

### Output files

| File | Description |
|---|---|
| `<write_fname>_ocn_dem.tif` | Synthetic elevation (float32). |
| `<write_fname>_ocn_drainage_cells.tif` | OCN drainage area in pixel counts. |
| `<write_fname>_ocn_summary.csv` | One row: `n_rows, n_cols, cell_size, gamma, n_iter, cooling_rate, temperature_scale, initial_no_cooling_fraction, seed, initial_energy, final_energy, routing_fidelity`, followed by the equilibration columns `equilibrated, eq_cycles_run, eq_converged, eq_limit_cycle_period, eq_divide_step_p90_m, eq_ocn_agreement`. When `run_equilibration` is false these read `0,0,0,0,-1,-1`. |
| `<write_fname>_ocn_equilibrated_dem.tif` | Equilibrated elevation, when `run_equilibration: true`. |
| `<write_fname>_ocn_snap_convergence.csv` | Per-cycle equilibration diagnostics, when `run_equilibration: true`. Same schema as the `lsdtt-lem` file of that name. |
| `<write_fname>_README.*` | Standard parameter README pair. |

**No flow-direction, receiver or network-topology product is written.** Downstream tools must re-route the DEM — which is exactly why `verify_routing` exists.

### Routing fidelity: read the number, do not assume it

`verify_routing` computes the fraction of interior, non-outlet land pixels whose D8 receiver on the synthesised DEM equals their OCN receiver, prints it as `D8 round-trip fidelity (interior)` and records it as `routing_fidelity` in the summary CSV. A failure inside the verification step warns and leaves `routing_fidelity` at `-1`; it does not stop the run.

This number is a **gate, not a diagnostic**. Fidelity is high for relaxed (low-temperature) networks and falls for random (high-temperature) ones, but it does not reach 1 at any temperature or anneal length: it plateaus near 0.935 even at more than four times the default `n_iter`, so several per cent of interior pixels always re-route. If your workflow re-routes the DEM — and every lsdtt3 tool downstream of this one does — then the network you analyse is not exactly the network you generated. The measured fidelity-versus-temperature-versus-`n_iter` table, and what to do about it, are in `docs/developer_notes/testing_and_benchmarking.md` §8b.1.

---

## Example parameter files

### A relaxed rectangular OCN

```
write_path:  /data/output/
write_fname: ocn_cold

n_rows:            128
n_cols:            128
cell_size_metres:  50

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
cell_size_metres:     30

domain_type:          noise_island
outlet_mode:          all_coast
island_fill_fraction: 0.6
coastline_blur_passes: 6

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
polygon_fname:  catchment_outline.gpkg
polygon_sizing: true_area
cell_size_metres: 30

energy_exponent_gamma: 0.5
temperature_scale:     0.0
random_seed:           1

steady_state_mode: relief
target_relief:     800

write_output_format: gtiff
```

---

## Equilibration: OCN → equilibrate in one invocation

Because the chi DEM never exactly carries the OCN graph (see *Routing fidelity*
above), the delivered surface is not self-consistent: re-route it and several
per cent of interior pixels pick a different receiver. `run_equilibration: true`
runs the landscape equilibration cycle — pit-fill, re-route, re-snap, repeat —
in the same invocation, so the DEM you hand downstream is one whose own D8
routing is stable.

Parameter names are taken **verbatim from `lsdtt-lem`**, so anything you know
there works here. See `cli_reference_lsdtt_lem.md` for the full descriptions;
this page lists what is available.

| Group | Parameters |
|---|---|
| switch | `run_equilibration` |
| cycles | `steady_state_n_cycles` (a **maximum** when convergence checking is on) |
| convergence | `steady_state_check_convergence`, `steady_state_tolerance_dz_over_relief`, `steady_state_tolerance_dz_metres`, `steady_state_tolerance_receiver_fraction`, `steady_state_convergence_patience` |
| divide diagnostics | `divide_merge_search_steps`, `divide_use_d8_diagonals`, `divide_gradient_target` |
| routing | `boundary_conditions`, `fill_min_slope`, `divide_area_by_dx`, `n_threads` |
| hillslopes | `run_hillslope_steady_state`, `hillslope_ss_model`, `incision_threshold_area_m2`, `hillslope_diffusivity`, `critical_hillslope_slope`, `enable_nonlinear_hillslopes`, `hillslope_ss_tolerance`, `hillslope_ss_max_iterations`, `hillslope_sor_omega`, `hillslope_every_cycle`, `channel_hysteresis_factor` |
| outputs | `write_equilibrated_dem_raster`, `write_snap_convergence_csv` |

The fluvial parameters the cycle snaps with (`uplift_rate`, `erodibility`,
`m_exponent`, `n_exponent`, `steady_state_mode`, `target_relief`,
`baselevel_elevation`) are the **same ones used to synthesise the OCN
elevations** — there is no second set to keep in step.

The un-equilibrated DEM is still written as `_ocn_dem.tif`, so the before/after
pair is always on disk.

### `eq_ocn_agreement` is drift, not an improved fidelity

After equilibration the CLI re-runs the same receiver comparison against the
**original OCN graph** and reports it as
`Equilibrated DEM vs the ORIGINAL OCN graph (drift, not fidelity)`.

Read it as drift. Equilibration makes the DEM self-consistent with its *own*
routing; it does not pull the surface back towards the OCN receivers, and
divides migrate while it works. **The number goes down, and that is expected.**
Measured on a 48×48, γ = 0.5, seed 42, `temperature_scale` 0 network:

| run | agreement with the OCN graph |
|---|---:|
| before equilibration (`routing_fidelity`) | 0.916 |
| after fluvial-only equilibration | 0.895 |
| after equilibration with Roering hillslopes | 0.290 |

The last row is not a bug. With a channel threshold that leaves only tens of
channel pixels, the hillslope solve reshapes almost the whole grid, so almost
nothing routes as the OCN did. If your study rests on the *branching statistics*
of the generated network, measure them after equilibration, not before — this
column is how you find out how much of the network survived.

### Two failure modes worth knowing

**A closed outlet edge is refused.** A `single` outlet OCN drains through the
south edge, so `boundary_conditions: nnnn` leaves no path to base level. The CLI
rejects this up front with a FATAL naming the setting, because the router on its
own can only report "could not assign a receiver at row R, col C" — true, but it
does not tell you the cause is a parameter you set.

**Watch the pre-flight line.** Before the first cycle the equilibration reports
predicted hillslope length, relief and maximum gradient. If it warns that
`predicted_hillslope_relief/computed_relief > 0.5`, the hillslope solve is being
asked for more relief than the landscape has, and it will usually diverge and
leave sinks that break the next cycle's routing. Lower `hillslope_diffusivity`
or raise `incision_threshold_area_m2` rather than pushing on.

### Example: generate and equilibrate

```
read_path:  ./
write_path: /data/output/
write_fname: ocn01

n_rows: 48
n_cols: 48
cell_size_metres: 100
energy_exponent_gamma: 0.5
random_seed: 42
temperature_scale: 0.0

steady_state_mode: relief
target_relief:     300
uplift_rate:       0.0001

run_equilibration:                  true
steady_state_n_cycles:              20
steady_state_check_convergence:     true
steady_state_convergence_patience:  2

run_hillslope_steady_state:  true
hillslope_ss_model:          roering
incision_threshold_area_m2:  20000
hillslope_diffusivity:       2.0
critical_hillslope_slope:    0.8
hillslope_sor_omega:         1.3

write_equilibrated_dem_raster: true
write_snap_convergence_csv:    true
write_output_format:           gtiff
```


---

## Why the energy is length-weighted, and what happens if it is not

Measured on a 400x200 grid at 30 m draining north and south (periodic sides,
gamma 0.5, `temperature_scale: 0`, fixed `n_iter`). Flow direction **weighted by
drainage area over channels larger than 2 km²** — that is, where the water
actually goes, not where the average headwater link points:

| | diagonal | pure N/S |
|---|---:|---:|
| `energy_length_weighted: false` (OCNet parity) | **0.937** | 0.062 |
| `energy_length_weighted: true` (default) | **0.056** | **0.944** |

With the length term dropped, 94% of the flow leaves its cell diagonally on a
strip whose only base levels are north and south. That is a lattice artefact,
not a landscape: real drainage follows steepest descent, where the diagonal
neighbour is sqrt(2) farther and its gradient is correspondingly lower, so
diagonal routing is penalised. The OCN energy only reproduces that penalty if
`L` is in it.

Two practical notes:

- **Use an area-weighted statistic to check this yourself.** Unweighted link
  counts are dominated by headwaters, which point everywhere in any network:
  the count-based diagonal fraction only moves 0.45 to 0.22 across the same
  change that moves the area-weighted one from 0.94 to 0.06.
- **Round-trip `routing_fidelity` is lower under the default** (0.837 vs 0.916
  on the 48x48 reference). That is expected: the default drives the graph toward
  straight cardinal trunks while the chi surface is re-routed by D8 steepest
  descent, which still favours diagonals, so the two disagree more. See
  *Routing fidelity* above and `docs/developer_notes/testing_and_benchmarking.md` 8b.1b.

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

### Synthetic domain

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_rows` | `int` | `128` | `param` | - | Number of rows in the synthetic model domain. | Used only when the domain is generated rather than read from a DEM. |
| `n_cols` | `int` | `128` | `param` | - | Number of columns in the synthetic model domain. | Used only when the domain is generated rather than read from a DEM. |
| `cell_size_metres` | `float` | `1.0` | `param` | - | Cell size of the synthetic model domain in metres. | Sets the map scale of the generated surface. Anything expressed as a pixel count (for example threshold_contributing_pixels) corresponds to a different ground area as this changes. |
| `random_seed` | `int` | `33` | `param` | - | Seed for every stochastic generator used to build the domain. | A fixed seed makes a run reproducible. Change it to draw a different realisation of the same statistical surface. |

### Lattice geometry

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `outlet_spacing` | `int` | `1` | `param` | - | Column spacing of baselevel outlets, for outlet_mode north_south. | This is the knob that sets BASIN WIDTH. At 1 every boundary pixel is its own baselevel, so the tile is partitioned into as many basins as there are boundary cells and the basin scale is imposed by the grid rather than emerging from the optimisation - annealing temperature cannot change it, because temperature controls how each basin BRANCHES and not how many there are. Larger values force the drainage to organise into wider basins: it is the OCN analogue of the diamond-square feature wavelength, and the continuous version of the single-outlet configuration used in much of the OCN literature (outlet_mode single). Ignored in the other outlet modes. |
| `outlet_col` | `int` | `-1` | `param` | - | Column of the single outlet on the south (bottom) edge. | -1 places the outlet at the centre column. Only used for outlet_mode single. |

### Domain (rectangle vs irregular island)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `domain_type` | `string` | `rectangle` | `param` | - | Shape of the domain the network is grown on. | rectangle: the full n_rows x n_cols grid. efd_island: a smooth random island from an Elliptical Fourier descriptor. noise_island: a radially-tapered red-noise island (fractal coast). polygon: a user polygon vector scaled to the grid (set polygon_fname). |
| `outlet_mode` | `string` | `auto` | `param` | - | Where drainage leaves the domain. | auto: islands use all_coast, the rectangle uses single. single: one outlet (whole domain drains through one network; most sensitive to boundary shape). all_coast: every coastal land pixel is a base level (an island draining to its whole coastline). north_south: every land pixel on the north and south edges, and nothing else -- the classic ridge-and-valley strip, and the mode to use with periodic sides (boundary_conditions: opop), where the main divide is free to sit anywhere between the two boundaries. |
| `island_fill_fraction` | `float` | `0.7` | `param` | - | Approximate fraction of the grid occupied by land (island generators). | Span (efd) or area fraction (noise) of the island within the grid. |
| `efd_n_harmonics` | `int` | `6` | `param` | - | Number of Elliptical Fourier harmonics for efd_island. | 1 = ellipse; more harmonics add lobes and bays. Try 4-12. |
| `efd_roughness` | `float` | `1.5` | `param` | - | Harmonic amplitude decay exponent for efd_island. | Higher = smoother coast (high harmonics suppressed); lower = rougher. |
| `coastline_blur_passes` | `int` | `6` | `param` | - | Smoothing passes for the noise_island coastline. | More passes = smoother, simpler coastline; fewer = rougher/fractal. |
| `polygon_fname` | `string` | `NULL` | `param` | - | Vector polygon file defining the domain (domain_type = polygon). | The first polygon's exterior ring is used as the island coastline. |
| `polygon_sizing` | `string` | `scaled` | `param` | - | How an imported polygon is fitted to a grid. | scaled: scale the polygon to fill the n_rows x n_cols grid (shape preserved, size set by the grid). true_area: derive the grid from the polygon's real extent and cell_size_metres, preserving real area and the polygon's CRS (n_rows/n_cols are computed, not taken from parameters). |

### Energy and branching

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `energy_length_weighted` | `bool` | `true` | `param` | - | Weight each link's energy by its length: E = sum(A^gamma * L). | This is the physically correct form and the default. Link energy expenditure is E_i ~ L_i Q_i S_i (Strong and Mudd 2022 PNAS eq. 4; Goren and Shelef 2024 Earth Surf. Dynam. 12:1347-1369 eq. 5), which reduces to A^gamma * L. On a 4-connected lattice L factors out, which is why the literature writes sum(A^gamma) -- but on the 8-connected lattice used here a diagonal link is sqrt(2) longer and dropping L under-charges it, so the annealer buys free lateral movement and the main channels come out running diagonally. It also selects a length-aware (Dijkstra) seed tree instead of the Chebyshev breadth-first one. Measured on a 400x200 strip draining north-south, area-weighted over channels &gt;2 km2: diagonal flow falls from 0.937 to 0.056 and north-south rises from 0.062 to 0.944. Set false to reproduce OCNet (Carraro et al. 2020) bit-for-bit; that is defensible if you want the flow GRAPH rather than a geometric network. |

### Energy / branching

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `energy_exponent_gamma` | `float` | `0.5` | `param` | - | Energy exponent gamma in E = sum(A^gamma) (OCNet expEnergy). | Must be in (0,1). Classical value 0.5. Smaller (more concave) gamma rewards aggregation more strongly -&gt; larger Tokunaga a,c and Horton ratios; gamma -&gt; 1 removes the aggregation incentive. |

### Annealing schedule

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_iter` | `int` | `0` | `param` | - | Number of annealing iterations (proposed flips). | 0 selects the OCNet default 40*n_rows*n_cols. |
| `cooling_rate` | `float` | `1.0` | `param` | - | Cooling rate (OCNet convention: larger = FASTER cooling). | Temperature after the no-cooling phase decays as T_k = T0*exp(-cooling_rate*k/N), N = n_rows*n_cols. Faster cooling freezes a more random network; slower cooling reaches lower-energy, more aggregated states. |
| `initial_no_cooling_fraction` | `float` | `0.0` | `param` | - | Fraction of iterations held at the initial temperature before cooling. | Lets the chain equilibrate at high T (more random) before cooling. |
| `initial_temperature` | `float` | `0.0` | `param` | - | Manual initial temperature; &lt;= 0 means auto. | When &lt;= 0, T0 = temperature_scale * E0 (E0 = initial-tree energy); temperature_scale = 1 reproduces OCNet's native T0 = E0. |
| `temperature_scale` | `float` | `1.0` | `param` | - | Multiplier on the auto initial temperature (T0 = scale * E0). | Continuous knob for sweeping Tokunaga statistics: raise toward random (c~2), lower toward 0 for the optimised low-T limit (greedy at 0). |

### Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `steady_state_mode` | `string` | `relief` | `param` | - | How the fluvial relief of the synthetic DEM is set. | relief: total relief is set to target_relief (K is back-calculated). K: erodibility is given directly and relief follows from uplift_rate/erodibility. Identical to lsdtt-lem's steady-state solver. |
| `target_relief` | `float` | `1000` | `param` | - | Target fluvial relief max(z)-baselevel in metres (relief mode). | Used when steady_state_mode = relief. |
| `uplift_rate` | `float` | `1e-4` | `param` | - | Uplift rate U (m/yr) for the steady-state stream-power solution. | U = K*A^m*S^n. In relief mode only the resulting steepness matters. |
| `erodibility` | `float` | `1e-3` | `param` | - | Erodibility K (K mode only). | Used when steady_state_mode = K; ignored in relief mode. |
| `m_exponent` | `float` | `0.45` | `param` | - | Drainage-area exponent m in the stream-power / chi law. | Concavity is m/n. Set the same m and n as the LEM run you compare to. |
| `n_exponent` | `float` | `1.0` | `param` | - | Slope exponent n in the stream-power / chi law. | Concavity is m/n. |
| `baselevel_elevation` | `float` | `0.0` | `param` | - | Elevation assigned to the outlet(s) (base level). | Chi is integrated upstream from this elevation. Fluvial relief only. |

### Georeference

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `centre_latitude` | `float` | `-9999` | `param` | - | Latitude of the synthetic grid centre in decimal degrees. | Must be supplied together with centre_longitude. -9999 keeps the historical default synthetic anchor. |
| `centre_longitude` | `float` | `-9999` | `param` | - | Longitude of the synthetic grid centre in decimal degrees. | Must be supplied together with centre_latitude. -9999 keeps the historical default synthetic anchor. |
| `projection_epsg` | `int` | `0` | `param` | - | Projected CRS EPSG for synthetic georeference. | 0 derives WGS84 UTM from centre_latitude/centre_longitude. Non-zero uses the supplied EPSG directly. |

### Equilibration (optional, run_equilibration)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_equilibration` | `bool` | `false` | `switch` | - | Equilibrate the generated DEM in the same invocation (OCN -&gt; equilibrate). | Runs the landscape equilibration cycle on the synthetic DEM: pit-fill, re-route, re-snap, repeat. Without it the DEM is delivered exactly as the OCN synthesis produced it, and its D8 graph will differ from the OCN graph on roughly 6-9 percent of interior pixels. |
| `steady_state_n_cycles` | `int` | `10` | `param` | `run_equilibration` | Number of equilibration cycles (a MAXIMUM when convergence checking is on). | Each cycle recomputes the flow graph from the previous cycle's elevations, letting divides migrate one pixel per cycle. |
| `steady_state_check_convergence` | `bool` | `false` | `param` | `run_equilibration` | Stop early once the landscape stops changing. | Default false runs exactly steady_state_n_cycles cycles. |
| `steady_state_tolerance_dz_over_relief` | `float` | `1e-3` | `param` | `steady_state_check_convergence` | Convergence tolerance on max\|dz\|/relief between cycles. | Used when steady_state_tolerance_dz_metres is 0. |
| `steady_state_tolerance_dz_metres` | `float` | `0` | `param` | `steady_state_check_convergence` | Absolute convergence tolerance on max\|dz\| (m) between cycles. | If &gt; 0 this overrides steady_state_tolerance_dz_over_relief. |
| `steady_state_tolerance_receiver_fraction` | `float` | `1e-3` | `param` | `steady_state_check_convergence` | Convergence tolerance on the fraction of pixels changing receiver. | Set &lt;= 0 to disable. A landscape can stop moving vertically while divides are still migrating, which is what this criterion catches. |
| `steady_state_convergence_patience` | `int` | `2` | `param` | `steady_state_check_convergence` | Consecutive converged cycles required before stopping. | Cycle 1 is never eligible: it has no predecessor to compare against. |
| `divide_merge_search_steps` | `int` | `5` | `param` | `run_equilibration` | Downstream search depth for the divide diagnostics. Must be &gt; 0. | Two adjacent pixels are a divide pair when their flow paths do not meet within this many receiver steps. |
| `divide_use_d8_diagonals` | `bool` | `true` | `param` | `run_equilibration` | Include diagonal neighbours when enumerating divide pairs. | Affects diagnostics and convergence checks only, not the routing. |
| `divide_gradient_target` | `float` | `0` | `param` | `run_equilibration` | Reporting-only divide-gradient acceptance target. | When &gt; 0 the final cycle reports whether the median divide gradient is at or below this value. It does not change stopping. |
| `divide_area_by_dx` | `bool` | `true` | `param` | `run_equilibration` | Divide drainage area by cell size in the chi integrand. | Match this to how the OCN elevations were synthesised. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions (north/east/south/west): o=open, n=no-flux, p=periodic, b=forced outlet. | Used for equilibration routing AND, for the 'p' flags, for the annealing lattice itself: 'p' on east/west makes column 0 and column n_cols-1 D8 neighbours during generation, so the generated network and the routed DEM see the same lattice. East and west must match. Pair 'opop' with outlet_mode: north_south for a ridge-and-valley strip. |
| `fill_min_slope` | `float` | `1e-4` | `param` | `run_equilibration` | Minimum slope imposed by the pit-fill between cycles. | - |
| `n_threads` | `int` | `0` | `param` | - | Number of threads (0 = all available). | - |

### Hillslopes inside the equilibration cycle

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_hillslope_steady_state` | `bool` | `false` | `param` | `run_equilibration` | Solve hillslope steady state inside the equilibration cycle. | Requires run_equilibration: true and incision_threshold_area_m2 &gt; 0; channel pixels become the Dirichlet set that pins the solve. |
| `hillslope_ss_model` | `string` | `roering` | `param` | `run_hillslope_steady_state` | Hillslope solver: linear, roering, cidre, or critical_slope. | linear/roering are elliptic solves that FIND the divide position and reach a fixed point. critical_slope IMPOSES a gradient by upslope integration and may never settle -- it delivers a far larger divide step. Prefer an elliptic model unless you specifically want the critical-slope construction. |
| `incision_threshold_area_m2` | `float` | `0` | `param` | `run_hillslope_steady_state` | Drainage area threshold (m^2) separating channels from hillslopes. | Required (&gt; 0) when run_hillslope_steady_state is true. |
| `hillslope_diffusivity` | `float` | `0.005` | `param` | `run_hillslope_steady_state` | Hillslope diffusivity D (m^2/yr). | - |
| `critical_hillslope_slope` | `float` | `0.8` | `param` | `run_hillslope_steady_state` | Critical hillslope gradient Sc (dimensionless). | - |
| `enable_nonlinear_hillslopes` | `bool` | `true` | `param` | `run_hillslope_steady_state` | Use the nonlinear transport law. | - |
| `hillslope_ss_tolerance` | `float` | `1e-4` | `param` | `run_hillslope_steady_state` | Convergence criterion for the hillslope solve: max\|delta_z\|. | - |
| `hillslope_ss_max_iterations` | `int` | `10000` | `param` | `run_hillslope_steady_state` | Maximum iterations for the hillslope solver. | - |
| `hillslope_sor_omega` | `float` | `0` | `param` | `run_hillslope_steady_state` | SOR relaxation factor omega (0 = auto-estimate). | The auto-estimate is tuned for the linear solver and has been observed to DIVERGE with hillslope_ss_model: roering. Set 1.3 if that happens. |
| `hillslope_every_cycle` | `bool` | `true` | `param` | `run_hillslope_steady_state` | Solve the hillslopes on every cycle, not just the last. | Costs one extra solve per cycle, but it is the only way to show the divide is stable under re-routing. |
| `channel_hysteresis_factor` | `float` | `0.5` | `param` | `run_hillslope_steady_state` | Hysteresis on the channel (Dirichlet) set between cycles. | A channel pixel stays a channel until its area falls below this factor times incision_threshold_area_m2. 1.0 disables the hysteresis. |

### Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_dem_raster` | `bool` | `true` | `switch` | - | Write the synthetic elevation (DEM) raster. | Output file: &lt;write_fname&gt;_ocn_dem.tif. |
| `write_drainage_area_raster` | `bool` | `true` | `switch` | - | Write the OCN drainage-area (pixel count) raster. | Enabled by default because drainage area is a primary OCN result used to interpret the generated DEM, whereas analysis CLIs treat it as an optional diagnostic. Output file: &lt;write_fname&gt;_ocn_drainage_cells.tif. |
| `verify_routing` | `bool` | `true` | `param` | - | Route the synthetic DEM through D8 and report graph fidelity. | Fidelity = fraction of interior pixels whose D8 receiver equals the OCN receiver. ~1 for relaxed (low-T) networks; lower for random (high-T). |
| `write_equilibrated_dem_raster` | `bool` | `true` | `switch` | - | Write the equilibrated DEM (when run_equilibration is true). | Output file: &lt;write_fname&gt;_ocn_equilibrated_dem.tif. The un-equilibrated DEM is still written separately as &lt;write_fname&gt;_ocn_dem.tif, so the before/after pair is always available. |
| `write_snap_convergence_csv` | `bool` | `true` | `switch` | - | Write per-cycle equilibration diagnostics (when run_equilibration is true). | Output file: &lt;write_fname&gt;_ocn_snap_convergence.csv. Read this when a run does not converge. |
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | cog (default), gtiff, or same_as_input. |

<!-- END GENERATED CLI PARAMETERS -->
