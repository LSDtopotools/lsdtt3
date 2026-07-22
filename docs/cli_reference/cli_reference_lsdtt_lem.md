# CLI reference: lsdtt-lem

Landscape evolution model combining implicit stream-power fluvial incision, hillslope sediment transport, and uniform tectonic uplift. Evolved from MuddPILE and DONKmodel. Supports both time-stepping and analytical steady-state modes, with optional hillslope steady-state solvers.

## Usage

```
lsdtt-lem <path> <parameter_file>
lsdtt-lem <parameter_file>         # path defaults to ./
lsdtt-lem -h                       # print help and exit
lsdtt-lem -v                       # print version and exit
lsdtt-lem -c                       # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

Two mutually exclusive modes are available. Select one by setting `run_steady_state`.

### Analytical steady-state mode (`run_steady_state: true`)

| Stage | Description |
|---|---|
| 1 | Build or load initial DEM |
| 2 | Pit-fill and build flow routing |
| 3 | Compute fluvial steady-state elevations (one or more snapping cycles) |
| 4 | Write `_steady_state_elevation` raster and optional hillshade PNG |
| 5 | Solve hillslope steady state (optional, requires `run_hillslope_steady_state`) |
| 6 | Write `_hillslope_steady_state_elevation` raster and optional hillshade PNG |

### Time-stepping mode (`run_steady_state: false`, default)

| Stage | Description |
|---|---|
| 1 | Build or load initial DEM |
| 2 | Apply surface conditioning (parabola, breach, noise, fill) |
| 3 | Write `_initial_topography` raster |
| 4 | Time loop: fluvial incision → hillslope transport → uplift |
| 5 | Write intermediate rasters at `frequency_of_output_writing` (optional) |
| 6 | Write `_modelled_topography` and `_model_elevation_stats.csv` |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename (without extension). |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Output raster format. `cog`, `gtiff`, `same_as_input`, or `png`. PNG produces 8-bit greyscale; float32 elevation is normalised to 1–255 (nodata→0). A worldfile (`.pgw`) and GDAL PAM sidecar (`.png.aux.xml`) are written alongside for georeferencing. Use `cog`/`gtiff` for analysis. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. |

---

## Initial surface

### Loading from file

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_initial_dem` | bool | `false` | Load initial topography from `read_fname` instead of building a synthetic surface. |

### Synthetic surface (when `read_initial_dem: false`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `nrows` | int | `200` | Number of rows in the model domain. |
| `ncols` | int | `400` | Number of columns in the model domain. |
| `grid_resolution` | float | `30` | Cell size in metres. |
| `initial_surface_type` | string | `flat` | Surface type: `flat`, `white_noise`, `red_noise`, `diamond_square`. `flat` = all zeros. |
| `max_initial_elevation` | float | `1.0` | Maximum elevation (m) of a `white_noise` surface. Minimum is 0. |
| `noise_target_stddev` | float | `1.0` | Target standard deviation (m) for `red_noise` and `diamond_square` surfaces. |
| `noise_blur_passes` | int | `6` | Number of blur passes for `red_noise` surface generation. |
| `diamond_square_roughness` | float | `0.55` | Roughness exponent (0–1) for `diamond_square`. Higher = rougher. |
| `random_seed` | int | `33` | Random seed for initial noise surface. |

### Surface conditioning (applied to synthetic surfaces only)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `parabola_relief` | float | `0` | Peak elevation (m) of a north–south parabola added to the initial surface. Biases drainage toward north and south boundaries. Set `> 0` to enable. |
| `post_breach_noise_relief` | float | `0` | Half-amplitude (m) of uniform noise added after depression breaching. Breaks up post-breach flats. Set `> 0` to enable. |
| `fill_initial_surface` | bool | `false` | Apply pit-filling to the initial surface. Applied after breaching and post-breach noise. |
| `fill_min_slope` | float | `0.0001` | Minimum slope imposed by pit-filling. |

---

## Process switches

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enable_d8_fluvial_erosion` | bool | `true` | Enable D8 stream-power fluvial incision (O'Callaghan & Mark 1984). |
| `enable_dinf_fluvial_erosion` | bool | `false` | Enable D-infinity stream-power fluvial incision (Tarboton 1997). |
| `enable_hillslope_transport` | bool | `false` | Enable hillslope sediment transport. |
| `enable_nonlinear_hillslopes` | bool | `false` | Use nonlinear (Roering 1999) transport within `roering` and `cidre` hillslope models. |
| `enable_uplift` | bool | `true` | Apply uniform tectonic uplift. |

---

## Fluvial parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `erodibility` | float | `1e-3` | Stream-power erodibility K (m^(1−2m) yr⁻¹). |
| `m_exponent` | float | `0.45` | Drainage-area exponent m. |
| `n_exponent` | float | `1.0` | Slope exponent n. n=1 uses a closed-form solver; n≠1 uses Newton–Raphson. |
| `incision_threshold_area` | float | `0` | Minimum drainage area (m²) for fluvial incision. `0` = no threshold. Also defines channel pixels for hillslope coupling. |
| `divide_area_by_dx` | bool | `true` | Divide drainage area by cell size to obtain specific area before comparing against `incision_threshold_area`. |
| `slope_clamp` | float | `1e-5` | Minimum slope used as Newton–Raphson fallback in the stream-power solver. |

---

## Hillslope transport parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `hillslope_model` | string | `cidre` | Hillslope transport model for the time-stepping LEM. `linear`: explicit FD of D·∇²z (no flow graph, stable at any typical LEM timestep). `roering`: explicit FD of the Roering (1999) nonlinear law; auto-subcycles when near Sc. `cidre`: Carretier et al. (2016) finite-volume forward sweep with D-infinity routing. `roering_implicit`: implicit (unconditionally stable) backward-Euler / TR-BDF2 solver for the Roering law with channels as Dirichlet sinks — no CFL collapse on fine, near-critical grids (requires Eigen; falls back to explicit `roering` if Eigen is absent). Requires `enable_hillslope_transport: true`. |
| `hillslope_diffusivity` | float | `1e-3` | Hillslope diffusivity κ (m² yr⁻¹). |
| `critical_hillslope_slope` | float | `0.8` | Critical slope Sc (m/m) for nonlinear hillslope transport. |
| `channel_mask_fname` | string | `` | Optional raster (in `read_path`) whose cells > 0.5 are fixed-elevation channel pixels for the hillslope Dirichlet set, overriding the drainage-area threshold. Must match the DEM grid. Empty = derive channels from drainage area. |

### `roering_implicit` solver knobs

These only apply when `hillslope_model: roering_implicit`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `hillslope_implicit_time_integrator` | int | `0` | Time integrator: `0` = backward Euler (1st order, robust), `1` = TR-BDF2 (2nd order, L-stable, embedded-error PI step control — use for time-accurate transients). |
| `hillslope_implicit_error_tol` | float | `0` | TR-BDF2 embedded-error tolerance (elevation units) for the PI step controller (`0` = auto = a small fraction of the relief). Only used when `time_integrator = 1`. |
| `hillslope_implicit_solver_choice` | int | `0` | Linear solver: `0` = SimplicialLDLT (robust, moderate grids), `1` = CG + IncompleteCholesky (scales better on very large grids). |
| `hillslope_implicit_relaxation_omega` | float | `0.7` | Damped-Picard relaxation factor ω (≤ 1) for the nonlinear iteration. |
| `hillslope_implicit_picard_tol` | float | `1e-4` | Picard convergence tolerance on max\|Δz\| between iterates. |
| `hillslope_implicit_picard_max_iter` | int | `50` | Maximum Picard iterations per sub-step. |
| `hillslope_implicit_slope_ratio_cap` | float | `0.99` | Regularised slope-ratio cap r_max (< 1) bounding D_eff near Sc. |
| `hillslope_implicit_dt_init` | float | `0` | Initial adaptive sub-step (`0` = auto = the LEM timestep). |
| `hillslope_implicit_dt_min` | float | `0` | Floor on the adaptive sub-step (`0` = auto). |
| `hillslope_implicit_dt_max` | float | `0` | Cap on the adaptive sub-step (`0` = the LEM timestep). |

---

## Uplift parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `uplift_rate` | float | `1e-4` | Tectonic uplift rate (m yr⁻¹). Applied uniformly to all valid pixels. |

---

## Routing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `boundary_conditions` | string | `oooo` | Flow-routing boundary conditions: four characters for north, east, south, west. `o` = open, `n` = no-flux, `p` = periodic, `b` = forced outlet. |

---

## Analytical steady-state mode

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_steady_state` | bool | `false` | Compute analytical fluvial steady-state topography. Exits after writing output. |
| `steady_state_mode` | string | `K` | `K`: erodibility is specified by the user. `relief`: target relief is specified and K is back-calculated. |
| `target_relief` | float | `1000` | Target maximum relief (m) for relief-mode. The highest point will be `baselevel_elevation + target_relief`. |
| `baselevel_elevation` | float | `0` | Elevation (m) assigned to outlet nodes in the steady-state solution. |
| `steady_state_n_cycles` | int | `1` | Number of snapping cycles. Each cycle recomputes flow routing from the previous cycle's elevations, allowing drainage divides to migrate. Iterating 5–20 cycles smooths sharp elevation breaks across divides. |

---

## Hillslope steady-state

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_hillslope_steady_state` | bool | `false` | Solve hillslope steady-state elevations after fluvial steady-state. Requires `run_steady_state: true` and `incision_threshold_area > 0`. Channel pixels become Dirichlet boundary conditions. |
| `hillslope_ss_model` | string | `roering` | Hillslope steady-state solver. `linear`: Red-Black SOR solving D·∇²z = −U. `roering`: Red-Black SOR solving the nonlinear Roering (1999) law. `cidre`: pseudo-transient Carretier (2016) iterated to steady state (rebuilds D-inf flow graph each iteration). At steady state the nonlinear Roering and Cidre solvers produce the same elevation field. |
| `hillslope_ss_tolerance` | float | `1e-4` | Convergence criterion: max|Δz| over hillslope pixels per iteration. |
| `hillslope_ss_max_iterations` | int | `10000` | Maximum iterations for the hillslope steady-state solver. |
| `hillslope_sor_omega` | float | `0` | SOR relaxation factor ω (used by `linear` and `roering` solvers). `0` = auto-estimate from grid size (capped at 1.95). Reduce toward 1.3 if the nonlinear solver diverges. |

---

## Hillshade output (steady-state only)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_hillshade_png` | bool | `false` | Write a Horn (1981) hillshade of the final steady-state surface as PNG. If `run_hillslope_steady_state` is also true, a second hillshade of the combined surface is written. PNG is always used regardless of `write_output_format`. |
| `hillshade_azimuth` | float | `315` | Sun azimuth, degrees clockwise from north. 315 = northwest (ESRI/QGIS default). |
| `hillshade_altitude` | float | `45` | Sun altitude above the horizon, degrees. |
| `hillshade_z_factor` | float | `1.0` | Vertical exaggeration applied before computing the hillshade. |

---

## Time-stepping

| Parameter | Type | Default | Description |
|---|---|---|---|
| `timestep` | float | `100` | Model timestep (yr). |
| `number_of_iterations` | int | `100` | Number of timesteps. Ignored when `stop_at_steady_state` is true. |
| `stop_at_steady_state` | bool | `false` | Stop when mean erosion rate is within 1 % of uplift rate for 100 consecutive steps, then run 20× that step count more before stopping. |
| `frequency_of_output_writing` | int | `0` | Write an intermediate raster every N iterations. `0` = write final only. |

---

## Output files

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_elevation_stats_csv` | bool | `true` | `<write_fname>_model_elevation_stats.csv` — per-timestep elevation statistics. |

### Steady-state outputs (when `run_steady_state: true`)

| File | Description |
|---|---|
| `<write_fname>_steady_state_elevation.tif` | Fluvial steady-state elevation raster. |
| `<write_fname>_steady_state_stats.csv` | Steady-state solver diagnostics (K, relief, mode). Written when `write_elevation_stats_csv: true`. |
| `<write_fname>_steady_state_hillshade.png` | Hillshade PNG of fluvial steady-state surface (when `write_hillshade_png: true`). |
| `<write_fname>_hillslope_steady_state_elevation.tif` | Combined fluvial + hillslope steady-state elevation (when `run_hillslope_steady_state: true`). |
| `<write_fname>_hillslope_steady_state_stats.csv` | Hillslope solver diagnostics. Written when `write_elevation_stats_csv: true`. |
| `<write_fname>_hillslope_steady_state_hillshade.png` | Hillshade PNG of combined surface (when `write_hillshade_png: true`). |

### Time-stepping outputs

| File | Description |
|---|---|
| `<write_fname>_initial_topography.tif` | Initial surface written before the time loop. |
| `<write_fname>_model_iteration_N.tif` | Intermediate snapshot written every `frequency_of_output_writing` steps. |
| `<write_fname>_modelled_topography.tif` | Final topography after the run completes. |
| `<write_fname>_model_elevation_stats.csv` | Per-timestep statistics: iteration, time, delta max/mean elevation, max/mean erosion rate. |

---

## Example parameter files

### Analytical fluvial steady state

```
read_path:  /data/dems/
read_fname: initial_dem
write_path: /data/output/
write_fname: run01

run_steady_state:       true
steady_state_mode:      K
steady_state_n_cycles:  10

enable_d8_fluvial_erosion: true
erodibility:            0.0005
m_exponent:             0.45
n_exponent:             1.0
uplift_rate:            0.0001

write_output_format:    cog
write_hillshade_png:    true
```

### Fluvial + nonlinear hillslope steady state

```
read_path:  /data/dems/
read_fname: initial_dem
write_path: /data/output/
write_fname: run02

run_steady_state:             true
run_hillslope_steady_state:   true
steady_state_n_cycles:        10

erodibility:                  0.0005
uplift_rate:                  0.0001
incision_threshold_area:      1000000

hillslope_ss_model:           roering
hillslope_diffusivity:        0.005
critical_hillslope_slope:     0.8
hillslope_ss_max_iterations:  20000
hillslope_ss_tolerance:       0.0001

write_hillshade_png:          true
```

### Time-stepping LEM with Cidre hillslope transport

```
write_path:  /data/output/
write_fname: lem_run

nrows:               200
ncols:               400
grid_resolution:     30
initial_surface_type: white_noise
max_initial_elevation: 1.0
parabola_relief:     5.0
post_breach_noise_relief: 0.05
random_seed:         42

enable_d8_fluvial_erosion:   true
enable_hillslope_transport:  true
enable_nonlinear_hillslopes: true
hillslope_model:             cidre
enable_uplift:               true

erodibility:                 0.001
m_exponent:                  0.45
n_exponent:                  1.0
hillslope_diffusivity:       0.005
critical_hillslope_slope:    0.8
uplift_rate:                 0.0001
incision_threshold_area:     500000
timestep:                    500
number_of_iterations:        2000
frequency_of_output_writing: 200

write_elevation_stats_csv:   true
```
