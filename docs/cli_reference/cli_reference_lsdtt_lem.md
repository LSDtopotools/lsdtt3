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

Stage 3 has two implementations. The default is the legacy snapping loop, which
runs exactly `steady_state_n_cycles` cycles. Setting
`use_landscape_equilibration_cycle: true` selects the equilibration cycle
instead, which adds per-cycle divide diagnostics, convergence-based early
stopping, and the option to re-solve the hillslopes on every cycle so the divide
is tested under re-routing. See *Landscape equilibration cycle* below.

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
| `n_rows` | int | `200` | Number of rows in the model domain. Former `nrows` is an unchanged-value alias. |
| `n_cols` | int | `400` | Number of columns in the model domain. Former `ncols` is an unchanged-value alias. |
| `cell_size_metres` | float | `30` | Cell size in metres. Former `grid_resolution` is an unchanged-value alias. |
| `noise_type` | string | `flat` | Surface type: `flat`, `white_noise`, `red_noise`, `diamond_square`. Former `initial_surface_type` is an unchanged-value alias. |
| `max_initial_elevation` | float | `1.0` | Maximum elevation (m) of a `white_noise` surface. Minimum is 0. |
| `noise_target_stddev_metres` | float | `1.0` | Target standard deviation for `red_noise` and `diamond_square` surfaces. Former `noise_target_stddev` is an unchanged-value alias. |
| `noise_blur_passes` | int | `6` | Number of blur passes for `red_noise` surface generation. |
| `diamond_square_roughness` | float | `0.55` | Roughness exponent (0–1) for `diamond_square`. Higher = rougher. |
| `random_seed` | int | `33` | Random seed for initial noise surface. |

### Surface conditioning (applied to synthetic surfaces only)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `initial_surface_slope` | float | `0.0001` | Regional slope added toward the nearest declared open (`o`/`b`) edge. Prevents unresolved synthetic flats and gives initial flow a boundary-consistent direction. Set to `0` to disable. Closed and periodic domains are unchanged. |
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
| `incision_threshold_specific_area` | float | `0` | **With `divide_area_by_dx: true` (default):** specific catchment area A/dx in metres. **With `false`:** total map area A in m². `0` = no threshold. Also defines channel pixels for hillslope coupling. The former `incision_threshold_area` spelling is an unchanged-value alias. |
| `divide_area_by_dx` | bool | `true` | Divide drainage area by cell size to obtain specific area before comparing against `incision_threshold_specific_area`. |
| `slope_clamp` | float | `1e-5` | Minimum slope used as Newton–Raphson fallback in the stream-power solver. |

---

## Hillslope transport parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `hillslope_model` | string | `cidre` | `cellular`: fast conservative soft-critical face transfers with optional nonlocal transport/tracing. `cidre`: graph-routed Carretier transport. `roering_implicit`: Newton backward Euler. `adi`: linear ADI. |
| `hillslope_diffusivity` | float | `1e-3` | Hillslope diffusivity κ (m² yr⁻¹). |
| `critical_hillslope_slope` | float | `0.8` | Critical slope Sc (m/m) for nonlinear hillslope transport. |
| `channel_mask_fname` | string | `` | Optional raster (in `read_path`) whose cells > 0.5 are fixed-elevation channel pixels for the hillslope Dirichlet set, overriding the drainage-area threshold. Must match the DEM grid. Empty = derive channels from drainage area. |
| `channel_mask_threshold_specific_area` | float | `0` | **With `divide_area_by_dx: true` (default):** A/dx in metres. **With `false`:** map area A in m². Threshold for a written channel mask; `0` uses `incision_threshold_specific_area`. The former `channel_mask_threshold_area` spelling is an unchanged-value alias. |

### `cellular` transport and tracing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `cellular_r_match` | float | `0.85` | Roering/finite-continuation matching slope ratio. |
| `cellular_cfl_safety` | float | `0.8` | Explicit stability safety factor `(0,1]`. |
| `cellular_nonlocal_fraction` | float | `0` | Fraction `[0,1]` of near-critical excess flux routed nonlocally. |
| `cellular_nonlocal_activation_ratio` | float | `0.8` | Slope ratio where nonlocal routing begins. |
| `cellular_nonlocal_deposition_probability` | float | `0.5` | Geometric deposition probability per cell. |
| `cellular_nonlocal_max_hops` | int | `4` | Maximum nonlocal distance in cells, at least 2 when enabled. |
| `enable_cellular_particles` | bool | `false` | Enable sparse Lagrangian superparticles. |
| `cellular_particle_inverse_density` | int | `100` | Active cells per seeded particle. |
| `cellular_particle_mixing_depth` | float | `1` | Well-mixed active-layer depth (m). |
| `cellular_particle_seed` | int | `1` | Deterministic sampling seed. |
| `cellular_particle_restart_fname` | string | `` | Restart CSV in `read_path`; empty seeds new particles. |

Tracing writes `<write_fname>_cellular_particles.csv`. Reusing it as a restart preserves the RNG
substep counter, provenance, ages and distances.

### `roering_implicit` solver knobs

These only apply when `hillslope_model: roering_implicit`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `hillslope_implicit_time_integrator` | int | `0` | Time integrator: `0` = backward Euler (1st order, robust), `1` = TR-BDF2 (2nd order, L-stable, embedded-error PI step control — use for time-accurate transients). |
| `hillslope_implicit_error_tolerance_metres` | float | `0` | TR-BDF2 embedded-error tolerance (elevation units) for the PI step controller (`0` = auto = a small fraction of the relief). Only used when `time_integrator = 1`. |
| `hillslope_implicit_solver_choice` | int | `0` | Linear solver: `0` = SimplicialLDLT (robust, moderate grids), `1` = CG + IncompleteCholesky (scales better on very large grids). |
| `hillslope_implicit_relaxation_omega` | float | `0.7` | Damped-Picard relaxation factor ω (≤ 1) for the nonlinear iteration. |
| `hillslope_implicit_picard_tolerance_metres` | float | `1e-4` | Picard convergence tolerance on max\|Δz\| between iterates. |
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
| `steady_state_n_cycles` | int | `1` | Number of snapping cycles. Each cycle recomputes flow routing from the previous cycle's elevations, allowing drainage divides to migrate. Iterating 5–20 cycles smooths sharp elevation breaks across divides. **When `steady_state_check_convergence: true` this becomes a *maximum*, not a fixed count** — the run may stop earlier. |

---

## Landscape equilibration cycle

An alternative implementation of the snapping loop, enabled with
`use_landscape_equilibration_cycle: true` (requires `run_steady_state: true`).
It adds per-cycle divide diagnostics, convergence-based early stopping, and the
option to solve hillslopes on every cycle rather than only the last.

A single-pass snap leaves a sharp elevation break along drainage divides, because
each side is built by integrating up from its own channel with nothing forcing
the two to meet. Cycling re-routes and re-snaps so the divide can migrate; the
diagnostics below are how you tell whether it actually settled.

> **D-infinity is refused, not silently downgraded.** The equilibration cycle
> routes D8 only. Setting both `enable_dinf_fluvial_erosion: true` and
> `use_landscape_equilibration_cycle: true` is a **fatal error** naming both
> parameters, rather than a quiet fall back to D8 — silently routing D8 when
> D-infinity was requested is the one failure here that could mislead a science
> result.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `use_landscape_equilibration_cycle` | bool | `false` | Route analytical steady state through the equilibration-cycle module. `false` keeps the legacy steady-state behaviour unchanged. |
| `divide_merge_search_steps` | int | `5` | Downstream search depth for the divide diagnostics. Two adjacent pixels are a divide pair when their flow paths do not meet within this many receiver steps. Must be `> 0`. |
| `divide_use_d8_diagonals` | bool | `true` | Include diagonal neighbours when enumerating divide pairs (D8 adjacency). Affects the diagnostics and convergence checks only, not the routing. |
| `divide_gradient_target` | float | `0` | Reporting-only acceptance target. When `> 0`, the final cycle reports whether the median divide gradient is at or below this value. **It does not change stopping.** |

### Convergence and early stopping

Off by default: with `steady_state_check_convergence: false` the run always
executes `steady_state_n_cycles` cycles.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `steady_state_check_convergence` | bool | `false` | Enable convergence-based early stopping. |
| `steady_state_tolerance_dz_over_relief` | float | `1e-3` | Convergence tolerance on max\|Δz\|/relief between cycles. Used when `steady_state_tolerance_dz_metres` is `0`. |
| `steady_state_tolerance_dz_metres` | float | `0` | Absolute tolerance on max\|Δz\| (m) between cycles. When `> 0` this **overrides** `steady_state_tolerance_dz_over_relief`. |
| `steady_state_tolerance_receiver_fraction` | float | `1e-3` | Convergence tolerance on the fraction of pixels changing receiver each cycle. Set `<= 0` to disable this criterion. |
| `steady_state_convergence_patience` | int | `2` | Consecutive converged cycles required before stopping. Cycle 1 is never eligible, since there is no previous cycle to compare against. |

A cycle counts as converged only when the elevation criterion **and** the
receiver-fraction criterion both hold; a criterion whose tolerance is disabled
(`<= 0`) is treated as satisfied. The count must be *consecutive* — any cycle
that fails resets the streak to zero.

The two criteria are not redundant: a landscape can stop moving vertically while
divides are still migrating, so max\|Δz\| alone can look converged while the flow
graph is still changing. That is the failure the receiver-fraction criterion
exists to catch. Disabling both leaves nothing but patience, which is why cycle 1
is excluded — it has no predecessor, `frac_receiver_changed` is forced to `0.0`,
and a patience-1 run would otherwise declare convergence having compared nothing.

---

## Hillslope steady-state

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_hillslope_steady_state` | bool | `false` | Solve hillslope steady-state elevations after fluvial steady-state. Requires `run_steady_state: true` and `incision_threshold_specific_area > 0`. Channel pixels become Dirichlet boundary conditions. |
| `hillslope_ss_model` | string | `roering` | Hillslope steady-state solver. `linear`: Red-Black SOR solving D·∇²z = −U. `roering`: Red-Black SOR solving the nonlinear Roering (1999) law. `cidre`: pseudo-transient Carretier (2016) iterated to steady state (rebuilds D-inf flow graph each iteration). At steady state the nonlinear Roering and Cidre solvers produce the same elevation field. |
| `hillslope_ss_tolerance` | float | `1e-4` | Convergence criterion: max|Δz| over hillslope pixels per iteration. |
| `hillslope_ss_max_iterations` | int | `10000` | Maximum iterations for the hillslope steady-state solver. |
| `hillslope_sor_omega` | float | `0` | SOR relaxation factor ω (used by `linear` and `roering` solvers). `0` = auto-estimate from grid size (capped at 1.95). Reduce toward 1.3 if the nonlinear solver diverges. |

### Inside the equilibration cycle

These apply only when `use_landscape_equilibration_cycle` and
`run_hillslope_steady_state` are both true.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `hillslope_every_cycle` | bool | `true` | Solve the hillslopes on every cycle, not just the last one. Costs one extra solve per cycle, but it is the only way to show the divide is stable **under re-routing** — with `false` the hillslopes are solved once on the delivered surface, so a divide that would move if the hillslopes were re-solved never gets the chance to. |
| `channel_hysteresis_factor` | float | `0.5` | Hysteresis on the channel (Dirichlet) set between cycles. A pixel already flagged as channel stays channel until its drainage area falls below `channel_hysteresis_factor × incision_threshold_specific_area`. Stops the Dirichlet set flickering from cycle to cycle. `1.0` disables the hysteresis. |

> **If the nonlinear solver diverges, set `hillslope_sor_omega` manually.** With
> `hillslope_ss_model: roering` the auto-omega estimate has been observed to
> diverge at every diffusivity and tolerance tried — running to the iteration cap
> with a residual of order hundreds of metres, and producing correspondingly
> meaningless divide diagnostics. `hillslope_sor_omega: 1.3` converges. The
> auto-estimate is tuned for the linear solver; it is not a safe default for the
> nonlinear one.

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

### Intermediate diagnostics

When `frequency_of_output_writing` is greater than zero, optional PNG
diagnostics can be written at the same checkpoints as the elevation snapshots.
The D-infinity diagnostic is `log10` drainage area with an inverted grayscale
mapping, so larger drainage areas are darker. The diagnostics are also written
for the initial and final modelled topographies. PNG georeferencing sidecars are
removed by default because these files are intended for visual inspection.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_intermediate_hillshade_png` | bool | `false` | Write an intermediate hillshade PNG at every output checkpoint. |
| `intermediate_hillshade_write_path` | string | empty | Optional directory for intermediate hillshade PNGs; empty uses `write_path`. |
| `write_intermediate_log10_dinf_drainage_area_png` | bool | `false` | Write an inverted-grayscale `log10` D-infinity drainage-area PNG at every output checkpoint. |
| `intermediate_drainage_area_write_path` | string | empty | Optional directory for intermediate D-infinity area PNGs; empty uses `write_path`. |
| `write_intermediate_png_sidecars` | bool | `false` | Keep PNG georeferencing sidecars (`.pgw` and `.aux.xml`) instead of removing them. |

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

### Equilibration-cycle outputs

Both require `run_steady_state` and `use_landscape_equilibration_cycle` to be true.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_snap_convergence_csv` | bool | `true` | `<write_fname>_snap_convergence.csv` — one row per cycle. This is the file to read when a run does not converge. |
| `write_divide_step_raster` | bool | `false` | `<write_fname>_divide_step.<ext>` (extension from `write_output_format`) — the per-pixel elevation contrast across divide pairs on the final DEM, for locating *where* a residual divide step sits rather than just its distribution. |

`_snap_convergence.csv` columns, grouped:

| Group | Columns |
|---|---|
| identity / cost | `cycle`, `seconds` |
| elevation change | `max_abs_dz_m`, `mean_abs_dz_m`, `rms_dz_m`, `max_dz_over_relief` |
| routing change | `n_receiver_changed`, `frac_receiver_changed`, `n_receiver_changed_hillslope`, `frac_receiver_changed_hillslope`, `n_hillslope_compared` |
| cross-divide drainage | `n_cross_divide`, `frac_cross_divide` |
| divide step | `divide_n_divide_pairs`, `divide_step_median_m`, `divide_step_p90_m`, `divide_step_max_m`, `K_used`, `computed_relief` |
| divide gradient target | `divide_gradient_median`, `divide_gradient_p90`, `divide_gradient_target_used`, `divide_gradient_target_margin`, `divide_gradient_target_met` |
| limit cycle | `receiver_hash`, `limit_cycle_detected`, `limit_cycle_period`, `limit_cycle_first` |

The limit-cycle columns matter because **a snapping run need not converge to a
fixed point at all** — it can settle into a repeating orbit of flow graphs.
`receiver_hash` identifies the flow graph each cycle; a repeat means the
iteration has closed an orbit, and `limit_cycle_period` is its length
(`1` = a genuine fixed point). Critical-slope hillslopes are the case to watch:
they can wander indefinitely without ever closing an orbit, in which case
nothing is detected and the run simply stops at `steady_state_n_cycles`. An
elliptic hillslope solve (`linear` or `roering`) reaches a period-1 fixed point
instead, which is the practical reason to prefer one.

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
incision_threshold_specific_area: 1000000

hillslope_ss_model:           roering
hillslope_diffusivity:        0.005
critical_hillslope_slope:     0.8
hillslope_ss_max_iterations:  20000
hillslope_ss_tolerance:       0.0001

write_hillshade_png:          true
```

### Equilibrated steady state with divide diagnostics

Cycles until the divides stop moving rather than for a fixed count, solving the
hillslopes every cycle so the divide position is tested under re-routing.
`steady_state_n_cycles` is the cap, not the target.

```
read_path:  /data/dems/
read_fname: initial_dem
write_path: /data/output/
write_fname: run03

run_steady_state:                          true
use_landscape_equilibration_cycle:         true
steady_state_n_cycles:                     40

erodibility:                               0.0005
uplift_rate:                               0.0001
incision_threshold_specific_area:          1000000

steady_state_check_convergence:            true
steady_state_tolerance_dz_over_relief:     0.001
steady_state_tolerance_receiver_fraction:  0.001
steady_state_convergence_patience:         3

run_hillslope_steady_state:                true
hillslope_ss_model:                        roering
hillslope_diffusivity:                     0.005
critical_hillslope_slope:                  0.8
hillslope_sor_omega:                       1.3
hillslope_every_cycle:                     true

divide_merge_search_steps:                 5
divide_gradient_target:                    0.5

write_snap_convergence_csv:                true
write_divide_step_raster:                  true
```

Note `hillslope_sor_omega: 1.3` rather than the `0` auto-estimate — see the
warning under *Hillslope steady-state*. Leave `enable_dinf_fluvial_erosion` off:
combining it with the equilibration cycle is a fatal error, not a fallback.

### Time-stepping LEM with Cidre hillslope transport

```
write_path:  /data/output/
write_fname: lem_run

n_rows:           200
n_cols:           400
cell_size_metres: 30
noise_type:       white_noise
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
incision_threshold_specific_area: 500000
timestep:                    500
number_of_iterations:        2000
frequency_of_output_writing: 200

write_elevation_stats_csv:   true
```

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

### DEM loading vs synthetic surface

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_initial_dem` | `bool` | `false` | `param` | - | Load initial topography from read_fname instead of building a synthetic surface. | If true, the read_path / read_fname / read_ext from the parameter file are used. |

### Synthetic domain

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_rows` | `int` | `200` | `param` | - | Number of rows in the synthetic model domain. | Used only when the domain is generated rather than read from a DEM. |
| `n_cols` | `int` | `400` | `param` | - | Number of columns in the synthetic model domain. | Used only when the domain is generated rather than read from a DEM. |
| `cell_size_metres` | `float` | `30` | `param` | - | Cell size of the synthetic model domain in metres. | Sets the map scale of the generated surface. Anything expressed as a pixel count (for example threshold_contributing_pixels) corresponds to a different ground area as this changes. |
| `random_seed` | `int` | `33` | `param` | - | Seed for every stochastic generator used to build the domain. | A fixed seed makes a run reproducible. Change it to draw a different realisation of the same statistical surface. |

### Grid (when not loading a DEM)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `synthetic_crs_wkt` | `string` | - | `param` | - | Coordinate Reference System WKT string for synthetic surfaces. | If left empty (the default), this automatically assigns EPSG:32601 (UTM Zone 1N) and places the lower-left corner at 500000E, 0N. |

### Synthetic noise

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `noise_type` | `string` | `flat` | `param` | - | Type of noise surface to generate: flat, white_noise, red_noise, diamond_square. | flat = all zeros, i.e. no noise; white_noise = uniform random, equal power at every wavelength; red_noise = low-frequency dominated; diamond_square = fractal terrain whose roughness is set by diamond_square_roughness. |
| `noise_target_stddev_metres` | `float` | `1.0` | `param` | - | Target standard deviation of the generated noise surface (metres). | The raw surface is rescaled to this standard deviation. Used when noise_type is red_noise or diamond_square. |

### Initial conditions

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `max_initial_elevation` | `float` | `1.0` | `param` | - | Max elevation of white-noise surface (m). Min is 0. | Used when noise_type is white_noise. |
| `noise_blur_passes` | `int` | `6` | `param` | - | Number of blur passes for red_noise surface generation. | - |
| `diamond_square_roughness` | `float` | `0.55` | `param` | - | Roughness exponent (0-1) for diamond-square surface. Higher = rougher. | - |
| `parabola_relief` | `float` | `0` | `param` | - | Peak elevation (m) of a N-S parabola added to random initial surfaces. | The parabola is zero at the north and south boundaries and peaks at the centre row. Adding it biases drainage toward both N/S boundaries, reducing internal basins. Set &gt; 0 to enable; typical values 2-5x the noise relief. |
| `initial_surface_slope` | `float` | `0.0001` | `param` | - | Small regional slope added to generated surfaces toward the nearest declared open boundary. Prevents unresolved flats and gives initial flow a boundary-consistent direction. Set to 0 to disable. | - |
| `post_breach_noise_relief` | `float` | `0` | `param` | - | Half-amplitude (m) of uniform noise added after breaching the initial surface. | Breaking up post-breach flats prevents unnaturally regular drainage. Set &gt; 0 to enable; typical values 0.01-0.1 m. Breaching is performed whenever parabola_relief &gt; 0 or post_breach_noise_relief &gt; 0. |
| `fill_initial_surface` | `bool` | `false` | `param` | - | Apply pit-filling to the initial surface before running. | Applied after breaching/noise conditioning when those are also enabled. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope imposed by pit filling. | Used when fill_initial_surface is true. |

### Sinuous Boundaries

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `enable_sinuous_n_boundary` | `bool` | `false` | `switch` | - | Apply a sinuous (sine wave) mask to the North domain boundary. | - |
| `sinuous_n_amplitude_fraction` | `float` | `0.25` | `param` | - | Amplitude of the North sinuous boundary as a fraction of total N-S domain length. | This can max be 0.25, as it create waves that span half the N-S domain. |
| `sinuous_n_n_waves` | `int` | `1` | `param` | - | Number of complete sine waves along the North boundary. | - |
| `sinuous_n_phase_offset` | `float` | `0.0` | `param` | - | Phase offset of the North sinuous boundary as a fractional value [0, 1]. | - |
| `enable_sinuous_s_boundary` | `bool` | `false` | `switch` | - | Apply a sinuous (sine wave) mask to the South domain boundary. | - |
| `sinuous_s_amplitude_fraction` | `float` | `0.25` | `param` | - | Amplitude of the South sinuous boundary as a fraction of total N-S domain length. | This can max be 0.25, as it create waves that span half the N-S domain. |
| `sinuous_s_n_waves` | `int` | `1` | `param` | - | Number of complete sine waves along the South boundary. | - |
| `sinuous_s_phase_offset` | `float` | `0.0` | `param` | - | Phase offset of the South sinuous boundary as a fractional value [0, 1]. | - |

### Process switches

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `enable_d8_fluvial_erosion` | `bool` | `true` | `switch` | - | Enable D8 stream-power fluvial incision. | D8 with n=1 always uses closed-form and is therefore very stable. For other n uses Newton-Raphson. |
| `enable_dinf_fluvial_erosion` | `bool` | `false` | `switch` | - | Enable D-infinity stream-power fluvial incision. | Uses stable closed-form for n=1 only on nodes that have only a single receiver, otherwise uses Newton-Raphson. |
| `enable_hillslope_transport` | `bool` | `false` | `switch` | - | Enable Cidre hillslope diffusion. | - |
| `enable_nonlinear_hillslopes` | `bool` | `false` | `switch` | - | Use nonlinear transport (Roering 1999 or Carretier et al. 2016 depending on hillslope_model). | Requires enable_hillslope_transport. |
| `critical_slope_clamp_only` | `bool` | `false` | `switch` | - | With hillslope_ss_model: critical_slope, LOWER cells above the target gradient instead of imposing it everywhere. | The default imposes z = z_receiver + dx*S_c on every hillslope cell, manufacturing a landscape of uniform critical slopes and destroying whatever form the input had. Clamping preserves the input wherever it is already admissible and cuts back only the super-critical parts. Intended for handing an equilibrated LINEAR landscape to a Roering run: linear diffusion has no critical slope, so its steady state routinely carries gradients above S_c where the Roering flux law is singular. NOTE this limits the gradient ALONG THE FLOW PATH, not the full gradient magnitude the flux law uses -- cross-slope can still exceed S_c, so use critical_slope_clamp_fraction to leave headroom. |
| `critical_slope_clamp_fraction` | `float` | `0` | `param` | - | Multiply S_c by this before clamping. 0 = use S_c unchanged. | Clamping to exactly S_c leaves the Roering flux law at its singularity, where alpha = 1/(1-(S/S_c)^2) diverges. Clamping to 0.9 S_c gives alpha ~ 5, and 0.8 gives alpha ~ 2.8 -- both comfortably solvable. It also leaves headroom for the cross-slope component that the along-flow clamp does not control. |
| `hillslope_force_sor` | `bool` | `false` | `switch` | - | Force the legacy SOR + Picard steady-state solver instead of Newton. | Steady-state hillslope solves default to Newton, which is faster and more accurate for the linear problem and converges near Sc where Picard limit-cycles. SOR is still used automatically when Newton is unavailable (a build without Eigen) or fails to converge; this switch forces it for comparison or reproduction of an older run. With SOR the nonlinear operator needs hillslope_sor_omega set low (0.3-0.5); the auto value is the LINEAR optimum and diverges. |
| `hillslope_model` | `string` | `cidre` | `param` | - | Hillslope transport model: cellular, roering_implicit, cidre, or adi. | cellular: conservative soft-critical face transfers with optional bounded nonlocal transport and sparse Lagrangian tracing; permits finite-lived slopes above Sc and is the recommended fast nonlinear transient model. roering_implicit: unconditionally-stable implicit solver for the Roering (1999) nonlinear law. Steady-state runs use Newton; transient runs use backward Euler with adaptive sub-stepping. Channels are Dirichlet sinks from the area threshold or channel_mask_fname; needs a non-empty channel set or the step is skipped. Set enable_nonlinear_hillslopes false for linear diffusion. cidre: Carretier et al. (2016) finite-volume forward sweep using D-inf routing; explicit, so CFL-limited, but stable as slopes approach Sc via slope-restoration. Works with an emergent channel network (incision_threshold_specific_area = 0); since 2026-08-27 roering_implicit does too, for transient runs. adi: Crank-Nicolson implicit linear diffusion using uncoupled ADI sweeps (transient only; steady-state runs use Newton regardless). NOTE: 'linear' and 'roering' select EXPLICIT solvers and are run as asked: they are CFL-limited and less stable than the implicit path, and are retained for testing and for builds without Eigen. They are NOT redirected to roering_implicit -- if you want the implicit solver, ask for it by name. Requires enable_hillslope_transport. |

### Timestep control

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `adaptive_timestep` | `bool` | `false` | `switch` | - | Adapt the LEM timestep to keep the coupled solve accurate. | Both halves of a coupled run can be unconditionally STABLE and the run still be wrong at a large timestep: stability describes one operator, accuracy of an operator-split PAIR is a different question. The controller keys on the achieved maximum local elevation change of each step against the relief across one cell (dx * median gradient), holding that ratio near timestep_target_ratio. It shrinks faster than it grows, because an oversized step corrupts the landscape while an undersized one only costs time. `timestep` becomes the starting value. |
| `timestep_target_ratio` | `float` | `0.7` | `param` | - | Target for (elevation change per step) / (relief across one cell). | Measured on a 200x200 grid at 5 m: ratios of 0.13-0.67 gave 4-5% grid-scale roughness and a stable hillslope length; 2.67 gave 12.2% roughness and a hillslope length wrong by 5.7x, while still converging and still reporting ~99% of the uplift rate. Used both by the pre-run check and, with adaptive_timestep, by the controller. |
| `steady_state_dz_tolerance` | `float` | `0.0001` | `param` | - | Stop early when no cell moves more than this in one step (m). | The run ends at whichever comes first: this steady-state criterion or run_until_time / number_of_iterations. Convergence is REFUSED while the timestep is below steady_state_min_dt_fraction of its initial value, because an absolute per-step tolerance is trivially met by a collapsing timestep while the landscape is still evolving as fast as ever. Set 0 to disable early stopping. |
| `steady_state_rate_tolerance` | `float` | `0` | `param` | - | Steady state as a RATE (m/yr), dt-independent. 0 = use the per-step tolerance alone. | The per-step tolerance is dt-dependent: halving the timestep halves every change without the landscape being any steadier (the same relaxed landscape measured 0.006 m/step at dt = 250 and 0.029 m/step at dt = 2000). A rate means the same thing at any timestep. The measured floor for a coupled run is 1.4e-5 to 2.8e-5 m/yr -- 3.6-7.0% of the uplift rate -- set by continual reorganisation of the drainage network, so a tolerance below that cannot be met however long the run goes. A sensible choice is a small fraction of uplift_rate. |
| `steady_state_patience` | `int` | `5` | `param` | - | Consecutive steps below tolerance before declaring steady state. | One quiet step can happen by accident. |
| `steady_state_min_dt_fraction` | `float` | `0.25` | `param` | - | Refuse to declare steady state while dt is below this fraction of its initial value. | Guards against a collapsing adaptive timestep faking convergence. |
| `stop_on_stall` | `bool` | `false` | `switch` | - | End the run when a stall or oscillation is detected, rather than warning. | Off by default: stall detection is heuristic, and a false positive throws away a healthy run while a false negative only costs the time the run was going to take anyway. Convergence is a definite test and always stops the run; a stall only warns unless this is set. |
| `stall_window` | `int` | `1000` | `param` | - | Steps over which a lack of progress is judged (0 disables). | If max\|dz\| has not fallen over this window while still above tolerance, the run is reported as STALLED, and as OSCILLATING when the step-to-step changes alternate in sign -- a period-2 cycle, which points at the numerics (timestep too large, or a limit cycle in the nonlinear solve) rather than at needing more model time. A stalled run is NOT a converged one and must not be reported as steady state. |
| `run_until_time` | `float` | `0` | `param` | - | Run until this much model time has elapsed (yr). 0 = use number_of_iterations instead. | With a FIXED timestep the two are equivalent, but with adaptive_timestep they are not: if the controller shrinks dt, a fixed iteration count covers proportionally less model time, and different runs in a sweep end at different times. Measured on a 7-run sweep intended to cover 6 Myr each: the runs actually reached 0.70, 1.12 and 1.26 Myr because dt fell by up to 15x. Set this instead of number_of_iterations whenever the timestep adapts. number_of_iterations then acts as a safety cap. |
| `timestep_min` | `float` | `0` | `param` | - | Floor for the adaptive timestep (yr). 0 = timestep / 100. | - |
| `timestep_max` | `float` | `0` | `param` | - | Ceiling for the adaptive timestep (yr). 0 = timestep * 100. | - |

### Cellular hillslopes

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `cellular_r_match` | `float` | `0.85` | `param` | `hillslope_model` | Slope ratio below which cellular transport exactly follows Roering. | Above this ratio the mobility follows a finite tangent continuation. Only used when hillslope_model is cellular. |
| `cellular_cfl_safety` | `float` | `0.8` | `param` | `hillslope_model` | Explicit cellular stability safety factor in (0,1]. | Only used when hillslope_model is cellular. |
| `cellular_nonlocal_fraction` | `float` | `0` | `param` | `hillslope_model` | Fraction [0,1] of near-critical excess flux routed nonlocally. | Zero gives local soft-critical transport. Only used when hillslope_model is cellular. |
| `cellular_nonlocal_activation_ratio` | `float` | `0.8` | `param` | `cellular_nonlocal_fraction` | Slope ratio at which nonlocal routing begins. | Only used when cellular_nonlocal_fraction is greater than 0. |
| `cellular_nonlocal_deposition_probability` | `float` | `0.5` | `param` | `cellular_nonlocal_fraction` | Per-cell geometric deposition probability for nonlocal volume. | Only used when cellular_nonlocal_fraction is greater than 0. |
| `cellular_nonlocal_max_hops` | `int` | `4` | `param` | `cellular_nonlocal_fraction` | Maximum nonlocal travel distance in cells (&gt;=2). | Only used when cellular_nonlocal_fraction is greater than 0. |

### Cellular particle tracing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `enable_cellular_particles` | `bool` | `false` | `switch` | `hillslope_model` | Track sparse Lagrangian sediment superparticles. | Only valid when hillslope_model is cellular. |
| `cellular_particle_inverse_density` | `int` | `100` | `param` | `enable_cellular_particles` | Active cells represented per seeded superparticle. | 100 is the validated standard; smaller values seed more particles. |
| `cellular_particle_mixing_depth` | `float` | `1` | `param` | `enable_cellular_particles` | Well-mixed active-layer depth represented by particles (m). | - |
| `cellular_particle_seed` | `int` | `1` | `param` | `enable_cellular_particles` | Deterministic particle sampling seed. | - |
| `cellular_particle_restart_fname` | `string` | - | `file` | `enable_cellular_particles` | Particle CSV to restart from, relative to read_path. | Empty seeds a new population. |

### roering_implicit solver knobs (used only when hillslope_model = roering_implicit)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `hillslope_implicit_relaxation_omega` | `float` | `0.7` | `param` | - | Damped-Picard relaxation factor omega (0&lt;omega&lt;=1) for roering_implicit. | Values below 1 stabilise the nonlinear iteration near Sc; the controller auto-tunes it down further on divergence. Ignored by other hillslope models. |
| `hillslope_implicit_picard_tolerance_metres` | `float` | `1e-4` | `param` | - | Convergence tolerance max\|delta_z\| for the Picard iteration (roering_implicit). | - |
| `hillslope_implicit_picard_max_iter` | `int` | `50` | `param` | - | Maximum Picard iterations per sub-step (roering_implicit). | - |
| `hillslope_implicit_slope_ratio_cap` | `float` | `0.99` | `param` | - | Regularised slope-ratio cap r_max (&lt;1) bounding D_eff near Sc (roering_implicit). | - |
| `hillslope_implicit_solver_choice` | `int` | `0` | `param` | - | Linear solver for roering_implicit: 0 = SimplicialLDLT, 1 = CG+IncompleteCholesky. | LDLT is robust for moderate grids; CG+IC scales better on very large grids. |
| `hillslope_implicit_time_integrator` | `int` | `0` | `param` | - | Time integrator for roering_implicit: 0 = backward Euler, 1 = TR-BDF2 (2nd order). | TR-BDF2 is L-stable and time-accurate with embedded PI error control; use it for time-accurate transients rather than just the approach to steady state. |
| `hillslope_implicit_error_tolerance_metres` | `float` | `0` | `param` | - | TR-BDF2 embedded-error tolerance (elevation units) for the PI step controller (0 = auto = a small fraction of the relief). Only used when time_integrator = 1. | - |
| `hillslope_implicit_dt_init` | `float` | `0` | `param` | - | Initial adaptive sub-step for roering_implicit (0 = auto = the LEM timestep). | - |
| `hillslope_implicit_dt_min` | `float` | `0` | `param` | - | Floor on the adaptive sub-step for roering_implicit (0 = auto). | - |
| `hillslope_implicit_dt_max` | `float` | `0` | `param` | - | Cap on the adaptive sub-step for roering_implicit (0 = the LEM timestep). | - |
| `channel_mask_fname` | `string` | - | `param` | - | Optional raster (in read_path) whose cells &gt; 0.5 are fixed-elevation channel pixels for the hillslope Dirichlet set, overriding the drainage-area threshold. | Must match the DEM grid. Empty (default) = derive channels from incision_threshold_specific_area. Useful for prescribed-channel / hillslope-strip runs. |
| `hold_channel_mask_fixed` | `bool` | `false` | `param` | - | Hold the channel_mask_fname pixels at fixed elevation through every step. | Requires channel_mask_fname. The masked pixels join the open-edge boundary pixels as fixed-elevation cells, so uplift does not lift the channel network away from the hillslopes relaxing toward it. Use for hillslopes-only runs on a frozen network; leave false for coupled runs where the channels must incise. |
| `write_channel_mask_raster` | `bool` | `false` | `switch` | - | Write the hillslope Dirichlet channel mask (1 = channel). | Output file: &lt;write_fname&gt;_channel_mask.tif for the final topography, plus &lt;write_fname&gt;_model_iteration_&lt;N&gt;_channel_mask.tif at every output frame when frequency_of_output_writing &gt; 0 (the migrating network of an animation). The mask is the one the hillslope step would use: the channel_mask_fname raster when supplied, else the pixels whose D-infinity drainage area reaches incision_threshold_specific_area (or channel_mask_threshold_specific_area). Feed it back in as channel_mask_fname to freeze a fluvially-evolved network for a hillslopes-only run. |
| `write_erosion_rate_raster` | `bool` | `false` | `switch` | - | Write the per-pixel erosion rate (m/yr, positive = erosion, uplift backed out). | Output file: &lt;write_fname&gt;_model_iteration_&lt;N&gt;_erosion_rate.tif at every output frame (needs frequency_of_output_writing &gt; 0), plus &lt;write_fname&gt;_erosion_rate.tif for the final step. This is the field in which a transient pulse of incision is visible: knickpoints appear as a band of high erosion rate migrating up the channel network. Deposition appears as negative values. |
| `channel_mask_threshold_specific_area` | `float` | `0` | `param` | - | A/dx threshold in metres when divide_area_by_dx=true; map area A in m^2 when false. Defines the written channel mask (0 = use incision_threshold_specific_area). | Only affects write_channel_mask_raster; the model itself is unchanged. Lets a run erode fluvially everywhere (incision_threshold_specific_area: 0) while still writing a channel network defined at a positive threshold. The run header resolves the threshold into specific area, map area, and contributing pixels. |
| `enable_uplift` | `bool` | `true` | `switch` | - | Apply uniform tectonic uplift. | - |

### Process parameters

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `uplift_rate` | `float` | `1e-4` | `param` | - | Tectonic uplift rate (m/yr). | Applied uniformly to all valid pixels. Constant unless overridden by forcing_history_fname. |
| `baselevel_fall_rate` | `float` | `0` | `param` | - | Rate at which base level falls (m/yr, positive = down). | The fixed-elevation boundary pixels are lowered at this rate instead of being held at their pre-step elevation, so a wave of incision propagates up the channel network. 0 (default) = static base level. Constant unless overridden by forcing_history_fname. |
| `forcing_history_fname` | `string` | - | `param` | - | Optional CSV (in read_path) making uplift_rate, baselevel_fall_rate, erodibility and hillslope_diffusivity functions of model time. | Columns: time_yr, uplift_rate, baselevel_fall_rate (the latter two optional, default 0; any column order), plus optional erodibility and hillslope_diffusivity, where an absent column or empty cell keeps the parameter-file value. Stepping K or D walks the channelization index; the applied K, D and CI are written to the stats CSV. Applied as a piecewise-constant step function of model time -- the row in force is the last one with time_yr &lt;= t -- so a pulse has sharp on/off edges. Overrides the constant uplift_rate and baselevel_fall_rate parameters for the whole run. Use it to impose a pulse of incision: a window of rapid base-level fall that sends knickpoints up the network, which the hillslopes then respond to. The forcing actually applied at each step is recorded in the uplift_rate / baselevel_fall_rate columns of &lt;write_fname&gt;_model_elevation_stats.csv. |
| `erodibility` | `float` | `1e-3` | `param` | - | Stream-power erodibility K (m^(1-2m) yr^-1). | Scales fluvial incision rate. |
| `m_exponent` | `float` | `0.45` | `param` | - | Drainage area exponent m in stream power law. | - |
| `n_exponent` | `float` | `1.0` | `param` | - | Slope exponent n in stream power law. | n=1 uses closed-form solver; n!=1 uses Newton-Raphson. |
| `hillslope_diffusivity` | `float` | `1e-3` | `param` | - | Hillslope diffusivity kappa (m^2/yr). | Used when enable_hillslope_transport is true. |
| `critical_hillslope_slope` | `float` | `0.8` | `param` | - | Critical slope Sc (m/m) for nonlinear hillslope transport. | - |
| `incision_threshold_specific_area` | `float` | `0` | `param` | - | A/dx channel threshold in metres when divide_area_by_dx=true; map area A in m^2 when false. 0 = no threshold. | With divide_area_by_dx true (the DEFAULT) the threshold is compared against specific catchment area a = A/dx, a length in metres; with it false, against map area A in m^2. Two consequences worth knowing before choosing a value: (1) under the default, threshold T selects cells with at least T/dx upstream pixels, so the same .param file carried from a 10 m DEM to a 1 m DEM asks for a ten times larger network in pixel terms; (2) a value picked as though it were m^2 is wrong by a factor of dx - on a 10 m grid '1000' means 100 pixels, not the 10 that 1000 m^2 would give. The run header prints the resolved threshold in all three forms. The specific-area convention follows Bonetti et al. (2020) and Litwin et al. (2022) and exists to remove grid dependence from the erosion law; only its documentation was ever wrong. Also defines channel pixels for hillslope coupling. |
| `hillslope_channel_coupling` | `string` | `auto` | `param` | - | Which channel/hillslope coupling to run: auto, coupled, or channels_as_boundaries. | lsdtt3 contains BOTH standard couplings, and which one you get has until now been an unreported side effect of incision_threshold_specific_area. channels_as_boundaries makes channels FIXED (Dirichlet) nodes for the hillslope solve, so drainage density is an INPUT - the MuddPILE configuration, D8 routing. coupled gives no Dirichlet set: both processes act everywhere and valley spacing, hillslope length and drainage density EMERGE - the DonkModel configuration, which needs D-infinity because D8 on a diverging hillslope manufactures channel heads. These are different experiments, not different settings. 'auto' reproduces the historical inference exactly and only reports it, and is the default because parameter files across this project already depend on the implicit behaviour. NOTE roering_implicit cannot run the coupled regime: it needs a Dirichlet set, so the most robust hillslope solver is the one that cannot do the emergent-drainage-density experiment. |
| `divide_area_by_dx` | `bool` | `true` | `param` | - | Divide drainage area by cell size to get specific area a = A/dx. | Applies to both fluvial incision and the hillslope coupling threshold, and changes the interpretation of incision_threshold_specific_area from map area in m^2 to specific area in metres. True is the Bonetti et al. (2020) / Litwin et al. (2022) convention (specific area, grid-independent); false is the Perron et al. (2008) total-area convention. The choice also changes the exponent in the landscape Peclet number, so it is not a free setting. |
| `slope_clamp` | `float` | `1e-5` | `param` | - | Minimum slope used as Newton-Raphson fallback in stream power solver. | - |

### Analytical steady-state

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_steady_state` | `bool` | `false` | `switch` | - | Compute analytical steady-state topography instead of running a time-stepping LEM. | Uses the flow algorithm selected by enable_d8_fluvial_erosion or enable_dinf_fluvial_erosion. Writes _steady_state_elevation.tif and returns. |
| `steady_state_mode` | `string` | `K` | `param` | `run_steady_state` | Steady-state solver mode: K or relief. | K-mode: K (erodibility) is specified by the user. Relief-mode: target_relief is specified and K is back-calculated. |
| `target_relief` | `float` | `1000` | `param` | `run_steady_state` | Target maximum relief (m) for relief-mode steady-state. | The highest point in the output raster will be baselevel_elevation + target_relief. |
| `baselevel_elevation` | `float` | `0` | `param` | `run_steady_state` | Elevation (m) assigned to outlet nodes in the steady-state solution. | - |
| `steady_state_n_cycles` | `int` | `1` | `param` | `run_steady_state` | Number of snapping cycles for steady-state computation. | Each cycle recomputes the flow graph from the previous cycle's elevation, allowing drainage divides to migrate one pixel per cycle. Iterating (typically 5-20 cycles) smooths the sharp elevation breaks across divides that appear in a single-pass steady-state snap. The DEM is pit-filled between cycles. |
| `use_landscape_equilibration_cycle` | `bool` | `false` | `param` | `run_steady_state` | Route analytical steady-state through the landscape equilibration cycle module. | Default false keeps legacy steady-state behavior unchanged. |
| `divide_merge_search_steps` | `int` | `5` | `param` | `use_landscape_equilibration_cycle` | Downstream search depth used by divide diagnostics in equilibration-cycle runs. | Must be &gt; 0 when use_landscape_equilibration_cycle is true. |
| `divide_use_d8_diagonals` | `bool` | `true` | `param` | `use_landscape_equilibration_cycle` | Include diagonal neighbors in divide-step diagnostics (D8 adjacency). | Only used by equilibration-cycle diagnostics outputs/convergence checks. |
| `steady_state_check_convergence` | `bool` | `false` | `param` | `use_landscape_equilibration_cycle` | Enable convergence-based early stopping in equilibration-cycle steady-state runs. | Default false keeps current behavior (always runs steady_state_n_cycles). |
| `steady_state_tolerance_dz_over_relief` | `float` | `1e-3` | `param` | `steady_state_check_convergence` | Convergence tolerance for max\|dz\|/relief between equilibration cycles. | Used when steady_state_check_convergence is true and steady_state_tolerance_dz_metres is 0. |
| `steady_state_tolerance_dz_metres` | `float` | `0` | `param` | `steady_state_check_convergence` | Absolute convergence tolerance for max\|dz\| (m) between equilibration cycles. | If &gt; 0 this overrides steady_state_tolerance_dz_over_relief. |
| `steady_state_tolerance_receiver_fraction` | `float` | `1e-3` | `param` | `steady_state_check_convergence` | Convergence tolerance for fraction of pixels changing receiver each cycle. | Set &lt;= 0 to disable this convergence criterion. |
| `divide_gradient_target` | `float` | `0` | `param` | `use_landscape_equilibration_cycle` | Optional divide-gradient acceptance target for equilibration diagnostics. | Reporting-only threshold: when &gt; 0, the final cycle reports whether divide gradient median &lt;= this target; it does not change stopping. |
| `steady_state_convergence_patience` | `int` | `2` | `param` | `steady_state_check_convergence` | Required number of consecutive converged cycles before early stopping. | Only used when steady_state_check_convergence is true. |

### Hillslope steady-state

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_hillslope_on_loaded_dem` | `bool` | `false` | `switch` | - | Solve hillslopes on the loaded DEM without re-snapping the fluvial field. | Requires read_initial_dem: true, run_hillslope_steady_state: true and incision_threshold_specific_area &gt; 0. Channels at or above the resolved threshold keep the elevations they arrived with and act as Dirichlet boundaries; only hillslope pixels are solved. Use this when the fluvial field is a TRANSIENT that must be preserved -- a step-change landscape carries a migrating knickpoint, and run_steady_state would snap the channels to a single ksn and erase it. Exits after writing &lt;write_fname&gt;_hillslope_on_loaded_dem. |
| `run_hillslope_steady_state` | `bool` | `false` | `switch` | - | Solve for hillslope steady-state elevations after fluvial steady-state. | Solver is selected by hillslope_ss_model: linear, roering, or cidre. Requires run_steady_state: true and incision_threshold_specific_area &gt; 0. Channel pixels at or above the resolved threshold become Dirichlet BCs. |
| `hillslope_ss_tolerance` | `float` | `1e-4` | `param` | `run_hillslope_steady_state` | Convergence criterion for hillslope steady-state: max\|delta_z\| over hillslope pixels. | - |
| `hillslope_ss_max_iterations` | `int` | `10000` | `param` | `run_hillslope_steady_state` | Maximum SOR iterations for hillslope steady-state solver. | - |
| `hillslope_ss_model` | `string` | `roering` | `param` | `run_hillslope_steady_state` | Hillslope steady-state solver: linear, roering, cidre, or critical_slope. | linear: Red-Black SOR solving D*grad^2z = -U (standard Laplacian diffusion). roering: Red-Black SOR solving D*grad*[gradz/(1-(\|gradz\|/Sc)^2)] = -U (nonlinear). Both SOR solvers converge in O(N) iterations with auto-tuned omega. critical_slope: imposes Sc by upslope integration from the channels -- not an elliptic solve, so it never finds a divide position, it inherits one; only available with run_hillslope_on_loaded_dem. cidre: pseudo-transient Cidre (Carretier 2016) iterated until steady state. Rebuilds the D-inf flow graph each iteration (Picard update of routing directions). The nonlinear flag for cidre is taken from enable_nonlinear_hillslopes. At steady state the nonlinear Cidre and Roering solvers produce the same elevation field (Carretier et al. 2016 Sec3.3). |
| `hillslope_every_cycle` | `bool` | `true` | `param` | `run_hillslope_steady_state` | Solve the hillslopes on every equilibration cycle, not just the last one. | Only used when use_landscape_equilibration_cycle and run_hillslope_steady_state are both true. Solving every cycle costs one extra solve per cycle but is the only way to show the divide is stable under re-routing; false solves once on the delivered surface. |
| `channel_hysteresis_factor` | `float` | `0.5` | `param` | `run_hillslope_steady_state` | Hysteresis on the channel (Dirichlet) set between equilibration cycles. | A pixel already flagged as channel stays channel until its drainage area falls below channel_hysteresis_factor * incision_threshold_specific_area. Stops the Dirichlet set flickering from cycle to cycle. 1.0 disables the hysteresis. |
| `hillslope_sor_omega` | `float` | `0` | `param` | `run_hillslope_steady_state` | SOR relaxation factor omega for hillslope steady-state (0 = auto). | Values in (1, 2) accelerate convergence. 0 = auto-estimate from grid size (capped at 1.95). Typical manual values: 1.5 (safe) to 1.9 (fast). Reduce toward 1.3 if the nonlinear solver diverges. |

### Hillshade output (steady-state runs only)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_hillshade_png` | `bool` | `false` | `switch` | - | Write a hillshade PNG of the final steady-state elevation. | Computes a standard Horn (1981) hillshade and saves it as a PNG alongside the elevation output. If run_hillslope_steady_state is also true, a second hillshade of the combined fluvial+hillslope surface is written. PNG is always used regardless of write_output_format. |
| `hillshade_azimuth` | `float` | `315` | `param` | `write_hillshade_png` | Sun azimuth for hillshade, degrees clockwise from north. | 315 = northwest (ESRI/QGIS default). |
| `hillshade_altitude` | `float` | `45` | `param` | `write_hillshade_png` | Sun altitude above the horizon for hillshade, degrees. | 45 is the standard default. |
| `hillshade_z_factor` | `float` | `1` | `param` | `write_hillshade_png` | Vertical exaggeration factor for hillshade. | Use values &gt; 1 to emphasise relief on low-gradient landscapes. |

### Timestepping

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `timestep` | `float` | `100` | `param` | - | Model timestep (yr). | Stability depends on erodibility, resolution and slope. |
| `number_of_iterations` | `int` | `100` | `param` | - | Number of timesteps to run. | Ignored when stop_at_steady_state is true. |
| `snap_initial_surface_to_steady_state` | `bool` | `true` | `switch` | - | Before time-stepping, snap the initial surface to the analytical fluvial steady state (chi-based, K-mode). | A near-flat start cannot reach steady state sooner than R/U, the time uplift needs to build the final relief R, so a run from noise spends at least that long just gaining height (0.8 Myr for 165 m at 0.2 mm/yr, and ~3 Myr in practice). The snap sets z = baselevel + (U/K)^(1/n) chi on the pit-filled initial surface, routed with the same flow method as the fluvial step (D-infinity when enable_dinf_fluvial_erosion, else D8), so time-stepping starts from the right relief and only has to equilibrate hillslopes and divides. Uses uplift_rate, erodibility, m_exponent, n_exponent, baselevel_elevation and divide_area_by_dx; with a forcing history the uplift rate at t = 0 is used. Skipped, with a message, when no fluvial process is enabled or the uplift rate is not positive. Set false to keep the initial surface as given -- REQUIRED for a transient experiment on a loaded DEM, whose topography this replaces. |
| `stop_at_steady_state` | `bool` | `false` | `param` | - | Stop when mean erosion rate is within 1% of uplift rate for 100 consecutive steps. | After detection, runs 20x that step count more before stopping. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `frequency_of_output_writing` | `int` | `0` | `param` | `write_elevation_stats_csv` | Write an intermediate raster every N iterations. 0 = write final only. | - |
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format: cog, gtiff, same_as_input, png. | png: 8-bit greyscale PNG for quick visualisation. Float32 elevation is normalised to 1-255 (nodata-&gt;0). A worldfile (.pgw) and GDAL PAM sidecar (.png.aux.xml) are written alongside for georeferencing. PNG cannot store float32 values or embed a CRS natively, so use cog/gtiff for analysis. |
| `write_elevation_stats_csv` | `bool` | `true` | `switch` | - | Write per-timestep elevation statistics to a CSV file. | - |
| `write_snap_convergence_csv` | `bool` | `true` | `switch` | - | Write per-cycle snap convergence diagnostics for the landscape equilibration path. | Only used when run_steady_state and use_landscape_equilibration_cycle are both true. |
| `write_divide_step_raster` | `bool` | `false` | `switch` | - | Write a divide-step raster from the final equilibration-cycle steady-state DEM. | Only used when run_steady_state and use_landscape_equilibration_cycle are both true. |

### Routing and threading

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Flow routing boundary conditions (north/east/south/west). | o=open, n=no-flux, p=periodic, b=forced outlet. |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. 0 = all available. | - |

<!-- END GENERATED CLI PARAMETERS -->
