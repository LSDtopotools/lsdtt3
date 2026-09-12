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
