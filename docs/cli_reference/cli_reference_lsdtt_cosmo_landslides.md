# CLI reference: lsdtt-cosmo-landslides

> **Renamed 2026-08-05.** This tool was `lsdtt-landslides`. There is no
> compatibility alias — parameter files are unaffected, but scripts invoking the
> old executable name must be updated. The rename marks it as a cosmogenic tool
> rather than a purely geomorphic one: it answers *"could landsliding be biasing
> my cosmogenic samples?"*, and is gaining full production scaling and shielding
> so the concentrations it reports are real rather than indicative.

Triggers landslides with an **imposed** size distribution. You prescribe *when*
(a trigger time, or a recurrence) and *how much* (an event volume, or a
long-term erosion rate); the tool works out how many landslides that implies,
what sizes they are, and where they go, then excavates them and writes the
fields the cosmogenic and sediment-routing tools consume.

The DEM must be in a projected CRS — areas and volumes come from cell size in
metres.

## Usage

```
lsdtt-cosmo-landslides <path> <parameter_file>
lsdtt-cosmo-landslides <parameter_file>         # path defaults to ./
lsdtt-cosmo-landslides -h                       # print help and exit
lsdtt-cosmo-landslides -v                       # print version and exit
```

---

## What this tool does and does not do

**Sizes are imposed, not emergent.** Landslide areas are sampled from a
prescribed distribution rather than derived from slope stability. That is
deliberate: the purpose is to impose a *known* landslide flux and ask what it
does to something else — most immediately to cosmogenic nuclide concentrations,
which is what Niemi et al. (2005) and Yanites et al. (2009) both study.

The alternative — deriving sizes from Culmann stability, as HyLands does
(Campforts et al. 2020) — is more fundamental and reproduces the observed
magnitude–frequency distribution *without being told to*. It is the right tool
if you want the distribution to be an outcome. It is the wrong tool if you want
the landslide flux to be a controlled input, because the flux then depends on
rock-strength calibration.

**Landslides come from hillslopes.** Initiation is restricted to pixels whose
contributing area is below the channel initiation threshold, and each scar is
built from that pixel's upslope set. Two consequences follow by construction, not
by filtering: no scar can reach into a channel, and **no landslide can be larger
than the hillslope it comes from**.

**No runout, no deposition.** Scars are excavated; where the material goes is the
sediment-routing tools' problem.

---

## Pipeline

| Step | Description |
|---|---|
| 1 | Load DEM; check projected CRS |
| 2 | Pit fill (skipped when `fill_before_routing: false`) |
| 3 | D8 flow routing |
| 4 | Build the discrete size distribution over pixel counts |
| 5 | Build the eligible-host list (hillslope pixels, inside the mask) |
| 6 | For each event: solve for the landslide count against the volume budget, place largest first, excavate |
| 7 | Write excavation depth, time of failure, inventory CSV, optional CRN concentration |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | `cog`, `gtiff`, or `same_as_input`. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fill_before_routing` | bool | `true` | Fill pits before building flow routing. |
| `min_slope_for_fill` | float | `0.0001` | Minimum imposed slope during filling (m/m). |
| `boundary_conditions` | string | `oooo` | D8 boundary conditions (N,E,S,W): `o` = open, `b` = forced outlet. |

---

## Placement

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold_contributing_pixels` | int | `1000` | Channel initiation threshold, in pixels. Landslides start only **below** this. It also caps the largest landslide. |
| `mask_fname` | string | `NULL` | Optional raster restricting where landslides may occur. Non-zero and valid = active, unless `mask_active_values` narrows it. |
| `mask_active_values` | string | `` | Comma-separated list of mask values to treat as active (e.g. `4` or `3,4`). Empty = historical rule: any non-zero valid pixel is active. |
| `ice_mask_fname` | string | `NULL` | Raster marking ice cover (non-zero = ice). No landslide is placed on ice. |

`ice_mask_fname` is **intersected** with `mask_fname` when both are given: a
pixel can host a landslide only if it is inside the active mask *and* not under
ice. The number of pixels excluded is reported, because it changes the **active
area** — and the active area sets the landslide count for a prescribed erosion
rate, so two runs differing only in an ice mask are not comparable by count.

**`mask_active_values`: selecting a subset of a classified raster.** By default,
any non-zero valid pixel in the mask raster is active. `mask_active_values` allows
you to specify which values count as active, so a pre-classified raster can be used
directly without reclassifying to 0/1 first. **Example use case:** the step-change
tag raster from `lsdtt-cosmo-prediction` has four classes: 1 = fluvial-old,
2 = hillslope-old, 3 = fluvial-new, 4 = hillslope-new. Landslides only on
ground that has already experienced the faster post-knickpoint erosion means
`mask_fname: _SC_TAG.tif` and `mask_active_values: 4`. Without
`mask_active_values`, all four classes would be active, placing landslides
upstream of the knickpoint — silently inverting the intended scenario. The values
are matched to the nearest integer, so this works for integer class/tag rasters,
not continuous ones. The realised active-pixel count is always reported.

Initiation pixels are drawn uniformly from the eligible set. This is what Niemi
and Yanites both do, and both flag it as a simplification — real triggered
inventories cluster strongly.

**A size-dependent bias to be aware of.** Hosting an *n*-pixel landslide requires
*n* pixels upslope, and contributing area is itself power-law distributed, so the
eligible set shrinks rapidly with size. The largest landslides can only initiate
just above channel heads and will occupy most of a zero-order basin. That is
defensible — large landslides really do evacuate whole zero-order basins — but it
is imposed rather than chosen.

---

## Size distribution

| Parameter | Type | Default | Description |
|---|---|---|---|
| `landslide_size_model` | string | `inverse_gamma` | `inverse_gamma` (Malamud et al. 2004) or `truncated_pareto` (Yanites et al. 2009). |
| `inverse_gamma_rho` | float | `1.40` | ρ. Asymptotic tail exponent is −(ρ+1) = −2.40. |
| `inverse_gamma_a_m2` | float | `1280` | a (m²). Sets the modal size. |
| `inverse_gamma_s_m2` | float | `-132` | s (m²), **negative**. Sets the small-area rollover. |
| `pareto_alpha_pdf` | float | `2.3` | Truncated Pareto exponent — **PDF, not cumulative**. |
| `min_landslide_pixels` | int | `1` | Smallest landslide, in pixels. |
| `max_landslide_area_m2` | float | `0` | Largest landslide. `0` = use the channel threshold. May not exceed it. |
| `area_volume_preset` | string | `larsen_global` | `V = αA^γ` preset. See below. |

The inverse gamma is the default because it has both a power-law tail *and* a
small-area rollover. Malamud argues the rollover is real rather than a mapping
artefact, and HyLands independently produces one from mechanics.

> ### ⚠ Exponent convention
>
> `pareto_alpha_pdf` is the **probability-density** exponent. Niemi's `b = 1.18`
> is a **cumulative** exponent; the two are related by `α_pdf = b_cum + 1`, so
> Niemi corresponds to `pareto_alpha_pdf: 2.18`. Valagussa et al. (2019) measure
> 1.75–3.11 across lithologies for coseismic inventories.

### Area-to-volume presets

From Larsen, Montgomery & Korup (2010) Supplementary Table S1. A in m², V in m³.

| Preset | γ | α | n |
|---|---|---|---|
| `larsen_global` (default) | 1.332 | 0.1459 | 4231 |
| `larsen_soil_scar` | 1.262 | 0.2244 | 1617 |
| `larsen_bedrock_scar` | 1.410 | 0.2344 | 168 |
| `larsen_soil_all` | 1.145 | 0.3631 | 2136 |
| `larsen_bedrock_all` | 1.350 | 0.1862 | 604 |
| `hovius_slab` | 1.5 | 0.05 | — |
| `hovius_parabolic` | 1.5 | 0.025 | — |

**Prefer the scar presets.** This tool excavates scars and does not model
deposits. The commonly quoted γ_soil = 1.145 and γ_bedrock = 1.35 are *combined*
scar + deposit fits; the scar-only rows give 1.262 and 1.410. Using
`larsen_soil_all` for a scar model understates large soil landslides by a factor
that grows with area.

The `hovius_*` presets exist to reproduce Niemi, Yanites and Malamud, all of
which use γ = 1.5. Larsen shows γ = 1.5 **overestimates total volume by at least
a factor of two** versus γ = 1.4, and since the slab form is itself twice the
parabolic form at equal ε, a careless choice can be 4× out.

### Resolution limit

With `min_landslide_pixels: 1` the smallest representable landslide is one cell.
At 30 m that is 900 m², **above** Malamud's modal area of ~400 m². At that
resolution the rollover cannot be represented and the sampled distribution is
effectively a truncated power law whatever family you select. The tool warns:

```
WARNING: cell size cannot resolve the small-area rollover.
```

---

## Triggering

| Parameter | Type | Default | Description |
|---|---|---|---|
| `landslide_mode` | string | `single_event` | `single_event` or `series`. |
| `event_time_yr` | float | `0` | Trigger time for `single_event`. |
| `event_volume_m3` | float | `0` | Total volume to remove per event. |
| `event_volume_per_area_metres` | float | `0` | Event volume per unit **active** area (m³/m² = metres of lowering). Used when `event_volume_m3` is 0. |
| `series_start_yr` | float | `0` | Series start. |
| `series_end_yr` | float | `100000` | Series end. |
| `mean_recurrence_yr` | float | `1000` | Mean interval between events. |
| `poisson_timing` | bool | `true` | Exponential inter-event intervals (Yanites Eq. 14). `false` = exactly regular. |
| `long_term_erosion_rate_mm_kyr` | float | `0` | Landslide erosion rate over the active area. The Niemi convention: prescribe the rate, let the model find the events. Series mode only. |
| `anneal_time_yr` | float | `0` | A failed pixel cannot fail again for this long. |
| `budget_tolerance` | float | `0.02` | Relative tolerance on the landslide-count solve. |
| `budget_max_iterations` | int | `50` | Iteration cap. |
| `random_seed` | int | `1` | Same seed reproduces a run exactly. |

"Active area" is the ground that can host a landslide, **not** the whole raster.
Using the raster area would understate the rate on a channel-dominated or masked
DEM.

`anneal_time_yr` stands in for regolith regeneration and is **uncalibrated** —
there is no strong literature constraint on it.

Poisson timing follows Yanites et al. (2009), who are explicit that it is a
simplification: real landslides are triggered by storms and earthquakes and are
neither independent nor uncorrelated.

---

## Three things that will bite you

These are properties of a heavy-tailed size distribution, not defects, and all
three are reported at run time.

### 1. The event volume must be much larger than the mean landslide volume

Solving for the landslide count conditions the draw on its total volume. With
many landslides that conditioning is weak. With few it is severe, and at one
landslide per event it is absolute — the tool is no longer sampling the
distribution, it is searching for the single landslide whose volume equals the
budget, so realised sizes cluster at `(V_target/α)^(1/γ)`.

The run reports `E[V]`; **keep the per-event volume at tens of times that.**
Below ten landslides per event the tool warns:

```
WARNING: only N landslide(s) in this event. Solving for N against a small
budget SELECTS sizes rather than sampling them...
```

### 2. `A_max` barely truncates the distribution but dominates the volume

At a 1000-pixel threshold on a 30 m DEM, `A_max` = 0.9 km² removes **0.008 % of
landslides by number** and **59 % of the volume**. The tail is `p(A) ∝ A^−2.40`
and volume goes as `A^1.332`, so volume density goes as `A^−1.068` — each
logarithmic decade of area contributes roughly equally to the total.

The practical consequence: **`E[V]`, and hence the landslide count for a fixed
budget, is strongly sensitive to `threshold_contributing_pixels`.** Halving it
changes the count by tens of per cent. Both truncated fractions are echoed in
the run output.

### 3. Budget tolerance is unreachable when one landslide can carry the budget

A single landslide at `A_max` carries `α·A_max^γ`. When that is a large fraction
of the event volume, no choice of landslide count can land within a tight
tolerance. Roughly, `budget_tolerance` is reachable when
`event_volume_m3 >> α·A_max^γ`. The alternative — trimming the last landslide to
fit — would buy exactness by destroying the distribution at the large end, so it
is not offered.

---

## Cosmogenic coupling

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_crn_concentration` | bool | `false` | Compute the post-landslide surface concentration raster. |
| `run_sample_point_ensemble` | bool | `false` | Run a multi-realisation sample-point ensemble and write `<write_fname>_LS_SAMPLE_EFFECT.csv`. |
| `sample_points_fname` | string | `NULL` | Input sample-point CSV for ensemble mode. |
| `sample_name_column` | string | `sample_name` | Name column in `sample_points_fname`. |
| `sample_latitude_column` | string | `latitude` | Latitude column in `sample_points_fname` (EPSG:4326 degrees). |
| `sample_longitude_column` | string | `longitude` | Longitude column in `sample_points_fname` (EPSG:4326 degrees). |
| `n_realisations` | int | `50` | Number of realisations for ensemble mode. |
| `ensemble_seed_stride` | int | `7919` | Seed increment between realisations (`random_seed + i*stride`). |
| `sample_snap_max_search_radius_pixels` | int | `25` | Maximum snap radius when mapping sample points to threshold channels. |
| `nuclide` | string | `Be10` | `Be10`, `Al26` or `C14`. |
| `muon_scheme` | string | `BraucherBorchers` | Muon production scheme. |
| `scaling_scheme` | string | `uniform` | Production scaling path for CRN coupling: `uniform`, `stone`, or `lsdn`. `uniform` keeps the legacy scalar path (`production_scaling` and `total_shielding`). `stone`/`lsdn` build a per-pixel `CosmoSite` from the DEM. |
| `atmosphere` | string | `ERA40` | Atmosphere model for site-based scaling (`stone`/`lsdn`): `ERA40`, `NCEP`, `standard`. |
| `crn_nc_fname` | string | `""` | Path to `crn_parameters.nc` for site-based scaling and atmosphere handling. Empty uses bundled lookup/fallback logic. |
| `sample_age_yr` | float | `0` | Sample exposure age (yr BP) for time-integrated production correction in site-based scaling. |
| `background_eff_erate_g_cm2_yr` | float | `0.001` | Background effective erosion rate before and after failure. |
| `production_scaling` | float | `1.0` | Uniform production scaling factor. Used directly when `scaling_scheme: uniform`. |
| `total_shielding` | float | `1.0` | Uniform total shielding factor. Used directly when `scaling_scheme: uniform`; for site-based scaling use `topo_shield_value`, `snow_shield_eff_thickness`, `self_shield_eff_thickness`. |
| `topo_shield_value` | float | `1.0` | Uniform topographic shielding scalar (0-1) for site-based scaling. |
| `snow_shield_eff_thickness` | float | `0` | Uniform snow shielding as effective depth (g/cm2) for site-based scaling. |
| `self_shield_eff_thickness` | float | `0` | Uniform self/sample shielding as effective depth (g/cm2) for site-based scaling. |
| `rock_density_kg_m3` | float | `2650` | Density of the **eroding bedrock column**. |
| `scar_density_kg_m3` | float | `-1` | Density of the material a landslide **excavates**. The negative sentinel deliberately inherits the active `area_volume_preset`; prediction from an already-generated excavation raster instead defaults explicitly to 1500 kg/m3. |
| `crn_report_time_yr` | float | `0` | Time at which concentrations are reported, on the same clock as event times. |

`cronusEmulatorSt` and `cronusEmulatorLSDn` are intentionally rejected in this
tool's landslide CRN path. The lower-level step-change/transient kernels used
here do not yet support the emulator parameterisation consistently.

### Two densities, and why they are not the same number

`rock_density_kg_m3` is the eroding column; `scar_density_kg_m3` is the material
that failed. They are different quantities and conflating them is a real error,
not a rounding matter.

The area-to-volume scaling `V = αA^γ` is fitted to **measured** scar volumes, so
the excavation depth it implies is a true geometric thickness of whatever
failed — and for a shallow landslide that is soil, at roughly 1500 kg/m³.
Carrying it at 2650 overstates the mass depth by a factor of **1.77**, and since
production falls exponentially with mass depth that changes production at the
scar base by **2.05**.

The density therefore travels with the preset, because choosing the preset is
already a statement about the material:

| Preset | γ | Scar density | Material-specific? |
|---|---|---|---|
| `larsen_soil_scar` | 1.262 | 1500 | yes |
| `larsen_bedrock_scar` | 1.410 | 2650 | yes |
| `larsen_soil_all` | 1.145 | 1500 | yes |
| `larsen_bedrock_all` | 1.350 | 2650 | yes |
| `larsen_global` *(default)* | 1.332 | 2650 *placeholder* | **no — warns** |
| `hovius_slab` / `hovius_parabolic` | 1.5 | 2650 *placeholder* | **no — warns** |

`larsen_global` and the Hovius rows are fitted across materials and have **no
defensible density**, so leaving `scar_density_kg_m3` negative with one of those
selected raises a warning. The default preset is one of them, which is
deliberate: that is the case that would otherwise pass silently.

No new cosmogenic physics is introduced. The model composes two existing
validated calls: **inheritance** — a landslide stripping effective depth *d*
exposes rock that was at depth *d*, so its concentration is the steady-state
value there — and **re-accumulation** since failure, using the transient
integrator's instantaneous *removal* event, which is precisely a landslide.

The time-of-failure raster is both the annealing clock and the
time-since-change field, so it feeds `lsdtt-cosmo-prediction` directly.

### Feeding the landslide field into a step change

`lsdtt-cosmo-prediction` reads `<stem>_LS_TFAIL` and `<stem>_LS_EXCAVATION` and
superimposes the scars on a step-change concentration field, in either of two
idealised scenarios — `scar` (the debris is long gone, what remains is a
stripped column) or `flush` (the debris is still in the channel, so the sample
contains the whole excavated slab), or `both` split by age. See
[the landslide coupling section](cli_reference_lsdtt_cosmo_prediction.md#landslide-coupling).

Two conventions to match across the two tools:

- **the failure-time clock.** `_LS_TFAIL` holds an absolute time on whatever
  clock `event_time_yr` / the series uses. `lsdtt-cosmo-prediction` converts it
  to elapsed years with `landslide_report_time_yr`, so that must be "now" on the
  same clock.
- **the scar density.** Set `scar_density_kg_m3` consistently in both, or the
  excavated slab is carried at one density here and another there.

---

## Output

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_excavation_depth_raster` | bool | `true` | `<write_fname>_LS_EXCAVATION.tif` — cumulative depth removed (m). |
| `write_time_of_failure_raster` | bool | `true` | `<write_fname>_LS_TFAIL.tif` — time each pixel last failed. |
| `write_post_landslide_dem_raster` | bool | `false` | `<write_fname>_LS_DEM.tif` — the DEM after excavation. |
| `write_landslide_csv` | bool | `true` | `<write_fname>_LS_INVENTORY.csv` — one row per landslide. |
| `run_sample_point_ensemble` | bool | `false` | `<write_fname>_LS_SAMPLE_EFFECT.csv` — per-sample control vs landslide apparent-rate ensemble statistics. |

The inventory CSV has `landslide_id`, `n_pixels`, `area_m2`, `volume_m3`,
`mean_depth_m`. It is the realised size distribution, for checking against the
requested one.

### Checking the realised distribution

```
python3 tests/scripts/landslide_distribution_acceptance.py \
  --inventory <write_path>/<write_fname>_LS_INVENTORY.csv \
  --cell-size 10 --out-csv acceptance.csv
```

> **Do not compare the fitted tail exponent with ρ+1.** The inverse gamma reaches
> its asymptotic slope only for A ≫ a: the local log-log slope is
> `−(ρ+1) + a/(A−s)`, which for the Malamud parameters is −0.38 at the mode,
> −2.15 at 5,000 m², and −2.40 only around 10⁶ m². Since `A_max` is capped by
> hillslope size, a realistic run has **no range over which the asymptote is
> measurable** — a correct inventory fits to about 1.73 from the mode. The script
> compares against the exponent the target distribution itself gives over the
> same fit range, which is the meaningful comparison.

---

## Example: one triggered event

```
read_path:  /data/dems/
read_fname: catchment
write_path: /data/output/
write_fname: storm_event

threshold_contributing_pixels: 1000

landslide_size_model: inverse_gamma
area_volume_preset:   larsen_bedrock_scar

landslide_mode:  single_event
event_time_yr:   0
event_volume_m3: 5000000
random_seed:     1
```

## Example: a long series at a prescribed erosion rate

The Niemi convention — prescribe the rate, let the model find the events.

```
read_path:  /data/dems/
read_fname: catchment
write_path: /data/output/
write_fname: series_50

threshold_contributing_pixels: 1000
area_volume_preset:            larsen_bedrock_scar

landslide_mode:                series
series_start_yr:               0
series_end_yr:                 100000
mean_recurrence_yr:            1000
poisson_timing:                true
long_term_erosion_rate_mm_kyr: 50
anneal_time_yr:                2000

calculate_crn_concentration: true
nuclide:                     Be10
production_scaling:          4.0
crn_report_time_yr:          100000
```

Validated against Niemi et al. (2005): across 20–200 mm/kyr the realised
long-term rate returns the prescribed rate to within 1.9 %, against the ~3 %
they report for their own model.

---

## Sources

- **Malamud, B.D., Turcotte, D.L., Guzzetti, F. & Reichenbach, P. (2004)**
  "Landslide inventories and their statistical properties", *ESPL* 29, 687–711.
  The inverse-gamma distribution and its parameters.
- **Larsen, I.J., Montgomery, D.R. & Korup, O. (2010)** "Landslide erosion
  controlled by hillslope material", *Nature Geoscience* 3, 247–251. Area–volume
  scaling, including the scar/deposit split.
- **Niemi, N.A. et al. (2005)** "Effects of bedrock landslides on cosmogenically
  determined erosion rates", *EPSL* 237, 480–498. Prescribed-rate convention.
- **Yanites, B.J., Tucker, G.E. & Anderson, R.S. (2009)** *JGR-ES* 114, F01007.
  Truncated Pareto, Poisson timing, the CRN question.
- **Valagussa, A. et al. (2019)** *EPSL* 506, 268–281. How the exponent varies
  with trigger intensity and lithology.
- **Campforts, B. et al. (2020)** "HyLands 1.0", *GMD* 13, 3863–3886. The
  emergent alternative.

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

**Status: experimental.** This CLI is not yet fully tested.

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_before_routing` | `bool` | `true` | `switch` | - | Choose whether to pit-fill the DEM before building the flow network. | ASKS: should the DEM be filled at all? Filling is OPTIONAL here and this is a genuine analytical choice. True (the default) routes over a pit-free surface. False builds the flow network on the UNFILLED surface, so flow paths terminate at sinks instead of being carried through them. Do not confuse this with raster_is_filled, which is an assertion that the input is already filled. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_slope_for_fill` | `float` | `0.0001` | `param` | - | Minimum imposed slope during pit filling (m/m). | - |
| `boundary_conditions` | `string` | `oooo` | `param` | - | D8 flow routing boundary conditions. | Four-character string (N,E,S,W): o=open, b=forced outlet. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `param` | - | Raster output format: cog, gtiff, or same_as_input. | - |
| `write_excavation_depth_raster` | `bool` | `true` | `switch` | - | Write cumulative excavation depth (m). | Filename: &lt;write_fname&gt;_LS_EXCAVATION.tif |
| `write_time_of_failure_raster` | `bool` | `true` | `switch` | - | Write the time each pixel last failed. | This is the annealing clock AND the cosmogenic time-since-change field. Filename: &lt;write_fname&gt;_LS_TFAIL.tif |
| `write_post_landslide_dem_raster` | `bool` | `false` | `switch` | - | Write the DEM after excavation. | Filename: &lt;write_fname&gt;_LS_DEM.tif |
| `write_landslide_csv` | `bool` | `true` | `switch` | - | Write one row per landslide (size in pixels, area, volume). | The realised size distribution, for comparison against the requested one. Filename: &lt;write_fname&gt;_LS_INVENTORY.csv |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Placement

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `mask_fname` | `string` | `NULL` | `param` | - | Optional raster restricting where landslides may occur. | Non-zero and valid = active, unless mask_active_values narrows it. NULL uses the whole DEM. |
| `mask_active_values` | `string` | - | `param` | - | Comma-separated list of mask values to treat as active, e.g. "4" or "3,4". Empty (the default) keeps the historical rule: any non-zero valid pixel is active. | This exists so a CLASSIFIED raster can be used directly as a mask instead of having to be reclassified to 0/1 first. The motivating case is the step-change tag raster from lsdtt-cosmo-prediction (&lt;stem&gt;_SC_TAG.tif: 1 = fluvial-old, 2 = hillslope-old, 3 = fluvial-new, 4 = hillslope-new), where landsliding only on ground that has already felt the faster post-knickpoint erosion means mask_active_values: 4. Passing that raster without this parameter would make ALL FOUR classes active, putting landslides upstream of the knickpoint and quietly inverting the result. Values are matched to the nearest integer, so this is for class/tag rasters, not continuous ones. The realised host-pixel count is always reported -- check it. |
| `ice_mask_fname` | `string` | `NULL` | `param` | - | Raster marking ice cover (non-zero = ice). No landslide is placed on ice. | Intersected with mask_fname if both are given: a pixel can host a landslide only if it is inside the active mask AND not under ice. The number of pixels excluded is reported, since it changes the active area and therefore the landslide count for a prescribed erosion rate. |

### Size distribution

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `landslide_size_model` | `string` | `inverse_gamma` | `param` | - | Size distribution family: inverse_gamma or truncated_pareto. | inverse_gamma (Malamud et al. 2004) has a power-law tail AND a small-area rollover. truncated_pareto (Yanites et al. 2009) is a pure power law, kept to reproduce published results. |
| `inverse_gamma_rho` | `float` | `1.40` | `param` | - | Inverse gamma rho. Tail exponent is -(rho+1), so 1.40 gives -2.40. | Malamud et al. (2004) maximum-likelihood fit. |
| `inverse_gamma_a_m2` | `float` | `1280` | `param` | - | Inverse gamma a (m2). Sets the modal landslide size. | Malamud et al. (2004): 1.28e-3 km2. |
| `inverse_gamma_s_m2` | `float` | `-132` | `param` | - | Inverse gamma s (m2), NEGATIVE. Sets the small-area rollover. | Malamud et al. (2004): -1.32e-4 km2. |
| `pareto_alpha_pdf` | `float` | `2.3` | `param` | - | Truncated Pareto exponent. This is the PDF exponent, NOT cumulative. | alpha_pdf = b_cumulative + 1. Valagussa et al. (2019) measure 1.75-3.11 across lithologies; Niemi's b = 1.18 is alpha_pdf = 2.18. |
| `min_landslide_pixels` | `int` | `1` | `param` | - | Smallest landslide, in pixels. | At coarse resolution this sets a floor above the modal size and the rollover becomes unrepresentable; the tool warns when that happens. |
| `max_landslide_area_m2` | `float` | `0` | `param` | - | Largest landslide (m2). 0 = use the channel threshold, which is the physical cap. | May not exceed threshold_contributing_pixels x cell area: a landslide cannot be bigger than its hillslope. |
| `area_volume_preset` | `string` | `larsen_global` | `param` | - | V = alpha*A^gamma preset (Larsen et al. 2010 Table S1). | larsen_global, larsen_soil_scar, larsen_bedrock_scar, larsen_soil_all, larsen_bedrock_all, hovius_slab, hovius_parabolic. PREFER THE SCAR PRESETS: this tool excavates scars, and the commonly quoted soil/bedrock gammas are combined scar+deposit fits. hovius_* reproduce Niemi/Yanites but overestimate volume by &gt;= 2x per Larsen. |

### Triggering

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `landslide_mode` | `string` | `single_event` | `param` | - | single_event or series. | single_event triggers once at event_time_yr. series runs repeatedly between series_start_yr and series_end_yr. |
| `event_time_yr` | `float` | `0` | `param` | - | Trigger time for single_event mode (yr, on the user's own clock). | Stamped into the time-of-failure raster, which is both the annealing clock and the cosmogenic time-since-change field. |
| `event_volume_m3` | `float` | `0` | `param` | - | Total volume to remove in one event (m3). | - |
| `event_volume_per_area_metres` | `float` | `0` | `param` | - | Event volume per unit ACTIVE area (m3/m2, i.e. metres of lowering). | Used when event_volume_m3 is 0. Active area is the ground that can host a landslide, not the whole raster. |
| `series_start_yr` | `float` | `0` | `param` | - | Series start time (yr). | - |
| `series_end_yr` | `float` | `100000` | `param` | - | Series end time (yr). | - |
| `mean_recurrence_yr` | `float` | `1000` | `param` | - | Mean interval between events in series mode (yr). | - |
| `poisson_timing` | `bool` | `true` | `switch` | - | Draw inter-event intervals from an exponential distribution. | Yanites et al. (2009) Eq. 14. They are explicit that this is a simplification: real landslides are triggered by storms and earthquakes and are neither independent nor uncorrelated. Set false for exactly regular events. |
| `long_term_erosion_rate_mm_kyr` | `float` | `0` | `param` | - | Landslide erosion rate over the active area (mm/kyr), series mode. | The Niemi et al. (2005) convention: prescribe the rate and let the model find the events. Overrides event_volume_m3 in series mode. |
| `anneal_time_yr` | `float` | `0` | `param` | - | A failed pixel cannot fail again for this long (yr). | Stands in for regolith regeneration. UNCALIBRATED -- there is no strong literature constraint on it. |
| `budget_tolerance` | `float` | `0.02` | `param` | - | Relative tolerance when solving for the number of landslides. | Reachable only when the event volume is large compared with the biggest permitted landslide; see the CLI reference. |
| `budget_max_iterations` | `int` | `50` | `param` | - | Iteration cap for the landslide-count solve. | - |
| `random_seed` | `int` | `1` | `param` | - | Seed. The same seed reproduces a run exactly. | - |

### Cosmogenics

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_sample_point_ensemble` | `bool` | `false` | `switch` | - | Run an ensemble at user sample points and write sample-effect CSV. | Computes apparent erosion rates at snapped sample outlets with and without landslides over n_realisations. Reports both the ratio to the no-landslide control and the ratio of the apparent rate to the basin's TRUE mean erosion rate (arithmetic area mean), the latter being the error a real detrital study would incur. Filename: &lt;write_fname&gt;_LS_SAMPLE_EFFECT.csv |
| `sample_points_fname` | `string` | `NULL` | `io` | - | CSV of sample points used by run_sample_point_ensemble. | Must contain latitude/longitude columns (defaults: latitude, longitude). |
| `sample_name_column` | `string` | `sample_name` | `param` | - | Sample-name column in sample_points_fname. | If missing or empty, names default to sample_1, sample_2, ... |
| `sample_latitude_column` | `string` | `latitude` | `param` | - | Latitude column name in sample_points_fname (EPSG:4326 degrees). | - |
| `sample_longitude_column` | `string` | `longitude` | `param` | - | Longitude column name in sample_points_fname (EPSG:4326 degrees). | - |
| `n_realisations` | `int` | `50` | `param` | - | Number of landslide realisations for run_sample_point_ensemble. | - |
| `ensemble_seed_stride` | `int` | `7919` | `param` | - | Seed increment between ensemble realisations. | Realisation i uses random_seed + i*ensemble_seed_stride. |
| `sample_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum snap radius (pixels) for sample points to threshold channel. | Used by expanding-square snapping in run_sample_point_ensemble. |
| `calculate_crn_concentration` | `bool` | `false` | `switch` | - | Compute the post-landslide surface nuclide concentration raster. | Composes existing validated CRN calls: steady state for the inherited profile, transient re-accumulation since failure. Filename: &lt;write_fname&gt;_LS_CONC.tif |
| `nuclide` | `string` | `Be10` | `param` | - | Nuclide: Be10, Al26 or C14. | - |
| `muon_scheme` | `string` | `BraucherBorchers` | `param` | - | Muon production scheme. | - |
| `scaling_scheme` | `string` | `uniform` | `param` | - | Production scaling for the CRN concentration path: uniform, stone or lsdn. | uniform preserves the legacy scalar path using production_scaling and total_shielding (no crn_parameters.nc needed). stone/lsdn build a real CosmoSite from the DEM and use per-pixel production and shielding. |
| `atmosphere` | `string` | `ERA40` | `param` | - | Atmosphere model for site-based production: ERA40, NCEP or standard. | Used only when scaling_scheme is stone or lsdn. |
| `crn_nc_fname` | `string` | - | `param` | - | Path to crn_parameters.nc for site-based production and atmosphere. | Used only when scaling_scheme is stone or lsdn. Empty = bundled nc beside the executable if found, otherwise standard atmosphere fallback. |
| `sample_age_yr` | `float` | `0` | `param` | - | Sample exposure age (yr BP) for the time-integrated production correction. | Used only when scaling_scheme is stone or lsdn and crn_parameters.nc loads. |
| `background_eff_erate_g_cm2_yr` | `float` | `0.001` | `param` | - | Background effective erosion rate (g/cm2/yr), before and after failure. | ~3.8 mm/kyr at rho = 2650. |
| `production_scaling` | `float` | `1.0` | `param` | - | Stone 2000 production scaling (1.0 at SLHL). | - |
| `total_shielding` | `float` | `1.0` | `param` | - | Topographic x snow x self shielding. | Applied through the Vermeesch redistribution, separately from production_scaling. Used directly when scaling_scheme is uniform. With site-based scaling, keep this at 1 and use topo/snow/self parameters. |
| `topo_shield_value` | `float` | `1` | `param` | - | Uniform topographic shielding factor for site-based scaling (0-1). | Used when scaling_scheme is stone or lsdn. |
| `snow_shield_eff_thickness` | `float` | `0` | `param` | - | Uniform snow shielding as effective depth (g/cm2) for site-based scaling. | Used when scaling_scheme is stone or lsdn. |
| `self_shield_eff_thickness` | `float` | `0` | `param` | - | Uniform self/sample shielding as effective depth (g/cm2) for site-based scaling. | Used when scaling_scheme is stone or lsdn. |
| `rock_density_kg_m3` | `float` | `2650` | `param` | - | Density of the ERODING BEDROCK COLUMN (kg/m3). | Converts the background erosion rate between m/yr and g/cm2/yr. This is NOT the density of the material a landslide removes -- see scar_density_kg_m3. |
| `scar_density_kg_m3` | `float` | `-1` | `param` | - | Density of the material a landslide EXCAVATES (kg/m3). Negative takes it from area_volume_preset. | This is a different quantity from rock_density_kg_m3 and conflating them is a real error. Larsen's V = alpha A^gamma is fitted to MEASURED scar volumes, so the excavation depth is a true geometric thickness of whatever failed -- soil, at roughly 1500 kg/m3, for a shallow landslide. Carrying that at 2650 overstates the mass depth by about a factor of 1.77, and since production falls exponentially with mass depth it changes the concentration of the delivered material by more than a factor of two. The material-specific presets carry a density: larsen_soil_scar 1500, larsen_bedrock_scar 2650. larsen_global and the hovius presets are fitted across materials and carry none, so leaving this negative with one of those selected raises a warning. |
| `crn_report_time_yr` | `float` | `0` | `param` | - | Time at which concentrations are reported, on the same clock as the event times. | - |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

<!-- END GENERATED CLI PARAMETERS -->
