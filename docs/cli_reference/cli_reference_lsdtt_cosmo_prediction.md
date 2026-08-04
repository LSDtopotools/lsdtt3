# CLI reference: lsdtt-cosmo-prediction

Predicts cosmogenic radionuclide (CRN) concentrations across a DEM under three scenarios: steady-state erosion, transient step-change erosion, and erosion-weighted downstream accumulation. Uses Stone (2000) spallation scaling with optional ERA40 atmospheric and palaeomagnetic corrections via a NetCDF file.

**Supported nuclides:** ¹⁰Be (`Be10`), ²⁶Al (`Al26`), ¹⁴C (`C14`).

**Units:** concentrations in atoms g⁻¹; effective erosion rates in g cm⁻² yr⁻¹.  
To convert: ε [g cm⁻² yr⁻¹] = ε [m yr⁻¹] × ρ [kg m⁻³] × 0.1.

## Usage

```
lsdtt-cosmo-prediction  <path> <parameter_file>
lsdtt-cosmo-prediction  <parameter_file>        # path defaults to ./
lsdtt-cosmo-prediction  -h                      # print help and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

At least one prediction mode must be activated; the program exits with a message if none are set.

---

## Pipeline overview

| Step | Description |
|---|---|
| 1. Load DEM | GDAL-supported format; must be in projected metres |
| 2. Generate the step change | *(optional)* `make_step_change_rasters` — derives erosion rate, time-since-change and topography from the analytical fluvial solution |
| 3. Production raster | Compute Stone 2000 / LSDn, or load from disk |
| 4. Shielding rasters | Topographic, snow, self (raster or uniform scalar) |
| 5. Erosion raster | Raster, uniform scalar, or the generated field |
| 6. Quartz fertility and ice | *(optional)* fertility field; the ice mask rewrites snow shielding, erosion rate and fertility |
| 7. Steady-state concentrations | `calculate_crn_concentration_raster` |
| 8. Step-change concentrations | `calculate_step_change_crn_concentration_raster` or the multi-isotope modes |
| 9. Landslide coupling | *(optional)* `landslide_coupling` — scars and/or flushed debris superimposed on step 8 |
| 10. Accumulated concentrations | `calculate_accumulated_crn_concentration`, mixing on erosion rate or on quartz |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | `./` | Directory containing the input DEM. |
| `read_fname` | string | _(empty)_ | Input DEM filename (including extension). |
| `write_path` | string | `./` | Directory for output files. |
| `write_fname` | string | _(same as read_fname stem)_ | Output filename stem (no extension). |
| `write_output_format` | string | `cog` | Output raster format. `cog`, `gtiff`, or `same_as_input`. |

---

## Nuclide and muon scheme

| Parameter | Type | Default | Description |
|---|---|---|---|
| `nuclide` | string | `Be10` | Nuclide to predict. `Be10`, `Al26`, or `C14`. |
| `muon_scheme` | string | `BraucherBorchers` | Muon production scheme. `BraucherBorchers`, `Braucher`, `Granger`, `Schaller`, `newCRONUS`, `cronusEmulatorSt`, `cronusEmulatorLSDn`. The `cronusEmulator*` schemes reproduce the CRONUS v3 online calculator; they are pressure-dependent (a pressure raster is computed automatically) and are **steady-state only** — the step-change and transient prediction modes reject them. |

---

## Prediction modes

At least one must be set to `true`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_crn_concentration_raster` | bool | `false` | Compute steady-state concentration at every pixel. Requires erosion rate field and shielding rasters. |
| `calculate_step_change_crn_concentration_raster` | bool | `false` | Compute transient concentration after an instantaneous step change in erosion rate. Requires `step_change_eff_erate_old/new` and `time_since_change`. |
| `calculate_multi_isotope_step_change` | bool | `false` | Carry Be-10, Al-26 and C-14 through the **same** step change at every pixel, and write the two isotope-ratio rasters. |
| `calculate_multi_isotope_transient` | bool | `false` | The same three nuclides through an **arbitrary erosion history** from `transient_history_file`. The general case of the step-change mode. |
| `calculate_accumulated_crn_concentration` | bool | `false` | Accumulate concentration downstream using D8 flow routing. Requires a concentration raster (computed in the same run or loaded separately). |

Three further capabilities are *modifiers* rather than modes, each opt-in and a
no-op when left alone:

| Section | Switch | What it does |
|---|---|---|
| [Generating the step change](#generating-the-step-change-from-a-fluvial-solution) | `make_step_change_rasters` | builds the time-since-change field instead of requiring one |
| [Landslide coupling](#landslide-coupling) | `landslide_coupling` | superimposes scars and/or flushed debris |
| [Quartz fertility and ice](#quartz-fertility-and-ice) | `quartz_fertility*`, `ice_mask_raster_prefix` | changes what the downstream mixture is weighted by |

### Multi-isotope parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `transient_history_file` | string | `NULL` | CSV erosion history for `calculate_multi_isotope_transient`: `time,rate,removal` (yr BP, m/yr, m). Rows run most-recent-first; the last row must be `time = -9999` carrying the background rate. Unordered records are rejected fatally. |
| `transient_rock_density_kg_m3` | float | `2650` | Density converting that history's metre-based rates and removals into effective depths. |
| `write_isotope_ratio_rasters` | bool | `true` | Write `_R2610_step.tif` (Al26/Be10) and `_R1410_step.tif` (C14/Be10). Needs `calculate_multi_isotope_step_change`. |

### Scaling scheme

| Parameter | Type | Default | Description |
|---|---|---|---|
| `scaling_scheme` | string | `stone` | `stone` = Stone (2000)/Lal, the default and the CAIRN/OCTOPUS scaling. `lsdn` = Lifton, Sato & Dunai (2014), interpolated from the LSDn2015 grids; needs those grids in `crn_parameters.nc` and supports `Be10`, `Al26`, `C14` only. Both are evaluated **per pixel**. `muon_scheme` is unaffected: LSDn is a spallation scaling and the muon pathways receive it through the same Vermeesch redistribution as Stone. |

---

## Production raster

By default the Stone 2000 spallation production rate is computed from DEM elevation. A pre-computed production raster can be loaded and optionally saved to disk.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `load_production_raster` | bool | `false` | Load production raster from disk instead of computing. |
| `production_raster_prefix` | string | `NULL` | Filename prefix (no extension) for a pre-computed production raster. Path: `read_path + prefix + .tif`. Used when `load_production_raster: true`. |
| `write_production_raster` | bool | `false` | Write the computed production raster to disk. Output: `<write_fname>_crn_prod.tif`. |
| `crn_nc_file` | string | _(empty)_ | Path to `crn_parameters.nc` for ERA40 atmospheric scaling and palaeomagnetic correction. Leave empty to use the `crn_parameters.nc` bundled with lsdtt3 (found automatically); if absent, standard-atmosphere Stone 2000. |
| `sample_age_yr` | float | `0` | Sample exposure age (yr) for palaeomagnetic correction. Only used when `crn_nc_file` is set and value > 0. |
| `atmosphere` | string | `ERA40` | `ERA40` (Balco/Lifton climatology + latitude-dependent lapse, the default), `NCEP` (the atmosphere CAIRN/OCTOPUS use; needs the NCEP grids), or `standard` (SLHL barometric formula). `ERA40`/`NCEP` fall back to `standard` if their grids are absent; an unknown name is fatal. **Match this to `lsdtt-cosmo-tool` when using the two together** — the atmosphere sets the pressure at every pixel, so forward-modelling under one and inverting under another silently biases the recovered erosion rate. |

---

## Shielding rasters

Each shielding type accepts either a raster (load by prefix) or a uniform scalar. If the prefix is `NULL` the scalar is used.

### Topographic shielding (dimensionless, 0–1; 1 = fully exposed)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `topographic_shielding_raster_prefix` | string | `NULL` | Prefix of a topographic shielding raster. Path: `read_path + prefix + .tif`. |
| `topo_shield_value` | float | `1.0` | Uniform topographic shielding value (0–1) used when prefix is `NULL`. |

### Snow shielding (effective depth, g cm⁻²)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `snow_shielding_raster_prefix` | string | `NULL` | Prefix of a snow shielding raster (effective depth, g/cm²). |
| `snow_shield_eff_thickness` | float | `0.0` | Uniform snow shielding effective depth (g/cm²). `0` = no snow. |

### Self / sample shielding (effective depth, g cm⁻²)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `self_shielding_raster_prefix` | string | `NULL` | Prefix of a self/sample shielding raster (effective depth, g/cm²). |
| `self_shield_eff_thickness` | float | `0.0` | Uniform self shielding effective depth (g/cm²). `0` = surface sample. |

---

## Erosion rate

Used for steady-state and accumulation modes.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `background_erosion_rate_raster_prefix` | string | `NULL` | Prefix of a raster giving effective erosion rate (g/cm²/yr) at each pixel. |
| `effective_erosion_rate` | float | `0.001` | Uniform effective erosion rate (g/cm²/yr). Used when prefix is `NULL`. ≈ 0.001 g/cm²/yr ≈ 3.8 mm/kyr for ρ = 2650 kg/m³. |

---

## Step-change parameters

Used only when `calculate_step_change_crn_concentration_raster: true`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `step_change_eff_erate_old_g_cm2_yr` | float | `0.001` | Pre-change effective erosion rate (g/cm²/yr). |
| `step_change_eff_erate_new_g_cm2_yr` | float | `0.01` | Post-change effective erosion rate (g/cm²/yr). |
| `time_since_change_raster_prefix` | string | `NULL` | Prefix of a raster giving time since step change (yr) at each pixel. Pixels with value 0 receive steady-state concentration at the new erosion rate. |
| `time_since_change_yr` | float | `10000` | Uniform time since step change (yr) used when prefix is `NULL`. |

---

## Generating the step change from a fluvial solution

The time-since-change field above is usually the hard input to supply, because
it is not a free parameter: it is **the arrival time of the incision wave**. Set
`make_step_change_rasters: true` and the tool derives it, along with the erosion
rate and a matching topography, from the analytical transient response of
detachment-limited stream power to a step in uplift (Royden & Perron 2013;
Mitchell & Yanites 2019). The generated fields are handed to the step-change
concentration kernel **in memory** — there is no round trip through disk.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `make_step_change_rasters` | bool | `false` | Generate the erosion-rate, time-since-change and topography fields instead of loading them. |
| `step_change_uplift_old_m_yr` | float | `0.0001` | Uplift (= steady erosion) rate before the step, **m/yr**. |
| `step_change_uplift_new_m_yr` | float | `0.001` | Uplift rate after the step, **m/yr**. |
| `step_change_target_relief_m` | float | `1000` | Relief the old steady-state landscape should have. `K` is back-calculated from it. |
| `step_change_age_yr` | float | `100000` | Time since the step reached the basin outlet. Every pixel gets its own, smaller, elapsed time. |
| `step_change_n` | float | `1` | Slope exponent `n` in `E = K A^m S^n`. |
| `step_change_m_over_n` | float | `0.45` | Concavity for the chi integral. |
| `step_change_A_0` | float | `1` | Reference drainage area for chi. Any positive value is handled exactly. |
| `step_change_channel_pixels` | int | `500` | Contributing-pixel count above which a pixel is a channel; below it a critical hillslope gradient is imposed. |
| `step_change_Sc_old` | float | `0.25` | Critical hillslope gradient in the unadjusted part of the landscape. |
| `step_change_Sc_new` | float | `0.4` | Critical hillslope gradient in the adjusted part. |
| `step_change_basin_outlet_csv` | string | `NULL` | CSV of `latitude,longitude` giving the basin outlet. The first point is used and snapped to the channel network. `NULL` uses the largest basin in the DEM. |
| `rock_density_kg_m3` | float | `2650` | Converts the uplift rates (m/yr) to effective erosion rates (g/cm²/yr). |
| `step_change_production_from_generated_dem` | bool | `true` | Compute production from the generated topography rather than the input DEM. |
| `write_step_change_rasters` | bool | `false` | Write `_SC_DEM`, `_SC_ERATE`, `_SC_STEPTIME`, `_SC_TAG` and `_SC_CHI`. |

**Read the derived numbers.** `K` is back-calculated from
`step_change_target_relief_m`, so the entire configuration — knickpoint
position, wave celerity, the whole time-since-change field — follows from a
relief you chose rather than a rate constant you typed. The run prints `K`, both
steepness indices, the knickpoint position and elevation, and the adjusted area.
Report them.

**The uplift rates override the effective-rate parameters.**
`step_change_eff_erate_old_g_cm2_yr` and `..._new_...` are ignored in this mode;
the rates come from `step_change_uplift_*_m_yr × rock_density_kg_m3 × 0.1`. The
run says so explicitly. Likewise `time_since_change_raster_prefix` is ignored,
with a note.

**What the wave does depends on `n`.** With `n ≥ 1` and rising uplift the
transition steepens into a migrating knickpoint (a shock) separating an adjusted
reach below from an unadjusted reach above. With `n < 1` it stretches into a
rarefaction fan, inside which the erosion rate varies continuously between the
old and new values — those pixels are tagged separately. A *falling* step swaps
which case applies, and that is handled: the choice is made on the wave
celerities, not on `n` alone.

**Use the tag raster, never chi, to ask whether a pixel has been rejuvenated.**
`_SC_TAG` is 1 = unadjusted channel, 2 = unadjusted hillslope, 3 = adjusted
channel, 4 = adjusted hillslope, 5 and 6 = the same two inside a fan. Chi
increases upslope, so a hillslope always has a larger chi than the channel it
drains to; a hillslope draining to an adjusted channel therefore carries the new
erosion rate while sitting *above* the knickpoint in chi. That is correct — a
hillslope belongs to its channel — but it means a chi comparison misclassifies
hillslopes.

**Scope and side effects.** The generated fields cover **one basin**, so the
concentration output is restricted to it. The generated topography is offset so
its outlet sits at the input DEM's outlet elevation, keeping base level real.
The generated per-pixel erosion rate is used for steady-state concentrations and
as the flux weight in downstream accumulation. Flow routing for accumulation
still uses the **input** DEM, not the generated one.

### Example: generate and consume in one run

```
read_path: ./
read_fname: my_dem.tif
write_path: ./
write_fname: sc

make_step_change_rasters: true
step_change_uplift_old_m_yr: 0.0001
step_change_uplift_new_m_yr: 0.0005
step_change_target_relief_m: 1200
step_change_age_yr: 500000
step_change_n: 1.0
step_change_channel_pixels: 500
write_step_change_rasters: true

calculate_multi_isotope_step_change: true
```

---

## Landslide coupling

Superimposes landslide scars on a step-change concentration field, using the
time-of-failure and excavation-depth rasters written by `lsdtt-landslides`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `landslide_coupling` | string | `none` | `none`, `scar`, `flush` or `both`. |
| `landslide_time_raster_prefix` | string | `NULL` | Prefix of `<stem>_LS_TFAIL` from `lsdtt-landslides`. |
| `landslide_excavation_raster_prefix` | string | `NULL` | Prefix of `<stem>_LS_EXCAVATION`. |
| `landslide_report_time_yr` | float | `0` | Present time on the same clock as the failure-time raster. |
| `flush_time_threshold_yr` | float | `100` | How long a scar keeps delivering debris (used by `both`). |
| `scar_density_kg_m3` | float | `1500` | Density of the material the landslide **excavated**. |
| `write_delivery_mass_raster` | bool | `false` | Write `<stem>_LS_DELIVERY.tif`, the downstream mixing weight. |

**The `scar` scenario.** The landslide material has long since left the system.
What remains is a scar eroding at the catchment rate from a column whose
inheritance was stripped when it failed — so the surface there is depleted
relative to unfailed ground, and recovers as it re-accumulates. Each failed
pixel gets a four-row erosion history: the old background rate, the change, an
instantaneous removal of the excavated depth at the failure time, then the new
rate to the present.

**The `flush` scenario.** The debris is still in the channel, so the sample
contains the **whole excavated slab** rather than the surface of the scar: the
depth average of the pre-failure profile over `[0, d_exc]`. `both` applies the
flush treatment inside `flush_time_threshold_yr` and the scar treatment outside
it, which is the physically complete case — one kernel with one number, so
`scar` is the zero-window limit and `flush` the unbounded one.

**The two do not dilute by the same amount, and the ordering is the opposite of
the intuitive one.** At short ages

```
background  >  flush  >  scar
```

because the flushed slab *contains* the pre-failure column's concentrated top,
while the scar is what is left behind — the material from the very base of the
slab, exposed and barely re-accumulated. In the limit of a fresh failure the
inequality `flush > scar` is exact: the mean of a decreasing profile over
`[0, d]` must exceed its value at `d`. As a scar re-accumulates it eventually
overtakes the flush value.

**The flush timescale is an age, following Niemi.** Niemi et al. (2005) have no
sediment storage at all — landslide material passes through within one 100 yr
timestep and mixes homogeneously — which is where the default comes from. Set it
to the event recurrence interval to reproduce that convention. How long
landslide sediment actually persists in a channel network is a separate research
question and is not modelled here.

### The delivery-mass raster

`write_delivery_mass_raster` writes the weight a downstream mixing calculation
should use, in **kg m⁻² yr⁻¹**:

```
flushing pixel : rho_scar * d_exc / flush_window
everything else: the background erosion rate, converted from g/cm2/yr
```

Both terms are the same quantity in the same units, which is the point of not
reusing the erosion-rate raster as a weight. **The erosion-rate raster is left
untouched**, so the *true* erosion rate stays available to compare an apparent
rate against — the legacy implementation overwrote it in place and destroyed
exactly that.

**Unfailed ground is bitwise unchanged.** The ordinary step-change kernel
computes the whole raster and only failed pixels are overwritten, so switching
the coupling on cannot perturb ground that never failed. That is asserted by a
regression test, because a coupling that shifted the whole raster would move
every basin average for reasons nothing to do with landslides.

**Failures older than the incision wave are discarded.** A pixel whose landslide
pre-dates the arrival of the erosion-rate change at *that pixel* cannot have
been caused by it; those keep their plain step-change value and are counted
separately in the run report.

**"Never failed" is not nodata.** The failure-time raster holds an absolute
time with nodata on ground that never failed; it is converted here to elapsed
years with "never" carried as a distinct state. Encoding it as nodata would drop
unfailed ground out of the raster, and hence out of any basin average taken from
it — replacing a small correction with a large error.

### Two densities, and why they are not the same number

| Parameter | Meaning | Typical |
|---|---|---|
| `rock_density_kg_m3` | the eroding **bedrock column** | 2650 |
| `scar_density_kg_m3` | the material a landslide **excavated** | ~1500 soil, 2650 rock |

The excavation depth from `lsdtt-landslides` is a **true geometric thickness**,
because Larsen's `V = αA^γ` is fitted to *measured* scar volumes. For a shallow
landslide the material that failed is soil, at roughly 1500 kg/m³. Carrying it at
the bedrock density overstates the mass depth by a factor of about 1.77 — and
since production falls exponentially with mass depth, that is more than a factor
of two in the concentration of the material delivered.

`lsdtt-landslides` takes its scar density from the `area_volume_preset` in play
(`larsen_soil_scar` → 1500, `larsen_bedrock_scar` → 2650) and **warns** if the
preset selected is one fitted across both materials, such as the default
`larsen_global`, which has no defensible density. Set the two consistently
between the two tools.

### Example

```
calculate_step_change_crn_concentration_raster: true
step_change_eff_erate_old_g_cm2_yr: 0.001
step_change_eff_erate_new_g_cm2_yr: 0.005
time_since_change_yr: 20000

landslide_coupling: both
landslide_time_raster_prefix: ls_LS_TFAIL
landslide_excavation_raster_prefix: ls_LS_EXCAVATION
landslide_report_time_yr: 5000
flush_time_threshold_yr: 100
scar_density_kg_m3: 1500
write_delivery_mass_raster: true
```

Not yet available with `calculate_multi_isotope_step_change`; use a single
nuclide.

---

## Quartz fertility and ice

A basin-averaged sample is a mixture of **quartz**, not of rock. Two pixels
eroding at the same rate contribute unequally if their rock carries different
quartz, and a quartz-free unit — limestone, basalt, many volcanics and mafics —
contributes sediment but no sample at all. In mixed-lithology terrain, routing
concentrations without accounting for that is **wrong**, not merely imprecise.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `quartz_fertility` | float | `1` | Uniform quartz mass fraction (0–1). `1` is a no-op. |
| `quartz_fertility_raster_prefix` | string | `NULL` | Per-pixel fertility, typically a reclassified lithology map. |
| `quartz_fertility_regolith` | float | `-1` | Two-layer model: regolith fraction. Negative = off. |
| `quartz_fertility_bedrock` | float | `-1` | Two-layer model: fresh-bedrock fraction. |
| `quartz_fertility_regolith_depth_m` | float | `0` | Regolith thickness separating the layers. |
| `ice_mask_raster_prefix` | string | `NULL` | Present-day ice cover (non-zero = ice). |
| `ice_thickness_m` | float | `50` | Uniform ice thickness. |
| `ice_density_kg_m3` | float | `917` | Density of ice; lower for firn. |
| `erate_under_ice_g_cm2_yr` | float | `0` | Erosion rate beneath the ice. |
| `quartz_fertility_under_ice` | float | `0` | Quartz fraction beneath the ice. |
| `write_quartz_weight_raster` | bool | `false` | Write `<stem>_QUARTZ_WEIGHT.tif`. |

**Everything here is opt-in and a no-op by default.** With fertility 1 and no
ice mask, the accumulation keeps mixing on the erosion-rate raster exactly as
before. Note also that a *uniform* fertility changes nothing either: a constant
weight cancels in `Σ(Cw)/Σ(w)`. Fertility only bites where it **varies**.

**A quartz-free unit is excluded however fast it erodes.** That is the case
that makes this a correctness issue: a basalt unit eroding ten times faster than
its granite neighbour dominates the sediment flux and contributes nothing to the
sample.

**An all-quartz-free domain is reported, not returned as nodata.** "This
catchment yields no quartz" is a meaningful answer, and a common one in volcanic
terrain. The run says so and stops rather than emitting an empty raster.

**A mis-scaled fertility raster is caught.** Values outside [0,1] — a lithology
map supplied as percentages, say — are clamped and the number clamped is
reported.

### The two-layer model

Chemical weathering removes soluble minerals and leaves quartz behind, so
regolith is enriched relative to the fresh bedrock a deep landslide scar
exposes. Setting the three `quartz_fertility_regolith*` parameters makes a
landslide deliver the **depth average** of the two layers over its excavated
slab:

```
f(d) = [ f_reg * min(d, D) + f_bed * max(0, d - D) ] / d
```

If that picture is right, landslide debris is *under-represented* in the quartz
separate relative to its mass, and the dilution a landslide causes is **weaker**
than a mass-weighted calculation predicts. This is the switch that tests it.

### Ice

Ice changes exactly three things, and lsdtt3 already has a field for each: the
ice column shields cosmic rays (**snow shielding**, as an effective depth), the
**erosion rate** differs or is zero, and little or no **quartz** reaches the
sample. So the ice mask is a helper that writes three existing rasters, not a
model of its own. Ten metres of ice is already 917 g cm⁻², at which 0.3 % of
spallation survives.

`erate_under_ice_g_cm2_yr: 0` is the cold-based idealisation and is numerically
safe: the steady-state kernel forms its denominator as `(ε + Γλ)` and never
divides by ε, so ε = 0 returns the saturation concentration `P/λ`.

`lsdtt-landslides` accepts the same mask as `ice_mask_fname` and places no
landslide on ice; the pixels excluded are reported, since they change the active
area and hence the landslide count for a prescribed erosion rate.

> **Limitation, stated rather than hidden.** This models **present-day ice
> only**. Concentrations depend on the shielding *history*, and real glaciated
> landscapes advanced and retreated repeatedly. lsdtt3 can express a per-pixel
> erosion history but not a per-pixel time-varying shielding, so a
> glacial–interglacial cycle cannot be represented at all, and concentrations
> under formerly glaciated ground will be overestimated.

---

## Accumulated concentration

Used only when `calculate_accumulated_crn_concentration: true`.

The concentration raster to accumulate is selected in priority order:
1. `concentration_raster_prefix` (if set)
2. Step-change concentration computed in the same run
3. Steady-state concentration computed in the same run

| Parameter | Type | Default | Description |
|---|---|---|---|
| `concentration_raster_prefix` | string | `NULL` | Prefix of a pre-computed concentration raster to use for accumulation. |
| `fill_before_routing` | bool | `true` | Fill DEM pits before building D8 flow routing for accumulation. |
| `min_slope_for_fill` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. |
| `boundary_conditions` | string | `oooo` | D8 boundary conditions (N, E, S, W). `o` = open, `b` = forced outlet. |

---

## Sampling the accumulated raster at points

Closes a forward → inverse round trip: the CSV written here is read directly by
`lsdtt-cosmo-tool`, so an imposed erosion rate can be pushed through the forward
model and recovered by the inversion.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sample_points_filename` | string | `NULL` | CSV of `latitude,longitude` at which to sample the accumulated raster. Needs `calculate_accumulated_crn_concentration: true`. |
| `sample_column_name` | string | `sample_name` | Optional sample-identifier column copied through from the points CSV. Omitted from the output if absent. |
| `threshold_contributing_pixels` | int | `1000` | Minimum contributing pixels for a node to count as a channel when snapping. Matches `lsdtt-cosmo-tool`. |
| `snap_search_radius_pixels` | int | `8` | Maximum search radius when snapping a point to the channel network. Matches `lsdtt-cosmo-tool`. |

Output is `<write_fname>_sampled_crn_conc.csv` with columns `latitude`,
`longitude`, `<nuclide>_CONC` and `snapped_contributing_pixels`. Snapping mirrors
`lsdtt-cosmo-tool` so both tools resolve the same node.

`threshold_contributing_pixels` is also used to snap
`step_change_basin_outlet_csv` when generating a step change.

---

## Output write toggles

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_crn_concentration_raster` | bool | `true` | Write the steady-state concentration raster. Output: `<write_fname>_<nuclide>_crn_conc.tif`. |
| `write_step_change_concentration_raster` | bool | `true` | Write the step-change concentration raster. Output: `<write_fname>_<nuclide>_step_crn_conc.tif`. |
| `write_accumulated_crn_raster` | bool | `true` | Write the accumulated concentration raster. Output: `<write_fname>_<nuclide>_crn_acc.tif`. |
| `write_production_raster` | bool | `false` | Write the Stone 2000 production raster. Output: `<write_fname>_crn_prod.tif`. |

---

## lsdviztools3 interface notes

`lsdtt-cosmo-prediction` is the forward side used by lsdtt3 comparison
workflows (including point sweeps) that are then plotted by `lsdviztools3`.

- For point/sweep workflows, sampled concentrations are emitted as
  `<write_fname>_sampled_crn_conc.csv` and paired with inverse runs from
  `lsdtt-cosmo-tool`.
- The Phase-M2 sweep script
  `tests/scripts/run_cosmo_single_point_scaling_scan.py` writes
  `single_point_scan_comparison.csv`.
- `lsdviztools3` command for sweep comparison figures:

```
lsdviz3-cosmo-compare \
  --scan-csv tests/expected/cosmo/generated/single_point_scan_comparison.csv \
  --outdir demo_output/cosmo
```

Optional `--sample-csv` can be used to include a sample-level residual figure in
the same run.

---

## Example parameter files

### Steady-state concentration map

```
read_path:   /data/dems/
read_fname:  catchment.tif

nuclide:     Be10
muon_scheme: BraucherBorchers

calculate_crn_concentration_raster: true

effective_erosion_rate: 0.003   # ~ 11 mm/kyr for rock density 2650 kg/m3
topo_shield_value:      1.0
snow_shield_eff_thickness: 0.0
self_shield_eff_thickness: 0.0

write_crn_concentration_raster: true
```

### Step-change transient concentration

```
read_path:   /data/dems/
read_fname:  catchment.tif

nuclide:     Be10

calculate_step_change_crn_concentration_raster: true

step_change_eff_erate_old_g_cm2_yr: 0.001
step_change_eff_erate_new_g_cm2_yr: 0.01
time_since_change_yr:               5000
```

### Downstream accumulated concentration

```
read_path:   /data/dems/
read_fname:  catchment.tif

nuclide:     Be10

# Compute steady-state first, then accumulate
calculate_crn_concentration_raster:      true
calculate_accumulated_crn_concentration: true

effective_erosion_rate:  0.002
fill_before_routing:     true
boundary_conditions:     oooo
```

---

## References

Stone, J.O. (2000). Air pressure and cosmogenic isotope production. *Journal of Geophysical Research*, 105(B10), 23753–23759.

Braucher, R., et al. (2011). Production rates of cosmogenic radionuclides at the atmosphere–land surface interface. *Earth and Planetary Science Letters*, 304(3–4), 417–428.


### Multi-isotope step change

Isotope **ratios** carry the transient signal. The half-lives differ by orders of
magnitude — C-14 5.7 kyr, Al-26 0.7 Myr, Be-10 1.39 Myr — so after a change in
erosion rate the three re-equilibrate at different speeds, and the pair dates the
perturbation in a way one nuclide alone cannot.

A ratio is only meaningful if both members saw an identical scenario, so
`calculate_multi_isotope_step_change` computes all three from **one particle
evolution per pixel** — same shielding, same rates, same elapsed time — rather
than from three separate runs. Each nuclide is bitwise identical to what
`calculate_step_change_crn_concentration_raster` returns for the same inputs, and
the whole thing costs little more than a single-nuclide run.

Inputs are the same as the single-nuclide step change: `time_since_change_yr` (or
a raster), `step_change_eff_erate_old_g_cm2_yr`, `step_change_eff_erate_new_g_cm2_yr`,
the production raster and topographic shielding. The `nuclide` parameter is
ignored — all three are always produced.

Outputs: `_Be10_step_crn_conc.tif`, `_Al26_step_crn_conc.tif`,
`_C14_step_crn_conc.tif`, plus `_R2610_step.tif` and `_R1410_step.tif`. Pixels
with zero Be-10 are written as nodata in the ratio rasters and the count is
reported.

Note that after an erosion-rate *increase* all concentrations fall, C-14 reaches
its new value first, and the C14/Be10 ratio is therefore **suppressed** and
approaches its endpoint from below — it does not overshoot.


### Multi-isotope transient (arbitrary history)

`calculate_multi_isotope_transient` generalises the step-change mode: instead of
one change in erosion rate, every pixel follows the history in
`transient_history_file`. A two-row history reproduces the step-change mode
**bitwise**, so this is a strict generalisation rather than a parallel
implementation — the per-pixel walk is the same code the column tool uses.

History file format — three columns, `time,rate,removal`:

```
time,rate,removal
3000,0.00113,0        <- most recent phase: 0.00113 m/yr for the last 3 kyr
30000,0.0000377,0     <- before that, 0.0000377 m/yr back to 30 kyr
-9999,0.0000755,0     <- long-term background rate (required last row)
```

Rows run from the **most recent** (index 0) to the oldest, and the final row must
carry `time = -9999` with the long-term background rate that initialises the
column. A row with `rate = 0` and a non-zero `removal` is an instantaneous
stripping event in metres, and must repeat the time of the row before it.

**An unordered history is rejected**, not silently accepted: reversed rows would
produce negative timesteps and a confidently wrong answer.

Outputs are the same five rasters as the step-change mode. Cost is
`O(pixels × history rows)`, so this is substantially slower — inherent, since an
arbitrary history has no closed form.
