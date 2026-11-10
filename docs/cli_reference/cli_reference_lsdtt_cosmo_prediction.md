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
| `calculate_multi_isotope_transient` | bool | `false` | The same three nuclides through an **arbitrary erosion history** from `transient_history_fname`. The general case of the step-change mode. |
| `calculate_accumulated_crn_concentration` | bool | `false` | Accumulate concentration downstream using D8 flow routing. In multi-isotope mode (with no `concentration_raster_prefix`), Be10/Al26/C14 are accumulated together in one pass and accumulated ratios can be written. |

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
| `transient_history_fname` | string | `NULL` | CSV erosion history for `calculate_multi_isotope_transient`. `NULL` is deliberate because transient prediction is optional; the file is required when that mode is enabled. Columns are `time,rate,removal` (yr BP, m/yr, m); rows run most-recent-first and end with `time = -9999`. |
| `transient_rock_density_kg_m3` | float | `2650` | Density converting that history's metre-based rates and removals into effective depths. |
| `write_isotope_ratio_rasters` | bool | `true` | In multi-isotope step-change/transient mode, write `_R2610_step.tif` (Al26/Be10) and `_R1410_step.tif` (C14/Be10). If accumulation is also on, also write `_R2610_acc.tif` and `_R1410_acc.tif` from accumulated basin means. |

### Scaling scheme

| Parameter | Type | Default | Description |
|---|---|---|---|
| `scaling_scheme` | string | `stone` | `stone` = Stone (2000)/Lal, the site-based raster-prediction default and the CAIRN/OCTOPUS scaling. The landslide generator instead defaults to `uniform` to preserve its legacy scalar path. `lsdn` = Lifton, Sato & Dunai (2014), interpolated from the LSDn2015 grids. Both site schemes are evaluated **per pixel**. |

---

## Production raster

By default the Stone 2000 spallation production rate is computed from DEM elevation. A pre-computed production raster can be loaded and optionally saved to disk.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `load_production_raster` | bool | `false` | Load production raster from disk instead of computing. |
| `production_raster_prefix` | string | `NULL` | Filename prefix (no extension) for a pre-computed production raster. Path: `read_path + prefix + .tif`. Used when `load_production_raster: true`. |
| `write_production_raster` | bool | `false` | Write the computed production raster to disk. Output: `<write_fname>_crn_prod.tif`. |
| `crn_nc_fname` | string | _(empty)_ | Path to `crn_parameters.nc` for ERA40 atmospheric scaling and palaeomagnetic correction. Leave empty to use the `crn_parameters.nc` bundled with lsdtt3 (found automatically); if absent, standard-atmosphere Stone 2000. |
| `sample_age_yr` | float | `0` | Sample exposure age (yr) for palaeomagnetic correction. Only used when `crn_nc_fname` is set and value > 0. |
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
| `step_change_target_relief_metres` | float | `1000` | Relief the old steady-state landscape should have. `K` is back-calculated from it. |
| `step_change_age_yr` | float | `100000` | Time since the step reached the basin outlet. Every pixel gets its own, smaller, elapsed time. |
| `step_change_n` | float | `1` | Slope exponent `n` in `E = K A^m S^n`. |
| `step_change_m_over_n` | float | `0.45` | Concavity for the chi integral. |
| `step_change_A_0` | float | `1` | Reference drainage area for chi. Any positive value is handled exactly. |
| `step_change_channel_pixels` | int | `500` | Contributing-pixel count above which a pixel is a channel; below it a critical hillslope gradient is imposed. |
| `step_change_Sc_old` | float | `0.25` | Critical hillslope gradient in the unadjusted part of the landscape. |
| `step_change_Sc_new` | float | `0.4` | Critical hillslope gradient in the adjusted part. |
| `step_change_basin_outlet_fname` | string | `NULL` | CSV of `latitude,longitude` giving the basin outlet. The first point is used and snapped to the channel network. `NULL` uses the largest basin in the DEM. |
| `rock_density_kg_m3` | float | `2650` | Converts the uplift rates (m/yr) to effective erosion rates (g/cm²/yr). |
| `step_change_production_from_generated_dem` | bool | `true` | Compute production from the generated topography rather than the input DEM. |
| `step_change_route_on_generated_dem` | bool | `true` | Build D8 flow routing on the generated topography for accumulation. Set `false` to route on the original DEM instead. |
| `write_step_change_rasters` | bool | `false` | Write `_SC_DEM`, `_SC_ERATE`, `_SC_STEPTIME`, `_SC_TAG` and `_SC_CHI`. |

**Read the derived numbers.** `K` is back-calculated from
`step_change_target_relief_metres`, so the entire configuration — knickpoint
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
as the flux weight in downstream accumulation.

**Flow routing: generated or original DEM.** Set `step_change_route_on_generated_dem`
to control which DEM is used for D8 flow accumulation. The default (`true`) rebuilds
routing on the generated topography, which matches the scenario's rates and times.
However, the generator snaps hillslopes to different critical gradients on each side
of the knickpoint, creating an elevation cliff along the drainage divide that can
cause cross-divide flow artefacts. Setting `false` routes on the original DEM,
avoiding the cliff but routing against a topography whose rates no longer match the
model. Production, elevation and the step-change fields still come from the generated
surface either way. Accumulation needs routing; set `step_change_route_on_generated_dem: false`
if cross-divide flow is biasing your basin averages.

### Example: generate and consume in one run

```
read_path: ./
read_fname: my_dem.tif
write_path: ./
write_fname: sc

make_step_change_rasters: true
step_change_uplift_old_m_yr: 0.0001
step_change_uplift_new_m_yr: 0.0005
step_change_target_relief_metres: 1200
step_change_age_yr: 500000
step_change_n: 1.0
step_change_channel_pixels: 500
write_step_change_rasters: true

calculate_multi_isotope_step_change: true
```

---

## Landslide coupling

Superimposes landslide scars on a step-change concentration field, using the
time-of-failure and excavation-depth rasters written by `lsdtt-cosmo-landslides`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `landslide_coupling` | string | `none` | `none`, `scar`, `flush` or `both`. |
| `landslide_time_raster_prefix` | string | `NULL` | Prefix of `<stem>_LS_TFAIL` from `lsdtt-cosmo-landslides`. |
| `landslide_excavation_raster_prefix` | string | `NULL` | Prefix of `<stem>_LS_EXCAVATION`. |
| `landslide_report_time_yr` | float | `0` | Present time on the same clock as the failure-time raster. |
| `flush_time_threshold_yr` | float | `100` | How long a scar keeps delivering debris (used by `both`). |
| `scar_density_kg_m3` | float | `1500` | Density of material already represented by an imported excavation raster; the explicit soil-like default is independent of any area-volume preset. The generating landslide CLI uses a negative sentinel to inherit its active preset instead. |
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

The excavation depth from `lsdtt-cosmo-landslides` is a **true geometric thickness**,
because Larsen's `V = αA^γ` is fitted to *measured* scar volumes. For a shallow
landslide the material that failed is soil, at roughly 1500 kg/m³. Carrying it at
the bedrock density overstates the mass depth by a factor of about 1.77 — and
since production falls exponentially with mass depth, that is more than a factor
of two in the concentration of the material delivered.

`lsdtt-cosmo-landslides` takes its scar density from the `area_volume_preset` in play
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

Available with both `calculate_step_change_crn_concentration_raster` and
`calculate_multi_isotope_step_change`.

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
| `quartz_fertility_regolith_depth_metres` | float | `0` | Regolith thickness separating the layers. |
| `ice_mask_raster_prefix` | string | `NULL` | Present-day ice cover (non-zero = ice). |
| `ice_thickness_metres` | float | `50` | Uniform ice thickness. |
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

`lsdtt-cosmo-landslides` accepts the same mask as `ice_mask_fname` and places no
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

When a multi-isotope concentration is present in memory and
`concentration_raster_prefix` is `NULL`, accumulation uses the in-memory
Be10/Al26/C14 rasters and accumulates all three in one shared pass over the
domain (shared routing, shared weight raster).

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
| `sample_points_fname` | string | `NULL` | CSV of `latitude,longitude` at which to sample the accumulated raster. Needs `calculate_accumulated_crn_concentration: true`. |
| `sample_column_name` | string | `sample_name` | Optional sample-identifier column copied through from the points CSV. Omitted from the output if absent. |
| `threshold_contributing_pixels` | int | `1000` | Minimum contributing pixels for a node to count as a channel when snapping. Matches `lsdtt-cosmo-tool`. |
| `sample_snap_max_search_radius_pixels` | int | `8` | Maximum search radius when snapping a point to the channel network. Matches `lsdtt-cosmo-tool`. |

Output is `<write_fname>_sampled_crn_conc.csv`.

- Single-isotope accumulation writes `latitude`, `longitude`, `<nuclide>_CONC`,
  `snapped_contributing_pixels`.
- Multi-isotope accumulation writes `latitude`, `longitude`, `Be10_CONC`,
  `Al26_CONC`, `C14_CONC`, `R2610_ACC`, `R1410_ACC`,
  `snapped_contributing_pixels`.

For apparent-erosion checks in `lsdtt-cosmo-tool`, multi-isotope sampling also
writes three direct bridge CSVs:

- `<write_fname>_sampled_crn_conc_Be10.csv`
- `<write_fname>_sampled_crn_conc_Al26.csv`
- `<write_fname>_sampled_crn_conc_C14.csv`

Each contains one concentration column (`Be10_CONC`, `Al26_CONC`, `C14_CONC`
respectively) plus coordinates/snapping, so the tool can invert each isotope
without editing columns or re-running the forward model.

Snapping mirrors `lsdtt-cosmo-tool` so both tools resolve the same node.

`threshold_contributing_pixels` is also used to snap
`step_change_basin_outlet_fname` when generating a step change.

---

## Output write toggles

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_crn_concentration_raster` | bool | `true` | Write the steady-state concentration raster. Output: `<write_fname>_<nuclide>_crn_conc.tif`. |
| `write_step_change_concentration_raster` | bool | `true` | Write the step-change concentration raster. Output: `<write_fname>_<nuclide>_step_crn_conc.tif`. |
| `write_accumulated_crn_raster` | bool | `true` | Write accumulated concentration raster(s). Single-isotope mode writes `<write_fname>_<nuclide>_crn_acc.tif`; multi-isotope mode writes all three: `_Be10_crn_acc.tif`, `_Al26_crn_acc.tif`, `_C14_crn_acc.tif`. |
| `write_production_raster` | bool | `false` | Write the Stone 2000 production raster. Output: `<write_fname>_crn_prod.tif`. |

---

## lsdviztools3 interface notes

`lsdtt-cosmo-prediction` is the forward side used by lsdtt3 comparison
workflows (including point sweeps) that are then plotted by [`lsdviztools3`](https://lsdtopotools.github.io/lsdviztools3/cosmogenic_nuclides/).

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

If `calculate_accumulated_crn_concentration: true` is also set, and
`concentration_raster_prefix` is `NULL`, the same run also writes accumulated
multi-isotope rasters (`_Be10_crn_acc.tif`, `_Al26_crn_acc.tif`,
`_C14_crn_acc.tif`) and accumulated ratios (`_R2610_acc.tif`, `_R1410_acc.tif`),
and sampled points include `R2610_ACC`/`R1410_ACC`.

Note that after an erosion-rate *increase* all concentrations fall, C-14 reaches
its new value first, and the C14/Be10 ratio is therefore **suppressed** and
approaches its endpoint from below — it does not overshoot.


### Multi-isotope transient (arbitrary history)

`calculate_multi_isotope_transient` generalises the step-change mode: instead of
one change in erosion rate, every pixel follows the history in
`transient_history_fname`. A two-row history reproduces the step-change mode
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
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format: cog, gtiff, or same_as_input. | - |

### Nuclide and muon scheme

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `nuclide` | `string` | `Be10` | `param` | - | Nuclide to predict: Be10, Al26, or C14. | - |
| `muon_scheme` | `string` | `BraucherBorchers` | `param` | - | Muon production scheme. | Braucher, BraucherBorchers, Granger, Schaller, newCRONUS, cronusEmulatorSt, cronusEmulatorLSDn. The cronusEmulator schemes fit lsdtt3's three muon pathways to the CRONUS v3 tabulated model-1A muon inventory and are pressure-dependent; a pressure raster is computed automatically. |
| `scaling_scheme` | `string` | `stone` | `param` | - | Production-rate scaling scheme. | stone = Stone (2000)/Lal, the lsdtt3 default and the CAIRN/OCTOPUS scheme. This raster prediction workflow is site-based by default; the landslide generator instead defaults to uniform to preserve its legacy scalar path. lsdn = Lifton, Sato &amp; Dunai (2014) nuclide-specific scaling, interpolated from the LSDn2015 scaling grids; requires the LSDn grids in crn_parameters.nc and supports Be10, Al26 and C14 only. Both schemes are evaluated per pixel, including when sample_age_yr &gt; 0. muon_scheme is unaffected: LSDn is a spallation scaling and the muon pathways receive it through the same Vermeesch redistribution as Stone. |

### Prediction modes

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `calculate_crn_concentration_raster` | `bool` | `false` | `switch` | - | Compute steady-state CRN concentration at every pixel. | Requires: erosion raster (or effective_erosion_rate), production raster, shielding. |
| `calculate_step_change_crn_concentration_raster` | `bool` | `false` | `switch` | - | Compute transient CRN concentration after a step change in erosion rate. | Requires: time_since_change raster, step_change_eff_erate_old/new, production raster, topo shielding. |
| `calculate_multi_isotope_step_change` | `bool` | `false` | `switch` | - | Carry Be-10, Al-26 and C-14 through the SAME step change at every pixel, and write the two isotope-ratio rasters. | The isotope RATIOS are what carry a transient signal: the half-lives differ by orders of magnitude (C-14 5.7 kyr, Al-26 0.7 Myr, Be-10 1.39 Myr), so a pair dates the perturbation in a way one nuclide alone cannot. A ratio is only meaningful if both members saw an identical scenario, so all three are computed from ONE particle evolution per pixel -- same shielding, same rates, same elapsed time -- rather than from three separate runs. Each nuclide is bitwise identical to what calculate_step_change_crn_concentration_raster returns for the same inputs, and costs little more than a single-nuclide run. Uses the same inputs as that mode: time_since_change raster, step_change_eff_erate_old/new, production raster, topographic shielding. The nuclide parameter is ignored here -- all three are always produced. |
| `calculate_multi_isotope_transient` | `bool` | `false` | `switch` | - | Carry Be-10, Al-26 and C-14 through an ARBITRARY erosion history at every pixel, and write the two isotope-ratio rasters. | The general case of calculate_multi_isotope_step_change: instead of one step, the pixels follow a user-supplied history from transient_history_fname. A two-row history reproduces the step-change mode exactly (bitwise), so this is a strict generalisation. COST: every pixel walks the whole history, so this is O(pixels x history rows) and is much slower than the step-change mode -- inherent, since an arbitrary history has no closed form. The nuclide parameter is ignored; all three are always produced. |
| `transient_rock_density_kg_m3` | `float` | `2650` | `param` | `calculate_multi_isotope_transient` | Rock density (kg/m3) used to convert the erosion history's metre-based rates and removals into effective depths. | Only used by calculate_multi_isotope_transient. The history file gives rates in m/yr and removals in m; the physics works in g/cm2, and this is the conversion. |
| `transient_history_fname` | `string` | `NULL` | `param` | `calculate_multi_isotope_transient` | CSV erosion history for calculate_multi_isotope_transient. | NULL is deliberate because transient prediction is an optional mode and the file becomes required only when that mode is enabled. Three columns: time,rate,removal -- time in yr BP, rate in m/yr, removal in m. Rows run from the MOST RECENT (index 0) to the oldest, and the LAST row must have time = -9999 carrying the long-term background rate. A row with rate = 0 and a non-zero removal is an instantaneous stripping event and must share the time of the row before it. An unordered record is rejected: it would produce negative timesteps and a silently wrong answer. |
| `calculate_accumulated_crn_concentration` | `bool` | `false` | `switch` | - | Accumulate CRN concentration downstream using D8 flow routing. | Requires a concentration raster (computed in the same run or loaded separately). |

### Production raster

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `load_production_raster` | `bool` | `false` | `param` | `write_production_raster` | Load the Stone 2000 production raster from disk instead of computing it. | File: read_path + production_raster_prefix + .tif |
| `production_raster_prefix` | `string` | `NULL` | `param` | `write_production_raster` | Filename prefix (no extension) for loading a pre-computed production raster. | Used when load_production_raster is true. |
| `write_production_raster` | `bool` | `false` | `switch` | - | Write the production raster to a GeoTIFF. | Output: write_path + write_fname + _crn_prod.tif |
| `crn_nc_fname` | `string` | - | `param` | `write_production_raster` | Path to crn_parameters.nc for ERA40 atmospheric scaling and palaeomagnetic correction. | Leave empty to use the crn_parameters.nc bundled with lsdtt3 (found automatically); if absent, standard-atmosphere Stone 2000 (no palaeomagnetic correction). |
| `sample_age_yr` | `float` | `0` | `param` | `write_production_raster` | Sample exposure age (yr BP) driving the time-integrated production correction; what it selects depends on scaling_scheme. | Only used when crn_nc_fname is set and &gt; 0. Under scaling_scheme=stone this is the Muscheler VDM palaeomagnetic correction to Stone 2000 scaling. Under scaling_scheme=lsdn it drives the duration-weighted average of LSDn production over the Pavon-Carrasco geomagnetic record (Holocene moving paleopole, then a GAD dipole record to ~20 Ma), making it longitude-dependent as well as latitude-dependent. NOTE: "now" is assumed to be 2010 (the packed record's reference year); lsdtt3 has no sample-collection-year input yet, so no (collection year - 2010) shift is applied - negligible (&lt;=~16 yr against 1e2-1e7 yr records) but an assumption, not an exact treatment. |
| `atmosphere` | `string` | `ERA40` | `param` | `write_production_raster` | Atmosphere model: ERA40 (default), NCEP, or standard. | ERA40 = Balco/Lifton ERA40 climatology + latitude-dependent lapse (the lsdtt3 default). NCEP = Balco NCEPatm_2, the atmosphere CAIRN/OCTOPUS use; requires the NCEP grids in crn_parameters.nc. standard = SLHL barometric formula, no climatology. ERA40/NCEP fall back to standard if the required grids are absent. MATCH THIS TO lsdtt-cosmo-tool when using the two together: the atmosphere sets the pressure at every pixel, so forward-modelling under one and inverting under another silently biases the recovered erosion rate. |

### Topographic shielding

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `topographic_shielding_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of topographic shielding raster (dimensionless, 0-1). | If NULL, uniform value topo_shield_value is used. |
| `topo_shield_value` | `float` | `1` | `param` | - | Uniform topographic shielding factor (dimensionless, 0-1; 1 = unshielded). | Used when topographic_shielding_raster_prefix is NULL. |

### Snow shielding

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `snow_shielding_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of snow shielding raster (effective depth, g/cm2). | If NULL, uniform value snow_shield_eff_thickness is used. |
| `snow_shield_eff_thickness` | `float` | `0` | `param` | - | Uniform snow shielding as effective depth (g/cm2; 0 = no snow). | Used when snow_shielding_raster_prefix is NULL. |

### Self shielding

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `self_shielding_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of self/sample shielding raster (effective depth, g/cm2). | If NULL, uniform value self_shield_eff_thickness is used. |
| `self_shield_eff_thickness` | `float` | `0` | `param` | - | Uniform self shielding as effective depth (g/cm2; 0 = surface sample). | Used when self_shielding_raster_prefix is NULL. |

### Background erosion rate

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `background_erosion_rate_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of raster giving effective erosion rate (g/cm2/yr) at each pixel. | If NULL, uniform value effective_erosion_rate is used. |
| `effective_erosion_rate` | `float` | `0.001` | `param` | - | Uniform effective erosion rate (g/cm2/yr) across the raster. | Used when background_erosion_rate_raster_prefix is NULL. ~0.001 ~= 3.8 mm/kyr for rho=2650. |

### Step change parameters

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `step_change_eff_erate_old_g_cm2_yr` | `float` | `0.001` | `param` | - | Old (pre-change) effective erosion rate for step change calculation (g/cm2/yr). | - |
| `step_change_eff_erate_new_g_cm2_yr` | `float` | `0.01` | `param` | - | New (post-change) effective erosion rate for step change calculation (g/cm2/yr). | - |
| `time_since_change_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of raster giving years since the step change at each pixel. | Pixels with value 0 receive SS concentration at the new erosion rate. |
| `time_since_change_yr` | `float` | `10000` | `param` | - | Uniform time since step change (yr) used when time_since_change_raster_prefix is NULL. | - |

### Generate the step change

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `make_step_change_rasters` | `bool` | `false` | `switch` | - | Generate the erosion-rate, time-since-change and topography rasters from an analytical step change in uplift, instead of loading them. | The time-since-change field is not a free input: it is the arrival time of the incision wave, from the detachment-limited stream-power solution in chi (Royden &amp; Perron 2013; Mitchell &amp; Yanites 2019). The generated rasters feed the step-change concentration modes IN MEMORY, so no round trip through disk is needed. K is back-calculated from step_change_target_relief_metres and is echoed at run time -- read it. The generated fields cover ONE basin, so the concentration output is restricted to that basin. |
| `step_change_uplift_old_m_yr` | `float` | `0.0001` | `param` | `make_step_change_rasters` | Uplift (= steady erosion) rate before the step, in METRES per year. | Metres per year, not mm/yr. Converted to an effective erosion rate with rock_density_kg_m3, and that conversion OVERRIDES step_change_eff_erate_old_g_cm2_yr when make_step_change_rasters is on. |
| `step_change_uplift_new_m_yr` | `float` | `0.001` | `param` | `make_step_change_rasters` | Uplift rate after the step, in METRES per year. | Overrides step_change_eff_erate_new_g_cm2_yr as above. |
| `step_change_target_relief_metres` | `float` | `1000` | `param` | `make_step_change_rasters` | Relief the OLD steady-state landscape should have (m). K is back-calculated from it. | This is how a user without a calibrated K gets a landscape of the right size. Everything else -- knickpoint position, wave celerity, the whole time-since-change field -- follows from the K it implies, so the derived K is printed and should be reported. |
| `step_change_age_yr` | `float` | `100000` | `param` | `make_step_change_rasters` | Time since the step change reached the basin outlet (yr). | Every pixel gets its own, smaller, time since the wave reached IT. |
| `step_change_n` | `float` | `1` | `param` | `make_step_change_rasters` | Slope exponent n in E = K A^m S^n. | n &gt;= 1 with rising uplift gives a migrating knickpoint (a shock); n &lt; 1 gives a rarefaction fan in which the erosion rate varies continuously between the old and new values. A FALLING step swaps the two, which is handled. n = 1 is the classical case. |
| `step_change_m_over_n` | `float` | `0.45` | `param` | `make_step_change_rasters` | Concavity m/n used for the chi integral. | - |
| `step_change_A_0` | `float` | `1` | `param` | `make_step_change_rasters` | Reference drainage area for the chi integral (m2). | Any positive value is handled exactly; the A_0 dependence is carried through the steepness and the wave celerity rather than assumed away. |
| `step_change_channel_pixels` | `int` | `500` | `param` | `make_step_change_rasters` | Contributing-pixel count above which a pixel is a channel. | Below it the pixel is a hillslope and is given a critical gradient instead of the fluvial solution. Hillslopes inherit the erosion rate and time-since-change of the channel they drain to. |
| `step_change_Sc_old` | `float` | `0.25` | `param` | `make_step_change_rasters` | Critical hillslope gradient in the unadjusted part of the landscape. | - |
| `step_change_Sc_new` | `float` | `0.4` | `param` | `make_step_change_rasters` | Critical hillslope gradient in the adjusted part of the landscape. | - |
| `step_change_basin_outlet_fname` | `string` | `NULL` | `param` | `make_step_change_rasters` | CSV of latitude,longitude giving the basin outlet to model. | The first point is used and any others ignored. It is snapped to the channel network using threshold_contributing_pixels and sample_snap_max_search_radius_pixels. If NULL the largest basin in the DEM is used. |
| `rock_density_kg_m3` | `float` | `2650` | `param` | `make_step_change_rasters` | Rock density (kg/m3) converting the generated uplift rates (m/yr) to effective erosion rates (g/cm2/yr). | eff = rate_m_yr * rho_kg_m3 * 0.1. Matches lsdtt-cosmo-tool's parameter of the same name. Distinct from transient_rock_density_kg_m3, which belongs to the erosion-history file read by calculate_multi_isotope_transient. |
| `step_change_production_from_generated_dem` | `bool` | `true` | `switch` | `make_step_change_rasters` | Compute the production raster from the GENERATED topography rather than the input DEM. | On by default because the generated surface is the landscape the scenario describes: using the real DEM's elevations with the scenario's erosion rates mixes two landscapes, and production depends on elevation through pressure. The generated surface is offset so its outlet sits at the input DEM's outlet elevation, keeping the base level real. Set false to reproduce the legacy behaviour, which built production from the input DEM. |
| `step_change_route_on_generated_dem` | `bool` | `true` | `switch` | `make_step_change_rasters` | Build the accumulation flow routing on the GENERATED topography. Set false to route on the ORIGINAL input DEM instead. | Only has an effect when step_change_production_from_generated_dem is true, since that is what replaces the working DEM. The default preserves historical behaviour, but false is the more defensible choice and here is why. The generator snaps hillslopes draining below the knickpoint to step_change_Sc_new and the rest to step_change_Sc_old. Each side is built up from its own channel and nothing reconciles them at the drainage divide, so the generated surface carries an elevation CLIFF along the divide -- measured on the Meghalaya 5x case at a median 56 m, with the steepened side higher in 92 % of boundary cells. Rebuilding D8 on that surface makes divide pixels drain ACROSS the divide into the neighbouring unadjusted catchment: 5.4 % of hillslope-new pixels did so there, carrying post-knickpoint material into basins that are supposed to be entirely pre-knickpoint. A real landscape would migrate the divide rather than hold a cliff on it; representing that needs the divides iterated as the canyons steepen, which the one-pass construction cannot do. Routing on the original DEM sidesteps the artefact, at the cost of a flow network that no longer matches the modelled topography. Production, elevation and the step-change fields still come from the generated surface either way. |
| `write_step_change_rasters` | `bool` | `false` | `switch` | `make_step_change_rasters` | Write the generated topography, erosion rate, time-since-change, tag and chi rasters. | Files: &lt;write_fname&gt;_SC_DEM, _SC_ERATE, _SC_STEPTIME, _SC_TAG, _SC_CHI. The tag raster is 1 = unadjusted channel, 2 = unadjusted hillslope, 3 = adjusted channel, 4 = adjusted hillslope, 5/6 = the same inside a rarefaction fan. USE THE TAG RASTER, not chi, to ask whether a pixel has been rejuvenated: a hillslope always has a larger chi than the channel it drains to, so an adjusted hillslope legitimately sits above the knickpoint in chi. |

### Landslide coupling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `landslide_coupling` | `string` | `none` | `param` | - | How landslides modify the step-change concentration: none, scar, flush or both. | scar = the landslide material has long since left the system, and what remains is a scar eroding at the catchment rate from a column whose inheritance was stripped when it failed. flush = the debris is still in the channel, so the sample contains the whole excavated slab. both = flush inside flush_time_threshold_yr and scar outside it, which is the physically complete case. Needs landslide_time_raster_prefix and landslide_excavation_raster_prefix, both written by lsdtt-cosmo-landslides. Requires a step-change mode. |
| `landslide_time_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of the time-of-failure raster (lsdtt-cosmo-landslides &lt;stem&gt;_LS_TFAIL). | This holds an ABSOLUTE failure time on the model clock, with nodata where nothing failed. It is converted here to elapsed years using landslide_report_time_yr; 'never failed' is carried as a distinct state, NOT as nodata, so unfailed ground stays in the raster and in any basin average taken from it. |
| `landslide_excavation_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of the excavation-depth raster (lsdtt-cosmo-landslides &lt;stem&gt;_LS_EXCAVATION). | Real metres of material removed, not a rock-equivalent thickness. It is converted to a mass depth with scar_density_kg_m3. |
| `landslide_report_time_yr` | `float` | `0` | `param` | - | Present time on the same clock as the time-of-failure raster (yr). | elapsed = landslide_report_time_yr - time_of_failure. A failure at or after this time has not happened yet and is treated as 'never'. |
| `flush_time_threshold_yr` | `float` | `100` | `param` | - | How long a scar keeps delivering its debris to the channel (yr). | Used by landslide_coupling: both. This is an AGE, following Niemi et al. (2005), whose model has no sediment storage at all -- landslide material passes through within one 100 yr timestep and mixes homogeneously, which is where this default comes from. Set it to the event recurrence interval to reproduce that convention. How long landslide sediment actually persists in a channel network is a separate research question and is not modelled here. |
| `write_delivery_mass_raster` | `bool` | `false` | `switch` | - | Write the delivery-mass raster: the weight a downstream mixing calculation should use. | File: &lt;write_fname&gt;_LS_DELIVERY.tif, in kg m^-2 yr^-1. A flushing pixel delivers its whole excavated slab over the flush window; everything else delivers its background erosion. The erosion-rate raster is left untouched, so the TRUE erosion rate stays available to compare an apparent rate against. |
| `scar_density_kg_m3` | `float` | `1500` | `param` | - | Density of the material a landslide EXCAVATED (kg/m3). | A DIFFERENT quantity from rock_density_kg_m3, which is the eroding bedrock column. The excavation depth is a true geometric thickness of whatever failed, and for a shallow landslide that is soil at roughly 1500 kg/m3. Carrying it at 2650 overstates the mass depth by about a factor of 1.77 and changes the delivered concentration by more than a factor of two. Set 2650 for bedrock failures. The default matches lsdtt-cosmo-landslides' larsen_soil_scar preset. |

### Quartz fertility and ice

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `quartz_fertility` | `float` | `1` | `param` | - | Uniform quartz mass fraction (0-1) used as the downstream mixing weight. | A basin-averaged sample is a mixture of QUARTZ, not of rock. Two pixels eroding at the same rate contribute unequally if their rock carries different quartz, and a quartz-free unit (limestone, basalt, many volcanics) contributes sediment but no sample at all. In mixed-lithology terrain, routing concentrations without this is simply wrong. The default of 1 makes the whole feature a no-op: nothing changes unless you opt in. |
| `quartz_fertility_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of a raster of quartz mass fraction (0-1), typically a reclassified lithology map. | Overrides quartz_fertility where present. Values outside [0,1] are clamped and the number clamped is reported, so a map supplied as percentages rather than fractions surfaces rather than being silently truncated. |
| `quartz_fertility_regolith` | `float` | `-1` | `param` | - | Quartz fraction of the regolith, for the two-layer model. Negative = off. | Chemical weathering removes soluble minerals and leaves quartz behind, so regolith is enriched relative to the fresh bedrock a deep landslide scar exposes. Set this, quartz_fertility_bedrock and quartz_fertility_regolith_depth_metres together to test whether that makes landslide debris under-represented in the quartz separate, and so weakens the dilution a landslide causes. |
| `quartz_fertility_bedrock` | `float` | `-1` | `param` | - | Quartz fraction of fresh bedrock, for the two-layer model. Negative = off. | - |
| `quartz_fertility_regolith_depth_metres` | `float` | `0` | `param` | - | Regolith thickness (m) separating the two fertility layers. 0 = off. | - |
| `ice_mask_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of a raster marking present-day ice cover (non-zero = ice). | Ice changes exactly three things and lsdtt3 already has a field for each: the ice column shields cosmic rays (SNOW SHIELDING, as an effective depth), the erosion rate differs or is zero, and little or no quartz reaches the sample. LIMITATION: this is PRESENT-DAY ice only. Concentrations depend on the shielding HISTORY and lsdtt3 cannot yet express a per-pixel time-varying shielding, so a glacial-interglacial cycle is not representable and concentrations under formerly glaciated ground are overestimated. |
| `ice_thickness_metres` | `float` | `50` | `param` | - | Uniform ice thickness (m) used where the mask marks ice. | Converted to a snow-shielding mass depth at ice_density_kg_m3. Ten metres of ice is already essentially total shielding (0.3 % of spallation survives). |
| `ice_density_kg_m3` | `float` | `917` | `param` | - | Density of ice (kg/m3). Lower for firn. | - |
| `erate_under_ice_g_cm2_yr` | `float` | `0` | `param` | - | Effective erosion rate beneath the ice (g/cm2/yr). | Zero is the cold-based idealisation and is numerically safe: the steady-state kernel forms its denominator as (eps + Gamma*lambda), never divides by eps, and returns the saturation concentration P/lambda. |
| `quartz_fertility_under_ice` | `float` | `0` | `param` | - | Quartz fraction beneath the ice. | Zero means the ice contributes nothing to the sample, which is the usual assumption when the erosion rate under it is also zero. |
| `write_quartz_weight_raster` | `bool` | `false` | `switch` | - | Write the quartz mixing-weight raster. | File: &lt;write_fname&gt;_QUARTZ_WEIGHT.tif, kg of quartz m^-2 yr^-1. This is what the accumulation mixes on when fertility or ice is active. |

### Concentration raster loading (for accumulation)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `concentration_raster_prefix` | `string` | `NULL` | `param` | - | Prefix of a pre-computed concentration raster to use for accumulation. | If NULL and a concentration raster was computed in the same run, that is used. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_before_routing` | `bool` | `true` | `switch` | - | Choose whether to pit-fill the DEM before building the flow network. | ASKS: should the DEM be filled at all? Filling is OPTIONAL here and this is a genuine analytical choice. True (the default) routes over a pit-free surface. False builds the flow network on the UNFILLED surface, so flow paths terminate at sinks instead of being carried through them. Do not confuse this with raster_is_filled, which is an assertion that the input is already filled. |

### Flow routing (for accumulation)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_slope_for_fill` | `float` | `0.0001` | `param` | - | Minimum imposed slope during pit filling (m/m). | - |
| `boundary_conditions` | `string` | `oooo` | `param` | - | D8 flow routing boundary conditions (N,E,S,W): o=open, b=forced outlet. | - |

### Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `sample_points_fname` | `string` | `NULL` | `param` | - | CSV of latitude,longitude points at which to sample the accumulated CRN raster. | When set (not NULL), each point is snapped to the channel network and the accumulated (basin-averaged) result at the snapped node is written to &lt;write_fname&gt;_sampled_crn_conc.csv. Single-isotope mode writes &lt;nuclide&gt;_CONC; multi-isotope mode writes Be10_CONC, Al26_CONC, C14_CONC and the accumulated ratios R2610_ACC and R1410_ACC, plus one per-nuclide bridge CSV (&lt;write_fname&gt;_sampled_crn_conc_Be10.csv etc.) for direct inversion with lsdtt-cosmo-tool. The CSV is directly usable by lsdtt-cosmo-tool via concentration_column_name (for example Be10_CONC). Requires calculate_accumulated_crn_concentration: true. |
| `sample_points_stride` | `int` | `1` | `param` | - | Keep every Nth point from sample_points_fname (1 = keep all). | For thinning a dense channel-node CSV to a manageable sample set. Points are thinned BEFORE snapping, in input order. Values &lt; 1 are treated as 1. |
| `sample_column_name` | `string` | `sample_name` | `param` | - | Optional sample identifier column copied through from the points CSV. | If absent from the points CSV, the column is omitted from the output. |
| `sample_snap_max_search_radius_pixels` | `int` | `8` | `param` | - | Maximum search radius (pixels) when snapping sample points to the channel network. | Matches lsdtt-cosmo-tool. Used only when sample_points_fname is set. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Output write toggles

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_crn_concentration_raster` | `bool` | `true` | `switch` | - | Write the computed CRN concentration raster. | Output: write_fname + _&lt;nuclide&gt;_crn_conc.tif |
| `write_step_change_concentration_raster` | `bool` | `true` | `switch` | - | Write the step-change CRN concentration raster. | Output: write_fname + _&lt;nuclide&gt;_step_crn_conc.tif |
| `write_isotope_ratio_rasters` | `bool` | `true` | `switch` | `calculate_multi_isotope_step_change` | Write the Al26/Be10 and C14/Be10 ratio rasters. | In multi-isotope step-change/transient mode writes _R2610_step.tif and _R1410_step.tif. If accumulation is also on, writes _R2610_acc.tif and _R1410_acc.tif from the accumulated basin-mean concentrations. |
| `write_accumulated_crn_raster` | `bool` | `true` | `switch` | - | Write the accumulated CRN concentration raster. | Output: write_fname + _&lt;nuclide&gt;_crn_acc.tif |
| `write_upstream_mean_erate_raster` | `bool` | `false` | `switch` | - | Write the arithmetic upstream mean erosion rate (g/cm2/yr) -- the basin's ACTUAL mean erosion rate above each pixel, which is what an erosion rate inverted from the accumulated concentration should be compared against. Sampled to CSV as upstream_mean_eff_erate_g_cm2_yr regardless of this switch. | Output: write_fname + _mean_erate_acc.tif |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

<!-- END GENERATED CLI PARAMETERS -->
