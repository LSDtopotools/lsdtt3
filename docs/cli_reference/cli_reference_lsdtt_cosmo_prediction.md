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
| 2. Production raster | Compute Stone 2000 or load from disk |
| 3. Shielding rasters | Topographic, snow, self (raster or uniform scalar) |
| 4. Erosion raster | Raster or uniform scalar |
| 1b. Step-change generation | `make_step_change_rasters` — derives erosion rate, time-since-change and topography from the analytical fluvial solution |
| 5. Steady-state concentrations | `calculate_crn_concentration_raster` |
| 6. Step-change concentrations | `calculate_step_change_crn_concentration_raster` |
| `calculate_multi_isotope_step_change` | bool | `false` | Carry Be-10, Al-26 and C-14 through the **same** step change at every pixel, and write the two isotope-ratio rasters. See below. |
| `calculate_multi_isotope_transient` | bool | `false` | The same three nuclides through an **arbitrary erosion history** from `transient_history_file`. The general case of the step-change mode. See below. |
| `transient_history_file` | string | `NULL` | CSV erosion history: `time,rate,removal` (yr BP, m/yr, m). Rows run most-recent-first; the last row must be `time = -9999` carrying the background rate. Unordered records are rejected. |
| `transient_rock_density_kg_m3` | float | `2650` | Density used to convert the history's metre-based rates and removals into effective depths. |
| `write_isotope_ratio_rasters` | bool | `true` | Write `_R2610_step.tif` (Al26/Be10) and `_R1410_step.tif` (C14/Be10). Needs `calculate_multi_isotope_step_change`. |
| 7. Accumulated concentrations | `calculate_accumulated_crn_concentration` |

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
| `calculate_accumulated_crn_concentration` | bool | `false` | Accumulate concentration downstream using D8 flow routing. Requires a concentration raster (computed in the same run or loaded separately). |

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
