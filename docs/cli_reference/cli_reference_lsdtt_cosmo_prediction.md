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
| `muon_scheme` | string | `BraucherBorchers` | Muon production scheme. `Braucher`, `BraucherBorchers`, `Granger`, `Schaller`, or `newCRONUS`. |

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
