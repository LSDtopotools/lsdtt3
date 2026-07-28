# CLI reference: lsdtt-cosmo-tool

Inverts measured cosmogenic nuclide concentrations to basin-averaged effective erosion rates using the production–accumulation–inversion approach. Loads a DEM, builds D8 flow routing, computes a Stone 2000 production raster, snaps sample points to the channel network, and writes a CSV of erosion rate estimates with optional uncertainty bounds.

The DEM must be in a projected CRS (e.g. UTM). At least one sample must be successfully snapped to a channel.

## Usage

```
lsdtt-cosmo-tool <path> <parameter_file>
lsdtt-cosmo-tool <parameter_file>         # path defaults to ./
lsdtt-cosmo-tool -h                       # print help and exit
lsdtt-cosmo-tool -v                       # print version and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

| Step | Description |
|---|---|
| 1 | Load DEM; check projected CRS |
| 2 | Pit fill (skipped when `fill_before_routing: false`) |
| 3 | D8 flow routing |
| 4 | Stone 2000 production raster |
| 5 | Load sample points from CSV; snap to channel |
| 6 | For each sample: apply standardisation, invert to erosion rate, compute uncertainty bounds |
| 7 | Write output CSV |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename (without extension). |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

The output CSV is always written as `<write_fname>_CN_ERATE.csv`.

---

## CRN physics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `nuclide` | string | `Be10` | Nuclide used for erosion rate inversion. Options: `Be10`, `Al26`, `C14`. |
| `muon_scheme` | string | `BraucherBorchers` | Muon production parameter scheme. Options: `BraucherBorchers`, `Braucher`, `Granger`, `Schaller`, `newCRONUS`. |
| `rock_density_kg_m3` | float | `2650` | Rock bulk density (kg/m³). Used to convert effective erosion rate (g/cm²/yr) to mm/kyr. |
| `scaling_scheme` | string | `stone` | Production-rate scaling. `stone` = Stone (2000)/Lal, the default and the CAIRN/OCTOPUS scheme. `lsdn` = Lifton, Sato & Dunai (2014) nuclide-specific scaling, interpolated from the LSDn2015 grids; needs the LSDn grids in `crn_parameters.nc` and supports `Be10`, `Al26`, `C14` only. Both are evaluated **per pixel**, including when `sample_age_yr > 0`. `muon_scheme` is unaffected. |

---

## Atmosphere and palaeomagnetic correction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `crn_nc_file` | string | `` | Path to `crn_parameters.nc` (ERA40 grids + Muscheler VDM record), for the ERA40 atmosphere and palaeomagnetic correction. Leave empty to use the `crn_parameters.nc` bundled with lsdtt3 (found automatically relative to the executable); if that is missing too, the standard (SLHL) atmosphere is used with no correction. |
| `sample_age_yr` | float | `0` | Sample exposure age (yr BP). Drives the time-integrated production correction; what it selects depends on `scaling_scheme`. Under `stone` it is the Muscheler VDM palaeomagnetic correction; under `lsdn` it is the duration-weighted average of LSDn production over the Pavón-Carrasco record, so it is longitude-dependent as well as latitude-dependent. `0` = modern field, no correction. Requires the record in `crn_parameters.nc`. "Now" is taken as 2010, the reference year of the packed record. |
| `atmosphere` | string | `ERA40` | `ERA40` (Balco/Lifton climatology + latitude-dependent lapse, the lsdtt3 default), `NCEP` (the atmosphere CAIRN/OCTOPUS use; needs the NCEP grids), or `standard` (SLHL barometric formula). `ERA40`/`NCEP` fall back to `standard` if their grids are absent. |
| `self_consistent_age` | bool | `false` | Solve the sample age and the erosion rate **together**, iterating to the erosional integration timescale `t = Γ_sp/ε` rather than holding the age at `sample_age_yr`. This is the CRONUScalc / Balco v3 convention, and it is how most published LSDn rates were computed — see below. |
| `self_consistent_age_tol` | float | `0.001` | Relative convergence tolerance on ε for the iteration. The default sits just inside the inversion solver's own round-trip tolerance, so tightening it buys precision the forward model cannot deliver. |
| `self_consistent_age_max_iter` | int | `10` | Iteration cap. A sample that has not converged is reported with a warning and flagged `self_consistent_converged = 0` in the output CSV; its last iterate is written, never silently accepted. |

### Self-consistent age

A sample's nuclide inventory accumulates over roughly one attenuation length of
erosion, `τ = Γ_sp/ε` (Γ_sp = 160 g/cm² in every muon scheme lsdtt3 ships). The
geomagnetic field it actually integrated is therefore the field over *that*
window — but τ depends on ε, the very thing being solved for. With
`self_consistent_age: true` the tool closes the loop: invert at `sample_age_yr`,
form `t = Γ_sp/ε`, rebuild production, re-invert, repeat until ε stops moving.
`sample_age_yr` becomes only the starting guess (`0` = start from the modern
field). Typically 2–4 iterations.

It is **off by default** for two reasons. It changes the meaning of the reported
rate, so it should be a deliberate choice; and it costs one production-raster
rebuild per iteration **per sample** — the age, and hence the whole raster,
differs sample to sample. With `scaling_scheme: lsdn` on a large DEM this
dominates the run time.

Turn it on when comparing against published LSDn or CRONUScalc rates, which are
computed this way. The effect is largest where the past field departed furthest
from modern — low latitudes, and slowly eroding samples with long τ. On the
bundled tiny fixture the converged rate sits ~9 % above the modern-field rate.

Three extra columns appear in the output CSV: `self_consistent_t_yr` (the
converged age), `self_consistent_n_iter`, and `self_consistent_converged` (1/0).
Uncertainty bounds are inverted against the *central* sample's converged
production raster, matching how the reference calculators propagate
concentration uncertainty with scaling held fixed.

---

## Sample points

| Parameter | Type | Default | Description |
|---|---|---|---|
| `points_filename` | string | `CRN_points.csv` | CSV file containing sample locations and concentrations. Must have `latitude` and `longitude` columns (EPSG:4326). |
| `sample_column_name` | string | `sample_name` | Column name for sample identifiers. Optional; omitted from output if not present. |
| `concentration_column_name` | string | `` | Column name for measured nuclide concentration (atoms/g). Defaults to `<nuclide>_CONC` if not set. |
| `uncertainty_column_name` | string | `` | Column name for 1-sigma concentration uncertainty (atoms/g). Defaults to `<nuclide>_uncertainty`. If absent, uncertainty bounds are skipped. |
| `standardisation_column_name` | string | `standardisation` | Column name for AMS standardisation name. If absent, no standardisation correction is applied. |

---

## Shielding

Shielding values are applied uniformly across the entire basin.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `self_shield_eff_thickness` | float | `0` | Sample self-shielding as effective depth (g/cm²). `0` = surface sample. |
| `snow_shield_eff_thickness` | float | `0` | Snow shielding as effective depth (g/cm²). `0` = no snow correction. |
| `topo_shield_value` | float | `1` | Topographic shielding factor (dimensionless, 0–1). `1` = no topographic shielding. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fill_before_routing` | bool | `true` | Fill DEM pits before building D8 flow routing. Set `false` only if the DEM is already pit-filled. |
| `min_slope_for_fill` | float | `0.0001` | Minimum imposed slope during pit filling (m/m). Used only when `fill_before_routing: true`. |
| `boundary_conditions` | string | `oooo` | D8 flow routing boundary conditions. Four-character string (N,E,S,W): `o` = open, `b` = forced outlet. |

---

## Channel snapping

Sample points are snapped to the nearest channel pixel meeting the contributing-pixel threshold using an expanding square search.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold_contributing_pixels` | int | `1000` | Minimum contributing pixels for a node to be considered a channel. |
| `snap_search_radius_pixels` | int | `8` | Maximum search radius (pixels) when snapping points to a channel. |

---

## Output

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_production_raster` | bool | `false` | `<write_fname>_PROD.tif` — Stone 2000 production scaling raster. |

The erosion rate CSV `<write_fname>_CN_ERATE.csv` is always written. Columns include latitude, longitude, measured concentration, predicted concentration, effective erosion rate (g/cm²/yr), and erosion rate (mm/kyr). Uncertainty bound columns are added when an uncertainty column is present in the input.

---

## lsdviztools3 interface notes

`lsdviztools3` can plot comparison figures from the standard cosmogenic CSV
workflows that include `lsdtt-cosmo-tool` outputs.

- Basin/sample comparison tiers (`run_cosmo_validation.py`,
  `run_cosmo_riversand_comparison.py`) produce summary CSVs under
  `tests/expected/cosmo/published/` such as:
  - `<dataset>_summary.csv`
  - `<dataset>_lsdn_summary.csv`
  - `<dataset>_riversand_summary.csv`
- `lsdviztools3` command:

```
lsdviz3-cosmo-published \
  --published-dir tests/expected/cosmo/published \
  --datasets dethier2014,clapp2002,dibiase2010 \
  --mode auto \
  --outdir demo_output/cosmo
```

For sweep-style lsdtt3-vs-reference figures, see the single-point scan script
`tests/scripts/run_cosmo_single_point_scaling_scan.py` and
`lsdviz3-cosmo-compare`.

---

## Example parameter file

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_samples

nuclide:       Be10
muon_scheme:   BraucherBorchers
rock_density_kg_m3: 2650

fill_before_routing:             true
min_slope_for_fill:              0.0001
boundary_conditions:             oooo
threshold_contributing_pixels:   1000
snap_search_radius_pixels:       8

points_filename:                 CRN_points.csv
concentration_column_name:       Be10_CONC
uncertainty_column_name:         Be10_uncertainty
standardisation_column_name:     standardisation
sample_column_name:              sample_name

write_production_raster: false
```
