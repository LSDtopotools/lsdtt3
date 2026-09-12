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
| 4 | Pressure raster (from `atmosphere`), then the production raster — Stone 2000 or LSDn per `scaling_scheme` |
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
| `muon_scheme` | string | `BraucherBorchers` | Muon production parameter scheme. Options: `BraucherBorchers`, `Braucher`, `Granger`, `Schaller`, `newCRONUS`, `cronusEmulatorSt`, `cronusEmulatorLSDn`. The two `cronusEmulator*` schemes fit lsdtt3's three muon pathways to the CRONUS v3 **tabulated** model-1A muon inventory — **use them when the goal is agreement with the CRONUS v3 online calculator** (measured mean −0.23 % on Be-10; see the benchmarking doc). They are pressure-dependent, so a pressure raster is computed automatically, and nuclide-specific. Pair `cronusEmulatorSt` with `scaling_scheme: stone` and `cronusEmulatorLSDn` with `scaling_scheme: lsdn`; the pairing is enforced by the spallation reference production, which differs between the two. `cronusEmulatorLSDn` carries a larger residual than the St variant for structural reasons documented in the benchmarking doc. |
| `rock_density_kg_m3` | float | `2650` | Rock bulk density (kg/m³). Used to convert effective erosion rate (g/cm²/yr) to mm/kyr. |
| `scaling_scheme` | string | `stone` | Production-rate scaling. `stone` = Stone (2000)/Lal, the default and the CAIRN/OCTOPUS scheme. `lsdn` = Lifton, Sato & Dunai (2014) nuclide-specific scaling, interpolated from the LSDn2015 grids; needs the LSDn grids in `crn_parameters.nc` and supports `Be10`, `Al26`, `C14` only. Both are evaluated **per pixel**, including when `sample_age_yr > 0`. `muon_scheme` is unaffected. |

---

## Atmosphere and palaeomagnetic correction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `crn_nc_fname` | string | `` | Path to `crn_parameters.nc` (ERA40 grids + Muscheler VDM record), for the ERA40 atmosphere and palaeomagnetic correction. Leave empty to use the `crn_parameters.nc` bundled with lsdtt3 (found automatically relative to the executable); if that is missing too, the standard (SLHL) atmosphere is used with no correction. |
| `sample_age_yr` | float | `0` | Sample exposure age (yr BP). Drives the time-integrated production correction; what it selects depends on `scaling_scheme`. Under `stone` it is the Muscheler VDM palaeomagnetic correction; under `lsdn` it is the duration-weighted average of LSDn production over the Pavón-Carrasco record, so it is longitude-dependent as well as latitude-dependent. `0` = modern field, no correction — the default, so **no time-series correction is applied unless you ask for it**; see below. Requires the record in `crn_parameters.nc`. "Now" is taken as 2010, the reference year of the packed record. |
| `atmosphere` | string | `ERA40` | `ERA40` (Balco/Lifton climatology + latitude-dependent lapse, the lsdtt3 default), `NCEP` (the atmosphere CAIRN/OCTOPUS use; needs the NCEP grids), or `standard` (SLHL barometric formula). `ERA40`/`NCEP` fall back to `standard` if their grids are absent. |
| `self_consistent_age` | bool | `false` | Solve the sample age and the erosion rate **together**, iterating to the erosional integration timescale `t = Γ_sp/ε` rather than holding the age at `sample_age_yr`. This is the CRONUScalc / Balco v3 convention, and it is how most published LSDn rates were computed — see below. Leaving it `false` while `muon_scheme: cronusEmulatorLSDn` raises an advisory warning, because that combination cannot reproduce CRONUS v3: measured bias at Sulawasi is +10.9 % at `sample_age_yr: 0` against +0.5 % with the iteration on. The warning is advisory — a fixed-age LSDn run is perfectly legitimate, it is just not the route to v3 agreement. |
| `self_consistent_age_tolerance` | float | `0.001` | Relative convergence tolerance on ε for the iteration. The default sits just inside the inversion solver's own round-trip tolerance, so tightening it buys precision the forward model cannot deliver. |
| `self_consistent_age_max_iter` | int | `10` | Iteration cap. A sample that has not converged is reported with a warning and flagged `self_consistent_converged = 0` in the output CSV; its last iterate is written, never silently accepted. |

### No time-series correction is applied by default

**Time-integrated production is opt-in, not the default.** It is switched on by
*either* `sample_age_yr > 0` *or* `self_consistent_age: true`, and both are off
out of the box. With the defaults, production is evaluated at the present-day
field under either scaling scheme — `stone` applies no Muscheler VDM correction,
and `lsdn` reads record step 0 with no time loop. Every run reports which it did:

```
Palaeomag     : off (modern field)
```

Check that line if you are unsure what a given parameter file produced.

### CAIRN-comparable run

CAIRN applied no palaeomagnetic time-series correction, so reproducing that
convention means leaving the time settings alone and matching the scaling and
atmosphere:

```
muon_scheme:         Braucher   # CAIRN's default muon scheme
scaling_scheme:      stone      # Stone (2000)/Lal — the CAIRN/OCTOPUS scaling
atmosphere:          NCEP       # the atmosphere CAIRN/OCTOPUS use
sample_age_yr:       0          # default — modern field, no correction
self_consistent_age: false      # default — no iteration to t = Γ_sp/ε
topo_shield_value:   1          # default — see the shielding caveat below
```

Two of these differ from the lsdtt3 defaults: `muon_scheme` (lsdtt3 defaults to
`BraucherBorchers`) and `atmosphere` (lsdtt3 defaults to `ERA40`). Both change
the answer on their own, so set both explicitly for a like-for-like comparison.
The atmosphere difference is a pressure difference, not a time-series one.

**Topographic shielding is the one convention lsdtt3 deliberately does not
follow.** CAIRN applies a topographic shielding correction; lsdtt3 defaults to
`topo_shield_value: 1`, i.e. none, following DiBiase (2018) — see
[Why topographic shielding defaults to off](#why-topographic-shielding-defaults-to-off).
This is a *deliberate methodological difference from CAIRN, not an omission*,
and it is why a CAIRN-convention run and an lsdtt3 default run will not agree
exactly even with everything above matched. Set `topo_shield_value` to a basin
shielding factor if you want CAIRN's convention.

Finally, this note fixes the *convention*, not the numerical agreement: for
validated per-sample targets against the CAIRN datasets (Bierman 2005, Dethier
2014, Kirchner 2001, Munack, …) see `docs/testing_and_benchmarking.md`, which
also records that CAIRN used SRTM 90 m throughout, so DEM resolution is a
further difference to control.

### Self-consistent age

**It applies to both scaling schemes, not just LSDn.** The iteration rebuilds the
production raster at the converged age whatever `scaling_scheme` is set to, and
what the age changes is scheme-specific — under `lsdn` it is the duration-weighted
LSDn average, under `stone` it is the Muscheler VDM palaeomagnetic correction.
Measured on the bundled tiny fixture, turning it on moves the rate **+9.0 % under
`stone`** (19.96 → 21.76 mm/kyr) and **+8.7 % under `lsdn`** (17.72 → 19.26
mm/kyr).

The one case where it does nothing is when `crn_parameters.nc` is unavailable:
the palaeomagnetic correction needs the VDM record, so with no NetCDF file the
age has no effect under `stone` and the iteration converges immediately at the
modern-field rate. (`lsdn` fails outright without the file, rather than falling
back.)

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
computed this way — but note that "computed this way" is a property of the
reference calculator, not of the scaling scheme, so the same argument applies to
a `stone` run being compared against CRONUScalc. The effect is largest where the
past field departed furthest from modern — low latitudes, and slowly eroding
samples with long τ.

Three extra columns appear in the output CSV: `self_consistent_t_yr` (the
converged age), `self_consistent_n_iter`, and `self_consistent_converged` (1/0).
Uncertainty bounds are inverted against the *central* sample's converged
production raster, matching how the reference calculators propagate
concentration uncertainty with scaling held fixed.

---

## Sample points

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sample_points_fname` | string | `CRN_points.csv` | CSV file containing sample locations and concentrations. Must have `latitude` and `longitude` columns (EPSG:4326). |
| `sample_column_name` | string | `sample_name` | Column name for sample identifiers. Optional; omitted from output if not present. |
| `concentration_column_name` | string | `` | Column name for measured nuclide concentration (atoms/g). Defaults to `<nuclide>_CONC` if not set. |
| `uncertainty_column_name` | string | `` | Column name for 1-sigma concentration uncertainty (atoms/g). Defaults to `<nuclide>_uncertainty`. If absent, uncertainty bounds are skipped. |
| `standardisation_column_name` | string | `standardisation` | Column name for AMS standardisation name. If absent, no standardisation correction is applied. Supply the concentration **as measured on its own standard** — see the warning below. |

> ### ⚠ Standardisation is applied by lsdtt3 — do not pre-normalise
>
> When `standardisation_column_name` names a column that is present and
> non-empty, lsdtt3 multiplies both the concentration and its uncertainty by the
> AMS standardisation factor for that standard (e.g. Be-10 `KNSTD` → `07KNSTD`,
> factor 0.9042). The `concentration_*` value written to the output CSV is the
> **normalised** one.
>
> Two consequences, and both have caused real errors:
>
> **1. Give the tool concentrations as your lab reported them.** If you
> normalise to 07KNSTD yourself *and* leave the standardisation column
> populated, the factor is applied twice. A 10 % error in concentration
> propagates straight into the erosion rate. To supply
> already-normalised values, either drop the column or set
> `standardisation_column_name` to a name that does not exist in the file.
>
> **2. Do not feed the output CSV's `concentration_*` back into another
> calculator together with the standard name.** It is already normalised;
> pairing it with the standard name asks the second calculator to normalise it
> again. This is exactly what happened in the OCTOPUS comparison tier — the
> concentrations arrived already on the 07KNSTD basis and were corrected a
> second time, which inverted the apparent sign of the residual (clapp2002
> read +6.6 % instead of −4.09 %) and sent the investigation after a
> non-existent physics discrepancy. See `docs/testing_and_benchmarking.md`.
>
> The CRONUS v3 input file that this tool writes handles the second case
> correctly: it carries the **as-measured** concentration with the standard
> name, so the calculator applies the factor exactly once.
>
> Known exception: `tests/scripts/compare_cosmo_tool_vs_cronus_v3.py` reads the
> *normalised* `concentration_*` from the output CSV, so passing it
> `--standardisation-col` reintroduces the double-count. Prefer the CLI's own
> `<write_fname>_CRONUS_V3_INPUT.txt`.

---

## Uncertainty propagation

**The internal/external split follows the CRONUS v3 convention.** CRONUS
reports an *internal* (analytical, measurement-only) uncertainty and an
*external* uncertainty that additionally carries the production-rate calibration
and decay-constant terms. lsdtt3 uses the same two-tier structure and the same
names, so `apparent_e_unc_int_mm_kyr` is comparable to the CRONUS internal
uncertainty and `apparent_e_unc_ext_mm_kyr` to the external one.

The consequence is that **lsdtt3's `unc_ext` is only comparable to a CRONUS
external uncertainty once you switch the extra terms on** — see the two
`propagate_*` parameters below. Left at their defaults, `unc_ext` equals
`unc_int` and both are concentration-only, so quoting `unc_ext` as an external
uncertainty would understate it. This is the single most common way to
under-report absolute uncertainty from this tool.

Uncertainty bounds are only produced when the input CSV carries the column named
by `uncertainty_column_name` (default `<nuclide>_uncertainty`). The concentration
sigma is scaled by the same AMS standardisation factor as the concentration
itself before it is used — so it must be supplied **as measured**, on the same
basis as the concentration. See the standardisation warning under
[Sample points](#sample-points).

The **internal** uncertainty is a re-inversion at each concentration bound:
`unc_int = ½ |ε(C+σ) − ε(C−σ)|`. Both bounds are inverted against the *central*
sample's production raster, with scaling held fixed — this matches how the
reference calculators propagate concentration uncertainty, and keeps the bounds
an uncertainty on the concentration rather than a second coupled solve.

The **external** uncertainty adds the terms below in quadrature:
`unc_ext = √(unc_int² + (ε·production_rel_sigma)² + (ε·decay_rel_sigma)²)`.

Both switches are **off by default**, so `unc_ext` equals `unc_int` unless you
turn them on. Reported uncertainties are then concentration-only and are an
underestimate of the absolute rate uncertainty — enable these for publishable
absolute uncertainties.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `propagate_production_uncertainty` | bool | `false` | Add a production-calibration term to the external erosion-rate uncertainty, using `production_rel_sigma`. |
| `production_rel_sigma` | float | `0` | Relative 1-sigma production uncertainty (dimensionless, e.g. `0.05` for 5 %). Used only when `propagate_production_uncertainty: true`. Must be ≥ 0 or the run exits fatally. |
| `propagate_decay_uncertainty` | bool | `false` | Add a decay-constant term to the external erosion-rate uncertainty, using `decay_rel_sigma`. |
| `decay_rel_sigma` | float | `0` | Relative 1-sigma decay uncertainty term (dimensionless). Used only when `propagate_decay_uncertainty: true`. Must be ≥ 0 or the run exits fatally. |

---

## Shielding

Shielding values are applied uniformly across the entire basin.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `self_shield_eff_thickness` | float | `0` | Sample self-shielding as effective depth (g/cm²). `0` = surface sample. |
| `snow_shield_eff_thickness` | float | `0` | Snow shielding as effective depth (g/cm²). `0` = no snow correction. |
| `topo_shield_value` | float | `1` | Topographic shielding factor (dimensionless, 0–1). `1` = no topographic shielding — the default, deliberately; see below. |

### Why topographic shielding defaults to off

lsdtt3 applies **no topographic shielding correction by default**, which is a
deliberate departure from CAIRN and from most older basin-averaged workflows.

The reasoning follows DiBiase (2018). The conventional correction accounts only
for the reduction of production *at the surface*, while assuming the effective
mass attenuation length is spatially uniform. It is not: attenuation length
increases on steep slopes, and the two effects largely offset. Correcting for
one and not the other introduces a **slope-dependent systematic error** — which
is worse than a bias, because slope correlates with erosion rate, so it can
manufacture spurious topography-erosion relationships. In steep landscapes the
surface-only correction reaches ~20 %.

DiBiase's conclusion is that topographic shielding corrections are inappropriate
for catchment-mean erosion rates in most cases, and are needed only for steep
catchments with a **non-uniform distribution of quartz and/or erosion rate**.
That is the case in which to set `topo_shield_value` yourself — along with
reproducing a published rate computed the older way, such as a CAIRN rate.

Two things to know if you do set it: it is applied as a **single uniform value
across the whole basin**, not per pixel; and in lsdtt3's published-dataset
comparison tier the topographic shielding term is decomposed at roughly
**+0.70 %** at basin scale (`docs/testing_and_benchmarking.md`).

> DiBiase, R. A. (2018). Short communication: Increasing vertical attenuation
> length of cosmogenic nuclide production on steep slopes negates topographic
> shielding corrections for catchment erosion rates. *Earth Surface Dynamics*,
> 6, 923–931. <https://esurf.copernicus.org/articles/6/923/2018/>

The same argument does not apply to `self_shield_eff_thickness` or
`snow_shield_eff_thickness`, which are ordinary depth-attenuation terms and
default to `0` simply because the common case is an unshielded surface sample.

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fill_before_routing` | bool | `true` | Fill DEM pits before building D8 flow routing. When `false`, route the input DEM unchanged; use this for an already filled DEM or when sinks must intentionally be preserved. |
| `min_slope_for_fill` | float | `0.0001` | Minimum imposed slope during pit filling (m/m). Used only when `fill_before_routing: true`. |
| `boundary_conditions` | string | `oooo` | D8 flow routing boundary conditions. Four-character string (N,E,S,W): `o` = open, `b` = forced outlet. |

---

## Channel snapping

Sample points are snapped to the nearest channel pixel meeting the contributing-pixel threshold using an expanding square search.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold_contributing_pixels` | int | `1000` | Minimum contributing pixels for a node to be considered a channel. |
| `sample_snap_max_search_radius_pixels` | int | `8` | Maximum search radius (pixels) when snapping points to a channel. |

---

## Output

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_production_raster` | bool | `false` | `<write_fname>_PROD.tif` — production scaling raster (Stone 2000 or LSDn, per `scaling_scheme`). |
| `write_basin_production_csv` | bool | `false` | `<write_fname>_BASIN_PROD.csv` — per-sample basin-mean production scaling and mean elevation. |
| `write_cronus_v3_input` | bool | **`true`** | `<write_fname>_CRONUS_V3_INPUT.txt` — a ready-to-paste CRONUS v3 batch request. See below. |

The erosion rate CSV `<write_fname>_CN_ERATE.csv` is always written. Columns include:

- sample metadata (`sample_name` when present, `nuclide`),
- snapped outlet location (`latitude`, `longitude`),
- basin comparison metadata (`effective_pressure_hpa`, `centroid_latitude_deg`),
- measured + predicted concentration,
- effective erosion rate (`apparent_eff_e_g_cm2_yr`),
- erosion rate in mm/kyr (`*_mm_kyr`).

For backward compatibility, duplicated legacy `*_mm_yr` aliases are still written and carry the same **mm/kyr** numeric values.

When an uncertainty column is present in the input, these are appended: the
lower and upper bounds (`concentration_minus_*`, `predicted_conc_minus_*`,
`apparent_eff_e_minus_g_cm2_yr`, `apparent_e_minus_mm_kyr` and the matching
`*_plus_*` set), and the summary uncertainties `apparent_eff_e_unc_int_g_cm2_yr`
/ `apparent_e_unc_int_mm_kyr` and `apparent_eff_e_unc_ext_g_cm2_yr` /
`apparent_e_unc_ext_mm_kyr`.

With `self_consistent_age: true`, three further columns appear:
`self_consistent_t_yr`, `self_consistent_n_iter`, `self_consistent_converged`.

### Basin production diagnostic

`write_basin_production_csv` is a **reporting-only diagnostic — it has no effect
on the inversion**. For each snapped sample it averages the production scaling
raster and the DEM elevation over the *exact* upslope pixel set the inversion
uses, so the basin-mean production can be checked against an external oracle
such as OCTOPUS `BEPROD`/`ELEV_AVE`. Columns: `sample_name`,
`snapped_latitude`, `snapped_longitude`, `n_upslope_pixels`, `mean_prod_scaling`,
`mean_elev_m`.

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

For **sample-set** comparison against CRONUS v3 (separate from synthetic point
scan), use:

```
python3 tests/scripts/compare_cosmo_tool_vs_cronus_v3.py \
  --cosmo-csv <write_path>/<write_fname>_CN_ERATE.csv \
  --points-csv <input_points.csv> \
  --uncertainty-col <sigma_column> \
  --standardisation-col <standardisation_column> \
  --scheme St \
  --out-csv <write_path>/<write_fname>_CRONUS_COMPARE.csv
```

This builds a CRONUS v3 batch request from each sample's nuclide,
concentration, basin-effective pressure, and basin-centroid latitude, then
writes per-sample lsdtt3-vs-CRONUS erosion-rate residuals.

---

## Building a CRONUS v3-compatible input file

**The CLI writes this file for you — no Python required.** `write_cronus_v3_input`
defaults to `true`, so every run produces
`<write_fname>_CRONUS_V3_INPUT.txt` alongside the erosion-rate CSV. Set it to
`false` to suppress it.

A basin-averaged sample has no single site elevation, so the correct way to
present it to a point calculator is as the basin-effective **pressure** together
with the basin **centroid latitude** — exactly what CRONUS v3's `pre`
(pre-calculated pressure) input mode takes. The tool has both numbers already;
they are also reported in the output CSV as `effective_pressure_hpa` and
`centroid_latitude_deg`.

The file contains two lines per sample:

```
S00001 0.000027 -7.488704 1002.3418 pre 0 2.65 1 0.0001 2010 ;
S00001 Be-10 quartz 100000.000000 3000.000000 07KNSTD ;
```

The first line is the sample/site record, the second the nuclide record:

| line | fields |
|---|---|
| site | `id`, latitude, longitude, **pressure (hPa)**, `pre` (flag: the preceding value is a pressure, not an elevation), sample thickness (`0`), density (`2.65`), **shielding (`1`)**, erosion-rate guess (`0.0001`), sample year (`2010`) |
| nuclide | `id`, nuclide (`Be-10` / `Al-26` / `C-14`), target mineral (`quartz`), concentration (atoms/g), 1-sigma uncertainty (atoms/g), AMS standard (`07KNSTD`, `KNSTD`, …) |

Details worth knowing:

- **Shielding is written from `topo_shield_value`**, so it is `1` (none) by
  default, matching the DiBiase (2018) reasoning above. Change
  `topo_shield_value` and the file follows; the run prints a note when it is
  not 1.
- **Density is written from `rock_density_kg_m3`**, converted to g/cm³.
- **Concentrations are written as measured, not standardisation-normalised**,
  and are paired with the sample's standardisation name so the calculator
  applies its own factor. lsdtt3 normalises internally, so writing the
  normalised value *and* the standard name would apply the correction twice —
  the error that once produced a spurious sign flip in the OCTOPUS comparison
  tier.
- **Sample names are sanitised**, since CRONUS ids must be whitespace-free and
  `;` terminates a record: anything outside `[A-Za-z0-9_-]` becomes `_`.
  Duplicate names gain a `_2`, `_3`, … suffix, because the calculator keys its
  results by id and would otherwise merge them.
- If the input has no uncertainty column, **sigma is written as `0`** and the
  run says so. It is never replaced by a guessed percentage — a fabricated
  uncertainty in a file you might publish is worse than an obvious zero.
- Samples whose pressure or centroid latitude did not resolve are omitted with a
  warning rather than written as `nan`, which the calculator would reject.

### Sending it to the online calculator

Paste the file's contents into the text box on the CRONUS v3 erosion-rate input
page (`hess.ess.washington.edu`, with a mirror at `stoneage.hzdr.de`) and
submit. The block is already in the calculator's native batch format, so no
reformatting is needed.

The online calculator applies its own scaling and muon treatment, so a matching
input file is necessary but not sufficient for agreement — see the CRONUS v3
example parameter file below for the lsdtt3 settings that target it.

### Automating the comparison

`tests/scripts/compare_cosmo_tool_vs_cronus_v3.py` submits a batch and writes
per-sample residuals, for when you want the round trip scripted rather than
manual:

```
python3 tests/scripts/compare_cosmo_tool_vs_cronus_v3.py \
  --cosmo-csv   <write_path>/<write_fname>_CN_ERATE.csv \
  --points-csv  <input_points.csv> \
  --uncertainty-col     <sigma_column> \
  --standardisation-col <standardisation_column> \
  --scheme St \
  --out-csv <write_path>/<write_fname>_CRONUS_COMPARE.csv
```

Match `--scheme` to `scaling_scheme` (`St` for `stone`, `LSDn` for `lsdn`). Add
`--no-submit --text-block-out <file>` to build a request without contacting the
calculator, though for that the CLI's own output file is usually what you want.

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
sample_snap_max_search_radius_pixels: 8

sample_points_fname:             CRN_points.csv
concentration_column_name:       Be10_CONC
uncertainty_column_name:         Be10_uncertainty
standardisation_column_name:     standardisation
sample_column_name:              sample_name

write_production_raster: false
```

## Example parameter file: agreement with the CRONUS v3 online calculator

Use this pairing when the goal is to reproduce CRONUS v3. Keep
`sample_age_yr: 0`, pair `cronusEmulatorSt` with `scaling_scheme: stone`, and
match `--scheme St` on `tests/scripts/compare_cosmo_tool_vs_cronus_v3.py`.

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_samples

nuclide:        Be10
muon_scheme:    cronusEmulatorSt
scaling_scheme: stone
atmosphere:     ERA40
sample_age_yr:  0

sample_points_fname:         CRN_points.csv
uncertainty_column_name:     Be10_uncertainty
standardisation_column_name: standardisation

propagate_production_uncertainty: true
production_rel_sigma:             0.05
```

For the LSDn route use `muon_scheme: cronusEmulatorLSDn` with
`scaling_scheme: lsdn` **and** `self_consistent_age: true` — see the
`self_consistent_age` row above for why the fixed-age variant will not
reproduce v3.

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

### CRN physics

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `nuclide` | `string` | `Be10` | `param` | - | Nuclide to use for erosion rate inversion. | Allowed values: Be10, Al26, C14. |
| `muon_scheme` | `string` | `BraucherBorchers` | `param` | - | Muon production parameter scheme. | Allowed values: BraucherBorchers, Braucher, Granger, Schaller, newCRONUS, cronusEmulatorSt, cronusEmulatorLSDn. The cronusEmulator schemes fit lsdtt3's three muon pathways to the CRONUS v3 tabulated model-1A muon inventory; use them when the goal is agreement with the CRONUS v3 online calculator. They are pressure-dependent (a pressure raster is computed automatically) and nuclide-specific. Match cronusEmulatorSt with scaling_scheme: stone and cronusEmulatorLSDn with scaling_scheme: lsdn. |
| `scaling_scheme` | `string` | `stone` | `param` | - | Production-rate scaling scheme. | stone = Stone (2000)/Lal, the lsdtt3 default and the CAIRN/OCTOPUS scheme. lsdn = Lifton, Sato &amp; Dunai (2014) nuclide-specific scaling, interpolated from the LSDn2015 scaling grids; requires the LSDn grids in crn_parameters.nc and supports Be10, Al26 and C14 only. Both schemes are evaluated per pixel, including when sample_age_yr &gt; 0. muon_scheme is unaffected: LSDn is a spallation scaling and the muon pathways receive it through the same Vermeesch redistribution as Stone. |
| `rock_density_kg_m3` | `float` | `2650` | `param` | - | Rock (bulk) density in kg/m3. | Used to convert effective erosion rate (g/cm2/yr) to mm/kyr and mm/yr. |

### Atmosphere and palaeomagnetic correction

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `crn_nc_fname` | `string` | - | `param` | - | Path to crn_parameters.nc (ERA40 grids + Muscheler VDM record). | For the ERA40 atmosphere and paleomagnetic correction. Leave empty to use the crn_parameters.nc bundled with lsdtt3 (found automatically); if that is not present either, the standard (SLHL) atmosphere is used with no correction. |
| `sample_age_yr` | `float` | `0` | `param` | - | Sample exposure age in years (yr BP). Drives the time-integrated production correction; what it selects depends on scaling_scheme. | Set &gt; 0 to apply a time-integrated correction; requires crn_nc_fname. Default 0 = no correction (modern field). Under scaling_scheme=stone this is the Muscheler VDM Stone/Lal paleomagnetic correction. Under scaling_scheme=lsdn it instead drives the duration-weighted average of LSDn production over the Pavon-Carrasco geomagnetic record (Holocene moving paleopole, then a GAD dipole record to ~20 Ma), so it is longitude-dependent as well as latitude-dependent. NOTE: "now" is assumed to be 2010 (the reference year of the packed record); lsdtt3 has no sample-collection-year input yet, so no (collection year - 2010) shift is applied. That offset is at most ~16 yr against records of 1e2-1e7 yr and is therefore negligible, but it IS an assumption, not an exact treatment. |
| `self_consistent_age` | `bool` | `false` | `param` | - | Iterate the sample age to the erosional integration timescale t = Gamma_sp/epsilon instead of holding it at sample_age_yr. This is the CRONUScalc / Balco v3 convention, and it is what published LSDn rates are computed with. | The nuclide inventory in an eroding column is accumulated over roughly one attenuation length of erosion, tau = Gamma_sp/epsilon (Gamma_sp = 160 g/cm2 in every muon scheme lsdtt3 ships), so the geomagnetic field the sample actually integrated is the field over that window, not the modern field. But tau depends on epsilon, which is the unknown, so the two are solved together: invert at sample_age_yr, form t = Gamma_sp/epsilon, rebuild production, re-invert, repeat until epsilon stops moving. Off by default, because it changes the meaning of the reported rate and because it costs one production-raster rebuild per iteration PER SAMPLE (the age, hence the whole raster, differs sample to sample). Typically converges in 2-4 iterations. With scaling_scheme=lsdn on a large DEM this is the dominant cost of a run. Matters most where the field departed furthest from modern, i.e. at low latitude and for slowly eroding samples with long tau; at 20 ka the effect on the inverted rate runs to ~10 %. sample_age_yr is then only the starting guess (0 = start from the modern field). |
| `self_consistent_age_tolerance` | `float` | `0.001` | `param` | - | Relative convergence tolerance on the erosion rate for self_consistent_age. | Iteration stops when \|eps_new - eps_old\|/eps_new falls below this. The default 1e-3 is well inside the inversion solver's own convergence tolerance (~1e-3 round-trip), so tightening it buys precision the forward model cannot deliver. Ignored unless self_consistent_age: true. |
| `self_consistent_age_max_iter` | `int` | `10` | `param` | - | Maximum self_consistent_age iterations before giving up on a sample. | A sample that has not converged by this many iterations is reported with a warning naming the sample and the last two erosion rates, and its LAST iterate is written out -- it is never silently accepted as converged. Ignored unless self_consistent_age: true. |
| `atmosphere` | `string` | `ERA40` | `param` | - | Atmosphere model: ERA40 (default), NCEP, or standard. | ERA40 = Balco/Lifton ERA40 climatology + latitude-dependent lapse (lsdtt3 default). NCEP = Balco NCEPatm_2 (the atmosphere CAIRN/OCTOPUS use; requires NCEP grids in the crn_parameters.nc). standard = SLHL barometric formula, no climatology. ERA40/NCEP fall back to standard if the required grids are absent. |

### Sample points

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `sample_points_fname` | `string` | `CRN_points.csv` | `param` | - | CSV file containing sample point locations and concentrations. | Must have latitude and longitude columns. |
| `sample_column_name` | `string` | `sample_name` | `param` | - | Name of the sample identifier column in the points CSV. | Optional; if absent the column is omitted from output. |
| `concentration_column_name` | `string` | - | `param` | - | Column name for the measured nuclide concentration (atoms/g). | Defaults to &lt;nuclide&gt;_CONC (e.g. Be10_CONC) if not set. |
| `uncertainty_column_name` | `string` | - | `param` | - | Column name for concentration measurement uncertainty (1-sigma, atoms/g). | Defaults to &lt;nuclide&gt;_uncertainty. If absent, uncertainty bounds are skipped. |
| `propagate_production_uncertainty` | `bool` | `false` | `switch` | - | If true, add a production-calibration term to external erosion-rate uncertainty. | Added in quadrature with concentration-derived uncertainty using production_rel_sigma. |
| `production_rel_sigma` | `float` | `0` | `param` | - | Relative 1-sigma production uncertainty (dimensionless, e.g. 0.05 for 5%). | Used only when propagate_production_uncertainty is true. |
| `propagate_decay_uncertainty` | `bool` | `false` | `switch` | - | If true, add a decay-constant term to external erosion-rate uncertainty. | Added in quadrature with concentration-derived uncertainty using decay_rel_sigma. |
| `decay_rel_sigma` | `float` | `0` | `param` | - | Relative 1-sigma decay uncertainty term (dimensionless). | Used only when propagate_decay_uncertainty is true. |
| `standardisation_column_name` | `string` | `standardisation` | `param` | - | Column name for the AMS standardisation name. | If absent, no standardisation correction is applied. |

### Shielding

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `self_shield_eff_thickness` | `float` | `0` | `param` | - | Sample self-shielding as effective depth (g/cm2). | 0 = surface sample. |
| `snow_shield_eff_thickness` | `float` | `0` | `param` | - | Snow shielding as effective depth (g/cm2). | 0 = no snow correction. Applied uniformly across the basin. |
| `topo_shield_value` | `float` | `1` | `param` | - | Topographic shielding factor (dimensionless, 0-1). | 1 = no topographic shielding. Applied uniformly across the basin. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_before_routing` | `bool` | `true` | `switch` | - | Choose whether to pit-fill the DEM before building the flow network. | ASKS: should the DEM be filled at all? Filling is OPTIONAL here and this is a genuine analytical choice. True (the default) routes over a pit-free surface. False builds the flow network on the UNFILLED surface, so flow paths terminate at sinks instead of being carried through them. Do not confuse this with raster_is_filled, which is an assertion that the input is already filled. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_slope_for_fill` | `float` | `0.0001` | `param` | - | Minimum imposed slope during pit filling (m/m). | Used only when fill_before_routing is true. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | D8 flow routing boundary conditions. | Four-character string (N,E,S,W): o=open, b=forced outlet. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Channel snapping

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `sample_snap_max_search_radius_pixels` | `int` | `8` | `param` | - | Maximum search radius (pixels) when snapping points to channel. | Expanding square search up to this radius. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_production_raster` | `bool` | `false` | `switch` | - | Write the production scaling raster to a GeoTIFF. | Stone 2000 or LSDn, per scaling_scheme. Filename: &lt;write_fname&gt;_PROD.tif |
| `write_basin_production_csv` | `bool` | `false` | `switch` | - | Write a per-sample basin-mean production-scaling + mean-elevation CSV (diagnostic; does not affect the inversion). | Filename: &lt;write_fname&gt;_BASIN_PROD.csv |
| `write_cronus_v3_input` | `bool` | `true` | `switch` | - | Write a CRONUS v3 input file that can be pasted into the online calculator. | Basin-effective pressure + centroid latitude, in the calculator's batch format. Filename: &lt;write_fname&gt;_CRONUS_V3_INPUT.txt |

<!-- END GENERATED CLI PARAMETERS -->
