# CLI reference: lsdtt-cosmo-column-tool

Simulates cosmogenic nuclide concentrations (¹⁰Be, ²⁶Al, ¹⁴C) in a bedrock or sediment column undergoing a known or hypothesised erosion history. No DEM or raster is required — this is a purely 1-D particle calculation.

Seven simulation modes are available, selected by `column_mode`.

## Usage

```
lsdtt-cosmo-column-tool <path> <parameter_file>
lsdtt-cosmo-column-tool <parameter_file>         # path defaults to ./
lsdtt-cosmo-column-tool -h                       # print help and exit
lsdtt-cosmo-column-tool -v                       # print version and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Simulation modes

| Mode | Description |
|---|---|
| `steady_state` | Single particle at a given effective depth under a constant erosion rate. |
| `step_change` | Single particle that was at one steady-state erosion rate for a long time, then changed to a new rate. |
| `full_transient` | Column of N evenly-spaced particles driven by a user-supplied erosion-history CSV. Writes a depth-profile CSV. |
| `burial_pair` | Two-nuclide burial/inheritance inversion for a single sample (fast/slow pair). |
| `step_pair` | Two-nuclide step-change inversion for a single sample (fast/slow pair). |
| `burial_pair_batch` | Batch burial/inheritance inversion from a CSV (one output row per input row). |
| `step_pair_batch` | Batch step-change inversion from a CSV (one output row per input row). |

---

## Mode and physics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `column_mode` | string | `steady_state` | Simulation mode. Options: `steady_state`, `step_change`, `full_transient`, `burial_pair`, `step_pair`, `burial_pair_batch`, `step_pair_batch`. |
| `nuclide` | string | `Be10` | Nuclide for output labelling. All three nuclides are always computed. Options: `Be10`, `Al26`, `C14`. |
| `muon_scheme` | string | `BraucherBorchers` | Muon production parameter scheme. Options: `BraucherBorchers`, `Braucher`, `Granger`, `Schaller`, `newCRONUS`. The pressure-dependent `cronusEmulatorSt` / `cronusEmulatorLSDn` schemes available in `lsdtt-cosmo-tool` are **not usable here** — see the note below. |

> **The `cronusEmulator*` schemes do not work in this tool.** Their muon terms
> depend on both site pressure and nuclide, and this tool works on a single
> column with no pressure raster, so there is no correct parameter set to
> install. The CLI's `muon_scheme` check validates against the shared
> whole-of-lsdtt3 scheme list, which *does* include them, so the name is
> accepted at parse time and the run then **aborts** (uncaught
> `std::invalid_argument` from `reset_scaling`, exit 134) rather than exiting
> with a clean error message. Use `lsdtt-cosmo-tool` if you need the emulator
> schemes. Note that `muon_scheme` is also a per-row column in both batch pair
> CSVs, and the same restriction applies there.
| `rock_density_kg_m3` | float | `2650` | Rock density (kg/m³). Used to convert depth to effective depth. |

---

## Production scaling

Supply either `production_scaling` directly, or let the tool compute it from site location. If `latitude` and `elevation` are both set (i.e. > -9998), Stone 2000 scaling is computed from those values and overrides `production_scaling`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `production_scaling` | float | `1.0` | Stone 2000 production scaling factor (dimensionless). Overridden when `latitude` and `elevation` are supplied. |
| `latitude` | float | `-9999` | Site latitude (degrees). |
| `longitude` | float | `-9999` | Site longitude (degrees). |
| `elevation` | float | `-9999` | Site elevation (m). |
| `total_shielding` | float | `1.0` | Total shielding factor (topographic × snow × self, dimensionless 0–1). Applied on top of `production_scaling`. `1` = fully exposed. |

---

## Parameters for `steady_state` and `step_change` modes

| Parameter | Type | Default | Description |
|---|---|---|---|
| `effective_depth_g_cm2` | float | `0` | Effective depth of sample (g/cm²). `0` = surface. |
| `effective_erosion_rate_g_cm2_yr` | float | `0.1` | Effective erosion rate (g/cm²/yr). Used as current rate in `steady_state` and as the new rate in `step_change`. |

### Additional parameters for `step_change` mode

| Parameter | Type | Default | Description |
|---|---|---|---|
| `effective_erosion_rate_old_g_cm2_yr` | float | `0.1` | Old (pre-step) effective erosion rate (g/cm²/yr). |
| `time_since_step_yr` | float | `10000` | Time since the erosion-rate step change (yr). |

---

## Parameters for pair inversions (`burial_pair` and `step_pair`)

These modes work on a **fast/slow nuclide pair**, so users can choose any two of `Be10`, `Al26`, `C14` (with `nuclide_fast` the shorter-lived member).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `nuclide_fast` | string | `C14` | Fast (shorter-lived) nuclide. Typical: `C14`/`Be10` or `Al26`/`Be10`. |
| `nuclide_slow` | string | `Be10` | Slow (longer-lived) nuclide. |
| `conc_fast_atoms_g` | float | `0` | Measured concentration of `nuclide_fast` (atoms/g). |
| `conc_slow_atoms_g` | float | `0` | Measured concentration of `nuclide_slow` (atoms/g). |
| `conc_fast_sigma_atoms_g` | float | `0` | 1-sigma concentration uncertainty of `nuclide_fast` (atoms/g). `0` disables concentration-corner uncertainty propagation. |
| `conc_slow_sigma_atoms_g` | float | `0` | 1-sigma concentration uncertainty of `nuclide_slow` (atoms/g). `0` disables concentration-corner uncertainty propagation. |

Additional parameter for `step_pair`:

| Parameter | Type | Default | Description |
|---|---|---|---|
| `effective_erosion_rate_old_g_cm2_yr` | float | `0.1` | Assumed pre-step erosion rate (g/cm²/yr). |

Legacy compatibility for `burial_pair`: `conc_26Al_atoms_g` and `conc_10Be_atoms_g` are still accepted as fallback inputs. Legacy fallback sigmas are also supported: `conc_26Al_sigma_atoms_g`, `conc_10Be_sigma_atoms_g`.

---

## Pair uncertainty controls

These controls apply to `burial_pair`, `step_pair`, and both pair batch modes.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `propagate_production_uncertainty` | bool | `false` | If `true`, adds a production term to external erosion uncertainty. |
| `production_rel_sigma` | float | `0` | Relative 1-sigma production uncertainty (e.g. `0.05` = 5%). Used only when `propagate_production_uncertainty: true`. |
| `propagate_decay_uncertainty` | bool | `false` | If `true`, adds a decay term to external erosion uncertainty. |
| `decay_rel_sigma` | float | `0` | Relative 1-sigma decay uncertainty term. Used only when `propagate_decay_uncertainty: true`. |

Uncertainty method for pair inversions:
- concentration-corner re-solves at `fast ± sigma_fast` and `slow ± sigma_slow`
- `unc_int` from half-range of converged corner solutions
- `unc_ext` from quadrature of `unc_int` plus enabled production/decay terms

---

## Parameters for `full_transient` mode

| Parameter | Type | Default | Description |
|---|---|---|---|
| `transient_history_fname` | string | `erosion_history.csv` | CSV file with erosion history. This history-first column mode defaults to a conventional filename; raster prediction instead uses `NULL` because its transient mode is optional. Required columns: `time` (yr BP), `rate` (m/yr), `removal` (m, used when `rate` = 0). **Row order matters - see below.** Last row: `time = -9999`, `rate` = background rate. |

> ### ⚠ Row order: the FIRST data row is the MOST RECENT
>
> `time` must **increase with row index**, so row 0 is the youngest (smallest
> yr BP) and the record runs backwards in time as you read down. A reversed
> record is rejected fatally:
>
> ```
> FATAL: lsdtt_crn_column: erosion history must be ordered with index 0 the
> MOST RECENT row and yr BP increasing with index...
> ```
>
> This is a hard failure by design. Before the check was added, a reversed
> record fed negative timesteps straight into the integrator and returned a
> confidently wrong number with no error, no warning and no NaN — measured
> 72,357 forward against 178,467 reversed.
>
> Each row's `time` is the yr BP at which that step *began*, so `time` on row 0
> is the duration of the most recent step. It must be greater than 0, or the
> most recent step has zero length.
>
> Two rows may share a `time` **only** when the second is a removal event
> (`rate = 0`, `removal > 0`), which is instantaneous.
| `n_samples` | int | `10` | Number of evenly-spaced depth samples in the column. |
| `sample_spacing_metres` | float | `0.05` | Spacing between depth samples (m). First sample is at the surface (0 m). |
| `use_initial_concentration` | bool | `false` | If true, particles start with a user-specified concentration instead of steady state. Useful for re-exposed surfaces with a known inherited inventory. |
| `initial_10Be_conc_atoms_g` | float | `0` | Initial ¹⁰Be concentration throughout the column (atoms/g). Used only when `use_initial_concentration: true`. |
| `initial_26Al_conc_atoms_g` | float | `0` | Initial ²⁶Al concentration (atoms/g). |
| `initial_14C_conc_atoms_g` | float | `0` | Initial ¹⁴C concentration (atoms/g). |

---

## Parameters for batch pair inversions

| Parameter | Type | Default | Description |
|---|---|---|---|
| `burial_pairs_fname` | string | `burial_pairs.csv` | Input CSV for `burial_pair_batch`. |
| `step_pairs_fname` | string | `step_pairs.csv` | Input CSV for `step_pair_batch`. |

`burial_pair_batch` required CSV columns:
- `conc_fast_atoms_g`, `conc_slow_atoms_g`

`burial_pair_batch` optional per-row columns:
- `sample_name`, `nuclide_fast`, `nuclide_slow`, `total_shielding`, `muon_scheme`, `conc_fast_sigma_atoms_g`, `conc_slow_sigma_atoms_g`

Legacy `burial_pair_batch` fallback columns (still accepted):
- `conc_26Al_atoms_g`, `conc_10Be_atoms_g`
- legacy fallback sigmas: `conc_26Al_sigma_atoms_g`, `conc_10Be_sigma_atoms_g`

`step_pair_batch` required CSV columns:
- `conc_fast_atoms_g`, `conc_slow_atoms_g`

`step_pair_batch` optional per-row columns:
- `sample_name`, `nuclide_fast`, `nuclide_slow`, `effective_erosion_rate_old_g_cm2_yr`, `total_shielding`, `muon_scheme`, `conc_fast_sigma_atoms_g`, `conc_slow_sigma_atoms_g`

---

## Output

All modes write CSV output to `output_file` in `write_path`.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `output_file` | string | `crn_column.csv` | Output CSV filename. `steady_state`/`step_change`/`full_transient` use columns: `depth_m`, `effective_depth_g_cm2`, `C_10Be`, `C_26Al`, `C_14C`; pair inversion modes write diagnostic rows including `status`/`status_code`, uncertainty fields (`*_unc_int_*`, `*_unc_ext_*`), and corner bounds (`*_lo_*`, `*_hi_*`); batch modes write one diagnostic row per input row. |
| `write_cronus_v3_input` | bool | **`true`** | Write `<output_file stem>_CRONUS_V3_INPUT.txt`, a ready-to-paste CRONUS v3 batch request. See below. |

### CRONUS v3 input file

The tool writes a CRONUS v3 input file by default, so a column result can be
checked against the online calculator without any scripting. Set
`write_cronus_v3_input: false` to suppress it. The format matches
`lsdtt-cosmo-tool`'s, and the file is named from the `output_file` stem — e.g.
`output_file: burial_pair_c14_be10.csv` gives
`burial_pair_c14_be10_CRONUS_V3_INPUT.txt`.

**It requires `latitude` and `elevation`.** A CRONUS site record is a real
location, so a run configured with `production_scaling` alone cannot produce
one; the tool says so once and carries on. Unlike `lsdtt-cosmo-tool` (which
supplies a basin-effective pressure with the `pre` flag), this tool writes the
elevation with the `std` flag and lets the calculator derive pressure from the
standard atmosphere.

What each mode contributes:

| Mode | Written? | Contents |
|---|---|---|
| `burial_pair`, `step_pair` | yes | One site, two nuclide records — the **measured** concentrations and sigmas. This is the natural CRONUS input. |
| `burial_pair_batch`, `step_pair_batch` | yes | One site per input row, named from the row's `sample_name`, honouring any per-row `total_shielding`. |
| `steady_state`, `step_change` | only at the surface | The **modelled** concentration for the selected `nuclide`, with sigma `0`. A round-trip check: lsdtt3 forward-models, CRONUS inverts. Skipped with a note when `effective_depth_g_cm2` is not `0`, because CRONUS would treat a shielded sample as exposed. |
| `full_transient` | no | A depth profile is not a point sample, and a transient history cannot be expressed in a CRONUS erosion-rate record. |

Other behaviour: shielding is written from `total_shielding` (note that in this
tool it is the *combined* topographic × snow × self factor, whereas the CRONUS
field is nominally topographic); density from `rock_density_kg_m3`; names are
sanitised to `[A-Za-z0-9_-]` and duplicates get a `_2`, `_3`, … suffix; a
concentration with no uncertainty is written as sigma `0` and reported, never
replaced by a guess.

> ### ⚠ AMS standardisation
>
> **This tool applies no AMS standardisation correction.** It has no
> standardisation parameter: concentrations are used exactly as given, which is
> the right behaviour for a forward model but means the *basis* of a measured
> concentration is yours to track.
>
> In the CRONUS v3 input file, each concentration is therefore written with the
> calculator's **default** standard for its nuclide — `07KNSTD` for Be-10,
> `KNSTD` for Al-26, none for C-14. If your measurement is on a different
> standard, edit that field in the generated file before submitting, or CRONUS
> will interpret it on the wrong basis.
>
> Conversely, do **not** normalise a concentration yourself and then leave the
> default standard in place — that applies the factor twice. The safe
> combination is an as-measured concentration with its own standard named in the
> file. `lsdtt-cosmo-tool` automates this via `standardisation_column_name`;
> this tool does not, because a single column has no sample metadata to read it
> from.

---

## lsdviztools3 interface notes

`lsdtt-cosmo-column-tool` is the lightest route for direct point-style CRN
experiments (no DEM). Those CSVs can be reshaped into comparison tables for
plotting with `lsdviztools3` chart helpers.

- For lsdtt3 standard comparison figure workflows, the preferred ready-made
  tables are generated by:
  - `tests/scripts/run_cosmo_single_point_scaling_scan.py` (sweeps)
  - `tests/scripts/run_cosmo_validation.py` /
    `tests/scripts/run_cosmo_riversand_comparison.py` (sample summaries)
- lsdviztools3 plotting entry points:
  - `lsdviz3-cosmo-compare` (sweep CSVs)
  - `lsdviz3-cosmo-published` (published summary CSVs)

---

## Example parameter file: steady state

```
column_mode:  steady_state
nuclide:      Be10
muon_scheme:  BraucherBorchers
rock_density_kg_m3: 2650

latitude:   45.0
elevation:  1500.0
total_shielding: 0.95

effective_depth_g_cm2:            0.0
effective_erosion_rate_g_cm2_yr:  0.05
```

## Example parameter file: full transient

```
column_mode:  full_transient
muon_scheme:  BraucherBorchers
rock_density_kg_m3: 2650

latitude:   45.0
elevation:  1500.0
total_shielding: 1.0

transient_history_fname: erosion_history.csv
n_samples:               20
sample_spacing_metres:   0.05

use_initial_concentration: false

output_file: transient_column_output.csv
```

## Example parameter file: single pair inversion with C14/Be10

```
column_mode:  burial_pair
muon_scheme:  BraucherBorchers
rock_density_kg_m3: 2650

nuclide_fast: C14
nuclide_slow: Be10
conc_fast_atoms_g: 62000
conc_slow_atoms_g: 410000
conc_fast_sigma_atoms_g: 1200
conc_slow_sigma_atoms_g: 8200
propagate_production_uncertainty: true
production_rel_sigma: 0.05
total_shielding: 1.0

output_file: burial_pair_c14_be10.csv
```

## Example parameter file: batch pair inversion

```
column_mode: step_pair_batch
step_pairs_fname: step_pairs.csv
output_file: step_pair_batch_results.csv
```

### Step-pair batch CSV format

```
sample_name,nuclide_fast,nuclide_slow,conc_fast_atoms_g,conc_slow_atoms_g,conc_fast_sigma_atoms_g,conc_slow_sigma_atoms_g,effective_erosion_rate_old_g_cm2_yr,total_shielding,muon_scheme
s1,C14,Be10,74000,355000,1500,7000,0.002,1.0,BraucherBorchers
s2,Al26,Be10,2100000,300000,42000,6000,0.002,1.0,Braucher
```

### Burial-pair batch CSV format

```
sample_name,nuclide_fast,nuclide_slow,conc_fast_atoms_g,conc_slow_atoms_g,conc_fast_sigma_atoms_g,conc_slow_sigma_atoms_g,total_shielding,muon_scheme
b1,C14,Be10,62000,410000,1200,8200,1.0,BraucherBorchers
b2,Al26,Be10,2300000,330000,46000,6600,0.8,Braucher
```

### Erosion history CSV format

Row 0 is the **most recent** step; `time` increases downwards. Both examples
below are verified to run.

```
time,rate,removal
10000,0.0005,0
50000,0.001,0
100000,0.0002,0
-9999,0.0002,0
```

Reading down: the last 10,000 yr eroded at 0.0005 m/yr; from 10,000 to 50,000 yr
BP at 0.001; from 50,000 to 100,000 yr BP at 0.0002; and 0.0002 m/yr as the
background before that.

**A landslide** is an instantaneous removal — `rate = 0` with `removal` in
metres, sharing the `time` of the row above it. Here 2 m is stripped at 5,000 yr
BP, with 0.0005 m/yr background either side:

```
time,rate,removal
5000,0.0005,0
5000,0,2.0
100000,0.0005,0
-9999,0.0005,0
```
