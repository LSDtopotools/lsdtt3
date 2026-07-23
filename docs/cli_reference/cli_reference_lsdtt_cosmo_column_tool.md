# CLI reference: lsdtt-cosmo-column-tool

Simulates cosmogenic nuclide concentrations (¹⁰Be, ²⁶Al, ¹⁴C) in a bedrock or sediment column undergoing a known or hypothesised erosion history. No DEM or raster is required — this is a purely 1-D particle calculation.

Three simulation modes are available, selected by `column_mode`.

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

---

## Mode and physics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `column_mode` | string | `steady_state` | Simulation mode. Options: `steady_state`, `step_change`, `full_transient`. |
| `nuclide` | string | `Be10` | Nuclide for output labelling. All three nuclides are always computed. Options: `Be10`, `Al26`, `C14`. |
| `muon_scheme` | string | `BraucherBorchers` | Muon production parameter scheme. Options: `BraucherBorchers`, `Braucher`, `Granger`, `Schaller`, `newCRONUS`. |
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

## Parameters for `full_transient` mode

| Parameter | Type | Default | Description |
|---|---|---|---|
| `transient_history_file` | string | `erosion_history.csv` | CSV file with erosion history. Required columns: `time` (yr BP, ascending), `rate` (m/yr), `removal` (m, used when `rate` = 0). Last row: `time = -9999`, `rate` = background rate. |
| `n_samples` | int | `10` | Number of evenly-spaced depth samples in the column. |
| `sample_spacing_m` | float | `0.05` | Spacing between depth samples (m). First sample is at the surface (0 m). |
| `use_initial_concentration` | bool | `false` | If true, particles start with a user-specified concentration instead of steady state. Useful for re-exposed surfaces with a known inherited inventory. |
| `initial_10Be_conc_atoms_g` | float | `0` | Initial ¹⁰Be concentration throughout the column (atoms/g). Used only when `use_initial_concentration: true`. |
| `initial_26Al_conc_atoms_g` | float | `0` | Initial ²⁶Al concentration (atoms/g). |
| `initial_14C_conc_atoms_g` | float | `0` | Initial ¹⁴C concentration (atoms/g). |

---

## Output

For `steady_state` and `step_change` modes, results are printed to standard output only.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `output_file` | string | `crn_column.csv` | Output CSV filename for column depth profile (`full_transient` mode only). Written to `write_path`. Columns: `depth_m`, `effective_depth_g_cm2`, `C_10Be`, `C_26Al`, `C_14C`. |

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

transient_history_file:  erosion_history.csv
n_samples:               20
sample_spacing_m:        0.05

use_initial_concentration: false

output_file: transient_column_output.csv
```

### Erosion history CSV format

```
time,rate,removal
100000,0.0002,0
50000,0.001,0
10000,0.0005,0
0,0.0005,0
-9999,0.0005,0
```
