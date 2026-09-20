# CLI reference: lsdtt-distributed-hydrology

Runs a time-evolving distributed hydrological model combining TOPMODEL saturation-excess routing with optional interception, infiltration excess, and channel routing. Simulates the full water balance from rainfall through groundwater and surface routing to outlet discharge.

This CLI is **hydrological** (generates discharge from rainfall), not hydraulic. For inundation mapping or hydraulic routing of a known discharge, see [`lsdtt-graphflood`](cli_reference_lsdtt_graphflood.md). For static drainage-related metrics (TWI, STI, connectivity), see [`lsdtt-drainage-analysis`](cli_reference_lsdtt_drainage_analysis.md).

See [Hydrological modelling in lsdtt3](../design_and_algorithms/hydrological_modelling.md) for a three-tool comparison and guidance on choosing between them.

## Usage

```
lsdtt-distributed-hydrology <path> <parameter_file>
lsdtt-distributed-hydrology <parameter_file>         # path defaults to ./
lsdtt-distributed-hydrology -h                       # print help and exit
lsdtt-distributed-hydrology -v                       # print version and exit
lsdtt-distributed-hydrology -c                       # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Overview

The model runs on a regular raster grid and follows a column-based approach:

1. **Rainfall** falls on the canopy (if interception is enabled) and throughfall reaches the soil.
2. **Infiltration** is partitioned into root-zone storage, unsaturated-zone storage, and infiltration excess (if enabled).
3. **Saturation** is computed pixel-by-pixel either via TOPMODEL (classic index mode) or with explicit subsurface routing (dynamic mode).
4. **Surface routing** conveys water through the flow network to the outlet, with travel times set by velocity and a width-function or dynamic channel routing.

The model combines **store-based hydrology** (deficits, depths) with **topographic routing** (drainage-area-based velocity, flow graph routing). Water is routed at every timestep; outputs include hydrographs, water-balance closures, and optional pixel-by-pixel diagnostics.

**Experimental status:** This CLI is version 0.1.0 and has not yet been published in peer-reviewed literature. Use with caution and report issues via GitHub.

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename (without extension). |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Output raster format. Options: `cog`, `gtiff`, `bil`, `png`, `same_as_input`. |
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | If true, skip pit-filling and use the input raster directly. Set to true only if the DEM has been pre-filled. |
| `fill_min_slope` | float | `0.0001` | Minimum gradient (m/m) imposed during pit-filling. Used only when `raster_is_filled: false`. |
| `boundary_conditions` | string | `oooo` | Boundary conditions for flow routing. Four-character string (N,E,S,W): `o` = open (flows may exit), `n` = no-flux (closed), `p` = periodic. |

---

## Flow routing and slopes

| Parameter | Type | Default | Description |
|---|---|---|---|
| `flow_method` | string | `dinf` | Flow routing method used for drainage area and all graph sweeps. **Allowed values: `dinf` (Tarboton 1997) or `mfd` (Quinn et al. 1991).** D8 is **not** available here: it sends all of a cell's area to one neighbour, which creates artificial stripes in the wetness index and all deficits derived from it. |
| `mfd_exponent` | string | `1.0` | Slope-weighting exponent for MFD flow partitioning. Quinn et al. (1991) use 1.0; Freeman (1991) suggest 1.1. Used only when `flow_method: mfd`. |
| `twi_window_radius_metres` | float | `100.0` | Polyfit window radius (metres) for the slope used in the Topographic Wetness Index and velocity calculations. Must be large enough to capture local slopes but small enough to capture hillslope-scale variation. |

---

## Model structure

| Parameter | Type | Default | Description |
|---|---|---|---|
| `saturated_zone_mode` | string | `index` | How the saturated zone is represented. **`index`**: each cell's saturation deficit is derived from the catchment-mean deficit and the topographic index (classic TOPMODEL; the configuration with an external reference model). **`dynamic`**: each cell carries its own deficit as a state variable and subsurface flux is routed down the flow graph. |
| `routing_method` | string | `width_function` | How runoff is routed to the outlet. **`width_function`** (default): pure translation using the precomputed travel-time distribution (the width function or GIUH; exact for a static velocity field). **`linear_reservoir`** adds a per-cell store, **`muskingum_cunge`** matches numerical diffusion to hydraulic diffusivity, and **`kinematic_wave`** uses an implicit four-point solve. |
| `velocity_mode` | string | `hydraulic_geometry` | How per-cell overland velocity is derived. **`hydraulic_geometry`** (default): velocity from Manning equation with hydraulic radius derived from drainage area (continuous ordering, low-order channels slower). **`roughness_raster`**: velocity from a single supplied `roughness_fname` raster (SD-TOPMODEL mode). **`strahler_table`**: explicit velocity per Strahler order (diagnostic override). **`constant`**: one channel and one hillslope velocity (classic TOPMODEL equivalent). |
| `channel_extent_mode` | string | `none` | How channel cells are identified for roughness or velocity override. **`none`** (default): geometry-only (roughness and velocity apply everywhere; smoothly varying hydraulic radius). **`area_threshold`**: cells with drainage area > `channel_extent_area_threshold_m2` are channels. **`area_slope_threshold`**: cells with `A * S^theta > const` are channels (steeper channel heads at smaller area). **`extraction`**: use a prior channel-extraction result. **`vector`**: snap a mapped river network (GIS vector) to the flow graph. All modes override velocity and roughness on channel cells only; the model runs correctly without any channel mask. |

---

## TOPMODEL core stores

All depths and deficits are in metres. All time-constant parameters must be specified.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `szm_metres` | float | `0.032` | Exponential store depth m (metres) for the saturated zone deficit exponential profile. Sets the decay length for transmissivity with depth. Typical range: 0.01–0.1 m. |
| `ln_t0` | float | `5.0` | Natural logarithm of the saturated transmissivity at zero deficit, T0 (m²/h). Equivalent to `log(T0_m2_per_hour)`. Classic TOPMODEL calibration parameter. Typical range (log scale): 3–8. **Note:** When `transmissivity_source: sti_table` is used, this value is overridden by the spatially distributed field derived from soil data; `ln_t0` is retained for fallback or diagnostics. |
| `td_hours_per_metre` | float | `50.0` | Unsaturated-zone time-delay parameter: time (hours) per unit deficit (metres) for drainage from unsaturated storage. Equation: `q_uz = S_uz / (D_i * td_hours_per_metre)`, capped at available storage. Typical range: 10–200 h/m. |
| `srmax_metres` | float | `0.05` | Maximum root-zone deficit (metres), i.e., plant-available water capacity. Storage is lost to actual ET when the root zone is not at capacity. Typical range: 0.01–0.2 m. |
| `q0_metres_per_step` | float | `3.28e-5` | Initial subsurface flow per unit area (metres per outer timestep) used as a baseline in saturated-zone routing. Typical range: 1e-6–1e-4 m/step. |
| `sr0_metres` | float | `0.002` | Initial root-zone deficit (metres) at the start of the run. Defines the starting wetness. Typical: 0.001–0.01 m. |
| `initial_deficit_metres` | float | `0.3` | Initial saturation deficit (metres) at every pixel in `saturated_zone_mode: dynamic`. Used only in dynamic mode and only at the first timestep. Typical range: 0.1–1.0 m. |
| `saturated_zone_substeps` | int | `1` | Number of sub-steps within one outer timestep for saturated-zone routing in dynamic mode. Set `0` to choose a power-of-two count automatically when the run reports stiff cells. |

---

## Transmissivity and soil data

Transmissivity (saturated conductivity times depth) comes either from a uniform value or a lookup table keyed on soil series and land-use class.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `transmissivity_source` | string | `uniform` | Where the saturated transmissivity T0 comes from. **`uniform`**: use the scalar `ln_t0` everywhere. **`sti_table`**: derive a spatial T0 field by looking up measured soil conductivity profiles keyed on soil series (polygon) and land-use class (raster). |
| `write_ln_t0_raster` | bool | `false` | Write the measured or assigned ln(T0) field as a raster for inspection. Useful for validating `sti_table` lookups. |
| `sti_soil_fname` | string | | Soil series polygon vector (GIS file, e.g., `.shp`, `.fgb`) used to key the transmissivity lookup. Required when `transmissivity_source: sti_table`. |
| `sti_soil_layer_name` | string | | Layer name within `sti_soil_fname` carrying the series code. Required when `transmissivity_source: sti_table`. |
| `sti_soil_field` | string | | Attribute field in `sti_soil_fname` holding the soil series code. Required when `transmissivity_source: sti_table`. |
| `sti_ksat_fname` | string | | CSV of measured saturated conductivity profiles, keyed on soil series and land-use code. Format: one row per profile; columns include series code, land-use code, and Ks values at depth. Required when `transmissivity_source: sti_table`. |
| `sti_depth_cap_centimetres` | float | `100.0` | Maximum profile depth (cm) used in the transmissivity integration. Horizons deeper than this are ignored. |
| `sti_min_valid_depth_centimetres` | float | `30.0` | Minimum profile depth (cm) accepted into the lookup table. Shallow profiles are rejected. |
| `sti_ks_unit_scale_to_metres_per_day` | float | `1.0` | Multiplier converting the CSV conductivity column to metres per day. E.g., if CSV is in cm/day, use `0.01`. |

---

## Interception (optional)

Canopy interception reduces throughfall and may increase actual evapotranspiration (wet-canopy evaporation). Disabled by default.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_interception` | bool | `false` | Enable the canopy interception store. When true, rainfall is intercepted before reaching the root zone; excess drains as throughfall. |
| `canopy_storage_mm_per_lai` | float | `0.2` | Canopy storage capacity (mm) per unit Leaf Area Index. Typical: 0.1–0.5 mm/LAI. |
| `lai_gap_coefficient` | float | `0.6` | Beer-Lambert extinction coefficient for the free-throughfall fraction. Typical: 0.4–0.7. |
| `lai_seasonal_factor` | float | `1.0` | Multiplier applied to every LAI value before building the canopy store. Use to represent seasonal leaf-off or phenological cycles. E.g., `0.5` for winter dormancy. |
| `interception_evap_factor` | float | `1.0` | Wet-canopy evaporation (from intercepted water) as a multiple of the supplied PET. Typical: 1.0–1.5 (wet canopy often evaporates faster than dry surfaces). |
| `landcover_fname` | string | | Land-cover class raster (integer codes). Used to look up LAI per class when `lai_fname` is not supplied. |
| `landcover_table` | string | | Land-cover class definitions. Format: semicolon-delimited class definitions, each with `name:code=N,lai=L,n=R,soil=S`. Example: `grazing:code=1,lai=1.5,n=0.30,soil=CULT;woodland:code=2,lai=5.0,n=0.60,soil=SEMI`. **Note: parameter values must not contain spaces.** |
| `lai_fname` | string | | Leaf Area Index raster (float, dimensionless). Overrides land-cover-based LAI if supplied. |

---

## Land-cover scenarios (optional)

Apply reclassification rules to the land-cover raster before the run, optionally limited to a spatial mask.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `landcover_scenario` | string | | Reclassification rules applied to `landcover_fname` before the run. Format: semicolon-delimited class-name rules, each `source->target`. Example: `grazing->woodland;arable->woodland`. `*` matches any source class; rules are applied in order. |
| `landcover_intervention_mask_fname` | string | | Raster limiting the spatial extent of `landcover_scenario`. Cells with values > 0.5 receive the reclassification; others are unchanged. Useful for NFM scenarios with localized treatments. |

---

## Velocity, manning roughness, and hydraulic geometry

Velocity is computed from the Manning equation with a flow-depth term (hydraulic radius). Roughness can come from a raster, land cover, or defaults.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `roughness_fname` | string | | Manning n raster for hillslope and floodplain cells (overrides land-cover-derived n if supplied). Used only in `velocity_mode: roughness_raster`. |
| `roughness_default` | float | `0.035` | Manning n used on cells with no assigned land cover (default fallback). Typical range: 0.02–0.1 (grass/short vegetation). |
| `hydraulic_radius_coeff` | float | `0.1` | Coefficient a_R in the hydraulic-radius law: `R = max(h_sheet, a_R * (A / A_ref)^c_R)`. Sets the strength of the downstream-deepening trend. |
| `hydraulic_radius_exponent` | float | `0.4` | Exponent c_R in the hydraulic-radius law. Typical: 0.3–0.5. Lower values weaken the relationship; higher values strengthen it. |
| `sheet_flow_min_depth_metres` | float | `0.001` | Floor on the hydraulic radius (metres), i.e., minimum sheet-flow depth. Prevents unrealistically shallow depths on hillslopes. |
| `velocity_scale` | float | `1.0` | Global multiplier on the entire velocity field. Calibration handle: scales the whole field uniformly so the area-weighted mean travel time can match a target. |
| `target_mean_traveltime_seconds` | float | `0.0` | If set to a positive value, `velocity_scale` is solved iteratively so the area-weighted mean travel time to outlet equals this target (seconds). If `0` or negative, `velocity_scale` is fixed at its specified value. Useful for rough calibration without detailed velocity data. |
| `v_min_metres_per_second` | float | `0.001` | Lower clamp on velocity (m/s). Prevents unrealistically slow flow (e.g., on flat floodplains). |
| `v_max_metres_per_second` | float | `5.0` | Upper clamp on velocity (m/s). Prevents unrealistically fast flow (numerical stability). |
| `v_channel_metres_per_second` | float | `1.0` | Fixed channel velocity (m/s) used only in `velocity_mode: constant`. Ignored in other modes. |
| `v_hillslope_metres_per_second` | float | `0.005` | Fixed hillslope velocity (m/s) used only in `velocity_mode: constant`. Ignored in other modes. |
| `strahler_velocity_table` | string | | Velocity per Strahler order, used only in `velocity_mode: strahler_table`. Format: semicolon-delimited `order:velocity`, e.g., `1:0.005;2:0.01;3:0.05;4:0.2;5:1.0` (m/s). |
| `strahler_order_fname` | string | | Strahler order raster, required when `velocity_mode: strahler_table`. |

---

## Channel extent (optional)

Define which cells are channels for the purpose of roughness or velocity override. If not set, geometry-based hydraulic radius applies everywhere.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `channel_extent_area_threshold_m2` | float | `100000.0` | Drainage area (m²) above which a cell is classified as channel. Used only when `channel_extent_mode: area_threshold`. Typical range: 1e4–1e6 m² (0.01–1 km²). |
| `channel_extent_area_slope_theta` | float | `2.0` | Exponent θ in the area-slope criterion: cells with `A * S^θ > threshold` are channels. Used only when `channel_extent_mode: area_slope_threshold`. Typical: 1.5–2.5. |
| `channel_extent_area_slope_threshold` | float | `100.0` | Threshold on `A * S^θ` defining a channel cell. Typical range: 10–1000 m^(1+θ). |
| `channel_extent_raster_fname` | string | | Channel raster from a prior extraction (binary: 0 = not channel, 1 = channel). Used only when `channel_extent_mode: extraction`. |
| `channel_extent_vector_fname` | string | | Mapped river network (GIS vector: `.shp`, `.fgb`, etc.) to snap to the flow graph as the channel extent. Used only when `channel_extent_mode: vector`. |
| `channel_extent_snap_max_search_radius_pixels` | int | `2` | Search radius (pixels) within which each vector pixel may snap to the flow graph. Prevents snapping across ridges. |
| `channel_extent_snap_min_area_m2` | float | `0.0` | Minimum drainage area (m²) a snapped pixel must reach to be kept as channel. Removes spurious snaps in headwater areas. |

---

## Channel roughness (Manning n variation)

By default, channel cells have a fixed Manning n. Optionally, n can vary with drainage area or channel steepness.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `channel_roughness_mode` | string | `constant` | How channel resistance is derived. **`constant`**: fixed `channel_roughness` (default). **`drainage_area`**: n decreases downstream with drainage area. **`ksn`** and **`slope_area`** use channel steepness. **`grain_size`** uses D50 through Strickler; **`jarrett`** and **`ferguson_vpe`** are depth-dependent high-gradient alternatives. |
| `channel_roughness` | float | `0.035` | Manning n in channels when `channel_roughness_mode: constant`, and the reference n for scaling modes. Typical: 0.025–0.05 (gravel/cobble). |
| `channel_n_area_exponent` | float | `0.08` | Exponent e_A in `n = n_ref * (A / A_ref)^(-e_A)`. Negative exponent means n decreases downstream (typical for channels). Typical: 0.05–0.15. |
| `channel_n_ksn_exponent` | float | `0.17` | Exponent e_k in `n = n_ref * (ksn / ksn_ref)^(e_k)`. Positive exponent means n increases with steepness. Typical: 0.1–0.3. |
| `channel_n_ref_area_m2` | float | `1e6` | Reference drainage area A_ref (m²) for area-based scaling. Typical: 1e5–1e7 m². |
| `channel_n_ref_ksn` | float | `0.0` | Reference ksn for ksn and slope-area roughness modes. Required only if `channel_roughness_mode` uses ksn. |
| `chi_m_over_n` | float | `0.45` | Concavity (m/n exponent) used to form `ksn = S * A^θ` in slope-area mode. Default assumes m/n=0.45 typical of fluvial channels. |
| `ksn_fname` | string | | Normalised channel steepness (ksn) raster, required when `channel_roughness_mode` uses ksn. |
| `d50_fname` | string | | Median grain size (D50) raster (metres), used in grain-size-based roughness. |
| `d84_fname` | string | | 84th-percentile grain size (D84) raster (metres), used in grain-size-based roughness. |
| `strickler_coefficient` | float | `0.047` | Coefficient a_str in the Strickler equation: `n = a_str * D50^(1/6)`. D50 in metres, n dimensionless. |

---

## Surface classification (optional)

Identify floodplain or terrace cells for resistance-partition reporting. At present both floodplain and hillslope resistance use land cover; the distinction prepares for future compound-channel conveyance.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `surface_class_fname` | string | | Surface class raster from `lsdtt-fluvial-surfaces` (floodplain, terrace, fan classifications). Optional; if omitted, all non-channel cells are treated as hillslope. |
| `surface_class_floodplain_code` | float | `1.0` | Value in `surface_class_fname` that denotes floodplain. Floodplain cells receive floodplain-specific roughness if assigned. |

---

## Infiltration excess (optional)

Exponential Green-Ampt infiltration capacity that decays as the wetting front advances. Disabled by default.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_infiltration_excess` | bool | `false` | Enable exponential-K Green-Ampt infiltration excess modelling. When true, rainfall exceeding the current infiltration capacity becomes surface runoff. |
| `xk0_metres_per_hour` | float | `1.0` | Surface saturated hydraulic conductivity Ks (metres/hour) at the wetting front. Typical range: 0.01–10 m/h depending on soil type. |
| `wetting_front_suction_metres` | float | `0.02` | Wetting-front suction head (metres). Typical: 0.01–0.1 m (fine sands to clays). |
| `delta_theta` | float | `0.1` | Volumetric water-content change across the wetting front (dimensionless fraction, 0–1). Typical: 0.05–0.3. |

---

## Forcing and time stepping

Rainfall and potential evapotranspiration are supplied as time series (one value per outer timestep).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `rainfall_fname` | string | | CSV file containing the rainfall time series. One row per timestep with a `rainfall_metres` or `rainfall_mm` column holding the depth during that timestep. Required. |
| `pet_fname` | string | | CSV file containing potential evapotranspiration. One row per timestep with a `pet_metres` or `pet_mm` column holding the depth during that timestep. Omit it for zero PET. |
| `dt_seconds` | float | `3600.0` | Outer timestep duration in seconds. Forcing values are depths for this timestep, not rates. Typical: 600–3600 s. |
| `n_timesteps` | int | `0` | Number of outer timesteps to run. If `0`, runs for as many steps as the rainfall CSV contains. |
| `routing_substeps` | int | `0` | Sub-steps within one outer timestep for surface routing. `0` chooses a deterministic power-of-two count from a Courant estimate; `1` disables sub-cycling. It controls accuracy, not stability. |

---

## Outputs

Select which rasters and tables to write.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_discharge_csv` | bool | `true` | `<write_fname>_discharge.csv` — outlet hydrograph and per-step rainfall, ET, runoff, deficit, contributing-area, and discharge fields. |
| `write_water_balance_csv` | bool | `true` | `<write_fname>_water_balance.csv` — catchment, canopy, and routing water-balance terms, including closure errors. |
| `write_saturation_raster` | bool | `false` | `<write_fname>_saturation.tif` — binary raster of saturated cells (1) vs. unsaturated (0) at model end. |
| `write_water_table_raster` | bool | `false` | `<write_fname>_water_table.tif` — local saturation deficit (metres) at model end. High values = dry, low values = wet. |
| `write_velocity_raster` | bool | `false` | `<write_fname>_velocity.tif` — computed flow velocity field (m/s) at model end. |
| `write_traveltime_raster` | bool | `false` | `<write_fname>_traveltime.tif` — travel time to outlet (seconds) for each cell. |
| `write_twi_raster` | bool | `false` | `<write_fname>_twi.tif` — Topographic Wetness Index: ln((A/pixel_area) / slope). |
| `write_iuh_csv` | bool | `false` | `<write_fname>_iuh.csv` — unit hydrograph (width function or GIUH), if using width-function routing. Columns: lag time (s), relative volume (fraction of peak). |

---

## Example: basic TOPMODEL run

```yaml
# Minimal parameter file for a 1-hour storm
read_fname: dem.tif
write_fname: test_run
n_timesteps: 12
dt_seconds: 300

# Rainfall and PET time series (in read_path)
rainfall_fname: rainfall.csv
pet_fname: pet.csv

# TOPMODEL stores (metres and hours)
srmax_metres: 0.05
szm_metres: 0.032
ln_t0: 5.0
td_hours_per_metre: 50.0

# Routing
routing_method: width_function
velocity_mode: hydraulic_geometry

# Outputs
write_discharge_csv: true
write_water_balance_csv: true
write_traveltime_raster: true
```

---

## Example: dynamic saturated-zone TOPMODEL with interception and channel routing

```yaml
read_fname: dem.tif
write_fname: dynamic_run
rainfall_fname: rain.csv
pet_fname: pet.csv
dt_seconds: 600
n_timesteps: 100

# TOPMODEL stores
saturated_zone_mode: dynamic
initial_deficit_metres: 0.5
srmax_metres: 0.1
szm_metres: 0.05
ln_t0: 5.5
td_hours_per_metre: 40.0

# Velocity
velocity_mode: hydraulic_geometry
hydraulic_radius_coeff: 0.15
hydraulic_radius_exponent: 0.4
velocity_scale: 1.2
target_mean_traveltime_seconds: 1800.0

# Routing
routing_method: muskingum_cunge
routing_substeps: 3

# Channel definition
channel_extent_mode: area_threshold
channel_extent_area_threshold_m2: 50000.0
channel_roughness_mode: drainage_area
channel_n_area_exponent: 0.1
channel_n_ref_area_m2: 1e6

# Interception
run_interception: true
landcover_fname: lulc.tif
landcover_table: "grazing:code=1,lai=1.5,n=0.30,soil=CULT;woodland:code=2,lai=5.0,n=0.60,soil=SEMI"
canopy_storage_mm_per_lai: 0.3
lai_seasonal_factor: 1.0

# Outputs
write_discharge_csv: true
write_water_balance_csv: true
write_water_table_raster: true
write_velocity_raster: true
write_traveltime_raster: true
```

---

## Known limitations and workarounds

1. **Experimental tool:** Version 0.1.0 has unit, smoke, water-balance, and reference-model tests, but it has not been calibrated or evaluated against an observed catchment hydrograph. Treat a study-site result as a model experiment until it has been independently evaluated.

2. **Index-mode assumption:** `saturated_zone_mode: index` follows TOPMODEL's quasi-steady topographic-index relation. Use `dynamic` when a per-cell prognostic deficit is needed, and test timestep sensitivity for the study site.

3. **Forcing is depth per timestep:** A rain-gauge value in millimetres belongs in `rainfall_mm`; do not pre-convert it to a rate. The model converts millimetres to metres at input.

---

## Related CLIs

- **[`lsdtt-drainage-analysis`](cli_reference_lsdtt_drainage_analysis.md)** — Computes static drainage metrics (TWI, STI, connectivity) without time stepping. Use this to inspect catchment properties before running a full simulation.
- **[`lsdtt-graphflood`](cli_reference_lsdtt_graphflood.md)** — Hydraulic routing model (Manning-Chezy) for inundation mapping. Solves depth-averaged flow and water depth on an evolving surface. Complementary to distributed-hydrology: graphflood *routes* a discharge; distributed-hydrology *generates* discharge from rainfall.

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

### Flow routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string for north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. A positive value is honoured verbatim and OVERRIDES OMP_NUM_THREADS, so setting 1 here runs serial on any machine. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling. | Passed to fill_pits() to gently enforce downstream descent across filled depressions and flats. Only used when filling actually happens. |

### Model structure

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `saturated_zone_mode` | `string` | `index` | `param` | - | How the saturated zone is represented. | index: the local deficit is derived every step from the catchment mean and the topographic index, D_i = S_bar + m(lambda - TWI_i). This is classic TOPMODEL and it is the only configuration with an external reference, so it is the validation mode. dynamic: each cell carries its own deficit as a state and the subsurface flux is routed down the flow graph, so cells wet and dry at different rates. |
| `routing_method` | `string` | `width_function` | `param` | - | How runoff is routed to the outlet. | width_function: pure translation through the travel-time distribution, i.e. the geomorphological instantaneous unit hydrograph. Exact for a static velocity field and the default. linear_reservoir: a store per cell with k = dx/v. muskingum_cunge: explicit, with numerical diffusion matched to the physical hydraulic diffusivity. kinematic_wave: implicit four-point, one scalar Newton solve per cell, unconditionally stable. |
| `velocity_mode` | `string` | `hydraulic_geometry` | `param` | - | How the per-cell velocity is derived. | hydraulic_geometry: v = scale (1/n) R^(2/3) S^(1/2) with R from downstream hydraulic geometry, so low-order channels come out slower without any order lookup. roughness_raster: as above but the supplied roughness_fname overrides the hillslope/floodplain/channel partition entirely. strahler_table: one velocity per Strahler order, from strahler_order_fname and strahler_velocity_table. constant: one channel and one hillslope velocity. |
| `channel_extent_mode` | `string` | `none` | `param` | - | How channel cells are identified. | The channel extent selects which RESISTANCE closure applies; it does not decide where water goes, and the model runs correctly without one. none: emergent from the hydraulic geometry, no parameter. area_threshold: drainage area above channel_extent_area_threshold_m2. area_slope_threshold: A S^theta above a threshold (Montgomery &amp; Dietrich). extraction: a channel raster from lsdtt-channel-extraction, named by channel_extent_raster_fname. vector: a mapped network burned in from channel_extent_vector_fname and snapped to the flow field. |
| `flow_method` | `string` | `dinf` | `param` | - | Flow routing method used for drainage area and all graph sweeps. | Allowed values: dinf (Tarboton 1997) or mfd (Quinn et al. 1991). D8 is NOT offered here: it sends all of a cell's area to one neighbour, which stripes the wetness index and therefore every deficit derived from it. |
| `mfd_exponent` | `float` | `1.0` | `param` | - | Slope-weighting exponent for MFD flow partitioning. | Quinn et al. (1991) use 1.0; Freeman (1991) suggest 1.1. Used only when flow_method is mfd. |
| `twi_window_radius_metres` | `float` | `100.0` | `param` | - | Polyfit window radius (metres) for the slope used by TWI and velocity. | Deliberately the same name and meaning as in lsdtt-drainage-analysis: the two CLIs share one TWI implementation and must be able to produce the same raster from the same settings. |

### TOPMODEL core

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `szm_metres` | `float` | `0.032` | `param` | - | The exponential store depth m (metres). | Controls how fast transmissivity falls as the saturation deficit grows, and therefore the recession shape. The single most sensitive TOPMODEL parameter. |
| `ln_t0` | `float` | `5.0` | `param` | - | ln of the saturated transmissivity at zero deficit, T0 in m^2/h. | Used directly when transmissivity_source is uniform. When it is sti_table this value is replaced by the areal mean of the measured ln(T0) field and is used only where that field has no value. |
| `td_hours_per_metre` | `float` | `50.0` | `param` | - | Unsaturated-zone time delay per unit deficit (hours per metre). | Sets the gravity drainage rate q_v = S_uz / (D_i * td), so it controls how long water takes to reach the water table. |
| `srmax_metres` | `float` | `0.05` | `param` | - | Maximum root zone deficit (metres), i.e. plant-available storage. | Evapotranspiration draws the root zone down to this deficit and no further. |
| `q0_metres_per_step` | `float` | `3.28e-5` | `param` | - | Initial subsurface flow per unit area (metres per timestep). | Used to initialise the catchment mean deficit, S_bar = -m ln(q0/szq). Must be greater than zero or that initialisation is undefined. |
| `sr0_metres` | `float` | `0.002` | `param` | - | Initial root zone deficit (metres). | Zero means the root zone starts at field capacity. |
| `initial_deficit_metres` | `float` | `0.3` | `param` | - | Initial per-cell saturation deficit (metres) in dynamic mode. | The prognostic saturated zone has no lambda to initialise from, so every cell starts here and relaxes towards its own equilibrium. Ignored when saturated_zone_mode is index. |
| `saturated_zone_substeps` | `int` | `1` | `param` | - | Sub-steps within one outer step for the saturated-zone routing. | Dynamic mode only. 0 lets the module choose from how hard the outflow is pushing against the water available, snapped to a power of two. Raise this if the run reports many stiff cells. |

### Transmissivity

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `transmissivity_source` | `string` | `uniform` | `param` | - | Where the saturated transmissivity T0 comes from. | uniform: one value for the whole catchment, ln_t0. sti_table: a measured, spatially distributed field, looked up per pixel on (soil series, land use) from the same profile database the Soil Topographic Index uses. Land use is an INDEX INTO MEASURED DATA there, not a free parameter, so a land-cover scenario re-keys the lookup rather than applying a multiplier. Requires sti_soil_fname, sti_ksat_fname, landcover_fname and landcover_table with soil= keys. |
| `sti_soil_fname` | `string` | - | `param` | - | Soil series polygon vector keying the transmissivity lookup. | Rasterised on the DEM grid and joined to the measured profiles. This is the same input lsdtt-drainage-analysis calls sti_soil_fname. |
| `sti_soil_layer_name` | `string` | - | `param` | - | Layer name within sti_soil_fname. | Required when transmissivity_source is sti_table. |
| `sti_soil_field` | `string` | - | `param` | - | Attribute field in sti_soil_fname holding the series code. | The code is matched against the (series, land use) keys of the measured profile table. |
| `sti_ksat_fname` | `string` | - | `param` | - | CSV of measured saturated conductivity profiles. | Conductivity is integrated over depth per profile, which IS transmissivity, so this file supplies ln(T0) by (series, land use) rather than a calibrated scalar. It is the same input lsdtt-drainage-analysis calls sti_ksat_fname. |
| `sti_depth_cap_centimetres` | `float` | `100.0` | `param` | - | Profile depth (cm) beyond which horizons are ignored. | Sets the depth over which conductivity is integrated to give transmissivity. Same name and meaning as in lsdtt-drainage-analysis. |
| `sti_min_valid_depth_centimetres` | `float` | `30.0` | `param` | - | Shortest profile (cm) accepted into the lookup table. | Profiles thinner than this are dropped, so a series described only by a topsoil horizon does not set the value for a whole polygon. |
| `sti_ks_unit_scale_to_metres_per_day` | `float` | `1.0` | `param` | - | Multiplier converting the CSV conductivity column to metres per day. | Transmissivity enters the model inside an exponential, so a wrong scale here is a multiplicative error in baseflow, not a shift. |
| `write_ln_t0_raster` | `bool` | `false` | `switch` | - | Write the measured ln(T0) field as a raster. | Output file: &lt;write_fname&gt;_ln_t0.tif. Only produced when transmissivity_source is sti_table. |

### Interception

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_interception` | `bool` | `false` | `switch` | - | Run the canopy interception store above the root zone. | Interception is not a skim off the top: wet-canopy evaporation runs at the potential rate and suppresses transpiration, so only the residual demand reaches the root zone. Needs LAI, from lai_fname or from landcover_fname plus landcover_table. |
| `landcover_fname` | `string` | - | `param` | - | Land-cover class raster. | Integer class codes matching the code= field of landcover_table. Drives LAI and the hillslope/floodplain roughness, and selects which measured soil profile supplies transmissivity. It does NOT set channel roughness, which is grain-size controlled. |
| `landcover_table` | `string` | - | `param` | - | Land-cover class definitions. | Format, with NO SPACES anywhere (a parameter value is a single token): 'grazing:code=1,lai=1.5,n=0.30,soil=CULT;woodland:code=2,lai=5.0,n=0.60,soil=SEMI'. code, lai and n are required; soil is the land-use key used against the measured profile table. A ks= or kv= field is a hard error: conductivity is not a land-cover property here. |
| `lai_fname` | `string` | - | `param` | - | Leaf area index raster. | Overrides the lai= values of landcover_table when supplied, so a measured or remotely sensed LAI can be swapped in by configuration rather than by code. |
| `canopy_storage_mm_per_lai` | `float` | `0.2` | `param` | `run_interception` | Canopy storage capacity per unit LAI (millimetres). | C = k_S * LAI, with k_S about 0.2 mm (Dickinson 1984; van Dijk &amp; Bruijnzeel 2001). |
| `lai_gap_coefficient` | `float` | `0.6` | `param` | `run_interception` | Beer-Lambert extinction coefficient for the free throughfall fraction. | p = exp(-kappa * LAI), the share of rain reaching the ground without touching the canopy (Aston 1979). Usually 0.5 to 0.7. |
| `lai_seasonal_factor` | `float` | `1.0` | `param` | `run_interception` | Multiplier applied to every LAI before the canopy is built. | Winter deciduous LAI is about 0.5 to 1 against about 5 in summer, an order of magnitude in interception capacity, and UK floods are winter dominated. Running a winter storm at annual-mean or summer LAI overstates the benefit of planting. Choose this deliberately. |
| `interception_evap_factor` | `float` | `1.0` | `param` | `run_interception` | Wet-canopy evaporation as a multiple of the supplied PET. | Greater than one for tall aerodynamically rough canopies, which are better ventilated than the grass reference surface, so total evapotranspiration can legitimately exceed PET. |

### Land-cover scenarios

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `landcover_scenario` | `string` | - | `param` | - | Reclassification rules applied before the run. | Format, with no spaces: 'grazing-&gt;woodland;arable-&gt;woodland'. A source of '*' matches any class and rules are applied in order. Roughness, LAI and transmissivity all move together, because they are consequences of one physical change. |
| `landcover_intervention_mask_fname` | `string` | - | `param` | - | Raster limiting where landcover_scenario applies. | Non-zero marks a cell the scenario may reclassify. Leave empty to apply the scenario everywhere. |

### Velocity and resistance

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `roughness_fname` | `string` | - | `param` | - | Manning n raster. | Required by velocity_mode roughness_raster, where it replaces the computed hillslope/floodplain/channel partition outright. |
| `roughness_default` | `float` | `0.035` | `param` | - | Manning n used on hillslope and floodplain cells with no land cover. | Applied where no land-cover raster is supplied, or where it is nodata. Channel cells take the channel_roughness_mode law instead. |
| `hydraulic_radius_coeff` | `float` | `0.1` | `param` | - | Coefficient a_R in R = max(h_sheet, a_R (A/A_ref)^c_R). | The hydraulic radius in metres at the reference drainage area. |
| `hydraulic_radius_exponent` | `float` | `0.4` | `param` | - | Exponent c_R in R = max(h_sheet, a_R (A/A_ref)^c_R). | Leopold &amp; Maddock (1953) give a downstream depth exponent near 0.4. This is what makes velocity rise along the network, and therefore what makes low-order channels slower, with no Strahler lookup. |
| `sheet_flow_min_depth_metres` | `float` | `0.001` | `param` | - | Floor on the hydraulic radius (metres), i.e. sheet flow depth. | This is what makes the hillslope-to-channel transition smooth: a channel is simply where hydraulic-geometry depth has overtaken sheet flow depth. There is no threshold and no discontinuity. |
| `velocity_scale` | `float` | `1.0` | `param` | - | One global multiplier on the whole velocity field. | A per-cell velocity field cannot be calibrated cell by cell, so the field sets the shape and this sets the magnitude (Maidment et al. 1996). Applied before the min/max clamp. |
| `target_mean_traveltime_seconds` | `float` | `0.0` | `param` | - | Area-weighted mean travel time to solve velocity_scale against. | Greater than zero replaces velocity_scale with the value that makes the mean travel time match this target. That is a closed form, not a search, because travel time scales exactly as 1/velocity. |
| `v_min_metres_per_second` | `float` | `0.001` | `param` | - | Lower clamp on cell velocity (metres per second). | Applied after velocity_scale, so it is a statement about the physical range rather than about the calibration. |
| `v_max_metres_per_second` | `float` | `5.0` | `param` | - | Upper clamp on cell velocity (metres per second). | Applied after velocity_scale. |
| `v_channel_metres_per_second` | `float` | `1.0` | `param` | - | Channel velocity (metres per second) in velocity_mode constant. | Cells at or above channel_extent_area_threshold_m2 take this value. |
| `v_hillslope_metres_per_second` | `float` | `0.005` | `param` | - | Hillslope velocity (metres per second) in velocity_mode constant. | Also used for cells with no Strahler order in strahler_table mode. |
| `strahler_velocity_table` | `string` | - | `param` | - | Velocity per Strahler order, in metres per second. | Format, with no spaces: '1:0.2,2:0.4,3:0.8,4:1.2,5:1.8'. Orders above the largest key take the largest key's velocity. Required by velocity_mode strahler_table, which is a diagnostic override: the ordering is emergent under hydraulic_geometry and does not need a table. |
| `strahler_order_fname` | `string` | - | `param` | - | Strahler order raster. | Required by velocity_mode strahler_table. Produce one with lsdtt-channel-network-analysis; this CLI does not extract a network. |

### Channel extent

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `channel_extent_area_threshold_m2` | `float` | `100000.0` | `param` | - | Drainage area (m^2) above which a cell is channel. | Used by channel_extent_mode area_threshold, and by velocity_mode constant to split channel from hillslope. |
| `channel_extent_area_slope_theta` | `float` | `2.0` | `param` | - | Exponent theta in the A S^theta channel-head criterion. | Montgomery &amp; Dietrich (1988, 1992): channel heads sit at smaller area where the ground is steeper. |
| `channel_extent_area_slope_threshold` | `float` | `100.0` | `param` | - | Threshold on A S^theta defining a channel cell. | Units depend on theta, so this is a calibrated number for a given landscape rather than a transferable one. |
| `channel_extent_raster_fname` | `string` | - | `param` | - | Channel raster from a prior extraction. | Any non-zero value marks a channel cell. Required by channel_extent_mode extraction; produce one with lsdtt-channel-extraction rather than reinventing the extraction here. |
| `channel_extent_vector_fname` | `string` | - | `param` | - | Mapped river network to burn in as the channel extent. | For example OS Open Rivers or the MasterMap Water Network. Required by channel_extent_mode vector. Burned pixels are snapped to the largest nearby drainage area, because a survey and a DEM routinely disagree by a pixel or two. |
| `channel_extent_snap_max_search_radius_pixels` | `int` | `2` | `param` | - | Radius (pixels) each burned vector pixel may be snapped within. | 0 takes the burn as drawn. Only used by channel_extent_mode vector. |
| `channel_extent_snap_min_area_m2` | `float` | `0.0` | `param` | - | Drainage area (m^2) a snapped pixel must reach to be kept. | Stops a mis-registered line from painting a hillslope as channel. Rejected pixels are counted and reported. 0 disables the check. |

### Channel roughness

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `channel_roughness_mode` | `string` | `constant` | `param` | - | How roughness varies within channels. | Channel n is set by the BED MATERIAL, never by land cover, so a planting scenario leaves it untouched. constant; drainage_area (downstream fining, n falls as A^-e_A); ksn or slope_area (steeper is coarser is rougher); grain_size (Strickler, from d50_fname); jarrett (high-gradient streams, depth dependent); ferguson_vpe (the variable power equation, for low relative submergence where Manning is the wrong shape, from d84_fname). |
| `channel_roughness` | `float` | `0.035` | `param` | - | Manning n in channels, and n_ref for the scaling modes. | The value the drainage_area, ksn and slope_area laws are anchored on at their reference area or steepness. |
| `channel_n_area_exponent` | `float` | `0.08` | `param` | - | Exponent e_A in n = n_ref (A/A_ref)^(-e_A). | Sternberg fining gives D50 ~ A^-b over b about 0.3 to 0.7, and Strickler gives n ~ D50^(1/6), so e_A lands near 0.05 to 0.12. |
| `channel_n_ksn_exponent` | `float` | `0.17` | `param` | - | Exponent e_k in n = n_ref (ksn/ksn_ref)^(e_k). | D50 scales with slope at fixed area and n ~ D50^(1/6), which gives about 0.17. |
| `channel_n_ref_area_m2` | `float` | `1e6` | `param` | - | Reference drainage area A_ref (m^2). | Anchors both the hydraulic-radius closure and the drainage_area roughness law, so channel_roughness is the value AT this area. |
| `channel_n_ref_ksn` | `float` | `0.0` | `param` | - | Reference ksn for the ksn and slope_area roughness modes. | 0 uses the median over channel cells, which centres the field on channel_roughness rather than on an arbitrary absolute value. |
| `chi_m_over_n` | `float` | `0.45` | `param` | - | Concavity used to form ksn = S A^theta in slope_area mode. | The usual reference concavity for bedrock rivers. |
| `ksn_fname` | `string` | - | `param` | - | Normalised channel steepness raster. | Required by channel_roughness_mode ksn. Produce one with the chi machinery in lsdtt-chi-analysis; use slope_area instead to get the same physics without chi segmentation. |
| `d50_fname` | `string` | - | `param` | - | Median grain size raster, in metres. | Required by channel_roughness_mode grain_size, which is the physically direct option: n = a_str D50^(1/6). |
| `d84_fname` | `string` | - | `param` | - | 84th-percentile grain size raster, in metres. | Required by channel_roughness_mode ferguson_vpe, whose resistance law is written in terms of relative submergence R/D84. |
| `strickler_coefficient` | `float` | `0.047` | `param` | - | Coefficient a_str in n = a_str D50^(1/6), D50 in metres. | Strickler (1923). Used by channel_roughness_mode grain_size. |

### Surface class

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `surface_class_fname` | `string` | - | `param` | - | Surface class raster from lsdtt-fluvial-surfaces. | The &lt;stem&gt;_surface_class.tif output, used to separate floodplain from hillslope. Both draw roughness from land cover today, so this changes the reported partition statistics and nothing else; it becomes load bearing when compound-channel conveyance arrives. |
| `surface_class_floodplain_code` | `float` | `1.0` | `param` | - | Value in surface_class_fname that means floodplain. | Cells carrying this value are counted as floodplain rather than hillslope. |

### Infiltration excess

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_infiltration_excess` | `bool` | `false` | `switch` | - | Run the exponential-K Green-Ampt infiltration excess store. | Beven (1984). Rain arriving faster than the surface can absorb becomes runoff without ever entering the column. Saturation excess is always on; this adds the Hortonian path. |
| `xk0_metres_per_hour` | `float` | `1.0` | `param` | `run_infiltration_excess` | Surface saturated hydraulic conductivity (metres per hour). | The K at the soil surface; it falls with depth on the same exponential scale szm_metres that the saturated zone uses. |
| `wetting_front_suction_metres` | `float` | `0.02` | `param` | `run_infiltration_excess` | Wetting front suction head (metres). | Enters only as the product with delta_theta, so the two are not separately identifiable from a hydrograph. |
| `delta_theta` | `float` | `0.1` | `param` | `run_infiltration_excess` | Water content change across the wetting front. | Enters only as the product with wetting_front_suction_metres. |

### Forcing and time

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `rainfall_fname` | `string` | - | `param` | - | CSV of the rainfall time series, one row per timestep. | Must carry a column named rainfall_metres or rainfall_mm, holding the depth falling during each step. Required: there is no default storm. |
| `pet_fname` | `string` | - | `param` | - | CSV of the potential evapotranspiration series, one row per timestep. | Must carry a column named pet_metres or pet_mm. Leave empty for zero PET, which is reasonable for a single winter storm and not for a seasonal run. |
| `dt_seconds` | `float` | `3600.0` | `param` | - | The outer timestep (seconds). | The forcing and the reporting both live on this step and it never moves: there is no adaptive stepping, because the routing schemes have no stability limit to dodge. Sub-cycling inside the routing operator is available separately. |
| `n_timesteps` | `int` | `0` | `param` | - | Number of timesteps to run. | 0 runs the length of the rainfall series. A larger value continues past the end of the series with zero rainfall, which is how you let a recession play out. |
| `routing_substeps` | `int` | `0` | `param` | - | Sub-steps within one outer step for the surface routing. | 0 chooses the count from a Courant estimate and snaps it to a power of two, which keeps it deterministic. 1 disables sub-cycling. This is an accuracy control, not a stability one. |

### Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_discharge_csv` | `bool` | `true` | `switch` | - | Write the outlet hydrograph and the per-step fluxes to a CSV table. | Output file: &lt;write_fname&gt;_discharge.csv. |
| `write_water_balance_csv` | `bool` | `true` | `switch` | - | Write the run water balance to a CSV table. | Output file: &lt;write_fname&gt;_water_balance.csv. The closure error is a small difference of large numbers and is the best single tripwire this model has for a structural mistake, so it is a first-class output and not a line of console text. |
| `write_saturation_raster` | `bool` | `false` | `switch` | - | Write a raster of the saturated cells at the end of the run. | Output file: &lt;write_fname&gt;_saturation.tif. 1 where the local saturation deficit has reached zero, 0 elsewhere. |
| `write_water_table_raster` | `bool` | `false` | `switch` | - | Write a raster of the local saturation deficit at the end of the run. | Output file: &lt;write_fname&gt;_water_table.tif, in metres below saturation. |
| `write_velocity_raster` | `bool` | `false` | `switch` | - | Write the velocity field as a raster. | Output file: &lt;write_fname&gt;_velocity.tif, in metres per second. |
| `write_traveltime_raster` | `bool` | `false` | `switch` | - | Write the travel time to the outlet as a raster. | Output file: &lt;write_fname&gt;_traveltime.tif, in seconds. |
| `write_twi_raster` | `bool` | `false` | `switch` | - | Write the Topographic Wetness Index as a raster. | Output file: &lt;write_fname&gt;_twi.tif. Deliberately the same switch name as in lsdtt-drainage-analysis: both CLIs call one shared TWI implementation, so for the same DEM and the same window radius they produce the same raster. |
| `write_iuh_csv` | `bool` | `false` | `switch` | - | Write the width-function unit hydrograph to a CSV table. | Output file: &lt;write_fname&gt;_iuh.csv. This is the geomorphological instantaneous unit hydrograph and it integrates to one. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format: cog, gtiff, png or same_as_input. | cog is a Cloud Optimised GeoTIFF and is the project default. It falls back to gtiff for small rasters, where the tiling overhead is not worth it. There is ONE raster format parameter per CLI: it governs every raster the tool writes. Vector outputs are governed separately by write_vector_format. |

<!-- END GENERATED CLI PARAMETERS -->
