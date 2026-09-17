# CLI reference: lsdtt-distributed-hydrology

Runs a time-evolving distributed hydrological model combining TOPMODEL saturation-excess routing with optional interception, infiltration excess, and channel routing. Simulates the full water balance from rainfall through groundwater and surface routing to outlet discharge.

This CLI is **hydrological** (generates discharge from rainfall), not hydraulic. For inundation mapping or hydraulic routing of a known discharge, see [`lsdtt-graphflood`](cli_reference_lsdtt_graphflood.md). For static drainage-related metrics (TWI, STI, connectivity), see [`lsdtt-drainage-analysis`](cli_reference_lsdtt_drainage_analysis.md).

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
| `saturated_zone_mode` | string | `index` | How the saturated zone is represented. **`index`**: each cell's saturation deficit is derived from the catchment-mean deficit and the topographic index (classic TOPMODEL; only configuration with an external reference oracle). **`dynamic`**: each cell carries its own deficit as a state variable and subsurface flux is routed down the flow graph (cells wet and dry at different times; more realistic but unvalidated). |
| `routing_method` | string | `width_function` | How runoff is routed to the outlet. **`width_function`** (default): pure translation using precomputed travel-time distribution (the width function or GIUH; exact for static velocity, one topological sweep to build). **`linear_reservoir`**, **`muskingum_cunge`**, **`kinematic_wave`**: explicit dynamic routing with per-cell storage or momentum terms. **Important limitation:** dynamic routers are only active where a polyfit slope exists, which insets the routed region from the raster edge by the slope window radius. If a catchment outlet sits at or near the DEM edge, discharge may be lost. Use `width_function` or inset the domain. |
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
| `saturated_zone_substeps` | int | `1` | Number of sub-steps within one outer timestep for saturated-zone routing in dynamic mode. Improves stability if subsurface wave speed exceeds the Courant limit. Used only in dynamic mode. Typical: 1–5. |

---

## Transmissivity and soil data

Transmissivity (saturated conductivity × depth) can come from a uniform value, a raster, or a lookup table keyed on soil series and land-use class.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `transmissivity_source` | string | `uniform` | Where the saturated transmissivity T0 comes from. **`uniform`**: use the scalar `ln_t0` everywhere. **`sti_table`**: look up T0 from a soil conductivity table keyed on soil series (polygon) and land-use class (raster). **Note:** `sti_table` mode is compiled and unit-tested but has not been run end-to-end via the CLI; use with caution and report issues. |
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
| `landcover_scenario` | string | | Reclassification rules applied to `landcover_fname` before the run. Format: semicolon-delimited rules, each `old_code:new_code`. Example: `1:2;3:2` reclassifies codes 1 and 3 to 2. Ignored if not set. |
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
| `channel_extent_area_threshold_m2` | float | `100000.0` | Drainage area (m²) above which a cell is classified as channel. Used only when `channel_extent_mode: area_threshold`. Typical range: 1e4–1e6 m² (1–100 km²). |
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
| `channel_roughness_mode` | string | `constant` | How Manning n varies within channels. **`constant`**: fixed `channel_roughness` (default). **`area`**: n varies inversely with drainage area. **`ksn`**: n varies with channel steepness (ksn). **`slope_area`**: n varies with both slope and area (compound exponent). |
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

Assign floodplain or terrace cells a separate roughness or velocity (usually lower than hillslope).

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
| `rainfall_fname` | string | | CSV file containing the rainfall time series. One row per timestep; first column is the rate (metres per timestep). Required unless `n_timesteps: 0`. |
| `pet_fname` | string | | CSV file containing the potential evapotranspiration (PET) time series. One row per timestep; first column is the rate (metres per timestep). Required unless `n_timesteps: 0`. |
| `dt_seconds` | float | `3600.0` | Outer timestep duration (seconds). All rainfall and PET values are interpreted as rates per this timestep. Typical: 600–3600 s (10 min – 1 hour). |
| `n_timesteps` | int | `0` | Number of outer timesteps to run. If `0`, runs for as many steps as the rainfall CSV contains. |
| `routing_substeps` | int | `0` | Sub-steps within one outer timestep for surface routing in dynamic modes. Improves stability if flow speeds are high. If `0`, no sub-stepping. Typical: 1–10. |

---

## Outputs

Select which rasters and tables to write.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_discharge_csv` | bool | `true` | `<write_fname>_discharge.csv` — outlet hydrograph and per-step fluxes. Columns: time (s), Q (m³/s), flux in, flux out, internal balance. |
| `write_water_balance_csv` | bool | `true` | `<write_fname>_water_balance.csv` — catchment-wide water balance summary. Columns: total input (rainfall), total output (discharge + ET), residuals, store changes. |
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
channel_roughness_mode: area
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

1. **Dynamic routing inset from raster edge:** The dynamic routers (`linear_reservoir`, `muskingum_cunge`, `kinematic_wave`) are only active where a polyfit slope exists, which insets the routed region from the raster margin by the slope window radius (typically 100 m). If a catchment outlet sits on or near the DEM edge, discharge may be lost. **Workaround:** Use `width_function` routing (the default), which is unaffected and is exact for static velocity fields. Alternatively, inset the catchment domain so the outlet is well within the interior.

2. **`transmissivity_source: sti_table` is not yet validated end-to-end:** This mode compiles and its components are unit-tested, but the CLI wiring has never been exercised in a full run. Use with caution and report issues.

3. **`saturated_zone_mode: dynamic` is unvalidated:** The dynamic deficit routing is physically motivated and numerically stable but has not been tested against observed hydrographs. Treat outputs as exploratory.

---

## Related CLIs

- **[`lsdtt-drainage-analysis`](cli_reference_lsdtt_drainage_analysis.md)** — Computes static drainage metrics (TWI, STI, connectivity) without time stepping. Use this to inspect catchment properties before running a full simulation.
- **[`lsdtt-graphflood`](cli_reference_lsdtt_graphflood.md)** — Hydraulic routing model (Manning-Chezy) for inundation mapping. Solves depth-averaged flow and water depth on an evolving surface. Complementary to distributed-hydrology: graphflood *routes* a discharge; distributed-hydrology *generates* discharge from rainfall.
