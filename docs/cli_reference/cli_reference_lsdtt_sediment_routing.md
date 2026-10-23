# CLI reference: lsdtt-sediment-routing

One-dimensional, reach-scale bedload routing with several grain-size classes
on a D8 channel network, built from a port of d-CASCADE (Battista et al.
2022). Reach length, slope and drainage area come from the DEM. Channel width
and the grain-size distribution are set in the parameter file. Discharge comes
from a drainage-area power law, or from a geomorphic instantaneous unit
hydrograph (GIUH) with `use_giuh_discharge: true`.

## Usage

```
lsdtt-sediment-routing  <path> <parameter_file>
lsdtt-sediment-routing  <parameter_file>        # path defaults to ./
lsdtt-sediment-routing  -h                      # print help and exit
lsdtt-sediment-routing  -v                      # print version and exit
lsdtt-sediment-routing  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

Keys use `snake_case`. All keys and their defaults are listed under
[Parameters](#parameters) below.

---

## Pipeline

1. Load the DEM and fill pits (skipped with `raster_is_filled: true`).
2. Route flow with D8 and extract the channel network
   (`threshold_contributing_pixels`, default 1000). Each network link is a
   reach.
3. Give each reach a width, `width_coefficient × A^width_exponent` (A in m²),
   and a discharge, either `discharge_coefficient × A^discharge_exponent` or,
   with `use_giuh_discharge: true`, a synthetic rainfall pulse convolved with
   the network's width-function GIUH (Rodríguez-Iturbe and Valdés 1979).
4. March through the network from upstream to downstream `n_timesteps` times,
   each step `timestep_s` seconds long (defaults: 365 daily steps). Transport
   capacity per grain-size class comes from `transport_formula_choice`:
   `wilcock_crowe` (default; Wilcock and Crowe 2003), `wong_parker` (Wong and
   Parker 2006, with the Parker 1982 hiding function) or `engelund_hansen`
   (Engelund and Hansen 1967, D50 only). The active layer is
   `active_layer_depth_metres` thick, or 2 × D90 per reach when that is ≤ 0
   (default).

The grain-size classes are `psi_mm` (class sizes in mm) with the fractions in
`gsd_fractions`; both are comma-separated lists of the same length.

## Outputs

| File | Written when | Contents |
|---|---|---|
| `<write_fname>_cumulative_deposition.tif` | `write_deposition_raster: true` (default) | Cumulative deposited volume over the run (m³). Every pixel of a reach carries that reach's total, so do not sum pixels. |
| `<write_fname>_d50_al.tif` | `write_d50_raster: true` (default) | Median grain size of the active layer at the end of the run (m) |

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

### Threading

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. Applies to pit filling and flow-routing construction. |

### DEM and flow routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for D8 flow routing. | Four-character string for north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope imposed by pit filling. | Used only when raster_is_filled is false. Passed to fill_pits() to gently enforce downstream descent. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Transport formula

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `transport_formula_choice` | `string` | `wilcock_crowe` | `param` | - | Bedload transport capacity formula. | wilcock_crowe: Wilcock &amp; Crowe (2003) mixed grain-size formula. Accounts for sand fraction and grain-size-dependent hiding. wong_parker: Wong &amp; Parker (2006) with Parker (1982) hiding function. Best for unimodal gravel beds. engelund_hansen: Engelund &amp; Hansen (1967) total-load formula using D50. Simple but does not resolve grain-size classes. |

### Routing time control

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_timesteps` | `int` | `365` | `param` | - | Number of routing timesteps. | The model marches upstream-to-downstream this many times. Use 365 with timestep_s = 86400 for one year of daily steps. |
| `timestep_s` | `float` | `86400` | `param` | - | Duration of each routing timestep in seconds. | Default is one day (86400 s). Must be consistent with the discharge input if use_giuh_discharge is true. |

### Sediment properties

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `porosity` | `float` | `0.4` | `param` | - | Bed sediment porosity (dimensionless). | Used to convert between bulk volume and solid volume when updating the active layer. Typical gravel beds: 0.35-0.45. |
| `sediment_density_kg_m3` | `float` | `2650` | `param` | - | Sediment grain density [kg/m^3]. | Quartz density is 2650 kg/m^3. Used in Shields stress and transport capacity calculations. |
| `active_layer_depth_metres` | `float` | `-1` | `param` | - | Fixed active-layer depth [m] applied to all reaches. | The active layer is the mobile surface sediment that participates in exchange each timestep. Values &lt;= 0 (default) use 2 x D90 per reach (Parker 1991), which scales with the coarsest grains. |

### Grain-size distribution

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `psi_mm` | `string` | `4,8,16,32,64` | `param` | - | Grain-class centre diameters [mm], comma-separated. | Defines the grain-size classes used in transport calculations. Values correspond to phi-scale class centres (4, 8, 16, 32, 64 mm span fine gravel to coarse gravel). Must have the same number of entries as gsd_fractions. |
| `gsd_fractions` | `string` | `0.10,0.25,0.35,0.20,0.10` | `param` | - | Mass fraction per grain class, comma-separated (must sum to 1). | Applied uniformly to all reaches unless overridden. Must have the same number of entries as psi_mm. |
| `d50_metres` | `float` | `0.016` | `param` | - | Median grain diameter [m] used in transport and hydraulic calculations. | Used directly by the Engelund-Hansen formula and as the reference grain size for hiding-function calculations in Wong-Parker. |
| `d84_metres` | `float` | `0.045` | `param` | - | D84 grain diameter [m]. | 84th percentile grain diameter. Used in Wilcock-Crowe hiding function. |
| `d90_metres` | `float` | `0.055` | `param` | - | D90 grain diameter [m]. | 90th percentile grain diameter. Sets the default active-layer depth (2 x D90) when active_layer_depth_metres &lt;= 0. |
| `mannings_n` | `float` | `0.04` | `param` | - | Manning's roughness coefficient for flow-depth estimation. | Used to compute normal flow depth from discharge, width, and slope via Manning's equation. Typical gravel-bed rivers: 0.03-0.07. |

### Channel width (hydraulic geometry)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `width_coefficient` | `float` | `2.76` | `param` | - | Hydraulic-geometry coefficient for channel width. | Width [m] = width_coefficient x drainage_area [m^2]^width_exponent. Default (2.76) is from Leopold &amp; Maddock (1953) for gravel-bed rivers. |
| `width_exponent` | `float` | `0.40` | `param` | - | Hydraulic-geometry exponent for channel width. | Typical range: 0.35-0.50 for natural rivers. |

### Discharge (hydraulic geometry fallback)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `discharge_coefficient` | `float` | `1e-5` | `param` | - | Hydraulic-geometry coefficient for mean discharge when GIUH is not used. | Q [m^3/s] = discharge_coefficient x drainage_area [m^2]^discharge_exponent. This is a simple power-law scaling. Use use_giuh_discharge for a physics-based estimate. |
| `discharge_exponent` | `float` | `0.80` | `param` | - | Hydraulic-geometry exponent for mean discharge. | Typical range: 0.7-0.9. Used only when use_giuh_discharge is false. |

### GIUH discharge estimation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `use_giuh_discharge` | `bool` | `false` | `switch` | - | Use the geomorphic IUH to estimate discharge for sediment routing. | When true, the width-function-based IUH (Rodriguez-Iturbe &amp; Valdes 1979) is computed from the channel network geometry. A synthetic rainfall pulse is convolved with the IUH to produce per-reach discharge at each timestep. When false, discharge is estimated from drainage area using the discharge_coefficient and discharge_exponent power law. |
| `giuh_v_ref_m_per_s` | `float` | `1.0` | `param` | `use_giuh_discharge` | Reference channel velocity [m/s] for GIUH travel-time estimation. | Travel velocity in each reach is scaled from this reference: v = v_ref x (slope / slope_ref)^velocity_exponent. Typical bankfull channel velocities are 0.5-2.0 m/s. |
| `giuh_slope_ref` | `float` | `0.01` | `param` | `use_giuh_discharge` | Reference slope corresponding to giuh_v_ref_m_per_s. | A reach with slope equal to giuh_slope_ref is assigned velocity giuh_v_ref_m_per_s. Steeper reaches are faster; shallower are slower. |
| `giuh_velocity_exponent` | `float` | `0.5` | `param` | `use_giuh_discharge` | Exponent beta in the kinematic-wave velocity scaling v proportional_to slope^beta. | beta = 0.5 is consistent with Manning's equation and the kinematic wave approximation. Typical range: 0.3-0.6. |
| `giuh_hillslope_length_metres` | `float` | `-1` | `param` | `use_giuh_discharge` | Mean hillslope length [m] added as overland-flow travel time per reach. | When &gt; 0, an overland-flow travel time = hillslope_length / v_hillslope is added to each reach, delaying the IUH peak. Values &lt;= 0 disable the hillslope contribution. |
| `giuh_v_hillslope_m_per_s` | `float` | `0.005` | `param` | `use_giuh_discharge` | Overland-flow velocity [m/s] for hillslope travel time. | Only used when giuh_hillslope_length_metres &gt; 0. Typical sheet-flow velocities: 0.001-0.01 m/s. |
| `giuh_n_bins` | `int` | `200` | `param` | `use_giuh_discharge` | Number of time bins used to discretise the IUH. | More bins give a smoother IUH but use more memory. 200 bins at 1-hour resolution covers 200 hours of response time. |
| `giuh_dt_bin_s` | `float` | `3600` | `param` | `use_giuh_discharge` | Width of each IUH time bin in seconds. | Set to match or be smaller than timestep_s for accurate convolution. Default is 1 hour (3600 s). |
| `giuh_runoff_coefficient` | `float` | `1.0` | `param` | `use_giuh_discharge` | Fraction of rainfall that becomes runoff (dimensionless, 0-1). | Applied to the rainfall rate before convolving with the IUH. 1.0 = all rain reaches the channel (no infiltration losses). Lower values reduce the discharge estimate. |
| `giuh_rainfall_rate_m_per_s` | `float` | `1e-5` | `param` | `use_giuh_discharge` | Rainfall rate [m/s] for the synthetic IUH pulse (~= 36 mm/hr). | A uniform pulse of this intensity lasting giuh_pulse_duration_s is convolved with the IUH to produce the discharge time series. Adjust to represent the design storm relevant to your setting. |
| `giuh_pulse_duration_s` | `float` | `86400` | `param` | `use_giuh_discharge` | Duration of the synthetic rainfall pulse in seconds. | The pulse is applied uniformly over this duration at giuh_rainfall_rate_m_per_s. Default is 86400 s (one day), matching the default timestep_s. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | cog: Cloud Optimised GeoTIFF (default). gtiff: plain compressed GeoTIFF. same_as_input: match the input raster format. |
| `write_deposition_raster` | `bool` | `true` | `switch` | - | Write cumulative sediment deposition [m^3] per reach as a raster. | Each channel pixel in a reach is assigned the cumulative deposition volume for that reach summed over all timesteps. Output file: &lt;write_fname&gt;_cumulative_deposition.tif. |
| `write_d50_raster` | `bool` | `true` | `switch` | - | Write the final active-layer D50 [m] per reach as a raster. | D50 is computed from the grain-size distribution of the active layer at the end of the last timestep. Output file: &lt;write_fname&gt;_d50_al.tif. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

<!-- END GENERATED CLI PARAMETERS -->
