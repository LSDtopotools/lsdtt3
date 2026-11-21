# lsdtt-slope-stability

!!! warning "Experimental"
    This CLI is under development. Its first component, the SINMAP stability
    index, is tested against analytic solutions; the user guide will follow once
    the remaining components are built and measured.

## Intent

`lsdtt-slope-stability` estimates where slopes are prone to failure, as maps
and statistics. It is not a tool for mapping individual landslides.

This version computes a stability index inspired by SINMAP (Pack, Tarboton,
Goodwin & Prasad 2005), which itself built on SHALSTAB (Montgomery & Dietrich
1994). It is an infinite-slope factor of safety with steady-state wetness, in
which cohesion, the friction coefficient and the ratio of transmissivity to
recharge are uncertain between bounds. Like SINMAP it applies to shallow
translational failure, not to deep-seated slides.

## How this differs from SINMAP

lsdtt3 keeps SINMAP's model, stability index and six classes, but it is not a
replica:

1. **Contour width.** SINMAP divides contributing area by the cell size
   whatever the flow direction. On a D-infinity grid that underestimates the
   specific catchment area by up to 29 % where flow is diagonal. By default
   lsdtt3 divides by the width of the D-infinity flow tube,
   cell size x max(|sin β|, |cos β|) for flow direction β, which is exact on
   planes. `contour_width_method: cell_size` gives SINMAP's choice.
2. **Routing surface.** Depressions are carved (least-cost breach) and then
   filled; SINMAP only fills. `breach_before_fill: false` gives SINMAP's choice.
3. **Probability of stability.** Computed from closed forms re-derived from the
   model and checked against Monte Carlo; several printed formulas in the
   SINMAP manual's appendix contain typesetting errors.
4. **Parameters.** Friction is given as a coefficient, tan φ, not an angle.
   Regions come from a raster and a CSV table; a cell whose region has no row
   is an error.
5. **Flat cells** are unconditionally stable (class 1), as in SINMAP, and
   nodata in the stability index and factor-of-safety rasters.
6. **Saturation** keeps SINMAP's coding, but writes the continuous wetness
   instead of splitting it into "partially wet" and "low moisture".
7. **Edge contamination** is optional and off by default.
8. **No interactive calibration** and, for now, no landslide-inventory
   statistics. The slope - specific catchment area boundaries are written as a
   table for plotting.

## Relation to the topographic wetness index

With the same slope and drainage area, SINMAP's wetness and the topographic
wetness index (TWI) written by `lsdtt-drainage-analysis` are the same index:
w = (R/T) × cell size × exp(TWI) / cos θ. The models differ in three ways,
each a deliberate choice in SINMAP:

- SINMAP uses sin θ, Darcy flow parallel to the slope; TOPMODEL uses tan θ,
  the small-slope form. They differ by cos θ: 3 % at a gradient of 0.25 and
  29 % at a gradient of 1.
- SINMAP assumes uniform hydraulic conductivity over a finite soil, so wetness
  rises linearly with a / sin θ until the soil is saturated. TOPMODEL's
  conductivity decays exponentially with depth, so its saturation deficit is
  linear in the logarithm of the index. Both rank cells in the same order.
- SINMAP's recharge is a steady rate for a critical wet period; TOPMODEL's
  deficit follows a catchment water balance through time.

The two rasters also differ through their inputs. `lsdtt-drainage-analysis`
computes TWI from a polynomial-fit slope (100 m window by default), whereas
this CLI uses the D-infinity facet slope on the carved-then-filled surface.

## Pipeline

1. Read the DEM, handle interior nodata, then carve (least-cost breach, by
   default) and fill it, as `lsdtt-slope-area` does.
2. Build a D-infinity flow graph on the routing surface: slope from the facet
   of steepest descent, and drainage area A. The specific catchment area is
   a = A / w, with w the contour width (see below).
3. Optionally read a region raster and a region table, giving one parameter
   set per region.
4. Compute, per cell, FS_min, FS_max, the probability that FS > 1, the
   stability index and class, and the saturation.
5. Write the area in each stability class, and the class boundaries of SINMAP's
   slope - specific catchment area plot.

## Usage

```text
lsdtt-slope-stability <path> <parameter_file>
lsdtt-slope-stability <parameter_file>
lsdtt-slope-stability -h
```

## Outputs

| File | Contents |
|---|---|
| `<write_fname>_sinmap_stability_index.tif` | SI: FS_min where FS_min > 1, P(FS > 1) where FS_min <= 1 <= FS_max, 0 where FS_max < 1 |
| `<write_fname>_sinmap_stability_class.tif` | 1 stable, 2 moderately stable, 3 quasi-stable, 4 lower threshold, 5 upper threshold, 6 defended |
| `<write_fname>_sinmap_probability_stable.tif` | P(FS > 1) |
| `<write_fname>_sinmap_fs_min.tif`, `_sinmap_fs_max.tif` | factor of safety at the most and least conservative parameter sets |
| `<write_fname>_sinmap_saturation.tif` | 3 always saturated, 2 saturated for part of the T/R range, otherwise the wettest possible wetness |
| `<write_fname>_specific_catchment_area.tif`, `_dinf_slope.tif` | a (m) and slope (m/m) |
| `<write_fname>_sinmap_class_areas.csv` | area by stability class and region |
| `<write_fname>_sinmap_slope_area_curves.csv` | class boundaries for the slope - specific catchment area plot |

Flat cells are class 1 (unconditionally stable) and nodata in the stability
index and factor-of-safety rasters.

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

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_before_routing` | `bool` | `true` | `switch` | - | Choose whether to pit-fill the DEM before building the flow network. | ASKS: should the DEM be filled at all? Filling is OPTIONAL here and this is a genuine analytical choice. True (the default) routes over a pit-free surface. False builds the flow network on the UNFILLED surface, so flow paths terminate at sinks instead of being carried through them. Do not confuse this with raster_is_filled, which is an assertion that the input is already filled. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope` | `float` | `0.0001` | `param` | `fill_before_routing` | Minimum slope imposed during pit filling. | Used when fill_before_routing is true. |
| `breach_before_fill` | `bool` | `true` | `switch` | - | Breach (carve) depressions before filling pits. | The carve comes first, then the fill (fill_before_routing) takes only the pits the carve could not solve, so large filled areas (and their artificially low slopes) are avoided. SINMAP itself only fills: false reproduces that. |
| `breach_algorithm_choice` | `string` | `least_cost` | `param` | `breach_before_fill` | Algorithm for depression breaching. | lindsay2016: Lindsay (2016) efficient breaching. least_cost: Lindsay (2020) pit-centric least-cost breaching. |
| `breach_lc_max_dist` | `int` | `1000 (outlet) / 100 (whitebox)` | `param` | `breach_before_fill` | Maximum breach path length in cells for least-cost breaching. | Counted in steps from the pit. Pits that need a longer path are left unbreached (and filled if fill_before_routing is true). 1 to 32767. If not set, 1000 with breach_lc_edge_policy outlet and 100 with whitebox. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_max_cost` | `float` | `unlimited` | `param` | `breach_before_fill` | Maximum cost of a least-cost breach path. | The summed depth each cell on the path is lowered (m; m x cells with breach_lc_minimize_dist). A pit whose cheapest path costs more is left unbreached. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_minimize_dist` | `bool` | `false` | `param` | `breach_before_fill` | Weight least-cost breach paths by their length. | false: a path costs the total depth lowered (m). true: each step is weighted by its length in cells. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_edge_policy` | `string` | `outlet` | `param` | `breach_before_fill` | Whether a least-cost breach may drain off the raster edge. | outlet: the raster edge is an outlet. whitebox: it is not (as in WhiteboxTools). Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_flat_increment` | `float` | `0.0` | `param` | `breach_before_fill` | Elevation drop (m) per cell along a least-cost breach channel. | 0 or less selects the automatic value (twice the float32 spacing at the largest absolute elevation). Only used when breach_algorithm_choice is least_cost. |
| `write_conditioned_pixels_raster` | `bool` | `false` | `switch` | - | Write the conditioned-pixels raster. | 0 untouched, 1 filled (the routing surface is above the input DEM), 2 carved. Their slope is the conditioning's, not the landscape's, which matters for stability. Output: &lt;write_fname&gt;_conditioned_pixels.tif. |

### Flow information

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four characters in north, east, south, west order. |

### Threading

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 uses the runtime default. Results do not depend on it. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

### Parameter regions

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `slope_stability_region_fname` | `string` | - | `param` | - | Raster of integer region ids, one parameter set per region. | Same grid as the DEM (for example a lithology raster). Cells whose id is nodata are left out (nodata in every output). Needs slope_stability_region_table_fname. Empty: the parameter-file values apply everywhere. A relative path is read from read_path. |
| `slope_stability_region_table_fname` | `string` | - | `param` | `slope_stability_region_fname` | CSV of parameters by region id. | Columns: region_id (integer, unique), optionally region_name, and any of the sinmap_* parameter names. A blank or absent column takes the parameter-file value. Any other column is an error. Every region id in the raster must have a row. A relative path is read from read_path. |

### SINMAP

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_sinmap` | `bool` | `true` | `switch` | - | Compute the SINMAP stability index (Pack et al. 2005). | Infinite slope, soil thickness perpendicular to the slope, steady-state wetness w = min(R a / (T sin theta), 1), with cohesion, friction coefficient and T/R uniform between bounds. SI = FS_min where FS_min &gt; 1; P(FS &gt; 1) where FS_min &lt;= 1 &lt;= FS_max; 0 where FS_max &lt; 1. Slope is the D-infinity facet slope (as SINMAP) and a = A / w (contour_width_method). Inspired by SINMAP, not a replica: see the CLI reference for the differences. For shallow translational failure only (SINMAP manual section 1.3). |
| `sinmap_dimensionless_cohesion_min` | `float` | `0.0` | `param` | `run_sinmap` | Lower bound of the dimensionless cohesion C = (C_root + C_soil) / (h rho_s g). | h is the soil thickness perpendicular to the slope. SINMAP default. |
| `sinmap_dimensionless_cohesion_max` | `float` | `0.25` | `param` | `run_sinmap` | Upper bound of the dimensionless cohesion. | SINMAP default. Equal bounds make cohesion certain. |
| `sinmap_friction_coefficient_min` | `float` | `0.5773503` | `param` | `run_sinmap` | Lower bound of the friction coefficient tan(phi) (dimensionless). | 0.577 = tan 30 degrees, SINMAP's default lower friction angle. |
| `sinmap_friction_coefficient_max` | `float` | `1.0` | `param` | `run_sinmap` | Upper bound of the friction coefficient tan(phi). | 1.0 = tan 45 degrees, SINMAP's default. |
| `sinmap_transmissivity_over_recharge_min_metres` | `float` | `2000` | `param` | `run_sinmap` | Lower bound of T/R (m): transmissivity over the effective recharge of a critical wet period. | The smaller T/R is the wetter bound. (T/R) sin(theta) is the length of planar hillslope needed to saturate. SINMAP default. |
| `sinmap_transmissivity_over_recharge_max_metres` | `float` | `3000` | `param` | `run_sinmap` | Upper bound of T/R (m). | SINMAP default. |
| `sinmap_water_to_soil_density_ratio` | `float` | `0.5` | `param` | `run_sinmap` | Water-to-soil density ratio r = rho_w / rho_s. | SINMAP treats it as constant: 1000 / 2000 kg/m3 by default. 0 to 1. |
| `contour_width_method` | `string` | `dinf_flow_direction` | `param` | `run_sinmap` | Contour width w in the specific catchment area a = A / w. | dinf_flow_direction (default): w = cell size x max(\|sin beta\|, \|cos beta\|) for each cell's D-infinity flow direction beta, the width of a D-infinity flow tube, which makes a exact on planes at any orientation. cell_size: w = cell size whatever the direction, as SINMAP does; on a D-infinity grid that underestimates a by up to 29 % where flow is diagonal (0.71 of the flow distance on a plane draining at 45 degrees). |
| `exclude_edge_contaminated_cells` | `bool` | `false` | `switch` | `run_sinmap` | Leave out cells whose contributing area reaches the raster edge or interior nodata. | Their drainage area may be truncated, so their wetness is a lower bound and their stability index an upper bound. Off by default: on a DEM clipped to a basin, the clip boundary counts as an edge and the whole basin would be excluded. Excluded cells are nodata in every output; the mask is written as &lt;write_fname&gt;_edge_contaminated.tif. |
| `write_sinmap_rasters` | `bool` | `true` | `switch` | `run_sinmap` | Write the SINMAP rasters. | &lt;write_fname&gt;_sinmap_stability_index.tif, _sinmap_stability_class.tif (1 stable, 2 moderately stable, 3 quasi-stable, 4 lower threshold, 5 upper threshold, 6 defended), _sinmap_probability_stable.tif, _sinmap_fs_min.tif, _sinmap_fs_max.tif, _sinmap_saturation.tif (3 always saturated, 2 saturated for part of the T/R range, else the wettest possible wetness in [0, 1)), _specific_catchment_area.tif (a, m) and _dinf_slope.tif (m/m). Flat cells are class 1 and nodata in the stability index and FS rasters. |
| `write_sinmap_class_area_table` | `bool` | `true` | `switch` | `run_sinmap` | Write the area in each stability class (by region). | &lt;write_fname&gt;_sinmap_class_areas.csv: region_id (0 without a region raster), stability_class, n_cells, area_m2, fraction_of_region. |
| `write_sinmap_slope_area_curves_table` | `bool` | `true` | `switch` | `run_sinmap` | Write the class boundaries of SINMAP's slope - specific catchment area plot. | &lt;write_fname&gt;_sinmap_slope_area_curves.csv: region_id, curve, slope (m/m), specific_catchment_area_m. Curves: always_saturated_above, never_saturated_below, fs_min_1.5, fs_min_1.0, fs_max_1.0 (manual Fig. 4, Appendix 1). Slope runs to the steepest cell, a from one cell size to the largest a. For plotting with the cell cloud. |

<!-- END GENERATED CLI PARAMETERS -->
