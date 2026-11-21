# lsdtt-slope-stability

!!! warning "Experimental"
    This CLI is under development. Its first two components, a SINMAP-style
    stability index and planar rupture surfaces after scaLr, are tested against
    analytic solutions; the user guide will follow once the remaining components
    are built and measured.

## Intent

`lsdtt-slope-stability` estimates where slopes are prone to failure, as maps
and statistics. It is not a tool for mapping individual landslides.

This version computes a stability index inspired by SINMAP (Pack, Tarboton,
Goodwin & Prasad 2005), which itself built on SHALSTAB (Montgomery & Dietrich
1994). It is an infinite-slope factor of safety with steady-state wetness, in
which cohesion, the friction coefficient and the ratio of transmissivity to
recharge are uncertain between bounds. Like SINMAP it applies to shallow
translational failure, not to deep-seated slides.

It also computes **planar rupture surfaces** after scaLr (Steer, Pelascini &
Longuevergne 2026). For each hillslope cell R, every cell D on its D8 path
down to the channel is a candidate point where the rupture surface reaches the
ground (the daylight point). A plane from depth h below R to D gives an
infinite-slope factor of safety at R. The depth minimising it has a closed form,
and the D with the lowest factor of safety is kept. The result is a potential
rupture depth and factor of safety for every hillslope cell. As in scaLr, these
rasters describe potential failure statistically; they do not map individual
landslides.

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

## Notes on the rupture planes

- The optimal depth is computed from a rearranged form of the scaLr paper's
  root (its eq. 6) that has no singular denominator. The paper notes that at
  very high cohesion the other root (its eq. 5) is the positive one. That root
  is indeed positive there, but it always lies deeper than the daylight point,
  so eq. 6 gives the optimal depth in every case. As printed, both roots divide
  by zero when c Δz = μ ρ g Δx²; the form used here does not.
- With zero cohesion the optimal depth is 0 and the factor of safety is
  μ Δx / Δz, the cohesionless infinite slope along the chord.
- Elevation differences come from the input DEM, after interior nodata is
  handled; paths follow the carved-then-filled routing surface. A candidate
  daylight point higher than its source cell is skipped.
- Channel cells are those with at least `threshold_contributing_pixels`
  contributing pixels. The walk stops before the first channel cell, and
  channel cells are not sources.
- Pore pressure is hydrostatic below a water table (`water_table_method`, see
  below). **Correction to the scaLr paper:** its printed eq. A2 is missing a
  factor 1/h. Developing its own eq. A1, with pore pressure
  ψ = ρ_w g (h − d) for h > d, gives

      F_wet(h) = F(h) − μ (ρ_w/ρ) ((h − d)/h) (Δx² + (Δz − h)²) / (Δx (Δz − h)),

  whereas A2 is printed with (h − d) in place of (h − d)/h. As printed, the
  pore-pressure term has units of metres while F is dimensionless, so A2 is
  dimensionally inconsistent. lsdtt3 uses the form derived from A1. With water, the minimum is no longer closed-form: it is
  found numerically on each side of the water table.

## Groundwater

The water table for pore pressure comes from a Dupuit - Boussinesq model of an
unconfined aquifer. It follows the Landlab GroundwaterDupuitPercolator (Litwin
et al. 2020): the aquifer lies on an impermeable base `regolith_thickness_metres`
below the surface, water flows down the water-table gradient, and water returns
to the surface where the regolith saturates.

- Channel cells are seepage outlets, with the water table at the surface. The
  raster edge and nodata cells are no-flux.
- `dupuit_steady` gives the steady state under a constant recharge.
- `dupuit_transient` starts from that steady state, runs a recharge series and
  evaluates the rupture planes at each output time. It writes a time series of
  storage, saturated cells and unstable cells.
- With `groundwater_transmissivity_profile: exponential`, conductivity decays
  with depth as TOPMODEL assumes. The steady state then follows TOPMODEL's
  index relation, with sin β cos β in place of tan β.
- With uniform conductivity and a steep, thin aquifer, the steady saturation
  equals SINMAP's steady-state wetness.

## Potential-landslide catalogues

A landslide is the set of unstable cells (factor of safety at most 1) that
share a daylight point. Like scaLr's code, lsdtt3 does two things the paper
does not describe:

- Every daylight point also joins its own landslide when its own factor of
  safety is at most 1.
- With `rupture_min_daylight_neighbours` (default 2), a daylight point is kept
  only if at least that many of its 8 neighbours share it, stable or not. This
  removes single lines of cells down one flow path.

Touching landslides are then merged (8-connected) into an amalgamated
catalogue, which is what a mapper would see. Volume is area × mean depth. The
landslide safety factor is scaLr's eq. 7, the summed resisting over the summed
driving stresses. Length and width are the axes of the ellipse with the same
second moments as the landslide's cells. These catalogues are statistics of
potential landslides; they are not maps of individual landslides.

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
6. Build a D8 flow graph on the routing surface. Mark channel cells with
   `threshold_contributing_pixels`, and find each hillslope cell's optimal
   planar rupture surface.
7. Group unstable cells into landslides by shared daylight point, and write
   the two catalogues.

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
| `<write_fname>_rupture_safety_factor.tif` | factor of safety of the optimal planar rupture surface |
| `<write_fname>_rupture_depth.tif` | its depth below the cell (m) |
| `<write_fname>_rupture_unstable.tif` | 1 where the factor of safety is at most 1, else 0 |
| `<write_fname>_rupture_length.tif` | horizontal distance to the daylight point (m) |
| `<write_fname>_rupture_slope.tif`, `_rupture_topographic_slope.tif` | gradient of the rupture plane and of the line from the cell to the daylight point (m/m) |
| `<write_fname>_landslide_catalogue.csv` | one row per daylight point: area, depth, volume, safety factors, centroid, length, width, perimeter |
| `<write_fname>_landslide_catalogue_amalgamated.csv` | the same for touching landslides merged |
| `<write_fname>_landslide_label.tif`, `_landslide_label_amalgamated.tif` | catalogue id of each landslide cell |
| `<write_fname>_water_table_depth.tif`, `_saturated_thickness.tif`, `_seepage.tif` | the groundwater state used for pore pressure |
| `<write_fname>_t<time>s_water_table_depth.tif`, `_t<time>s_rupture_safety_factor.tif` | the same at each transient output time |
| `<write_fname>_groundwater_time_series.csv` | storage, saturated cells and unstable cells at each output time |

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

### Rupture planes

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_rupture_planes` | `bool` | `true` | `switch` | - | Compute planar rupture surfaces after scaLr (Steer et al. 2026). | For each hillslope cell R, every cell D down its D8 path (stopping before the first channel cell) is a candidate daylight point; a plane from depth h below R to D gives an infinite-slope safety factor at R, minimised in closed form over h, and the D with the lowest factor is kept. Outputs the potential rupture depth and its safety factor per cell. Elevation differences come from the input DEM (after interior nodata handling); paths follow the routing surface. Dry: no pore pressure yet. |
| `rupture_cohesion_pa` | `float` | `1000` | `param` | `run_rupture_planes` | Cohesion c of the slope material (Pa). | scaLr's reference value is 1 kPa, deliberately low as an apparent strength during a triggering event (Steer et al. 2026, section 3). Rupture depth rises with cohesion. |
| `rupture_friction_coefficient` | `float` | `0.7002075` | `param` | `run_rupture_planes` | Friction coefficient mu = tan(phi) of the slope material (dimensionless). | 0.700 = tan 35 degrees, scaLr's reference value. |
| `rupture_density_kg_m3` | `float` | `2600` | `param` | `run_rupture_planes` | Bulk density rho of the slope material (kg/m3). | As in scaLr. Gravity is 9.81 m/s2. |
| `write_rupture_rasters` | `bool` | `true` | `switch` | `run_rupture_planes` | Write the rupture-plane rasters. | &lt;write_fname&gt;_rupture_safety_factor.tif (F at the daylight point), _rupture_depth.tif (optimal depth h, m), _rupture_unstable.tif (1 where F &lt;= 1), _rupture_length.tif (horizontal distance R to D, m), _rupture_slope.tif ((dz - h)/dx) and _rupture_topographic_slope.tif (dz/dx). Channel cells and cells with no admissible D are nodata. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Landslide catalogues

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `rupture_min_daylight_neighbours` | `int` | `2` | `param` | `run_rupture_planes` | Fewest neighbours of a daylight point that must share it, or its landslide is dropped. | Removes "1-D landslides", single lines of cells down one flow path, as the scaLr code does (its iopt = 2; the paper does not describe it). The 8 neighbours of the daylight point D are counted if their daylight point is D, stable or not. 0 keeps every landslide. 0 to 8. |
| `write_landslide_catalogue` | `bool` | `true` | `switch` | `run_rupture_planes` | Write the potential-landslide catalogues as CSV. | A landslide is the unstable cells (F &lt;= 1) that share a daylight point, as in scaLr (every daylight point also joins its own landslide if its own F &lt;= 1, as the scaLr code does). &lt;write_fname&gt;_landslide_catalogue.csv: one row per daylight point; _landslide_catalogue_amalgamated.csv: touching landslides merged (8-connected), as a mapper would see them. Columns: area, mean and max depth, volume (area x mean depth), mean F, the landslide safety factor (scaLr eq. 7), centroid, length and width (axes of the ellipse with the same second moments), perimeter. These are statistics of potential landslides, not maps of individual ones (Steer et al. 2026, section 4.1). |
| `write_landslide_label_rasters` | `bool` | `true` | `switch` | `run_rupture_planes` | Write the landslide id rasters behind the catalogues. | &lt;write_fname&gt;_landslide_label.tif and _landslide_label_amalgamated.tif: the catalogue id of each landslide cell, nodata elsewhere. For checking the statistics. |

### Groundwater

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `water_table_method` | `string` | `none` | `param` | - | Where the water table for pore pressure in the rupture planes comes from. | none: dry. dupuit_steady: the steady state of a Dupuit - Boussinesq aquifer on an impermeable base (regolith_thickness_metres below the surface) under groundwater_recharge_m_per_s, after the Landlab GroundwaterDupuitPercolator (Litwin et al. 2020); channel cells (threshold_contributing_pixels) are seepage outlets with the water table at the surface, the raster edge and nodata are no-flux. dupuit_transient: start from that steady state and run the recharge series in groundwater_recharge_series_fname, evaluating the rupture planes at each of groundwater_output_times_s. Pore pressure is hydrostatic below the water table (scaLr Appendix A, corrected). SINMAP keeps its own steady-state wetness. |
| `regolith_thickness_metres` | `float` | `2.0` | `param` | `water_table_method` | Thickness (m) of the permeable regolith above the impermeable base. | Measured vertically. Used where regolith_thickness_fname is empty. |
| `regolith_thickness_fname` | `string` | - | `param` | `water_table_method` | Raster of regolith thickness (m), same grid as the DEM. | Cells with nodata or a thickness &lt;= 0 are outside the aquifer. A relative path is read from read_path. |
| `groundwater_hydraulic_conductivity_m_per_s` | `float` | `0.001` | `param` | `water_table_method` | Saturated hydraulic conductivity K (m/s); at the surface for the exponential profile. | Landlab's default. Sands are about 1e-5 to 1e-3, silts 1e-8 to 1e-6. |
| `groundwater_drainable_porosity` | `float` | `0.2` | `param` | `water_table_method` | Drainable porosity n (0 to 1]. | Landlab's default. |
| `groundwater_transmissivity_profile` | `string` | `uniform` | `param` | `water_table_method` | How conductivity varies with depth: uniform or exponential. | uniform: T = K h (Landlab). exponential: K = K0 exp(-depth / z_m), as TOPMODEL assumes; with it the steady state follows TOPMODEL's index relation with sin cos in place of tan. |
| `groundwater_conductivity_decay_depth_metres` | `float` | `1.0` | `param` | `groundwater_transmissivity_profile=exponential` | Decay depth z_m (m) of the exponential conductivity profile. | TOPMODEL's m = drainable porosity x z_m. |
| `groundwater_recharge_m_per_s` | `float` | `1e-08` | `param` | `water_table_method` | Recharge f (m/s) for the steady state (dupuit_steady), and the antecedent state of dupuit_transient. | 1e-8 m/s is about 0.9 mm/day. For a design storm use the effective recharge of the critical wet period. |
| `groundwater_recharge_series_fname` | `string` | - | `param` | `water_table_method=dupuit_transient` | CSV recharge series for dupuit_transient: columns time_s, recharge_m_per_s. | Piecewise constant: each rate applies from its time to the next; the first time must be 0, times increasing; the last rate continues to the last output time. A relative path is read from read_path. |
| `groundwater_output_times_s` | `string` | - | `param` | `water_table_method=dupuit_transient` | Comma-separated times (s) at which dupuit_transient writes the water table and rupture safety factor. | &lt;write_fname&gt;_t&lt;time&gt;s_water_table_depth.tif and, with run_rupture_planes, _t&lt;time&gt;s_rupture_safety_factor.tif, plus one row each in &lt;write_fname&gt;_groundwater_time_series.csv. The full outputs are written for the last time. |
| `groundwater_steady_tolerance` | `float` | `1e-06` | `param` | `water_table_method` | Steady state reached when max \|n dh/dt\| &lt;= this x max recharge. | Relative to the recharge rate. |
| `groundwater_max_substeps` | `int` | `10000000` | `param` | `water_table_method` | Most substeps in the march to steady state before giving up. | Substeps are limited by Courant and von Neumann conditions; fine grids and conductive regolith need many. |
| `write_groundwater_rasters` | `bool` | `true` | `switch` | `water_table_method` | Write the groundwater rasters. | &lt;write_fname&gt;_water_table_depth.tif (m below the surface), _saturated_thickness.tif (m above the base) and _seepage.tif (return flow, m/s). |

<!-- END GENERATED CLI PARAMETERS -->
