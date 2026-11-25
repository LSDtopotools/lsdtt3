# lsdtt-slope-stability

!!! warning "Experimental"
    This CLI is under development. Its components so far (a SINMAP-style
    stability index, planar rupture surfaces after scaLr, a Dupuit - Boussinesq
    water table, a probability of failure, 2-D slip surfaces by the method of
    slices and 3-D ellipsoids by the method of columns) are tested against
    analytic solutions. How to choose between them, what the results mean and
    how they were checked are in
    [Slope stability](../design_and_algorithms/slope_stability.md); the
    [Python guide](../python/slope_stability.md) shows the same steps in Python.

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

For deep-seated slides it can integrate the forces over the sliding body: by
the method of slices on scaLr's plane, circular arcs, compound (bilinear) or
structurally controlled surfaces, or in 3-D by the method of columns on
ellipsoids (see "Slip surfaces" below).

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
  `groundwater_steady_method` chooses how it is found:
  - `implicit` (the default) solves the steady equations directly by Newton's
    method. Each cell is either saturated (h = d, seeping), unsaturated
    (inflow balances outflow), or dry. It uses pseudo-transient continuation
    with an active-set Newton step, so it takes tens of iterations whatever the
    grid's stability limit.
  - `explicit` marches the time-dependent model until it stops changing. Its
    time steps are limited by the Courant and von Neumann conditions, so fine
    grids and conductive regolith need many.
  - Both solve the same discrete equations and agree to
    `groundwater_steady_tolerance`.
  - The implicit solve reports its budget as steady rates (m³/s: recharge, flow
    to channels, seepage). The explicit march and the transient run report a
    cumulative budget (m³).
- `dupuit_transient` starts from that steady state, runs a recharge series and
  evaluates the rupture planes at each output time. It writes a time series of
  storage, saturated cells and unstable cells.
  `groundwater_transient_method` chooses the stepper:
  - `explicit` (the default) takes substeps at the Courant and von Neumann
    limits. Its substep is set by stability, not accuracy, so on coarse grids
    it can be hours long.
  - `implicit` uses backward Euler. Each step is solved by the same Newton
    method as the implicit steady state. Steps are at most
    `groundwater_implicit_max_dt_s` and land on every recharge breakpoint and
    output time.
  - Backward Euler is first order in the step and damps fast transients, so the
    maximum step sets how well a storm is resolved.
    `groundwater_implicit_dt_tolerance_metres` adds step-doubling error control.
  - The explicit stepper regularises seepage within about 1 % of saturation,
    while the implicit one holds saturated cells exactly at the surface. Their
    counts of saturated cells therefore differ slightly near the seepage faces.
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
catalogue, which is what a mapper would see.

**Frequency – area.** For each catalogue the CLI also writes the
log-binned probability density p(A) (`frequency_area_bins_per_decade` bins per
decade) and fits its tail with a power law, p(A) ~ A^η.

- **The fit:** η = −α, where α is found by maximum likelihood. A_min is chosen
  to minimise the Kolmogorov–Smirnov distance between the data and the fitted
  tail (Clauset, Shalizi & Newman 2009, SIAM Review 51:661).
- **Discrete areas:** areas are whole numbers of cells, so the discrete form of
  the estimator is used.
- **Standard error:** (α − 1)/√n_tail. The tail must hold at least
  `frequency_area_min_tail` landslides.
- **Binned slope:** the least-squares slope of the binned density above A_min
  is written as well, for comparison with studies that fit binned densities. Volume is area × mean depth. The
landslide safety factor is scaLr's eq. 7, the summed resisting over the summed
driving stresses. Length and width are the axes of the ellipse with the same
second moments as the landslide's cells. These catalogues are statistics of
potential landslides; they are not maps of individual landslides.

## Slip surfaces for deep-seated landslides

scaLr judges a rupture plane by the infinite-slope formula at its crown R
alone. That formula sees only the depth there and the plane's angle, so bending
the surface would change nothing. With `slip_surface_family`, the forces are
instead summed over the whole body between R and the toe.

**2-D families (method of slices)** on a section along R's D8 path, unrolled
along the flow distance, with the ground linear between cell centres:

- `plane_integrated` keeps scaLr's plane, so it shows the effect of
  integration alone. On a planar slope its ordinary-method factor of safety is
  scaLr's formula with the cohesion term doubled, because the body's mean
  thickness is half the depth at R.
- `circular` uses arcs from depth h below R to the toe D, with a sagitta s (the
  depth of the arc below its chord); s = 0 is the plane. These represent
  rotational slides.
- `bilinear` (compound slides) has a steep back segment from depth h below R,
  at the Rankine active angle 45° + φ/2 (φ = atan μ), down to a flatter basal
  plane through D. A basal plane at the chord's gradient gives the plane.
- `structural` (dip-slope slides) follows a plane from depth h below R at the
  apparent dip of a bedding or weak layer: tan(dip) × cos(bearing of R to D −
  dip direction). The toe is wherever that plane meets the ground. A plane
  that cannot meet the ground before the channel is rejected.

A path that bends in plan stretches the section, so its sinuosity (flow
distance over straight distance) is written. Surfaces that rise above the
ground are rejected.

- **Back scarp** (`slip_surface_back_scarp`):
  - `vertical` (the default) runs from the ground at R straight down to the
    surface, with no resistance, like a tension crack.
  - `rankine` runs from the ground at R down at the Rankine angle until it
    meets the surface, and carries cohesion and friction.
  - The Rankine scarp also cuts off the wedge above it, so it can raise or
    lower F.
  - `bilinear` already has a Rankine back segment and takes only `vertical`.
- **Methods** (`slip_surface_method`):
  - ordinary (Fellenius): F = Σ[c l + μ (W cos α − u l)] / Σ W sin α;
  - simplified Bishop: F = Σ[(c b + μ (W − u b)) / m_α] / Σ W sin α, with
    m_α = cos α + μ sin α / F, iterated;
  - simplified Janbu: force equilibrium, Σ[(c b + μ (W − u b)) / (cos α m_α)]
    / Σ W tan α, times Janbu's correction factor
    f0 = 1 + b1 (d/L − 1.4 (d/L)²) (`slip_surface_janbu_correction`);
  - Spencer: force and moment equilibrium, with interslice forces at a constant
    inclination. It is slow when used to search.
  - Surfaces with m_α < 0.2 are rejected.
  - `auto` uses ordinary for planes, Bishop for circles, and Janbu for
    `bilinear` and `structural`.
  - `slip_surface_spencer_check` instead evaluates Spencer only on each
    critical surface found by the family's own method, as a check of the
    simplified methods. In the unit tests Spencer is within 0.4 % of Bishop on
    circles. On bilinear surfaces simplified Janbu is below Spencer, and its
    correction moves it closer.
- **Pore pressure** u = ρ_w g (z_w − z_b) is the *vertical* hydrostatic head of
  the water table above the base, as in the rupture planes and scaLr's eq. A1.
  It is not the slope-parallel seepage form, ρ_w g (z_w − z_b) cos²θ, often
  used for infinite slopes with flow parallel to the surface. The vertical form
  gives the larger pore pressure, by 1/cos²θ (about 6 % at a gradient of 0.25,
  and twice as much at a gradient of 1), and so the lower factor of safety.
- **Parameters:** c, μ and ρ come from the crown's region for the whole body.
- **Search** (`slip_surface_search`):
  - `seeded`, the fast default, fixes the toe at the rupture plane's daylight
    point. `full` searches every toe down the path and is never above the
    seeded result. For `structural` the toe follows from the plane.
  - Only crowns whose rupture-plane factor of safety is at most
    `slip_surface_screen_safety_factor` are evaluated. The rest are written as
    not evaluated, not as stable.
  - The shape is searched on a grid and refined with Brent's method. Planes
    belong to the `circular` and `bilinear` families, so their searches are
    never above `plane_integrated` with the same method.

**3-D family (method of columns): `ellipsoid`.** These are ellipsoids through
R and the rupture plane's daylight point D, after r.slope.stability (Mergili
et al. 2014).

- **Geometry:**
  - one half-axis lies along the chord R–D (half its length);
  - the width half-axis is horizontal across it, and the depth half-axis is
    normal to the chord;
  - the depth and width are searched;
  - the slip surface is the ellipsoid's lower surface, and each cell above it
    and below the ground is a column.
- **Factor of safety:**
  F = Σ f [c A + μ (G cos β_c − u A)] / Σ f G sin β_m, with no inter-column
  forces, where:
  - A is the column's slip-surface area;
  - β_c and β_m are the true dip, and the apparent dip in the sliding
    direction R to D;
  - G is the column's weight;
  - f = cos β_m for `xie_revised` (the default), or 1 for `hovland`
    (`slip_surface_column_method`).
- **Pore pressure** is the vertical head, as in 2-D.
- **Small ellipsoids:** crowns whose R–D half-length is below
  `slip_surface_ellipsoid_min_half_axis_pixels` cell sizes are not evaluated,
  because small ellipsoids are poorly resolved by the columns.
- **Layers:** with `slip_surface_layer_depth_fnames`, each ellipsoid is also
  evaluated truncated at each layer base, and the lowest factor of safety is
  kept. A layer therefore never raises F.
- **Search:** `slip_surface_ellipsoid_search: random` replaces the grid search
  with reproducible random draws, which also let the centre move up or down.
  Use it to check the seeded search.
- **Per-cell outputs:** `_slip_cover_safety_factor.tif` gives the lowest
  factor of safety of the critical ellipsoids covering each cell, and with
  uncertainty `_slip_cover_failure_probability.tif` the highest probability.

**Bodies:** an unstable crown's body covers its path cells (2-D) or its
columns (3-D).

- **Depth at a cell:**
  - 2-D: the mean thickness of the body over the cell's span along the section,
    so a thin wedge between cell centres still counts;
  - 3-D: the column's thickness.
- `_slip_body_depth.tif` is the deepest unstable body over each cell, so a
  volume is Σ depth × cell area.
- The by-toe catalogue groups unstable crowns sharing a toe. Its groups can
  overlap.
- The amalgamated catalogue takes the connected parts of the body-depth raster,
  which do not overlap.

**Probability of failure:** with `parameter_uncertainty`, Pf is computed on
each crown's critical surface, found with the central values. It is therefore
a lower bound on the probability that some surface fails. Screened-out crowns
get Pf = 0, and `_slip_evaluated.tif` marks them.

**Two dimensions are not three:** a section has no side resistance. For slides
of low length-to-width ratio, the 3-D factor of safety is typically higher than
the 2-D one; the `ellipsoid` family includes it.

## Probability of failure

With `parameter_uncertainty: equal_density`, cohesion and the friction
coefficient are uncertain, each uniform (between `_min` and `_max`) or normal.
A normal parameter has its mean at the central value and is truncated at
mean ± 3 SD and at 0.

- Each parameter takes `uncertainty_samples_per_parameter` values at equal
  intervals of its distribution's CDF, with trapezoid weights. This is the
  equal-density strategy of r.slope.stability (Mergili et al. 2014), with exact
  quantiles. Every combination of the two is used, and the weights multiply.
- For each cell and each combination, the rupture-plane search is repeated in
  full, so the daylight point may differ between combinations.
- The probability of failure is the weighted fraction of combinations with a
  factor of safety of at most 1. The weighted mean and SD of the factor of
  safety are written as well.
- Density and the water table are fixed. With `dupuit_transient`, the water
  table is the one at the last output time.
- The deterministic rasters and the landslide catalogues use the central
  values, `rupture_cohesion_pa` and `rupture_friction_coefficient`.
- A region table may set the bounds and SDs per region: the
  `rupture_cohesion_pa_min`, `_max` and `_sd` and the
  `rupture_friction_coefficient_min`, `_max` and `_sd` columns.
- The probability is a trapezoid rule on a step function, so its error is at
  most about 1/(n − 1) for n values per parameter. The default of 5 values
  gives 25 combinations and an error of up to about 0.25. Use more values where
  the probability itself matters.
- The cost is roughly the number of combinations times the deterministic
  search. The downslope path is walked once per cell.

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

## Potential landslides from a DEM of difference

With `run_dem_change_landslides: true` the tool finds candidate areas of change in
a DEM of difference, with their uncertainty. `read_fname` is the **newer** DEM on
the same grid. The method and its limits are in
[Slope stability](../design_and_algorithms/slope_stability.md#potential-landslides-from-a-dem-of-difference).

```
read_fname: newer_dem.tif
run_sinmap: false
run_rupture_planes: false
run_dem_change_landslides: true
dem_difference_fname: newer_minus_older.tif
dem_difference_stable_mask_fname: stable_terrain.tif
dem_change_sign: both
dem_change_min_slope: 0.1
```

Make the difference with
[`lsdtt-raster-calculator`](cli_reference_lsdtt_raster_calculator.md), with the
newer DEM as the first raster. If the difference is older minus newer, set
`dem_difference_is_newer_minus_older: false`.

| Output | Content |
|---|---|
| `_dem_change_uncertainty.csv` | Bias, NMAD, RMSE and the variogram fit |
| `_dem_change_uncertainty_by_slope.csv` | Bias and NMAD by slope |
| `_dem_change_variogram.csv` | Empirical and fitted semivariance |
| `_dem_change_patches.csv` | Each patch: area, shape, change, volume, uncertainties, significance |
| `_dem_change_coverage.csv` | Loss and gain: patches, area, share of the assessed area, volume |
| `_dem_change_sensitivity.csv` | Patches and area at other thresholds |
| `_dem_change_false_detection.csv` | Patches and area from noise alone |
| `_dem_change_false_detection_by_area.csv` | Patches found and from noise alone, by size |
| `_dem_change_frequency_area.csv`, `_fit.csv` | Binned p(A) and its power-law tail |
| `_dem_change_volume_area_fit.csv` | \|V\| = αA^γ for loss and gain, OLS and RMA (bootstrap intervals with `scaling_fit_bootstrap_samples`) |
| `_dem_change_depth_area_fit.csv` | \|mean change\| = αA^γ, as above |
| `_dem_change_frequency_volume.csv`, `_fit.csv` | Binned p(V) and its power-law tail (continuous maximum likelihood) |
| `_dem_change_corrected`, `_sigma`, `_z` | Bias-corrected change, error, and change / error |
| `_dem_change_class`, `_label` | -1 loss, 1 gain, 0 none; patch id |
| `_dem_change_exclusion.tif` | With `write_dem_change_exclusion_raster`: 1 excluded (masks and channels) |
| `_dem_change_patches.fgb` or `.gpkg` | With `write_dem_change_patch_polygons`: one feature per 4-connected part of a patch, with holes and the patch table's fields |

**Exclusions.** `dem_change_exclusion_mask_fname` takes one or more rasters,
comma separated (canopy from both surveys, water, buildings); their union is
removed from the detection and the error model. Channel beds aggrade and incise
without landsliding: `dem_change_channel_specific_catchment_metres` removes cells
whose D-infinity specific catchment area on the carved and filled DEM is at least
that, grown by `dem_change_channel_buffer_pixels`. Detection itself stays on the
input DEM. `dem_change_min_patch_mean_change_metres` drops shallow patches.

```
dem_change_exclusion_mask_fname: canopy_older.tif, canopy_newer.tif, water.tif
dem_change_channel_specific_catchment_metres: 5000
dem_change_channel_buffer_pixels: 2
write_dem_change_patch_polygons: true
scaling_fit_bootstrap_samples: 1000
```

Draw them with `lsdviz3-slope-stability` from lsdviztools3.

## Coseismic hazard area

`run_coseismic_hazard_area: true` computes the SHALRUN-EQ hazard area of Milledge et
al. (2019). It is the area of cells steep enough to initiate a landslide whose debris
can reach each cell by slope-weighted multiple flow directions, without crossing a
cell gentler than the stopping slope. It is a relative index of exposure.

```yaml
read_fname: my_dem
write_fname: my_site
run_sinmap: false
run_rupture_planes: false
run_coseismic_hazard_area: true
coseismic_hazard_initiation_slope: 0.839
coseismic_hazard_stopping_slope: 0.176
```

Both thresholds are gradients (m/m). The defaults come from 30 m DEMs and depend on
resolution. The method, its limits and its differences from the paper are in
[Slope stability](../design_and_algorithms/slope_stability.md#coseismic-hazard-area).

## Scar volumes from a landslide inventory

`run_scar_drape_volumes: true` estimates the volume of each mapped scar from a
post-event DEM and an inventory of scar polygons, when there is no pre-event
DEM. For each scar a cloth is held at the ground on a ring of cells just
outside the scar and is free inside it, but cannot pass below the ground. Its
resting shape estimates the pre-event surface, and the cloth minus the DEM is
the depth removed.

```yaml
read_fname: post_event_dem.tif
write_fname: my_event
run_sinmap: false
run_rupture_planes: false
run_scar_drape_volumes: true
landslide_inventory_fname: scars.gpkg
landslide_inventory_id_field_name: slide_id
scar_drape_rim_buffer_pixels: 1
```

The inventory can be any vector format GDAL reads (shapefile, GeoPackage,
FlatGeobuf), in any CRS. Draw each polygon around the source area only, out
to the crown, so that the ring outside it is undisturbed ground.

**Testing against a pre-event DEM.** Where a pre-event DEM exists, give it as
`scar_drape_reference_dem_fname` (on the same grid; `lsdtt-raster-calculator`
with `write_second_coregistered_raster: true` gives one). Each scar then also
gets the volume removed between the two DEMs, and the drape is repeated on the
pre-event DEM as a control. A perfect method gives zero there, so the control
measures concavity that existed before the event and that the cloth counts
wrongly. Each scar also gets rim descriptors: the relief of the rim, the
gradient of the plane through it, and how far the rim departs from that plane
(`rim_plane_residual_m`), which is the curvature a taut cloth cannot follow.

| Output | Content |
|---|---|
| `_scar_drape_volumes.csv` | Per scar: id, status, mapped area and volume, mean and maximum depth, area and volume including the buffer, centroid, cell counts (rim, rim off the DEM or nodata, nodata inside, cells where the cloth rests on the ground), flags and convergence |
| `_scar_drape_volume_area_fit.csv` | V = α A^γ over the scars with status `ok`, by OLS and by reduced major axis (RMA), with bootstrap 95 % intervals when `scaling_fit_bootstrap_samples` > 0; with a reference DEM, the same for the reference volumes |
| `_scar_drape_validation.csv` | With `scar_drape_reference_dem_fname`: totals, drape / reference ratio, median and quartiles of the per-scar ratio, correlation of log volumes, and the control as a fraction of the reference |
| `_scar_drape_surface.tif`, `_depth.tif`, `_label.tif` | Draped surface, depth (m) and scar label (position in the layer + 1); where scars overlap, the deeper cloth |

**What it gets right and wrong.** A taut cloth spans a scar linearly between
its rims. On a planar slope it restores the slope exactly. Across a convex
nose it lies too low, so the volume is too small. Across a hollow, or a scar cut
into an older gully, it lies too high, so the volume is too large. The rim
buffer trades two errors against each other. With no buffer, a polygon drawn
inside the crown puts the ring on disturbed ground and the volume is too small.
A wide buffer takes in more of the surrounding slope's curvature.
`scar_drape_sag_curvature_per_metre` lets the cloth sag, which lowers volumes.
The method and a test against repeat lidar are in
[Slope stability](../design_and_algorithms/slope_stability.md#scar-volumes-from-a-draped-cloth).

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
7. With `parameter_uncertainty`, repeat the search for each (c, μ)
   combination and write the probability of failure.
8. With `slip_surface_family`, find each screened crown's critical slip
   surface (method of slices in 2-D, of columns in 3-D), build the bodies of
   the unstable ones, and write their rasters and catalogues.
9. Group unstable cells into landslides by shared daylight point, and write
   the two catalogues.

With `run_coseismic_hazard_area`, the hazard area is computed on the routing
surface after step 1. `run_scar_drape_volumes` and `run_dem_change_landslides`
use the input DEM after step 1's nodata handling, before carving and filling.

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
| `<write_fname>_rupture_failure_probability.tif` | probability that the factor of safety is at most 1 over the (c, μ) combinations (`parameter_uncertainty`) |
| `<write_fname>_rupture_safety_factor_mean.tif`, `_rupture_safety_factor_sd.tif` | weighted mean and SD of the factor of safety over the combinations |
| `<write_fname>_slip_safety_factor.tif`, `_slip_unstable.tif` | factor of safety of the critical slip surface, and 1 where it is at most 1 (`slip_surface_family`) |
| `<write_fname>_slip_crown_depth.tif`, `_slip_max_depth.tif`, `_slip_sagitta.tif` | depth of the surface below the crown (ellipsoid: the depth half-axis), deepest body thickness, and sagitta (circular) (m) |
| `<write_fname>_slip_basal_gradient.tif` | bilinear: gradient of the basal plane; structural: the apparent dip (m/m) |
| `<write_fname>_slip_half_width.tif`, `_slip_cover_safety_factor.tif`, `_slip_cover_failure_probability.tif` | ellipsoid: the lateral half-axis (m); the lowest F and the highest Pf of the critical ellipsoids covering each cell |
| `<write_fname>_slip_spencer_safety_factor.tif` | Spencer's F on each critical 2-D surface (`slip_surface_spencer_check`) |
| `<write_fname>_slip_length.tif`, `_slip_sinuosity.tif` | flow distance to the toe (ellipsoid: the straight R–D distance) (m), and flow over straight distance (2-D only) |
| `<write_fname>_slip_evaluated.tif` | 1 evaluated, 0 screened out by `slip_surface_screen_safety_factor` |
| `<write_fname>_slip_body_depth.tif`, `_slip_body_label_amalgamated.tif` | deepest unstable body over each cell (m), and the amalgamated body id |
| `<write_fname>_slip_body_catalogue.csv`, `_slip_body_catalogue_amalgamated.csv` | unstable bodies by toe, and merged: crowns, cells, area, depth, volume, safety factors, centroid, length, width |
| `<write_fname>_slip_failure_probability.tif`, `_slip_safety_factor_mean.tif`, `_slip_safety_factor_sd.tif` | probability of failure on the fixed critical surface (`parameter_uncertainty`) |
| `<write_fname>_landslide_catalogue.csv` | one row per daylight point: area, depth, volume, safety factors, centroid, length, width, perimeter |
| `<write_fname>_landslide_catalogue_amalgamated.csv` | the same for touching landslides merged |
| `<write_fname>_landslide_frequency_area.csv`, `_landslide_frequency_area_fit.csv` | log-binned probability density p(A) of each catalogue, and its power-law tail p(A) ~ A^η: η by maximum likelihood with A_min chosen by the Kolmogorov–Smirnov distance (discrete in cells), its standard error, and the least-squares slope of the binned density above A_min |
| `<write_fname>_slip_body_frequency_area.csv`, `_slip_body_frequency_area_fit.csv` | the same for the slip-body catalogues |
| `<write_fname>_landslide_label.tif`, `_landslide_label_amalgamated.tif` | catalogue id of each landslide cell |
| `<write_fname>_hazard_area.tif` | coseismic hazard area A_h (m²/m); 0 where no source can reach |
| `<write_fname>_hazard_area_steepest_slope.tif` | steepest downslope gradient used by the thresholds (m/m) |
| `<write_fname>_hazard_area_zone.tif` | 0 transit, 1 source, 2 stopping |
| `<write_fname>_hazard_area_edge_influenced.tif` | 1 where debris from a DEM-edge cell arrives, so A_h may be a lower bound |
| `<write_fname>_hazard_area_summary.csv` | the thresholds, cell counts and source area |
| `<write_fname>_scar_drape_volumes.csv`, `_scar_drape_volume_area_fit.csv` | per-scar volumes from a draped cloth, and the volume - area fit |
| `<write_fname>_scar_drape_surface.tif`, `_scar_drape_depth.tif`, `_scar_drape_label.tif` | draped surface, depth removed (m) and scar label |
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

### DEM of difference

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_dem_change_landslides` | `bool` | `false` | `switch` | - | Find potential landslides in a DEM of difference, with their uncertainty. | Needs dem_difference_fname. read_fname is the newer DEM on the same grid, used for slope. The results are CANDIDATE areas of change: vegetation, snow, water, radar penetration and DEM artefacts also change elevation between two products. |
| `dem_difference_fname` | `string` | - | `param` | `run_dem_change_landslides` | DEM of difference (metres), for example from lsdtt-raster-calculator. | On the grid of read_fname. Read from read_path if the path is relative. |
| `dem_difference_is_newer_minus_older` | `bool` | `true` | `param` | `run_dem_change_landslides` | The difference is newer minus older, so a loss of elevation is negative. | false: older minus newer, which is converted. lsdtt-raster-calculator writes first minus second. |
| `dem_difference_stable_mask_fname` | `string` | - | `param` | `run_dem_change_landslides` | Raster of stable terrain (non-zero = stable), for the error model. | On the grid of read_fname. Without it every assessable cell within 4 NMAD of the median is treated as stable, which is reasonable only where change is rare. |
| `dem_change_sign` | `string` | `loss` | `param` | `run_dem_change_landslides` | Which change to detect: loss, gain or both. | loss finds source areas (scars). Gain finds deposits, and also vegetation growth. |
| `dem_change_detection_sigma` | `float` | `1.96` | `param` | `run_dem_change_landslides` | Multiple of the slope-dependent error at which a cell can belong to a patch. | - |
| `dem_change_seed_sigma` | `float` | `3.0` | `param` | `run_dem_change_landslides` | A connected group of candidate cells needs one cell at this multiple. | This hysteresis keeps weak cells attached to a strong signal and drops isolated weak groups. |
| `dem_change_patch_significance` | `float` | `1.96` | `param` | `run_dem_change_landslides` | Minimum \|mean change\| / its standard error for a patch to be kept. | The standard error allows for spatially correlated DEM error. |
| `dem_change_min_patch_pixels` | `int` | `5` | `param` | `run_dem_change_landslides` | Smallest patch, in pixels. | - |
| `dem_change_aggregation_radius_pixels` | `int` | `2` | `param` | `run_dem_change_landslides` | Join changed cells with up to this many unchanged cells between them. | 0 joins only cells that touch. Unchanged cells between them are not added to the patch's area or volume. Larger values merge neighbouring landslides. |
| `dem_change_min_slope` | `float` | `0.0` | `param` | `run_dem_change_landslides` | Ignore cells with a gradient (m/m) below this. | Useful for excluding valley floors, reservoirs and farmland. |
| `dem_change_n_slope_bins` | `int` | `10` | `param` | `run_dem_change_landslides` | Quantile bins of slope for the error model. | - |
| `dem_change_n_null_simulations` | `int` | `20` | `param` | `run_dem_change_landslides` | Noise simulations used to estimate false detections. | Each runs the detection on Gaussian noise with the fitted error and spatial correlation. Real DEM errors have heavier tails, so this is a lower bound. 0 turns it off. |
| `dem_change_random_seed` | `int` | `33` | `param` | `run_dem_change_landslides` | Seed for the variogram sample and the noise simulations. | - |
| `dem_change_exclusion_mask_fname` | `string` | - | `param` | `run_dem_change_landslides` | Rasters of cells to exclude (non-zero = excluded), comma separated. | Their union is used: canopy masks from both surveys, water, buildings. Excluded cells are never detected and are not in the error model. On the grid of read_fname; read from read_path if relative. |
| `dem_change_channel_specific_catchment_metres` | `float` | `0.0` | `param` | `run_dem_change_landslides` | Exclude channels: cells whose D-infinity specific catchment area (m) is at least this. | Computed on the carved and filled DEM (breach_before_fill, fill_before_routing); detection stays on the input DEM. Channel beds aggrade and incise, which is not landsliding. 0 turns it off. |
| `dem_change_channel_buffer_pixels` | `int` | `1` | `param` | `dem_change_channel_specific_catchment_metres` | Grow the channel exclusion by this many pixels (Euclidean). | - |
| `dem_change_min_patch_mean_change_metres` | `float` | `0.0` | `param` | `run_dem_change_landslides` | Drop patches whose \|mean change\| (m) is below this. | Applied to the magnitude; the sign is set by dem_change_sign. Also applied in the sensitivity runs and the noise simulations. |
| `write_dem_change_exclusion_raster` | `bool` | `false` | `switch` | `run_dem_change_landslides` | Write the union of the exclusions (masks and channels) as a raster. | &lt;write_fname&gt;_dem_change_exclusion.tif: 1 excluded, 0 not. |
| `write_dem_change_patch_polygons` | `bool` | `false` | `switch` | `run_dem_change_landslides` | Write the patches as polygons, with the fields of the patch table. | &lt;write_fname&gt;_dem_change_patches.fgb or .gpkg (write_vector_format). One feature per 4-connected part, with holes; a patch joined across unchanged cells has several parts, each carrying the patch's fields. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

### Coseismic hazard area

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_coseismic_hazard_area` | `bool` | `false` | `switch` | - | Compute the coseismic landslide hazard area (SHALRUN-EQ, Milledge et al. 2019). | The area of cells steep enough to initiate whose debris can reach each cell by slope-weighted multiple flow directions without crossing a cell gentler than the stopping slope. A relative index of exposure. Runs on the routing surface. Outputs: &lt;write_fname&gt;_hazard_area.tif (m2/m), _hazard_area_steepest_slope.tif, _hazard_area_zone.tif (0 transit, 1 source, 2 stopping), _hazard_area_edge_influenced.tif and _hazard_area_summary.csv. |
| `coseismic_hazard_initiation_slope` | `float` | `0.839` | `param` | `run_coseismic_hazard_area` | Gradient (m/m) at or above which a cell can initiate a landslide. | Default from the six-earthquake average of Milledge et al. (2019) on 30 m DEMs. An effective value that depends on DEM resolution, not a material property. |
| `coseismic_hazard_stopping_slope` | `float` | `0.176` | `param` | `run_coseismic_hazard_area` | Gradient (m/m) below which run-out stops. | A gentler cell is still reached and counted, but passes nothing on. Must be above 0 and below coseismic_hazard_initiation_slope. Default from Milledge et al. (2019) on 30 m DEMs; depends on DEM resolution. |

### Scar drape volumes

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_scar_drape_volumes` | `bool` | `false` | `switch` | - | Estimate the volume of each mapped landslide scar with a cloth draped over it. | read_fname is the POST-event DEM; landslide_inventory_fname holds the scar polygons. For each scar a cloth is held at the ground on a ring just outside the scar (scar_drape_rim_buffer_pixels) and is free inside it, but may not pass below the ground. Its resting shape estimates the pre-event surface; the volume is the cloth minus the DEM. A taut cloth interpolates linearly across the scar, so it is exact on a planar slope, too shallow across a convex nose and too deep across a hollow. Uses the input DEM (no carving or filling). Outputs: &lt;write_fname&gt;_scar_drape_volumes.csv, _scar_drape_volume_area_fit.csv and (write_scar_drape_rasters) _scar_drape_surface.tif, _scar_drape_depth.tif, _scar_drape_label.tif. |
| `landslide_inventory_fname` | `string` | - | `param` | `run_scar_drape_volumes` | Landslide inventory: polygons of the scars (shapefile, GeoPackage, FlatGeobuf, ...). | Any OGR vector format. Reprojected to the DEM's CRS. Each polygon is one scar; draw it around the source area only, not the runout, and to the crown, so that the ring outside it is undisturbed ground. Read from read_path if the path is relative. |
| `landslide_inventory_layer_name` | `string` | - | `param` | `run_scar_drape_volumes` | Layer of landslide_inventory_fname to read. | Empty: the first layer. |
| `landslide_inventory_id_field_name` | `string` | - | `param` | `run_scar_drape_volumes` | Attribute that identifies each scar in the outputs. | Empty: the feature's position in the layer (0-based). The run stops, listing the fields, if the layer has no such field. |
| `scar_drape_rim_buffer_pixels` | `int` | `1` | `param` | `run_scar_drape_volumes` | Cells added around each scar polygon before the cloth is held at the ground. | The cloth is free over the polygon grown by this many cells (8-connected) and held on the ring of cells just outside that. 1 tolerates a polygon drawn slightly inside the crown; larger values move the ring further onto stable ground but take in more of the surrounding slope. Volumes are reported over the polygon and over the free area. |
| `scar_drape_sag_curvature_per_metre` | `float` | `0.0` | `param` | `run_scar_drape_volumes` | Curvature (1/m) the cloth takes where it hangs clear of the ground. | 0: a taut cloth (harmonic surface between the rim and the ground it rests on). A positive value lets it sag: a span of L metres sags by about value x L^2 / 8, so 0.01 sags a 20 m scar by 0.5 m, giving smaller volumes. In physical units, so the result does not depend on the grid. |
| `scar_drape_tolerance_metres` | `float` | `1e-05` | `param` | `run_scar_drape_volumes` | Stop when the cloth moves less than this in one sweep (m). | Scars that reach scar_drape_max_iterations first are reported with status not_converged. |
| `scar_drape_max_iterations` | `int` | `100000` | `param` | `run_scar_drape_volumes` | Largest number of relaxation sweeps per scar. | Over-relaxation needs a few times the scar's width in cells. |
| `scar_drape_reference_dem_fname` | `string` | - | `param` | `run_scar_drape_volumes` | Pre-event DEM on the post-event DEM's grid, to validate the drape volumes. | When given, each scar also gets the volume removed between this DEM and read_fname over its mapped cells (net, loss and gain), and the same drape is run on this DEM as a control: a perfect method would give 0 there, so the control measures older concavity the cloth wrongly counts. Writes &lt;write_fname&gt;_scar_drape_validation.csv. Make the grids match with lsdtt-raster-calculator (write_second_coregistered_raster). Read from read_path if relative. |
| `scar_drape_min_reference_volume_m3` | `float` | `1.0` | `param` | `scar_drape_reference_dem_fname` | Smallest net reference volume (m3) for a scar to enter the per-scar ratios. | Scars with less, or with net gain, are left out of the ratio statistics but kept in the totals. |
| `write_scar_drape_rasters` | `bool` | `true` | `switch` | `run_scar_drape_volumes` | Write the draped surface, the depth and the scar labels as rasters. | &lt;write_fname&gt;_scar_drape_surface.tif (cloth elevation), _scar_drape_depth.tif (cloth minus DEM, m) and _scar_drape_label.tif (feature position + 1, 0 off scars). Where scars overlap the deeper cloth is kept. |

### Scaling fits

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `scaling_fit_bootstrap_samples` | `int` | `0` | `param` | - | Bootstrap resamples for the 95 % intervals of the scaling fits (0 = none). | Applies to the volume - area and depth - area fits, OLS and RMA, of the scar drape and of the DEM of difference. 1000 is a reasonable choice. |

### Rupture planes

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_rupture_planes` | `bool` | `true` | `switch` | - | Compute planar rupture surfaces after scaLr (Steer et al. 2026). | For each hillslope cell R, every cell D down its D8 path (stopping before the first channel cell) is a candidate daylight point; a plane from depth h below R to D gives an infinite-slope safety factor at R, minimised in closed form over h, and the D with the lowest factor is kept. Outputs the potential rupture depth and its safety factor per cell. Elevation differences come from the input DEM (after interior nodata handling); paths follow the routing surface. Pore pressure from water_table_method. |
| `rupture_cohesion_pa` | `float` | `1000` | `param` | `run_rupture_planes` | Cohesion c of the slope material (Pa). | scaLr's reference value is 1 kPa, deliberately low as an apparent strength during a triggering event (Steer et al. 2026, section 3). Rupture depth rises with cohesion. |
| `rupture_friction_coefficient` | `float` | `0.7002075` | `param` | `run_rupture_planes` | Friction coefficient mu = tan(phi) of the slope material (dimensionless). | 0.700 = tan 35 degrees, scaLr's reference value. |
| `rupture_density_kg_m3` | `float` | `2600` | `param` | `run_rupture_planes` | Bulk density rho of the slope material (kg/m3). | As in scaLr. Gravity is 9.81 m/s2. |
| `rupture_peak_ground_acceleration_fname` | `string` | - | `param` | `run_rupture_planes` | Peak ground acceleration raster, in g, for a pseudo-static earthquake load. | Must be on the DEM's grid; nodata cells are left out. Adds a horizontal force k W down the dip of each plane, k = rupture_seismic_coefficient_ratio x PGA. Empty means static. Cannot be combined with water_table_method or slip_surface_family. |
| `rupture_seismic_coefficient_ratio` | `float` | `0.5` | `param` | `rupture_peak_ground_acceleration_fname` | Pseudo-static coefficient k as a fraction of PGA (in g). | 0.5 follows Hynes-Griffin and Franklin (1984). It is an empirical choice, and F depends on it directly, so test a range. |
| `write_rupture_rasters` | `bool` | `true` | `switch` | `run_rupture_planes` | Write the rupture-plane rasters. | &lt;write_fname&gt;_rupture_safety_factor.tif (F at the daylight point), _rupture_depth.tif (optimal depth h, m), _rupture_unstable.tif (1 where F &lt;= 1), _rupture_length.tif (horizontal distance R to D, m), _rupture_slope.tif ((dz - h)/dx) and _rupture_topographic_slope.tif (dz/dx). Channel cells and cells with no admissible D are nodata. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Rupture uncertainty

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `parameter_uncertainty` | `string` | `none` | `param` | `run_rupture_planes` | Probability of failure of the rupture planes from uncertain cohesion and friction. | none: deterministic only. equal_density: uncertainty_samples_per_parameter values of cohesion and of the friction coefficient at equal intervals of each distribution's CDF, with trapezoid weights (the equal-density strategy of r.slope.stability, Mergili et al. 2014), all combinations. For each cell and each combination the rupture-plane search is repeated exactly (the daylight point may differ between combinations). Writes &lt;write_fname&gt;_rupture_failure_probability.tif (the weighted fraction of combinations with F &lt;= 1), _rupture_safety_factor_mean.tif and _rupture_safety_factor_sd.tif (weighted mean and SD of F). Density and the water table are fixed (with dupuit_transient, the water table at the last output time). The deterministic rasters and the landslide catalogues use the central values (rupture_cohesion_pa, rupture_friction_coefficient). Cost: about the number of combinations times the deterministic search. |
| `rupture_cohesion_distribution` | `string` | `uniform` | `param` | `parameter_uncertainty` | Distribution of cohesion for parameter_uncertainty. | uniform: between rupture_cohesion_pa_min and _max. normal: mean rupture_cohesion_pa, SD rupture_cohesion_pa_sd, truncated to mean +- 3 SD and at 0. Per region from the region table if it has these columns. |
| `rupture_cohesion_pa_min` | `float` | `500` | `param` | `rupture_cohesion_distribution=uniform` | Lower bound of uniform cohesion (Pa). | - |
| `rupture_cohesion_pa_max` | `float` | `1500` | `param` | `rupture_cohesion_distribution=uniform` | Upper bound of uniform cohesion (Pa). | - |
| `rupture_cohesion_pa_sd` | `float` | `250` | `param` | `rupture_cohesion_distribution=normal` | SD of normal cohesion (Pa); the mean is rupture_cohesion_pa. | 0 makes cohesion certain. |
| `rupture_friction_coefficient_distribution` | `string` | `uniform` | `param` | `parameter_uncertainty` | Distribution of the friction coefficient for parameter_uncertainty. | uniform: between rupture_friction_coefficient_min and _max. normal: mean rupture_friction_coefficient, SD rupture_friction_coefficient_sd, truncated to mean +- 3 SD and at 0. Per region from the region table if it has these columns. |
| `rupture_friction_coefficient_min` | `float` | `0.6` | `param` | `rupture_friction_coefficient_distribution=uniform` | Lower bound of the uniform friction coefficient (dimensionless). | - |
| `rupture_friction_coefficient_max` | `float` | `0.8` | `param` | `rupture_friction_coefficient_distribution=uniform` | Upper bound of the uniform friction coefficient (dimensionless). | - |
| `rupture_friction_coefficient_sd` | `float` | `0.05` | `param` | `rupture_friction_coefficient_distribution=normal` | SD of the normal friction coefficient; the mean is rupture_friction_coefficient. | 0 makes friction certain. |
| `uncertainty_samples_per_parameter` | `int` | `5` | `param` | `parameter_uncertainty` | Number of equal-density values per uncertain parameter. | 5 gives 25 combinations. The error in the probability of failure is at most about 1/(n - 1) (the trapezoid rule on a step); 1 uses the median only. |

### Slip surfaces

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `slip_surface_family` | `string` | `none` | `param` | `run_rupture_planes` | Slip surfaces integrated over the sliding body (deep-seated landslides). | none: off. On a 2-D section (R's D8 path, unrolled along the flow distance), by the method of slices: plane_integrated, scaLr's plane but with the forces summed over the body from the crown R to the toe D (isolates the effect of integration); circular, arcs from depth h below R to D with a sagitta s (rotational slides); bilinear, a steep back segment at the Rankine angle 45 deg + phi/2 from depth h below R down to a flatter basal plane through D (compound slides); structural, a plane from depth h below R at the apparent dip of a given bedding or weak layer, ending where it daylights (dip-slope slides). In 3-D: ellipsoid, ellipsoids through R and D by the method of columns (Hovland / Xie, as r.slope.stability), which include side resistance. A surface above the ground anywhere is rejected. Needs run_rupture_planes (the seeds and screen). The 2-D families have no side resistance, so for slides of low length-to-width ratio the 3-D factor of safety is typically higher. |
| `slip_surface_method` | `string` | `auto` | `param` | `slip_surface_family` | Limit-equilibrium method for the slip surfaces. | ordinary (Fellenius): F = sum[c l + mu (W cos a - u l)] / sum W sin a. bishop (simplified): F = sum[(c b + mu (W - u b)) / m_a] / sum W sin a, m_a = cos a + mu sin a / F, iterated; surfaces with m_a &lt; 0.2 are rejected. janbu (simplified, force equilibrium, with Janbu's correction factor f0). spencer (force and moment equilibrium, constant interslice inclination; slow, for checking). auto: ordinary for plane_integrated, bishop for circular, janbu for bilinear and structural. The ellipsoid family uses slip_surface_column_method instead. |
| `slip_surface_search` | `string` | `seeded` | `param` | `slip_surface_family` | Which toes are searched. | seeded (fast): the toe is the rupture plane's daylight point D* (plus or minus slip_surface_toe_window_pixels along the path). full: every toe down the path (up to slip_surface_max_length_metres); never above the seeded factor of safety, at roughly path-length times the cost. |
| `slip_surface_spencer_check` | `bool` | `false` | `switch` | `slip_surface_family` | Also write Spencer's factor of safety on each critical 2-D surface. | &lt;write_fname&gt;_slip_spencer_safety_factor.tif: Spencer's method (force and moment equilibrium) evaluated on the surface the family's own method found, as a check of the simplified methods. Much cheaper than slip_surface_method: spencer, which searches with it. Nodata where Spencer has no solution. |
| `slip_surface_back_scarp` | `string` | `vertical` | `param` | `slip_surface_family` | The back scarp at the crown of a 2-D slip surface. | vertical: from the ground at R straight down to the surface, with no resistance (a tension crack). rankine: from the ground at R down at the Rankine active angle 45 deg + phi/2 (phi = atan of the crown's friction coefficient) until it meets the surface; the scarp carries cohesion and friction. Not for bilinear, whose back segment is already at that angle. |
| `slip_surface_janbu_correction` | `bool` | `true` | `switch` | `slip_surface_method=janbu` | Apply Janbu's correction factor f0 = 1 + b1 (d/L - 1.4 (d/L)^2). | d is the deepest point of the surface below the line from its upper end to the toe, L that line's length; b1 = 0.50 (c and phi), 0.31 (phi only), 0.69 (c only). |
| `slip_surface_structural_dip_gradient` | `float` | `0.3` | `param` | `slip_surface_family=structural` | structural: the dip of the bedding or weak layer, as a gradient tan(dip) (m/m). | The apparent dip along each crown's section is this times cos(bearing of R to D* - dip direction). |
| `slip_surface_structural_dip_direction_degrees` | `float` | `0` | `param` | `slip_surface_family=structural` | structural: the direction the layer dips towards, degrees clockwise from north. | - |
| `slip_surface_structural_dip_gradient_fname` | `string` | - | `param` | `slip_surface_family=structural` | structural: a raster of dip gradient, in place of the constant; read at each crown. | Same grid as the DEM. |
| `slip_surface_structural_dip_direction_fname` | `string` | - | `param` | `slip_surface_family=structural` | structural: a raster of dip direction (degrees clockwise from north), in place of the constant. | Same grid as the DEM. |
| `slip_surface_column_method` | `string` | `xie_revised` | `param` | `slip_surface_family=ellipsoid` | ellipsoid: the column sum. | F = sum f [c A + mu (G cos beta_c - u A)] / sum f G sin beta_m over the columns, with f = cos beta_m (xie_revised) or 1 (hovland). No inter-column forces. |
| `slip_surface_ellipsoid_search` | `string` | `seeded` | `param` | `slip_surface_family=ellipsoid` | ellipsoid: how the depth and width half-axes are searched. | seeded: a slip_surface_search_points grid over the depth half-axis c (up to half the R-D distance a) and the lateral half-axis b (up to 2a), refined with Brent; the ellipsoid passes through R and D*. random: slip_surface_ellipsoid_random_draws draws, also letting the centre move up or down by up to c/2 (reproducible from slip_surface_ellipsoid_random_seed). |
| `slip_surface_ellipsoid_random_draws` | `int` | `200` | `param` | `slip_surface_ellipsoid_search=random` | ellipsoid, random search: draws per crown. | - |
| `slip_surface_ellipsoid_random_seed` | `int` | `1` | `param` | `slip_surface_ellipsoid_search=random` | ellipsoid, random search: the seed. | - |
| `slip_surface_ellipsoid_min_half_axis_pixels` | `float` | `3` | `param` | `slip_surface_family=ellipsoid` | ellipsoid: the smallest half-axis, in cell sizes. | Smaller ellipsoids are poorly resolved by the columns (Mergili et al. 2014); crowns whose R-D half-length is below this are not evaluated. |
| `slip_surface_layer_depth_fnames` | `string` | - | `param` | `slip_surface_family=ellipsoid` | ellipsoid: comma-separated rasters of layer depth below the ground (m). | Each ellipsoid is also evaluated truncated at each layer (the slip surface follows the layer base where the ellipsoid goes deeper), and the lowest F kept. |
| `slip_surface_screen_safety_factor` | `float` | `1.5` | `param` | `slip_surface_family` | Only crowns whose rupture-plane factor of safety is at most this are evaluated. | Others are written as not evaluated (_slip_evaluated.tif = 0, nodata factor of safety), not as stable. Set it very large to evaluate every crown. |
| `slip_surface_toe_window_pixels` | `int` | `0` | `param` | `slip_surface_search=seeded` | Seeded search: toe candidates either side of D* along the path. | - |
| `slip_surface_max_length_metres` | `float` | `0` | `param` | `slip_surface_search=full` | Full search: the longest flow distance from crown to toe (m); 0 = no cap. | - |
| `slip_surface_slices_per_cell` | `int` | `4` | `param` | `slip_surface_family` | Slices per path segment between cell centres. | The sums converge to the arc integrals at second order in the slice width. |
| `slip_surface_search_points` | `int` | `12` | `param` | `slip_surface_family` | Grid points per shape parameter (h, and s / chord for circles) before Brent refinement. | - |
| `write_slip_surface_rasters` | `bool` | `true` | `switch` | `slip_surface_family` | Write the slip-surface rasters. | &lt;write_fname&gt;_slip_safety_factor.tif, _slip_crown_depth.tif (h, m), _slip_max_depth.tif (deepest body thickness at a cell centre, m), _slip_sagitta.tif (circular, m), _slip_length.tif (flow distance to the toe, m), _slip_sinuosity.tif (flow over straight distance), _slip_unstable.tif, _slip_evaluated.tif (1 evaluated, 0 screened out) and _slip_body_depth.tif (the deepest unstable body over each cell, m). Also _slip_basal_gradient.tif (bilinear: the basal plane; structural: the apparent dip). ellipsoid: _slip_crown_depth.tif is the depth half-axis, _slip_half_width.tif the lateral half-axis, _slip_length.tif the R-D distance, and _slip_cover_safety_factor.tif the lowest factor of safety of the critical ellipsoids covering each cell; no sinuosity. |
| `write_slip_body_catalogue` | `bool` | `true` | `switch` | `slip_surface_family` | Write the catalogues of unstable slip bodies as CSV. | &lt;write_fname&gt;_slip_body_catalogue.csv: unstable crowns sharing a toe; the footprint is the union of their bodies, the depth at a cell the deepest of them, volume = sum(depth) x cell area (groups can overlap). _slip_body_catalogue_amalgamated.csv: 8-connected components of the body depth raster, which do not overlap. Statistics of potential landslides, not maps of individual ones. |

### Landslide catalogues

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `rupture_min_daylight_neighbours` | `int` | `2` | `param` | `run_rupture_planes` | Fewest neighbours of a daylight point that must share it, or its landslide is dropped. | Removes "1-D landslides", single lines of cells down one flow path, as the scaLr code does (its iopt = 2; the paper does not describe it). The 8 neighbours of the daylight point D are counted if their daylight point is D, stable or not. 0 keeps every landslide. 0 to 8. |
| `write_landslide_catalogue` | `bool` | `true` | `switch` | `run_rupture_planes` | Write the potential-landslide catalogues as CSV. | A landslide is the unstable cells (F &lt;= 1) that share a daylight point, as in scaLr (every daylight point also joins its own landslide if its own F &lt;= 1, as the scaLr code does). &lt;write_fname&gt;_landslide_catalogue.csv: one row per daylight point; _landslide_catalogue_amalgamated.csv: touching landslides merged (8-connected), as a mapper would see them. Columns: area, mean and max depth, volume (area x mean depth), mean F, the landslide safety factor (scaLr eq. 7), centroid, length and width (axes of the ellipse with the same second moments), perimeter. These are statistics of potential landslides, not maps of individual ones (Steer et al. 2026, section 4.1). |
| `frequency_area_bins_per_decade` | `int` | `10` | `param` | `write_landslide_catalogue` | Logarithmic bins per decade of area for the frequency - area tables. | With each catalogue, &lt;write_fname&gt;_landslide_frequency_area.csv (and _slip_body_frequency_area.csv) holds the log-binned probability density p(A) per m2, and _frequency_area_fit.csv the power-law tail p(A) ~ A^eta: eta by maximum likelihood with A_min chosen by the Kolmogorov-Smirnov distance (Clauset et al. 2009, discrete in cells), with its standard error, and the least-squares slope of the binned density above A_min. |
| `frequency_area_min_tail` | `int` | `50` | `param` | `write_landslide_catalogue` | Fewest landslides in the power-law tail of the frequency - area fit. | - |
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
| `groundwater_steady_method` | `string` | `implicit` | `param` | `water_table_method` | How the steady water table is found (dupuit_steady, and the antecedent state of dupuit_transient). | implicit: solved directly by Newton's method. Pseudo-transient continuation with a semismooth (active-set) Newton step for the saturated (h = d) and dry (h = 0) cells; sparse direct or BiCGSTAB solves. Few iterations, independent of the stability limit. explicit: marched in time with the explicit scheme until it stops changing; substeps are limited by the Courant and von Neumann conditions, so fine grids and conductive regolith need many. Both solve the same discrete equations and agree to the tolerance. |
| `groundwater_steady_tolerance` | `float` | `1e-06` | `param` | `water_table_method` | Steady state reached when the residual is at most this x max recharge. | explicit: max \|n dh/dt\|. implicit: max \|net inflow\| on unsaturated cells (and its sign on saturated or dry ones). Both relative to the recharge rate. |
| `groundwater_max_substeps` | `int` | `10000000` | `param` | `groundwater_steady_method=explicit` | explicit steady method: most substeps in the march to steady state before giving up. | Substeps are limited by Courant and von Neumann conditions; fine grids and conductive regolith need many. |
| `groundwater_transient_method` | `string` | `explicit` | `param` | `water_table_method=dupuit_transient` | How dupuit_transient steps the recharge series: explicit or implicit. | explicit: substeps at the Courant and von Neumann limits, the semi-analytic thickness update with regularised seepage. implicit: backward Euler, each step solved by the same Newton method as the implicit steady state, with steps of at most groundwater_implicit_max_dt_s landing on every recharge breakpoint and output time. Backward Euler is first order in the step and damps fast transients, so the maximum step sets how well a storm is resolved. The explicit substep is limited by stability, not accuracy: on coarse grids it can be hours long. |
| `groundwater_implicit_max_dt_s` | `float` | `3600` | `param` | `groundwater_transient_method=implicit` | implicit transient method: longest backward-Euler step (s). | - |
| `groundwater_implicit_dt_tolerance_metres` | `float` | `0` | `param` | `groundwater_transient_method=implicit` | implicit transient method: step-doubling error tolerance on the water table (m); 0 = off. | Each step is also taken as two half steps; where they differ by more than this the step is halved, otherwise the half-step result is kept. The largest accepted difference is reported. |
| `groundwater_max_newton_iterations` | `int` | `200` | `param` | `groundwater_steady_method=implicit` | implicit steady method: most Newton (pseudo-time) steps before giving up. | - |
| `write_groundwater_rasters` | `bool` | `true` | `switch` | `water_table_method` | Write the groundwater rasters. | &lt;write_fname&gt;_water_table_depth.tif (m below the surface), _saturated_thickness.tif (m above the base) and _seepage.tif (return flow, m/s). |

<!-- END GENERATED CLI PARAMETERS -->
