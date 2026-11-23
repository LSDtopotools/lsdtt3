# Slope stability

`lsdtt-slope-stability` and its Python routines estimate where slopes are
prone to failure, as maps and statistics. They are **not** tools for mapping
individual landslides. Each model answers a different question, and they
share one pipeline: a DEM, its routing surface, a channel network, a water
table and material parameters.

This page explains what each model assumes, how to choose between them, what
the results mean, and how they were checked. The parameters, outputs and
equations in full are in the
[CLI reference](../cli_reference/cli_reference_lsdtt_slope_stability.md);
run `lsdtt-slope-stability -h` for the same list. The
[Python guide](../python/slope_stability.md) shows the same steps in Python.

## Which model answers which question

| Model | Question | Failure it represents | Main output |
|---|---|---|---|
| SINMAP stability index (Pack et al. 1998) | Which cells are stable, and how confidently, given ranges of parameters? | Shallow, translational, on an infinite slope | Stability index and its six classes; probability of stability |
| Rupture planes (scaLr; Steer et al. 2026) | How deep and how stable is the most dangerous plane from each hillslope cell to the ground below it? | Planar slides from a cell to a daylight point on its flow path | Factor of safety F, rupture depth, potential-landslide catalogues |
| Slip surfaces, 2-D (method of slices) | Does integrating the forces over the sliding body, or bending the surface, change the answer? | Deep-seated slides: planar, rotational, compound, dip-slope | F of the critical surface per crown, body depth, body catalogues |
| Ellipsoids, 3-D (method of columns; after Mergili et al. 2014) | How much does side resistance add? | Slides of low length-to-width ratio | F per crown and the lowest F covering each cell |
| Dupuit groundwater (after Litwin et al. 2020) | Where is the water table, and how does a storm raise it? | (Supplies pore pressure to the others) | Water-table depth, saturated thickness, seepage |
| Probability of failure | How likely is failure given uncertain cohesion and friction? | (Applies to the rupture planes, slip surfaces and ellipsoids) | P(F ≤ 1), with the mean and SD of F |

The models build on each other. The rupture planes seed the slip surfaces and
ellipsoids, and screen which crowns they evaluate. The groundwater model
gives the rupture planes, slip surfaces and ellipsoids their pore pressure.
SINMAP keeps its own steady wetness. Slopes are gradients (m/m) throughout,
and friction is a coefficient μ = tan φ, never an angle.

## The models in brief

**SINMAP.** The infinite-slope factor of safety with a steady wetness
w = min(1, R a / (T sin θ)), where a is the specific catchment area, evaluated
over uniform ranges of dimensionless cohesion, the friction coefficient and
T/R. The stability index is FS_min where that exceeds 1, the probability that
FS > 1 where the range straddles 1, and 0 where even FS_max is below 1.
lsdtt3 follows SINMAP's model and classes but differs in eight places, all
listed in the CLI reference. The largest is the contour width: dividing the
D-infinity area by the cell size, as SINMAP does, underestimates the specific
catchment area by up to 29 % where flow is diagonal, so lsdtt3 divides by the
width of the D-infinity flow tube instead.

**Rupture planes.** From each hillslope cell R, every cell D on its D8 path to
the channel is a candidate daylight point. A plane from depth h below R to D
has a factor of safety F(h) with a closed-form minimum h\*. The candidate with
the lowest F is kept. The result describes potential failure statistically,
as in scaLr.

**Slip surfaces.** scaLr judges a plane by the infinite-slope formula at the
crown alone, so bending the surface would change nothing. The slip surfaces
instead sum the forces over the body between the crown and the toe, on a
section along the crown's D8 path, by the ordinary, Bishop, Janbu or Spencer
method.

**Ellipsoids.** The same idea in 3-D: ellipsoids through R and D, with the
forces summed over columns (Hovland, or the revised method of Xie et al.).

**Groundwater.** An unconfined aquifer on an impermeable base at the regolith
depth. Water flows down the water-table gradient (the Dupuit–Boussinesq
equation), channel cells are seepage outlets with the water table at the
surface, and water seeps out where the regolith saturates.

## Choosing a slip-surface family

- **Start with the rupture planes.** They are fast, they give every hillslope
  cell a depth and F, and the slip surfaces need them as seeds.
- **`plane_integrated`** keeps scaLr's plane but integrates the forces over
  the body. Use it to see the effect of integration alone. On a planar slope
  its ordinary-method F is scaLr's with the cohesion term doubled, because the
  body's mean thickness is half the depth at the crown.
- **`circular`** (Bishop by default) represents rotational slides in roughly
  homogeneous material.
- **`bilinear`** (Janbu by default) represents compound slides: a steep back
  segment at the Rankine active angle onto a flatter basal plane, as over a
  weak layer.
- **`structural`** (Janbu by default) follows the apparent dip of a bedding
  plane or weak layer, for dip-slope slides. Give the dip as a gradient and the
  dip direction as a bearing, constant or as rasters.
- **`ellipsoid`** adds side resistance. Use it where slides are short and wide,
  where a 2-D section overstates instability.
- **Methods.** `auto` picks ordinary for planes, Bishop for circles and Janbu
  (with its correction factor) for compound and structural surfaces. Spencer
  satisfies force and moment equilibrium but is about 400 times slower to
  search with; `slip_surface_spencer_check` evaluates it only on each critical
  surface, as a check. In the unit tests Spencer was within 0.33 % of Bishop
  on random circles, and simplified Janbu was 0.78–0.94 of Spencer on bilinear
  surfaces (0.84–0.97 with its correction).
- **Search.** `seeded` fixes the toe at the rupture plane's daylight point;
  `full` tries every toe down the path, costs roughly the path length times
  more, and is never above the seeded F. Use `full` to check the seeded
  results on a sample.

## What the results mean

- **Potential, not mapped, landslides.** The catalogues group unstable cells
  into potential landslides (by shared daylight point, then merged where they
  touch). They are statistics of where failure is possible, not an inventory.
  Their frequency – area exponent, area – volume scaling and shapes depend on
  the DEM, its resolution, the channel threshold and the 1-D rule
  (`rupture_min_daylight_neighbours`), so compare them with other studies only
  on matched choices.
- **Screened is not stable.** The slip surfaces and ellipsoids are evaluated
  only at crowns whose rupture-plane F is at most
  `slip_surface_screen_safety_factor`. The rest are marked "not evaluated"
  (`_slip_evaluated.tif`, `SlipSurfaces.evaluated`), not stable.
- **The slip-surface probability of failure is a lower bound.** It is
  evaluated on each crown's critical surface found with the central parameter
  values, held fixed. Another surface could fail under other values. The
  rupture-plane probability has no such limit: its search is repeated in full
  for every parameter combination.
- **The probability's resolution.** It is a trapezoid rule over equal-density
  samples of each parameter, so its error is at most about 1/(n − 1) for n
  samples per parameter: up to about 0.25 at the default of 5. Use more samples
  where the probability itself matters.
- **Two dimensions are not three.** A section has no side resistance, so for
  slides of low length-to-width ratio the 3-D factor of safety is typically
  higher than the 2-D one.
- **One material per body.** Cohesion, friction and density are the crown's
  region values for the whole body, not each slice's own.
- **Pore pressure is the vertical head.** On the rupture planes, slip surfaces
  and ellipsoids, u = ρ_w g (z_w − z_b), as in scaLr's eq. A1. The
  slope-parallel seepage form used for many infinite-slope analyses,
  ρ_w g (z_w − z_b) cos²θ, is smaller by cos²θ: about 6 % at a gradient of
  0.25, and half at a gradient of 1. The vertical head therefore gives the
  lower factor of safety. The slope-parallel form is not offered.

## The water table

- **Steady state.** `groundwater_steady_method: implicit` (the default) solves
  the steady equations directly. Each cell is saturated (the water table at
  the surface, seeping), unsaturated (inflow balancing outflow) or dry, and
  Newton's method with an active-set step finds which, in tens of iterations
  whatever the grid. `explicit` marches the time-dependent model until it stops
  changing; its steps are limited by the Courant and von Neumann conditions, so
  fine grids and conductive regolith need many. Both solve the same discrete
  equations and agree to the tolerance. On coarse grids the explicit march can
  be the faster of the two.
- **Storms.** `dupuit_transient` starts from the steady state and runs a
  recharge series. The default explicit stepper's substep is set by stability,
  not accuracy, and on coarse grids it can be hours long; in one test it
  crossed two days in 7 substeps and ended 3.5 cm from the converged answer.
  The implicit stepper (backward Euler) takes steps of at most
  `groundwater_implicit_max_dt_s`, landing on every breakpoint and output time;
  it is first order in the step and damps fast transients. Step-doubling error
  control (`groundwater_implicit_dt_tolerance_metres`) bounds its error per
  step.
- **Seepage.** The explicit scheme regularises seepage within about 1 % of
  saturation; the implicit one holds saturated cells exactly at the surface.
  Their counts of saturated cells differ slightly near seepage faces.
- **Transient probability.** In transient runs, the probability of failure and
  the slip surfaces use the water table at the last output time; the rupture
  F is written at every output time.

## Relation to TOPMODEL and the wetness index

With the same slope and drainage area, SINMAP's unsaturated wetness and the
topographic wetness index (TWI) written by `lsdtt-drainage-analysis` are the
same index: w = (R/T) b exp(TWI) / cos θ, with TWI = ln((A/b²)/tan θ) and b
the cell size. The models differ in three deliberate ways:

- **sin θ versus tan θ.** SINMAP's flux is T sin θ, with the soil thickness
  measured perpendicular to the slope. That is Darcy flow parallel to the
  slope, and it equals the Dupuit flux. TOPMODEL's T tan β is the small-slope
  form. tan β exceeds sin β by 1/cos β: 3 % at a gradient of 0.25, 41 % at a
  gradient of 1. lsdtt3's TWI keeps ln(a / tan β), for continuity with
  published parameter sets and earlier LSDTopoTools; see
  [Hydrological modelling](hydrological_modelling.md).
- **Transmissivity.** SINMAP has a uniform conductivity over a finite soil, so
  w rises linearly with a / sin θ until the soil saturates. TOPMODEL's
  conductivity decays exponentially with depth, so its deficit is linear in
  the logarithm of the index. Both rank cells in the same order. The Dupuit
  model offers both (`groundwater_transmissivity_profile`); with the
  exponential profile its steady state follows TOPMODEL's index relation with
  sin β cos β in place of tan β.
- **Forcing.** SINMAP uses a steady recharge for a critical wet period;
  TOPMODEL a catchment water balance through time.

In practice the two rasters differ more through their inputs than through the
model. `lsdtt-drainage-analysis` computes TWI from a polynomial-fit slope over
a 100 m window on a filled surface, while this tool uses the D-infinity facet
slope on a carved-then-filled surface. On a 30 m test DEM the Spearman
correlation between the indices was 0.9995 on shared inputs but 0.69 between
the two rasters as each tool writes them.

## Corrections to the published equations

- **scaLr eq. A2 is missing a factor 1/h.** Developing the paper's own eq. A1,
  with pore pressure ψ = ρ_w g (h − d) below a water table at depth d, gives

      F_wet(h) = F(h) − μ (ρ_w/ρ) ((h − d)/h) (Δx² + (Δz − h)²) / (Δx (Δz − h)),   h > d.

  A2 is printed with (h − d) in place of (h − d)/h. As printed, the
  pore-pressure term has units of metres while F is dimensionless, so it is
  dimensionally inconsistent. lsdtt3 uses the form above. F_wet has a kink at
  h = d, so its minimum is found separately on each side of the water table:
  in closed form above it, by Brent's method below.
- **scaLr's optimal depth.** The paper gives two roots (its eqs. 5 and 6) and
  notes that at high cohesion eq. 5 is the positive one. It is positive there,
  but it always lies deeper than the daylight point, so eq. 6 gives the optimal
  depth in every case. As printed, both divide by zero when c Δz = μ ρ g Δx²;
  lsdtt3 uses a rearranged form of eq. 6 that does not.
- **SINMAP manual.** The printed probability formulas A-19, A-33 and A-38
  have typesetting losses (the friction bounds printed as (t1, t1), and r
  dropped from the unsaturated term). lsdtt3 re-derives them, and the unit
  test checks them against Monte Carlo.

## Potential landslides from a DEM of difference

`run_dem_change_landslides` looks for clusters of elevation change in a DEM of
difference from
[`lsdtt-raster-calculator`](../cli_reference/cli_reference_lsdtt_raster_calculator.md).
Its results are **candidates**. Forestry, vegetation growth, snow, water, radar
penetration, layover and DEM artefacts also change elevation between two
products, so a patch is an area that changed by more than the DEM error, not a
mapped landslide.

**Error model.** On stable terrain the tool measures the bias (median), the
NMAD in quantile bins of slope, and a variogram of the standardised error, fitted
with a nugget and one spherical model (Rolstad et al. 2009; Hugonnet et al. 2022).
Cells more than 4 NMAD from the median are left out, so real change does not
inflate the error. Without a stable-terrain mask every assessable cell is treated
as stable, which holds only where change is rare. In a dense landslide cluster,
give a mask.

**Detection.** A cell is a candidate if its bias-corrected change exceeds
`dem_change_detection_sigma` times the error at its slope. A connected group is
kept only if one cell exceeds `dem_change_seed_sigma`. Changed cells with up to
`dem_change_aggregation_radius_pixels` unchanged cells between them form one
patch. The unchanged cells are not added to its area or volume.

**Patch uncertainty.** Errors are correlated, so averaging over a patch reduces
them less than for independent cells. The variance of a patch mean is

$$
\sigma_{\bar h}^2 = \overline{\sigma^2}\left[ n \frac{a}{A} + (1 - n)\, C(r, L) \right]
$$

with $n$ the nugget fraction, $a$ the cell area, $A$ the patch area, and $C$ the
mean spherical correlation (range $L$) between two points in a disc of area $A$.
A patch is kept if its mean change is at least `dem_change_patch_significance`
standard errors. Area uncertainty is half a cell along the perimeter. Regional
volume uncertainty treats patches as independent, so it is a lower bound where
patches are closer together than the correlation range.

**False detections.** The detection is rerun on simulated noise with the fitted
error, nugget and range. This gives the number and area of patches that noise
alone produces, overall and by size. Small patches are usually dominated by noise,
so read a frequency - area rollover and any power-law fit only above the size at
which noise becomes negligible. The simulated noise is Gaussian and real DEM
errors have heavier tails, so this is a **lower bound** on false detections. It
does not include systematic change such as forestry.

Run the same analysis on a nearby area with no known landslide event. Patches
found there show how much of the signal comes from land-cover change and sensor
differences.

## Design decisions

These are the choices that shape the results:

- **No landslide outlines.** The catalogues are statistics, with label rasters
  for checking; mapping individual landslides is out of scope.
- **Friction is a coefficient**, μ = tan φ, everywhere. The only angle among
  the parameters is the structural dip direction, because it is an azimuth.
- **Channels** are cells with at least `threshold_contributing_pixels`
  contributing pixels, as in the other tools. In the groundwater model they are
  seepage outlets with the water table at the surface.
- **Paths are D8 flow paths** on the routing surface, measured on the input
  DEM. A 2-D section is the crown's D8 path unrolled along the flow distance;
  its sinuosity is written, because a path that bends in plan stretches the
  section.
- **1-D landslides are removed** by default: a daylight point is kept only if
  at least `rupture_min_daylight_neighbours` (2) of its 8 neighbours share it.
  The rule is not applied to slip bodies, which are one cell wide by
  construction.
- **Curved and 3-D searches are seeded** from the rupture planes (toe, depth
  and screen) to save time; the full and random searches remain as checks.
- **Body depth** is the deepest unstable body over each cell (the envelope), so
  a body volume is the sum of depth times cell area.
- **Bilinear surfaces** use a fixed Rankine angle, 45° + φ/2 from the crown's
  μ, not a searched one. The basal plane passes through the toe and is no
  steeper than the crown–toe chord.
- **Structural surfaces** end where the plane meets the ground; the apparent
  dip uses the bearing of the straight line from crown to daylight point.
- **Ellipsoids** have one half-axis along the crown–daylight chord, a
  horizontal width and a depth normal to the chord; the centre lies on the
  chord in the seeded search; the sliding direction is the crown-to-daylight
  bearing, not the terrain aspect. Crowns shorter than
  `slip_surface_ellipsoid_min_half_axis_pixels` (3) are not evaluated. Column
  pore pressure is u A on the slip surface with the vertical head, so 2-D and
  3-D share one convention and F stays linear in (c, μ). The cover rasters
  take each evaluated crown's critical ellipsoid only.
- **Uncertainty** uses equal-density sampling of cohesion and friction
  (uniform or truncated normal), with density and the water table fixed.

## How it was checked

No original code (ArcMap, GRASS or MATLAB) is run. The components are tested
against analytic solutions on analytic surfaces, in the CTest compartment
`ctest -L slope_stability -LE heavy`:

- **SINMAP** (`unit_slope_stability_sinmap`): FS_min and FS_max against
  SINMAP's eqs. 3 and 8; the probability of stability against a 10⁶-draw Monte
  Carlo; the D-infinity specific catchment area on planes and cones, exact with
  the contour width at every flow angle; the identity with TWI.
- **Rupture planes** (`unit_slope_stability_rupture`): the closed-form depth
  against the minimum of F (300 random cases); the cohesionless limit; the
  wet minimum against a 2×10⁵-point search; a water table below the plane
  giving the dry result exactly; no interior local maximum of F_wet over a
  wide sweep of parameters; a brute-force search on a curved profile; thread
  independence.
- **Probability** (`unit_slope_stability_rupture_probability`,
  `unit_equal_density_sampling`): Pf against the exact probability on a plane,
  within the trapezoid bound; a single sample reproduces the deterministic run.
- **Slip surfaces** (`unit_slope_stability_slip_surface`,
  `unit_slope_stability_slip_compound`): the infinite-slope limit of a long
  slope-parallel base for every method; the plane-integrated F against its
  closed form; the slice sums converging to the arc integrals at second order;
  Spencer within 0.33 % of Bishop on 40 random circles; the Rankine wedge,
  bilinear and structural surfaces against their closed forms; Janbu's
  correction; the seeded search never below the full one.
- **Ellipsoids** (`unit_slope_stability_columns`): the infinite-slope limit
  of the column sums (Hovland and Xie, dry and wet); the ellipsoid volume
  converging to its closed form; truncation by layers never raising F; F
  linear in (c, μ) on a fixed ellipsoid.
- **Catalogues** (`unit_slope_stability_catalogue`, `unit_power_law_tail_fit`):
  the 1-D rule on constructed cases; area, depth, volume, axes, perimeter and
  F_L of a 3 × 7 rectangle; the power-law exponent recovered from synthetic
  samples.
- **Groundwater** (`unit_slope_stability_groundwater_dupuit`): the Dupuit
  mound and divide, the kinematic limit on a steep plane (equal to SINMAP's
  wetness), TOPMODEL's relation with the exponential profile, the decay of a
  linearised mound (and backward Euler's exact amplification), the water
  budget to round-off, implicit against explicit, and thread independence.
- **The CLI**: smoke tests of every family and option, and refusals of invalid
  parameter files.

Results on real landscapes depend on the DEM and on parameters that are rarely
known well. Treat them as relative measures of where failure is more likely,
unless the parameters have been calibrated locally.

## References

Hugonnet, R., et al. (2022). Uncertainty analysis of digital elevation models by
spatial inference from stable terrain. *IEEE JSTARS* 15, 6456-6472.

Rolstad, C., Haug, T. and Denby, B. (2009). Spatially integrated geodetic glacier
mass balance and its uncertainty based on geostatistical analysis. *Journal of
Glaciology* 55(192), 666-680.

- Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). Power-law
  distributions in empirical data. *SIAM Review*, 51(4), 661–703.
- Hovland, H. J. (1977). Three-dimensional slope stability analysis method.
  *Journal of the Geotechnical Engineering Division, ASCE*, 103(9), 971–986.
- Litwin, D. G., Tucker, G. E., Barnhart, K. R., & Harman, C. J. (2020).
  GroundwaterDupuitPercolator: a Landlab component for groundwater flow.
  *Journal of Open Source Software*, 5(46), 1935.
  https://doi.org/10.21105/joss.01935
- Mergili, M., Marchesini, I., Rossi, M., Guzzetti, F., & Fellin, W. (2014).
  Spatially distributed three-dimensional slope stability modelling in a raster
  GIS. *Geomorphology*, 206, 178–195.
- Montgomery, D. R., & Dietrich, W. E. (1994). A physically based model for the
  topographic control on shallow landsliding. *Water Resources Research*, 30(4),
  1153–1171.
- Pack, R. T., Tarboton, D. G., & Goodwin, C. N. (1998). The SINMAP approach to
  terrain stability mapping. *8th Congress of the International Association of
  Engineering Geology*, Vancouver.
- Pack, R. T., Tarboton, D. G., Goodwin, C. N., & Prasad, A. (2005). SINMAP 2:
  a stability index approach to terrain stability hazard mapping, user's
  manual. Utah State University.
- Steer, P., Pelascini, L., & Longuevergne, L. (2026). A physics-based
  morphometric model to explain the emergence of landslide rupture geometry and
  to identify hillslope transience. *EGUsphere* preprint,
  https://doi.org/10.5194/egusphere-2026-4142.
- Xie, M., Esaki, T., Zhou, G., & Mitani, Y. (2003). Three-dimensional
  stability evaluation of landslides and a sliding process simulation using a
  new geographic information systems component. *Environmental Geology*, 43,
  503–512.
