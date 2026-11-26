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
| Coseismic hazard area (SHALRUN-EQ; Milledge et al. 2019) | Which cells could debris from steep ground reach during an earthquake? | Initiation and run-out, from slope thresholds alone | Hazard area $A_h$ (m²/m) |
| Scar drape volumes | How much material left each mapped scar, given only a post-event DEM? | (Measures failures that happened) | Volume, mean and maximum depth per scar; V = αA^γ |
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
as in scaLr. An optional [pseudo-static earthquake load](#pseudo-static-shaking-on-the-rupture-planes)
adds shaking from a peak ground acceleration raster.

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

**Exclusions.** Cells in any exclusion mask, and channel cells, are neither
detected nor used for the error model, exactly like cells gentler than
`dem_change_min_slope`. Channels are cells whose D-infinity specific catchment
area $a = A/b$ ($b$ the cell size) on the carved and filled DEM reaches
`dem_change_channel_specific_catchment_metres`, grown by a Euclidean buffer. The
carved surface is used only to route flow; the change is measured on the input
DEM. Without this, aggradation and incision along channels form long patches
that are not landslides.

**Polygons and scaling.** Patches are written as polygons along cell edges, one
feature per 4-connected part, with holes. A patch joined across unchanged cells
has several parts, and each carries the patch's fields, so summing `area_m2`
over a patch's parts gives its area. For loss and gain separately the tool fits
$|V| = \alpha A^\gamma$ and $|\bar h| = \alpha A^\gamma$ by ordinary least
squares and reduced major axis on the logarithms (both are reported because
detection truncates small, shallow patches and areas carry error, which biases
OLS low; see `src/utils/lsdtt_scaling_fit.hpp`), and the frequency - volume tail
$p(V) \propto V^{\eta}$ by the continuous maximum likelihood of Clauset et al.
(2009), with $V_{min}$ chosen by the Kolmogorov - Smirnov distance. Areas are
whole numbers of cells, so the frequency - area fit is discrete; volumes are not.
A DEM-of-difference inventory gives the observed scaling against which the
models' potential-landslide catalogues are compared.

**A mapped inventory against the DoD.** With `landslide_inventory_fname`, each
mapped polygon is matched to the detected patches through shared cells
(`match_inventory_to_dem_change`). Per polygon: the area detected as loss and
gain, the DoD's net, loss and gain volume over all its cells (detected or not),
and the patches it touches; it is matched if it touches a loss patch. Per
patch: its cells inside polygons and the polygons it overlaps. The summary
gives two area fractions: of the mapped area (union of polygons) detected as
loss, and of the detected loss inside the polygons. Overlapping polygons each
count the shared cells; the union counts them once. On East San Lorenzo, 39 of
81 mapped source polygons touch a loss patch, 24 % of the mapped area is
detected and 0.7 % of the detected loss is mapped: at an NMAD of 0.33 m most
shallow scars are below detection, and most detected loss is noise, channels
not removed, or unmapped change.

Run the same analysis on a nearby area with no known landslide event. Patches
found there show how much of the signal comes from land-cover change and sensor
differences.

## Pseudo-static shaking on the rupture planes

`rupture_peak_ground_acceleration_fname` adds an earthquake load to the rupture
planes. This is an addition to scaLr. A horizontal force $kW$ acts on each column,
directed down the dip of the plane, which is its least stable direction, with

$$
k = r \, \mathrm{PGA}
$$

where PGA is in g and $r$ is `rupture_seismic_coefficient_ratio` (default 0.5,
after Hynes-Griffin and Franklin 1984). With $s = \Delta z - h$ and
$C = c / \rho g$,

$$
F_k(h) = \frac{C (\Delta x^2 + s^2) / (h \Delta x) + \mu \max(\Delta x - k s, 0)}{s + k \Delta x}.
$$

Friction is zero where the inertial force removes the normal stress ($k s > \Delta x$).
Away from that limit, $dF_k/dh = 0$ is a quadratic, because the cubic terms cancel. Its
smaller positive root, with $E = \Delta z + k \Delta x$ and $L^2 = \Delta x^2 + \Delta z^2$, is

$$
h^* = \frac{C L^2 E}{C L^2 + \Delta x \sqrt{(1 + k^2) \, C L^2 (C + \mu E)}}.
$$

This is the static optimum at $k = 0$. Unlike the static case, $F_k$ stays finite as
$h \to \Delta z$, so lsdtt3 takes the lowest of $h^*$, that limit and, where $h^*$
falls in the unloaded part, the frictionless optimum and the point $k s = \Delta x$.

- **The coefficient is empirical.** $F$ depends directly on $r$, so run a range of values.
- **Treat $F$ as a ranking under strong shaking.** A pseudo-static load applied everywhere
  can bring most hillslope cells to $F \le 1$, so the unstable mask alone says little.
- **The load is steady and one-directional.** It has no duration, topographic
  amplification, or dependence on the direction of shaking.
- **Not combined with a water table or slip surfaces.** Those runs are refused.

## Coseismic hazard area

`run_coseismic_hazard_area` computes the hazard area of Milledge et al. (2019), the
SHALRUN-EQ index. It needs only a DEM, and it answers a different question from a
factor of safety: how much steep ground lies upslope of a cell with a viable
run-out path to it.

A cell $i$ is a source if its steepest downslope gradient $S_i$ is at least the
initiation gradient $S_m$. Debris moves by multiple flow directions, with each
receiver weighted by its gradient (Quinn et al. 1991). A cell with $S < S_s$, the
stopping gradient, is reached but passes nothing on. The hazard area of cell $j$ is

$$
A_{h,j} = \sum_i \frac{a_i}{l_j}\, P_{m,i}\, f_{ij}
$$

with $a_i$ the cell area, $l_j$ the cell width, $P_{m,i}$ equal to 1 on sources,
and $f_{ij}$ the fraction of cell $i$'s material that arrives at $j$ along paths
that never pass through a stopping cell. One sweep in topological order computes
it for every cell.

**Parameters.** The defaults, 0.839 and 0.176, are the gradients of the paper's
six-earthquake average (40° and 10°) on 30 m DEMs. The best values differed between
earthquakes, and skill depends on the DEM and its resolution. They are effective
values and do not measure a material property.

**Departures from the paper's text.**

- Eq. 7 writes the connection as binary, but the paper reports hazard areas near
  0.1 m²/m on a 30 m grid. A binary count cannot be below one cell width, so lsdtt3
  routes flow fractions. The two readings give the same cells with $A_h = 0$, but
  rank the other cells differently.
- Fig. 2e shows zero at two cells that its own slopes and eqs. 6–9 give as 30 m²/m.
  lsdtt3 follows the equations.
- The paper filled pits. The CLI carves and then fills by default. Set
  `breach_before_fill: false` to use filling alone.

**Limits.** $A_h$ is a relative index. It ignores shaking, rock strength, wetness,
landslide size and run-out across reverse slopes. Where debris from a DEM-edge cell
arrives, sources may lie outside the DEM and $A_h$ may be too low; the
`_hazard_area_edge_influenced.tif` raster marks these cells. Published skill scores
depend on the DEM, its resolution and the study area, so they do not transfer
directly to another site.

## Scar volumes from a draped cloth

Landslide volumes are best measured by differencing DEMs from before and after
an event, but usually only a post-event DEM and a mapped inventory exist.
`run_scar_drape_volumes` reconstructs the pre-event surface over each scar with
a cloth and integrates the gap.

**The problem solved.** For one scar, let $\Omega$ be the polygon's cells grown
by `scar_drape_rim_buffer_pixels` (8-connected), and $\partial\Omega$ the cells
4-adjacent to $\Omega$ outside it. The cloth $z$ satisfies

$$
\nabla^2 z = s \ \text{where } z > z_{\text{DEM}}, \qquad z \ge z_{\text{DEM}} \ \text{in } \Omega, \qquad z = z_{\text{DEM}} \ \text{on } \partial\Omega,
$$

an obstacle problem: the shape of a membrane under tension, held on the rim and
resting on any ground that rises into it. $s \ge 0$ is the sag curvature
(`scar_drape_sag_curvature_per_metre`, 1/m). $s = 0$ is a taut cloth, harmonic
where it is clear of the ground. It is discretised with the 5-point Laplacian
and solved by projected successive over-relaxation, with
$\omega = 2/(1 + \sin(\pi/N))$ for a window $N$ cells across. The problem is
convex, so the solution is unique. Cells off the DEM or with nodata carry no
flux and are counted. Scars are independent and run in parallel. The volume is
$\sum (z - z_{\text{DEM}})\,\Delta x^2$ over the mapped cells, and also over
$\Omega$.

**Relation to the cloth in `lsdtt-surface-metrics`.** That cloth (after Zhang et
al. 2016) falls under constant gravity and rests where its spring force
balances gravity. Away from the ground each cell then sits gravity/rigidity
metres below its neighbours' mean, so the sag depends on the grid. Here the
resting shape is solved directly, and the sag is a curvature in physical units.
With the cloth pinned outside the scar and gravity set to zero, the two give the
same surface.

**What it can and cannot recover.** Linear interpolation across the scar is
exact on a planar slope. Across a convex nose it lies too low (volume too
small), and across a hollow, or a scar cut into an older gully, too high (volume
too large). Deposits left in the scar, or a polygon drawn inside the crown, put
the ring or the floor too high and the volume too low.

**Fitted surfaces** (`scar_drape_surface_type`). Two alternatives to the
membrane, both kept no lower than the ground ($z = \max(z_{\text{DEM}}, f)$):
- `rim_plane`: $f$ is the least-squares plane through the band of cells up to
  `scar_drape_rim_fit_width_pixels` (8-connected rings) outside $\Omega$. Exact
  on a planar slope, like the membrane.
- `rim_quadric`: $f = a + bx + cy + dx^2 + exy + fy^2$ by least squares through
  the same band. Exact wherever the pre-event surface had uniform curvature: a
  bowl, a hollow or a nose, where the membrane is too deep or too shallow. The
  band must be at least two cells wide. On a closed ring the six coefficients
  are barely determined: on a circle $x^2 + y^2$ is constant, so its
  coefficient (the mean curvature, which sets the depth at the centre) is not
  fixed by the ring at all. A band of width $w$ fixes it through the radial
  change across the band.

Coordinates are centred and scaled before the normal equations are solved
(Gaussian elimination with partial pivoting); a band that cannot determine the
surface (a pivot below $10^{-10}$ of the largest diagonal) gives status
`rim_fit_failed`. `surface_fit_n_cells` and `surface_fit_residual_m` (RMS misfit
in the band) are reported; a large residual means the surface is not a good
description of the ground around the scar.

**Uncertainty.** Two kinds are reported separately and not added:
- *DEM error* (`scar_drape_uncertainty_samples`): Monte Carlo. Each sample adds
  Gaussian error to the post-event DEM and re-drapes every scar.
  `volume_sigma_m3` is the sample standard deviation of each scar's mapped
  volume, and `volume_noise_bias_m3` the mean sample volume minus the volume.
  Because the cloth is re-solved, the error of the rim is carried through the
  surface, and error correlated over the scar and its rim cancels (a uniform
  offset moves cloth and ground together). For a linear surface (`rim_plane`,
  no contact) with white error of standard deviation $\sigma$,
  $\mathrm{var}\,V = a^2\sigma^2\,[\,n + s^\top (X^\top X)^{-1} s\,]$, with $n$
  mapped cells of area $a$, $X$ the band cells' $(1, x, y)$ and $s$ the sum of
  $(1, x, y)$ over the mapped cells; the Monte Carlo reproduces it. The noise
  bias is positive for the membrane: a cloth over a noisy surface rests on the
  highs, so it measures a noise floor. The error is either
  `scar_drape_dem_sigma_metres` with `scar_drape_error_nugget_fraction` and
  `scar_drape_error_correlation_range_metres`, or, when it is 0, the error
  model of the reference DEM of difference: its slope-binned NMAD divided by
  √2 (equal, independent errors in both DEMs), nugget and range, estimated off
  the scars. The noise generator is the one behind the DoD false-detection
  simulations; each sample's seed is a splitmix64 hash of the seed and sample
  number, and results do not depend on the thread count.
- *Rim placement* (`scar_drape_buffer_spread_max_pixels`): the mapped volume
  with rim buffers 0 to that number; `volume_buffer_spread_m3` is max − min.

On East San Lorenzo (81 scars, reference = 2015 DTM), the drape / reference
total is 1.71 (membrane), 1.88 (`rim_plane`) and 1.40 (`rim_quadric`), and the
control 1.01, 0.95 and 0.52 of the reference. With the DoD error model (NMAD
0.32 m, nugget 0.35, range 15 m) and 50 samples, the total volume sigma is
48 m³ and the noise bias +262 m³, while buffers 0–2 span 1,980–4,485 m³: where
the rim is placed matters far more than DEM error.

**Validation outputs.** With `scar_drape_reference_dem_fname` the tool computes
the reference volume of each scar (net, loss and gain), the control drape on the
reference DEM, and a summary (`_scar_drape_validation.csv`). Volume - area fits
use the shared scaling fit (`src/utils/lsdtt_scaling_fit`): OLS and reduced
major axis, with optional bootstrap intervals. With error in area, OLS γ is
attenuated by λ = var(log A) / (var(log A) + var(error)) and RMA by √λ, so both
are reported.

**Test against repeat lidar.** In Gleann Easan Biorach (Arran) the 2012
Phase 2 and 2025 National LiDAR Programme DTMs (1 m) give 115 scars of at least
10 m² with lowering of at least 0.3 m. Draping the 2025 DTM over their outlines,
with the outlines taken from the difference itself:

| Rim buffer (cells) | Total vs. differencing | Median per-scar ratio (IQR) | r, log volumes |
|---|---|---|---|
| 0 | 0.76 | 0.38 (0.16–0.77) | 0.88 |
| 1 | 1.27 | 0.91 (0.40–1.43) | 0.86 |
| 2 | 1.68 | 1.18 (0.65–1.80) | 0.85 |

Draping the 2012 DTM over the same outlines, where the answer should be zero,
gives 21 %, 50 % and 84 % of the differenced volume. These scars enlarge
existing gullies, so the cloth also spans the older gully. With no buffer the
ring lies on ground that was itself lowered by less than the 0.3 m threshold.
The drape ranks scars by volume well, but the total depends on where the rim is
placed to within a factor of about two. This site is a hard case. Discrete
shallow slides on planar slopes, with polygons drawn to the crown, should do
better, and should be checked against repeat surveys where they exist.

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
- **Pseudo-static shaking** (`unit_slope_stability_rupture`, `smoke_slope_stability_rupture_pga`): $F_k$ against the
  stresses, the closed-form optimum against a dense search over 600 random cases, and
  $k = 0$ bit-identical to the static model.
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
- **Coseismic hazard area** (`unit_slope_stability_hazard_area`): the worked example
  of Milledge et al. (2019) Fig. 2, with its two inconsistent cells corrected; an exact
  profile; ordinary MFD accumulation when both thresholds are off; conservation of
  source area; hazard never rising with either threshold; the $A_h = 0$ cells against
  a brute-force search from every source.
- **Scar drape** (`unit_slope_stability_scar_drape`): a pit dug into a planar
  slope restored exactly (volume to 2 × 10⁻⁶) for rim buffers of 0, 1 and 3; the
  Poisson square, centre at −0.0736713 s L², within 0.03 % on a 2 m grid and
  0.008 % on a 1 m grid; the cloth never below the ground and resting on a
  pillar; DEM edges, nodata rings, overlapping scars; the volume - area fit;
  polygons read from a GeoPackage, with ids, the cell-centre rule, small
  polygons, overlap and a missing field.
- **The CLI**: smoke tests of every family and option, and refusals of invalid
  parameter files.

Results on real landscapes depend on the DEM and on parameters that are rarely
known well. Treat them as relative measures of where failure is more likely,
unless the parameters have been calibrated locally.

## References

Hynes-Griffin, M. E. and Franklin, A. G. (1984). Rationalizing the seismic coefficient
method. Miscellaneous Paper GL-84-13, US Army Engineer Waterways Experiment Station.

Hugonnet, R., et al. (2022). Uncertainty analysis of digital elevation models by
spatial inference from stable terrain. *IEEE JSTARS* 15, 6456-6472.

Milledge, D. G., Densmore, A. L., Bellugi, D., Rosser, N. J., Watt, J., Li, G. and
Oven, K. J. (2019). Simple rules to minimise exposure to coseismic landslide hazard.
*Natural Hazards and Earth System Sciences* 19, 837-856.

Quinn, P., Beven, K., Chevallier, P. and Planchon, O. (1991). The prediction of
hillslope flow paths for distributed hydrological modelling using digital terrain
models. *Hydrological Processes* 5, 59-79.

Rolstad, C., Haug, T. and Denby, B. (2009). Spatially integrated geodetic glacier
mass balance and its uncertainty based on geostatistical analysis. *Journal of
Glaciology* 55(192), 666-680.

Zhang, W., Qi, J., Wan, P., Wang, H., Xie, D., Wang, X. and Yan, G. (2016). An
easy-to-use airborne LiDAR data filtering method based on cloth simulation.
*Remote Sensing* 8(6), 501. https://doi.org/10.3390/rs8060501

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
