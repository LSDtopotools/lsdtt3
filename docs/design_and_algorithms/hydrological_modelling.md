# Hydrological Modelling in lsdtt3

This page explains the three hydrology tools available in lsdtt3: what each computes, when you should use each one, where the theory comes from, and how they differ fundamentally. It is not a CLI reference — for parameter names and default values see the CLI reference pages for [lsdtt-drainage-analysis](../cli_reference/cli_reference_lsdtt_drainage_analysis.md), [lsdtt-distributed-hydrology](../cli_reference/cli_reference_lsdtt_distributed_hydrology.md), and [lsdtt-graphflood](../cli_reference/cli_reference_lsdtt_graphflood.md). Instead, this page is for understanding which tool answers your research question.

If you are new to lsdtt3 hydrology, start with [The three tools](#the-three-tools) to understand the distinction, then move to [When to use each tool](#when-to-use-each-tool) to match your question to a method. If you already know which tool you want, turn to the corresponding CLI reference.

---

## Overview

lsdtt3 offers three hydrological tools that approach the problem from different angles:

1. **lsdtt-drainage-analysis**: Computes **static topographic indices** (drainage area, Topographic Wetness Index, Soil Topographic Index, connectivity). No time, no rainfall, no water balance. This is what the catchment *looks like*.

2. **lsdtt-distributed-hydrology**: A **time-evolving rainfall-driven hydrological simulation** with a full water balance, built on TOPMODEL. Simulates how water moves vertically (interception → infiltration → root zone → unsaturated zone → saturated zone) and laterally (surface and subsurface routing). This is what happens *during a storm*.

3. **lsdtt-graphflood**: A **stationary graph-based hydraulic routing** solver with Manning friction. Routes a supplied discharge over an evolving water surface (z + h, where z is bed elevation and h is water depth). This is where the water *goes and how deep* it is.

**Critical distinction: hydrological versus hydraulic routing.** The second tool *generates* discharge from rainfall and routes it through time; the third tool takes a *steady* water input and solves for the depths it produces. They are complementary, not alternatives.

| Question | `lsdtt-drainage-analysis` | `lsdtt-distributed-hydrology` | `lsdtt-graphflood` |
|---|---|---|---|
| What is supplied? | DEM and optional soil data | DEM, rainfall, PET, and optional soil/land-cover data | DEM, roughness field, and a rainfall-rate field (optionally a point inflow) |
| What is computed? | Static drainage and wetness metrics | Rainfall-runoff response and outlet hydrograph | Flow depth and inundation for the supplied water |
| Does it have time? | No | Yes | No (stationary solver) |
| Does it generate discharge? | No | Yes | No |
| Best first question | Where is the catchment prone to wetness? | How does this storm produce runoff? | Where does this steady water input inundate? |

---

## The three tools

### lsdtt-drainage-analysis: static topographic indices

This tool computes **metrics of the catchment's structure**, using only the DEM and optionally a soil data table. It produces no time series and has no weather. The metrics are:

- **Drainage area** (`A`): the accumulated upslope area contributing to each cell. D8, D-infinity and MFD are all available (`write_d8_drainage_area_raster` gives the D8 field). The D8 restriction applies only to the indices built on specific area: `twi_method` and `sti_method` do not offer it, because sending all of a cell's area to one neighbour stripes the resulting index.
- **Topographic Wetness Index** (TWI, also called the wetness index): `ln(a / tan β)`, where `a` is the specific contributing area (area per unit contour width) and `β` is the local slope. High TWI means a place where water tends to gather (convergent topography + large upslope area); low TWI means a ridge or hillslope. TWI is the foundation of TOPMODEL and is also used in lsdtt-distributed-hydrology to diagnose saturation.
- **Soil Topographic Index** (STI): `ln(a / tan β)` weighted by soil transmissivity. Beven & Kirkby (1979) showed that saturated areas group along contours of constant ln(a/tanβ), so the Index is a proxy for the areal pattern of saturation. This tool computes STI from measured soil data (soil series and optional land-use table) rather than a generic proxy.
- **Index of Connectivity** (IC): how readily sediment and water move from a cell to a downslope target, after Cavalli et al. (2013). It is `IC = log10(D_up / D_down)`, where the upslope component `D_up = W̄ S̄ √A` combines mean weighting factor, mean slope and contributing area over the upslope region, and the downslope component `D_down = Σ d_i / (W_i S_i)` accumulates path length along the flow route, penalised by local weighting and slope. The weighting factor is a roughness measure (`connectivity_roughness_method`, default `madk2`), so rough surfaces impede connection and smooth ones promote it. Controlled by `write_connectivity_raster`, `connectivity_drainage_area_method`, `connectivity_area_threshold_m2` (default 1e6), `connectivity_window_radius_metres` (default 100), `connectivity_slope_method` (default `polyfit_direct`) and `connectivity_roughness_method`.

**Output**: rasters, no time dimension. No calibration needed; results depend only on the DEM and (optionally) the soil data table.

### lsdtt-distributed-hydrology: TOPMODEL with distributed stores

This tool simulates **what happens to water during a rainfall event**, hour by hour (or at whatever timestep you choose). It is a raster implementation of TOPMODEL (Beven & Kirkby 1979), extended with:

- **Stores per cell**: optional canopy interception, root zone, unsaturated zone, saturated zone, and surface-routing storage.
- **Optional infiltration-excess routing**: where rainfall intensity exceeds the soil's ability to accept water, excess runs overland immediately.
- **Saturation-excess runoff**: when the saturated zone fills (deficit reaches zero), water spills as overland flow.
- **Dynamic velocity field**: flow speed depends on water depth, slope, and surface roughness (Manning's equation). Low-order channels are slower than high-order ones because they are shallower and have coarser bed material.
- **Two saturated-zone modes**: a quasi-steady `index` mode (the TOPMODEL-compatible validation mode) and a fully prognostic `dynamic` mode (one saturated-zone store per cell).
- **Scenario-driven land-cover change**: the tool reads measured soil profiles keyed on `(soil series, land use)`, so converting a pixel from grazing to woodland re-keys the lookup to the measured transmissivity for that soil and cover, rather than applying a universal multiplier.

**Input**: DEM, rainfall time series, PET time series, soil data table, optional LAI (Leaf Area Index) scenarios, optional vector channel network.

**Output**: outlet streamflow time series plus rasters of saturation deficit, saturation, velocity, travel time, and TWI. Uniform-transmissivity runs require calibration of `ln_t0` (default 5.0) and the unsaturated-zone delay `td_hours_per_metre` (default 50.0); STI-table runs instead derive spatial transmissivity from the supplied measured profiles.

### lsdtt-graphflood: graph-based hydraulic routing

This tool solves for the **steady-state flow depths** produced by a supplied water input, using graph-based hydraulic routing with Manning friction. It is not a rainfall-runoff simulator: it does not generate the water, and it does not evolve it through time. The model computes:

- **Water depth** at each cell as flow accumulates downslope, using Manning friction.
- **Inundation extent and depth** across the landscape, including floodplains and out-of-bank regions.
- **Flow direction** at each cell based on the water surface (z + h), which can differ from the bed elevation surface.

**Input**: DEM, a **rainfall-rate raster** (`rainfall_rate_fname`, in m/s) and a **friction-coefficient raster** (`friction_coefficient_fname`) — both required; the tool exits if either is missing. A single point inflow can be added with `use_point_inflow` and `point_inflow_discharge_m3_per_s`. Routing is MFD by default (`flow_routing_method: mfd`, `mfd_exponent: 1.0`) or D-infinity.

**Output**: rasters of bed elevation, water depth, water surface elevation, `qin`, `qout`, `dhw`, hydraulic slope, specific discharge and total discharge. Note that the hydraulic slope is the water-surface slope; there is no separate bed-slope output.

!!! warning "There is no discharge time series"

    A common misreading is that graphflood accepts a hydrograph. It does not:
    the water input is a **rate field**, and the solver is **stationary**.
    `n_steps` (default 100) is the maximum number of stationary hydraulic
    iterations, and `hydraulic_dt` (default 1e-3 s) is a relaxation parameter
    that controls how fast the iteration approaches its fixed point — neither
    is physical time. There is consequently no way to feed an
    `lsdtt-distributed-hydrology` hydrograph into this tool, and questions
    about flood propagation, wave celerity or the timing of a flood peak are
    outside what it can answer.

The implementation is conceptually inspired by the graph-based hydraulic-routing approach described by Gailleton et al. (2024), with Manning friction and settable depth and slope exponents to accommodate different friction laws (Manning, Chézy, etc.).

---

## When to use each tool

### Use lsdtt-drainage-analysis if:

- You want to understand the **structure of a catchment** without any rainfall data or time dimension.
- You need static metrics like TWI or STI for a paper, or as input to other analyses.
- You are exploring a new study site and want to see where water tends to gather.
- You want to **validate** your model by comparing its TWI to observations (e.g., mapped saturated areas in Google Earth imagery or from field surveys).

### Use lsdtt-distributed-hydrology if:

- You have **rainfall and PET time series** and want to simulate the **hydrological response** (streamflow, groundwater levels, saturation patterns) hour by hour.
- You want to test **land-cover change scenarios** — e.g., "what if we planted trees here?" — and see how that changes the hydrograph.
- You want to study **infiltration-excess flooding** (surface water), **saturation-excess flooding** (groundwater-driven), or both.
- You need to model **interception** (canopy storage of rainfall), which matters for NFM (Natural Flood Management) studies with newly planted woodland.
- You want a physically grounded model with **measured soil data** rather than a black-box emulator.

**Limitations**: The tool has not yet been calibrated against observed catchments. It is an event-scale, grid-based hydrological model, not a full 2D groundwater model. In `index` mode, the quasi-steady topographic-index assumption can break down during rapid transients; use `dynamic` mode and timestep-sensitivity checks when that assumption is unsuitable.

### Use lsdtt-graphflood if:

- You have a **steady water input** — a design rainfall rate, or a constant point inflow — and want to **map the inundation depth and extent it produces**.
- You want to compare inundation between scenarios: different roughness fields, a modified DEM, a different input rate.
- You need to include **out-of-bank flow** and **floodplain routing** with physically justified friction.

**Limitations**: You must supply the water input — it does not compute it from rainfall through a soil column. It is a **hydraulic** model (solves depths given a water input), not a **hydrological** model (generates discharge from weather). Manning's equation is only appropriate for subcritical flow with low Froude number; it degrades in steep cascades or hydraulic jumps.

Most importantly, **it is stationary**: it converges to the depth field a sustained input would produce, and has no time axis. Flood propagation, wave celerity and the timing of a peak are not questions it can answer. If you need those, you need an unsteady hydraulic model, which lsdtt3 does not provide.

### Neither model alone captures modern flood forecasting

Real flood forecasting needs both a hydrological model that generates the discharge and a hydraulic model that routes it — and, for timing, an unsteady one.

**lsdtt3 does not currently connect the two.** `lsdtt-graphflood` takes a rainfall-rate field, not a hydrograph, so there is no parameter through which an `lsdtt-distributed-hydrology` outlet time series can be passed to it. The two tools answer adjacent questions on the same DEM; they are not a pipeline. Coupling them would need either an unsteady hydraulic solver or a convention for turning a hydrograph into a sequence of stationary solves, and neither exists yet.

---

## Provenance and validation

Understanding where each model comes from matters for knowing how far to trust it.

### lsdtt-distributed-hydrology: TOPMODEL lineage

The core is **TOPMODEL**, introduced by Beven & Kirkby (1979) as a semi-distributed rainfall-runoff model built around the Topographic Wetness Index. The index itself — `ln(a/tanβ)` — predicts which parts of a catchment will saturate during a storm, and Beven & Kirkby showed that the model could reproduce observed hydrographs using this structure plus a simple subsurface-routing assumption (water table parallel to a recharge-steady solution).

The column equations follow **TOPMODEL 95.02** (Beven et al. 1995). Eight CTest
suites cover the hydrology modules — `unit_topmodel_column`,
`unit_topmodel_dynamic`, `unit_infiltration_excess`,
`unit_canopy_interception`, `unit_hydro_roughness`, `unit_hydro_routing`,
`unit_hydro_velocity_field` and `unit_giuh_discharge` — checking mass
conservation, closed-form solutions (the exponential baseflow recession, the
Green-Ampt infiltration capacity, the index relation as the steady state of
the dynamic saturated zone), limiting behaviours and convergence.

Two caveats. These are implementation checks: they do not replace calibration
and evaluation against observed catchment data, and **the tool has not been
calibrated against any observed catchment**. And unlike the flow-routing and
cosmogenic work, these suites are not yet written up in
`docs/developer_notes/testing_and_benchmarking.md`, so there is no durable
record of what they measured beyond the test sources themselves.

The distributed form is informed by:

- **Dynamic TOPMODEL** (Beven & Freer 2001): a fully prognostic saturated zone where `D_i` (saturation deficit) is a state variable at every cell, not derived from a catchment-mean.
- **SD-TOPMODEL** (Gao, Holden & Kirkby 2015): spatially distributed roughness and velocity fields built on Darcy-Weisbach resistance, where velocity is dynamic (depth-dependent). Zhu et al. (2025) report calibrations of this approach on UK sites.
- **Routing numerics** from Li, Simons & Stevens (1975) (implicit kinematic wave), Cunge (1969) and Ponce & Yevjevich (1978) (Muskingum-Cunge matching diffusion), and Ponce, Li & Simons (1978) on when kinematic routing is valid.

**Why this matters for your model**: a TOPMODEL parameter set published for a different river basin is not a substitute for calibration. Parameters such as uniform `ln_t0` and the recession scale `szm_metres` (default 0.032 m) are catchment-specific. The index mode reproduces the classic TOPMODEL formulation; it does not establish predictive skill for a new site.

### lsdtt-graphflood: Gailleton et al. (2024)

The implementation is conceptually inspired by graph-based hydraulic-routing ideas described by **Gailleton et al. (2024)**, with **Manning friction** over a routing graph derived from the water surface.

**The routing graph is rebuilt, not fixed.** Flow directions come from `z + h`, so as depth changes the ordering can change too. The solver tracks this with dirty tiles: `mark_dirty()` flags a tile whose depths moved, `dirty_topology_still_valid()` checks whether the ordering within it survived, and `rebuild_topology()` re-derives it when it did not. `reuse_topology`, `dirty_tile_size` (default 64) and `full_validation_every` (default 128) control how aggressively this is amortised, and the run reports `topology_rebuilds` so you can see how often it happened. This matters for backwater: the flow path *does* respond to the rising water surface.

**The depth update is explicit, not implicit.** Each iteration applies forward Euler, `h ← max(0, h + (Δt/A)(q_in − q_out))`. Stability does not come from implicitness; it comes from a flux limiter that caps outflow at the water physically present, `q_available = q_in + h·A/Δt`, so a cell can never discharge more than it holds. This is why `hydraulic_dt` can be small without the run being unstable, and why it should be read as a relaxation parameter rather than a physical timestep.

The merit of the Gailleton approach is that it treats routing as a network problem (a directed acyclic graph) rather than a grid stencil problem, making it naturally parallel within each topological level and avoiding the checkerboard artifacts of some stencil schemes.

**Why this matters**: graphflood produces physically grounded depths using industry-standard Manning friction, so depths are directly comparable to observations and to regulatory design-flood models. The limitation to keep in view is not the graph — it is that the solution is a steady state, so nothing about the timing of a flood can be read from it.

### lsdtt-drainage-analysis: Beven & Kirkby indices

The TWI and STI are built on foundations laid by **Beven & Kirkby (1979)** and the theory of **variable-source-area hydrology**. The premise is simple: saturated areas are not fixed; they expand and contract through a storm as the water table rises and falls. Where will they be? Where topography is convergent (negative curvature) and upslope area is large, the water table is shallowest — so saturation starts there. The ln(a/tanβ) index captures this pattern, and a histogram of the index over a catchment predicts the shape of the saturation-area curve.

The STI layer adds **measured soil transmissivity** (from the HOST database or your own field surveys), keyed on soil series and land use. This is more faithful to the physics than a single basin-wide `ln_t0` parameter, because soils differ in structure and permeability even within the same catchment.

**Why this matters**: lsdtt-drainage-analysis produces no discharge or floods by itself, so there is no validation against gauge data. But the TWI pattern *is* observable in the field — look for seeps, spring lines, and saturated-area patches during a rainfall event or after snow melt. If your TWI does not match the observed saturation pattern, your DEM may be wrong (e.g., inverted hillshade, undetected levees, or feet-on-metre-grid), or your routing method choice (D8 vs D-infinity) may be inaccurate.

---

## Cross-comparison: the broader landscape

### Hydrological vs hydraulic models

The distinction between the second and third tools is worth belaboring because they sound similar but do opposite jobs:

| Aspect | lsdtt-distributed-hydrology | lsdtt-graphflood |
|---|---|---|
| **Input forcing** | Rainfall + PET (generates discharge internally) | Discharge (must be supplied) |
| **Vertical process** | Full water balance: interception, infiltration, storage, groundwater | None; flow enters at specified points or from upstream |
| **Horizontal process** | Routes the generated runoff over hillslopes and channels | Routes supplied discharge with friction |
| **Time scale** | Event scale (hours to days for a storm) | Any time scale; discharge is prescribed |
| **Calibration** | Requires streamflow or groundwater observations | Requires channel friction calibration |
| **Output** | Hydrograph, saturation pattern, water-table depth | Inundation map, flow depth, water-surface slope |

**Use both for modern flood forecasting**: run distributed-hydrology to compute the hydrograph (rain → discharge), then feed that discharge into graphflood to compute depths and inundation. That is the complete workflow.

### Land-cover change: measured soil data vs multipliers

A philosophical difference worth understanding, because it affects how you parameterize land-cover change.

**SD-TOPMODEL** (and many hydrological models) represents land-cover change as a multiplier: convert grazing to woodland, and multiply the conductivity `Ks` by 1.5 to account for reduced soil compaction. Transmissivity `T0` is then `Ks` × profile depth, and the multiplier propagates.

**lsdtt-distributed-hydrology** does this differently. The transmissivity `ln_t0` comes from a lookup table keyed on `(soil series, land use)`. When you convert a pixel from grazing to woodland, the model does not multiply an existing value; it re-keys the lookup. The measured profile for that soil under grazing is replaced by the measured profile for the same soil under woodland. **The two profiles can differ in transmissivity and in sign** — e.g., a clay soil might have lower permeability under woodland (root reinforcement, less tillage), while a sandy soil gains higher permeability (less compaction). The model captures this variation rather than asserting a universal multiplier.

**Why this matters**: if you are modelling NFM (land-cover change for flood mitigation), using measured soil profiles is more faithful to the physics than a multiplier. But you need measured data: either field observations from your site, or a soil-survey database (like the UK's HOST, or USDA SSURGO). A model run with made-up multipliers is no better than a model run with no physics at all.

---

## The process modules

`lsdtt-distributed-hydrology` is an assembly of library modules rather than a
monolithic solver. Each carries its own physics and its own provenance, and
each is separately testable. This section covers the ones a user has to make
decisions about.

### Roughness: three regimes, and where Manning stops working

Roughness is partitioned into three regimes, because they are set by different
things. Hillslope and floodplain `n` come from **land cover**. Channel `n` does
not — it is a **grain-size property**.

That distinction has a consequence worth stating plainly: a land-cover
scenario, such as planting a field or aerating a soil, must leave channel `n`
untouched. Changing what grows on the bank does not change the bed material,
and `rigorous_roughness_partition` asserts exactly that.

Channel roughness follows Strickler, `n ~ D50^(1/6)`, combined with Sternberg
(1875) downstream fining, `D50 ~ A^-b` with `b ~ 0.3-0.7`. Channel `n`
therefore *falls* downstream as `A^(-0.05..-0.12)`. This reinforces the
"low-order channels are slower" behaviour through a second, independent
mechanism: headwater channels are both shallower **and** rougher.

At low relative submergence Manning is not merely imprecise — it is the wrong
shape. In headwater channels `h/D84` is often 1-3, where resistance rises
steeply and no constant `n` describes it. Two modes address this:

| Mode | Form | Note |
|---|---|---|
| `jarrett` | `n = 0.39 S^0.38 R^(-0.16)` | Fitted to high-gradient streams. Depth **and** slope dependent, so `n` is not a constant field. |
| `ferguson_vpe` | Variable power equation | Abandons Manning's form entirely, then reports the **equivalent** `n` that would give the same velocity. |

`ferguson_vpe` keeps one currency — everything downstream still consumes a
roughness raster — while the physics underneath is not Manning's. Both modes
are depth-dependent, which is the same point the velocity field makes:
roughness alone is not a resistance model.

### Infiltration excess: Green-Ampt with declining conductivity

Infiltration-excess generation is a Green-Ampt wetting front in a soil whose
saturated conductivity falls exponentially with depth, after Beven (1984). It
is a port of TOPMODEL 95.02's `EXPINF`.

It was ported **from the FORTRAN, not from the C**. `tmod9502.c`'s translation
of `EXPINF` is structurally broken: `IROF`, `DF` and `CUMF` are passed by value
where FORTRAN passes them by reference as outputs, so the ponding state and the
infiltration increment never reach the caller; `DTH` is re-declared as an
uninitialised local that shadows the parameter, making `cd = hf*dth` garbage;
`fx` is a typo for `fc` in the series and is never assigned; and `rint`, a
rainfall rate, is declared `int`, truncating any realistic storm to zero. The C
author's own header records that it was only validated with `infex = 0`, and
the shipped `params.dat` sets exactly that — so that path was never exercised.

Two further defects are latent in the FORTRAN itself, and lsdtt3 departs from
it on both:

1. `TP` (time to ponding) and `CONST` (the definite-integral constant) are
   local variables that the next call reads back when ponding persists, so
   correctness depended on the compiler giving locals static storage. Many F77
   compilers did; the standard does not require it. Here they are explicit
   state on `infiltration_state`.
2. The final test is `IF (F .LT. CUMF + RINT)`, comparing a depth against a
   depth plus a **rate**. Dimensionally it wants `CUMF + RINT*DT`. The
   published demo runs at `DT` = 1 hour, where the two coincide numerically,
   which is presumably why it survived. `use_dimensional_ponding_test` selects
   the corrected form; **the default reproduces the original**, so a run that
   wants dimensional consistency has to ask for it.

### Canopy interception

A Rutter-type interception store, with both canopy parameters driven by LAI:

```
C = k_S · LAI          canopy storage capacity, k_S ~ 0.2 mm
p = exp(−κ · LAI)      free throughfall fraction, κ ~ 0.5–0.7
```

`C` follows Dickinson (1984) and the `S ~ S_L·LAI` form of van Dijk &
Bruijnzeel (2001); the Beer-Lambert gap fraction for `p` follows Aston (1979).

Wet-canopy evaporation is `E_canopy = min(S, PET · evap_factor)`, and the
residual demand passed on is `PET_residual = PET − E_canopy/evap_factor`. **The
division by `evap_factor` is the part that is easy to drop.** `evap_factor`
exists because wet-canopy evaporation can exceed grass-reference PET for tall,
aerodynamically rough vegetation — the canopy is better ventilated than the
reference surface. With `evap_factor = 2`, one millimetre of canopy evaporation
consumes only half a millimetre of transpiration demand, so total
evapotranspiration legitimately exceeds PET. Subtracting `E_canopy` directly
would silently cap total ET at PET and understate forest water use.

!!! warning "Seasonality is decisive, not decorative"

    Winter deciduous LAI is ~0.5–1 against ~5 in summer — an order of magnitude
    in interception capacity. UK flood events are winter-dominated, so newly
    planted deciduous woodland intercepts **least** exactly when it matters
    most. A study that runs winter storms with annual-mean or summer LAI will
    overstate the benefit of planting. `lai_seasonal_factor` multiplies LAI at
    construction, and it is a parameter you must think about rather than a
    default that quietly applies summer conditions to a January storm.

### Velocity field and routing

The velocity field converts the roughness field into travel times, and the
routing module moves water along them. Three routing schemes are available —
a linear reservoir, Muskingum-Cunge, and a kinematic wave — with the numerics
following Li, Simons & Stevens (1975), Cunge (1969), Ponce & Yevjevich (1978),
and Ponce, Li & Simons (1978) on when kinematic routing is valid.

The point the velocity field makes, and which the roughness modes above
reinforce, is that **roughness alone is not a resistance model**: velocity
depends on depth, so a static `n` raster does not by itself determine travel
time.

### GIUH

The geomorphological instantaneous unit hydrograph is computed from the width
function — the distribution of flow-path lengths to the outlet — with
velocities scaled by slope, `v_i = v_ref (S_i/S_ref)^β` with `β = 0.5`, giving
travel time `τ = L/v`. It follows Rodriguez-Iturbe & Valdes (1979) and Mesa &
Mifflin (1986).

Note that the GIUH lives in `lsdtt-sediment-routing`, not in any of the three
tools above, so it is a fourth route into catchment hydrology in lsdtt3.

---

## Known limitations

### Index vs dynamic mode (distributed-hydrology)

The **index mode** computes the saturation deficit at each cell from a single catchment-mean value plus the TWI. This is cheap and reproduces TOPMODEL's published parameter sets exactly, but it assumes the water table is everywhere parallel to a steady-recharge solution. During rapid transients (sudden heavy rain or snowmelt), different parts of the catchment dry at different rates, and the assumption breaks down.

The **dynamic mode** makes the saturation deficit prognostic (a state variable) at every cell. It is more expensive and should be checked for timestep sensitivity. Its lateral flux follows the flow graph; it is not a full 2D Dupuit-Boussinesq groundwater solution.

### Manning's assumptions (graphflood)

Manning's equation is a **semi-empirical correlation** tuned for subcritical flow in natural channels and floodplains (Froude number < 1, shallow submergence, moderate slopes). It is reasonable for floods but:

- **Breaks down in steep cascades** (Fr > 1, hydraulic jumps, shooting flow).
- **Degrades at very low relative submergence** (depth comparable to grain size). In headwater channels where depth is millimetres to centimetres, Manning's `n` is not constant; Ferguson's (2007) variable power equation is more accurate but requires `D84` grain size rather than just roughness.
- **Is insensitive to air entrainment and aeration**, which matter in highly turbulent cascades.

For UK upland applications (the design drivers for lsdtt3), Manning is reasonable in most channel and floodplain cells, though headwater channels may need special treatment.

---

## References

**TOPMODEL core**

- Beven, K.J., & Kirkby, M.J. (1979). *A physically based, variable contributing area model of basin hydrology.* Hydrological Sciences Bulletin 24(1):43–69.
- Beven, K.J., Lamb, R., Quinn, P., Romanowicz, R., & Freer, J. (1995). TOPMODEL. Computer models of watershed hydrology. Water Resources Publications, pp. 627–668.
- Beven, K.J., Kirkby, M.J., Freer, J.E., & Lamb, R. (2021). *A history of TOPMODEL.* Hydrology and Earth System Sciences 25(2):527–549. doi:10.5194/hess-25-527-2021

**Distributed TOPMODEL**

- Gao, H., Holden, J., & Kirkby, M. (2015). *A distributed TOPMODEL for modelling impacts of land-cover change on river flow in upland peatland catchments.* Hydrological Processes 29(13):2867–2879. doi:10.1002/hyp.10408
- Beven, K.J., & Freer, J. (2001). *A dynamic TOPMODEL.* Hydrological Processes 15:1993–2011.
- Zhu, S., Klaar, M.J., Willis, T., & Holden, J. (2025). *Use of spatially distributed TOPMODEL to assess the effectiveness of diverse natural flood management techniques in a UK catchment.* Hydrological Processes 39(4):e70122. doi:10.1002/hyp.70122

**Routing and hydrology numerics**

- Li, R.M., Simons, D.B., & Stevens, M.A. (1975). *Nonlinear kinematic wave approximation for water routing.* Water Resources Research 11(2):245–252.
- Cunge, J.A. (1969). *On the subject of a flood propagation computation method (Muskingum method).* Journal of Hydraulic Research 7(2):205–230.
- Ponce, V.M., & Yevjevich, V. (1978). *Muskingum-Cunge method with variable parameters.* Journal of the Hydraulics Division, ASCE 104(12):1663–1667.
- Ponce, V.M., Li, R.M., & Simons, D.B. (1978). *Applicability of kinematic and diffusion models.* Journal of the Hydraulics Division, ASCE 104(3):353–360.

**Hydraulic routing**

- Gailleton, B., Steer, P., Davy, P., Schwanghart, W., & Bernard, T. (2024). *GraphFlood 1.0: an efficient algorithm to approximate 2D hydrodynamics for landscape evolution models.* Earth Surface Dynamics 12(6):1295–1313. doi:10.5194/esurf-12-1295-2024

**Soil transmissivity and green-ampt infiltration**

- Beven, K.J. (1984). *Infiltration into a class of vertically non-uniform soils.* Hydrological Sciences Journal 29(4):425–434.

**Interception and vegetation**

- Dickinson, R.E. (1984). *Modeling evapotranspiration for three-dimensional global climate models.* Climate Processes and Climate Sensitivity, Geophysical Monograph 29:58–72.
- van Dijk, A.I.J.M., & Bruijnzeel, L.A. (2001). *Modelling rainfall interception by vegetation of variable density using an adapted analytical model. Part 1. Model description.* Journal of Hydrology 247:230–238.
- Aston, A.R. (1979). *Rainfall interception by eight small trees.* Journal of Hydrology 42:383–396.

**Roughness and flow resistance**

- Strickler, A. (1923). *Beiträge zur Frage der Geschwindigkeitsformel und der Rauhigkeitszahlen für Ströme, Kanäle und geschlossene Leitungen.* Mitteilungen des Eidgenössischen Amtes für Wasserwirtschaft 16, Bern.
- Sternberg, H. (1875). *Untersuchungen über Längen- und Querprofil geschiebeführender Flüsse.* Zeitschrift für Bauwesen 25:483–506.
- Jarrett, R.D. (1984). *Hydraulics of high-gradient streams.* Journal of Hydraulic Engineering 110(11):1519–1539. doi:10.1061/(ASCE)0733-9429(1984)110:11(1519)

**Unit hydrograph theory**

- Rodríguez-Iturbe, I., & Valdés, J.B. (1979). *The geomorphologic structure of hydrologic response.* Water Resources Research 15(6):1409–1420. doi:10.1029/WR015i006p01409
- Mesa, O.J., & Mifflin, E.R. (1986). *On the relative role of hillslope and network geometry in hydrologic response.* In: Gupta, V.K., Rodríguez-Iturbe, I. & Wood, E.F. (eds), Scale Problems in Hydrology. Reidel, Dordrecht. <!-- page range not independently verified -->

**Sediment and hydrological connectivity**

- Cavalli, M., Trevisani, S., Comiti, F., & Marchi, L. (2013). *Geomorphometric assessment of spatial sediment connectivity in small Alpine catchments.* Geomorphology 188:31–41. doi:10.1016/j.geomorph.2012.05.007

**Indices and catchment structure**

- Leopold, L.B., & Maddock, T. (1953). *The hydraulic geometry of stream channels and some physiographic implications.* U.S. Geological Survey Professional Paper 252. doi:10.3133/pp252
- Rickenmann, D., & Recking, A. (2011). *Evaluation of flow resistance in gravel-bed rivers through a large field data set.* Water Resources Research 47:W07538. doi:10.1029/2010WR009793
- Ferguson, R. (2007). *Flow resistance equations for gravel- and boulder-bed streams.* Water Resources Research 43(5):W05427. doi:10.1029/2006WR005422
