# Hydrological Modelling in lsdtt3

This page explains the three hydrology tools available in lsdtt3: what each computes, when you should use each one, where the theory comes from, and how they differ fundamentally. It is not a CLI reference — for parameter names and default values see the CLI reference pages for [lsdtt-drainage-analysis](../cli_reference/cli_reference_lsdtt_drainage_analysis.md), [lsdtt-distributed-hydrology](../cli_reference/cli_reference_lsdtt_distributed_hydrology.md), and [lsdtt-graphflood](../cli_reference/cli_reference_lsdtt_graphflood.md). Instead, this page is for understanding which tool answers your research question.

If you are new to lsdtt3 hydrology, start with [The three tools](#the-three-tools) to understand the distinction, then move to [When to use each tool](#when-to-use-each-tool) to match your question to a method. If you already know which tool you want, turn to the corresponding CLI reference.

---

## Overview

lsdtt3 offers three hydrological tools that approach the problem from different angles:

1. **lsdtt-drainage-analysis**: Computes **static topographic indices** (drainage area, Topographic Wetness Index, Soil Topographic Index, connectivity). No time, no rainfall, no water balance. This is what the catchment *looks like*.

2. **lsdtt-distributed-hydrology**: A **time-evolving rainfall-driven hydrological simulation** with a full water balance, built on TOPMODEL. Simulates how water moves vertically (interception → infiltration → root zone → unsaturated zone → saturated zone) and laterally (surface and subsurface routing). This is what happens *during a storm*.

3. **lsdtt-graphflood**: A **stationary graph-based hydraulic routing** solver with Manning friction. Routes a supplied discharge over an evolving water surface (z + h, where z is bed elevation and h is water depth). This is where the water *goes and how deep* it is.

**Critical distinction: hydrological versus hydraulic routing.** The second tool *generates* discharge from rainfall and routes it; the third tool *takes* a discharge as input and routes it while computing depths. They are complementary, not alternatives.

| Question | `lsdtt-drainage-analysis` | `lsdtt-distributed-hydrology` | `lsdtt-graphflood` |
|---|---|---|---|
| What is supplied? | DEM and optional soil data | DEM, rainfall, PET, and optional soil/land-cover data | DEM, roughness, and prescribed discharge or inflow |
| What is computed? | Static drainage and wetness metrics | Rainfall-runoff response and outlet hydrograph | Flow depth and inundation for the supplied water |
| Does it have time? | No | Yes | Yes |
| Does it generate discharge? | No | Yes | No |
| Best first question | Where is the catchment prone to wetness? | How does this storm produce runoff? | Where does this discharge inundate? |

---

## The three tools

### lsdtt-drainage-analysis: static topographic indices

This tool computes **metrics of the catchment's structure**, using only the DEM and optionally a soil data table. It produces no time series and has no weather. The metrics are:

- **Drainage area** (`A`): the accumulated upslope area contributing to each cell. Computed using D-infinity or MFD routing (D8 is not offered because it produces stripe artifacts in the Wetness Index).
- **Topographic Wetness Index** (TWI, also called the wetness index): `ln(a / tan β)`, where `a` is the specific contributing area (area per unit contour width) and `β` is the local slope. High TWI means a place where water tends to gather (convergent topography + large upslope area); low TWI means a ridge or hillslope. TWI is the foundation of TOPMODEL and is also used in lsdtt-distributed-hydrology to diagnose saturation.
- **Soil Topographic Index** (STI): `ln(a / tan β)` weighted by soil transmissivity. Beven & Kirkby (1979) showed that saturated areas group along contours of constant ln(a/tanβ), so the Index is a proxy for the areal pattern of saturation. This tool computes STI from measured soil data (soil series and optional land-use table) rather than a generic proxy.
- **Connectivity**: how easily water moves downslope. This is a placeholder for future work.

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

**Output**: outlet streamflow time series plus rasters of saturation deficit, saturation, velocity, travel time, and TWI. Uniform-transmissivity runs require calibration of `ln_t0` and the delay parameter `td`; STI-table runs instead derive spatial transmissivity from the supplied measured profiles.

### lsdtt-graphflood: graph-based hydraulic routing

This tool routes a **supplied discharge** over the DEM using stationary graph-based hydraulic routing (Manning's equation). It is not a rainfall-runoff simulator — you must supply the water budget (discharge time series, or a constant source). The model computes:

- **Water depth** at each cell as flow accumulates downslope, using Manning friction.
- **Inundation extent and depth** across the landscape, including floodplains and out-of-bank regions.
- **Flow direction** at each cell based on the water surface (z + h), which can differ from the bed elevation surface.

**Input**: DEM, discharge time series (or point inflows), Manning roughness field, channel extent (optional).

**Output**: rasters of water depth, water surface elevation, specific discharge, bed and water-surface slopes.

The implementation is conceptually inspired by the graph-based hydraulic-routing approach described by Gailleton et al. (2024), using a fixed flow graph (no dynamic rerouting) and Manning friction with settable depth and slope exponents to accommodate different friction laws (Manning, Chézy, etc.).

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

- You have a **known discharge time series** (from a gauge, or from lsdtt-distributed-hydrology, or from a regional model) and want to **map inundation depth and extent**.
- You want to route a **design flood** or test a specific discharge scenario on your landscape.
- You need to include **out-of-bank flow** and **floodplain routing** with physically justified friction.
- You are studying **flood propagation**, wave celerity, or the effect of channel geometry on timing.

**Limitations**: You must supply the discharge — it does not compute it from rainfall. It is a **hydraulic** model (routes water given a discharge), not a **hydrological** model (generates discharge from rainfall). Manning's equation is only appropriate for subcritical flow with low Froude number; it degrades in steep cascades or hydraulic jumps. The flow graph is static (does not recompute routing as the water surface changes), which is a reasonable approximation for transient flows but may not capture backwater in very flat reaches.

### Neither model alone captures modern flood forecasting

Real flood forecasting often needs both: run lsdtt-distributed-hydrology to get the hydrograph, then feed that into lsdtt-graphflood to map depths and inundation. The two are complementary.

---

## Provenance and validation

Understanding where each model comes from matters for knowing how far to trust it.

### lsdtt-distributed-hydrology: TOPMODEL lineage

The core is **TOPMODEL**, introduced by Beven & Kirkby (1979) as a semi-distributed rainfall-runoff model built around the Topographic Wetness Index. The index itself — `ln(a/tanβ)` — predicts which parts of a catchment will saturate during a storm, and Beven & Kirkby showed that the model could reproduce observed hydrographs using this structure plus a simple subsurface-routing assumption (water table parallel to a recharge-steady solution).

lsdtt-distributed-hydrology is validated against **TOPMODEL 95.02**, the published Fortran implementation (Beven et al. 1995). Its column physics is tested against the reference output, and its infiltration-excess routine is tested against a frozen trace of the original `EXPINF` subroutine. These implementation checks do not replace calibration and evaluation against observed catchment data.

The distributed form is informed by:

- **Dynamic TOPMODEL** (Beven & Freer 2001): a fully prognostic saturated zone where `D_i` (saturation deficit) is a state variable at every cell, not derived from a catchment-mean.
- **SD-TOPMODEL** (Gao, Holden & Kirkby 2015): spatially distributed roughness and velocity fields built on Darcy-Weisbach resistance, where velocity is dynamic (depth-dependent). Zhu et al. (2025) report calibrations of this approach on UK sites.
- **Routing numerics** from Li, Simons & Stevens (1975) (implicit kinematic wave), Cunge (1969) and Ponce & Yevjevich (1978) (Muskingum-Cunge matching diffusion), and Ponce, Li & Simons (1978) on when kinematic routing is valid.

**Why this matters for your model**: a TOPMODEL parameter set published for a different river basin is not a substitute for calibration. Parameters such as uniform `ln_t0` and `m` are catchment-specific. The index mode is a **validation mode** that reproduces the TOPMODEL reference case as a tripwire; it does not establish predictive skill for a new site.

### lsdtt-graphflood: Gailleton et al. (2024)

The implementation is conceptually inspired by graph-based hydraulic-routing ideas described by **Gailleton et al. (2024)**. The solver uses a **fixed routing graph** (copied from the flow routing at the start and not updated), with **Manning friction** and **implicit numerical schemes** that are unconditionally stable, so the timestep is not restricted by a Courant criterion.

The merit of the Gailleton approach is that it treats routing as a network problem (directed acyclic graph) rather than a grid stencil problem, making it naturally parallel (within each topological level) and avoiding the checkerboard artifacts of some stencil schemes. The fixed-graph assumption is reasonable for transient floods (minutes to hours) but may miss backwater effects in very flat reaches that take longer to equilibrate.

**Why this matters**: graphflood produces physically grounded depths using industry-standard Manning friction, so depths are directly comparable to observations and to regulatory design-flood models. However, it does not include dynamic rerouting (the flow path does not change as the water surface rises), which is a simplification but a fast one.

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
- Beven, K.J., Buytaert, W., & Pfister, L. (2021). *Fifty years of TOPMODEL: a review.* Hydrology and Earth System Sciences 25:527–549.

**Distributed TOPMODEL**

- Gao, H., Holden, J., & Kirkby, M. (2015). *A distributed model for simulating water table rise and surface saturation area with implications for wetland flood management.* Hydrological Processes 29:2867–2879.
- Beven, K.J., & Freer, J. (2001). *A dynamic TOPMODEL.* Hydrological Processes 15:1993–2011.
- Zhu, S., Peng, J., Gao, X., Liao, J., & Lü, N. (2025). *SD-TOPMODEL coupling with a distributed flood routing module: application to flood management in UK catchments.* Hydrological Processes 39(1):e70122.

**Routing and hydrology numerics**

- Li, R.M., Simons, D.B., & Stevens, M.A. (1975). *Nonlinear kinematic wave approximation for water routing.* Water Resources Research 11(2):245–252.
- Cunge, J.A. (1969). *On the subject of a flood propagation computation method (Muskingum method).* Journal of Hydraulic Research 7(2):205–230.
- Ponce, V.M., & Yevjevich, V. (1978). *Muskingum-Cunge method with variable parameters.* Journal of the Hydraulics Division, ASCE 104(12):1663–1667.
- Ponce, V.M., Li, R.M., & Simons, D.B. (1978). *Applicability of kinematic and diffusion models.* Journal of the Hydraulics Division, ASCE 104(3):353–360.

**Hydraulic routing**

- Gailleton, B., Carretier, S., Durand, S., Mudd, S.M., & Jiao, R. (2024). *Tracing water and sediment routes: a graph-based approach for routing on rasters.* Environmental Modelling & Software 182:105927.

**Soil transmissivity and green-ampt infiltration**

- Beven, K.J. (1984). *Infiltration into a class of vertically non-uniform soils.* Hydrological Sciences Journal 29(4):425–434.

**Interception and vegetation**

- Dickinson, R.E. (1984). *Modeling evapotranspiration for three-dimensional global climate models.* Climate Processes and Climate Sensitivity, Geophysical Monograph 29:58–72.
- van Dijk, A.I.J.M., & Bruijnzeel, L.A. (2001). *Forest patch size and LAI affect throughfall and stemflow rates in the lutz spruce forest in Germany.* Journal of Hydrology 247:230–238.
- Aston, A.R. (1979). *Rainfall interception by eight small trees.* Journal of Hydrology 42:383–396.

**Indices and catchment structure**

- Leopold, L.B., & Maddock, T. (1953). *The hydraulic geometry of stream channels and some physiographic implications.* U.S. Geological Survey Professional Paper 252.
- Rickenmann, D., & Recking, A. (2011). *Evaluation of flow resistance equations using a large field database.* Water Resources Research 47:W07538.
- Ferguson, R. (2007). *Flow resistance in gravel-bed rivers at bankfull stage.* Journal of Hydraulic Engineering 133(6):671–678.
