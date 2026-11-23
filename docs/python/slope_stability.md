# Slope stability

These functions and classes run the routines of
[`lsdtt-slope-stability`](../cli_reference/cli_reference_lsdtt_slope_stability.md)
on rasters in memory. What each model assumes, how to choose between them and
what the results mean are in
[Slope stability](../design_and_algorithms/slope_stability.md).

The options have the CLI's names without the `sinmap_`, `rupture_`,
`slip_surface_` and `groundwater_` prefixes, and the tables have the CLI's
columns. Python and the CLI call the same C++ routines. The CLI holds its
parameters in single precision: its default friction coefficient is
tan 35° rounded to a float, 0.7002075, while Python's is the double. Pass the
CLI's values, rounded the same way, to reproduce a CLI run bit for bit; the
default differs from it by about 2 parts in 10⁷.

!!! note "Experimental"
    These bindings are new. `helper_for_llm` marks them `[experimental]` until
    they are promoted.

| Python | Returns | CLI equivalent |
|---|---|---|
| `lsdtt3.sinmap(slope, drainage_area, ...)` | dict of Rasters, counts and the class-area table | `run_sinmap` |
| `lsdtt3.sinmap_slope_area_curves(slope_max, area_min, area_max, ...)` | DataFrame | `_sinmap_slope_area_curves.csv` |
| `FlowGraph.dinf_contour_width()`, `FlowGraph.dinf_edge_contaminated_mask()` | Raster | `contour_width_method`, `exclude_edge_contaminated_cells` |
| `FlowInfo.rupture_planes(geometry, ...)` | `RupturePlanes` | `run_rupture_planes` |
| `RupturePlanes.failure_probability(...)` | dict | `parameter_uncertainty: equal_density` |
| `RupturePlanes.landslide_catalogues(...)` | dict of DataFrames and label Rasters | `write_landslide_catalogue`, `write_landslide_label_rasters` |
| `RupturePlanes.slip_surfaces(family=..., ...)` | `SlipSurfaces` | `slip_surface_family` |
| `SlipSurfaces.failure_probability()`, `.spencer_check()`, `.body_catalogues()` | dicts | the slip-surface probability, Spencer check and body catalogues |
| `lsdtt3.landslide_frequency_area(areas_m2, cell_area_m2, ...)` | (DataFrame, dict) | `_frequency_area.csv`, `_frequency_area_fit.csv` |
| `lsdtt3.GroundwaterDupuit(surface, fixed_head, ...)` | the aquifer | `water_table_method` |
| `lsdtt3.dem_change_landslides(difference, dem, ...)` | dict of tables and Rasters | `run_dem_change_landslides` |

## Routing surface, slope and area

The CLI carves, then fills, then routes. SINMAP uses the D-infinity slope and
area of that surface; the rupture planes follow its D8 paths but measure
heights on the input DEM.

```python
import lsdtt3

dem = lsdtt3.Raster.read("my_site_DEM.tif")
routing = lsdtt3.fill_pits(lsdtt3.breach_least_cost(dem, fill_remainder=False))

graph = lsdtt3.FlowGraph(routing, method="dinf")
flow = lsdtt3.FlowInfo(routing)
```

## SINMAP

```python
s = lsdtt3.sinmap(graph.slope(), graph.drainage_area(),
                  contour_width=graph.dinf_contour_width(),   # the CLI's default
                  dimensionless_cohesion_max=0.25,
                  friction_coefficient_min=0.6, friction_coefficient_max=1.0)
s["stability_class"].write("my_site_sinmap_stability_class.tif")
s["class_areas"]          # region_id, stability_class, n_cells, area_m2, fraction_of_region
```

The result also holds `stability_index`, `probability_stable`, `fs_min`,
`fs_max`, `saturation` and `specific_catchment_area`. Pass
`exclude=graph.dinf_edge_contaminated_mask()` to leave out cells whose
contributing area may be truncated by the raster edge.

## Rupture planes and potential landslides

```python
rp = flow.rupture_planes(dem, threshold_contributing_pixels=1000,
                         cohesion_pa=1000.0, friction_coefficient=0.7)
rp.safety_factor, rp.depth, rp.unstable      # Rasters
rp.n_unstable

cat = rp.landslide_catalogues(min_daylight_neighbours=2)
cat["amalgamated"]                           # the CLI's _landslide_catalogue_amalgamated.csv columns

density, fit = lsdtt3.landslide_frequency_area(cat["amalgamated"]["area_m2"],
                                               dem.metadata.cell_size ** 2)
fit["eta"], fit["eta_standard_error"], fit["area_min_m2"]
```

`channel=` takes a channel mask in place of the contributing-pixel threshold.

## Probability of failure

```python
pf = rp.failure_probability(cohesion_pa_min=500.0, cohesion_pa_max=1500.0,
                            friction_coefficient_min=0.6, friction_coefficient_max=0.8,
                            uncertainty_samples_per_parameter=5)
pf["failure_probability"]                     # P(F <= 1)
```

`cohesion_distribution="normal"` (and `friction_coefficient_distribution`)
uses a normal distribution about the central value, with `cohesion_pa_sd`,
truncated at ±3 SD and at 0.

## Slip surfaces and ellipsoids

```python
ss = rp.slip_surfaces(family="circular")      # method="auto" -> Bishop
ss.safety_factor, ss.body_depth, ss.evaluated
check = ss.spencer_check()                     # Spencer on each critical surface
bodies = ss.body_catalogues()                  # by_toe, amalgamated, label_amalgamated
pf_slip = ss.failure_probability()             # a lower bound (fixed critical surface)

ell = rp.slip_surfaces(family="ellipsoid", column_method="xie_revised")
ell.cover_safety_factor
```

The other families are `plane_integrated`, `bilinear` and `structural`; give
`structural_dip_gradient` and `structural_dip_direction_degrees` as numbers or
Rasters. Crowns whose rupture F is above `screen_safety_factor` (1.5) are not
evaluated: `ss.evaluated` is 0 there.

## Water table

```python
aquifer = lsdtt3.GroundwaterDupuit(dem, rp.channel, regolith_thickness_metres=2.0,
                                   hydraulic_conductivity_m_per_s=1e-5,
                                   drainable_porosity=0.2)
aquifer.run_to_steady(1e-8)                    # implicit by default
aquifer.steady_rates()                         # m3/s: recharge, to channels, seepage

# a storm, then the planes with that water table
aquifer.run(12 * 3600.0, 5e-7)                  # explicit by default
wet = flow.rupture_planes(dem, water_table_depth=aquifer.water_table_depth())
aquifer.budget()                               # cumulative m3
```

`run` advances one constant-recharge segment; call it once per segment of a
recharge series, as the CLI does. `transient_method="implicit"` uses backward
Euler with steps of at most `implicit_max_dt_s` (3600 s).

## Potential landslides from a DEM of difference

`dem_change_landslides` is `run_dem_change_landslides`. The keywords are the CLI's
without the `dem_change_` and `dem_difference_` prefixes. The results are candidate
areas of change, not mapped landslides.

```python
out = lsdtt3.dem_change_landslides(
    difference, newer_dem, stable_mask=stable, change_sign="both", min_slope=0.1
)
print(out["uncertainty"]["nmad_metres"], out["false_detection"]["mean_false_area_m2"])
patches = out["patches"]
loss = patches[patches.change == "loss"]
density, fit = lsdtt3.landslide_frequency_area(loss.area_m2, newer_dem.metadata.cell_size ** 2)
```

`coregister_dem(reference, moving, stable_mask=stable)` fits the shift on stable
terrain only.

## Regions

A region raster and a dict of parameters per region replace the CLI's region
raster and table. The keys are the table's column names, so a CLI region table
can be passed straight in:

```python
import pandas as pd

table = pd.read_csv("my_site_regions.csv").set_index("region_id")
params = table.to_dict("index")              # {region_id: {"rupture_cohesion_pa": ..., ...}}
region = lsdtt3.Raster.read("my_site_regions.tif")

rp = flow.rupture_planes(dem, region=region, region_parameters=params)
s = lsdtt3.sinmap(graph.slope(), graph.drainage_area(), region=region,
                  region_parameters=params)
```

A name missing from a region takes the keyword value, `sinmap_*` names are
ignored by the rupture planes (and `rupture_*` by SINMAP), and an unknown name
is an error. The per-region bounds (`rupture_cohesion_pa_min`, ...) are used by
`failure_probability`.

Figures of the CLI's outputs (the slope - area plot, catalogues,
frequency - area with its fit, and maps) are drawn by
`lsdviz3-slope-stability` in lsdviztools3.
