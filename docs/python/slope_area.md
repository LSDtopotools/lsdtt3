# Slope-area and low-concavity reaches

These functions run the routines of
[`lsdtt-slope-area`](../cli_reference/cli_reference_lsdtt_slope_area.md) on
rasters in memory. The method is explained in
[Slope-area analysis and low-concavity reaches](../design_and_algorithms/slope_area_and_low_concavity.md).
The options have the CLI's names without the `low_concavity_` prefix, and the
tables have the CLI's columns.

| Function | Returns | CLI output |
|---|---|---|
| `lsdtt3.conditioned_pixel_mask(input_dem, routing_dem)` | Raster: 0 untouched, 1 filled, 2 carved | `_conditioned_pixels.tif` |
| `lsdtt3.slope_area_bins(drainage_area, slope, n_log_bins=20, exclude=None)` | DataFrame of binned quantiles | `_lambda_XXm.csv`, `_d8window_XXm.csv` |
| `FlowInfo.d8_flow_path_gradient(routing_dem, drainage_area, window_metres=150, exclude=None)` | Raster of gradient (m/m) | `_d8window_XXm_slope.tif` |
| `FlowInfo.low_concavity(drainage_area, gradient, ...)` | DataFrame, one row per channel pixel | `_low_concavity_points.csv` |
| `FlowInfo.low_concavity_segments(drainage_area, gradient, ...)` | GeoDataFrame of segments | `_low_concavity_segments.fgb` |

## The CLI's default run

The CLI carves, then fills, then routes, and leaves the carved and filled
pixels out. In Python:

```python
import lsdtt3

dem = lsdtt3.Raster.read("my_site_DEM.tif")

# carve then fill (breach_before_fill + fill_before_routing)
carved = lsdtt3.breach_least_cost(dem, fill_remainder=False)
routing = lsdtt3.fill_pits(carved)
conditioned = lsdtt3.conditioned_pixel_mask(dem, routing)   # the exclusion mask

# drainage area: D-infinity (the default) or D8
flow = lsdtt3.FlowInfo(routing)
area = lsdtt3.FlowGraph(routing, method="dinf").drainage_area()

# polyfit gradient on the INPUT DEM, binned by area
polyfit = lsdtt3.polyfit_metrics(dem, ["slope"], lambda_metres=150.0)["slope"]
bins = lsdtt3.slope_area_bins(area, polyfit, n_log_bins=30, exclude=conditioned)

# or the D8 flow-path gradient, on the routing surface
d8 = flow.d8_flow_path_gradient(routing, area, window_metres=150.0, exclude=conditioned)
d8_bins = lsdtt3.slope_area_bins(area, d8, n_log_bins=30, exclude=conditioned)
```

`bins` has the columns `drainage_area`, `d05_slope`, `d25_slope`, `d50_slope`,
`d75_slope`, `d95_slope` and `number_of_pixels`.

- `exclude` takes any raster on the same grid; nonzero pixels are left out.
- Without `exclude` every pixel is used, as with
  `exclude_carved_and_filled_pixels: false`.
- For D8 area use `flow.drainage_area()`.

## Low-concavity reaches

The finder always uses the D8 flow-path gradient:

```python
gradient = flow.d8_flow_path_gradient(routing, area, window_metres=150.0,
                                      exclude=conditioned)

points = flow.low_concavity(area, gradient,
                            threshold_contributing_pixels=100,
                            window_pixels=10, min_pixels=5,
                            max_abs_concavity=0.15, exclude=conditioned)

segments = flow.low_concavity_segments(area, gradient,
                                       threshold_contributing_pixels=100,
                                       min_run_pixels=3, exclude=conditioned)
segments.to_file("my_site_low_concavity_segments.fgb")
```

`points` has one row per channel pixel, with these columns:

- `flat`, `x`, `y`;
- `drainage_area_m2`, `slope`;
- `theta`, `theta_se`, `r2`, `n_pixels`, `log10_area_range`;
- `n_junctions_in_window`;
- `low_concavity`: 1, 0, or -1 for no valid fit.

`segments` has the attributes of the CLI's segments file. Pass
`include=` (nonzero = kept) to keep only the pixels inside a basin, as
the CLI's basin clip does.

The CLI and Python run the same C++ routines. With the same settings,
the two give the same flags and segments.

Figures of these outputs (binned plots, the theta histogram, the slope-area
cloud and the segments map) are drawn by `lsdviz3-slope-area` in lsdviztools3
from the CLI's files.
