# Ridgeline profiles

A `Ridgeline` is an ordered line of pixels along a ridge or divide. It holds
the ridge's elevation profile, and its methods fit and detrend that profile
as `lsdtt-ridgeline-analysis` does. The pixels can come from a divide you
traced yourself, a basin outline, or any other source.

```python
import numpy as np
import lsdtt3

dem = lsdtt3.Raster.read("dem.tif")

# N x 2 (row, col) pixels along the ridge, in order, four-connected.
ridge = lsdtt3.Ridgeline(dem, ridge_pixels)
ridge.fill_elevation(dem)

ridge.segment_snout(smoothing_window_fraction=0.1)
ridge.fit_profile(detrend_moving_average_window_pixels=[11, 31])

print(ridge.snout_index, ridge.slope_fit, ridge.r_squared_fit)
fitted, detrended = ridge.moving_average_fits[31]
```

Call the methods in this order:

1. `fill_elevation(dem)` reads the profile from the DEM. Before it, the
   elevations are empty.
2. `segment_snout` finds the snout, the steep end where the ridge drops into
   a valley, following Robinson et al. (2025). The profile is LOWESS
   smoothed over `smoothing_window_fraction` of its length, and the snout is
   cut at the curvature maximum. The CLI's default fraction
   (`ridgeline_smoothing_window_fraction`) is 0.05. `snout_index` is the
   first node after the snout. Every fit below uses only the nodes from
   `snout_index` on. Without this step `snout_index` is 0 and the whole
   profile is used.
3. `fit_profile` fits a straight line to elevation against distance along
   the ridge. `slope_fit` (m/m) and `r_squared_fit` always come from that
   line. With `use_polynomial=True` the profile is detrended against a
   quadratic instead, and with `use_constant=True` against its mean (the
   CLI's `use_constant_fit_for_mdd_ridgeline`, for main drainage divides).
   The statistics stay those of the line.

`slope_fit` is signed in the order of the pixels: it is negative when
elevation falls along the ridge. `reverse()` flips the order. Distances
along the ridge are in `distance_metres`, and the total in `length_metres`.

## Moving-average detrending

`detrend_moving_average_window_pixels` adds one detrended profile per window
width, each against a moving average of that many pixels instead of the
fitted line. Widths must be positive, odd and given once each. This is the
CLI parameter `ridgeline_detrend_moving_average_window_pixels`. A short
window keeps only the bumps narrower than itself; a long one keeps broader
undulations.

`moving_average_fits` maps each width to a pair of arrays, the moving
average and the elevation minus it, one value per node. Nodes before
`snout_index` hold the DEM's nodata value, so mask them before plotting:

```python
fitted, detrended = ridge.moving_average_fits[31]
fitted = np.asarray(fitted)
valid = np.arange(len(fitted)) >= ridge.snout_index
distance = np.asarray(ridge.distance_metres)[valid]
```

## Ridgeline steepness

`detrend(smoothing_window_pixels)` smooths the detrended profile, as the
CLI does before it looks for maxima and minima. `steepness(area, concavity)`
normalises the fitted slope by the drainage area of the basin the ridge
bounds: `slope_fit / area**(-concavity)`. The area is in m^2, and the
concavity is the same theta used for ksn (0.45 is common).

## Ridgelines from basin outlines

`basin_ridgeline_extrema(basins, network, dem)` splits each basin outline
into ridgelines and runs the snout cut, the fit and the smoothing on each.
It returns the maxima and minima along every ridge, with their fit and
topological-distance diagnostics, as a GeoDataFrame. See
[Flow, channels, and basins](flow_networks.md) for making `basins` and
`network`.

The [lsdtt-ridgeline-analysis reference](../cli_reference/cli_reference_lsdtt_ridgeline_analysis.md)
covers the CLI's ridge types and its outputs. The
[API reference](../python_api/sediment.md) lists every `Ridgeline` method.
