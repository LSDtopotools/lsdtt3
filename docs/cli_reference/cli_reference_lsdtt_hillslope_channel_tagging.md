# lsdtt-hillslope-channel-tagging

!!! warning "Experimental"
    This CLI is under development. Its routines are tested against analytic
    and hand-built cases, but its options and output names may still change.
    The method is in
    [Hillslope-channel tagging](../design_and_algorithms/hillslope_channel_tagging.md),
    and the [Python guide](../python/hillslope_channel_tagging.md) shows the same
    steps in Python.

## Intent

`lsdtt-hillslope-channel-tagging` transposes channel information onto hillslope
pixels. Each hillslope pixel is linked to the first channel pixel on its D8
path, and that link is the same whichever ridge method is chosen. From it the
tool reports:

- **HAND**, the height above the nearest drainage (Nobre et al. 2011);
- **ridge relief**, the elevation of the pixel's ridge above its channel pixel;
- **relative height**, HAND divided by ridge relief;
- **tags**, the value of the channel pixel copied to its hillslope pixels: chi,
  k_sn, or an imported attribute such as sinuosity or planform curvature.

## Usage

```bash
lsdtt-hillslope-channel-tagging <path> <parameter_file>
lsdtt-hillslope-channel-tagging <parameter_file>
lsdtt-hillslope-channel-tagging -h
```

## Choosing a ridge

`hillslope_ridge_method` selects how the ridge is found:

| Method | Ridge |
|---|---|
| `farthest` (default) | The pixel with the longest D8 path to the channel, in the pixel's group |
| `highest` | The highest pixel in the pixel's group |
| `aspect_ray` | Where a straight ray, running upslope along the pixel's aspect, stops at a drainage divide |

**Groups.** About half of the hillslope pixels of a typical DEM receive no flow
from any neighbour, so a ridge searched only among a pixel's own donors would
often be the pixel itself, with a relative height of 1. `farthest` and
`highest` therefore share a ridge within a group, set by
`hillslope_relief_grouping`: `channel_pixel` (default), `entry_pixel` (pixels
entering the channel through the same last hillslope pixel) or `pixel`.

**Relative height above 1.** With `farthest` and a shared group, the farthest
pixel need not be the highest, so a pixel on a shorter, steeper branch can
stand above its ridge. These values are reported as they are, and counted in the
run summary.

**Aspect rays.** The divides are the boundaries between the catchments of the
channel links, with no filtering. A ray stops at the first divide that is at or
above its starting pixel and is a crest along the ray, meaning the next cell is
not higher. If the ray reaches a channel, nodata or the grid edge first, the
ridge is the highest cell on the ray. A pixel with no aspect has no ridge. The
crest test looks one cell ahead, so a small bump on a rising slope can stop a
ray early.

## Chi and k_sn tags

`tag_chi` and `tag_k_sn` need a true baselevel, and a DEM's nodata cannot show
whether a basin ends at the coast or is cut off by the edge of the data. Declare
the valid basins with exactly one of:

- `basin_outlet_fname`: these basins;
- `chi_max_outlet_elevation_metres`: every baselevel basin whose outlet is at or
  below this elevation. On a DEM with the sea masked, a few metres keeps the
  basins that reach the coast.

Basins with a pixel on the grid edge are always excluded. Hillslopes of excluded
basins are nodata and counted. Check this count, because the limit can exclude a
large share of a DEM.

Chi is fitted on the network given by `threshold_contributing_pixels`. If you
ingest a different network with `channel_network_fname`, channel pixels outside
the fitted network have no chi or k_sn, and their hillslopes are counted as
having no value.

## Imported channel attributes

`channel_attributes_fname` is a CSV or vector file of channel points or lines.
`channel_attribute_columns` lists the numeric columns to tag.

| Input | How it joins |
|---|---|
| `row`, `col` and `flat` columns | Exactly, as written by `write_channel_network_points` |
| CSV with `latitude` and `longitude` | Snapped to the nearest channel pixel |
| Points in a vector file | Reprojected, then snapped |
| Lines in a vector file | Reprojected, sampled every half cell, and each sample snapped |

- Records further than `channel_attribute_max_snap_distance_metres` from every
  channel pixel are dropped and counted.
- Several records on one channel pixel give their median, ignoring missing
  values. A line counts once per pixel.
- A channel pixel with no record takes the nearest record on its own
  junction-to-junction link, up to
  `channel_attribute_max_gap_fill_distance_metres` along the channel. Filling
  never crosses a junction, and equidistant records give their median.
- Empty cells, NaN and `channel_attribute_nodata_value` are missing values.

## Outputs

| File | Written by | Content |
|---|---|---|
| `<write_fname>_HAND` | `write_hand_raster` | Height above the nearest drainage (m) |
| `<write_fname>_ridge_relief` | `write_ridge_relief_raster` | Ridge elevation above the channel pixel (m) |
| `<write_fname>_relative_height` | `write_relative_height_raster` | HAND / ridge relief |
| `<write_fname>_ridge_distance` | `write_ridge_distance_raster` | D8 path from ridge to channel, or for `aspect_ray` the straight line from pixel to ridge (m) |
| `<write_fname>_ridge_status` | `write_ridge_status_raster` | 0 channel, 1 no channel, 2 D8, 3 crest on a divide, 4 ray reached a channel, 5 ray reached nodata or the edge, 6 no aspect |
| `<write_fname>_tag_chi`, `_tag_k_sn` | `tag_chi`, `tag_k_sn` | Chi and k_sn of the channel pixel |
| `<write_fname>_tag_<column>` | `channel_attribute_columns` | The imported attribute of the channel pixel |
| `<write_fname>_channel_pixel_summary.csv` | `write_channel_pixel_summary` | One row per channel pixel: hillslope pixels and median relative height |

Channel pixels carry their own value in the tag rasters. Undefined values are
nodata, never 0.

## Example parameter file

```yaml
read_fname: my_dem.tif
write_fname: my_run
threshold_contributing_pixels: 500

hillslope_ridge_method: aspect_ray
polyfit_lambda_metres: 70

tag_chi: true
tag_k_sn: true
m_over_n: 0.45
chi_max_outlet_elevation_metres: 5

channel_attributes_fname: my_channel_planform.fgb
channel_attribute_columns: sinuosity,planform_curvature
channel_attribute_max_snap_distance_metres: 60
channel_attribute_max_gap_fill_distance_metres: 150
```

## Reference

Nobre, A. D., et al. (2011). Height Above the Nearest Drainage: a hydrologically
relevant new terrain model. *Journal of Hydrology* 404, 13-29.

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

### Flow routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string for north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. A positive value is honoured verbatim and OVERRIDES OMP_NUM_THREADS, so setting 1 here runs serial on any machine. |

### Input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_elevation` | `float` | `-9999` | `param` | - | Cells at or below this elevation (m) become nodata before any analysis. -9999 disables. | For removing the sea, which many DEMs store as valid cells at about 0 m, and nodata stored as 0 (a block of zeros beside real terrain is a cliff that captures flow paths). Removed cells that touch the DEM edge act as an open boundary, so a river reaching the coast drains off the DEM there. Negative values apply as given: where land lies below sea level (fens, polders, Death Valley -86 m, the Dead Sea -430 m) set it just below the lowest land, or below the sea surface of the DEM, so that flow is still routed over the sub-zero land. -9999 (or lower) disables. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `breach_before_fill` | `bool` | `true` | `switch` | - | Breach (carve) depressions before filling pits. | The carve comes first, then the fill (fill_before_routing) takes only the pits the carve could not solve, so large filled areas (and their artificially low slopes) are avoided. SINMAP itself only fills: false reproduces that. |
| `breach_algorithm_choice` | `string` | `least_cost` | `param` | `breach_before_fill` | Algorithm for depression breaching. | lindsay2016: Lindsay (2016) efficient breaching. least_cost: Lindsay (2020) pit-centric least-cost breaching. |
| `breach_lc_max_dist` | `int` | `1000 (outlet) / 100 (whitebox)` | `param` | `breach_before_fill` | Maximum breach path length in cells for least-cost breaching. | Counted in steps from the pit. Pits that need a longer path are left unbreached (and filled if fill_before_routing is true). 1 to 32767. If not set, 1000 with breach_lc_edge_policy outlet and 100 with whitebox. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_max_cost` | `float` | `unlimited` | `param` | `breach_before_fill` | Maximum cost of a least-cost breach path. | The summed depth each cell on the path is lowered (m; m x cells with breach_lc_minimize_dist). A pit whose cheapest path costs more is left unbreached. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_minimize_dist` | `bool` | `false` | `param` | `breach_before_fill` | Weight least-cost breach paths by their length. | false: a path costs the total depth lowered (m). true: each step is weighted by its length in cells. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_edge_policy` | `string` | `outlet` | `param` | `breach_before_fill` | Whether a least-cost breach may drain off the raster edge. | outlet: the raster edge is an outlet. whitebox: it is not (as in WhiteboxTools). Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_flat_increment` | `float` | `0.0` | `param` | `breach_before_fill` | Elevation drop (m) per cell along a least-cost breach channel. | 0 or less selects the automatic value (twice the float32 spacing at the largest absolute elevation). Only used when breach_algorithm_choice is least_cost. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_before_routing` | `bool` | `true` | `switch` | - | Choose whether to pit-fill the DEM before building the flow network. | ASKS: should the DEM be filled at all? Filling is OPTIONAL here and this is a genuine analytical choice. True (the default) routes over a pit-free surface. False builds the flow network on the UNFILLED surface, so flow paths terminate at sinks instead of being carried through them. Do not confuse this with raster_is_filled, which is an assertion that the input is already filled. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill across filled depressions and flats, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |
| `channel_network_fname` | `string` | - | `param` | - | Optional channel network to ingest instead of using the threshold. | A point file written by lsdtt-channel-extraction (write_channel_network_points), with row and col. It must have been extracted on this DEM's grid. |

### Ridges

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `hillslope_ridge_method` | `string` | `farthest` | `param` | - | How each hillslope pixel's ridge is chosen: farthest, highest or aspect_ray. | farthest: the pixel with the longest D8 path to the channel in the pixel's group. highest: the group's highest pixel. aspect_ray: a straight ray from the pixel, upslope along its aspect, stopping at the first drainage divide that is at or above the pixel and is a crest along the ray; a ray reaching a channel, nodata or the grid edge first uses the highest cell on the ray. The channel a pixel is linked to is the D8 one for all three. With farthest and a shared group, relative height can exceed 1: it is not clipped. |
| `hillslope_relief_grouping` | `string` | `channel_pixel` | `param` | `hillslope_ridge_method` | The pixels sharing a ridge, for farthest and highest. | channel_pixel: all pixels draining to one channel pixel. entry_pixel: all pixels entering the channel through the same last hillslope pixel. pixel: the pixel and those draining through it; about half of the hillslope pixels of a DEM have no donors, so these are their own ridge with relative height 1. Not used by aspect_ray. |
| `hillslope_aspect_method` | `string` | `polyfit` | `param` | `hillslope_ridge_method` | Aspect used for aspect_ray: polyfit or dinf. | polyfit: from the fitted polynomial surface, at polyfit_lambda_metres. dinf: the D-infinity facet direction, which is noisier. |
| `polyfit_lambda_metres` | `float` | `70.0` | `param` | `hillslope_ridge_method` | Polyfit length scale (half-wavelength) for the aspect of aspect_ray (m). | The fitting window follows from this half-wavelength, as in lsdtt-surface-metrics. 0 or less uses three cells. |
| `polyfit_on_routing_surface` | `bool` | `false` | `switch` | `hillslope_ridge_method` | Fit the aspect to the routing surface rather than the input DEM. | The routing surface has carved and filled areas whose slopes come from the conditioning rather than the landscape. |

### Chi and k_sn tags

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `tag_chi` | `bool` | `false` | `switch` | - | Tag hillslope pixels with the chi of their channel pixel. | Needs basin_outlet_fname or chi_max_outlet_elevation_metres, because chi needs a true baselevel. Output: &lt;write_fname&gt;_tag_chi. |
| `tag_k_sn` | `bool` | `false` | `switch` | - | Tag hillslope pixels with the k_sn of their channel pixel. | Needs basin_outlet_fname or chi_max_outlet_elevation_metres, because chi needs a true baselevel. Output: &lt;write_fname&gt;_tag_k_sn. |
| `m_over_n` | `float` | `0.45` | `param` | `tag_chi` | Concavity index for chi. | - |
| `chi_A0` | `float` | `1.0` | `param` | `tag_chi` | Reference drainage area A0 for chi (m2). | Leave at 1 so that the slope of a chi profile is k_sn. |
| `fit_critical_divergence_metres` | `float` | `5.0` | `param` | `tag_chi` | Critical elevation divergence (m) used to split chi-elevation branches. | - |
| `chi_max_outlet_elevation_metres` | `float` | `not set` | `param` | `tag_chi` | Use every baselevel basin with its outlet at or below this elevation (m). | An alternative to basin_outlet_fname; give one of them. On a DEM with the sea masked, a few metres keeps the basins that reach the coast. A basin that a clipped DEM cuts off has a high outlet and is left out. Basins touching the grid edge are always left out. Their hillslopes are nodata and counted. |

### Basin outlets

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_fname` | `string` | `NULL` | `param` | - | File of requested basin outlet points (CSV, FlatGeobuf or GeoPackage). | Format is auto-detected from the file extension: *.csv is read as a text file with latitude and longitude columns (WGS84 / EPSG:4326 decimal degrees); *.fgb or *.gpkg is read as a point vector layer in its own CRS, which is then transformed to the raster CRS. For CSV inputs the file must contain columns named exactly 'latitude' and 'longitude'. For vector inputs all features must be point geometry. |
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | How outlet points are snapped onto the channel network. | Allowed values: none, expanding_square, downstream, largest_area. none: use the mapped raster pixel directly. expanding_square: search outward ring-by-ring until a pixel reaching basin_outlet_contributing_pixels_threshold is found. downstream: follow the D8 receiver chain until a qualifying pixel is found. largest_area: take the pixel with the most contributing pixels in the square of half-width basin_outlet_snap_max_search_radius_pixels, if it reaches basin_outlet_contributing_pixels_threshold (a large radius can jump onto a bigger neighbouring river). |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Contributing-pixel threshold a pixel must reach (&gt;=) to accept an outlet. | This is the total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Used only when basin_outlet_snapping_method is not none. It is a separate quantity from threshold_contributing_pixels, which sets where the channel network begins. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | - | Minimum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. 0 means the mapped input pixel is tested first. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square or largest_area. For expanding_square the search stops at the first radius containing one or more qualifying pixels; for largest_area it is the half-width of the square searched. |

### Imported channel attributes

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `channel_attributes_fname` | `string` | - | `param` | - | CSV or vector file of channel points or lines with numeric attributes. | With row, col and flat columns the join is exact. Otherwise CSV needs latitude and longitude (EPSG:4326), and points or lines in a vector file are reprojected and snapped to the nearest channel pixel. |
| `channel_attribute_columns` | `string` | - | `param` | `channel_attributes_fname` | Comma-separated numeric columns to tag onto hillslopes. | Each gives one raster, &lt;write_fname&gt;_tag_&lt;column&gt;. A missing or text column stops the run. Several records on a channel pixel give their median. |
| `channel_attribute_max_snap_distance_metres` | `float` | `60.0` | `param` | `channel_attributes_fname` | Furthest a record can be from the channel pixel it joins (m). | Records further from every channel pixel are dropped and counted. |
| `channel_attribute_max_gap_fill_distance_metres` | `float` | `150.0` | `param` | `channel_attributes_fname` | Furthest along the channel a record's value can fill a gap (m). | A channel pixel with no record takes the nearest record on its own junction-to-junction link by D8 channel distance. Filling never crosses a junction, and equidistant records give their median. 0 turns filling off. |
| `channel_attribute_nodata_value` | `float` | `-9999` | `param` | `channel_attributes_fname` | Attribute value meaning 'no data'. | Empty and NaN values are always nodata as well. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format: cog, gtiff, png or same_as_input. | cog is a Cloud Optimised GeoTIFF and is the project default. It falls back to gtiff for small rasters, where the tiling overhead is not worth it. There is ONE raster format parameter per CLI: it governs every raster the tool writes. Vector outputs are governed separately by write_vector_format. |
| `write_hand_raster` | `bool` | `false` | `switch` | - | Write height above the nearest drainage (m). | Output: &lt;write_fname&gt;_HAND. |
| `write_ridge_relief_raster` | `bool` | `true` | `switch` | - | Write relief from channel to ridge (m). | Output: &lt;write_fname&gt;_ridge_relief. |
| `write_relative_height_raster` | `bool` | `true` | `switch` | - | Write HAND / ridge relief. | Output: &lt;write_fname&gt;_relative_height. |
| `write_ridge_distance_raster` | `bool` | `false` | `switch` | - | Write ridge distance (m): the D8 path from ridge to channel, or for aspect_ray the straight line from pixel to ridge. | Output: &lt;write_fname&gt;_ridge_distance. |
| `write_ridge_status_raster` | `bool` | `false` | `switch` | - | Write how each ridge was found (0 channel, 1 no channel, 2 D8, 3 crest on a divide, 4 ray reached a channel, 5 ray reached nodata or the edge, 6 no aspect). | Output: &lt;write_fname&gt;_ridge_status. |
| `write_channel_pixel_summary` | `bool` | `false` | `switch` | - | Write a CSV with one row per channel pixel: hillslope pixels and median relative height. | Output: &lt;write_fname&gt;_channel_pixel_summary.csv. |

<!-- END GENERATED CLI PARAMETERS -->
