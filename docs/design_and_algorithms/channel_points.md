# Field points on the channel network

Field data such as grain-size counts, sample sites or gauging stations sit
near a channel but rarely exactly on a DEM's channel cell.
`lsdtt-channel-extraction` can attach channel properties to such points:

- `run_channel_point_attributes` moves each point to the nearest channel
  cell. It then adds that cell's drainage area, flow distance, elevation and
  along-channel slope, and keeps every attribute the point already had.
- `write_channel_slope_points` writes every channel cell with the same slope
  and smoothed elevation. It gives the whole network that the field points
  are compared against.

Both use one definition of the channel slope, the same one
[`lsdtt-bankfull-width`](../cli_reference/cli_reference_lsdtt_bankfull_width.md)
uses.

---

## A parameter file

```text
# channel_points.param
read_path: /data/my_site/
read_fname: my_site_dem.tif
write_path: /data/my_site/results/
write_fname: my_site

# Channel cells drain at least this many pixels
threshold_contributing_pixels: 500

# The field points: csv with latitude/longitude, or a gpkg/fgb/shp/geojson
run_channel_point_attributes: true
channel_points_fname: field_sites.csv
channel_points_snap_max_search_radius_pixels: 10

# The slope window: +/- this distance along the channel
slope_window_half_length_metres: 150

# Every channel cell too, for comparison
write_channel_slope_points: true

write_vector_format: csv
```

Run it with:

```bash
lsdtt-channel-extraction channel_points.param
```

`channel_points_fname` is relative to `read_path` unless it is an absolute
path. A csv needs `latitude` and `longitude` columns in WGS84 degrees, or
`easting` and `northing` in the DEM's coordinate system. Any other format
is read through GDAL and reprojected to the DEM's coordinate system. For a
GeoPackage with several layers, choose one with `channel_points_layer`.

## How a point is moved

Each point moves to the channel cell whose centre is closest to it in a
straight line. A channel cell is one with at least
`threshold_contributing_pixels` contributing pixels. The search goes up to
`channel_points_snap_max_search_radius_pixels` cells from the point.

A point does **not** follow its D8 flow path. A site on a channel bank is
put on that channel, not wherever its hillslope drains. A hillslope can
drain down a tributary and join the channel some distance downstream of the
site.

A point with no channel cell within the radius keeps its own position. Its
`snap_status` says why:

- `threshold_not_found_in_search_window`: no channel cell within the radius;
- `input_out_of_bounds`: the point is outside the DEM.

Its added values are left empty.

## What is added

Every input attribute is kept, in its original order. In a csv output, each
value is written exactly as it was read, so an identifier such as `01` stays
`01`. The added fields follow:

| Field | Meaning |
|---|---|
| `original_x`, `original_y` | the point's position before snapping, in the DEM's coordinate system |
| `snap_distance_m` | how far it moved (m) |
| `snap_status` | `success`, or why it did not snap |
| `snapped_flat`, `contributing_pixels` | the channel cell and its upslope pixel count |
| `drainage_area_m2` | drainage area of that cell (m²) |
| `flow_distance_m` | distance along the flow path to the outlet (m) |
| `elevation_m` | elevation of the pit-filled DEM at the cell (m) |
| `slope_m_per_m` | along-channel slope (m/m) |
| `smoothed_elevation_m` | the fitted line's elevation at the cell (m) |
| `slope_valid` | 1 if the slope fit was usable, else 0 |
| `slope_window_length_m` | the span along the channel the fit used (m) |

If an input column already has an added field's name, for example
`slope_m_per_m`, the added field gets the suffix `_lsdtt` and the log says
so. The input column is not changed.

The geometry of the output is the centre of the snapped cell. The output is
`<write_fname>_channel_points.<ext>`. Its format is
`channel_points_write_vector_format`, or `write_vector_format` if that is
not set: csv (latitude/longitude), fgb or gpkg.

## The channel slope

The slope is the least-squares slope of elevation against distance along
the channel, over ± `slope_window_half_length_metres` around the cell:

- **Downstream**, the window follows the flow path.
- **Upstream**, at each junction it follows the tributary with the larger
  drainage area. It only enters cells that are themselves channel cells.
- **The smoothed elevation** is the fitted line's value at the cell.

Near a channel head or an outlet the window is one-sided. It is still used
if it spans at least the half-length. Otherwise `slope_valid` is 0 and the
slope is empty.

The half-length must be at least two cell sizes. Roughly 5–10 channel widths
smooths out DEM noise without averaging across real changes in gradient.

The elevations come from the pit-filled DEM used for flow routing. Where
filling has flattened the channel, for example behind a road embankment or
in a lake, the slope there is the small gradient filling imposes on flats
(`fill_min_slope`, 0.0001 by default), not the true channel gradient. A
cluster of slopes at exactly that value marks such reaches.

## Plotting

[lsdviztools3](https://lsdtopotools.github.io/lsdviztools3/channel_networks/) draws both
outputs together:

```bash
lsdviz3-channel-slopes --slope-points my_site_channel_slope_points.csv \
    --field-points my_site_channel_points.csv --field-color-by d50_mm \
    --dem my_site_dem.tif --out-dir figures
```

It draws four figures:

- long profiles, raw and smoothed;
- slope against flow distance;
- slope against drainage area;
- a map coloured by slope.

The field points are drawn on top, coloured by one of their own attributes.

## From Python

`FlowInfo.channel_point_attributes` calls the same C++ routine, with the
same parameter names. It takes and returns a GeoDataFrame:

```python
import geopandas as gpd
import lsdtt3

dem = lsdtt3.Raster.read("my_site_dem.tif")
filled = lsdtt3.fill_pits(dem)
flow = lsdtt3.FlowInfo(filled)
sites = gpd.read_file("field_sites.gpkg")
out = flow.channel_point_attributes(
    filled, sites,
    threshold_contributing_pixels=500,
    slope_window_half_length_metres=150.0,
    channel_points_snap_max_search_radius_pixels=10,
)
```

The result is in the DEM's coordinate system, with the columns above. Empty
values are `NaN`, or `<NA>` for the integer columns.
