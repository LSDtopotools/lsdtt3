# Analysing along a feature: oblique frames

A fault, a trunk river, a coastline or a mountain range rarely runs north–south
or east–west. Analysed on a UTM grid, a long diagonal strip sits in a bounding
box that is mostly not the strip: a 100 km × 10 km strip at 45° fills a
78 km × 78 km box, six times the cells. Orientations come out relative to
north, when what you usually want is relative to the feature.

lsdtt3 can instead put the DEM in an **oblique frame**: a grid whose x axis runs
along the feature. The DEM is still an ordinary north-up GeoTIFF; only its CRS
is different (a Hotine Oblique Mercator projection centred on the feature), so
every lsdtt3 tool reads and analyses it unchanged. When the analysis is done,
products such as basins, channel masks and slope are rotated back onto a UTM
grid. How the frame is built and how well results agree between frames and
UTM is in [Oblique analysis frames](../design_and_algorithms/oblique_frames.md).

This page works through two examples with the command-line tools:

1. a strip along the San Andreas Fault in the Carrizo Plain, California, fetched
   straight into its frame, analysed, rotated back and plotted;
2. the Wind River Range, Wyoming: an existing UTM DEM rotated into the frame of
   a polygon drawn around the range, then profiled along its axis.

The same steps from Python are in [Oblique frames in Python](../python/oblique_frames.md).

---

## What a frame looks like

- **+x runs along the feature** and the frame centre is (0, 0). The axis
  azimuth is stored folded into (−90°, 90°], so for a feature trending
  north-west to south-east (azimuth −45°) +x points north-west and the map has
  its north-west end on the right.
- **Up is not north.** In that −45° frame, grid up points south-west, and true
  north is 135° clockwise from grid up. Tools that write bearings print this at
  start-up:

  ```
  Grid north is 135.4 deg anticlockwise of true north at the grid centre;
  bearings are relative to grid north (bearings_relative_to_true_north is false).
  ```

- **Bearings are relative to the grid** unless you ask otherwise: aspect,
  channel-link bearings, basin orientation and the other bearing outputs are
  clockwise from grid up, so 90° means "along the feature". Set
  `bearings_relative_to_true_north: true` in any tool that writes bearings to
  get them from true north instead; each value is converted at its own
  location. Angles between two directions (junction angles, deflections) are
  the same either way, and hillshade illumination always stays relative to the
  grid (315° lights the map from its top left).
- **A rotation frame file** is written next to every DEM put into a frame
  (`<name>_rotation_frame.txt`). It records the frame and the UTM grid to
  rotate products back onto. Keep it with the DEM: rotate-back reads it.

---

## Example 1: a fault strip

### The feature file

The fault is a line in any vector format and CRS. Here it is an approximate
trace of the San Andreas Fault across the Carrizo Plain, digitised for this
example, saved as `carrizo_fault.geojson`:

```json
{"type": "FeatureCollection", "name": "carrizo_fault",
 "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},
 "features": [{"type": "Feature", "properties": {"name": "San Andreas Fault, Carrizo Plain (approximate)"},
   "geometry": {"type": "LineString", "coordinates": [[-120.02, 35.40], [-119.80, 35.22], [-119.55, 35.02]]}}]}
```

A line (or a multi-part line) is framed by its two furthest end points; a file
with two points uses them as the axis ends; a polygon is framed by its
minimum-area rectangle, with the longest side as the axis (Example 2).

### 1. Fetch the strip into its frame

```
# fetch.param: an 8 km wide strip along the fault, straight into its frame.
fetch_from_linear_feature: true
fetch_feature_fname: carrizo_fault.geojson
linear_feature_fetch_width: 8000
grid_spacing: 30
source: cop30_aws
write_fname: carrizo
```

```bash
lsdtt-fetch-raster fetch.param
```

```
Centre  : -119.7844525, 35.21023104
Azimuth : -45.42894765 deg
Size    : 60720 m x 8040 m (2024 x 268 cells)
Wrote: ./carrizo_oblique_DEM.tif
```

The source tiles are warped straight into the frame, so the DEM is resampled
once. The strip is the line's length plus an end buffer (`linear_feature_fetch_end_buffer`,
default ten cells at each end) by the width. Next to the DEM is
`carrizo_oblique_DEM_rotation_frame.txt`, whose rotate-back target is UTM
zone 11N (EPSG:32611), the zone of the frame centre, on the 30 m lattice a
plain fetch uses. `rotate_after_fetch: false` would instead fetch the strip's
bounding box on that UTM grid.

The frame holds 2024 × 268 = 542,432 cells. The same strip rotated back onto
UTM 11N needs a 1597 × 1648 grid, 2,631,856 cells.

### 2. Analyse in the frame

Any tool runs on the frame DEM as on any projected DEM. Here, complete
second-order basins (in an 8 km strip, third-order basins were all cut by the
strip's edges), and slope and aspect:

```
# basins.param
read_fname: carrizo_oblique_DEM.tif
write_fname: carrizo
threshold_contributing_pixels: 500
select_all_basins_of_stream_order: true
selected_stream_order: 2
filter_incomplete_basins: true
write_basins_raster: true
write_basin_outlines_as_polygon: true
```

```
# metrics.param: aspect is clockwise from grid up, 90 = along the fault.
read_fname: carrizo_oblique_DEM.tif
write_fname: carrizo
polyfit_metrics: slope,aspect_degrees
polyfit_lambda_metres: 90
write_hillshade_raster: true
```

```bash
lsdtt-basin-extraction basins.param      # 65 basins: carrizo_junction_basins.tif + _polygons.fgb
lsdtt-surface-metrics metrics.param      # carrizo_direct_polyfit_slope.tif, _aspect_degrees.tif
```

### 3. Rotate the products back

Rotate-back is a separate mode of `lsdtt-raster-preprocessing`. List the rasters
one per line, with `near` for anything categorical (basin IDs, masks, junction
numbers, lithology); the default method is `bilinear`:

```
# rotate_back_list.txt
carrizo_junction_basins.tif, near
carrizo_direct_polyfit_slope.tif
```

```
# rotate_back.param
run_rotate_back: true
rotation_frame_fname: carrizo_oblique_DEM_rotation_frame.txt
rotate_back_file_list: rotate_back_list.txt
write_fname: carrizo
```

```bash
lsdtt-raster-preprocessing rotate_back.param
```

Each raster is written as `<name>_rotated_back.tif` on EPSG:32611 at 30 m, here
1597 × 1648 cells, with all 65 basin IDs. Vector outputs such as
`carrizo_junction_basin_polygons.fgb` need no rotating back: they carry their
CRS, and QGIS, geopandas and lsdviztools3 reproject them.

**Rotate back products, never the DEM.** Rotating a DEM into a frame and back
resamples it twice: for the 30 m Copernicus DEM of Example 2 (cubic both ways,
5.6 million cells), the round trip changed elevations by 0.36 m RMS, 1.3 m at
the 99th percentile and up to 18 m on cliffs. Keep the original DEM for
anything that needs elevations on the UTM grid.

### 4. Figures

[lsdviztools3](https://lsdtopotools.github.io/lsdviztools3/) draws a frame DEM
as a rectangle in frame km, with axes labelled along and across the feature, a
title giving the frame azimuth, the UTM grid drawn over it as dashed lines with
their eastings and northings, and a north arrow that points to true north at
its own position:

```python
from lsdviztools3.io.raster import get_crs, load_raster, to_coordinates, to_gmt_grid
from lsdviztools3.render.base import grid_to_display
from lsdviztools3.render.decorations import add_north_arrow
from lsdviztools3.workflows import hillshade_with_vectors


def map_with_arrow(base, layers, png, arrow_size=0.07, **kwargs):
    """A base raster + vector layers + a north arrow pointing to true north."""
    fig = hillshade_with_vectors(base, layers, legend=False, **kwargs)
    grid = to_gmt_grid(to_coordinates(load_raster(base), "local"))
    crs = get_crs(grid)
    _, extent = grid_to_display(grid, crs)
    add_north_arrow(fig, extent, crs, anchor="BR", label_halo="white",
                    size_frac=arrow_size)
    fig.savefig(png, dpi=200)


fault = {"path": "carrizo_fault.geojson", "kind": "line", "color": "red", "width": "1.5p"}
basins = {"path": "carrizo_junction_basin_polygons.fgb", "kind": "polygon",
          "color": "blue", "width": "0.6p"}

# The strip in its frame: the fault runs along +x.
map_with_arrow("carrizo_oblique_DEM.tif", [basins, fault], "carrizo_frame_map.png",
               arrow_size=0.25)

# Slope rotated back onto UTM 11N; the basin polygons, still in the frame CRS,
# are reprojected onto the map.
map_with_arrow("carrizo_direct_polyfit_slope_rotated_back.tif", [basins, fault],
               "carrizo_utm_map.png", shade=False, base_cmap="viridis")
```

The graticule is on by default for frame maps and off for all others;
`MapStyle(graticule=False)` turns it off, and `graticule_interval_km` and
`reference_crs` set its spacing and CRS. The north arrow is turned by the grid
convergence on every projected map: by up to about 3° on UTM, by 135° in this
frame.

---

## Example 2: an elongated range from an existing DEM

When the DEM is already on disk, `lsdtt-raster-preprocessing` rotates it into a
frame in place of its usual reprojection step. Here the DEM is a plain UTM fetch
of the Wind River Range, and the frame comes from a polygon drawn around the
range (124 km × 40 km, saved as `wind_river_range.geojson`; any vector format
and CRS):

```json
{"type": "FeatureCollection", "name": "wind_river_range",
 "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},
 "features": [{"type": "Feature", "properties": {"name": "Wind River Range (approximate)"},
   "geometry": {"type": "Polygon", "coordinates": [[[-109.751, 43.507], [-110.148, 43.292],
     [-109.245, 42.392], [-108.854, 42.607], [-109.751, 43.507]]]}}]}
```

### 1. The UTM DEM

```
# fetch.param: an ordinary UTM DEM covering the polygon's bounding box.
bounds_fname: wind_river_range.geojson
grid_spacing: 30
source: cop30_aws
write_fname: windriver
```

```bash
lsdtt-fetch-raster fetch.param    # windriver_DEM.tif: 3593 x 4191 cells, EPSG:32612
```

### 2. Rotate it into the range's frame

```
# rotate.param
read_fname: windriver_DEM.tif
write_fname: windriver
rotate_to_linear_feature: true
rotation_feature_fname: wind_river_range.geojson
write_fill_raster: true
write_hillshade_raster: true
```

```bash
lsdtt-raster-preprocessing rotate.param
```

```
Centre  : -109.4957472, 42.94983547
Azimuth : -36.30276092 deg
Size    : 124920 m x 40680 m (4164 x 1356 cells of 30 m)
Rotate back : EPSG:32612 at 30, lattice origin (568830, 4693650)
```

The rotated DEM is `windriver_oblique.tif` (4164 × 1356 = 5,646,384 cells,
against 15,058,263 for the UTM DEM), and every later step in the same run
(here the fill and the hillshade) uses it. Because the source DEM is projected
in metres, the frame keeps its 30 m cells and the rotate-back target is its own
grid: products rotated back land on the source DEM's pixels. The frame is on
the source DEM's own datum, so a NAD83 DEM gets a NAD83 frame.

Other ways to define the frame:

- `rotate_to_linear_feature` with a line needs `linear_feature_rotation_width`
  (metres); a polygon takes its width from the polygon.
- `rotate_to_data_footprint: true` frames the DEM's valid cells (the
  minimum-area rectangle around them) instead of a feature file, for a DEM
  already clipped to an elongated area. A footprint with no clear long axis is
  warned.
- `linear_feature_rotation_end_buffer` adds length beyond the axis ends (all
  four sides for a polygon or a footprint); the default, −1, is ten cells.

### 3. A swath along the range

In the frame the range axis is the x axis, so a swath baseline along the range
is the line y = 0. `lsdtt-swath-profile` reads the baseline's coordinates as
frame metres:

```json
{"type": "FeatureCollection", "name": "range_axis",
 "features": [{"type": "Feature", "properties": {},
   "geometry": {"type": "LineString", "coordinates": [[-62000.0, 0.0], [62000.0, 0.0]]}}]}
```

```
# swath.param
read_fname: windriver_oblique.tif
write_fname: windriver
baseline_vector_fname: range_axis.geojson
half_width_metres: 15000
bin_width_metres: 2000
write_distance_rasters: false
```

```bash
lsdtt-swath-profile swath.param    # windriver_swath_profile.csv
```

and plot it, with the map of the range in its frame:

```python
from lsdviztools3.charts import plot_swath
from lsdviztools3.io.swath import load_swath

fig, ax = plot_swath(load_swath("windriver_swath_profile.csv"))
ax.set_xlabel("Distance along the range axis (km)")
fig.savefig("windriver_swath.png", dpi=200, bbox_inches="tight")
```

---

## Rules of thumb

- **Analyse in the frame; bring back products, never the DEM.**
- **Bearings are from grid up** unless `bearings_relative_to_true_north: true`.
- **Do not set `force_reproject` on a frame DEM**: it would warp it back to UTM.
  Preprocessing refuses `force_reproject` together with a rotation.
- **Vector inputs in a frame CRS go in FlatGeobuf or GeoPackage.** GeoJSON
  cannot record a custom CRS, so a GeoJSON written in a frame CRS is read as
  WGS84. Tools whose help asks for a vector in the DEM's CRS (a crossing line
  for `lsdtt-basin-extraction`) need it in the frame CRS;
  `lsdtt-swath-profile` reads its baseline's coordinates as-is, which is why
  the GeoJSON axis above works.
- **The strip's long sides are grid edges**, so boundary conditions can be set
  per side, for example no flux across the strip and open at its ends.
- **Nodata wedges need no special handling.** Where the source DEM is smaller
  than the frame, the corners are nodata, and the fills drain cells beside
  nodata as at any edge.
- **Lost the frame file?** `run_rotate_back` without `rotation_frame_fname`
  rebuilds the frame from each raster's own CRS and rotates back to the UTM
  zone of the frame centre on the lattice anchored at (0, 0), with a warning: a
  source DEM origin off that lattice cannot be recovered without the file.

All the parameters are in the [lsdtt-fetch-raster](../cli_reference/cli_reference_lsdtt_fetch_raster.md)
and [lsdtt-raster-preprocessing](../cli_reference/cli_reference_lsdtt_raster_preprocessing.md)
references.
