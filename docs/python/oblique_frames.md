# Oblique frames in Python

An oblique frame is a grid whose x axis runs along a feature (a fault, a river,
a range), so a long diagonal strip is analysed on a grid that holds the strip
and not its bounding box. The concepts, and the same work with the
command-line tools, are in
[Analysing along a feature](../design_and_algorithms/oblique_frames.md), and
the method at the end of that page. The
Python calls run the same C++ routines as the command-line tools, so the same
frame and settings give the same rasters.

## A fault strip, start to finish

The feature is an approximate trace of the San Andreas Fault across the
Carrizo Plain, California, saved as `carrizo_fault.geojson` (any vector format
and CRS; the file is listed on the command-line page).

```python
import lsdtt3

# 1. The frame of an 8 km wide strip along the fault, and the DEM fetched
#    straight into it (one resampling pass from the source tiles).
frame = lsdtt3.RotationFrame.from_feature("carrizo_fault.geojson", width=8000.0,
                                          grid_spacing=30.0)
print(frame)
# RotationFrame(line, azimuth -45.428948 deg, 2024 x 268 cells of 30.000000 m)
print(frame.target_epsg)          # 32611: products rotate back onto UTM 11N
frame.write("carrizo_rotation_frame.txt")

dem = lsdtt3.fetch_dem(frame=frame, grid_spacing=30.0)

# 2. Analysis in the frame: the fault runs along +x.
filled = lsdtt3.fill_pits(dem)
flow = lsdtt3.FlowInfo(filled)
network = flow.channel_network(flow.channel_scope(threshold_contributing_pixels=500))
links = network.link_bearings()          # bearing_deg clockwise from grid up
links_true = network.link_bearings(bearings_relative_to_true_north=True)

junctions = network.junctions()
outlets = junctions[junctions.stream_order == 2][["row", "col"]].to_numpy()
basins = flow.basins(outlets)
labels = basins.labels                   # basin_id + 1, 0 outside every basin
slope = lsdtt3.slope(filled)

# 3. Bring the products back onto UTM: IDs with "near", values with "bilinear".
labels_utm = lsdtt3.rotate_back(labels, frame, resampling="near")
slope_utm = lsdtt3.rotate_back(slope, frame)
print(labels_utm.shape)                  # (1648, 1597) cells of 30 m on EPSG:32611
labels_utm.write("carrizo_basins_utm.tif")
```

`fetch_dem` always returns a `Raster`; keep the `RotationFrame` to rotate
products back. `fetch_dem(fetch_from_linear_feature=True,
fetch_feature_fname=..., linear_feature_fetch_width=...)` builds the same frame
from the command-line tool's keywords, but does not return it, so build it with
`RotationFrame.from_feature` when you will rotate back.

## Rotating a DEM you already have

```python
dem = lsdtt3.Raster.read("windriver_DEM.tif")       # a UTM DEM

frame = lsdtt3.RotationFrame.from_feature("wind_river_range.geojson", source=dem)
rotated = lsdtt3.rotate_to_frame(dem, frame)         # resampling="cubic" by default

# ... analysis on `rotated` ...

back = lsdtt3.rotate_back(product, frame, resampling="near")
```

With `source=`, the frame is built as `lsdtt-raster-preprocessing` builds it:

- on the DEM's own geographic datum (a NAD83 DEM gets a NAD83 frame);
- at the DEM's cell size when it is projected in metres (`grid_spacing=0`, the
  default), or at the cell size reprojection would choose for a geographic DEM;
- with the DEM's own grid as the rotate-back target, so products land on its
  pixels.

Without `source=`, the frame is for `fetch_dem`: on WGS84, at
`grid_spacing` (30 m when 0), rotating back to the UTM (or polar stereographic
beyond 75°) zone of the frame centre, or to `target_epsg`.

`RotationFrame.from_valid_data(dem)` frames the minimum-area rectangle of the
DEM's valid cells instead of a feature file.

## Frame files

`RotationFrame.write` and `RotationFrame.read` use the same
`<name>_rotation_frame.txt` file as the command-line tools, so a frame made in
Python can rotate products back with `lsdtt-raster-preprocessing`
(`run_rotate_back`) and the other way round. If the file is lost,
`RotationFrame.from_raster(raster)` rebuilds the frame from a frame raster's own
CRS and warns that the rotate-back target falls back to the UTM zone of the
frame centre on the lattice anchored at (0, 0).

## Properties

| Property | Meaning |
|---|---|
| `azimuth_degrees` | Axis azimuth at the centre, clockwise from true north, in (−90, 90] |
| `centre_longitude`, `centre_latitude` | Frame centre, in the base geographic CRS |
| `axis_length_metres`, `width_metres`, `end_buffer_metres` | Frame dimensions |
| `cell_size`, `n_cols`, `n_rows`, `extent` | The frame grid; `extent` is (x_min, x_max, y_min, y_max) in frame metres |
| `crs_wkt`, `proj_string`, `base_crs_wkt` | The frame CRS and its geographic CRS |
| `frame_source` | `line`, `two_points`, `polygon`, `data_footprint` or `crs` |
| `target_epsg`, `target_crs_wkt`, `target_cell_size`, `target_origin` | The rotate-back grid |

`frame.grid_convergence_degrees(x, y)` gives the angle of true north clockwise
from grid up at a point (frame metres): a grid bearing g there is the true
bearing g − c.

## Bearings

Every bound routine that returns a bearing, aspect or orientation takes
`bearings_relative_to_true_north` (default `False`): `aspect`,
`polyfit_metrics`, `ChannelNetwork.link_bearings`,
`ChannelNetwork.basin_com_bearings`, `BasinCollection.statistics`
(`orientation_degrees`), `basin_surface_metrics`, `find_roads`,
`find_road_flow_discrimination` and `detect_directional_periodicity`. With
`False` the values are clockwise from grid up, so in a frame 90° is along the
feature; with `True` each value is converted with the grid convergence at its
own location. Angles between two directions (`junction_angles`) are unchanged.
