# Flow, channels, and basins

For D8 channel and basin work, create one `FlowInfo` and reuse it. This keeps
flow directions, drainage area, channel scope, and network indexing coherent.

```python
import lsdtt3

dem = lsdtt3.fill_pits(lsdtt3.Raster.read("dem.tif"))
flow = lsdtt3.FlowInfo(dem)

# A simple area-threshold channel definition.
scope = flow.channel_scope(threshold_contributing_pixels=500)
network = flow.channel_network(scope)

area = flow.drainage_area()
channels = network.lines_geodataframe()
basins = lsdtt3.junction_basins(network)
```

The threshold definition counts every traced pixel above the threshold as a
channel pixel. To keep only the pixels of Strahler order `min_stream_order` or
higher, pass `channel_definition_algorithm_choice="stream_order"`, the same
choice that `lsdtt-fluvial-surfaces` offers. The threshold still decides which
first-order channels exist:

```python
scope = flow.channel_scope(threshold_contributing_pixels=500,
                           channel_definition_algorithm_choice="stream_order",
                           min_stream_order=3)
```

The default `min_stream_order` is 3, the order Clubb et al. (2017) use.

`junction_basins` selects basins at channel junctions. To get the basins
upstream of outlets you choose, pass their pixels to `flow.basins` as an N by 2
array of `(row, col)`:

```python
import numpy as np

outlet_pixels = np.array([[1358, 483], [802, 158]])
basins = flow.basins(outlet_pixels)

basins.table      # basin_id, outlet row/col, n_pixels, area_metres2, complete, ...
basins.polygons   # basin outlines as a GeoDataFrame
basins.labels     # Raster: basin_id + 1 inside each basin, 0 outside
```

Check the `complete` column: a basin cut by the edge of the DEM or by nodata
(such as the sea) is incomplete, and its chi profile or basin-average values
would be wrong.

An outlet is used exactly where it is given. A point taken from a map is
rarely on the D8 channel pixel, so snap it onto the channel first with
`flow.snap_points`, which also converts latitude and longitude (or x and y)
to pixels:

```python
outlet_pixels = flow.snap_points(latitude=[37.94234], longitude=[15.79008],
                                 method="downstream",
                                 threshold_contributing_pixels=500)
```

`method="downstream"` follows the D8 receivers until the contributing-pixel
count reaches the threshold. The methods (`downstream`, `expanding_square`,
`largest_area`) and their rules are those of the command-line programs'
`basin_outlet_snapping_method`, and the same C++ routines do the work. The
`python_binding_examples` notebook in lsdtt3_notebooks has a complete
example.

Use `extract_channel_heads` when a simple accumulation threshold is not the
right source-detection method. The resulting heads can be passed to
`flow.scope_from_heads`, then into the same `channel_network` workflow.

`FlowGraph` is complementary: use it for distributed MFD or D-infinity
products such as TWI or sediment connectivity, rather than D8 channel trees.

See [Flow, channels, and basins API](../python_api/flow.md).
