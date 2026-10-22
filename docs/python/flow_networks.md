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
rarely on the D8 channel pixel, so move it downstream onto the channel first:
follow `flow.receiver_of(row, col)` until `flow.contributing_pixels()` at the
pixel reaches the channel threshold. That is what the command-line
programs' `basin_outlet_snapping_method: downstream` does. The
`python_binding_examples` notebook in lsdtt3_notebooks has a complete version,
including the conversion from latitude and longitude to a pixel.

Use `extract_channel_heads` when a simple accumulation threshold is not the
right source-detection method. The resulting heads can be passed to
`flow.scope_from_heads`, then into the same `channel_network` workflow.

`FlowGraph` is complementary: use it for distributed MFD or D-infinity
products such as TWI or sediment connectivity, rather than D8 channel trees.

See [Flow, channels, and basins API](../python_api/flow.md).
