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

Use `extract_channel_heads` when a simple accumulation threshold is not the
right source-detection method. The resulting heads can be passed to
`flow.scope_from_heads`, then into the same `channel_network` workflow.

`FlowGraph` is complementary: use it for distributed MFD or D-infinity
products such as TWI or sediment connectivity, rather than D8 channel trees.

See [Flow, channels, and basins API](../python_api/flow.md).
