# Chi and ksn

Chi analysis builds on the D8 channel network. The normal sequence is DEM,
`FlowInfo`, a channel scope/network, chi coordinates, and a piecewise-linear
fit.

```python
import lsdtt3

dem = lsdtt3.fill_pits(lsdtt3.Raster.read("dem.tif"))
flow = lsdtt3.FlowInfo(dem)
network = flow.channel_network(flow.channel_scope(threshold_pixels=500))

chi = network.chi(m_over_n=0.45)
fit = chi.fit(dem, critical_divergence_metres=0.1)

nodes = fit.nodes       # GeoDataFrame / table of fitted channel nodes
segments = fit.segments # one row per fitted chi segment, including ksn
```

Use the segment table for mapped `k_sn`; join or assign its values to the node
table for point-based plotting. See the repository example
`python_bindings/examples/chi_ksn.py` for a map plus chi-elevation panel.

See [Chi and ksn API](../python_api/chi.md).
