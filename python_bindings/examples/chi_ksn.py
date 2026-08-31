"""Segment ksn map and chi-elevation diagnostic on the Green River DEM."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


dem = lsdtt3.fill_pits(lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif"))
flow = lsdtt3.FlowInfo(dem)
scope = flow.channel_scope(max(100, dem.valid_pixel_count // 250))
network = flow.channel_network(scope)

fit = network.chi(m_over_n=0.45).fit(
    dem, critical_divergence_metres=0.1
)
nodes = fit.nodes.copy()
segments = fit.segments
meta = dem.metadata

# Assign every fitted node the ksn of its containing piecewise-linear segment.
nodes["ksn"] = np.nan
for segment in segments.itertuples():
    in_segment = (
        (nodes.basin_id == segment.basin_id)
        & (nodes.branch_id == segment.branch_id)
        & (nodes.chi >= segment.downstream_chi - 1e-6)
        & (nodes.chi <= segment.upstream_chi + 1e-6)
    )
    nodes.loc[in_segment, "ksn"] = segment.ksn

nodes["row"] = nodes.flat // meta.n_cols
nodes["col"] = nodes.flat % meta.n_cols
hillshade = lsdtt3.multidirectional_hillshade(dem)

fig, (map_axis, profile_axis) = plt.subplots(
    1, 2, figsize=(12, 5), constrained_layout=True
)

map_axis.imshow(
    np.where(hillshade.mask_numpy(), hillshade.to_numpy(), np.nan),
    cmap="gray", origin="upper",
)
points = map_axis.scatter(
    nodes.col, nodes.row, c=nodes.ksn, s=5, cmap="magma", linewidths=0, vmin = np.nanpercentile(nodes.ksn, 20), vmax = np.nanpercentile(nodes.ksn, 80)
)
map_axis.set_axis_off()
map_axis.set_title("Fitted channel-network $k_{sn}$")
fig.colorbar(points, ax=map_axis, label="$k_{sn}$")

profile_axis.scatter(
    nodes.chi, nodes.observed_elevation, c=nodes.ksn,
    s=4, cmap="magma", linewidths=0, alpha=0.7, vmin = np.nanpercentile(nodes.ksn, 20), vmax = np.nanpercentile(nodes.ksn, 80)
)
for segment in segments.itertuples():
    profile_axis.plot(
        [segment.downstream_chi, segment.upstream_chi],
        [segment.downstream_elevation, segment.upstream_elevation],
        color="black", linewidth=1.0,
    )
profile_axis.set(
    xlabel="chi (m)", ylabel="elevation (m)",
    title="Chi-elevation profiles and fitted segments",
)

plt.show()
