"""Compact channel-network analysis on the Green River DEM."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


dem = lsdtt3.fill_pits(lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif"))
flow = lsdtt3.FlowInfo(dem)
scope = flow.channel_scope(max(100, dem.valid_pixel_count // 250))
network = flow.channel_network(scope)

area = flow.drainage_area()
junctions = network.junctions()
outlets = junctions.loc[junctions.is_outlet, ["row", "col"]].to_numpy()
chi = network.chi(outlets=outlets).raster()
basins = flow.basins()
links = network.link_bearings()


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
axes[0].imshow(np.log10(as_nan(area)), cmap="viridis")
axes[0].contour(as_nan(network.channel_mask), levels=[0.5], colors="white", linewidths=0.4)
axes[0].set_title("Channel network")

axes[1].imshow(as_nan(chi), cmap="magma")
axes[1].set_title("Channel chi")

axes[2].imshow(as_nan(basins.labels), cmap="tab20")
axes[2].contour(as_nan(network.channel_mask), levels=[0.5], colors="black", linewidths=0.35)
axes[2].set_title("Outlet basins")

for axis in axes:
    axis.set_axis_off()

print(links[["upstream_junction", "downstream_junction", "bearing_deg"]].head())
print(network.junction_angles().head())
plt.show()
