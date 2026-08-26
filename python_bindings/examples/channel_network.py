"""Minimal channel scope and junction-network example on Green River."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


dem = lsdtt3.fill_pits(lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif"))
flow = lsdtt3.FlowInfo(dem)
scope = flow.channel_scope(max(100, dem.valid_pixel_count // 250))
network = flow.channel_network(scope)
junctions = network.junctions()

valid = dem.mask_numpy()
area = np.where(valid, flow.drainage_area().to_numpy(), np.nan)

fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
axes[0].imshow(np.log10(area), cmap="viridis")
axes[0].contour(scope.mask.to_numpy(), levels=[0.5], colors="white", linewidths=0.5)
axes[0].set_title("Resolved channel scope")

axes[1].imshow(np.log10(area), cmap="viridis")
axes[1].contour(network.channel_mask.to_numpy(), levels=[0.5], colors="white", linewidths=0.5)
axes[1].scatter(junctions.col, junctions.row, s=8, c=junctions.stream_order, cmap="magma")
axes[1].set_title("Junction network")

for axis in axes:
    axis.set_axis_off()

print(network.links().head())
plt.show()
