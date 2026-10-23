"""Minimal channel scope and junction-network example on Green River."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


# Fetch a small 30 m Copernicus DEM of the Green River canyon beside
# Canyonlands National Park, Utah (projected to UTM zone 12N). The tile is
# cached in dem_cache/ so later runs do not download it again.
dem = lsdtt3.fetch_dem(
    west=-110.10, south=38.36, east=-109.96, north=38.47,
    grid_spacing=30.0, cache_dir="dem_cache",
)
dem = lsdtt3.fill_pits(dem)
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
