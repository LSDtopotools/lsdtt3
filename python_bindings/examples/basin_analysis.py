"""Minimal basin selection and multiscale surface-metrics example."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


dem = lsdtt3.fill_pits(lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif"))
flow = lsdtt3.FlowInfo(dem)
scope = flow.channel_scope(max(100, dem.valid_pixel_count // 250))
network = flow.channel_network(scope)

# Pick complete second-order junction basins, then retain a compact subset.
selected = lsdtt3.junction_basins(network).filter_incomplete().filter_stream_order(2, 2)
table = selected.table.head(6)
outlets = table[["outlet_row", "outlet_col"]].to_numpy()
basins = flow.basins(outlets)

# Conventional summaries and divide/interior polyfit metrics are DataFrames.
print(basins.statistics(dem))
metrics = lsdtt3.basin_surface_metrics(
    basins, dem, scales=[30.0, 90.0], metrics=["slope", "maximum_curvature"]
)

hillshade = as_nan(lsdtt3.multidirectional_hillshade(dem))
labels = as_nan(basins.labels)
print(metrics["curvature_flips"])
curvature = as_nan(lsdtt3.polyfit_metrics(dem, ["maximum_curvature"], lambda_metres=90.0)["maximum_curvature"])
curvature[~np.isfinite(labels)] = np.nan

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

axes[0].imshow(hillshade, cmap="gray")
axes[0].contour(network.channel_mask.to_numpy(), levels=[0.5], colors="deepskyblue", linewidths=0.4)
axes[0].imshow(labels, cmap="tab20", alpha=0.35, interpolation="nearest")
axes[0].scatter(table.outlet_col, table.outlet_row, marker="x", color="crimson", s=35)
axes[0].set_title("Selected junction basins")

axes[1].imshow(hillshade, cmap="gray")
image = axes[1].imshow(curvature, cmap="coolwarm", vmin=-0.02, vmax=0.02, alpha=0.8)
axes[1].contour(labels, levels=np.arange(0.5, basins.size + 0.5), colors="black", linewidths=0.4)
axes[1].set_title("Maximum curvature, 90 m scale")
fig.colorbar(image, ax=axes[1], shrink=0.8, label="curvature")

for axis in axes:
    axis.set_axis_off()
plt.show()
