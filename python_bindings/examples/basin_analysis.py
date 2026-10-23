"""Minimal basin selection and multiscale surface-metrics example."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


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

# Pick complete junction basins, retaining a compact robust subset.
selected = lsdtt3.junction_basins(network).filter_incomplete()
table = selected.table.nlargest(6, "area_metres2")
outlets = table[["outlet_row", "outlet_col"]].to_numpy()
basins = flow.basins(outlets)

# Conventional summaries and divide/interior polyfit metrics are DataFrames.
print(basins.statistics(dem, include_shape_metrics=True))
channel_lengths = lsdtt3.basin_channel_lengths(basins, network)
print(channel_lengths)
metrics = lsdtt3.basin_surface_metrics(
    basins, dem, scales=[30.0, 90.0], metrics=["slope", "maximum_curvature"]
)

hillshade = as_nan(lsdtt3.multidirectional_hillshade(dem))
labels = as_nan(basins.labels)
print(metrics["curvature_flips"])
print(metrics["slope_area_summary"])
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

# The tabular products are ready for ordinary plotting or modelling.
slope_area = metrics["slope_area"]
fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), constrained_layout=True)
for basin_id, group in slope_area.groupby("basin_id"):
    axes[0].plot(group.area_metres2, group.median_slope, marker=".", lw=1, label=f"basin {basin_id}")
axes[0].set(xscale="log", xlabel="drainage area (m²)", ylabel="median slope", title="Slope–area")
axes[0].legend(fontsize=7)

length_columns = [column for column in channel_lengths if column.startswith("channel_length_order_")]
axes[1].bar(range(1, len(length_columns) + 1), channel_lengths.loc[0, length_columns])
axes[1].set(xlabel="Strahler order", ylabel="channel length (m)", title="Channel lengths: basin 0")
plt.show()
