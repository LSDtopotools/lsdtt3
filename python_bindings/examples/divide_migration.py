"""Compare drainage-divide migration proxies on the Green River DEM."""

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
network = flow.channel_network(flow.channel_scope(300))

# Use independent, complete channel basins as the two/three comparison units.
choices = lsdtt3.junction_basins(network).filter_incomplete().filter_nested()
outlets = choices.table.nlargest(3, "area_metres2")[["outlet_row", "outlet_col"]]
result = lsdtt3.divide_migration(
    flow, dem, outlets.to_numpy(), source_pixels=300, outlet_control="none"
)

hillshade_raster = lsdtt3.multidirectional_hillshade(dem)
hillshade = np.where(
    hillshade_raster.mask_numpy(), hillshade_raster.to_numpy(), np.nan
)
meta = dem.metadata
extent = (
    meta.x_min,
    meta.x_min + meta.n_cols * meta.cell_size,
    meta.y_min,
    meta.y_min + meta.n_rows * meta.cell_size,
)
fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
axes[0].imshow(hillshade, cmap="gray", extent=extent, origin="upper")
result.divides.plot(ax=axes[0], color="orange", linewidth=1.5)
result.channel_heads.plot(ax=axes[0], column="chi", cmap="viridis", markersize=8)
axes[0].set_title("Divides and channel heads")
axes[0].set_axis_off()

axes[1].imshow(hillshade, cmap="gray", extent=extent, origin="upper")
heads = result.channel_heads
image = axes[1].scatter(
    heads.head_x, heads.head_y, c=heads.gradient_to_divide,
    cmap="coolwarm", s=12,
)
result.divides.plot(ax=axes[1], color="black", linewidth=1)
fig.colorbar(image, ax=axes[1], label="head-to-divide gradient")
axes[1].set_title("Migration gradient proxy")
axes[1].set_axis_off()
print(result.pair_summary)
plt.show()
