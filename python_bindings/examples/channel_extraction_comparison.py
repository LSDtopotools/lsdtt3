"""Compare threshold and Wiener channel-head extraction on Green River."""

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
hillshade_raster = lsdtt3.multidirectional_hillshade(dem)
hillshade = np.where(
    hillshade_raster.mask_numpy(), hillshade_raster.to_numpy(), np.nan
)

methods = [
    ("Threshold: 100 pixels", dict(method="threshold", threshold_contributing_pixels=100)),
    ("Threshold: 500 pixels", dict(method="threshold", threshold_contributing_pixels=500)),
    # Wiener lengths and areas must span several 30 m cells (one cell is 900 m²).
    ("Wiener: 0.1 km²", dict(method="wiener", pruning_drainage_area_m2=1.0e5,
                             surface_fitting_radius_metres=60.0,
                             connected_components_threshold_pixels=10)),
    ("Wiener: 0.5 km²", dict(method="wiener", pruning_drainage_area_m2=5.0e5,
                             surface_fitting_radius_metres=60.0,
                             connected_components_threshold_pixels=10)),
]

fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
for ax, (title, options) in zip(axes.flat, methods):
    heads = lsdtt3.extract_channel_heads(dem, flow, **options)
    network = flow.channel_network(flow.scope_from_heads(heads))
    channels = network.points_geodataframe()
    source_table = heads.table

    ax.imshow(hillshade, cmap="gray", origin="upper")
    ax.scatter(channels.col, channels.row, s=1, c="dodgerblue", alpha=0.7)
    ax.scatter(source_table.col, source_table.row, marker="x", s=42,
               c="crimson", linewidths=1.4)
    ax.set(title=f"{title} ({heads.n_heads} heads)", xticks=[], yticks=[])

plt.show()
