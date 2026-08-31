import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


size = 64
meta = lsdtt3.RasterMeta.from_shape(size, size, cell_size=1.0)
row, col = np.indices((size, size), dtype=np.float32)
rng = np.random.default_rng(8)
terrain = 0.08 * row + 0.04 * (col - size / 2) ** 2 + 0.03 * rng.standard_normal((size, size))
dem = lsdtt3.Raster.from_numpy(terrain, np.ones((size, size), dtype=bool), meta)

flats = lsdtt3.find_flat_surfaces(
    dem, madk2_window_radius_metres=3.0, use_residual_iqr=False,
    min_component_pixels=1, relief_window_pixels=5, relief_threshold_metres=1.0,
)
roads = lsdtt3.find_roads(dem, smoothing_window_metres=3.0)
road_flow = lsdtt3.find_road_flow_discrimination(
    dem, regional_smoothing_metres=5.0, flow_accumulation_threshold_pixels=5,
    intermediate_window_metres=3.0, roughness_radius_metres=1.0,
)
# Marsh extraction is demonstrated on a separate, low-relief synthetic tile.
marsh_dem = lsdtt3.Raster.from_numpy(
    0.01 * row + 0.005 * col + 0.03 * rng.standard_normal((size, size)),
    np.ones((size, size), dtype=bool), meta,
)
marsh = lsdtt3.find_marsh_platform(
    marsh_dem, max_scarp_growth_iterations=2, max_platform_growth_iterations=2,
)

channel = np.column_stack((np.arange(4, size - 4), np.full(size - 8, size // 2)))
valleys = lsdtt3.extract_valley_centrelines(
    dem, [channel], corridor_width_metres=15.0, fallback_valley_hand_metres=3.0,
    centreline_loops=1,
)
ridge = lsdtt3.Ridgeline(dem, np.column_stack((np.full(30, 8), np.arange(16, 46))))
ridge.fill_elevation(dem)
ridge.fit_profile()

def image(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
axes[0, 0].imshow(image(dem), cmap="terrain")
axes[0, 0].plot(ridge.pixels[:, 1], ridge.pixels[:, 0], color="white", lw=1)
axes[0, 0].set_title("DEM and Ridgeline")
axes[0, 1].imshow(image(flats.labelled_flat_mask), cmap="tab20")
axes[0, 1].set_title("Flat surfaces")
axes[0, 2].imshow(image(roads.mask), cmap="gray_r", vmin=0, vmax=1)
axes[0, 2].set_title("Road/bench mask")
axes[1, 0].imshow(image(road_flow.flow_path_divergence), cmap="magma")
axes[1, 0].set_title("Road-flow divergence")
axes[1, 1].imshow(image(valleys.valley_mask), cmap="Blues")
centreline = valleys.centreline_pixels(0)
axes[1, 1].plot(centreline[:, 1], centreline[:, 0], color="orange", lw=1)
axes[1, 1].set_title("Valley and centreline")
axes[1, 2].imshow(image(marsh.platform_mask), cmap="viridis", vmin=0, vmax=1)
axes[1, 2].set_title("Marsh-platform mask")
for axis in axes.flat:
    axis.set_axis_off()
plt.show()
