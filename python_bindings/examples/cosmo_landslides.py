"""Place a stochastic landslide event and inspect its Be10 response."""

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
event = lsdtt3.trigger_landslides(
    dem, flow, target_volume_m3=2.0e5,
    channel_threshold_pixels=max(200, dem.valid_pixel_count // 400), seed=42,
)
be10 = lsdtt3.landslide_crn_concentration(
    event.excavation_depth, event.failure_time,
    background_erosion_rate_g_cm2_yr=0.001,
)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
for ax, raster, title, cmap in zip(
    axes,
    (event.elevation, event.excavation_depth, be10),
    ("Post-event elevation", "Excavation depth", "Post-landslide Be10"),
    ("terrain", "magma", "viridis"),
):
    image = ax.imshow(np.where(raster.mask_numpy(), raster.to_numpy(), np.nan), cmap=cmap)
    ax.set_title(title)
    ax.set_axis_off()
    fig.colorbar(image, ax=ax)

print(event.summary)
plt.show()
