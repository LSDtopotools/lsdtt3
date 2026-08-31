"""Place a stochastic landslide event and inspect its Be10 response."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3

dem = lsdtt3.fill_pits(lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif"))
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
