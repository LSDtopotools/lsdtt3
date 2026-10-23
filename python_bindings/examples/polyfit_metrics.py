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
metrics = lsdtt3.polyfit_metrics(
    dem,
    ["slope_degrees", "mean_curvature", "plan_curvature"],
    lambda_metres=90.0,
)

fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
for axis, (name, raster) in zip(axes, metrics.items()):
    image = axis.imshow(as_nan(raster), cmap="terrain")
    axis.set_title(name.replace("_", " "))
    axis.set_axis_off()
    fig.colorbar(image, ax=axis, shrink=0.8)

plt.show()
