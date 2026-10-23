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
shape = lsdtt3.classify_surface_shape(dem, lambda_metres=90.0, kt=1e-4)

plt.imshow(as_nan(shape["shape_class"]), cmap="RdBu", vmin=-4, vmax=3)
plt.colorbar(label="surface-shape class")
plt.axis("off")
plt.show()
