"""Fetch and display a roughly 100 km DEM in the Romanian Carpathians."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


# Approximate centre near the Eastern Carpathians, Romania. These degree spans
# are about 100 km at this latitude; the fetched raster is projected to UTM.
centre_lon, centre_lat = 25.5, 47.1
half_width_deg = 0.66

dem = lsdtt3.fetch_dem(
    west=centre_lon - half_width_deg,
    south=centre_lat - 0.45,
    east=centre_lon + half_width_deg,
    north=centre_lat + 0.45,
    grid_spacing=30.0,
    cache_dir="dem_cache",
)
meta = dem.metadata
elevation = np.where(dem.mask_numpy(), dem.to_numpy(), np.nan)

hillshade_raster = lsdtt3.multidirectional_hillshade(dem)
hillshade = np.where(
    hillshade_raster.mask_numpy(), hillshade_raster.to_numpy(), np.nan
)

extent = [meta.x_min, meta.x_max, meta.y_min, meta.y_max]
fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)

image = ax.imshow(elevation, extent=extent, origin="upper", cmap="terrain")
ax.imshow(hillshade, extent=extent, origin="upper", cmap="gray", alpha=0.38)

ax.set_aspect("equal")
ax.set(title="Romanian Carpathians — UTM zone 35N",
       xlabel="easting (m)", ylabel="northing (m)")
fig.colorbar(image, ax=ax, label="elevation (m)")
plt.show()
