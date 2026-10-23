import matplotlib.pyplot as plt

import lsdtt3


# Fetch a small 30 m Copernicus DEM of the Green River canyon beside
# Canyonlands National Park, Utah (projected to UTM zone 12N). The tile is
# cached in dem_cache/ so later runs do not download it again.
dem = lsdtt3.fetch_dem(
    west=-110.10, south=38.36, east=-109.96, north=38.47,
    grid_spacing=30.0, cache_dir="dem_cache",
)
result = lsdtt3.spectral_analysis(dem, n_background_iterations=20)

plt.loglog(result.binned["wavelength"], result.binned["psd"])
plt.axvline(result.rollover_wavelength, color="black", linestyle="--")
plt.xlabel("wavelength (m)")
plt.ylabel("power spectral density")
plt.show()
