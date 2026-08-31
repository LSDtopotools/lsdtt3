import matplotlib.pyplot as plt

import lsdtt3


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
result = lsdtt3.spectral_analysis(dem, n_background_iterations=20)

plt.loglog(result.binned["wavelength"], result.binned["psd"])
plt.axvline(result.rollover_wavelength, color="black", linestyle="--")
plt.xlabel("wavelength (m)")
plt.ylabel("power spectral density")
plt.show()
