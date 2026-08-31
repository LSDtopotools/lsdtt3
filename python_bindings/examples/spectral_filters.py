import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
rasters = {
    "DEM": dem,
    "Low-pass": lsdtt3.lowpass_filter(dem, cutoff_metres=500.0),
    "High-pass": lsdtt3.highpass_filter(dem, cutoff_metres=500.0),
    "Band-pass": lsdtt3.bandpass_filter(dem, short_metres=30.0, long_metres=500.0),
    "Wiener": lsdtt3.wiener_filter(dem, signal_wavelength_long_m=500.0,
                                    signal_wavelength_short_m=30.0),
}

fig, axes = plt.subplots(1, len(rasters), figsize=(16, 4), constrained_layout=True)
for axis, (name, raster) in zip(axes, rasters.items()):
    image = axis.imshow(as_nan(raster), cmap="terrain")
    axis.set_title(name)
    axis.set_axis_off()
    fig.colorbar(image, ax=axis, shrink=0.75)

plt.show()
