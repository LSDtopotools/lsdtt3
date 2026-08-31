"""Fetch a Carpathian DEM, make a synthetic 10Be outlet sample, and invert it."""

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


# Roughly 100 x 100 km in the Romanian Carpathians.  Ninety metres keeps this
# downloaded demonstration light; use 30 m for an actual study.
dem = lsdtt3.fetch_dem(
    west=25.0, south=46.70, east=26.35, north=47.60,
    grid_spacing=90.0, cache_dir="dem_cache",
)
dem = lsdtt3.fill_pits(dem)
flow = lsdtt3.FlowInfo(dem)
meta = dem.metadata

# The maximum contributing-area pixel is the largest D8 outlet represented by
# this clipped DEM.  It is a useful synthetic end-to-end target, not a claim
# that the real river mouth is inside this bounding box.
area = as_nan(flow.drainage_area())
outlet_row, outlet_col = np.unravel_index(np.nanargmax(area), area.shape)
outlet_x = meta.x_min + (outlet_col + 0.5) * meta.cell_size
outlet_y = meta.y_max - (outlet_row + 0.5) * meta.cell_size

# Forward-model a known steady erosion rate, then use its outlet concentration
# as a fake field measurement. The inversion should recover ~0.01 g cm-2 yr-1.
production = lsdtt3.crn_production_scaling(dem, atmosphere="standard")
known_erosion = lsdtt3.crn_uniform_raster(dem, 0.01)
local_be10 = lsdtt3.crn_steady_state(known_erosion, production, nuclide="Be10")
fake_concentration = lsdtt3.crn_outlet_concentration(
    local_be10, known_erosion, flow, int(outlet_row), int(outlet_col),
)
samples = gpd.GeoDataFrame(
    {"sample": ["synthetic-largest-outlet"], "concentration": [fake_concentration]},
    geometry=[Point(outlet_x, outlet_y)], crs=meta.crs_wkt,
)
result = lsdtt3.crn_invert_basin_erosion(
    flow, samples, production,
    channel_threshold_pixels=1_000,
    snap_radius_pixels=10,
)
print(result[["sample", "crn_status", "erosion_rate_g_cm2_yr", "erosion_rate_mm_kyr"]])

extent = (meta.x_min, meta.x_max, meta.y_min, meta.y_max)
fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
axes[0].imshow(as_nan(dem), cmap="terrain", extent=extent)
shade = as_nan(lsdtt3.multidirectional_hillshade(dem))
axes[0].imshow(shade, cmap="gray", alpha=0.35, extent=extent)
axes[0].scatter(outlet_x, outlet_y, marker="x", s=80, color="white", label="synthetic sample")
axes[0].scatter(result.snapped_x, result.snapped_y, s=30, color="crimson", label="snapped outlet")
axes[0].legend()
axes[0].set_title("Carpathian DEM and largest D8 outlet")

image = axes[1].imshow(area / 1e6, cmap="viridis", extent=extent)
axes[1].scatter(result.snapped_x, result.snapped_y, s=35, color="crimson")
axes[1].set_title("D8 drainage area (km²)")
fig.colorbar(image, ax=axes[1], label="km²")
plt.show()
