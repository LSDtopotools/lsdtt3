"""Inspect a closed depression and retain it as a deliberate drainage exit."""

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

# Detect nested closed depressions and immediately obtain GIS-ready products.
depressions = lsdtt3.find_internal_drainage(
    dem,
    method="hierarchy",
    min_pixels=1,
    min_depth_metres=0.5,
)
print(depressions.candidates[["sink_id", "depth_m", "area_m2", "volume_m3"]])

# Keep a selected subset as drainage exits during filling/routing.
selected = depressions.candidates.nlargest(5, "depth_m")
if selected.empty:
    raise RuntimeError(
        "No depression passed these example filters. Lower min_pixels/min_depth_metres "
        "or choose another DEM."
    )
pits = selected[["pit_row", "pit_col"]].to_numpy(dtype=np.int64)
conditioned = lsdtt3.fill_pits_with_sinks(dem, pits, min_slope=1e-4)
flow = lsdtt3.FlowInfo(conditioned, forced_outlets=pits, allow_pits_as_outlets=True)

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
hillshade = as_nan(lsdtt3.multidirectional_hillshade(dem))
axes[0].imshow(as_nan(dem), cmap="terrain")
axes[0].imshow(hillshade, cmap="gray", alpha=0.35)
if not depressions.polygons.empty:
    depressions.polygons.boundary.plot(ax=axes[0], color="deepskyblue", linewidth=0.8)
depressions.candidates.plot(ax=axes[0], color="crimson", marker="x", markersize=30)
axes[0].set(title="Detected depressions")

fill_change = as_nan(lsdtt3.subtract(conditioned, dem))
axes[1].imshow(fill_change, cmap="magma")
axes[1].set(title=f"Sink-preserving fill ({flow.n_outlets} outlets)")
for axis in axes:
    axis.set_axis_off()
plt.show()
