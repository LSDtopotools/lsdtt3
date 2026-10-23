"""Run repeated landslides and map their final cosmogenic response."""

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
series = lsdtt3.run_landslide_series(
    dem, flow,
    end_time_yr=50_000.0,
    mean_recurrence_yr=2_500.0,
    long_term_erosion_rate_mm_kyr=0.2,
    channel_threshold_pixels=max(200, dem.valid_pixel_count // 400),
    anneal_time_yr=10_000.0,
    seed=42,
)
be10 = lsdtt3.landslide_crn_concentration(
    series.excavation_depth, series.failure_time,
    report_time_yr=50_000.0,
)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
for ax, raster, title, cmap in zip(
    axes,
    (series.elevation, series.excavation_depth, be10),
    ("Final elevation", "Latest excavation depth", "Final Be10"),
    ("terrain", "magma", "viridis"),
):
    image = ax.imshow(np.where(raster.mask_numpy(), raster.to_numpy(), np.nan), cmap=cmap)
    ax.set_title(title)
    ax.set_axis_off()
    fig.colorbar(image, ax=ax)

print(series.summary)
print("event times (yr):", series.event_times_yr)
plt.show()
