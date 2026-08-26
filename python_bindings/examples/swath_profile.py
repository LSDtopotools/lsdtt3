import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString

import lsdtt3


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
meta = dem.metadata
baseline = gpd.GeoDataFrame(
    geometry=[LineString([
        (meta.x_min + 0.10 * (meta.x_max - meta.x_min), meta.y_min + 0.35 * (meta.y_max - meta.y_min)),
        (meta.x_min + 0.28 * (meta.x_max - meta.x_min), meta.y_min + 0.62 * (meta.y_max - meta.y_min)),
        (meta.x_min + 0.48 * (meta.x_max - meta.x_min), meta.y_min + 0.42 * (meta.y_max - meta.y_min)),
        (meta.x_min + 0.70 * (meta.x_max - meta.x_min), meta.y_min + 0.68 * (meta.y_max - meta.y_min)),
        (meta.x_min + 0.90 * (meta.x_max - meta.x_min), meta.y_min + 0.48 * (meta.y_max - meta.y_min)),
    ])],
    crs=meta.crs_wkt,
)
profile = lsdtt3.swath_profile(dem, baseline, half_width_metres=500.0, bin_width_metres=100.0)

table = profile.table
fig, (map_axis, profile_axis) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

image = map_axis.imshow(
    np.where(dem.mask_numpy(), dem.to_numpy(), np.nan),
    extent=[meta.x_min, meta.x_max, meta.y_min, meta.y_max],
    origin="upper",
    cmap="terrain",
)
map_axis.imshow(
    np.where(profile.perpendicular_distance.mask_numpy(), 1.0, np.nan),
    extent=[meta.x_min, meta.x_max, meta.y_min, meta.y_max],
    origin="upper",
    cmap="Blues",
    alpha=0.18,
)
baseline.plot(ax=map_axis, color="black", linewidth=2)
map_axis.set_aspect("equal")
map_axis.set_title("DEM and swath baseline")
fig.colorbar(image, ax=map_axis, label="elevation")

profile_axis.fill_between(
    table["along_axis_centre"], table["p5"], table["p95"], alpha=0.25
)
profile_axis.plot(table["along_axis_centre"], table["mean"], color="black")
profile_axis.set_xlabel("distance along baseline (m)")
profile_axis.set_ylabel("elevation")
plt.show()
