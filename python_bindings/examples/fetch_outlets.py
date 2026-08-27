"""Fetch grouped in-memory DEMs around outlet points."""

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point

import lsdtt3

# WGS84 outlets in the Romanian Carpathians. Nearby points share one DEM.
outlets = gpd.GeoDataFrame(
    {"name": ["upper", "lower"]},
    geometry=[Point(25.50, 47.10), Point(25.64, 47.02)],
    crs="EPSG:4326",
)

# These demonstration coordinates are approximate rather than digitised river
# cells, so search locally for a sufficiently contributing channel pixel.
result = lsdtt3.fetch_dem_for_outlets(
    outlets,
    initial_radius_km=15.0,
    snapping_method="expanding_square",
    max_search_radius_pixels=100,
)
print(result.groups[["group_id", "outlet_ids", "area_km2"]])

for dem, basin_labels, group in zip(
    result.rasters, result.basin_labels, result.groups.itertuples()
):
    dem.write(f"carpathians_group_{group.group_id}.tif")

    meta = dem.metadata
    elevation = dem.to_numpy()
    elevation[~dem.mask_numpy()] = np.nan
    labels = basin_labels.to_numpy()
    labels[~basin_labels.mask_numpy()] = np.nan
    projected_outlets = outlets.to_crs(meta.crs_wkt)
    snapped = result.snapped_outlets[
        result.snapped_outlets.group_id == group.group_id
    ]

    fig, ax = plt.subplots(figsize=(8, 6))
    image = ax.imshow(
        elevation,
        cmap="terrain",
        origin="upper",
        extent=(meta.x_min, meta.x_max, meta.y_min, meta.y_max),
    )
    ax.imshow(
        labels,
        cmap="tab20",
        alpha=0.35,
        origin="upper",
        interpolation="none",
        extent=(meta.x_min, meta.x_max, meta.y_min, meta.y_max),
    )
    ax.scatter(projected_outlets.geometry.x, projected_outlets.geometry.y,
               marker="x", color="black", s=60, label="input outlets")
    ax.scatter(snapped.x, snapped.y, marker="o", facecolors="none",
               edgecolors="crimson", s=70, label="snapped outlets")
    ax.set_title(f"Fetched outlet DEM — group {group.group_id}")
    ax.set_aspect("equal")
    ax.legend()
    fig.colorbar(image, ax=ax, label="elevation (m)")

plt.show()
