"""Cellular hillslopes, depression topology, and CRN-history products."""

import numpy as np

import lsdtt3


# Fetch a small 30 m Copernicus DEM of the Green River canyon beside
# Canyonlands National Park, Utah (projected to UTM zone 12N). The tile is
# cached in dem_cache/ so later runs do not download it again.
dem = lsdtt3.fetch_dem(
    west=-110.10, south=38.36, east=-109.96, north=38.47,
    grid_spacing=30.0, cache_dir="dem_cache",
)
channels = lsdtt3.crn_uniform_raster(dem, 0.0)

# Conservative cellular transport can also run independently of fluvial incision.
cellular = lsdtt3.hillslope_transport_cellular(
    dem, channels, duration_years=1_000.0, diffusivity=1.0e-2
)

# The stateful LEM accepts the same cellular process model.
model = lsdtt3.LandscapeModel(
    dem, timestep_years=1_000.0, hillslope_model="cellular",
    fluvial_method="dinf", coupling_regime="coupled", uplift_rate=0.0,
)
final_surface = model.run(10).elevation
print(model.coupling_regime, cellular.n_face_events)

# Inspect all nested depressions rather than only candidate sinks.
hierarchy = lsdtt3.depression_hierarchy(final_surface)
print(hierarchy.nodes[["depression_id", "parent_id", "volume_m3"]].head())

# One present-to-past erosion history drives Be-10, Al-26, and C-14 together.
column = lsdtt3.CRNColumn()
column.set_erosion_history(
    [5_000.0, 20_000.0, -9_999.0],
    [8.0e-4, 2.0e-4, 2.0e-4],
    [0.0, 0.0, 0.0],
)
production = lsdtt3.crn_uniform_raster(final_surface, 1.0)
concentrations = lsdtt3.crn_transient_history(column, production)
be10 = concentrations["be10"]
