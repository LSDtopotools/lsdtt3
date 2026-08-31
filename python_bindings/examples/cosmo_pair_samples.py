"""Apply two-nuclide burial diagnostics directly to a sample GeoDataFrame."""

import geopandas as gpd
from shapely.geometry import Point

import lsdtt3


samples = gpd.GeoDataFrame(
    {
        "sample": ["GR-01", "GR-02"],
        "Al26_CONC": [1.5e6, 1.2e6],
        "Be10_CONC": [2.5e5, 2.2e5],
    },
    geometry=[Point(500_000, 4_400_000), Point(501_000, 4_401_000)],
    crs="EPSG:32633",
)

burial = lsdtt3.crn_pair_diagnostics(
    samples,
    mode="burial",
    nuclide_fast="Al26", nuclide_slow="Be10",
)
print(burial[[
    "sample", "crn_pair_status", "crn_pair_burial_required",
    "crn_pair_erosion_rate_g_cm2_yr", "crn_pair_time_yr",
]])

# A known-pre-change-rate step-change inversion uses the same input table.
step = lsdtt3.crn_pair_diagnostics(
    samples,
    mode="step_change",
    nuclide_fast="Al26", nuclide_slow="Be10",
    effective_erosion_rate_old_g_cm2_yr=0.01,
)
print(step[["sample", "crn_pair_status", "crn_pair_erosion_rate_g_cm2_yr", "crn_pair_time_yr"]])
