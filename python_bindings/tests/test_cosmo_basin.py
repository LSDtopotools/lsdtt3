import geopandas as gpd
import numpy as np
from shapely.geometry import Point

import lsdtt3


def test_basin_inversion_returns_annotated_geodataframe():
    values = np.array([[100.0, 99.0, 98.0], [101.0, 100.0, 97.0], [102.0, 101.0, 96.0]], dtype=np.float32)
    meta = lsdtt3.RasterMeta.from_shape(3, 3, cell_size=30.0, x_min=500_000.0, y_min=4_400_000.0)
    dem = lsdtt3.Raster.from_numpy(values, np.ones(values.shape, dtype=bool), meta)
    flow = lsdtt3.FlowInfo(dem)
    production = lsdtt3.crn_uniform_raster(dem, 1.0)
    samples = gpd.GeoDataFrame(
        {"concentration": [200_000.0]},
        geometry=[Point(500_045.0, 4_400_045.0)],
    )

    result = lsdtt3.crn_invert_basin_erosion(
        flow, samples, production, channel_threshold_pixels=1, snap_radius_pixels=1,
    )

    assert result.loc[0, "crn_status"] == "ok"
    assert result.loc[0, "basin_pixels"] >= 1
    assert result.loc[0, "erosion_rate_g_cm2_yr"] > 0
    assert result.loc[0, "predicted_concentration_atoms_g"] > 0


def test_pair_diagnostics_returns_annotated_geodataframe():
    column = lsdtt3.CRNColumn()
    steady = column.steady_state(
        effective_erosion_rate_g_cm2_yr=0.01,
        production_scaling=1.0,
    )
    samples = gpd.GeoDataFrame(
        {"Al26_CONC": [steady.concentration_al26], "Be10_CONC": [steady.concentration_be10]},
        geometry=[Point(0.0, 0.0)],
    )
    result = lsdtt3.crn_pair_diagnostics(samples, mode="burial")
    assert result.loc[0, "crn_pair_status"] == "Converged"
    assert result.loc[0, "crn_pair_erosion_rate_g_cm2_yr"] > 0
    assert result.loc[0, "crn_pair_time_yr"] >= 0
