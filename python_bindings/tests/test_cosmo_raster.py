import numpy as np

import lsdtt3


def _raster(values):
    values = np.asarray(values, dtype=np.float32)
    meta = lsdtt3.RasterMeta.from_shape(*values.shape, cell_size=30.0)
    return lsdtt3.Raster.from_numpy(values, np.ones(values.shape, dtype=bool), meta)


def test_steady_state_ratio_and_accumulation():
    dem = _raster([[1000, 990, 980], [1010, 1000, 970], [1020, 1010, 960]])
    erosion = lsdtt3.crn_uniform_raster(dem, 0.01)
    production = lsdtt3.crn_uniform_raster(dem, 1.0)
    be10 = lsdtt3.crn_steady_state(erosion, production, "Be10")
    al26 = lsdtt3.crn_steady_state(erosion, production, "Al26")
    ratio = lsdtt3.crn_isotope_ratio(al26, be10)
    assert ratio["n_valid"] == 9
    assert np.all(ratio["raster"].to_numpy() > 0)
    flow = lsdtt3.FlowInfo(dem)
    accumulated = lsdtt3.crn_accumulated_concentration(be10, erosion, flow)
    assert accumulated.shape == dem.shape


def test_arbitrary_history_returns_all_three_nuclides():
    dem = _raster([[1000, 990], [1010, 980]])
    column = lsdtt3.CRNColumn()
    column.set_erosion_history(
        [1_000.0, -9_999.0], [1.0e-3, 1.0e-3], removal_metres=[0.0, 0.0]
    )
    result = lsdtt3.crn_transient_history(column, lsdtt3.crn_uniform_raster(dem, 1.0))
    assert set(result) == {"be10", "al26", "c14"}
    assert all(raster.shape == dem.shape for raster in result.values())
