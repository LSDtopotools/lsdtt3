"""lsdtt-slope-area's routines in Python (plan_slope_area_gradient_and_low_concavity.md,
step 7): conditioned_pixel_mask, slope_area_bins, FlowInfo.d8_flow_path_gradient,
FlowInfo.low_concavity and FlowInfo.low_concavity_segments. The cases mirror the
C++ driver tests/cpp/test_slope_area.cpp; the column names are the CLI's."""

import numpy as np
import pytest
from pyproj import CRS

import lsdtt3

pytestmark = [pytest.mark.unit]

WKT = CRS.from_epsg(32633).to_wkt()


def _raster(values, cell_size=10.0):
    values = np.asarray(values, dtype=np.float32)
    meta = lsdtt3.RasterMeta.from_shape(values.shape[0], values.shape[1], x_min=0.0, y_min=0.0,
                                        cell_size=cell_size, crs_wkt=WKT, epsg_code=32633)
    return lsdtt3.Raster.from_numpy(values, np.ones(values.shape, dtype=bool), meta)


def _east_plane(nr=9, nc=30, cell_size=10.0, gradient=0.05):
    z = 300.0 - gradient * cell_size * np.arange(nc)[None, :] * np.ones((nr, 1))
    return _raster(z, cell_size)


def _valid(raster):
    return raster.mask_numpy().astype(bool)


@pytest.mark.exact_solution
def test_d8_flow_path_gradient_on_a_plane():
    dem = _east_plane()
    flow = lsdtt3.FlowInfo(dem)
    g = flow.d8_flow_path_gradient(dem, flow.drainage_area(), window_metres=60.0)
    row = g.to_numpy()[4, 1:-1]
    assert _valid(g)[4, 1:-1].all()
    np.testing.assert_allclose(row, 0.05, atol=1e-6)
    with pytest.raises(ValueError, match="four cell sizes"):
        flow.d8_flow_path_gradient(dem, flow.drainage_area(), window_metres=39.0)


@pytest.mark.exact_solution
def test_one_sided_windows_on_30_m_cells():
    """Plan section 10 F1: the divide pixel's window reaches 60 m of the 75 m
    half-window; the span rule (half less one diagonal step) accepts it."""
    dem = _east_plane(nc=20, cell_size=30.0)
    flow = lsdtt3.FlowInfo(dem)
    g = flow.d8_flow_path_gradient(dem, flow.drainage_area(), window_metres=150.0)
    assert _valid(g)[4, 0] and _valid(g)[4, -1]
    assert g.to_numpy()[4, 0] == pytest.approx(0.05, abs=1e-6)


@pytest.mark.exact_solution
def test_conditioned_pixel_mask_codes_and_exclusion():
    z = 100.0 - 0.5 * np.arange(15)[None, :] * np.ones((15, 1))
    z[6:9, 6:9] -= 5.0                                    # a depression
    dem = _raster(z)
    filled = lsdtt3.fill_pits(dem)
    codes = lsdtt3.conditioned_pixel_mask(dem, filled).to_numpy()
    raised = filled.to_numpy() > dem.to_numpy()
    assert np.array_equal(codes == 1, raised) and raised[6:9, 6:9].all()
    assert not (codes == 2).any()
    carved = lsdtt3.breach_least_cost(dem)
    ccodes = lsdtt3.conditioned_pixel_mask(dem, carved).to_numpy()
    assert np.array_equal(ccodes == 2, carved.to_numpy() < dem.to_numpy()) and (ccodes == 2).any()
    assert not lsdtt3.conditioned_pixel_mask(dem, dem).to_numpy().any()


@pytest.mark.property
@pytest.mark.negative
def test_slope_area_bins_columns_counts_and_exclusion():
    rng = np.random.default_rng(3)
    area = _raster(np.exp(rng.uniform(np.log(900.0), np.log(1e7), size=(40, 40))))
    slope = _raster(rng.uniform(0.01, 0.6, size=(40, 40)))
    bins = lsdtt3.slope_area_bins(area, slope, n_log_bins=10)
    assert list(bins.columns) == ["drainage_area", "d05_slope", "d25_slope", "d50_slope", "d75_slope",
                                  "d95_slope", "number_of_pixels"]
    assert bins.number_of_pixels.sum() == 1600
    assert (np.diff(bins.drainage_area) > 0).all()
    q = bins[["d05_slope", "d25_slope", "d50_slope", "d75_slope", "d95_slope"]].to_numpy()
    assert (np.diff(q, axis=1) >= 0).all()
    ex = np.zeros((40, 40), dtype=np.float32)
    ex[:5, :] = 2.0                                       # 200 "carved" pixels
    fewer = lsdtt3.slope_area_bins(area, slope, n_log_bins=10, exclude=_raster(ex))
    assert fewer.number_of_pixels.sum() == 1400
    with pytest.raises(ValueError, match="n_log_bins"):
        lsdtt3.slope_area_bins(area, slope, n_log_bins=0)


def _power_law_rasters(theta, nr=5, nc=30):
    a = 1000.0 * 1.3 ** np.arange(nc)[None, :] * np.ones((nr, 1))
    return _raster(a), _raster(0.5 * a ** (-theta))


@pytest.mark.exact_solution
@pytest.mark.parametrize("theta, flag", [(0.5, 0), (0.1, 1), (-0.1, 1), (-0.3, 0)])
def test_low_concavity_flag_rule(theta, flag):
    dem = _east_plane(nr=5, nc=30, gradient=0.1)
    flow = lsdtt3.FlowInfo(dem)
    area, gradient = _power_law_rasters(theta)
    t = flow.low_concavity(area, gradient, threshold_contributing_pixels=1)
    assert list(t.columns) == ["flat", "x", "y", "drainage_area_m2", "slope", "theta", "theta_se", "r2",
                               "n_pixels", "log10_area_range", "n_junctions_in_window", "low_concavity"]
    assert len(t) == 5 * 30
    row = t[t.flat == 2 * 30 + 5].iloc[0]
    assert row.theta == pytest.approx(theta, abs=1e-5)
    assert row.theta_se == pytest.approx(0.0, abs=1e-5)
    assert row.low_concavity == flag
    near_outlet = t[t.flat == 2 * 30 + 27].iloc[0]          # 3 pixels left before the outlet
    assert near_outlet.low_concavity == -1


@pytest.mark.exact_solution
@pytest.mark.negative
def test_low_concavity_exclusion_and_refusals():
    dem = _east_plane(nr=5, nc=30, gradient=0.1)
    flow = lsdtt3.FlowInfo(dem)
    area, gradient = _power_law_rasters(0.1)
    ex = np.zeros((5, 30), dtype=np.float32)
    ex[2, 5] = 2.0
    ex[2, 8] = 1.0
    t = flow.low_concavity(area, gradient, threshold_contributing_pixels=1, exclude=_raster(ex))
    assert t[t.flat == 2 * 30 + 5].iloc[0].low_concavity == -1
    r = t[t.flat == 2 * 30 + 4].iloc[0]
    assert r.n_pixels == 8 and r.theta == pytest.approx(0.1, abs=1e-5)
    with pytest.raises(ValueError):
        flow.low_concavity(area, gradient, threshold_contributing_pixels=1, window_pixels=1)
    with pytest.raises(ValueError):
        flow.low_concavity(area, gradient, threshold_contributing_pixels=0)


@pytest.mark.property
def test_low_concavity_segments_partition_the_channel_pixels():
    dem = _east_plane(nr=5, nc=30, gradient=0.1)
    flow = lsdtt3.FlowInfo(dem)
    area, gradient = _power_law_rasters(0.1)
    t = flow.low_concavity(area, gradient, threshold_contributing_pixels=1)
    seg = flow.low_concavity_segments(area, gradient, threshold_contributing_pixels=1)
    assert seg.n_pixels.sum() == len(t)
    for flag in (-1, 0, 1):
        assert seg.loc[seg.low_concavity == flag, "n_pixels"].sum() == (t.low_concavity == flag).sum()
    assert seg.crs is not None and (seg.geom_type == "LineString").all()
    assert {"median_theta_se", "stream_order", "link_id"} <= set(seg.columns)
    # five rows, each one link: the flag runs 1 (cols 0-25) then -1 (the last 4)
    assert seg.link_id.nunique() == 5
