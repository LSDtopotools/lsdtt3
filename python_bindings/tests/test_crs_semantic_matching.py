"""Rasters in one CRS whose WKT text differs are accepted together.

An lsdtt3-written DEM and a mask made by gdalwarp in the same oblique Mercator
frame carry the same CRS in different WKT dialects. Every binding grid check
compares CRSs by meaning; a frame with a different azimuth is still refused.
"""

import numpy as np
import pytest
from pyproj import CRS

import lsdtt3


def _omerc(alpha):
    return CRS.from_proj4(
        f"+proj=omerc +lat_0=38.4 +lonc=-110.0 +alpha={alpha} +gamma=90 +k=1 "
        "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
    )


WKT1 = _omerc(30.0).to_wkt("WKT1_GDAL")
WKT2 = _omerc(30.0).to_wkt("WKT2_2019")
ROTATED = _omerc(31.0).to_wkt("WKT2_2019")
N = 21


def _raster(values, wkt):
    meta = lsdtt3.RasterMeta.from_shape(
        N, N, x_min=-105.0, y_min=-105.0, cell_size=10.0, crs_wkt=wkt, epsg_code=0
    )
    return lsdtt3.Raster.from_numpy(
        np.asarray(values, dtype=np.float32), np.ones((N, N), dtype=bool), meta
    )


@pytest.fixture(scope="module")
def frame():
    rows, cols = np.indices((N, N))
    elevation = 200.0 - 4.0 * rows + 2.0 * np.abs(cols - 10)
    dem = _raster(elevation, WKT1)
    flow = lsdtt3.FlowInfo(
        dem, boundary_conditions="nnnn", forced_outlets=np.array([[N - 1, 10]]),
        allow_pits_as_outlets=True,
    )
    return elevation, dem, flow


def test_wkt_dialects_differ_as_text():
    assert WKT1 != WKT2


CASES = {
    "flow_accumulate": lambda e, dem, flow, wkt: flow.accumulate(
        _raster(np.ones_like(e), wkt)),
    "flow_graph_accumulate": lambda e, dem, flow, wkt: lsdtt3.FlowGraph(
        dem, method="dinf", boundary_conditions="nnnn").accumulate(
        _raster(np.ones_like(e), wkt)),
    "topographic_wetness_index": lambda e, dem, flow, wkt: lsdtt3.topographic_wetness_index(
        flow.drainage_area(), _raster(np.full_like(e, 0.1), wkt)),
    "extract_channel_heads": lambda e, dem, flow, wkt: lsdtt3.extract_channel_heads(
        _raster(e, wkt), flow, threshold_contributing_pixels=5),
    "crn_transient_history": lambda e, dem, flow, wkt: _crn(dem, wkt),
    "fill_hydraulic_surface": lambda e, dem, flow, wkt: _hydraulic(e, dem, wkt),
    "basin_outline_ridgelines": lambda e, dem, flow, wkt: flow.basins(
        np.array([[N - 1, 10]])).outline_ridgelines(_raster(e, wkt)),
    "basin_ridgeline_extrema": lambda e, dem, flow, wkt: lsdtt3.basin_ridgeline_extrema(
        flow.basins(np.array([[N - 1, 10]])),
        flow.channel_network(flow.channel_scope(5)), _raster(e, wkt),
        detrend_smoothing_pixels=1),
}


def _crn(dem, wkt):
    column = lsdtt3.CRNColumn()
    column.set_erosion_history([1_000.0, -9_999.0], [1.0e-3, 1.0e-3], [0.0, 0.0])
    return lsdtt3.crn_transient_history(
        column, lsdtt3.crn_uniform_raster(dem, 1.0),
        topographic_shielding=_raster(np.ones((N, N)), wkt),
    )


def _hydraulic(elevation, dem, wkt):
    states = np.ones((N, N), dtype=np.float32)
    states[N - 1, 10] = 3
    return lsdtt3.fill_hydraulic_surface(dem, _raster(states, wkt))


@pytest.mark.parametrize("name", sorted(CASES))
def test_same_frame_in_another_wkt_dialect_is_accepted(frame, name):
    CASES[name](*frame, WKT2)


@pytest.mark.parametrize("name", sorted(CASES))
def test_frame_with_another_azimuth_is_refused(frame, name):
    with pytest.raises(ValueError):
        CASES[name](*frame, ROTATED)
