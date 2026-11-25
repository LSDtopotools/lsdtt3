"""lsdtt-hillslope-channel-tagging's routines in Python: FlowInfo.hillslope_relief,
hillslope_chi_tags, channel_attributes_to_raster, tag_hillslopes_from_channel and
summarise_hillslopes_by_channel_pixel. The numerics are tested in the C++ driver
(tests/cpp/test_hillslope_channel_tagging.cpp); these carry the plumbing."""

import numpy as np
import pandas as pd
import pytest
from pyproj import CRS

import lsdtt3

pytestmark = [pytest.mark.unit]

WKT = CRS.from_epsg(32633).to_wkt()


def _raster(values, cell_size=10.0, mask=None):
    values = np.asarray(values, dtype=np.float32)
    meta = lsdtt3.RasterMeta.from_shape(values.shape[0], values.shape[1], x_min=500000.0, y_min=4000000.0,
                                        cell_size=cell_size, crs_wkt=WKT, epsg_code=32633)
    if mask is None:
        mask = np.ones(values.shape, dtype=bool)
    return lsdtt3.Raster.from_numpy(values, mask, meta)


def _valley(nr=21, nc=15, c0=7, gradient=0.5, cell_size=10.0):
    """A V-valley running south, with a nodata collar, and its channel in column c0."""
    r, c = np.mgrid[0:nr, 0:nc]
    z = 1000.0 + gradient * cell_size * np.abs(c - c0) + 0.01 * cell_size * (nr - 1 - r)
    mask = np.zeros((nr, nc), dtype=bool)
    mask[1:-1, 1:-1] = True
    z = np.where(mask, z, -9999.0)
    channel = np.where((c == c0) & mask, 1.0, 0.0)
    dem = _raster(z, cell_size, mask)
    return dem, _raster(channel, cell_size), lsdtt3.FlowInfo(dem, boundary_conditions="oooo"), (nr, nc, c0)


@pytest.mark.exact_solution
def test_relief_on_a_planar_valley():
    dem, channels, flow, (nr, nc, c0) = _valley()
    r = flow.hillslope_relief(dem, channels, ridge_method="highest")
    hand = r["hand"].to_numpy()
    row = nr // 2
    assert hand[row, c0 + 3] == pytest.approx(15.0, abs=0.05)
    relative = r["relative_height"].to_numpy()
    assert relative[row, c0 + 3] == pytest.approx(3.0 / 6.0, abs=0.02)
    assert r["n_hillslope"] + r["n_no_channel"] + r["n_channel"] == int(dem.mask_numpy().sum())
    assert r["channel_flat"].ndim == 1 and r["channel_flat"].dtype == np.int64
    assert r["channel_flat"][row * nc + c0 + 3] == row * nc + c0


@pytest.mark.property
def test_aspect_ray_needs_a_surface_and_bounds_ratio():
    dem, channels, flow, _ = _valley()
    with pytest.raises(ValueError, match="aspect_surface"):
        flow.hillslope_relief(dem, channels, ridge_method="aspect_ray")
    r = flow.hillslope_relief(dem, channels, ridge_method="aspect_ray", aspect_surface=dem, polyfit_lambda_metres=30.0)
    ratio = r["relative_height"].to_numpy()[r["relative_height"].mask_numpy().astype(bool)]
    assert ratio.size > 0 and ratio.min() >= 0.0 and ratio.max() <= 1.0
    assert r["n_crest_divide"] + r["n_ray_channel"] + r["n_ray_boundary"] + r["n_no_aspect"] == r["n_hillslope"]


@pytest.mark.exact_solution
def test_tag_recovers_the_channel_link():
    dem, channels, flow, (nr, nc, _) = _valley()
    r = flow.hillslope_relief(dem, channels)
    ramp = np.arange(nr * nc, dtype=np.float32)
    values = _raster(ramp.reshape(nr, nc), mask=channels.to_numpy() != 0)
    out = lsdtt3.tag_hillslopes_from_channel(r["channel_flat"], values, dem)
    tag = out["tag"].to_numpy().ravel()
    tagged = out["tag"].mask_numpy().astype(bool).ravel()
    assert tagged.sum() == out["n_hillslope_tagged"] + out["n_channel_tagged"]
    assert np.array_equal(tag[tagged], r["channel_flat"][tagged].astype(np.float32))


@pytest.mark.negative
def test_chi_tags_need_exactly_one_baselevel_setting():
    dem, channels, flow, (nr, nc, c0) = _valley()
    r = flow.hillslope_relief(dem, channels)
    with pytest.raises(ValueError):
        flow.hillslope_chi_tags(dem, r["channel_flat"], threshold_contributing_pixels=5)
    with pytest.raises(ValueError):
        flow.hillslope_chi_tags(dem, r["channel_flat"], outlet_flats=[(nr - 2) * nc + c0],
                                chi_max_outlet_elevation_metres=2000.0, threshold_contributing_pixels=5)


@pytest.mark.property
def test_chi_tags_equal_their_channel_pixel():
    dem, channels, flow, (nr, nc, c0) = _valley()
    r = flow.hillslope_relief(dem, channels)
    out = flow.hillslope_chi_tags(dem, r["channel_flat"], outlet_flats=[(nr - 2) * nc + c0],
                                  threshold_contributing_pixels=5)
    assert out["n_fitted_nodes"] > 5
    chi = out["chi"]["tag"].to_numpy().ravel()
    tagged = out["chi"]["tag"].mask_numpy().astype(bool).ravel()
    assert tagged.sum() > 50
    assert np.array_equal(chi[tagged], chi[r["channel_flat"][tagged]])


@pytest.mark.exact_solution
def test_attributes_from_a_table_and_from_a_file_are_identical(tmp_path):
    dem, channels, flow, (nr, nc, c0) = _valley()
    table = pd.DataFrame({"row": [5, 5, 12], "col": [c0, c0, c0], "flat": [5 * nc + c0, 5 * nc + c0, 12 * nc + c0],
                          "v": [1.0, 3.0, -9999.0]})
    path = tmp_path / "attributes.csv"
    table.to_csv(path, index=False)
    a = flow.channel_attributes_to_raster(table, channels, ["v"], max_gap_fill_distance_metres=0.0)
    b = flow.channel_attributes_to_raster(path, channels, ["v"], max_gap_fill_distance_metres=0.0)
    va, vb = a["values"]["v"], b["values"]["v"]
    assert np.array_equal(va.to_numpy(), vb.to_numpy()) and np.array_equal(va.mask_numpy(), vb.mask_numpy())
    assert va.to_numpy()[5, c0] == 2.0
    assert not va.mask_numpy()[12, c0]
    assert a["counts"]["v"]["n_multiple_records"] == 1


@pytest.mark.property
def test_attributes_from_a_geodataframe():
    gpd = pytest.importorskip("geopandas")
    from shapely.geometry import LineString

    dem, channels, flow, (nr, nc, c0) = _valley()
    x = 500000.0 + (c0 + 0.5) * 10.0
    y_top = 4000000.0 + (nr - 3 - 0.5 + 1.0) * 10.0
    y_bottom = y_top - 60.0
    frame = gpd.GeoDataFrame({"v": [4.0]}, geometry=[LineString([(x, y_top), (x, y_bottom)])], crs=WKT)
    out = flow.channel_attributes_to_raster(frame, channels, ["v"], max_snap_distance_metres=5.0,
                                            max_gap_fill_distance_metres=0.0)
    assert out["n_records_matched"] == 1 and out["counts"]["v"]["n_direct"] >= 5
    assert np.all(out["values"]["v"].to_numpy()[out["values"]["v"].mask_numpy().astype(bool)] == 4.0)


@pytest.mark.exact_solution
def test_summary_row_per_channel_pixel():
    dem, channels, flow, (nr, nc, c0) = _valley()
    r = flow.hillslope_relief(dem, channels, ridge_method="highest")
    frame = lsdtt3.summarise_hillslopes_by_channel_pixel(r["channel_flat"], r["relative_height"])
    assert len(frame) == r["n_channel"]
    assert frame["n_hillslope_pixels"].sum() == r["n_hillslope"]
    assert list(frame.columns) == ["flat", "row", "col", "latitude", "longitude", "n_hillslope_pixels",
                                   "n_with_value", "median"]
