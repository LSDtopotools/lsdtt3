# SPDX-License-Identifier: AGPL-3.0-or-later
"""Point snapping (FlowInfo.snap_points) and RasterMeta coordinate helpers.

The synthetic DEM is a V-shaped valley that drains south: every hillside
pixel flows diagonally onto the centre column, which carries all of the
drainage to the southern (base-level) edge. Drainage area therefore grows
down the centre column and is one pixel on the valley sides' outer edges.

Authors
-------
S. M. Mudd
"""

import inspect
import pathlib
import re

import geopandas as gpd
import numpy as np
import pytest
from pyproj import CRS, Transformer
from shapely.geometry import LineString, Point

import lsdtt3

N_ROWS, N_COLS, CELL = 21, 21, 10.0
X_MIN, Y_MIN = 500000.0, 4000000.0
CENTRE = N_COLS // 2
EPSG = 32630
UNDESCRIBED = "Not yet described."


def _meta(crs=True):
    return lsdtt3.RasterMeta.from_shape(
        N_ROWS, N_COLS, x_min=X_MIN, y_min=Y_MIN, cell_size=CELL,
        crs_wkt=CRS.from_epsg(EPSG).to_wkt() if crs else "", epsg_code=EPSG if crs else 0,
    )


def _valley(mask=None):
    rows, cols = np.indices((N_ROWS, N_COLS), dtype=np.float32)
    z = (N_ROWS - rows) * 1.0 + 0.5 * np.abs(cols - CENTRE)
    valid = np.ones((N_ROWS, N_COLS), dtype=bool) if mask is None else mask
    return lsdtt3.Raster.from_numpy(z, valid, _meta())


@pytest.fixture(scope="module")
def flow():
    return lsdtt3.FlowInfo(_valley(), boundary_conditions="nnbn", n_threads=1)


@pytest.fixture(scope="module")
def contributing(flow):
    return flow.contributing_pixels().to_numpy()


def _centre_xy(row, col):
    return X_MIN + (col + 0.5) * CELL, Y_MIN + N_ROWS * CELL - (row + 0.5) * CELL


def test_valley_drains_onto_the_centre_column(contributing):
    """The fixture is what the snapping tests assume."""
    assert contributing[0, 0] == 1
    assert np.all(np.diff(contributing[:, CENTRE]) > 0)
    assert contributing[-1, CENTRE] == contributing[:, CENTRE].max()


def _walk_downstream(flow, contributing, row, col, threshold):
    """The CLIs' downstream rule, written out: follow receiver_of until the
    contributing-pixel count reaches the threshold (>=)."""
    while contributing[row, col] < threshold:
        receiver = flow.receiver_of(row, col)
        if receiver is None:
            return None
        row, col = receiver
    return row, col


def test_downstream_matches_the_d8_receiver_walk(flow, contributing):
    """'downstream' is the CLI basin_outlet_snapping_method: downstream rule."""
    starts = [(r, c) for r in range(0, N_ROWS, 3) for c in range(0, N_COLS, 2)]
    n_reached = n_outlet = 0
    for threshold in (1, 7, 20, int(contributing[-1, CENTRE])):
        for r, c in starts:
            expected = _walk_downstream(flow, contributing, r, c, threshold)
            if expected is None:
                # The walk reached a base-level outlet first.
                n_outlet += 1
                with pytest.raises(ValueError, match="outlet_reached_before_threshold"):
                    flow.snap_points(*_centre_xy(r, c), threshold_contributing_pixels=threshold)
                continue
            n_reached += 1
            snapped = flow.snap_points(*_centre_xy(r, c), threshold_contributing_pixels=threshold)
            assert snapped.dtype == np.int64
            assert tuple(snapped[0]) == expected, (r, c, threshold)
    assert n_reached > 100 and n_outlet > 0
    # The default method is downstream, and the walk stops AT the threshold.
    row, col = 3, CENTRE
    exact = int(contributing[row, col])
    assert flow.snap_points(*_centre_xy(0, CENTRE), threshold_contributing_pixels=exact).tolist() \
        == [[row, col]]


def test_threshold_area_converts_exactly_to_contributing_pixels(flow, contributing):
    """N pixels qualify for an area A when N * pixel_area >= A."""
    x, y = _centre_xy(0, CENTRE)
    pixel_area = CELL * CELL
    for pixels in (5, 12):
        by_count = flow.snap_points(x, y, threshold_contributing_pixels=pixels)
        assert flow.snap_points(x, y, threshold_drainage_area_m2=pixels * pixel_area).tolist() \
            == by_count.tolist()
        # Just above pixels * pixel_area needs one more pixel.
        above = flow.snap_points(x, y, threshold_drainage_area_m2=pixels * pixel_area + 1e-6)
        assert above.tolist() == flow.snap_points(
            x, y, threshold_contributing_pixels=pixels + 1).tolist()
        assert contributing[tuple(above[0])] >= pixels + 1
    with pytest.raises(ValueError, match="not both"):
        flow.snap_points(x, y, threshold_contributing_pixels=5, threshold_drainage_area_m2=500.0)
    with pytest.raises(ValueError, match=">= 1"):
        flow.snap_points(x, y, threshold_contributing_pixels=0)
    with pytest.raises(ValueError, match="must be positive"):
        flow.snap_points(x, y, threshold_drainage_area_m2=-1.0)


CLI_DEFAULT_THRESHOLD = 1000


def test_cli_default_threshold_is_the_one_the_cli_declares():
    """The constant this file checks against is the CLI's own default."""
    header = (pathlib.Path(__file__).resolve().parents[2] / "src" / "core"
              / "lsdtt_common_parameters.hpp")
    if not header.exists():
        pytest.skip("needs the lsdtt3 source tree")
    found = re.search(r"default_basin_outlet_contributing_pixels_threshold\s*=\s*(\d+)",
                      header.read_text())
    assert found and int(found.group(1)) == CLI_DEFAULT_THRESHOLD


@pytest.mark.parametrize("method", ["downstream", "expanding_square", "largest_area"])
def test_without_a_threshold_every_method_uses_the_cli_default(method):
    """No threshold means basin_outlet_contributing_pixels_threshold's default."""
    size = 41   # the lower centre column drains more than 1000 pixels
    rows, cols = np.indices((size, size), dtype=np.float32)
    z = (size - rows) * 1.0 + 0.5 * np.abs(cols - size // 2)
    meta = lsdtt3.RasterMeta.from_shape(size, size, x_min=X_MIN, y_min=Y_MIN, cell_size=CELL,
                                        crs_wkt=CRS.from_epsg(EPSG).to_wkt(), epsg_code=EPSG)
    big = lsdtt3.FlowInfo(lsdtt3.Raster.from_numpy(z, np.ones(z.shape, bool), meta),
                          boundary_conditions="nnbn", n_threads=1)
    counts = big.contributing_pixels().to_numpy()
    start = 31   # 892 contributing pixels; row 34 is the first with >= 1000
    assert counts[33, size // 2] < CLI_DEFAULT_THRESHOLD <= counts[34, size // 2]
    x = X_MIN + (size // 2 + 0.5) * CELL
    y = Y_MIN + size * CELL - (start + 0.5) * CELL
    radius = dict(max_search_radius_pixels=5)
    default = big.snap_points(x, y, method=method, **radius)
    explicit = big.snap_points(x, y, method=method, **radius,
                               threshold_contributing_pixels=CLI_DEFAULT_THRESHOLD)
    assert default.tolist() == explicit.tolist()
    assert counts[tuple(default[0])] >= CLI_DEFAULT_THRESHOLD


def test_downstream_reports_the_status(flow, contributing):
    x, y = _centre_xy(0, 0)
    # The fixture has 441 pixels, so the default 1000 cannot be reached.
    with pytest.raises(ValueError, match="outlet_reached_before_threshold.*threshold 1000"):
        flow.snap_points(x, y)
    too_big = int(contributing.max()) + 1
    with pytest.raises(ValueError, match="outlet_reached_before_threshold"):
        flow.snap_points(x, y, threshold_contributing_pixels=too_big)


def test_expanding_square_lands_on_the_nearest_channel_pixel(flow, contributing):
    threshold_pixels = 20
    x, y = _centre_xy(5, CENTRE + 4)
    snapped = flow.snap_points(
        x, y, method="expanding_square", max_search_radius_pixels=6,
        threshold_drainage_area_m2=threshold_pixels * CELL * CELL,
    )
    row, col = snapped[0]
    assert contributing[row, col] >= threshold_pixels
    channel = np.argwhere(contributing >= threshold_pixels)
    rings = np.max(np.abs(channel - [5, CENTRE + 4]), axis=1)
    assert max(abs(row - 5), abs(col - CENTRE - 4)) == rings.min()
    # The input pixel itself is not a channel pixel.
    assert contributing[5, CENTRE + 4] < threshold_pixels
    # min_search_radius_pixels skips the inner rings, as in the CLIs.
    ring = max(abs(row - 5), abs(col - CENTRE - 4))
    outer = flow.snap_points(x, y, method="expanding_square", min_search_radius_pixels=ring + 1,
                             max_search_radius_pixels=10,
                             threshold_contributing_pixels=threshold_pixels)
    assert max(abs(outer[0, 0] - 5), abs(outer[0, 1] - CENTRE - 4)) == ring + 1


def test_expanding_square_requires_a_reachable_channel(flow):
    x, y = _centre_xy(0, 0)
    with pytest.raises(ValueError, match="threshold_not_found_in_search_window"):
        flow.snap_points(x, y, method="expanding_square", max_search_radius_pixels=2,
                         threshold_drainage_area_m2=50 * CELL * CELL)
    with pytest.raises(ValueError, match="min_search_radius_pixels"):
        flow.snap_points(x, y, method="expanding_square", min_search_radius_pixels=3,
                         max_search_radius_pixels=2, threshold_contributing_pixels=5)


def test_largest_area_takes_the_largest_area_in_the_window(flow, contributing):
    x, y = _centre_xy(10, CENTRE + 2)
    snapped = flow.snap_points(x, y, method="largest_area", max_search_radius_pixels=3,
                               threshold_contributing_pixels=1)
    row, col = snapped[0]
    window = contributing[7:14, CENTRE - 1:CENTRE + 6]
    assert contributing[row, col] == window.max()
    assert (row, col) == (13, CENTRE)
    # Radius 0 keeps the point in its own pixel.
    assert flow.snap_points(x, y, method="largest_area", max_search_radius_pixels=0,
                            threshold_contributing_pixels=1).tolist() == [[10, CENTRE + 2]]
    # The threshold accepts or rejects the largest pixel; it does not move it.
    largest = int(window.max())
    assert flow.snap_points(x, y, method="largest_area", max_search_radius_pixels=3,
                            threshold_contributing_pixels=largest).tolist() == [[13, CENTRE]]
    with pytest.raises(ValueError, match="threshold_not_found_in_search_window"):
        flow.snap_points(x, y, method="largest_area", max_search_radius_pixels=3,
                         threshold_contributing_pixels=largest + 1)


def test_largest_area_ties_go_to_the_nearer_pixel(flow, contributing):
    """Hillside pixels of one row have equal areas; the nearer one wins."""
    x, y = _centre_xy(0, 2)
    snapped = flow.snap_points(x, y, method="largest_area", max_search_radius_pixels=1,
                               threshold_contributing_pixels=1)
    window = contributing[0:2, 1:4]
    assert np.all(window[1] == window.max())   # three tied pixels in row 1
    assert snapped.tolist() == [[1, 2]]


def test_max_search_radius_metres_rounds_up_to_whole_pixels(flow):
    x, y = _centre_xy(10, CENTRE + 2)
    in_metres = flow.snap_points(x, y, method="largest_area", max_search_radius_metres=25.0,
                                 threshold_contributing_pixels=1)
    in_pixels = flow.snap_points(x, y, method="largest_area", max_search_radius_pixels=3,
                                 threshold_contributing_pixels=1)
    np.testing.assert_array_equal(in_metres, in_pixels)


def test_off_grid_and_nodata_points_raise_naming_each_point(flow):
    inside = _centre_xy(10, CENTRE)
    with pytest.raises(ValueError, match=r"1 of 2 points .*point 1 .*input_out_of_bounds"):
        flow.snap_points([inside[0], X_MIN - 5.0], [inside[1], inside[1]],
                         threshold_contributing_pixels=1)
    # The eastern boundary line belongs to no pixel.
    with pytest.raises(ValueError, match="input_out_of_bounds"):
        flow.snap_points(X_MIN + N_COLS * CELL, inside[1], threshold_contributing_pixels=1)

    mask = np.ones((N_ROWS, N_COLS), dtype=bool)
    mask[0, 0] = False
    holed = lsdtt3.FlowInfo(_valley(mask), boundary_conditions="nnbn", n_threads=1)
    # downstream needs a valid starting pixel; the window methods do not.
    with pytest.raises(ValueError, match=r"point 0 at pixel \(0, 0\): input_invalid_flow_cell"):
        holed.snap_points(*_centre_xy(0, 0), threshold_contributing_pixels=1)
    moved = holed.snap_points(*_centre_xy(0, 0), method="largest_area", max_search_radius_pixels=1,
                              threshold_contributing_pixels=1)
    assert moved.tolist() != [[0, 0]]


def test_input_forms_and_argument_errors(flow):
    x, y = _centre_xy(12, CENTRE)
    to_lonlat = Transformer.from_crs(EPSG, 4326, always_xy=True)
    longitude, latitude = to_lonlat.transform(x, y)
    expected = [[12, CENTRE]]
    keep = dict(method="largest_area", max_search_radius_pixels=0, threshold_contributing_pixels=1)
    assert flow.snap_points(x, y, **keep).tolist() == expected
    assert flow.snap_points(longitude=longitude, latitude=latitude, **keep).tolist() == expected
    # x/y in another CRS, and a GeoDataFrame with its own CRS.
    assert flow.snap_points(longitude, latitude, crs="EPSG:4326", **keep).tolist() == expected
    points = gpd.GeoDataFrame(geometry=[Point(longitude, latitude)], crs="EPSG:4326")
    assert flow.snap_points(points=points, **keep).tolist() == expected
    assert flow.snap_points(points=points.geometry, **keep).tolist() == expected

    with pytest.raises(ValueError, match="exactly one of"):
        flow.snap_points(x, y, longitude=longitude, latitude=latitude, **keep)
    with pytest.raises(ValueError, match="exactly one of"):
        flow.snap_points(**keep)
    with pytest.raises(ValueError, match="both x and y"):
        flow.snap_points(x, **keep)
    with pytest.raises(ValueError, match="method must be"):
        flow.snap_points(x, y, method="nearest_channel")
    with pytest.raises(ValueError, match="max_search_radius_pixels"):
        flow.snap_points(x, y, method="largest_area", max_search_radius_pixels=-1)
    with pytest.raises(ValueError, match="min_search_radius_pixels"):
        flow.snap_points(x, y, method="largest_area", min_search_radius_pixels=-1)
    lines = gpd.GeoDataFrame(geometry=[LineString([(x, y), (x + CELL, y)])], crs=EPSG)
    with pytest.raises(ValueError, match="must be a Point"):
        flow.snap_points(points=lines, **keep)
    with pytest.raises(ValueError, match="crs is not used with points"):
        flow.snap_points(points=points, crs="EPSG:4326", **keep)
    assert flow.snap_points([], [], threshold_contributing_pixels=1).shape == (0, 2)


def test_snapped_outlets_feed_basins(flow, contributing):
    xs, ys = zip(*(_centre_xy(r, c) for r, c in ((4, CENTRE + 1), (15, CENTRE - 2))))
    outlets = flow.snap_points(list(xs), list(ys), method="largest_area", max_search_radius_pixels=2,
                               threshold_contributing_pixels=1)
    basins = flow.basins(outlets)
    table = basins.table
    assert table["basin_id"].tolist() == [0, 1]
    assert table["outlet_row"].tolist() == outlets[:, 0].tolist()
    assert table["n_pixels"].tolist() == [contributing[r, c] for r, c in outlets]


def test_xy_row_col_round_trip_and_north_edge_orientation():
    meta = _meta()
    rows, cols = np.indices((N_ROWS, N_COLS))
    x, y = meta.row_col_to_xy(rows, cols)
    assert x.shape == rows.shape and x.dtype == np.float64
    back_rows, back_cols = meta.xy_to_row_col(x, y)
    np.testing.assert_array_equal(back_rows, rows)
    np.testing.assert_array_equal(back_cols, cols)

    # Row 0 is the northern edge: y falls as row grows.
    assert meta.row_col_to_xy(0, 0) == (X_MIN + 0.5 * CELL, meta.y_max - 0.5 * CELL)
    assert np.all(np.diff(y[:, 0]) < 0)
    assert meta.xy_to_row_col(X_MIN + 1.0, meta.y_max - 0.1) == (0, 0)
    assert meta.xy_to_row_col(X_MIN + 1.0, meta.y_min + 0.1) == (N_ROWS - 1, 0)
    row, col = meta.xy_to_row_col(X_MIN + 1.0, meta.y_max - 0.1)
    assert type(row) is int and type(col) is int


def test_xy_to_row_col_marks_off_grid_points():
    meta = _meta()
    rows, cols = meta.xy_to_row_col(
        [X_MIN - 1.0, meta.x_max, X_MIN + 1.0, X_MIN + 1.0, np.nan],
        [meta.y_max - 1.0, meta.y_max - 1.0, meta.y_max + 1.0, meta.y_min, 0.0],
    )
    assert rows.tolist() == [-1] * 5 and cols.tolist() == [-1] * 5
    # A scalar broadcasts against an array.
    rows, cols = meta.xy_to_row_col(X_MIN + 1.0, [meta.y_max - 1.0, meta.y_max - 11.0])
    assert rows.tolist() == [0, 1] and cols.tolist() == [0, 0]


def test_lonlat_to_row_col_uses_the_raster_crs():
    meta = _meta()
    x, y = meta.row_col_to_xy([3, 17], [5, 9])
    longitude, latitude = Transformer.from_crs(EPSG, 4326, always_xy=True).transform(x, y)
    rows, cols = meta.lonlat_to_row_col(longitude, latitude)
    assert rows.tolist() == [3, 17] and cols.tolist() == [5, 9]
    assert meta.lonlat_to_row_col(0.0, 0.0) == (-1, -1)
    with pytest.raises(ValueError, match="no CRS"):
        _meta(crs=False).lonlat_to_row_col(longitude, latitude)
    with pytest.raises(ValueError, match="not a CRS"):
        meta.lonlat_to_row_col(longitude, latitude, crs="not a crs")


@pytest.mark.parametrize("owner, member", [
    ("FlowInfo", "snap_points"),
    ("RasterMeta", "xy_to_row_col"),
    ("RasterMeta", "row_col_to_xy"),
    ("RasterMeta", "lonlat_to_row_col"),
])
def test_new_helpers_are_fully_documented(owner, member):
    value = getattr(getattr(lsdtt3, owner), member)
    doc = inspect.cleandoc(inspect.unwrap(value).__doc__ or "")
    assert UNDESCRIBED not in doc
    for heading in ("Parameters", "Returns", "Examples"):
        assert f"\n{heading}\n{'-' * len(heading)}\n" in doc, heading
    parameters = doc.split("\nParameters\n", 1)[1].split("\nReturns\n", 1)[0]
    definitions = [definition for definition, _doc, _defaults in inspect.unwrap(value).__nb_signature__]
    for name in re.findall(r"(?:\(|, )([A-Za-z_]\w*)(?=:)", definitions[0]):
        if name != "self":
            assert re.search(rf"(?<![\w]){name}\b", parameters), name
    # The class docstring lists the new member.
    assert re.search(rf"^{member}\(", inspect.getdoc(getattr(lsdtt3, owner)), re.MULTILINE)


def test_no_export_has_an_undescribed_entry():
    """Overloads merge their per-form Parameters, so no placeholder is left."""
    found = [
        name for name in lsdtt3.__all__
        if UNDESCRIBED in (inspect.getdoc(getattr(lsdtt3, name)) or "")
    ]
    assert not found, found
    doc = inspect.getdoc(lsdtt3.fill_pits_with_sinks)
    assert "sinks : Raster or array-like of int" in doc
    assert "Call form 2: N x 2 array" in doc
    assert inspect.getdoc(lsdtt3.classify_surface_shape).count("\nParameters\n----------\n") == 1
