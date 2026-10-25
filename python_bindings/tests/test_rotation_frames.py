"""Oblique analysis frames from Python: RotationFrame, rotate_to_frame,
rotate_back, fetch_dem(frame=...) argument checks, and
bearings_relative_to_true_north on the bound bearing outputs.

The fixtures are the ones the CLI tests use (tests/data/rotation and
tests/data/true_north). lakeland_line_oblique.tif and its frame file were
written by lsdtt-raster-preprocessing (rotate_to_linear_feature, the line,
width 1200 m), so rotate_to_frame on the same inputs must reproduce that
raster exactly: a .param run and the equivalent Python call give identical
rasters. The planar DEMs dip to TRUE azimuth 200 deg at their centre, with
the grid bearing set from pyproj's meridian convergence.

Authors
-------
S. M. Mudd
"""

from pathlib import Path

import numpy as np
import pytest

import lsdtt3

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
ROTATION = DATA / "rotation"
TRUE_NORTH = DATA / "true_north"

pytestmark = pytest.mark.skipif(
    not (ROTATION / "lakeland_line_oblique.tif").exists(),
    reason="repository test data not available",
)


@pytest.fixture(scope="module")
def dem():
    return lsdtt3.Raster.read(str(DATA / "flat_resolution_lakeland_crop_dem.tif"))


@pytest.fixture(scope="module")
def line_frame(dem):
    return lsdtt3.RotationFrame.from_feature(
        str(ROTATION / "lakeland_line.geojson"), width=1200.0, source=dem)


@pytest.fixture(scope="module")
def cli_rotated():
    return lsdtt3.Raster.read(str(ROTATION / "lakeland_line_oblique.tif"))


def _same_crs(a, b):
    from pyproj import CRS
    return CRS.from_user_input(a) == CRS.from_user_input(b)


# --- frames ------------------------------------------------------------------

def test_line_frame_built_for_a_source_dem(line_frame):
    f = line_frame
    assert f.frame_source == "line"
    assert (f.n_cols, f.n_rows) == (144, 40)
    assert f.cell_size == 30.0
    assert f.azimuth_degrees == pytest.approx(50.127289171, abs=1e-8)
    assert f.extent == (-2160.0, 2160.0, -600.0, 600.0)
    assert f.width_metres == 1200.0
    assert f.end_buffer_metres == 300.0          # -1 = 10 cells
    # Rotate-back target: the source DEM's own grid (UTM 15N, 30 m, its corner).
    assert f.target_epsg == 32615
    assert f.target_cell_size == 30.0
    assert f.target_origin == (564930.0, 5600760.0)
    assert f.feature_fname == "lakeland_line.geojson"
    assert f.grid_convergence_degrees(0.0, 0.0) == pytest.approx(90.0 - f.azimuth_degrees, abs=1e-6)
    assert "RotationFrame(line" in repr(f)


def test_frame_file_from_the_cli_matches_the_python_frame(line_frame):
    f = lsdtt3.RotationFrame.read(str(ROTATION / "lakeland_line_oblique_rotation_frame.txt"))
    assert f.azimuth_degrees == pytest.approx(line_frame.azimuth_degrees, abs=1e-12)
    assert f.extent == line_frame.extent
    assert f.centre_longitude == pytest.approx(line_frame.centre_longitude, abs=1e-12)
    assert f.target_epsg == line_frame.target_epsg
    assert f.target_origin == line_frame.target_origin
    assert f.source_dem_fname == "flat_resolution_lakeland_crop_dem.tif"


def test_frame_file_round_trip(line_frame, tmp_path):
    path = tmp_path / "x_rotation_frame.txt"
    line_frame.write(str(path))
    back = lsdtt3.RotationFrame.read(str(path))
    assert back.azimuth_degrees == line_frame.azimuth_degrees
    assert back.extent == line_frame.extent
    assert back.target_origin == line_frame.target_origin
    assert back.frame_source == "line"
    assert "created_by: 'lsdtt3 (Python)'" in path.read_text()


def test_fetch_frame_defaults_without_a_source():
    f = lsdtt3.RotationFrame.from_feature(str(ROTATION / "lakeland_line.geojson"), width=1200.0)
    assert f.cell_size == 30.0                   # the fetch_dem default
    assert f.target_epsg == 32615                # UTM zone of the frame centre
    assert f.target_origin == (0.0, 0.0)         # the fetch (-tap) lattice
    g = lsdtt3.RotationFrame.from_feature(
        str(ROTATION / "lakeland_line.geojson"), width=1200.0, grid_spacing=10.0,
        target_epsg=32616)
    assert (g.cell_size, g.target_epsg) == (10.0, 32616)


def test_polygon_frame(dem):
    f = lsdtt3.RotationFrame.from_feature(str(ROTATION / "lakeland_polygon.geojson"), source=dem)
    assert f.frame_source == "polygon"
    assert f.buffer_all_sides
    assert (f.n_cols, f.n_rows) == (122, 52)


def test_footprint_frame_of_a_square_dem_warns(dem):
    with pytest.warns(UserWarning, match="ill-defined"):
        f = lsdtt3.RotationFrame.from_valid_data(dem)
    assert f.frame_source == "data_footprint"
    assert (f.n_cols, f.n_rows) == (150, 150)


def test_frame_rebuilt_from_a_frame_raster_warns(cli_rotated, line_frame):
    with pytest.warns(UserWarning):
        f = lsdtt3.RotationFrame.from_raster(cli_rotated)
    assert f.frame_source == "crs"
    assert f.azimuth_degrees == pytest.approx(line_frame.azimuth_degrees, abs=1e-9)
    assert f.extent == line_frame.extent


def test_frame_refusals(dem):
    with pytest.raises(ValueError):
        lsdtt3.RotationFrame.from_feature(str(ROTATION / "missing.geojson"), width=100.0)
    with pytest.raises(ValueError, match="target_epsg"):
        lsdtt3.RotationFrame.from_feature(
            str(ROTATION / "lakeland_line.geojson"), width=1200.0, source=dem, target_epsg=32615)
    with pytest.raises(ValueError):
        lsdtt3.RotationFrame.from_raster(dem)    # UTM, not a frame


# --- rotation ----------------------------------------------------------------

def test_python_rotation_equals_the_cli_raster(dem, line_frame, cli_rotated):
    rotated = lsdtt3.rotate_to_frame(dem, line_frame)
    a, b = rotated.metadata, cli_rotated.metadata
    assert rotated.shape == cli_rotated.shape == (40, 144)
    assert (a.x_min, a.y_min, a.cell_size) == (b.x_min, b.y_min, b.cell_size)
    assert _same_crs(a.crs_wkt, b.crs_wkt)
    assert np.array_equal(rotated.mask_numpy(), cli_rotated.mask_numpy())
    m = rotated.mask_numpy()
    assert np.array_equal(rotated.to_numpy()[m], cli_rotated.to_numpy()[m])


def test_rotate_back_lands_on_the_source_lattice(cli_rotated, line_frame):
    back = lsdtt3.rotate_back(cli_rotated, line_frame, resampling="near")
    meta = back.metadata
    assert back.shape == (125, 136)              # as smoke_raster_preprocessing_rotate_back
    assert meta.cell_size == 30.0
    for offset in ((meta.x_min - 564930.0) / 30.0, (meta.y_min - 5600760.0) / 30.0):
        assert offset == pytest.approx(round(offset), abs=1e-6)


def test_ids_survive_the_round_trip(dem, line_frame):
    rows, cols = np.indices(dem.shape)
    ids = (rows // 8) * 100 + cols // 8
    id_raster = lsdtt3.Raster.from_numpy(ids.astype(np.float32), dem.mask_numpy(), dem.metadata)
    rotated = lsdtt3.rotate_to_frame(id_raster, line_frame, resampling="near")
    back = lsdtt3.rotate_back(rotated, line_frame, resampling="near")
    values = set(np.unique(back.to_numpy()[back.mask_numpy()]).astype(int))
    assert values <= set(np.unique(ids))
    assert len(values) > 50


def test_rotate_back_refuses_a_raster_not_in_the_frame(dem, line_frame):
    with pytest.raises(ValueError, match="frame CRS"):
        lsdtt3.rotate_back(dem, line_frame)


# --- fetch_dem argument checks (no network: all refused first) ----------------

def test_fetch_dem_frame_argument_checks(line_frame):
    fetch_frame = lsdtt3.RotationFrame.from_feature(
        str(ROTATION / "lakeland_line.geojson"), width=1200.0)
    with pytest.raises(ValueError, match="west/south/east/north"):
        lsdtt3.fetch_dem(-92.1, 50.5, -92.0, 50.6, frame=fetch_frame)
    with pytest.raises(ValueError, match="target_epsg"):
        lsdtt3.fetch_dem(frame=fetch_frame, target_epsg=32615)
    with pytest.raises(ValueError, match="cell_size"):
        lsdtt3.fetch_dem(frame=fetch_frame, grid_spacing=10.0)
    with pytest.raises(ValueError, match="size_pixels"):
        lsdtt3.fetch_dem(frame=fetch_frame, size_pixels=10)
    with pytest.raises(ValueError, match="cannot be combined"):
        lsdtt3.fetch_dem(frame=fetch_frame, fetch_from_linear_feature=True)
    with pytest.raises(ValueError, match="fetch_feature_fname"):
        lsdtt3.fetch_dem(fetch_from_linear_feature=True)
    with pytest.raises(ValueError, match="give west"):
        lsdtt3.fetch_dem()


# --- bearings_relative_to_true_north ------------------------------------------

def _plane(name):
    return lsdtt3.Raster.read(str(TRUE_NORTH / name))


def _centre(raster):
    a = raster.to_numpy()
    return float(a[a.shape[0] // 2, a.shape[1] // 2])


@pytest.mark.parametrize("name, grid_aspect", [
    ("plane_frame_true200.tif", 239.8727),
    ("plane_utm15n_true200.tif", 199.2717),
])
def test_aspect_of_a_plane_dipping_to_true_200(name, grid_aspect):
    dem = _plane(name)
    grid = lsdtt3.aspect(dem, lambda_metres=90.0)
    true = lsdtt3.aspect(dem, lambda_metres=90.0, bearings_relative_to_true_north=True)
    g = grid.to_numpy()[grid.mask_numpy()]
    t = true.to_numpy()[true.mask_numpy()]
    assert np.all(np.abs(g - grid_aspect) < 1e-3)
    assert np.all(np.abs(t - 200.0) < 0.02)
    metrics = lsdtt3.polyfit_metrics(dem, ["aspect_radians", "slope"], lambda_metres=90.0,
                                     bearings_relative_to_true_north=True)
    assert np.degrees(_centre(metrics["aspect_radians"])) == pytest.approx(200.0, abs=0.001)
    plain = lsdtt3.polyfit_metrics(dem, ["slope"], lambda_metres=90.0)
    assert np.array_equal(metrics["slope"].to_numpy(), plain["slope"].to_numpy())


def test_true_north_needs_a_crs():
    meta = lsdtt3.RasterMeta.from_shape(10, 10, cell_size=10.0)
    raster = lsdtt3.Raster.from_numpy(
        np.arange(100, dtype=np.float32).reshape(10, 10), np.ones((10, 10), bool), meta)
    with pytest.raises(ValueError, match="no CRS"):
        lsdtt3.aspect(raster, bearings_relative_to_true_north=True)


@pytest.fixture(scope="module")
def frame_network(cli_rotated):
    filled = lsdtt3.fill_pits(cli_rotated)
    flow = lsdtt3.FlowInfo(filled)
    network = flow.channel_network(flow.channel_scope(threshold_contributing_pixels=15))
    return filled, flow, network


def test_link_and_com_bearings_convert_at_their_own_location(frame_network, line_frame):
    _filled, _flow, network = frame_network
    grid = network.link_bearings()
    true = network.link_bearings(bearings_relative_to_true_north=True)
    assert len(grid) > 50
    c = np.array([line_frame.grid_convergence_degrees(x, y)
                  for x, y in zip(grid.representative_x, grid.representative_y)])
    diff = np.mod(true.bearing_deg - (grid.bearing_deg - c) + 180.0, 360.0) - 180.0
    assert np.nanmax(np.abs(diff)) < 1e-6
    assert np.array_equal(true.direction_x, grid.direction_x)

    grid = network.basin_com_bearings()
    true = network.basin_com_bearings(bearings_relative_to_true_north=True)
    c = np.array([line_frame.grid_convergence_degrees((ox + cx) / 2, (oy + cy) / 2)
                  for ox, oy, cx, cy in zip(grid.outlet_x, grid.outlet_y,
                                            grid.centre_of_mass_x, grid.centre_of_mass_y)])
    diff = np.mod(true.bearing_deg - (grid.bearing_deg - c) + 180.0, 360.0) - 180.0
    assert np.nanmax(np.abs(diff)) < 1e-6


def test_basin_orientation_converts_at_the_centroid(frame_network, line_frame):
    filled, flow, network = frame_network
    heads = network.junctions()
    outlets = heads[heads.stream_order >= 2][["row", "col"]].to_numpy()[:5]
    basins = flow.basins(outlets)
    grid = basins.statistics(filled, include_shape_metrics=True)
    true = basins.statistics(filled, include_shape_metrics=True,
                             bearings_relative_to_true_north=True)
    c = np.array([line_frame.grid_convergence_degrees(x, y)
                  for x, y in zip(grid.centroid_x, grid.centroid_y)])
    diff = np.mod(true.orientation_degrees - (grid.orientation_degrees - c) + 90.0, 180.0) - 90.0
    assert np.max(np.abs(diff)) < 1e-6
    assert np.array_equal(true.solidity, grid.solidity)


def test_basin_surface_aspect_and_raster_orientations_shift(frame_network):
    filled, flow, network = frame_network
    heads = network.junctions()
    basins = flow.basins(heads[heads.stream_order >= 2][["row", "col"]].to_numpy()[:3])
    grid = lsdtt3.basin_surface_metrics(basins, filled, [90.0], ["aspect_degrees"],
                                        include_outline=False)["interior"]
    true = lsdtt3.basin_surface_metrics(basins, filled, [90.0], ["aspect_degrees"],
                                        include_outline=False,
                                        bearings_relative_to_true_north=True)["interior"]
    assert not np.allclose(grid["median"], true["median"])

    plane = _plane("plane_frame_true200.tif")
    roads = lsdtt3.find_roads(plane, smoothing_window_metres=90.0)
    roads_true = lsdtt3.find_roads(plane, smoothing_window_metres=90.0,
                                   bearings_relative_to_true_north=True)
    shift = _centre(roads.feature_orientation_degrees) - _centre(roads_true.feature_orientation_degrees)
    assert np.mod(shift, 180.0) == pytest.approx(39.87, abs=0.01)
