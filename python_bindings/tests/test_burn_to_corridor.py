"""burn_to_corridor and fetch_river_lines: the CLIs' routines from Python.

burn_to_corridor runs lsdtt-stream-burning's burn_to_corridor through the same
C++ routines (the option checks, the burn, the final conditioning, the checks,
the profile rows and the burn depth). On a synthetic valley with a dam and
noise pits it must give exactly the CLI's DEM and the CLI's profile table, when
the CLI binary is available (the packaged driver, or the repository's build
directory); the rest of the tests need only the binding.

fetch_river_lines reads the committed OS Open Rivers extract (no network) with
the parameters of the smoke_fetch_river_lines_local CTest, so it must give that
test's 24 lines in EPSG:32630, oriented downstream.
"""

import math
import os
import subprocess
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from pyproj import CRS
from shapely.geometry import LineString

import lsdtt3

pytestmark = [pytest.mark.unit]

REPO = Path(__file__).resolve().parents[2]
NR, NC, CS = 40, 80, 10.0
X0, Y0 = 500000.0, 4200000.0
FLOOR = 20
WKT = CRS.from_epsg(32630).to_wkt()


def _valley_dem():
    rows, cols = np.indices((NR, NC))
    z = 100.0 - 0.05 * CS * cols + 1.0 * np.abs(rows - FLOOR)
    z[:, 50] += np.where(np.abs(rows[:, 50] - FLOOR) <= 3, 4.0, 0.0)    # a 4 m dam across the floor
    z[FLOOR, 20] -= 1.5                                                    # noise pits on the floor
    z[FLOOR, 35] -= 2.0
    meta = lsdtt3.RasterMeta.from_shape(NR, NC, x_min=X0, y_min=Y0, cell_size=CS, crs_wkt=WKT, epsg_code=32630)
    return lsdtt3.Raster.from_numpy(z.astype(np.float32), np.ones((NR, NC), dtype=bool), meta)


def _lines():
    y = Y0 + (NR - FLOOR - 0.5) * CS
    return gpd.GeoDataFrame({"name": ["Floor Beck"]},
                            geometry=[LineString([(X0 + 1.5 * CS, y), (X0 + (NC - 0.5) * CS, y)])], crs="EPSG:32630")


def _cli_binary():
    try:
        from lsdtt3._cli import _driver_path
        return _driver_path("lsdtt-stream-burning")
    except Exception:  # noqa: BLE001 - an editable build has no packaged drivers
        p = REPO / "build" / "lsdtt-stream-burning"
        return p if p.is_file() else None


@pytest.fixture(scope="module")
def run():
    return lsdtt3.burn_to_corridor(_valley_dem(), _lines(), vector_name_field="name")


@pytest.mark.property
def test_burn_drains_and_reports(run):
    r = run.report
    assert r["n_fitted"] == 1 and r["n_skipped"] == 0 and r["n_paths"] == 1
    assert r["n_paths_rising_after_burn"] == 0 and r["n_paths_rising_final"] == 0
    assert r["final_line_cells_changed"] == 0 and r["final_other_cells_changed"] == 0
    z = run.profiles["z_final"].to_numpy()
    assert np.all(np.diff(z) < 0)                                  # falls downstream
    assert run.profiles["name"].eq("Floor Beck").all()
    assert {"path_id", "z_original", "z_profile", "z_final", "pit_filled", "starts_on_path_id"} <= set(run.profiles)
    assert len(run.channels) == 1 and run.channels.crs.equals(CRS.from_wkt(WKT))
    dam = run.profiles.loc[(run.profiles["x"] - X0) // CS == 50]
    assert (dam["z_original"] - dam["z_final"]).iloc[0] > 3.0      # the dam is cut
    depth = run.burn_depth.to_numpy()
    assert np.nanmax(depth) > 3.0


@pytest.mark.metamorphic
def test_lines_in_another_crs_are_reprojected(run):
    other = lsdtt3.burn_to_corridor(_valley_dem(), _lines().to_crs("EPSG:4326"), vector_name_field="name")
    assert np.array_equal(other.dem.to_numpy(), run.dem.to_numpy())


@pytest.mark.negative
def test_refusals():
    dem, lines = _valley_dem(), _lines()
    with pytest.raises(ValueError, match="burn_profile_fit_algorithm_choice"):
        lsdtt3.burn_to_corridor(dem, lines, burn_profile_fit_algorithm_choice="cubic")
    with pytest.raises(ValueError, match="burn_profile_max_pit_pixels"):
        lsdtt3.burn_to_corridor(dem, lines, burn_profile_max_pit_pixels=-1)
    with pytest.raises(ValueError, match="no_such_field"):
        lsdtt3.burn_to_corridor(dem, lines, vector_name_field="no_such_field")
    with pytest.raises(ValueError, match="breach_algorithm_choice"):
        lsdtt3.burn_to_corridor(dem, lines, breach_algorithm_choice="bogus")


@pytest.mark.metamorphic
def test_same_dem_and_profiles_as_the_cli(run, tmp_path):
    exe = _cli_binary()
    if exe is None:
        pytest.skip("lsdtt-stream-burning binary not available")
    _valley_dem().write(str(tmp_path / "valley.tif"))
    _lines().to_file(tmp_path / "lines.fgb", driver="FlatGeobuf")
    (tmp_path / "burn.param").write_text(
        f"read_path: {tmp_path}/\nread_fname: valley.tif\nwrite_path: {tmp_path}/\nwrite_fname: cli\n"
        f"stream_vector_path: {tmp_path}/lines.fgb\nstream_burn_mode: burn_to_corridor\n"
        "run_stream_burn: true\nvector_name_field: name\nwrite_burned_outputs: true\n"
        "write_burned_profiles_csv: true\nn_threads: 1\n")
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GDAL_", "PROJ_"))}
    try:   # the environment the wheel's own console scripts use (DLL dirs on Windows)
        from lsdtt3._cli import _child_environment
        env = _child_environment(env)
    except ImportError:
        pass
    subprocess.run([str(exe), str(tmp_path / "burn.param")], check=True, capture_output=True, env=env)
    cli_dem = lsdtt3.Raster.read(str(tmp_path / "cli_burned.tif")).to_numpy()
    assert np.array_equal(cli_dem, run.dem.to_numpy())
    cli = pd.read_csv(tmp_path / "cli_burned_profiles.csv")
    py = run.profiles
    assert list(cli.columns) == list(py.columns) and len(cli) == len(py)
    for c in ("z_original", "z_cleaned", "z_profile", "z_final", "distance_m", "x", "y"):
        assert np.allclose(cli[c].to_numpy(), py[c].to_numpy(), rtol=0, atol=1e-6), c
    for c in ("path_id", "reach_id", "pit_filled", "is_record_low", "is_rejected_low"):
        assert np.array_equal(cli[c].to_numpy(), py[c].to_numpy()), c


@pytest.mark.characterisation
def test_fetch_river_lines_local_extract():
    data = REPO / "tests" / "data" / "river_lines" / "os_open_rivers_extract" / "Data" / "oprvrs_gb.gpkg"
    if not data.is_file():
        pytest.skip("committed Open Rivers extract not found")
    r = lsdtt3.fetch_river_lines(
        river_line_source="local_file", west=-2.52, south=53.01, east=-2.32, north=53.12,
        river_lines_fname=str(data), river_lines_layer="watercourse_link", feature_class_field_name="form",
        direction_field_name="flow_direction", direction_field_mapping="in_direction=1,in_opposite_direction=-1",
        start_node_field_name="start_node", end_node_field_name="end_node")
    assert len(r.lines) == 24 and r.lines.crs.to_epsg() == 32630
    assert {"line_id", "feature_class", "from_node", "to_node", "strahler_order", "length_m"} <= set(r.lines)
    assert len(r.nodes) > 0 and len(r.topology) >= 1
    assert "local file" in r.attribution and "oprvrs_gb.gpkg" in r.attribution
    with pytest.raises(ValueError, match="river_line_source"):
        lsdtt3.fetch_river_lines(river_line_source="")
    with pytest.raises(ValueError, match="area"):
        lsdtt3.fetch_river_lines(river_line_source="local_file", river_lines_fname=str(data))
