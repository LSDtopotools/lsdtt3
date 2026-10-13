"""Flow-benchmark DEM fixture registry and fetch (stdlib only).

Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3

Copyright (C) 2026 Simon M. Mudd
SPDX-License-Identifier: AGPL-3.0-or-later

The alpine / lakeland x small / medium / large fixtures of
run_flow_benchmark.py, and the fetch that caches them. Kept apart from that
driver so a stdlib-only test (check_flat_resolution_real_terrain_routing.py)
can fetch the same fixtures without importing a module that needs numpy and
GDAL's Python bindings. run_flow_benchmark.py re-exports everything here.
"""

from __future__ import annotations

import pathlib
import sys

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))  # python -I does not add it
import lsdtt_test_env  # noqa: E402

SKIP_EXIT = 77


def die_skip(msg: str):
    print(f"SKIP: {msg}", file=sys.stderr)
    sys.exit(SKIP_EXIT)


# -----------------------------------------------------------------------------
# Fixture registry
#
# Fixtures are fetched by the lsdtt3 CLI itself (`lsdtt-fetch-raster`), not by a
# Python helper -- the Python fetcher was retired once the C++ path was shown to
# produce byte-identical output. Two terrain regimes, because they
# diverge very differently: see docs/developer_notes/testing_and_benchmarking.md
#
# Each fixture is a centre point plus a pixel count. The bounding box is widened
# by a margin so the exact centre crop cannot run off the edge onto nodata after
# reprojection.
# -----------------------------------------------------------------------------
FIXTURES = {
    "alpine":   {"centre": (10.40, 46.50), "epsg": 32632,
                 "description": "high relief, Italian Alps (Val Venosta / Ortler)"},
    "lakeland": {"centre": (25.40, 61.50), "epsg": 32635,
                 "description": "low relief, Finnish Lakeland (many closed basins)"},
}

SIZES = {"small": 256, "medium": 720, "large": 2400}

GRID_SPACING_M = 30.0


def fixture_bbox(name: str, size_px: int):
    """WGS84 bbox for a fixture, wide enough to survive reprojection."""
    import math
    lon, lat = FIXTURES[name]["centre"]
    half_m = 0.5 * size_px * GRID_SPACING_M * 1.25
    dlat = half_m / 111_320.0
    dlon = dlat / max(math.cos(math.radians(lat)), 1e-6)
    return (lon - dlon, lat - dlat, lon + dlon, lat + dlat)


def ensure_fixture(fixture_key: str, fetch_exe: pathlib.Path,
                   fixture_dir: pathlib.Path) -> pathlib.Path:
    """Fetch a fixture into `fixture_dir` with lsdtt-fetch-raster if it is not
    already cached there."""
    out_path = fixture_dir / f"{fixture_key}.tif"
    if out_path.exists():
        return out_path

    if "_" not in fixture_key:
        die_skip(f"fixture name {fixture_key!r} should be <name>_<size>")
    name, size_name = fixture_key.rsplit("_", 1)

    if name not in FIXTURES or size_name not in SIZES:
        die_skip(f"unknown fixture {fixture_key!r}; "
                 f"names={sorted(FIXTURES)} sizes={sorted(SIZES)}")

    if not fetch_exe or not pathlib.Path(fetch_exe).exists():
        die_skip(f"fixture {out_path} is missing and lsdtt-fetch-raster was not "
                 f"found (looked for {fetch_exe}). Build it, or pass --fetch-exe.")

    size_px = SIZES[size_name]
    west, south, east, north = fixture_bbox(name, size_px)

    fixture_dir.mkdir(parents=True, exist_ok=True)
    param_dir = fixture_dir / "_fetch_params"
    param_dir.mkdir(parents=True, exist_ok=True)
    param_path = param_dir / f"{fixture_key}.param"

    param_path.write_text(
        f"write_path: {fixture_dir}\n"
        f"write_fname: {fixture_key}_raw\n"
        f"west: {west:.10f}\n"
        f"south: {south:.10f}\n"
        f"east: {east:.10f}\n"
        f"north: {north:.10f}\n"
        f"dem_source: cop30_aws\n"
        f"target_epsg: {FIXTURES[name]['epsg']}\n"
        f"grid_spacing: {GRID_SPACING_M}\n"
        f"size_pixels: {size_px}\n"
    )

    print(f"fetching fixture {fixture_key} ({FIXTURES[name]['description']})")
    proc = lsdtt_test_env.run([str(fetch_exe), str(param_dir), param_path.name],
                          capture_output=True, text=True, env=lsdtt_test_env.system_env())

    if proc.returncode != 0:
        # An unreachable source is a missing dependency, not a test failure.
        die_skip(f"could not fetch fixture {fixture_key}:\n{proc.stdout[-1500:]}\n"
                 f"{proc.stderr[-1500:]}")

    written = fixture_dir / f"{fixture_key}_raw_DEM.tif"
    if not written.exists():
        die_skip(f"lsdtt-fetch-raster reported success but {written} is missing")
    written.rename(out_path)
    return out_path
