#!/usr/bin/env python3
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Real-terrain routing gate for fill_flat_resolution_choice (network tier).

Runs `lsdtt-drainage-analysis` (raster_is_filled: false, fill_min_slope 1e-4,
D8 drainage area) on fetched Copernicus GLO-30 fixtures with each of the three
flat-handling choices -- inline_gradient (the default, a control),
barnes2014b_float_step and barnes2014b_gradient -- and asserts that every run
exits 0. For the two barnes2014b runs it also asserts the final
`Rim cells without descent: 0` line, and prints the rounds and the per-round
counts of cells without descent, and the final priority-flood pass.

Why: a single Barnes et al. (2014b) pass (RichDEM ResolveFlats) leaves rim
cells and low edges without a lower neighbour on lake-rich real terrain, and
the D8 router then refuses the surface with "could not assign a receiver".
The synthetic fixtures of unit_flat_resolution did not show it. The committed
crop (unit_flat_resolution_real_terrain and the
rigorous_drainage_analysis_flat_resolution_real_terrain_* tests) is the
hermetic gate; this script runs the same command on whole fetched fixtures.

Fixtures, in two groups:

  * the flow-benchmark registry (flow_benchmark_fixtures.py): alpine_small,
    lakeland_small, alpine_medium, lakeland_medium (Italian Alps / Finnish
    Lakeland). Fetched with lsdtt-fetch-raster when not cached (unless
    --no-fetch); a fixture that cannot be obtained is a SKIP (exit 77).
  * the benchmark-suite ladder (run_benchmark_suite.py): alpine_500,
    lakeland_500, alpine_1000, lakeland_1000 (Himalaya / Ontario lakeland),
    which are much richer in depressions. They are cropped from a source
    raster by the suite (with osgeo), so this stdlib-only script cannot make
    them: they are used when cached in the fixture directory and listed as
    "not cached" otherwise.

Exit codes: 0 all runs routed; 1 a run failed; 77 a fixture could not be
obtained, or none was available.

Usage:
    python3 tests/scripts/check_flat_resolution_real_terrain_routing.py \\
        --bin-dir build --work-dir /tmp/frr
    python3 tests/scripts/check_flat_resolution_real_terrain_routing.py \\
        --bin-dir build --work-dir /tmp/frr \\
        --fixture-dir /tmp/lsdtt3_benchmark_fixtures --no-fetch \\
        --fixtures alpine_small,lakeland_small,alpine_500,lakeland_500
"""

from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import sys

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))  # python -I does not add it
import lsdtt_test_env  # noqa: E402
import flow_benchmark_fixtures as fbf  # noqa: E402  (stdlib only)

# The cache run_flow_benchmark.py and run_benchmark_suite.py share.
DEFAULT_FIXTURE_DIR = SCRIPT_DIR.parent / "data" / "dem_cache"

FETCHABLE = ["alpine_small", "lakeland_small", "alpine_medium", "lakeland_medium"]
LADDER = ["alpine_500", "lakeland_500", "alpine_1000", "lakeland_1000"]

CHOICES = ["inline_gradient", "barnes2014b_float_step", "barnes2014b_gradient"]

PARAM_TEMPLATE = """read_path: {read_path}
read_fname: {read_fname}
write_path: {write_path}
write_fname: {write_fname}

raster_is_filled: false
fill_flat_resolution_choice: {choice}
fill_min_slope: 0.0001
write_output_format: gtiff

write_d8_drainage_area_raster: true
write_dinf_drainage_area_raster: false
write_mfd_drainage_area_raster: false
write_twi_raster: false
write_connectivity_raster: false
write_height_above_channel_raster: false
"""


def cli_env() -> dict[str, str]:
    """Environment for an lsdtt3 CLI started from a Pixi/conda Python.

    lsdtt_test_env.system_env() drops the environment's GDAL/PROJ library and
    data paths. PROJ_NETWORK is dropped as well: an activated Pixi environment
    sets PROJ_NETWORK=ON, and a system-GDAL CLI that also picked up the
    environment's GDAL plugins ran about 20x slower per call when launched
    from Pixi Python. Nothing here needs PROJ's network grids.
    """
    env = lsdtt_test_env.system_env()
    env.pop("PROJ_NETWORK", None)
    for key in [k for k in env if k.startswith("CONDA_")]:
        env.pop(key, None)
    return env


def line_value(stdout: str, prefix: str) -> str | None:
    for line in stdout.splitlines():
        s = line.strip()
        if s.startswith(prefix):
            return s[len(prefix):].strip()
    return None


def run_one(exe: pathlib.Path, fixture: pathlib.Path, choice: str,
            work_dir: pathlib.Path) -> tuple[bool, str]:
    out_dir = work_dir / f"{fixture.stem}_{choice}"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    stem = "frr"
    (out_dir / "run.param").write_text(PARAM_TEMPLATE.format(
        read_path=fixture.parent, read_fname=fixture.name,
        write_path=out_dir, write_fname=stem, choice=choice))
    proc = lsdtt_test_env.run([str(exe), str(out_dir), "run.param"],
                              capture_output=True, text=True, env=cli_env(),
                              cwd=str(out_dir))
    (out_dir / "stdout.txt").write_text(proc.stdout)
    (out_dir / "stderr.txt").write_text(proc.stderr)

    problems = []
    if proc.returncode != 0:
        err_lines = [l for l in proc.stderr.splitlines() if l.strip()]
        problems.append(f"exit {proc.returncode}"
                        + (f": {err_lines[-1]}" if err_lines else ""))
    if not (out_dir / f"{stem}_d8_drainage_area.tif").exists():
        problems.append("no D8 drainage area raster written")

    summary = ""
    if choice != "inline_gradient":
        rim = line_value(proc.stdout, "Rim cells without descent:")
        rounds = line_value(proc.stdout, "Rounds:")
        final = line_value(proc.stdout, "Final priority-flood pass:")
        if rim != "0":
            problems.append(f"'Rim cells without descent:' is {rim!r}, want '0'")
        if rounds is None:
            problems.append("no 'Rounds:' line")
        summary = f"rounds {rounds}; final pass {final}"
        m = re.search(r"Compute time: ([0-9.eE+-]+) s", proc.stdout)
        if m:
            summary += f"; fill + resolve {float(m.group(1)):.3f} s"
    ok = not problems
    return ok, (summary if ok else "; ".join(problems) + (f" ({summary})" if summary else ""))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", required=True,
                        help="directory holding lsdtt-drainage-analysis (and lsdtt-fetch-raster)")
    parser.add_argument("--work-dir", required=True, help="scratch directory for the runs")
    parser.add_argument("--fixture-dir", default=None,
                        help=f"where fixtures are cached (default {DEFAULT_FIXTURE_DIR})")
    parser.add_argument("--fixtures", default=None,
                        help="comma-separated fixture keys (default: "
                             + ",".join(FETCHABLE + LADDER) + ")")
    parser.add_argument("--no-fetch", action="store_true",
                        help="never download; a missing registry fixture is a skip")
    args = parser.parse_args()

    bin_dir = pathlib.Path(args.bin_dir).resolve()
    exe = bin_dir / "lsdtt-drainage-analysis"
    if not exe.exists():
        print(f"FAIL: {exe} not found", file=sys.stderr)
        return 1
    fetch_exe = bin_dir / "lsdtt-fetch-raster"
    fixture_dir = (pathlib.Path(args.fixture_dir).resolve() if args.fixture_dir
                   else DEFAULT_FIXTURE_DIR)
    work_dir = pathlib.Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    keys = args.fixtures.split(",") if args.fixtures else FETCHABLE + LADDER
    fixtures: list[pathlib.Path] = []
    for key in keys:
        cached = fixture_dir / f"{key}.tif"
        if cached.exists():
            fixtures.append(cached)
        elif key in LADDER:
            print(f"  {key}: not cached in {fixture_dir} (made by run_benchmark_suite.py); not run")
        elif args.no_fetch:
            fbf.die_skip(f"fixture {cached} is missing and --no-fetch was given")
        else:
            # exits 77 when the fixture cannot be fetched
            fixtures.append(fbf.ensure_fixture(key, fetch_exe, fixture_dir))
    if not fixtures:
        fbf.die_skip(f"no fixture available in {fixture_dir}")

    n_fail = 0
    for fixture in fixtures:
        for choice in CHOICES:
            ok, msg = run_one(exe, fixture, choice, work_dir)
            print(f"{'PASS' if ok else 'FAIL'}  {fixture.stem:15s} {choice:24s} {msg}")
            n_fail += 0 if ok else 1

    n_runs = len(fixtures) * len(CHOICES)
    print(f"\n{n_runs} runs on {len(fixtures)} fixtures, {n_fail} failed")
    print("RESULT: " + ("PASS" if n_fail == 0 else "FAIL"))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
