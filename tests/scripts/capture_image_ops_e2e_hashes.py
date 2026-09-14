#!/usr/bin/env python3
# capture_image_ops_e2e_hashes.py
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# -----------------------------------------------------------------------------
# End-to-end half of the image-ops oracle (OpenCV removal, Phase 0).
#
# Runs the two CLIs whose outputs depend on the OpenCV image operations, on the
# fluvial-surfaces fixture DEM, and records a SHA-256 of every file they write:
#
#   lsdtt-flat-surface-finder  slope-break outputs (always run; OpenCV adaptive
#                              threshold inside find_vec_skeleton)
#   lsdtt-surface-metrics      Otsu, Hough lines, VEC adaptive threshold and
#                              VEC adaptive skeleton rasters
#
# Each CLI is run in two OpenCV modes, each twice (to detect files that are not
# byte-deterministic, e.g. embedded timestamps):
#
#   default    OpenCV as built (IPP + dispatched SIMD)
#   reference  OPENCV_IPP=disabled and OPENCV_CPU_DISABLE=<all dispatched
#              features>; baseline SIMD cannot be disabled at run time
#
# NOT part of CTest: run by hand while OpenCV is linked.
#
#   python3 tests/scripts/capture_image_ops_e2e_hashes.py \
#       --bin-dir build --fixture-dir build/tests/fluvial_surfaces_fixture \
#       --work-dir <scratch> --out tests/expected/image_ops/e2e_hashes.csv
# -----------------------------------------------------------------------------

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

REFERENCE_ENV = {
    "OPENCV_IPP": "disabled",
    "OPENCV_CPU_DISABLE": "SSE4_1,SSE4_2,POPCNT,AVX,FP16,AVX2,FMA3,AVX512_SKX",
}

FLAT_SURFACE_FINDER_PARAMS = """\
read_path: {read_path}
read_fname: fluvsurf_dem.tif
write_path: {write_path}
write_fname: fsf

n_threads: 1
write_output_format: gtiff

relief_window_side_pixels: 5
relief_threshold_metres: 2.0
adaptive_relief_threshold: false
flat_closing_radius_metres: 0
flat_opening_radius_metres: 0
min_component_pixels: 100
run_plane_analysis: false

write_combined_flat_mask_raster: true
write_metric_flat_mask_rasters: false
write_texture_rasters: false
"""

SURFACE_METRICS_PARAMS = """\
read_path: {read_path}
read_fname: fluvsurf_dem.tif
write_path: {write_path}
write_fname: sm

n_threads: 1
write_output_format: gtiff

write_otsu_raster: true
write_hough_lines_raster: true
write_vec_adaptive_threshold_rasters: true
write_vec_adaptive_skeleton_rasters: true
"""

CASES = [
    ("lsdtt-flat-surface-finder", FLAT_SURFACE_FINDER_PARAMS),
    ("lsdtt-surface-metrics", SURFACE_METRICS_PARAMS),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_case(exe: Path, template: str, read_path: Path, work: Path, env_extra: dict):
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    param = work / "run.param"
    param.write_text(template.format(read_path=str(read_path) + "/",
                                     write_path=str(work) + "/"))
    env = dict(os.environ)
    env.update(env_extra)
    proc = subprocess.run([str(exe), str(work) + "/", "run.param"],
                          cwd=work, env=env, capture_output=True, text=True)
    (work / "stdout.txt").write_text(proc.stdout)
    (work / "stderr.txt").write_text(proc.stderr)
    if proc.returncode != 0:
        sys.exit(f"{exe.name} failed ({proc.returncode}); see {work}/stderr.txt")
    # README and ingested-parameter files embed the run's write path, so they
    # differ between work directories and carry no image-ops information.
    skip = {"run.param", "stdout.txt", "stderr.txt"}
    skip_suffixes = ("_README.csv", "_README.html", "_ingested_params.param")
    return {p.name: sha256(p) for p in sorted(work.iterdir())
            if p.is_file() and p.name not in skip
            and not p.name.endswith(skip_suffixes)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bin-dir", required=True, type=Path)
    ap.add_argument("--fixture-dir", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rows = []
    for exe_name, template in CASES:
        exe = (args.bin_dir / exe_name).resolve()
        if not exe.exists():
            sys.exit(f"missing executable: {exe}")
        runs = {}
        for mode, env_extra in (("default", {}), ("reference", REFERENCE_ENV)):
            for rep in (1, 2):
                work = args.work_dir / f"{exe_name}_{mode}_{rep}"
                runs[(mode, rep)] = run_case(exe, template, args.fixture_dir.resolve(),
                                             work, env_extra)
        names = sorted(set().union(*[r.keys() for r in runs.values()]))
        for name in names:
            ref1 = runs[("reference", 1)].get(name, "")
            ref2 = runs[("reference", 2)].get(name, "")
            def1 = runs[("default", 1)].get(name, "")
            def2 = runs[("default", 2)].get(name, "")
            rows.append({
                "cli": exe_name,
                "file": name,
                "reference_sha256": ref1,
                "deterministic": "yes" if (ref1 == ref2 and def1 == def2) else "no",
                "default_matches_reference": "yes" if def1 == ref1 else "no",
            })
        print(f"{exe_name}: {len(names)} output files")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    n_nondet = sum(r["deterministic"] == "no" for r in rows)
    n_mode_diff = sum(r["default_matches_reference"] == "no" for r in rows)
    print(f"files: {len(rows)}; not byte-deterministic: {n_nondet}; "
          f"default differs from reference: {n_mode_diff}")


if __name__ == "__main__":
    main()
