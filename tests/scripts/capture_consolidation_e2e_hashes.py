#!/usr/bin/env python3
# capture_consolidation_e2e_hashes.py
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# -----------------------------------------------------------------------------
# End-to-end numbers gate for refactors that must not change any output (the
# src/utils + src/cpu duplicate-code consolidation is the first user).
#
# Most CLI tests check that outputs exist, have the right shape, or sit inside a
# tolerance, so a refactor that moves results in the last bits passes them. This
# script runs the CLI-level CTest tiers and hashes EVERY file those runs write,
# so any change in any output is seen:
#
#   rasters   decoded pixel data per band, plus data type, size, nodata value
#             and geotransform. File bytes are NOT hashed, so GDAL metadata
#             differences cannot raise false alarms.
#   vectors   every feature in file order: field values and geometry as WKB.
#   CSV       file bytes, except that a wall-clock column (WALL_CLOCK_CSV_COLUMNS,
#             e.g. "seconds" in *_snap_convergence.csv) is removed first, so the
#             numbers beside it still gate.
#   other     file bytes (text, JSON, ...).
#
# Two commands:
#
#   capture   run the selected tests N times (default 2) with the default
#             thread count and once with OMP_NUM_THREADS=1, and write the
#             manifest. A file whose hash differs between the repeated default
#             runs is marked deterministic=no and is not used as a gate.
#   check     run once per thread setting and compare against the manifest.
#             Exit status 1 if any deterministic file changed, vanished or
#             appeared.
#
#   python3 tests/scripts/capture_consolidation_e2e_hashes.py capture \
#       --build-dir build --manifest build/consolidation_e2e/manifest.csv
#   python3 tests/scripts/capture_consolidation_e2e_hashes.py check \
#       --build-dir build --manifest build/consolidation_e2e/manifest.csv
#
# The manifest belongs in the build directory, not the repository: it depends
# on the GDAL and compiler on the machine that made it, so it is a gate for a
# piece of work, not a CTest (the same reasoning as
# capture_image_ops_e2e_hashes.py).
#
# Files are found by modification time: everything under --build-dir written
# after the run started, outside CMakeFiles/ and Testing/. Nothing may be built
# while the script runs.
# -----------------------------------------------------------------------------

import argparse
import csv
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path

import lsdtt_test_env

DEFAULT_LABELS = "^cli$"
DEFAULT_EXCLUDE = "^(heavy|needs_network)$"

RASTER_SUFFIXES = {".tif", ".tiff", ".bil", ".flt", ".asc", ".img", ".nc", ".vrt"}
VECTOR_SUFFIXES = {".shp", ".gpkg", ".geojson", ".fgb", ".gml", ".kml"}
# Sidecars are covered by the file they belong to.
SIDECAR_SUFFIXES = {".hdr", ".aux", ".xml", ".prj", ".dbf", ".shx", ".cpg", ".qix",
                    ".ovr", ".msk", ".pyc"}
# Columns holding elapsed wall-clock time. Matched exactly against the header;
# they are dropped before a CSV is hashed and nothing else in the file is.
WALL_CLOCK_CSV_COLUMNS = {"seconds"}
SKIP_DIRS = {"CMakeFiles", "Testing", "__pycache__", ".pytest_cache"}

THREAD_PASSES = {
    "default": {},
    "omp1": {"OMP_NUM_THREADS": "1"},
}


def hash_raster(path: Path):
    from osgeo import gdal
    ds = gdal.Open(str(path))
    if ds is None:
        return None
    h = hashlib.sha256()
    h.update(f"{ds.RasterXSize}x{ds.RasterYSize}x{ds.RasterCount}".encode())
    h.update(repr(tuple(ds.GetGeoTransform(can_return_null=True) or ())).encode())
    for b in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(b)
        h.update(f"|{gdal.GetDataTypeName(band.DataType)}|{band.GetNoDataValue()!r}|".encode())
        h.update(band.ReadRaster())
    return h.hexdigest()


def hash_vector(path: Path):
    from osgeo import ogr
    ds = ogr.Open(str(path))
    if ds is None:
        return None
    h = hashlib.sha256()
    for li in range(ds.GetLayerCount()):
        layer = ds.GetLayerByIndex(li)
        defn = layer.GetLayerDefn()
        names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
        h.update(("layer|" + layer.GetName() + "|" + ",".join(names)).encode())
        layer.ResetReading()
        for feat in layer:
            for i in range(len(names)):
                h.update(b"\x1f")
                if feat.IsFieldSetAndNotNull(i):
                    h.update(str(feat.GetField(i)).encode())
            geom = feat.GetGeometryRef()
            h.update(b"\x1e" + (bytes(geom.ExportToWkb()) if geom is not None else b""))
    return h.hexdigest()


def hash_bytes(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_csv(path: Path):
    with path.open("r", newline="", errors="surrogateescape") as f:
        rows = list(csv.reader(f))
    if not rows:
        return hash_bytes(path)
    drop = {i for i, name in enumerate(rows[0]) if name.strip() in WALL_CLOCK_CSV_COLUMNS}
    if not drop:
        return hash_bytes(path)
    h = hashlib.sha256()
    for row in rows:
        kept = [v for i, v in enumerate(row) if i not in drop]
        h.update(("\x1f".join(kept) + "\x1e").encode("utf-8", "surrogateescape"))
    return h.hexdigest()


def hash_file(path: Path):
    suffix = path.suffix.lower()
    if suffix in RASTER_SUFFIXES:
        d = hash_raster(path)
        if d is not None:
            return "raster", d
    if suffix in VECTOR_SUFFIXES:
        d = hash_vector(path)
        if d is not None:
            return "vector", d
    if suffix == ".csv":
        return "csv", hash_csv(path)
    return "bytes", hash_bytes(path)


def written_since(build_dir: Path, t0: float):
    for dirpath, dirnames, filenames in os.walk(build_dir):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for name in filenames:
            p = Path(dirpath) / name
            if p.suffix.lower() in SIDECAR_SUFFIXES:
                continue
            try:
                if p.stat().st_mtime >= t0 and p.is_file():
                    yield p
            except OSError:
                continue


def run_suite(build_dir: Path, labels: str, exclude: str, jobs: int, env_extra: dict,
              log: Path, own_dir: Path):
    """Run the selected tests; return {relative path: (kind, digest)} and failures."""
    env = dict(os.environ)
    env.update(env_extra)
    # Drop the sub-second part so a file written in the same second as the start
    # is still counted on file systems with coarse timestamps.
    t0 = float(int(time.time()))
    time.sleep(1.1)
    cmd = ["ctest", "--test-dir", str(build_dir), "-L", labels, "-LE", exclude,
           "-j", str(jobs), "--output-on-failure", "--timeout", "1800"]
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as f:
        proc = lsdtt_test_env.run(cmd, env=env, stdout=f, stderr=subprocess.STDOUT)
    failed = parse_failed(log)
    outputs = {}
    for p in written_since(build_dir, t0):
        if p.is_relative_to(own_dir):   # this script's manifest and logs
            continue
        outputs[p.relative_to(build_dir).as_posix()] = hash_file(p)
    return outputs, failed, proc.returncode


def parse_failed(log: Path):
    failed = []
    in_list = False
    for line in log.read_text(errors="ignore").splitlines():
        if line.startswith("The following tests FAILED"):
            in_list = True
            continue
        if in_list:
            parts = line.split()
            if len(parts) >= 3 and parts[1] == "-":
                failed.append(parts[2])
            elif line.strip() and not line.startswith("\t"):
                in_list = False
    return failed


def common_args(ap):
    ap.add_argument("--build-dir", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--labels", default=DEFAULT_LABELS,
                    help=f"ctest -L regex (default: {DEFAULT_LABELS})")
    ap.add_argument("--exclude-labels", default=DEFAULT_EXCLUDE,
                    help=f"ctest -LE regex (default: {DEFAULT_EXCLUDE})")
    ap.add_argument("-j", "--jobs", type=int, default=4)
    ap.add_argument("--passes", default="default,omp1",
                    help="comma list of thread passes: default, omp1")


def cmd_capture(args):
    build = args.build_dir.resolve()
    logs = args.manifest.resolve().parent / "logs"
    passes = args.passes.split(",")
    runs = {}
    failures = {}
    for tp in passes:
        reps = args.repeats if tp == "default" else 1
        for rep in range(1, reps + 1):
            key = (tp, rep)
            print(f"[capture] pass {tp} repeat {rep} ...", flush=True)
            t = time.time()
            runs[key], failures[key], rc = run_suite(
                build, args.labels, args.exclude_labels, args.jobs,
                THREAD_PASSES[tp], logs / f"capture_{tp}_{rep}.log", logs.parent)
            print(f"[capture]   {len(runs[key])} files, {len(failures[key])} failed tests, "
                  f"{time.time() - t:.0f} s", flush=True)

    rows = []
    for tp in passes:
        reps = [runs[(tp, r)] for r in range(1, (args.repeats if tp == "default" else 1) + 1)]
        names = sorted(set().union(*[r.keys() for r in reps]))
        for name in names:
            got = [r.get(name) for r in reps]
            first = got[0] or ("", "")
            det = all(g == got[0] for g in got) if len(got) > 1 else None
            rows.append({
                "pass": tp,
                "file": name,
                "kind": first[0],
                "sha256": first[1],
                "deterministic": "unknown" if det is None else ("yes" if det else "no"),
            })
    # Determinism measured on the default pass carries over to the other passes.
    det_default = {r["file"]: r["deterministic"] for r in rows if r["pass"] == "default"}
    for r in rows:
        if r["deterministic"] == "unknown":
            r["deterministic"] = det_default.get(r["file"], "unknown")

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with (args.manifest.parent / "capture_failed_tests.txt").open("w") as f:
        for key, names in failures.items():
            for n in names:
                f.write(f"{key[0]}\t{key[1]}\t{n}\n")

    n_nondet = sum(r["deterministic"] == "no" for r in rows if r["pass"] == "default")
    n_default = sum(r["pass"] == "default" for r in rows)
    print(f"[capture] manifest: {args.manifest}")
    print(f"[capture] files (default pass): {n_default}; not deterministic: {n_nondet}")
    for key, names in failures.items():
        if names:
            print(f"[capture] failed tests in {key[0]}/{key[1]}: {' '.join(names)}")
    if "default" in passes and "omp1" in passes:
        d = {r["file"]: r["sha256"] for r in rows if r["pass"] == "default"}
        o = {r["file"]: r["sha256"] for r in rows if r["pass"] == "omp1"}
        n_thread_diff = sum(1 for k in d.keys() & o.keys()
                            if d[k] != o[k] and det_default.get(k) == "yes")
        print(f"[capture] deterministic files whose output differs between "
              f"default threads and OMP_NUM_THREADS=1: {n_thread_diff}")


def cmd_check(args):
    build = args.build_dir.resolve()
    logs = args.manifest.resolve().parent / "logs"
    want = {}
    for r in csv.DictReader(args.manifest.open()):
        want[(r["pass"], r["file"])] = r
    passes = sorted({p for (p, _) in want} & set(args.passes.split(",")))
    baseline_failed = set()
    failed_file = args.manifest.parent / "capture_failed_tests.txt"
    if failed_file.exists():
        for line in failed_file.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) == 3:
                baseline_failed.add((parts[0], parts[2]))
    n_same = n_bad = n_ignored = 0
    report = []
    for tp in passes:
        print(f"[check] pass {tp} ...", flush=True)
        got, failed, rc = run_suite(build, args.labels, args.exclude_labels, args.jobs,
                                    THREAD_PASSES[tp], logs / f"check_{tp}.log", logs.parent)
        names = sorted({f for (p, f) in want if p == tp} | set(got))
        for name in names:
            w = want.get((tp, name))
            g = got.get(name)
            if w is not None and w["deterministic"] != "yes":
                n_ignored += 1
                continue
            if w is None:
                status = "NEW"
            elif g is None:
                status = "MISSING"
            elif g[1] != w["sha256"]:
                status = "DIFFERS"
            else:
                status = "same"
            if status == "same":
                n_same += 1
            else:
                n_bad += 1
                report.append(f"{status:8s} {tp:8s} {name}")
        new_failed = [t for t in failed if (tp, t) not in baseline_failed]
        old_failed = [t for t in failed if (tp, t) in baseline_failed]
        if new_failed:
            n_bad += len(new_failed)
            report.append(f"NEWLY-FAILING {tp}: {' '.join(new_failed)}")
        if old_failed:
            report.append(f"(already failing at capture) {tp}: {' '.join(old_failed)}")
    for line in report:
        print(line)
    print(f"[check] identical: {n_same}; changed/missing/new: {n_bad}; "
          f"ignored (not deterministic): {n_ignored}")
    sys.exit(0 if n_bad == 0 else 1)


def main():
    ap = argparse.ArgumentParser(
        description="End-to-end output hashes for no-numbers-change refactors.")
    sub = ap.add_subparsers(dest="command", required=True)
    cap = sub.add_parser("capture", help="run the tests and write the manifest")
    common_args(cap)
    cap.add_argument("--repeats", type=int, default=2,
                     help="default-thread repeats used to detect nondeterminism")
    chk = sub.add_parser("check", help="run the tests and compare with the manifest")
    common_args(chk)
    args = ap.parse_args()
    if args.command == "capture":
        cmd_capture(args)
    else:
        cmd_check(args)


if __name__ == "__main__":
    main()
