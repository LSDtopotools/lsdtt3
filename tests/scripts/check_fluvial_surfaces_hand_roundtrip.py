#!/usr/bin/env python3
"""Pin the hand_fname composability guarantee of lsdtt-fluvial-surfaces.

HAND is the expensive stage of lsdtt-fluvial-surfaces and the level count is
what a user iterates on, so the tool lets you compute HAND once, write it, and
feed it back with hand_fname to try different ladder settings for the cost of
the ladder alone. That is only worth anything if the two routes agree exactly -
a supplied HAND that gave slightly different bands would make every iteration a
different experiment from the final run.

So this runs the CLI twice on the same DEM:

  arm 1  computes HAND and writes it (write_hand_raster)
  arm 2  reads that raster back as hand_fname

and asserts the band index and the banded labels are IDENTICAL by GDAL
checksum, not merely similar.

Both arms use a global ladder (ladder_reach_window_metres < 0), and that is
load-bearing rather than incidental: reach-relative windowing needs each
pixel's HAND donor node, which a raster cannot carry, so a supplied HAND is
necessarily banded globally. Comparing a windowed computed run against a
supplied one would compare two different ladders and fail for a reason that has
nothing to do with the guarantee.
"""

import argparse
import pathlib
import re
import shutil
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify hand_fname reproduces the computed HAND path")
    parser.add_argument("--exe", required=True)
    parser.add_argument("--param-template", required=True)
    parser.add_argument("--read-path", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--stem", required=True,
                        help="write_fname stem used by the parameter template")
    return parser.parse_args()


def raster_checksum(path: pathlib.Path) -> int:
    result = subprocess.run(
        ["gdalinfo", "-checksum", str(path)],
        capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"gdalinfo failed for {path}: {result.stderr}")
    match = re.search(r"Checksum=(-?\d+)", result.stdout)
    if match is None:
        raise RuntimeError(f"no checksum reported for {path}")
    return int(match.group(1))


def run_cli(exe: pathlib.Path, work_dir: pathlib.Path, out_dir: pathlib.Path,
            param_text: str, tag: str) -> None:
    param_name = f"run_{tag}.txt"
    (work_dir / param_name).write_text(param_text, encoding="utf-8")
    result = subprocess.run(
        [str(exe), str(work_dir), param_name],
        cwd=str(work_dir), capture_output=True, text=True, check=False)
    (work_dir / f"stdout_{tag}.log").write_text(result.stdout, encoding="utf-8")
    (work_dir / f"stderr_{tag}.log").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        print(f"Validation failed: {tag} arm returned {result.returncode}",
              file=sys.stderr)
        print(f"Stdout: {work_dir / f'stdout_{tag}.log'}", file=sys.stderr)
        print(f"Stderr: {work_dir / f'stderr_{tag}.log'}", file=sys.stderr)
        sys.exit(result.returncode)
    _ = out_dir


def main() -> int:
    args = parse_args()
    exe = pathlib.Path(args.exe).resolve()
    template = pathlib.Path(args.param_template).resolve()
    read_path = pathlib.Path(args.read_path).resolve()
    work_dir = pathlib.Path(args.work_dir).resolve()
    out_dir = work_dir / "out"
    stem = args.stem

    if work_dir.exists():
        shutil.rmtree(work_dir)
    out_dir.mkdir(parents=True)

    base = (template.read_text(encoding="utf-8")
            .replace("__READ_PATH__", str(read_path))
            .replace("__WRITE_PATH__", str(out_dir)))

    # Arm 1: compute HAND and write it.
    run_cli(exe, work_dir, out_dir, base, "computed")

    hand_path = out_dir / f"{stem}_hand.tif"
    if not hand_path.exists():
        print(f"Validation failed: {hand_path.name} was not written, so the "
              "second arm has nothing to read", file=sys.stderr)
        return 1

    # Arm 2: read that raster back. read_path becomes the output directory so
    # hand_fname resolves against it, and the DEM is copied alongside.
    dem_name = None
    for line in base.splitlines():
        if line.strip().startswith("read_fname:"):
            dem_name = line.split(":", 1)[1].strip()
    if dem_name is None:
        print("Validation failed: the template has no read_fname",
              file=sys.stderr)
        return 1
    shutil.copy(read_path / dem_name, out_dir / dem_name)

    supplied = []
    for line in base.splitlines():
        if line.strip().startswith("read_path:"):
            supplied.append(f"read_path: {out_dir}")
        elif line.strip().startswith("write_fname:"):
            supplied.append(f"write_fname: {stem}_supplied")
        else:
            supplied.append(line)
    supplied.append(f"hand_fname: {stem}_hand.tif")
    run_cli(exe, work_dir, out_dir, "\n".join(supplied) + "\n", "supplied")

    stdout_supplied = (work_dir / "stdout_supplied.log").read_text(
        encoding="utf-8")
    if "Using the supplied HAND raster" not in stdout_supplied:
        print("Validation failed: the second arm did not report using the "
              "supplied HAND raster, so it may have recomputed it",
              file=sys.stderr)
        return 1

    errors = 0
    for suffix in ("_band_index.tif", "_banded_labels.tif"):
        a = out_dir / f"{stem}{suffix}"
        b = out_dir / f"{stem}_supplied{suffix}"
        for path in (a, b):
            if not path.exists():
                print(f"Validation failed: missing {path.name}",
                      file=sys.stderr)
                errors += 1
        if errors:
            continue
        ca, cb = raster_checksum(a), raster_checksum(b)
        if ca != cb:
            print(f"Validation failed: {suffix} differs between the computed "
                  f"and supplied HAND paths (checksum {ca} vs {cb})",
                  file=sys.stderr)
            errors += 1
        else:
            print(f"  OK  {suffix}: identical, checksum {ca}")

    if errors:
        print(f"FAIL: {errors} check(s) failed", file=sys.stderr)
        return 1
    print("PASS: hand_fname reproduces the computed HAND path exactly")
    return 0


if __name__ == "__main__":
    sys.exit(main())
