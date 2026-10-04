#!/usr/bin/env python3
"""Run the two-stage fluvial-surfaces pipeline and check it hands over.

After phase 5 of the fluvial-surfaces split the valley is
identified by lsdtt-fluvial-surfaces and traced by lsdtt-valley-centreline, so
neither tool can be smoke-tested end to end on its own: CLI 3 will not start
without a mask, and the only supported source of that mask is CLI 2.

What this checks is the HANDOVER, which is the part a single-CLI test cannot
see. Both stages must run, and CLI 3 must actually consume what CLI 2 wrote -
so the run is required to say that identification was supplied, which is the
one line that distinguishes "read the mask" from "quietly derived its own".
"""

import argparse
import pathlib
import shutil
import subprocess
import sys
import lsdtt_test_env


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the fluvial pipeline")
    parser.add_argument("--surfaces-exe", required=True)
    parser.add_argument("--centreline-exe", required=True)
    parser.add_argument("--surfaces-param", required=True)
    parser.add_argument("--centreline-param", required=True)
    parser.add_argument("--read-path", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--expect", action="append", default=[],
                        help="File that must exist and be non-empty after stage 2")
    parser.add_argument("--stage2-stdout-contains", action="append", default=[],
                        help="Required substring in stage 2's stdout")
    parser.add_argument("--stage2-stderr-contains", action="append", default=[],
                        help="Required substring in stage 2's stderr")
    return parser.parse_args()


def run_stage(exe: pathlib.Path, work_dir: pathlib.Path, out_dir: pathlib.Path,
              template: pathlib.Path, read_path: pathlib.Path,
              tag: str) -> tuple[str, str]:
    rendered = (template.read_text(encoding="utf-8")
                .replace("__READ_PATH__", str(read_path))
                .replace("__WRITE_PATH__", str(out_dir)))
    param_name = f"{tag}.param"
    (work_dir / param_name).write_text(rendered, encoding="utf-8")

    result = lsdtt_test_env.run(
        [str(exe), str(work_dir), param_name],
        cwd=str(work_dir), capture_output=True, text=True, check=False)
    (work_dir / f"stdout_{tag}.log").write_text(result.stdout, encoding="utf-8")
    (work_dir / f"stderr_{tag}.log").write_text(result.stderr, encoding="utf-8")

    if result.returncode != 0:
        print(f"Validation failed: {tag} returned {result.returncode}",
              file=sys.stderr)
        print(f"Stdout: {work_dir / f'stdout_{tag}.log'}", file=sys.stderr)
        print(f"Stderr: {work_dir / f'stderr_{tag}.log'}", file=sys.stderr)
        sys.exit(result.returncode)
    return result.stdout, result.stderr


def main() -> int:
    args = parse_args()
    read_path = pathlib.Path(args.read_path).resolve()
    work_dir = pathlib.Path(args.work_dir).resolve()
    out_dir = work_dir / "out"

    if work_dir.exists():
        shutil.rmtree(work_dir)
    out_dir.mkdir(parents=True)

    # Stage 2 reads the mask stage 1 wrote, so its read path is the output
    # directory. The DEM and the source points are copied alongside it, which is
    # also what a user doing this by hand ends up with.
    run_stage(pathlib.Path(args.surfaces_exe).resolve(), work_dir, out_dir,
              pathlib.Path(args.surfaces_param).resolve(), read_path, "surfaces")

    for name in ("fluvsurf_dem.tif", "fluvsurf_sources.fgb"):
        src = read_path / name
        if src.exists():
            shutil.copy(src, out_dir / name)

    stage2_stdout, stage2_stderr = run_stage(
        pathlib.Path(args.centreline_exe).resolve(), work_dir,
        out_dir, pathlib.Path(args.centreline_param).resolve(),
        out_dir, "centreline")

    errors = 0
    for needle in args.stage2_stdout_contains:
        if needle not in stage2_stdout:
            print(f"Validation failed: stage 2 stdout lacks: {needle}",
                  file=sys.stderr)
            errors += 1
        else:
            print(f"  OK  stage 2 reported: {needle}")

    for needle in args.stage2_stderr_contains:
        if needle not in stage2_stderr:
            print(f"Validation failed: stage 2 stderr lacks: {needle}",
                  file=sys.stderr)
            errors += 1
        else:
            print(f"  OK  stage 2 warning reported: {needle}")

    for rel in args.expect:
        path = out_dir / rel
        if not path.exists():
            print(f"Validation failed: missing {rel}", file=sys.stderr)
            errors += 1
        elif path.stat().st_size == 0:
            print(f"Validation failed: {rel} is empty", file=sys.stderr)
            errors += 1
        else:
            print(f"  OK  {rel} ({path.stat().st_size} bytes)")

    if errors:
        print(f"FAIL: {errors} check(s) failed", file=sys.stderr)
        return 1
    print("PASS: the pipeline runs and stage 2 consumed stage 1's valley mask")
    return 0


if __name__ == "__main__":
    sys.exit(main())
