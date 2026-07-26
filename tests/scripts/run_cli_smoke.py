#!/usr/bin/env python3

import argparse
import hashlib
import pathlib
import re
import shutil
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one lsdtt CLI smoke test")
    parser.add_argument("--exe", required=True, help="Path to CLI executable")
    parser.add_argument("--param-template", required=True, help="Template parameter file")
    parser.add_argument("--read-path", required=True, help="Input data directory")
    parser.add_argument("--work-dir", required=True, help="Working directory for this test")
    parser.add_argument("--expect", action="append", required=True,
                        help="Expected output file path relative to the test output directory")
    parser.add_argument(
        "--expect-failure",
        action="store_true",
        help="Expect the CLI to fail (non-zero exit code)",
    )
    parser.add_argument(
        "--stderr-contains",
        action="append",
        default=[],
        help="Required substring in stderr.log",
    )
    parser.add_argument(
        "--stdout-contains",
        action="append",
        default=[],
        help="Required substring in stdout.log (used to pin run-header reporting)",
    )
    parser.add_argument(
        "--expect-raster-min",
        action="append",
        default=[],
        help=(
            "Expected minimum band-1 value as rel_path:MIN; the raster must also "
            "be 100%% valid, so a NaN or nodata pixel fails. "
            "For example out_PROD.tif:0.0"
        ),
    )
    parser.add_argument(
        "--expect-valid-percent-100",
        action="append",
        default=[],
        help=(
            "Expected raster output (relative path) that must report "
            "STATISTICS_VALID_PERCENT=100 from gdalinfo -stats"
        ),
    )
    parser.add_argument(
        "--expect-valid-percent",
        action="append",
        default=[],
        help=(
            "Expected raster valid percent as rel_path:PERCENT, "
            "for example tiny_nodata_filled.tif:96"
        ),
    )
    parser.add_argument(
        "--expect-raster-size",
        action="append",
        default=[],
        help=(
            "Expected raster size as rel_path:WIDTHxHEIGHT, "
            "for example tiny_fill.tif:5x5"
        ),
    )
    parser.add_argument(
        "--expect-vector-feature-min",
        action="append",
        default=[],
        help=(
            "Expected minimum feature count as rel_path:MIN, "
            "for example tiny_channel_network_lines.fgb:1"
        ),
    )
    parser.add_argument(
        "--expect-raster-checksum",
        action="append",
        default=[],
        help=(
            "Expected band-1 GDAL checksum as rel_path:CHECKSUM, "
            "for example tiny_difference.tif:12345"
        ),
    )
    parser.add_argument(
        "--expect-sha256",
        action="append",
        default=[],
        help=(
            "Expected SHA256 as rel_path:HEX, for example out.csv:abcd..."
        ),
    )
    return parser.parse_args()


def assert_valid_percent_100(raster_path: pathlib.Path) -> None:
    cmd = ["gdalinfo", "-stats", str(raster_path)]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "gdalinfo failed for "
            f"{raster_path} (exit {result.returncode}). "
            "Ensure GDAL command-line tools are installed."
        )

    if "STATISTICS_VALID_PERCENT=100" not in result.stdout:
        raise AssertionError(
            f"Expected STATISTICS_VALID_PERCENT=100 in {raster_path}, "
            "but it was not found in gdalinfo output."
        )


def assert_raster_min(raster_path: pathlib.Path, minimum: float) -> None:
    """Band-1 minimum must be >= minimum, and every pixel must be valid.

    gdalinfo -stats computes its statistics over valid pixels only, so the
    valid-percent check is what makes a NaN or nodata pixel fail rather than
    quietly drop out of the minimum.
    """
    cmd = ["gdalinfo", "-stats", str(raster_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"gdalinfo -stats failed for {raster_path} (exit {result.returncode})."
        )

    match = re.search(r"STATISTICS_MINIMUM=([-0-9.eE+]+)", result.stdout)
    if not match:
        raise AssertionError(f"Could not parse STATISTICS_MINIMUM for {raster_path}.")
    got = float(match.group(1))
    if not got >= minimum:
        raise AssertionError(
            f"Minimum too low for {raster_path}: expected >= {minimum}, got {got}."
        )
    if "STATISTICS_VALID_PERCENT=100" not in result.stdout:
        raise AssertionError(
            f"Expected every pixel of {raster_path} to be valid and finite "
            "(STATISTICS_VALID_PERCENT=100)."
        )


def read_valid_percent(raster_path: pathlib.Path) -> str:
    cmd = ["gdalinfo", "-stats", str(raster_path)]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"gdalinfo -stats failed for {raster_path} (exit {result.returncode})."
        )

    match = re.search(r"STATISTICS_VALID_PERCENT=([0-9]+(?:\.[0-9]+)?)", result.stdout)
    if not match:
        raise AssertionError(
            f"Could not parse STATISTICS_VALID_PERCENT for {raster_path}."
        )
    return match.group(1)


def assert_raster_size(raster_path: pathlib.Path, width: int, height: int) -> None:
    cmd = ["gdalinfo", str(raster_path)]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"gdalinfo failed for {raster_path} (exit {result.returncode})."
        )

    match = re.search(r"Size is\s+(\d+),\s*(\d+)", result.stdout)
    if not match:
        raise AssertionError(f"Could not parse raster size for {raster_path}.")

    got_w = int(match.group(1))
    got_h = int(match.group(2))
    if got_w != width or got_h != height:
        raise AssertionError(
            f"Raster size mismatch for {raster_path}: "
            f"expected {width}x{height}, got {got_w}x{got_h}."
        )


def assert_vector_feature_min(vector_path: pathlib.Path, minimum: int) -> None:
    cmd = ["ogrinfo", "-so", "-al", str(vector_path)]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"ogrinfo failed for {vector_path} (exit {result.returncode})."
        )

    match = re.search(r"Feature Count:\s*(\d+)", result.stdout)
    if not match:
        raise AssertionError(
            f"Could not parse feature count for {vector_path}."
        )

    count = int(match.group(1))
    if count < minimum:
        raise AssertionError(
            f"Feature count too low for {vector_path}: "
            f"expected >= {minimum}, got {count}."
        )


def assert_raster_checksum(raster_path: pathlib.Path, expected: int) -> None:
    cmd = ["gdalinfo", "-checksum", str(raster_path)]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"gdalinfo -checksum failed for {raster_path} (exit {result.returncode})."
        )

    match = re.search(r"Checksum=([0-9]+)", result.stdout)
    if not match:
        raise AssertionError(f"Could not parse checksum for {raster_path}.")

    got = int(match.group(1))
    if got != expected:
        raise AssertionError(
            f"Checksum mismatch for {raster_path}: expected {expected}, got {got}."
        )


def assert_sha256(file_path: pathlib.Path, expected_hex: str) -> None:
    h = hashlib.sha256()
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    got = h.hexdigest()
    if got != expected_hex.lower():
        raise AssertionError(
            f"SHA256 mismatch for {file_path}: expected {expected_hex}, got {got}."
        )


def main() -> int:
    args = parse_args()

    exe = pathlib.Path(args.exe).resolve()
    template_path = pathlib.Path(args.param_template).resolve()
    read_path = pathlib.Path(args.read_path).resolve()
    work_dir = pathlib.Path(args.work_dir).resolve()
    out_dir = work_dir / "out"

    if not exe.exists():
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not template_path.exists():
        raise FileNotFoundError(f"Template not found: {template_path}")
    if not read_path.exists():
        raise FileNotFoundError(f"Read path not found: {read_path}")

    if work_dir.exists():
        shutil.rmtree(work_dir)
    out_dir.mkdir(parents=True)

    template_text = template_path.read_text(encoding="utf-8")
    rendered = (template_text
                .replace("__READ_PATH__", str(read_path))
                .replace("__WRITE_PATH__", str(out_dir)))

    param_name = "run_params.txt"
    param_file = work_dir / param_name
    param_file.write_text(rendered, encoding="utf-8")

    cmd = [str(exe), str(work_dir), param_name]
    result = subprocess.run(
        cmd,
        cwd=str(work_dir),
        capture_output=True,
        text=True,
        check=False,
    )

    (work_dir / "stdout.log").write_text(result.stdout, encoding="utf-8")
    (work_dir / "stderr.log").write_text(result.stderr, encoding="utf-8")

    if args.expect_failure:
        if result.returncode == 0:
            print("Validation failed: CLI was expected to fail but returned 0", file=sys.stderr)
            print(f"Command: {' '.join(cmd)}", file=sys.stderr)
            print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
            print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
            return 4
    else:
        if result.returncode != 0:
            print("Smoke test failed: CLI returned non-zero exit code", file=sys.stderr)
            print(f"Command: {' '.join(cmd)}", file=sys.stderr)
            print(f"Exit code: {result.returncode}", file=sys.stderr)
            print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
            print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
            return result.returncode

    for needle in args.stdout_contains:
        if needle not in result.stdout:
            print(
                f"Validation failed: stdout does not contain required text: {needle}",
                file=sys.stderr,
            )
            print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
            print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
            return 5

    stderr_text = result.stderr
    for needle in args.stderr_contains:
        if needle not in stderr_text:
            print(
                f"Validation failed: stderr does not contain required text: {needle}",
                file=sys.stderr,
            )
            print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
            print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
            return 5

    if args.expect_failure:
        print(f"PASS (expected failure): {exe.name}")
        return 0

    missing = []
    empty = []
    for expected in args.expect:
        p = out_dir / expected
        if not p.exists():
            missing.append(str(p))
            continue
        if p.stat().st_size == 0:
            empty.append(str(p))

    if missing or empty:
        if missing:
            print("Missing outputs:", file=sys.stderr)
            for m in missing:
                print(f"  {m}", file=sys.stderr)
        if empty:
            print("Empty outputs:", file=sys.stderr)
            for e in empty:
                print(f"  {e}", file=sys.stderr)
        print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
        print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
        return 2

    try:
        for rel_path in args.expect_valid_percent_100:
            assert_valid_percent_100(out_dir / rel_path)

        for spec in args.expect_valid_percent:
            rel_path, pct = spec.rsplit(":", 1)
            got = read_valid_percent(out_dir / rel_path)
            if got != pct:
                raise AssertionError(
                    f"Valid percent mismatch for {rel_path}: expected {pct}, got {got}."
                )

        for spec in args.expect_raster_min:
            rel_path, min_spec = spec.rsplit(":", 1)
            assert_raster_min(out_dir / rel_path, float(min_spec))

        for spec in args.expect_raster_size:
            rel_path, size_spec = spec.rsplit(":", 1)
            width_str, height_str = size_spec.lower().split("x", 1)
            assert_raster_size(out_dir / rel_path, int(width_str), int(height_str))

        for spec in args.expect_vector_feature_min:
            rel_path, min_spec = spec.rsplit(":", 1)
            assert_vector_feature_min(out_dir / rel_path, int(min_spec))

        for spec in args.expect_raster_checksum:
            rel_path, checksum_spec = spec.rsplit(":", 1)
            assert_raster_checksum(out_dir / rel_path, int(checksum_spec))

        for spec in args.expect_sha256:
            rel_path, hash_hex = spec.rsplit(":", 1)
            assert_sha256(out_dir / rel_path, hash_hex)
    except Exception as exc:
        print(f"Validation check failed: {exc}", file=sys.stderr)
        print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
        print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
        return 3

    print(f"PASS: {exe.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
