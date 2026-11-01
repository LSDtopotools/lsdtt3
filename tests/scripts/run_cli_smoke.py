#!/usr/bin/env python3

import argparse
import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import lsdtt_test_env


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one lsdtt CLI smoke test")
    parser.add_argument("--exe", required=True, help="Path to CLI executable")
    parser.add_argument(
        "--cli-src-dir",
        default=None,
        help=(
            "The repository's cli/ directory. When given, the executable is "
            "checked against cli/<name>.cpp and the test FAILS if the binary is "
            "older than its source. See assert_binary_is_current()."
        ),
    )
    parser.add_argument("--param-template", required=True, help="Template parameter file")
    parser.add_argument("--read-path", required=True, help="Input data directory")
    parser.add_argument("--work-dir", required=True, help="Working directory for this test")
    parser.add_argument("--expect", action="append", required=True,
                        help="Expected output file path relative to the test output directory")
    parser.add_argument(
        "--expect-absent",
        action="append",
        default=[],
        help="Output file path, relative to the test output directory, that must "
             "NOT exist (e.g. an output whose switch is off)",
    )
    parser.add_argument(
        "--expect-failure",
        action="store_true",
        help="Expect the CLI to fail (non-zero exit code)",
    )
    parser.add_argument(
        "--expect-exit-code",
        type=int,
        default=None,
        help="Expect this exact exit code. A non-zero value implies --expect-failure",
    )
    parser.add_argument(
        "--stderr-contains",
        action="append",
        default=[],
        help="Required substring in stderr.log",
    )
    parser.add_argument(
        "--stderr-not-contains",
        action="append",
        default=[],
        help="Substring that must not appear in stderr.log",
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
        "--expect-raster-stat",
        action="append",
        default=[],
        help=(
            "Expected band-1 statistic within a band, as rel_path:STAT:LO:HI "
            "with STAT one of min, max, mean. Two-sided, so it pins a value "
            "rather than a direction. For example "
            "lemstrip_modelled_topography.tif:max:4.90:4.99"
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
        "--expect-raster-values",
        action="append",
        default=[],
        help=(
            "Expected exact set of non-negative band-1 values as "
            "rel_path:V1,V2,...; nodata and negative sentinels are ignored"
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
    parser.add_argument(
        "--expect-file-contains",
        action="append",
        default=[],
        help=(
            "Substring an output file must contain, as rel_path:SUBSTRING, "
            "for example out.csv:self_consistent_t_yr. Splits on the FIRST "
            "colon, so the substring may itself contain colons."
        ),
    )
    parser.add_argument(
        "--expect-file-not-contains",
        action="append",
        default=[],
        help=(
            "Substring an output file must NOT contain, as rel_path:SUBSTRING. "
            "The negative form is what catches a nodata sentinel or a NaN "
            "reaching a data column: a positive assertion cannot, because the "
            "value is still present and still numeric. Splits on the FIRST "
            "colon."
        ),
    )
    return parser.parse_args()


def assert_binary_is_current(exe: pathlib.Path,
                             cli_src_dir: pathlib.Path | None) -> None:
    """Fail if the executable is older than the CLI source it was built from.

    WHY THIS IS A TEST FAILURE AND NOT A WARNING.

    A build system will happily run yesterday's binary. That makes "the test
    passed" and "the code passed" different statements, and the difference is
    invisible in the result: a stale binary produces a confident green.

    Both halves of that bit this project on 2026-09-10. The parameter-convention
    baseline was generated against a tree with ten stale binaries and one
    missing, so it recorded two violations the source had already fixed and
    missed nine that were really there. And smoke_feature_finder passed against
    a lsdtt-feature-finder built before a signature change, while a full build
    was at that moment FAILING on two call sites in that very file.

    So the check is fatal, and it names the rebuild. A source with no binary at
    all is caught by the existing FileNotFoundError on `exe`.
    """
    if cli_src_dir is None:
        return

    source = cli_src_dir / (exe.name + ".cpp")
    if not source.is_file():
        # Not a cli/ program (a test driver, say). Nothing to compare against.
        return

    if exe.stat().st_mtime < source.stat().st_mtime:
        raise AssertionError(
            f"{exe.name} is STALE: the binary predates {source.name}, so this "
            f"test would report on code that is no longer in the tree.\n"
            f"  binary: {exe}\n"
            f"  source: {source}\n"
            f"Rebuild and re-run:  cmake --build <build dir> --target {exe.name}"
        )


def assert_file_contains(file_path: pathlib.Path, needle: str) -> None:
    if not file_path.is_file():
        raise AssertionError(f"Expected file {file_path} does not exist")
    text = file_path.read_text(encoding="utf-8", errors="replace")
    if needle not in text:
        raise AssertionError(
            f"Expected {file_path} to contain {needle!r}, but it does not. "
            f"First line was: {text.splitlines()[0] if text else '(empty)'!r}"
        )


def assert_file_not_contains(file_path: pathlib.Path, needle: str) -> None:
    if not file_path.is_file():
        raise AssertionError(f"Expected file {file_path} does not exist")
    text = file_path.read_text(encoding="utf-8", errors="replace")
    if needle in text:
        offending = next(
            (ln for ln in text.splitlines() if needle in ln), ""
        )
        raise AssertionError(
            f"Expected {file_path} NOT to contain {needle!r}, but it does. "
            f"Offending line: {offending[:200]!r}"
        )


def assert_valid_percent_100(raster_path: pathlib.Path) -> None:
    cmd = ["gdalinfo", "-stats", str(raster_path)]
    result = lsdtt_test_env.run(
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
    result = lsdtt_test_env.run(cmd, capture_output=True, text=True, check=False)
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


def assert_raster_stat(
    raster_path: pathlib.Path, stat: str, low: float, high: float
) -> None:
    """Band-1 min / max / mean must lie in [low, high].

    Two-sided on purpose: a physical steady state is pinned by a band around a
    known value, and a one-sided bound would be satisfied by a run that simply
    did not evolve.
    """
    keys = {
        "min": "STATISTICS_MINIMUM",
        "max": "STATISTICS_MAXIMUM",
        "mean": "STATISTICS_MEAN",
    }
    if stat not in keys:
        raise AssertionError(f"Unknown raster stat {stat!r}; use min, max or mean.")

    cmd = ["gdalinfo", "-stats", str(raster_path)]
    result = lsdtt_test_env.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"gdalinfo -stats failed for {raster_path} (exit {result.returncode})."
        )

    match = re.search(rf"{keys[stat]}=([-0-9.eE+]+)", result.stdout)
    if not match:
        raise AssertionError(f"Could not parse {keys[stat]} for {raster_path}.")

    got = float(match.group(1))
    if not (low <= got <= high):
        raise AssertionError(
            f"{keys[stat]} out of band for {raster_path}: expected "
            f"[{low}, {high}], got {got}."
        )


def read_valid_percent(raster_path: pathlib.Path) -> str:
    cmd = ["gdalinfo", "-stats", str(raster_path)]
    result = lsdtt_test_env.run(
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
    result = lsdtt_test_env.run(
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
    result = lsdtt_test_env.run(
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
    result = lsdtt_test_env.run(
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


def assert_raster_values(raster_path: pathlib.Path, expected: set[float]) -> None:
    cmd = ["gdal_translate", "-q", "-of", "XYZ", str(raster_path), "/vsistdout/"]
    result = lsdtt_test_env.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"gdal_translate failed for {raster_path} (exit {result.returncode})."
        )
    got = set()
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        value = float(fields[2])
        if value >= 0.0:
            got.add(value)
    if got != expected:
        raise AssertionError(
            f"Raster values mismatch for {raster_path}: expected {sorted(expected)}, "
            f"got {sorted(got)}."
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


# Files that legitimately differ between a run and its re-run: they name the
# write_path, or record the run itself.
RERUN_DIFFERENT_BY_DESIGN = ("_README.csv", "_README.html",
                             "_ingested_params.param", "_run_summary.json",
                             # lsdtt-fluvial-surfaces provenance: its source
                             # column says "preset" on the first run and
                             # "parameter file" on the rerun, which reads the
                             # resolved values back from the ingested file
                             "_method_record.csv",
                             # GDAL statistics sidecars, written by this
                             # script's own raster checks after the run
                             ".aux.xml")


def rerun_from_ingested(exe: pathlib.Path, work_dir: pathlib.Path,
                        out_dir: pathlib.Path, stdout: str) -> str:
    """Re-run the tool from its own <write_fname>_ingested_params.param.

    The ingested file records every value the run used, read or default, so
    running the tool on it (with write_path moved to a fresh directory) must
    write the same files with the same bytes. A default that does not survive
    the round trip (e.g. a float written with too few digits) shows up here.
    Also checks the re-run's ingested file equals the first one apart from
    write_path, i.e. the record is a fixed point. Returns "" or the problem.
    """
    m = re.search(r"summary: (.+)$", stdout.strip().splitlines()[-1])
    summary = json.loads(pathlib.Path(m.group(1)).read_text(encoding="utf-8"))
    write_fname = summary["write_fname"]
    ingested = out_dir / f"{write_fname}_ingested_params.param"
    if not ingested.is_file():
        return f"no {ingested.name} was written"

    rerun_dir = work_dir / "rerun"
    out2 = rerun_dir / "out"
    if rerun_dir.exists():
        shutil.rmtree(rerun_dir)
    out2.mkdir(parents=True)

    def moved(text: str) -> list[str]:
        return [f"write_path: {out2}" if l.startswith("write_path:") else l
                for l in text.splitlines()]

    def values(text: str) -> dict[str, str]:
        # key -> value, ignoring comments, sections and write_path: a value
        # read from the file and the same value as a default are the same.
        out = {}
        for line in text.splitlines():
            if not line.strip() or line.lstrip().startswith("#") or ":" not in line:
                continue
            key, value = (part.strip() for part in line.split(":", 1))
            if key != "write_path":
                out[key] = value
        return out

    first_ingested = ingested.read_text(encoding="utf-8")
    (rerun_dir / "rerun_params.txt").write_text("\n".join(moved(first_ingested)) + "\n",
                                                 encoding="utf-8")
    result = lsdtt_test_env.run_logged([str(exe), str(rerun_dir), "rerun_params.txt"],
                                       rerun_dir, cwd=rerun_dir)
    if result.returncode != 0:
        return f"the re-run exited {result.returncode}"

    def outputs(root: pathlib.Path) -> dict[str, pathlib.Path]:
        found = {}
        for entry in root.iterdir():
            if not entry.name.startswith(write_fname + "_"):
                continue
            files = [entry] if entry.is_file() else [p for p in entry.rglob("*") if p.is_file()]
            for p in files:
                rel = p.relative_to(root).as_posix()
                if not rel.endswith(RERUN_DIFFERENT_BY_DESIGN):
                    found[rel] = p
        return found

    a, b = outputs(out_dir), outputs(out2)
    if set(a) != set(b):
        return (f"different output files: only in the first run "
                f"{sorted(set(a) - set(b))[:5]}, only in the re-run "
                f"{sorted(set(b) - set(a))[:5]}")
    differ = [rel for rel in sorted(a) if a[rel].read_bytes() != b[rel].read_bytes()]
    if differ:
        return f"outputs with different bytes: {differ[:5]}"

    second = out2 / ingested.name
    if not second.is_file():
        return f"the re-run wrote no {ingested.name}"
    v1, v2 = values(first_ingested), values(second.read_text(encoding="utf-8"))
    if v1 != v2:
        diff = sorted(k for k in set(v1) | set(v2) if v1.get(k) != v2.get(k))
        return ("the re-run's ingested parameters differ: "
                + ", ".join(f"{k}: {v1.get(k)!r} -> {v2.get(k)!r}" for k in diff[:5]))
    return ""


def check_run_record(result, out_dir: pathlib.Path, expect_failure: bool) -> str:
    """Check how the run ended (lsdtt_run_record). Returns "" or the problem.

    A failed run shows the boxed RUN FAILED message on stdout. A successful
    run ends with "lsdtt3: <tool> finished OK; <n> output files; summary:
    <path>", and that summary says status ok and lists every file the run
    wrote under write_path in an entry named <write_fname>_...
    """
    if expect_failure:
        if "|| RUN FAILED: " not in result.stdout:
            return "a failed run did not print the RUN FAILED box on stdout"
        return ""

    lines = [l for l in result.stdout.splitlines() if l.strip()]
    last = lines[-1] if lines else ""
    m = re.match(r"lsdtt3: \S+ finished OK; (\d+) output files?; summary: (.+)$", last)
    if not m:
        return f"the last stdout line is not the finished-OK line: {last!r}"
    summary_path = pathlib.Path(m.group(2))
    if not summary_path.is_file():
        return f"the run summary named on the last line does not exist: {summary_path}"
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return f"the run summary is not valid JSON: {exc}"
    if summary.get("status") != "ok":
        return f"the run summary status is {summary.get('status')!r}, not 'ok'"
    listed = {o["path"] for o in summary.get("outputs", [])}
    if int(m.group(1)) != len(listed):
        return (f"the last line counts {m.group(1)} outputs, the summary lists "
                f"{len(listed)}")
    write_path = pathlib.Path(summary["write_path"])
    prefix = summary["write_fname"] + "_"
    written = set()
    if write_path.resolve() == out_dir.resolve():
        for entry in out_dir.iterdir():
            if not entry.name.startswith(prefix):
                continue
            files = [entry] if entry.is_file() else [p for p in entry.rglob("*") if p.is_file()]
            for p in files:
                rel = p.relative_to(out_dir).as_posix()
                if rel != summary_path.name:
                    written.add(rel)
    missing = sorted(written - listed)
    if missing:
        return f"the run summary does not list these outputs: {missing[:5]}"
    return ""


def main() -> int:
    args = parse_args()

    exe = pathlib.Path(args.exe).resolve()
    template_path = pathlib.Path(args.param_template).resolve()
    read_path = pathlib.Path(args.read_path).resolve()
    work_dir = pathlib.Path(args.work_dir).resolve()
    out_dir = work_dir / "out"

    if not exe.exists():
        raise FileNotFoundError(f"Executable not found: {exe}")
    assert_binary_is_current(
        exe,
        pathlib.Path(args.cli_src_dir).resolve() if args.cli_src_dir else None,
    )
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

    # Streamed to stdout.log / stderr.log while the CLI runs, with command.txt
    # written first, so a hung or killed test leaves its evidence on disk.
    cmd = [str(exe), str(work_dir), param_name]
    result = lsdtt_test_env.run_logged(cmd, work_dir, cwd=work_dir)

    if args.expect_exit_code is not None and args.expect_exit_code != 0:
        args.expect_failure = True
    if args.expect_exit_code is not None and result.returncode != args.expect_exit_code:
        print(
            f"Validation failed: expected exit code {args.expect_exit_code}, "
            f"got {result.returncode}",
            file=sys.stderr,
        )
        print(f"Command: {' '.join(cmd)}", file=sys.stderr)
        print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
        print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
        return 4

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

    for needle in args.stderr_not_contains:
        if needle in stderr_text:
            print(
                f"Validation failed: stderr contains forbidden text: {needle}",
                file=sys.stderr,
            )
            print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
            print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
            return 5

    record_problem = check_run_record(result, out_dir, args.expect_failure)
    if record_problem:
        print(f"Validation failed: {record_problem}", file=sys.stderr)
        print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
        print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
        return 6

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
    present = [str(out_dir / name) for name in args.expect_absent
               if (out_dir / name).exists()]

    if missing or empty or present:
        if missing:
            print("Missing outputs:", file=sys.stderr)
            for m in missing:
                print(f"  {m}", file=sys.stderr)
        if empty:
            print("Empty outputs:", file=sys.stderr)
            for e in empty:
                print(f"  {e}", file=sys.stderr)
        if present:
            print("Outputs that should not have been written:", file=sys.stderr)
            for p in present:
                print(f"  {p}", file=sys.stderr)
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

        for spec in args.expect_raster_stat:
            rel_path, stat, low_spec, high_spec = spec.rsplit(":", 3)
            assert_raster_stat(out_dir / rel_path, stat.lower(),
                               float(low_spec), float(high_spec))

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

        for spec in args.expect_raster_values:
            rel_path, values_spec = spec.rsplit(":", 1)
            expected = {float(v) for v in values_spec.split(",") if v}
            assert_raster_values(out_dir / rel_path, expected)

        for spec in args.expect_sha256:
            rel_path, hash_hex = spec.rsplit(":", 1)
            assert_sha256(out_dir / rel_path, hash_hex)

        for spec in args.expect_file_contains:
            rel_path, needle = spec.split(":", 1)
            assert_file_contains(out_dir / rel_path, needle)

        for spec in args.expect_file_not_contains:
            rel_path, needle = spec.split(":", 1)
            assert_file_not_contains(out_dir / rel_path, needle)
    except Exception as exc:
        print(f"Validation check failed: {exc}", file=sys.stderr)
        print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
        print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
        return 3

    if os.environ.get("LSDTT_RERUN_FROM_INGESTED") == "1":
        problem = rerun_from_ingested(exe, work_dir, out_dir, result.stdout)
        if problem:
            print(f"Rerun from ingested parameters failed: {problem}", file=sys.stderr)
            print(f"Rerun logs: {work_dir / 'rerun'}", file=sys.stderr)
            return 7
        print(f"Rerun from ingested parameters reproduced the outputs: {exe.name}")

    print(f"PASS: {exe.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
