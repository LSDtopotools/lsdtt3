#!/usr/bin/env python3
"""Real-landscape cosmogenic validation harness (cosmo plan Phase 6 / §5b).

"Download -> run -> delete" tier. For each dataset it:

  1. reads the small, tracked manifest + sample table under
     tests/expected/cosmo/published/;
  2. fetches the DEM KEYLESS via the lsdtt-fetch-raster CLI (Copernicus GLO-30
     from AWS Open Data / Microsoft Planetary Computer -- NO API key), clipped to
     the manifest bounding box and warped to the manifest's projected CRS, into a
     TEMPORARY working directory;
  3. runs lsdtt-cosmo-tool to invert the measured concentrations to
     basin-averaged erosion rates;
  4. compares against the published rates and APPENDS a small summary CSV;
  5. **deletes the working directory -- including on failure (cleanup trap)** so
     no multi-hundred-MB DEM survives a crashed run.

The DEM and all intermediate rasters are NEVER committed. Only the manifest,
sample table, and frozen summary (a few kB each) are tracked.

**NOT in CI.** This needs network and time. The DEM fetch is keyless (Copernicus
GLO-30 / Planetary Computer, the same pathway as lsdtt-fetch-raster), so no
OpenTopography key is required. When the network is unavailable the harness exits
77 (CTest "skip"), never a red failure. The CTest label is `cosmo_realworld`,
kept OUT of the default `unit;cosmo` suite.

Usage:
  python3 tests/scripts/run_cosmo_validation.py <dataset|all> \
      --exe build/lsdtt-cosmo-tool --fetch-exe build/lsdtt-fetch-raster \
      [--keep] [--published-dir DIR]

  python3 tests/scripts/run_cosmo_validation.py --self-test
      # offline; exercises the cleanup-trap + manifest/schema logic with a mock
      # workdir and proves the trap fires on failure. Needs no network or key.

See docs/plan_docs/cosmo_cli_testing_plan.md §5b and tests/expected/cosmo/
published/README.md.
"""

import argparse
import csv
import datetime
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Exit codes: 0 pass, 1 real failure, 77 skip (network/key/dep unavailable),
# matching the CTest SKIP_RETURN_CODE convention used elsewhere in the suite.
EXIT_OK = 0
EXIT_FAIL = 1
EXIT_SKIP = 77

TODO_SENTINEL = "TODO"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_published_dir() -> Path:
    return repo_root() / "tests" / "expected" / "cosmo" / "published"


def git_sha() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(repo_root()),
            stderr=subprocess.DEVNULL,
        )
        return out.decode().strip()
    except Exception:  # noqa: BLE001 - metadata only
        return "unknown"


# --------------------------------------------------------------------------- #
# Manifest + sample table (tracked, small)
# --------------------------------------------------------------------------- #
def load_manifest(dataset: str, published_dir: Path) -> dict:
    path = published_dir / f"{dataset}_manifest.json"
    if not path.is_file():
        raise FileNotFoundError(f"no manifest for dataset '{dataset}': {path}")
    with open(path) as f:
        # Strip nothing; the manifest is strict JSON (comments live in _comment
        # keys, not // lines, so a plain parser suffices and stays dependency-free).
        manifest = json.load(f)
    for key in ("dataset", "dem", "sample_table", "summary"):
        if key not in manifest:
            raise ValueError(f"manifest {path} missing required key '{key}'")
    return manifest


def _row_has_todo(row: dict) -> bool:
    return any(str(v).strip() == TODO_SENTINEL for v in row.values())


def read_sample_table(path: Path) -> list:
    """Read the per-sample table, skipping '#' comment lines.

    Raises if the table still carries the TODO sentinel: an un-populated
    template must never silently drive a comparison (plan §5b, README).
    """
    if not path.is_file():
        raise FileNotFoundError(f"sample table not found: {path}")
    with open(path, newline="") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    rows = list(reader)
    if not rows:
        raise ValueError(f"sample table has no data rows: {path}")
    if any(_row_has_todo(r) for r in rows):
        raise ValueError(
            f"sample table {path} is still a TEMPLATE (contains the "
            f"'{TODO_SENTINEL}' sentinel). Populate it from the data repository "
            f"named in the manifest before running validation -- do NOT invent "
            f"or paper-scrape the numbers (plan §5b)."
        )
    return rows


def manifest_has_unresolved_todo(manifest: dict) -> bool:
    """True if any manifest value still carries the TODO sentinel."""
    def walk(v):
        if isinstance(v, dict):
            return any(walk(x) for k, x in v.items() if k != "_comment")
        if isinstance(v, list):
            return any(walk(x) for x in v)
        return isinstance(v, str) and v.strip().startswith(TODO_SENTINEL)

    return walk(manifest)


# --------------------------------------------------------------------------- #
# DEM fetch (keyless, large, temporary, never committed)
# --------------------------------------------------------------------------- #
class SkipError(RuntimeError):
    """Raised when a prerequisite (network / fetch) is missing -> exit 77, not
    fail."""


def write_fetch_param(manifest: dict, workdir: Path) -> Path:
    """Write an lsdtt-fetch-raster parameter file that clips the manifest bbox
    keyless (Copernicus GLO-30 / Planetary Computer, no API key)."""
    dem = manifest["dem"]
    bbox = dem["bounding_box"]
    write_fname = f"{manifest['dataset']}"
    lines = [
        f"read_path: {workdir}",
        f"write_path: {workdir}",
        f"write_fname: {write_fname}",
        f"west: {bbox['west']}",
        f"south: {bbox['south']}",
        f"east: {bbox['east']}",
        f"north: {bbox['north']}",
        f"dem_source: {dem.get('dem_source', 'cop30_aws')}",
        f"grid_spacing: {dem.get('grid_spacing', 30.0)}",
        f"resampling_method: {dem.get('resampling_method', 'cubic')}",
    ]
    if dem.get("target_epsg"):
        lines.append(f"target_epsg: {dem['target_epsg']}")
    if dem.get("pc_collection"):
        lines.append(f"pc_collection: {dem['pc_collection']}")
    if dem.get("size_pixels"):
        lines.append(f"size_pixels: {dem['size_pixels']}")
    param = workdir / "fetch.param"
    param.write_text("\n".join(lines) + "\n")
    return param


def fetch_dem(fetch_exe: Path, manifest: dict, workdir: Path) -> Path:
    """Drive lsdtt-fetch-raster to fetch+clip the DEM into workdir. A network or
    fetch failure is a SKIP (exit 77), not a red failure -- this tier needs the
    network by design."""
    param = write_fetch_param(manifest, workdir)
    try:
        subprocess.run([str(fetch_exe), str(workdir), param.name],
                       check=True, cwd=str(workdir))
    except subprocess.CalledProcessError as exc:
        raise SkipError(
            f"lsdtt-fetch-raster failed (rc={exc.returncode}); the DEM source "
            f"is likely unreachable. Skipping this dataset."
        )
    out_path = workdir / f"{manifest['dataset']}_DEM.tif"
    if not out_path.is_file() or out_path.stat().st_size == 0:
        raise SkipError("DEM fetch produced no/empty raster; skipping.")
    return out_path


# --------------------------------------------------------------------------- #
# Run + compare
# --------------------------------------------------------------------------- #
def run_cosmo_tool(exe: Path, dem_path: Path, sample_table: Path,
                   workdir: Path) -> Path:
    """Invert concentrations with lsdtt-cosmo-tool; return the _CN_ERATE.csv."""
    # The CLI reads a parameter file; write a minimal one into the workdir.
    points_csv = workdir / "points.csv"
    shutil.copyfile(sample_table, points_csv)
    write_base = "validation"
    param = workdir / "cosmo_validation.param"
    param.write_text(
        f"read path: {workdir}\n"
        f"read fname: {dem_path.stem}\n"
        f"write path: {workdir}\n"
        f"write fname: {write_base}\n"
        f"nuclide: Be10\n"
        f"points_filename: {points_csv.name}\n"
    )
    subprocess.run([str(exe), str(workdir), param.name], check=True, cwd=str(workdir))
    out = workdir / f"{write_base}_CN_ERATE.csv"
    if not out.is_file():
        raise RuntimeError(f"cosmo-tool produced no output: {out}")
    return out


def compare(lsdtt_csv: Path, samples: list, manifest: dict) -> list:
    """Join lsdtt results to published rates by row order; build summary rows."""
    with open(lsdtt_csv, newline="") as f:
        results = list(csv.DictReader(f))
    if len(results) != len(samples):
        raise RuntimeError(
            f"result count ({len(results)}) != sample count ({len(samples)})"
        )
    dem = manifest["dem"]
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    sha = git_sha()
    out_rows = []
    for res, samp in zip(results, samples):
        lsdtt_e = float(res["apparent_eff_e_g_cm2_yr"])
        pub_e = float(samp["published_erate_gcm2yr"])
        pct = 100.0 * (lsdtt_e - pub_e) / pub_e if pub_e else float("nan")
        out_rows.append({
            "sample_name": samp.get("sample_name", ""),
            "nuclide": "Be-10",
            "standardisation": samp.get("standardisation", ""),
            "lsdtt_erate_gcm2yr": f"{lsdtt_e:.8g}",
            "published_erate_gcm2yr": f"{pub_e:.8g}",
            "pct_diff": f"{pct:.4g}",
            "published_calculator": samp.get("published_calculator", ""),
            "dem_source": f"{dem.get('provider')}:{dem.get('dataset_id')}",
            "resolution_m": dem.get("resolution_m", ""),
            "run_date_utc": now,
            "git_sha": sha,
        })
    return out_rows


SUMMARY_FIELDS = [
    "sample_name", "nuclide", "standardisation", "lsdtt_erate_gcm2yr",
    "published_erate_gcm2yr", "pct_diff", "published_calculator",
    "dem_source", "resolution_m", "run_date_utc", "git_sha",
]


def append_summary(summary_path: Path, rows: list) -> None:
    """Append rows to the tracked summary CSV, preserving its comment header and
    the field order frozen in the skeleton."""
    with open(summary_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        for r in rows:
            w.writerow({k: r.get(k, "") for k in SUMMARY_FIELDS})


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def validate_one(dataset: str, args) -> int:
    published_dir = args.published_dir
    manifest = load_manifest(dataset, published_dir)
    if manifest_has_unresolved_todo(manifest):
        raise SkipError(
            f"manifest for '{dataset}' still has unresolved TODO fields "
            f"(DEM footprint / data repo not confirmed); populate it first."
        )
    sample_table = published_dir / manifest["sample_table"]
    samples = read_sample_table(sample_table)  # raises on TODO template
    summary_path = published_dir / manifest["summary"]

    workdir = Path(tempfile.mkdtemp(prefix=f"cosmo_val_{dataset}_"))
    print(f"[{dataset}] working directory: {workdir}")
    try:
        if args.exe is None:
            raise ValueError("--exe (path to lsdtt-cosmo-tool) is required")
        if args.fetch_exe is None:
            raise ValueError("--fetch-exe (path to lsdtt-fetch-raster) is required")
        dem_path = fetch_dem(Path(args.fetch_exe), manifest, workdir)
        print(f"[{dataset}] DEM: {dem_path} ({dem_path.stat().st_size} bytes)")
        lsdtt_csv = run_cosmo_tool(Path(args.exe), dem_path, sample_table, workdir)
        rows = compare(lsdtt_csv, samples, manifest)
        append_summary(summary_path, rows)
        for r in rows:
            print(f"[{dataset}] {r['sample_name']}: lsdtt="
                  f"{r['lsdtt_erate_gcm2yr']} published="
                  f"{r['published_erate_gcm2yr']} pct_diff={r['pct_diff']}")
        print(f"[{dataset}] appended {len(rows)} rows to {summary_path}")
        return EXIT_OK
    finally:
        # CLEANUP TRAP: always fires, including on exception, unless --keep.
        if args.keep:
            print(f"[{dataset}] --keep set: leaving {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)
            print(f"[{dataset}] cleaned up {workdir} (exists={workdir.exists()})")


def cleanup_trap_fires(keep: bool, force_failure: bool, seen: list) -> None:
    """Reusable copy of the try/finally cleanup contract, isolated so the
    self-test can exercise it without any network. Records the workdir path in
    `seen` BEFORE anything can fail, so the caller can inspect it afterwards
    (it must be gone unless keep=True) even when force_failure raises."""
    workdir = Path(tempfile.mkdtemp(prefix="cosmo_val_selftest_"))
    seen.append(workdir)
    (workdir / "dem.tif").write_bytes(b"x" * 1024)  # stand-in for a huge DEM
    try:
        if force_failure:
            raise RuntimeError("simulated mid-run failure (network/CLI crash)")
    finally:
        if not keep:
            shutil.rmtree(workdir, ignore_errors=True)


def self_test() -> int:
    """Offline: prove the cleanup trap + manifest/schema guards behave. Prints
    PASS/FAIL per check, like the C++ cosmo driver."""
    published_dir = default_published_dir()
    checks = []

    def check(name, ok):
        checks.append(ok)
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")

    # 1. Cleanup trap fires on a FAILED run: the simulated failure must
    #    propagate AND the workdir (containing a stand-in DEM) must be gone.
    seen = []
    try:
        cleanup_trap_fires(keep=False, force_failure=True, seen=seen)
        raised = False
    except RuntimeError:
        raised = True
    failed_wd = seen[-1]
    check("cleanup trap: exception propagates from a failed run", raised)
    check("cleanup trap: workdir deleted after a failed run",
          not failed_wd.exists())

    # 2. --keep preserves the workdir, then we tidy it ourselves.
    seen = []
    cleanup_trap_fires(keep=True, force_failure=False, seen=seen)
    kept = seen[-1]
    check("cleanup trap: --keep leaves the workdir in place", kept.exists())
    shutil.rmtree(kept, ignore_errors=True)

    # 3. Manifest loads and validates for the pilot dataset.
    try:
        manifest = load_manifest("dethier2014", published_dir)
        loaded = True
    except Exception as exc:  # noqa: BLE001
        loaded = False
        print(f"       manifest load error: {exc}")
    check("manifest: dethier2014 loads with required keys", loaded)

    # 4. The pilot manifest is correctly flagged as still having TODOs (so the
    #    harness will SKIP rather than run against an unfinished manifest).
    if loaded:
        check("manifest: unresolved-TODO guard trips on the pilot template",
              manifest_has_unresolved_todo(manifest))

    # 5. The sample-table TEMPLATE is rejected (TODO sentinel), so a bogus
    #    comparison can never be produced from an un-populated table.
    rejected = False
    try:
        read_sample_table(published_dir / "dethier2014_samples.csv")
    except ValueError:
        rejected = True
    check("sample table: TODO template is refused (no fabricated comparison)",
          rejected)

    # 6. Summary append round-trips against a COPY of the frozen skeleton, and
    #    leaves the tracked file untouched.
    tracked_summary = published_dir / "dethier2014_summary.csv"
    before = tracked_summary.read_text()
    with tempfile.TemporaryDirectory() as td:
        tmp_summary = Path(td) / "summary.csv"
        tmp_summary.write_text(before)
        append_summary(tmp_summary, [{
            "sample_name": "MOCK-1", "nuclide": "Be-10",
            "standardisation": "07KNSTD", "lsdtt_erate_gcm2yr": "1.23e-4",
            "published_erate_gcm2yr": "1.20e-4", "pct_diff": "2.5",
            "published_calculator": "CRONUS", "dem_source": "OpenTopography:SRTMGL3",
            "resolution_m": 90, "run_date_utc": "2026-07-24T00:00:00Z",
            "git_sha": "deadbee",
        }])
        with open(tmp_summary, newline="") as f:
            data_rows = [ln for ln in f if ln.strip()
                         and not ln.startswith("#") and "sample_name" not in ln]
        appended_ok = len(data_rows) == 1 and "MOCK-1" in data_rows[0]
    untouched = tracked_summary.read_text() == before
    check("summary: append writes one schema-correct row to a copy", appended_ok)
    check("summary: tracked frozen skeleton left byte-unchanged", untouched)

    # 7. End-to-end validate_one cleanup: a POPULATED (no-TODO) manifest + a
    #    fetch-exe that fails mid-run must (a) SKIP, not fail, and (b) leave no
    #    workdir behind -- exercising validate_one's own try/finally, not just
    #    the mirrored helper above.
    import glob
    with tempfile.TemporaryDirectory() as td:
        pdir = Path(td)
        (pdir / "mock_manifest.json").write_text(json.dumps({
            "dataset": "mock",
            "dem": {"dem_source": "cop30_aws", "provider": "COP30",
                    "dataset_id": "COP30", "resolution_m": 30,
                    "bounding_box": {"west": 0.0, "south": 0.0,
                                     "east": 0.1, "north": 0.1}},
            "sample_table": "mock_samples.csv",
            "summary": "mock_summary.csv",
        }))
        (pdir / "mock_samples.csv").write_text(
            "sample_name,latitude,longitude,Be10_CONC,Be10_uncertainty,"
            "standardisation,published_erate_gcm2yr,published_erate_mm_kyr,"
            "published_calculator\n"
            "M1,0.05,0.05,1.0e5,3.0e3,07KNSTD,1.0e-4,3.77,CRONUS\n")
        (pdir / "mock_summary.csv").write_text(",".join(SUMMARY_FIELDS) + "\n")

        false_exe = shutil.which("false") or "/bin/false"
        mock_args = argparse.Namespace(
            published_dir=pdir, exe=false_exe, fetch_exe=false_exe, keep=False)
        pattern = str(Path(tempfile.gettempdir()) / "cosmo_val_mock_*")
        before = set(glob.glob(pattern))
        skipped = False
        try:
            validate_one("mock", mock_args)
        except SkipError:
            skipped = True
        except Exception:  # noqa: BLE001
            pass
        leaked = set(glob.glob(pattern)) - before
        check("validate_one: failing fetch is a SKIP, not a hard fail", skipped)
        check("validate_one: cleanup trap removed the real workdir",
              not any(Path(p).exists() for p in leaked))

    ok = all(checks)
    print(f"\nRESULT: {'PASS' if ok else 'FAIL'} "
          f"({sum(checks)}/{len(checks)} checks)")
    return EXIT_OK if ok else EXIT_FAIL


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dataset", nargs="?",
                    help="dataset name (e.g. dethier2014) or 'all'")
    ap.add_argument("--exe", default=os.environ.get("LSDTT_COSMO_TOOL"),
                    help="path to lsdtt-cosmo-tool")
    ap.add_argument("--fetch-exe", default=os.environ.get("LSDTT_FETCH_RASTER"),
                    help="path to lsdtt-fetch-raster (keyless DEM fetch)")
    ap.add_argument("--published-dir", type=Path, default=default_published_dir())
    ap.add_argument("--keep", action="store_true",
                    help="keep the working directory for debugging (default: delete)")
    ap.add_argument("--self-test", action="store_true",
                    help="offline harness self-test (cleanup trap + schema guards)")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    if not args.dataset:
        ap.error("a dataset name (or 'all', or --self-test) is required")

    if args.dataset == "all":
        datasets = sorted(
            p.name[: -len("_manifest.json")]
            for p in args.published_dir.glob("*_manifest.json")
        )
    else:
        datasets = [args.dataset]

    rc = EXIT_OK
    for ds in datasets:
        try:
            r = validate_one(ds, args)
            rc = rc or r
        except SkipError as exc:
            print(f"[{ds}] SKIP: {exc}")
            if rc == EXIT_OK:
                rc = EXIT_SKIP
        except Exception as exc:  # noqa: BLE001 - report and fail this dataset
            print(f"[{ds}] FAIL: {exc}")
            rc = EXIT_FAIL
    return rc


if __name__ == "__main__":
    sys.exit(main())
