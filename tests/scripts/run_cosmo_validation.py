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
def write_crn_input(samples: list, workdir: Path) -> Path:
    """Translate our tracked sample table into the points CSV lsdtt-cosmo-tool
    reads. The CLI requires `latitude`, `longitude`, a `sample_name` column, and
    (for nuclide Be10, by default) `Be10_CONC` / `Be10_uncertainty` /
    `standardisation` columns. We emit exactly those so the CLI's defaults apply --
    we do NOT blindly copy the tracked table, whose extra provenance columns
    (elevation, published rate) and comment header are ours, not the CLI's."""
    points_csv = workdir / "points.csv"
    fields = ["sample_name", "latitude", "longitude", "Be10_CONC",
              "Be10_uncertainty", "standardisation"]
    with open(points_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for s in samples:
            w.writerow({
                "sample_name": s.get("sample_name", ""),
                "latitude": s.get("latitude", ""),
                "longitude": s.get("longitude", ""),
                "Be10_CONC": s.get("Be10_CONC", ""),
                "Be10_uncertainty": s.get("Be10_uncertainty", ""),
                "standardisation": s.get("standardisation", ""),
            })
    return points_csv


def run_cosmo_tool(exe: Path, dem_path: Path, samples: list,
                   workdir: Path, run_params: dict = None) -> Path:
    """Invert concentrations with lsdtt-cosmo-tool; return the _CN_ERATE.csv.

    The lsdtt param parser uses UNDERSCORE keys (read_fname/read_path/...); with
    a space ("read fname") the key is ignored and read_fname silently defaults to
    the parameter-file stem -- which is exactly the DEM-not-found trap this hit on
    the first networked run."""
    run_params = run_params or {}
    points_csv = write_crn_input(samples, workdir)
    write_base = "validation"
    param = workdir / "cosmo_validation.param"
    lines = [
        f"read_path: {workdir}",
        f"read_fname: {dem_path.name}",       # parser strips the .tif extension
        f"write_path: {workdir}",
        f"write_fname: {write_base}",
        "nuclide: Be10",
        f"points_filename: {points_csv.name}",
    ]
    # Optional per-dataset run knobs (e.g. a channel-snap threshold suited to the
    # basin sizes -- the 1000-pixel default is wrong for sub-km2 catchments).
    for key, val in run_params.items():
        if key == "_comment":
            continue
        lines.append(f"{key}: {val}")
    param.write_text("\n".join(lines) + "\n")
    subprocess.run([str(exe), str(workdir), param.name], check=True, cwd=str(workdir))
    out = workdir / f"{write_base}_CN_ERATE.csv"
    if not out.is_file():
        raise RuntimeError(f"cosmo-tool produced no output: {out}")
    return out


def compare(lsdtt_csv: Path, samples: list, manifest: dict) -> list:
    """Join lsdtt results to published rates BY sample_name (channel snapping may
    drop or reorder samples, so a positional zip is unsafe); build summary rows."""
    with open(lsdtt_csv, newline="") as f:
        results = {r["sample_name"]: r for r in csv.DictReader(f)}
    missing = [s["sample_name"] for s in samples
               if s["sample_name"] not in results]
    if missing:
        print(f"WARNING: {len(missing)} sample(s) absent from cosmo-tool output "
              f"(snapping/extent): {', '.join(missing)}")
    dem = manifest["dem"]
    # Basins whose true (OCTOPUS) area is below the snap-threshold footprint cannot
    # be delineated faithfully -- the outlet snaps onto a larger channel, so the
    # comparison for them reflects delineation mismatch, not the numerics. Flag
    # them so the frozen band can be read on the reliable subset.
    res_m = float(dem.get("resolution_m", 30) or 30)
    thr_px = float((manifest.get("run_params") or {}).get(
        "threshold_contributing_pixels", 1000) or 1000)
    min_reliable_km2 = thr_px * res_m * res_m / 1.0e6
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    sha = git_sha()
    out_rows = []
    for samp in samples:
        res = results.get(samp["sample_name"])
        if res is None:
            continue
        lsdtt_e = float(res["apparent_eff_e_g_cm2_yr"])
        pub_e = float(samp["published_erate_gcm2yr"])
        pct = 100.0 * (lsdtt_e - pub_e) / pub_e if pub_e else float("nan")
        area = samp.get("octopus_area_km2", "")
        try:
            reliable = "yes" if float(area) >= min_reliable_km2 else "no"
        except (TypeError, ValueError):
            reliable = ""
        out_rows.append({
            "sample_name": samp.get("sample_name", ""),
            "nuclide": "Be-10",
            "standardisation": samp.get("standardisation", ""),
            # elevation carried through so the residual-vs-elevation trend is
            # visible in the tracked summary (Simon 2026-07-24: residual is
            # expected to be elevation-dependent, see manifest residual_expectation).
            "elevation_ave_m": samp.get("elevation_ave_m", ""),
            "octopus_area_km2": area,
            # "no" => basin smaller than the snap-threshold footprint; its rate
            # reflects delineation mismatch, exclude from the residual band.
            "delineation_reliable": reliable,
            "lsdtt_erate_gcm2yr": f"{lsdtt_e:.8g}",
            "published_erate_gcm2yr": f"{pub_e:.8g}",
            "published_erate_mm_kyr": samp.get("published_erate_mm_kyr", ""),
            "pct_diff": f"{pct:.4g}",
            "abs_pct_diff": f"{abs(pct):.4g}",
            "published_calculator": samp.get("published_calculator", ""),
            "dem_source": f"{dem.get('provider')}:{dem.get('dataset_id')}",
            "resolution_m": dem.get("resolution_m", ""),
            "run_date_utc": now,
            "git_sha": sha,
        })
    return out_rows


SUMMARY_FIELDS = [
    "sample_name", "nuclide", "standardisation", "elevation_ave_m",
    "octopus_area_km2", "delineation_reliable",
    "lsdtt_erate_gcm2yr", "published_erate_gcm2yr", "published_erate_mm_kyr",
    "pct_diff", "abs_pct_diff", "published_calculator",
    "dem_source", "resolution_m", "run_date_utc", "git_sha",
]


def read_frozen_summary(summary_path: Path) -> dict:
    """Read the tracked (frozen) summary's data rows, keyed by sample_name."""
    with open(summary_path, newline="") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    return {r["sample_name"]: r for r in csv.DictReader(lines)}


# Regression tolerance: a rerun on the same code + the deterministic COP30 fetch
# should reproduce the frozen inverted rate closely. Allow a little slack for
# fetch/interpolation non-determinism, not a physics band.
CHECK_REL_TOL = 0.02


def check_against_frozen(new_rows: list, summary_path: Path) -> int:
    """Compare a fresh run against the frozen tracked summary WITHOUT modifying it
    (the plan's cheap tracked regression check). Returns EXIT_OK / EXIT_FAIL."""
    frozen = read_frozen_summary(summary_path)
    if not frozen:
        print("CHECK: no frozen rows to compare against; run once without --check "
              "to freeze a reference first.")
        return EXIT_FAIL
    ok = True
    for r in new_rows:
        name = r["sample_name"]
        ref = frozen.get(name)
        if ref is None:
            print(f"[CHECK][FAIL] {name}: not in frozen summary")
            ok = False
            continue
        new_e = float(r["lsdtt_erate_gcm2yr"])
        ref_e = float(ref["lsdtt_erate_gcm2yr"])
        rel = abs(new_e - ref_e) / ref_e if ref_e else float("inf")
        status = "PASS" if rel <= CHECK_REL_TOL else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"[CHECK][{status}] {name}: lsdtt={new_e:.6g} frozen={ref_e:.6g} "
              f"rel={rel:.3g}")
    missing = set(frozen) - {r["sample_name"] for r in new_rows}
    if missing:
        print(f"[CHECK][FAIL] frozen samples absent from this run: "
              f"{', '.join(sorted(missing))}")
        ok = False
    print(f"CHECK: {'PASS' if ok else 'FAIL'} "
          f"(tol {CHECK_REL_TOL:.0%}, {len(new_rows)} samples)")
    return EXIT_OK if ok else EXIT_FAIL


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
        lsdtt_csv = run_cosmo_tool(Path(args.exe), dem_path, samples, workdir,
                                   manifest.get("run_params"))
        rows = compare(lsdtt_csv, samples, manifest)
        for r in rows:
            print(f"[{dataset}] {r['sample_name']}: lsdtt="
                  f"{r['lsdtt_erate_gcm2yr']} published="
                  f"{r['published_erate_gcm2yr']} pct_diff={r['pct_diff']}")
        if getattr(args, "check", False):
            # Regression mode: compare against the frozen summary, do NOT touch it.
            return check_against_frozen(rows, summary_path)
        append_summary(summary_path, rows)
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

    # 4a. The pilot manifest is now fully populated (OCTOPUS-harvested), so the
    #     unresolved-TODO guard must be CLEAN -- the harness will run it, not skip.
    if loaded:
        check("manifest: pilot dethier2014 has no unresolved TODO (ready to run)",
              not manifest_has_unresolved_todo(manifest))

    # 4b. But the guard must still TRIP on an un-finished (TODO) manifest, so an
    #     unpopulated future dataset can never silently run.
    check("manifest: unresolved-TODO guard still trips on a TODO template",
          manifest_has_unresolved_todo({"dem": {"target_epsg": "TODO: set me"}}))

    # 5a. The pilot sample table is now populated (real OCTOPUS rows) and reads
    #     cleanly -- no TODO sentinel, at least one data row.
    real_rows = []
    try:
        real_rows = read_sample_table(published_dir / "dethier2014_samples.csv")
        table_ok = len(real_rows) > 0
    except Exception as exc:  # noqa: BLE001
        table_ok = False
        print(f"       sample table error: {exc}")
    check("sample table: pilot dethier2014 reads real rows (no TODO sentinel)",
          table_ok)

    # 5b. A TODO template is still refused, so a bogus comparison can never be
    #     produced from an un-populated table.
    rejected = False
    with tempfile.TemporaryDirectory() as td:
        tmpl = Path(td) / "tmpl.csv"
        tmpl.write_text(
            "sample_name,latitude,longitude,Be10_CONC,Be10_uncertainty,"
            "standardisation,published_erate_gcm2yr\n"
            "TODO,TODO,TODO,TODO,TODO,TODO,TODO\n")
        try:
            read_sample_table(tmpl)
        except ValueError:
            rejected = True
    check("sample table: a TODO template is refused (no fabricated comparison)",
          rejected)

    # 6. Summary append round-trips against a COPY of the tracked summary (which
    #    now holds a frozen accepted run), adding EXACTLY one row and leaving the
    #    tracked file untouched.
    tracked_summary = published_dir / "dethier2014_summary.csv"
    before = tracked_summary.read_text()

    def data_rows_of(text):
        return [ln for ln in text.splitlines() if ln.strip()
                and not ln.startswith("#") and not ln.startswith("sample_name,")]

    n_before = len(data_rows_of(before))
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
        after_rows = data_rows_of(tmp_summary.read_text())
        appended_ok = (len(after_rows) == n_before + 1
                       and "MOCK-1" in after_rows[-1])
    untouched = tracked_summary.read_text() == before
    check("summary: append adds exactly one schema-correct row to a copy",
          appended_ok)
    check("summary: tracked frozen summary left byte-unchanged", untouched)

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
    ap.add_argument("--check", action="store_true",
                    help="regression mode: compare a fresh run against the frozen "
                         "summary WITHOUT modifying it (used by the ctest wiring)")
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
