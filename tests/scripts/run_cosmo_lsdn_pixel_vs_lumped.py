#!/usr/bin/env python3
"""L13 generation step: per-pixel basin-mean production scalings (LSDn and Stone).

WHAT THIS IS (lsdn plan Phase 6 deliverable (ii); cosmo plan §6f row L13)
-------------------------------------------------------------------------
L13 asks how lsdtt3's PRODUCT -- the PER-PIXEL LSDn basin-mean production -- differs
from the LUMPED LSDn value OCTOPUS/RIVERSAND uses, i.e. a single point evaluation at
(CENTR_LAT, ATM_PRESS).  That difference is a SCIENTIFIC RESULT, not an error term:
per-pixel evaluation is lsdtt3's default and its whole point (lsdn plan §3.0), and the
per-pixel number must NEVER be tuned toward the lumped one or toward OCTOPUS.

The work splits in two, deliberately:

  * THIS SCRIPT (expensive, networked, label `cosmo_realworld`) fetches each dataset's
    DEM, runs lsdtt-cosmo-tool over the tracked samples, and records the basin-mean
    production SCALING that `write_basin_production_csv` averages over the EXACT
    upslope node set the inversion uses.  That is quantity (a).  It also records the
    per-pixel STONE basin mean, which is precisely the quantity OCTOPUS inverts to
    obtain ATM_PRESS -- so the frozen artefact carries everything needed to reproduce
    their reduction as well.
  * THE C++ DRIVER (cheap, offline, default `unit;cosmo` suite) reads the frozen
    artefact together with the frozen OCTOPUS point reference, evaluates the LUMPED
    LSDn kernel in process, emits `lsdn_pixel_vs_lumped.csv`, and regression-checks
    the result.  See phase7_L13_pixel_vs_lumped() in tests/cpp/test_cosmo_physics.cpp.

That two-step is the point: the expensive half runs manually and is frozen; the cheap
half runs in CI on every build.

SAMPLE AGE.  lsdn_scaling() is time dependent, and `sample_age_yr` is a single scalar
per tool run, so the per-pixel side is recorded at THREE ages -- 0, 5000 and 20000 yr.
The C++ driver evaluates the lumped side at exactly the same ages, so the (a)-vs-(b)
ratio is always taken with t matched on both sides and the time-averaging factor is
common to both.  What remains is the SPATIAL term, which is what L13 is about.  The
three ages then bound how much of the answer depends on that choice; t = Gamma_sp/eps
for these basins is ~3-16 kyr, so 5000 yr is the headline and the other two bracket it.

MATCHED CONFIGURATION.  muon_scheme Braucher, atmosphere NCEP, 07KNSTD, and the
per-dataset threshold/snap knobs from the manifest -- identical to the §6e isolation
run, so the per-pixel Stone mean here is directly comparable with the BEPROD-matched
result already established there (-0.19 %).

FROZEN-ARTEFACT DISCIPLINE.  tests/expected/cosmo/lsdn_pixel_vs_lumped_perpixel.csv is
TRACKED and FROZEN, so the DEFAULT mode is --check: compare and report, never write.
Pass --freeze to regenerate deliberately, then review the diff -- a numeric change
there is a change in the per-pixel production result.  The '# Run: ... ; git <sha>'
provenance line is neither compared nor rewritten.

Exit 77 = SKIP (network/DEM source unavailable), like the rest of the real-landscape
tier.  DEMs live in a temp workdir with a `finally` cleanup and are NEVER committed.

Usage:
  # regression check (what `ctest -L cosmo_realworld` runs; artefact untouched)
  python3 tests/scripts/run_cosmo_lsdn_pixel_vs_lumped.py all \
      --exe build/lsdtt-cosmo-tool --fetch-exe build/lsdtt-fetch-raster [--keep]

  # deliberate re-freeze (review the diff before committing)
  python3 tests/scripts/run_cosmo_lsdn_pixel_vs_lumped.py all \
      --exe build/lsdtt-cosmo-tool --fetch-exe build/lsdtt-fetch-raster --freeze
"""

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_cosmo_validation as rcv  # reuse fetch/cleanup/sample-table plumbing

EXIT_OK, EXIT_FAIL, EXIT_SKIP = 0, 1, 77

# The three tracked datasets of the real-landscape tier.  L13 covers ALL of them
# rather than only dethier2014 (which is all the §6e isolation script handles),
# because the frozen OCTOPUS point reference L13 is compared against carries all
# 50 basins across the three, and restricting to 14 would leave two thirds of the
# lumped-side oracle unused while adding only two more modest DEM fetches.
DATASETS = ["dethier2014", "clapp2002", "dibiase2010"]

# Sample ages at which the per-pixel LSDn basin mean is recorded.  See the SAMPLE
# AGE note in the module docstring; 5000 is the headline, 0 and 20000 bracket it.
SAMPLE_AGES_YR = [0, 5000, 20000]
NUCLIDES = ["Be10", "Al26", "C14"]

ARTEFACT_NAME = "lsdn_pixel_vs_lumped_perpixel.csv"


def expected_cosmo_dir() -> Path:
    return rcv.repo_root() / "tests" / "expected" / "cosmo"


# --------------------------------------------------------------------------- #
# Run the tool once per (scheme, age) and read back the basin-production CSV
# --------------------------------------------------------------------------- #
def run_tool(exe: Path, dem_path: Path, points_csv: Path, workdir: Path,
             manifest: dict, scheme: str, age_yr: int, nuclide: str) -> dict:
    """One lsdtt-cosmo-tool run; returns {sample_name: basin-prod row}.

    Only `write_basin_production_csv` is switched on -- the production RASTER is
    not needed and writing it would cost time and disk for nothing.
    """
    write_base = f"l13_{scheme}_{nuclide}_{age_yr}"
    param = workdir / f"{write_base}.param"
    run_params = manifest.get("run_params") or {}
    points_for_run = points_csv
    if nuclide != "Be10":
        # The tracked validation sample tables are Be10-based. For this L13
        # generation step we need Al26/C14 production scalings too; those do not
        # depend on the measured concentrations, but the CLI requires a nuclide-
        # specific concentration column to parse points. Build a run-local copy
        # by duplicating Be10 concentrations into the requested nuclide columns.
        points_for_run = workdir / f"points_{nuclide}.csv"
        with open(points_csv, newline="") as f_in, open(points_for_run, "w", newline="") as f_out:
            rows = list(csv.DictReader(f_in))
            fields = [
                "sample_name", "latitude", "longitude",
                f"{nuclide}_CONC", f"{nuclide}_uncertainty", "standardisation",
            ]
            w = csv.DictWriter(f_out, fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow({
                    "sample_name": r.get("sample_name", ""),
                    "latitude": r.get("latitude", ""),
                    "longitude": r.get("longitude", ""),
                    f"{nuclide}_CONC": r.get("Be10_CONC", ""),
                    f"{nuclide}_uncertainty": r.get("Be10_uncertainty", ""),
                    "standardisation": "KNSTD",
                })

    lines = [
        f"read_path: {workdir}",
        f"read_fname: {dem_path.name}",
        f"write_path: {workdir}",
        f"write_fname: {write_base}",
        f"nuclide: {nuclide}",
        f"points_filename: {points_for_run.name}",
        # --- matched to OCTOPUS, identical to the §6e isolation run ---
        "muon_scheme: Braucher",
        "atmosphere: NCEP",
        f"scaling_scheme: {scheme}",
        f"sample_age_yr: {age_yr}",
        "write_basin_production_csv: true",
    ]
    for key, val in run_params.items():
        if key == "_comment":
            continue
        lines.append(f"{key}: {val}")
    param.write_text("\n".join(lines) + "\n")
    subprocess.run([str(exe), str(workdir), param.name], check=True, cwd=str(workdir))
    bp = workdir / f"{write_base}_BASIN_PROD.csv"
    if not bp.is_file():
        raise RuntimeError(f"cosmo-tool produced no basin-prod CSV: {bp}")
    with open(bp, newline="") as f:
        return {r["sample_name"]: r for r in csv.DictReader(f)}


def rows_for_dataset(dataset: str, args) -> list:
    pdir = Path(args.published_dir) if args.published_dir else rcv.default_published_dir()
    manifest = rcv.load_manifest(dataset, pdir)
    samples = rcv.read_sample_table(pdir / manifest["sample_table"])

    workdir = Path(tempfile.mkdtemp(prefix=f"cosmo_l13_{dataset}_"))
    print(f"[{dataset}] working directory: {workdir}")
    try:
        dem_path = rcv.fetch_dem(Path(args.fetch_exe), manifest, workdir)
        print(f"[{dataset}] DEM: {dem_path} ({dem_path.stat().st_size} bytes)")
        points_csv = rcv.write_crn_input(samples, workdir)

        stone = run_tool(Path(args.exe), dem_path, points_csv, workdir,
                         manifest, "stone", 0, "Be10")
        lsdn = {
            nuc: {
                age: run_tool(Path(args.exe), dem_path, points_csv, workdir,
                              manifest, "lsdn", age, nuc)
                for age in SAMPLE_AGES_YR
            }
            for nuc in NUCLIDES
        }
    finally:
        if args.keep:
            print(f"--keep set: leaving {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)
            print(f"[{dataset}] cleaned up {workdir} (exists={workdir.exists()})")

    # Delineation reliability, exactly as run_cosmo_validation.py / the §6e isolation
    # script define it: a basin whose true (OCTOPUS) area is below the snap-threshold
    # footprint cannot be delineated faithfully -- the outlet snaps onto a much larger
    # channel, so its per-pixel mean describes a different landscape.  Flagged, not
    # dropped, so the frozen artefact keeps the evidence and the C++ driver can report
    # the reliable subset separately.
    res_m = float(manifest["dem"].get("resolution_m", 30) or 30)
    thr_px = float((manifest.get("run_params") or {}).get(
        "threshold_contributing_pixels", 1000) or 1000)
    min_reliable_km2 = thr_px * res_m * res_m / 1.0e6

    rows = []
    for s in samples:
        name = s["sample_name"]
        if name not in stone:
            print(f"[{dataset}] WARNING: {name} absent from tool output "
                  f"(snap/extent) -- skipped")
            continue
        st = stone[name]
        try:
            reliable = float(s.get("octopus_area_km2", "")) >= min_reliable_km2
        except (TypeError, ValueError):
            reliable = False
        row = {
            "dataset": dataset,
            "basin": name,
            "delineation_reliable": "yes" if reliable else "no",
            "n_upslope_pixels": st["n_upslope_pixels"],
            "snapped_latitude": f"{float(st['snapped_latitude']):.6f}",
            "snapped_longitude": f"{float(st['snapped_longitude']):.6f}",
            "mean_elev_m": f"{float(st['mean_elev_m']):.4f}",
            "mean_scaling_stone": f"{float(st['mean_prod_scaling']):.8f}",
        }
        ok = True
        for nuc in NUCLIDES:
            for age in SAMPLE_AGES_YR:
                r = lsdn[nuc][age].get(name)
                if r is None:
                    ok = False
                    break
                row[f"mean_scaling_lsdn_{nuc}_t{age}"] = f"{float(r['mean_prod_scaling']):.8f}"
                # The upslope set must not depend on nuclide, scaling scheme or
                # sample age -- only the production values do.
                if r["n_upslope_pixels"] != st["n_upslope_pixels"]:
                    print(f"[{dataset}] FAIL: {name} upslope pixel count changed with "
                          f"nuclide/scheme/age ({st['n_upslope_pixels']} stone vs "
                          f"{r['n_upslope_pixels']} lsdn {nuc} t={age})")
                    ok = False
            if not ok:
                break
        if ok:
            rows.append(row)
            print(f"  [{name}] n={row['n_upslope_pixels']} "
                  f"S_stone={row['mean_scaling_stone']} "
                  f"S_lsdn_Be10(5k)={row['mean_scaling_lsdn_Be10_t5000']} "
                  f"S_lsdn_Al26(5k)={row['mean_scaling_lsdn_Al26_t5000']} "
                  f"S_lsdn_C14(5k)={row['mean_scaling_lsdn_C14_t5000']}")
    return rows


FIELDS = (["dataset", "basin", "delineation_reliable", "n_upslope_pixels",
           "snapped_latitude", "snapped_longitude", "mean_elev_m",
           "mean_scaling_stone"]
          + [f"mean_scaling_lsdn_{n}_t{a}" for n in NUCLIDES for a in SAMPLE_AGES_YR])

# Relative tolerance for the frozen-artefact regression check.  These quantities are
# deterministic given the same DEM and executable; the tolerance absorbs COP30
# tile-server and float-formatting noise, not real drift.  Same value and rationale
# as run_cosmo_production_isolation.py's ISOLATION_CHECK_REL_TOL.
CHECK_REL_TOL = 1e-3
CHECK_COLS = ["mean_elev_m", "mean_scaling_stone"] + \
             [f"mean_scaling_lsdn_{n}_t{a}" for n in NUCLIDES for a in SAMPLE_AGES_YR]
CHECK_EXACT_COLS = ["n_upslope_pixels"]


def read_frozen(path: Path) -> dict:
    if not path.is_file():
        return {}
    with open(path, newline="") as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    return {r["basin"]: r for r in csv.DictReader(lines) if r.get("basin")}


def check_against_frozen(rows: list, out_csv: Path, datasets: list) -> int:
    """Compare a fresh run against the tracked artefact WITHOUT modifying it.

    The '# Run: ... ; git <sha>' provenance line is deliberately NOT compared and
    NOT rewritten: rewriting it dirties a tracked file on every suite run.
    """
    frozen = read_frozen(out_csv)
    if not frozen:
        print(f"CHECK: no frozen artefact at {out_csv}; run once with --freeze to "
              "create a reference first.")
        return EXIT_FAIL
    ok = True
    for r in rows:
        name = r["basin"]
        ref = frozen.get(name)
        if ref is None:
            print(f"[CHECK][FAIL] {name}: not in the frozen artefact")
            ok = False
            continue
        for col in CHECK_COLS:
            try:
                new_v, ref_v = float(r.get(col)), float(ref.get(col))
            except (TypeError, ValueError):
                print(f"[CHECK][FAIL] {name}.{col}: non-numeric "
                      f"(new={r.get(col)!r} frozen={ref.get(col)!r})")
                ok = False
                continue
            denom = abs(ref_v)
            rel = abs(new_v - ref_v) / denom if denom else abs(new_v - ref_v)
            if rel > CHECK_REL_TOL:
                print(f"[CHECK][FAIL] {name}.{col}: new={new_v:.8g} "
                      f"frozen={ref_v:.8g} rel={rel:.3g}")
                ok = False
        for col in CHECK_EXACT_COLS:
            if str(r.get(col)) != str(ref.get(col)):
                print(f"[CHECK][FAIL] {name}.{col}: new={r.get(col)} "
                      f"frozen={ref.get(col)} (must match exactly)")
                ok = False
    # A basin that the frozen artefact has but this run did not produce is a real
    # failure -- silently dropping one would let a delineation regression pass.
    # Restricted to the datasets ACTUALLY RUN, so a single-dataset check (useful
    # when triaging one landscape) is not swamped by the other two.
    ran = set(datasets)
    expected = {b for b, r in frozen.items() if r.get("dataset") in ran}
    missing = expected - {r["basin"] for r in rows}
    if missing:
        print(f"[CHECK][FAIL] frozen basins absent from this run: "
              f"{', '.join(sorted(missing))}")
        ok = False
    print(f"CHECK: {'PASS' if ok else 'FAIL'} (rel tol {CHECK_REL_TOL:g}, "
          f"{len(rows)} basins, artefact NOT modified)")
    if not ok:
        print("  If the change is intended, re-run with --freeze and review the diff.")
    return EXIT_OK if ok else EXIT_FAIL


def freeze(rows: list, out_csv: Path) -> int:
    sha = rcv.git_sha()
    with open(out_csv, "w", newline="") as f:
        f.write("# L13 PER-PIXEL side -- FROZEN artefact (lsdn plan Phase 6(ii); cosmo "
                "plan 6f row L13).\n")
        f.write("# Generated by tests/scripts/run_cosmo_lsdn_pixel_vs_lumped.py "
                "--freeze under the cosmo_realworld\n# tier (needs network: keyless "
                "COP30 fetch per dataset).  NOT fabricated.  The DEMs live in a temp\n"
                "# workdir and are deleted; only this few-kB table is tracked.\n")
        f.write("# mean_scaling_* = basin-mean PRODUCTION SCALING from lsdtt-cosmo-tool's "
                "_BASIN_PROD.csv, i.e.\n# the mean of calculate_production_raster over the "
                "EXACT upslope node set the inversion uses.\n# NO P0 multiply: it is a "
                "P0-normalised scaling, in the same units as OCTOPUS BEPROD.\n")
        f.write("# mean_scaling_stone is the quantity OCTOPUS iterates ATM_PRESS to "
                "reproduce, so it drives the\n# ATM_PRESS-recovery half of L13.  The lsdn "
                "columns are lsdtt3's PRODUCT -- the per-pixel LSDn\n# basin mean -- at "
                "three sample ages; the C++ driver evaluates the LUMPED LSDn at the SAME\n"
                "# age, so the per-pixel-vs-lumped ratio is always taken with t matched.\n")
        f.write("# Config: muon_scheme=Braucher atmosphere=NCEP standardisation=07KNSTD, "
                "per-manifest\n# threshold_contributing_pixels/snap_search_radius_pixels "
                "(matched to OCTOPUS, as 6e).\n")
        f.write("# THIS FILE IS FROZEN AND TRACKED.  The default mode of the generator is "
                "--check; regenerate\n# only with --freeze, and review the diff: a numeric "
                "change here is a change in the per-pixel\n# production result, not "
                "bookkeeping.\n")
        f.write(f"# Run: datasets={','.join(DATASETS)} ages={SAMPLE_AGES_YR} ; git {sha}\n")
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})
    print(f"\nREFROZE artefact (--freeze): {out_csv}")
    print("  This file is tracked.  Review the diff before committing.")
    return EXIT_OK


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dataset", nargs="?", default="all",
                    help="dataset name or 'all' (default)")
    ap.add_argument("--exe", required=True, help="path to lsdtt-cosmo-tool")
    ap.add_argument("--fetch-exe", required=True, help="path to lsdtt-fetch-raster")
    ap.add_argument("--published-dir")
    ap.add_argument("--keep", action="store_true", help="keep the workdir(s)")
    ap.add_argument("--freeze", action="store_true",
                    help="REGENERATE the tracked frozen artefact instead of checking "
                         "against it.  Deliberate, reviewed act; default is check-only.")
    args = ap.parse_args()

    datasets = DATASETS if args.dataset == "all" else [args.dataset]
    # A partial re-freeze would silently drop the datasets not run, so refuse it.
    if args.freeze and set(datasets) != set(DATASETS):
        print("FAIL: --freeze regenerates the whole artefact; run it with 'all'.")
        return EXIT_FAIL

    rows = []
    try:
        for d in datasets:
            rows.extend(rows_for_dataset(d, args))
    except rcv.SkipError as exc:
        print(f"SKIP: {exc}")
        return EXIT_SKIP
    if not rows:
        print("FAIL: no basins produced")
        return EXIT_FAIL

    out_csv = expected_cosmo_dir() / ARTEFACT_NAME
    if args.freeze:
        return freeze(rows, out_csv)
    return check_against_frozen(rows, out_csv, datasets)


if __name__ == "__main__":
    sys.exit(main())
