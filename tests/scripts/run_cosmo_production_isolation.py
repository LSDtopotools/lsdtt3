#!/usr/bin/env python3
"""Production-vs-BEPROD isolation test (cosmo plan Phase 6 / §6e).

DECISIVE test for whether the ~+4% lsdtt3-vs-OCTOPUS erosion-rate residual on the
Dethier (2014) Front Range basins is GEOMETRIC (COP30 DEM hypsometry / basin-mean
production) rather than tabular. At fixed measured concentration the inferred
erosion rate is proportional to basin-mean production, so if lsdtt3's basin-mean
production is ~+4% above OCTOPUS's stored CAIRN basin-mean production (BEPROD), the
whole residual is production/DEM and the inversion is exonerated.

Procedure (plan §6e):
  1. Fetch the keyless COP30 DEM (lsdtt-fetch-raster) into a temp workdir, warped to
     the manifest target_epsg -- deleted on exit (cleanup trap).
  2. Run lsdtt-cosmo-tool MATCHED to OCTOPUS: muon_scheme Braucher, atmosphere NCEP,
     sample_age_yr 0 (paleomag off), threshold_contributing_pixels 250,
     snap_search_radius_pixels 20, write_production_raster + write_basin_production_csv.
  3. The tool's _BASIN_PROD.csv reports, per snapped sample, the mean Stone-2000
     production SCALING and mean elevation over the EXACT upslope pixel set the
     inversion averages. P_lsdtt = P0_ref x mean_scaling, P0_ref = Braucher total
     SLHL 10Be production 4.30 atoms/g/yr (set_braucher_parameters).
  4. Compare per basin: (a) P_lsdtt vs BEPROD; (b) P_lsdtt x BETOPO vs BETOTS;
     (c) lsdtt mean elevation vs ELEV_AVE, against the frozen artefact.

FROZEN-ARTEFACT DISCIPLINE (fixed 2026-07-26). dethier2014_production_isolation.csv
is TRACKED and described as frozen, so the DEFAULT mode is now --check: compare and
report, never write. Previously every run rewrote the file -- in practice only its
'# Run: ... ; git <sha>' provenance line moved, so a clean tree went dirty purely by
running the suite, and the numbers were never actually regression-checked despite
the file being called frozen. Pass --freeze to regenerate deliberately, then review
the diff: a numeric change there is a change in the production/DEM result.

This is NOT in CI (needs network + time). Exit 77 = SKIP (network down), like the
rest of the real-landscape tier. --harvest pulls BEPROD/BETOTS/BETOPO/ELEV_AVE from
the OCTOPUS WFS (STUDYID S169); --run does the DEM fetch + tool run + compare.

Usage:
  # regression check (what CTest runs; does not modify the artefact)
  python3 tests/scripts/run_cosmo_production_isolation.py \
      --exe build/lsdtt-cosmo-tool --fetch-exe build/lsdtt-fetch-raster [--keep]

  # deliberate re-freeze (review the diff before committing)
  python3 tests/scripts/run_cosmo_production_isolation.py \
      --exe build/lsdtt-cosmo-tool --fetch-exe build/lsdtt-fetch-raster --freeze
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_cosmo_validation as rcv  # reuse fetch/cleanup plumbing

EXIT_OK, EXIT_FAIL, EXIT_SKIP = 0, 1, 77

# P0_ref: Braucher total SLHL 10Be production (atoms/g/yr), read EXACTLY from
# set_braucher_parameters() in src/cpu/lsdtt_crn_parameters.cpp (P0_10Be_ = 4.30).
# Reported for provenance only -- see note below on why it is NOT applied.
P0_REF_BRAUCHER_10BE = 4.30

# UNITS NOTE (tester, 2026-07-24): the plan §6e specifies P_lsdtt = P0_ref x
# mean(scaling). Empirically that OVER-scales by ~4.3x (= P0) and gives an absurd
# +330% "residual". lsdtt-cosmo-tool's PROD raster stores stone2000sp(lat,P,Fsp=1),
# i.e. the Stone-2000 SPALLATION scaling S (~6.3 at 2500 m), and its basin mean
# equals OCTOPUS BEPROD to ~1% DIRECTLY (6.205 vs 6.221; 11.16 vs 11.07; ...).
# OCTOPUS BEPROD is therefore a P0-normalised production SCALING, in the same units
# as lsdtt's PROD raster -- so the comparable quantity is mean(scaling) itself, with
# NO P0 multiply. We report P_lsdtt = mean(scaling) and flag the plan's units error.

PILOT_OBSIDS = [
    "S169WTS036", "S169WTS035", "S169WTS032", "S169WTS029", "S169WTS005",
    "S169WTS030", "S169WTS031", "S169WTS033", "S169WTS001", "S169WTS020",
    "S169WTS021", "S169WTS022", "S169WTS014", "S169WTS024",
]

WFS = "http://geoserver.octopusdata.org/geoserver/wfs"
BASINS_TYPE = "be10-denude:crn_int_basins"


def published_dir() -> Path:
    return rcv.repo_root() / "tests" / "expected" / "cosmo" / "published"


# --------------------------------------------------------------------------- #
# 1. Harvest OCTOPUS BEPROD/BETOTS/BETOPO/ELEV_AVE (never fabricated)
# --------------------------------------------------------------------------- #
def harvest_octopus(out_path: Path, studyid: str = "S169") -> int:
    params = {
        "service": "WFS", "version": "2.0.0", "request": "GetFeature",
        "typeName": BASINS_TYPE, "outputFormat": "application/json",
        "srsName": "EPSG:4326", "CQL_FILTER": f"STUDYID='{studyid}'",
    }
    url = WFS + "?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=180) as r:
            feats = json.load(r)["features"]
    except Exception as exc:  # noqa: BLE001
        print(f"HARVEST SKIP: OCTOPUS WFS unreachable ({exc})")
        return EXIT_SKIP
    by_oid = {f["properties"]["OBSID1"]: f["properties"] for f in feats}
    missing = [o for o in PILOT_OBSIDS if o not in by_oid]
    if missing:
        print(f"HARVEST FAIL: OBSIDs absent from S169: {missing}")
        return EXIT_FAIL
    with open(out_path, "w", newline="") as f:
        f.write("# OCTOPUS BEPROD/BETOTS harvest for the §6e production isolation "
                "test (STUDYID S169).\n")
        f.write("# HARVESTED from the OCTOPUS WFS (be10-denude:crn_int_basins), NOT "
                "fabricated. BEPROD = basin-mean\n# production before shielding "
                "(atoms/g/yr); BETOTS = BEPROD*BETOPO*BESELF*BESNOW (after shielding).\n")
        w = csv.writer(f)
        w.writerow(["octopus_obsid", "BEPROD", "BETOTS", "BETOPO",
                    "BESELF", "BESNOW", "ELEV_AVE"])
        for oid in PILOT_OBSIDS:
            p = by_oid[oid]
            w.writerow([oid, p["BEPROD"], p["BETOTS"], p["BETOPO"],
                        p["BESELF"], p["BESNOW"], p["ELEV_AVE"]])
    print(f"harvested {len(PILOT_OBSIDS)} basins -> {out_path}")
    return EXIT_OK


def read_octopus(path: Path) -> dict:
    with open(path, newline="") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    return {r["octopus_obsid"]: r for r in csv.DictReader(lines)}


# --------------------------------------------------------------------------- #
# 2/3. Run lsdtt-cosmo-tool matched to OCTOPUS, read the basin-prod CSV
# --------------------------------------------------------------------------- #
def run_tool(exe: Path, dem_path: Path, samples: list, workdir: Path) -> Path:
    points_csv = rcv.write_crn_input(samples, workdir)
    write_base = "isolation"
    param = workdir / "isolation.param"
    lines = [
        f"read_path: {workdir}",
        f"read_fname: {dem_path.name}",
        f"write_path: {workdir}",
        f"write_fname: {write_base}",
        "nuclide: Be10",
        f"points_filename: {points_csv.name}",
        # --- matched to OCTOPUS (plan §6e) ---
        "muon_scheme: Braucher",
        "atmosphere: NCEP",
        "sample_age_yr: 0",
        "threshold_contributing_pixels: 250",
        "snap_search_radius_pixels: 20",
        "write_production_raster: true",
        "write_basin_production_csv: true",
    ]
    param.write_text("\n".join(lines) + "\n")
    subprocess.run([str(exe), str(workdir), param.name], check=True, cwd=str(workdir))
    bp = workdir / f"{write_base}_BASIN_PROD.csv"
    if not bp.is_file():
        raise RuntimeError(f"cosmo-tool produced no basin-prod CSV: {bp}")
    return bp


def read_basin_prod(path: Path) -> dict:
    with open(path, newline="") as f:
        return {r["sample_name"]: r for r in csv.DictReader(f)}


# --------------------------------------------------------------------------- #
# 4/5. Compare + emit frozen artefact
# --------------------------------------------------------------------------- #
def run_isolation(args) -> int:
    pdir = published_dir()
    manifest = rcv.load_manifest("dethier2014", pdir)
    samples = rcv.read_sample_table(pdir / manifest["sample_table"])

    octo_path = pdir / "dethier2014_octopus_production.csv"
    if args.harvest or not octo_path.is_file():
        rc = harvest_octopus(octo_path)
        if rc != EXIT_OK:
            return rc
    octo = read_octopus(octo_path)

    workdir = Path(tempfile.mkdtemp(prefix="cosmo_iso_dethier2014_"))
    print(f"working directory: {workdir}")
    try:
        dem_path = rcv.fetch_dem(Path(args.fetch_exe), manifest, workdir)
        print(f"DEM: {dem_path} ({dem_path.stat().st_size} bytes)")
        bp_csv = run_tool(Path(args.exe), dem_path, samples, workdir)
        bp = read_basin_prod(bp_csv)
        erate_csv = bp_csv.parent / "isolation_CN_ERATE.csv"
        erate = {}
        if erate_csv.is_file():
            with open(erate_csv, newline="") as f:
                erate = {r["sample_name"]: r for r in csv.DictReader(f)}
    except rcv.SkipError as exc:
        print(f"SKIP: {exc}")
        return EXIT_SKIP
    finally:
        if args.keep:
            print(f"--keep set: leaving {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)
            print(f"cleaned up {workdir} (exists={workdir.exists()})")

    # min reliable basin area for the snap threshold footprint (250 px @ 30 m).
    min_reliable_km2 = 250 * 30.0 * 30.0 / 1.0e6
    samp_by_name = {s["sample_name"]: s for s in samples}

    rows, a_pcts, b_pcts, d_elevs, e_pcts = [], [], [], [], []
    for oid in PILOT_OBSIDS:
        if oid not in bp:
            print(f"WARNING: {oid} absent from tool output (snap/extent) -- skipped")
            continue
        r = bp[oid]
        o = octo[oid]
        mean_scal = float(r["mean_prod_scaling"])
        mean_elev = float(r["mean_elev_m"])
        # P_lsdtt = basin-mean Stone-2000 production scaling (SAME units as OCTOPUS
        # BEPROD; NO P0 multiply -- see UNITS NOTE at top of file).
        p_lsdtt = mean_scal
        beprod = float(o["BEPROD"])
        betots = float(o["BETOTS"])
        betopo = float(o["BETOPO"])
        elev_ave = float(o["ELEV_AVE"])
        p_lsdtt_x_topo = p_lsdtt * betopo
        pct_a = 100.0 * (p_lsdtt - beprod) / beprod
        pct_b = 100.0 * (p_lsdtt_x_topo - betots) / betots
        d_elev = mean_elev - elev_ave

        # erate residual (lsdtt vs OCTOPUS EBE_GCMYR) for the same basin -- lets us
        # say whether production matching coincides with erate matching.
        pct_e = ""
        er = erate.get(oid)
        sm = samp_by_name.get(oid)
        if er is not None and sm is not None and sm.get("published_erate_gcm2yr"):
            le = float(er["apparent_eff_e_g_cm2_yr"])
            pe = float(sm["published_erate_gcm2yr"])
            if pe:
                pct_e = 100.0 * (le - pe) / pe

        # delineation reliability: OCTOPUS area vs the snap-threshold footprint.
        try:
            area = float(sm.get("octopus_area_km2", "")) if sm else float("nan")
            reliable = area >= min_reliable_km2
        except (TypeError, ValueError):
            reliable = False

        if reliable:
            a_pcts.append(pct_a)
            b_pcts.append(pct_b)
            d_elevs.append(d_elev)
            if pct_e != "":
                e_pcts.append(pct_e)
        rows.append({
            "basin": oid,
            "n_upslope_pixels": r["n_upslope_pixels"],
            "delineation_reliable": "yes" if reliable else "no",
            "P_lsdtt": f"{p_lsdtt:.5f}",
            "BEPROD": f"{beprod:.5f}",
            "pct_a": f"{pct_a:.3f}",
            "P_lsdtt_x_BETOPO": f"{p_lsdtt_x_topo:.5f}",
            "BETOTS": f"{betots:.5f}",
            "pct_b": f"{pct_b:.3f}",
            "lsdtt_mean_elev_m": f"{mean_elev:.2f}",
            "ELEV_AVE": f"{elev_ave:.2f}",
            "d_elev_m": f"{d_elev:.2f}",
            "erate_pct": f"{pct_e:.3f}" if pct_e != "" else "",
        })
        print(f"[{oid}] {'OK ' if reliable else 'UNREL'} "
              f"P_lsdtt={p_lsdtt:.3f} BEPROD={beprod:.3f} pct_a={pct_a:+.2f}%  "
              f"d_elev={d_elev:+.1f} m  erate_pct="
              f"{(f'{pct_e:+.2f}%' if pct_e != '' else 'NA')}  n={r['n_upslope_pixels']}")

    if not rows:
        print("FAIL: no basins compared")
        return EXIT_FAIL

    def mean(xs):
        return sum(xs) / len(xs) if xs else float("nan")

    print("\n=== MEANS over {} DELINEATION-RELIABLE basins ===".format(len(a_pcts)))
    print(f"  (a) mean pct P_lsdtt vs BEPROD        : {mean(a_pcts):+.2f}%")
    print(f"  (b) mean pct P_lsdtt*BETOPO vs BETOTS : {mean(b_pcts):+.2f}%")
    print(f"  (c) mean d_elev (lsdtt - ELEV_AVE)    : {mean(d_elevs):+.1f} m")
    print(f"  (+) mean pct lsdtt erate vs OCTOPUS   : {mean(e_pcts):+.2f}% "
          f"(n={len(e_pcts)})")

    out_csv = pdir / "dethier2014_production_isolation.csv"
    fields = ["basin", "delineation_reliable", "P_lsdtt", "BEPROD", "pct_a",
              "P_lsdtt_x_BETOPO", "BETOTS", "pct_b", "lsdtt_mean_elev_m",
              "ELEV_AVE", "d_elev_m", "erate_pct", "n_upslope_pixels"]

    # DEFAULT IS --check: compare against the tracked artefact WITHOUT touching it.
    # This file is described as FROZEN and is tracked by git, so a plain suite run
    # must not rewrite it.  Regenerating is an explicit, reviewed act (--freeze),
    # exactly as for the sibling summaries in run_cosmo_validation.py.
    if not getattr(args, "freeze", False):
        return check_against_frozen_isolation(rows, out_csv)

    sha = rcv.git_sha()
    with open(out_csv, "w", newline="") as f:
        f.write("# Dethier (2014) §6e production-vs-BEPROD isolation -- FROZEN artefact.\n")
        f.write("# P_lsdtt = basin-mean Stone-2000 production SCALING from "
                "lsdtt-cosmo-tool _BASIN_PROD.csv\n# (stone2000sp, Fsp=1), averaged over "
                "the EXACT upslope set the inversion uses. NO P0 multiply:\n# OCTOPUS "
                "BEPROD is itself a P0-normalised scaling in the same units (lsdtt "
                "mean-S == BEPROD\n# to ~1%); the plan's 'x P0_ref' over-scales by ~4.3x "
                "and is a units error (tester 2026-07-24).\n")
        f.write("# BEPROD/BETOTS/BETOPO/ELEV_AVE harvested from OCTOPUS WFS (STUDYID "
                "S169), NOT fabricated.\n")
        f.write("# pct_a = 100*(P_lsdtt-BEPROD)/BEPROD ; pct_b = "
                "100*(P_lsdtt*BETOPO-BETOTS)/BETOTS ; d_elev_m = lsdtt_mean - ELEV_AVE ;\n"
                "# erate_pct = 100*(lsdtt EBE - OCTOPUS EBE_GCMYR)/EBE_GCMYR (inversion "
                "residual, same basin).\n")
        f.write(f"# Run: muon_scheme=Braucher atmosphere=NCEP sample_age_yr=0 "
                f"threshold_px=250 snap_px=20 ; git {sha}\n")
        f.write(f"# MEANS over {len(a_pcts)} delineation-reliable basins: "
                f"pct_a={mean(a_pcts):+.2f}% pct_b={mean(b_pcts):+.2f}% "
                f"d_elev={mean(d_elevs):+.1f} m erate_pct={mean(e_pcts):+.2f}%.\n")
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    print(f"\nREFROZE artefact (--freeze): {out_csv}")
    print("  This file is tracked. Review the diff before committing: a numeric "
          "change here is a change in the production/DEM result, not bookkeeping.")
    return EXIT_OK


# Relative tolerance for the frozen-artefact regression check.  The quantities are
# deterministic given the same DEM and executable, so this is tight on purpose; it
# exists to absorb COP30 tile-server and float-formatting noise, not real drift.
ISOLATION_CHECK_REL_TOL = 1e-3

# Columns compared numerically by --check.  n_upslope_pixels is an integer pixel
# count and is compared exactly: if the delineation moves, the test must say so
# rather than average it away.
ISOLATION_CHECK_COLS = ["P_lsdtt", "BEPROD", "pct_a", "P_lsdtt_x_BETOPO", "BETOTS",
                        "pct_b", "lsdtt_mean_elev_m", "ELEV_AVE", "d_elev_m",
                        "erate_pct"]
ISOLATION_CHECK_EXACT_COLS = ["n_upslope_pixels"]


def read_frozen_isolation(path: Path) -> dict:
    """Read the tracked artefact, skipping its '#' comment header.  Keyed by basin."""
    if not path.is_file():
        return {}
    with open(path, newline="") as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    return {r["basin"]: r for r in csv.DictReader(lines) if r.get("basin")}


def check_against_frozen_isolation(rows: list, out_csv: Path) -> int:
    """Compare a fresh run against the tracked frozen artefact WITHOUT modifying it.

    The artefact carries a '# Run: ... ; git <sha>' provenance line.  That line is
    deliberately NOT compared and NOT rewritten here: it changes on every commit and
    rewriting it dirtied a tracked file on every suite run, which is the defect this
    mode fixes.  Only the numbers are regression-checked.
    """
    frozen = read_frozen_isolation(out_csv)
    if not frozen:
        print(f"CHECK: no frozen artefact at {out_csv}; run once with --freeze to "
              "create a reference first.")
        return EXIT_FAIL

    ok = True
    for r in rows:
        name = r.get("basin")
        ref = frozen.get(name)
        if ref is None:
            print(f"[CHECK][FAIL] {name}: not in the frozen artefact")
            ok = False
            continue
        for col in ISOLATION_CHECK_COLS:
            try:
                new_v, ref_v = float(r.get(col)), float(ref.get(col))
            except (TypeError, ValueError):
                print(f"[CHECK][FAIL] {name}.{col}: non-numeric "
                      f"(new={r.get(col)!r} frozen={ref.get(col)!r})")
                ok = False
                continue
            denom = abs(ref_v)
            rel = abs(new_v - ref_v) / denom if denom else abs(new_v - ref_v)
            if rel > ISOLATION_CHECK_REL_TOL:
                print(f"[CHECK][FAIL] {name}.{col}: new={new_v:.6g} "
                      f"frozen={ref_v:.6g} rel={rel:.3g}")
                ok = False
        for col in ISOLATION_CHECK_EXACT_COLS:
            if str(r.get(col)) != str(ref.get(col)):
                print(f"[CHECK][FAIL] {name}.{col}: new={r.get(col)} "
                      f"frozen={ref.get(col)} (must match exactly)")
                ok = False

    missing = set(frozen) - {r.get("basin") for r in rows}
    if missing:
        print(f"[CHECK][FAIL] frozen basins absent from this run: "
              f"{', '.join(sorted(missing))}")
        ok = False

    print(f"CHECK: {'PASS' if ok else 'FAIL'} (rel tol {ISOLATION_CHECK_REL_TOL:g}, "
          f"{len(rows)} basins, artefact NOT modified)")
    if not ok:
        print("  If the change is intended, re-run with --freeze and review the diff.")
    return EXIT_OK if ok else EXIT_FAIL


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", help="path to lsdtt-cosmo-tool")
    ap.add_argument("--fetch-exe", help="path to lsdtt-fetch-raster")
    ap.add_argument("--harvest", action="store_true",
                    help="re-harvest OCTOPUS BEPROD/BETOTS even if cached")
    ap.add_argument("--harvest-only", action="store_true",
                    help="only harvest the OCTOPUS fields, do not run the tool")
    ap.add_argument("--keep", action="store_true", help="keep the workdir")
    ap.add_argument("--freeze", action="store_true",
                    help="REGENERATE the tracked frozen artefact instead of "
                         "checking against it. Deliberate, reviewed act: the "
                         "artefact is tracked by git and a numeric change in it is "
                         "a change in the result. Default is check-only.")
    args = ap.parse_args()

    if args.harvest_only:
        return harvest_octopus(published_dir() / "dethier2014_octopus_production.csv")
    if not args.exe or not args.fetch_exe:
        ap.error("--exe and --fetch-exe are required (unless --harvest-only)")
    try:
        return run_isolation(args)
    except rcv.SkipError as exc:
        print(f"SKIP: {exc}")
        return EXIT_SKIP


if __name__ == "__main__":
    sys.exit(main())
