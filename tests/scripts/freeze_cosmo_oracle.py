#!/usr/bin/env python3
"""Freeze and check the cosmogenic CLIs' outputs -- the refactor oracle.

WHY THIS EXISTS
---------------
The cosmogenic CLI family is about to be refactored: the production, pressure
and shielding setup is CLI-local and duplicated, and it is being extracted into
a shared component. **The oracle for a refactor is the pre-refactor binary, and
that binary stops existing the moment the refactor lands.** So it has to be
captured first, deliberately, as a tracked artefact.

Every phase of that work is then one of exactly two kinds, and must say which:

  behaviour-PRESERVING   the frozen set must come back UNCHANGED. That is the
                         proof the refactor was safe.
  behaviour-CHANGING     the frozen set is re-frozen with --freeze, and the
                         commit message states which artefacts moved and why.

A phase that is unclear about which kind it is has not been thought through.

WHAT IS HASHED, AND WHY IT IS EXACT
-----------------------------------
Rasters are converted to ENVI (raw, uncompressed, no header in the data file)
with gdal_translate and the raw bytes are sha256'd. That is **bitwise** on the
pixel values -- verified round-tripping through numpy against the source band.

Two rejected alternatives, recorded so they are not retried:

  * `gdalinfo -checksum` alone. GDAL's band checksum is a 16-bit value, and an
    all-zero band (a landslide excavation raster, say) checksums to 0 -- so it
    cannot distinguish "all zeros" from "all zeros in a different shape".
  * `gdalinfo -stats`. It truncates STATISTICS_MEAN and _MAXIMUM to two decimal
    places, which is nowhere near enough to support a bitwise claim, and it
    writes a .aux.xml sidecar that then caches stale statistics.

CSV outputs are hashed after stripping absolute paths, so a run from a different
directory compares equal.

**stdout is deliberately NOT hashed.** It carries absolute paths and free text
and is not a stable contract; hashing it would produce failures that mean
nothing. Everything stdout reports that actually matters -- derived K, the
selected atmosphere, the scaling scheme -- is expressed in the rasters, so a
change there moves a hash. Exit codes ARE recorded, because a CLI that starts
failing (or starts succeeding) is a real change.

USAGE
-----
  # check against the tracked oracle (the default, and what CTest runs)
  freeze_cosmo_oracle.py --exe lsdtt-cosmo-prediction=build/lsdtt-cosmo-prediction \\
      ... --read-path tests/data --work-dir /tmp/oracle --oracle tests/data/cosmo_oracle.csv

  # or, for a single-config build tree, let it find them
  freeze_cosmo_oracle.py --build-dir build ...

  # (re)generate the tracked oracle -- only ever deliberately
  freeze_cosmo_oracle.py ... --freeze

  # narrow while developing
  freeze_cosmo_oracle.py ... --only P3 T2

A missing oracle file in check mode is a **failure**, never a skip: the file is
tracked, so its absence means a broken checkout rather than an optional step.
That is the same rule the cosmo_realworld tier already follows.
"""

import argparse
import csv
import hashlib
import pathlib
import re
import shutil
import subprocess
import sys
import lsdtt_test_env


# ---------------------------------------------------------------------------
# The matrix.
#
# tier "cheap" runs in the default suite. tier "heavy" needs the LSDn scaling
# grids and/or the CRONUS emulator schemes, which are slower and depend on
# crn_parameters.nc carrying the LSDn grids.
#
# Each case is (id, exe, tier, param_body). read_path/write_path/write_fname
# are filled in per run. Keep write_fname unique per case -- outputs are
# collected by scanning the work directory.
# ---------------------------------------------------------------------------

PREDICTION = "lsdtt-cosmo-prediction"
TOOL       = "lsdtt-cosmo-tool"
LANDSLIDES = "lsdtt-cosmo-landslides"      # renamed to lsdtt-cosmo-landslides in R2a
COLUMN     = "lsdtt-cosmo-column-tool"

CASES = [

    # --- lsdtt-cosmo-prediction ------------------------------------------
    ("P1", PREDICTION, "cheap", """
calculate_crn_concentration_raster: true
effective_erosion_rate: 0.002
nuclide: Be10
muon_scheme: BraucherBorchers
scaling_scheme: stone
write_crn_concentration_raster: true
write_production_raster: true
"""),

    ("P2", PREDICTION, "cheap", """
calculate_step_change_crn_concentration_raster: true
step_change_eff_erate_old_g_cm2_yr: 0.002
step_change_eff_erate_new_g_cm2_yr: 0.02
time_since_change_yr: 8000
nuclide: Be10
scaling_scheme: stone
write_step_change_concentration_raster: true
"""),

    # Multi-isotope: this is the case that carries the RATIO rasters, which
    # stage SA1 promotes to a shared helper. It must be frozen before that.
    ("P3", PREDICTION, "cheap", """
calculate_multi_isotope_step_change: true
step_change_eff_erate_old_g_cm2_yr: 0.002
step_change_eff_erate_new_g_cm2_yr: 0.02
time_since_change_yr: 8000
scaling_scheme: stone
write_isotope_ratio_rasters: true
"""),

    # Non-unit shielding, separately from scaling. This is the combination a
    # scaling/shielding merge would break, and the codebase has already paid
    # 11.4 % for that mistake once.
    ("P4", PREDICTION, "cheap", """
calculate_crn_concentration_raster: true
effective_erosion_rate: 0.002
nuclide: Be10
scaling_scheme: stone
topo_shield_value: 0.85
snow_shield_eff_thickness: 12.0
self_shield_eff_thickness: 5.0
write_crn_concentration_raster: true
write_production_raster: true
"""),

    # Accumulation + sample points: the path goal 1 reuses.
    ("P5", PREDICTION, "cheap", """
calculate_crn_concentration_raster: true
calculate_accumulated_crn_concentration: true
effective_erosion_rate: 0.002
nuclide: Be10
scaling_scheme: stone
threshold_contributing_pixels: 2
sample_snap_max_search_radius_pixels: 3
sample_points_fname: cosmo_tiny_points.csv
write_crn_concentration_raster: true
write_accumulated_crn_raster: true
"""),

    ("P6", PREDICTION, "heavy", """
calculate_multi_isotope_step_change: true
step_change_eff_erate_old_g_cm2_yr: 0.002
step_change_eff_erate_new_g_cm2_yr: 0.02
time_since_change_yr: 8000
scaling_scheme: lsdn
sample_age_yr: 5000
write_isotope_ratio_rasters: true
"""),

    ("P7", PREDICTION, "heavy", """
calculate_crn_concentration_raster: true
effective_erosion_rate: 0.002
nuclide: Be10
muon_scheme: cronusEmulatorSt
scaling_scheme: stone
write_crn_concentration_raster: true
write_production_raster: true
"""),

    ("P8", PREDICTION, "heavy", """
calculate_crn_concentration_raster: true
effective_erosion_rate: 0.002
nuclide: Be10
muon_scheme: cronusEmulatorLSDn
scaling_scheme: lsdn
sample_age_yr: 5000
write_crn_concentration_raster: true
write_production_raster: true
"""),

    # --- lsdtt-cosmo-tool -------------------------------------------------
    ("T1", TOOL, "cheap", """
sample_points_fname: cosmo_tiny_points.csv
nuclide: Be10
muon_scheme: BraucherBorchers
scaling_scheme: stone
threshold_contributing_pixels: 2
sample_snap_max_search_radius_pixels: 3
write_production_raster: true
"""),

    # Non-unit shielding through the INVERSION, which reaches the scalar
    # shielding path this CLI has instead of shielding rasters.
    ("T2", TOOL, "cheap", """
sample_points_fname: cosmo_tiny_points.csv
nuclide: Be10
scaling_scheme: stone
topo_shield_value: 0.85
snow_shield_eff_thickness: 12.0
threshold_contributing_pixels: 2
sample_snap_max_search_radius_pixels: 3
"""),

    ("T3", TOOL, "heavy", """
sample_points_fname: cosmo_tiny_points.csv
nuclide: Be10
scaling_scheme: lsdn
threshold_contributing_pixels: 2
sample_snap_max_search_radius_pixels: 3
write_production_raster: true
"""),

    # THE HOT PATH OF E2. self_consistent_age re-evaluates production over one
    # basin's pixels per iterate, via the sc_eval_set overload. If the site
    # extraction forces whole-raster evaluation the VALUES stay identical and
    # only the cost changes -- so this case pins correctness, and the separate
    # performance guard pins the cost. Both are needed.
    ("T4", TOOL, "heavy", """
sample_points_fname: cosmo_tiny_points.csv
nuclide: Be10
scaling_scheme: lsdn
threshold_contributing_pixels: 2
sample_snap_max_search_radius_pixels: 3
self_consistent_age: true
"""),

    # --- lsdtt-cosmo-landslides (becomes lsdtt-cosmo-landslides at R2a) ---------
    # Its concentration output is EXPECTED to move at R1 (the scalar path was
    # never the reference) and again at R2c (real production and shielding).
    # Frozen here so that when it moves, the size of the move is known rather
    # than discovered.
    ("L1", LANDSLIDES, "cheap", """
landslide_mode: single_event
threshold_contributing_pixels: 20
event_volume_m3: 20
event_time_yr: 0
random_seed: 7
area_volume_preset: larsen_soil_scar
calculate_crn_concentration: true
background_eff_erate_g_cm2_yr: 0.002
production_scaling: 4.0
total_shielding: 1.0
nuclide: Be10
muon_scheme: BraucherBorchers
"""),

    # --- lsdtt-cosmo-column-tool -----------------------------------------
    ("C1", COLUMN, "cheap", """
column_mode: steady_state
nuclide: Be10
muon_scheme: BraucherBorchers
effective_erosion_rate_g_cm2_yr: 0.002
production_scaling: 4.0
total_shielding: 1.0
"""),

    ("C2", COLUMN, "cheap", """
column_mode: step_change
nuclide: Be10
muon_scheme: BraucherBorchers
effective_erosion_rate_old_g_cm2_yr: 0.002
effective_erosion_rate_g_cm2_yr: 0.02
time_since_step_yr: 8000
production_scaling: 4.0
"""),
]


RASTER_EXT = {".tif", ".tiff", ".bil", ".asc"}
TEXT_EXT   = {".csv"}
# README/help artefacts embed absolute paths and the parameter dump; they are
# not results and would fail for reasons that mean nothing.
SKIP_PATTERNS = (re.compile(r"_README\."), re.compile(r"\.aux\.xml$"),
                 re.compile(r"\.prj$"), re.compile(r"\.hdr$"))


# re.MULTILINE is load-bearing: without it `^` anchors to the start of the
# whole capture and findall returns nothing, which is why an empty parse is a
# hard error below rather than an empty "nothing to check" set.
_DECLARED_RE = re.compile(r"^ +([a-z_0-9]+) \((?:string|bool|float|int)",
                          re.MULTILINE)

_declared_cache = {}


def declared_parameters(exe: pathlib.Path) -> set:
    """Every parameter name a CLI declares, from its own -h output.

    Help is rendered from the declarations, so this is authoritative rather
    than a second list to keep in step.
    """
    key = str(exe)
    if key not in _declared_cache:
        r = lsdtt_test_env.run([str(exe), "-h"], capture_output=True, text=True,
                           check=False)
        names = set(_DECLARED_RE.findall(r.stdout))
        if not names:
            raise SystemExit(f"could not parse declared parameters from "
                             f"{exe} -h; the help format may have changed.")
        _declared_cache[key] = names
    return _declared_cache[key]


def check_case_parameters(case_id: str, exe: pathlib.Path, body: str) -> list:
    """Names in a case body that the CLI does not declare.

    AN UNRECOGNISED PARAMETER IS SILENTLY IGNORED by the parameter parser, so a
    case with a typo freezes the CLI's DEFAULT instead of what it meant to
    exercise -- an oracle that looks fine and guards nothing. This was not
    hypothetical: the first draft of case C1 wrote `effective_erosion_rate`
    where lsdtt-cosmo-column-tool declares `effective_erosion_rate_g_cm2_yr`,
    and froze a run at the default 0.1 g/cm2/yr instead of the intended 0.002.
    """
    declared = declared_parameters(exe)
    bad = []
    for line in body.strip().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or ":" not in line:
            continue
        name = line.split(":", 1)[0].strip()
        if name not in declared:
            bad.append(name)
    return bad


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def hash_raster(path: pathlib.Path, scratch: pathlib.Path) -> str:
    """Bitwise hash of a raster's pixel values, via a raw ENVI dump."""
    out = scratch / (path.name + ".envi.img")
    for stale in (out, out.with_suffix(".hdr"), pathlib.Path(str(out) + ".aux.xml")):
        if stale.exists():
            stale.unlink()
    r = lsdtt_test_env.run(
        ["gdal_translate", "-of", "ENVI", "-q", str(path), str(out)],
        capture_output=True, text=True, check=False)
    if r.returncode != 0:
        raise SystemExit(f"gdal_translate failed for {path}:\n{r.stderr}")
    return sha256_bytes(out.read_bytes())


def hash_text(path: pathlib.Path, work: pathlib.Path) -> str:
    """Hash a CSV with absolute paths stripped, so the run location is irrelevant."""
    txt = path.read_text(errors="replace")
    txt = txt.replace(str(work), "<WORK>")
    txt = re.sub(r"/[^\s,]*/(?=[^\s,/]+\.(?:tif|asc|csv))", "", txt)
    return sha256_bytes(txt.encode())


def collect_artefacts(work: pathlib.Path, scratch: pathlib.Path):
    """(name, kind, hash) for every result file in the work dir, sorted."""
    rows = []
    for p in sorted(work.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(work).as_posix()
        if any(pat.search(rel) for pat in SKIP_PATTERNS):
            continue
        ext = p.suffix.lower()
        if ext in RASTER_EXT:
            rows.append((rel, "raster", hash_raster(p, scratch)))
        elif ext in TEXT_EXT:
            rows.append((rel, "csv", hash_text(p, work)))
    return rows


def run_case(case, exe_map: dict, read_path: pathlib.Path,
             work_root: pathlib.Path):
    """Run one case; return (rows, exit_code) with rows as (artefact, kind, hash)."""
    case_id, exe_name, _tier, body = case
    exe = exe_map.get(exe_name)
    if exe is None or not exe.exists():
        raise SystemExit(
            f"{case_id}: executable not found for {exe_name}: {exe}\n"
            f"       Pass --exe {exe_name}=<path>, or --build-dir for a "
            f"single-config tree.")

    unknown = check_case_parameters(case_id, exe, body)
    if unknown:
        raise SystemExit(
            f"{case_id}: parameter(s) not declared by {exe_name}: "
            f"{', '.join(unknown)}\n"
            f"       The parser ignores unknown names silently, so this case "
            f"would freeze a DEFAULT\n"
            f"       rather than the behaviour it names. Fix the case.")

    work = work_root / case_id
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    scratch = work_root / "_scratch" / case_id
    scratch.mkdir(parents=True, exist_ok=True)

    param = work / f"{case_id}.param"
    param.write_text(
        f"read_path: {read_path}\n"
        f"read_fname: tiny_dem.asc\n"
        f"write_path: {work}\n"
        f"write_fname: {case_id.lower()}\n"
        f"write_output_format: gtiff\n"
        + body.strip() + "\n")

    proc = lsdtt_test_env.run([str(exe), param.name], cwd=str(work),
                          capture_output=True, text=True, check=False)
    (work / "stdout.log").write_text(proc.stdout)
    (work / "stderr.log").write_text(proc.stderr)

    rows = collect_artefacts(work, scratch)
    # The param file and the logs are inputs/diagnostics, not results.
    rows = [r for r in rows if not r[0].endswith(".param")]
    return rows, proc.returncode


def read_oracle(path: pathlib.Path):
    rec = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            rec.setdefault(row["case"], {})[row["artefact"]] = {
                "kind": row["kind"], "sha256": row["sha256"],
                "exit_code": row["exit_code"]}
    return rec


def write_oracle(path: pathlib.Path, results):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "exe", "tier", "artefact", "kind", "sha256",
                    "exit_code"])
        for case_id, exe, tier, rows, code in results:
            if not rows:
                w.writerow([case_id, exe, tier, "<no artefacts>", "none", "",
                            code])
            for artefact, kind, digest in rows:
                w.writerow([case_id, exe, tier, artefact, kind, digest, code])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build-dir", type=pathlib.Path, default=None,
                    help="single-config build tree; executables looked up by "
                         "name inside it. Prefer --exe from CMake, which is "
                         "generator-safe.")
    ap.add_argument("--exe", action="append", default=[], metavar="NAME=PATH",
                    help="explicit path for one executable; repeatable. "
                         "CMake passes $<TARGET_FILE:...> this way.")
    ap.add_argument("--read-path", type=pathlib.Path, required=True)
    ap.add_argument("--work-dir",  type=pathlib.Path, required=True)
    ap.add_argument("--oracle",    type=pathlib.Path, required=True)
    ap.add_argument("--tier", choices=["cheap", "all"], default="cheap")
    ap.add_argument("--freeze", action="store_true",
                    help="rewrite the tracked oracle. Deliberate acts only.")
    ap.add_argument("--only", nargs="*", default=None, metavar="CASE_ID")
    args = ap.parse_args(argv)

    exe_map = {}
    for spec in args.exe:
        if "=" not in spec:
            raise SystemExit(f"--exe expects NAME=PATH, got {spec!r}")
        name, path = spec.split("=", 1)
        exe_map[name] = pathlib.Path(path).resolve()
    if args.build_dir is not None:
        for name in (PREDICTION, TOOL, LANDSLIDES, COLUMN):
            exe_map.setdefault(name, (args.build_dir / name).resolve())
    if not exe_map:
        raise SystemExit("give --exe NAME=PATH (repeatable) or --build-dir")

    read_path = args.read_path.resolve()
    work_root = args.work_dir.resolve()
    oracle    = args.oracle.resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    if not shutil.which("gdal_translate"):
        raise SystemExit("gdal_translate not on PATH; it is required to hash "
                         "raster outputs exactly.")

    cases = [c for c in CASES if args.tier == "all" or c[2] == "cheap"]
    if args.only:
        wanted = set(args.only)
        cases = [c for c in cases if c[0] in wanted]
        missing = wanted - {c[0] for c in cases}
        if missing:
            raise SystemExit(f"--only named unknown/filtered cases: "
                             f"{sorted(missing)}")
    if not cases:
        raise SystemExit("no cases selected")

    print("=" * 68)
    print(" cosmogenic CLI oracle -- "
          + ("FREEZING" if args.freeze else "checking"))
    print("=" * 68)
    print(f"  executables : {len(exe_map)} declared")
    print(f"  oracle      : {oracle}")
    print(f"  tier        : {args.tier}  ({len(cases)} cases)")
    print()

    results = []
    for case in cases:
        rows, code = run_case(case, exe_map, read_path, work_root)
        results.append((case[0], case[1], case[2], rows, code))
        print(f"  {case[0]:<4} {case[1]:<26} exit={code}  "
              f"{len(rows)} artefact(s)")

    if args.freeze:
        write_oracle(oracle, results)
        print(f"\n  FROZE {sum(len(r[3]) for r in results)} artefacts across "
              f"{len(results)} cases -> {oracle}")
        print("  This is a deliberate act. State in the commit message which "
              "artefacts moved and why.")
        return 0

    # ---- check mode ----
    if not oracle.exists():
        print(f"\nFAIL: no oracle at {oracle}.\n"
              f"      The oracle is a TRACKED artefact, so its absence means a "
              f"broken checkout,\n"
              f"      not an optional step. Regenerate deliberately with "
              f"--freeze.", file=sys.stderr)
        return 1

    ref = read_oracle(oracle)
    n_fail = 0
    n_ok = 0
    for case_id, _exe, _tier, rows, code in results:
        if case_id not in ref:
            print(f"  FAIL {case_id}: not in the oracle. Freeze it "
                  f"deliberately if it is new.")
            n_fail += 1
            continue
        want = ref[case_id]
        got = {a: {"kind": k, "sha256": h} for a, k, h in rows}

        # Exit code first, and unconditionally: every row of a case carries the
        # same one, and a case with NO artefacts still has an exit code worth
        # checking -- a CLI that starts failing is the loudest change there is.
        want_code = next(iter(want.values()))["exit_code"]
        if str(code) != want_code:
            print(f"  FAIL {case_id}: exit code {code}, oracle says {want_code}")
            n_fail += 1

        for artefact, w in sorted(want.items()):
            if artefact == "<no artefacts>":
                continue
            g = got.get(artefact)
            if g is None:
                print(f"  FAIL {case_id}/{artefact}: not produced by this run")
                n_fail += 1
            elif g["sha256"] != w["sha256"]:
                print(f"  FAIL {case_id}/{artefact}: hash differs\n"
                      f"        oracle {w['sha256'][:16]}  "
                      f"got {g['sha256'][:16]}")
                n_fail += 1
            else:
                n_ok += 1

        extra = set(got) - set(want)
        if extra:
            print(f"  FAIL {case_id}: NEW artefacts not in the oracle: "
                  f"{sorted(extra)}")
            n_fail += len(extra)

    missing_cases = set(ref) - {r[0] for r in results}
    if missing_cases and not args.only and args.tier == "all":
        print(f"  FAIL: oracle cases absent from this run: "
              f"{sorted(missing_cases)}")
        n_fail += len(missing_cases)

    print()
    print(f"  matched {n_ok}, failed {n_fail}")
    if n_fail:
        print("\n  A failure here means a cosmogenic CLI's output MOVED.\n"
              "  If that was intended, re-freeze with --freeze and say so in "
              "the commit message.\n"
              "  If it was not, the change is a regression -- do not re-freeze "
              "to make it pass.", file=sys.stderr)
        return 1
    print("  all artefacts match the frozen oracle")
    return 0


if __name__ == "__main__":
    sys.exit(main())
