#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Generate the HYPE oracle files for the lsdtt3 HYPE-column tests.

WHAT THIS IS
------------
The lsdtt3 HYPE column is a reimplementation of the HYPE hydrological model's
process equations (SMHI, LGPLv3, https://sourceforge.net/projects/hype/). Its
CTest tests compare lsdtt3 against HYPE's own output, but CTest NEVER downloads
or compiles HYPE. Instead this script, run by hand, does that once and freezes
the results as small text files under tests/expected/hype/. The tests read only
those files.

Regenerating is a reviewable event, like re-fetching the CRONUS reference:
commit the regenerated files and say why in the commit message.

WHAT IT DOES
------------
1. Finds the latest HYPE release that ships source on SourceForge (or uses
   --hype-archive), downloads hype_X_Y_Z_src.tgz and checks its SHA-256 against
   the version pinned in tests/expected/hype/manifest.json. A different version
   or hash stops the script unless --accept-new-version is given.
2. Builds HYPE twice with gfortran in a temporary directory outside the repo:
     * "single": the shipped makefile unchanged (default REAL = single
       precision). This is the released program and is the ORACLE.
     * "double": the same flags plus -fdefault-real-8 -fdefault-double-8. Used
       ONLY to measure how far single-precision round-off moves each output;
       its values are not committed.
3. For every case definition in tests/expected/hype/cases/*.json, writes two
   views of the same set-up:
     HYPE input (work directory only; regenerated from the JSON every run):
                          info.txt GeoData.txt GeoClass.txt par.txt Pobs.txt Tobs.txt
     <case>/lsdtt_input/  general_parameters.csv landuse_parameters.csv
                          soil_parameters.csv classes.csv subbasins.csv
                          subbasin_classes.csv forcing.csv
   runs both builds on the HYPE input (in the temporary work directory, so
   HYPE's log files never land in the repo), and writes the single-precision
   results, at HYPE's 9 significant figures in shortest form, to
     <case>/expected/class_<var>.csv   columns: date, <subid>_<classid> ...
     <case>/expected/basin_<var>.csv   columns: date, <subid> ...
4. Writes tests/expected/hype/manifest.json (provenance and, per case and
   variable, the measured single-vs-double difference) and README.md.

--check regenerates everything into a temporary directory and compares it with
the committed files byte for byte, ignoring only the provenance fields that
legitimately change between runs (generation time, host, script revision).

Usage:
  python3 tests/scripts/generate_hype_oracle.py [--accept-new-version]
  python3 tests/scripts/generate_hype_oracle.py --check
  python3 tests/scripts/generate_hype_oracle.py --hype-archive /path/hype_5_36_0_src.tgz

Needs: python3 (standard library only), gfortran, make, network access unless
--hype-archive is given.
"""

import argparse
import csv
import datetime
import filecmp
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

import lsdtt_test_env  # non-interactive process launching (tests/scripts rule)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO_ROOT / "tests" / "expected" / "hype"
SF_FILES = "https://sourceforge.net/projects/hype/files/"
SF_GIT = "https://git.code.sf.net/p/hype/code"

# HYPE prints with ES format; above 9 significant figures the single build's
# format statement breaks (Fortran runtime error, ES16.;). 9 is the maximum.
SIGNIFICANT_FIGURES = 9

# Flags of the shipped makefile's gfortran release configuration, minus -Wall,
# plus the precision promotion. Used only for the "double" measurement build.
DOUBLE_SWITCH = ("-O2 -fbacktrace -ffree-line-length-0 -fno-automatic -std=gnu "
                 "-cpp -fdefault-real-8 -fdefault-double-8")

# Provenance fields excluded from --check.
VOLATILE_KEYS = {"generated_utc", "host", "script_git_revision"}


# ---------------------------------------------------------------------------
# HYPE source: discovery, download, build
# ---------------------------------------------------------------------------
def http_get(url, binary=False):
    # SourceForge serves its HTML "your download will start" page, not the file,
    # to browser-like User-Agents on /download links. A plain tool User-Agent
    # gets redirected to the mirror and receives the archive itself.
    ua = "lsdtt3-oracle/1.0" if binary else "Mozilla/5.0 lsdtt3-oracle"
    req = urllib.request.Request(url, headers={"User-Agent": ua})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    if binary and data[:2] != b"\x1f\x8b":
        sys.exit(f"ERROR: {url} did not return a gzip archive (got {data[:15]!r})")
    return data if binary else data.decode("utf-8", errors="replace")


def latest_release_with_source():
    """Highest release_hype_X_Y_Z folder on SourceForge holding hype_X_Y_Z_src.tgz."""
    html = http_get(SF_FILES)
    versions = sorted({tuple(int(x) for x in m) for m in
                       re.findall(r"release_hype_(\d+)_(\d+)_(\d+)", html)}, reverse=True)
    for v in versions:
        tag = "_".join(str(x) for x in v)
        folder = f"{SF_FILES}release_hype_{tag}/"
        if f"hype_{tag}_src.tgz" in http_get(folder):
            return ".".join(str(x) for x in v), (
                f"{folder}hype_{tag}_src.tgz/download")
    sys.exit("ERROR: no HYPE release with a source archive found on SourceForge")


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit_of_release(version):
    """Commit of the release in the SourceForge git repo, if master points at it.
    Informational only; returns '' when it cannot be determined."""
    try:
        out = lsdtt_test_env.run(["git", "ls-remote", SF_GIT, "refs/heads/master"],
                                 capture_output=True, text=True, timeout=60).stdout
        return out.split()[0] if out.strip() else ""
    except (OSError, subprocess.TimeoutExpired):
        return ""


def build_hype(archive, workdir, label, extra_switch=None):
    src_parent = workdir / f"build_{label}"
    src_parent.mkdir()
    with tarfile.open(archive) as tf:
        tf.extractall(src_parent)
    srcdirs = [p for p in src_parent.iterdir() if p.is_dir()]
    if len(srcdirs) != 1:
        sys.exit(f"ERROR: unexpected archive layout in {archive}")
    cmd = ["make", "hype"]
    if extra_switch:
        cmd.append(f"switch={extra_switch}")
    log = srcdirs[0] / "make.log"
    with open(log, "w") as f:
        rc = lsdtt_test_env.run(cmd, cwd=srcdirs[0], stdout=f, stderr=subprocess.STDOUT).returncode
    exe = srcdirs[0] / "hype"
    if rc != 0 or not exe.exists():
        sys.exit(f"ERROR: building HYPE ({label}) failed; see {log}")
    return exe


def gfortran_version():
    out = lsdtt_test_env.run(["gfortran", "--version"], capture_output=True, text=True).stdout
    return out.splitlines()[0].strip() if out else "unknown"


# ---------------------------------------------------------------------------
# Case definition -> forcing
# ---------------------------------------------------------------------------
def dates_of(case):
    start = datetime.date.fromisoformat(case["dates"]["start"])
    end = datetime.date.fromisoformat(case["dates"]["end"])
    n = (end - start).days + 1
    return [start + datetime.timedelta(days=i) for i in range(n)]


def forcing_of(case, dates):
    """Deterministic forcing, written to file with fixed formatting.

    precipitation: precip_pattern_mm repeated cyclically, times precip_scale.
    temperature:   temp_mean_c + temp_amp_c * sin(2 pi (step - temp_phase_steps) / temp_period_steps)
                   + temp_trend_c_per_day * step, rounded to 2 decimals.
    A short temp_period_steps puts a whole seasonal cycle into a short case,
    which keeps the committed files small.
    """
    f = case["forcing"]
    pat = f["precip_pattern_mm"]
    scale = f.get("precip_scale", 1.0)
    prec, temp = [], []
    for i, d in enumerate(dates):
        prec.append(round(pat[i % len(pat)] * scale, 3))
        t = (f["temp_mean_c"]
             + f.get("temp_amp_c", 0.0) * math.sin(
                 2 * math.pi * (i - f.get("temp_phase_steps", 0)) / f.get("temp_period_steps", 365))
             + f.get("temp_trend_c_per_day", 0.0) * i)
        temp.append(round(t, 2))
    return prec, temp


def fmt(x):
    """Stable text form for input files: integers stay integers."""
    if isinstance(x, int) or (isinstance(x, float) and x.is_integer() and abs(x) < 1e15):
        return str(int(x))
    return repr(float(x))


# ---------------------------------------------------------------------------
# Case definition -> HYPE input
# ---------------------------------------------------------------------------
def class_ids(case):
    return sorted(int(k) for k in case["classes"])


def write_hype_inputs(case, dates, prec, temp, d):
    d.mkdir(parents=True, exist_ok=True)
    subs = case["subbasins"]
    cids = class_ids(case)
    nl = (lambda s: s + "\n")

    # info.txt. Model options not listed are HYPE's defaults (0).
    lines = [
        f"bdate {dates[0].isoformat()}",
        f"cdate {dates[0].isoformat()}",
        f"edate {dates[-1].isoformat()}",
        "resultdir ./results/",
        "readobsid n",
    ]
    for opt, val in sorted(case.get("modeloptions", {}).items()):
        lines.append(f"modeloption {opt} {val}")
    if case.get("class_outputs"):
        lines += [
            "classoutput 1 variable " + " ".join(case["class_outputs"]),
            "classoutput 1 allbasin",
            "classoutput 1 definegroup allclass",
            "classoutput 1 meanperiod 1",
            f"classoutput 1 signfigures {SIGNIFICANT_FIGURES}",
        ]
    if case.get("basin_outputs"):
        lines += [
            "basinoutput variable " + " ".join(case["basin_outputs"]),
            "basinoutput allbasin",
            "basinoutput meanperiod 1",
            f"basinoutput signfigures {SIGNIFICANT_FIGURES}",
        ]
    (d / "info.txt").write_text("".join(nl(l) for l in lines))

    # GeoData.txt
    head = ["subid", "maindown", "area", "elev_mean", "slope_mean", "loc_rivlen", "rivlen"]
    head += [f"slc_{c}" for c in cids] + [f"dhslc_{c}" for c in cids]
    rows = ["\t".join(head)]
    for s in subs:
        cl = s["classes"]
        r = [s["subid"], s["maindown"], s["area_m2"], s["elev_mean_m"], s["slope_mean_percent"],
             s["loc_rivlen_m"], s["rivlen_m"]]
        r += [cl.get(str(c), {}).get("fraction", 0.0) for c in cids]
        r += [cl.get(str(c), {}).get("dh_m", 0.0) for c in cids]
        rows.append("\t".join(fmt(x) for x in r))
    (d / "GeoData.txt").write_text("".join(nl(r) for r in rows))

    # GeoClass.txt: class landuse soil crop1 crop2 rotation vegtype special
    #               tiledepth streamdepth nlayers depth1..3 (layer bottoms, m)
    rows = ["!class\tlanduse\tsoil\tcrop1\tcrop2\trotation\tvegtype\tspecial\t"
            "tiledepth\tstreamdepth\tnlayers\tdepth1\tdepth2\tdepth3"]
    for c in cids:
        k = case["classes"][str(c)]
        depths = k["layer_bottoms_m"]
        r = [c, k["landuse"], k["soil"], 0, 0, 0, k.get("vegtype", 1), 0,
             k.get("tile_depth_m", 0.0), k["stream_depth_m"], len(depths)] + depths
        rows.append("\t".join(fmt(x) for x in r))
    (d / "GeoClass.txt").write_text("".join(nl(r) for r in rows))

    # par.txt: general values, then per-land-use and per-soil vectors.
    rows = []
    for k, v in sorted(case["general"].items()):
        rows.append(f"{k}\t{fmt(v)}")
    for table in ("landuse", "soil"):
        codes = sorted(int(c) for c in case[table])
        n = max(codes)
        names = sorted({p for c in case[table].values() for p in c})
        for p in names:
            vals = [case[table].get(str(i), {}).get(p, 0.0) for i in range(1, n + 1)]
            rows.append(p + "\t" + "\t".join(fmt(v) for v in vals))
    (d / "par.txt").write_text("".join(nl(r) for r in rows))

    # Pobs.txt / Tobs.txt: one identical column per sub-basin (first version
    # applies one series to the whole domain).
    ids = [str(s["subid"]) for s in subs]
    for name, series in (("Pobs.txt", prec), ("Tobs.txt", temp)):
        rows = ["DATE\t" + "\t".join(ids)]
        for dt, v in zip(dates, series):
            rows.append(dt.isoformat() + "\t" + "\t".join(fmt(v) for _ in ids))
        (d / name).write_text("".join(nl(r) for r in rows))


# ---------------------------------------------------------------------------
# Case definition -> lsdtt3 input (the view the C++ tests read)
# ---------------------------------------------------------------------------
def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([fmt(x) if isinstance(x, (int, float)) else x for x in r])


def write_lsdtt_inputs(case, dates, prec, temp, d):
    d.mkdir(parents=True, exist_ok=True)
    write_csv(d / "general_parameters.csv", ["name", "value"],
              sorted(case["general"].items()))
    for table, fname in (("landuse", "landuse_parameters.csv"), ("soil", "soil_parameters.csv")):
        names = sorted({p for c in case[table].values() for p in c})
        rows = [[int(code)] + [case[table][code].get(p, 0.0) for p in names]
                for code in sorted(case[table], key=int)]
        write_csv(d / fname, ["code"] + names, rows)
    rows = []
    for c in class_ids(case):
        k = case["classes"][str(c)]
        depths = list(k["layer_bottoms_m"]) + [""] * (3 - len(k["layer_bottoms_m"]))
        rows.append([c, k["landuse"], k["soil"], len(k["layer_bottoms_m"])] + depths +
                    [k["stream_depth_m"], k.get("tile_depth_m", 0.0)])
    write_csv(d / "classes.csv",
              ["class_id", "landuse", "soil", "n_layers", "layer1_bottom_m", "layer2_bottom_m",
               "layer3_bottom_m", "stream_depth_m", "tile_depth_m"], rows)
    write_csv(d / "subbasins.csv",
              ["subid", "maindown", "area_m2", "elev_mean_m", "slope_mean_percent",
               "loc_rivlen_m", "rivlen_m"],
              [[s["subid"], s["maindown"], s["area_m2"], s["elev_mean_m"],
                s["slope_mean_percent"], s["loc_rivlen_m"], s["rivlen_m"]]
               for s in case["subbasins"]])
    rows = []
    for s in case["subbasins"]:
        for c in sorted(s["classes"], key=int):
            k = s["classes"][c]
            rows.append([s["subid"], int(c), k["fraction"], k.get("dh_m", 0.0)])
    write_csv(d / "subbasin_classes.csv", ["subid", "class_id", "fraction", "dh_m"], rows)
    write_csv(d / "forcing.csv", ["date", "precipitation_mm", "temperature_c"],
              [[dt.isoformat(), p, t] for dt, p, t in zip(dates, prec, temp)])


# ---------------------------------------------------------------------------
# Run HYPE, read its output
# ---------------------------------------------------------------------------
def run_hype(exe, hype_input, workdir, label):
    run_dir = workdir / label
    shutil.copytree(hype_input, run_dir)
    (run_dir / "results").mkdir()
    proc = lsdtt_test_env.run([str(exe), str(run_dir) + "/"], cwd=run_dir,
                              capture_output=True, text=True)
    # HYPE prints its final status line to stderr.
    ok = "Simulation will halt: successfully" in (proc.stdout + proc.stderr)
    if proc.returncode != 0 or not ok:
        sys.exit(f"ERROR: HYPE ({label}) failed (exit {proc.returncode}) in {run_dir}:\n"
                 f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
    for log in run_dir.glob("hyss_*.log"):
        text = log.read_text(errors="replace")
        if "WARNING" in text:
            warn = [l for l in text.splitlines() if "WARNING" in l]
            sys.exit(f"ERROR: HYPE ({label}) warned; cases must run clean:\n" + "\n".join(warn))
    return run_dir / "results"


def read_table(path):
    """HYPE time file -> (column ids, rows of (date, [float...]))."""
    lines = [l.rstrip("\r\n") for l in open(path) if l.strip()]
    lines = [l for l in lines if not l.startswith("!!")]
    header = lines[0].split("\t")
    body = [l for l in lines[1:] if not l.startswith("UNITS")]
    rows = []
    for l in body:
        parts = l.split("\t")
        rows.append((parts[0], [float(x) for x in parts[1:]]))
    return header[1:], rows


def collect_outputs(case, results):
    """Return {('class'|'basin', var): (header, rows)} in a stable layout."""
    out = {}
    cids = class_ids(case)
    for var in case.get("class_outputs", []):
        merged_head, merged = [], None
        for g, cid in enumerate(cids, start=1):
            f = results / f"time{var.upper()}_DD_CG{g}.txt"
            if not f.exists():
                sys.exit(f"ERROR: HYPE wrote no class output for '{var}' ({f.name}); "
                         "the variable may not exist at class level")
            subs, rows = read_table(f)
            merged_head += [f"{s}_{cid}" for s in subs]
            if merged is None:
                merged = [[d, list(v)] for d, v in rows]
            else:
                for m, (d, v) in zip(merged, rows):
                    assert m[0] == d
                    m[1].extend(v)
        out[("class", var)] = (merged_head, merged)
    for var in case.get("basin_outputs", []):
        head, merged = [], None
        for s in case["subbasins"]:
            f = results / f"{s['subid']:07d}.txt"
            cols, rows = read_table(f)
            j = cols.index(var)
            head.append(str(s["subid"]))
            if merged is None:
                merged = [[d, [v[j]]] for d, v in rows]
            else:
                for m, (d, v) in zip(merged, rows):
                    m[1].append(v[j])
        out[("basin", var)] = (head, merged)
    return out


def compact(v):
    """HYPE printed 9 significant figures. Keep exactly those digits, in the
    shortest text: '0' for zero, plain decimals where that is shorter. Lossless
    with respect to what HYPE wrote, and roughly halves the file size."""
    if v == 0.0:
        return "0"
    return f"{v:.{SIGNIFICANT_FIGURES}g}"


def write_expected(outputs, d):
    d.mkdir(parents=True, exist_ok=True)
    for (kind, var), (head, rows) in sorted(outputs.items()):
        with open(d / f"{kind}_{var}.csv", "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["date"] + head)
            for dt, vals in rows:
                w.writerow([dt] + [compact(v) for v in vals])


def precision_floor(single, double):
    """Per variable: max |single - double| and max |double| over all columns and steps."""
    res = {}
    for key, (head, srows) in single.items():
        _, drows = double[key]
        dmax, amax = 0.0, 0.0
        for (_, sv), (_, dv) in zip(srows, drows):
            for a, b in zip(sv, dv):
                dmax = max(dmax, abs(a - b))
                amax = max(amax, abs(b))
        res[f"{key[0]}_{key[1]}"] = {"max_abs_diff": float(f"{dmax:.3e}"),
                                     "max_abs_value": float(f"{amax:.6e}")}
    return res


# ---------------------------------------------------------------------------
# Manifest and README
# ---------------------------------------------------------------------------
def script_revision():
    try:
        rev = lsdtt_test_env.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
                                 capture_output=True, text=True).stdout.strip()
        dirty = lsdtt_test_env.run(["git", "status", "--porcelain", "--", str(Path(__file__))],
                                   cwd=REPO_ROOT, capture_output=True, text=True).stdout.strip()
        return rev + ("+uncommitted" if dirty else "")
    except OSError:
        return "unknown"


README_TEMPLATE = """# HYPE oracle files

Reference output from the HYPE hydrological model, used by the lsdtt3 HYPE-column
tests. CTest reads these files only; it never downloads or builds HYPE.

## Provenance

| | |
|---|---|
| HYPE version | {version} |
| Source archive | `{archive_name}` |
| Download URL | {url} |
| Archive SHA-256 | `{sha256}` |
| SourceForge git master at generation | `{git_commit}` |
| Compiler | {gfortran} |
| Oracle build | shipped makefile, `make hype` (default REAL = single precision) |
| Precision-floor build | `make hype switch="{double_switch}"` |
| Output precision | {signif} significant figures (`signfigures`); 9 is the most HYPE's output format accepts |

HYPE is developed by SMHI and released under the GNU Lesser General Public
Licence v3. These files are HYPE *outputs* for input sets written for lsdtt3;
no HYPE source is included.

## Regenerating

    python3 tests/scripts/generate_hype_oracle.py            # regenerate
    python3 tests/scripts/generate_hype_oracle.py --check    # verify committed files

The script uses the latest HYPE release that ships source. If that differs from
the version above, it stops until rerun with `--accept-new-version`.
Regenerating is a reviewable event: commit the changed files and state why.

## Layout

    cases/<case>.json             case definitions (the single source for both views)
    <case>/lsdtt_input/           the set-up as lsdtt3 tables (the HYPE input files are
                                  regenerated from cases/<case>.json; --keep-work keeps them)
    <case>/expected/class_<v>.csv per-class output; columns <subid>_<classid>
    <case>/expected/basin_<v>.csv per-sub-basin output; columns <subid>
    manifest.json                 provenance and precision floors

Values are HYPE's printed 9 significant figures, written in shortest form
('0', '12.3456789'). Forcing uses a compressed seasonal cycle
(temp_period_steps) so each case is short.

All runs use daily steps, one precipitation and one temperature series for the
whole domain, and HYPE's default model options unless a case lists others.
Initial state is HYPE's default: soil at field capacity, no snow, rivers at
their dead volume.

## Precision floor (tolerance basis)

The oracle is HYPE's released, single-precision build. The same cases were run
with a double-precision build; the largest difference per variable is the
floor below which a double-precision reimplementation cannot be judged against
the single-precision oracle.

| Case | Variable | Max abs diff single vs double | Max abs value |
|---|---|---|---|
{floor_rows}

## Cases

{case_rows}
"""


def write_readme(out, manifest):
    h = manifest["hype"]
    floor_rows = []
    case_rows = []
    for name, c in sorted(manifest["cases"].items()):
        for var, f in sorted(c["precision_floor"].items()):
            floor_rows.append(f"| {name} | {var} | {f['max_abs_diff']:.3e} | {f['max_abs_value']:.6g} |")
        case_rows.append(f"- **{name}** ({c['n_steps']} daily steps): {c['description']}")
    (out / "README.md").write_text(README_TEMPLATE.format(
        version=h["version"], archive_name=h["archive"], url=h["url"], sha256=h["sha256"],
        git_commit=h.get("git_master_commit") or "not determined",
        gfortran=manifest["build"]["gfortran"], double_switch=DOUBLE_SWITCH,
        signif=SIGNIFICANT_FIGURES, floor_rows="\n".join(floor_rows),
        case_rows="\n".join(case_rows)))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def generate(out, archive_arg, accept_new, keep_work, check_mode):
    cases_dir = (DEFAULT_OUT if check_mode else out) / "cases"
    case_files = sorted(cases_dir.glob("*.json"))
    if not case_files:
        sys.exit(f"ERROR: no case definitions in {cases_dir}")
    old_manifest_path = DEFAULT_OUT / "manifest.json"
    pinned = json.loads(old_manifest_path.read_text())["hype"] if old_manifest_path.exists() else None

    work = Path(tempfile.mkdtemp(prefix="hype_oracle_"))
    try:
        if archive_arg:
            archive = Path(archive_arg).resolve()
            m = re.search(r"hype_(\d+)_(\d+)_(\d+)_src\.tgz$", archive.name)
            if not m:
                sys.exit("ERROR: --hype-archive must be named hype_X_Y_Z_src.tgz")
            version = ".".join(m.groups())
            url = (f"{SF_FILES}release_hype_{'_'.join(m.groups())}/{archive.name}/download")
        else:
            version, url = latest_release_with_source()
            archive = work / f"hype_{version.replace('.', '_')}_src.tgz"
            print(f"Latest HYPE release with source: {version}\nDownloading {url}")
            archive.write_bytes(http_get(url, binary=True))
        sha = sha256_of(archive)
        print(f"SHA-256 {sha}")
        if pinned and (pinned["version"] != version or pinned["sha256"] != sha) and not accept_new:
            sys.exit(f"ERROR: HYPE {version} (sha256 {sha}) differs from the pinned "
                     f"{pinned['version']} (sha256 {pinned['sha256']}). Rerun with "
                     "--accept-new-version to regenerate against it.")
        if not pinned and not accept_new and not check_mode:
            sys.exit("ERROR: no pinned version yet (no manifest.json). Rerun with "
                     "--accept-new-version to pin this one.")

        print("Building HYPE (single precision, shipped makefile) ...")
        exe_s = build_hype(archive, work, "single")
        print("Building HYPE (double precision, measurement only) ...")
        exe_d = build_hype(archive, work, "double", DOUBLE_SWITCH)

        manifest = {
            "hype": {"version": version, "archive": archive.name, "url": url, "sha256": sha,
                     "git_master_commit": git_commit_of_release(version)},
            "build": {"gfortran": gfortran_version(), "single": "make hype",
                      "double": f'make hype switch="{DOUBLE_SWITCH}"',
                      "significant_figures": SIGNIFICANT_FIGURES},
            "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "host": f"{platform.system()} {platform.release()} {platform.machine()}",
            "script_git_revision": script_revision(),
            "cases": {},
        }
        if not check_mode:
            out.mkdir(parents=True, exist_ok=True)
        else:
            (out / "cases").mkdir(parents=True)
            for cf in case_files:
                shutil.copy2(cf, out / "cases" / cf.name)

        for cf in case_files:
            case = json.loads(cf.read_text())
            name = case["name"]
            if name != cf.stem:
                sys.exit(f"ERROR: {cf.name}: 'name' must equal the file stem")
            print(f"Case {name} ...")
            dates = dates_of(case)
            prec, temp = forcing_of(case, dates)
            cdir = out / name
            if cdir.exists():
                shutil.rmtree(cdir)
            # The HYPE view is regenerated from the case JSON every run, so it
            # lives only in the work directory (--keep-work to inspect it).
            hype_in = work / f"{name}_hype_input"
            write_hype_inputs(case, dates, prec, temp, hype_in)
            write_lsdtt_inputs(case, dates, prec, temp, cdir / "lsdtt_input")
            res_s = collect_outputs(case, run_hype(exe_s, hype_in, work, f"{name}_single"))
            res_d = collect_outputs(case, run_hype(exe_d, hype_in, work, f"{name}_double"))
            write_expected(res_s, cdir / "expected")
            manifest["cases"][name] = {
                "description": case["description"],
                "n_steps": len(dates),
                "modeloptions": case.get("modeloptions", {}),
                "precision_floor": precision_floor(res_s, res_d),
            }

        (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        write_readme(out, manifest)
    finally:
        if keep_work:
            print(f"Work directory kept: {work}")
        else:
            shutil.rmtree(work, ignore_errors=True)


def compare_trees(new, old):
    problems = []
    for p in sorted(new.rglob("*")):
        if p.is_dir():
            continue
        rel = p.relative_to(new)
        q = old / rel
        if not q.exists():
            problems.append(f"missing in committed set: {rel}")
        elif rel.name == "manifest.json":
            a = json.loads(p.read_text())
            b = json.loads(q.read_text())
            for k in VOLATILE_KEYS:
                a.pop(k, None)
                b.pop(k, None)
            if a != b:
                problems.append("manifest.json differs (beyond provenance fields)")
        elif not filecmp.cmp(p, q, shallow=False):
            problems.append(f"differs: {rel}")
    for q in sorted(old.rglob("*")):
        if q.is_file() and not (new / q.relative_to(old)).exists():
            problems.append(f"not regenerated: {q.relative_to(old)}")
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT))
    ap.add_argument("--hype-archive", help="local hype_X_Y_Z_src.tgz instead of downloading")
    ap.add_argument("--accept-new-version", action="store_true")
    ap.add_argument("--check", action="store_true",
                    help="regenerate into a temporary directory and compare with the committed files")
    ap.add_argument("--keep-work", action="store_true", help="keep the temporary build/run directory")
    a = ap.parse_args()

    if shutil.which("gfortran") is None or shutil.which("make") is None:
        sys.exit("ERROR: gfortran and make are required")

    if a.check:
        tmp = Path(tempfile.mkdtemp(prefix="hype_oracle_check_"))
        try:
            generate(tmp / "hype", a.hype_archive, a.accept_new_version, a.keep_work, True)
            problems = compare_trees(tmp / "hype", Path(a.out_dir))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        if problems:
            print("CHECK FAILED:\n  " + "\n  ".join(problems))
            sys.exit(1)
        print("CHECK PASSED: regenerated files are identical to the committed set.")
    else:
        generate(Path(a.out_dir), a.hype_archive, a.accept_new_version, a.keep_work, False)
        print(f"Done. Oracle written to {a.out_dir}")


if __name__ == "__main__":
    main()
