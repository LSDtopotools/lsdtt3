#!/usr/bin/env python3
"""Single-point scaling scan dataset generator (Phase M2 wired).

Generates point-scan datasets across latitude/elevation/erosion-rate space for
Be10/Al26/C14, runs lsdtt3 forward->inverse at each point, and optionally joins
CRONUS v3 references.

This script is data-workflow infrastructure for:
  docs/plan_docs/cosmo_visualisation_scaling_plan.md

Notes:
  - Uses synthetic 3x3 ASCII DEMs centred at each scan point.
  - lsdtt3 concentration comes from lsdtt-cosmo-prediction sampled points CSV.
  - lsdtt3 erosion inversion comes from lsdtt-cosmo-tool.
  - CRONUS v3 is optional (network POST).
"""

from __future__ import annotations

import argparse
import csv
import datetime
import math
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

from pyproj import Transformer

EXIT_OK = 0
EXIT_FAIL = 1
EXIT_SKIP = 77

NUCLIDES = ("Be10", "Al26", "C14")
SCALING_VALUES = ("stone", "lsdn")
MUON_VALUES = ("braucherborchers", "tunedexp", "emulator")

MUON_TO_LSDTT = {
    "braucherborchers": "BraucherBorchers",
    "tunedexp": "newCRONUS",
    # The CRONUS emulator comes in a St and an LSDn flavour that must be paired
    # with the matching scaling_scheme: its spallation reference production is
    # scheme-specific (refP_St vs refP_LSDn).  Resolved per (scaling, muon)
    # combination by lsdtt_muon_name(); "emulator" is the scan-level alias.
    "emulator": None,
}


def lsdtt_muon_name(muon_scheme: str, scaling: str) -> str:
    """Scan muon alias + scaling scheme -> lsdtt3 muon_scheme value."""
    if muon_scheme == "emulator":
        return "cronusEmulatorLSDn" if scaling == "lsdn" else "cronusEmulatorSt"
    return MUON_TO_LSDTT[muon_scheme]

FILE_INPUTS = "single_point_scan_inputs.csv"
FILE_LSDTT = "single_point_scan_predictions_lsdtt3.csv"
FILE_CRONUS = "single_point_scan_reference_cronus_v3.csv"
FILE_COMPARISON = "single_point_scan_comparison.csv"

DEFAULT_CRONUS_ENDPOINTS = [
    "https://hess.ess.washington.edu/cgi-bin/matweb",
    "http://stoneage.hzdr.de/cgi/matweb",
]

INPUT_FIELDS = [
    "run_id",
    "lat_deg",
    "lon_deg",
    "pressure_hpa",
    "elevation_m",
    "erate_gcm2yr",
    "sample_age_yr",
    "atmosphere_source",
    "timestamp_utc",
    "git_sha",
]

MODEL_FIELDS = [
    "run_id",
    "source",
    "nuclide",
    "scaling_scheme",
    "muon_scheme",
    "lat_deg",
    "lon_deg",
    "pressure_hpa",
    "elevation_m",
    "erate_gcm2yr",
    "conc_atoms_g",
    "conc_sigma_atoms_g",
    "erate_est_gcm2yr",
    "sample_age_yr",
    "atmosphere_source",
    "timestamp_utc",
    "git_sha",
]

COMPARISON_FIELDS = [
    "run_id",
    "nuclide",
    "scaling_scheme",
    "muon_scheme",
    "lat_deg",
    "lon_deg",
    "pressure_hpa",
    "elevation_m",
    "erate_gcm2yr",
    "conc_atoms_g",
    "conc_sigma_atoms_g",
    "erate_lsdtt3_gcm2yr",
    "erate_cronus_v3_gcm2yr",
    "cronus_v3_status",
    "cronus_v3_zero_erate_flag",
    "pct_diff_vs_cronus_v3",
    "pct_diff_lsdtt3_minus_reference",
    "sample_age_yr",
    "atmosphere_source",
    "timestamp_utc",
    "git_sha",
]

IGNORE_CHECK_FIELDS = {"timestamp_utc", "git_sha"}


@dataclass(frozen=True)
class GridPoint:
    run_id: str
    lat_deg: float
    lon_deg: float
    elevation_m: float
    pressure_hpa: float
    erate_gcm2yr: float


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_expected_dir() -> Path:
    return repo_root() / "tests" / "expected" / "cosmo"


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


def std_pressure_hpa(elevation_m: float) -> float:
    t = 1.0 - 2.25577e-5 * elevation_m
    return 1013.25 * max(t, 1e-4) ** 5.25588


def parse_choice(value: str, valid_values: tuple[str, ...]) -> list[str]:
    if value == "both":
        return list(valid_values)
    if value not in valid_values:
        raise ValueError(f"invalid choice '{value}', expected one of {valid_values} or 'both'")
    return [value]


def build_grid(preset: str) -> list[GridPoint]:
    if preset == "mini":
        lats = [-60.0, -30.0, 0.0, 30.0, 60.0]
        elevations = [0.0, 1000.0, 2000.0, 3000.0]
        erates = [0.001, 0.003, 0.01, 0.03, 0.1]
    elif preset == "latsweep":
        # Plan Phase 4 validation grid.  The erosion-rate ceiling is 1.0, not
        # 0.1: the mini/full presets stopped at 0.1, which is a large part of
        # why the balcoCRONUSv3PD defect went unnoticed on real samples.
        lats = [-60.0, -30.0, -1.0, 30.0, 60.0]
        elevations = [0.0, 500.0, 1000.0, 1500.0, 2500.0, 3500.0]
        erates = [0.001, 0.01, 0.05, 0.1, 0.3, 0.55, 1.0]
    elif preset == "full":
        lats = [-75.0, -60.0, -45.0, -30.0, -15.0, 0.0, 15.0, 30.0, 45.0, 60.0, 75.0]
        elevations = [0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0]
        erates = [0.0007, 0.001, 0.002, 0.003, 0.007, 0.01, 0.02, 0.03, 0.07, 0.1]
    else:
        raise ValueError(f"unknown preset {preset}")

    points: list[GridPoint] = []
    idx = 0
    lon = 0.0
    for lat in lats:
        for elev in elevations:
            p = std_pressure_hpa(elev)
            for e in erates:
                idx += 1
                points.append(
                    GridPoint(
                        run_id=f"P{idx:05d}",
                        lat_deg=lat,
                        lon_deg=lon,
                        elevation_m=elev,
                        pressure_hpa=p,
                        erate_gcm2yr=e,
                    )
                )
    return points


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def is_float_like(s: str) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False


def compare_csv_rows(new_rows: list[dict], old_rows: list[dict], key: str) -> tuple[bool, str]:
    if len(new_rows) != len(old_rows):
        return False, f"row count mismatch new={len(new_rows)} frozen={len(old_rows)}"
    if not new_rows:
        return False, "no rows"
    old_by_key = {r[key]: r for r in old_rows}
    for r in new_rows:
        k = r[key]
        if k not in old_by_key:
            return False, f"missing key '{k}' in frozen"
        o = old_by_key[k]
        for col, v_new in r.items():
            if col in IGNORE_CHECK_FIELDS:
                continue
            v_old = o.get(col, "")
            if is_float_like(v_new) and is_float_like(v_old):
                a = float(v_new)
                b = float(v_old)
                if math.isnan(a) and math.isnan(b):
                    continue
                denom = abs(b) if b != 0.0 else 1.0
                if abs(a - b) / denom > 1e-10:
                    return False, f"numeric mismatch key={k} col={col} new={a} frozen={b}"
            elif v_new != v_old:
                return False, f"value mismatch key={k} col={col} new='{v_new}' frozen='{v_old}'"
    return True, "match"


def run_cli(exe: Path, workdir: Path, param_file: Path, label: str) -> None:
    proc = subprocess.run(
        [str(exe), str(workdir), param_file.name],
        cwd=str(workdir),
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"{label} failed (rc={proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-3000:]}\n"
            f"--- stderr ---\n{proc.stderr[-3000:]}"
        )


def write_point_and_dem(workdir: Path, point: GridPoint) -> tuple[Path, Path]:
    dem_asc = workdir / f"dem_{point.run_id}.asc"
    # 3x3 projected DEM (EPSG:3857) centred on transformed (lon,lat),
    # constant elevation. Cellsize 30 m gives a compact valid projected raster.
    tx = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    x_c, y_c = tx.transform(point.lon_deg, point.lat_deg)
    cell = 30.0
    xll = x_c - 1.5 * cell
    yll = y_c - 1.5 * cell
    z = f"{point.elevation_m:.3f}"
    dem_asc.write_text(
        "\n".join(
            [
                "ncols         3",
                "nrows         3",
                f"xllcorner     {xll:.8f}",
                f"yllcorner     {yll:.8f}",
                f"cellsize      {cell:.6f}",
                "NODATA_value  -9999",
                f"{z} {z} {z}",
                f"{z} {z} {z}",
                f"{z} {z} {z}",
            ]
        )
        + "\n"
    )

    dem_tif = workdir / f"dem_{point.run_id}.tif"
    proc = subprocess.run(
        [
            "gdal_translate",
            "-of",
            "GTiff",
            "-a_srs",
            "EPSG:3857",
            str(dem_asc),
            str(dem_tif),
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not dem_tif.is_file():
        raise RuntimeError(
            f"failed to create projected DEM TIFF for {point.run_id}\n"
            f"stdout: {proc.stdout[-1000:]}\n"
            f"stderr: {proc.stderr[-1000:]}"
        )

    pts = workdir / f"points_{point.run_id}.csv"
    with open(pts, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_name", "latitude", "longitude"])
        w.writerow([point.run_id, f"{point.lat_deg:.8f}", f"{point.lon_deg:.8f}"])
    return dem_tif, pts


def lsdtt_forward_concentration(
    pred_exe: Path,
    workdir: Path,
    point: GridPoint,
    nuclide: str,
    scaling: str,
    muon_scheme: str,
    sample_age_yr: float,
    atmosphere_source: str,
) -> float:
    dem, pts = write_point_and_dem(workdir, point)
    stem = f"fwd_{point.run_id}_{nuclide}_{scaling}_{muon_scheme}"
    param = workdir / f"{stem}.param"
    lines = [
        f"read_path: {workdir}",
        f"read_fname: {dem.name}",
        f"write_path: {workdir}",
        f"write_fname: {stem}",
        "calculate_crn_concentration_raster: true",
        "calculate_accumulated_crn_concentration: true",
        f"effective_erosion_rate: {point.erate_gcm2yr}",
        f"nuclide: {nuclide}",
        f"scaling_scheme: {scaling}",
        f"muon_scheme: {lsdtt_muon_name(muon_scheme, scaling)}",
        f"sample_age_yr: {sample_age_yr}",
        f"atmosphere: {atmosphere_source}",
        f"sample_points_filename: {pts.name}",
        "threshold_contributing_pixels: 1",
        "snap_search_radius_pixels: 2",
        "fill_before_routing: true",
    ]
    param.write_text("\n".join(lines) + "\n")
    run_cli(pred_exe, workdir, param, f"forward[{stem}]")

    sampled = workdir / f"{stem}_sampled_crn_conc.csv"
    rows = list(csv.DictReader(sampled.open(newline="")))
    if len(rows) != 1:
        raise RuntimeError(f"forward expected 1 sampled row, got {len(rows)}")
    conc_col = next((c for c in rows[0] if c.endswith("_CONC")), None)
    if conc_col is None:
        raise RuntimeError(f"forward output has no *_CONC column: {sampled}")
    return float(rows[0][conc_col])


def lsdtt_inverse_erate(
    tool_exe: Path,
    workdir: Path,
    point: GridPoint,
    nuclide: str,
    scaling: str,
    muon_scheme: str,
    concentration: float,
    conc_sigma: float,
    sample_age_yr: float,
    atmosphere_source: str,
) -> float:
    dem, _ = write_point_and_dem(workdir, point)

    points = workdir / f"inv_points_{point.run_id}_{nuclide}_{scaling}_{muon_scheme}.csv"
    with open(points, "w", newline="") as f:
        if nuclide == "Be10":
            fields = ["sample_name", "latitude", "longitude", "Be10_CONC", "Be10_uncertainty", "standardisation"]
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerow(
                {
                    "sample_name": point.run_id,
                    "latitude": f"{point.lat_deg:.8f}",
                    "longitude": f"{point.lon_deg:.8f}",
                    "Be10_CONC": f"{concentration:.10g}",
                    "Be10_uncertainty": f"{conc_sigma:.10g}",
                    "standardisation": "07KNSTD",
                }
            )
        elif nuclide == "Al26":
            fields = ["sample_name", "latitude", "longitude", "Al26_CONC", "Al26_uncertainty", "standardisation"]
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerow(
                {
                    "sample_name": point.run_id,
                    "latitude": f"{point.lat_deg:.8f}",
                    "longitude": f"{point.lon_deg:.8f}",
                    "Al26_CONC": f"{concentration:.10g}",
                    "Al26_uncertainty": f"{conc_sigma:.10g}",
                    "standardisation": "KNSTD",
                }
            )
        else:
            fields = ["sample_name", "latitude", "longitude", "C14_CONC", "C14_uncertainty"]
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerow(
                {
                    "sample_name": point.run_id,
                    "latitude": f"{point.lat_deg:.8f}",
                    "longitude": f"{point.lon_deg:.8f}",
                    "C14_CONC": f"{concentration:.10g}",
                    "C14_uncertainty": f"{conc_sigma:.10g}",
                }
            )

    stem = f"inv_{point.run_id}_{nuclide}_{scaling}_{muon_scheme}"
    param = workdir / f"{stem}.param"
    lines = [
        f"read_path: {workdir}",
        f"read_fname: {dem.name}",
        f"write_path: {workdir}",
        f"write_fname: {stem}",
        f"points_filename: {points.name}",
        f"nuclide: {nuclide}",
        f"scaling_scheme: {scaling}",
        f"muon_scheme: {lsdtt_muon_name(muon_scheme, scaling)}",
        f"sample_age_yr: {sample_age_yr}",
        f"atmosphere: {atmosphere_source}",
        "threshold_contributing_pixels: 1",
        "snap_search_radius_pixels: 2",
    ]
    param.write_text("\n".join(lines) + "\n")
    run_cli(tool_exe, workdir, param, f"inverse[{stem}]")

    out = workdir / f"{stem}_CN_ERATE.csv"
    rows = list(csv.DictReader(out.open(newline="")))
    if len(rows) != 1:
        raise RuntimeError(f"inverse expected 1 row, got {len(rows)}")
    return float(rows[0]["apparent_eff_e_g_cm2_yr"])


def build_cronus_text_block(obs_rows: list[dict]) -> tuple[str, dict[str, str]]:
    lines: list[str] = []
    idmap: dict[str, str] = {}
    for i, r in enumerate(obs_rows):
        # CRONUS v3 requires sample names <= 32 chars.
        name = f"SP{i:06d}"
        idmap[name] = r["run_id"]
        lat = float(r["lat_deg"])
        lon = float(r["lon_deg"])
        p = float(r["pressure_hpa"])
        conc = float(r["conc_atoms_g"])
        sigma = float(r["conc_sigma_atoms_g"])
        nuclide = r["nuclide"]
        lines.append(f"{name} {lat:.5f} {lon:.5f} {p:.4f} pre 0 2.65 1 0.0001 2010 ;")
        if nuclide == "Be10":
            lines.append(f"{name} Be-10 quartz {conc:.6f} {sigma:.6f} 07KNSTD ;")
        elif nuclide == "Al26":
            lines.append(f"{name} Al-26 quartz {conc:.6f} {sigma:.6f} KNSTD ;")
        else:
            lines.append(f"{name} C-14 quartz {conc:.6f} {sigma:.6f} ;")
    return "\n".join(lines), idmap


def post_cronus(endpoint: str, text_block: str) -> str:
    data = urllib.parse.urlencode(
        {
            "mlmfile": "erosion_input_v3",
            "reportType": "XML",
            "resultType": "long",
            "summary": "no",
            "plotFlag": "no",
            "text_block": text_block,
        }
    ).encode()
    req = urllib.request.Request(endpoint, data=data)
    with urllib.request.urlopen(req, timeout=300) as resp:
        return resp.read().decode(errors="replace")


def parse_cronus_xml(raw: str) -> dict[tuple[str, str], float]:
    start = raw.find("<")
    if start < 0:
        raise RuntimeError("CRONUS response contains no XML")
    root = ET.fromstring(raw[start:])

    out: dict[tuple[str, str], float] = {}
    for res in root.iter():
        if not res.tag.endswith("erosionRateResult"):
            continue
        name = (res.findtext("sample_name") or "").strip()
        for nres in res.findall("nuclide_result"):
            n = (nres.findtext("nuclide") or "").strip()
            if "10" in n:
                nuclide = "Be10"
            elif "26" in n:
                nuclide = "Al26"
            elif "14" in n:
                nuclide = "C14"
            else:
                continue
            for scheme_tag, scheme_name in (("St", "stone"), ("LSDn", "lsdn")):
                t = nres.findtext(f"E_gcm2_{scheme_tag}")
                if t is None:
                    continue
                out[(name, f"{nuclide}|{scheme_name}")] = float(t)
    return out


def fetch_cronus_v3(obs_rows: list[dict], endpoints: list[str], raw_xml_cache: Path | None) -> dict[str, float]:
    text_block, idmap = build_cronus_text_block(obs_rows)
    raw = None
    last_err = None
    for ep in endpoints:
        try:
            raw = post_cronus(ep, text_block)
            break
        except Exception as exc:  # noqa: BLE001
            last_err = exc
    if raw is None:
        raise RuntimeError(f"all CRONUS endpoints failed; last error: {last_err}")
    if raw_xml_cache is not None:
        raw_xml_cache.write_text(raw)
    parsed = parse_cronus_xml(raw)
    refs: dict[str, float] = {}
    for sid, rid in idmap.items():
        # scaling is already encoded in run_id; nuclide/scaling also in row data.
        # Map by scanning obs_rows for this run_id entry.
        row = next((x for x in obs_rows if x["run_id"] == rid), None)
        if row is None:
            continue
        key = (sid, f"{row['nuclide']}|{row['scaling_scheme']}")
        if key in parsed:
            refs[rid] = parsed[key]
    return refs


def pct_diff(a: float, b: float) -> float:
    return 100.0 * (a - b) / b if b != 0.0 else float("nan")


def generate_rows(
    points: list[GridPoint],
    scaling_values: list[str],
    muon_values: list[str],
    sample_age_yr: float,
    atmosphere_source: str,
    pred_exe: Path,
    tool_exe: Path,
    with_cronus_v3: bool,
    cronus_endpoints: list[str],
    cronus_raw_xml_cache: Path | None,
) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    sha = git_sha()

    inputs: list[dict] = []
    lsdtt_rows: list[dict] = []
    obs_rows_for_cronus: list[dict] = []

    workdir = Path(tempfile.mkdtemp(prefix="single_point_scan_m2_"))
    try:
        for p in points:
            inputs.append(
                {
                    "run_id": p.run_id,
                    "lat_deg": f"{p.lat_deg:.6f}",
                    "lon_deg": f"{p.lon_deg:.6f}",
                    "pressure_hpa": f"{p.pressure_hpa:.6f}",
                    "elevation_m": f"{p.elevation_m:.2f}",
                    "erate_gcm2yr": f"{p.erate_gcm2yr:.10g}",
                    "sample_age_yr": f"{sample_age_yr:.1f}",
                    "atmosphere_source": atmosphere_source,
                    "timestamp_utc": now,
                    "git_sha": sha,
                }
            )

            for n in NUCLIDES:
                for s in scaling_values:
                    for m in muon_values:
                        combo_id = f"{p.run_id}_{n}_{s}_{m}"
                        conc = lsdtt_forward_concentration(
                            pred_exe=pred_exe,
                            workdir=workdir,
                            point=p,
                            nuclide=n,
                            scaling=s,
                            muon_scheme=m,
                            sample_age_yr=sample_age_yr,
                            atmosphere_source=atmosphere_source,
                        )
                        sigma = 0.03 * conc
                        erate = lsdtt_inverse_erate(
                            tool_exe=tool_exe,
                            workdir=workdir,
                            point=p,
                            nuclide=n,
                            scaling=s,
                            muon_scheme=m,
                            concentration=conc,
                            conc_sigma=sigma,
                            sample_age_yr=sample_age_yr,
                            atmosphere_source=atmosphere_source,
                        )

                        row_base = {
                            "run_id": combo_id,
                            "nuclide": n,
                            "scaling_scheme": s,
                            "muon_scheme": m,
                            "lat_deg": f"{p.lat_deg:.6f}",
                            "lon_deg": f"{p.lon_deg:.6f}",
                            "pressure_hpa": f"{p.pressure_hpa:.6f}",
                            "elevation_m": f"{p.elevation_m:.2f}",
                            "erate_gcm2yr": f"{p.erate_gcm2yr:.10g}",
                            "conc_atoms_g": f"{conc:.10g}",
                            "conc_sigma_atoms_g": f"{sigma:.10g}",
                            "sample_age_yr": f"{sample_age_yr:.1f}",
                            "atmosphere_source": atmosphere_source,
                            "timestamp_utc": now,
                            "git_sha": sha,
                        }

                        lsdtt_rows.append({**row_base, "source": "lsdtt3", "erate_est_gcm2yr": f"{erate:.10g}"})
                        obs_rows_for_cronus.append({**row_base})
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    cronus_refs: dict[str, float] = {}
    if with_cronus_v3:
        cronus_refs = fetch_cronus_v3(obs_rows_for_cronus, cronus_endpoints, cronus_raw_xml_cache)

    cronus_rows: list[dict] = []
    comparison_rows: list[dict] = []
    by_id_lsdtt = {r["run_id"]: r for r in lsdtt_rows}

    for rid, lr in by_id_lsdtt.items():
        e_lsdtt = float(lr["erate_est_gcm2yr"])
        e_cv3 = cronus_refs.get(rid)
        zero_reference = e_cv3 is not None and e_cv3 == 0.0
        if e_cv3 is None:
            cronus_v3_status = "missing_reference"
        elif zero_reference:
            cronus_v3_status = "zero_reference_erate"
        else:
            cronus_v3_status = "ok"

        if e_cv3 is not None:
            cronus_rows.append({**lr, "source": "cronus_v3", "erate_est_gcm2yr": f"{e_cv3:.10g}"})

        d_cv3 = pct_diff(e_lsdtt, e_cv3) if e_cv3 is not None else float("nan")
        d_mix = d_cv3

        comparison_rows.append(
            {
                "run_id": rid,
                "nuclide": lr["nuclide"],
                "scaling_scheme": lr["scaling_scheme"],
                "muon_scheme": lr["muon_scheme"],
                "lat_deg": lr["lat_deg"],
                "lon_deg": lr["lon_deg"],
                "pressure_hpa": lr["pressure_hpa"],
                "elevation_m": lr["elevation_m"],
                "erate_gcm2yr": lr["erate_gcm2yr"],
                "conc_atoms_g": lr["conc_atoms_g"],
                "conc_sigma_atoms_g": lr["conc_sigma_atoms_g"],
                "erate_lsdtt3_gcm2yr": lr["erate_est_gcm2yr"],
                "erate_cronus_v3_gcm2yr": "" if e_cv3 is None else f"{e_cv3:.10g}",
                "cronus_v3_status": cronus_v3_status,
                "cronus_v3_zero_erate_flag": "" if e_cv3 is None else ("1" if zero_reference else "0"),
                "pct_diff_vs_cronus_v3": "" if math.isnan(d_cv3) else f"{d_cv3:.8g}",
                "pct_diff_lsdtt3_minus_reference": "" if math.isnan(d_mix) else f"{d_mix:.8g}",
                "sample_age_yr": lr["sample_age_yr"],
                "atmosphere_source": lr["atmosphere_source"],
                "timestamp_utc": lr["timestamp_utc"],
                "git_sha": lr["git_sha"],
            }
        )

    return inputs, lsdtt_rows, cronus_rows, comparison_rows


def run_self_test() -> int:
    mini = build_grid("mini")
    if len(mini) != 100:
        print(f"SELFTEST FAIL: expected 100 mini points, got {len(mini)}")
        return EXIT_FAIL
    if MUON_TO_LSDTT["tunedexp"] != "newCRONUS":
        print("SELFTEST FAIL: tunedexp map")
        return EXIT_FAIL
    p = mini[0]
    ptest = std_pressure_hpa(p.elevation_m)
    if abs(ptest - p.pressure_hpa) > 1e-12:
        print("SELFTEST FAIL: pressure consistency")
        return EXIT_FAIL
    print("SELFTEST: PASS")
    return EXIT_OK


def copy_file(src: Path, dst: Path) -> None:
    dst.write_text(src.read_text())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="mini", choices=["mini", "full", "latsweep"])
    ap.add_argument("--max-points", type=int, default=0,
                    help="optional cap on number of grid points (0 = no cap)")
    ap.add_argument("--scaling", default="both", help="stone|lsdn|both")
    ap.add_argument("--muon-scheme", default="both", help="braucherborchers|tunedexp|emulator|both")
    ap.add_argument("--sample-age-yr", type=float, default=20000.0)
    ap.add_argument("--atmosphere-source", default="standard",
                    choices=["standard", "ERA40", "NCEP"])

    ap.add_argument("--pred-exe", default="/mnt/c/LSDTopoTools/lsdtt3/build/lsdtt-cosmo-prediction")
    ap.add_argument("--tool-exe", default="/mnt/c/LSDTopoTools/lsdtt3/build/lsdtt-cosmo-tool")

    ap.add_argument("--with-cronus-v3", action="store_true",
                    help="query CRONUS v3 online calculator")
    ap.add_argument("--cronus-endpoint", action="append", default=None,
                    help="CRONUS endpoint (can repeat); defaults to known mirrors")
    ap.add_argument("--cronus-raw-xml-cache", help="optional path to cache raw CRONUS XML")

    ap.add_argument("--outdir", help="directory for generated run CSVs")
    ap.add_argument("--expected-dir", help="tracked expected directory (freeze/check target)")
    ap.add_argument("--freeze", action="store_true", help="write tracked expected CSVs")
    ap.add_argument("--check", action="store_true", help="check generated vs tracked expected CSVs")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return run_self_test()

    if args.freeze and args.check:
        print("FAIL: choose only one of --freeze or --check")
        return EXIT_FAIL

    pred_exe = Path(args.pred_exe)
    tool_exe = Path(args.tool_exe)
    if not pred_exe.is_file() or not tool_exe.is_file():
        print(f"SKIP: missing executable(s): pred={pred_exe} tool={tool_exe}")
        return EXIT_SKIP

    try:
        scaling_values = parse_choice(args.scaling, SCALING_VALUES)
        muon_values = parse_choice(args.muon_scheme, MUON_VALUES)
    except ValueError as exc:
        print(f"FAIL: {exc}")
        return EXIT_FAIL

    points = build_grid(args.preset)
    if args.max_points > 0:
        points = points[: args.max_points]
    if not points:
        print("FAIL: no points selected")
        return EXIT_FAIL

    expected_dir = Path(args.expected_dir) if args.expected_dir else default_expected_dir()
    outdir = Path(args.outdir) if args.outdir else expected_dir / "generated"
    outdir.mkdir(parents=True, exist_ok=True)

    endpoints = args.cronus_endpoint or DEFAULT_CRONUS_ENDPOINTS
    raw_cache = Path(args.cronus_raw_xml_cache) if args.cronus_raw_xml_cache else None

    try:
        inputs, lsdtt_rows, cronus_rows, comparison_rows = generate_rows(
            points=points,
            scaling_values=scaling_values,
            muon_values=muon_values,
            sample_age_yr=args.sample_age_yr,
            atmosphere_source=args.atmosphere_source,
            pred_exe=pred_exe,
            tool_exe=tool_exe,
            with_cronus_v3=args.with_cronus_v3,
            cronus_endpoints=endpoints,
            cronus_raw_xml_cache=raw_cache,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"FAIL: generation error: {exc}")
        return EXIT_FAIL

    paths = {
        FILE_INPUTS: outdir / FILE_INPUTS,
        FILE_LSDTT: outdir / FILE_LSDTT,
        FILE_CRONUS: outdir / FILE_CRONUS,
        FILE_COMPARISON: outdir / FILE_COMPARISON,
    }
    write_csv(paths[FILE_INPUTS], INPUT_FIELDS, inputs)
    write_csv(paths[FILE_LSDTT], MODEL_FIELDS, lsdtt_rows)
    write_csv(paths[FILE_CRONUS], MODEL_FIELDS, cronus_rows)
    write_csv(paths[FILE_COMPARISON], COMPARISON_FIELDS, comparison_rows)

    print(f"WROTE: {paths[FILE_INPUTS]}")
    print(f"WROTE: {paths[FILE_LSDTT]}")
    print(f"WROTE: {paths[FILE_CRONUS]}")
    print(f"WROTE: {paths[FILE_COMPARISON]}")
    print(
        "ROWS: "
        f"points={len(inputs)} "
        f"lsdtt={len(lsdtt_rows)} "
        f"cronus_v3={len(cronus_rows)}"
    )

    if args.freeze:
        expected_dir.mkdir(parents=True, exist_ok=True)
        for name, src in paths.items():
            dst = expected_dir / name
            copy_file(src, dst)
            print(f"FROZE: {dst}")
        return EXIT_OK

    if args.check:
        checks = [
            (FILE_INPUTS, "run_id"),
            (FILE_LSDTT, "run_id"),
            (FILE_CRONUS, "run_id"),
            (FILE_COMPARISON, "run_id"),
        ]
        for fname, key in checks:
            fr = expected_dir / fname
            if not fr.is_file():
                print(f"CHECK FAIL: missing frozen file {fr}")
                return EXIT_FAIL
            ok, msg = compare_csv_rows(read_csv(paths[fname]), read_csv(fr), key)
            print(f"CHECK {fname}: {'PASS' if ok else 'FAIL'} ({msg})")
            if not ok:
                return EXIT_FAIL
        return EXIT_OK

    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
