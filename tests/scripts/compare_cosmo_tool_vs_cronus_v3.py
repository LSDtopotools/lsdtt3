#!/usr/bin/env python3
"""Compare lsdtt-cosmo-tool sample erosion rates against CRONUS v3.

Reads ``<write_fname>_CN_ERATE.csv`` from lsdtt-cosmo-tool, builds a CRONUS v3
batch request using per-sample nuclide, concentration, latitude and effective
pressure, submits once, and writes a comparison CSV.
"""

from __future__ import annotations

import argparse
import csv
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

DEFAULT_ENDPOINTS = [
    "https://hess.ess.washington.edu/cgi-bin/matweb",
    "http://stoneage.hzdr.de/cgi/matweb",
]


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _infer_nuclide(row: dict) -> str:
    if row.get("nuclide"):
        return row["nuclide"]
    for key in row:
        if key.startswith("concentration_"):
            return key.split("concentration_", 1)[1]
    raise ValueError("Could not infer nuclide from cosmo CSV row")


def _conc_field(nuclide: str) -> str:
    return f"concentration_{nuclide}"


def _default_std(nuclide: str) -> str:
    if nuclide == "Be10":
        return "07KNSTD"
    if nuclide == "Al26":
        return "KNSTD"
    return ""


def _build_text_block(rows: list[dict], points_by_name: dict[str, dict], uncertainty_col: str | None, standardisation_col: str | None) -> tuple[str, dict[str, tuple[str, str]]]:
    lines: list[str] = []
    idmap: dict[str, tuple[str, str]] = {}
    for i, r in enumerate(rows, start=1):
        sample = r.get("sample_name") or f"sample_{i}"
        nuclide = _infer_nuclide(r)
        conc = float(r[_conc_field(nuclide)])
        lat = float(r.get("centroid_latitude_deg") or r["latitude"])
        lon = float(r["longitude"])
        pressure = float(r["effective_pressure_hpa"])
        sid = f"S{i:05d}"
        idmap[sid] = (sample, nuclide)

        sigma = 0.03 * conc
        std = _default_std(nuclide)
        p = points_by_name.get(sample)
        if p is not None:
            if uncertainty_col and p.get(uncertainty_col, ""):
                sigma = float(p[uncertainty_col])
            if standardisation_col and p.get(standardisation_col, ""):
                std = p[standardisation_col]

        lines.append(
            f"{sid} {lat:.6f} {lon:.6f} {pressure:.4f} pre 0 2.65 1 0.0001 2010 ;"
        )
        if nuclide == "C14":
            lines.append(f"{sid} C-14 quartz {conc:.6f} {sigma:.6f} ;")
        elif nuclide == "Al26":
            lines.append(f"{sid} Al-26 quartz {conc:.6f} {sigma:.6f} {std or 'KNSTD'} ;")
        else:
            lines.append(f"{sid} Be-10 quartz {conc:.6f} {sigma:.6f} {std or '07KNSTD'} ;")

    return "\n".join(lines), idmap


def _post(endpoint: str, text_block: str) -> str:
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


def _parse_xml(raw: str, scheme: str) -> dict[tuple[str, str], float]:
    start = raw.find("<")
    if start < 0:
        raise RuntimeError("CRONUS response did not contain XML")
    root = ET.fromstring(raw[start:])
    out: dict[tuple[str, str], float] = {}
    for res in root.iter():
        if not res.tag.endswith("erosionRateResult"):
            continue
        sid = (res.findtext("sample_name") or "").strip()
        for nres in res.findall("nuclide_result"):
            nu = (nres.findtext("nuclide") or "").strip()
            if "10" in nu:
                nu = "Be10"
            elif "26" in nu:
                nu = "Al26"
            elif "14" in nu:
                nu = "C14"
            key = f"E_gcm2_{scheme}"
            e = nres.findtext(key)
            if e is not None:
                out[(sid, nu)] = float(e)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cosmo-csv", required=True, type=Path)
    ap.add_argument("--points-csv", type=Path, default=None)
    ap.add_argument("--uncertainty-col", default=None)
    ap.add_argument("--standardisation-col", default="standardisation")
    ap.add_argument("--scheme", default="St", choices=["St", "Lm", "LSDn"])
    ap.add_argument("--endpoint", action="append", default=None)
    ap.add_argument("--raw-xml-out", type=Path, default=None)
    ap.add_argument("--out-csv", required=True, type=Path)
    args = ap.parse_args()

    rows = _read_csv(args.cosmo_csv)
    if not rows:
        raise SystemExit("No rows in cosmo csv")

    if "effective_pressure_hpa" not in rows[0]:
        raise SystemExit("cosmo csv missing effective_pressure_hpa; rerun lsdtt-cosmo-tool with updated build")

    points_by_name: dict[str, dict] = {}
    if args.points_csv is not None:
        points_rows = _read_csv(args.points_csv)
        for r in points_rows:
            name = r.get("sample_name")
            if name:
                points_by_name[name] = r

    text_block, idmap = _build_text_block(
        rows,
        points_by_name,
        args.uncertainty_col,
        args.standardisation_col,
    )

    endpoints = args.endpoint or DEFAULT_ENDPOINTS
    raw = None
    last_err = None
    for ep in endpoints:
        try:
            raw = _post(ep, text_block)
            break
        except Exception as exc:  # noqa: BLE001
            last_err = exc
    if raw is None:
        raise SystemExit(f"All CRONUS endpoints failed: {last_err}")

    if args.raw_xml_out is not None:
        args.raw_xml_out.parent.mkdir(parents=True, exist_ok=True)
        args.raw_xml_out.write_text(raw)

    refs = _parse_xml(raw, args.scheme)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "sample_name",
            "nuclide",
            "latitude_deg",
            "longitude_deg",
            "effective_pressure_hpa",
            "lsdtt_erate_gcm2yr",
            "cronus_v3_erate_gcm2yr",
            "pct_diff_lsdtt_minus_cronus",
        ])
        for i, r in enumerate(rows, start=1):
            sid = f"S{i:05d}"
            sample, nuclide = idmap[sid]
            lsdtt_e = float(r["apparent_eff_e_g_cm2_yr"])
            cv3 = refs.get((sid, nuclide))
            pct = ""
            if cv3 is not None and cv3 != 0.0:
                pct = f"{100.0 * (lsdtt_e - cv3) / cv3:.8g}"
            w.writerow([
                sample,
                nuclide,
                r.get("centroid_latitude_deg") or r.get("latitude"),
                r.get("longitude"),
                r.get("effective_pressure_hpa"),
                f"{lsdtt_e:.10g}",
                "" if cv3 is None else f"{cv3:.10g}",
                pct,
            ])

    print(f"WROTE: {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
