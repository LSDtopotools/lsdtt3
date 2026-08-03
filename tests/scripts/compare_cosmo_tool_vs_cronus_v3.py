#!/usr/bin/env python3
"""Compare lsdtt-cosmo-tool sample erosion rates against CRONUS v3.

Reads ``<write_fname>_CN_ERATE.csv`` from lsdtt-cosmo-tool, builds a CRONUS v3
batch request using per-sample nuclide, concentration, latitude and effective
pressure, submits once, and writes a comparison CSV.

WARNING -- AMS standardisation double-count with --standardisation-col.
The ``concentration_*`` column in ``_CN_ERATE.csv`` is ALREADY normalised by
lsdtt-cosmo-tool (e.g. KNSTD -> 07KNSTD, x0.9042).  Passing
``--standardisation-col`` sends that normalised value together with the sample's
standard name, so CRONUS applies the factor a SECOND time.  This is the error
that inverted the apparent sign of the OCTOPUS-tier residual (clapp2002 read
+6.6 % rather than -4.09 %); see docs/testing_and_benchmarking.md.

Omit ``--standardisation-col`` unless every sample is already on the calculator's
default standard, in which case it is a no-op.  For a correct request built from
as-measured concentrations, prefer the file lsdtt-cosmo-tool now writes itself:
``<write_fname>_CRONUS_V3_INPUT.txt`` (parameter ``write_cronus_v3_input``,
default true).
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
    # Writing the request out is useful on its own: it is exactly the text that
    # the CRONUS v3 erosion-rate input page accepts, so a user can paste it into
    # the web form instead of submitting from here (no network, and the online
    # calculator's own reporting).  With --no-submit it is the only output, so
    # --out-csv stops being required.
    ap.add_argument("--text-block-out", type=Path, default=None,
                    help="Write the CRONUS v3 input text block to this file.")
    ap.add_argument("--no-submit", action="store_true",
                    help="Build the input text block but do not contact CRONUS. "
                         "Requires --text-block-out.")
    ap.add_argument("--out-csv", type=Path, default=None)
    args = ap.parse_args()

    if args.no_submit and args.text_block_out is None:
        raise SystemExit("--no-submit requires --text-block-out")
    if not args.no_submit and args.out_csv is None:
        raise SystemExit("--out-csv is required unless --no-submit is given")

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

    # The docstring warning only reaches someone who runs --help; this reaches
    # the person who is actually about to submit a double-counted request.
    # Warn on ANY non-empty standard rather than trying to exclude "default"
    # ones: the default is nuclide-dependent (07KNSTD for Be-10, KNSTD for
    # Al-26), so KNSTD is a no-op on Al-26 but a 0.9042 double-count on Be-10.
    # A spurious warning costs a line of output; a missed one corrupts a rate.
    if args.standardisation_col and points_by_name:
        seen_standards = {
            r.get(args.standardisation_col, "")
            for r in points_by_name.values()
            if r.get(args.standardisation_col, "")
        }
        if seen_standards:
            print("WARNING: --standardisation-col double-counts AMS standardisation.\n"
                  "         concentration_* in the cosmo CSV is already normalised by\n"
                  "         lsdtt-cosmo-tool; sending it with a standard name applies the\n"
                  "         factor twice. Standards seen: "
                  + ", ".join(sorted(seen_standards)) + "\n"
                  "         Prefer <write_fname>_CRONUS_V3_INPUT.txt from the CLI.")

    text_block, idmap = _build_text_block(
        rows,
        points_by_name,
        args.uncertainty_col,
        args.standardisation_col,
    )

    if args.text_block_out is not None:
        args.text_block_out.parent.mkdir(parents=True, exist_ok=True)
        args.text_block_out.write_text(text_block + "\n")
        print(f"Wrote CRONUS v3 input text block: {args.text_block_out} "
              f"({len(rows)} sample(s))")

    if args.no_submit:
        return 0

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
