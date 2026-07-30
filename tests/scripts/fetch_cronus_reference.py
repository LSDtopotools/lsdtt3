#!/usr/bin/env python3
"""Fetch the frozen CRONUS v3 erosion-rate reference for the lsdtt3 cosmo
comparison tests (see docs/testing_and_benchmarking.md).

Reads the fixture sample table written by the C++ test driver
(cronus_fixture_samples.csv), submits ONE batch request to the CRONUS v3
online calculator (the closed-source v3 accepts an HTTP POST of the v2/v3
plain-text input block - the same mechanism the open-source `riversand`
package uses), parses the XML response, and writes
tests/expected/cosmo/cronus_v3_reference.csv.

This is a MANUAL, one-off step: CI never calls the network. Refreshing the
reference is a reviewable event - commit the regenerated CSV and record why.
Be polite to the CRONUS server: this script sends a single batched request
and caches the raw XML next to the output.

Usage:
  python3 tests/scripts/fetch_cronus_reference.py \
      [--samples build/tests/cosmo_test_output/cronus_fixture_samples.csv] \
      [--out-dir tests/expected/cosmo] \
      [--endpoint https://hess.ess.washington.edu/cgi-bin/matweb]

Physics notes:
- Pressure is submitted directly ('pre' flag) using OUR standard-atmosphere
  pressures from the fixture CSV, so atmosphere-model differences (CRONUS
  'std' = NCEP) are excluded from the comparison by construction.
- Thickness 0, density 2.65, shielding 1, standards 07KNSTD / KNSTD.
- C-14 is included when `conc_C14_C14STD` is present in the samples CSV.
  CRONUS v3 expects C-14 lines without an AMS standard token.
- CRONUS reports erosion in m/Myr per scheme (St/Lm/LSDn) plus g/cm2/yr;
  we keep g/cm2/yr as primary (density-free, matches lsdtt3).
"""

import argparse
import csv
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

DEFAULT_ENDPOINTS = [
    "https://hess.ess.washington.edu/cgi-bin/matweb",
    "http://stoneage.hzdr.de/cgi/matweb",
]

# CRONUS v3 erosion output reports m/Myr; it also reports gcm2/yr in the
# long XML. Density used to build the input lines (thickness is 0 so it only
# matters for CRONUS's internal m/Myr conversion).
DENSITY = 2.65
UNC_FRACTION = 0.03  # 3 percent one-sigma on the synthetic concentrations


def read_samples(path: Path):
    samples = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            samples.append(row)
    if not samples:
        sys.exit(f"no samples found in {path}")
    return samples


def build_text_block(samples):
    """v3 input: sample line (10 fields) + one nuclide line per measurement."""
    lines = []
    has_c14 = "conc_C14_C14STD" in samples[0]
    for s in samples:
        name = s["name"]
        lines.append(
            f"{name} {float(s['lat']):.5f} {float(s['lon']):.5f} "
            f"{float(s['pressure_hPa']):.4f} pre 0 {DENSITY} 1 0.0001 2010 ;"
        )
        c10 = float(s["conc_Be10_07KNSTD"])
        c26 = float(s["conc_Al26_KNSTD"])
        lines.append(
            f"{name} Be-10 quartz {c10:.4f} {c10 * UNC_FRACTION:.4f} 07KNSTD ;"
        )
        lines.append(
            f"{name} Al-26 quartz {c26:.4f} {c26 * UNC_FRACTION:.4f} KNSTD ;"
        )
        if has_c14:
            c14 = float(s["conc_C14_C14STD"])
            lines.append(
                f"{name} C-14 quartz {c14:.4f} {c14 * UNC_FRACTION:.4f} ;"
            )
    return "\n".join(lines)


def post_request(endpoint: str, text_block: str) -> str:
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


def parse_xml(raw: str):
    """Yields (name, nuclide, scheme, erate_gcm2yr, unc_int, unc_ext).

    Observed v3 XML (2026-07-16):
      <calcs_v3_erosion_data>
        <erosionRateResult>
          <sample_name>..</sample_name>
          <nuclide_result>
            <nuclide>N10quartz</nuclide>
            <E_gcm2_St>..</E_gcm2_St>
            <delE_gcm2_St>internal</delE_gcm2_St>
            <delE_gcm2_St>external</delE_gcm2_St>
            ... _Lm, _LSDn ...
          </nuclide_result>
          ...
    Erosion rates are g/cm2/yr; the duplicated delE tags are internal then
    external one-sigma uncertainty.
    """
    start = raw.find("<")
    if start < 0:
        sys.exit("response contains no XML; raw response saved for inspection")
    root = ET.fromstring(raw[start:])

    rows = []
    for res in root.iter():
        if not res.tag.endswith("erosionRateResult"):
            continue
        name = (res.findtext("sample_name") or "").strip()
        for nres in res.findall("nuclide_result"):
            nuclide = (nres.findtext("nuclide") or "").strip()
            if "10" in nuclide:
                nuclide = "Be-10"
            elif "26" in nuclide:
                nuclide = "Al-26"
            elif "14" in nuclide:
                nuclide = "C-14"
            for scheme in ("St", "Lm", "LSDn"):
                e_text = nres.findtext(f"E_gcm2_{scheme}")
                if e_text is None:
                    continue
                dels = [
                    float(d.text)
                    for d in nres.findall(f"delE_gcm2_{scheme}")
                    if d.text
                ]
                unc_int = dels[0] if len(dels) > 0 else None
                unc_ext = dels[1] if len(dels) > 1 else None
                rows.append(
                    (name, nuclide, scheme, float(e_text), unc_int, unc_ext)
                )
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    repo = Path(__file__).resolve().parents[2]
    ap.add_argument(
        "--samples",
        type=Path,
        default=repo / "build/tests/cosmo_test_output/cronus_fixture_samples.csv",
    )
    ap.add_argument("--out-dir", type=Path, default=repo / "tests/expected/cosmo")
    ap.add_argument("--endpoint", action="append", default=None)
    args = ap.parse_args()

    samples = read_samples(args.samples)
    text_block = build_text_block(samples)
    print(f"submitting {len(samples)} samples (one batched request)")

    raw = None
    endpoints = args.endpoint or DEFAULT_ENDPOINTS
    for ep in endpoints:
        try:
            print(f"POST {ep} ...")
            raw = post_request(ep, text_block)
            break
        except Exception as exc:  # noqa: BLE001 - report and try the mirror
            print(f"  failed: {exc}")
    if raw is None:
        sys.exit("all endpoints failed")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.out_dir / "cronus_v3_raw.xml"
    raw_path.write_text(raw)
    print(f"raw response cached at {raw_path} ({len(raw)} bytes)")

    rows = parse_xml(raw)
    if not rows:
        sys.exit(
            "no erosion results parsed - inspect the cached raw XML; the "
            "server may have returned an error page"
        )

    out_path = args.out_dir / "cronus_v3_reference.csv"
    with open(out_path, "w", newline="") as f:
        f.write(
            "# CRONUS v3 online-calculator erosion rates for the lsdtt3\n"
            "# comparison fixture - see docs/testing_and_benchmarking.md\n"
            "# FROZEN reference - refreshing it is a reviewable event.\n"
        )
        w = csv.writer(f)
        w.writerow(
            ["name", "nuclide", "scheme", "erate_gcm2yr",
             "unc_int_gcm2yr", "unc_ext_gcm2yr"]
        )
        for r in rows:
            w.writerow(["" if v is None else v for v in r])
    print(f"wrote {len(rows)} reference rows to {out_path}")


if __name__ == "__main__":
    main()
