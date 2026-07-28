#!/usr/bin/env python3
"""Summarize muon-isolation decomposition CSVs for quick review diffs.

Reads `muon_isolation_comparison.csv` and prints stable, grouped stats with
explicit solve-status counts and split residual summaries:

- all solved rows (`ok` + bound-limited)
- interior rows only (`solve_status == ok`)

This is reporting only; it does not modify files.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def _fmt(x: float) -> str:
    return format(x, ".6g")


def _stats(values: list[float]) -> tuple[float, float, float]:
    if not values:
        return (math.nan, math.nan, math.nan)
    return (min(values), max(values), sum(values) / float(len(values)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "csv_path",
        nargs="?",
        default="build/tests/cosmo_test_output/muon_isolation_comparison.csv",
        help="Path to muon_isolation_comparison.csv",
    )
    args = ap.parse_args()

    path = Path(args.csv_path)
    if not path.is_file():
        raise FileNotFoundError(path)

    groups: dict[tuple[str, str], dict[str, object]] = {}
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            variant = row["variant"]
            nuclide = row["nuclide"]
            status = row.get("solve_status", "") or "unknown"

            key = (variant, nuclide)
            g = groups.setdefault(
                key,
                {
                    "n_total": 0,
                    "n_ok": 0,
                    "n_low": 0,
                    "n_high": 0,
                    "n_none": 0,
                    "all_pct": [],
                    "ok_pct": [],
                },
            )

            g["n_total"] = int(g["n_total"]) + 1
            if status == "ok":
                g["n_ok"] = int(g["n_ok"]) + 1
            elif status == "bound_limited_low":
                g["n_low"] = int(g["n_low"]) + 1
            elif status == "bound_limited_high":
                g["n_high"] = int(g["n_high"]) + 1
            elif status == "no_solution":
                g["n_none"] = int(g["n_none"]) + 1

            try:
                pct = float(row["pct_diff"])
            except ValueError:
                continue
            if not math.isfinite(pct):
                continue

            cast_all = g["all_pct"]
            assert isinstance(cast_all, list)
            cast_all.append(pct)
            if status == "ok":
                cast_ok = g["ok_pct"]
                assert isinstance(cast_ok, list)
                cast_ok.append(pct)

    print(f"source: {path}")
    print(
        "variant,nuclide,n_total,n_ok,n_low,n_high,n_none,"
        "all_min,all_max,all_mean,ok_min,ok_max,ok_mean"
    )

    for key in sorted(groups.keys()):
        variant, nuclide = key
        g = groups[key]
        all_min, all_max, all_mean = _stats(g["all_pct"])  # type: ignore[arg-type]
        ok_min, ok_max, ok_mean = _stats(g["ok_pct"])      # type: ignore[arg-type]
        print(
            f"{variant},{nuclide},"
            f"{g['n_total']},{g['n_ok']},{g['n_low']},{g['n_high']},{g['n_none']},"
            f"{_fmt(all_min)},{_fmt(all_max)},{_fmt(all_mean)},"
            f"{_fmt(ok_min)},{_fmt(ok_max)},{_fmt(ok_mean)}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
