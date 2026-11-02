#!/usr/bin/env python3
"""Gather lsdtt-valley-centreline width tables from several DEMs.

Writes one table with every station of every DEM (a `dem` column first), and
one summary row per DEM: station counts by `excluded_reason`, and the
quartiles of `valley_width_m` over the kept stations (`excluded_reason` none),
and how many kept stations have a width of 0.
"""

import argparse
import csv
import statistics


def quantile(sorted_values, q):
    if not sorted_values:
        return ""
    if len(sorted_values) == 1:
        return f"{sorted_values[0]:.3f}"
    return f"{statistics.quantiles(sorted_values, n=4, method='inclusive')[q - 1]:.3f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--widths", nargs="+", required=True)
    ap.add_argument("--table-out", required=True)
    ap.add_argument("--summary-out", required=True)
    args = ap.parse_args()
    if len(args.names) != len(args.widths):
        raise SystemExit("--names and --widths must have the same length")

    reasons = set()
    per_dem = []
    header = None
    with open(args.table_out, "w", newline="", encoding="utf-8") as fout:
        writer = None
        for name, path in zip(args.names, args.widths):
            with open(path, newline="", encoding="utf-8") as fin:
                reader = csv.DictReader(fin)
                if writer is None:
                    header = ["dem"] + list(reader.fieldnames)
                    writer = csv.DictWriter(fout, fieldnames=header)
                    writer.writeheader()
                elif ["dem"] + list(reader.fieldnames) != header:
                    raise SystemExit(f"{path}: columns differ from the first table")
                counts = {}
                kept = []
                segments = set()
                for row in reader:
                    reason = row["excluded_reason"]
                    counts[reason] = counts.get(reason, 0) + 1
                    reasons.add(reason)
                    if reason == "none":
                        kept.append(float(row["valley_width_m"]))
                        segments.add(row["segment_id"])
                    writer.writerow({"dem": name, **row})
            per_dem.append((name, counts, sorted(kept), len(segments)))

    excluded = sorted(r for r in reasons if r != "none")
    with open(args.summary_out, "w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout)
        writer.writerow(["dem", "n_stations", "n_kept", *[f"n_{r}" for r in excluded],
                         "n_kept_zero_width", "n_segments_with_kept_stations",
                         "kept_width_p25_m", "kept_width_median_m", "kept_width_p75_m"])
        for name, counts, kept, n_segments in per_dem:
            writer.writerow([name, sum(counts.values()), counts.get("none", 0),
                             *[counts.get(r, 0) for r in excluded],
                             sum(1 for w in kept if w == 0.0), n_segments,
                             quantile(kept, 1), quantile(kept, 2), quantile(kept, 3)])


if __name__ == "__main__":
    main()
