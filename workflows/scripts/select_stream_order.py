#!/usr/bin/env python3

import argparse
import csv
import pathlib
import re
import subprocess
import sys


def run_for_order(exe: pathlib.Path, dem_path: pathlib.Path, work_dir: pathlib.Path,
                  source_threshold: int, stream_order: int) -> int:
    param_path = work_dir / f"probe_order_{stream_order}.param"
    text = (
        f"read_path: {dem_path.parent}\n"
        f"read_fname: {dem_path.name}\n"
        f"write_path: {work_dir}\n"
        f"write_fname: probe_order_{stream_order}\n\n"
        "raster_is_filled: true\n"
        f"source_pixels_threshold: {source_threshold}\n\n"
        "select_all_basins_of_stream_order: true\n"
        f"selected_stream_order: {stream_order}\n\n"
        "write_basin_outlines: false\n"
        "write_channel_network_lines: false\n"
        "write_channel_network_points: false\n"
        "write_channel_network_junctions: false\n"
        "write_basins_raster: false\n"
        "calculate_spatial_index: false\n"
    )
    param_path.write_text(text, encoding="utf-8")

    proc = subprocess.run(
        [str(exe), str(work_dir), param_path.name],
        cwd=str(work_dir),
        capture_output=True,
        text=True,
        check=False,
    )

    # Empty selection currently exits nonzero ("basin selection produced no basins").
    if proc.returncode != 0:
        if "produced no basins" in proc.stderr:
            return 0
        raise RuntimeError(
            f"probe order {stream_order} failed with exit {proc.returncode}\n"
            f"--- stdout ---\n{proc.stdout[-4000:]}\n"
            f"--- stderr ---\n{proc.stderr[-4000:]}"
        )

    m = re.search(r"Selected basins:\s*([0-9]+)", proc.stdout)
    if not m:
        raise RuntimeError(
            f"could not parse selected-basin count for stream order {stream_order}\n"
            f"--- stdout ---\n{proc.stdout[-4000:]}"
        )
    return int(m.group(1))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--dem", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--source-threshold", type=int, required=True)
    ap.add_argument("--min-basins", type=int, required=True)
    ap.add_argument("--max-order", type=int, required=True)
    ap.add_argument("--selected-order-out", required=True)
    ap.add_argument("--scan-report-out", required=True)
    args = ap.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    dem = pathlib.Path(args.dem).resolve()
    work_dir = pathlib.Path(args.work_dir).resolve()
    selected_out = pathlib.Path(args.selected_order_out).resolve()
    report_out = pathlib.Path(args.scan_report_out).resolve()

    if not exe.is_file():
        raise SystemExit(f"lsdtt-basin-extraction executable not found: {exe}")
    if not dem.is_file():
        raise SystemExit(f"DEM not found: {dem}")

    work_dir.mkdir(parents=True, exist_ok=True)
    report_out.parent.mkdir(parents=True, exist_ok=True)
    selected_out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    chosen = None
    for order in range(args.max_order, 0, -1):
        n_basins = run_for_order(exe, dem, work_dir, args.source_threshold, order)
        rows.append((order, n_basins))
        if n_basins >= args.min_basins:
            chosen = order
            break

    with report_out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["stream_order", "n_basins"])
        for order, n_basins in rows:
            w.writerow([order, n_basins])

    if chosen is None:
        raise SystemExit(
            f"no stream order in 1..{args.max_order} had at least {args.min_basins} basins"
        )

    selected_out.write_text(f"{chosen}\n", encoding="utf-8")
    print(
        f"Selected stream order {chosen} (largest with >= {args.min_basins} basins)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
