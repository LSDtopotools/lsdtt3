#!/usr/bin/env python3
"""Phase 5 sweep: every CLI that ingests vector data, against multipart input.

Until 2026-08-26 lsdtt3 had no branch for ANY multipart geometry and discarded
it silently, with a zero exit code. The tables in
docs/plan_docs/plan_multilinestring_channel_outputs.md audit the CALL SITES;
this script audits the BEHAVIOUR, which is not the same thing.

THE FAILURE MODE IS SILENCE, so every case here asserts a positive result --
a feature count, a changed raster, a non-empty table -- and never a clean exit.
Both encodings of every fixture exited 0 before the fix.

Each case runs the same CLI twice on geometrically identical input, once
single-part and once multipart, and requires the multipart run to produce a
result comparable to the single-part one.

Usage: multipart_sweep.py <build_dir> <work_dir>
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import lsdtt_test_env

# Before the third-party imports, so a missing package is a clear message.
lsdtt_test_env.require_modules(["geopandas", "numpy", "rasterio", "shapely"])

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import (LineString, MultiLineString, MultiPoint,
                              MultiPolygon, Point, Polygon)

NR = NC = 80
CELL = 10.0

results: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str) -> None:
    results.append((name, ok, detail))
    print(f"  {'ok  ' if ok else 'FAIL'}: {name} -- {detail}")


def build_fixtures(work: Path) -> None:
    r, c = np.indices((NR, NC))
    z = (200.0 - 1.2 * r + 0.5 * np.abs(c - NC / 2) + 0.01 * c).astype("float32")
    with rasterio.open(work / "dem.tif", "w", driver="GTiff", height=NR, width=NC,
                       count=1, dtype="float32", crs="EPSG:32630",
                       transform=from_origin(0, NR * CELL, CELL, CELL),
                       nodata=-9999) as d:
        d.write(z, 1)

    def X(col): return (col + 0.5) * CELL
    def Y(row): return (NR - row - 0.5) * CELL

    a = [(X(NC // 2), Y(i)) for i in range(5, 40)]
    b = [(X(NC // 2), Y(i)) for i in range(42, 75)]
    pts = [(X(20), Y(20)), (X(40), Y(40)), (X(60), Y(60))]
    sq1 = Polygon([(X(10), Y(10)), (X(30), Y(10)), (X(30), Y(30)), (X(10), Y(30))])
    sq2 = Polygon([(X(50), Y(50)), (X(70), Y(50)), (X(70), Y(70)), (X(50), Y(70))])

    out = {
        # Geometrically identical content, two encodings.
        "line_single": ([LineString(a + b)], "FlatGeobuf", "fgb"),
        "line_multi":  ([MultiLineString([a, b])], "FlatGeobuf", "fgb"),
        # Shapefile is the only format that preserves LineString next to
        # MultiLineString in one layer, so it is where a mixed layer is built.
        # lsdtt3 never writes shapefile; this is purely ingest.
        # The multipart half must have TWO parts. A one-part MultiLineString
        # round-trips through shapefile as a plain LineString -- the format has
        # no distinction -- so a fixture built that way is not mixed at all and
        # passes with or without the fix.
        "line_mixed":  ([LineString(a),
                         MultiLineString([b, [(X(10), Y(10)), (X(10), Y(35))]])],
                        "ESRI Shapefile", "shp"),
        "pt_single":   ([Point(p) for p in pts], "FlatGeobuf", "fgb"),
        "pt_multi":    ([MultiPoint(pts)], "FlatGeobuf", "fgb"),
        "poly_single": ([sq1, sq2], "FlatGeobuf", "fgb"),
        "poly_multi":  ([MultiPolygon([sq1, sq2])], "FlatGeobuf", "fgb"),
    }
    for name, (geoms, drv, ext) in out.items():
        gpd.GeoDataFrame({"id": range(len(geoms)), "geometry": geoms},
                         crs="EPSG:32630").to_file(work / f"{name}.{ext}", driver=drv)


def run(build: Path, work: Path, tool: str, params: str, tag: str) -> int:
    (work / f"{tag}.param").write_text(params)
    with open(work / f"{tag}.log", "w") as log:
        return lsdtt_test_env.run([str(build / tool), "./", f"{tag}.param"],
                              cwd=work, stdout=log, stderr=subprocess.STDOUT).returncode


def n_features(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    try:
        return len(gpd.read_file(path))
    except Exception:
        return 0


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

def case_linear_stats(build, work):
    got = {}
    for k in ("single", "multi", "mixed"):
        ext = "shp" if k == "mixed" else "fgb"
        run(build, work, "lsdtt-linear-stats",
            f"read_path: ./\nread_fname: dem.tif\nwrite_path: ./\n"
            f"write_fname: ls_{k}\nvector_fname: line_{k}.{ext}\n", f"ls_{k}")
        got[k] = n_features(work / f"ls_{k}_linear_stats.fgb")
    # mixed: 1 single-part + 1 two-part feature -> 3 records after exploding.
    record("lsdtt-linear-stats (line)",
           got["multi"] >= 2 and got["single"] >= 1 and got["mixed"] >= 3,
           f"single={got['single']} multi={got['multi']} mixed={got['mixed']}")


def case_stream_burning(build, work):
    got = {}
    for k in ("single", "multi"):
        run(build, work, "lsdtt-stream-burning",
            f"read_path: ./\nread_fname: dem.tif\nwrite_path: ./\n"
            f"write_fname: sb_{k}\nwrite_output_format: gtiff\n"
            f"stream_vector_path: line_{k}.fgb\nstream_burn_mode: burn_in_place\n"
            f"run_stream_burn: true\n", f"sb_{k}")
        # `_vector_edges.fgb` is the WRONG probe: it is produced by
        # lsdtt_vector_to_junction_network, which has handled MultiLineString
        # since before this work (plan section 2.2), so it passes either way and
        # proves nothing about lsdtt_burn_streams_cpu.
        #
        # `_burned_profiles.csv` comes from the burn itself, which is what
        # exercises the seven get_if<lsdtt_line_geom> sites.
        prof = work / f"sb_{k}_burned_profiles.csv"
        got[k] = sum(1 for _ in prof.open()) - 1 if prof.exists() else 0
    record("lsdtt-stream-burning (burn itself, 7 get_if sites)",
           got["multi"] > 0 and got["single"] > 0,
           f"burned profile rows: single={got['single']} multi={got['multi']}")


def case_channel_extraction(build, work):
    got = {}
    for k in ("single", "multi"):
        run(build, work, "lsdtt-channel-extraction",
            f"read_path: ./\nread_fname: dem.tif\nwrite_path: ./\n"
            f"write_fname: ce_{k}\nwrite_output_format: gtiff\n"
            f"source_extraction_algorithm_choice: threshold\n"
            f"threshold_contributing_pixels: 60\n"
            f"write_channel_network_points: true\n"
            f"channel_network_points_write_vector_format: flatgeobuf\n"
            f"stream_vector_path: line_{k}.fgb\ncalculate_spatial_index: false\n"
            f"write_burned_stream_profiles: true\n",
            f"ce_{k}")
        # The channel POINTS come from the drainage-area threshold, not from
        # the supplied vector, so counting them proves nothing about the vector
        # ingest -- they were identical before and after the fix. The burned
        # network is what the vector drives.
        # write_burned_stream_profiles is what routes the supplied vector
        # through lsdtt_burn_streams_cpu here; without it stream_vector_path is
        # read but never used, and the case exercises nothing.
        prof = list(work.glob(f"ce_{k}*burned*"))
        got[k] = sum(1 for _ in prof[0].open()) - 1 if prof else 0
    record("lsdtt-channel-extraction (stream_vector_path)",
           got["multi"] > 0 and got["single"] > 0,
           f"burned profile rows: single={got['single']} multi={got['multi']}")


def case_raster_preprocessing_points(build, work):
    (work / "rasters.csv").write_text("raster,attribute\ndem.tif,elev\n")
    got = {}
    for k in ("single", "multi"):
        run(build, work, "lsdtt-raster-preprocessing",
            f"read_path: ./\nread_fname: dem.tif\nwrite_path: ./\n"
            f"write_fname: rp_{k}\npoint_data_fname: pt_{k}.fgb\n"
            f"sample_rasters_to_points: true\nraster_attribute_fname: rasters.csv\n",
            f"rp_{k}")
        got[k] = n_features(work / f"rp_{k}_sampled_points.fgb")
    record("lsdtt-raster-preprocessing (point_data_fname)",
           got["multi"] == got["single"] == 3,
           f"sampled points: single={got['single']} multi={got['multi']}")


def case_swath_profile(build, work):
    got = {}
    for k in ("single", "multi"):
        run(build, work, "lsdtt-swath-profile",
            f"read_path: ./\nread_fname: dem.tif\nwrite_path: ./\n"
            f"write_fname: sw_{k}\nbaseline_vector_fname: line_{k}.fgb\n",
            f"sw_{k}")
        hits = list(work.glob(f"sw_{k}*.csv")) + list(work.glob(f"sw_{k}*.fgb"))
        got[k] = sum(1 for h in hits if h.stat().st_size > 0)
    record("lsdtt-swath-profile (baseline_vector_fname)",
           got["multi"] > 0 and got["single"] > 0,
           f"non-empty outputs: single={got['single']} multi={got['multi']}")


CASES = [
    case_linear_stats,
    case_stream_burning,
    case_channel_extraction,
    case_raster_preprocessing_points,
    case_swath_profile,
]


def main() -> int:
    if len(sys.argv) < 3:
        print("usage: multipart_sweep.py <build_dir> <work_dir>", file=sys.stderr)
        return 2
    build = Path(sys.argv[1]).resolve()
    work = Path(sys.argv[2]).resolve()
    work.mkdir(parents=True, exist_ok=True)

    print("=" * 55)
    print(" Phase 5 sweep: CLIs against multipart geometry")
    print("=" * 55)
    build_fixtures(work)

    for case in CASES:
        try:
            case(build, work)
        except Exception as exc:  # noqa: BLE001
            record(case.__name__, False, f"raised {type(exc).__name__}: {exc}")

    n_ok = sum(1 for _, ok, _ in results if ok)
    print("=" * 55)
    print(f" {n_ok} / {len(results)} cases passed")
    print("=" * 55)
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
