#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""lsdtt-fetch-river-lines -> lsdtt-stream-burning, end to end (plan phase 1b).

Four stream-burning runs on the stream-burning fixture DEM:

  raw_clean    the committed network as it is (the reference);
  prep_clean   the same network after lsdtt-fetch-river-lines;
  prep_messy   a deliberately messy copy after lsdtt-fetch-river-lines: every
               line reversed, one junction end moved 0.6 m off its junction,
               and two lines merged into one multipart feature;
  raw_messy    the messy copy given to stream burning directly.

What must hold: every run's vector network summary (junctions, edges,
outlets, loops, sources, sources mapped, edge sources carved) equals the
reference's, and every burned DEM has the sha256 given by
--expect-burned-sha256 (the frozen oracle of smoke_stream_burning). So
preparing the lines, clean or messy, does not change what stream burning
does. raw_messy is reported, and must also match: stream burning's own
snapping copes with this mess, so this fixture pins compatibility, not an
improvement.
"""

import argparse
import hashlib
import pathlib
import re
import shutil
import sys

import lsdtt_test_env

SUMMARY_KEYS = ("Vector junctions", "Vector edges", "Vector outlets",
                "Loops detected / resolved", "Vector source junctions",
                "Source junctions mapped", "Edge sources carved inward")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--fetch-exe", required=True)
    p.add_argument("--burn-exe", required=True)
    p.add_argument("--burn-param-template", required=True)
    p.add_argument("--fixture-dir", required=True,
                   help="stream-burning fixture: sb_diamond_square.tif + stream_burning_network.fgb")
    p.add_argument("--work-dir", required=True)
    p.add_argument("--expect-burned-sha256", required=True)
    return p.parse_args()


def fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def run(exe: str, work: pathlib.Path, param_name: str, tag: str) -> str:
    r = lsdtt_test_env.run([exe, str(work), param_name], cwd=str(work),
                           capture_output=True, text=True, check=False)
    (work / f"stdout_{tag}.log").write_text(r.stdout, encoding="utf-8")
    (work / f"stderr_{tag}.log").write_text(r.stderr, encoding="utf-8")
    if r.returncode != 0:
        fail(f"{tag}: {exe} returned {r.returncode}; see {work / f'stderr_{tag}.log'}")
    return r.stdout


def summary(stdout: str) -> dict:
    out = {}
    for key in SUMMARY_KEYS:
        m = re.search(r"^" + re.escape(key) + r"\s*:\s*(.+)$", stdout, re.M)
        out[key] = re.sub(r"\s+", " ", m.group(1)).strip() if m else None
    return out


def make_messy(clean: pathlib.Path, messy: pathlib.Path) -> None:
    from osgeo import ogr
    src = ogr.Open(str(clean))
    layer = src.GetLayer(0)
    pts = {}
    for f in layer:
        g = f.GetGeometryRef()
        pts[f.GetField("id")] = [g.GetPoint_2D(k) for k in range(g.GetPointCount())]
    if sorted(pts) != [1, 2, 3, 4, 5]:
        fail(f"unexpected fixture ids {sorted(pts)}")
    out = ogr.GetDriverByName("FlatGeobuf").CreateDataSource(str(messy))
    lo = out.CreateLayer("stream_network", layer.GetSpatialRef(), ogr.wkbUnknown)
    lo.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    def line(p):
        g = ogr.Geometry(ogr.wkbLineString)
        for x, y in p:
            g.AddPoint_2D(x, y)
        return g

    def put(i, geom):
        f = ogr.Feature(lo.GetLayerDefn())
        f.SetField("id", i)
        f.SetGeometry(geom)
        lo.CreateFeature(f)

    for i in (1, 2, 5):
        p = pts[i][::-1]
        if i == 2:
            # line 2 ends where line 5 ends; reversed, that end comes first
            p[0] = (p[0][0] + 0.6, p[0][1])
        put(i, line(p))
    m = ogr.Geometry(ogr.wkbMultiLineString)
    m.AddGeometry(line(pts[3][::-1]))
    m.AddGeometry(line(pts[4][::-1]))
    put(34, m)
    out = None


def main() -> int:
    a = parse_args()
    work = pathlib.Path(a.work_dir)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    fixture = pathlib.Path(a.fixture_dir)
    shutil.copy(fixture / "sb_diamond_square.tif", work)
    shutil.copy(fixture / "stream_burning_network.fgb", work / "clean.fgb")
    make_messy(work / "clean.fgb", work / "messy.fgb")

    for v in ("clean", "messy"):
        (work / f"prep_{v}.param").write_text(
            "read_path: ./\nwrite_path: ./\n"
            f"write_fname: prep_{v}\n"
            "river_line_source: local_file\n"
            f"river_lines_fname: {v}.fgb\n"
            "bounds_fname: sb_diamond_square.tif\n", encoding="utf-8")
        out = run(a.fetch_exe, work, f"prep_{v}.param", f"prep_{v}")
        if f"stream_vector_path: prep_{v}_river_lines.fgb" not in out:
            fail(f"prep_{v}: no stream_vector_path line in the log")

    template = pathlib.Path(a.burn_param_template).read_text(encoding="utf-8")
    runs = {"raw_clean": "clean.fgb", "prep_clean": "prep_clean_river_lines.fgb",
            "prep_messy": "prep_messy_river_lines.fgb", "raw_messy": "messy.fgb"}
    results = {}
    for tag, vec in runs.items():
        text = (template.replace("__READ_PATH__/stream_burning_network.fgb", vec)
                .replace("__READ_PATH__", "./").replace("__WRITE_PATH__", "./"))
        text = re.sub(r"^write_fname: .*$", f"write_fname: {tag}", text, flags=re.M)
        (work / f"{tag}.param").write_text(text, encoding="utf-8")
        s = summary(run(a.burn_exe, work, f"{tag}.param", tag))
        sha = hashlib.sha256((work / f"{tag}_burned.tif").read_bytes()).hexdigest()
        results[tag] = (s, sha)
        print(f"{tag:11s} sha256 {sha[:16]}  " + "; ".join(f"{k}: {v}" for k, v in s.items()))

    ref = results["raw_clean"][0]
    if any(v is None for v in ref.values()):
        fail(f"reference summary incomplete: {ref}")
    ok = True
    for tag, (s, sha) in results.items():
        if s != ref:
            print(f"FAIL: {tag} summary differs from raw_clean:\n  {s}\n  {ref}", file=sys.stderr)
            ok = False
        if sha != a.expect_burned_sha256:
            print(f"FAIL: {tag} burned DEM sha256 {sha} != {a.expect_burned_sha256}", file=sys.stderr)
            ok = False
    if not ok:
        return 1
    print("PASS: prepared clean and messy networks give the reference stream-burning result")
    return 0


if __name__ == "__main__":
    sys.exit(main())
