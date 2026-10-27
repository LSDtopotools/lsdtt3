#!/usr/bin/env python3
"""Check segmentation_algorithm_choice whole_mask against within_bands.

Fluvial-surface methods plan, phase 6. whole_mask is the riser segmentation:
the admitted ground minus a riser mask is labelled once, each patch taking its
majority level. It is kept so it can be compared with band-then-label, and the
comparison is the test. On the dissected staircase (4 treads, 3 graded gaps,
known tread map), admitted by gradient so the gap ramps connect the treads,
banded on the generator's own levels so only the segmentation differs, and
scored on tread cells (the tread map's id 0 - gaps, wall, channel - is not a
tread):

  * within_bands: every patch is at least 95 % one tread (the baseline);
  * whole_mask with a CLOSED riser mask - every boundary of the tread map,
    gaps included - does as well: no patch spans two treads;
  * whole_mask with the find_risers centreline (lsdtt-feature-finder, no gap
    closing), which stops at each graded gap: treads join through the gaps, so
    at least one patch holds two treads with at least 20 % of its cells each.
    This is the documented failure of riser segmentation, asserted so the
    option's behaviour stays honest.
"""

from __future__ import annotations

import argparse
import pathlib
import struct
import sys
from collections import Counter, defaultdict

import lsdtt_test_env


def run(exe: str, template: pathlib.Path, read_path: pathlib.Path,
        out_dir: pathlib.Path, work_dir: pathlib.Path, tag: str,
        extra: dict[str, str]) -> str:
    """Render a template, run the CLI on it, and return its stdout."""
    text = (template.read_text(encoding="utf-8")
            .replace("__READ_PATH__", str(read_path))
            .replace("__WRITE_PATH__", str(out_dir)))
    for k, v in extra.items():
        text = text.replace(k, v)
    param = work_dir / f"{tag}.param"
    param.write_text(text, encoding="utf-8")
    proc = lsdtt_test_env.run([exe, str(work_dir), param.name],
                              cwd=str(work_dir), capture_output=True,
                              text=True, check=False)
    (work_dir / f"stdout_{tag}.log").write_text(proc.stdout, encoding="utf-8")
    (work_dir / f"stderr_{tag}.log").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout[-3000:] + proc.stderr[-3000:])
        raise SystemExit(f"run {tag} exited {proc.returncode}")
    return proc.stdout


def read_raster(path: pathlib.Path):
    """Band 1 as a flat list of floats, plus the dataset (ReadRaster, not
    ReadAsArray: gdal_array is tied to the NumPy it was built against)."""
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(str(path))
    w, h = ds.RasterXSize, ds.RasterYSize
    raw = ds.GetRasterBand(1).ReadRaster(0, 0, w, h, buf_type=gdal.GDT_Float32)
    return list(struct.unpack(f"<{w * h}f", raw)), ds


def write_closed_riser_mask(tread_path: pathlib.Path, out: pathlib.Path) -> int:
    """Every cell whose tread id differs from its east or south neighbour: the
    tread map's boundaries. The gaps carry tread id 0, so this also cuts them
    off the treads up- and down-valley, closing every tread."""
    from osgeo import gdal
    tread, ds = read_raster(tread_path)
    w, h = ds.RasterXSize, ds.RasterYSize
    mask = [0.0] * (w * h)
    for r in range(h):
        for c in range(w):
            i = r * w + c
            if c + 1 < w and int(tread[i]) != int(tread[i + 1]):
                mask[i] = 1.0
            if r + 1 < h and int(tread[i]) != int(tread[i + w]):
                mask[i] = 1.0
    drv = gdal.GetDriverByName("GTiff")
    o = drv.Create(str(out), w, h, 1, gdal.GDT_Float32)
    o.SetGeoTransform(ds.GetGeoTransform())
    o.SetProjection(ds.GetProjection())
    o.GetRasterBand(1).WriteRaster(0, 0, w, h, struct.pack(f"<{w * h}f", *mask))
    o.FlushCache()
    return int(sum(mask))


def patch_treads(label_path: pathlib.Path, tread_path: pathlib.Path):
    """Per patch label: Counter of tread ids over its tread cells. Id 0 is
    the gaps, the valley wall and the channel - not a tread - and is not
    scored."""
    lab, _ = read_raster(label_path)
    tread, _ = read_raster(tread_path)
    per = defaultdict(Counter)
    for l, t in zip(lab, tread):
        if l > 0.5 and int(t) > 0:
            per[int(l)][int(t)] += 1
    return per


def main() -> int:
    lsdtt_test_env.require_modules(["osgeo.gdal"])
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaces-exe", required=True)
    ap.add_argument("--feature-finder-exe", required=True)
    ap.add_argument("--param-dir", required=True)
    ap.add_argument("--read-path", required=True)
    ap.add_argument("--work-dir", required=True)
    args = ap.parse_args()

    params = pathlib.Path(args.param_dir)
    read_path = pathlib.Path(args.read_path)
    work = pathlib.Path(args.work_dir)
    out = work / "out"
    out.mkdir(parents=True, exist_ok=True)
    tread_path = read_path / "staircase_tread_id.tif"

    n_closed = write_closed_riser_mask(tread_path, out / "closed_risers.tif")
    print(f"  closed riser mask: {n_closed} cells")
    run(args.feature_finder_exe, params / "feature_finder_risers_open.param.in",
        read_path, out, work, "risers", {})

    template = params / "fluvial_surfaces_whole_mask.param.in"
    cases = {
        "wmb": ("within_bands", ""),
        "wmc": ("whole_mask", f"riser_mask_fname: {out / 'closed_risers.tif'}"),
        "wmo": ("whole_mask", f"riser_mask_fname: {out / 'wmr_riser_centreline.tif'}"),
    }
    logs = {}
    for stem, (seg, riser) in cases.items():
        logs[stem] = run(args.surfaces_exe, template, read_path, out, work, stem,
                         {"__SEGMENTATION__": seg, "__RISER_MASK__": riser,
                          "__STEM__": stem})

    failures = []

    def check(ok: bool, what: str) -> None:
        print(("  ok    " if ok else "  FAIL  ") + what)
        if not ok:
            failures.append(what)

    def summary(stem: str):
        per = patch_treads(out / f"{stem}_banded_labels.tif", tread_path)
        impure = 0
        joined = 0
        for counts in per.values():
            n = sum(counts.values())
            if n == 0:
                continue
            top = counts.most_common(2)
            if top[0][1] / n < 0.95:
                impure += 1
            if len(top) > 1 and top[1][1] / n >= 0.20:
                joined += 1
        print(f"  {stem}: {len(per)} patches, {impure} below 95 % one tread, "
              f"{joined} with a second tread >= 20 %")
        return len(per), impure, joined

    _, impure_b, _ = summary("wmb")
    _, impure_c, _ = summary("wmc")
    _, _, joined_o = summary("wmo")

    check(impure_b == 0, "within_bands: every patch at least 95 % one tread")
    check(impure_c == 0,
          "whole_mask, closed risers: every patch at least 95 % one tread")
    check("whole_mask:" in logs["wmc"] and "removed as riser" in logs["wmc"],
          "whole_mask reports the riser cells it removed")
    check(joined_o >= 1,
          "whole_mask, risers open at the gaps: treads join into one patch "
          "(the documented failure)")

    if failures:
        print(f"FAILED: {len(failures)} check(s)")
        return 1
    print("all whole_mask segmentation checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
