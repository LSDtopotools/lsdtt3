#!/usr/bin/env python3
"""Relational check on the band-then-label outputs of lsdtt-fluvial-surfaces.

Phase 1 of the fluvial-surfaces split makes the patch
label the join key for everything downstream: _banded_labels.tif carries it per
pixel, _band_index.tif carries the band, and _surface_patches.csv is the table
those two are meant to join to. The band-packed
label id specifically because that raster plus that table make packing
unnecessary - which is only true if the join actually holds. Nothing in a
per-file schema or checksum assertion can see a join, so this script exists.

It runs the CLI once and then asserts, across the emitted rasters and tables:

  1. the labels in _banded_labels.tif are a contiguous run 1..N and are exactly
     the labels in _surface_patches.csv, one row each;
  2. the label -> band map read off the two rasters is SINGLE-VALUED (no label
     appears in two bands) and agrees with the band column of the patch table;
  3. patch n_pixels equals the pixel count of that label in the raster;
  4. _surface_class.csv is 1:1 with the patch table on label, and the class,
     source_id and window centre columns agree row for row;
  5. every row of _reference_field_profile.csv joins to a patch by label, and its band and
     class columns agree with the rasters at that row's own coordinates.

With --expect-empty-tables the assertion is inverted: the run must emit zero
patches and still write both CSVs with their header rows, so a degenerate
result is an empty table rather than a missing file.
"""

import argparse
import csv
import math
import pathlib
import shutil
import struct
import subprocess
import sys

import lsdtt_test_env

# Before the osgeo import, so a missing package is a clear message.
lsdtt_test_env.require_modules(["osgeo.gdal"])

from osgeo import gdal


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify the band-then-label output joins")
    parser.add_argument("--exe", required=True, help="Path to the CLI executable")
    parser.add_argument("--param-template", required=True,
                        help="Template parameter file")
    parser.add_argument("--read-path", required=True, help="Input data directory")
    parser.add_argument("--work-dir", required=True,
                        help="Working directory for this test")
    parser.add_argument("--stem", required=True,
                        help="write_fname stem used by the parameter template")
    parser.add_argument(
        "--expect-empty-tables",
        action="store_true",
        help=(
            "Assert the run produced no patches and that the patch and ladder "
            "CSVs still exist carrying only their header rows"
        ),
    )
    return parser.parse_args()


def read_band(path: pathlib.Path):
    """Return (values, nodata, n_cols, geotransform) for band 1 as a flat list.

    Reads through ReadRaster + struct rather than ReadAsArray so the check does
    not depend on the GDAL numpy bindings, which are ABI-coupled to the numpy
    build and are not a dependency any other test in this suite takes. Every
    band is pulled as float64 regardless of its on-disk type, so a UInt32 label
    and an Int16 band index come back through the same path.
    """
    ds = gdal.Open(str(path))
    if ds is None:
        raise RuntimeError(f"could not open raster: {path}")
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    n_cols = ds.RasterXSize
    n_rows = ds.RasterYSize
    gt = ds.GetGeoTransform()
    raw = band.ReadRaster(0, 0, n_cols, n_rows, buf_type=gdal.GDT_Float64)
    data = list(struct.unpack(f"<{n_rows * n_cols}d", raw))
    ds = None
    return data, nodata, n_cols, gt


def read_csv_rows(path: pathlib.Path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def rowcol_from_xy(gt, x: float, y: float):
    """Invert the geotransform for a north-up grid (all lsdtt3 output is)."""
    col = int(math.floor((x - gt[0]) / gt[1]))
    row = int(math.floor((y - gt[3]) / gt[5]))
    return row, col


def fail(message: str) -> None:
    print(f"Validation failed: {message}", file=sys.stderr)


def run_cli(args) -> pathlib.Path:
    exe = pathlib.Path(args.exe).resolve()
    template_path = pathlib.Path(args.param_template).resolve()
    read_path = pathlib.Path(args.read_path).resolve()
    work_dir = pathlib.Path(args.work_dir).resolve()
    out_dir = work_dir / "out"

    if not exe.exists():
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not template_path.exists():
        raise FileNotFoundError(f"Template not found: {template_path}")
    if not read_path.exists():
        raise FileNotFoundError(f"Read path not found: {read_path}")

    if work_dir.exists():
        shutil.rmtree(work_dir)
    out_dir.mkdir(parents=True)

    rendered = (template_path.read_text(encoding="utf-8")
                .replace("__READ_PATH__", str(read_path))
                .replace("__WRITE_PATH__", str(out_dir)))
    param_name = "run_params.txt"
    (work_dir / param_name).write_text(rendered, encoding="utf-8")

    result = lsdtt_test_env.run(
        [str(exe), str(work_dir), param_name],
        cwd=str(work_dir), capture_output=True, text=True, check=False)

    (work_dir / "stdout.log").write_text(result.stdout, encoding="utf-8")
    (work_dir / "stderr.log").write_text(result.stderr, encoding="utf-8")

    if result.returncode != 0:
        fail(f"CLI returned {result.returncode}")
        print(f"Stdout: {work_dir / 'stdout.log'}", file=sys.stderr)
        print(f"Stderr: {work_dir / 'stderr.log'}", file=sys.stderr)
        sys.exit(result.returncode)

    return out_dir


PATCH_HEADERS = [
    "label", "band", "class", "source_id", "n_pixels",
    "area_m2", "median_field_m", "aspect_degrees",
    "terrain_aspect_degrees", "aspect_r_bar", "window_centre_m"]

CLASS_HEADERS = [
    "label", "class", "source_id", "n_pixels",
    "median_hand_m", "mean_hand_m", "min_hand_m", "max_hand_m",
    "min_dist_to_channel_m", "adjacent_to_channel", "leak_split",
    "split_hand_m", "hand_band", "band_split", "hand_spread_m",
    "spread_guard_failed", "median_channel_distance_m",
    "ladder_window_centre_m"]

LADDER_HEADERS = [
    "source_id", "window_centre_m", "n_pixels", "k_selected",
    "status", "used_global_breaks", "breaks_m", "band_centres_m",
    "min_band_centre_gap_m", "required_gap_m",
    "k_bootstrap_mode", "k_stability", "k_ci_lo", "k_ci_hi", "n_bootstrap",
    "vertical_accuracy_m", "median_gradient", "goodness_of_variance_fit",
    "separation_ratio"]


def check_empty_tables(out_dir: pathlib.Path, stem: str) -> int:
    """A result with no patches must be an empty table, not a missing file.

    The writers used to suppress a CSV with no rows, which is the one output
    shape a schema assertion cannot distinguish from a crash: a downstream
    reader sees a missing file either way. So the contract is that the header
    row is always written.

    The ladder-window table is deliberately NOT required to be header-only. The
    ladder is upstream of the size filter that empties the patch table, so it
    has a real verdict to report and section 4.5 makes that verdict an output in
    its own right. What is asserted is that its schema is intact and that every
    row is well formed, not that it is empty.
    """
    errors = 0

    for name, headers in ((f"{stem}_surface_patches.csv", PATCH_HEADERS),
                          (f"{stem}_surface_class.csv", CLASS_HEADERS)):
        path = out_dir / name
        if not path.exists():
            fail(f"{name} was not written; an empty result must still emit "
                 "the header row")
            errors += 1
            continue
        fields, rows = read_csv_rows(path)
        if fields != headers:
            fail(f"{name} header is {fields}, expected {headers}")
            errors += 1
        if rows:
            fail(f"{name} has {len(rows)} data row(s); this run was expected "
                 "to produce none")
            errors += 1
        else:
            print(f"  OK  {name}: header only, {len(fields)} column(s)")

    ladder_path = out_dir / f"{stem}_ladder_windows.csv"
    if not ladder_path.exists():
        fail(f"{ladder_path.name} was not written; an empty result must still "
             "emit the header row")
        errors += 1
    else:
        fields, rows = read_csv_rows(ladder_path)
        if fields != LADDER_HEADERS:
            fail(f"{ladder_path.name} header is {fields}, expected "
                 f"{LADDER_HEADERS}")
            errors += 1
        valid_status = {"resolved", "unresolved", "forced", "degenerate"}
        bad = [r for r in rows
               if r["status"] not in valid_status or int(r["k_selected"]) < 1]
        if bad:
            fail(f"{ladder_path.name} has {len(bad)} malformed row(s); "
                 f"e.g. status={bad[0]['status']}, "
                 f"k_selected={bad[0]['k_selected']}")
            errors += 1
        else:
            print(f"  OK  {ladder_path.name}: header plus {len(rows)} "
                  "well-formed window(s)")

    # A labelled raster over an empty admitted set must carry no label at all.
    labels_path = out_dir / f"{stem}_banded_labels.tif"
    if labels_path.exists():
        label_data, label_nd, _, _ = read_band(labels_path)
        n_labelled = sum(
            1 for v in label_data
            if not (label_nd is not None and v == label_nd) and int(round(v)) > 0)
        if n_labelled:
            fail(f"{labels_path.name} carries {n_labelled} labelled pixel(s) "
                 "but the patch table is empty")
            errors += 1
        else:
            print(f"  OK  {labels_path.name}: no labelled pixel, consistent "
                  "with an empty patch table")

    return errors


def check_joins(out_dir: pathlib.Path, stem: str) -> int:
    errors = 0

    labels_path = out_dir / f"{stem}_banded_labels.tif"
    band_path = out_dir / f"{stem}_band_index.tif"
    class_raster_path = out_dir / f"{stem}_surface_class.tif"
    patches_path = out_dir / f"{stem}_surface_patches.csv"
    class_path = out_dir / f"{stem}_surface_class.csv"
    profile_path = out_dir / f"{stem}_reference_field_profile.csv"

    for path in (labels_path, band_path, class_raster_path,
                 patches_path, class_path, profile_path):
        if not path.exists():
            fail(f"missing expected output: {path.name}")
            errors += 1
    if errors:
        return errors

    label_data, label_nd, n_cols, gt = read_band(labels_path)
    band_data, band_nd, _, _ = read_band(band_path)
    class_data, _, _, _ = read_band(class_raster_path)

    # The schemas the joins below are written against. Pinned here as well as in
    # the smoke test so a column rename shows up as a schema failure rather than
    # as a KeyError halfway through a join.
    for path, headers in ((patches_path, PATCH_HEADERS),
                          (class_path, CLASS_HEADERS)):
        fields, _ = read_csv_rows(path)
        if fields != headers:
            fail(f"{path.name} header is {fields}, expected {headers}")
            errors += 1
    if errors:
        return errors

    # ---- 1 / 2 / 3: the rasters against each other and against the patch table
    raster_bands: dict[int, set[int]] = {}
    raster_counts: dict[int, int] = {}
    for i, lv in enumerate(label_data):
        if label_nd is not None and lv == label_nd:
            continue
        label = int(round(lv))
        if label <= 0:
            continue
        bv = band_data[i]
        if band_nd is not None and bv == band_nd:
            fail(f"pixel {i} carries label {label} but nodata in the band index")
            errors += 1
            continue
        band = int(round(bv))
        if band < 0:
            fail(f"pixel {i} carries label {label} but band {band}: a labelled "
                 "pixel must sit in a band")
            errors += 1
            continue
        raster_bands.setdefault(label, set()).add(band)
        raster_counts[label] = raster_counts.get(label, 0) + 1

    multi = {lab: bands for lab, bands in raster_bands.items() if len(bands) > 1}
    if multi:
        sample = list(multi.items())[:5]
        fail(f"the label -> band map is multi-valued for {len(multi)} label(s); "
             f"e.g. {sample}")
        errors += 1
    else:
        print(f"  OK  label -> band map is single-valued over "
              f"{len(raster_bands)} label(s)")

    patch_fields, patch_rows = read_csv_rows(patches_path)
    patch_by_label = {}
    for row in patch_rows:
        label = int(row["label"])
        if label in patch_by_label:
            fail(f"_surface_patches.csv has a duplicate row for label {label}")
            errors += 1
        patch_by_label[label] = row

    raster_labels = set(raster_bands)
    patch_labels = set(patch_by_label)
    if raster_labels != patch_labels:
        only_raster = sorted(raster_labels - patch_labels)[:10]
        only_table = sorted(patch_labels - raster_labels)[:10]
        fail(f"label sets differ: {len(raster_labels - patch_labels)} only in "
             f"the raster (e.g. {only_raster}), "
             f"{len(patch_labels - raster_labels)} only in the patch table "
             f"(e.g. {only_table})")
        errors += 1
    else:
        print(f"  OK  _banded_labels.tif and _surface_patches.csv carry the "
              f"same {len(patch_labels)} label(s)")

    if patch_labels and sorted(patch_labels) != list(range(1, len(patch_labels) + 1)):
        fail(f"labels are not a contiguous run 1..{len(patch_labels)}: "
             f"min {min(patch_labels)}, max {max(patch_labels)}")
        errors += 1
    else:
        print(f"  OK  labels are a contiguous run 1..{len(patch_labels)}")

    band_mismatch = 0
    count_mismatch = 0
    for label, row in patch_by_label.items():
        if label not in raster_bands:
            continue
        table_band = int(row["band"])
        raster_band = next(iter(raster_bands[label]))
        if table_band != raster_band:
            if band_mismatch < 5:
                fail(f"label {label}: patch table says band {table_band}, "
                     f"_band_index.tif says {raster_band}")
            band_mismatch += 1
        table_n = int(row["n_pixels"])
        if table_n != raster_counts[label]:
            if count_mismatch < 5:
                fail(f"label {label}: patch table says {table_n} pixel(s), "
                     f"the raster has {raster_counts[label]}")
            count_mismatch += 1
    if band_mismatch:
        fail(f"{band_mismatch} label(s) disagree on band between the raster "
             "and the patch table")
        errors += 1
    else:
        print("  OK  patch band agrees with _band_index.tif for every label")
    if count_mismatch:
        fail(f"{count_mismatch} label(s) disagree on n_pixels")
        errors += 1
    else:
        print("  OK  patch n_pixels agrees with the labelled raster")

    # ---- 4: the class table against the patch table
    class_fields, class_rows = read_csv_rows(class_path)
    class_by_label = {}
    for row in class_rows:
        label = int(row["label"])
        if label in class_by_label:
            fail(f"_surface_class.csv has a duplicate row for label {label}")
            errors += 1
        class_by_label[label] = row

    if set(class_by_label) != patch_labels:
        only_class = sorted(set(class_by_label) - patch_labels)[:10]
        only_patch = sorted(patch_labels - set(class_by_label))[:10]
        fail(f"_surface_class.csv is not 1:1 with _surface_patches.csv: "
             f"{len(set(class_by_label) - patch_labels)} only in the class "
             f"table (e.g. {only_class}), "
             f"{len(patch_labels - set(class_by_label))} only in the patch "
             f"table (e.g. {only_patch})")
        errors += 1
    else:
        print(f"  OK  _surface_class.csv joins 1:1 to _surface_patches.csv on "
              f"{len(class_by_label)} label(s)")

    col_mismatch = 0
    for label, prow in patch_by_label.items():
        crow = class_by_label.get(label)
        if crow is None:
            continue
        for patch_col, class_col in (("class", "class"),
                                     ("source_id", "source_id"),
                                     ("window_centre_m", "ladder_window_centre_m")):
            if prow[patch_col] != crow[class_col]:
                if col_mismatch < 5:
                    fail(f"label {label}: patch {patch_col}={prow[patch_col]} "
                         f"but class {class_col}={crow[class_col]}")
                col_mismatch += 1
        if int(prow["n_pixels"]) != int(crow["n_pixels"]):
            if col_mismatch < 5:
                fail(f"label {label}: patch n_pixels={prow['n_pixels']} but "
                     f"class n_pixels={crow['n_pixels']}")
            col_mismatch += 1
        if int(prow["band"]) != int(crow["hand_band"]):
            if col_mismatch < 5:
                fail(f"label {label}: patch band={prow['band']} but class "
                     f"hand_band={crow['hand_band']}")
            col_mismatch += 1
    if col_mismatch:
        fail(f"{col_mismatch} shared column value(s) disagree between the "
             "patch and class tables")
        errors += 1
    else:
        print("  OK  shared columns agree between the patch and class tables")

    # ---- 5: the HAND profile against the rasters and the patch table
    profile_fields, profile_rows = read_csv_rows(profile_path)
    if not profile_rows:
        fail("_reference_field_profile.csv has no data rows, so it proves no join")
        errors += 1

    orphan = 0
    band_disagree = 0
    label_disagree = 0
    class_disagree = 0
    for row in profile_rows:
        label = int(row["label"])
        prow = patch_by_label.get(label)
        if prow is None:
            if orphan < 5:
                fail(f"_reference_field_profile.csv row carries label {label}, which has "
                     "no row in _surface_patches.csv")
            orphan += 1
            continue
        if int(row["band"]) != int(prow["band"]):
            if band_disagree < 5:
                fail(f"_reference_field_profile.csv label {label}: band="
                     f"{row['band']} but patch band={prow['band']}")
            band_disagree += 1

        r, c = rowcol_from_xy(gt, float(row["x"]), float(row["y"]))
        idx = r * n_cols + c
        if idx < 0 or idx >= len(label_data):
            fail(f"_reference_field_profile.csv row maps outside the grid: x={row['x']}, "
                 f"y={row['y']}")
            errors += 1
            continue
        if int(round(label_data[idx])) != label:
            if label_disagree < 5:
                fail(f"_reference_field_profile.csv row at ({row['x']}, {row['y']}) says "
                     f"label {label} but _banded_labels.tif says "
                     f"{int(round(label_data[idx]))}")
            label_disagree += 1
        if int(round(class_data[idx])) != int(row["class"]):
            if class_disagree < 5:
                fail(f"_reference_field_profile.csv row at ({row['x']}, {row['y']}) says "
                     f"class {row['class']} but _surface_class.tif says "
                     f"{int(round(class_data[idx]))}")
            class_disagree += 1

    for count, message in (
            (orphan, "HAND-profile row(s) carry a label with no patch row"),
            (band_disagree, "HAND-profile row(s) disagree with the patch band"),
            (label_disagree,
             "HAND-profile row(s) disagree with _banded_labels.tif"),
            (class_disagree,
             "HAND-profile row(s) disagree with _surface_class.tif")):
        if count:
            fail(f"{count} {message}")
            errors += 1
    if not (orphan or band_disagree or label_disagree or class_disagree):
        print(f"  OK  all {len(profile_rows)} _reference_field_profile.csv row(s) join to "
              "a patch and agree with both rasters")

    return errors


def main() -> int:
    args = parse_args()
    out_dir = run_cli(args)

    print(f"Checking fluvial-surface output joins in {out_dir}")
    if args.expect_empty_tables:
        errors = check_empty_tables(out_dir, args.stem)
    else:
        errors = check_joins(out_dir, args.stem)

    if errors:
        print(f"FAIL: {errors} join check(s) failed", file=sys.stderr)
        return 1
    print("PASS: all join checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
