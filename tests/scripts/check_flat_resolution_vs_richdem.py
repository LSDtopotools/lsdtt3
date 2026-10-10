#!/usr/bin/env python3
"""Oracle O1: lsdtt3's barnes2014b_float_step flat resolution vs RichDEM.

Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3

Copyright (C) 2026 Simon M. Mudd
SPDX-License-Identifier: AGPL-3.0-or-later

WHAT IT PROVES

  The FIRST ROUND of `fill_flat_resolution_choice: barnes2014b_float_step`
  is meant to be bit-for-bit RichDEM's `ResolveFlats` on a float32 DEM:
  Barnes, Lehman and Mulla (2014b), "An efficient assignment of drainage
  direction over flat surfaces in raster digital elevation models",
  Computers & Geosciences 62, https://doi.org/10.1016/j.cageo.2013.01.009.
  Each flat cell is raised by `mask` successive float32 steps (nextafter
  towards +inf), where the mask combines the distance to the flat's low edges
  and away from its high edges.

  That single pass can leave cells with no lower neighbour (a rim cell a few
  ULPs above a flat that the raise overtops, or a flat's low edge overtopped
  by a lower flat's raise), and a D8 router refuses them. lsdtt3 therefore
  refills and resolves again (rounds 2-4) and, if cells are still stranded,
  runs a final strict-rise priority-flood pass. Those later changes are the
  ONLY permitted differences from RichDEM, and lsdtt3 reports how many cells
  they changed ("Changed after round 1: <n> cell(s)").

  For every fixture this script:

    1. fills the DEM with lsdtt-raster-preprocessing (vincent1993,
       fill_min_slope 0), giving F, a surface with TRUE flats;
    2. runs RichDEM ResolveFlats on a float32 rdarray of F, giving R, and
       captures the pair count RichDEM prints;
    3. runs lsdtt-raster-preprocessing on F with barnes2014b_float_step
       (refilling F is a no-op), giving L, and parses its pair count;
    4. asserts L and R have the same bits (uint32 view) on every valid cell,
       except on near-zero-guarded flats (below) and on the cells lsdtt3
       changed after round 1: the number of differing unguarded cells must
       equal lsdtt3's "Changed after round 1" count (n_unexplained == 0), and
       every one of them must be strictly higher in L than in R (ordered
       float32 keys; no lsdtt3 cell is ever below RichDEM's);
    5. asserts lsdtt3's ROUND-1 pair count ("Round 1 (RichDEM-comparable):
       <p> pair(s), <c> cell(s) ...") equals the pair / cell recount from R,
       and RichDEM's "Cells inappropriately raised above surrounding terrain"
       too wherever that count is order independent (otherwise it must lie in
       the range its in-place loop allows; see richdem_printed_pair_bound);
       and that lsdtt3's final "Raised above surroundings" line equals the
       recount from L;
    5b. counts the cells without a lower neighbour (valid, off the border, no
       nodata neighbour) in R and in L: R's count must equal lsdtt3's round-1
       count, and L's must be 0 -- lsdtt3 always returns a routable surface;
    6. asserts both outputs have the same nodata cells, and that F, L and the
       fixture share raster size and geotransform.

  The count equality in 4 cannot see a round-1 bug confined to a cell the
  loop later changed; the C++ driver (unit_flat_resolution, F21-F23 and F29:
  round 1 equals max_rounds = 1, unchanged cells bitwise) covers that.

  Fixtures are generated here (no network): hand-built valleys, terraces, a
  mesa, an interior nodata hole, a big flat whose rim is 1 ULP above it, a
  valley with a 2-ULP bump (pins RichDEM's `>=` pair rule: it prints 8, where
  a `>` rule would give 7; RichDEM leaves the bump itself without a lower
  neighbour, so lsdtt3's second round raises it and its final count is 7),
  valley floors at +0, -0, -50 and 4000 m, and two fixed-seed 256 x 256
  integer-metre fractal surfaces around 300 m and 3000 m.

  Fixtures on which RichDEM itself strands cells (asserted: R has cells
  without a lower neighbour, L has none): a knob 1 ULP above the valley floor,
  two terraces 2 and 3 ULPs apart, a cascade of five terraces 2 ULPs apart,
  and the committed Copernicus GLO-30 lakeland crop
  (tests/data/flat_resolution_lakeland_crop_dem.tif).

WHY ELEVATIONS, NOT D8 DIRECTIONS

  The comparison is on the resolved ELEVATIONS. RichDEM's D8 is not a valid
  reference for flow direction: it picks the lowest neighbour rather than the
  steepest descent (see tests/scripts/run_flow_benchmark.py, the accum_d8
  notes), so comparing routed directions would score that rule difference as
  lsdtt3 error. Equal elevations, bit for bit, imply that any steepest-descent
  router sees the same surface.

THE NEAR-ZERO GUARD EXCEPTION

  On a flat with |z| < 1 m lsdtt3 does NOT take native float32 steps: near 0
  they are subnormal (1.4e-45 m at z = 0), which float32 slope arithmetic
  rounds to zero. It uses a fixed 2^-23 m per mask unit instead,
  z' = float(double(z) + mask * 2^-23). On those cells, and only those, L
  differs from R by design. The script asserts that the differing cells are
  EXACTLY the guarded flat cells (raised, non-border cells of flats with
  |z| < 1 m), that each guarded L value is the guard formula applied to
  RichDEM's own mask (its raise in ULPs), and hence that L orders the cells
  exactly as R does.

RUN

  pixi run -e benchmark python -I tests/scripts/check_flat_resolution_vs_richdem.py \\
      --bin-dir build --work-dir <tmp>

  Exit 0 pass, 1 fail, 77 skip (RichDEM, GDAL or numpy not importable).
  Registered in CTest as benchmark_flat_resolution_richdem (needs
  -DLSDTT_ENABLE_BENCHMARKS=ON).
"""

from __future__ import annotations

import argparse
import os
import pathlib
import re
import shutil
import sys

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))  # python -I does not add it
import lsdtt_test_env  # noqa: E402

SKIP_EXIT = 77

try:
    import numpy as np
except ImportError as exc:
    print(f"SKIP  numpy is not importable ({exc}). Run in the pixi benchmark "
          f"environment: pixi run -e benchmark python -I {__file__} ...")
    sys.exit(SKIP_EXIT)
try:
    from osgeo import gdal, osr
except ImportError as exc:
    print(f"SKIP  GDAL (osgeo) is not importable ({exc}). Run in the pixi "
          f"benchmark environment: pixi run -e benchmark python -I {__file__} ...")
    sys.exit(SKIP_EXIT)
try:
    import richdem as rd
except ImportError as exc:
    print(f"SKIP  richdem is not importable ({exc}). It comes from conda-forge "
          f"in the pixi benchmark environment: pixi run -e benchmark python -I "
          f"{__file__} ...")
    sys.exit(SKIP_EXIT)

gdal.UseExceptions()

ND = np.float32(-9999.0)
EPSG = 32630                 # WGS 84 / UTM zone 30N; any projected CRS will do
CELL = 10.0                  # metres
GUARD_STEP = 2.0 ** -23      # lsdtt3's near-zero guard step, metres per mask unit
GUARD_LIMIT = 1.0            # flats with |z| < 1 m are guarded
EXE_NAME = "lsdtt-raster-preprocessing"
# Committed Copernicus GLO-30 crop (provenance in tests/README.md).
REAL_TERRAIN_CROP = SCRIPT_DIR.parent / "data" / "flat_resolution_lakeland_crop_dem.tif"
N8 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


# ---------------------------------------------------------------------------
# float32 helpers
# ---------------------------------------------------------------------------
def ordered_key(a: np.ndarray) -> np.ndarray:
    """float32 -> int64 key, monotone in value; -0.0 and +0.0 both map to 0,
    so key(b) - key(a) is the number of nextafter steps from a up to b."""
    i = np.asarray(a, np.float32).view(np.int32).astype(np.int64)
    return np.where(i < 0, -(i & 0x7FFFFFFF), i)


def shifted(a: np.ndarray, dr: int, dc: int, fill) -> np.ndarray:
    """out[r, c] = a[r + dr, c + dc], `fill` where that is off the raster."""
    out = np.full_like(a, fill)
    R, C = a.shape
    rs, re_ = max(0, -dr), min(R, R - dr)
    cs, ce = max(0, -dc), min(C, C - dc)
    out[rs:re_, cs:ce] = a[rs + dr:re_ + dr, cs + dc:ce + dc]
    return out


def count_raised_above(before: np.ndarray, after: np.ndarray,
                       valid: np.ndarray) -> tuple[int, int]:
    """RichDEM's "raised above surrounding terrain" count, which counts PAIRS:
    ordered 8-neighbour pairs (c, n), both valid, with before[n] > before[c]
    and after[c] >= after[n] (the `>=` rule, settled in Step 0). Also returns
    the number of distinct cells c involved."""
    pairs = 0
    cells = np.zeros(before.shape, bool)
    for dr, dc in N8:
        nb = shifted(before, dr, dc, np.float32(np.nan))
        na = shifted(after, dr, dc, np.float32(np.nan))
        nv = shifted(valid, dr, dc, False)
        hit = valid & nv & (nb > before) & (after >= na)
        pairs += int(hit.sum())
        cells |= hit
    return pairs, int(cells.sum())


def richdem_printed_pair_bound(before: np.ndarray, after: np.ndarray,
                               valid: np.ndarray) -> int:
    """Upper bound on the pair count RichDEM PRINTS for a given output.

    RichDEM (flats/Barnes2014.hpp, ResolveFlatsEpsilon_Barnes2014) raises the
    cells IN PLACE, in a parallel loop, and compares each raised cell with its
    neighbours' CURRENT values: a neighbour on another flat that is processed
    earlier has already been raised, one processed later has not. So its count
    for a pair (c, n) uses v = after[n] or before[n], depending on the order:
    hit if before[c] < v <= after[c]. With v = after[n] for every pair this is
    count_raised_above() (the lower bound, which lsdtt3 reports); this is the
    count when either value may be used. The two agree unless a raised cell
    borders a raised cell of another flat (the lakeland crop: 364 vs 370)."""
    pairs = 0
    for dr, dc in N8:
        nb = shifted(before, dr, dc, np.float32(np.nan))
        na = shifted(after, dr, dc, np.float32(np.nan))
        nv = shifted(valid, dr, dc, False)
        hit_b = (before < nb) & (nb <= after)
        hit_a = (before < na) & (na <= after) & (nb > before)
        pairs += int((valid & nv & (hit_b | hit_a)).sum())
    return pairs


def flat_cells(z: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """The cells Barnes (2014b) raises: NO_FLOW cells (no lower neighbour,
    not draining) plus low edges (draining cells next to a NO_FLOW cell of
    equal elevation). Border and nodata-adjacent cells drain."""
    R, C = z.shape
    drains = np.zeros((R, C), bool)
    drains[0, :] = drains[-1, :] = drains[:, 0] = drains[:, -1] = True
    for dr, dc in N8:
        nv = shifted(valid, dr, dc, True)
        nz = shifted(z, dr, dc, np.float32(np.inf))
        drains |= ~nv | (nv & (nz < z))
    noflow = valid & ~drains
    low_edge = np.zeros((R, C), bool)
    for dr, dc in N8:
        nnf = shifted(noflow, dr, dc, False)
        nz = shifted(z, dr, dc, np.float32(np.nan))
        low_edge |= valid & drains & nnf & (nz == z)
    return noflow | low_edge


def without_descent(z: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Cells a D8 router cannot give a receiver: valid, off the border, no
    nodata neighbour, and no strictly lower neighbour (the predicate of
    lsdtt3's "cells without descent" count)."""
    R, C = z.shape
    drains = np.zeros((R, C), bool)
    drains[0, :] = drains[-1, :] = drains[:, 0] = drains[:, -1] = True
    for dr, dc in N8:
        nv = shifted(valid, dr, dc, True)
        nz = shifted(z, dr, dc, np.float32(np.inf))
        drains |= ~nv | (nz < z)
    return valid & ~drains


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------
def valley(floor: float = 100.0, R: int = 11, C: int = 24) -> np.ndarray:
    """11 x 24 valley: flat floor, walls 10 m higher on all four sides, outlet
    1 m below the floor at row R//2, col 0."""
    z = np.full((R, C), np.float32(floor + 10.0), np.float32)
    z[1:R - 1, 1:C - 1] = np.float32(floor)
    z[R // 2, 0] = np.float32(floor - 1.0)
    return z


def fractal(seed: int, base: float, n: int = 256, relief: float = 120.0) -> np.ndarray:
    """Fixed-seed spectral-synthesis surface (power spectrum ~ k^-3.2),
    rounded to integer metres. Rounding makes plateaus and the fill turns the
    closed depressions into flats: plenty of flats of every shape."""
    rng = np.random.default_rng(seed)
    ky = np.fft.fftfreq(n)[:, None]
    kx = np.fft.rfftfreq(n)[None, :]
    k = np.hypot(kx, ky)
    k[0, 0] = 1.0
    amp = k ** (-3.2 / 2.0)
    amp[0, 0] = 0.0
    spec = amp * (rng.standard_normal(k.shape) + 1j * rng.standard_normal(k.shape))
    f = np.fft.irfft2(spec, s=(n, n))
    f = (f - f.min()) / (f.max() - f.min())
    return np.round(base + relief * (f - 0.5)).astype(np.float32)


def build_fixtures() -> dict[str, np.ndarray]:
    F: dict[str, np.ndarray] = {}
    F["valley"] = valley()
    z = valley(); z[5, 0] = np.float32(110.0); z[4:7, 0] = np.float32(100.0)
    F["open_flat_edge"] = z
    z = valley(); z[5, 12] = ND
    F["interior_nodata"] = z
    z = np.full((11, 24), np.float32(120.0), np.float32)
    z[1:10, 1:12] = np.float32(105.0)          # upper terrace
    z[1:10, 12:23] = np.float32(100.0)         # lower terrace
    z[5, 23] = np.float32(99.0)                # outlet on the east edge
    F["two_terraces"] = z
    z = np.full((11, 24), np.float32(90.0), np.float32)
    z[2:9, 3:21] = np.float32(100.0)           # mesa top: no high edges
    F["mesa_no_high_edges"] = z
    z = valley(R=41, C=81)
    rim = np.nextafter(np.float32(100.0), np.float32(200.0))
    z[0, :] = z[-1, :] = z[:, -1] = z[:, 0] = rim
    z[20, 0] = np.float32(99.0)
    F["big_flat_rim_1ulp"] = z
    z = valley()
    bump = np.float32(100.0)
    for _ in range(2):
        bump = np.nextafter(bump, np.float32(200.0))
    z[3, 2] = bump                              # 2 ULPs: one pair ends level
    F["bump_2ulp"] = z
    for name, fl in (("floor_+0", 0.0), ("floor_-0", -0.0),
                     ("floor_-50", -50.0), ("floor_+4000", 4000.0)):
        F[f"valley_{name}"] = valley(floor=fl)
    F["fractal_300m"] = fractal(seed=20141, base=300.0)
    F["fractal_3000m"] = fractal(seed=20142, base=3000.0)
    return F


def step_up(z: float, k: int) -> np.float32:
    """k float32 steps above z (nextafter towards +inf)."""
    v = np.float32(z)
    for _ in range(k):
        v = np.nextafter(v, np.float32(np.inf))
    return v


def build_stranding_fixtures() -> dict[str, np.ndarray | pathlib.Path]:
    """Fixtures on which RichDEM ResolveFlats leaves cells without a lower
    neighbour. The geometries are those of the C++ driver's F21-F23."""
    F: dict[str, np.ndarray | pathlib.Path] = {}
    z = valley()
    z[5, 12] = step_up(100.0, 1)               # the knob, inside the floor
    F["knob_1ulp"] = z
    for gap in (2, 3):                         # upper terrace gap ULPs up
        z = np.full((11, 24), np.float32(120.0), np.float32)
        z[1:10, 1:12] = step_up(100.0, gap)
        z[1:10, 12:23] = np.float32(100.0)
        z[5, 23] = np.float32(99.0)
        F[f"terraces_{gap}ulp"] = z
    k, w = 5, 6                                # five terraces, 2 ULPs apart
    z = np.full((11, k * w + 2), np.float32(120.0), np.float32)
    for i in range(k):
        z[1:10, 1 + i * w:(i + 1) * w + 1] = step_up(100.0, 2 * (k - 1 - i))
    z[5, k * w + 1] = np.float32(99.0)
    F["cascade_5"] = z
    F["lakeland_crop"] = REAL_TERRAIN_CROP
    return F


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def write_tif(path: pathlib.Path, z: np.ndarray) -> None:
    R, C = z.shape
    ds = gdal.GetDriverByName("GTiff").Create(str(path), C, R, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((500000.0, CELL, 0.0, 5000000.0 + R * CELL, 0.0, -CELL))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(float(ND))
    band.WriteArray(z)
    band.FlushCache()
    ds = None


def read_tif(path: pathlib.Path):
    ds = gdal.Open(str(path))
    band = ds.GetRasterBand(1)
    if band.DataType != gdal.GDT_Float32:
        raise RuntimeError(f"{path.name}: band type "
                           f"{gdal.GetDataTypeName(band.DataType)}, not Float32")
    z = band.ReadAsArray().astype(np.float32, copy=False)
    nd = band.GetNoDataValue()
    valid = ~np.isnan(z)
    if nd is not None:
        valid &= z != np.float32(nd)
    info = (ds.RasterXSize, ds.RasterYSize, tuple(ds.GetGeoTransform()))
    ds = None
    return z, valid, info


def run_cli(exe: pathlib.Path, run_dir: pathlib.Path, src: pathlib.Path,
            stem: str, extra: dict[str, str]) -> str:
    """Run lsdtt-raster-preprocessing; return its stdout."""
    params = {
        "read_path": str(src.parent),
        "read_fname": src.name,
        "write_path": str(run_dir),
        "write_fname": stem,
        "enforce_projection": "false",
        "write_output_format": "gtiff",
        "fill_algorithm_choice": "vincent1993",
        "fill_min_slope": "0",
        "write_fill_raster": "true",
        # The default (0) masks every cell at or below 0 m to nodata, which
        # would delete the floor of the +0, -0 and -50 m valleys. Negative
        # disables it.
        "min_elevation": "-1",
    }
    params.update(extra)
    pfile = run_dir / f"{stem}.param"
    pfile.write_text("".join(f"{k}: {v}\n" for k, v in params.items()))
    # The binary is built against the system GDAL; strip the pixi env's
    # library and data paths (see lsdtt_test_env.system_env).
    proc = lsdtt_test_env.run([str(exe), str(run_dir), pfile.name],
                              capture_output=True, text=True,
                              env=lsdtt_test_env.system_env())
    if proc.returncode != 0:
        hint = lsdtt_test_env.explain_failure(proc.stderr)
        raise RuntimeError(f"{exe.name} failed (exit {proc.returncode}) on "
                           f"{pfile}\n--- stdout ---\n{proc.stdout[-2000:]}\n"
                           f"--- stderr ---\n{proc.stderr[-2000:]}\n{hint}")
    return proc.stdout


def run_richdem(z: np.ndarray, log_path: pathlib.Path) -> tuple[np.ndarray, int]:
    """rd.ResolveFlats on a float32 rdarray. RichDEM prints its counts from
    C++, so capture at the file-descriptor level, not via sys.stdout."""
    a = rd.rdarray(z.astype(np.float32, copy=True), no_data=float(ND))
    sys.stdout.flush(); sys.stderr.flush()
    saved = os.dup(1), os.dup(2)
    fd = os.open(str(log_path), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    os.dup2(fd, 1); os.dup2(fd, 2)
    try:
        out = rd.ResolveFlats(a)
    finally:
        sys.stdout.flush(); sys.stderr.flush()
        os.dup2(saved[0], 1); os.dup2(saved[1], 2)
        os.close(saved[0]); os.close(saved[1]); os.close(fd)
    log = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "",
                 log_path.read_text(errors="replace"))
    m = re.search(r"Cells inappropriately raised above surrounding terrain\s*=\s*(\d+)", log)
    if not m:   # printed even when there are no flats ("There were no flats!")
        raise RuntimeError(f"RichDEM pair count not found in its output "
                           f"({log_path}):\n{log[-2000:]}")
    out = np.asarray(out)
    if out.dtype != np.float32:
        raise RuntimeError(f"RichDEM returned {out.dtype}, not float32")
    return out, int(m.group(1))


# ---------------------------------------------------------------------------
# one fixture
# ---------------------------------------------------------------------------
PAIR_RE = re.compile(r"Raised above surroundings:\s*(\d+) pair\(s\),\s*(\d+) cell\(s\)")
ROUND1_RE = re.compile(r"Round 1 \(RichDEM-comparable\):\s*(\d+) pair\(s\),\s*(\d+) cell\(s\) "
                       r"raised above surroundings,\s*(\d+) cell\(s\) without descent")
CHANGED_RE = re.compile(r"Changed after round 1:\s*(\d+) cell\(s\)")
ROUNDS_RE = re.compile(r"Rounds:\s*(\d+) \(cells without descent after each round: ([^)]*)\)")
FINAL_PASS_RE = re.compile(r"Final priority-flood pass:\s*(?:(\d+) cell\(s\)|not needed)")


def check_fixture(name: str, z0: np.ndarray | pathlib.Path, exe: pathlib.Path,
                  root: pathlib.Path, richdem_strands: bool = False) -> tuple[bool, str]:
    """One fixture. `z0` is an array (written here) or a GeoTIFF path (copied).
    With `richdem_strands` the fixture must be one on which RichDEM leaves
    cells without a lower neighbour."""
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    errors: list[str] = []

    src = d / f"{name}.tif"
    if isinstance(z0, pathlib.Path):
        if not z0.is_file():
            return False, f"{name}: fixture {z0} not found"
        shutil.copyfile(z0, src)
    else:
        write_tif(src, z0)
    _, v0, info0 = read_tif(src)

    # 1. F: vincent1993 at fill_min_slope 0 keeps true flats.
    run_cli(exe, d, src, f"{name}_F", {})
    F, vF, infoF = read_tif(d / f"{name}_F_fill.tif")

    # 2. R: RichDEM on F.
    R, rd_pairs = run_richdem(np.where(vF, F, ND).astype(np.float32),
                              d / "richdem_stdout.txt")
    vR = R != ND

    # 3. L: lsdtt3 float step on F.
    out = run_cli(exe, d, d / f"{name}_F_fill.tif", f"{name}_L",
                  {"fill_flat_resolution_choice": "barnes2014b_float_step"})
    (d / "lsdtt_stdout.txt").write_text(out)
    L, vL, infoL = read_tif(d / f"{name}_L_fill.tif")
    m = PAIR_RE.search(out)
    m1 = ROUND1_RE.search(out)
    mc = CHANGED_RE.search(out)
    mr = ROUNDS_RE.search(out)
    mf = FINAL_PASS_RE.search(out)
    if not (m and m1 and mc and mr and mf):
        return False, (f"{name}: a 'Raised above surroundings', 'Round 1', 'Changed after "
                       f"round 1', 'Rounds' or 'Final priority-flood pass' line is missing "
                       f"from lsdtt3 stdout")
    ls_pairs, ls_cells = int(m.group(1)), int(m.group(2))
    r1_pairs, r1_cells, r1_stranded = int(m1.group(1)), int(m1.group(2)), int(m1.group(3))
    n_changed_r1 = int(mc.group(1))
    n_rounds = int(mr.group(1))
    n_final_pass = int(mf.group(1)) if mf.group(1) else 0

    # 6. Geometry and validity.
    if not (info0 == infoF == infoL):
        errors.append(f"size/geotransform differ: fixture {info0}, F {infoF}, L {infoL}")
    if R.shape != F.shape:
        errors.append(f"RichDEM shape {R.shape} != F shape {F.shape}")
        return False, f"{name}: " + "; ".join(errors)
    if not np.array_equal(v0, vF):
        errors.append(f"the fill changed validity: fixture {int(v0.sum())}, "
                      f"F {int(vF.sum())} valid cells")
    if not np.array_equal(vF, vL) or not np.array_equal(vF, vR):
        errors.append(f"validity differs: F {int(vF.sum())}, L {int(vL.sum())}, "
                      f"R {int(vR.sum())} valid cells")
    valid = vF & vL & vR

    # 4. Bits, outside the guard: the differing cells are exactly the cells
    # lsdtt3 changed after round 1 (by count), each strictly above RichDEM.
    border = np.zeros(F.shape, bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    flats = flat_cells(F, vF) & ~border
    guarded = flats & (np.abs(F) < GUARD_LIMIT)
    differ = valid & (L.view(np.uint32) != R.view(np.uint32))
    unguarded_diff = differ & ~guarded
    n_unguarded_diff = int(unguarded_diff.sum())
    n_unexplained = n_unguarded_diff - n_changed_r1
    below = unguarded_diff & (ordered_key(L) < ordered_key(R))
    if guarded.any() and n_changed_r1 > 0:
        errors.append("fixture mixes guarded and loop-changed cells; extend the check")
    if n_unexplained != 0:
        first = ""
        if unguarded_diff.any():
            rr, cc = np.argwhere(unguarded_diff)[0]
            first = (f"; first differing ({rr},{cc}): F {F[rr, cc]!r} L {L[rr, cc]!r} "
                     f"R {R[rr, cc]!r}")
        errors.append(f"{n_unguarded_diff} cell(s) differ outside the guard but lsdtt3 "
                      f"changed {n_changed_r1} after round 1 (unexplained "
                      f"{n_unexplained}){first}")
    if below.any():
        rr, cc = np.argwhere(below)[0]
        errors.append(f"{int(below.sum())} cell(s) are BELOW RichDEM; first ({rr},{cc}): "
                      f"L {L[rr, cc]!r} R {R[rr, cc]!r}")

    # Guarded cells: must ALL differ (both sides raise them), and L must be the
    # guard formula applied to RichDEM's mask, so the ordering is RichDEM's.
    mask_R = ordered_key(R) - ordered_key(F)
    if guarded.any():
        missing = guarded & valid & ~differ
        if missing.any():
            errors.append(f"{int(missing.sum())} guarded cell(s) do not differ "
                          f"from RichDEM")
        g = guarded & valid
        expect = (F[g].astype(np.float64) + mask_R[g] * GUARD_STEP).astype(np.float32)
        bad = expect.view(np.uint32) != L[g].view(np.uint32)
        if bad.any():
            errors.append(f"{int(bad.sum())} guarded cell(s) are not "
                          f"float(z + mask_R * 2^-23)")
        if (mask_R[g] <= 0).any():
            errors.append("a guarded cell is not raised by RichDEM")
        for dr, dc in N8:     # same ordering on every guarded neighbour pair
            ng = shifted(g, dr, dc, False)
            both = g & ng
            sR = np.sign(mask_R - shifted(mask_R, dr, dc, 0))[both]
            Ld = L.astype(np.float64)
            sL = np.sign(Ld - shifted(Ld, dr, dc, 0.0))[both]
            if (sR != sL).any():
                errors.append("guarded cells ordered differently from RichDEM")
                break

    # 5. Pair counts. lsdtt3's ROUND-1 count (round 1 is the RichDEM-equivalent
    # pass) must equal the recount from R, and its final line the recount from
    # L. RichDEM's PRINTED count is order dependent where two raised flats
    # touch (richdem_printed_pair_bound), so it must lie between the recount
    # from R and that bound; where they coincide (every synthetic fixture
    # here) that is equality. Across the two tools they must
    # agree except where the guard could change it: the guard raises by up to
    # mask * 2^-23 m, much more than RichDEM's subnormal steps, so a guarded
    # flat could overtop a rim RichDEM does not. The guarded fixtures here have
    # rims 10 m above the floor and max masks below 100 (raise < 1.2e-5 m), so
    # no pair can end level on either side and both counts are 0: they ARE
    # compared.
    py_R, py_R_cells = count_raised_above(F, R, valid)
    py_R_hi = richdem_printed_pair_bound(F, R, valid)
    py_L, py_L_cells = count_raised_above(F, L, valid)
    if not py_R <= rd_pairs <= py_R_hi:
        errors.append(f"RichDEM printed {rd_pairs} pairs, outside its recount range "
                      f"[{py_R}, {py_R_hi}]")
    if (py_L, py_L_cells) != (ls_pairs, ls_cells):
        errors.append(f"lsdtt3 printed {ls_pairs} pairs/{ls_cells} cells, "
                      f"recount gives {py_L}/{py_L_cells}")
    if (r1_pairs, r1_cells) != (py_R, py_R_cells):
        errors.append(f"lsdtt3 round 1 printed {r1_pairs} pairs/{r1_cells} cells, the "
                      f"recount from RichDEM's output gives {py_R} pairs/{py_R_cells} cells")

    # 5b. Cells without a lower neighbour: RichDEM's must be lsdtt3's round-1
    # count, and lsdtt3's final surface must have none.
    s_R = int(without_descent(R, vR).sum())
    s_L = int(without_descent(L, vL).sum())
    if s_R != r1_stranded:
        errors.append(f"RichDEM leaves {s_R} cell(s) without descent, lsdtt3 round 1 "
                      f"reports {r1_stranded}")
    if s_L != 0:
        errors.append(f"lsdtt3 leaves {s_L} cell(s) without descent")
    if richdem_strands and s_R == 0:
        errors.append("RichDEM leaves no cell without descent on a fixture chosen "
                      "because it does")
    if s_R == 0 and (n_rounds != 1 or n_changed_r1 != 0):
        errors.append(f"nothing stranded in round 1, yet {n_rounds} rounds and "
                      f"{n_changed_r1} cells changed after round 1")

    max_ulp_R = int(mask_R[valid].max()) if valid.any() else 0
    max_ulp_L = int((ordered_key(L) - ordered_key(F))[valid].max()) if valid.any() else 0
    line = (f"{name:22s} compared {int(valid.sum()):6d}  flat {int(flats.sum()):5d}  "
            f"differing {int(differ.sum()):4d} (guarded {int(guarded.sum()):4d}, "
            f"changed after r1 {n_changed_r1:4d})  without descent RichDEM {s_R:4d} "
            f"lsdtt3 {s_L}  rounds {n_rounds} final pass {n_final_pass:4d}  "
            f"pairs RichDEM {rd_pairs:4d}"
            + (f" (recount {py_R}-{py_R_hi})" if py_R_hi != py_R else "")
            + f" lsdtt3 r1 {r1_pairs:4d} final {ls_pairs:4d}  "
            f"max ULP raise RichDEM {max_ulp_R} lsdtt3 {max_ulp_L}")
    if errors:
        return False, line + "\n    FAIL: " + "\n    FAIL: ".join(errors)
    return True, line


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin-dir", required=True, type=pathlib.Path,
                    help=f"directory containing {EXE_NAME}")
    ap.add_argument("--work-dir", required=True, type=pathlib.Path,
                    help="scratch directory for fixtures and outputs (created)")
    ap.add_argument("--keep", action="store_true",
                    help="keep the rasters and logs after a pass")
    args = ap.parse_args()

    exe = (args.bin_dir / EXE_NAME).resolve()
    if not exe.is_file():
        print(f"FAIL  {exe} not found (build the {EXE_NAME} target)")
        return 1
    root = args.work_dir.resolve() / "flat_resolution_o1"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)

    ok_all = True
    cases = [(n, z, False) for n, z in build_fixtures().items()]
    cases += [(n, z, True) for n, z in build_stranding_fixtures().items()]
    for name, z, strands in cases:
        try:
            ok, line = check_fixture(name, z, exe, root, richdem_strands=strands)
        except RuntimeError as exc:
            ok, line = False, f"{name}: ERROR {exc}"
        ok_all &= ok
        print(line, flush=True)

    if ok_all:
        print("PASS  barnes2014b_float_step matches RichDEM ResolveFlats bit for "
              "bit outside near-zero-guarded flats and the cells lsdtt3 changed "
              "after round 1 (all above RichDEM); lsdtt3 leaves no cell without "
              "descent")
        if not args.keep:
            shutil.rmtree(root)
        return 0
    print(f"FAIL  see above; outputs kept in {root}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
