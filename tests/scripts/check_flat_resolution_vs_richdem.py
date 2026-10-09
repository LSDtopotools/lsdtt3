#!/usr/bin/env python3
"""Oracle O1: lsdtt3's barnes2014b_float_step flat resolution vs RichDEM.

Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3

Copyright (C) 2026 Simon M. Mudd
SPDX-License-Identifier: AGPL-3.0-or-later

WHAT IT PROVES

  `fill_flat_resolution_choice: barnes2014b_float_step` is meant to be
  bit-for-bit RichDEM's `ResolveFlats` on a float32 DEM: Barnes, Lehman and
  Mulla (2014b), "An efficient assignment of drainage direction over flat
  surfaces in raster digital elevation models", Computers & Geosciences 62,
  https://doi.org/10.1016/j.cageo.2013.01.009. Each flat cell is raised by
  `mask` successive float32 steps (nextafter towards +inf), where the mask
  combines the distance to the flat's low edges and away from its high edges.

  For every fixture this script:

    1. fills the DEM with lsdtt-raster-preprocessing (vincent1993,
       fill_min_slope 0), giving F, a surface with TRUE flats;
    2. runs RichDEM ResolveFlats on a float32 rdarray of F, giving R, and
       captures the pair count RichDEM prints;
    3. runs lsdtt-raster-preprocessing on F with barnes2014b_float_step
       (refilling F is a no-op), giving L, and parses its pair count;
    4. asserts L and R have the same bits (uint32 view) on every valid cell,
       except on near-zero-guarded flats (below);
    5. asserts lsdtt3's "Raised above surroundings: <n> pair(s)" equals
       RichDEM's "Cells inappropriately raised above surrounding terrain",
       and that each side's printed count equals a count recomputed here
       from its own output;
    6. asserts both outputs have the same nodata cells, and that F, L and the
       fixture share raster size and geotransform.

  Fixtures are generated here (no network, no committed rasters): hand-built
  valleys, terraces, a mesa, an interior nodata hole, a big flat whose rim is
  1 ULP above it, a valley with a 2-ULP bump (pins RichDEM's `>=` pair rule:
  it prints 8, where a `>` rule would give 7), valley floors at +0, -0, -50
  and 4000 m, and two fixed-seed 256 x 256 integer-metre fractal surfaces
  around 300 m and 3000 m.

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


def check_fixture(name: str, z0: np.ndarray, exe: pathlib.Path,
                  root: pathlib.Path) -> tuple[bool, str]:
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    errors: list[str] = []

    src = d / f"{name}.tif"
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
    if not m:
        return False, f"{name}: 'Raised above surroundings' line not in lsdtt3 stdout"
    ls_pairs, ls_cells = int(m.group(1)), int(m.group(2))

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

    # 4. Bits, outside the guard.
    border = np.zeros(F.shape, bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    flats = flat_cells(F, vF) & ~border
    guarded = flats & (np.abs(F) < GUARD_LIMIT)
    differ = valid & (L.view(np.uint32) != R.view(np.uint32))
    unexpected = differ & ~guarded
    if unexpected.any():
        rr, cc = np.argwhere(unexpected)[0]
        errors.append(f"{int(unexpected.sum())} cell(s) differ outside the guard; "
                      f"first ({rr},{cc}): F {F[rr, cc]!r} L {L[rr, cc]!r} "
                      f"R {R[rr, cc]!r}")

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

    # 5. Pair counts. Each printed count must equal the count recomputed from
    # that side's own output. Across the two tools they must agree except where
    # the guard could change it: the guard raises by up to mask * 2^-23 m, much
    # more than RichDEM's subnormal steps, so a guarded flat could overtop a
    # rim RichDEM does not. The guarded fixtures here have rims 10 m above the
    # floor and max masks below 100 (raise < 1.2e-5 m), so no pair can end
    # level on either side and both counts are 0: they ARE compared.
    py_R, _ = count_raised_above(F, R, valid)
    py_L, py_L_cells = count_raised_above(F, L, valid)
    if py_R != rd_pairs:
        errors.append(f"RichDEM printed {rd_pairs} pairs, recount gives {py_R}")
    if (py_L, py_L_cells) != (ls_pairs, ls_cells):
        errors.append(f"lsdtt3 printed {ls_pairs} pairs/{ls_cells} cells, "
                      f"recount gives {py_L}/{py_L_cells}")
    if ls_pairs != rd_pairs:
        errors.append(f"pair counts differ: lsdtt3 {ls_pairs}, RichDEM {rd_pairs}")

    max_ulp_R = int(mask_R[valid].max()) if valid.any() else 0
    max_ulp_L = int((ordered_key(L) - ordered_key(F))[valid].max()) if valid.any() else 0
    line = (f"{name:22s} compared {int(valid.sum()):6d}  flat {int(flats.sum()):5d}  differing {int(differ.sum()):4d} "
            f"(guarded {int(guarded.sum()):4d})  pairs lsdtt3 {ls_pairs:4d} "
            f"RichDEM {rd_pairs:4d}  max ULP raise RichDEM {max_ulp_R} lsdtt3 {max_ulp_L}")
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
    for name, z in build_fixtures().items():
        try:
            ok, line = check_fixture(name, z, exe, root)
        except RuntimeError as exc:
            ok, line = False, f"{name}: ERROR {exc}"
        ok_all &= ok
        print(line, flush=True)

    if ok_all:
        print("PASS  barnes2014b_float_step matches RichDEM ResolveFlats bit for "
              "bit outside near-zero-guarded flats")
        if not args.keep:
            shutil.rmtree(root)
        return 0
    print(f"FAIL  see above; outputs kept in {root}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
