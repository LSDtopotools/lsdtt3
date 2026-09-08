#!/usr/bin/env python3
"""Determine which D8 rule Whitebox Workflows uses, against analytic oracles.

WHY THIS EXISTS

  lsdtt3's D8 picks the neighbour of steepest DESCENT -- drop divided by
  distance, so a diagonal must be sqrt(2) times the drop of a cardinal to win.
  RichDEM instead picks the neighbour of lowest ELEVATION
  (richdem/flowmet/d8_flowdirs.hpp, d8_FlowDir), which is a different rule, and
  that is why `run_flow_benchmark.py` records `accum_d8@richdem` as report-only
  rather than asserting a tolerance against it.

  TopoToolbox uses steepest descent and IS asserted.  Whitebox has never been
  characterised either way, so `accum_d8@whitebox` is a timing row only and the
  D8 correctness claim rests on a single independent implementation.  This
  script settles which rule Whitebox uses, so that row can be promoted to a
  value comparison or documented as one that never can be.

  It is a one-shot characterisation, not part of any suite.  It asserts nothing
  about lsdtt3 and cannot fail a build.

WHAT IT MEASURES

  Fixture 1, PLANAR RAMP -- the P8 oracle from tests/cpp/test_flow_properties
  .cpp.  A surface tilted purely north-south: every cell drains down its own
  column, so row r accumulates r+1 cells.  Exact analytic answer, no reference
  implementation needed.

  *** THE RAMP ALONE DOES NOT IDENTIFY THE RULE. ***  On a pure ramp the S, SE
  and SW neighbours are all at the SAME elevation, so a lowest-elevation rule
  faces a tie, and how it breaks that tie decides the answer.  RichDEM breaks
  toward cardinals and therefore PASSES the ramp despite not being steepest
  descent -- the correction recorded at test_flow_properties.cpp:996.  A tool
  passing fixture 1 has been shown to be sane, not to be steepest descent.

  Fixture 2, DISCRIMINATOR -- a surface where the two rules genuinely disagree,
  with no tie to break.  One cell is given a diagonal neighbour that is LOWER
  and a cardinal neighbour that is STEEPER:

      cardinal S  : drop  8 m over 30.0 m  ->  gradient 0.2667   <- steeper
      diagonal SE : drop 10 m over 42.43 m ->  gradient 0.2357   <- lower

  Steepest descent routes S; lowest elevation routes SE.  The cell sits at the
  foot of a long feeder column, so whichever neighbour it chooses receives a
  large accumulation and the other does not.  The difference is unmissable.

  Rather than hand-deriving what every other cell on that surface does, the
  script computes the WHOLE expected accumulation raster twice with its own
  reference D8 -- once per rule -- and classifies Whitebox's output as matching
  one, the other, or neither.  Hand-derivation is where this kind of check
  usually goes wrong.

CALIBRATION, NOT ASSUMPTION

  Whitebox's `out_type="cells"` may or may not count the cell itself.  The
  script does not assume: it reads the offset off fixture 1, where the answer
  is known exactly, and applies that same offset to fixture 2.  If the ramp
  matches neither convention the script says so and stops rather than reporting
  a rule it inferred from a misread raster.

USAGE

    pixi run -e benchmark python tests/scripts/check_whitebox_d8_rule.py

  --keep writes the fixtures and Whitebox's output to a directory for
  inspection.  --threads pins Whitebox's pool (irrelevant to the answer; the
  rule is not thread-dependent, but a pinned run is reproducible).
"""

from __future__ import annotations

import argparse
import contextlib
import io
import math
import os
import pathlib
import sys
import tempfile

import numpy as np

CELL = 30.0

# Direction offsets (dr, dc) and whether each is diagonal.
NEIGHBOURS = [(-1, 0), (-1, 1), (0, 1), (1, 1),
              (1, 0), (1, -1), (0, -1), (-1, -1)]


# -----------------------------------------------------------------------------
# Reference D8, both rules
#
# Deliberately simple and slow. This is an oracle, not a kernel: it exists to
# be obviously correct by inspection, so it never uses a clever traversal.
# -----------------------------------------------------------------------------
def d8_receivers(z: np.ndarray, rule: str) -> np.ndarray:
    """Receiver index per cell, -1 for a sink. rule: 'steepest' or 'lowest'.

    'lowest' reproduces the minimum-elevation family, breaking ties toward
    CARDINALS -- which is what RichDEM does, and is the harder case for this
    check because it passes the ramp.
    """
    n_rows, n_cols = z.shape
    recv = np.full(z.size, -1, dtype=np.int64)

    for r in range(n_rows):
        for c in range(n_cols):
            here = z[r, c]
            best_flat = -1
            best_score = 0.0
            best_is_diagonal = True

            for dr, dc in NEIGHBOURS:
                rr, cc = r + dr, c + dc
                if not (0 <= rr < n_rows and 0 <= cc < n_cols):
                    continue
                drop = here - z[rr, cc]
                if drop <= 0.0:
                    continue

                diagonal = (dr != 0 and dc != 0)
                if rule == "steepest":
                    score = drop / (CELL * math.sqrt(2.0) if diagonal else CELL)
                elif rule == "lowest":
                    score = drop          # elevation difference alone
                else:
                    raise ValueError(rule)

                # Strictly greater takes it; on an exact tie a cardinal wins.
                if (score > best_score
                        or (score == best_score and best_is_diagonal
                            and not diagonal)):
                    best_score = score
                    best_flat = rr * n_cols + cc
                    best_is_diagonal = diagonal

            recv[r * n_cols + c] = best_flat

    return recv


def d8_accumulate(recv: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Cell counts including the cell itself, by draining in elevation order."""
    n = recv.size
    accum = np.ones(n, dtype=np.int64)

    # Kahn's algorithm on the receiver graph: repeatedly drain cells with no
    # remaining donors. O(n) and needs no elevation sort.
    n_donors = np.zeros(n, dtype=np.int64)
    for i in range(n):
        if recv[i] >= 0:
            n_donors[recv[i]] += 1

    stack = [i for i in range(n) if n_donors[i] == 0]
    processed = 0
    while stack:
        i = stack.pop()
        processed += 1
        j = recv[i]
        if j >= 0:
            accum[j] += accum[i]
            n_donors[j] -= 1
            if n_donors[j] == 0:
                stack.append(j)

    if processed != n:
        raise RuntimeError(
            f"receiver graph is not acyclic: {n - processed} cells in a cycle")

    return accum.reshape(shape)


# -----------------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------------
def ramp_surface(n: int = 50) -> np.ndarray:
    """P8: tilted purely north-south. Row r accumulates r+1 cells."""
    z = np.zeros((n, n), dtype=np.float64)
    for r in range(n):
        z[r, :] = 100.0 - r
    return z


def discriminator_surface(n: int = 21) -> tuple[np.ndarray, tuple[int, int]]:
    """A surface on which steepest-descent and lowest-elevation disagree.

    Returns the surface and the (row, col) of the cell where they disagree.

    The base is the same north-south tilt as the ramp, at 2 m per row, so every
    column drains down itself and the focus cell inherits a long feeder column.
    Two cells below the focus are then carved:

        S  neighbour  ->  8 m below the focus   (steeper: 8/30   = 0.2667)
        SE neighbour  -> 10 m below the focus   (lower:  10/42.4 = 0.2357)

    Both carvings are deeper than the 2 m/row base, so neither creates a pit and
    the surface still drains to the southern edge under either rule.
    """
    z = np.zeros((n, n), dtype=np.float64)
    for r in range(n):
        z[r, :] = 1000.0 - 2.0 * r

    fr, fc = n // 2, n // 2
    z_focus = z[fr, fc]
    z[fr + 1, fc] = z_focus - 8.0        # cardinal south: steeper
    z[fr + 1, fc + 1] = z_focus - 10.0   # diagonal south-east: lower

    # Keep the two carved cells draining onward rather than becoming pits.
    for r in range(fr + 2, n):
        z[r, fc] = min(z[r, fc], z[fr + 1, fc] - 2.0 * (r - fr - 1))
        z[r, fc + 1] = min(z[r, fc + 1], z[fr + 1, fc + 1] - 2.0 * (r - fr - 1))

    return z, (fr, fc)


def write_geotiff(z: np.ndarray, path: str) -> None:
    from osgeo import gdal, osr
    gdal.UseExceptions()
    n_rows, n_cols = z.shape
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, n_cols, n_rows, 1, gdal.GDT_Float32)
    # North-up: origin at the top-left, negative north-south pixel size.
    ds.SetGeoTransform((0.0, CELL, 0.0, n_rows * CELL, 0.0, -CELL))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32645)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.WriteArray(z.astype(np.float32))
    band.SetNoDataValue(-9999.0)
    ds.FlushCache()
    ds = None


@contextlib.contextmanager
def suppressed_output():
    """Whitebox prints a banner per call regardless of wbe.verbose."""
    buf_out, buf_err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(buf_out), contextlib.redirect_stderr(buf_err):
        yield


def whitebox_pointer(path: str, threads: int) -> np.ndarray:
    """Whitebox's D8 receiver direction per cell.

    A direct reading of the routing decision, which the accumulation raster can
    only imply. Kept as a SECOND, independent confirmation rather than the sole
    evidence, because `d8_flow_accum` is the routine the benchmark actually
    times and a pointer that disagreed with it would itself be a finding.

    The direction CODES are not assumed. They are calibrated on the ramp, where
    every interior cell must route due south, so whichever code covers the ramp
    IS south. That avoids depending on which of the two published Whitebox
    pointer conventions this build uses.
    """
    import whitebox_workflows as wb
    wbe = wb.WbEnvironment()
    wbe.verbose = False
    wbe.max_procs = int(threads)
    wbe.working_directory = os.path.dirname(os.path.abspath(path)) or "."
    with suppressed_output():
        dem = wbe.read_raster(os.path.basename(path))
        ptr = wbe.hydrology.d8_pointer(dem=dem)
    return np.asarray(ptr.to_numpy(), dtype=np.float64)


def whitebox_accum(path: str, threads: int, out_path: str | None) -> np.ndarray:
    import whitebox_workflows as wb
    wbe = wb.WbEnvironment()
    wbe.verbose = False
    wbe.max_procs = int(threads)
    wbe.working_directory = os.path.dirname(os.path.abspath(path)) or "."
    with suppressed_output():
        dem = wbe.read_raster(os.path.basename(path))
        # out_type="cells", exactly as whitebox_reference.py uses it: specific
        # catchment area would fold in the cell size and hide a routing
        # difference behind a unit conversion.
        acc = wbe.hydrology.d8_flow_accum(input=dem, out_type="cells")
        if out_path:
            wbe.write_raster(acc, os.path.abspath(out_path))
    return np.asarray(acc.to_numpy(), dtype=np.float64)


# -----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--keep", default="",
                    help="directory to keep the fixtures and outputs in")
    args = ap.parse_args()

    try:
        import whitebox_workflows  # noqa: F401
        from osgeo import gdal  # noqa: F401
    except Exception as exc:
        print(f"SKIP: needs whitebox_workflows and GDAL ({exc})")
        print("      pixi run -e benchmark python tests/scripts/"
              "check_whitebox_d8_rule.py")
        return 77

    workdir = args.keep or tempfile.mkdtemp(prefix="wb_d8_rule_")
    pathlib.Path(workdir).mkdir(parents=True, exist_ok=True)
    print(f"working directory: {workdir}\n")

    # -- Fixture 1: the ramp, and the counting-convention calibration --------
    print("=" * 72)
    print("FIXTURE 1 -- planar ramp (analytic oracle)")
    print("=" * 72)
    n = 50
    z = ramp_surface(n)
    ramp_path = os.path.join(workdir, "ramp.tif")
    write_geotiff(z, ramp_path)

    wb_acc = whitebox_accum(ramp_path, args.threads,
                            os.path.join(workdir, "ramp_wb_accum.tif")
                            if args.keep else None)

    # Analytic: interior columns only. The outermost columns have their upslope
    # set clipped by the edge, exactly as the C++ P8 test excludes them.
    rows = np.arange(n).reshape(-1, 1)
    expected = np.repeat(rows + 1, n, axis=1).astype(np.float64)
    interior = np.zeros((n, n), dtype=bool)
    interior[:, 1:-1] = True

    offset = None
    for cand in (0, -1):
        if np.allclose(wb_acc[interior], (expected + cand)[interior], atol=1e-6):
            offset = cand
            break

    if offset is None:
        diff = wb_acc[interior] - expected[interior]
        print("  RESULT: Whitebox does NOT reproduce the analytic ramp under "
              "either counting convention.")
        print(f"    observed minus expected: min {diff.min():.3f} "
              f"max {diff.max():.3f} mean {diff.mean():.3f}")
        print("    Nothing further can be concluded; the discriminator would "
              "be read against a raster this script does not understand.")
        return 1

    convention = ("includes the cell itself" if offset == 0
                  else "counts UPSLOPE cells only, excluding itself")
    print(f"  ramp reproduced exactly.  out_type='cells' {convention}.")
    print("  NOTE: passing the ramp does NOT establish steepest descent -- on a")
    print("        ramp the S/SE/SW neighbours are tied in elevation, and a")
    print("        lowest-elevation rule that breaks ties toward cardinals")
    print("        passes it too (RichDEM does exactly that).")
    print()

    # -- Fixture 2: the discriminator ---------------------------------------
    print("=" * 72)
    print("FIXTURE 2 -- discriminator (no tie: lower diagonal vs steeper cardinal)")
    print("=" * 72)
    m = 21
    zd, (fr, fc) = discriminator_surface(m)
    disc_path = os.path.join(workdir, "discriminator.tif")
    write_geotiff(zd, disc_path)

    cardinal_grad = 8.0 / CELL
    diagonal_grad = 10.0 / (CELL * math.sqrt(2.0))
    print(f"  focus cell            : row {fr}, col {fc}")
    print(f"  cardinal S  drop  8 m : gradient {cardinal_grad:.4f}  <- steeper")
    print(f"  diagonal SE drop 10 m : gradient {diagonal_grad:.4f}  <- lower")
    print()

    exp_steep = d8_accumulate(d8_receivers(zd, "steepest"), zd.shape).astype(float)
    exp_lowest = d8_accumulate(d8_receivers(zd, "lowest"), zd.shape).astype(float)

    if np.array_equal(exp_steep, exp_lowest):
        print("  ERROR: the two rules produce the same raster on this fixture, "
              "so it discriminates nothing. The fixture is wrong, not the tool.")
        return 1

    s_cell = (fr + 1, fc)
    se_cell = (fr + 1, fc + 1)
    print("  reference expectations at the two candidate receivers "
          "(cells, self included):")
    print(f"    steepest descent : S {exp_steep[s_cell]:6.0f}   "
          f"SE {exp_steep[se_cell]:6.0f}")
    print(f"    lowest elevation : S {exp_lowest[s_cell]:6.0f}   "
          f"SE {exp_lowest[se_cell]:6.0f}")
    print(f"    the two rasters differ on "
          f"{int((exp_steep != exp_lowest).sum())} of {zd.size} cells")
    print()

    wb_disc = whitebox_accum(disc_path, args.threads,
                             os.path.join(workdir, "discriminator_wb_accum.tif")
                             if args.keep else None)

    # Compare on the interior only: edge behaviour is a boundary convention,
    # not a routing rule, and the two references share this script's convention
    # rather than Whitebox's.
    inner = np.zeros(zd.shape, dtype=bool)
    inner[1:-1, 1:-1] = True

    match_steep = np.allclose(wb_disc[inner], (exp_steep + offset)[inner], atol=1e-6)
    match_lowest = np.allclose(wb_disc[inner], (exp_lowest + offset)[inner], atol=1e-6)

    print(f"  Whitebox at S  {wb_disc[s_cell]:6.0f}     "
          f"at SE {wb_disc[se_cell]:6.0f}")

    # Direct confirmation from the pointer, with the codes calibrated on the
    # ramp rather than assumed.
    pointer_verdict = ""
    try:
        ptr_ramp = whitebox_pointer(ramp_path, args.threads)
        codes, counts = np.unique(ptr_ramp[1:-1, 1:-1], return_counts=True)
        if len(codes) == 1:
            south = int(codes[0])
            ptr_disc = whitebox_pointer(disc_path, args.threads)
            focus_code = int(ptr_disc[fr, fc])
            went_south = (focus_code == south)
            print()
            print(f"  d8_pointer: due south = code {south} "
                  f"(all {int(counts[0])} interior ramp cells agree)")
            print(f"              focus cell code {focus_code} -> "
                  f"{'SOUTH, the cardinal' if went_south else 'NOT south'}")
            pointer_verdict = "steepest" if went_south else "lowest"
        else:
            print("\n  d8_pointer: ramp did not route uniformly south "
                  f"(codes {dict(zip(codes.astype(int), counts))}); "
                  "codes could not be calibrated, so the pointer is not used.")
    except Exception as exc:
        print(f"\n  d8_pointer unavailable ({type(exc).__name__}); "
              "relying on the accumulation comparison alone.")
    print()
    print("=" * 72)
    if pointer_verdict and (
            (pointer_verdict == "steepest") != bool(match_steep and not match_lowest)):
        print("RESULT: INCONCLUSIVE -- the accumulation raster and the pointer")
        print("        raster imply DIFFERENT rules. That is itself a finding")
        print("        and must be understood before either is used as a")
        print("        reference. Do not promote this row.")
        print("=" * 72)
        return 1

    if match_steep and not match_lowest:
        print("RESULT: Whitebox uses STEEPEST DESCENT (drop / distance).")
        if pointer_verdict == "steepest":
            print("        Confirmed twice: whole-raster accumulation AND the")
            print("        d8_pointer reading at the focus cell.")
        print()
        print("  Same rule as lsdtt3 and TopoToolbox. `accum_d8@whitebox` can")
        print("  be promoted from a timing row to a VALUE comparison in")
        print("  run_flow_benchmark.py's COMPARISONS, giving the D8 correctness")
        print("  claim a second independent reference.")
        rc = 0
    elif match_lowest and not match_steep:
        print("RESULT: Whitebox uses LOWEST ELEVATION, not steepest descent.")
        print()
        print("  Same family as RichDEM. `accum_d8@whitebox` must stay a TIMING")
        print("  row: pinning a tolerance against it would enshrine a different")
        print("  routing rule as if it were lsdtt3 error. Say so in the paper.")
        rc = 0
    else:
        n_vs_steep = int((~np.isclose(wb_disc[inner],
                                      (exp_steep + offset)[inner])).sum())
        n_vs_lowest = int((~np.isclose(wb_disc[inner],
                                       (exp_lowest + offset)[inner])).sum())
        print("RESULT: Whitebox matches NEITHER reference rule exactly.")
        print(f"  cells differing from steepest descent : {n_vs_steep}")
        print(f"  cells differing from lowest elevation : {n_vs_lowest}")
        print()
        print("  Most likely a third tie-breaking or flat-handling convention.")
        print("  It stays a timing row, and the paper should say the rule was")
        print("  tested and found to match neither, which is a stronger")
        print("  statement than leaving it uncharacterised.")
        rc = 0
    print("=" * 72)

    if args.keep:
        print(f"\nfixtures and outputs kept in {workdir}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
