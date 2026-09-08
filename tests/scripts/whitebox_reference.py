#!/usr/bin/env python3
"""Whitebox Workflows reference implementations for the benchmark suite.

Third sibling of `richdem_reference.py` and `topotoolbox_reference.py`, and the
only one of the three that is a compiled, multi-threaded implementation across
BOTH flow routing and surface metrics.  That is why it was added: the surface
metrics comparison against TopoToolbox is scipy-serial vs OpenMP C++, which
measures the language as much as the algorithm.  Whitebox is Rust with a real
thread pool, so it is a fair speed opponent rather than a strawman.

WHAT IT IS VALID FOR

  Flow routing   -- TIMING, plus VALUE comparison for fill (see below).
  Surface metrics-- TIMING ONLY.  Whitebox's slope/curvature are 3x3 finite
                    difference operators (Horn / Zevenbergen-Thorne lineage);
                    lsdtt3 fits a 6-term quadratic over a user-set window.  At
                    window half-width 1 the work is comparable, but the sign
                    conventions and curvature definitions are not identical and
                    matching them is out of scope.  Never turn a surface row
                    into a value assertion.
  Depressions    -- TIMING ONLY, and read the note on `depression_storage`.

  fill IS value-comparable, and is a better epsilon reference than RichDEM.
  Whitebox's `fill_depressions(fix_flats=True, flat_increment=e)` imposes a
  CONSTANT increment across flats, which is what lsdtt3's `fill_min_slope`
  does.  RichDEM's `epsilon=True` instead raises each cell by one ULP, which
  scales with elevation magnitude -- the reason the pinned `fill_epsilon`
  tolerances in run_flow_benchmark.py differ ~40x between the alpine and
  lakeland fixtures.  Against Whitebox that comparison is apples to apples.

  OBSERVED, UNRESOLVED -- check this before quoting a `fill` value agreement.
  On tests/data/dem_cache/alpine_small.tif (no nodata, 256^2) the raw minimum
  1177.1334 sits at row 254, an INTERIOR cell, and the minimum along the raster
  border is 1178.3820.  A depression fill should therefore bring the global
  minimum to the border minimum.  `fill_depressions_wang_and_liu` does exactly
  that (1178.3820).  `fill_depressions` -- the Barnes tool -- returns 1184.8361,
  6.45 m HIGHER, and the two fills differ by up to 9.78 m.  lsdtt3 is bit-exact
  against both RichDEM and TopoToolbox on this same operation, so the odd one
  out is Whitebox's Barnes tool, most likely a boundary convention (edge cells
  treated as closed).  Until that is understood, treat `fill@whitebox` as a
  TIMING row and let `depression_fill_wl` carry the value comparison.

LICENSING -- READ BEFORE ADDING A ROUTINE

  This script constructs `wb.WbEnvironment()`, the "open tier", which needs no
  licence key.  733 of Whitebox's 775 tools are open tier and every routine
  used here is one of them (verified against 2.0.6).  Do NOT add a routine that
  requires a Pro entitlement, and never call `from_floating_license_id` or
  `from_signed_entitlement_*` here: a benchmark a peer reviewer cannot run
  without buying a licence is not a benchmark.

THREADING -- THE REASON THIS SCRIPT TAKES --threads

  `WbEnvironment.max_procs` defaults to -1, meaning "use every core".  A
  timing table that compares an all-cores Whitebox against a single-threaded
  lsdtt3 is meaningless, so --threads is REQUIRED to be echoed into the JSON
  and is set explicitly on the environment.  The caller pins the same number
  across all four implementations.

  Note the asymmetry that remains and cannot be fixed: RichDEM and TopoToolbox
  are effectively serial for these routines, so their rows are only ever
  1-thread rows.  The scaling sweep is lsdtt3 vs Whitebox.

API NOTES (verified against whitebox-workflows 2.0.6, 2026-08-31)

  These all CHANGED from the 1.x API that most published examples use, and each
  one fails loudly but unhelpfully if you get it wrong:

    * Tools are NOT flat methods on the environment.  `wbe.fill_depressions(...)`
      raises "Flat WbEnvironment tool methods are disabled in this build".
      They live under category namespaces: `wbe.hydrology.fill_depressions`,
      `wbe.terrain.slope`.
    * Category tool calls are KEYWORD-ONLY.  A positional argument raises
      "category tool calls only accept keyword arguments".
    * The first parameter is named `dem` for some tools and `input` for others
      -- `fill_depressions(dem=)` but `d8_flow_accum(input=)` -- so they cannot
      be called through one generic wrapper.  The table below spells each out.
    * Raster metadata moved from `.configs` to `.metadata`.
    * `Raster.to_numpy()` exists, which is how values reach the comparison.

Usage:
    pixi run -e benchmark python tests/scripts/whitebox_reference.py \
        --routine fill --in dem.tif --out ref.tif --threads 1 --repeats 5
"""

import argparse
import contextlib
import json
import os
import shutil
import tempfile
import resource
import statistics
import sys
import time

# -----------------------------------------------------------------------------
# Routine registry
#
# Each entry is (category, tool, kwargs-builder).  The builder takes the loaded
# raster and the CLI args and returns the keyword dict, because the input
# parameter is spelled `dem` for some tools and `input` for others.
# -----------------------------------------------------------------------------

FLOW_ROUTINES = (
    "fill",             # plain fill, true flats -- the tight comparison
    "fill_epsilon",     # constant increment across flats, like lsdtt3
    "breach",           # Lindsay 2016 least-cost breaching
    "accum_d8",
    "accum_dinf",
    "accum_mfd",
)

SURFACE_ROUTINES = (
    "slope",
    "aspect",
    "curvature_profc",
    "curvature_planc",
    "curvature_all2",   # profile + plan back to back, one shared timing
    "hillshade",
)

DEPRESSION_ROUTINES = (
    "depression_sink",      # sink: labels depressions
    "depression_fill_wl",   # Wang & Liu fill -- a second fill lineage
    "depression_fill_wl_epsilon",  # ... with a constant flat gradient
    "basins_from_seeds",    # nested basins, one pass per nesting level
)

ROUTINES = FLOW_ROUTINES + SURFACE_ROUTINES + DEPRESSION_ROUTINES

# Routines whose timing may be RECORDED but whose ratio must never be quoted as
# an lsdtt3 speed claim, because the two tools are not doing the same work.
#
#   aspect  -- Whitebox uses a 3x3 Horn operator; lsdtt3 derives aspect from a
#              fitted quadratic surface.  Different amounts of work by design.
#   hillshade -- a rendering operation, included because users time it, not
#              because it is an algorithmic comparison.
#   depression_sink -- Whitebox labels depressions; lsdtt3's depression
#              hierarchy builds a nested tree with far more structure.  The
#              hierarchy is strictly more work, so a ratio flatters Whitebox.
REPORT_ONLY_ROUTINES = frozenset({
    "aspect", "hillshade", "depression_sink",
})

# All surface routines are timing-only; kept separate from REPORT_ONLY because
# a timing comparison IS meaningful for slope and curvature (same 3x3-scale
# work), it is only the VALUES that must never be asserted.
VALUE_COMPARABLE = frozenset({"fill", "depression_fill_wl",
                              "depression_fill_wl_epsilon"})
# `fill_epsilon` is NOT here: it calls the Barnes tool, which is the one whose
# boundary convention is unresolved (see the module docstring). Use
# `depression_fill_wl_epsilon` for a constant-gradient value comparison.


def peak_rss_mb() -> float:
    """Peak resident set size of this process, in MiB.

    ru_maxrss is kibibytes on Linux and bytes on macOS.  The Linux reading is
    used because that is where the benchmark is run; on macOS the number reads
    ~1024x high, which is obvious rather than subtle.

    This is the peak for the WHOLE process, so it includes the loaded DEM and
    the interpreter, not just the routine.  The caller subtracts the baseline
    reported by `--routine` on a no-op to get the routine's own footprint; see
    `rss_after_load_mb` below, which is measured before the timing loop.
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


@contextlib.contextmanager
def suppressed_output():
    """Silence Whitebox's progress output.

    `wbe.verbose = False` handles most of it, but the Rust layer writes to the
    C-level file descriptors, which would write past a Python-level capture and
    corrupt the single JSON line this script emits.  Redirect the descriptors.
    """
    saved = os.dup(1), os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        os.close(devnull)
        os.close(saved[0])
        os.close(saved[1])


def _nodata_value(raster):
    """Whitebox's nodata sentinel for a raster, or None if it cannot be read.

    Note `metadata` is a METHOD returning a `RasterConfigs`, not a property --
    `raster.metadata.nodata` silently yields a bound-method attribute error
    rather than the value.  (In 1.x this object was reached as `.configs`.)
    The sentinel is -9999.0 by default, NOT NaN, so masking on `np.isfinite`
    alone leaves it in the data.

    Returning None means "do not mask", which is the safe direction: it can
    leave a sentinel in the summary, but it can never silently discard real
    terrain.
    """
    try:
        return float(raster.metadata().nodata)
    except Exception:
        return None


def _make_call(wbe, routine: str, dem, args):
    """Return a zero-argument callable performing one routine.

    *** RESOLVE THE TOOL HANDLE OUTSIDE THE LAMBDA. THIS IS NOT STYLE. ***

    Attribute access on a category namespace -- `wbe.hydrology.sink` -- is not a
    plain attribute read.  It performs a full tool-registry resolution through
    the Rust layer, measured at ~190 ms PER ACCESS (1000 accesses took 191.4 s
    on the development machine, whitebox-workflows 2.0.6).  The namespace itself
    (`wbe.hydrology`) is free: 1000 accesses took 0.8 ms.

    If the lookup happens inside the timed lambda, every Whitebox measurement
    carries a fixed ~190 ms tax that has nothing to do with the algorithm.  On
    the 256^2 fixture that tax IS the entire measurement:

        fill on 256^2, lookup inside the loop   0.207 s
        fill on 256^2, handle resolved once     0.019 s     <- 10.9x

    An earlier version of this file made exactly that mistake, and it inflated
    Whitebox's cost by an order of magnitude in lsdtt3's favour on every small
    fixture.  Binding the handle first is what makes the number the algorithm.
    """
    H = wbe.hydrology
    T = wbe.terrain

    # -- flow --------------------------------------------------------------
    if routine == "fill":
        # fix_flats=False gives TRUE flats, matching RichDEM epsilon=False and
        # lsdtt3 fill_min_slope=0.  This is the comparison that should be
        # bit-exact if all three implement Barnes 2014 faithfully.
        fn = H.fill_depressions
        return lambda: fn(dem=dem, fix_flats=False)
    if routine == "fill_epsilon":
        fn = H.fill_depressions
        return lambda: fn(dem=dem, fix_flats=True,
                          flat_increment=args.flat_increment)
    if routine == "breach":
        # max_dist defaults to 100 cells, and an unbounded search is a
        # different algorithm with a different cost.  Surfaced as a CLI flag so
        # the caller states what it used rather than inheriting it silently.
        fn = H.breach_depressions_least_cost
        return lambda: fn(dem=dem, max_dist=args.breach_max_dist,
                          fill_deps=True)
    if routine == "accum_d8":
        # out_type="cells" not the default "sca": lsdtt3's drainage area raster
        # is compared in cell counts, and specific catchment area would fold in
        # the cell size and hide a routing difference behind a unit conversion.
        fn = H.d8_flow_accum
        return lambda: fn(input=dem, out_type="cells")
    if routine == "accum_dinf":
        fn = H.dinf_flow_accum
        return lambda: fn(input=dem, out_type="cells")
    if routine == "accum_mfd":
        # Freeman FD8 with an exponent, the same family as RichDEM's
        # method="Freeman".  lsdtt3's MFD is Quinn (1991); the three do not
        # agree on proportioning and this row is a timing row.
        fn = H.fd8_flow_accum
        return lambda: fn(dem=dem, out_type="cells",
                          exponent=args.mfd_exponent)

    # -- surface -----------------------------------------------------------
    if routine == "slope":
        fn = T.slope
        return lambda: fn(input=dem, units="degrees")
    if routine == "aspect":
        fn = T.aspect
        return lambda: fn(input=dem)
    if routine == "curvature_profc":
        fn = T.profile_curvature
        return lambda: fn(input=dem)
    if routine == "curvature_planc":
        fn = T.plan_curvature
        return lambda: fn(input=dem)
    if routine == "curvature_all2":
        # What a user actually pays to get both curvatures.  Whitebox, like
        # TopoToolbox, recomputes its derivatives per curvature type; one
        # lsdtt3 polyfit yields every coefficient at once.  Two calls, not
        # four, because Whitebox exposes plan and profile as separate tools but
        # has no direct tangential/mean equivalents in the same family.
        prof = T.profile_curvature
        plan = T.plan_curvature
        return lambda: [prof(input=dem), plan(input=dem)]
    if routine == "hillshade":
        fn = T.hillshade
        return lambda: fn(input=dem, azimuth=args.azimuth,
                          altitude=args.altitude)

    # -- depressions -------------------------------------------------------
    if routine == "depression_sink":
        fn = H.sink
        return lambda: fn(dem=dem)
    if routine == "depression_fill_wl":
        fn = H.fill_depressions_wang_and_liu
        return lambda: fn(dem=dem, fix_flats=False)
    if routine == "depression_fill_wl_epsilon":
        # Wang & Liu WITH a constant gradient across flats -- the same operation
        # lsdtt3's fill_min_slope performs.
        #
        # Routed through wang_and_liu rather than through `fill_epsilon`, which
        # uses Whitebox's Barnes tool. Measured on alpine_small, that tool
        # disagrees with lsdtt3 by 9.78 m -- exactly the offset by which it
        # disagrees with Whitebox's OWN wang_and_liu -- so `fill_epsilon` is
        # contaminated by that tool's boundary convention rather than by any
        # difference in the increment scheme. See VALUE_COMPARABLE below.
        fn = H.fill_depressions_wang_and_liu
        return lambda: fn(dem=dem, fix_flats=True,
                          flat_increment=args.flat_increment)

    raise ValueError(f"unknown routine {routine!r}")


# -----------------------------------------------------------------------------
# WHITEBOX CANNOT READ A DEFLATE-PREDICTOR GeoTIFF.
#
# Measured against whitebox-workflows 2.0.6: a Float32 tiled DEFLATE GeoTIFF
# written with PREDICTOR=2 *or* PREDICTOR=3 is read back as ALL NaN, silently.
# The same data written with no predictor reads correctly, tiled or striped.
# This is a limitation of Whitebox's own TIFF reader, not of the file: GDAL,
# TopoToolbox and RichDEM all read those files without complaint.
#
# It fails in the worst possible way -- no error, no warning, just NaN -- and
# what surfaces downstream is a plausible-looking *result*. It was found because
# `accum_d8@whitebox` reported a maximum flow accumulation of 36 cells against
# lsdtt3's 41,437: Whitebox was faithfully routing a raster of NaN, and it read
# as a routing disagreement rather than as a failure to load the DEM.
#
# lsdtt3 writes tiled DEFLATE with a predictor by design (it is 18% smaller than
# no predictor on float terrain), and it should not give that up for one
# reader's limitation. So the fix belongs here: if the input carries a
# predictor, transcode it to a predictor-free temporary and hand Whitebox that.
# The pixel values are bit-identical -- only the compression of the bytes on
# disk changes -- so nothing about the measurement is affected.
#
# Staging happens before the timing starts, and read time is reported separately
# from compute in any case, so this costs the benchmark nothing.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Nested basins, one pass per nesting level
#
# `unnest_basins` is the literal analogue -- it exists precisely because
# Whitebox cannot put nested basins in ONE raster and emits one per nesting
# level. `watershed` is used here instead because it returns a single
# predictable raster per call, so calling it once per level measures the same
# per-level cost without a variable-length output to handle. Either way the
# shape is the point: cost linear in DEPTH, where lsdtt3 reads a sub-interval.
#
# Do not "optimise" this into a single call with all seeds. That would measure
# the disjoint-basin case, which every tool does in one pass and which
# therefore distinguishes nothing. See run_svector_benchmark.py.
# -----------------------------------------------------------------------------
def _seeds_to_vector(seeds_csv: str, dem_path: str, out_path: str) -> int:
    """Reproject lat/lon seeds into the DEM's CRS and write a point shapefile."""
    import csv as _csv
    from osgeo import gdal, ogr, osr
    gdal.UseExceptions(); osr.UseExceptions()

    ds = gdal.Open(dem_path)
    dst = osr.SpatialReference(); dst.ImportFromWkt(ds.GetProjection())
    ds = None
    src = osr.SpatialReference(); src.ImportFromEPSG(4326)
    for srs in (src, dst):
        try:
            srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except Exception:
            pass
    tr = osr.CoordinateTransformation(src, dst)

    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_path):
        drv.DeleteDataSource(out_path)
    vds = drv.CreateDataSource(out_path)
    layer = vds.CreateLayer("seeds", dst, ogr.wkbPoint)
    layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    n = 0
    with open(seeds_csv) as fh:
        for i, row in enumerate(_csv.DictReader(fh)):
            x, y, *_ = tr.TransformPoint(float(row["longitude"]),
                                         float(row["latitude"]))
            pt = ogr.Geometry(ogr.wkbPoint); pt.AddPoint(x, y)
            feat = ogr.Feature(layer.GetLayerDefn())
            feat.SetGeometry(pt); feat.SetField("id", i)
            layer.CreateFeature(feat); feat = None
            n += 1
    vds = None
    return n


def run_basins_from_seeds(wbe, dem, dem_path: str, seeds_csv: str,
                          levels: int, repeats: int):
    import numpy as np

    stage = tempfile.mkdtemp(prefix="wb_seeds_")
    try:
        shp = os.path.join(stage, "seeds.shp")
        n_seeds = _seeds_to_vector(seeds_csv, dem_path, shp)

        with suppressed_output():
            pour = wbe.read_vector(shp)

        def one_pass():
            # The D8 pointer is built INSIDE the timed region: lsdtt-basin-
            # extraction builds its flow routing on every invocation, so both
            # sides must include it or the comparison is between different
            # operations.
            pntr = wbe.hydrology.d8_pointer(dem=dem)
            out = None
            for _ in range(max(1, levels)):
                out = wbe.hydrology.watershed(d8_pntr=pntr, pour_pts=pour)
            return out

        with suppressed_output():
            one_pass()                              # warm-up, discarded
        times = []
        for _ in range(max(1, repeats)):
            t0 = time.perf_counter()
            with suppressed_output():
                result = one_pass()
            times.append(time.perf_counter() - t0)

        arr = np.asarray(result.to_numpy())
        assert_readable(arr, dem_path)
        return {"reference": "whitebox", "routine": "basins_from_seeds",
                "seconds": min(times), "seconds_all": times,
                "seconds_median": statistics.median(times),
                "levels": levels, "n_seeds": n_seeds,
                "threads": 1, "timing_scope": "in_process_routine",
                "peak_rss_mb": peak_rss_mb()}
    finally:
        shutil.rmtree(stage, ignore_errors=True)


def stage_predictor_free(in_path: str, tmp_dir: str) -> tuple[str, bool]:
    """Return a path Whitebox can read, and whether a copy had to be made."""
    try:
        from osgeo import gdal
        gdal.UseExceptions()
    except Exception:
        # No GDAL here: hand back the original. If it does carry a predictor
        # the guard below still catches the all-NaN read.
        return in_path, False

    ds = gdal.Open(in_path)
    if ds is None:
        return in_path, False
    predictor = (ds.GetMetadata("IMAGE_STRUCTURE") or {}).get("PREDICTOR")
    if predictor in (None, "", "1"):
        ds = None
        return in_path, False

    out = os.path.join(tmp_dir, "wb_nopredictor_" + os.path.basename(in_path))
    gdal.Translate(out, ds, creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])
    ds = None
    return out, True


def assert_readable(arr, in_path: str) -> None:
    """Refuse to time a routine on a raster that did not load.

    An all-NaN read is indistinguishable from a real result once it reaches a
    summary table, so it has to be caught here, at the only point where the
    cause is still visible.
    """
    import numpy as np
    a = np.asarray(arr)
    if a.size and np.isnan(a).all():
        raise RuntimeError(
            f"Whitebox read {in_path} as ALL NaN. This is the DEFLATE-predictor "
            f"limitation described above; the predictor-free staging should have "
            f"prevented it. Do NOT report timings from this run."
        )


def run(routine: str, in_path: str, out_path, args):
    """Time one routine, best-of-N after a discarded warm-up.

    Best-of-N is the headline, matching topotoolbox_reference.py: the fastest
    run is the one least polluted by scheduler noise and page faults, and a
    mean that includes a cold first call measures the allocator.  The median
    and the full list are reported alongside so a single lucky run cannot be
    quoted on its own, and so the caller can show spread.

    I/O is measured SEPARATELY and never folded into `seconds`.  The four
    implementations have very different I/O layers (GDAL, libtopotoolbox,
    Whitebox's own reader) and bundling read time into compute time compares
    file formats rather than algorithms.
    """
    import whitebox_workflows as wb

    wbe = wb.WbEnvironment()
    wbe.verbose = False
    # Pin the thread pool.  Default is -1 = all cores; see the module docstring.
    wbe.max_procs = int(args.threads)

    in_dir = os.path.dirname(os.path.abspath(in_path)) or "."

    # See stage_predictor_free(): a predictor-compressed input reads as all NaN.
    #
    # The staging directory must outlive every TOOL CALL, not just the read.
    # Whitebox's Raster object keeps a handle to the file on disk and the tools
    # re-read it, so deleting the staged copy after read_raster() returns fails
    # later with "GeoTIFF read error: No such file or directory" from inside the
    # timing loop. Hence mkdtemp plus a finally, rather than a with-block around
    # the read.
    stage_dir = tempfile.mkdtemp(prefix="wb_stage_")
    try:
        read_path, staged = stage_predictor_free(in_path, stage_dir)
        wbe.working_directory = (os.path.dirname(os.path.abspath(read_path))
                                 or in_dir)

        t0 = time.perf_counter()
        with suppressed_output():
            dem = wbe.read_raster(os.path.basename(read_path))
        seconds_read = time.perf_counter() - t0
        assert_readable(dem.to_numpy(), in_path)

        # Returns its own dict: it produces no single comparable raster and the
        # shared timing path below assumes one.
        if routine == "basins_from_seeds":
            out = run_basins_from_seeds(wbe, dem, read_path, args.seeds,
                                        args.levels, args.repeats)
            out["seconds_read"] = seconds_read
            out["staged_predictor_free"] = bool(staged)
            return out

        return _run_timed(wbe, routine, dem, in_path, out_path, args,
                          seconds_read, staged)
    finally:
        shutil.rmtree(stage_dir, ignore_errors=True)


def _run_timed(wbe, routine: str, dem, in_path: str, out_path, args,
               seconds_read: float, staged: bool):
    """The timing loop, split out so the staged input stays alive around it."""

    rss_after_load = peak_rss_mb()

    call = _make_call(wbe, routine, dem, args)

    with suppressed_output():
        result = call()                       # warm-up, discarded

    timings = []
    for _ in range(max(1, args.repeats)):
        t0 = time.perf_counter()
        with suppressed_output():
            result = call()
        timings.append(time.perf_counter() - t0)

    rss_peak = peak_rss_mb()

    # A list result (curvature_all2) has no single raster to write; take the
    # first, which is enough for a shape/extent sanity check.  These routines
    # are timing-only anyway.
    primary = result[0] if isinstance(result, list) else result

    seconds_write = None
    if out_path:
        t0 = time.perf_counter()
        with suppressed_output():
            wbe.write_raster(primary, os.path.abspath(out_path))
        seconds_write = time.perf_counter() - t0

    info = {
        "reference": "whitebox",
        # True when the input carried a DEFLATE predictor and had to be
        # transcoded before Whitebox could read it. Recorded so a campaign can
        # be audited for it rather than it being invisible.
        "staged_predictor_free": bool(staged),
        "routine": routine,
        "seconds": min(timings),
        "seconds_median": statistics.median(timings),
        "seconds_all": timings,
        "seconds_stdev": statistics.stdev(timings) if len(timings) > 1 else 0.0,
        "repeats": len(timings),
        "threads": int(args.threads),
        "seconds_read": seconds_read,
        "seconds_write": seconds_write,
        "peak_rss_mb": rss_peak,
        "rss_after_load_mb": rss_after_load,
        # The routine's own peak footprint above the loaded DEM.  Can read 0.0
        # for a cheap routine, which is honest: ru_maxrss only moves when the
        # high-water mark moves.
        "rss_routine_mb": max(0.0, rss_peak - rss_after_load),
        "report_only": routine in REPORT_ONLY_ROUTINES,
        "value_comparable": routine in VALUE_COMPARABLE,
    }

    try:
        import numpy as np
        arr = primary.to_numpy()
        # Mask nodata EXPLICITLY.  Whitebox writes a sentinel value (-9999 by
        # default), not NaN, so `np.isfinite` alone leaves it in and the
        # reported minimum becomes the sentinel rather than the terrain.  That
        # is not cosmetic: `min` is what a reader uses to sanity-check that a
        # fill did not move the global minimum, and -9999 destroys that check.
        nodata = _nodata_value(primary)
        valid = np.isfinite(arr)
        if nodata is not None:
            valid &= (arr != nodata)
        finite = arr[valid]
        info["shape"] = list(arr.shape)
        info["nodata"] = nodata
        info["n_valid"] = int(finite.size)
        info["min"] = float(finite.min()) if finite.size else None
        info["max"] = float(finite.max()) if finite.size else None
    except Exception as exc:
        # A missing numpy view is not a benchmark failure -- the timing is
        # still valid -- so record why rather than aborting the run.
        info["value_summary_error"] = str(exc)

    return info


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--routine", required=True, choices=ROUTINES)
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--out", dest="out_path", default=None,
                    help="optional output raster; omit for a timing-only run")
    ap.add_argument("--threads", type=int, default=1,
                    help="pinned thread count (WbEnvironment.max_procs). "
                         "Default 1, NOT Whitebox's own -1 = all cores, so an "
                         "unconfigured run is comparable rather than flattering.")
    ap.add_argument("--seeds", default="",
                    help="basins_from_seeds: CSV of latitude,longitude seeds")
    ap.add_argument("--levels", type=int, default=1,
                    help="basins_from_seeds: nesting levels to resolve")
    ap.add_argument("--repeats", type=int, default=1,
                    help="timed repeats after a discarded warm-up")
    ap.add_argument("--mfd-exponent", type=float, default=1.0)
    ap.add_argument("--flat-increment", type=float, default=0.0001,
                    help="constant increment across flats for fill_epsilon; "
                         "match lsdtt3's fill_min_slope")
    ap.add_argument("--breach-max-dist", type=int, default=100,
                    help="least-cost breach search radius in cells")
    ap.add_argument("--azimuth", type=float, default=315.0)
    ap.add_argument("--altitude", type=float, default=60.0)
    ap.add_argument("--list", action="store_true",
                    help="print the routine registry as JSON and exit")
    args = ap.parse_args()

    if args.list:
        print(json.dumps({
            "flow": list(FLOW_ROUTINES),
            "surface": list(SURFACE_ROUTINES),
            "depression": list(DEPRESSION_ROUTINES),
            "report_only": sorted(REPORT_ONLY_ROUTINES),
            "value_comparable": sorted(VALUE_COMPARABLE),
        }))
        return 0

    try:
        import whitebox_workflows  # noqa: F401
    except ImportError:
        # Exit 77 is CTest's skip code; the driver propagates it.  An optional
        # reference that is not installed must never turn a run red.
        print("SKIP: whitebox_workflows not importable. Install via "
              "`pixi install -e benchmark` (the `whitebox` feature). "
              "No wheel exists for osx-64; use the benchmark-lite environment "
              "there.", file=sys.stderr)
        return 77

    info = run(args.routine, args.in_path, args.out_path, args)
    # Single JSON line on stdout: the driver parses this, and suppressed_output()
    # guarantees the Rust layer has not polluted the stream.
    print(json.dumps(info))
    return 0


if __name__ == "__main__":
    sys.exit(main())
