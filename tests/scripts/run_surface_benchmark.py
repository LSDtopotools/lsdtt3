#!/usr/bin/env python3
"""Surface-metrics SPEED benchmark driver: lsdtt-surface-metrics vs TopoToolbox.

Phases 3 and 4 of the surface-metrics testing. Sibling of
run_flow_benchmark.py, and deliberately built the same way: --fixture,
--comparison, --check, --work-dir, and exit 77 (CTest's skip code) whenever a
dependency or the network is missing, so an optional benchmark can never turn
CI red.

WHAT THIS IS AND IS NOT
-----------------------
**TopoToolbox is a PERFORMANCE reference here, not a correctness reference.**
Correctness for the polyfit chain is handled by the analytic Taylor oracle in
`lsdtt_raster_synthetic_surfaces_cpu.hpp` (Tier A), which is exact, offline
and covers all 27 derived metrics -- 23 of which have no TopoToolbox
counterpart at all.

So there are NO value-accuracy assertions in this driver. Sign conventions,
aspect handedness and flat-ground masking deliberately differ between the two
tools and matching them is explicitly out of scope. The one correctness-shaped
number computed here is a rank correlation, and it exists purely as a GUARD:
near-zero correlation means the two tools computed unrelated things and the
timing is between unlike computations, so it must not be published. A high
correlation with the opposite sign is fine and expected.

EQUIVALENCE OF WORK -- a timing number that violates these is invalid
-------------------------------------------------------------------------------
1. Same window support. lsdtt3's half-width is
   `W = max(1, round(polyfit_lambda_metres / (2 * cell_size)))`, giving a
   (2W+1)^2 window; the Evans kernels TopoToolbox uses are fixed 3x3. This
   driver sets `polyfit_lambda_metres = 2 * cell_size` explicitly to force
   W = 1 rather than relying on the default 70 m happening to round to 1 on
   30 m data, and reports the W it derived in every row.
2. Same number of output quantities. lsdtt3 is asked for exactly the matching
   `polyfit_metrics`, never `all`.
3. `modified=False` on evansslope (see topotoolbox_reference.py).
4. Valid-pixel count recorded on both sides. A large asymmetry means one tool
   did less work; it is reported, not averaged away.
5. Float precision recorded on both sides.

THREADING -- READ BEFORE QUOTING ANY NUMBER
-------------------------------------------
`evansslope` and `curvature` are single-threaded `scipy.ndimage.convolve`
calls. Presenting an OpenMP-threaded lsdtt3 run against them as an algorithmic
win would be misleading, so this driver measures lsdtt3 TWICE:

    lsdtt_seconds_1t    n_threads = 1   <-- PRIMARY, the honest comparison
    lsdtt_seconds_mt    n_threads = 0   <-- secondary, what a user actually gets

Both are written to the CSV and both are printed, each labelled. `speedup_1t`
is the number to quote; `speedup_mt` is a throughput statement about the tool,
not about the algorithm.

WHAT IS ACTUALLY TIMED
----------------------
Not process start-up, not GDAL open, not the raster write.
`lsdtt-surface-metrics` already prints an in-process timing for the fit itself
("Running direct fit... done (X s)", cli/lsdtt-surface-metrics.cpp), which
brackets exactly `lsdtt::direct_polyfit` and nothing else. This driver parses
that line, so no C++ harness is needed and the comparison is compute-to-compute.

The one exception is hillshade, which the CLI does not time internally. There
the number is whole-process wall time minus a baseline run that opens the same
DEM and writes no raster, so it still excludes start-up and GDAL open but DOES
include one GeoTIFF write -- which on this fixture costs an order of magnitude
more than the reference's entire compute. Those rows are therefore REPORT-ONLY:
a ratio built on them would be a filesystem measurement. They carry an
`lsdtt_est_compute_1t` column that subtracts a separately measured per-write
cost, and that is an estimate, not a timing. Adding a `std::chrono` bracket
around `compute_hillshade` in cli/lsdtt-surface-metrics.cpp -- matching the one
already around `direct_polyfit` -- would promote them to real numbers; it is
the only comparison in this file blocked on a source change.

Rows carry a `timing_basis` column saying which basis they used; never compare
a number across bases.

Usage:
    pixi run -e benchmark python tests/scripts/run_surface_benchmark.py \\
        --exe build/lsdtt-surface-metrics --fixture alpine_medium
    python tests/scripts/run_surface_benchmark.py --list
"""

import argparse
import csv
import datetime
import hashlib
import json
import pathlib
import re
import shutil
import statistics
import subprocess
import sys
import time
import lsdtt_test_env

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

# Reuse the flow driver's fixture registry, fetch logic, raster loader and
# environment scrubbing rather than forking a second copy. Importing it is
# side-effect free (everything is behind `if __name__ == "__main__"`), and it
# means a fixture fetched for one benchmark is shared with the other.
import run_flow_benchmark as flow  # noqa: E402

SKIP_EXIT = 77
die_skip = flow.die_skip
system_env = flow.system_env
load_raster = flow.load_raster

REFERENCE_SCRIPT = SCRIPT_DIR / "topotoolbox_reference.py"

# The polyfit is what the CLI times internally, and it is the expensive step.
FIT_TIMING_RE = re.compile(
    r"(?:Running (direct|wavelet) fit|Computed (hillshade))\.\.\. done "
    r"\(([-+0-9.eE]+) s\)")

# Correlation below which a timing comparison is declared meaningless: the two
# tools were not computing the same kind of quantity, so the stopwatch is
# measuring unlike work. Deliberately loose -- this is a sanity floor, not an
# accuracy tolerance. Sign is discarded before comparing (|rho|), because
# opposite curvature sign conventions are expected and harmless.
CORRELATION_GUARD = 0.5

# All 27 derived metrics, in the order the CLI reference lists them. Used for
# the marginal-cost-per-extra-metric measurement, which is the number most
# likely to matter to a user: one lsdtt3 polyfit yields every one of these,
# while TopoToolbox recomputes five convolutions for each curvature type.
ALL_METRICS = [
    "fitted_surface", "dzdx", "dzdy", "d2zdx2", "d2zdxdy", "d2zdy2",
    "slope", "slope_degrees", "aspect_radians", "aspect_degrees",
    "laplacian", "hessian_determinant", "mean_curvature", "gaussian_curvature",
    "maximum_curvature", "minimum_curvature", "curvedness", "shape_index",
    "profile_curvature", "plan_curvature", "tangential_curvature",
    "difference_curvature", "horizontal_excess_curvature",
    "vertical_excess_curvature", "accumulation_curvature", "ring_curvature",
    "rotor",
]

FOUR_METRICS = ["slope", "profile_curvature", "plan_curvature",
                "tangential_curvature"]

# -----------------------------------------------------------------------------
# Comparisons.
#
# `note` is printed and stored, because a bare ratio with no statement of what
# was compared is exactly the kind of benchmark number that gets misquoted.
#
# ASYMMETRY THAT FAVOURS TOPOTOOLBOX AND MUST BE STATED: lsdtt3's direct polyfit
# always solves for all SIX quadratic coefficients, whatever single metric is
# requested. TopoToolbox's `evansslope` runs two convolutions, `curvature` five.
# So on the single-metric rows lsdtt3 is doing at least as much work as the
# reference and usually more; the ratio is therefore conservative, not flattered.
# The `curvature_all4` row is where that structure pays off and is the fair
# whole-task comparison.
# -----------------------------------------------------------------------------
COMPARISONS = {
    "slope": {
        "routine": "evansslope",
        "metrics": ["slope"],
        "suffix": "_direct_polyfit_slope.tif",
        "note": ("Evans 3x3 quadratic fit both sides -- the closest "
                 "equivalence of work in the API"),
    },
    "dzdx": {
        "routine": "evansslope_partials",
        "metrics": ["dzdx"],
        "suffix": "_direct_polyfit_dzdx.tif",
        "note": "first partial derivative; evansslope(partial_derivatives=True)",
    },
    "profile_curvature": {
        "routine": "curvature_profc",
        "metrics": ["profile_curvature"],
        "suffix": "_direct_polyfit_profile_curvature.tif",
        "note": "Evans 3x3 second derivatives; opposite sign convention expected",
    },
    "plan_curvature": {
        "routine": "curvature_planc",
        "metrics": ["plan_curvature"],
        "suffix": "_direct_polyfit_plan_curvature.tif",
        "note": "Evans 3x3 second derivatives",
    },
    "tangential_curvature": {
        "routine": "curvature_tangc",
        "metrics": ["tangential_curvature"],
        "suffix": "_direct_polyfit_tangential_curvature.tif",
        "note": "Evans 3x3 second derivatives",
    },
    "mean_curvature": {
        "routine": "curvature_meanc",
        "metrics": ["mean_curvature"],
        "suffix": "_direct_polyfit_mean_curvature.tif",
        "note": "Evans 3x3 second derivatives",
    },
    "curvature_all4": {
        "routine": "curvature_all4",
        "metrics": ["profile_curvature", "plan_curvature",
                    "tangential_curvature", "mean_curvature"],
        "suffix": "_direct_polyfit_mean_curvature.tif",
        "note": ("four curvature types: ONE lsdtt3 polyfit vs FOUR TopoToolbox "
                 "calls (20 convolutions) -- the structural comparison"),
    },
    "aspect": {
        "routine": "aspect",
        "metrics": ["aspect_degrees"],
        "suffix": "_direct_polyfit_aspect_degrees.tif",
        "report_only": True,
        # Two independent reasons this row can never be quoted as a ratio:
        # (a) np.gradient central differences are strictly less work than a
        #     quadratic fit, and (b) aspect is a CIRCULAR quantity, so the rank
        #     correlation guard is meaningless on it -- 359 deg and 1 deg are
        #     adjacent directions but maximally distant ranks. The ~0.5 rho
        #     this row reports is an artefact of the wraparound, not evidence
        #     of disagreement.
        "note": ("REPORT ONLY: TopoToolbox aspect is np.gradient central "
                 "differences, strictly less work than a quadratic fit; and "
                 "the correlation guard is meaningless on a circular quantity"),
    },
    # -------------------------------------------------------------------------
    # HILLSHADE AZIMUTH HANDEDNESS -- measured, not assumed.
    #
    # lsdtt3 at azimuth 315 has Spearman rho = +1.0000 against TopoToolbox at
    # azimuth 45, and only +0.13 against TopoToolbox at 315. Sweeping the
    # reference azimuth gives the whole picture:
    #
    #   ttb azimuth   0     45     90    135    180    225    270    315
    #   rho        +0.75  +1.000 +0.73  +0.06  -0.55  -0.76  -0.49  +0.13
    #
    # A second sweep with lsdtt3 at azimuth 45 peaks at ttb 135, also with
    # rho = 1.0000. Two points, and they discriminate between the two obvious
    # explanations: a mirror about north would have put lsdtt3 45 at ttb 315,
    # and it does not. The relation is a ROTATION:
    #
    #     ttb_azimuth = lsdtt3_azimuth + 90  (mod 360)
    #
    # So lsdtt3 keeps the CLI default 315 and the reference is asked for 45,
    # which is the same physical light direction and makes the comparison
    # like-for-like.
    #
    # This is a by-product of the timing work, not a correctness assertion --
    # but it is exactly the trap a user comparing the two tools would fall
    # into, so it belongs in docs/developer_notes/testing_and_benchmarking.md.
    #
    # REPORT-ONLY, and not because of handedness. lsdtt-surface-metrics does
    # not time its hillshade step internally the way it times the polyfit, so
    # the only number obtainable without touching cli/ is wall time minus a
    # no-output baseline. That still contains one GeoTIFF write, which on this
    # fixture costs an order of magnitude MORE than the reference's entire
    # compute. A ratio built on it would be a filesystem measurement wearing a
    # benchmark's clothes. `lsdtt_seconds_1t_est_compute` subtracts a measured
    # per-write cost to give an indicative figure; treat it as an estimate.
    #
    # AZIMUTH: both sides now take the SAME azimuth. Until 2026-07-22 this
    # spec passed 45 to TopoToolbox against lsdtt3's 315, to compensate for a
    # measured `ttb_az = lsdtt3_az + 90` offset that was recorded as a
    # convention difference. It was not - it was a defect in lsdtt3's
    # compute_hillshade (`A = az - pi/2` where the physics gives `A = az - pi`),
    # since fixed. The compensation would now introduce the error it used to
    # cancel, so it is gone.
    #
    # RESOLVED 2026-07-22: cli/lsdtt-surface-metrics.cpp now prints
    # "Computed hillshade... done (X s)" around compute_hillshade, matching the
    # polyfit timer, so both sides are in-process compute-only and these rows
    # are no longer report-only. This is the only compiled-vs-compiled
    # comparison available (libtopotoolbox fused hillshade), so it carries the
    # most weight of any row here.
    # -------------------------------------------------------------------------
    "hillshade": {
        "routine": "hillshade_fused",
        "metrics": [],
        "params": {"write_hillshade_raster": "true",
                   "hillshade_azimuth": "315", "hillshade_altitude": "60",
                   "hillshade_z_factor": "1.0"},
        "ref_args": ["--azimuth", "315", "--altitude", "60"],
        "suffix": "_hillshade.tif",
        "timing_basis": "in_process",
        "report_only": False,
        "note": ("compiled vs compiled (libtopotoolbox fused hillshade), but "
                 "in-process on both sides since the CLI gained a timer "
                 "around compute_hillshade (2026-07-22). "
                 "same azimuth on both sides since the lsdtt3 sun-azimuth defect was fixed 2026-07-22"),
    },
    "hillshade_unfused": {
        "routine": "hillshade_unfused",
        "metrics": [],
        "params": {"write_hillshade_raster": "true",
                   "hillshade_azimuth": "315", "hillshade_altitude": "60",
                   "hillshade_z_factor": "1.0"},
        "ref_args": ["--azimuth", "315", "--altitude", "60"],
        "suffix": "_hillshade.tif",
        "timing_basis": "in_process",
        "report_only": False,
        "note": ("numpy hillshade path, trades memory for speed vs fused; "
                 "in-process on both sides, as `hillshade`"),
    },
}

# Window half-widths for the O(n) claim in the CLI reference: the `direct` fit
# is advertised as a separable rolling-moment algorithm whose cost is
# independent of window size. W = 40 is an 81x81 window, 6561 cells -- a naive
# implementation would be ~4000x slower than W = 1 there, so the claim is not
# subtle and a flat curve is real evidence.
WINDOW_HALF_WIDTHS = (1, 5, 15, 40)

CSV_FIELDS = [
    "timestamp", "git_sha", "fixture", "fixture_sha", "n_cells",
    "comparison", "reference", "routine", "report_only", "timing_basis",
    "window_half_width", "lambda_metres", "cell_size", "n_metrics",
    "lsdtt_seconds_1t", "lsdtt_median_1t", "lsdtt_seconds_mt", "lsdtt_threads_mt",
    "lsdtt_median_mt", "ref_seconds", "ref_median", "ref_seconds_nan_fill",
    "ref_seconds_comparable", "ref_n_convolutions",
    "write_cost_seconds", "lsdtt_est_compute_1t",
    "speedup_1t", "speedup_mt",
    "abs_correlation", "n_valid_lsdtt", "n_valid_ref",
    "lsdtt_dtype", "ref_dtype",
    "lsdtt_peak_rss_mb", "ref_peak_rss_mb", "note",
]


def git_sha() -> str:
    """Short SHA of the working tree, suffixed -dirty when it is.

    Every row carries it because a timing table without a commit is a number
    nobody can reproduce or bisect.
    """
    try:
        sha = lsdtt_test_env.run(["git", "-C", str(REPO_ROOT), "rev-parse",
                              "--short", "HEAD"],
                             capture_output=True, text=True).stdout.strip()
        dirty = lsdtt_test_env.run(["git", "-C", str(REPO_ROOT), "status",
                                "--porcelain"],
                               capture_output=True, text=True).stdout.strip()
        return f"{sha}-dirty" if dirty else sha
    except Exception:
        return "unknown"


def raster_shape(path: pathlib.Path):
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(str(path))
    return ds.RasterYSize, ds.RasterXSize


def band_dtype(path: pathlib.Path) -> str:
    """On-disk band type. NOT what load_raster returns.

    load_raster casts everything to float64 for comparison arithmetic, so
    asking it for the dtype would report float64 for every raster and quietly
    hide a genuine precision difference between the two tools -- which must
    be stated next to the timing.
    """
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(str(path))
    return gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)


def cell_size_of(path: pathlib.Path) -> float:
    from osgeo import gdal
    gdal.UseExceptions()
    gt = gdal.Open(str(path)).GetGeoTransform()
    return abs(gt[1])


def window_half_width(lambda_m: float, cell_size: float) -> int:
    """Mirror of lsdtt_raster_surface_metrics_cpu.cpp:2499.

    Duplicated here on purpose so the driver can REPORT the window it believes
    lsdtt3 used. If the C++ formula ever changes, this comparison silently
    stops being like-for-like, so it needs a separate unit test
    asserting W == 1 at 30 m / lambda 70 rather than relying on this copy.
    """
    return max(1, int(round(lambda_m / (2.0 * cell_size))))


# -----------------------------------------------------------------------------
# Running lsdtt3
# -----------------------------------------------------------------------------

def run_lsdtt_once(exe, fixture, out_dir, params):
    """One CLI invocation. Returns (wall_seconds, fit_seconds_or_None, rss_mb).

    `fit_seconds` is the CLI's own in-process timing of the polynomial fit and
    is the number the benchmark actually uses -- wall time here includes
    start-up, GDAL open and every raster write, none of which is the compute
    being compared.

    Peak RSS comes from /usr/bin/time rather than getrusage(RUSAGE_CHILDREN),
    which reports a running maximum over ALL children ever reaped and would
    silently attribute the largest earlier run's footprint to this one.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    base = {
        "read_path": str(fixture.parent),
        "read_fname": fixture.name,
        "write_path": str(out_dir),
        "write_fname": fixture.stem,
        # The fixture is already in metric UTM. Letting the CLI reproject would
        # score the two tools on different rasters.
        "enforce_projection": "false",
        "write_output_format": "gtiff",
    }
    base.update(params)

    param_path = out_dir / "surface_benchmark.param"
    param_path.write_text("".join(f"{k}: {v}\n" for k, v in base.items()))

    cmd = [str(exe), str(out_dir), param_path.name]
    rss_mb = float("nan")
    time_exe = pathlib.Path("/usr/bin/time")
    if time_exe.exists():
        cmd = [str(time_exe), "-f", "LSDTT_MAXRSS_KB %M"] + cmd

    t0 = time.perf_counter()
    proc = lsdtt_test_env.run(cmd, capture_output=True, text=True, env=system_env())
    wall = time.perf_counter() - t0

    if proc.returncode != 0:
        raise RuntimeError(
            f"{pathlib.Path(exe).name} failed (exit {proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-2000:]}\n"
            f"--- stderr ---\n{proc.stderr[-2000:]}")

    m = re.search(r"LSDTT_MAXRSS_KB (\d+)", proc.stderr)
    if m:
        rss_mb = int(m.group(1)) / 1024.0

    # Group 1 = polyfit kind, group 2 = "hillshade", group 3 = seconds.
    fits = {(a or b): float(secs)
            for a, b, secs in FIT_TIMING_RE.findall(proc.stdout)}
    # A hillshade-only run has no polyfit line and vice versa, so preferring
    # "direct" and falling back to "hillshade" picks the right timer without
    # the caller having to say which step it asked for.
    return wall, fits.get("direct", fits.get("hillshade")), rss_mb


def time_lsdtt(exe, fixture, out_dir, params, repeats):
    """Best-of-N in-process fit timing after a discarded warm-up.

    Returns dict with best/median fit seconds, best/median wall seconds and
    peak RSS. Best-of-N because the fastest run is the one least polluted by
    the scheduler and the page cache; the median is carried too so a single
    lucky run cannot be quoted alone.
    """
    run_lsdtt_once(exe, fixture, out_dir, params)      # warm-up, discarded

    walls, fits, rss = [], [], []
    for _ in range(max(1, repeats)):
        w, f, r = run_lsdtt_once(exe, fixture, out_dir, params)
        walls.append(w)
        if f is not None:
            fits.append(f)
        if r == r:
            rss.append(r)

    return {
        "wall_best": min(walls),
        "wall_median": statistics.median(walls),
        "fit_best": min(fits) if fits else None,
        "fit_median": statistics.median(fits) if fits else None,
        "peak_rss_mb": max(rss) if rss else float("nan"),
    }


def assert_metrics_written(out_dir, stem, metrics):
    """Every requested metric must have produced a raster.

    Guards the equivalence of work directly rather than trusting the CLI: an
    unrecognised metric name is dropped SILENTLY (no warning, exit 0), so
    without this check a typo -- or the whitespace bug worked around in
    polyfit_params -- would show up as lsdtt3 being impressively fast at
    computing four metrics while it had actually computed one.
    """
    missing = [m for m in metrics
               if not (out_dir / f"{stem}_direct_polyfit_{m}.tif").exists()]
    if missing:
        raise RuntimeError(
            f"lsdtt3 was asked for {len(metrics)} metrics but wrote no raster "
            f"for {missing}. It exits 0 and prints no warning when a metric "
            f"name is not recognised, so this would otherwise have been "
            f"measured as lsdtt3 doing less work than the reference.")


def polyfit_params(metrics, lambda_m, n_threads, extra=None):
    p = {
        "polyfit_method": "direct",
        "polyfit_lambda_metres": f"{lambda_m:g}",
        "n_threads": str(n_threads),
    }
    if metrics:
        # NO SPACE after the comma. lsdtt-surface-metrics does not trim
        # whitespace when it splits `polyfit_metrics`, so "a, b, c" silently
        # yields only `a` -- every name after the first arrives with a leading
        # space, matches nothing, and is dropped without a warning. The CLI
        # reference's own worked example uses ", " and would therefore compute
        # one metric instead of four. Reported separately; worked around here
        # so the benchmark cannot be flattered by lsdtt3 quietly doing less.
        p["polyfit_metrics"] = ",".join(metrics)
    if extra:
        p.update(extra)
    return p


# -----------------------------------------------------------------------------
# Running the reference
# -----------------------------------------------------------------------------

def run_reference(routine, fixture, out_path, repeats, extra_args=None):
    """Invoke topotoolbox_reference.py as a subprocess.

    A subprocess, not an import, so a missing optional package is a clean
    skip (77) rather than an ImportError inside the driver -- and so the
    reference's peak RSS is that of a process doing only the reference work.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(REFERENCE_SCRIPT), "--routine", routine,
           "--in", str(fixture), "--out", str(out_path),
           "--repeats", str(repeats)] + list(extra_args or [])
    proc = lsdtt_test_env.run(cmd, capture_output=True, text=True)

    if proc.returncode == SKIP_EXIT:
        die_skip("topotoolbox not available. Use `pixi run -e benchmark ...`")
    if proc.returncode != 0:
        raise RuntimeError(f"topotoolbox_reference.py {routine} failed:\n"
                           f"{proc.stderr[-2000:]}")
    return json.loads(proc.stdout.strip().splitlines()[-1])


# -----------------------------------------------------------------------------
# The correlation GUARD (not an accuracy check -- see module docstring)
# -----------------------------------------------------------------------------

def correlation_guard(lsdtt_path, ref_path, max_samples=200_000):
    """|Spearman rho| between the two outputs on commonly-valid cells.

    Returns (abs_rho, n_valid_lsdtt, n_valid_ref, lsdtt_dtype).

    Sign is discarded: lsdtt3 and TopoToolbox use opposite curvature sign
    conventions by design, and aspect handedness differs, so a rho of -0.99 is
    a pass. The only thing this can fail on is rho near zero, which would mean
    the stopwatch was timing two unrelated computations.

    Interior cells only. Both tools handle the boundary differently (lsdtt3
    masks it, TopoToolbox extends with mode='nearest'), and an edge disagreement
    says nothing about whether the interiors are the same quantity.
    """
    import numpy as np

    a, va = load_raster(lsdtt_path)
    b, vb = load_raster(ref_path)
    dtype = band_dtype(lsdtt_path)
    if a.shape != b.shape:
        return float("nan"), int(va.sum()), int(vb.sum()), dtype

    n_a, n_b = int(va.sum()), int(vb.sum())

    interior = np.zeros_like(va)
    pad = 2
    interior[pad:-pad, pad:-pad] = True
    both = va & vb & interior
    if int(both.sum()) < 100:
        return float("nan"), n_a, n_b, dtype

    x, y = a[both], b[both]
    if x.size > max_samples:
        idx = np.linspace(0, x.size - 1, max_samples).astype("int64")
        x, y = x[idx], y[idx]

    # Constant input makes the correlation undefined rather than zero; report
    # NaN so it cannot masquerade as a guard failure.
    if np.ptp(x) == 0.0 or np.ptp(y) == 0.0:
        return float("nan"), n_a, n_b, dtype

    try:
        from scipy.stats import spearmanr
        rho = float(spearmanr(x, y).statistic)
    except Exception:
        rx = np.argsort(np.argsort(x)).astype("float64")
        ry = np.argsort(np.argsort(y)).astype("float64")
        rho = float(np.corrcoef(rx, ry)[0, 1])

    return abs(rho), n_a, n_b, dtype


# -----------------------------------------------------------------------------
# Sub-studies
# -----------------------------------------------------------------------------

def measure_write_cost(exe, fixture, work_dir, lambda_m, repeats):
    """Seconds to write one extra output GeoTIFF, measured not guessed.

    Difference in whole-process wall time between asking for one metric and
    asking for two. Both runs perform exactly one polyfit and open the DEM
    once, so the fit, the start-up and the GDAL open cancel; what remains is
    one extra raster write plus the (near-free) algebra for a second metric.

    Needed because the comparisons the CLI does not time internally --
    currently only hillshade -- can otherwise report a filesystem measurement
    as though it were compute. On a WSL2/NTFS working directory this number is
    large enough to swamp everything else, which is precisely why it is
    measured on the same filesystem as the benchmark rather than assumed small.
    """
    d1 = work_dir / "writecost_1"
    d2 = work_dir / "writecost_2"
    for d in (d1, d2):
        if d.exists():
            shutil.rmtree(d)
    t1 = time_lsdtt(exe, fixture, d1,
                    polyfit_params(["slope"], lambda_m, 1), repeats)
    t2 = time_lsdtt(exe, fixture, d2,
                    polyfit_params(["slope", "slope_degrees"], lambda_m, 1),
                    repeats)
    return max(t2["wall_best"] - t1["wall_best"], 0.0)


def study_marginal_metric_cost(exe, fixture, work_dir, lambda_m, repeats, row0,
                               write_cost=0.0):
    """Marginal cost of an extra derived metric on the lsdtt3 side.

    The structural claim being tested: ONE polyfit yields all 27 metrics, so
    asking for 27 instead of 4 should cost only the per-metric algebra and the
    per-metric GeoTIFF write, NOT another fit. TopoToolbox has no equivalent --
    every extra curvature type is another five convolutions.

    Reported honestly: the wall-time delta includes 23 extra raster WRITES,
    which on a large fixture dominate the algebra. The fit-only time is
    reported alongside and should be unchanged between the two, which is the
    actual evidence for the claim.
    """
    rows, out = [], {}
    for label, metrics in (("marginal_metrics_4", FOUR_METRICS),
                           ("marginal_metrics_27", ALL_METRICS)):
        d = work_dir / label
        if d.exists():
            shutil.rmtree(d)
        t = time_lsdtt(exe, fixture, d,
                       polyfit_params(metrics, lambda_m, 1), repeats)
        assert_metrics_written(d, fixture.stem, metrics)
        out[label] = t
        rows.append({**row0, "comparison": label, "reference": "",
                     "routine": "", "timing_basis": "in_process_fit",
                     "n_metrics": len(metrics),
                     "lambda_metres": lambda_m,
                     "lsdtt_seconds_1t": t["fit_best"],
                     "lsdtt_median_1t": t["fit_median"],
                     "lsdtt_peak_rss_mb": t["peak_rss_mb"],
                     "note": f"wall_best={t['wall_best']:.4f}s "
                             f"(fit + {len(metrics)} raster writes)"})
    return rows, out


def study_window_scaling(exe, fixture, work_dir, cell_size, repeats, row0):
    """Is the `direct` fit really O(n), independent of window size?

    The CLI reference asserts a separable rolling-moment fit whose cost does
    not depend on the window. W = 40 is an 81x81 window; a naive per-pixel fit
    would be thousands of times slower there than at W = 1, so if the curve is
    anywhere near flat the claim holds, and if it is not, the documentation is
    wrong. Reported, never asserted.
    """
    rows = []
    for W in WINDOW_HALF_WIDTHS:
        lam = 2.0 * W * cell_size
        d = work_dir / f"window_W{W}"
        if d.exists():
            shutil.rmtree(d)
        t = time_lsdtt(exe, fixture, d,
                       polyfit_params(["slope"], lam, 1), repeats)
        rows.append({**row0, "comparison": f"window_W{W}", "reference": "",
                     "routine": "", "timing_basis": "in_process_fit",
                     "window_half_width": W, "lambda_metres": lam,
                     "cell_size": cell_size, "n_metrics": 1,
                     "lsdtt_seconds_1t": t["fit_best"],
                     "lsdtt_median_1t": t["fit_median"],
                     "lsdtt_peak_rss_mb": t["peak_rss_mb"],
                     "note": f"{2 * W + 1}x{2 * W + 1} window "
                             f"({(2 * W + 1) ** 2} cells)"})
    return rows


# -----------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", help="path to lsdtt-surface-metrics")
    ap.add_argument("--fetch-exe",
                    help="path to lsdtt-fetch-raster (fetches missing "
                         "fixtures; defaults to build/lsdtt-fetch-raster)")
    ap.add_argument("--fixture", default="alpine_medium",
                    help="fixture basename in tests/data/dem_cache "
                         "(default alpine_medium; 256x256 'small' fixtures are "
                         "too fast to time reliably)")
    ap.add_argument("--comparison", choices=list(COMPARISONS) + ["all"],
                    default="all")
    ap.add_argument("--repeats", type=int, default=5,
                    help="timed repetitions after a discarded warm-up (>=5)")
    ap.add_argument("--threads", type=int, default=0,
                    help="n_threads for the SECONDARY multi-threaded column; "
                         "0 = all hardware threads")
    ap.add_argument("--work-dir", default=None)
    ap.add_argument("--check", action="store_true",
                    help="fail if the correlation GUARD trips, i.e. if a "
                         "timing comparison turns out to be between unrelated "
                         "computations. Timings themselves are never asserted.")
    ap.add_argument("--skip-studies", action="store_true",
                    help="skip the marginal-cost and window-scaling studies")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    if args.list:
        for name, spec in COMPARISONS.items():
            flag = " [report-only]" if spec.get("report_only") else ""
            print(f"{name:22s}{flag} {spec['note']}")
        return 0

    try:
        import numpy  # noqa: F401
        from osgeo import gdal  # noqa: F401
    except ImportError:
        die_skip("numpy/GDAL not importable. Use `pixi run -e benchmark ...`")

    if not args.exe:
        die_skip("--exe not given (path to lsdtt-surface-metrics)")
    exe = pathlib.Path(args.exe).resolve()
    if not exe.exists():
        die_skip(f"executable {exe} not found; build lsdtt3 first")

    fetch_exe = pathlib.Path(args.fetch_exe).resolve() if args.fetch_exe else \
        (REPO_ROOT / "build" / "lsdtt-fetch-raster")
    fixture = flow.ensure_fixture(args.fixture, fetch_exe)

    work_dir = pathlib.Path(args.work_dir).resolve() if args.work_dir else \
        REPO_ROOT / "build" / "surface_benchmark"
    work_dir.mkdir(parents=True, exist_ok=True)

    fixture_sha = hashlib.sha256(fixture.read_bytes()).hexdigest()[:12]
    cell_size = cell_size_of(fixture)
    rows_shape = raster_shape(fixture)
    n_cells = rows_shape[0] * rows_shape[1]

    # Force W = 1 explicitly rather than relying on the default 70 m rounding
    # to 1 on 30 m data -- on any other cell size the default would give lsdtt3
    # a bigger window than the fixed 3x3 Evans kernels and invalidate every
    # timing in this file.
    lambda_m = 2.0 * cell_size
    W = window_half_width(lambda_m, cell_size)
    if W != 1:
        die_skip(f"could not force a 3x3 window: lambda={lambda_m} at "
                 f"cell_size={cell_size} gives W={W}. The window formula in "
                 f"lsdtt_raster_surface_metrics_cpu.cpp has probably changed; "
                 f"every Tier B timing depends on it.")

    row0 = {f: "" for f in CSV_FIELDS}
    row0.update({
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
        "git_sha": git_sha(),
        "fixture": args.fixture, "fixture_sha": fixture_sha,
        "n_cells": n_cells, "cell_size": cell_size,
        "window_half_width": W, "lambda_metres": lambda_m,
        "report_only": False,
    })

    print(f"fixture      : {args.fixture}  {rows_shape[1]}x{rows_shape[0]} "
          f"({n_cells} cells), cell {cell_size:g} m, sha {fixture_sha}")
    print(f"git          : {row0['git_sha']}")
    print(f"window       : lambda {lambda_m:g} m -> W = {W} "
          f"({2 * W + 1}x{2 * W + 1}), matching the fixed Evans 3x3 kernels")
    print(f"repeats      : {args.repeats} timed, 1 discarded warm-up, best-of-N")
    print("threading    : lsdtt3 measured at n_threads=1 (PRIMARY) and "
          f"n_threads={args.threads} (secondary). TopoToolbox's evansslope and")
    print("               curvature are single-threaded scipy convolutions, so "
          "only the n_threads=1 ratio is an algorithmic statement.")
    if args.repeats < 5:
        print("WARNING: --repeats below 5; N >= 5 is required for a stable median.")

    names = list(COMPARISONS) if args.comparison == "all" else [args.comparison]

    rows, failures = [], []

    # The baseline and the per-write cost each take 2 * (repeats + 1) CLI runs,
    # so they are measured LAZILY. Selecting one in-process-timed comparison
    # (which is what CI does) would otherwise pay ~24 process launches for two
    # numbers it never uses.
    _cache = {}

    def baseline_run():
        if "baseline" not in _cache:
            d = work_dir / "baseline"
            if d.exists():
                shutil.rmtree(d)
            b = time_lsdtt(exe, fixture, d, {}, args.repeats)
            print(f"\nbaseline (open DEM, no compute, no raster write): "
                  f"{b['wall_best']:.4f} s")
            _cache["baseline"] = b
        return _cache["baseline"]

    def write_cost_s():
        if "write_cost" not in _cache:
            c = measure_write_cost(exe, fixture, work_dir, lambda_m,
                                   args.repeats)
            print(f"cost of one extra output GeoTIFF on this filesystem: "
                  f"{c * 1e3:.1f} ms (measured on the same filesystem, because "
                  f"on WSL2/NTFS it")
            print("  dwarfs the compute and would otherwise be reported as one)")
            _cache["write_cost"] = c
        return _cache["write_cost"]

    for name in names:
        spec = COMPARISONS[name]
        print(f"\n=== {name} on {args.fixture} ===")
        print(f"    {spec['note']}")

        case_dir = work_dir / f"{args.fixture}_{name}"
        if case_dir.exists():
            shutil.rmtree(case_dir)

        basis = spec.get("timing_basis", "in_process_fit")
        report_only = bool(spec.get("report_only"))
        est_compute = float("nan")

        def measure(n_threads):
            d = case_dir / f"t{n_threads}"
            params = polyfit_params(spec["metrics"], lambda_m, n_threads,
                                    spec.get("params"))
            return time_lsdtt(exe, fixture, d, params, args.repeats)

        t1 = measure(1)
        tm = measure(args.threads)

        if basis in ("in_process_fit", "in_process"):
            l1, l1m = t1["fit_best"], t1["fit_median"]
            lm, lmm = tm["fit_best"], tm["fit_median"]
            if l1 is None or lm is None:
                raise RuntimeError(
                    f"{name}: the CLI printed no compute-timing line "
                    f"('Running direct fit... done' or 'Computed hillshade... "
                    f"done'), so there is no compute-only timing to report. "
                    f"The CLI's timing output has probably changed; fix "
                    f"FIT_TIMING_RE rather than falling back to wall time.")
        else:
            # Legacy path: wall minus a no-compute baseline. No comparison uses
            # it any more (hillshade moved to in_process on 2026-07-22), and it
            # cannot be used as-is: a run requesting no outputs is now a hard
            # error, because for a user it is always a mistake. Kept only so a
            # future wall-clock-only comparison has a documented starting point.
            raise RuntimeError(
                f"{name}: timing_basis 'wall_minus_baseline' needs a run with "
                f"no outputs enabled, which lsdtt-surface-metrics now rejects. "
                f"Add an in-process timer to the step and use 'in_process'.")
            baseline = baseline_run()
            l1 = max(t1["wall_best"] - baseline["wall_best"], 0.0)
            l1m = max(t1["wall_median"] - baseline["wall_median"], 0.0)
            lm = max(tm["wall_best"] - baseline["wall_best"], 0.0)
            lmm = max(tm["wall_median"] - baseline["wall_median"], 0.0)
            est_compute = max(l1 - write_cost_s(), 0.0)

        assert_metrics_written(case_dir / "t1", fixture.stem, spec["metrics"])
        lsdtt_out = case_dir / "t1" / f"{fixture.stem}{spec['suffix']}"
        if not lsdtt_out.exists():
            raise RuntimeError(f"expected lsdtt3 output {lsdtt_out} not written")

        # The cache key includes any extra reference arguments: without them a
        # re-run at a different hillshade azimuth would silently score the
        # correlation guard against the previous azimuth's raster.
        ref_key = "_".join([spec["routine"]] + list(spec.get("ref_args", [])))
        ref_key = re.sub(r"[^A-Za-z0-9_.-]", "", ref_key)
        ref_out = (work_dir / "topotoolbox_cache" /
                   f"{args.fixture}_{fixture_sha}_{ref_key}.tif")
        ref = run_reference(spec["routine"], fixture, ref_out, args.repeats,
                            spec.get("ref_args"))

        rho, n_a, n_b, l_dtype = correlation_guard(lsdtt_out, ref_out)

        # evansslope spends a large fraction of its time in
        # distance_transform_edt in-filling NaNs before it convolves anything.
        # That is real work, but it is not work lsdtt3 does, so quoting the raw
        # ratio would credit lsdtt3 for skipping a step rather than for being
        # faster at the fit. Discount it and quote the conservative number.
        ref_adj = ref["seconds"] - ref.get("seconds_nan_fill", 0.0)
        ref_adj = max(ref_adj, 1e-9)

        # If the discount ate most of the measurement, what is left is
        # difference-of-two-noisy-numbers and any ratio built on it is
        # meaningless -- it can read 0.00x purely from sampling. Demote the row
        # rather than publish a figure that looks like a 100x result.
        adj_fraction = ref_adj / ref["seconds"] if ref["seconds"] else 1.0
        discount_dominates = ("seconds_nan_fill" in ref) and adj_fraction < 0.25
        if discount_dominates:
            report_only = True

        sp1 = ref_adj / l1 if l1 else float("nan")
        spm = ref_adj / lm if lm else float("nan")

        mt_label = "all" if args.threads == 0 else str(args.threads)
        print(f"    lsdtt3, 1 thread    : {l1 * 1e3:8.3f} ms  (median "
              f"{l1m * 1e3:.3f} ms)   <-- PRIMARY")
        print(f"    lsdtt3, {mt_label} threads : {lm * 1e3:8.3f} ms  (median "
              f"{lmm * 1e3:.3f} ms)   secondary")
        print(f"    topotoolbox         : {ref['seconds'] * 1e3:8.3f} ms  (median "
              f"{ref['seconds_median'] * 1e3:.3f} ms, "
              f"{ref.get('n_convolutions')} convolutions, single-threaded)")
        if "seconds_nan_fill" in ref:
            print(f"      of which NaN in-fill (distance_transform_edt, work "
                  f"lsdtt3 does not do): {ref['seconds_nan_fill'] * 1e3:.3f} ms")
            print(f"      -> {ref_adj * 1e3:.3f} ms of comparable work; "
                  f"ratios below use this, not the raw number")
            if discount_dominates:
                print(f"      DEMOTED TO REPORT-ONLY: the NaN in-fill is "
                      f"{100 * (1 - adj_fraction):.0f}% of the reference "
                      f"timing, so the remainder is a difference of two noisy")
                print("      numbers and no ratio can be built on it. Use a "
                      "larger fixture if this row matters.")
        if basis == "wall_minus_baseline":
            print(f"      lsdtt3 figure above is wall-minus-baseline and still "
                  f"contains one GeoTIFF write ({write_cost_s() * 1e3:.1f} ms)")
            if est_compute > 0.0:
                print(f"      -> indicative lsdtt3 COMPUTE, 1 thread: "
                      f"{est_compute * 1e3:.1f} ms (ESTIMATE, not a measured "
                      f"compute-only timing)")
            else:
                # The write alone costs more than the whole measured delta, so
                # the compute is somewhere inside the write's own noise. Saying
                # "0.0 ms" would be a fabricated number; say what is actually
                # known instead.
                print("      -> lsdtt3 COMPUTE is smaller than the run-to-run "
                      "spread of the write itself, so it cannot be separated "
                      "here at all.")
                print("         An in-process timer in the CLI is the only way "
                      "to get this number.")
        if report_only:
            print("    ratio               : NOT QUOTED -- report-only; see the "
                  "note above for why the two sides are not like-for-like")
        else:
            print(f"    speedup, 1 thread   : {sp1:6.2f}x   <-- the honest number")
            print(f"    speedup, {mt_label} threads : {spm:6.2f}x   "
                  f"(threaded C++ vs serial scipy; throughput, not algorithm)")
        print(f"    valid cells         : lsdtt3 {n_a}  topotoolbox {n_b}")
        print(f"    band dtype          : lsdtt3 {l_dtype}  topotoolbox "
              f"{ref.get('out_dtype')}  (input DEM {ref.get('dem_dtype')})")
        print(f"    peak RSS            : lsdtt3 {t1['peak_rss_mb']:.1f} MiB  "
              f"topotoolbox {ref.get('peak_rss_mb', float('nan')):.1f} MiB")
        print(f"    |rank correlation|  : {rho:.4f}  (GUARD only, not accuracy)")

        rows.append({**row0, "comparison": name, "reference": "topotoolbox",
                     "routine": spec["routine"], "report_only": report_only,
                     "timing_basis": basis, "n_metrics": len(spec["metrics"]),
                     "lsdtt_seconds_1t": l1, "lsdtt_median_1t": l1m,
                     "lsdtt_seconds_mt": lm, "lsdtt_median_mt": lmm,
                     "lsdtt_threads_mt": args.threads,
                     "ref_seconds": ref["seconds"],
                     "ref_median": ref["seconds_median"],
                     "ref_seconds_nan_fill": ref.get("seconds_nan_fill", ""),
                     "ref_seconds_comparable": ref_adj,
                     "ref_n_convolutions": ref.get("n_convolutions"),
                     "speedup_1t": "" if report_only else sp1,
                     "speedup_mt": "" if report_only else spm,
                     "write_cost_seconds": (write_cost_s()
                                            if basis != "in_process_fit" else ""),
                     "lsdtt_est_compute_1t": (est_compute
                                              if basis != "in_process_fit" else ""),
                     "abs_correlation": rho,
                     "n_valid_lsdtt": n_a, "n_valid_ref": n_b,
                     "lsdtt_dtype": l_dtype, "ref_dtype": ref.get("out_dtype"),
                     "lsdtt_peak_rss_mb": t1["peak_rss_mb"],
                     "ref_peak_rss_mb": ref.get("peak_rss_mb"),
                     "note": spec["note"]})

        # The only thing that can fail: a timing comparison between two tools
        # that turn out not to be computing the same kind of quantity.
        if args.check and not report_only and rho == rho and rho < CORRELATION_GUARD:
            failures.append(
                f"{name}: |rank correlation| = {rho:.4f} below the guard "
                f"{CORRELATION_GUARD}. The two tools are not computing the same "
                f"kind of quantity, so this timing comparison is meaningless "
                f"and must not be published.")

        # A tool that produced far fewer valid cells did less work.
        if n_a and n_b and min(n_a, n_b) / max(n_a, n_b) < 0.9:
            print(f"    NOTE: valid-cell counts differ by more than 10% "
                  f"({n_a} vs {n_b}); the two tools did not compute over the "
                  f"same domain and this timing is not like-for-like.")

    if not args.skip_studies:
        print("\n=== marginal cost per extra derived metric (lsdtt3 only) ===")
        mrows, m = study_marginal_metric_cost(exe, fixture, work_dir, lambda_m,
                                              args.repeats, row0,
                                              write_cost_s())
        rows.extend(mrows)
        a, b = m["marginal_metrics_4"], m["marginal_metrics_27"]
        n_extra = len(ALL_METRICS) - len(FOUR_METRICS)
        d_wall = b["wall_best"] - a["wall_best"]
        print(f"    4 metrics : fit {a['fit_best'] * 1e3:.3f} ms, "
              f"whole run {a['wall_best']:.4f} s")
        print(f"    27 metrics: fit {b['fit_best'] * 1e3:.3f} ms, "
              f"whole run {b['wall_best']:.4f} s")
        print(f"    the FIT -- the expensive part -- is paid ONCE either way: "
              f"{a['fit_best'] * 1e3:.3f} ms vs {b['fit_best'] * 1e3:.3f} ms")
        per_metric = d_wall / n_extra
        wc = write_cost_s()
        print(f"    marginal cost per extra metric: {per_metric * 1e3:.2f} ms, "
              f"of which {wc * 1e3:.2f} ms is the GeoTIFF write")
        algebra = per_metric - wc
        algebra_txt = (f"{algebra * 1e3:.2f} ms" if algebra > 0.0 else
                       "below the resolution of this measurement (the write "
                       "alone accounts for the whole marginal cost)")
        print(f"    -> ALGEBRA-ONLY marginal cost per extra metric: "
              f"{algebra_txt}")
        print(f"       Compare {a['fit_best'] * 1e3:.1f} ms for the fit itself. "
              f"Adding a metric is essentially free;")
        print("       what you pay for is writing the raster. TopoToolbox has "
              "no equivalent: each extra curvature type")
        print("       costs another five convolutions there -- see the "
              "curvature_all4 row above for that cost measured.")

        print("\n=== window-size scaling: is `direct` really O(n)? ===")
        wrows = study_window_scaling(exe, fixture, work_dir, cell_size,
                                     args.repeats, row0)
        rows.extend(wrows)
        secs = [r["lsdtt_seconds_1t"] for r in wrows]
        for W, s in zip(WINDOW_HALF_WIDTHS, secs):
            print(f"    W = {W:2d}  ({2 * W + 1:3d}x{2 * W + 1:<3d}, "
                  f"{(2 * W + 1) ** 2:5d} cells): {s * 1e3:8.3f} ms")
        if all(s > 0 for s in secs):
            spread = max(secs) / min(secs)
            naive = ((2 * WINDOW_HALF_WIDTHS[-1] + 1) ** 2 /
                     (2 * WINDOW_HALF_WIDTHS[0] + 1) ** 2)
            print(f"    max/min = {spread:.2f}x across a {naive:.0f}x increase "
                  f"in window cells")
            print(f"    -> the O(n) claim in the CLI reference "
                  f"{'HOLDS' if spread < 3.0 else 'DOES NOT HOLD'} "
                  f"(flat curve threshold 3x, reported not asserted)")

    csv_path = work_dir / "surface_benchmark.csv"
    write_header = not csv_path.exists()
    with open(csv_path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if write_header:
            w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {csv_path}")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("\n(timings are REPORTED, never asserted -- no wall-clock threshold "
          "is enforced in CI)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
