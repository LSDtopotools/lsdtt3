#!/usr/bin/env python3
"""Flow-routing benchmark driver: lsdtt3 vs RichDEM.

docs/developer_notes/testing_and_benchmarking.md.

Runs an lsdtt3 CLI and the corresponding RichDEM reference over the same fetched
fixture, computes the comparison metrics of docs/developer_notes/testing_and_benchmarking.md, and writes one row per
comparison to benchmark_results.csv.

Phase 2 scope is the harness: metrics are computed and REPORTED.  Tolerances are
deliberately not asserted yet -- docs/developer_notes/testing_and_benchmarking.md requires pinning them just above
observed divergence once real numbers exist, otherwise a tolerance asserts
nothing.  Pass --check to enforce whatever tolerances have been pinned in
TOLERANCES so far.

Skips cleanly (exit 77, CTest's skip code) when RichDEM or the fixtures are
unavailable, so an optional dependency never turns CI red.

Usage:
    pixi run -e benchmark python tests/scripts/run_flow_benchmark.py \
        --exe build/lsdtt-raster-preprocessing
    python tests/scripts/run_flow_benchmark.py --list
"""

import argparse
import csv
import hashlib
import json
import pathlib
import shutil
import subprocess
import sys
import time

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))  # python -I does not add it
import lsdtt_test_env  # noqa: E402

SKIP_EXIT = 77

# Comparisons that are recorded but NOT asserted, keyed "comparison@reference".
# RichDEM's accumulation is here because its D8 rule differs from lsdtt3's by
# design (see above), so pinning a tolerance against it would enshrine that
# difference as if it were lsdtt3 error.
REPORT_ONLY = {
    "accum_d8@richdem", "accum_dinf@richdem", "accum_mfd@richdem",
}

# Reference implementations. Each is a standalone script invoked as a
# subprocess so that a missing optional package is a clean skip rather than an
# import failure in the driver.
#
# richdem     -- valid for depression handling; NOT valid for D8 flow direction.
#                It selects the LOWEST-ELEVATION neighbour rather than the
#                steepest descent (drop / distance), so it takes a diagonal
#                wherever the diagonal is lower but a cardinal is steeper.
#                (Corrected 2026-08-21: its tie-break prefers CARDINALS, and it
#                passes a planar-ramp oracle -- the divergence is the selection
#                rule, not ties. See docs/developer_notes/testing_and_benchmarking.md 5.1.)
# topotoolbox -- true steepest descent, same rule as lsdtt3.
# whitebox    -- ALSO true steepest descent. MEASURED, not assumed:
#                tests/scripts/check_whitebox_d8_rule.py builds a cell sitting
#                between a LOWER diagonal (10 m drop over 42.43 m) and a STEEPER
#                cardinal (8 m over 30 m), so the two rules MUST disagree there,
#                and scores the whole accumulation raster against a reference
#                implementation of each rule. Whitebox matched steepest descent
#                on every cell; a d8_pointer read, with the direction codes
#                calibrated on a planar ramp rather than assumed, confirmed it
#                independently.
#
#                A planar ramp on its own would NOT have established this, which
#                is why the discriminating fixture exists: on a ramp the S/SE/SW
#                neighbours are TIED in elevation, so a lowest-elevation rule
#                that breaks ties toward cardinals passes it. RichDEM does
#                exactly that.
#
# D8 accumulation is therefore asserted against TWO independent implementations,
# and RichDEM's exclusion is a statement about its selection rule rather than a
# hole in the comparison.
#
# See docs/developer_notes/testing_and_benchmarking.md.
REFERENCE_SCRIPTS = {
    "richdem":     "richdem_reference.py",
    "topotoolbox": "topotoolbox_reference.py",
    "whitebox":    "whitebox_reference.py",
}

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "dem_cache"

# Fixture registry and fetch: flow_benchmark_fixtures.py (stdlib only, so a
# stdlib-tier test can share it). FIXTURE_DIR above is where this driver caches
# them; main() rebinds it for --fixture-dir, and ensure_fixture() reads it at
# call time.
from flow_benchmark_fixtures import (  # noqa: E402
    FIXTURES, GRID_SPACING_M, SIZES, fixture_bbox)
import flow_benchmark_fixtures  # noqa: E402


def ensure_fixture(fixture_key: str, fetch_exe: pathlib.Path) -> pathlib.Path:
    """Fetch a fixture with lsdtt-fetch-raster if it is not already cached."""
    return flow_benchmark_fixtures.ensure_fixture(fixture_key, fetch_exe, FIXTURE_DIR)

# float32(0.001) as a double: the flat increment the lsdtt3 CLI actually uses
# when given 0.001 (see breach_least_cost below). Computed with the struct
# module so the comparison table does not need numpy at import time.
_F32_0P001 = __import__("struct").unpack("f", __import__("struct").pack("f", 0.001))[0]

# Comparison definitions.  `richdem_routine` names a routine in
# richdem_reference.py; `lsdtt_suffix` is the output file lsdtt3 writes.
#
# NOTE on fill: lsdtt3's fill_min_slope defaults to 1e-4, which
# imposes an epsilon gradient across flats, while RichDEM's epsilon=False gives
# true flats.  These SHOULD disagree.  The tight comparison is min_slope=0 vs
# epsilon=False; the epsilon pair is a loose shape-agreement check only.
#
# The two epsilon variants are NOT the same operation: RichDEM
# raises each cell by one ULP (std::nextafter), which scales with elevation
# magnitude and dtype, while lsdtt3 applies a constant slope.  They happen to
# nearly coincide at Alpine elevations in float32 (ULP = 1.22e-4 vs 1e-4) and
# diverge by ~40x in volume moved at Lakeland elevations (ULP = 7.6e-6).  That
# is why the epsilon tolerances below differ so wildly between fixtures.
COMPARISONS = {
    "fill": {
        # Whitebox enters this comparison through `depression_fill_wl`
        # (Wang & Liu, fix_flats=False), NOT through its Barnes tool. On
        # alpine_small the Barnes tool returns a global minimum 6.45 m above the
        # border minimum where wang_and_liu returns exactly the border minimum,
        # and the two Whitebox fills differ by up to 9.78 m. lsdtt3 is bit-exact
        # against both RichDEM and TopoToolbox on this same operation, so the
        # odd one out is that tool -- most likely an edge-cell boundary
        # convention. Until that is understood it stays out of the value
        # comparison; `fill@whitebox` remains a TIMING row in the campaign.
        "references": {"richdem": "fill", "topotoolbox": "fill",
                       "whitebox": "depression_fill_wl"},
        "lsdtt_suffix": "_fill.tif",
        "params": {"write_fill_raster": "true", "fill_min_slope": "0"},
        "kind": "elevation",
        "note": "Barnes 2014, no epsilon -- the tight comparison",
    },
    "fill_epsilon": {
        # WHITEBOX IS THE BETTER REFERENCE HERE, and the reason is the whole
        # point of splitting this case out. lsdtt3 imposes a CONSTANT
        # fill_min_slope across flats. RichDEM's epsilon=True raises each cell
        # by one ULP of std::nextafter, which scales with elevation magnitude --
        # a different operation, which is why the pinned richdem tolerances
        # differ ~40x between the two fixtures (volume ratio 1.005 vs 45.0).
        # Whitebox's fill_depressions(fix_flats=True, flat_increment=e) is a
        # constant increment, i.e. the same operation lsdtt3 performs, and
        # whitebox_reference.py defaults --flat-increment to 1e-4 to match
        # fill_min_slope below.
        "references": {"richdem": "fill_epsilon",
                       "whitebox": "depression_fill_wl_epsilon"},
        "lsdtt_suffix": "_fill.tif",
        "params": {"write_fill_raster": "true", "fill_min_slope": "0.0001"},
        "kind": "elevation",
        "note": ("lsdtt3 steps by a constant min_slope, RichDEM by one ULP "
                 " -- shape agreement only"),
    },
    "breach": {
        "references": {"richdem": "breach"},
        "lsdtt_suffix": "_breached.tif",
        "params": {"write_breached_raster": "true",
                   "breach_algorithm_choice": "lindsay2016"},
        "kind": "elevation",
        "note": "both implementations pit-flatten before carving",
    },
    "breach_least_cost": {
        # Least-cost breaching, Lindsay (2020): lsdtt3's port of Whitebox
        # BreachDepressionsLeastCost against whitebox_workflows 2.0.6 itself.
        #
        # ASSERTED BITWISE: float32(Whitebox's f64 output) == lsdtt3's float32
        # output on every valid cell, and the nodata masks are identical.
        # Whitebox writes Float64, so the comparison rounds it once to float32
        # (`metric: float32_exact`) -- the same single rounding lsdtt3 does at
        # the end of its double-precision core.
        #
        # Matched parameters, and why each one:
        #   breach_lc_edge_policy whitebox  lsdtt3's default (outlet) treats
        #       the raster edge as an outlet; Whitebox does not.
        #   fill_remainder_after_breach false / fill_deps False  the two
        #       remainder fills are different algorithms.
        #   max_dist 100, minimize_dist false  Whitebox's defaults.
        #   flat increment float32(0.001)  the lsdtt3 CLI reads the parameter
        #       as float32, so Whitebox is handed that exact value
        #       (0.0010000000474974513), not the double 0.001; the channel
        #       targets pit_z - L * increment would otherwise differ. Explicit
        #       on both sides: Whitebox's own auto increment (~1e-10 m) is
        #       below float32 resolution and lsdtt3's auto rule is different
        #       by design.
        "references": {"whitebox": "breach"},
        "lsdtt_suffix": "_breached.tif",
        "params": {"write_breached_raster": "true",
                   "breach_algorithm_choice": "least_cost",
                   "breach_lc_max_dist": "100",
                   "breach_lc_minimize_dist": "false",
                   "breach_lc_flat_increment": "0.001",
                   "breach_lc_edge_policy": "whitebox",
                   "fill_remainder_after_breach": "false"},
        "reference_args": {"whitebox": [
            "--breach-max-dist", "100",
            "--breach-flat-increment", repr(float(_F32_0P001))]},
        "kind": "elevation",
        "metric": "float32_exact",
        "note": "Lindsay 2020 least-cost breach vs Whitebox 2.0.6 -- bitwise",
    },
    "flats_float_step": {
        # Fill, then Barnes, Lehman and Mulla (2014b) flat resolution with
        # float32 nextafter steps, against RichDEM FillDepressions
        # (epsilon=False) + ResolveFlats on the same float32 raster. The
        # synthetic oracle check_flat_resolution_vs_richdem.py proves the two
        # are equal BIT FOR BIT; this extends that to real terrain, whose flats
        # (lake surfaces, quantised plateaus, filled depressions) are messier
        # than anything generated there.
        #
        # ASSERTED BITWISE, zero tolerance, on ROUND 1. lsdtt3's first round
        # is RichDEM's single pass; where that pass leaves cells without a
        # lower neighbour (a D8 router refuses them), lsdtt3 refills and
        # resolves again and, if needed, runs a final strict-rise pass, and
        # prints "Changed after round 1: <n> cell(s)". Those cells are the
        # only permitted differences: `n_unexplained_differences` (cells whose
        # float32 bits differ, minus n) is pinned at 0, and so are
        # `n_below_richdem` (cells where lsdtt3 is LOWER than RichDEM, by
        # ordered float32 key) and `n_stranded_lsdtt` (cells left without a
        # lower neighbour). `n_differing_unguarded` and `n_stranded_richdem`
        # are reported. The one other exemption is lsdtt3's near-zero guard:
        # on a flat with |z| < 1 m a float32 step is subnormal and rounds away
        # in slope arithmetic, so lsdtt3 steps by 2^-23 m per mask unit there
        # instead, by design. Those cells are counted separately (`n_guarded`)
        # and not asserted. Neither fixture has a cell below 1 m, so on these
        # the exemption is empty.
        #
        # Flats are found on the true-flat fill (`fill` params above, itself
        # pinned bit-exact against RichDEM), run as a second lsdtt3 stage.
        #
        # ELEVATIONS, NOT D8 DIRECTIONS: RichDEM's D8 picks the lowest
        # neighbour rather than the steepest descent (see accum_d8 below), so
        # comparing directions would score that rule difference as lsdtt3
        # error. Identical elevations mean any steepest-descent router sees the
        # same surface.
        #
        # min_elevation stays at its default (0, which masks cells at or below
        # 0 m): both fixtures' minima are well above it (alpine 1177 m,
        # lakeland 74.9 m), and the always-enforced mask check would catch any
        # cell one side dropped.
        "references": {"richdem": "fill_resolve_flats"},
        "lsdtt_suffix": "_fill.tif",
        "params": {"write_fill_raster": "true", "fill_min_slope": "0",
                   "fill_flat_resolution_choice": "barnes2014b_float_step"},
        "flat_baseline_params": {"write_fill_raster": "true",
                                 "fill_min_slope": "0"},
        "kind": "elevation",
        "metric": "float32_bits",
        "note": "Barnes 2014b float32 steps vs RichDEM ResolveFlats -- bitwise",
    },

    # -------------------------------------------------------------------------
    # Flow accumulation (docs/developer_notes/testing_and_benchmarking.md, Phase 4)
    #
    # Run on the SAME filled DEM that lsdtt3 produced, with raster_is_filled so
    # the drainage CLI does not fill again. That isolates the accumulation
    # algorithm: any difference here cannot be blamed on the fill, which the
    # `fill` comparison has already shown to be bit-exact.
    #
    # THREE references, and they are NOT interchangeable.
    #
    #   topotoolbox  ASSERTED. Steepest descent, same rule as lsdtt3.
    #   whitebox     ASSERTED. Steepest descent, MEASURED not assumed -- see
    #                REFERENCE_SCRIPTS above and check_whitebox_d8_rule.py.
    #   richdem      REPORT-ONLY. Different selection rule; see below.
    #
    # *** RICHDEM IS REPORT-ONLY, NOT ASSERTED ***
    # RichDEM is not a valid reference for D8 flow direction: its D8 picks the
    # neighbour with the LOWEST ELEVATION (richdem/flowmet/d8_flowdirs.hpp:63-71)
    # rather than the steepest DESCENT (drop / distance), so it takes a diagonal
    # wherever the diagonal is lower but a cardinal is steeper. Pinning
    # tolerances against it would enshrine that as if it were lsdtt3 error.
    #
    # Note it is NOT ties that separate them -- RichDEM's tie-break prefers
    # cardinals and it passes a planar-ramp oracle. The divergence needs a case
    # where the lowest neighbour is diagonal AND a cardinal is steeper.
    #
    # The reference-free assertions live in tests/cpp/test_flow_properties.cpp:
    # P4 mass conservation, P5 monotonicity, P8 exact planar-ramp oracle. Those
    # remain the primary evidence; agreement with two independent
    # implementations is corroboration on real terrain, not a substitute.
    # -------------------------------------------------------------------------
    "accum_d8": {
        "references": {"richdem": "accum_d8", "topotoolbox": "accum_d8",
                       "whitebox": "accum_d8"},
        "lsdtt_suffix": "_d8_drainage_area.tif",
        "params": {"write_d8_drainage_area_raster": "true",
                   "raster_is_filled": "true"},
        "kind": "accumulation",
        "stage": "drainage",
        "note": "D8, O'Callaghan & Mark (1984)",
    },
    "accum_dinf": {
        "references": {"richdem": "accum_dinf"},
        "lsdtt_suffix": "_dinf_drainage_area.tif",
        "params": {"write_dinf_drainage_area_raster": "true",
                   "raster_is_filled": "true"},
        "kind": "accumulation",
        "stage": "drainage",
        "note": "D-infinity, Tarboton (1997) -- proportioning rules diverge more",
    },
    "accum_mfd": {
        "references": {"richdem": "accum_mfd"},
        "lsdtt_suffix": "_mfd_drainage_area.tif",
        "params": {"write_mfd_drainage_area_raster": "true",
                   "raster_is_filled": "true", "mfd_exponent": "1.0"},
        "kind": "accumulation",
        "stage": "drainage",
        "note": "MFD, Quinn (1991) vs RichDEM's Freeman -- loosest comparison",
    },
}

# -----------------------------------------------------------------------------
# Pinned tolerances (docs/developer_notes/testing_and_benchmarking.md, Phase 3)
#
# Pinned just above OBSERVED divergence, per docs/developer_notes/testing_and_benchmarking.md: a tolerance set to a
# round number far above reality asserts nothing.
#
# Keyed by FIXTURE as well as comparison, because breach divergence turned out
# to be fixture-dependent -- alpine breach measured 4.54 m on one fixture build
# and 0.00 m on another. A single global tolerance would either be too loose for
# alpine or too tight for lakeland.
#
# Each entry:
#   max_abs          maximum |lsdtt3 - richdem| over valid cells, metres
#   frac_differing   fraction of valid cells differing by more than 1 cm
#   volume_ratio_max upper bound on max(r, 1/r) where r is the ratio of material
#                    each implementation moved. 1.0 means they moved exactly the
#                    same volume; 40 means one moved 40x the other.
#
# The exact-zero pins are deliberate, not laziness. lsdtt3's Barnes 2014 fill is
# bit-exact against RichDEM on both terrain regimes, and both routines are
# serial, so the result is deterministic. Pinning 0.0 makes any future
# divergence -- however small -- a failure, which is exactly what we want from a
# port that is supposed to reproduce a reference implementation. It is only safe
# because EXPECTED_FIXTURE_SHA below guarantees we are scoring the same raster.
#
# The volume bounds tell a story the cell-wise metrics do not: in LOW relief
# lsdtt3's epsilon fill moves ~40x more material than RichDEM's, and its breach
# carves ~3.7x more, while in HIGH relief the same comparisons are 1.002x and
# 1.000x. The epsilon schemes are genuinely different operations (lsdtt3 imposes
# min_slope across every flat; RichDEM adds a minimal float increment) and that
# only becomes visible where flats dominate. These bounds are pinned to hold the
# CURRENT behaviour, not to bless it -- see docs/developer_notes/testing_and_benchmarking.md
# -----------------------------------------------------------------------------
TOLERANCES = {
    "alpine_small": {
        # Barnes 2014 vs RichDEM epsilon=False: bit-exact.
        "fill@richdem":         {"max_abs": 0.0,  "frac_differing": 0.0,
                                 "volume_ratio_max": 1.0},
        # ... and bit-exact against TopoToolbox's fillsinks too, which is a
        # genuinely independent implementation rather than a sibling port.
        "fill@topotoolbox":     {"max_abs": 0.0,  "frac_differing": 0.0,
                                 "volume_ratio_max": 1.0},
        # ... and bit-exact against Whitebox's wang_and_liu as well. THREE
        # independent implementations reproducing lsdtt3's fill cell-for-cell.
        # Measured 2026-09-08.
        "fill@whitebox":        {"max_abs": 0.0,  "frac_differing": 0.0,
                                 "volume_ratio_max": 1.0},
        "fill_epsilon@richdem": {"max_abs": 0.12, "frac_differing": 0.006,
                                 "volume_ratio_max": 1.005},
        # Whitebox, with flat_increment matched to fill_min_slope * cell_size,
        # is 4x closer than RichDEM on every metric (0.0246 m vs 0.1052 m;
        # 0.105 % vs 0.523 % of cells; volume 1.00037 vs 1.00230). Measured
        # 2026-09-08. Pinned just above observed.
        "fill_epsilon@whitebox": {"max_abs": 0.026, "frac_differing": 0.0012,
                                  "volume_ratio_max": 1.001},
        "breach@richdem":       {"max_abs": 0.0,  "frac_differing": 0.0,
                                 "volume_ratio_max": 1.0},
        # Float32 bits identical on every valid cell outside the near-zero
        # guard and the cells lsdtt3 changed after round 1 (see COMPARISONS).
        # Zero is the claim, not a measurement margin. Measured 2026-10-10:
        # RichDEM leaves 1 cell without a lower neighbour; lsdtt3 takes 2
        # rounds and changes 1 cell after round 1, so 1 cell differs (it
        # was `n_differing_unguarded: 0` before the round loop).
        "flats_float_step@richdem": {"n_unexplained_differences": 0,
                                     "n_below_richdem": 0,
                                     "n_stranded_lsdtt": 0},
        # Least-cost breach vs Whitebox 2.0.6: float32 bits identical on every
        # valid cell. Zero is the claim (a port), not a measurement margin.
        "breach_least_cost@whitebox": {"n_differing": 0, "max_abs": 0.0},
        # D8 accumulation is IDENTICAL to TopoToolbox on this fixture -- same
        # steepest-descent rule, and steep terrain leaves almost no ties to
        # break. Pinned tight so any divergence fails.
        "accum_d8@topotoolbox": {"rank_correlation_min": 0.9999,
                                 "top_channels_rel": 1e-6, "total_rel": 1e-6},
        # Whitebox is EXACT here too, on all three metrics: rank correlation
        # 1.0000000000, top-channel and total agreement 0.0. Measured
        # 2026-09-08. Two independent implementations reproducing lsdtt3's D8
        # cell-for-cell on high relief is the strongest statement this
        # comparison can make, so it is pinned as tight as topotoolbox's.
        "accum_d8@whitebox":    {"rank_correlation_min": 0.9999,
                                 "top_channels_rel": 1e-6, "total_rel": 1e-6},
    },
    "lakeland_small": {
        "fill@richdem":         {"max_abs": 0.0,  "frac_differing": 0.0,
                                 "volume_ratio_max": 1.0},
        "fill@topotoolbox":     {"max_abs": 0.0,  "frac_differing": 0.0,
                                 "volume_ratio_max": 1.0},
        # Bit-exact on the low-relief fixture too, where the flats that make
        # every other comparison diverge are at their worst. Measured 2026-09-08.
        "fill@whitebox":        {"max_abs": 0.0,  "frac_differing": 0.0,
                                 "volume_ratio_max": 1.0},
        "fill_epsilon@richdem": {"max_abs": 0.55, "frac_differing": 0.92,
                                 "volume_ratio_max": 45.0},
        # THE ROW THAT JUSTIFIES SPLITTING THIS CASE. On flat-dominated terrain
        # RichDEM's one-ULP scheme moves 40.5x the material lsdtt3 does; with a
        # matched constant increment Whitebox moves 1.39x. Same fixture, same
        # lsdtt3 output -- the 40x was the epsilon SCHEME, not lsdtt3.
        #
        # 1.39 is still not 1.0, and should not be forced to be: the two impose
        # the gradient differently along a flat (lsdtt3 scales by distance, so a
        # diagonal step rises by sqrt(2) more). That is a real operational
        # difference, so the pin sits above what was measured rather than at a
        # number that would assert the two schemes are identical.
        # Measured 2026-09-08.
        "fill_epsilon@whitebox": {"max_abs": 0.15, "frac_differing": 0.82,
                                  "volume_ratio_max": 1.40},
        "breach@richdem":       {"max_abs": 0.95, "frac_differing": 0.032,
                                 "volume_ratio_max": 4.0},
        # As alpine_small. Measured 2026-10-10: RichDEM leaves 41 cells without
        # a lower neighbour; lsdtt3 takes 2 rounds and changes 51 cells after
        # round 1, which are exactly the 51 cells that differ.
        "flats_float_step@richdem": {"n_unexplained_differences": 0,
                                     "n_below_richdem": 0,
                                     "n_stranded_lsdtt": 0},
        # As alpine_small: bitwise.
        "breach_least_cost@whitebox": {"n_differing": 0, "max_abs": 0.0},
        # Low relief does NOT reproduce exactly, and that is expected rather
        # than alarming: after epsilon filling, flat ground still leaves many
        # near-ties, and D8 must pick one neighbour. Tie-breaking moves WHICH
        # cell carries a value without changing the network, which is precisely
        # what rank correlation tolerates and a cell-wise diff would not. The
        # binding constraint here is the rank correlation, not the channel
        # agreement.
        "accum_d8@topotoolbox": {"rank_correlation_min": 0.93,
                                 "top_channels_rel": 0.17, "total_rel": 0.002},
        # Whitebox agrees with lsdtt3 MORE CLOSELY than TopoToolbox does on low
        # relief, which is worth stating because it is the opposite of what the
        # alpine row would lead you to expect. Measured 2026-09-08:
        #
        #                  rank corr    top 1% channels    total
        #   whitebox       0.99968          0.60 %        0.0051 %
        #   topotoolbox    0.94018         15.11 %        0.1466 %
        #
        # Same D8 rule in all three, so the spread is tie-breaking among the
        # near-ties that survive an epsilon fill on flat ground -- and Whitebox
        # happens to break them the way lsdtt3 does. That is a coincidence of
        # implementation, not a correctness result, so the pins are set just
        # above what was measured rather than at the round numbers the alpine
        # row can afford.
        "accum_d8@whitebox":    {"rank_correlation_min": 0.999,
                                 "top_channels_rel": 0.0065,
                                 "total_rel": 0.00006},
    },
}

# SHA256 (first 12 hex chars) of the fixture each tolerance set was pinned
# against. Fixture content depends on the GDAL build that fetched it, and a
# regenerated fixture is a different raster -- applying pinned tolerances to it
# would be meaningless. On mismatch the run SKIPS with an explanation rather
# than failing, because a different GDAL is an environment difference, not a
# defect in lsdtt3.
EXPECTED_FIXTURE_SHA = {
    "alpine_small":   "899cc313e27c",
    "lakeland_small": "2af105584b68",
}


def die_skip(msg: str):
    print(f"SKIP: {msg}", file=sys.stderr)
    sys.exit(SKIP_EXIT)


def load_raster(path: pathlib.Path):
    """Return (array, nodata_mask) as float64 / bool."""
    import numpy as np
    from osgeo import gdal
    gdal.UseExceptions()

    ds = gdal.Open(str(path))
    if ds is None:
        raise RuntimeError(f"could not open {path}")
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype("float64")
    nodata = band.GetNoDataValue()
    valid = np.isfinite(arr)
    if nodata is not None:
        valid &= (arr != nodata)
    return arr, valid


def elevation_metrics(a_path, b_path, original_path=None, tol=0.01):
    """Comparison metrics for two elevation rasters.

    `original_path` is the unmodified input DEM. It is needed for the volume
    check: docs/developer_notes/testing_and_benchmarking.md asks for a "total fill volume ratio", i.e. how much
    material each implementation moved. Comparing sums of ELEVATIONS instead
    would be datum-dependent -- adding 1000 m to every cell would change the
    ratio without changing anything about the algorithms -- and in low relief it
    is dominated by the datum rather than by the fill. Comparing sums of
    |change from the original| is datum-independent and is what the metric was
    meant to be.
    """
    import numpy as np

    a, va = load_raster(a_path)
    b, vb = load_raster(b_path)

    if a.shape != b.shape:
        raise RuntimeError(f"shape mismatch {a.shape} vs {b.shape}")

    # The nodata masks must match exactly.  A mismatch is a bug, never
    # something to absorb into a tolerance.
    mask_identical = bool((va == vb).all())
    both = va & vb
    n = int(both.sum())
    if n == 0:
        raise RuntimeError("no overlapping valid cells")

    d = a[both] - b[both]
    ad = np.abs(d)

    # Volume of material moved by each implementation, relative to the original
    # DEM. Integral check: catches a systematic offset that a cell-wise maximum
    # would miss.
    volume_ratio = float("nan")
    if original_path is not None:
        original, vo = load_raster(original_path)
        common = both & vo
        if common.any():
            moved_a = float(np.abs(a[common] - original[common]).sum())
            moved_b = float(np.abs(b[common] - original[common]).sum())
            if moved_b != 0.0:
                volume_ratio = moved_a / moved_b
            elif moved_a == 0.0:
                # Neither implementation changed anything: agreement is perfect.
                volume_ratio = 1.0

    return {
        "mask_identical": mask_identical,
        "n_valid": n,
        "max_abs": float(ad.max()),
        "mean_abs": float(ad.mean()),
        "frac_differing": float((ad > tol).sum() / n),
        "signed_min": float(d.min()),
        "signed_max": float(d.max()),
        "volume_ratio": volume_ratio,
    }


# lsdtt3's near-zero guard for barnes2014b_float_step: flats with |z| below
# this are raised by 2^-23 m per mask unit instead of by float32 steps.
FLAT_GUARD_LIMIT_M = 1.0


def load_raster_float32(path: pathlib.Path):
    """Return (array, valid_mask) as float32 / bool, refusing any other band
    type: a bitwise comparison of a raster that went through another dtype
    would compare the conversion, not the algorithm."""
    import numpy as np
    from osgeo import gdal
    gdal.UseExceptions()

    ds = gdal.Open(str(path))
    if ds is None:
        raise RuntimeError(f"could not open {path}")
    band = ds.GetRasterBand(1)
    if band.DataType != gdal.GDT_Float32:
        raise RuntimeError(f"{path.name}: band type "
                           f"{gdal.GetDataTypeName(band.DataType)}, not Float32")
    arr = band.ReadAsArray().astype(np.float32, copy=False)
    nodata = band.GetNoDataValue()
    valid = ~np.isnan(arr)
    if nodata is not None:
        valid &= (arr != np.float32(nodata))
    return arr, valid


def _ordered_key(a):
    """float32 -> int64, monotone in value, with -0.0 and +0.0 both 0, so
    key(b) - key(a) is the number of nextafter steps from a up to b."""
    import numpy as np
    i = np.asarray(a, np.float32).view(np.int32).astype(np.int64)
    return np.where(i < 0, -(i & 0x7FFFFFFF), i)


def _shifted(a, dr: int, dc: int, fill):
    """out[r, c] = a[r + dr, c + dc], `fill` where that is off the raster."""
    import numpy as np
    out = np.full_like(a, fill)
    R, C = a.shape
    rs, re_ = max(0, -dr), min(R, R - dr)
    cs, ce = max(0, -dc), min(C, C - dc)
    out[rs:re_, cs:ce] = a[rs + dr:re_ + dr, cs + dc:ce + dc]
    return out


def _flat_cells(z, valid):
    """Cells Barnes (2014b) raises: NO_FLOW cells (no lower neighbour, not
    draining) plus low edges (draining cells next to a NO_FLOW cell of equal
    elevation). Raster-border and nodata-adjacent cells drain. Same definition
    as check_flat_resolution_vs_richdem.py."""
    import numpy as np
    R, C = z.shape
    n8 = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if dr or dc]
    drains = np.zeros((R, C), bool)
    drains[0, :] = drains[-1, :] = drains[:, 0] = drains[:, -1] = True
    for dr, dc in n8:
        nv = _shifted(valid, dr, dc, True)
        nz = _shifted(z, dr, dc, np.float32(np.inf))
        drains |= ~nv | (nv & (nz < z))
    noflow = valid & ~drains
    low_edge = np.zeros((R, C), bool)
    for dr, dc in n8:
        nnf = _shifted(noflow, dr, dc, False)
        nz = _shifted(z, dr, dc, np.float32(np.nan))
        low_edge |= valid & drains & nnf & (nz == z)
    return noflow | low_edge


def _without_descent(z, valid):
    """Cells a D8 router cannot give a receiver: valid, off the border, no
    nodata neighbour, no strictly lower neighbour. Same predicate as
    check_flat_resolution_vs_richdem.py and lsdtt3's own count."""
    import numpy as np
    R, C = z.shape
    n8 = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if dr or dc]
    drains = np.zeros((R, C), bool)
    drains[0, :] = drains[-1, :] = drains[:, 0] = drains[:, -1] = True
    for dr, dc in n8:
        nv = _shifted(valid, dr, dc, True)
        nz = _shifted(z, dr, dc, np.float32(np.inf))
        drains |= ~nv | (nz < z)
    return valid & ~drains


CHANGED_AFTER_ROUND1_RE = r"Changed after round 1:\s*(\d+) cell\(s\)"


def changed_after_round1(stdout_path: pathlib.Path) -> int:
    """lsdtt3's "Changed after round 1: <n> cell(s)" from a saved stdout."""
    import re
    m = re.search(CHANGED_AFTER_ROUND1_RE, stdout_path.read_text())
    if not m:
        raise RuntimeError(f"'Changed after round 1' line not found in {stdout_path}")
    return int(m.group(1))


def float32_bit_metrics(lsdtt_path, ref_path, flat_baseline_path,
                        n_changed_after_round1: int = 0):
    """Bitwise comparison of two flat-resolved float32 DEMs.

    `flat_baseline_path` is the true-flat fill both resolved; it locates the
    flats (for the near-zero guard) and is the zero point for the raise in
    ULPs. A cell is GUARDED when it lies on a flat (off the raster border,
    which Barnes 2014b never raises) with |z| < FLAT_GUARD_LIMIT_M; every
    other valid cell must carry the same bits in both outputs, except the
    `n_changed_after_round1` cells lsdtt3 reports changing after its first
    (RichDEM-equivalent) round: `n_unexplained_differences` is the differing
    unguarded count minus that, and must be 0.
    """
    import numpy as np

    L, vL = load_raster_float32(lsdtt_path)
    R, vR = load_raster_float32(ref_path)
    F, vF = load_raster_float32(flat_baseline_path)
    if not (L.shape == R.shape == F.shape):
        raise RuntimeError(f"shape mismatch: lsdtt3 {L.shape}, reference "
                           f"{R.shape}, flat baseline {F.shape}")

    # Validity must agree three ways: the two outputs, and the fill both
    # started from. Folded into mask_identical, which is always enforced.
    mask_identical = bool(np.array_equal(vL, vR) and np.array_equal(vL, vF))
    valid = vL & vR & vF
    n = int(valid.sum())
    if n == 0:
        raise RuntimeError("no overlapping valid cells")

    border = np.zeros(F.shape, bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    flats = _flat_cells(F, vF) & ~border & valid
    guarded = flats & (np.abs(F) < FLAT_GUARD_LIMIT_M)
    differ = valid & (L.view(np.uint32) != R.view(np.uint32))
    unguarded = differ & ~guarded

    below = unguarded & (_ordered_key(L) < _ordered_key(R))
    raise_L = (_ordered_key(L) - _ordered_key(F))[valid]
    raise_R = (_ordered_key(R) - _ordered_key(F))[valid]

    first = ""
    if unguarded.any():
        r, c = (int(x) for x in np.argwhere(unguarded)[0])

        def show(z):    # 9 significant digits round-trip a float32
            return f"{float(z[r, c]):.9g} (0x{int(z.view(np.uint32)[r, c]):08x})"
        first = (f"row {r} col {c}: fill {show(F)}, lsdtt3 {show(L)}, "
                 f"reference {show(R)}")

    return {
        "mask_identical": mask_identical,
        "n_valid": n,
        "n_flat": int(flats.sum()),
        "n_differing": int(differ.sum()),
        "n_guarded": int(guarded.sum()),
        "n_differing_unguarded": int(unguarded.sum()),
        "n_changed_after_round1": int(n_changed_after_round1),
        "n_unexplained_differences": int(unguarded.sum()) - int(n_changed_after_round1),
        "n_below_richdem": int(below.sum()),
        "n_stranded_lsdtt": int(_without_descent(L, vL).sum()),
        "n_stranded_richdem": int(_without_descent(R, vR).sum()),
        "max_ulp_raise_lsdtt": int(raise_L.max()),
        "max_ulp_raise_reference": int(raise_R.max()),
        "max_abs": float(np.abs(L[valid].astype(np.float64)
                                - R[valid].astype(np.float64)).max()),
        "first_unguarded_cell": first,
    }


def float32_exact_metrics(lsdtt_path, ref_path, original_path):
    """Bitwise comparison of two elevation rasters at float32.

    The reference is rounded ONCE to float32 (Whitebox writes Float64; lsdtt3
    rounds its double working copy once at the end), then the float32 bit
    patterns are compared on every valid cell. Also reports which cells each
    side changed relative to the input, and the Jaccard index of those sets,
    so a failure says whether the two disagree on WHERE they carved or only on
    the values.
    """
    import numpy as np
    from osgeo import gdal
    gdal.UseExceptions()

    def read(path):
        ds = gdal.Open(str(path))
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()
        nodata = band.GetNoDataValue()
        valid = np.isfinite(arr)
        if nodata is not None:
            valid &= (arr != nodata)
        return arr, valid

    a, va = read(lsdtt_path)
    b, vb = read(ref_path)
    o, vo = read(original_path)
    if a.shape != b.shape or a.shape != o.shape:
        raise RuntimeError(f"shape mismatch {a.shape} / {b.shape} / {o.shape}")
    a32 = a.astype(np.float32)
    b32 = b.astype(np.float32)
    o32 = o.astype(np.float32)

    both = va & vb
    differ = both & (a32.view(np.uint32) != b32.view(np.uint32))
    diff = np.abs(a32.astype(np.float64) - b32.astype(np.float64))
    mod_a = va & vo & (a32.view(np.uint32) != o32.view(np.uint32))
    mod_b = vb & vo & (b32.view(np.uint32) != o32.view(np.uint32))
    union = int((mod_a | mod_b).sum())
    first = ""
    if differ.any():
        r, c = np.argwhere(differ)[0]
        first = f"({r},{c}) lsdtt3 {a32[r, c]!r} vs ref {b32[r, c]!r}"
    return {
        "mask_identical": bool(np.array_equal(va, vb)),
        "n_valid": int(both.sum()),
        "n_differing": int(differ.sum()),
        "max_abs": float(diff[both].max()) if both.any() else 0.0,
        "n_modified_lsdtt": int(mod_a.sum()),
        "n_modified_reference": int(mod_b.sum()),
        "modified_jaccard": (float((mod_a & mod_b).sum()) / union) if union else 1.0,
        "first_differing_cell": first,
    }


def system_env():
    """Environment for running lsdtt3 binaries from inside the pixi env.

    The driver runs under `pixi run -e benchmark` so that RichDEM is importable,
    but the lsdtt3 CLIs are C++ binaries built against the SYSTEM GDAL. Left
    alone, they inherit pixi's LD_LIBRARY_PATH and load pixi's libgdal.so.38,
    which needs a newer libstdc++ than the system one they were linked with:

        libstdc++.so.6: version `CXXABI_1.3.15' not found (required by
        .../.pixi/envs/benchmark/lib/libgdal.so.38)

    Stripping the environment-injected library and data paths makes each binary
    resolve the libraries it was actually built against. This bit
    lsdtt-drainage-analysis and not lsdtt-raster-preprocessing, purely because
    of which libraries each happens to pull in -- so it is worth removing for
    all of them rather than waiting for the next one to trip over it.
    """
    return lsdtt_test_env.system_env()


def cell_area_of(path: pathlib.Path) -> float:
    """Pixel area in projection units, from the geotransform."""
    from osgeo import gdal
    gdal.UseExceptions()
    gt = gdal.Open(str(path)).GetGeoTransform()
    return abs(gt[1] * gt[5])


def cell_size_of(path: pathlib.Path) -> float:
    """Cell size in projection units (the N-S and E-W pitch are equal here)."""
    from osgeo import gdal
    gdal.UseExceptions()
    gt = gdal.Open(str(path)).GetGeoTransform()
    return abs(gt[1])


def accumulation_metrics(lsdtt_path, richdem_path, cell_area, top_fraction=0.01):
    """Comparison metrics for two flow-accumulation rasters.

    UNIT CONVENTION: lsdtt3 reports drainage AREA in m^2 (its accumulate()
    weights each cell by cell_size^2); RichDEM reports a CELL COUNT, with
    sources = 1. Dividing the lsdtt3 raster by the cell area puts both in cell
    counts. Getting this wrong would show up as a constant 900x discrepancy on
    a 30 m grid rather than as anything subtle, but it is worth stating.
    """
    import numpy as np

    a, va = load_raster(lsdtt_path)
    b, vb = load_raster(richdem_path)

    if a.shape != b.shape:
        raise RuntimeError(f"shape mismatch {a.shape} vs {b.shape}")

    both = va & vb
    n = int(both.sum())
    if n == 0:
        raise RuntimeError("no overlapping valid cells")

    a_cells = a[both] / cell_area      # m^2 -> cell counts
    b_cells = b[both]

    # Rank correlation over all cells. Deliberately the primary structural
    # metric: D8 tie-breaking on flats differs between implementations, which
    # moves WHICH cell carries a value without changing the network. A
    # cell-wise diff punishes that; a rank correlation tolerates it while still
    # catching a genuine routing difference.
    try:
        from scipy.stats import spearmanr
        rho = float(spearmanr(a_cells, b_cells).statistic)
    except Exception:
        # Fall back to Pearson on ranks, which is the same quantity.
        ra = np.argsort(np.argsort(a_cells)).astype("float64")
        rb = np.argsort(np.argsort(b_cells)).astype("float64")
        rho = float(np.corrcoef(ra, rb)[0, 1])

    # Main channels: the top fraction of cells by accumulation. These carry the
    # signal most analyses actually depend on, and agreement there matters far
    # more than agreement on hillslope cells with accumulation of 1-3.
    k = max(1, int(round(top_fraction * n)))
    top_a = np.sort(a_cells)[-k:]
    top_b = np.sort(b_cells)[-k:]
    denom = np.where(top_b == 0.0, 1.0, top_b)
    top_rel = float(np.abs(top_a - top_b).max() / np.abs(denom).max())

    sum_a, sum_b = float(a_cells.sum()), float(b_cells.sum())
    total_rel = abs(sum_a - sum_b) / sum_b if sum_b else float("nan")

    return {
        "mask_identical": bool((va == vb).all()),
        "n_valid": n,
        "rank_correlation": rho,
        "top_channels_rel": top_rel,
        "total_rel": total_rel,
        "max_cells_lsdtt": float(a_cells.max()),
        "max_cells_richdem": float(b_cells.max()),
    }


def count_interior_minima(path: pathlib.Path) -> int:
    """Interior cells strictly lower than all eight neighbours.

    This is the most valuable line in the benchmark: it says the
    filled DEM is genuinely depression-free regardless of what RichDEM thinks,
    so it still has teeth if both implementations are wrong in the same way.

    The Phase 1 property tests assert this on synthetic surfaces; running it
    here extends the same guarantee to real terrain, where flats and nodata
    edges are far messier than anything a fractal generator produces.
    """
    import numpy as np

    arr, valid = load_raster(path)

    # Restrict to cells whose full 8-neighbourhood is in bounds and valid. Edge
    # cells and cells touching nodata are legitimate drainage exits under open
    # boundary conditions, so they can never be spurious minima.
    interior = valid.copy()
    interior[0, :] = interior[-1, :] = False
    interior[:, 0] = interior[:, -1] = False
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            interior[1:-1, 1:-1] &= valid[1 + dr:arr.shape[0] - 1 + dr,
                                          1 + dc:arr.shape[1] - 1 + dc]

    has_non_higher = np.zeros_like(interior)
    centre = arr[1:-1, 1:-1]
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            neighbour = arr[1 + dr:arr.shape[0] - 1 + dr,
                            1 + dc:arr.shape[1] - 1 + dc]
            has_non_higher[1:-1, 1:-1] |= (neighbour <= centre)

    return int((interior & ~has_non_higher).sum())


def ensure_filled_dem(exe: pathlib.Path, fixture: pathlib.Path,
                      work_dir: pathlib.Path, digest: str) -> pathlib.Path:
    """Fill the fixture once and reuse it for every accumulation comparison.

    Uses the DEFAULT fill_min_slope (1e-4), not 0, even though min_slope=0 is
    the variant proven bit-exact against RichDEM. Reason: lsdtt3's D8 router
    cannot route a truly flat surface --

        lsdtt_flowinfo: could not assign a receiver at row 12, col 129.
        The DEM likely contains an unresolved sink or flat.

    -- so a min_slope=0 fill produces a DEM that lsdtt-drainage-analysis
    refuses outright. That is by design: the imposed gradient is what makes
    flats routable (see docs/developer_notes/testing_and_benchmarking.md).

    Using the epsilon-filled DEM does not weaken the comparison. Both tools
    receive the SAME raster, so the accumulation algorithms are still compared
    in isolation; and because that raster has no flats, flat-handling is
    removed as a confounder for RichDEM too.
    """
    out_dir = work_dir / f"filled_{digest}"
    filled = out_dir / f"{fixture.stem}_fill.tif"
    if filled.exists():
        return filled
    run_lsdtt(exe, fixture, out_dir,
              {"write_fill_raster": "true", "fill_min_slope": "0.0001"})
    if not filled.exists():
        raise RuntimeError(f"expected filled DEM {filled} not written")
    return filled


def run_lsdtt(exe: pathlib.Path, fixture: pathlib.Path, out_dir: pathlib.Path,
              params: dict) -> float:
    """Run the lsdtt3 CLI on the fixture. Returns wall time in seconds."""
    out_dir.mkdir(parents=True, exist_ok=True)

    base = {
        "read_path": str(fixture.parent),
        "read_fname": fixture.name,
        "write_path": str(out_dir),
        "write_fname": fixture.stem,
        # The fixture is already in metric UTM (fetched that way). Letting
        # the CLI reproject would score the two tools on different rasters.
        "enforce_projection": "false",
        "write_output_format": "gtiff",
    }
    base.update(params)

    param_path = out_dir / "benchmark.param"
    param_path.write_text("".join(f"{k}: {v}\n" for k, v in base.items()))

    # The CLI takes the parameter file as <directory> <filename>, not a single
    # path (see tests/scripts/run_cli_smoke.py).
    t0 = time.perf_counter()
    proc = lsdtt_test_env.run([str(exe), str(out_dir), param_path.name],
                          capture_output=True, text=True, env=system_env())
    elapsed = time.perf_counter() - t0
    (out_dir / "lsdtt_stdout.txt").write_text(proc.stdout)

    if proc.returncode != 0:
        raise RuntimeError(
            f"{exe.name} failed (exit {proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-2000:]}\n"
            f"--- stderr ---\n{proc.stderr[-2000:]}")
    return elapsed


def run_reference(reference: str, routine: str, fixture: pathlib.Path,
                  out_path: pathlib.Path, cache: bool = True,
                  extra_args: list | None = None):
    """Run one reference implementation. Returns (seconds, cached_bool)."""
    if cache and out_path.exists():
        meta_path = out_path.with_suffix(".json")
        if meta_path.exists():
            return json.loads(meta_path.read_text())["seconds"], True

    out_path.parent.mkdir(parents=True, exist_ok=True)
    proc = lsdtt_test_env.run(
        [sys.executable, str(SCRIPT_DIR / REFERENCE_SCRIPTS[reference]),
         "--routine", routine, "--in", str(fixture), "--out", str(out_path)]
        + list(extra_args or []),
        capture_output=True, text=True)

    if proc.returncode == SKIP_EXIT:
        die_skip(f"{reference} not available. Use `pixi run -e benchmark ...`")
    if proc.returncode != 0:
        raise RuntimeError(f"{REFERENCE_SCRIPTS[reference]} failed:\n"
                           f"{proc.stderr[-2000:]}")

    info = json.loads(proc.stdout.strip().splitlines()[-1])
    out_path.with_suffix(".json").write_text(json.dumps(info))
    return info["seconds"], False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", help="path to lsdtt-raster-preprocessing")
    ap.add_argument("--drainage-exe",
                    help="path to lsdtt-drainage-analysis (needed for the "
                         "accum_* comparisons)")
    ap.add_argument("--fetch-exe",
                    help="path to lsdtt-fetch-raster (used to fetch missing "
                         "fixtures; defaults to build/lsdtt-fetch-raster)")
    ap.add_argument("--fixture", default="alpine_small",
                    help="fixture basename in tests/data/dem_cache (default alpine_small)")
    ap.add_argument("--comparison", choices=list(COMPARISONS) + ["all"], default="all")
    # --out-dir is the name the other two drivers use; --work-dir is the
    # original spelling, kept working so existing scripts and job files do not
    # break. Both write to the same place. This directory holds the RESULTS
    # (benchmark_results.csv) AND the scratch: staged filled DEMs and the
    # per-reference output caches, which run to hundreds of MB on the larger
    # fixtures. Point it somewhere with room.
    ap.add_argument("--out-dir", "--work-dir", dest="out_dir", default=None,
                    help="results and scratch directory. Default "
                         "build/flow_benchmark inside the repo. Give the "
                         "campaign directory to keep everything together.")
    ap.add_argument("--fixture-dir", default="",
                    help="where fixture DEMs are cached and fetched to. "
                         "Default tests/data/dem_cache inside the repo, which "
                         "puts downloads in your working tree. Point it at the "
                         "campaign's benchmark_fixtures directory to share the "
                         "cache with run_benchmark_suite.py and avoid "
                         "re-downloading.")
    ap.add_argument("--check", action="store_true",
                    help="fail on pinned tolerance violations (Phase 3)")
    ap.add_argument("--list", action="store_true", help="list comparisons and exit")
    args = ap.parse_args()

    if args.list:
        for name, spec in COMPARISONS.items():
            print(f"{name:18s} {spec['note']}")
        return 0

    try:
        import numpy  # noqa: F401
        from osgeo import gdal  # noqa: F401
    except ImportError:
        die_skip("numpy/GDAL not importable. Use `pixi run -e benchmark ...`")

    # Rebind the module global before anything can fetch: ensure_fixture()
    # reads it at call time, so this redirects both the lookup and the
    # download.
    if args.fixture_dir:
        global FIXTURE_DIR
        FIXTURE_DIR = pathlib.Path(args.fixture_dir).resolve()

    fetch_exe = pathlib.Path(args.fetch_exe).resolve() if args.fetch_exe else \
        (REPO_ROOT / "build" / "lsdtt-fetch-raster")
    fixture = ensure_fixture(args.fixture, fetch_exe)

    # Short content hash of the fixture, used in the RichDEM cache key below.
    fixture_digest = hashlib.sha256(fixture.read_bytes()).hexdigest()[:12]

    # Tolerances are pinned against a specific fixture build. A different GDAL
    # produces a slightly different raster (a one-row registration shift was
    # measured between GDAL 3.4.1 and 3.12.3), and pinned numbers do not
    # transfer. Skip rather than fail: a different GDAL is an environment
    # difference, not a defect.
    expected_sha = EXPECTED_FIXTURE_SHA.get(args.fixture)
    if args.check and expected_sha and expected_sha != fixture_digest:
        die_skip(
            f"fixture {args.fixture} has sha {fixture_digest}, but the pinned "
            f"tolerances were measured against {expected_sha}. This normally "
            f"means a different GDAL build fetched the fixture. Re-measure the "
            f"tolerances (run without --check) and update TOLERANCES and "
            f"EXPECTED_FIXTURE_SHA together before trusting a pass.")

    if not args.exe:
        die_skip("--exe not given")
    exe = pathlib.Path(args.exe).resolve()
    if not exe.exists():
        die_skip(f"executable {exe} not found; build lsdtt3 first")

    work_dir = pathlib.Path(args.out_dir).resolve() if args.out_dir else \
        REPO_ROOT / "build" / "flow_benchmark"
    work_dir.mkdir(parents=True, exist_ok=True)

    drainage_exe = pathlib.Path(args.drainage_exe).resolve() if args.drainage_exe \
        else (REPO_ROOT / "build" / "lsdtt-drainage-analysis")

    names = list(COMPARISONS) if args.comparison == "all" else [args.comparison]

    rows, failures = [], []
    for name in names:
        spec = COMPARISONS[name]
        print(f"\n=== {name} on {args.fixture} ===")
        print(f"    {spec['note']}")

        case_dir = work_dir / f"{args.fixture}_{name}"
        if case_dir.exists():
            shutil.rmtree(case_dir)

        if spec.get("stage") == "drainage":
            if drainage_exe is None or not drainage_exe.exists():
                die_skip(f"{name} needs lsdtt-drainage-analysis "
                         f"(looked for {drainage_exe}); pass --drainage-exe")
            stage_exe   = drainage_exe
            stage_input = ensure_filled_dem(exe, fixture, work_dir, fixture_digest)
        else:
            stage_exe   = exe
            stage_input = fixture

        lsdtt_s = run_lsdtt(stage_exe, stage_input, case_dir, spec["params"])
        lsdtt_out = case_dir / f"{stage_input.stem}{spec['lsdtt_suffix']}"
        if not lsdtt_out.exists():
            raise RuntimeError(f"expected lsdtt3 output {lsdtt_out} not written")

        # Depression-free property is a property of the lsdtt3 output alone,
        # so it is computed once, not per reference.
        # The true-flat fill the flats were resolved from, for the bitwise
        # comparison's guard and ULP count. Same lsdtt3 executable and input.
        flat_baseline = None
        if spec.get("metric") == "float32_bits":
            base_dir = case_dir / "flat_baseline"
            run_lsdtt(stage_exe, stage_input, base_dir,
                      spec["flat_baseline_params"])
            flat_baseline = base_dir / f"{stage_input.stem}_fill.tif"
            if not flat_baseline.exists():
                raise RuntimeError(f"expected flat baseline {flat_baseline} "
                                   f"not written")

        minima = None
        if spec["kind"] == "elevation" and name in ("fill", "fill_epsilon"):
            minima = count_interior_minima(lsdtt_out)
            print(f"    interior minima: {minima}  (must be 0)")
        if minima is not None and minima > 0:
            failures.append(
                f"{name}: filled DEM still contains {minima} interior minima "
                f"-- it is not depression-free")

        # One comparison per available reference. Both are scored on the SAME
        # lsdtt3 output, so any difference between the two rows is a difference
        # between the references, not between lsdtt3 runs.
        for reference, routine in spec["references"].items():
            key = f"{name}@{reference}"

            # Key the cache on the fixture's CONTENT, not just its name. Keying
            # on the name alone silently reuses a reference computed from a
            # different raster whenever the fixture is regenerated -- which
            # happened when fixture generation moved from pixi's GDAL to the C++
            # tool's, and produced a bogus "99.9 % of cells differ" instead of
            # the true 0 %.
            # lsdtt3's `fill_min_slope` is a GRADIENT: the rise imposed across a
            # flat is `min_slope * distance`, so on a 30 m grid a min_slope of
            # 1e-4 raises a cardinal neighbour by 3e-3 m. Whitebox's
            # `flat_increment` is metres per cell directly. Handing it the raw
            # min_slope compares two different increments -- measured at 30x
            # apart, which showed up as a volume ratio of 25 on the flat-dominated
            # lakeland fixture. Convert.
            #
            # RichDEM is deliberately NOT converted: its epsilon=True steps by
            # one ULP of std::nextafter, which is not an increment in metres at
            # all and has no matching value.
            extra: list[str] = list(spec.get("reference_args", {}).get(reference, []))
            if (reference == "whitebox"
                    and "fill_min_slope" in spec.get("params", {})):
                min_slope = float(spec["params"]["fill_min_slope"])
                if min_slope > 0.0:
                    cell = cell_size_of(stage_input)
                    extra += ["--flat-increment", repr(min_slope * cell)]

            # Cache key must include the extra args, or a run with a different
            # increment silently reuses the previous raster.
            suffix = ("_" + "_".join(extra).replace("-", "").replace(".", "p")
                      if extra else "")
            ref_out = (work_dir / f"{reference}_cache" /
                       f"{args.fixture}_{fixture_digest}_{routine}{suffix}.tif")
            ref_s, cached = run_reference(reference, routine, stage_input,
                                          ref_out, extra_args=extra)

            if spec["kind"] == "accumulation":
                m = accumulation_metrics(lsdtt_out, ref_out,
                                         cell_area_of(stage_input))
            elif spec.get("metric") == "float32_exact":
                m = float32_exact_metrics(lsdtt_out, ref_out, stage_input)
            elif spec.get("metric") == "float32_bits":
                m = float32_bit_metrics(
                    lsdtt_out, ref_out, flat_baseline,
                    changed_after_round1(case_dir / "lsdtt_stdout.txt"))
            else:
                m = elevation_metrics(lsdtt_out, ref_out, original_path=fixture)

            report_only = key in REPORT_ONLY
            print(f"    -- vs {reference}{'  (report-only)' if report_only else ''}")
            print(f"       mask identical: {m['mask_identical']}")
            if spec["kind"] == "accumulation":
                print(f"       rank correlation: {m['rank_correlation']:.6f}")
                print(f"       top 1% channels:  {100.0 * m['top_channels_rel']:.4f} % rel")
                print(f"       total accum:      {100.0 * m['total_rel']:.6f} % rel")
                print(f"       max cells:        lsdtt3 {m['max_cells_lsdtt']:.0f}"
                      f" / {reference} {m['max_cells_richdem']:.0f}")
            elif spec.get("metric") == "float32_exact":
                print(f"       cells compared:   {m['n_valid']}")
                print(f"       bits differing:   {m['n_differing']}  (must be 0)")
                print(f"       max |diff|:       {m['max_abs']:.9g} m")
                print(f"       cells modified:   lsdtt3 {m['n_modified_lsdtt']}"
                      f" / {reference} {m['n_modified_reference']}"
                      f"  (Jaccard {m['modified_jaccard']:.6f})")
                if m["first_differing_cell"]:
                    print(f"       first differing:  {m['first_differing_cell']}")
            elif spec.get("metric") == "float32_bits":
                print(f"       cells compared:   {m['n_valid']}")
                print(f"       flat cells:       {m['n_flat']}")
                print(f"       bits differing:   {m['n_differing']}"
                      f"  (guarded |z| < {FLAT_GUARD_LIMIT_M:g} m: "
                      f"{m['n_guarded']} cells on flats)")
                print(f"       differing outside guard: "
                      f"{m['n_differing_unguarded']}  (lsdtt3 changed "
                      f"{m['n_changed_after_round1']} after round 1; unexplained "
                      f"{m['n_unexplained_differences']}, must be 0)")
                print(f"       lsdtt3 below {reference}: {m['n_below_richdem']}"
                      f"  (must be 0)")
                print(f"       without descent:  lsdtt3 {m['n_stranded_lsdtt']}"
                      f" (must be 0) / {reference} {m['n_stranded_richdem']}")
                if m["n_guarded"] and m["n_changed_after_round1"]:
                    failures.append(f"{key}: fixture mixes guarded and loop-changed "
                                    f"cells; extend the check")
                print(f"       max ULP raise:    lsdtt3 {m['max_ulp_raise_lsdtt']}"
                      f" / {reference} {m['max_ulp_raise_reference']}")
                print(f"       max |diff|:       {m['max_abs']:.9g} m")
                if m["first_unguarded_cell"]:
                    print(f"       first unguarded:  {m['first_unguarded_cell']}")
            else:
                print(f"       max |diff|:       {m['max_abs']:.4f} m")
                print(f"       cells >1cm:       {100.0 * m['frac_differing']:.4f} %")
                print(f"       volume ratio:     {m['volume_ratio']:.9f}")
            print(f"       lsdtt3 {lsdtt_s:.3f} s   {reference} {ref_s:.3f} s"
                  f"{'  (cached)' if cached else ''}")

            rows.append({"fixture": args.fixture, "fixture_sha": fixture_digest,
                         "comparison": name, "reference": reference,
                         "lsdtt_seconds": lsdtt_s, "reference_seconds": ref_s,
                         "interior_minima": "" if minima is None else minima,
                         **m})

            # Always enforced: a nodata mismatch means the two rasters do not
            # describe the same cells, so every other metric is meaningless.
            if not m["mask_identical"]:
                failures.append(f"{key}: nodata masks differ")

            if not args.check:
                continue
            if report_only:
                print("       (report-only: not asserted)")
                continue

            tol = TOLERANCES.get(args.fixture, {}).get(key)
            if tol is None:
                failures.append(
                    f"{key}: --check requested but no tolerances are pinned "
                    f"for fixture '{args.fixture}'. Run without --check, then "
                    f"pin TOLERANCES and EXPECTED_FIXTURE_SHA together.")
                continue

            for metric, limit in tol.items():
                if limit is None:
                    continue
                if metric == "volume_ratio_max":
                    ratio = m["volume_ratio"]
                    if ratio != ratio or ratio <= 0.0:
                        continue
                    value = max(ratio, 1.0 / ratio)
                elif metric == "rank_correlation_min":
                    # A minimum, not a maximum: flag when correlation FALLS below.
                    if m["rank_correlation"] < limit:
                        failures.append(
                            f"{key}: rank_correlation = "
                            f"{m['rank_correlation']:.6g} below pinned {limit:g}")
                    continue
                elif metric == "n_unexplained_differences":
                    # Signed in the row; either sign is a failure.
                    value = abs(m[metric])
                else:
                    value = m[metric]

                if value > limit:
                    first = m.get("first_unguarded_cell") or m.get("first_differing_cell")
                    detail = f"; first at {first}" if first else ""
                    failures.append(
                        f"{key}: {metric} = {value:.6g} exceeds pinned "
                        f"{limit:g}{detail}")

    csv_path = work_dir / "benchmark_results.csv"
    write_header = not csv_path.exists()
    # Union of every row's keys, in first-seen order -- NOT rows[0].keys().
    # Elevation comparisons (fill, breach) and accumulation comparisons
    # (accum_*) emit different metric sets, so taking the header from the first
    # row alone made any run containing both abort at the very end with
    # "dict contains fields not in fieldnames: rank_correlation, ...", after
    # all the work was done and the results printed. restval fills the gaps for
    # rows that lack a column.
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with open(csv_path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, restval="")
        if write_header:
            w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {csv_path}")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        return 1

    if not args.check:
        print("\n(reporting only -- tolerances not pinned yet, docs/developer_notes/testing_and_benchmarking.md)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
