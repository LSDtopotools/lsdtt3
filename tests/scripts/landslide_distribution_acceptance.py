#!/usr/bin/env python3
"""Compare a realised landslide inventory against the distribution that was asked for.

Reads ``<write_fname>_LS_INVENTORY.csv`` from lsdtt-landslides and reports the
realised frequency-area distribution beside the requested one: the tail
exponent, the mean and modal areas, and where the two diverge.

WHY THIS EXISTS. The size distribution in this module is IMPOSED, so the
realised inventory should reproduce it. It does not always, and the ways it
fails are informative rather than random:

  * A per-event volume that is small compared with E[V] makes the budget solve
    SELECT sizes rather than sample them. At one landslide per event the bias is
    total -- realised areas cluster at (V_target/alpha)^(1/gamma). The CLI warns,
    but this script shows the damage.
  * The A_max cutoff removes almost nothing by number and most of the volume, so
    a realised distribution can look perfect while the volume budget is
    dominated by the truncation.
  * A cell size coarser than the modal area (~400 m2 for the Malamud fit) makes
    the small-area rollover unrepresentable, so an inverse-gamma request comes
    back looking like a pure power law. That is a resolution limit, not an error.

The tail exponent is fit by maximum likelihood (Clauset et al. 2009) rather than
by least squares on a log-log histogram, which is biased by binning.

DO NOT COMPARE THE FITTED EXPONENT WITH rho+1. The inverse gamma only reaches
its asymptotic slope for A >> a: the local log-log slope is

    d ln p / d ln A  =  -(rho+1) + a/(A-s)

which for the Malamud parameters is -0.38 at the mode, -2.15 at 5,000 m2, -2.36
at 30,000 m2, and reaches -2.40 only around 10^6 m2. Because A_max is capped by
hillslope size -- a landslide cannot exceed the hillslope it comes from -- a
realistic run has NO range over which the asymptote is measurable. Fitting from
the mode over a realistic range returns about 1.73, and that is CORRECT
behaviour, not a bias.

This script therefore compares the fitted exponent with the exponent the TARGET
distribution would itself produce over the same fit range, computed by numerical
integration. That comparison is meaningful; the one against rho+1 is not.

Usage:
  landslide_distribution_acceptance.py --inventory run_LS_INVENTORY.csv \\
      [--model inverse_gamma|truncated_pareto] [--rho 1.40] [--a 1280] \\
      [--s -132] [--alpha-pdf 2.3] [--cell-size 10] [--out-csv report.csv]

No plotting and no third-party dependencies beyond the standard library, so it
runs anywhere the tests do.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_inventory(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def mle_pareto_exponent(areas: list[float], a_min: float) -> tuple[float, int]:
    """Maximum-likelihood PDF exponent for a power law above a_min.

    Clauset et al. (2009) Eq. 3.1: alpha = 1 + n / sum(ln(x_i / x_min)).
    Returns (alpha_pdf, n_used). Least-squares on a binned log-log histogram is
    the common alternative and is biased; this is not.
    """
    tail = [a for a in areas if a >= a_min]
    if len(tail) < 10:
        return float("nan"), len(tail)
    acc = sum(math.log(a / a_min) for a in tail)
    if acc <= 0.0:
        return float("nan"), len(tail)
    return 1.0 + len(tail) / acc, len(tail)


def inverse_gamma_mean(rho: float, a: float, s: float) -> float:
    """Malamud et al. (2004) Eq. 14. Diverges for rho <= 1."""
    if rho <= 1.0:
        return float("inf")
    return a / (rho - 1.0) + s


def inverse_gamma_pdf(area: float, rho: float, a: float, s: float) -> float:
    x = area - s
    if x <= 0.0:
        return 0.0
    return math.exp(-math.log(a) - math.lgamma(rho)
                    + (rho + 1.0) * (math.log(a) - math.log(x)) - a / x)


def expected_mle_exponent(a_min: float, a_max: float, model: str,
                          rho: float, a: float, s: float,
                          alpha_pdf: float, n_steps: int = 4000) -> float:
    """The exponent an MLE fit would return FROM THIS DISTRIBUTION over [a_min, a_max].

    The population MLE satisfies alpha = 1 + 1/E[ln(A/a_min)], with the
    expectation taken over the target density restricted to the fit range. This
    is the number the realised inventory should reproduce -- not the asymptotic
    tail exponent, which a truncated range cannot reach.
    """
    if not (a_max > a_min > 0.0):
        return float("nan")

    # Integrated numerically for both families. A closed form exists for the
    # bounded Pareto but is easy to get subtly wrong, and this runs once.
    lo, hi = math.log(a_min), math.log(a_max)
    h = (hi - lo) / n_steps
    norm = 0.0
    acc = 0.0
    for i in range(n_steps + 1):
        u = lo + i * h
        A = math.exp(u)
        if model == "inverse_gamma":
            dens = inverse_gamma_pdf(A, rho, a, s)
        else:
            dens = A ** (-alpha_pdf)
        wgt = 1.0 if i in (0, n_steps) else (4.0 if i % 2 else 2.0)
        # integrating in log space: dA = A du
        norm += wgt * dens * A
        acc += wgt * dens * A * math.log(A / a_min)
    if norm <= 0.0 or acc <= 0.0:
        return float("nan")
    return 1.0 + norm / acc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inventory", required=True, type=Path,
                    help="<write_fname>_LS_INVENTORY.csv from lsdtt-landslides")
    ap.add_argument("--model", default="inverse_gamma",
                    choices=["inverse_gamma", "truncated_pareto"])
    ap.add_argument("--rho", type=float, default=1.40)
    ap.add_argument("--a", type=float, default=1280.0, help="inverse gamma a (m2)")
    ap.add_argument("--s", type=float, default=-132.0, help="inverse gamma s (m2)")
    ap.add_argument("--alpha-pdf", type=float, default=2.3)
    ap.add_argument("--cell-size", type=float, default=None,
                    help="metres; inferred from the smallest area if omitted")
    ap.add_argument("--out-csv", type=Path, default=None)
    args = ap.parse_args()

    rows = read_inventory(args.inventory)
    if not rows:
        raise SystemExit(f"No landslides in {args.inventory}")

    areas = [float(r["area_m2"]) for r in rows]
    pixels = [int(r["n_pixels"]) for r in rows]
    volumes = [float(r["volume_m3"]) for r in rows]
    areas.sort()

    cell = args.cell_size
    if cell is None:
        cell = math.sqrt(min(areas) / max(1, min(pixels)))

    n = len(areas)
    mean_area = sum(areas) / n
    total_vol = sum(volumes)

    # Tail fit above the modal scale. For the inverse gamma the mode is near
    # a/(rho+1)+s; below it the distribution rolls over and is not a power law,
    # so fitting through it would report an exponent that means nothing.
    if args.model == "inverse_gamma":
        mode = args.a / (args.rho + 1.0) + args.s
        expected_tail = args.rho + 1.0
        expected_mean = inverse_gamma_mean(args.rho, args.a, args.s)
    else:
        mode = min(areas)
        expected_tail = args.alpha_pdf
        expected_mean = float("nan")

    a_min_fit = max(mode, min(areas))
    fitted, n_tail = mle_pareto_exponent(areas, a_min_fit)
    expected_over_range = expected_mle_exponent(
        a_min_fit, max(areas), args.model, args.rho, args.a, args.s,
        args.alpha_pdf)

    print(f"Inventory      : {args.inventory}")
    print(f"Landslides     : {n}")
    print(f"Cell size      : {cell:g} m  ({cell*cell:g} m2 per pixel)")
    print(f"Area range     : {min(areas):g} - {max(areas):g} m2")
    print(f"Mean area      : {mean_area:.1f} m2")
    print(f"Total volume   : {total_vol:.1f} m3")
    print()
    print(f"Requested model: {args.model}")
    if args.model == "inverse_gamma":
        print(f"  rho={args.rho}  a={args.a} m2  s={args.s} m2")
        print(f"  modal area   : {mode:.1f} m2")
        print(f"  expected mean: {expected_mean:.1f} m2 (Malamud Eq. 14, UNtruncated)")
    else:
        print(f"  alpha_pdf={args.alpha_pdf}")
    print(f"  ASYMPTOTIC tail exponent: {expected_tail:.3f}"
          f"   (reached only for A >> a; usually NOT measurable)")
    print(f"  exponent this distribution gives over the fit range: "
          f"{expected_over_range:.3f}")
    print()
    print(f"Realised tail (MLE above {a_min_fit:.1f} m2, n={n_tail}):")
    print(f"  fitted exponent: {fitted:.3f}")
    if not math.isnan(fitted) and not math.isnan(expected_over_range):
        print(f"  vs this distribution over the same range: "
              f"{fitted - expected_over_range:+.3f}   <-- THE MEANINGFUL COMPARISON")
        print(f"  vs the asymptote                        : "
              f"{fitted - expected_tail:+.3f}   (expected to be large and negative)")
    print()

    # --- Diagnoses, in the order they usually bite -------------------------
    problems = []

    if cell * cell >= mode and args.model == "inverse_gamma":
        problems.append(
            f"RESOLUTION: one pixel ({cell*cell:g} m2) is at or above the modal "
            f"area ({mode:.0f} m2), so the rollover cannot be represented. The "
            f"realised inventory will look like a pure power law however the "
            f"model was configured.")

    if (not math.isnan(fitted) and not math.isnan(expected_over_range)
            and abs(fitted - expected_over_range) > 0.25):
        problems.append(
            f"TAIL EXPONENT off by {fitted - expected_over_range:+.2f} from what "
            f"this distribution gives over the same range ({expected_over_range:.2f}). "
            f"This is the comparison that matters -- a large gap here means the "
            f"realised sizes are NOT the requested distribution.")

    frac_at_min = sum(1 for p in pixels if p == min(pixels)) / n
    if frac_at_min > 0.5:
        problems.append(
            f"CLUSTERING AT THE FLOOR: {100*frac_at_min:.0f}% of landslides are "
            f"the minimum size. Either the distribution is unresolved at this "
            f"cell size, or the per-event budget is too small.")

    distinct = len(set(pixels))
    if distinct < 10 and n > 50:
        problems.append(
            f"SELECTION, NOT SAMPLING: only {distinct} distinct sizes across "
            f"{n} landslides. A per-event volume comparable to E[V] makes the "
            f"budget solve pick sizes rather than sample them. Raise the "
            f"per-event volume to many times E[V].")

    if args.model == "inverse_gamma" and math.isfinite(expected_mean):
        rel = abs(mean_area - expected_mean) / expected_mean
        if rel > 0.5:
            problems.append(
                f"MEAN AREA off by {100*rel:.0f}% from the untruncated "
                f"expectation. Some of this is legitimate -- A_max truncation "
                f"lowers the realised mean -- but a large gap alongside few "
                f"distinct sizes points at budget selection.")

    if problems:
        print("DIAGNOSIS")
        for p in problems:
            print("  - " + p)
    else:
        print("DIAGNOSIS: realised inventory is consistent with the request.")

    if args.out_csv is not None:
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.out_csv.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["metric", "value"])
            w.writerow(["n_landslides", n])
            w.writerow(["cell_size_m", cell])
            w.writerow(["min_area_m2", min(areas)])
            w.writerow(["max_area_m2", max(areas)])
            w.writerow(["mean_area_m2", mean_area])
            w.writerow(["total_volume_m3", total_vol])
            w.writerow(["asymptotic_tail_exponent", expected_tail])
            w.writerow(["expected_exponent_over_fit_range", expected_over_range])
            w.writerow(["fitted_tail_exponent", fitted])
            w.writerow(["n_tail", n_tail])
            w.writerow(["distinct_sizes", distinct])
            w.writerow(["fraction_at_min_size", frac_at_min])
            w.writerow(["n_problems", len(problems)])
        print(f"\nWrote {args.out_csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
