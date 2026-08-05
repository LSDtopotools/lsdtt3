#!/usr/bin/env python3
"""Check a realised Niemi ensemble against the claims Niemi et al. (2005) make.

Reads the basin CSVs the Niemi ensemble emits and evaluates, one by one, the
statements in the paper that lsdtt3 claims to reproduce. Each is reported with
the number that supports or refutes it, so a run that fails says *which* claim
failed rather than "acceptance failed".

WHICH RATE YOU DIVIDE BY IS THE WHOLE PROBLEM
---------------------------------------------
A landsliding landscape has more than one defensible "true" erosion rate and the
answer to "is the cosmogenic estimate biased?" depends entirely on the choice:

  imposed      the long-term rate the scenario was asked to deliver. No spread --
               it is a number someone chose.
  volumetric   E_v, the mass actually leaving the basin AT THE SAMPLING INSTANT:
               background erosion plus whatever debris is in transit. Large
               spread.
  apparent     E_c, what an analyst gets -- the flux-weighted basin-averaged
               concentration inverted through the steady-state assumption.

**Niemi's comparison is E_c against E_v, and both are snapshots.** His sec. 2.4.2
defines the volumetric rate as "the total depth of material removed from the
watershed divided by the watershed area multiplied by the model TIME STEP" --
the same per-timestep depth array that weights his Eq. 6 mixture. His Fig. 5
shows E_v with a LARGER spread than the cosmogenic rates at every catchment
scale, which a long-term average could not have.

Dividing by the imposed rate instead and reading the shortfall as a model bias
is the mistake this script exists to stop anyone repeating: it cost lsdtt3 a
fictitious "unexplained factor of 1.7". Niemi reports the shortfall himself, for
his own model, in sec. 3.1 -- and reports it for the VOLUMETRIC rate too.

THE CLAIMS CHECKED
------------------
  C1  the ensemble is large enough and spans enough basin area to say anything
  C2  (sec. 3.1) the median of BOTH the cosmogenic and the volumetric rate falls
      below the imposed rate. Checked on the volumetric rate, which contains no
      cosmogenics at all -- so this is a statement about episodic erosion, not
      about the CRN model.
  C3  against the volumetric rate the apparent rate is close to unbiased: the
      interquartile range brackets unity and the median is near it
  C4  (sec. 3.1) "At small catchment scales, CRN-derived rates are typically
      higher than volumetric rates ... but converge at larger catchment scales"
  C5  (sec. 3.1) the scatter falls with basin area
  C6  (sec. 4.3) the median across basins is closer to the truth than a typical
      single sample -- his practical conclusion
  C7  (Eq. 7) the area above which samples cluster, reported and compared with
      A_avg[km2] = 100 / E[km/Myr]. REPORTED ONLY on a narrow-area fixture,
      where the criterion cannot bind.

Usage:
  niemi_ensemble_acceptance.py --basins niemi_basins.csv \\
      --volumetric niemi_basins_volumetric.csv \\
      [--imposed-erate-g-cm2-yr 0.00575] [--out-csv report.csv] [--strict]

--strict makes an unmet claim exit non-zero. Without it every claim is reported
and the exit status is 0 unless a file is unreadable, which is the right default
for a diagnostic run on a new landscape.

No plotting and no third-party dependencies beyond the standard library, so it
runs anywhere the tests do.
"""

import argparse
import csv
import math
import pathlib
import sys


# ---------------------------------------------------------------------------
# statistics, written out rather than imported, so this runs with a bare python
# ---------------------------------------------------------------------------

def quantile(sorted_vals, q):
    """Linear-interpolated quantile of an already-sorted list."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def summarise(vals):
    s = sorted(vals)
    return {
        "n": len(s),
        "median": quantile(s, 0.5),
        "p25": quantile(s, 0.25),
        "p75": quantile(s, 0.75),
        "p05": quantile(s, 0.05),
        "p95": quantile(s, 0.95),
        "min": s[0] if s else float("nan"),
        "max": s[-1] if s else float("nan"),
    }


def read_basins(path):
    """Rows of a niemi_basins CSV, valid ones only, numbers coerced."""
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if str(r.get("valid", "1")).strip() not in ("1", "true", "True"):
                continue
            try:
                rows.append({
                    "realisation": int(float(r["realisation"])),
                    "area_km2": float(r["area_km2"]),
                    "fraction_failed": float(r["fraction_failed"]),
                    "apparent": float(r["apparent_erate_g_cm2_yr"]),
                    "true": float(r["true_erate_g_cm2_yr"]),
                    "ratio": float(r["ratio"]),
                })
            except (KeyError, ValueError) as exc:
                raise SystemExit(f"{path}: unreadable row: {exc}")
    if not rows:
        raise SystemExit(f"{path}: no valid basin rows")
    return rows


def log_area_bands(rows, bands_per_decade=2):
    """Group rows into log-spaced area bands, smallest first."""
    areas = [r["area_km2"] for r in rows if r["area_km2"] > 0.0]
    lo, hi = min(areas), max(areas)
    if not (hi > lo):
        return []
    width = 1.0 / bands_per_decade
    n_bands = max(1, int(math.ceil((math.log10(hi) - math.log10(lo)) / width)))
    out = []
    for i in range(n_bands):
        a_lo = lo * 10.0 ** (i * width)
        a_hi = lo * 10.0 ** ((i + 1) * width)
        sel = [r for r in rows
               if r["area_km2"] >= a_lo and (r["area_km2"] < a_hi
                                             or (i == n_bands - 1
                                                 and r["area_km2"] <= a_hi))]
        if len(sel) >= 5:
            st = summarise([r["ratio"] for r in sel])
            st["area_lo_km2"] = a_lo
            st["area_hi_km2"] = a_hi
            st["iqr"] = st["p75"] - st["p25"]
            out.append(st)
    return out


# ---------------------------------------------------------------------------
# the claims
# ---------------------------------------------------------------------------

class Report:
    def __init__(self):
        self.claims = []      # (id, met, headline, detail)

    def claim(self, cid, met, headline, detail=""):
        self.claims.append((cid, met, headline, detail))

    def note(self, cid, headline, detail=""):
        self.claims.append((cid, None, headline, detail))

    def emit(self):
        print()
        for cid, met, headline, detail in self.claims:
            mark = "ok  " if met else ("FAIL" if met is False else "note")
            print(f"  {mark}: {cid}  {headline}")
            for line in detail.splitlines():
                if line.strip():
                    print(f"          {line}")
        n_fail = sum(1 for _, met, _, _ in self.claims if met is False)
        print()
        print(f"  claims checked: "
              f"{sum(1 for _, m, _, _ in self.claims if m is not None)}, "
              f"unmet: {n_fail}")
        return n_fail


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Check a Niemi ensemble against the paper's claims.")
    ap.add_argument("--basins", required=True, type=pathlib.Path,
                    help="niemi_basins.csv -- apparent against the IMPOSED rate")
    ap.add_argument("--volumetric", type=pathlib.Path, default=None,
                    help="niemi_basins_volumetric.csv -- apparent against the "
                         "SNAPSHOT volumetric rate. This is Niemi's comparison; "
                         "without it C2-C4 cannot be evaluated.")
    ap.add_argument("--imposed-erate-g-cm2-yr", type=float, default=None,
                    help="the long-term rate the scenario imposed. Defaults to "
                         "the median 'true' column of --basins, which IS the "
                         "imposed rate when that file was written with the "
                         "imposed denominator.")
    ap.add_argument("--out-csv", type=pathlib.Path, default=None)
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any claim is unmet")
    args = ap.parse_args(argv)

    rep = Report()

    imposed_rows = read_basins(args.basins)
    vol_rows = read_basins(args.volumetric) if args.volumetric else None

    imposed = args.imposed_erate_g_cm2_yr
    if imposed is None:
        imposed = summarise([r["true"] for r in imposed_rows])["median"]

    print("=" * 62)
    print(" Niemi et al. (2005) ensemble acceptance")
    print("=" * 62)
    print(f"  basins (imposed denominator): {args.basins}")
    if args.volumetric:
        print(f"  basins (volumetric denom.) : {args.volumetric}")
    print(f"  imposed erosion rate        : {imposed:.6g} g/cm2/yr")

    # -- C1 -----------------------------------------------------------------
    areas = sorted(r["area_km2"] for r in imposed_rows)
    decades = math.log10(areas[-1] / areas[0]) if areas[0] > 0 else 0.0
    n_real = len(set(r["realisation"] for r in imposed_rows))
    rep.claim("C1", len(imposed_rows) >= 100 and decades >= 1.0,
              "the ensemble can support a distributional claim",
              f"{len(imposed_rows)} basin samples over {n_real} realisations; "
              f"area {areas[0]:.4g}-{areas[-1]:.4g} km2 ({decades:.2f} decades)")

    st_imposed = summarise([r["ratio"] for r in imposed_rows])
    print(f"\n  apparent / IMPOSED : median {st_imposed['median']:.4g}"
          f"  (p25 {st_imposed['p25']:.4g}, p75 {st_imposed['p75']:.4g})")

    if vol_rows is None:
        rep.note("C2-C4", "SKIPPED -- no --volumetric file supplied",
                 "Niemi's comparison is against the snapshot volumetric rate.\n"
                 "Without it the only ratio available is against the imposed\n"
                 "rate, whose median falls below 1 in his model too, so a low\n"
                 "value here is NOT evidence of a biased cosmogenic estimate.")
    else:
        st_vol = summarise([r["ratio"] for r in vol_rows])
        ev = summarise([r["true"] for r in vol_rows])
        print(f"  apparent / VOLUMETRIC: median {st_vol['median']:.4g}"
              f"  (p25 {st_vol['p25']:.4g}, p75 {st_vol['p75']:.4g})")
        print(f"  volumetric rate      : median {ev['median']:.4g} g/cm2/yr"
              f"  (E_v/E_imposed = {ev['median'] / imposed:.4g})")

        # -- C2 ---------------------------------------------------------------
        rep.claim("C2", ev["median"] < 0.9 * imposed,
                  "the VOLUMETRIC rate alone median-under-reads the imposed "
                  "rate (sec. 3.1)",
                  f"median E_v/E_imposed = {ev['median'] / imposed:.4g}. "
                  f"No cosmogenics enter this number.\n"
                  f"It is the whole reason apparent/imposed "
                  f"({st_imposed['median']:.4g}) sits below 1.")

        # -- C3 ---------------------------------------------------------------
        brackets = st_vol["p25"] < 1.0 < st_vol["p75"]
        near = 0.85 < st_vol["median"] < 1.15
        rep.claim("C3", brackets and near,
                  "against the volumetric rate the apparent rate is close to "
                  "unbiased",
                  f"IQR [{st_vol['p25']:.4g}, {st_vol['p75']:.4g}] "
                  f"{'brackets' if brackets else 'does NOT bracket'} unity; "
                  f"median {st_vol['median']:.4g} "
                  f"{'within' if near else 'OUTSIDE'} 15 % of unity")

        # -- C4 ---------------------------------------------------------------
        bands = log_area_bands(vol_rows)
        if len(bands) >= 3:
            print("\n  E_c/E_v by area band")
            print("    area (km2)              n    median       IQR")
            for b in bands:
                print(f"    {b['area_lo_km2']:8.4g} - {b['area_hi_km2']:8.4g}"
                      f"  {b['n']:5d}  {b['median']:8.4g}  {b['iqr']:8.4g}")
            declines = all(bands[i]["median"] <= bands[i - 1]["median"]
                           for i in range(1, len(bands)))
            starts_high = bands[0]["median"] > 1.0
            rep.claim("C4", declines and starts_high,
                      "small basins read high against volumetric and converge "
                      "with area (sec. 3.1)",
                      f"smallest band {bands[0]['median']:.4g}, largest "
                      f"{bands[-1]['median']:.4g}; "
                      f"{'monotonic' if declines else 'NOT monotonic'}")
        else:
            rep.note("C4", "too few populated area bands to test the trend",
                     f"{len(bands)} bands with >= 5 samples")

    # -- C5 -----------------------------------------------------------------
    ibands = log_area_bands(imposed_rows)
    if len(ibands) >= 2:
        rep.claim("C5", ibands[-1]["iqr"] <= ibands[0]["iqr"],
                  "the scatter falls with basin area (sec. 3.1)",
                  f"IQR {ibands[0]['iqr']:.4g} in the smallest band, "
                  f"{ibands[-1]['iqr']:.4g} in the largest")
    else:
        rep.note("C5", "too few area bands to compare scatter")

    # -- C6 -----------------------------------------------------------------
    #
    # DO NOT compare |median(r) - 1| with median(|r - 1|). When every ratio
    # falls below 1 -- which is the normal case here, and is Niemi's own
    # sec. 3.1 result -- those two are ALGEBRAICALLY IDENTICAL:
    #
    #     r < 1 for all r  =>  |r - 1| = 1 - r
    #                      =>  median(|r-1|) = 1 - median(r) = |median(r) - 1|
    #
    # so the comparison can only be won by the handful of samples above 1 and
    # says nothing about averaging. It is a degenerate check dressed as a
    # result.
    #
    # Niemi's sec. 4.3 claim is about SAMPLING: "measuring CRN-derived erosion
    # rates from several third- or fourth-order catchments should allow" a
    # representative answer. The testable content is that the median of a
    # handful of basins is a much more REPEATABLE estimator than one basin --
    # i.e. the spread of per-realisation medians is far smaller than the spread
    # of individual samples. That is what makes his advice good advice.
    ref6 = vol_rows if vol_rows else imposed_rows
    per_real = {}
    for r in ref6:
        per_real.setdefault(r["realisation"], []).append(r["ratio"])
    med_by_real = [summarise(v)["median"] for v in per_real.values()
                   if len(v) >= 5]
    if len(med_by_real) >= 5:
        st_meds = summarise(med_by_real)
        st_singles = summarise([r["ratio"] for r in ref6])
        iqr_meds = st_meds["p75"] - st_meds["p25"]
        iqr_singles = st_singles["p75"] - st_singles["p25"]
        rep.claim("C6", iqr_meds < 0.9 * iqr_singles,
                  "the median over several basins is a more repeatable "
                  "estimator than one basin (sec. 4.3)",
                  f"IQR of per-realisation medians {iqr_meds:.4g} against "
                  f"{iqr_singles:.4g} for individual samples "
                  f"({iqr_meds / iqr_singles:.3g}x), over "
                  f"{len(med_by_real)} realisations.\n"
                  f"THE IMPROVEMENT IS MODEST, AND THAT IS A RESULT. Sampled "
                  f"basins are NESTED --\n"
                  f"a large one contains the small ones and they share one "
                  f"landslide field -- so they\n"
                  f"are far from independent and the median does not average "
                  f"down like sqrt(n).\n"
                  f"Niemi's 46 sub-basins are nested too, so this caveat "
                  f"applies to his advice as\n"
                  f"much as to ours: 'sample several catchments' buys less "
                  f"than it appears to\n"
                  f"unless the catchments are genuinely independent.")
    else:
        rep.note("C6", "too few realisations to compare medians with singles",
                 f"{len(med_by_real)} realisations with >= 5 basins")

    # -- C7 -----------------------------------------------------------------
    # Niemi Eq. 7: A_avg [km2] = 100 [km3/Myr] / E [km/Myr].
    # E in km/Myr is (g/cm2/yr) / rho * 10 with rho in g/cm3; at 2.65 that is
    # e_eff * 3.7736 km/Myr.
    e_km_myr = imposed / 2.65 * 10.0
    a_avg = 100.0 / e_km_myr if e_km_myr > 0 else float("inf")
    ref = vol_rows if vol_rows else imposed_rows
    ref_median = summarise([r["ratio"] for r in ref])["median"]
    bands = log_area_bands(ref)
    a_meet = None
    for b in bands:
        within = [1 for b2 in [b]
                  if abs(b["median"] - ref_median) < 0.2 * abs(ref_median)]
        if within and b["iqr"] < 0.4 * abs(ref_median):
            a_meet = b["area_lo_km2"]
            break
    rep.note("C7",
             "Niemi Eq. 7 averaging area -- REPORTED, not asserted",
             f"Eq. 7 predicts A_avg = {a_avg:.4g} km2 at E = {e_km_myr:.4g} "
             f"km/Myr.\n"
             + (f"Observed clustering sets in by {a_meet:.4g} km2."
                if a_meet is not None else
                "The criterion never binds within the sampled area range, so "
                "the scaling\ncannot be measured here. It needs a DEM spanning "
                "more decades of area.")
             + f"\nSampled range is {areas[0]:.4g}-{areas[-1]:.4g} km2; Eq. 7 "
               f"is only testable if\nthat range brackets A_avg.")

    n_fail = rep.emit()

    if args.out_csv:
        with open(args.out_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["claim", "met", "headline"])
            for cid, met, headline, _ in rep.claims:
                w.writerow([cid,
                            "" if met is None else ("1" if met else "0"),
                            headline])
        print(f"  wrote {args.out_csv}")

    if args.strict and n_fail:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
