#!/usr/bin/env python3
"""Seed-ensemble statistics for the cosmogenic transient landslide experiment.

Workflow port of p5_analyse.py from the hand-run experiment
(data/cosmo_transient_test/lsdtt3_version/ls_phase3/p5_analyse.py). The
statistics are unchanged; the run layout, the ensemble axes and the output
paths are now supplied on the command line instead of being hardcoded.

For each (magnitude, seed, scenario) it reads the joined sample table and the
matching per-magnitude control, forms the PAIRED delta

    d = ratio_scenario - ratio_control

per sample point -- matched on (latitude, longitude), which is stable because
every run in a magnitude samples the same channel-network point set -- takes
the MEDIAN of those paired deltas as the per-seed statistic, and then
summarises the distribution of those medians ACROSS SEEDS with a median and an
IQR.

Validation emitted alongside the statistics:

  * anchor invariance. Samples whose control sc_steptime_yr == 0 lie entirely
    upstream of the knickpoint, so no landslide can be upstream of them and
    their delta must be zero to floating-point. The count and the worst
    |delta| are written per combination, and the run FAILS if any exceeds
    --anchor-alarm.
  * pairing integrity. Any control sample with no scenario partner (or vice
    versa) is counted and reported rather than silently dropped.
  * scenario equality assertions, e.g. S1 == S1f, which must hold exactly.

Data wrangling and descriptive statistics only.
"""
import argparse
import csv
import sys
from pathlib import Path


def load(path):
    if not path.exists():
        sys.exit(f"FATAL: missing {path}")
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {(r["latitude"], r["longitude"]): r for r in rows}


def fnum(r, k):
    try:
        return float(r[k])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def median(v):
    v = sorted(v)
    n = len(v)
    if n == 0:
        return float("nan")
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def quartiles(v):
    """Medians of the split halves, matching the hand-run report."""
    v = sorted(v)
    n = len(v)
    lo = v[: n // 2]
    hi = v[(n + 1) // 2:]
    return median(lo), median(hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs-root", required=True)
    ap.add_argument("--magnitudes", required=True, help="comma separated, e.g. 5x,10x")
    ap.add_argument("--seeds", required=True, help="comma separated, e.g. s01,s02")
    ap.add_argument("--scenarios", required=True, help="comma separated, e.g. S1,S1f,S2,S3")
    ap.add_argument("--control", required=True, help="control scenario name, e.g. S0")
    ap.add_argument("--ratios", required=True, help="comma separated ratio columns")
    ap.add_argument("--anchor-eps", type=float, default=1e-6)
    ap.add_argument("--anchor-alarm", type=float, default=1e-4)
    ap.add_argument("--expected-anchors", type=int, default=0,
                    help="expected anchor sample count; 0 disables the check")
    ap.add_argument("--assert-equal", default="",
                    help="semicolon separated pairs of scenarios that must be "
                         "identical, e.g. 'S1=S1f'")
    ap.add_argument("--out-per-seed", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-anchors", required=True)
    ap.add_argument("--out-report", required=True)
    args = ap.parse_args()

    root = Path(args.runs_root).resolve()
    runs = [x for x in args.magnitudes.split(",") if x]
    seeds = [x for x in args.seeds.split(",") if x]
    scenarios = [x for x in args.scenarios.split(",") if x]
    ratios = [x for x in args.ratios.split(",") if x]
    equal_pairs = [tuple(p.split("=", 1))
                   for p in args.assert_equal.split(";") if "=" in p]

    per_seed = []
    anchor_rows = []
    mismatch = []
    pairing = []
    counts = {}
    lines = []

    def say(msg=""):
        print(msg)
        lines.append(msg)

    for run in runs:
        ctrl = load(root / run / "CTRL" / f"{args.control}_joined_samples.csv")
        for seed in seeds:
            inv = root / run / seed / "ls_LS_INVENTORY.csv"
            with open(inv) as f:
                counts[(run, seed)] = sum(1 for _ in f) - 1
            tabs = {sc: load(root / run / seed / f"{sc}_joined_samples.csv")
                    for sc in scenarios}

            for a, b in equal_pairs:
                if a not in tabs or b not in tabs:
                    continue
                ndiff = 0
                for key, r in tabs[a].items():
                    o = tabs[b].get(key)
                    if o is None or any(r[c] != o[c] for c in r):
                        ndiff += 1
                if ndiff:
                    mismatch.append((run, seed, a, b, ndiff))

            for sc in scenarios:
                tab = tabs[sc]
                # Pairing integrity: the comparison is per sample point against
                # the matching control, so an unmatched key on either side is a
                # defect, not something to average over.
                missing_in_scen = sum(1 for k in ctrl if k not in tab)
                extra_in_scen = sum(1 for k in tab if k not in ctrl)
                if missing_in_scen or extra_in_scen:
                    pairing.append((run, seed, sc, missing_in_scen, extra_in_scen))

                for ratio in ratios:
                    d, anch = [], []
                    for key, c in ctrl.items():
                        s = tab.get(key)
                        if s is None:
                            continue
                        a, b = fnum(s, ratio), fnum(c, ratio)
                        if a != a or b != b:
                            continue
                        dd = a - b
                        d.append(dd)
                        if fnum(c, "sc_steptime_yr") == 0.0:
                            anch.append(abs(dd))
                    if not d:
                        continue
                    worst_anchor = max(anch) if anch else 0.0
                    per_seed.append({
                        "run": run, "scenario": sc, "ratio": ratio,
                        "seed": seed, "n_samples": len(d),
                        "n_landslides": counts[(run, seed)],
                        "median_delta": median(d),
                        "n_anchor": len(anch),
                        "max_abs_anchor_delta": worst_anchor,
                    })
                    anchor_rows.append({
                        "run": run, "seed": seed, "scenario": sc, "ratio": ratio,
                        "n_anchor": len(anch),
                        "max_abs_anchor_delta": f"{worst_anchor:.6g}",
                        "status": ("ALARM" if worst_anchor > args.anchor_alarm
                                   else "note" if worst_anchor > args.anchor_eps
                                   else "pass"),
                    })

    if not per_seed:
        sys.exit("FATAL: no per-seed statistics produced.")

    Path(args.out_per_seed).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_per_seed, "w", newline="") as f:
        cols = ["run", "scenario", "ratio", "seed", "n_samples",
                "n_landslides", "median_delta", "n_anchor",
                "max_abs_anchor_delta"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(per_seed)

    with open(args.out_anchors, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(anchor_rows[0].keys()))
        w.writeheader()
        w.writerows(anchor_rows)

    summary = []
    for run in runs:
        for sc in scenarios:
            for ratio in ratios:
                m = [r["median_delta"] for r in per_seed
                     if r["run"] == run and r["scenario"] == sc
                     and r["ratio"] == ratio]
                if not m:
                    continue
                q1, q3 = quartiles(m)
                med = median(m)
                if q1 > 0 and q3 > 0:
                    verdict = "AMPLIFIES"
                elif q1 < 0 and q3 < 0:
                    verdict = "DAMPS"
                else:
                    verdict = "NOISE_ONLY"
                summary.append({
                    "run": run, "scenario": sc, "ratio": ratio,
                    "n_seeds": len(m),
                    "ensemble_median": f"{med:.6g}",
                    "q1": f"{q1:.6g}", "q3": f"{q3:.6g}",
                    "iqr": f"{q3 - q1:.6g}",
                    "min": f"{min(m):.6g}", "max": f"{max(m):.6g}",
                    "frac_seeds_positive":
                        f"{sum(1 for x in m if x > 0) / len(m):.3f}",
                    "iqr_excludes_zero": "yes" if verdict != "NOISE_ONLY" else "no",
                    "verdict": verdict,
                })

    with open(args.out_summary, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    say("=== ENSEMBLE SUMMARY ===")
    say(f"{'run':>4} {'scen':>4} {'ratio':>10} {'median':>11} {'IQR':>10} "
        f"{'q1':>11} {'q3':>11} {'min':>11} {'max':>11} {'f>0':>5}  verdict")
    for r in summary:
        say(f"{r['run']:>4} {r['scenario']:>4} "
            f"{r['ratio'].replace('app_erate_ratio_',''):>10} "
            f"{float(r['ensemble_median']):>11.6f} {float(r['iqr']):>10.6f} "
            f"{float(r['q1']):>11.6f} {float(r['q3']):>11.6f} "
            f"{float(r['min']):>11.6f} {float(r['max']):>11.6f} "
            f"{float(r['frac_seeds_positive']):>5.2f}  {r['verdict']}")

    failures = []

    say("")
    say("=== ANCHOR INVARIANCE ===")
    worst = max(anchor_rows, key=lambda t: float(t["max_abs_anchor_delta"]))
    say(f"max |anchor delta| over all {len(anchor_rows)} "
        f"(run,seed,scenario,ratio) combinations = "
        f"{float(worst['max_abs_anchor_delta']):.3e} at "
        f"{worst['run']} {worst['seed']} {worst['scenario']} {worst['ratio']}")
    anchor_counts = sorted({int(t["n_anchor"]) for t in anchor_rows})
    say(f"anchor sample count per combination: {anchor_counts}")
    if args.expected_anchors and anchor_counts != [args.expected_anchors]:
        msg = (f"*** FAIL: expected {args.expected_anchors} anchor samples in "
               f"every combination, saw {anchor_counts}")
        say(msg)
        failures.append(msg)
    bad = [t for t in anchor_rows if t["status"] == "ALARM"]
    warn = [t for t in anchor_rows if t["status"] == "note"]
    if bad:
        msg = (f"*** ALARM: {len(bad)} combinations exceed {args.anchor_alarm:g} "
               f"-- DIVIDE ARTEFACT MAY BE BACK ***")
        say(msg)
        for t in bad[:20]:
            say(f"    {t['run']} {t['seed']} {t['scenario']} {t['ratio']}: "
                f"{t['max_abs_anchor_delta']}")
        failures.append(msg)
    elif warn:
        say(f"note: {len(warn)} combinations between {args.anchor_eps:g} and "
            f"{args.anchor_alarm:g} (above float32 epsilon, far below artefact scale)")
    else:
        say(f"all anchor deltas <= {args.anchor_eps:g} -- PASS")

    say("")
    say("=== SAMPLE PAIRING ===")
    if pairing:
        msg = f"*** FAIL: {len(pairing)} scenario tables do not pair 1:1 with the control"
        say(msg)
        for t in pairing[:20]:
            say(f"    {t[0]} {t[1]} {t[2]}: {t[3]} control-only, {t[4]} scenario-only")
        failures.append(msg)
    else:
        say("every scenario sample pairs 1:1 with the control -- PASS")

    if equal_pairs:
        say("")
        say("=== SCENARIO EQUALITY ASSERTIONS ===")
        if mismatch:
            msg = "*** FAIL: scenarios asserted equal are not equal"
            say(msg)
            for t in mismatch:
                say(f"    {t[0]} {t[1]}: {t[2]} != {t[3]} in {t[4]} rows")
            failures.append(msg)
        else:
            for a, b in equal_pairs:
                say(f"{a} == {b} exactly (all columns, all rows) for every "
                    f"run/seed -- PASS")

    say("")
    say("=== REALISED LANDSLIDE COUNT PER SEED ===")
    for run in runs:
        c = [counts[(run, s)] for s in seeds]
        say(f"  {run}: " + " ".join(f"{s}={counts[(run, s)]}" for s in seeds))
        say(f"        median={median(c):.0f} min={min(c)} max={max(c)}")

    Path(args.out_report).write_text("\n".join(lines) + "\n", encoding="utf-8")

    if failures:
        sys.exit("ensemble validation failed:\n  " + "\n  ".join(failures))


if __name__ == "__main__":
    main()
