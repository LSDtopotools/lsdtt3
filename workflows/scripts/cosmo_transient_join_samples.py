#!/usr/bin/env python3
"""Join the sampled-concentration table to the three per-nuclide inversions.

Workflow port of join_samples.py from the hand-run experiment
(data/cosmo_transient_test/lsdtt3_version/ls_phase3/join_samples.py). The
logic is unchanged; the one behavioural difference is deliberate:

    the original resolved every path relative to __file__, so it only ever
    worked when it sat inside the run directory. This version takes the run
    directory as an argument.

Data wrangling only: every column written here is either copied verbatim from
an lsdtt3 output or is one of the two apparent-erosion-rate ratios, which are
plain divisions of emitted columns. No physics is computed in this script.

Reads   <run-dir>/<stem>_conc_sampled_crn_conc.csv
        <run-dir>/<stem>_inv_{Be10,Al26,C14}_CN_ERATE.csv
        <run-dir>/<stem>_conc.param      (for rock_density_kg_m3)
Writes  the --out CSV.

The join is NOT on latitude/longitude: lsdtt-cosmo-tool reports the
coordinates of the channel node it snapped the sample onto, which is generally
a different cell from the one the sample was taken at. Rows are instead matched
on the concentration value, which the inversion echoes back verbatim in
`concentration_<nuclide>` -- a key that is unique per sample and survives the
snap. Ties and misses are reported rather than silently dropped.
"""
import argparse
import csv
import sys
from pathlib import Path

NUCLIDES = ("Be10", "Al26", "C14")


def load(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def fnum(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def read_density(param_path):
    """Read rock_density_kg_m3 from a param file rather than hardcoding it.

    This is the ROCK density (2650). It is not the scar density (1500), which
    lives under scar_density_kg_m3 and is used for the geometric scar volume.
    """
    try:
        with open(param_path) as f:
            for line in f:
                key, _, val = line.partition(":")
                if key.strip() == "rock_density_kg_m3":
                    return float(val.split("#")[0].strip())
    except OSError:
        pass
    sys.exit(f"FATAL: could not read rock_density_kg_m3 from {param_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", required=True,
                    help="directory holding the run's outputs and param files")
    ap.add_argument("--stem", required=True,
                    help="scenario stem, e.g. S1 (files are <stem>_conc_*)")
    ap.add_argument("--out", required=True, help="output joined CSV")
    # The inversion prints concentrations at lower precision than the sampler,
    # so match on a relative tolerance rather than equality.
    ap.add_argument("--rel-tol", type=float, default=1e-6)
    args = ap.parse_args()

    here = Path(args.run_dir).resolve()
    stem = args.stem
    rel_tol = args.rel_tol

    rho_kg_m3 = read_density(here / f"{stem}_conc.param")
    print(f"rock density: {rho_kg_m3:g} kg/m3 (from {stem}_conc.param)")
    sampled_path = here / f"{stem}_conc_sampled_crn_conc.csv"
    sampled = load(sampled_path)
    print(f"sampled rows: {len(sampled)}  ({sampled_path.name})")

    # nuclide -> list of (concentration, inversion row).
    inv = {}
    for nuc in NUCLIDES:
        p = here / f"{stem}_inv_{nuc}_CN_ERATE.csv"
        if not p.exists():
            sys.exit(f"FATAL: missing {p.name} -- run the {nuc} inversion first.")
        rows = load(p)
        inv[nuc] = [(fnum(r, f"concentration_{nuc}"), r) for r in rows]
        print(f"  {nuc}: {len(rows)} inverted rows ({p.name})")

    out_rows = []
    unmatched = {n: 0 for n in NUCLIDES}
    for s in sampled:
        rec = {
            "latitude": s["latitude"],
            "longitude": s["longitude"],
            "snapped_contributing_pixels": s.get("snapped_contributing_pixels", ""),
            # The basin's ACTUAL mean erosion rate above this point (arithmetic
            # upslope mean). This is the ground truth the apparent rates are
            # scored against.
            "upstream_mean_eff_erate_g_cm2_yr":
                s.get("upstream_mean_eff_erate_g_cm2_yr", ""),
            "sc_steptime_yr": s.get("sc_steptime_yr", ""),
            "sc_erate_m_yr": s.get("sc_erate_m_yr", ""),
            "sc_eff_erate_g_cm2_yr": s.get("sc_eff_erate_g_cm2_yr", ""),
            "sc_elevation_m": s.get("sc_elevation_m", ""),
            "Be10_CONC": s.get("Be10_CONC", ""),
            "Al26_CONC": s.get("Al26_CONC", ""),
            "C14_CONC": s.get("C14_CONC", ""),
            "R2610_ACC": s.get("R2610_ACC", ""),
            "R1410_ACC": s.get("R1410_ACC", ""),
        }
        app = {}
        for nuc in NUCLIDES:
            target = fnum(s, f"{nuc}_CONC")
            best, best_err = None, float("inf")
            for conc, r in inv[nuc]:
                if conc == 0.0:
                    continue
                err = abs(conc - target) / abs(target) if target else abs(conc)
                if err < best_err:
                    best, best_err = r, err
            if best is None or best_err > rel_tol:
                unmatched[nuc] += 1
                rec[f"app_eff_erate_{nuc}_g_cm2_yr"] = ""
                rec[f"predicted_conc_{nuc}"] = ""
                app[nuc] = float("nan")
            else:
                v = fnum(best, "apparent_eff_e_g_cm2_yr")
                rec[f"app_eff_erate_{nuc}_g_cm2_yr"] = f"{v:.10g}"
                rec[f"predicted_conc_{nuc}"] = best.get(f"predicted_conc_{nuc}", "")
                app[nuc] = v

        # The figures colour by kyr, and the plotting layer must only draw
        # columns it is handed.
        t_yr = fnum(s, "sc_steptime_yr")
        rec["sc_steptime_kyr"] = f"{t_yr / 1000.0:.10g}" if t_yr == t_yr else ""

        # Unit conversion, not a computation: g/cm2/yr -> mm/yr at the rock
        # density pinned in the param file.
        eff = fnum(s, "upstream_mean_eff_erate_g_cm2_yr")
        rec["upstream_mean_erate_mm_yr"] = (
            f"{eff * 10.0 / (rho_kg_m3 / 1000.0):.10g}" if eff == eff else "")

        # The only derived quantities, and only divisions of emitted columns.
        for num, den, name in (("C14", "Be10", "app_erate_ratio_C14_Be10"),
                               ("C14", "Al26", "app_erate_ratio_C14_Al26")):
            a, b = app.get(num, float("nan")), app.get(den, float("nan"))
            rec[name] = f"{a / b:.10g}" if b and b == b and a == a and b != 0 else ""
        out_rows.append(rec)

    cols = ["latitude", "longitude", "snapped_contributing_pixels",
            "upstream_mean_eff_erate_g_cm2_yr", "upstream_mean_erate_mm_yr",
            "sc_steptime_yr", "sc_steptime_kyr",
            "sc_erate_m_yr", "sc_eff_erate_g_cm2_yr",
            "sc_elevation_m",
            "Be10_CONC", "Al26_CONC", "C14_CONC", "R2610_ACC", "R1410_ACC",
            "predicted_conc_Be10", "predicted_conc_Al26", "predicted_conc_C14",
            "app_eff_erate_Be10_g_cm2_yr", "app_eff_erate_Al26_g_cm2_yr",
            "app_eff_erate_C14_g_cm2_yr",
            "app_erate_ratio_C14_Be10", "app_erate_ratio_C14_Al26"]

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(out_rows)

    print(f"\nwrote {out_path.name}: {len(out_rows)} rows")
    for nuc in NUCLIDES:
        if unmatched[nuc]:
            print(f"  WARNING: {unmatched[nuc]} sample(s) had no {nuc} match "
                  f"within rel tol {rel_tol:g}")

    # Round-trip check: the inversion's forward prediction must reproduce the
    # concentration it was given.
    print("\nround-trip (predicted_conc vs input concentration):")
    for nuc in NUCLIDES:
        worst, n = 0.0, 0
        for r in out_rows:
            c, p = fnum(r, f"{nuc}_CONC"), fnum(r, f"predicted_conc_{nuc}")
            if c == c and p == p and c:
                worst = max(worst, abs(p - c) / abs(c))
                n += 1
        print(f"  {nuc}: n={n}  worst relative error = {worst:.3e}")

    # Recovery check. The apparent rate is inverted from a BASIN-AVERAGED
    # concentration, so the only meaningful comparison is against the basin's
    # mean erosion rate. The local pixel rate is shown alongside purely to make
    # the difference between the two visible.
    print("\nrecovered apparent rate vs TRUE basin mean (local rate for contrast):")
    for label, pred in (("t == 0 (upstream, old rate)",
                         lambda r: fnum(r, "sc_steptime_yr") == 0.0),
                        ("t > 200000 (long adjusted, new rate)",
                         lambda r: fnum(r, "sc_steptime_yr") > 200000.0)):
        sel = [r for r in out_rows if pred(r)]
        if not sel:
            continue
        local = [fnum(r, "sc_eff_erate_g_cm2_yr") for r in sel]
        true_ = [fnum(r, "upstream_mean_eff_erate_g_cm2_yr") for r in sel]
        true_mean = sum(true_) / len(true_)
        print(f"  {label}: n={len(sel)}")
        print(f"      TRUE basin mean = {true_mean:.4f}   "
              f"(local pixel rate = {sum(local)/len(local):.4f})")
        for nuc in NUCLIDES:
            v = [fnum(r, f"app_eff_erate_{nuc}_g_cm2_yr") for r in sel]
            v = [x for x in v if x == x]
            if v:
                m = sum(v) / len(v)
                print(f"      {nuc}: apparent = {m:.4f}  "
                      f"apparent/true = {m / true_mean:.4f}")


if __name__ == "__main__":
    main()
