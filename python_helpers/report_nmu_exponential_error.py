#!/usr/bin/env python3
"""Freeze the diagnosis behind the removal of balcoCRONUSv3/balcoCRONUSv3PD.

Measures, over the whole Balco CRONUS v3 tabulated model-1A muon inventory:

  1. the error in v3's shipped single-exponential muon approximation
     (make_consts_v3.m:193-211) -- the model the removed lsdtt3 schemes
     implemented;
  2. the surface-production limit (eps -> 0), which isolates Pmu0 / L_mu_atm
     from the attenuation length;
  3. the effective attenuation length the table actually implies, against
     Balco's Lambda_mu(p) = 1/(a + b p);
  4. how many exponentials of lsdtt3's own kernel form are needed to reproduce
     the table -- the feasibility number behind the replacement scheme.

Reads the FROZEN extraction (tests/expected/cosmo/balco_v3_nmu_reference.csv),
not the GPL reference tree, so it is reproducible without crn_reference/.

Usage:
  python3 python_helpers/report_nmu_exponential_error.py \
      --reference tests/expected/cosmo/balco_v3_nmu_reference.csv \
      --out tests/expected/cosmo/nmu_exponential_error_report.csv
"""

from __future__ import annotations

import argparse
import csv
import datetime
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

# make_consts_v3.m:208-211.  Index 4 = Be-10 quartz, 5 = C-14 quartz,
# 7 = Al-26 quartz in v3's 8-element nuclide vectors.
PMU0 = {"Be10": 0.07345, "Al26": 0.67644, "C14": 3.06690}
L_MU_ATM = {"Be10": 299.2, "Al26": 288.0, "C14": 267.8}
LMU_A, LMU_B = 0.01036, -9.697e-6

# make_consts_v3.m:134,142,147
LAMBDA = {
    "Be10": math.log(2.0) / 1.387e6,
    "Al26": 9.83e-7,
    "C14": math.log(2.0) / 5730.0,
}

NUCLIDES = ("Be10", "Al26", "C14")


def read_reference(path: Path):
    """Return (meta_lines, {nuclide: (p_axis, log10e_axis, table[ie, ip])})."""
    meta = [l.rstrip("\n") for l in path.open() if l.startswith("#")]
    rows = [l for l in path.open() if not l.startswith("#")]
    recs = list(csv.DictReader(rows))

    out = {}
    for nuc in NUCLIDES:
        sub = [r for r in recs if r["nuclide"] == nuc]
        p_ax = sorted({float(r["p_hPa"]) for r in sub})
        e_ax = sorted({float(r["log10_erate"]) for r in sub})
        tab = np.full((len(e_ax), len(p_ax)), np.nan)
        pi = {v: i for i, v in enumerate(p_ax)}
        ei = {v: i for i, v in enumerate(e_ax)}
        for r in sub:
            tab[ei[float(r["log10_erate"])], pi[float(r["p_hPa"])]] = \
                float(r["N_mu_atoms_g"])
        if np.isnan(tab).any():
            raise ValueError(f"{nuc}: reference grid has holes")
        out[nuc] = (np.array(p_ax), np.array(e_ax), tab)
    return meta, out


def v3_exponential(nuc, pressure, erate):
    """v3's shipped fast approximation -- what the removed schemes computed."""
    p_mu = PMU0[nuc] * np.exp((1013.25 - pressure) / L_MU_ATM[nuc])
    lmu = 1.0 / (LMU_A + LMU_B * pressure)
    return p_mu * lmu / (LAMBDA[nuc] * lmu + erate)


def fit_k_exponentials(k, n_col, erates, lam, seeds=6):
    """Fit N = sum_k P_k L_k /(eps + L_k lambda) to one pressure slice."""
    def residual(q):
        P, L = np.exp(q[:k]), np.exp(q[k:])
        model = sum(P[i] * L[i] / (erates + L[i] * lam) for i in range(k))
        return np.log(model) - np.log(n_col)

    best = None
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        q0 = np.concatenate([
            np.log(np.full(k, 0.05)) + rng.normal(0, 0.5, k),
            np.log(np.geomspace(500.0, 8000.0, k)) + rng.normal(0, 0.3, k),
        ])
        r = least_squares(residual, q0, method="lm", max_nfev=20000)
        if best is None or r.cost < best.cost:
            best = r
    P, L = np.exp(best.x[:k]), np.exp(best.x[k:])
    model = sum(P[i] * L[i] / (erates + L[i] * lam) for i in range(k))
    return float(np.max(np.abs(model / n_col - 1.0)) * 100.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference",
                    default="tests/expected/cosmo/balco_v3_nmu_reference.csv")
    ap.add_argument("--out",
                    default="tests/expected/cosmo/nmu_exponential_error_report.csv")
    ap.add_argument("--no-fit", action="store_true",
                    help="skip the K-exponential feasibility summary")
    args = ap.parse_args()

    ref_path = Path(args.reference)
    meta, grids = read_reference(ref_path)

    fit_summary = {}
    if not args.no_fit:
        for nuc in NUCLIDES:
            p_ax, e_ax, tab = grids[nuc]
            erates = 10.0 ** e_ax
            for k in (1, 2, 3):
                errs = [fit_k_exponentials(k, tab[:, j], erates, LAMBDA[nuc])
                        for j in range(len(p_ax))]
                fit_summary[(nuc, k)] = (max(errs), float(np.median(errs)))

    rows = []
    for nuc in NUCLIDES:
        p_ax, e_ax, tab = grids[nuc]
        lam = LAMBDA[nuc]
        for j, pr in enumerate(p_ax):
            p_mu = PMU0[nuc] * np.exp((1013.25 - pr) / L_MU_ATM[nuc])
            lmu_balco = 1.0 / (LMU_A + LMU_B * pr)
            for i, le in enumerate(e_ax):
                eps = 10.0 ** le
                n_tab = tab[i, j]
                n_exp = v3_exponential(nuc, pr, eps)
                denom = p_mu - lam * n_tab
                l_implied = eps * n_tab / denom if denom > 0 else float("nan")
                rows.append({
                    "nuclide": nuc,
                    "p_hPa": f"{pr:.6g}",
                    "log10_erate": f"{le:.17g}",
                    "n_mu_tabulated": f"{n_tab:.10g}",
                    "n_mu_v3_exponential": f"{n_exp:.10g}",
                    "ratio_exp_over_tab": f"{n_exp / n_tab:.10g}",
                    "lambda_mu_balco_g_cm2": f"{lmu_balco:.10g}",
                    "lambda_mu_implied_g_cm2": f"{l_implied:.10g}",
                })

    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    with Path(args.out).open("w", newline="") as f:
        f.write("# Balco CRONUS v3 single-exponential muon approximation: error report\n")
        f.write("# FROZEN reference - refreshing it is a reviewable event.\n")
        f.write(f"# generated_utc: {now}\n")
        f.write(f"# source_reference: {ref_path.as_posix()}\n")
        for line in meta:
            if line.startswith(("# source_path", "# source_sha256", "# Nmu.version")):
                f.write(f"# upstream{line[1:]}\n")
        f.write("#\n")
        f.write("# WHY THIS EXISTS: lsdtt3's balcoCRONUSv3 / balcoCRONUSv3PD schemes\n")
        f.write("# implemented the model in column n_mu_v3_exponential.  CRONUS v3's\n")
        f.write("# erosion-rate root-find (get_erates_v3.m:155 St, :268 LSDn) uses the\n")
        f.write("# TABULATED inventory instead; the exponential is used only for ages,\n")
        f.write("# saturation checks and plotting.  Those schemes were removed\n")
        f.write("# 2026-08-01.  ratio_exp_over_tab is how much muon inventory they\n")
        f.write("# actually delivered.\n")
        f.write("#\n")
        f.write("# lambda_mu_balco_g_cm2   = 1/(0.01036 - 9.697e-6 p), make_consts_v3.m:211\n")
        f.write("# lambda_mu_implied_g_cm2 = eps N /(P_mu - lambda N), i.e. the effective\n")
        f.write("#                           attenuation the tabulated inventory implies\n")
        if fit_summary:
            f.write("#\n")
            f.write("# K-exponential feasibility (max |rel error| % per pressure slice,\n")
            f.write("# fitting N = sum_k P_k L_k/(eps + L_k lambda), lsdtt3's kernel form):\n")
            for nuc in NUCLIDES:
                for k in (1, 2, 3):
                    worst, med = fit_summary[(nuc, k)]
                    f.write(f"#   {nuc:5s} K={k}: worst {worst:8.3f}%  median {med:8.3f}%\n")
            f.write("#   K=1 is what the removed schemes had; lsdtt3's Gamma_[1..3]\n")
            f.write("#   provide K=3 with no kernel change.\n")
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"WROTE: {args.out}  ({len(rows)} rows)")
    for nuc in NUCLIDES:
        p_ax, e_ax, tab = grids[nuc]
        r = [v3_exponential(nuc, pr, 10.0 ** le) / tab[i, j]
             for j, pr in enumerate(p_ax) for i, le in enumerate(e_ax)]
        print(f"  {nuc}: ratio_exp_over_tab min {min(r):.3f} max {max(r):.3f}")


if __name__ == "__main__":
    main()
