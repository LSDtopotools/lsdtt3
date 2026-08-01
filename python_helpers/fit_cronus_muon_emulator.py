#!/usr/bin/env python3
"""Fit lsdtt3's three-exponential muon model to the CRONUS v3 tabulated
model-1A muon inventory.

Phase 1 of docs/plan_docs/cronus_to_lsdtt3_emulator_plan.md.

Model (lsdtt3's own kernel form, so no change to lsdtt_crn_particle):

    N_mu(eps, p) = sum_{k=1..3}  P_k(p) * L_k(p) / (eps + L_k(p) * lambda)

with smooth pressure dependence, x = (p - 1013.25)/100:

    log P_k(x) = a_k0 + a_k1 x + a_k2 x^2
    log D_k(x) = b_k0 + b_k1 x + b_k2 x^2
        L_1 = D_1;  L_2 = L_1 + D_2;  L_3 = L_2 + D_3

The cumulative form for L_k guarantees L_1 < L_2 < L_3 at EVERY pressure, so the
pathway labelling stays continuous in p and the frozen coefficients are
diffable.  18 coefficients per nuclide.

Fitting is done on log(N), so relative error is minimised uniformly across the
five decades of erosion rate the table spans.

Two variants are fitted and compared, resolving the open question in plan
section 4.3:
  * free     -- no constraint on the surface production
  * anchored -- soft penalty pulling sum_k P_k(p) onto Balco's
                P_mu,0(p) = Pmu0 exp((1013.25 - p)/L_mu_atm), which section 2.3
                measured as accurate to 1-3%

Reads the FROZEN extraction, not the GPL reference tree.

Usage:
  python3 python_helpers/fit_cronus_muon_emulator.py
"""

from __future__ import annotations

import argparse
import csv
import datetime
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

# make_consts_v3.m:208-211 (index 4 = Be-10 quartz, 5 = C-14, 7 = Al-26)
PMU0 = {"Be10": 0.07345, "Al26": 0.67644, "C14": 3.06690}
L_MU_ATM = {"Be10": 299.2, "Al26": 288.0, "C14": 267.8}

# make_consts_v3.m:134,142,147
LAMBDA = {
    "Be10": math.log(2.0) / 1.387e6,
    "Al26": 9.83e-7,
    "C14": math.log(2.0) / 5730.0,
}

NUCLIDES = ("Be10", "Al26", "C14")
N_EXP = 3
P_REF = 1013.25
P_SCALE = 100.0

# Plan section 5, Phase 1 exit criterion.
EXIT_MAX_REL_PCT = {"Be10": 1.0, "Al26": 2.0, "C14": 2.0}


def read_reference(path: Path):
    meta = [l.rstrip("\n") for l in path.open() if l.startswith("#")]
    recs = list(csv.DictReader(l for l in path.open() if not l.startswith("#")))
    out = {}
    for nuc in NUCLIDES:
        sub = [r for r in recs if r["nuclide"] == nuc]
        p_ax = np.array(sorted({float(r["p_hPa"]) for r in sub}))
        e_ax = np.array(sorted({float(r["log10_erate"]) for r in sub}))
        tab = np.full((len(e_ax), len(p_ax)), np.nan)
        pi = {v: i for i, v in enumerate(p_ax)}
        ei = {v: i for i, v in enumerate(e_ax)}
        for r in sub:
            tab[ei[float(r["log10_erate"])], pi[float(r["p_hPa"])]] = \
                float(r["N_mu_atoms_g"])
        if np.isnan(tab).any():
            raise ValueError(f"{nuc}: reference grid has holes")
        out[nuc] = (p_ax, e_ax, tab)
    return meta, out


def unpack(theta):
    """theta -> (a[3,3], b[3,3]) coefficient blocks."""
    a = theta[: N_EXP * 3].reshape(N_EXP, 3)
    b = theta[N_EXP * 3:].reshape(N_EXP, 3)
    return a, b


def eval_pathways(theta, p):
    """Return P_k(p) and L_k(p), each shape (3, len(p)), ordered L_1<L_2<L_3."""
    a, b = unpack(theta)
    x = (np.asarray(p, float) - P_REF) / P_SCALE
    powers = np.vstack([np.ones_like(x), x, x * x])          # (3, np)

    with np.errstate(over="ignore"):
        P = np.exp(a @ powers)                                # (3, np)
        inc = np.exp(b @ powers)                              # (3, np)
    L = np.empty_like(inc)
    L[0] = inc[0]
    L[1] = L[0] + inc[1]
    L[2] = L[1] + inc[2]
    return P, L


def model_inventory(theta, p_ax, erates, lam):
    """N_mu on the (erate, pressure) grid -> shape (n_e, n_p)."""
    P, L = eval_pathways(theta, p_ax)
    e = erates[:, None]
    out = np.zeros((len(erates), len(p_ax)))
    # Wild seeds can overflow; least_squares rejects the resulting inf/nan cost
    # on its own, so the warnings are noise.
    with np.errstate(over="ignore", invalid="ignore"):
        for k in range(N_EXP):
            out += P[k][None, :] * L[k][None, :] / (e + L[k][None, :] * lam)
    return out


def surface_production(theta, p_ax):
    P, _ = eval_pathways(theta, p_ax)
    return P.sum(axis=0)


def balco_surface_production(nuc, p_ax):
    return PMU0[nuc] * np.exp((P_REF - np.asarray(p_ax, float)) / L_MU_ATM[nuc])


def table_surface_vs_balco(nuc, p_ax, e_ax, tab):
    """How far the TABLE's own eps->0 surface production sits from Balco's
    P_mu,0(p) formula.  This is the number that decides the anchor question: if
    the table disagrees with the formula by this much, anchoring the fit to the
    formula can only pull it away from the ground truth."""
    n_lo = tab[0, :]                       # smallest erosion rate in the grid
    implied = n_lo * LAMBDA[nuc]
    return float(np.max(np.abs(implied / balco_surface_production(nuc, p_ax) - 1.0)) * 100.0)


def fit_nuclide(nuc, p_ax, e_ax, tab, anchor_weight, seeds=40):
    lam = LAMBDA[nuc]
    erates = 10.0 ** e_ax
    log_tab = np.log(tab)
    p_anchor = balco_surface_production(nuc, p_ax)

    def residual(theta):
        model = model_inventory(theta, p_ax, erates, lam)
        r = (np.log(model) - log_tab).ravel()
        if anchor_weight > 0.0:
            s = surface_production(theta, p_ax)
            r = np.concatenate([r, anchor_weight * (np.log(s) - np.log(p_anchor))])
        return r

    # Seed: split Balco's SLHL surface production across three pathways with
    # attenuation lengths spanning the range the table implies (roughly
    # 300-40000 g/cm2), flat in pressure, then perturb.
    best = None
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        a0 = np.zeros((N_EXP, 3))
        b0 = np.zeros((N_EXP, 3))
        frac = np.array([0.75, 0.22, 0.03])
        for k in range(N_EXP):
            a0[k, 0] = np.log(PMU0[nuc] * frac[k])
        b0[0, 0] = np.log(1200.0)
        b0[1, 0] = np.log(4000.0)     # increment, not absolute
        b0[2, 0] = np.log(20000.0)
        theta0 = np.concatenate([a0.ravel(), b0.ravel()])
        if seed > 0:
            theta0 = theta0 + rng.normal(0.0, 0.45, theta0.shape)
        try:
            r = least_squares(residual, theta0, method="lm", max_nfev=200000)
        except Exception:  # noqa: BLE001 - bad seed, try the next one
            continue
        if best is None or r.cost < best.cost:
            best = r

    if best is None:
        raise RuntimeError(f"{nuc}: all seeds failed")

    theta = best.x
    model = model_inventory(theta, p_ax, erates, lam)
    rel = model / tab - 1.0
    return theta, rel


def stats(rel, p_ax=None, e_ax=None):
    out = {
        "max_abs_pct": float(np.max(np.abs(rel)) * 100.0),
        "rms_pct": float(np.sqrt(np.mean(rel ** 2)) * 100.0),
        "mean_pct": float(np.mean(rel) * 100.0),
        "min_pct": float(np.min(rel) * 100.0),
        "max_pct": float(np.max(rel) * 100.0),
    }
    if p_ax is not None:
        i, j = np.unravel_index(np.argmax(np.abs(rel)), rel.shape)
        out["worst_p_hPa"] = float(p_ax[j])
        out["worst_log10_erate"] = float(e_ax[i])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference",
                    default="tests/expected/cosmo/balco_v3_nmu_reference.csv")
    ap.add_argument("--out-coeffs",
                    default="tests/expected/cosmo/cronus_muon_emulator_coeffs.csv")
    ap.add_argument("--out-report",
                    default="tests/expected/cosmo/cronus_muon_emulator_fit_report.csv")
    ap.add_argument("--anchor-weight", type=float, default=3.0,
                    help="soft-anchor weight on sum_k P_k(p) vs Balco P_mu,0(p)")
    ap.add_argument("--seeds", type=int, default=40)
    args = ap.parse_args()

    ref_path = Path(args.reference)
    meta, grids = read_reference(ref_path)

    results = {}
    for nuc in NUCLIDES:
        p_ax, e_ax, tab = grids[nuc]
        for variant, w in (("free", 0.0), ("anchored", args.anchor_weight)):
            theta, rel = fit_nuclide(nuc, p_ax, e_ax, tab, w, seeds=args.seeds)
            s = stats(rel, p_ax, e_ax)
            surf = surface_production(theta, p_ax)
            anchor = balco_surface_production(nuc, p_ax)
            s["surface_vs_balco_max_pct"] = float(
                np.max(np.abs(surf / anchor - 1.0)) * 100.0)
            results[(nuc, variant)] = (theta, rel, s)
            print(f"{nuc:5s} {variant:9s}  max|rel| {s['max_abs_pct']:6.3f}%   "
                  f"rms {s['rms_pct']:6.3f}%   "
                  f"surface vs Balco max {s['surface_vs_balco_max_pct']:6.3f}%")

    # ANCHOR QUESTION RESOLVED (plan section 4.3, measured 2026-08-01): reject
    # the anchor, freeze 'free' for every nuclide.
    #
    # The tabulated inventory is the ground truth, and the table's OWN eps->0
    # surface production already disagrees with Balco's P_mu,0(p) formula by up
    # to 3.00% (Be10), 2.19% (Al26), 0.55% (C14) -- see the
    # table_surface_vs_balco_max_pct column.  The free fit's deviation from the
    # formula tracks those figures almost exactly, i.e. it is reproducing the
    # table, not drifting.  Anchoring to the formula therefore pulls the fit
    # AWAY from ground truth, and measurably so: Be10 max|rel| degrades from
    # 0.77% to 1.57%, straight through the 1% gate.
    chosen = {nuc: "free" for nuc in NUCLIDES}

    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")

    def provenance(f):
        f.write("# FROZEN reference - refreshing it is a reviewable event.\n")
        f.write(f"# generated_utc: {now}\n")
        f.write("# generator: python_helpers/fit_cronus_muon_emulator.py\n")
        f.write(f"# source_reference: {ref_path.as_posix()}\n")
        for line in meta:
            if line.startswith(("# source_path", "# source_sha256", "# Nmu.version")):
                f.write(f"# upstream{line[1:]}\n")

    with Path(args.out_coeffs).open("w", newline="") as f:
        f.write("# lsdtt3 CRONUS v3 muon emulator coefficients\n")
        provenance(f)
        f.write("#\n")
        f.write("# Model, with x = (p_hPa - 1013.25)/100:\n")
        f.write("#   log P_k(x) = a_k0 + a_k1 x + a_k2 x^2        [atoms g^-1 yr^-1]\n")
        f.write("#   log D_k(x) = b_k0 + b_k1 x + b_k2 x^2        [g cm^-2]\n")
        f.write("#   L_1 = D_1;  L_2 = L_1 + D_2;  L_3 = L_2 + D_3   (so L_1<L_2<L_3)\n")
        f.write("#   N_mu(eps,p) = sum_k P_k(p) L_k(p) / (eps + L_k(p) lambda)\n")
        f.write("#\n")
        f.write("# Clamp p to [450, 1020] hPa before evaluating: CRONUS v3's intNmu.m\n")
        f.write("# clamps to its grid bounds, and these polynomials must not be\n")
        f.write("# extrapolated beyond the fitted range.\n")
        f.write("#\n")
        f.write("# variant: 'anchored' = sum_k P_k(p) softly tied to Balco's\n")
        f.write("#   P_mu,0(p) = Pmu0 exp((1013.25-p)/L_mu_atm); 'free' = unconstrained.\n")
        for nuc in NUCLIDES:
            f.write(f"# chosen_variant[{nuc}]: {chosen[nuc]}\n")
        w = csv.writer(f)
        w.writerow(["nuclide", "variant", "chosen", "lambda_yr", "pathway",
                    "a0", "a1", "a2", "b0", "b1", "b2"])
        for nuc in NUCLIDES:
            for variant in ("free", "anchored"):
                theta = results[(nuc, variant)][0]
                a, b = unpack(theta)
                for k in range(N_EXP):
                    # %.17g round-trips a double exactly and stays plain ASCII;
                    # repr() on a numpy scalar emits "np.float64(...)", which
                    # would make this frozen file unparseable from C++.
                    g = lambda v: f"{float(v):.17g}"
                    w.writerow([nuc, variant,
                                "yes" if chosen[nuc] == variant else "no",
                                g(LAMBDA[nuc]), k + 1,
                                g(a[k, 0]), g(a[k, 1]), g(a[k, 2]),
                                g(b[k, 0]), g(b[k, 1]), g(b[k, 2])])
    print(f"WROTE: {args.out_coeffs}")

    with Path(args.out_report).open("w", newline="") as f:
        f.write("# lsdtt3 CRONUS v3 muon emulator: fit quality\n")
        provenance(f)
        f.write("#\n")
        f.write("# Plan Phase 1 exit criterion: max |rel| < 1% (Be10), < 2% (Al26, C14).\n")
        f.write("# rel = N_fit/N_tab - 1 over the full 20 pressure x 28 erosion-rate grid.\n")
        f.write("#\n")
        f.write("# ANCHOR QUESTION (plan 4.3) RESOLVED: 'free' is frozen for every\n")
        f.write("# nuclide.  Compare fit_surface_vs_balco_max_pct against\n")
        f.write("# table_surface_vs_balco_max_pct: the TABLE itself already departs from\n")
        f.write("# Balco's P_mu,0(p) formula by that much, and the free fit departs by the\n")
        f.write("# same amount, i.e. it is tracking ground truth.  Anchoring to the formula\n")
        f.write("# costs real accuracy (Be10 0.77% -> 1.57%, through the gate) to chase a\n")
        f.write("# quantity the reference does not itself honour.\n")
        w = csv.writer(f)
        w.writerow(["nuclide", "variant", "chosen", "n_points",
                    "max_abs_pct", "rms_pct", "mean_pct", "min_pct", "max_pct",
                    "worst_p_hPa", "worst_log10_erate",
                    "fit_surface_vs_balco_max_pct",
                    "table_surface_vs_balco_max_pct",
                    "exit_gate_pct", "passes_gate"])
        for nuc in NUCLIDES:
            p_ax, e_ax, tab = grids[nuc]
            for variant in ("free", "anchored"):
                s = results[(nuc, variant)][2]
                gate = EXIT_MAX_REL_PCT[nuc]
                w.writerow([nuc, variant, "yes" if chosen[nuc] == variant else "no",
                            tab.size,
                            f"{s['max_abs_pct']:.6f}", f"{s['rms_pct']:.6f}",
                            f"{s['mean_pct']:.6f}", f"{s['min_pct']:.6f}",
                            f"{s['max_pct']:.6f}",
                            f"{s['worst_p_hPa']:.6g}",
                            f"{s['worst_log10_erate']:.6g}",
                            f"{s['surface_vs_balco_max_pct']:.6f}",
                            f"{table_surface_vs_balco(nuc, p_ax, e_ax, tab):.6f}",
                            f"{gate:.1f}",
                            "yes" if s["max_abs_pct"] < gate else "no"])
    print(f"WROTE: {args.out_report}")

    print("\nPhase 1 exit criterion:")
    all_pass = True
    for nuc in NUCLIDES:
        s = results[(nuc, chosen[nuc])][2]
        gate = EXIT_MAX_REL_PCT[nuc]
        ok = s["max_abs_pct"] < gate
        all_pass &= ok
        print(f"  {nuc:5s} [{chosen[nuc]:8s}] max|rel| = {s['max_abs_pct']:6.3f}% "
              f"vs gate {gate:.1f}%  -> {'PASS' if ok else 'FAIL'}")
    print("PHASE 1: " + ("PASS" if all_pass else "FAIL"))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
