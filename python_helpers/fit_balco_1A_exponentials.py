#!/usr/bin/env python3
"""Fit a 3-exponential lsdtt-style muon approximation to Balco v3 tabulated 1A.

This is an offline analysis helper for the muon-scheme isolation work.
It reads the frozen `balco_v3_nmu_reference.csv`, fits shared attenuation lengths
with a shared muon split ratio, and writes a frozen fit report CSV.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from datetime import datetime, timezone

import numpy as np
from scipy.optimize import least_squares


NUCLIDES = ("Be10", "Al26", "C14")


@dataclass(frozen=True)
class NuclideConst:
    lambda_yr: float
    p0: float
    pmu0_slhl: float
    l_mu_atm: float


CONST = {
    "Be10": NuclideConst(
        lambda_yr=math.log(2.0) / 1.387e6,
        p0=4.086 + 0.07345,
        pmu0_slhl=0.07345,
        l_mu_atm=299.2,
    ),
    "Al26": NuclideConst(
        lambda_yr=9.83e-7,
        p0=28.535 + 0.67644,
        pmu0_slhl=0.67644,
        l_mu_atm=288.0,
    ),
    "C14": NuclideConst(
        lambda_yr=math.log(2.0) / 5730.0,
        p0=12.1 + 3.06690,
        pmu0_slhl=3.06690,
        l_mu_atm=267.8,
    ),
}


def _split_csv(line: str) -> list[str]:
    return [f.strip() for f in line.rstrip("\r\n").split(",")]


def _format_g(x: float) -> str:
    return format(float(x), ".17g")


class NmuTable:
    def __init__(self) -> None:
        self.p_axis: np.ndarray | None = None
        self.loge_axis: np.ndarray | None = None
        self.grids: dict[str, np.ndarray] = {}

    @staticmethod
    def from_csv(path: str) -> "NmuTable":
        tab = NmuTable()
        rows: list[tuple[str, float, float, float]] = []
        with open(path, newline="", encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                if line.startswith("#"):
                    if line.startswith("# p_hPa_axis_csv:"):
                        tab.p_axis = np.array(
                            [float(x) for x in _split_csv(line.split(":", 1)[1])],
                            dtype=np.float64,
                        )
                    elif line.startswith("# log10_erate_axis_csv:"):
                        tab.loge_axis = np.array(
                            [float(x) for x in _split_csv(line.split(":", 1)[1])],
                            dtype=np.float64,
                        )
                    continue
                if line.startswith("nuclide,"):
                    continue
                f = _split_csv(line)
                if len(f) != 4:
                    continue
                rows.append((f[0], float(f[1]), float(f[2]), float(f[3])))

        if tab.p_axis is None or tab.loge_axis is None:
            raise RuntimeError("missing axis comments in reference CSV")

        np_ = tab.p_axis.size
        ne_ = tab.loge_axis.size
        p0 = float(tab.p_axis[0])
        l0 = float(tab.loge_axis[0])
        dp = float((tab.p_axis[-1] - tab.p_axis[0]) / (np_ - 1))
        dl = float((tab.loge_axis[-1] - tab.loge_axis[0]) / (ne_ - 1))

        for nuc in NUCLIDES:
            tab.grids[nuc] = np.full((ne_, np_), np.nan, dtype=np.float64)

        for nuc, p, le, n in rows:
            if nuc not in tab.grids:
                continue
            ix = int(round((p - p0) / dp))
            iy = int(round((le - l0) / dl))
            if ix < 0 or ix >= np_ or iy < 0 or iy >= ne_:
                raise RuntimeError("row out of axis bounds")
            tab.grids[nuc][iy, ix] = n

        for nuc in NUCLIDES:
            if not np.isfinite(tab.grids[nuc]).all():
                raise RuntimeError(f"grid has holes for {nuc}")
        return tab

    def interp(self, nuc: str, p_hpa: float, eps: np.ndarray) -> np.ndarray:
        assert self.p_axis is not None and self.loge_axis is not None
        g = self.grids[nuc]
        loge = np.log10(eps)

        def frac(ax: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            step = (ax[-1] - ax[0]) / (ax.size - 1)
            t = (v - ax[0]) / step
            t = np.clip(t, 0.0, ax.size - 1.0)
            i = np.minimum(np.floor(t).astype(np.int64), ax.size - 2)
            f = t - i
            return i, f

        ix, fx = frac(self.p_axis, np.full_like(loge, p_hpa))
        iy, fy = frac(self.loge_axis, loge)
        g00 = g[iy, ix]
        g01 = g[iy, ix + 1]
        g10 = g[iy + 1, ix]
        g11 = g[iy + 1, ix + 1]
        return (1.0 - fy) * ((1.0 - fx) * g00 + fx * g01) + fy * ((1.0 - fx) * g10 + fx * g11)


def balco_one_exp(nuc: str, p_hpa: float, eps: np.ndarray) -> np.ndarray:
    c = CONST[nuc]
    p_mu = c.pmu0_slhl * np.exp((1013.25 - p_hpa) / c.l_mu_atm)
    lambda_mu = 1.0 / (0.01036 - 9.697e-6 * p_hpa)
    return p_mu / (c.lambda_yr + eps / lambda_mu)


def fit_model(nuc: str, p_hpa: float, eps: np.ndarray,
              logit_r: float, log_g1: float, log_dg: float) -> np.ndarray:
    c = CONST[nuc]
    fmu = c.pmu0_slhl / c.p0
    r = 1.0 / (1.0 + np.exp(-logit_r))
    f1 = fmu * r
    f2 = fmu - f1
    g1 = np.exp(log_g1)
    g2 = g1 + np.exp(log_dg)
    s1 = np.exp((1013.25 - p_hpa) / c.l_mu_atm)
    term1 = s1 * f1 * g1 / (eps + g1 * c.lambda_yr)
    term2 = f2 * g2 / (eps + g2 * c.lambda_yr)
    return c.p0 * (term1 + term2)


def rel_stats(pred: np.ndarray, ref: np.ndarray) -> tuple[float, float, float, float]:
    rel = pred / ref - 1.0
    return (
        float(np.min(rel)),
        float(np.max(rel)),
        float(np.mean(rel)),
        float(np.sqrt(np.mean(rel * rel))),
    )


def main() -> int:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--reference",
        default=os.path.join(repo_root, "tests", "expected", "cosmo", "balco_v3_nmu_reference.csv"),
    )
    ap.add_argument(
        "--output",
        default=os.path.join(repo_root, "tests", "expected", "cosmo", "balco_1A_exp_fit_report.csv"),
    )
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out = os.path.abspath(args.output)
    if os.path.exists(out) and not args.force:
        raise FileExistsError(f"{out} exists; pass --force")

    tab = NmuTable.from_csv(os.path.abspath(args.reference))
    eps = np.logspace(-5.0, -1.0, 200)
    p_fit = 1013.25

    def residuals(x: np.ndarray) -> np.ndarray:
        rr = []
        for nuc in NUCLIDES:
            ref = tab.interp(nuc, p_fit, eps)
            pred = fit_model(nuc, p_fit, eps, x[0], x[1], x[2])
            rr.append(pred / ref - 1.0)
        return np.concatenate(rr)

    x0 = np.array([0.0, np.log(1800.0), np.log(8000.0)], dtype=np.float64)
    lsq = least_squares(residuals, x0, method="trf")
    x = lsq.x
    r = 1.0 / (1.0 + np.exp(-x[0]))
    g1 = float(np.exp(x[1]))
    g2 = float(g1 + np.exp(x[2]))

    rows: list[list[str]] = []
    for p_hpa in (1013.25, 700.0):
        all_fit = []
        all_one = []
        for nuc in NUCLIDES:
            ref = tab.interp(nuc, p_hpa, eps)
            pred_fit = fit_model(nuc, p_hpa, eps, x[0], x[1], x[2])
            pred_one = balco_one_exp(nuc, p_hpa, eps)
            all_fit.append(pred_fit / ref - 1.0)
            all_one.append(pred_one / ref - 1.0)

            mn, mx, mean, rms = rel_stats(pred_fit, ref)
            rows.append([
                nuc,
                _format_g(p_hpa),
                "fit_3exp",
                _format_g(mn),
                _format_g(mx),
                _format_g(mean),
                _format_g(rms),
                _format_g(max(abs(mn), abs(mx))),
            ])

            mn, mx, mean, rms = rel_stats(pred_one, ref)
            rows.append([
                nuc,
                _format_g(p_hpa),
                "runB_oneexp",
                _format_g(mn),
                _format_g(mx),
                _format_g(mean),
                _format_g(rms),
                _format_g(max(abs(mn), abs(mx))),
            ])

        all_fit_vec = np.concatenate(all_fit)
        all_one_vec = np.concatenate(all_one)
        for model, vec in (("fit_3exp", all_fit_vec), ("runB_oneexp", all_one_vec)):
            mn = float(np.min(vec))
            mx = float(np.max(vec))
            mean = float(np.mean(vec))
            rms = float(np.sqrt(np.mean(vec * vec)))
            rows.append([
                "ALL",
                _format_g(p_hpa),
                model,
                _format_g(mn),
                _format_g(mx),
                _format_g(mean),
                _format_g(rms),
                _format_g(max(abs(mn), abs(mx))),
            ])

    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        fh.write("# Balco v3 1A three-exponential fit report for lsdtt3\n")
        fh.write("# FROZEN reference - refreshing it is a reviewable event.\n")
        fh.write(f"# generated_utc: {ts}\n")
        fh.write(f"# source_reference: {os.path.abspath(args.reference)}\n")
        fh.write("# fit_target: relative residual N_fit/N_tab - 1 at p=1013.25 over eps in [1e-5,1e-1], 200 log-spaced points\n")
        fh.write(f"# fit_param_r (muon split to Gamma1): {_format_g(r)}\n")
        fh.write(f"# fit_param_Gamma1_g_cm2: {_format_g(g1)}\n")
        fh.write(f"# fit_param_Gamma2_g_cm2: {_format_g(g2)}\n")
        fh.write("nuclide,p_hPa,model,min_rel,max_rel,mean_rel,rms_rel,max_abs_rel\n")
        w = csv.writer(fh)
        w.writerows(rows)

    print(f"wrote: {out}")
    print(f"r={r:.17g} Gamma1={g1:.17g} Gamma2={g2:.17g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
