#!/usr/bin/env python3
"""Functional coverage for lsdtt-cosmo-column-tool.

Until now this CLI had exactly one test: `-h`. It was the only cosmo tool with
no functional coverage at all, which matters more than it looks, because it is
the natural place to compare lsdtt3 against point calculators (CRONUS v3,
CRONUScalc): no DEM, no flow routing, no basin averaging, so nothing can hide a
discrepancy in the scaling or the muon terms.

The checks below need NO external reference. They are LIMIT IDENTITIES of the
column model itself -- cases where the answer is forced by the physics, so a
wrong answer is unambiguous:

  * step_change at t = 0 must equal steady state at the OLD rate;
  * step_change at t -> infinity must equal steady state at the NEW rate;
  * a full_transient run at a CONSTANT rate must relax to steady state at that
    rate;
  * concentration must fall monotonically with depth, and with erosion rate;
  * production_scaling and total_shielding must be INTERCHANGEABLE, because as
    of 2026-07-28 both go through the Vermeesch redistribution as a single
    product, matching the raster pipeline; and neither may be a flat P0 multiply,
    since the redistribution attenuates the pathways unequally;
  * C-14 must re-equilibrate far faster than Be-10 after a step change (the
    short half-life is the whole reason C-14 is interesting).

Sign conventions and magnitudes are pinned elsewhere (the sentinel, the CRONUS
comparison). This file pins that the CLI *plumbs its parameters into the engine
correctly* -- which no library test can do.

Exit codes: 0 pass, 1 fail, 77 skip (missing executable).
"""

import argparse
import csv
import math
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

EXIT_OK, EXIT_FAIL, EXIT_SKIP = 0, 1, 77

n_pass = n_fail = 0


def check(name, ok, detail=""):
    global n_pass, n_fail
    if ok:
        n_pass += 1
    else:
        n_fail += 1
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


def rel(a, b):
    return abs(a - b) / abs(b) if b else (0.0 if a == b else math.inf)


def run(exe, workdir, params, tag, numeric=True):
    """Run one column-tool case; return CSV rows as dicts."""
    param = workdir / f"{tag}.param"
    lines = [f"write_path: {workdir}", f"output_file: {tag}.csv"]
    lines += [f"{k}: {v}" for k, v in params.items()]
    param.write_text("\n".join(lines) + "\n")
    proc = subprocess.run([str(exe), str(workdir), param.name],
                          cwd=str(workdir), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"{tag} failed (exit {proc.returncode})\n"
                           f"--- stderr ---\n{proc.stderr[-2000:]}")
    out = workdir / f"{tag}.csv"
    if not out.is_file():
        raise RuntimeError(f"{tag} wrote no CSV: {out}")
    rows = list(csv.DictReader(out.open(newline="")))
    if not rows:
        raise RuntimeError(f"{tag} wrote an empty CSV")
    if numeric:
        rows = [{k: float(v) for k, v in r.items()} for r in rows]
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=os.environ.get("LSDTT_COSMO_COLUMN_TOOL"))
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    if not args.exe or not pathlib.Path(args.exe).is_file():
        print(f"SKIP: lsdtt-cosmo-column-tool not found ({args.exe})")
        return EXIT_SKIP
    exe = pathlib.Path(args.exe)

    workdir = pathlib.Path(tempfile.mkdtemp(prefix="cosmo_column_"))
    try:
        E_OLD, E_NEW = 0.002, 0.02

        def ss(erate, depth=0.0, scaling=1.0, shield=1.0, tag=None):
            return run(exe, workdir, {
                "column_mode": "steady_state",
                "effective_depth_g_cm2": depth,
                "effective_erosion_rate_g_cm2_yr": erate,
                "production_scaling": scaling,
                "total_shielding": shield,
            }, tag or f"ss_{erate}_{depth}_{scaling}_{shield}")[0]

        def step(t, depth=0.0, shield=1.0, scaling=1.0, tag=None):
            return run(exe, workdir, {
                "column_mode": "step_change",
                "effective_depth_g_cm2": depth,
                "effective_erosion_rate_old_g_cm2_yr": E_OLD,
                "effective_erosion_rate_g_cm2_yr": E_NEW,
                "time_since_step_yr": t,
                "total_shielding": shield,
                "production_scaling": scaling,
            }, tag or f"step_{t}_{depth}_{shield}_{scaling}")[0]

        def burial_pair(n_fast, n_slow, shield=1.0,
                        fast="Al26", slow="Be10", tag=None):
            return run(exe, workdir, {
                "column_mode": "burial_pair",
                "nuclide_fast": fast,
                "nuclide_slow": slow,
                "conc_fast_atoms_g": n_fast,
                "conc_slow_atoms_g": n_slow,
                "total_shielding": shield,
            }, tag or f"burial_{fast}_{slow}_{n_fast}_{n_slow}_{shield}",
                numeric=False)[0]

        def step_pair(n_fast, n_slow, e_old=E_OLD,
                      fast="C14", slow="Be10", shield=1.0, tag=None):
            return run(exe, workdir, {
                "column_mode": "step_pair",
                "nuclide_fast": fast,
                "nuclide_slow": slow,
                "conc_fast_atoms_g": n_fast,
                "conc_slow_atoms_g": n_slow,
                "effective_erosion_rate_old_g_cm2_yr": e_old,
                "total_shielding": shield,
            }, tag or f"step_pair_{fast}_{slow}_{n_fast}_{n_slow}_{shield}",
                numeric=False)[0]

        # --- 1. every mode writes a machine-readable CSV ---------------------
        base = ss(E_NEW, tag="base")
        check("steady_state writes a CSV with the documented columns",
              set(base) >= {"depth_m", "effective_depth_g_cm2",
                            "C_10Be", "C_26Al", "C_14C"})
        check("steady_state concentrations are positive",
              all(base[c] > 0 for c in ("C_10Be", "C_26Al", "C_14C")),
              str({c: base[c] for c in ("C_10Be", "C_26Al", "C_14C")}))

        # --- 2. step_change limits ------------------------------------------
        # These are the sharpest checks here: at t=0 nothing has happened yet,
        # so the answer MUST be steady state at the old rate; after long enough
        # the memory of the old rate is gone. A CLI that dropped, swapped or
        # ignored either rate fails one of the two.
        ss_old, ss_new = ss(E_OLD, tag="ss_old"), ss(E_NEW, tag="ss_new")
        s0 = step(0.0, tag="step_t0")
        for nuc in ("C_10Be", "C_26Al", "C_14C"):
            check(f"step_change at t=0 == steady state at the OLD rate ({nuc})",
                  rel(s0[nuc], ss_old[nuc]) < 1e-9,
                  f"{s0[nuc]:.10g} vs {ss_old[nuc]:.10g}")
        s_inf = step(3.0e6, tag="step_tinf")
        for nuc in ("C_10Be", "C_26Al", "C_14C"):
            check(f"step_change at t=3 Myr -> steady state at the NEW rate ({nuc})",
                  rel(s_inf[nuc], ss_new[nuc]) < 1e-6,
                  f"{s_inf[nuc]:.10g} vs {ss_new[nuc]:.10g}")

        # The old and new states must be genuinely different, or the two checks
        # above are trivially satisfied and prove nothing.
        check("the step is a real step (old and new steady states differ)",
              rel(ss_old["C_10Be"], ss_new["C_10Be"]) > 0.5,
              f"old {ss_old['C_10Be']:.6g} vs new {ss_new['C_10Be']:.6g}")

        # --- 3. C-14 relaxes far faster than Be-10 --------------------------
        # The physical reason C-14 is used for transients. Measured at a time
        # where the contrast is largest; asserted as an ORDERING plus a wide
        # factor, not a pinned number.
        s_mid = step(6000.0, tag="step_mid")
        f_be = ((s_mid["C_10Be"] - ss_old["C_10Be"])
                / (ss_new["C_10Be"] - ss_old["C_10Be"]))
        f_c14 = ((s_mid["C_14C"] - ss_old["C_14C"])
                 / (ss_new["C_14C"] - ss_old["C_14C"]))
        check("C-14 re-equilibrates faster than Be-10 after a step change",
              f_c14 > f_be, f"C-14 {f_c14:.3f} vs Be-10 {f_be:.3f} of the way")
        # Measured 1.62 at 6 kyr (C-14 76.0 % of the way, Be-10 46.9 %); frozen
        # at 1.3, the house convention. NOT a predicted value -- measured first,
        # then backed off, so it guards the ordering without pinning noise.
        check("C-14 is >= 1.3x further relaxed than Be-10 at 6 kyr",
              f_c14 > 1.3 * f_be, f"ratio {f_c14 / f_be:.2f}")

        # --- 3b. THE SAME LIMIT IDENTITY, UNDER SHIELDING AND SCALING -------
        # This is the check that was missing, and its absence let a real defect
        # sit undetected: the two point modes applied total_shielding by
        # DIFFERENT mechanisms (steady_state through the Vermeesch
        # redistribution, step_change as a flat P0 multiply), so they disagreed
        # by 11.4 % at S = 0.5 while agreeing perfectly at S = 1. Asserting the
        # identity only at S = 1 -- as this file originally did -- cannot see
        # that class of fault at all. Fixed 2026-07-27 by giving the column an
        # additive overload that takes the two factors separately.
        for shield, scaling in ((0.5, 1.0), (1.0, 0.5), (0.7, 1.3), (0.5, 0.5)):
            sg = step(0.0, shield=shield, scaling=scaling,
                      tag=f"step_t0_s{shield}_p{scaling}")
            sr = ss(E_OLD, shield=shield, scaling=scaling,
                    tag=f"ss_old_s{shield}_p{scaling}")
            for nuc in ("C_10Be", "C_26Al", "C_14C"):
                check(f"step_change t=0 == steady state OLD rate at "
                      f"shielding {shield}, scaling {scaling} ({nuc})",
                      rel(sg[nuc], sr[nuc]) < 1e-9,
                      f"{sg[nuc]:.10g} vs {sr[nuc]:.10g}")

        # --- 4. monotonicity in depth and in erosion rate -------------------
        depths = [0.0, 50.0, 150.0, 400.0]
        cbd = [ss(E_NEW, depth=d, tag=f"depth_{int(d)}")["C_10Be"] for d in depths]
        check("concentration falls monotonically with depth",
              all(x > y for x, y in zip(cbd, cbd[1:])),
              " > ".join(f"{c:.4g}" for c in cbd))
        rates = [0.0005, 0.002, 0.01, 0.05]
        cbr = [ss(e, tag=f"rate_{e}")["C_10Be"] for e in rates]
        check("concentration falls monotonically with erosion rate",
              all(x > y for x, y in zip(cbr, cbr[1:])),
              " > ".join(f"{c:.4g}" for c in cbr))

        # --- 5. production_scaling and total_shielding are DIFFERENT ---------
        # They are deliberately not the same operation (lsdtt_crn_column.cpp
        # :126-128): production_scaling multiplies P0 uniformly and may exceed 1
        # (high-altitude enhancement), while total_shielding goes through the
        # Vermeesch virtual-attenuation redistribution in scale_F_values, which
        # attenuates the pathways unequally. So S=0.5 must NOT equal
        # scaling=0.5 -- asserting equality here would be asserting a bug.
        # Both must still MOVE the answer, or the knob is not wired up.
        a = ss(E_NEW, scaling=1.0, shield=0.5, tag="prod_shield")
        b = ss(E_NEW, scaling=0.5, shield=1.0, tag="prod_scaled")
        check("production_scaling and total_shielding are INTERCHANGEABLE "
              "(one convention: only their product reaches the physics)",
              rel(a["C_10Be"], b["C_10Be"]) < 1e-9,
              f"S=0.5 {a['C_10Be']:.10g} vs scaling=0.5 {b['C_10Be']:.10g}")
        check("halving the production scaling is NOT a no-op",
              rel(b["C_10Be"], base["C_10Be"]) > 0.1,
              f"{b['C_10Be']:.6g} vs {base['C_10Be']:.6g}")
        check("halving the shielding is NOT a no-op",
              rel(a["C_10Be"], base["C_10Be"]) > 0.1,
              f"{a['C_10Be']:.6g} vs {base['C_10Be']:.6g}")
        # Both knobs go through the Vermeesch virtual-attenuation
        # redistribution, which attenuates spallation and the muon pathways
        # UNEQUALLY -- so halving either must NOT halve the concentration. A
        # flat P0 multiply would be exactly linear, and that is precisely the
        # convention removed on 2026-07-28 to match the raster pipeline.
        check("halving the scaling is NOT a flat P0 multiply "
              "(Vermeesch redistribution, so not exactly linear)",
              rel(b["C_10Be"], 0.5 * base["C_10Be"]) > 0.05,
              f"{b['C_10Be']:.10g} vs a flat half {0.5 * base['C_10Be']:.10g}")

        # --- 6. full_transient at a constant rate relaxes to steady state ---
        # Drives the third mode, its history-file reader, and the multi-sample
        # profile writer. A long constant history must reproduce the analytic
        # steady state the first mode computes -- two independent code paths
        # meeting at the same number.
        hist = workdir / "const_history.csv"
        density, erate_m_yr = 2650.0, E_NEW / (2650.0 * 0.1)
        # Rows run from the MOST RECENT (index 0) to the oldest, with a final
        # row at time = -9999 carrying the long-term background rate. This was
        # originally written 2000000/0/-9999, i.e. reversed; it ran and produced
        # negative timesteps, and only passed because the rate was constant. The
        # validator rejects that ordering as of 2026-07-28.
        hist.write_text("time,rate,removal\n"
                        + f"2000000,{erate_m_yr},0\n"
                        + f"-9999,{erate_m_yr},0\n")
        prof = run(exe, workdir, {
            "column_mode": "full_transient",
            "transient_history_file": str(hist),
            "n_samples": 4,
            "sample_spacing_m": 0.5,
            "production_scaling": 1.0,
            "rock_density_kg_m3": density,
        }, "transient_const")
        check("full_transient writes one row per sample", len(prof) == 4,
              f"got {len(prof)}")
        surf = min(prof, key=lambda r: r["effective_depth_g_cm2"])
        ss_at_surf = ss(E_NEW, depth=surf["effective_depth_g_cm2"],
                        tag="ss_transient_cmp")
        check("full_transient at a constant rate relaxes to steady state",
              rel(surf["C_10Be"], ss_at_surf["C_10Be"]) < 1e-3,
              f"{surf['C_10Be']:.8g} vs {ss_at_surf['C_10Be']:.8g}")

        # --- 7. two-nuclide burial diagnostics are exposed on the CLI --------
        # Build a synthetic buried sample from a known pre-burial steady state,
        # then require the CLI inversion to recover both the rate and age.
        E_TRUE, T_BURIAL = 0.005, 1.0e6
        L10, L26 = 500.0e-9, 980.0e-9
        ss_true = ss(E_TRUE, tag="burial_ss_true")
        n10 = ss_true["C_10Be"] * math.exp(-L10 * T_BURIAL)
        n26 = ss_true["C_26Al"] * math.exp(-L26 * T_BURIAL)
        b = burial_pair(n26, n10, fast="Al26", slow="Be10",
                        tag="burial_pair_roundtrip")
        check("burial_pair writes and surfaces solve status",
              b.get("status") == "converged" and "status_code" in b,
              f"status={b.get('status')} code={b.get('status_code')}")
        check("burial_pair recovers pre-burial erosion rate",
              rel(float(b["erate_preburial_g_cm2_yr"]), E_TRUE) < 1e-5,
              f"{b['erate_preburial_g_cm2_yr']} vs {E_TRUE}")
        check("burial_pair recovers burial age",
              rel(float(b["burial_age_yr"]), T_BURIAL) < 1e-5,
              f"{b['burial_age_yr']} vs {T_BURIAL}")
        check("burial_pair reports burial_required when ratio is sub-SS",
              int(float(b["burial_required"])) == 1,
              f"burial_required={b['burial_required']}")

        # Uncertainty fields are present and positive when sigma inputs are set.
        b_unc = run(exe, workdir, {
            "column_mode": "burial_pair",
            "nuclide_fast": "Al26",
            "nuclide_slow": "Be10",
            "conc_fast_atoms_g": n26,
            "conc_slow_atoms_g": n10,
            "conc_fast_sigma_atoms_g": 0.02 * n26,
            "conc_slow_sigma_atoms_g": 0.02 * n10,
            "propagate_production_uncertainty": "true",
            "production_rel_sigma": 0.05,
        }, "burial_pair_unc", numeric=False)[0]
        check("burial_pair writes internal/external uncertainty fields",
              "erate_preburial_unc_int_g_cm2_yr" in b_unc
              and "erate_preburial_unc_ext_g_cm2_yr" in b_unc,
              str(sorted(b_unc.keys())[:8]))
        check("burial_pair uncertainty is positive when sigmas are provided",
              float(b_unc["erate_preburial_unc_int_g_cm2_yr"]) > 0.0
              and float(b_unc["erate_preburial_unc_ext_g_cm2_yr"]) > 0.0,
              "int=" + b_unc["erate_preburial_unc_int_g_cm2_yr"]
              + " ext=" + b_unc["erate_preburial_unc_ext_g_cm2_yr"])

        # Same calculation with a C14/Be10 pair: this is the new capability.
        T_BURIAL_C14 = 5000.0
        L14 = math.log(2.0) / 5730.0
        n10_c14 = ss_true["C_10Be"] * math.exp(-L10 * T_BURIAL_C14)
        n14 = ss_true["C_14C"] * math.exp(-L14 * T_BURIAL_C14)
        b14 = burial_pair(n14, n10_c14, fast="C14", slow="Be10",
                          tag="burial_pair_roundtrip_c14")
        check("burial_pair supports C14 as the fast nuclide",
              b14.get("status") == "converged"
              and rel(float(b14["erate_preburial_g_cm2_yr"]), E_TRUE) < 1e-4,
              "status=" + str(b14.get("status"))
              + " erate=" + b14["erate_preburial_g_cm2_yr"]
              + " age=" + b14["burial_age_yr"])

        # --- 8. two-nuclide step-pair inversion is exposed on the CLI --------
        # Use a forward modelled C14/Be10 pair from step_change and recover the
        # imposed (new rate, time) in step_pair.
        T_TRUE = 8000.0
        step_true = step(T_TRUE, tag="step_pair_forward")
        p = step_pair(step_true["C_14C"], step_true["C_10Be"],
                      tag="step_pair_roundtrip")
        check("step_pair writes and surfaces solve status",
              p.get("status") == "converged" and "status_code" in p,
              f"status={p.get('status')} code={p.get('status_code')}")
        check("step_pair recovers new erosion rate",
              rel(float(p["erate_new_g_cm2_yr"]), E_NEW) < 1e-3,
              f"{p['erate_new_g_cm2_yr']} vs {E_NEW}")
        check("step_pair recovers time since change",
              rel(float(p["time_since_change_yr"]), T_TRUE) < 1e-3,
              f"{p['time_since_change_yr']} vs {T_TRUE}")

        p_unc = run(exe, workdir, {
            "column_mode": "step_pair",
            "nuclide_fast": "C14",
            "nuclide_slow": "Be10",
            "conc_fast_atoms_g": step_true["C_14C"],
            "conc_slow_atoms_g": step_true["C_10Be"],
            "effective_erosion_rate_old_g_cm2_yr": E_OLD,
            "conc_fast_sigma_atoms_g": 0.02 * step_true["C_14C"],
            "conc_slow_sigma_atoms_g": 0.02 * step_true["C_10Be"],
        }, "step_pair_unc", numeric=False)[0]
        check("step_pair writes erate/time uncertainty fields",
              "erate_new_unc_int_g_cm2_yr" in p_unc
              and "time_since_change_unc_int_yr" in p_unc,
              str(sorted(p_unc.keys())[:8]))
        check("step_pair uncertainty is positive when sigmas are provided",
              float(p_unc["erate_new_unc_int_g_cm2_yr"]) > 0.0
              and float(p_unc["time_since_change_unc_int_yr"]) > 0.0,
              "e=" + p_unc["erate_new_unc_int_g_cm2_yr"]
              + " t=" + p_unc["time_since_change_unc_int_yr"])

        # --- 9. non-identifiable pairs are reported as status, not hidden -----
        # Al26/Be10 cannot uniquely date this short transient; status must be
        # not unique, with >1 brackets reported.
        step_amb = step(6000.0, tag="step_pair_amb_forward")
        p_amb = step_pair(step_amb["C_26Al"], step_amb["C_10Be"],
                          fast="Al26", slow="Be10",
                          tag="step_pair_not_unique")
        check("step_pair reports NotUnique for an unidentifiable Al26/Be10 pair",
              p_amb.get("status") == "not unique" and int(float(p_amb["n_brackets"])) > 1,
              f"status={p_amb.get('status')} brackets={p_amb.get('n_brackets')}")

        # --- 10. burial_pair_batch processes multiple rows --------------------
        burial_in = workdir / "burial_pairs.csv"
        burial_in.write_text(
            "sample_name,nuclide_fast,nuclide_slow,conc_fast_atoms_g,conc_slow_atoms_g,total_shielding\n"
            + f"buried_al26,Al26,Be10,{n26},{n10},1.0\n"
            + f"buried_c14,C14,Be10,{n14},{n10_c14},1.0\n"
            + "bad,Al26,Be10,0,0,1.0\n"
        )
        burial_rows = run(exe, workdir, {
            "column_mode": "burial_pair_batch",
            "burial_pairs_csv_file": str(burial_in),
        }, "burial_pair_batch", numeric=False)
        check("burial_pair_batch writes one result row per input row",
              len(burial_rows) == 3, f"got {len(burial_rows)}")
        check("burial_pair_batch converged row recovers erosion rate",
              burial_rows[0]["status"] == "converged"
              and rel(float(burial_rows[0]["erate_preburial_g_cm2_yr"]), E_TRUE) < 1e-5,
              f"status={burial_rows[0]['status']} erate={burial_rows[0]['erate_preburial_g_cm2_yr']}")
        check("burial_pair_batch supports C14 rows",
              burial_rows[1]["status"] == "converged"
              and burial_rows[1]["nuclide_fast"] == "C14",
              f"status={burial_rows[1]['status']} fast={burial_rows[1]['nuclide_fast']}")
        check("burial_pair_batch emits uncertainty columns",
              "erate_preburial_unc_int_g_cm2_yr" in burial_rows[0]
              and "erate_preburial_unc_ext_g_cm2_yr" in burial_rows[0])
        check("burial_pair_batch reports impossible rows as no solution",
              burial_rows[2]["status"] == "no solution",
              f"status={burial_rows[2]['status']}")

        # --- 11. step_pair_batch supports per-row nuclide/status output -------
        step_in = workdir / "step_pairs.csv"
        step_in.write_text(
            "sample_name,nuclide_fast,nuclide_slow,conc_fast_atoms_g,conc_slow_atoms_g,effective_erosion_rate_old_g_cm2_yr\n"
            + f"resolvable,C14,Be10,{step_true['C_14C']},{step_true['C_10Be']},{E_OLD}\n"
            + f"ambiguous,Al26,Be10,{step_amb['C_26Al']},{step_amb['C_10Be']},{E_OLD}\n"
        )
        step_rows = run(exe, workdir, {
            "column_mode": "step_pair_batch",
            "step_pairs_csv_file": str(step_in),
        }, "step_pair_batch", numeric=False)
        check("step_pair_batch writes one result row per input row",
              len(step_rows) == 2, f"got {len(step_rows)}")
        check("step_pair_batch resolvable row converges and recovers new rate",
              step_rows[0]["status"] == "converged"
              and rel(float(step_rows[0]["erate_new_g_cm2_yr"]), E_NEW) < 1e-3,
              f"status={step_rows[0]['status']} erate={step_rows[0]['erate_new_g_cm2_yr']}")
        check("step_pair_batch emits uncertainty columns",
              "erate_new_unc_int_g_cm2_yr" in step_rows[0]
              and "time_since_change_unc_int_yr" in step_rows[0])
        check("step_pair_batch ambiguous row reports not unique with >1 bracket",
              step_rows[1]["status"] == "not unique"
              and int(float(step_rows[1]["n_brackets"])) > 1,
              f"status={step_rows[1]['status']} brackets={step_rows[1]['n_brackets']}")

        # --- 12. an unknown mode is a clean failure, not a silent no-op ------
        param = workdir / "badmode.param"
        param.write_text(f"write_path: {workdir}\ncolumn_mode: not_a_mode\n")
        proc = subprocess.run([str(exe), str(workdir), param.name],
                              cwd=str(workdir), capture_output=True, text=True)
        check("an unknown column_mode exits non-zero and lists the valid modes",
              proc.returncode != 0 and "steady_state" in proc.stderr,
              f"exit {proc.returncode}")
    except Exception as exc:                                   # noqa: BLE001
        check("harness ran to completion", False, str(exc))
    finally:
        if args.keep:
            print(f"kept working directory: {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)

    print(f"\n{n_pass}/{n_pass + n_fail} column-tool checks passed")
    return EXIT_OK if n_fail == 0 else EXIT_FAIL


if __name__ == "__main__":
    sys.exit(main())
