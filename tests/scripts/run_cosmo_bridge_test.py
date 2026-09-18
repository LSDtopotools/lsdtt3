#!/usr/bin/env python3
"""Forward -> inverse CLI bridge, plus knob coverage (cosmo CLI test plan Phase 4).

THE POINT OF THIS TEST. Every other cosmo test drives the LIBRARY. This one
drives the two CLIs against each other: lsdtt-cosmo-prediction forward-models a
concentration for a KNOWN erosion rate, writes the points CSV it advertises as
ingestible, and lsdtt-cosmo-tool inverts that CSV back to an erosion rate. The
imposed rate must come back. That closes a loop no library test can close,
because the failure modes here are plumbing -- a column renamed on one side, a
knob that reaches the forward model but not the inverse, a shielding term
applied twice -- and every one of them produces perfectly plausible numbers.

The knob cases exist for exactly that reason. Each one changes a parameter on
BOTH sides; if a knob is honoured asymmetrically the round trip stops closing
even though neither CLI errors. A case that passes only because a knob is
ignored on both sides would be worthless, so the knobs are chosen to MOVE the
concentration substantially (verified: see --report), and a self-test asserts
that they do.

Not a physics test: it asserts self-consistency of the pair, so it cannot
validate the physics itself. The sentinel, the CRONUS comparison and the
real-landscape tiers do that.

Exit codes: 0 pass, 1 fail, 77 skip (missing executables).
"""

import argparse
import csv
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import lsdtt_test_env

EXIT_OK, EXIT_FAIL, EXIT_SKIP = 0, 1, 77

# The bridge closes to a few 1e-4 relative. The floor is not the inversion
# solver (it converges to ~1e-3 on the concentration, which maps to far less on
# the rate) but the fact that the forward model hands over the accumulated
# concentration through a FLOAT32 raster: 55120.117188 above, ~7 significant
# digits. Measured worst case at first run was 6.2e-4; frozen at 5e-3, the
# house ~8x convention, so ordinary float32 noise cannot make this flap.
RECOVERY_REL_TOL = 5.0e-3

# A knob must actually change the forward concentration, or the case it guards
# proves nothing. Measured smallest real effect at first run was ~2.2 % (aged
# Stone); frozen an order of magnitude below that.
KNOB_EFFECT_MIN_REL = 2.0e-3

# Per-knob overrides, each with its measured value and the physical reason the
# general floor is wrong for it. An override is a statement about the FIXTURE,
# never a concession to a knob that might not be wired up.
KNOB_EFFECT_MIN_OVERRIDE = {
    # ERA40 vs NCEP differ by only 0.057 % HERE because tiny_dem.asc sits at
    # ~sea level on the equator, where the two climatologies very nearly agree.
    # The floor is set just under the measured value: it still distinguishes
    # "NCEP is honoured" from "NCEP is silently ignored" (which would be exactly
    # 0), which is the only thing this check is for. Raising the effect would
    # mean a high-altitude, high-latitude fixture, not a looser test.
    "atm_NCEP": 4.0e-4,
}


def q(v):
    """Render a Python value as the lsdtt param parser expects it."""
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


def run_cli(exe, workdir, param_name, label):
    proc = lsdtt_test_env.run([str(exe), str(workdir), param_name],
                          cwd=str(workdir), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{label} failed (exit {proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-3000:]}\n"
            f"--- stderr ---\n{proc.stderr[-3000:]}")
    return proc


def forward(pred_exe, workdir, data_dir, erate, knobs, tag):
    """Forward-model the accumulated concentration; return (conc, csv_path)."""
    base = f"fwd_{tag}"
    lines = [
        f"read_path: {data_dir}",
        "read_fname: tiny_dem.asc",
        f"write_path: {workdir}",
        f"write_fname: {base}",
        "calculate_crn_concentration_raster: true",
        "calculate_accumulated_crn_concentration: true",
        f"effective_erosion_rate: {erate}",
        "sample_points_fname: cosmo_tiny_points.csv",
        "threshold_contributing_pixels: 2",
        "sample_snap_max_search_radius_pixels: 3",
    ]
    lines += [f"{k}: {q(v)}" for k, v in knobs.items()]
    param = workdir / f"{base}.param"
    param.write_text("\n".join(lines) + "\n")
    run_cli(pred_exe, workdir, param.name, f"forward[{tag}]")

    out = workdir / f"{base}_sampled_crn_conc.csv"
    if not out.is_file():
        raise RuntimeError(f"forward[{tag}] wrote no points CSV: {out}")
    rows = list(csv.DictReader(out.open(newline="")))
    if len(rows) != 1:
        raise RuntimeError(f"forward[{tag}] sampled {len(rows)} points, expected 1")
    conc_col = next(c for c in rows[0] if c.endswith("_CONC"))
    return float(rows[0][conc_col]), out


def inverse(tool_exe, workdir, data_dir, points_csv, knobs, tag):
    """Invert the forward CSV; return the full result row."""
    base = f"inv_{tag}"
    lines = [
        f"read_path: {data_dir}",
        "read_fname: tiny_dem.asc",
        f"write_path: {workdir}",
        f"write_fname: {base}",
        f"sample_points_fname: {points_csv}",
        "threshold_contributing_pixels: 2",
        "sample_snap_max_search_radius_pixels: 3",
    ]
    lines += [f"{k}: {q(v)}" for k, v in knobs.items()]
    param = workdir / f"{base}.param"
    param.write_text("\n".join(lines) + "\n")
    run_cli(tool_exe, workdir, param.name, f"inverse[{tag}]")

    out = workdir / f"{base}_CN_ERATE.csv"
    if not out.is_file():
        raise RuntimeError(f"inverse[{tag}] wrote no result CSV: {out}")
    rows = list(csv.DictReader(out.open(newline="")))
    if len(rows) != 1:
        raise RuntimeError(f"inverse[{tag}] returned {len(rows)} rows, expected 1")
    return rows[0]


# --------------------------------------------------------------------------- #
# Cases.  Each is (tag, erate, shared_knobs) applied to BOTH CLIs.
# --------------------------------------------------------------------------- #
CASES = [
    ("baseline",      0.01,   {}),
    # Rate sweep: slow samples are muon-heavy, fast ones spallation-dominated,
    # so this exercises different parts of the same forward model.
    ("slow",          0.0005, {}),
    ("fast",          0.05,   {}),
    # Nuclide: different P0, F, decay constant and standardisation map.
    ("nuclide_Al26",  0.01,   {"nuclide": "Al26"}),
    ("nuclide_C14",   0.01,   {"nuclide": "C14"}),
    # Muon scheme: changes F and Gamma, which is where an asymmetric knob shows.
    ("muon_Braucher", 0.01,   {"muon_scheme": "Braucher"}),
    ("muon_Granger",  0.01,   {"muon_scheme": "Granger"}),
    ("muon_Schaller", 0.01,   {"muon_scheme": "Schaller"}),
    ("muon_newCRONUS", 0.01,  {"muon_scheme": "newCRONUS"}),
    # CRONUS v3 emulator (plan Phase 3).  These take a DIFFERENT code path from
    # every case above -- a dedicated per-pixel branch in the forward raster and
    # a per-pixel basin precompute in the inversion, because their muon terms
    # are pressure-dependent.  The generic reset_scaling()/scale_F_values() path
    # is not used at all, so nothing else here covers this plumbing.  Swept
    # across erosion rate as well: the muon share of the inventory falls with
    # rate, so a muon term wired into the forward but not the inverse (or
    # vice versa) shows up at one end of the sweep even if it hides at the other.
    ("emu_st_slow",   0.001,  {"muon_scheme": "cronusEmulatorSt"}),
    ("emu_st",        0.01,   {"muon_scheme": "cronusEmulatorSt"}),
    ("emu_st_fast",   0.3,    {"muon_scheme": "cronusEmulatorSt"}),
    ("emu_st_Al26",   0.01,   {"muon_scheme": "cronusEmulatorSt",
                               "nuclide": "Al26"}),
    ("emu_st_C14",    0.01,   {"muon_scheme": "cronusEmulatorSt",
                               "nuclide": "C14"}),
    # Shielding matters especially here: the emulator applies the site scaling
    # and topographic shielding to the SPALLATION pathway only, so an
    # accidental application to the muon pathways breaks the round trip.
    ("emu_st_shield", 0.01,   {"muon_scheme": "cronusEmulatorSt",
                               "topo_shield_value": 0.8,
                               "snow_shield_eff_thickness": 20.0,
                               "self_shield_eff_thickness": 15.0}),
    ("emu_lsdn",      0.01,   {"muon_scheme": "cronusEmulatorLSDn",
                               "scaling_scheme": "lsdn"}),
    ("emu_lsdn_fast", 0.1,    {"muon_scheme": "cronusEmulatorLSDn",
                               "scaling_scheme": "lsdn"}),
    # Production scaling, both schemes, and LSDn with a time-averaged field.
    ("scaling_lsdn",  0.01,   {"scaling_scheme": "lsdn"}),
    ("lsdn_aged",     0.01,   {"scaling_scheme": "lsdn", "sample_age_yr": 20000}),
    ("stone_aged",    0.01,   {"sample_age_yr": 20000}),
    # Shielding, one term at a time then together -- a term applied on one side
    # only, or twice, breaks the round trip while both CLIs still succeed.
    ("shield_topo",   0.01,   {"topo_shield_value": 0.8}),
    ("shield_snow",   0.01,   {"snow_shield_eff_thickness": 20.0}),
    ("shield_self",   0.01,   {"self_shield_eff_thickness": 15.0}),
    ("shield_all",    0.01,   {"topo_shield_value": 0.8,
                               "snow_shield_eff_thickness": 20.0,
                               "self_shield_eff_thickness": 15.0}),
    # Atmosphere. Covered on BOTH sides since `atmosphere` was added to
    # lsdtt-cosmo-prediction (2026-07-27); before that it was tool-only, so no
    # matching atmosphere could be imposed and setting it on one side alone
    # closed to 0.38 % on this fixture -- a pass that established nothing.
    ("atm_NCEP",      0.01,   {"atmosphere": "NCEP"}),
    ("atm_standard",  0.01,   {"atmosphere": "standard"}),
    ("atm_ERA40",     0.01,   {"atmosphere": "ERA40"}),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-exe", default=os.environ.get("LSDTT_COSMO_PREDICTION"))
    ap.add_argument("--tool-exe", default=os.environ.get("LSDTT_COSMO_TOOL"))
    ap.add_argument("--data-dir", type=pathlib.Path,
                    default=pathlib.Path(__file__).resolve().parents[1] / "data")
    ap.add_argument("--keep", action="store_true",
                    help="keep the working directory for debugging")
    ap.add_argument("--report", action="store_true",
                    help="print the forward concentration of every case")
    args = ap.parse_args()

    if not args.pred_exe or not args.tool_exe:
        print("SKIP: set --pred-exe/--tool-exe (or LSDTT_COSMO_PREDICTION / "
              "LSDTT_COSMO_TOOL); nothing to bridge.")
        return EXIT_SKIP
    pred_exe, tool_exe = pathlib.Path(args.pred_exe), pathlib.Path(args.tool_exe)
    for e in (pred_exe, tool_exe):
        if not e.is_file():
            print(f"SKIP: executable not found: {e}")
            return EXIT_SKIP
    for f in ("tiny_dem.asc", "cosmo_tiny_points.csv"):
        if not (args.data_dir / f).is_file():
            print(f"FAIL: fixture missing: {args.data_dir / f}")
            return EXIT_FAIL

    workdir = pathlib.Path(tempfile.mkdtemp(prefix="cosmo_bridge_"))
    n_pass = n_fail = 0
    concs = {}
    try:
        # The fixture points CSV lives beside the DEM; the forward CLI resolves
        # it against read_path, so nothing is copied.
        for tag, erate, knobs in CASES:
            try:
                conc, points_csv = forward(pred_exe, workdir, args.data_dir,
                                           erate, knobs, tag)
                row = inverse(tool_exe, workdir, args.data_dir, points_csv,
                              knobs, tag)
                got = float(row["apparent_eff_e_g_cm2_yr"])
                rel = abs(got - erate) / erate
                ok = rel <= RECOVERY_REL_TOL
            except Exception as exc:                       # noqa: BLE001
                print(f"[FAIL] {tag}: {exc}")
                n_fail += 1
                continue
            concs[tag] = conc
            print(f"[{'PASS' if ok else 'FAIL'}] {tag}: imposed {erate:g} "
                  f"recovered {got:.6g} rel={rel:.3g} (conc {conc:.6g})")
            n_pass, n_fail = (n_pass + 1, n_fail) if ok else (n_pass, n_fail + 1)

        # ---- the knobs must MATTER, or the cases above prove nothing --------
        base = concs.get("baseline")
        if base is None:
            print("[FAIL] knob-effect check: baseline concentration missing")
            n_fail += 1
        else:
            for tag in ("nuclide_Al26", "muon_Braucher", "scaling_lsdn",
                        "shield_topo", "shield_snow", "shield_self",
                        "lsdn_aged", "stone_aged", "slow", "fast",
                        "atm_NCEP", "atm_standard"):
                if tag not in concs:
                    continue
                rel = abs(concs[tag] - base) / base
                floor = KNOB_EFFECT_MIN_OVERRIDE.get(tag, KNOB_EFFECT_MIN_REL)
                ok = rel >= floor
                print(f"[{'PASS' if ok else 'FAIL'}] knob-effect {tag}: "
                      f"moved the forward concentration by {rel:.4%} "
                      f"(need >= {floor:.4%})")
                n_pass, n_fail = (n_pass + 1, n_fail) if ok else (n_pass, n_fail + 1)

        # ---- density is a pure unit conversion, so assert the identity ------
        for rho in (2650.0, 2000.0):
            tag = f"density_{int(rho)}"
            try:
                conc, points_csv = forward(pred_exe, workdir, args.data_dir,
                                           0.01, {}, tag)
                row = inverse(tool_exe, workdir, args.data_dir, points_csv,
                              {"rock_density_kg_m3": rho}, tag)
                eff = float(row["apparent_eff_e_g_cm2_yr"])
                # BOTH conversions are asserted. This check previously read
                # apparent_e_mm_yr and compared it to the mm/KYR formula, which
                # passed only because the column was written with the mm/kyr
                # factor -- the test was pinning the bug in place. Checking the
                # two columns against their own formulas is what makes a
                # repeat of that impossible.
                for col, factor, unit in (
                        ("apparent_e_mm_kyr", 1.0e7, "mm/kyr"),
                        ("apparent_e_mm_yr", 1.0e4, "mm/yr")):
                    got = float(row[col])
                    expect = eff * factor / rho
                    rel = abs(got - expect) / expect
                    ok = rel <= 1.0e-6  # a printed-precision check, not physics
                    print(f"[{'PASS' if ok else 'FAIL'}] {tag}: {unit} "
                          f"{got:.6g} vs eff*{factor:.0e}/rho {expect:.6g} "
                          f"rel={rel:.3g}")
                    n_pass, n_fail = ((n_pass + 1, n_fail) if ok
                                      else (n_pass, n_fail + 1))
                # The two must also differ by exactly 1000, which is the
                # relationship the bug destroyed.
                ratio = (float(row["apparent_e_mm_kyr"])
                         / float(row["apparent_e_mm_yr"]))
                ok = abs(ratio - 1000.0) / 1000.0 <= 1.0e-6
                print(f"[{'PASS' if ok else 'FAIL'}] {tag}: mm/kyr / mm/yr "
                      f"= {ratio:.6g} (must be 1000)")
                n_pass, n_fail = (n_pass + 1, n_fail) if ok else (n_pass, n_fail + 1)
            except Exception as exc:                       # noqa: BLE001
                print(f"[FAIL] {tag}: {exc}")
                n_fail += 1

        # ---- standardisation: a non-unity factor must be undone -------------
        # KNSTD -> 07KNSTD is 0.9042. Feed conc/0.9042 labelled KNSTD; the tool
        # multiplies by the factor, recovering the true concentration and hence
        # the imposed rate. If the factor were ignored, or applied the wrong
        # way, the rate would be off by ~10 % -- far outside RECOVERY_REL_TOL.
        try:
            factor, erate = 0.9042, 0.01
            conc, _ = forward(pred_exe, workdir, args.data_dir, erate, {}, "std")
            std_csv = workdir / "std_points.csv"
            with std_csv.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["sample_name", "latitude", "longitude",
                            "Be10_CONC", "Be10_standardisation"])
                src = list(csv.DictReader(
                    (workdir / "fwd_std_sampled_crn_conc.csv").open(newline="")))[0]
                w.writerow([src["sample_name"], src["latitude"], src["longitude"],
                            conc / factor, "KNSTD"])
            row = inverse(tool_exe, workdir, args.data_dir, std_csv,
                          {"standardisation_column_name": "Be10_standardisation"},
                          "std")
            got = float(row["apparent_eff_e_g_cm2_yr"])
            rel = abs(got - erate) / erate
            ok = rel <= RECOVERY_REL_TOL
            print(f"[{'PASS' if ok else 'FAIL'}] standardisation KNSTD: "
                  f"imposed {erate:g} recovered {got:.6g} rel={rel:.3g}")
            n_pass, n_fail = (n_pass + 1, n_fail) if ok else (n_pass, n_fail + 1)
        except Exception as exc:                           # noqa: BLE001
            print(f"[FAIL] standardisation: {exc}")
            n_fail += 1

        # ---- +/-1 sigma bounds must BRACKET and be correctly ordered --------
        # More nuclide means slower erosion, so the +sigma rate is the LOWER
        # one. Getting this backwards is a classic sign error that no
        # single-value test can see.
        try:
            erate = 0.01
            conc, _ = forward(pred_exe, workdir, args.data_dir, erate, {}, "sig")
            sig_csv = workdir / "sig_points.csv"
            with sig_csv.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["sample_name", "latitude", "longitude",
                            "Be10_CONC", "Be10_uncertainty"])
                src = list(csv.DictReader(
                    (workdir / "fwd_sig_sampled_crn_conc.csv").open(newline="")))[0]
                w.writerow([src["sample_name"], src["latitude"], src["longitude"],
                            conc, conc * 0.05])
            row = inverse(tool_exe, workdir, args.data_dir, sig_csv, {}, "sig")
            e_mid = float(row["apparent_eff_e_g_cm2_yr"])
            e_lo_conc = float(row["apparent_eff_e_minus_g_cm2_yr"])   # conc - sigma
            e_hi_conc = float(row["apparent_eff_e_plus_g_cm2_yr"])    # conc + sigma
            ok = e_hi_conc < e_mid < e_lo_conc
            print(f"[{'PASS' if ok else 'FAIL'}] +/-1sigma bracketing: "
                  f"conc+sigma {e_hi_conc:.6g} < mid {e_mid:.6g} < "
                  f"conc-sigma {e_lo_conc:.6g}")
            n_pass, n_fail = (n_pass + 1, n_fail) if ok else (n_pass, n_fail + 1)
        except Exception as exc:                           # noqa: BLE001
            print(f"[FAIL] +/-1sigma: {exc}")
            n_fail += 1

        if args.report:
            print("\nForward concentrations (atoms/g):")
            for tag, c in concs.items():
                print(f"  {tag:18s} {c:.8g}")
    finally:
        if args.keep:
            print(f"kept working directory: {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)

    total = n_pass + n_fail
    print(f"\n{n_pass}/{total} bridge checks passed")
    return EXIT_OK if n_fail == 0 else EXIT_FAIL


if __name__ == "__main__":
    sys.exit(main())
