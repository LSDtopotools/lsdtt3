#!/usr/bin/env python3
"""Append the LSDn (Lifton-Sato-Dunai 2014) scaling data to crn_parameters.nc.

lsdtt3 implements LSDn production scaling by bilinearly interpolating the
precomputed `b`/`m` scaling grids that Balco's v3 online calculator itself
interpolates (`get_LSDnSF.m`), rather than by porting the Sato/PARMA spectral
integration.  This script packs those grids, together with the Lifton solar
modulation record and the Pavon-Carrasco paleomagnetic records needed for the
time-dependent cutoff rigidity, into the bundled parameter file.

It mirrors `add_ncep_to_crn_nc.py`: a non-destructive, idempotent append.  It
creates dimensions and variables that do not exist and overwrites values in
place for those that do.  It NEVER touches a pre-existing variable -- that is
verified, not assumed: every existing variable's bytes are read before the
write and compared after it, and a mismatch is a hard failure.

Input is the intermediate `.npz` written by `extract_lsdn_reference.py`, which
carries the sha256 of every source `.mat` alongside the arrays; those hashes
are copied into a global `lsdn_provenance` attribute so the packed data can be
traced back to exactly the files that were checked in Phase 0.  Alternatively
`--data-dir` re-runs the extraction (and its C1-C7 checks) in-process.

Encoding.  The six 380x201 grids are float32 + zlib complevel 9 + shuffle.
That combination was measured, not guessed: float32 round-trips the float64
source to ~5.9e-8 relative (far below bilinear-interpolation error, and the
same order as the existing NCEP packing), and `shuffle` is worth ~0.49 MiB of
the compressed payload -- do not drop it.  The small f8 records are stored
uncompressed.  Total added payload is ~1.03 MiB.

Storage order for the grids is `(lsdn_p, lsdn_rc)` -- p slow, Rc fast -- so the
C++ flat index is `ip * 201 + irc` and the MATLAB 380x201 shape needs no
transpose at read time.  The shape is asserted before writing.

LICENCE NOTE.  The v3 reference implementation is GPL-v2.  What is packed here
is DATA (grid axes, coefficient tables, paleomagnetic and solar records) --
facts, not expression -- cited as reference data on the same footing as the
existing Lifton lapse-rate polynomial and Knudsen VDM record.  No MATLAB source
text is copied into this repository.

Usage:
    pixi run python python_helpers/add_lsdn_to_crn_nc.py \
        --npz /path/to/scratch/lsdn/lsdn_reference.npz \
        --nc data/crn_parameters.nc --verify
"""
import argparse
import json
import os
import sys

import numpy as np
import netCDF4 as nc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extract_lsdn_reference as ref  # noqa: E402  (sibling helper, path set above)

N_RC, N_P = ref.RC_N, ref.P_N
NUCLIDES = ("Be10", "Al26", "C14")

CITE_LSDN = "Lifton, Sato & Dunai (2014) EPSL 386:149, doi:10.1016/j.epsl.2013.10.052"
CITE_SATO = "Sato et al. (2008) J. Nucl. Sci. Technol. 45:392, doi:10.1080/18811248.2008.9711443"
CITE_PAVON = ("Pavon-Carrasco et al. (2014) EPSL 388:98, doi:10.1016/j.epsl.2013.11.046; "
              "compiled for Balco v3 by N. Lifton")
CITE_SPHI = ("Usoskin et al. (2011) JGR 116:A02104, doi:10.1029/2010JA016105; "
             "as used by Lifton et al. (2014)")

# The single constraint that must travel with the data, repeated verbatim in
# pavon_t2's attributes because that is where a future reader will look.
NOTE_T2 = ("get_DipRc.m branch threshold ONLY (elseif in.t <= m.t2(end)). "
           "NEVER use as an averaging weight -- the duration weighting uses "
           "pavon_t2min/pavon_t2max. Conversely NEVER use pavon_t2max as the "
           "branch threshold: t2[-1] = 2.00074e7 yr but t2max[-1] = 2.02044e7 yr.")


def snapshot_existing(path):
    """Read every existing variable's raw bytes + attributes, for the
    non-destructiveness check.  Done before the append so the comparison is
    against the real prior state rather than against an assumption."""
    snap = {}
    with nc.Dataset(path, "r") as ds:
        for name, v in ds.variables.items():
            snap[name] = (np.asarray(v[:]).tobytes(),
                          {a: repr(v.getncattr(a)) for a in v.ncattrs()})
        gattrs = {a: repr(ds.getncattr(a)) for a in ds.ncattrs()}
    return snap, gattrs


def verify_unchanged(path, snap, gattrs):
    """Confirm every pre-existing variable is byte-identical afterwards."""
    bad = []
    with nc.Dataset(path, "r") as ds:
        for name, (raw, attrs) in snap.items():
            if name not in ds.variables:
                bad.append(f"{name}: DISAPPEARED")
                continue
            v = ds.variables[name]
            if np.asarray(v[:]).tobytes() != raw:
                bad.append(f"{name}: VALUES CHANGED")
            now = {a: repr(v.getncattr(a)) for a in v.ncattrs()}
            if now != attrs:
                bad.append(f"{name}: attributes changed {attrs} -> {now}")
        for a, val in gattrs.items():
            if a not in ds.ncattrs() or repr(ds.getncattr(a)) != val:
                bad.append(f"global attr {a}: changed or removed")
    return bad


def load_arrays(args):
    """Return (arrays, provenance) from the .npz, or by re-extracting."""
    if args.npz:
        z = np.load(args.npz)
        arrays = {k: z[k] for k in z.files if k != "provenance_json"}
        prov = json.loads(str(z["provenance_json"])) if "provenance_json" in z.files else {}
        print(f"loaded {len(arrays)} arrays from {args.npz}")
        return arrays, prov
    log = ref.CheckLog()
    arrays, prov = ref.extract(args.data_dir, log)
    if not log.all_ok:
        raise SystemExit("extraction checks C1-C7 FAILED -- refusing to pack")
    return arrays, prov


def pack(ds, arrays, prov):
    """Create/overwrite the LSDn dimensions and variables.  Returns the list of
    variable names written."""
    for dim, n in (("lsdn_rc", N_RC), ("lsdn_p", N_P),
                   ("pavon_short", int(arrays["pavon_mm01"].size)),
                   ("pavon_long", int(arrays["pavon_mm02"].size))):
        if dim not in ds.dimensions:
            ds.createDimension(dim, n)
        elif len(ds.dimensions[dim]) != n:
            raise ValueError(f"dimension {dim} exists with length "
                             f"{len(ds.dimensions[dim])}, need {n}")

    grid_kw = dict(zlib=True, complevel=9, shuffle=True)
    written = []

    def put(name, dtype, dims, arr, units, longname, source, **kw):
        a = np.asarray(arr)
        want = tuple(len(ds.dimensions[d]) for d in dims)
        if a.shape != want:
            raise ValueError(f"{name}: array shape {a.shape} != dims {dims} {want}")
        if name not in ds.variables:
            v = ds.createVariable(name, dtype, dims, **kw)
        else:
            v = ds.variables[name]
        v[...] = a.astype(dtype)
        v.units = units
        v.long_name = longname
        v.source = source
        written.append(name)

    # -- axes -----------------------------------------------------------------
    put("lsdn_rc", "f8", ("lsdn_rc",), arrays["lsdn_rc"], "GV",
        "LSDn cutoff-rigidity axis, ascending 0 to 20 step 0.1",
        f"LSDn2015_grid.mat (sfgrid.Rc); {CITE_LSDN}")
    put("lsdn_p", "f8", ("lsdn_p",), arrays["lsdn_p"], "hPa",
        "LSDn pressure axis, ascending 300 to 1058 step 2",
        f"LSDn2015_grid.mat (sfgrid.p); {CITE_LSDN}")

    # -- nuclide grids --------------------------------------------------------
    # b is the solar-independent term, m the solar-modulation slope, so the
    # raw scaling for solar modulation S (MV) is b(Rc,p) + S*m(Rc,p).
    mat_of = {n: ref.NUCLIDE_MATS[n][0] for n in NUCLIDES}
    struct_of = {n: ref.NUCLIDE_MATS[n][1] for n in NUCLIDES}
    label = {"Be10": "Be-10", "Al26": "Al-26", "C14": "C-14"}
    for n in NUCLIDES:
        src = (f"{mat_of[n]} ({struct_of[n]}.%s); {CITE_LSDN}; "
               f"spectra from {CITE_SATO}")
        put(f"lsdn_b_{n}", "f4", ("lsdn_p", "lsdn_rc"), arrays[f"lsdn_b_{n}"],
            "atoms g^-1 yr^-1",
            f"LSDn scaling grid, solar-independent term b, {label[n]} in quartz "
            f"(uncalibrated; raw scaling = b + S*m)", src % "b", **grid_kw)
        put(f"lsdn_m_{n}", "f4", ("lsdn_p", "lsdn_rc"), arrays[f"lsdn_m_{n}"],
            "atoms g^-1 yr^-1 MV^-1",
            f"LSDn scaling grid, solar-modulation slope m, {label[n]} in quartz",
            src % "m", **grid_kw)

    # -- solar modulation record ---------------------------------------------
    put("lsdn_sphi", "f8", ("pavon_short",), arrays["lsdn_sphi"], "MV",
        "Solar modulation potential per Holocene time step",
        f"LiftonS2015.mat (SPhi); {CITE_SPHI}")
    put("lsdn_sphi_inf", "f8", (), arrays["lsdn_sphi_inf"], "MV",
        "Long-term mean solar modulation potential, used for the long record",
        f"LiftonS2015.mat (SPhiInf); {CITE_SPHI}")

    # -- Pavon paleomagnetic records ------------------------------------------
    src_p = f"PavonDip.mat (m.%s); {CITE_PAVON}"
    put("pavon_mm01", "f8", ("pavon_short",), arrays["pavon_mm01"], "1",
        "Dipole moment ratio M_t/M_0 for the Holocene (short) record",
        src_p % "MM01")
    put("pavon_t1min", "f8", ("pavon_short",), arrays["pavon_t1min"],
        "yr before pavon_tzero",
        "Lower bound of each Holocene time step", src_p % "t1min")
    put("pavon_t1max", "f8", ("pavon_short",), arrays["pavon_t1max"],
        "yr before pavon_tzero",
        "Upper bound of each Holocene time step "
        "(t1max[-1] abuts pavon_t2min[0] exactly)", src_p % "t1max")
    put("pavon_pole_lat", "f8", ("pavon_short",), arrays["pavon_pole_lat"], "rad",
        "Paleopole latitude per Holocene time step", src_p % "latr_pp")
    put("pavon_pole_lon", "f8", ("pavon_short",), arrays["pavon_pole_lon"], "rad",
        "Paleopole longitude per Holocene time step", src_p % "lonr_pp")
    put("pavon_mm02", "f8", ("pavon_long",), arrays["pavon_mm02"], "1",
        "Dipole moment ratio M_t/M_0 for the long (GAD) record", src_p % "MM02")
    put("pavon_t2min", "f8", ("pavon_long",), arrays["pavon_t2min"],
        "yr before pavon_tzero",
        "Lower bound of each long-record time step "
        "(t2min[0] abuts pavon_t1max[-1] exactly)", src_p % "t2min")
    put("pavon_t2max", "f8", ("pavon_long",), arrays["pavon_t2max"],
        "yr before pavon_tzero",
        "Upper bound of each long-record time step", src_p % "t2max")
    put("pavon_t2", "f8", ("pavon_long",), arrays["pavon_t2"],
        "yr before pavon_tzero",
        "Long-record time-step CENTRES. " + NOTE_T2, (src_p % "t2") + ". " + NOTE_T2)
    put("pavon_tzero", "f8", (), arrays["pavon_tzero"], "yr CE",
        "Reference year of the paleomagnetic records", src_p % "tzero")

    # Provenance travels with the file: the sha256 of each source .mat, so the
    # packed values can be traced to exactly the files Phase 0 checked.
    ds.lsdn_provenance = json.dumps({
        "source_sha256": prov.get("sha256", {}),
        "extracted_utc": prov.get("extracted_utc", "unknown"),
        "extractor": "python_helpers/extract_lsdn_reference.py",
        "packer": "python_helpers/add_lsdn_to_crn_nc.py",
        "encoding": "b/m grids f4 zlib-9 shuffle; records f8 uncompressed",
        "citations": [CITE_LSDN, CITE_SATO, CITE_PAVON, CITE_SPHI],
    }, sort_keys=True)
    return written


def verify(path, arrays):
    """Re-read the packed file and reproduce C1/C2/C3 from the stored arrays,
    plus the float32 round-trip deviation against the float64 source."""
    print("\n--- verify: re-reading the packed file ---")
    ok = True
    with nc.Dataset(path, "r") as ds:
        rc = np.asarray(ds.variables["lsdn_rc"][:], dtype=np.float64)
        p = np.asarray(ds.variables["lsdn_p"][:], dtype=np.float64)
        b10 = np.asarray(ds.variables["lsdn_b_Be10"][:], dtype=np.float64)
        m10 = np.asarray(ds.variables["lsdn_m_Be10"][:], dtype=np.float64)

        ok_rc, det_rc = ref.axis_report("Rc [GV]", rc, ref.RC_MIN, ref.RC_STEP, ref.RC_N)
        ok_p, det_p = ref.axis_report("p [hPa]", p, ref.P_MIN, ref.P_STEP, ref.P_N)
        print(f"C1 {'PASS' if ok_rc and ok_p else 'FAIL'}: {det_rc}; {det_p}")
        ok &= ok_rc and ok_p

        got = (float(b10.min()), float(b10.max()), float(m10.min()), float(m10.max()))
        # Same tolerances as the extractor's C2 -- the float32 packing costs
        # ~6e-8, nothing like enough to need looser bounds here.
        ok_c2 = (ref.close(got[0], ref.C2_B_RANGE[0], ref.TOL_QUOTED_REL_LOOSE)
                 and ref.close(got[1], ref.C2_B_RANGE[1], ref.TOL_QUOTED_REL)
                 and ref.close(got[2], ref.C2_M_RANGE[0], ref.TOL_QUOTED_REL_LOOSE)
                 and ref.close(got[3], ref.C2_M_RANGE[1], ref.TOL_QUOTED_REL_LOOSE))
        print(f"C2 {'PASS' if ok_c2 else 'FAIL'}: b {got[0]:.6g}..{got[1]:.6g}, "
              f"m {got[2]:.6g}..{got[3]:.6g} (expected {ref.C2_B_RANGE}, {ref.C2_M_RANGE})")
        ok &= ok_c2

        v = (ref.bilinear_regular(b10, rc, p, ref.C3_RC_GV, ref.C3_P_HPA)
             + ref.C3_S_MV * ref.bilinear_regular(m10, rc, p, ref.C3_RC_GV, ref.C3_P_HPA))
        # float32 storage costs ~6e-8 relative, so C3 is checked at 1e-4 (the
        # precision the expected value is quoted at), not at machine epsilon.
        ok_c3 = ref.close(v, ref.C3_EXPECTED, ref.TOL_QUOTED_REL)
        print(f"C3 {'PASS' if ok_c3 else 'FAIL'}: {v:.6f} (expected {ref.C3_EXPECTED}, "
              f"off-node in p -- also confirms axis orientation survived packing)")
        ok &= ok_c3

        worst = 0.0
        for n in NUCLIDES:
            for kind in ("b", "m"):
                name = f"lsdn_{kind}_{n}"
                if name not in ds.variables:
                    continue
                src = np.asarray(arrays[name], dtype=np.float64)
                back = np.asarray(ds.variables[name][:], dtype=np.float64)
                d = float(np.max(np.abs(back - src) / np.maximum(np.abs(src), 1e-300)))
                print(f"  {name}: max float32 round-trip rel dev {d:.3e}")
                worst = max(worst, d)
        print(f"worst float32 round-trip relative deviation: {worst:.3e} "
              f"(expected ~5.9e-8; an input to the L2 tolerance, not a substitute for it)")
        # A gross encoding error (wrong dtype, wrong scaling) would show up here
        # as orders of magnitude, so gate on that rather than on the exact figure.
        if worst > 1e-6:
            print("  FAIL: round-trip deviation far above the float32 expectation")
            ok = False
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--npz", help="lsdn_reference.npz from extract_lsdn_reference.py")
    src.add_argument("--data-dir", help="v3 data dir; re-runs the extraction + C1-C7")
    ap.add_argument("--nc", default="data/crn_parameters.nc",
                    help="crn_parameters.nc to append to (in place)")
    ap.add_argument("--verify", action="store_true",
                    help="re-read the packed file and reproduce C1/C2/C3")
    args = ap.parse_args()

    arrays, prov = load_arrays(args)
    if prov.get("sha256"):
        print("source .mat sha256 recorded in the nc provenance:")
        for fn, h in sorted(prov["sha256"].items()):
            print(f"  {h}  {fn}")

    size_before = os.path.getsize(args.nc)
    snap, gattrs = snapshot_existing(args.nc)
    print(f"\n{args.nc}: {size_before} B before, "
          f"{len(snap)} existing variables snapshotted")

    ds = nc.Dataset(args.nc, "a")
    try:
        already = "lsdn_b_Be10" in ds.variables
        if already:
            print("LSDn variables already present; overwriting values in place.")
        written = pack(ds, arrays, prov)
    finally:
        ds.close()

    size_after = os.path.getsize(args.nc)
    print(f"wrote {len(written)} variables: {', '.join(written)}")
    print(f"size: {size_before} B -> {size_after} B "
          f"(added {size_after - size_before} B, "
          f"{(size_after - size_before) / 1048576:.3f} MiB; "
          f"total {size_after / 1048576:.3f} MiB)")

    bad = verify_unchanged(args.nc, snap, gattrs)
    if bad:
        print("\nPRE-EXISTING DATA WAS MODIFIED -- this append is not "
              "non-destructive:", file=sys.stderr)
        for b in bad:
            print("  " + b, file=sys.stderr)
        return 1
    print(f"non-destructive: all {len(snap)} pre-existing variables byte-identical, "
          f"attributes unchanged")

    if args.verify and not verify(args.nc, arrays):
        print("\nVERIFY FAILED", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
