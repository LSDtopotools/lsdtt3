#!/usr/bin/env python3
"""Extract and verify the LSDn (Lifton-Sato-Dunai 2014) reference data.

lsdtt3 implements LSDn production scaling by interpolating the precomputed
`b`/`m` scaling grids that Balco's v3 online calculator itself interpolates
(`get_LSDnSF.m`), rather than by porting the Sato/PARMA spectral integration.
This script is the one-off extraction step: it reads the v3 reference `.mat`
files, verifies them against a fixed set of checks, and writes an intermediate
`.npz` plus a human-readable dump that the nc packer (Phase 1) consumes.

This is ORCHESTRATION ONLY.  No production physics lives here and none ever
will -- the numerics are C++.  The script exists so that the numbers that end
up in `data/crn_parameters.nc` are traceable to a checksummed source file and
so that a silent transcription or axis-orientation error is caught once, here,
instead of surfacing as a plausible-but-wrong production surface later.

LICENCE NOTE.  The v3 reference implementation is GPL-v2.  What is extracted
here is DATA (grid axes, coefficient tables, paleomagnetic and solar records) --
facts, not expression -- and no MATLAB source text is copied into this
repository.  The `.mat` files themselves are read from an untracked directory
outside the repo and are never vendored.

Source files (all in the v3 `data/` directory):
    LSDn2015_grid.mat        sfgrid.Rc, sfgrid.p        -- full 380x201 meshes
    LSDn2015_N10quartz.mat   sfN10quartz.b, .m          -- Be-10 in quartz
    LSDn2015_N26quartz.mat   sfN26quartz.b, .m          -- Al-26 in quartz
    LSDn2015_N14quartz.mat   sfN14quartz.b, .m          -- C-14 in quartz
    LiftonS2015.mat          SPhi (91), SPhiInf         -- solar modulation
    PavonDip.mat             m.MM01, m.t1min, m.t1max, m.latr_pp, m.lonr_pp,
                             m.MM02, m.t2min, m.t2max, m.t2, m.tzero

Checks (all reported PASS/FAIL; a failure is a hard stop, not something to be
tuned away):
    C1  axis regularity and extent      C5  MM01[0] ~ 1, tzero = 2010
    C2  Be-10 b/m ranges                C6  SPhi[0], SPhiInf, record lengths
    C3  SLHL Be-10 value, OFF-NODE in p C7  Pavon record time bounds
    C4  Rc < 1 GV columns flat

C3 is the important one: p = 1013.25 hPa lies between the 1012 and 1014 nodes,
so it is a genuine bilinear evaluation and only reproduces 4.8809 if the axis
orientation (p slow, Rc fast) is right.  A transposed grid gives a plausible
surface and a wrong answer.

Usage (pixi environment; see pixi.toml):
    pixi run python python_helpers/extract_lsdn_reference.py \
        --data-dir /mnt/c/LSDTopoTools/crn_reference/v3/data \
        --out-dir  /path/to/scratch/lsdn

Outputs (written to --out-dir, which must NOT be inside the repository):
    lsdn_reference.npz        arrays + a provenance JSON blob
    lsdn_reference_dump.txt   the same content in human-readable form

Exit status is non-zero if any check fails.
"""
import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
import scipy.io as sio

# ---------------------------------------------------------------------------
# Source layout.  `nuclide -> (mat file, struct name)`; the struct name is
# verified against the file rather than trusted.
# ---------------------------------------------------------------------------
GRID_MAT = "LSDn2015_grid.mat"
SOLAR_MAT = "LiftonS2015.mat"
PAVON_MAT = "PavonDip.mat"
NUCLIDE_MATS = {
    "Be10": ("LSDn2015_N10quartz.mat", "sfN10quartz"),
    "Al26": ("LSDn2015_N26quartz.mat", "sfN26quartz"),
    "C14": ("LSDn2015_N14quartz.mat", "sfN14quartz"),
}
SOURCE_FILES = [GRID_MAT, SOLAR_MAT, PAVON_MAT] + [v[0] for v in NUCLIDE_MATS.values()]

# Expected axis definition (C1).
RC_MIN, RC_STEP, RC_N = 0.0, 0.1, 201
P_MIN, P_STEP, P_N = 300.0, 2.0, 380

# Expected reference values.  These are independently measured numbers being
# used as transcription checks -- if one disagrees, report it, do not widen it.
C2_B_RANGE = (1.369, 1673.87)
C2_M_RANGE = (-1.0688, -4.96e-5)
C3_P_HPA = 1013.25          # deliberately off-node (between 1012 and 1014)
C3_RC_GV = 0.0
C3_S_MV = 624.5718          # SPhi[0]
C3_EXPECTED = 4.8809
C6_SPHI0 = 624.5718
C6_SPHIINF = 462.5406
C6_N_SHORT = 91
C5_TZERO = 2010
# C7 checks the step-BOUND arrays (m.t1min/t1max, m.t2min/t2max).  PavonDip.mat
# also carries the step-CENTRE arrays m.t1/m.t2, whose values are 13860 / 14060 /
# 2.00074e7 -- close enough to be mistaken for these and the source of an earlier
# transcription error.  Bounds are what the duration-weighted averaging uses, so
# bounds are what is checked and packed; m.t2 is packed too, but only as the
# get_DipRc.m branch threshold (see NOTE_T2 below).
C7_T1MAX_LAST = 13960.0
C7_T2MIN_FIRST = 13960.0
C7_T2MAX_LAST = 2.02044e7

# Tolerances.  The expected values above are quoted to 5-6 significant figures
# in the research note, so they are compared at the precision they were quoted
# at; they are not free parameters.
TOL_QUOTED_REL = 1e-4       # C2/C3/C6 against 5-6 s.f. quoted values
TOL_QUOTED_REL_LOOSE = 1e-3 # C2 lower bounds, quoted to 4 s.f. only
TOL_AXIS_ABS = 1e-9         # C1 axis regularity
TOL_FLAT_REL = 1e-12        # C4 flatness of the sub-1 GV columns
TOL_TIME_REL = 1e-3         # C7 record bounds


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def sha256_of(path):
    """sha256 of a file, read in chunks (the grids are ~1.1 MB each)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_struct(path, expect_name):
    """Load a .mat holding a single struct; return (struct, actual name)."""
    d = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    keys = [k for k in d if not k.startswith("__")]
    if len(keys) != 1:
        raise ValueError(f"{path}: expected one variable, found {keys}")
    if keys[0] != expect_name:
        print(f"  NOTE: {os.path.basename(path)} holds '{keys[0]}', "
              f"plan says '{expect_name}'")
    return d[keys[0]], keys[0]


def axis_from_mesh(mesh, along):
    """Pull a 1-D axis out of a full 2-D mesh.

    `along` is the axis the values VARY along (0 = down rows, 1 = across
    columns).  Asserts the mesh is genuinely constant in the other direction,
    which is what makes the orientation determination trustworthy rather than
    an assumption.
    """
    a = np.asarray(mesh, dtype=np.float64)
    if a.ndim != 2:
        raise ValueError(f"expected a 2-D mesh, got shape {a.shape}")
    other = 1 - along
    spread = np.ptp(a, axis=other)
    if np.max(spread) > TOL_AXIS_ABS:
        raise ValueError(f"mesh is not constant along axis {other} "
                         f"(max spread {np.max(spread):g}) -- orientation wrong?")
    return a[:, 0].copy() if along == 0 else a[0, :].copy()


def axis_report(name, ax, x0, step, n):
    """Regularity/extent report for one axis; returns (ok, detail string)."""
    dif = np.diff(ax)
    max_step_err = float(np.max(np.abs(dif - step))) if dif.size else 0.0
    ok = (ax.size == n
          and abs(ax[0] - x0) <= TOL_AXIS_ABS
          and abs(ax[-1] - (x0 + step * (n - 1))) <= TOL_AXIS_ABS
          and max_step_err <= TOL_AXIS_ABS)
    detail = (f"{name}: n={ax.size} (exp {n}), {ax[0]:g}..{ax[-1]:g} "
              f"(exp {x0:g}..{x0 + step * (n - 1):g}), "
              f"max |d-{step:g}| = {max_step_err:.3e}")
    return ok, detail


def bilinear_regular(grid, x_ax, y_ax, x, y):
    """Bilinear interpolation on a regular grid, grid[iy, ix].

    Mirrors what the C++ `bilinear_regular()` will do (index by arithmetic,
    no search), so that C3 exercises the same arithmetic path the runtime
    will take.  Coordinates are clamped to the grid, as the runtime clamps.
    """
    def frac(ax, v):
        step = (ax[-1] - ax[0]) / (ax.size - 1)
        t = (v - ax[0]) / step
        t = min(max(t, 0.0), ax.size - 1.0)
        i = min(int(np.floor(t)), ax.size - 2)
        return i, t - i

    ix, fx = frac(x_ax, x)
    iy, fy = frac(y_ax, y)
    g = np.asarray(grid, dtype=np.float64)
    return float((1 - fy) * ((1 - fx) * g[iy, ix] + fx * g[iy, ix + 1])
                 + fy * ((1 - fx) * g[iy + 1, ix] + fx * g[iy + 1, ix + 1]))


def close(actual, expected, rel):
    return abs(actual - expected) <= rel * max(abs(expected), 1e-300)


class CheckLog:
    """Accumulates PASS/FAIL lines so the summary table can be printed twice
    (stdout and the dump file) and so the exit status is derived, not set."""

    def __init__(self):
        self.rows = []

    def add(self, tag, ok, what, actual, expected):
        self.rows.append((tag, bool(ok), what, actual, expected))
        print(f"  [{'PASS' if ok else 'FAIL'}] {tag}  {what}")
        print(f"         actual   : {actual}")
        print(f"         expected : {expected}")

    @property
    def all_ok(self):
        return all(r[1] for r in self.rows)

    def table(self):
        w = max(len(r[2]) for r in self.rows)
        out = [f"{'chk':<4} {'result':<7} {'what':<{w}}  actual | expected", ]
        out.append("-" * (len(out[0]) + 8))
        for tag, ok, what, actual, expected in self.rows:
            out.append(f"{tag:<4} {'PASS' if ok else 'FAIL':<7} {what:<{w}}  "
                       f"{actual} | {expected}")
        return "\n".join(out)


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------
def extract(data_dir, log):
    """Read every source .mat, run C1-C7, return (arrays dict, provenance)."""
    prov = {"extracted_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "data_dir": os.path.abspath(data_dir),
            "sha256": {}, "field_names": {}, "shapes": {}}

    print("sha256 of source files:")
    for fn in sorted(SOURCE_FILES):
        p = os.path.join(data_dir, fn)
        if not os.path.isfile(p):
            raise FileNotFoundError(p)
        prov["sha256"][fn] = sha256_of(p)
        print(f"  {prov['sha256'][fn]}  {fn}")

    arrays = {}

    # -- axes -----------------------------------------------------------------
    sfgrid, gname = load_struct(os.path.join(data_dir, GRID_MAT), "sfgrid")
    prov["field_names"][GRID_MAT] = [gname + "." + f for f in sfgrid._fieldnames]
    rc_mesh, p_mesh = np.asarray(sfgrid.Rc), np.asarray(sfgrid.p)
    prov["shapes"]["sfgrid.Rc"] = list(rc_mesh.shape)
    prov["shapes"]["sfgrid.p"] = list(p_mesh.shape)
    print(f"\nmeshes: sfgrid.Rc {rc_mesh.shape} {rc_mesh.dtype}, "
          f"sfgrid.p {p_mesh.shape} {p_mesh.dtype}")

    # Determine orientation from the data rather than assuming it: Rc must
    # vary across columns and p down rows for the (p slow, Rc fast) layout.
    rc_ax = axis_from_mesh(rc_mesh, along=1)
    p_ax = axis_from_mesh(p_mesh, along=0)
    print("orientation: Rc varies across COLUMNS, p varies down ROWS "
          "-> storage (p, Rc), p slow")

    ok_rc, det_rc = axis_report("Rc [GV]", rc_ax, RC_MIN, RC_STEP, RC_N)
    ok_p, det_p = axis_report("p [hPa]", p_ax, P_MIN, P_STEP, P_N)
    print("\nchecks:")
    log.add("C1", ok_rc and ok_p, "axes regular and at expected extent",
            det_rc + "; " + det_p,
            f"Rc {RC_MIN}:{RC_STEP}:{RC_MIN + RC_STEP * (RC_N - 1)} ({RC_N} pts), "
            f"p {P_MIN}:{P_STEP}:{P_MIN + P_STEP * (P_N - 1)} ({P_N} pts), both exact")
    arrays["lsdn_rc"] = rc_ax
    arrays["lsdn_p"] = p_ax

    # -- nuclide grids --------------------------------------------------------
    for nuc, (fn, sname) in NUCLIDE_MATS.items():
        st, actual = load_struct(os.path.join(data_dir, fn), sname)
        prov["field_names"][fn] = [actual + "." + f for f in st._fieldnames]
        b = np.asarray(st.b, dtype=np.float64)
        m = np.asarray(st.m, dtype=np.float64)
        if b.shape != (P_N, RC_N) or m.shape != (P_N, RC_N):
            raise ValueError(f"{fn}: b/m shapes {b.shape}/{m.shape}, "
                             f"expected ({P_N}, {RC_N})")
        prov["shapes"][f"{actual}.b"] = list(b.shape)
        prov["shapes"][f"{actual}.m"] = list(m.shape)
        arrays[f"lsdn_b_{nuc}"] = b
        arrays[f"lsdn_m_{nuc}"] = m

    b10 = arrays["lsdn_b_Be10"]
    m10 = arrays["lsdn_m_Be10"]

    # C2 -- transcription check on the Be-10 grid ranges.
    got = (float(b10.min()), float(b10.max()), float(m10.min()), float(m10.max()))
    ok_c2 = (close(got[0], C2_B_RANGE[0], TOL_QUOTED_REL_LOOSE)
             and close(got[1], C2_B_RANGE[1], TOL_QUOTED_REL)
             and close(got[2], C2_M_RANGE[0], TOL_QUOTED_REL_LOOSE)
             and close(got[3], C2_M_RANGE[1], TOL_QUOTED_REL_LOOSE))
    log.add("C2", ok_c2, "Be-10 b and m grid ranges",
            f"b {got[0]:.6g}..{got[1]:.6g}, m {got[2]:.6g}..{got[3]:.6g}",
            f"b {C2_B_RANGE[0]}..{C2_B_RANGE[1]}, m {C2_M_RANGE[0]}..{C2_M_RANGE[1]}")

    # C3 -- SLHL Be-10 value, evaluated OFF-NODE in p.  This is the axis
    # orientation check as much as it is a value check.
    b_slhl = bilinear_regular(b10, rc_ax, p_ax, C3_RC_GV, C3_P_HPA)
    m_slhl = bilinear_regular(m10, rc_ax, p_ax, C3_RC_GV, C3_P_HPA)
    v_slhl = b_slhl + C3_S_MV * m_slhl
    log.add("C3", close(v_slhl, C3_EXPECTED, TOL_QUOTED_REL),
            f"b + {C3_S_MV}*m at Rc={C3_RC_GV} GV, p={C3_P_HPA} hPa (off-node)",
            f"{v_slhl:.6f}  (b={b_slhl:.6f}, m={m_slhl:.8f})",
            f"{C3_EXPECTED}")

    # C4 -- the generator (NeutronsXS.m) floors Rc at 1 GV, so every column
    # below 1 GV should equal the 1 GV column.  If it does not, the runtime
    # needs an explicit floor.  Checked for every nuclide, since the runtime
    # floor decision applies to all of them.
    i1 = int(round((1.0 - RC_MIN) / RC_STEP))
    dev = []
    for nuc in NUCLIDE_MATS:
        for name in ("b", "m"):
            g = arrays[f"lsdn_{name}_{nuc}"]
            ref = g[:, i1:i1 + 1]
            d = np.max(np.abs(g[:, :i1 + 1] - ref) / np.maximum(np.abs(ref), 1e-300))
            dev.append((f"{nuc}.{name}", float(d)))
    max_dev = max(d for _, d in dev)
    log.add("C4", max_dev <= TOL_FLAT_REL,
            f"columns Rc <= 1 GV flat (cols 0..{i1}), all nuclides",
            "max rel dev " + ", ".join(f"{n}={d:.3e}" for n, d in dev),
            f"<= {TOL_FLAT_REL:g} (flat; else add a runtime 1 GV floor)")
    prov["c4_flat_columns"] = {"i_1GV": i1, "max_rel_dev": max_dev}

    # -- solar record ---------------------------------------------------------
    sol = sio.loadmat(os.path.join(data_dir, SOLAR_MAT),
                      squeeze_me=True, struct_as_record=False)
    prov["field_names"][SOLAR_MAT] = sorted(k for k in sol if not k.startswith("__"))
    sphi = np.asarray(sol["SPhi"], dtype=np.float64)
    sphi_inf = float(sol["SPhiInf"])
    prov["shapes"]["SPhi"] = list(sphi.shape)
    arrays["lsdn_sphi"] = sphi
    arrays["lsdn_sphi_inf"] = np.float64(sphi_inf)

    # -- Pavon paleomagnetic record ------------------------------------------
    pv, pvname = load_struct(os.path.join(data_dir, PAVON_MAT), "m")
    prov["field_names"][PAVON_MAT] = [pvname + "." + f for f in pv._fieldnames]
    for fld in ("MM01", "t1min", "t1max", "latr_pp", "lonr_pp",
                "MM02", "t2min", "t2max", "t2", "tzero"):
        if fld not in pv._fieldnames:
            raise ValueError(f"{PAVON_MAT}: missing field '{pvname}.{fld}' "
                             f"(present: {pv._fieldnames})")
    # The .mat stores several of the short-record arrays as uint16; promote to
    # f8 here so the packer never has to think about integer time axes.
    mm01 = np.asarray(pv.MM01, dtype=np.float64)
    t1min = np.asarray(pv.t1min, dtype=np.float64)
    t1max = np.asarray(pv.t1max, dtype=np.float64)
    mm02 = np.asarray(pv.MM02, dtype=np.float64)
    t2min = np.asarray(pv.t2min, dtype=np.float64)
    t2max = np.asarray(pv.t2max, dtype=np.float64)
    # Long-record step CENTRES.  Extracted because get_DipRc.m uses m.t2(end)
    # as the "exposure longer than the long record" branch threshold.  It is a
    # threshold, never an averaging weight -- the averaging uses the bounds.
    t2c_arr = np.asarray(pv.t2, dtype=np.float64)
    tzero = float(pv.tzero)
    for nm, a in (("MM01", mm01), ("t1min", t1min), ("t1max", t1max),
                  ("latr_pp", pv.latr_pp), ("lonr_pp", pv.lonr_pp),
                  ("MM02", mm02), ("t2min", t2min), ("t2max", t2max),
                  ("t2", t2c_arr)):
        prov["shapes"][f"m.{nm}"] = list(np.shape(a))
    arrays["pavon_mm01"] = mm01
    arrays["pavon_t1min"] = t1min
    arrays["pavon_t1max"] = t1max
    arrays["pavon_pole_lat"] = np.asarray(pv.latr_pp, dtype=np.float64)
    arrays["pavon_pole_lon"] = np.asarray(pv.lonr_pp, dtype=np.float64)
    arrays["pavon_mm02"] = mm02
    arrays["pavon_t2min"] = t2min
    arrays["pavon_t2max"] = t2max
    arrays["pavon_t2"] = t2c_arr
    arrays["pavon_tzero"] = np.float64(tzero)

    # C5 -- the short record starts at the present-day field, at tzero.
    log.add("C5", close(float(mm01[0]), 1.0, 1e-3) and tzero == C5_TZERO,
            "MM01[0] ~ 1 and tzero = 2010",
            f"MM01[0] = {mm01[0]:.6f}, tzero = {tzero:g}",
            f"1 (+-1e-3), {C5_TZERO}")

    # C6 -- solar record pins and the SPhi/MM01 length coupling that
    # get_DipRc.m relies on ("assumes SPhi has same numel as MM01").
    ok_c6 = (close(float(sphi[0]), C6_SPHI0, TOL_QUOTED_REL)
             and close(sphi_inf, C6_SPHIINF, TOL_QUOTED_REL)
             and sphi.size == mm01.size == C6_N_SHORT)
    log.add("C6", ok_c6, "SPhi[0], SPhiInf, len(SPhi) == len(MM01) == 91",
            f"SPhi[0] = {sphi[0]:.6f}, SPhiInf = {sphi_inf:.6f}, "
            f"len(SPhi) = {sphi.size}, len(MM01) = {mm01.size}",
            f"{C6_SPHI0}, {C6_SPHIINF}, 91, 91")

    # C7 -- record time bounds, plus the exact contiguity that get_DipRc.m
    # relies on when it overwrites m.t2min(1) = temptmax(end): the short and
    # long records must abut with zero gap, and the long record's own steps
    # must be gapless.  The step-CENTRE arrays are reported alongside so that
    # a future reader can see at a glance which array a number came from.
    t1c = np.asarray(pv.t1, dtype=np.float64) if "t1" in pv._fieldnames else None
    t2c = np.asarray(pv.t2, dtype=np.float64) if "t2" in pv._fieldnames else None
    abut = float(t1max[-1] - t2min[0])
    gap2 = float(np.max(np.abs(t2min[1:] - t2max[:-1])))
    ok_c7 = (close(float(t1max[-1]), C7_T1MAX_LAST, TOL_TIME_REL)
             and close(float(t2min[0]), C7_T2MIN_FIRST, TOL_TIME_REL)
             and close(float(t2max[-1]), C7_T2MAX_LAST, TOL_TIME_REL)
             and abut == 0.0 and gap2 == 0.0)
    centres = ("" if t1c is None or t2c is None else
               f"; centres t1[-1] = {t1c[-1]:g}, t2[0] = {t2c[0]:g}, "
               f"t2[-1] = {t2c[-1]:g}")
    log.add("C7", ok_c7, "Pavon record bounds (yr) + exact contiguity",
            f"t1max[-1] = {t1max[-1]:g}, t2min[0] = {t2min[0]:g}, "
            f"t2max[-1] = {t2max[-1]:.6g}, t1max[-1]-t2min[0] = {abut:g}, "
            f"max|t2min[1:]-t2max[:-1]| = {gap2:g}{centres}",
            f"{C7_T1MAX_LAST:g}, {C7_T2MIN_FIRST:g}, {C7_T2MAX_LAST:.6g}, 0, 0")

    prov["n_pavon_short"] = int(mm01.size)
    prov["n_pavon_long"] = int(mm02.size)
    return arrays, prov


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def emit_reference_csv(path, arrays, prov):
    """Write the FROZEN L2 grid reference (lsdn_implementation_plan Sec5, L2).

    The values are evaluated from the float64 `.mat` grids with scipy's
    RegularGridInterpolator (method='linear'), which is the same 4-point formula
    MATLAB `interp2` linear uses.  lsdtt3 must reproduce these from the FLOAT32
    grids packed into `crn_parameters.nc`, so the residual is the float32
    round-trip, not a formula difference.  At grid NODES the formula is the
    identity, so those rows are a pure data-fidelity test.

    Regenerating this file is a MANUAL, REVIEWABLE event: it re-derives a frozen
    reference from the untracked v3 `.mat` files and must be signed off before
    it is committed.
    """
    from scipy.interpolate import RegularGridInterpolator

    rc = np.asarray(arrays["lsdn_rc"], dtype=np.float64)
    p = np.asarray(arrays["lsdn_p"], dtype=np.float64)
    sphi0 = float(np.asarray(arrays["lsdn_sphi"])[0])
    rc_max = float(rc[-1])

    # -- the case list --------------------------------------------------------
    # (a) the 15 CRONUS fixture samples: 5 latitudes x 3 elevations, standard
    #     atmosphere, present-day field and paleopole, longitude 0.
    cases = []  # (case, Rc_GV, p_hPa, S_MV)
    mm01_0 = float(np.asarray(arrays["pavon_mm01"])[0])
    pole_lat = float(np.asarray(arrays["pavon_pole_lat"])[0])
    pole_lon = float(np.asarray(arrays["pavon_pole_lon"])[0])
    for lat_deg in (0.0, 20.0, 40.0, 60.0, 70.0):
        for elev_m in (250.0, 1500.0, 3000.0):
            # Balco's barometric equation with SLHL values (lsdtt3's
            # atmospheric_pressure_standard).
            press = 1013.25 * np.exp(
                (-0.03417 / 0.0065)
                * (np.log(288.15) - np.log(288.15 - 0.0065 * elev_m)))
            lat = np.radians(lat_deg)
            cos_theta = (np.cos(lat) * np.cos(pole_lat) * np.cos(0.0 - pole_lon)
                         + np.sin(lat) * np.sin(pole_lat))
            gm_lat = abs(0.5 * np.pi - np.arccos(min(1.0, max(-1.0, cos_theta))))
            c = np.cos(gm_lat)
            poly = mm01_0 * (-448.004 * c**6 + 1189.18 * c**5 - 1152.15 * c**4
                             + 522.061 * c**3 - 103.241 * c**2 + 6.89901 * c)
            cases.append(("fixture", abs(float(poly)), float(press), sphi0))

    # (b) exact grid nodes: the four corners plus interior nodes.
    for irc, ip in [(0, 0), (0, P_N - 1), (RC_N - 1, 0), (RC_N - 1, P_N - 1),
                    (1, 1), (50, 100), (100, 190), (150, 250), (200, 379),
                    (10, 5), (75, 300), (199, 378)]:
        cases.append(("node", float(rc[irc]), float(p[ip]), sphi0))

    # (c) deliberate off-node points: cell midpoints and points 1e-3 inside
    #     each edge of the grid.
    for irc, ip in [(0, 0), (50, 100), (100, 190), (150, 250), (199, 378)]:
        cases.append(("midcell", float(rc[irc] + 0.05), float(p[ip] + 1.0), sphi0))
    cases += [
        ("edge_in", float(rc[0] + 1e-3), float(p[0] + 1e-3), sphi0),
        ("edge_in", float(rc[-1] - 1e-3), float(p[0] + 1e-3), sphi0),
        ("edge_in", float(rc[0] + 1e-3), float(p[-1] - 1e-3), sphi0),
        ("edge_in", float(rc[-1] - 1e-3), float(p[-1] - 1e-3), sphi0),
        ("offnode", 3.7317, 913.7, sphi0),
        ("offnode", 11.2481, 657.3, sphi0),
        ("offnode", 7.5, 1013.25, 462.54064664694187),
    ]

    # (d) the SLHL point - off-node in p, which is the ORIENTATION check: it
    #     only comes out at 4.8809 for Be-10 if p is slow and Rc fast.
    cases.append(("slhl", 0.0, 1013.25, sphi0))

    # (e) the Rc clamp (25 GV -> 20 GV, as the runtime and v3 both clamp) and a
    #     point below the generator's 1 GV floor.
    cases.append(("rc_clamp", 25.0, 850.0, sphi0))
    cases.append(("rc_sub1", 0.5, 850.0, sphi0))

    interps = {}
    for nuc in NUCLIDE_MATS:
        interps[nuc] = (
            RegularGridInterpolator((p, rc), arrays[f"lsdn_b_{nuc}"],
                                    method="linear", bounds_error=True),
            RegularGridInterpolator((p, rc), arrays[f"lsdn_m_{nuc}"],
                                    method="linear", bounds_error=True))

    with open(path, "w") as fh:
        fh.write(
            "# LSDn v3 GRID REFERENCE (frozen) - lsdn_implementation_plan Sec5 L2\n"
            "# Evaluated from the FLOAT64 v3 .mat grids with scipy\n"
            "# RegularGridInterpolator(method='linear'); lsdtt3 must reproduce\n"
            "# these from the float32 grids packed in data/crn_parameters.nc.\n"
            "# Regenerating this file is a MANUAL, REVIEWABLE event and needs\n"
            "# sign-off by Simon Mudd before it is committed.\n"
            f"# generated_utc: {prov['extracted_utc']}\n"
            f"# generating_script: python_helpers/extract_lsdn_reference.py "
            "--emit-reference\n")
        for fn, h in sorted(prov["sha256"].items()):
            fh.write(f"# source_sha256 {h}  {fn}\n")
        fh.write("# Rc_GV is the REQUESTED rigidity; values are evaluated at\n"
                 "# min(Rc, 20) because both v3 and lsdtt3 clamp there.\n")
        fh.write("nuclide,case,Rc_GV,p_hPa,S_MV,value\n")
        for nuc, (bi, mi) in interps.items():
            for case, rcv, pv_, s in cases:
                rc_eval = min(rcv, rc_max)
                bval = float(bi([[pv_, rc_eval]])[0])
                mval = float(mi([[pv_, rc_eval]])[0])
                fh.write(f"{nuc},{case},{rcv!r},{pv_!r},{s!r},"
                         f"{bval + s * mval!r}\n")
    print(f"wrote {path} ({len(cases) * len(interps)} rows)")


def write_dump(path, arrays, prov, log):
    """Human-readable companion to the .npz -- what was read, from what, and
    whether it checked out."""
    lines = ["LSDn reference extraction (Lifton-Sato-Dunai 2014 / Balco v3)",
             f"extracted    : {prov['extracted_utc']}",
             f"source dir   : {prov['data_dir']}",
             "",
             "sha256 of source .mat files",
             "---------------------------"]
    for fn, h in sorted(prov["sha256"].items()):
        lines.append(f"{h}  {fn}")

    lines += ["", "field names as found", "--------------------"]
    for fn, flds in sorted(prov["field_names"].items()):
        lines.append(f"{fn}: {', '.join(flds)}")

    lines += ["", "shapes", "------"]
    for k, v in prov["shapes"].items():
        lines.append(f"{k}: {tuple(v)}")

    lines += ["", "checks", "------", log.table(),
              "", f"OVERALL: {'PASS' if log.all_ok else 'FAIL'}",
              "", "extracted arrays", "----------------"]
    for k in sorted(arrays):
        a = np.asarray(arrays[k])
        if a.ndim == 0:
            lines.append(f"{k}: scalar {float(a):.10g}")
        else:
            lines.append(f"{k}: shape {a.shape} dtype {a.dtype} "
                         f"range {float(a.min()):.6g} .. {float(a.max()):.6g}")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default="/mnt/c/LSDTopoTools/crn_reference/v3/data",
                    help="directory holding the v3 reference .mat files")
    ap.add_argument("--emit-reference", metavar="CSV",
                    help="ALSO write the frozen L2 grid reference CSV to this "
                         "path (a manual, reviewable regeneration event; the "
                         "normal path is to leave the committed file alone)")
    ap.add_argument("--out-dir", required=True,
                    help="scratch directory for the .npz and the dump "
                         "(deliberately NOT inside the repository)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    log = CheckLog()
    arrays, prov = extract(args.data_dir, log)

    npz_path = os.path.join(args.out_dir, "lsdn_reference.npz")
    dump_path = os.path.join(args.out_dir, "lsdn_reference_dump.txt")
    np.savez_compressed(npz_path, provenance_json=np.array(json.dumps(prov, indent=2)),
                        **arrays)
    write_dump(dump_path, arrays, prov, log)

    if args.emit_reference:
        emit_reference_csv(args.emit_reference, arrays, prov)

    print("\n" + log.table())
    print(f"\nwrote {npz_path} ({os.path.getsize(npz_path)} B)")
    print(f"wrote {dump_path}")
    if not log.all_ok:
        print("\nONE OR MORE CHECKS FAILED -- stop and report; do not relax the "
              "expected values.", file=sys.stderr)
        return 1
    print("all checks PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
