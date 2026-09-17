#!/usr/bin/env python3
"""Regenerate tests/data/expinf_95_02_reference.csv.

MAINTAINER TOOL. You do not need this to build or test lsdtt3, and lsdtt3 has
no Fortran dependency of any kind: the test suite reads the committed CSV.

The CSV is a frozen trace of TOPMODEL 95.02's subroutine EXPINF. This script
holds that subroutine verbatim (see FORTRAN_SOURCE below, transcribed from
refs/original_code_fortran/TMOD9502.f in https://github.com/NOAA-OWP/topmodel),
writes it to a temporary directory, compiles it, runs it over the reference
storms, and writes the CSV. Nothing is left behind.

Run it only if the storms or parameters in the fixture need to change. Requires
gfortran at that moment and at no other time.

    python3 tests/scripts/regenerate_expinf_fixture.py \
        --out tests/data/expinf_95_02_reference.csv

    # or, to confirm the committed fixture still matches its source:
    python3 tests/scripts/regenerate_expinf_fixture.py \
        --out tests/data/expinf_95_02_reference.csv --check

TWO CHANGES TO THE ORIGINAL, both mechanical and both documented here because
they are the only places this differs from TMOD9502.f:

  * COMMON blocks become explicit arguments. The original reads XK0, HF, DTH,
    SZM and DT out of COMMON/PARAM/ and COMMON/FLOW/, and there is no model here
    to fill them.

  * TP, CONST and F become arguments rather than locals. The original RELIES on
    all three surviving between calls while ponding persists: TP and CONST are
    read by the Newton solve on a later call, and F carries the "initial
    estimate for next time step". F77 compilers commonly gave locals static
    storage; the standard does not require it. Making them arguments preserves
    the behaviour the original depended on instead of depending on the compiler
    to supply it.

    This matters. A first version carried TP and CONST but dropped the
    F = F + DF line, and the resulting trace disagreed with the port by 3.3e-6 m.

The arithmetic is otherwise untouched, including the ten-term series with its
POW and running factorial, the twice-evaluated DLOG(FC), and the dimensionally
odd final test F .LT. CUMF + RINT.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

# Reference storms. `dt4` exists because dt = 4 h is where the original's
# `F < CUMF + RINT` test differs from the dimensionally correct form.
STORMS = {
    "demo": dict(xk0=1.0, hf=0.02, dth=0.1, szm=0.032, dt=1.0,
                 rain=[0.002 * (i + 1) for i in range(12)] + [0.08] * 20 +
                      [0.08 / (i + 2) for i in range(12)]),
    "lowk": dict(xk0=0.02, hf=0.05, dth=0.3, szm=0.05, dt=1.0,
                 rain=[0.001] * 5 + [0.03] * 25 + [0.005] * 10),
    "dt4":  dict(xk0=1.0, hf=0.02, dth=0.1, szm=0.032, dt=4.0,
                 rain=[0.05] * 30),
}

HEADER = [
    "# Frozen trace of TOPMODEL 95.02 subroutine EXPINF (Beven, HSJ, 1984).",
    "# Source: refs/original_code_fortran/TMOD9502.f in NOAA-OWP/topmodel.",
    "#",
    "# THIS FILE IS THE ORACLE. lsdtt3 has no Fortran dependency: the test",
    "# suite reads this CSV. Regenerate with",
    "#   python3 tests/scripts/regenerate_expinf_fixture.py --out <this file>",
    "# which carries the subroutine verbatim and needs gfortran only then.",
]

FORTRAN_SOURCE = r"""
      PROGRAM EXPINF_DRIVER
      IMPLICIT NONE
      DOUBLE PRECISION XK0, HF, DTH, SZM, DT
      DOUBLE PRECISION RINT(10000)
      DOUBLE PRECISION DF, CUMF, TP, CONST, F
      INTEGER NSTEP, IT, IROF, I

      READ(5,*) XK0, HF, DTH, SZM, DT, NSTEP
      READ(5,*) (RINT(I), I = 1, NSTEP)

      IROF  = 0
      CUMF  = 0.0D0
      DF    = 0.0D0
      TP    = 0.0D0
      CONST = 0.0D0
      F     = 0.0D0

      DO 100 IT = 1, NSTEP
         IF (RINT(IT) .GT. 0.0D0) THEN
            CALL EXPINF(IROF, IT, RINT(IT), DF, CUMF, TP, CONST, F,
     1                  XK0, HF, DTH, SZM, DT)
         ELSE
            IROF = 0
            DF   = 0.0D0
         ENDIF
         WRITE(6,900) IT, IROF, DF, CUMF
  900    FORMAT(I6,1X,I2,1X,E22.15,1X,E22.15)
  100 CONTINUE

      END

C ----------------------------------------------------------------------
      SUBROUTINE EXPINF(IROF, IT, RINT, DF, CUMF, TP, CONST, F,
     1                  XK0, HF, DTH, SZM, DT)
      IMPLICIT NONE
      INTEGER IROF, IT
      DOUBLE PRECISION RINT, DF, CUMF, TP, CONST, F
      DOUBLE PRECISION XK0, HF, DTH, SZM, DT

      DOUBLE PRECISION SUM, FC, FUNC, CD, SZF, XKF
      DOUBLE PRECISION E, F1, F2, R2, FAC, ADD, DFUNC
      INTEGER I, J
      DATA E /0.00001D0/

C  Note that HF and DTH only appear in product CD
      CD  = HF * DTH
      SZF = 1.0D0 / SZM
      XKF = XK0
      IF (IROF .EQ. 1) GO TO 10
C  PONDING HAS ALREADY OCCURRED - GO TO EXCESS CALCULATION
      IF (CUMF .EQ. 0.0D0) GOTO 7
C  FIRST TIME STEP, OVERFLOW IF CUMF=0, GO DIRECT TO F2 CALCULATION
C  INITIAL ESTIMATE OF TIME TO PONDING
      F1 = CUMF
      R2 = -XKF * SZF * (CD + F1) / (1.0D0 - DEXP(SZF * F1))
      IF (R2 .LT. RINT) THEN
C  PONDING STARTS AT BEGINNING OF TIME STEP
         TP   = (IT - 1.0D0) * DT
         IROF = 1
         F    = CUMF
         GO TO 8
      ENDIF
    7 F2 = CUMF + DT * RINT
      IF (F2 .EQ. 0.0D0) GO TO 20
      R2 = -XKF * SZF * (CD + F2) / (1.0D0 - DEXP(SZF * F2))
      IF (R2 .GT. RINT) GO TO 20
      F = CUMF + R2 * DT
      DO 9 I = 1, 20
         R2 = -XKF * SZF * (CD + F) / (1.0D0 - DEXP(SZF * F))
         IF (R2 .GT. RINT) THEN
            F1 = F
            F  = (F2 + F) * 0.5D0
            IF (DABS(F - F1) .LT. E) GO TO 11
         ELSE
            F2 = F
            F  = (F1 + F) * 0.5D0
            IF (DABS(F - F2) .LT. E) GO TO 11
         ENDIF
    9 CONTINUE
   11 CONTINUE
      TP = (IT - 1) * DT + (F - CUMF) / RINT
      IF (TP .GT. IT * DT) GO TO 20
C
C  SET UP DEFINITE INTEGRAL CONSTANT USING FP
C
    8 CONST = 0.0D0
      FAC   = 1.0D0
      FC    = (F + CD)
      DO 12 J = 1, 10
         FAC   = FAC * J
         ADD   = (FC * SZF) ** J / (J * FAC)
         CONST = CONST + ADD
   12 CONTINUE
      CONST = DLOG(FC) - (DLOG(FC) + CONST) / DEXP(SZF * CD)
C
      IROF = 1
      F    = F + 0.5D0 * RINT * (IT * DT - TP)
   10 CONTINUE
C
C  NEWTON-RAPHSON SOLUTION FOR F(T)
      DO 14 I = 1, 20
         FC  = (F + CD)
         SUM = 0.0D0
         FAC = 1.0D0
         DO 13 J = 1, 10
            FAC = FAC * J
            ADD = (FC * SZF) ** J / (J * FAC)
            SUM = SUM + ADD
   13    CONTINUE
         FUNC = -(DLOG(FC) - (DLOG(FC) + SUM) / DEXP(SZF * CD) - CONST)
     1          / (XKF * SZF) - (IT * DT - TP)
         DFUNC = (DEXP(SZF * F) - 1.0D0) / (XKF * SZF * FC)
         DF    = -FUNC / DFUNC
         F     = F + DF
         IF (DABS(DF) .LE. E) GO TO 15
   14 CONTINUE
   15 CONTINUE
      IF (F .LT. CUMF + RINT) THEN
         DF   = F - CUMF
         CUMF = F
C  SET UP INITIAL ESTIMATE FOR NEXT TIME STEP
         F    = F + DF
         RETURN
      ENDIF
   20 CONTINUE
C  THERE IS NO PONDING IN THIS TIME STEP
      IROF = 0
      DF   = RINT * DT
      CUMF = CUMF + DF
      RETURN
      END
"""


def build_and_run(workdir):
    src = os.path.join(workdir, "expinf_reference.f")
    exe = os.path.join(workdir, "expinf_reference")
    with open(src, "w") as f:
        f.write(FORTRAN_SOURCE)

    fc = shutil.which("gfortran") or shutil.which("f77") or shutil.which("g77")
    if fc is None:
        sys.stderr.write(
            "No Fortran compiler found. One is needed ONLY to regenerate this\n"
            "fixture; building and testing lsdtt3 does not require one.\n")
        raise SystemExit(2)

    subprocess.run([fc, "-O2", "-o", exe, src], check=True)

    rows = []
    for name, s in STORMS.items():
        stdin = "%g %g %g %g %g %d\n" % (s["xk0"], s["hf"], s["dth"],
                                         s["szm"], s["dt"], len(s["rain"]))
        stdin += " ".join(repr(r) for r in s["rain"]) + "\n"
        out = subprocess.run([exe], input=stdin, capture_output=True,
                             text=True, check=True).stdout
        for line in out.strip().splitlines():
            it, irof, df, cumf = line.split()
            rows.append((name, s["xk0"], s["hf"], s["dth"], s["szm"], s["dt"],
                         int(it), s["rain"][int(it) - 1], int(irof),
                         df.replace("E", "e"), cumf.replace("E", "e")))
    return rows


def main():
    ap = argparse.ArgumentParser(
        description="Regenerate the EXPINF reference fixture (maintainer tool).")
    ap.add_argument("--out", required=True)
    ap.add_argument("--check", action="store_true",
                    help="compare against the existing file instead of writing")
    a = ap.parse_args()

    with tempfile.TemporaryDirectory(prefix="expinf_") as tmp:
        rows = build_and_run(tmp)

    body_lines = ["storm,xk0,hf,dth,szm,dt,it,rint,irof,df,cumf"]
    body_lines += [",".join(str(x) for x in r) for r in rows]

    if a.check:
        have = [l for l in open(a.out).read().splitlines()
                if l and not l.startswith("#")]
        if have != body_lines:
            sys.stderr.write("FIXTURE DRIFT: %s no longer matches EXPINF\n" % a.out)
            return 1
        print("fixture matches EXPINF (%d rows)" % len(rows))
        return 0

    with open(a.out, "w") as f:
        f.write("\n".join(HEADER) + "\n" + "\n".join(body_lines) + "\n")
    print("wrote %d rows to %s" % (len(rows), a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
