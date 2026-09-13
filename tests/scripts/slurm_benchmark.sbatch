#!/bin/bash
#
# Slurm submission wrapper for the lsdtt3 benchmark suite.
#
# ---------------------------------------------------------------------------
# THIS FILE IS OPTIONAL AND CONTAINS NO BENCHMARK LOGIC.
#
# It exists only to arrange for a machine with a scheduler to run exactly the
# same measurement a machine without one runs.  Nothing in the benchmark knows
# or cares that Slurm exists.  If you are on a workstation, a laptop, a cloud
# VM, or a cluster with PBS/LSF/SGE instead, ignore this file entirely and run
# the suite directly -- you will get an identical campaign directory.  Do not
# add measurement logic here, or results from Slurm sites and non-Slurm sites
# stop being comparable.
# ---------------------------------------------------------------------------
#
# NOTHING BELOW NEEDS EDITING.  Everything site-specific is an environment
# variable or an sbatch flag, so this file can stay exactly as committed.
#
#   sbatch --account=PROJECT tests/scripts/slurm_benchmark.sbatch
#
# Slurm's own settings: pass them to sbatch, which overrides the #SBATCH lines
# below, or export Slurm's native SBATCH_* variables.  There is deliberately NO
# `#SBATCH --account` line here -- a placeholder would silently become the
# default the one time you forgot to override it.
#
#   sbatch --account=PROJ --partition=compute --cpus-per-task=32 \
#          --time=08:00:00 tests/scripts/slurm_benchmark.sbatch
#   export SBATCH_ACCOUNT=PROJ SBATCH_PARTITION=compute   # same effect
#
# This script's settings: export any of these before submitting.
#
#   BENCH_SETUP_CMD    shell to eval before anything runs. This is where module
#                      loads and conda activation go, e.g.
#                        BENCH_SETUP_CMD="module load anaconda; conda activate lsdtt3"
#   BENCH_SETUP_SCRIPT a file to source instead, if the setup is long-lived.
#                      Both run if both are set: script first, then command.
#   BENCH_RUNNER       "pixi" (default) or "python". See RUNNERS below.
#   BENCH_PYTHON       interpreter for BENCH_RUNNER=python. Default "python".
#   BENCH_BUILD_DIR    directory holding the built lsdtt-* binaries.
#                      Default <repo>/build.
#   BENCH_REPO         the lsdtt3 checkout. Default: derived from this file's
#                      own path.
#   BENCH_OUT_DIR      campaign output. Default <submit dir>/benchmark_campaigns.
#   BENCH_FIXTURE_DIR  fixture cache. Default <submit dir>/benchmark_fixtures.
#   BENCH_PROFILE      "server" (default) or "laptop", for a quick shakedown.
#
# Anything after the script name is forwarded to the suite unchanged:
#
#   sbatch --account=PROJ tests/scripts/slurm_benchmark.sbatch --threads 1,2,4,8,16,32
#
# A full example: conda environment, own account, fixtures on scratch, results
# in the directory you submit from.
#
#   cd /path/to/my/data
#   export BENCH_SETUP_CMD="module load anaconda; conda activate lsdtt3"
#   export BENCH_RUNNER=python
#   export BENCH_BUILD_DIR=/path/to/lsdtt3/build
#   export BENCH_FIXTURE_DIR=/scratch/$USER/lsdtt3_fixtures
#   sbatch --account=PROJECT --cpus-per-task=32 \
#          /path/to/lsdtt3/tests/scripts/slurm_benchmark.sbatch --threads 1,2,4,8,16,32
#
#SBATCH --job-name=lsdtt3-bench
#
# Logs go in the SUBMIT DIRECTORY ITSELF, not a subdirectory of it.
#
# Slurm resolves a relative --output against the submit directory and opens the
# file BEFORE the job script runs. A path like `slurm_logs/...` therefore fails
# instantly, and silently, whenever that directory does not already exist --
# the job dies before it can create it, and the error has nowhere to be
# written. Submitting from a fresh data directory is exactly when that
# happens. Keep these paths free of any directory component.
#SBATCH --output=lsdtt3-bench-%j.out
#SBATCH --error=lsdtt3-bench-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
# --cpus-per-task MUST be at least the largest thread count in the profile.
# The `server` profile sweeps 1/2/4/8, so 8 is the floor; the guard below
# refuses to run rather than silently producing an 8-thread row that Slurm
# confined to fewer cores. Raise this and --threads together to sweep higher.
#SBATCH --cpus-per-task=8
#
# One node, exclusive. A benchmark sharing a node with another job is measuring
# the co-tenant's memory bandwidth as much as its own. This is the single most
# important line in the file for peer review; drop it and the numbers become
# unreproducible noise.
#SBATCH --exclusive
#
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${SUBMIT_DIR}"

echo "=== lsdtt3 benchmark campaign ==="
echo "job              : ${SLURM_JOB_ID:-<none>}"
echo "node             : $(hostname)"
echo "account          : ${SLURM_JOB_ACCOUNT:-<none>}"
echo "cpus-per-task    : ${SLURM_CPUS_PER_TASK:-<unset>}"
echo "submit directory : ${SUBMIT_DIR}"
echo "started          : $(date -Is)"
echo

# ---------------------------------------------------------------------------
# Site setup: modules, conda, anything else
#
# `set -u` is lifted around this. conda's own profile scripts, and many site
# module systems, reference unset variables internally and abort under -u --
# which would fail the job in a way that looks like a benchmark problem rather
# than an environment one.
# ---------------------------------------------------------------------------
if [ -n "${BENCH_SETUP_SCRIPT:-}" ]; then
    echo "sourcing BENCH_SETUP_SCRIPT: ${BENCH_SETUP_SCRIPT}"
    set +u
    # shellcheck disable=SC1090
    source "${BENCH_SETUP_SCRIPT}"
    set -u
fi

if [ -n "${BENCH_SETUP_CMD:-}" ]; then
    echo "running BENCH_SETUP_CMD: ${BENCH_SETUP_CMD}"
    set +u
    eval "${BENCH_SETUP_CMD}"
    set -u
fi

if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "conda environment: ${CONDA_DEFAULT_ENV}"
fi
echo

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# Derived from this file's own location: <repo>/tests/scripts/<this file>.
# Slurm copies the batch script to a spool directory on some configurations, so
# the derivation is verified below and BENCH_REPO overrides it.
_self="${BASH_SOURCE[0]:-$0}"
_derived_repo="$(cd "$(dirname "${_self}")/../.." 2>/dev/null && pwd || true)"
BENCH_REPO="${BENCH_REPO:-${_derived_repo}}"

SUITE="${BENCH_REPO}/tests/scripts/run_benchmark_suite.py"
if [ ! -f "${SUITE}" ]; then
    echo "ERROR: cannot find the benchmark suite at ${SUITE}" >&2
    echo "       Slurm may have copied this script out of the repository." >&2
    echo "       Export BENCH_REPO=/path/to/lsdtt3 and resubmit." >&2
    exit 1
fi

BENCH_BUILD_DIR="${BENCH_BUILD_DIR:-${BENCH_REPO}/build}"
BENCH_PROFILE="${BENCH_PROFILE:-server}"

# Both default into the SUBMIT directory rather than the repository: the repo
# usually sits on a quota'd home filesystem and the server profile stores
# ~300 MB of region sources. Submitting from a data directory therefore puts
# the campaign and the fixtures in that data directory.
#
# The fixture cache is worth putting somewhere persistent and shared -- it is
# reused by every later campaign and re-downloading it is the slowest part of a
# cold run. The campaign directory is small and is the thing you keep.
BENCH_OUT_DIR="${BENCH_OUT_DIR:-${SUBMIT_DIR}/benchmark_campaigns}"
BENCH_FIXTURE_DIR="${BENCH_FIXTURE_DIR:-${SUBMIT_DIR}/benchmark_fixtures}"
mkdir -p "${BENCH_OUT_DIR}" "${BENCH_FIXTURE_DIR}"

if [ ! -x "${BENCH_BUILD_DIR}/lsdtt-raster-preprocessing" ]; then
    echo "ERROR: no lsdtt3 binaries in ${BENCH_BUILD_DIR}" >&2
    echo "       (looked for lsdtt-raster-preprocessing)" >&2
    echo "       Export BENCH_BUILD_DIR=/path/to/your/build and resubmit." >&2
    exit 1
fi

# Submitting from inside the checkout dumps the campaign AND the fixture cache
# into the source tree -- ~300 MB of DEMs under, say, tests/scripts/. It is
# allowed, because someone may genuinely want it, but it is almost never what
# was intended and it is invisible until the repository grows.
_repo_real="$(cd "${BENCH_REPO}" && pwd -P)"
_submit_real="$(cd "${SUBMIT_DIR}" && pwd -P)"
case "${_submit_real}/" in
    "${_repo_real}"/*)
        echo "WARNING: submitted from inside the lsdtt3 checkout." >&2
        echo "         Campaign output and ~300 MB of cached DEMs will be" >&2
        echo "         written into the source tree:" >&2
        echo "           ${BENCH_OUT_DIR}" >&2
        echo "           ${BENCH_FIXTURE_DIR}" >&2
        echo "         Submit from a data directory instead, or set" >&2
        echo "         BENCH_OUT_DIR and BENCH_FIXTURE_DIR explicitly." >&2
        echo >&2
        ;;
esac

echo "repository      : ${BENCH_REPO}"
echo "binaries        : ${BENCH_BUILD_DIR}"
echo "profile         : ${BENCH_PROFILE}"
echo "campaign output : ${BENCH_OUT_DIR}"
echo "fixture cache   : ${BENCH_FIXTURE_DIR}"
echo "scratch (staged): ${TMPDIR:-/tmp}"
echo

# The suite pins OMP_NUM_THREADS per measurement, so anything inherited from
# the site profile or a conda activation would be overwritten anyway. Unsetting
# it makes that explicit rather than leaving a stale value for a reader to
# wonder about.
unset OMP_NUM_THREADS || true

# Refuse to run under-provisioned rather than produce mislabelled rows. An
# "8 thread" measurement on a 4-core allocation is not a slow result, it is a
# wrong one, and it is indistinguishable from a genuine scaling ceiling.
MAX_THREADS=8
if [ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ "${SLURM_CPUS_PER_TASK}" -lt "${MAX_THREADS}" ]; then
    echo "ERROR: --cpus-per-task=${SLURM_CPUS_PER_TASK} is below the ${MAX_THREADS}" >&2
    echo "       threads the server profile sweeps. Resubmit with" >&2
    echo "       --cpus-per-task=${MAX_THREADS} or higher, or pass a smaller" >&2
    echo "       --threads list to this script." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# RUNNERS
#
# pixi   the reference path. The `benchmark` environment holds RichDEM,
#        TopoToolbox and Whitebox at the pinned versions the published numbers
#        were measured against, so this is what makes a campaign comparable
#        with anyone else's.
#
# python any interpreter you have already put on PATH -- typically a conda
#        environment activated in BENCH_SETUP_CMD. Use this where pixi is not
#        available or the site requires conda.
#
#        *** The reference implementations must be importable in that
#        environment or you will get an lsdtt3-only campaign. *** The suite
#        skips a reference it cannot import (exit 77) rather than failing, so a
#        missing package is silent in the run and shows up only as absent rows.
#        The check below reports what is actually importable before the
#        campaign starts, so that is visible in the log up front.
# ---------------------------------------------------------------------------
# Report which reference implementations are actually importable, BEFORE the
# campaign starts. The suite skips a reference it cannot import (exit 77)
# rather than failing, so a missing package is otherwise silent in the run and
# shows up only as absent rows in a table nobody re-counts. Both runners call
# this for the same reason.
check_references() {
    echo "reference implementations importable in this environment:"
    "$@" - <<'PYEOF' || true
for name in ("richdem", "topotoolbox", "whitebox_workflows",
             "osgeo.gdal", "numpy", "scipy"):
    try:
        __import__(name)
        print(f"    {name:20s} yes")
    except Exception as exc:
        print(f"    {name:20s} NO  ({type(exc).__name__})")

# PROJ database check. An interpreter used WITHOUT properly activating its
# environment has GDAL but no PROJ_DATA, and GDAL then silently falls back to
# CRS parameters read from the GeoTIFF keys instead of the EPSG registry. The
# fixture crops carry a CRS and the fetch reprojects, so a broken PROJ is a
# correctness problem, not merely noisy warnings. It is easy to hit: it is
# exactly what happens if you point BENCH_PYTHON at an environment's python
# binary directly instead of activating the environment first.
try:
    from osgeo import osr
    osr.UseExceptions()
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32645)
    print("    PROJ database        yes")
except Exception as exc:
    print(f"    PROJ database        NO  ({type(exc).__name__})")
    print("    *** GDAL cannot read the EPSG registry. ACTIVATE the")
    print("        environment (conda activate <env>) rather than calling its")
    print("        python binary directly, or export PROJ_DATA/PROJ_LIB. ***")

print("    (a missing reference is SKIPPED, not an error -- it produces no")
print("     rows, so check this list before trusting the tables)")
PYEOF
    echo
}


BENCH_RUNNER="${BENCH_RUNNER:-pixi}"
BENCH_PYTHON="${BENCH_PYTHON:-python}"

case "${BENCH_RUNNER}" in
    pixi)
        if ! command -v pixi >/dev/null 2>&1; then
            echo "ERROR: BENCH_RUNNER=pixi but pixi is not on PATH on this node." >&2
            echo "       Compute nodes often have a different PATH than the" >&2
            echo "       login node. Either add it:" >&2
            echo "         export PATH=\"\$HOME/.pixi/bin:\$PATH\"" >&2
            echo "       or switch to a conda environment:" >&2
            echo "         export BENCH_RUNNER=python" >&2
            echo "         export BENCH_SETUP_CMD=\"conda activate <env>\"" >&2
            exit 1
        fi
        # Install before timing starts: a first-run solve during the campaign
        # would land inside a measurement.
        pixi install --manifest-path "${BENCH_REPO}/pixi.toml" -e benchmark
        # `-e benchmark` is NOT optional. Without it pixi runs in the DEFAULT
        # environment, which has no RichDEM, TopoToolbox or Whitebox -- and
        # because the suite treats a missing reference as a clean skip rather
        # than an error, the campaign completes "successfully" with only
        # lsdtt3 rows in it. That happened during development of this file and
        # is invisible unless you count the rows.
        RUN_CMD=(pixi run --manifest-path "${BENCH_REPO}/pixi.toml"
                 -e benchmark python "${SUITE}")
        check_references pixi run --manifest-path "${BENCH_REPO}/pixi.toml" \
                                  -e benchmark python
        ;;
    python)
        if ! command -v "${BENCH_PYTHON}" >/dev/null 2>&1; then
            echo "ERROR: BENCH_PYTHON=${BENCH_PYTHON} is not on PATH." >&2
            echo "       Did BENCH_SETUP_CMD activate the environment?" >&2
            exit 1
        fi
        echo "interpreter     : $(command -v "${BENCH_PYTHON}")"
        check_references "${BENCH_PYTHON}"
        RUN_CMD=("${BENCH_PYTHON}" "${SUITE}")
        ;;
    *)
        echo "ERROR: BENCH_RUNNER must be 'pixi' or 'python', got '${BENCH_RUNNER}'" >&2
        exit 1
        ;;
esac
# Absolute paths throughout, because pixi runs its tasks with the working
# directory set to the manifest root -- a relative --out-dir would resolve
# inside the repository no matter where this job was submitted from.
#
# "$@" forwards any extra sbatch arguments to the suite unchanged.
"${RUN_CMD[@]}" \
    --profile "${BENCH_PROFILE}" \
    --build-dir "${BENCH_BUILD_DIR}" \
    --out-dir "${BENCH_OUT_DIR}" \
    --fixture-dir "${BENCH_FIXTURE_DIR}" \
    "$@"

echo
echo "finished : $(date -Is)"
echo
echo "The campaign directory is under ${BENCH_OUT_DIR} and is"
echo "self-contained: tar it and send it as-is."
