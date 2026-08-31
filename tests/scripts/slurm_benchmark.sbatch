#!/bin/bash
#
# Slurm submission wrapper for the lsdtt3 benchmark suite.
#
# ---------------------------------------------------------------------------
# THIS FILE IS OPTIONAL AND CONTAINS NO BENCHMARK LOGIC.
#
# It exists only to arrange for a machine with a scheduler to run exactly the
# same command a machine without one runs:
#
#     pixi run bench-server
#
# Nothing in the benchmark knows or cares that Slurm exists.  If you are on a
# workstation, a laptop, a cloud VM, or a cluster with PBS/LSF/SGE instead,
# ignore this file entirely and run that command directly -- you will get an
# identical campaign directory.  Do not add measurement logic here, or results
# from Slurm sites and non-Slurm sites stop being comparable.
# ---------------------------------------------------------------------------
#
# Submit with:
#
#     sbatch tests/scripts/slurm_benchmark.sbatch
#
# Override any of the resource requests on the command line, which takes
# precedence over the #SBATCH lines below:
#
#     sbatch --cpus-per-task=32 --time=08:00:00 tests/scripts/slurm_benchmark.sbatch
#
# Pass benchmark arguments after the script name; they reach the suite
# unchanged:
#
#     sbatch tests/scripts/slurm_benchmark.sbatch --threads 1,2,4,8,16,32
#
#SBATCH --job-name=lsdtt3-bench
#SBATCH --output=slurm_logs/lsdtt3-bench-%j.out
#SBATCH --error=slurm_logs/lsdtt3-bench-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
# --cpus-per-task MUST be at least the largest thread count in the profile.
# The `server` profile sweeps 1/2/4/8, so 8 is the floor; the guard below
# refuses to run rather than silently producing an 8-thread row that Slurm
# confined to fewer cores. Raise both this and --threads together to sweep
# higher.
#SBATCH --cpus-per-task=8
#
# One node, exclusive. A benchmark sharing a node with another job is
# measuring that job's memory bandwidth as much as its own. This is the single
# most important line in the file for peer review; drop it and the numbers
# become unreproducible noise.
#SBATCH --exclusive
#
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -euo pipefail

mkdir -p slurm_logs

echo "=== lsdtt3 benchmark campaign ==="
echo "job              : ${SLURM_JOB_ID:-<none>}"
echo "node             : $(hostname)"
echo "cpus-per-task    : ${SLURM_CPUS_PER_TASK:-<unset>}"
echo "submit directory : ${SLURM_SUBMIT_DIR:-$PWD}"
echo "started          : $(date -Is)"
echo

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# The suite pins OMP_NUM_THREADS per measurement, so anything inherited from
# the site profile would be overwritten anyway. Unsetting it here makes that
# explicit rather than leaving a stale value in the environment for a reader to
# wonder about.
unset OMP_NUM_THREADS || true

# Refuse to run under-provisioned rather than produce mislabelled rows. An
# "8 thread" measurement on a 4-core allocation is not a slow result, it is a
# wrong one, and it would be indistinguishable from a genuine scaling ceiling.
MAX_THREADS=8
if [ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ "${SLURM_CPUS_PER_TASK}" -lt "${MAX_THREADS}" ]; then
    echo "ERROR: --cpus-per-task=${SLURM_CPUS_PER_TASK} is below the ${MAX_THREADS}" >&2
    echo "       threads the server profile sweeps. Resubmit with" >&2
    echo "       --cpus-per-task=${MAX_THREADS} or higher, or pass a smaller" >&2
    echo "       --threads list to this script." >&2
    exit 1
fi

if ! command -v pixi >/dev/null 2>&1; then
    echo "ERROR: pixi is not on PATH on this node." >&2
    echo "       Cluster nodes often have a different PATH than the login" >&2
    echo "       node. Load your module or add pixi's install directory," >&2
    echo "       e.g.  export PATH=\"\$HOME/.pixi/bin:\$PATH\"" >&2
    exit 1
fi

# Install into the compute node's view of the filesystem before timing starts.
# A first-run solve during the campaign would land inside a measurement.
pixi install -e benchmark

# "$@" forwards any extra sbatch arguments to the suite unchanged.
pixi run bench-server "$@"

echo
echo "finished : $(date -Is)"
echo
echo "The campaign directory is under benchmark_campaigns/ and is"
echo "self-contained: tar it and send it as-is."
