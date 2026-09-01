# Running benchmarks on your own system

lsdtt3 ships a benchmark suite that measures how it compares against three
other topographic-analysis implementations — **RichDEM**, **TopoToolbox** and
**Whitebox Workflows** — on *your* hardware, on real terrain.

It is deliberately **not** part of the test suite. `ctest` never runs it, it
asserts nothing, and it cannot fail a build. Its only product is a
self-contained directory of measurements you can read, keep, or send to
somebody else.

If you want the reasoning behind the design — why particular comparisons are
valid and others are not — see
[Testing and benchmarking](testing_and_benchmarking.md). This page is just how
to run it.

---

## Quick start

```sh
pixi run bench-laptop
```

That is the whole thing. It takes a few minutes, downloads ~8 MB of DEM the
first time, and writes a directory containing a `SUMMARY.md` you can open
straight away.

To see what it *would* do without measuring anything:

```sh
pixi run bench-laptop --dry-run
```

---

## What you need

**The Python environment.** The three reference implementations live in pixi's
`benchmark` environment:

```sh
pixi install -e benchmark
```

This pulls RichDEM and TopoToolbox and Whitebox Workflows at pinned versions.
It is a large environment (~1.3 GB) because it includes GDAL. All three
references are open-source and none needs a licence key.

!!! note "No Whitebox wheel for Intel Macs"
    Whitebox Workflows ships no `osx-64` wheel. On an Intel Mac use the
    `benchmark-lite` environment instead, which has RichDEM and TopoToolbox
    only:
    ```sh
    pixi run -e benchmark-lite python tests/scripts/run_benchmark_suite.py --profile laptop
    ```
    You will simply get fewer rows.

**The lsdtt3 binaries**, built in **Release** mode. A Debug build measures the
compiler's inlining decisions rather than the algorithms, and the summary will
warn you if it detects one:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

**Network access, once.** The first run downloads a DEM per region. After that
everything is cached and the suite runs offline.

---

## The two profiles

| | `laptop` | `server` |
|---|---|---|
| DEM sizes | 500², 1000² | 1000², 2000², 4000², 8000² |
| Threads | 1, 2, 4 | 1, 2, 4, 8 |
| Repeats | 3 | 5 |
| Download (first run) | ~4 MB per region | ~256 MB per region |
| Runtime | minutes | hours |

```sh
pixi run bench-laptop      # a real baseline, and the shakedown run
pixi run bench-server      # the numbers worth quoting
```

**Run `bench-laptop` first, always.** It exercises every code path in a few
minutes, so a missing binary or an unreachable download shows up immediately
rather than after a queue wait. It is genuinely useful in its own right, but
its DEMs are small enough that fixed overheads matter — see
[Reading the results](#reading-the-results).

---

## The terrain

Two regions, chosen as end-members because depression handling is the thing
that separates these tools most sharply:

| region | location | relief | depression cells |
|---|---|---|---|
| `alpine` | High Himalaya, Nepal (86.9 E, 27.9 N) | 1662 m | 6.4 % |
| `lakeland` | Canadian Shield, NW Ontario (92.09 W, 50.57 N) | 94 m | 18.5 % |

The suite downloads **one large DEM per region** and crops every size in the
ladder out of it, so terrain stays as consistent as possible as the grid grows
and you only ever download once.

This matters more than it sounds. In one measured campaign, Whitebox's fill was
**5.5× faster than lsdtt3 on high-relief terrain and 44× slower on low-relief
terrain** — same code, same grid size, only the landscape changed. A benchmark
on one terrain type would have told you something confidently wrong. Always
read both regions.

---

## What you get

One directory per run, `<profile>_<host>_<timestamp>/`, containing:

| file | what it is |
|---|---|
| `SUMMARY.md` | rendered tables — start here |
| `benchmark_raw.csv` | one row per measurement, every column |
| `run_manifest.json` | profile, fixture hashes, terrain stats, exact command line |
| `hardware_manifest.json` | CPU, cores, RAM, compiler, git SHA, reference versions |

The directory is self-contained: tar it and send it as-is. Everything needed to
interpret or reproduce the numbers is inside it.

### Where it gets written

```sh
pixi run bench-laptop                                    # into the repo
python tests/scripts/run_benchmark_suite.py --out-dir .  # into the current directory
```

A relative `--out-dir` resolves against the **current directory**.

!!! warning "pixi runs tasks from the repository root"
    `pixi run bench-server --out-dir results` puts results in the *repo*, not
    where you invoked it, because pixi sets the working directory to the
    manifest root. Pass an **absolute** path to collect results elsewhere.

Three separate locations, with three different lifetimes:

| what | default | override |
|---|---|---|
| campaign output | `./benchmark_campaigns/` | `--out-dir` |
| DEM cache | `<repo>/tests/data/dem_cache/` | `--fixture-dir` |
| working scratch | `$TMPDIR` | `--scratch-dir` |

On a shared machine or a cluster, move the DEM cache off the repository with
`--fixture-dir` — the `server` profile stores ~512 MB there, and it is reused by
every later run, so put it somewhere persistent.

---

## Reading the results

`SUMMARY.md` is written to be read top to bottom. Four things decide whether a
number means anything.

### 1. Check the header

```
- build type: `Release`, compiler `c++ (Ubuntu 11.4.0) 11.4.0`
- measured on: `ext2/ext3` at `/tmp/lsdtt3_benchmark_fixtures` (staged from v9fs)
```

A non-Release build, or a filesystem flagged as slow, and the summary says so
in bold. Both invalidate the run.

### 2. Only compare like with like

Every row records a `timing_scope`, and rows are only ever tabled against
others of the same scope:

- `cli_reported_compute` — the compute time an lsdtt3 CLI printed itself
- `in_process_routine` — a reference script timing the routine in-process
- `whole_process` — wall clock for an entire CLI run, **never** put in a ratio

This is not pedantry. lsdtt3's fill on a small DEM is ~7.6 ms of compute inside
a ~279 ms process — whole-process wall clock was 97 % startup and file I/O.

### 3. Respect the flags

- **`overhead_dominated`** — the measurement is within 3× of that engine's own
  measured fixed cost, so it is mostly overhead. Excluded from the ratio table
  and listed separately. Expect many of these on the `laptop` profile.
- **`report_only`** — the engines are not doing equivalent work, so a ratio
  would be a claim about the definitions rather than the implementations.
  Absolute times only. All surface-metric and multi-direction-flow rows are
  report-only.

RichDEM in particular carries ~10 ms of fixed cost per call, from a per-call
array construction plus the citation block it prints on every invocation. On a
small DEM that is most of what you measure. The summary states the figure.

### 4. Scaling curves need a ceiling, and the big DEMs

The `laptop` profile's speedups sit between about 0.8× and 1.4×. There is not
enough work in a 500² grid to amortise thread setup, so use the `server`
profile for anything about parallel scaling.

More importantly, **a scaling number means nothing without knowing what the
machine can do**. lsdtt3's surface-metric kernels are memory-bandwidth-bound,
not compute-bound, so they scale like a streaming kernel rather than like an
arithmetic one. On the development machine a compute-bound loop reached 3.09×
on 4 threads while a streaming loop managed only 1.14× — and 0.27× on 8, where
hyperthreads contend for the same load/store ports.

Two practical consequences:

- **Cap `--threads` at your physical core count.** Including hyperthreads makes
  bandwidth-bound work several times slower, which reads as a scaling collapse
  but is just core sharing.
- **Flow routing and depression filling are not swept at all.** Filling,
  breaching and flow-graph construction are ordered traversals with no
  parallelism to measure, so they run once, at one thread, and are excluded
  from the scaling table. `SUMMARY.md` lists which cases those are. Only the
  surface metrics are genuinely parallel code, so only they are swept.

The measured breakdown — which pass fails to scale and why, and a NUMA
first-touch hypothesis for multi-socket machines — is in
[Testing and benchmarking](testing_and_benchmarking.md).

---

## Running on a server or cluster

Everything is a flag or an environment variable; no file needs editing.

```sh
pixi run bench-server --out-dir "$PWD" --fixture-dir /scratch/$USER/lsdtt3_fixtures
```

### Under Slurm

`tests/scripts/slurm_benchmark.sbatch` contains **no measurement logic**. It
just arranges for a scheduler node to run the same command everyone else runs,
so results from Slurm and non-Slurm sites stay comparable. If your cluster uses
PBS, LSF or SGE, ignore the file and run the command above in your own job
script — you get an identical campaign directory.

```sh
cd /path/to/my/data
sbatch --account=PROJECT tests/scripts/slurm_benchmark.sbatch
```

Slurm settings go to `sbatch` (which overrides the `#SBATCH` lines) or through
Slurm's own `SBATCH_*` variables. There is deliberately **no `#SBATCH
--account`** line, so you cannot accidentally submit under a stale default.

The script's own settings are environment variables:

| variable | purpose |
|---|---|
| `BENCH_SETUP_CMD` | shell run first — where `module load` / `conda activate` go |
| `BENCH_SETUP_SCRIPT` | a file to source instead |
| `BENCH_RUNNER` | `pixi` (default) or `python` |
| `BENCH_PYTHON` | interpreter, when `BENCH_RUNNER=python` |
| `BENCH_BUILD_DIR` | where the `lsdtt-*` binaries are |
| `BENCH_REPO` | the lsdtt3 checkout |
| `BENCH_OUT_DIR` | campaign output (default: the submit directory) |
| `BENCH_FIXTURE_DIR` | DEM cache (default: the submit directory) |
| `BENCH_PROFILE` | `server` (default) or `laptop` |

Arguments after the script name reach the suite unchanged. A conda-based site,
with results where you submit from and DEMs on scratch:

```sh
cd /path/to/my/data
export BENCH_SETUP_CMD="module load anaconda; conda activate lsdtt3"
export BENCH_RUNNER=python
export BENCH_BUILD_DIR=/path/to/lsdtt3/build
export BENCH_FIXTURE_DIR=/scratch/$USER/lsdtt3_fixtures
sbatch --account=PROJECT --cpus-per-task=32 \
       /path/to/lsdtt3/tests/scripts/slurm_benchmark.sbatch --threads 1,2,4,8,16,32
```

The job requests an exclusive node. A benchmark sharing a node measures the
co-tenant's memory bandwidth as much as its own, so do not remove that.

If you raise `--threads`, raise `--cpus-per-task` to match. The script refuses
to start otherwise: an "8 thread" row on a 4-core allocation is not a slow
result, it is a wrong one, and it looks exactly like a real scaling ceiling.

---

## Customising a run

```sh
python tests/scripts/run_benchmark_suite.py --profile server \
    --regions alpine --ladder 1000,2000 --threads 1,4 --cases fill,breach
```

| flag | effect |
|---|---|
| `--regions` | which terrain (`alpine`, `lakeland`) |
| `--ladder` | crop sizes in pixels; the source is fetched at the largest |
| `--threads` | thread counts to sweep |
| `--repeats` | timed repeats per measurement |
| `--suites` | `flow`, `surface`, `mfd`, `depression` |
| `--cases` | individual cases — see `--list` |
| `--timeout-seconds` | per-measurement deadline (default 1800) |
| `--no-fetch` | never download; skip anything missing |
| `--dry-run` | print the plan and exit |
| `--list` | print every case and profile as JSON |

---

## Troubleshooting

**"No row was measured on a comparable basis"** — normal on the `laptop`
profile. The DEMs are small enough that most measurements are overhead-bound.
Run the `server` profile.

**Only lsdtt3 rows appear.** The suite *skips* any reference it cannot import
rather than failing, so a campaign can complete "successfully" with a third of
the data. Check you are in the `benchmark` environment:

```sh
pixi run -e benchmark python -c "import richdem, topotoolbox, whitebox_workflows; print('all three ok')"
```

**`timeout` rows.** A measurement exceeded `--timeout-seconds` and was killed;
the campaign continued. That is a real result — the engine was working, just
too slow — not a failure. Raise the limit if you want the number.

**Warnings about PROJ or the EPSG registry.** Your GDAL cannot find its PROJ
database. Activate the environment (`conda activate <env>`, or use
`pixi run`) rather than calling an interpreter's binary by its full path.

**"--threads goes up to N but this machine reports M logical cores."** Rows
above the core count are oversubscribed and their scaling figures are
meaningless. Pass a smaller `--threads`.

**A download fails.** The suite reports it and carries on with whatever it has.
Re-run when the network is back; anything already cached is not fetched again.

---

## Sharing results

Send the whole campaign directory. It carries the hardware, the compiler, the
lsdtt3 commit, the reference versions, the DEM hashes and the exact command
line, which is everything needed to interpret or repeat it.

Two things to state alongside any number you quote:

- **which region** — the alpine and lakeland results genuinely disagree, and a
  figure without its terrain is not interpretable;
- **that value correctness is not tested here.** This suite measures speed and
  memory. Correctness lives in the CTest suite; see
  [Testing and benchmarking](testing_and_benchmarking.md).
