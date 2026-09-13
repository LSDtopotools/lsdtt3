# Running benchmarks on your own system

lsdtt3 includes a benchmark suite that measures how it compares against three
other topographic analysis implementations — **RichDEM**, **TopoToolbox** and
**Whitebox Workflows** on *your* hardware, on real terrain.

It is deliberately not part of the test suite: the test suite is used to making sure the clis and routines do not break with updates and have the correct implementation. This is a different excercise. `ctest` never runs it  and it cannot fail a build. Its only product is a
self-contained directory of measurements you can use to see how fast and how memory thirsty the difference implemetaions are. 

If you want the reasoning behind the design — why particular comparisons are
valid and others are not — see
[Testing and benchmarking](../developer_notes/testing_and_benchmarking.md). This page is just how
to run it.

It should be noted here that the comparisons are not all like-for-like. Details about different implemetations and choices made are available in the report produced by these benchmarking tests. 

---

## Quick start

You need to install a bunch of stuff to get Whitebox Workflows, RichDEM and TopoToolbox3 working, we have wrapped all that up in a pixi environment. Once it is installed (see below) you run: 

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

!!! danger "`pixi install` creates the environment. `pixi run` is what uses it."
    **Always invoke the suite through `pixi run -e benchmark python …`, never
    through a bare `python`.** `pixi install -e benchmark` only builds the
    environment; it does not put you in it, and it does not change what `python`
    means in your shell. If you have a conda environment active — an `(lsdtt3)`
    prompt, say — a bare `python` is *that* interpreter, which has none of the
    reference packages.

    The failure is quiet and expensive. Every reference row comes back
    `skipped`, the lsdtt3 rows all succeed, and the campaign **exits
    successfully** with a third of the data:

    ```
    [  19/ 216] accum_d8  topotoolbox  alpine_500  1t  skipped  SKIP: topotoolbox not importable...
    [  20/ 216] accum_d8  whitebox     alpine_500  1t  skipped  SKIP: whitebox_workflows not importable...
    ```

    A reference that is installed but *broken* is worse still — it reports
    `error exit 1` rather than `skipped`. The usual cause is a pip-installed
    RichDEM: version 0.3.4 pins `numpy<2` and breaks against a newer NumPy.
    Use the conda-forge build that pixi provides.

    Check before committing to a long run:

    ```sh
    pixi run -e benchmark python -c "import richdem, topotoolbox, whitebox_workflows; print('all three ok')"
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

Only run the `server` version is you have powerful compute available. 
Note that if you are running on a personal computer (e.g., your laptop), the compute will become thermally limited (this is when you fan starts going crazy). In that event, your computer will throttle the CPU speed and the speed reporting will not be reliable. 

| | `laptop` | `server` |
|---|---|---|
| DEM sizes | 500², 1000² | 1000², 2000², 4000² |
| Threads | 1, 2, 4 | 1, 2, 4, 8 (Whitebox: 1, 2) |
| Repeats | 3 | 5 |
| Download (first run) | ~4 MB per region | ~64 MB per region |
| Runtime | minutes | hours |

The `server` ladder stopped at 8000² on 2026-09-13. On the first full server
campaign that rung took 4.2 h of 5.5 h, and for the surface metrics, D8 and
D-infinity every engine's throughput at 8000² was within 0.95–1.04× of its
4000² figure. Add it back with `--ladder 1000,2000,4000,8000` if you want a
size study. If the 8000 px source is already in your fixture cache it is reused
rather than downloading a new one, and existing crops keep their hashes.

Whitebox is measured at 1 and 2 threads only: on that campaign it scaled
0.96–1.04× at every thread count on every surface case, so higher counts
repeated the same number.

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
 
The different implementations do better or worse in different terrain types. 
In one measured campaign, Whitebox's fill was
**5.5× faster than lsdtt3 on high-relief terrain and 44× slower on low-relief
terrain**; this with the same code, same grid size, and only a changed landscape. 
We have given some end member landscapes to highlight this behaviour in the benchmarking.  

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

### Which profile you get

**`--profile` defaults to `laptop`.** So a bare invocation is a laptop run, and
the `server` profile always has to be asked for:

```sh
pixi run bench-laptop                                          # laptop, explicitly
pixi run bench-server                                          # server, explicitly

# the same two, spelled out
pixi run -e benchmark python tests/scripts/run_benchmark_suite.py
pixi run -e benchmark python tests/scripts/run_benchmark_suite.py --profile server
```

The two pixi tasks are just those calls with `--profile` filled in. If you are
not sure which one a directory came from, `run_manifest.json` records the
profile and the exact command line.

### Where it gets written

```sh
pixi run bench-laptop                                       # into the repo
pixi run -e benchmark python tests/scripts/run_benchmark_suite.py \
    --profile laptop --out-dir .                            # into the current directory
```

A relative `--out-dir` resolves against the **current directory**.

!!! note "`--out-dir` and `--build-dir` resolve against *different* places"
    This catches people, so it is worth stating plainly:

    | flag | a relative value resolves against |
    |---|---|
    | `--out-dir` | the **current directory** |
    | `--fixture-dir` | the **current directory** |
    | `--build-dir` | the **repository root**, derived from the script's own path |

    So the default `--build-dir build` means `<repo>/build` **no matter where
    you run from**. The binaries do *not* need to be on `$PATH`, and you do not
    need to be in the repository — running the script by its full path from
    a data directory finds them. Point `--build-dir` somewhere else only if you
    built out of tree.

    If a binary is missing the suite says so before measuring anything:

    ```
    WARNING: not found in /path/to/lsdtt3/build: lsdtt-chi-analysis
             lsdtt3 rows for those cases will be 'unavailable'.
             Build them, or pass --build-dir.
    ```

    Read that line. A campaign whose lsdtt3 rows are all `unavailable` still
    exits successfully.

!!! warning "pixi runs tasks from the repository root"
    `pixi run bench-server --out-dir results` puts results in the *repo*, not
    where you invoked it, because pixi sets the working directory to the
    manifest root. Pass an **absolute** path to collect results elsewhere.

Three separate locations, with three different lifetimes:

| what | default (direct run) | default (via sbatch) | override |
|---|---|---|---|
| campaign output | `./benchmark_campaigns/` | `<submit dir>/benchmark_campaigns/` | `--out-dir` / `BENCH_OUT_DIR` |
| DEM cache | **`<repo>/tests/data/dem_cache/`** | `<submit dir>/benchmark_fixtures/` | `--fixture-dir` / `BENCH_FIXTURE_DIR` |
| working scratch | `$TMPDIR` | `$TMPDIR` | `--scratch-dir` |

!!! warning "The DEM cache and the campaign do not default to the same place"
    On a **direct run** the campaign follows your current directory but the DEM
    cache goes **inside the repository**. So after
    `pixi run bench-laptop --out-dir /my/data` the DEMs are in
    `<repo>/tests/data/dem_cache/`, not in `/my/data/`. Only the sbatch wrapper
    puts the two side by side.

    The suite prints both paths in its first three lines. If you want them
    together on a direct run, say so:
    ```sh
    pixi run bench-laptop --out-dir /my/data --fixture-dir /my/data/fixtures
    ```

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

We do this for a reason. lsdtt3's fill on a small DEM is ~7.6 ms of compute inside
a ~279 ms process: spinup costs time but we wish to measure compute. 

### 3. Respect the flags

- **`overhead_dominated`**: the measurement is within 3× of that engine's own
  measured fixed cost, so it is mostly overhead. Excluded from the ratio table
  and listed separately. Expect many of these on the `laptop` profile.
- **`report_only`**: the engines are not doing equivalent work, so a ratio
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

!!! note "Only sweep what is actually swept"

    The `surface` suite is the only one swept over threads, and it can be run on
    its own — see [Running only the thread sweep](#running-only-the-thread-sweep).
    Any scaling figure quoted from a campaign generated before 2026-09-09 is
    invalid; that section explains why and what to re-run.

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
[Testing and benchmarking](../developer_notes/testing_and_benchmarking.md).

---

## Running on a server or cluster

Everything is a flag or an environment variable; no file needs editing.

```sh
pixi run bench-server --out-dir "$PWD" --fixture-dir /scratch/$USER/lsdtt3_fixtures
```

### Set up once (both routes)

**1. Clone and install the environment:**

```sh
git clone https://github.com/LSDtopotools/lsdtt3.git
cd lsdtt3
pixi install -e benchmark
```

**2. Build in Release mode:**

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Binaries land in `/path/to/lsdtt3/build/`, which is where the benchmark looks
by default. A Debug build measures the compiler's inlining decisions rather
than the algorithms, and the summary says so in bold if it detects one.

There are two routes, and they run the *same* measurement:

| | when |
|---|---|
| [**A plain bash session**](#route-a-a-plain-bash-session) | you have a shell on the machine and no scheduler, or you want a quick run now |
| [**Slurm**](#route-b-slurm) | the machine has a scheduler, and you want an **exclusive** node — required for anything quotable |

The Slurm wrapper contains no measurement logic. It exists only so a scheduler
site runs the identical command a non-scheduler site runs.

### How much memory the run needs

Not previously stated here, and it decides which ladder rungs you can use.

lsdtt3's flow object is about **43 bytes per cell**, and the suite runs one
engine at a time, so the peak is roughly one engine plus the staged DEM:

| rung | cells | lsdtt3 flow object | comfortable in |
|---|---:|---:|---|
| 1000² | 1 M | ~43 MB | anything |
| 2000² | 4 M | ~172 MB | anything |
| 4000² | 16 M | ~688 MB | **8 GB** |
| 8000² | 64 M | ~2.8 GB | **16 GB**, with room to watch |

Add the DEM itself (256 MB at 8000² in float32), the filled copy, and the
outputs. The Python references are the ones to watch rather than lsdtt3: they
work in NumPy and a 64 M-cell float64 array is 512 MB each, and RichDEM
constructs arrays per call.

**On 8 GB the `server` ladder should fit, but confirm rather than assume.**
Every row records `peak_rss_mb`, so after a first run:

```sh
python - <<'EOF'
import csv
rows = list(csv.DictReader(open('benchmark_campaigns/<run>/benchmark_raw.csv')))
peak = max(float(r['peak_rss_mb']) for r in rows if r.get('peak_rss_mb'))
print(f"peak observed: {peak:.0f} MB")
EOF
```

If a rung is too big the engine is killed, the suite records the failure and
carries on — so you get a campaign with holes rather than a crash. Check for
`error` values in the CSV before trusting a table. To stay clear of it, drop the
top rung:

```sh
--ladder 1000,2000
```

!!! warning "Memory columns from before 2026-09-13 are invalid"

    On Linux a child process inherits its parent's **peak** memory, and the
    kernel reports that inherited figure as the child's own. The suite used to
    launch every engine straight from its driver, which had already reached a
    large peak computing fixture statistics, so almost every row reported the
    driver's peak instead of the engine's. On the first server campaign 441 of
    504 rows read exactly 3557 MiB. Only values above that floor were real.

    Every command now starts through `tests/scripts/rusage_launcher.py`, a
    small separate process. Each campaign measures the remaining floor by
    launching `/bin/true` the same way, and prints it above the memory table.
    If that floor exceeds 64 MiB, the summary says the memory column must not
    be quoted. Cells marked `<` are close enough to the floor that they may be
    reporting it. The CSV also gains `cpu_seconds_process`: CPU time well above
    wall time means a command used more cores than its thread label says.

### Route A: a plain bash session

The honest caveat first, because it is the whole reason Slurm is preferred:

!!! warning "A shared machine measures its co-tenants"
    Without an exclusive allocation you are measuring whatever else is running,
    especially for the surface metrics, which are **memory-bandwidth-bound** and
    degrade under a co-tenant long before CPU contention shows up. A campaign
    on a busy machine produced a 0.64× "scaling" figure that took a separate
    investigation to attribute, and it was probably contention all along.

    Check the machine is quiet first (`uptime`, `top`), and record in your notes
    that it was not exclusive. Value comparisons — chi, fill, D8 agreement — are
    deterministic and unaffected; only timings are at risk.

With the [one-off setup](#set-up-once-both-routes) done:

```sh
cd /path/to/my/data
nohup pixi run --manifest-path /path/to/lsdtt3/pixi.toml -e benchmark \
    python /path/to/lsdtt3/tests/scripts/run_benchmark_suite.py \
        --profile server \
        --build-dir /path/to/lsdtt3/build \
        --out-dir "$PWD/benchmark_campaigns" \
        --fixture-dir /scratch/$USER/lsdtt3_fixtures \
        --threads 1,2,4,8 \
    > bench.log 2>&1 &
```

`nohup ... &` so it survives the session dropping — the `server` profile runs
for hours. `tmux` or `screen` do the same job and let you watch it; either is
fine.

Four things differ from the Slurm route, and all four are on you rather than the
scheduler:

- **Absolute paths.** `pixi run` sets the working directory to the manifest
  root, so a relative `--out-dir` lands in the repository rather than where you
  invoked it.
- **Cap `--threads` at the PHYSICAL core count.** `nproc` counts
  hyperthreads. Including them makes bandwidth-bound work several times slower
  and reads as a scaling collapse. `lscpu | grep -E '^Core|^Socket'` gives the
  real number.
- **Nothing enforces the thread/allocation match.** The sbatch wrapper refuses
  to start if `--cpus-per-task` is below the sweep; here there is no such
  guard, so an oversubscribed row is yours to avoid.
- **Watch progress** without the log:

  ```sh
  tail -f bench.log
  ls -lt benchmark_campaigns/          # hardware_manifest.json appears at once
  ls -l  /scratch/$USER/lsdtt3_fixtures  # grows during the first download
  ```

Start with the shakedown, exactly as under Slurm — it costs ten minutes and
catches a missing binary or an unreachable download before a multi-hour run:

```sh
pixi run --manifest-path /path/to/lsdtt3/pixi.toml -e benchmark \
    python /path/to/lsdtt3/tests/scripts/run_benchmark_suite.py \
        --profile laptop --build-dir /path/to/lsdtt3/build \
        --out-dir "$PWD/benchmark_campaigns"
```

Read the first lines of its output before anything else: all three references
must say `yes`. A missing one is skipped silently and you get a campaign with a
third of the rows.

### Route B: Slurm

The route for anything quotable: `#SBATCH --exclusive` is already in the
wrapper, and an exclusive node is the only way to keep a co-tenant out of your
timings.

**Submit from your data directory.** Everything else is derived:

| | derived from | value |
|---|---|---|
| repository | the script's own path | `/path/to/lsdtt3` |
| binaries | `<repo>/build` | `/path/to/lsdtt3/build` |
| pixi environment | `<repo>/pixi.toml`, `-e benchmark` | the one from step 1 |
| campaign output | the **submit** directory | `<submit>/benchmark_campaigns/` |
| DEM cache | the **submit** directory | `<submit>/benchmark_fixtures/` |
| job logs | the **submit** directory | `lsdtt3-bench-<jobid>.out` / `.err` |

Nothing is written into the clone.

### Three Slurm runs, smallest first

Do these in order. Each one proves something the next depends on, and the first
two cost minutes rather than hours.

#### a. The plumbing test — one case, one region, ~2 minutes

The smallest submission that exercises the entire chain: your account, the
scheduler, the repository detection, the pixi environment, the binaries, the
DEM download and the output paths. It measures one operation on one small DEM,
so the numbers are not interesting — **the point is that it finishes.**

```sh
cd /path/to/my/data
BENCH_PROFILE=laptop sbatch --account=PROJECT \
    /path/to/lsdtt3/tests/scripts/slurm_benchmark.sbatch \
    --cases fill --regions alpine --ladder 500 --threads 1 --repeats 2
```

Arguments after the script name pass straight through to the suite, so any of
`--cases`, `--regions`, `--ladder`, `--threads`, `--repeats` can be narrowed
this way.

When it finishes, read the **top** of `lsdtt3-bench-<jobid>.out` before
anything else:

```
repository      : /path/to/lsdtt3
binaries        : /path/to/lsdtt3/build
campaign output : /path/to/my/data/benchmark_campaigns
fixture cache   : /path/to/my/data/benchmark_fixtures

reference implementations importable in this environment:
    richdem              yes
    topotoolbox          yes
    whitebox_workflows   yes
    PROJ database        yes
```

Every path should be where you expect, and all three references should say
`yes`. A `NO` here is the one failure that looks like success — see
[Troubleshooting](#troubleshooting).

#### b. The laptop profile — ~10 minutes, ~8 MB

The full `laptop` profile: both regions, 500² and 1000² DEMs, threads 1/2/4,
every case. This is a real baseline you can keep, and it confirms the whole
suite runs under the scheduler before you commit to a long job.

```sh
cd /path/to/my/data
BENCH_PROFILE=laptop sbatch --account=PROJECT \
    /path/to/lsdtt3/tests/scripts/slurm_benchmark.sbatch
```

Expect most rows to be flagged `overhead_dominated` — the DEMs are small. That
is correct behaviour, not a problem. See
[Reading the results](#reading-the-results).

#### c. The server profile — hours, ~512 MB

The quotable run. `server` is the default, so there is nothing to set:

```sh
cd /path/to/my/data
sbatch --account=PROJECT /path/to/lsdtt3/tests/scripts/slurm_benchmark.sbatch
```

Raising the thread sweep means raising the allocation to match, or the job
refuses to start:

```sh
sbatch --account=PROJECT --cpus-per-task=32 \
    /path/to/lsdtt3/tests/scripts/slurm_benchmark.sbatch --threads 1,2,4,8,16,32
```

!!! tip "Point the DEM cache at scratch"
    ~512 MB, reused by every later campaign, so somewhere persistent beats a
    quota'd home directory:
    ```sh
    export BENCH_FIXTURE_DIR=/scratch/$USER/lsdtt3_fixtures
    ```

!!! warning "Directories you create are not used unless you say so"
    Making a `bench_results/` directory does not make the campaign go there.
    The defaults are `<submit>/benchmark_campaigns` and
    `<submit>/benchmark_fixtures`, and your own directories are simply ignored.
    To use them, name them:
    ```sh
    export BENCH_OUT_DIR=$PWD/bench_results
    export BENCH_FIXTURE_DIR=$PWD/lsdtt3_fixtures
    ```

### Where results go, exactly

**Results follow the directory you submit from, not where the script lives.**
Slurm sets `SLURM_SUBMIT_DIR` to your working directory at `sbatch` time, and
both the campaign and the DEM cache default there.

| you run | results land in |
|---|---|
| `cd /path/to/my/data && sbatch /path/to/lsdtt3/.../slurm_benchmark.sbatch` | `/path/to/my/data/` |
| `cd /path/to/lsdtt3/tests/scripts && sbatch slurm_benchmark.sbatch` | `/path/to/lsdtt3/tests/scripts/` |

!!! warning "Submitting from the script's own directory writes into the repository"
    The second row puts `benchmark_campaigns/` **and ~512 MB of cached DEMs**
    inside `tests/scripts/`, in your source tree. The job warns when it detects
    this but does not stop. Submit from a data directory, or set `BENCH_OUT_DIR`
    and `BENCH_FIXTURE_DIR` explicitly.

The same rule applies without Slurm: `$PWD` stands in for the submit directory,
so `cd` to where you want the results and run the suite directly.

**This rule covers `run_benchmark_suite.py` only.** The other two drivers
default into `build/` regardless of where you run them from, and must be given
`--out-dir` and `--fixture-dir` explicitly — see
[the other two drivers](#putting-their-output-in-your-campaign-directory).

### Under Slurm

`tests/scripts/slurm_benchmark.sbatch` contains **no measurement logic**. It
just arranges for a scheduler node to run the same command everyone else runs,
so results from Slurm and non-Slurm sites stay comparable. If your cluster uses
PBS, LSF or SGE, ignore the file and run the command above in your own job
script — you get an identical campaign directory.

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
pixi run -e benchmark python tests/scripts/run_benchmark_suite.py --profile server \
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

## Running only the thread sweep

The thread sweep is a small part of a campaign but it sits inside the longest
one. If you only need the scaling numbers — because a full `server` run already
gave you good single-thread timings, or because a previous sweep was invalid
(see the warning below) — run the `surface` suite on its own:

```sh
pixi run -e benchmark python tests/scripts/run_benchmark_suite.py \
    --profile server --suites surface
```

**Why that is the whole thread experiment.** Every case outside the `surface`
suite is declared `sequential` and runs once, at one thread, by design: filling,
breaching and flow-graph construction are ordered traversals with no
parallelism to measure. `lsdtt-surface-metrics` is the only genuinely parallel
CLI in the campaign, so `--suites surface` covers 100% of what the sweep
measures and skips the hours spent on cases that are pinned at one thread
anyway.

Narrow it further if you want a specific answer:

```sh
# one case, the full thread ladder, one region -- minutes, not hours
pixi run -e benchmark python tests/scripts/run_benchmark_suite.py \
    --profile server --cases polyfit_all --regions alpine

# just the top of the ladder, where the scaling signal is strongest
pixi run -e benchmark python tests/scripts/run_benchmark_suite.py \
    --profile server --suites surface --ladder 2000,4000

# check the plan before committing the time
pixi run -e benchmark python tests/scripts/run_benchmark_suite.py \
    --profile server --suites surface --dry-run
```

Use `--threads` to override the ladder — but cap it at your **physical** core
count, for the reason in [§4 above](#4-scaling-curves-need-a-ceiling-and-the-big-dems).

### Pick a rung with enough compute in it

Whole-process wall time includes GDAL read and write, and the surface-metric
kernels are fast. On a 2400² fixture the compute is a small fraction of the
run, so I/O dominates and the scaling signal is buried no matter how many
threads you give it. That is why the `server` ladder goes to 4000 px: at 16M
cells there is enough work to measure, and on the first server campaign the
surface-metric scaling curves at 8000² matched those at 4000². If you shorten
the ladder for a quick answer, shorten it from the *bottom*
(`--ladder 2000,4000`), never to the small rungs.

The per-run `Compute time` and `Write time` the CLI prints are captured in the
results, so you can check the ratio for yourself rather than guessing.

### If your campaign predates the n_threads fix

!!! warning "Thread sweeps from before 2026-09-09 are invalid"

    The harness set thread count **only** through the `OMP_NUM_THREADS`
    environment variable. That is not sufficient. A CLI that declares
    `n_threads` passes it into the library, and `resolve_threads()`
    (`src/cpu/lsdtt_raster_surface_fitting_cpu.cpp`) honours any positive value
    verbatim, applying it as an explicit `num_threads(n)` clause that
    **overrides the environment**. Only `n_threads: 0` falls through to
    `omp_get_max_threads()`, which is what actually reads `OMP_NUM_THREADS`.

    `lsdtt-surface-metrics` defaulted to `n_threads: 1` at the time. Since the harness never
    wrote `n_threads` into the generated `.param` file, every surface-metrics
    measurement ran **single-threaded at every thread count**, and the resulting
    scaling curves are flat by construction rather than by measurement.

    Because the `surface` suite is the only suite that is swept, this invalidated
    the entire thread experiment of any earlier campaign while leaving every
    single-thread timing untouched and still valid.

    The harness now writes `n_threads` into every generated parameter file. If
    you have a campaign from before the fix, you do **not** need to repeat it:
    re-run `--suites surface` as above and substitute the new scaling table. The
    non-surface numbers in your existing `SUMMARY.md` are unaffected.

Current `lsdtt-surface-metrics` and `lsdtt-spectral-analysis` defaults use
`n_threads: 0` (all hardware threads). Benchmark and shared-host runs should
still state an explicit physical-core limit rather than relying on that default.

---

## The other two drivers

The campaign above is the cross-package speed comparison. Two smaller drivers
answer questions it cannot, and both run the same way on a server as anywhere
else.

### Putting their output in your campaign directory

**Read this before running either one.** Unlike `run_benchmark_suite.py`, these
two do **not** create a timestamped campaign directory for you, and their
defaults write *inside the repository*:

| driver | results default to | fixtures default to |
|---|---|---|
| `run_benchmark_suite.py` | `./benchmark_campaigns/<profile>_<host>_<stamp>/` | `<out-dir>/benchmark_fixtures/` |
| `run_flow_benchmark.py` | `build/flow_benchmark/` | `tests/data/dem_cache/` |
| `run_svector_benchmark.py` | `build/svector_benchmark/<fixture>_<stamp>/` | `tests/data/dem_cache/` |

Both take `--out-dir` and `--fixture-dir` to override that. Pointing them at
the directory the suite already made keeps one run in one place, and — the part
that actually matters on a server — makes all three drivers **share one fixture
cache**, so the DEMs are downloaded once rather than once per driver.

Say your campaign is the one the suite just printed:

```sh
CAMPAIGN=$PWD/benchmark_campaigns/server_myhost_20260909_101147
FIXTURES=$PWD/benchmark_campaigns/benchmark_fixtures
```

Then:

```sh
# correctness -> $CAMPAIGN/flow_benchmark/
pixi run -e benchmark python tests/scripts/run_flow_benchmark.py \
    --exe build/lsdtt-raster-preprocessing \
    --drainage-exe build/lsdtt-drainage-analysis \
    --fetch-exe build/lsdtt-fetch-raster \
    --out-dir "$CAMPAIGN/flow_benchmark" \
    --fixture-dir "$FIXTURES" \
    --check

# s-vector payoff -> $CAMPAIGN/svector_benchmark/
pixi run -e benchmark python tests/scripts/run_svector_benchmark.py \
    --build-dir build \
    --out-dir "$CAMPAIGN/svector_benchmark" \
    --fixture-dir "$FIXTURES" \
    --depths 1,2,3,4,5
```

Afterwards the campaign holds everything:

```text
benchmark_campaigns/
├── benchmark_fixtures/              # shared by all three drivers
└── server_myhost_20260909_101147/
    ├── SUMMARY.md                   # from run_benchmark_suite.py
    ├── benchmark_raw.csv
    ├── flow_benchmark/
    │   ├── benchmark_results.csv    # value comparisons
    │   └── <fixture>_<case>/        # staged rasters and reference caches
    └── svector_benchmark/
        ├── svector_benchmark.csv
        └── run_manifest.json
```

Three things to know about these flags:

- `--out-dir` on `run_flow_benchmark.py` is **results and scratch together**.
  It holds `benchmark_results.csv` alongside staged filled DEMs and one output
  cache per reference package, which runs to hundreds of MB on the larger
  fixtures. Give it room. `--work-dir` is the original spelling of the same
  flag and still works.
- `--out-dir` on `run_svector_benchmark.py` is used **verbatim** — it does not
  append a timestamp the way the default does. Re-running with the same path
  appends to the existing `svector_benchmark.csv` rather than replacing it. Use
  a fresh path per run, or delete the CSV first.
- `--fixture-dir` on `run_flow_benchmark.py` is where fixtures are **fetched
  to** as well as read from, so it needs to be writable. On
  `run_svector_benchmark.py` it is read-only: that driver never downloads. Its
  `basins` and `chi` cases need the fixture to exist already — run the flow
  driver once first — and its `chi_value` case generates its own synthetic
  island and ignores `--fixture-dir` entirely.

### Correctness, not speed

```sh
pixi run -e benchmark python tests/scripts/run_flow_benchmark.py \
    --exe build/lsdtt-raster-preprocessing \
    --drainage-exe build/lsdtt-drainage-analysis --check
```

Compares lsdtt3's **values** against RichDEM, TopoToolbox and Whitebox: fill,
breach and D8 accumulation. `--check` enforces tolerances pinned just above
observed divergence, so it exits non-zero on a real disagreement.

Deterministic, so it does not need an idle machine — run it anywhere. It skips
(exit 77) if the fixture hashes do not match the ones the tolerances were pinned
against, because applying a pinned tolerance to a different raster asserts
nothing.

### What the s-vector buys

```sh
pixi run -e benchmark python tests/scripts/run_svector_benchmark.py \
    --build-dir build --depths 1,2,3,4,5
```

(Add `--out-dir` and `--fixture-dir` as above to keep the output with the rest
of the campaign; without them this writes to `build/svector_benchmark/`.)

lsdtt3 stores each cell's upslope set as a contiguous interval. That costs about
57 % of the D8 build and ~24 of the ~43 bytes per cell; the campaign measures
the cost. This driver measures the return, in the two places it is collected:

- **`basins_nested`** — nested basins swept over **nesting depth**. Not basin
  count: every tool does N *disjoint* basins in one pass, so that axis shows
  nothing. Depth is what the references pay for repeatedly (Whitebox's
  `unnest_basins` emits one raster per level), and what lsdtt3 reads as a
  sub-interval. **The result is the slope against depth, not any single time.**
- **`chi_network`** — DEM to chi, end-to-end and decomposed by phase. lsdtt3
  has no separate accumulation phase because the s-vector already produced
  drainage area; TopoToolbox must run `flow_accumulation` before
  `chitransform` can be called at all. That empty cell is the finding.
- **`chi_value`** — a VALUE comparison of chi against TopoToolbox, with
  `--check`. Deterministic and wired into CTest as
  `benchmark_chi_value_topotoolbox`, so it needs no idle machine.

```sh
pixi run -e benchmark python tests/scripts/run_svector_benchmark.py \
    --build-dir build --cases chi_value --check
```

It generates its own fixture — a synthetic island — and verifies the DEM's
sha256 before applying any tolerance. It cannot use a cropped real DEM at all:
lsdtt3 drops any basin whose upslope area touches invalid data, which on a crop
is every basin, so there would be no chi to compare.

### Enabling these under CTest

```sh
cmake -S . -B build -DLSDTT_ENABLE_BENCHMARKS=ON
ctest --test-dir build -L benchmark
```

OFF by default: they need network access for fixtures and the pixi `benchmark`
environment, neither of which may turn CI red. Anything unavailable exits 77 and
CTest records a skip rather than a failure.

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

**The log is empty and it looks hung.** Almost always stdout buffering, not a
stall. The suite line-buffers its own output, but if you are on a version that
predates that, or something else in the chain buffers, force it:

```sh
PYTHONUNBUFFERED=1 nohup pixi run bench-laptop --out-dir /abs/path &
```

To check progress without the log, look at the campaign directory — the
`hardware_manifest.json` appears almost immediately, and the DEM cache grows
during the first run:

```sh
ls -lt /abs/path/benchmark_campaigns/
ls -l  /abs/path/benchmark_fixtures/
```

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
  [Testing and benchmarking](../developer_notes/testing_and_benchmarking.md).
