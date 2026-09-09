# Running benchmarks on your own system

lsdtt3 includes a benchmark suite that measures how it compares against three
other topographic analysis implementations — **RichDEM**, **TopoToolbox** and
**Whitebox Workflows** on *your* hardware, on real terrain.

It is deliberately not part of the test suite: the test suite is used to making sure the clis and routines do not break with updates and have the correct implementation. This is a different excercise. `ctest` never runs it  and it cannot fail a build. Its only product is a
self-contained directory of measurements you can use to see how fast and how memory thirsty the difference implemetaions are. 

If you want the reasoning behind the design — why particular comparisons are
valid and others are not — see
[Testing and benchmarking](testing_and_benchmarking.md). This page is just how
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

**On 16 GB the full `server` ladder should fit, but confirm rather than assume.**
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
--ladder 1000,2000,4000
```

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

## The other two drivers

The campaign above is the cross-package speed comparison. Two smaller drivers
answer questions it cannot, and both run the same way on a server as anywhere
else.

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
  [Testing and benchmarking](testing_and_benchmarking.md).
