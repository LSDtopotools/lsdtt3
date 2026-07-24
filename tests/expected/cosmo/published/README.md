# Real-landscape cosmogenic validation tier (plan §5b / Phase 6)

This directory holds the **small, tracked** artefacts of the "download → run →
delete" real-landscape validation tier for the cosmogenic-nuclide CLIs. See
`docs/plan_docs/cosmo_cli_testing_plan.md` §5b and §6c(B).

**This tier is NOT in CI.** It needs network, time, and possibly an
OpenTopography API key. It is invoked manually via
`tests/scripts/run_cosmo_validation.py <dataset|all>` (CTest label
`cosmo_realworld`, deliberately kept out of the default `unit;cosmo` suite).

## What lives here (small, few kB, reviewed in diffs)

For each dataset `<name>`:

- `<name>_manifest.json` — names the DEM source (provider, dataset id, bounding
  box, resolution) and where the authoritative per-sample regression targets
  live (the LSDTopoTools GitHub data repositories, per CAIRN §9 — **not** the
  paper tables). The harness reads this to know what to download and where to
  find the sample table.
- `<name>_samples.csv` — the per-sample table: lat/lon, measured nuclide
  concentration, 1σ uncertainty, AMS standard, and the **published**
  basin-averaged erosion rate plus the metadata of the calculator that produced
  it (original author, calculator name/version). This is the comparison target.
- `<name>_summary.csv` — the **frozen** summary of the last accepted validation
  run: per-sample lsdtt rate vs published rate and the percent difference.
  Comparing a fresh run against this frozen summary is itself a cheap tracked
  regression check; it only needs regenerating when the code or the manifest
  changes (a reviewable event, like the CRONUS reference).

## What NEVER lives here (large)

DEMs (hundreds of MB to GB) and any intermediate rasters. The harness downloads
the DEM into a **temporary working directory**, runs the comparison, appends the
summary, then **deletes the working directory — including on failure (cleanup
trap)**. A `--keep` flag exists for debugging a single dataset; the default is
always delete. Nothing large ever lands in the repo tree or in version control.

## Populating a sample table (IMPORTANT — do not fabricate)

The pilot `dethier2014_samples.csv` is a **documented template**: it carries the
exact schema and a `TODO` marker in place of real data. The authoritative
lat/lon, concentrations, standards and published rates must be copied from the
LSDTopoTools data repository named in the manifest (`data_repo` field), NOT
invented and NOT read off the paper figure. The harness refuses to run against a
table still carrying the `TODO` sentinel, so an un-populated template cannot
silently produce a bogus comparison.

## Status

Pilot dataset: **dethier2014** (Dethier et al. 2014, small low-relief Front
Range catchments — chosen per plan §5b/§6c(B) to shake out the harness before
the GB-scale Himalayan datasets). Authored as scaffolding; the sample table is a
template pending population from the data repo, and no accepted run exists yet
(no network/key in the authoring environment). See the plan §9b findings log.
