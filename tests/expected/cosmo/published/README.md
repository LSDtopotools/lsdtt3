# Real-landscape cosmogenic validation tier (plan §5b / Phase 6)

This directory holds the **small, tracked** artefacts of the "download → run →
delete" real-landscape validation tier for the cosmogenic-nuclide CLIs. See
`docs/plan_docs/cosmo_cli_testing_plan.md` §5b and §6c(B).

**This tier is NOT in CI.** It needs network and time (the DEM fetch is
**keyless** — Copernicus GLO-30 via lsdtt-fetch-raster — so no API key). It is
invoked manually via `tests/scripts/run_cosmo_validation.py <dataset|all>`
(CTest label `cosmo_realworld`, deliberately kept out of the default `unit;cosmo`
suite).

## What lives here (small, few kB, reviewed in diffs)

For each dataset `<name>`:

- `<name>_manifest.json` — names the keyless DEM source (provider, bounding box,
  resolution, target UTM), the **OCTOPUS WFS query** the regression targets were
  harvested from, per-dataset run knobs (e.g. the channel-snap threshold), and
  the residual expectation. The harness reads this to know what to download and
  how to run.
- `<name>_samples.csv` — the per-basin table: outlet lat/lon, measured nuclide
  concentration, 1σ uncertainty, AMS standard, mean elevation, basin area, and
  the **published** basin-averaged erosion rate. This is the comparison target.
- `<name>_summary.csv` — the **frozen** summary of the last accepted validation
  run: per-basin lsdtt rate vs published rate, percent difference, elevation, and
  a `delineation_reliable` flag. Comparing a fresh run against this frozen summary
  (`--check`) is itself a cheap tracked regression check; it only needs
  regenerating when the code or the manifest changes (a reviewable event, like
  the CRONUS reference).

## What NEVER lives here (large)

DEMs (hundreds of MB to GB) and any intermediate rasters. The harness downloads
the DEM into a **temporary working directory**, runs the comparison, then
**deletes the working directory — including on failure (cleanup trap)**. A
`--keep` flag exists for debugging; the default is always delete. Nothing large
ever lands in the repo tree or in version control.

## Where the regression targets come from (do NOT fabricate)

The published rates and concentrations are **harvested from the OCTOPUS database**
(v2.3) via its WFS — not scraped from the paper and not invented. OCTOPUS
recalculates every rate with **CAIRN** (Mudd et al. 2016), the same algorithm
`lsdtt-cosmo-tool` implements, so these are a **near-CAIRN oracle, not bit-exact**
(lsdtt3 has since diverged — updated production parameters + a high-elevation
scaling fix this suite found). The exact WFS endpoint, `typeName`, `CQL_FILTER`,
and OCTOPUS DB version are recorded in each manifest's `octopus_source` block and
the sample-table header. The harness still refuses to run against a table
carrying the `TODO` sentinel, so an un-populated template can never silently
produce a bogus comparison.

## Running / re-freezing

```
# freeze / re-freeze the summary from a fresh accepted run:
python3 tests/scripts/run_cosmo_validation.py dethier2014 \
    --exe build/lsdtt-cosmo-tool --fetch-exe build/lsdtt-fetch-raster

# regression check against the frozen summary WITHOUT modifying it (ctest mode):
python3 tests/scripts/run_cosmo_validation.py dethier2014 --check \
    --exe build/lsdtt-cosmo-tool --fetch-exe build/lsdtt-fetch-raster

# offline harness self-test (cleanup trap + schema guards; in the default suite):
python3 tests/scripts/run_cosmo_validation.py --self-test
```

## Status

Pilot dataset: **dethier2014** (Dethier et al. 2014, Front Range). Pilot subset =
a compact Boulder Creek cluster of **14 small catchments (0.07–10.11 km²) spanning
1962–3480 m**, chosen per plan §5b/§6c(B) to shake out the harness and expose any
elevation-dependent residual before the GB-scale Himalayan datasets.

**First accepted run frozen 2026-07-24** (keyless COP30, St scaling, snap
threshold 250 px). The 12 delineation-reliable basins agree with OCTOPUS/CAIRN to
**|Δ| ≤ 9.3 % (mean ~5 %)** — tight for a near-CAIRN oracle, with no strong
elevation trend at COP30. One sub-threshold basin (`S169WTS036`, 0.07 km²) is
flagged `delineation_reliable=no` (+71 %: the outlet snaps onto a larger channel
than its true catchment — a delineation limit, not a numeric bug).

**Why not exact — three config differences vs OCTOPUS's CAIRN** (quantified in the
manifest `calculator_differences` block, Simon-directed 2026-07-24): (1) muon
scheme — OCTOPUS Braucher vs lsdtt3 **BraucherBorchers** ⇒ lsdtt3 ~5.65 % lower
production; (2) **topographic shielding** — OCTOPUS applies it (per-basin
`octopus_topo_shield`), lsdtt3 ignores it ⇒ lsdtt3 ~+0.7 % (up to +6 % steep);
(3) atmosphere — OCTOPUS **NCEP** vs lsdtt3 **ERA40**, not yet quantified. (1) and
(2) oppose and partly cancel; matching both overshoots to +3.9 %, and the per-basin
scatter is spatial (DEM/atmosphere), not these offsets. See the plan §9b findings log.
