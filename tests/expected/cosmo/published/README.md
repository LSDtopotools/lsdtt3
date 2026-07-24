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

## Adding a dataset (fully data-driven, no code edits)

A new dataset needs only (1) a harvest and (2) a networked freeze — no changes to
`run_cosmo_validation.py`:

```
# 1. harvest an OCTOPUS study's basins+outlets into the tracked artefacts
python3 tests/scripts/harvest_octopus.py --study S055 --dataset dibiase2010 \
    --out tests/expected/cosmo/published \
    --max-area-km2 8 --max-rel-uncert 0.5 --subset-bbox=-118.06,34.27,-117.95,34.39
#   -> writes <name>_samples.csv (real OCTOPUS rows + provenance header),
#      <name>_manifest.json (skeleton; hand-edit description + citations),
#      <name>_summary.csv (empty frozen skeleton)

# 2. freeze the summary from a real accepted run, then wire into CTest
#    (add_test cosmo_realworld_<name>, label cosmo_realworld, --check, skip-77)
```

The harness reads `<name>_manifest.json` for the DEM bbox, target UTM, and
per-dataset run knobs (`run_params`, e.g. `threshold_contributing_pixels`, or
`atmosphere: NCEP` to match OCTOPUS's atmosphere — default ERA40). `all` mode
(`run_cosmo_validation.py all`) discovers every `*_manifest.json`.

## Status — three datasets, chosen to span terrain and elevation

| dataset | OCTOPUS | region | basins | area (km²) | elev (m) | frozen band vs OCTOPUS |
|---|---|---|---|---|---|---|
| **dethier2014** (pilot) | S169 | Front Range, CO | 14 | 0.07–10.11 | 1962–3480 | \|Δ\| ≤ 9.3 % (mean ~5 %); 12 reliable |
| **clapp2002** | S014 | Yuma Wash, AZ (arid, low-relief) | 21 | 0.48–29.1 | 230–436 | **+5.1 to +7.3 % (mean +6.6 %)**; all 21 reliable |
| **dibiase2010** | S055 | San Gabriel Mts, CA (steep) | 15 | 0.08–6.4 | 1533–2079 | **−11.5 to +3.1 % (mean ~−5 %)**; 14 reliable |

All frozen **2026-07-24** (keyless COP30, St scaling, snap threshold 250 px).

**The two new datasets deliberately CONTRAST the pilot** — `clapp2002` is
low-relief / near-sea-level arid desert (topo shielding ~1.0, high-elevation
S_tot scaling absent), `dibiase2010` is steep high-relief with the largest topo-
shielding departures and erosion rates an order of magnitude faster. The residual
**flips sign** between them (clapp `+6.6 %`, dibiase `~−5 %`), a genuine
terrain/elevation-dependent divergence that feeds the parallel numerics
investigation. Agreement is near-CAIRN, not exact, for the same three documented
config differences as the pilot (muon scheme, topographic shielding, atmosphere;
see each manifest's `residual_expectation` and the pilot `calculator_differences`
block). Sub-snap-threshold basins are flagged `delineation_reliable=no` (pilot
`S169WTS036` +71 %, dibiase `S055WTS044`) and excluded from the band — a
delineation limit, not a numeric bug.

**Why not exact — three config differences vs OCTOPUS's CAIRN** (quantified in the
manifest `calculator_differences` block, Simon-directed 2026-07-24): (1) muon
scheme — OCTOPUS Braucher vs lsdtt3 **BraucherBorchers** ⇒ lsdtt3 ~5.65 % lower
production; (2) **topographic shielding** — OCTOPUS applies it (per-basin
`octopus_topo_shield`), lsdtt3 ignores it ⇒ lsdtt3 ~+0.7 % (up to +6 % steep);
(3) atmosphere — OCTOPUS **NCEP** vs lsdtt3 **ERA40**, not yet quantified. (1) and
(2) oppose and partly cancel; matching both overshoots to +3.9 %, and the per-basin
scatter is spatial (DEM/atmosphere), not these offsets. See the plan §9b findings log.
