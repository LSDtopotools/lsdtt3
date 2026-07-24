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
| **dethier2014** (pilot) | S169 | Front Range, CO | 14 | 0.07–10.11 | 1962–3480 | **−6.98 to −2.27 % (mean −5.59 %)**; 12 reliable |
| **clapp2002** | S014 | Yuma Wash, AZ (arid, low-relief) | 21 | 0.48–29.1 | 230–436 | **−5.49 to −3.57 % (mean −4.09 %)**; all 21 reliable |
| **dibiase2010** | S055 | San Gabriel Mts, CA (steep) | 15 | 0.08–6.4 | 1533–2079 | **−11.45 to +3.08 % (mean −5.25 %)**; 14 reliable |

All frozen **2026-07-24** (keyless COP30, St scaling, snap threshold 250 px), AFTER
the standardisation-double-count fix (below).

**All three datasets now share ONE uniform, single-signed offset: lsdtt3 ~4–6 %
LOWER than OCTOPUS.** That offset is the single *intentional* default difference vs
OCTOPUS's CAIRN — the **muon scheme** (lsdtt3 `BraucherBorchers` vs OCTOPUS
`Braucher`, ~5.65 % lower production). When lsdtt3 is **matched** to OCTOPUS
(`muon_scheme: Braucher`, `atmosphere: NCEP`, 07KNSTD) it reproduces OCTOPUS/CAIRN
erate to **+0.38 %** and basin-mean production to **−0.19 %** vs `BEPROD` — essentially
exact (the §6e isolation test). `clapp2002` (low-relief) shows the muon offset nearly
full (~−4 %, topo shielding ≈ 1); `dibiase2010` (steep) has larger per-basin scatter
from ignored per-pixel topo shielding. Sub-snap-threshold basins are flagged
`delineation_reliable=no` (`S169WTS036` +71 %, `S169WTS032`, `S055WTS044`) and
excluded from the band — a delineation limit, not a numeric bug.

**A standardisation double-count (harness bug) was found and fixed 2026-07-24.** The
earlier bands (`dethier` |Δ|≤9.3 %, `clapp` **+6.6 %**, `dibiase` ~−5 %) and the
apparent **sign flip** between `clapp` and `dibiase` were an artefact: OCTOPUS's
`BE10NC` is *already renormalised to 07KNSTD*, while `BESTND` is the original
source-paper standard kept as provenance (OCTOPUS docs; Codilean 2018/2022 ESSD). The
harness was harvesting the already-07KNSTD `BE10NC` **and** passing the native `BESTND`
label, so lsdtt3 (correctly, like CAIRN `LSDCosmoData.cpp:313`) applied the AMS factor a
*second* time, inflating erate by exactly `1/factor − 1` on every non-07KNSTD basin
(KNSTD 0.9042 ⇒ +10.6 %, LLNL3000 0.8644 ⇒ +15.7 %). The fix is harness-only —
`run_cosmo_validation.py write_crn_input()` now passes `07KNSTD`; the numeric-core
`be10_std_map_` is unchanged. Removing the double-count collapsed the sign-dependent
scatter to the clean uniform muon offset above. The three *calculator* differences
(muon 5.65 %, topo −0.7 %, atmosphere +0.64 %) remain real but small; see each
manifest's `residual_expectation`, the pilot `calculator_differences` +
`standardisation_convention` + `decisive_result` blocks, and the plan §6e / §9b log.
