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

## The LSDn tier — `dethier2014_lsdn`, `clapp2002_lsdn`, `dibiase2010_lsdn`

Added **2026-07-26** (LSDn plan Phase 6 deliverable (iii)): the same three landscapes
run with `scaling_scheme: lsdn` and compared against OCTOPUS's **`EBEL`** (LSDn)
instead of `EBE_*` (Stone/Lal).

**No data is duplicated.** Each `<name>_lsdn_manifest.json` points `sample_table` at
its St counterpart's **existing** `<name>_samples.csv` — the positions,
concentrations and uncertainties lsdtt-cosmo-tool consumes are identical between the
two tiers, so there is exactly one copy of every measured concentration in the
repository. Only the published *target* changes, supplied by the manifest's
`published_erate_override`, which reads `EBEL` (mm/kyr) from the tracked
`../octopus_lsdn_point_reference.csv` and converts with OCTOPUS's own density
`EBE_GCMYR/(EBE_MMKYR·1e-4)` = 2.650 g/cm³. A missing override file, a missing row, a
non-positive rate or an unexpected implied density is a **hard failure** — an absent
tracked input can never make the comparison report green.

| dataset | basins | frozen band vs `EBEL`, `delineation_reliable` only |
|---|---|---|
| **dethier2014_lsdn** | 14 | **+2.47 to +9.41 % (mean +4.51, median +4.28)**; 12 reliable |
| **clapp2002_lsdn** | 21 | **+2.61 to +4.66 % (mean +3.89, median +3.82)**; all 21 reliable |
| **dibiase2010_lsdn** | 15 | **+5.82 to +22.85 % (mean +11.45, median +10.14)**; 14 reliable |

**This is a MONITORING BAND, not an oracle.** It mixes the LSDn kernel residual
(pinned DEM-free by test **L12**), the per-pixel-vs-lumped difference (pinned by
**L13**, and a scientific *result*, not an error) and the DEM/delineation difference —
and it runs at the CLI default `sample_age_yr: 0` (modern 2010 field) while `EBEL` is
Balco v3 time-averaged over Γ/ε ≈ 12–38 kyr. **Mind the sign:** this tier reports
`(lsdtt3 − published)/published`, L12/L13 report the opposite. On this tier's sign the
matched-age expectation is ≈ **−6.6 %** (kernel −6.14 % on lsdtt3's own `P0`,
per-pixel/lumped −0.53 %, ERA40-vs-NCEP −0.64 %, ignored topo shielding +0.70 %), and
the ~+11-point offset to the frozen numbers is the time term — **measured**, not
argued: repeating the dethier run at `sample_age_yr` 5000 / 15000 / 20000 gives mean
−16.1 / −9.3 / −7.0 % on the same 12 basins. Never tune anything to shrink this band.

Independent cross-check: lsdtt3's per-basin LSDn/Stone *erate* ratio reproduces the
separately frozen scaling ratio (`mean_scaling_lsdn_t0 / mean_scaling_stone` in
`../lsdn_pixel_vs_lumped_perpixel.csv`) to better than one point on every basin.

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
