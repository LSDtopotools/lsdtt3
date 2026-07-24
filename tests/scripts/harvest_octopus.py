#!/usr/bin/env python3
"""Harvest an OCTOPUS study's basins + outlets into a tracked validation dataset.

Real-landscape cosmogenic validation tier (cosmo plan Phase 6 / §5b). This is the
small, committed helper that turns an OCTOPUS ``STUDYID`` into the two tracked
artefacts the download->run->delete harness (``run_cosmo_validation.py``) reads:

  * ``<dataset>_samples.csv`` -- the per-basin regression table (outlet lat/lon,
    measured 10Be concentration + 1sigma, AMS standard, snow/topographic
    shielding, mean elevation, atmospheric pressure, basin area, and the
    **published** basin-averaged erosion rate) HARVESTED from OCTOPUS, never
    invented and never paper-scraped;
  * ``<dataset>_manifest.json`` -- a skeleton naming the keyless DEM to fetch
    (bounding box = union of the selected basins' polygon bounds, warped to the
    OCTOPUS ``PROJEPSGID``), the OCTOPUS WFS provenance (endpoint, typeNames,
    ``CQL_FILTER``, DB version, retrieval date), and per-dataset run knobs.

OCTOPUS recalculates every rate with **CAIRN** (Mudd et al. 2016) -- the same
algorithm ``lsdtt-cosmo-tool`` implements -- so the harvested rates are a
**near-CAIRN oracle**, not a bit-exact target (see the manifest
``residual_expectation`` block). Adding a new dataset therefore needs only a run
of this helper plus a networked freeze of its summary; no per-dataset code edits.

Harvest specifics VERIFIED against the OCTOPUS server (2026-07-24):
  * WFS ``http://geoserver.octopusdata.org/geoserver/wfs``, ``version=2.0.0``,
    ``outputFormat=application/json``, ``srsName=EPSG:4326``.
  * Select a study with ``CQL_FILTER=STUDYID='SXXX'`` -- **NOT** the WFS ``bbox``
    parameter, which throws a ProjectionException on this server.
  * Basins in ``be10-denude:crn_int_basins`` (polygon + all fields); outlets in
    ``be10-denude:crn_int_outlets`` (Point), joined on ``OBSID1``.

Usage:
  python3 tests/scripts/harvest_octopus.py --study S014 --dataset clapp2002 \
      --out tests/expected/cosmo/published \
      [--max-area-km2 30] [--min-area-km2 0] \
      [--subset-bbox W,S,E,N] [--pad-deg 0.02] \
      [--target-epsg 32611] [--threshold-px 250] [--snap-radius 20]

The written manifest is a SKELETON: hand-edit its ``description`` and any
narrative caveats before freezing. The numeric fields (bbox, target_epsg,
octopus_source) are complete and correct as written.
"""

import argparse
import csv
import datetime
import json
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

WFS_ENDPOINT = "http://geoserver.octopusdata.org/geoserver/wfs"
BASINS_TYPE = "be10-denude:crn_int_basins"
OUTLETS_TYPE = "be10-denude:crn_int_outlets"

# Sample-table column order -- MUST match the pilot dethier2014_samples.csv so the
# harness (run_cosmo_validation.py) reads every dataset identically.
SAMPLE_FIELDS = [
    "sample_name", "latitude", "longitude", "Be10_CONC", "Be10_uncertainty",
    "standardisation", "snow_shielding", "octopus_topo_shield",
    "published_erate_gcm2yr", "published_erate_mm_kyr", "published_erate_mm_kyr_err",
    "elevation_ave_m", "atm_pressure_hpa", "octopus_area_km2",
    "published_calculator", "octopus_obsid",
]


def wfs_features(typename: str, studyid: str, timeout: int = 180) -> list:
    """GetFeature (GeoJSON) for one STUDYID via CQL_FILTER."""
    params = {
        "service": "WFS", "version": "2.0.0", "request": "GetFeature",
        "typeName": typename, "outputFormat": "application/json",
        "srsName": "EPSG:4326", "CQL_FILTER": f"STUDYID='{studyid}'",
    }
    url = WFS_ENDPOINT + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.load(r)["features"]


def _iter_coords(geom):
    """Yield (lon, lat) pairs from any GeoJSON geometry (Polygon/MultiPolygon)."""
    def walk(x):
        if (isinstance(x, list) and len(x) >= 2
                and all(isinstance(v, (int, float)) for v in x[:2])):
            yield (x[0], x[1])
        elif isinstance(x, list):
            for e in x:
                yield from walk(e)
    yield from walk(geom.get("coordinates", []))


def polygon_bounds(geom):
    lons = [c[0] for c in _iter_coords(geom)]
    lats = [c[1] for c in _iter_coords(geom)]
    return min(lons), min(lats), max(lons), max(lats)


def harvest(args) -> int:
    print(f"Harvesting OCTOPUS study {args.study} -> dataset '{args.dataset}'")
    basins = wfs_features(BASINS_TYPE, args.study)
    outlets = {f["properties"]["OBSID1"]: f
               for f in wfs_features(OUTLETS_TYPE, args.study)}
    if not basins:
        print(f"ERROR: no basins for STUDYID='{args.study}'", file=sys.stderr)
        return 1
    n_study = len(basins)

    sub = None
    if args.subset_bbox:
        sub = [float(v) for v in args.subset_bbox.split(",")]
        if len(sub) != 4:
            print("ERROR: --subset-bbox must be W,S,E,N", file=sys.stderr)
            return 1

    rows, poly_bounds, epsgs, dbvers = [], [], [], []
    for f in basins:
        p = f["properties"]
        oid = p["OBSID1"]
        if oid not in outlets:
            continue
        area = p.get("AREA")
        if area is None or area < args.min_area_km2 or area > args.max_area_km2:
            continue
        # Data-quality cut: drop samples whose analytical 1sigma is a large
        # fraction of the concentration (an unmeasurably noisy target, not a
        # meaningful regression point). Principled + documented, NOT cherry-picking.
        conc, cerr = p.get("BE10NC"), p.get("BE10NC_ERR")
        if (args.max_rel_uncert is not None and conc and cerr
                and cerr / conc > args.max_rel_uncert):
            print(f"  drop {oid}: rel. uncertainty {cerr/conc:.0%} > "
                  f"{args.max_rel_uncert:.0%}")
            continue
        lon, lat = outlets[oid]["geometry"]["coordinates"][:2]
        if sub and not (sub[0] <= lon <= sub[2] and sub[1] <= lat <= sub[3]):
            continue
        poly_bounds.append(polygon_bounds(f["geometry"]))
        epsgs.append(int(p["PROJEPSGID"]))
        dbvers.append(p.get("DBVERNO"))
        rows.append({
            "sample_name": oid,
            "latitude": f"{lat:.5f}",
            "longitude": f"{lon:.5f}",
            "Be10_CONC": p.get("BE10NC", ""),
            "Be10_uncertainty": p.get("BE10NC_ERR", ""),
            "standardisation": p.get("BESTND", ""),
            "snow_shielding": p.get("BESNOW", ""),
            "octopus_topo_shield": p.get("BETOPO", ""),
            "published_erate_gcm2yr": p.get("EBE_GCMYR", ""),
            "published_erate_mm_kyr": p.get("EBE_MMKYR", ""),
            "published_erate_mm_kyr_err": p.get("EBE_ERR", ""),
            "elevation_ave_m": round(p["ELEV_AVE"]) if p.get("ELEV_AVE") is not None else "",
            "atm_pressure_hpa": p.get("ATM_PRESS", ""),
            "octopus_area_km2": area,
            "published_calculator": "",  # filled after DB version resolved
            "octopus_obsid": oid,
        })
    if not rows:
        print("ERROR: no basins passed the area/bbox filters", file=sys.stderr)
        return 1

    # Sort by elevation so the residual-vs-elevation trend is easy to read.
    rows.sort(key=lambda r: (float(r["elevation_ave_m"] or 0), r["sample_name"]))

    db_version = Counter(dbvers).most_common(1)[0][0] or "unknown"
    calc = f"OCTOPUS_{db_version}_CAIRN_St_Braucher"
    for r in rows:
        r["published_calculator"] = calc

    target_epsg = args.target_epsg or Counter(epsgs).most_common(1)[0][0]

    pad = args.pad_deg
    west = min(b[0] for b in poly_bounds) - pad
    south = min(b[1] for b in poly_bounds) - pad
    east = max(b[2] for b in poly_bounds) + pad
    north = max(b[3] for b in poly_bounds) + pad

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    harvested = datetime.date.today().isoformat()

    samples_name = f"{args.dataset}_samples.csv"
    samples_path = out_dir / samples_name
    header = [
        f"# OCTOPUS study {args.study} 10Be basins -- REAL regression targets "
        f"(cosmo plan Phase 6).",
        "#",
        f"# HARVESTED {harvested} from the OCTOPUS database via WFS by "
        f"tests/scripts/harvest_octopus.py -- NOT scraped from the paper and NOT",
        "# invented. OCTOPUS recalculates every rate with CAIRN (Mudd et al. 2016) "
        "-- the",
        "# same algorithm lsdtt-cosmo-tool implements -- so these are a NEAR-CAIRN "
        "oracle,",
        "# not bit-exact (lsdtt3 has since diverged: updated production params + a "
        "high-",
        "# elevation S_tot>=1 scaling fix this suite found, plan §9b).",
        f"#   WFS: {WFS_ENDPOINT}",
        f"#        typeName={BASINS_TYPE} (+ {OUTLETS_TYPE}), "
        f"CQL_FILTER STUDYID='{args.study}'",
        f"#   OCTOPUS DBVERNO={db_version}   STUDYID={args.study}   "
        f"harvested {harvested}",
        "#",
        "# Scaling: EBE_* is OCTOPUS's default Stone/Lal (St) recalculation; run "
        "lsdtt-",
        "# cosmo-tool in the matching (St) scheme for a like-for-like compare.",
        "# Columns: latitude/longitude = OCTOPUS basin OUTLET point "
        "(crn_int_outlets);",
        "# published_erate_gcm2yr = EBE_GCMYR (the compared field); "
        "published_erate_mm_kyr",
        "# = EBE_MMKYR (+/- published_erate_mm_kyr_err = EBE_ERR); "
        "octopus_topo_shield = BETOPO;",
        "# snow_shielding = BESNOW; octopus_area_km2 = AREA; "
        "atm_pressure_hpa = ATM_PRESS.",
    ]
    with open(samples_path, "w", newline="") as f:
        for h in header:
            f.write(h + "\n")
        w = csv.DictWriter(f, fieldnames=SAMPLE_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {samples_path} ({len(rows)} basins of {n_study} in study)")

    summary_name = f"{args.dataset}_summary.csv"
    manifest = {
        "_comment": [
            "Real-landscape validation manifest (cosmo plan Phase 6 / §5b).",
            "SKELETON written by tests/scripts/harvest_octopus.py -- hand-edit the",
            "'description' and any narrative caveats; the numeric fields (bbox,",
            "target_epsg, octopus_source, run_params) are complete as written.",
            "The DEM is fetched KEYLESS (Copernicus GLO-30) into a temp workdir by",
            "run_cosmo_validation.py, compared, then deleted. No large data is",
            "ever committed.",
        ],
        "dataset": args.dataset,
        "description": (
            f"OCTOPUS study {args.study}: {len(rows)} small 10Be basins "
            f"(<= {args.max_area_km2} km2). TODO: describe region/relief/elevation "
            f"and the contrast with the pilot Front Range dataset."
        ),
        "published_reference": {
            "citation": "TODO: paper citation (author, year, journal).",
            "doi": "TODO",
            "calculator": (
                f"OCTOPUS {db_version} recalculation with CAIRN (Mudd et al. 2016), "
                "default Stone/Lal 'St' scaling (EBE_* fields)."
            ),
        },
        "octopus_source": {
            "_comment": [
                "Regression targets harvested from OCTOPUS, NOT scraped or invented.",
                "OCTOPUS recalculates with CAIRN = the algorithm lsdtt-cosmo-tool",
                "implements, so a NEAR-CAIRN oracle (see residual_expectation).",
            ],
            "wfs_endpoint": WFS_ENDPOINT,
            "typenames": [BASINS_TYPE, OUTLETS_TYPE],
            "cql_filter": f"STUDYID='{args.study}'",
            "cql_note": ("Filter with CQL_FILTER, not bbox: the OCTOPUS WFS 2.0.0 "
                         "bbox parameter throws a ProjectionException."),
            "output_format": "application/json",
            "srs": "EPSG:4326",
            "db_version": db_version,
            "study_id": args.study,
            "harvested_utc": harvested,
            "basins_in_study": n_study,
            "basins_in_subset": len(rows),
            "subset_bbox_deg": args.subset_bbox or "none (whole study, area-filtered)",
            "area_filter_km2": {"min": args.min_area_km2, "max": args.max_area_km2},
            "max_rel_uncertainty": args.max_rel_uncert,
            "scaling_scheme": ("St (Stone/Lal) -- EBE_MMKYR / EBE_GCMYR. OCTOPUS "
                               "also stores EBEL/EBELNS (LSDn)."),
        },
        "residual_expectation": {
            "_comment": [
                "Agreement will NOT be exact. lsdtt3 differs from OCTOPUS's CAIRN in",
                "muon scheme (BraucherBorchers vs Braucher), topographic shielding",
                "(lsdtt3 ignores it by default), and atmosphere (ERA40 vs NCEP), plus",
                "the high-elevation scaling fix this suite found. See the pilot",
                "dethier2014_manifest calculator_differences block for the quantified",
                "three-way breakdown.",
            ],
            "expect_zero": False,
            "band": ("measured-then-frozen on the FIRST accepted network run; do NOT "
                     "assert a tight zero-residual band."),
            "attribute_residual_to": ("muon-scheme + topographic-shielding + "
                                      "atmosphere + high-elevation scaling + DEM "
                                      "product (COP30 vs original) + delineation -- "
                                      "NOT bugs."),
        },
        "dem": {
            "provider": "Copernicus GLO-30 (keyless)",
            "dem_source": "cop30_aws",
            "dem_source_note": ("Keyless. lsdtt-fetch-raster retrieves+clips "
                                "Copernicus GLO-30 from AWS Open Data (cop30_aws) "
                                "or Microsoft Planetary Computer (planetary_computer) "
                                "with NO API key."),
            "pc_collection": "cop-dem-glo-30",
            "dataset_id": "COP30",
            "resolution_m": 30,
            "grid_spacing": 30.0,
            "resampling_method": "cubic",
            "target_epsg": target_epsg,
            "target_epsg_note": (f"UTM zone from OCTOPUS PROJEPSGID={target_epsg}. "
                                 "The cosmo CLIs require a PROJECTED CRS, so the DEM "
                                 "is warped to UTM at fetch time."),
            "resolution_choice_note": ("COP30 (30 m) vs the DEM in the OCTOPUS/CAIRN "
                                       "recalculation is a DOCUMENTED, attributable "
                                       "validation difference, not a bug."),
            "size_pixels": 0,
            "bounding_box": {
                "_comment": (f"Union of the {len(rows)} selected basins' OCTOPUS "
                             f"polygon bounds, padded {pad} deg."),
                "west": round(west, 4),
                "south": round(south, 4),
                "east": round(east, 4),
                "north": round(north, 4),
            },
        },
        "run_params": {
            "_comment": [
                "lsdtt-cosmo-tool run knobs. threshold_contributing_pixels lowered",
                "from the 1000 default (= 0.9 km2 at 30 m) so sub-km2 basins",
                "delineate near their true size; basins below this*cellsize^2 remain",
                "unreliable (outlet snaps onto a larger channel) and are flagged in",
                "the summary. To match OCTOPUS's atmosphere add 'atmosphere: NCEP'",
                "(default ERA40).",
            ],
            "threshold_contributing_pixels": args.threshold_px,
            "snap_search_radius_pixels": args.snap_radius,
        },
        "sample_table": samples_name,
        "summary": summary_name,
    }
    manifest_path = out_dir / f"{args.dataset}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {manifest_path} (target_epsg={target_epsg}, "
          f"bbox=[{west:.4f},{south:.4f},{east:.4f},{north:.4f}])")

    # Frozen-summary skeleton: header + column row only. The first accepted run of
    # run_cosmo_validation.py appends the data rows (a reviewable freeze).
    summary_path = out_dir / summary_name
    if not summary_path.exists():
        summary_fields = [
            "sample_name", "nuclide", "standardisation", "elevation_ave_m",
            "octopus_area_km2", "delineation_reliable", "lsdtt_erate_gcm2yr",
            "published_erate_gcm2yr", "published_erate_mm_kyr", "pct_diff",
            "abs_pct_diff", "published_calculator", "dem_source", "resolution_m",
            "run_date_utc", "git_sha",
        ]
        with open(summary_path, "w", newline="") as f:
            f.write(f"# {args.dataset} validation SUMMARY -- FROZEN skeleton "
                    f"(cosmo plan Phase 6).\n")
            f.write("# Written by run_cosmo_validation.py; first accepted networked "
                    "run appends one\n# row per basin (a reviewable freeze). "
                    "Difference is EXPECTED, not zero.\n")
            f.write(",".join(summary_fields) + "\n")
        print(f"wrote {summary_path} (empty frozen skeleton)")
    else:
        print(f"kept existing {summary_path} (not overwriting a frozen summary)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", required=True, help="OCTOPUS STUDYID, e.g. S014")
    ap.add_argument("--dataset", required=True,
                    help="dataset name for the tracked files, e.g. clapp2002")
    ap.add_argument("--out", default="tests/expected/cosmo/published",
                    help="output directory for the tracked artefacts")
    ap.add_argument("--max-area-km2", type=float, default=30.0)
    ap.add_argument("--min-area-km2", type=float, default=0.0)
    ap.add_argument("--max-rel-uncert", type=float, default=None,
                    help="drop samples with BE10NC_ERR/BE10NC above this fraction "
                         "(e.g. 0.5 = 50%%); an unmeasurably noisy target")
    ap.add_argument("--subset-bbox", default=None,
                    help="W,S,E,N outlet filter to select a compact cluster")
    ap.add_argument("--pad-deg", type=float, default=0.02)
    ap.add_argument("--target-epsg", type=int, default=None,
                    help="override the OCTOPUS PROJEPSGID for the DEM warp")
    ap.add_argument("--threshold-px", type=int, default=250)
    ap.add_argument("--snap-radius", type=int, default=20)
    args = ap.parse_args()
    try:
        return harvest(args)
    except Exception as exc:  # noqa: BLE001
        print(f"HARVEST FAILED: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
