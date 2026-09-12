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
      [--scaling St|LSDn|LSDnNS] \
      [--max-area-km2 30] [--min-area-km2 0] \
      [--subset-bbox W,S,E,N] [--pad-deg 0.02] \
      [--target-epsg 32611] [--threshold-px 250] [--snap-radius 20]

The written manifest is a SKELETON: hand-edit its ``description`` and any
narrative caveats before freezing. The numeric fields (bbox, target_epsg,
octopus_source) are complete and correct as written.

SECOND MODE -- ``--point-reference`` (lsdn plan Phase 6 deliverable (i), test L12)
---------------------------------------------------------------------------------
OCTOPUS does **not** compute LSDn from the DEM. Per Saktura et al. (2025, GChron
7:113) CAIRN reduces each basin to two scalars -- ``CENTR_LAT`` (centroid
latitude) and ``ATM_PRESS`` (the effective basin-averaged pressure, itself
iterated to match the spatially averaged **Stone--Lal** production) -- and feeds
them to the Balco v3 point calculator via RIVERSAND. ``EBEL`` is that point
result. So a *matched* LSDn comparison is DEM-FREE: evaluate our LSDn kernel at
exactly (``CENTR_LAT``, ``ATM_PRESS``), invert, compare to ``EBEL``.

One longitude is needed too (Holocene Rc is a distance to a MOVING paleopole),
and OCTOPUS publishes none. ``CENTR_LON`` in the frozen artefact is therefore
DERIVED BY US -- the area-weighted shoelace centroid of OCTOPUS's own published
basin polygon -- and is not a harvested field. The same computation's latitude is
written as ``CENTR_LAT_DERIVED``, and it does NOT reproduce the published
``CENTR_LAT`` (max 0.00863 deg, RMS 0.00200 deg over the 50 tracked basins; on
``S055WTS044`` the published 34.3838 lies outside the polygon's own latitude span
[34.37989, 34.38226], which no geometric centroid can do). The geometry is not at
fault -- those rings reproduce the published ``AREA`` to 9.2e-5 relative on
``S014WTS036`` -- so the leading, UNCONFIRMED hypothesis is that ``CENTR_LAT`` is
an *effective* latitude iterated against spatially averaged Stone--Lal
production, exactly as ``ATM_PRESS`` is. A separate research task is settling
that; ``CENTR_LAT_DERIVED`` is recorded so the discrepancy is frozen, and is
never asserted against ``CENTR_LAT`` anywhere.

``--point-reference`` harvests exactly that, for the basins already selected by
the tracked ``<dataset>_samples.csv`` files, into the small frozen
``tests/expected/cosmo/octopus_lsdn_point_reference.csv``:

  python3 tests/scripts/harvest_octopus.py --point-reference \
      [--datasets dethier2014,clapp2002,dibiase2010] \
      [--published-dir tests/expected/cosmo/published] \
      [--point-reference-out tests/expected/cosmo/octopus_lsdn_point_reference.csv] \
      [--freeze]

FROZEN-ARTEFACT DISCIPLINE: the default is CHECK-ONLY -- a fresh harvest is
compared against the tracked file and never writes it; the ``# HARVESTED ...``
provenance line is neither compared nor rewritten. ``--freeze`` regenerates
deliberately, for review. Same pattern as run_cosmo_production_isolation.py.

FIELD NAMES CONFIRMED 2026-07-26 against the live WFS ``DescribeFeatureType`` on
``be10-denude:crn_int_basins`` (the published field dictionary is incomplete; the
WFS schema is authoritative): ``EBEL``, ``EBEL_INT``, ``EBEL_EXT``, ``EBELNS``,
``EBELNS_INT``, ``EBELNS_EXT`` for 10Be, and ``EALL*`` / ``EALLNS*`` for 26Al.
``EBEL`` is in **mm/kyr**, the same units as ``EBE_MMKYR`` (median ratio 1.050
over S169's 46 basins); ``_INT``/``_EXT`` are Balco v3's internal/external
uncertainty pair. ``EBELNS`` is the **no-TOPOGRAPHIC-shielding** control:
``EBELNS >= EBEL`` always.

NS = NO TOPOGRAPHIC SHIELDING -- CORRECTED 2026-07-26. An earlier revision of
this file said "no snow shielding", inferred from ``EBELNS == EBEL`` on some
basins. That inference was wrong, and the data refute it: ``BESNOW`` is
identically 1.0000 on all 50 harvested basins, so it cannot explain any
difference at all. What the ratio actually tracks is topographic shielding --
``EBELNS/EBEL`` vs ``1/BETOPO`` agree to mean 4e-4 (max 4.5e-3, consistent with
EBEL being published to 3 significant figures) over the 50 basins. This also
matches Codilean & Munack (2025, GChron 7:113), which describes the pair as with
and without **topographic** shielding. The basins where ``EBELNS == EBEL`` are
simply those whose ``BETOPO >= 0.996``, where the difference falls below the
published precision.
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

# --scaling: which OCTOPUS 10Be rate fields populate published_erate_*.
#   St     -- EBE_*, OCTOPUS's own CAIRN recalculation with time-independent
#             Stone/Lal. The ONLY scheme CAIRN implements (Saktura et al. 2025,
#             GChron 7:113); the default, so the tracked St artefacts never move.
#   LSDn   -- EBEL, the Balco v3 point result at (CENTR_LAT, ATM_PRESS) obtained
#             through RIVERSAND, WITH topographic shielding.
#   LSDnNS -- EBELNS, the same WITHOUT topographic shielding. EBELNS >= EBEL
#             always, and EBELNS/EBEL tracks 1/BETOPO to mean 4e-4 over the 50
#             harvested basins. NOT snow shielding -- see the module docstring;
#             BESNOW is identically 1.0000 across every basin harvested here.
# EBEL/EBELNS are published in mm/kyr only; there is no EBEL_GCMYR field, so the
# g cm^-2 yr^-1 column is derived per basin with the density OCTOPUS itself used,
# read off its own EBE_GCMYR/EBE_MMKYR pair (see octopus_rock_density).
SCALING_FIELDS = {
    "St":     {"mmkyr": "EBE_MMKYR", "err": "EBE_ERR",     "gcmyr": "EBE_GCMYR"},
    "LSDn":   {"mmkyr": "EBEL",      "err": "EBEL_EXT",    "gcmyr": None},
    "LSDnNS": {"mmkyr": "EBELNS",    "err": "EBELNS_EXT",  "gcmyr": None},
}

# Fields written by --point-reference (lsdn plan Phase 6 deliverable (i), L12).
POINT_REFERENCE_FIELDS = [
    # keys
    "sample", "dataset", "study_id",
    # the two PUBLISHED scalars OCTOPUS reduces the basin to, and which L12
    # evaluates at
    "CENTR_LAT", "ATM_PRESS",
    # CENTR_LON is NOT an OCTOPUS field. The schema publishes CENTR_LAT with no
    # centroid longitude, and lsdn_scaling() needs a longitude because Holocene Rc
    # is the great-circle distance to a MOVING paleopole. So we DERIVE one here:
    # the area-weighted shoelace centroid of OCTOPUS's own basin polygon from
    # crn_int_basins (exterior rings positive, interior rings subtracted,
    # MultiPolygon parts area-weighted). Derived by us, clearly labelled, never
    # presented as harvested.
    "CENTR_LON",
    # CENTR_LAT_DERIVED is the latitude that falls out of the SAME computation.
    # It is RECORDED, NOT CHECKED, and it DOES NOT REPRODUCE the published
    # CENTR_LAT: over these 50 basins CENTR_LAT_DERIVED - CENTR_LAT reaches
    # 0.00863 deg with an RMS of 0.00200 deg, and on S055WTS044 the published
    # CENTR_LAT (34.3838) lies OUTSIDE the polygon's own latitude range
    # [34.37989, 34.38226] - a geometric centroid cannot do that, so the two
    # quantities are not the same quantity. Nothing anywhere asserts that they
    # agree; the column exists so the discrepancy is frozen and regresses.
    # The geometry itself is demonstrably right: the same rings reproduce the
    # published AREA to 9.2e-5 relative on S014WTS036 (22.3779 km2 computed vs
    # 22.38 published), and to better than 8e-4 on every basin large enough that
    # AREA's own 0.01 km2 quantisation does not dominate.
    # LEADING HYPOTHESIS, UNCONFIRMED: CENTR_LAT is an EFFECTIVE latitude,
    # iterated against the spatially averaged Stone-Lal production exactly as
    # ATM_PRESS is (Saktura et al. 2025, GChron 7:113, describes ATM_PRESS that
    # way), rather than a geometric centroid. A separate research task is running
    # to settle it; until it does, this stays a recorded discrepancy.
    "CENTR_LAT_DERIVED",
    # OUTLET_LON, the basin outlet longitude from crn_int_outlets, is kept as a
    # DIAGNOSTIC: it was the longitude L12 used before CENTR_LON existed, and
    # CENTR_LON - OUTLET_LON (max 0.0329 deg, RMS 0.0129 deg over these 50
    # basins) is the size of that former approximation.
    "OUTLET_LON",
    # the oracle: Balco v3 LSDn point rates (mm/kyr) + internal/external sigma
    "EBEL", "EBEL_INT", "EBEL_EXT", "EBELNS", "EBELNS_INT", "EBELNS_EXT",
    # what L12 inverts: the measured concentration (already renormalised, §9b)
    "BE10NC", "BE10NC_ERR", "BESTND",
    # shielding: EBEL carries TOPOGRAPHIC shielding, EBELNS does not (corrected
    # 2026-07-26 -- see the module docstring; BESNOW is 1.0000 on every basin
    # here, so it explains nothing, while EBELNS/EBEL tracks 1/BETOPO to 4e-4)
    "BETOPO", "BESELF", "BESNOW",
    # the St pair, which pins the mm/kyr -> g cm^-2 yr^-1 density conversion
    "EBE_GCMYR", "EBE_MMKYR",
    # context only, never entering the inversion: residual vs elevation/area
    "ELEV_AVE", "AREA",
]

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


def _ring_shoelace(ring):
    """(2*signed area, 6*Cx*signed area, 6*Cy*signed area) for one closed ring.

    Plain planar shoelace in (lon, lat) degrees.  Over basins of a few tens of
    km2 the degree-space distortion shifts the centroid by orders of magnitude
    less than the LSDn grid's Rc resolution, so no projection is warranted; the
    sign convention (CCW positive) is what lets interior rings subtract
    themselves automatically.
    """
    a2 = cx6 = cy6 = 0.0
    for (x0, y0), (x1, y1) in zip(ring, ring[1:]):
        cross = x0 * y1 - x1 * y0
        a2 += cross
        cx6 += (x0 + x1) * cross
        cy6 += (y0 + y1) * cross
    return a2, cx6, cy6


def polygon_centroid(geom):
    """Area-weighted geometric centroid (lon, lat) of a GeoJSON Polygon/MultiPolygon.

    Exterior rings contribute positively and interior rings (holes) subtract,
    which the signed shoelace does on its own provided the rings carry GeoJSON's
    winding; MultiPolygon parts combine by signed area, so the result is the
    centroid of the whole basin rather than of its largest part.  Degenerate
    (zero-area) geometries fall back to the mean vertex position.

    This is DERIVED BY US from OCTOPUS's own published basin polygon.  It is not
    a published OCTOPUS field, and it does NOT reproduce their published
    CENTR_LAT -- see the POINT_REFERENCE_FIELDS comment for the measured
    discrepancy and the leading hypothesis.
    """
    a2 = cx6 = cy6 = 0.0
    for poly in _iter_rings(geom):
        for ring in poly:
            if len(ring) < 4:
                continue
            r = ring if ring[0] == ring[-1] else list(ring) + [ring[0]]
            da2, dcx6, dcy6 = _ring_shoelace(r)
            a2 += da2
            cx6 += dcx6
            cy6 += dcy6
    if a2 == 0.0:
        pts = list(_iter_coords(geom))
        if not pts:
            return None, None
        return (sum(p[0] for p in pts) / len(pts),
                sum(p[1] for p in pts) / len(pts))
    return cx6 / (3.0 * a2), cy6 / (3.0 * a2)


def _iter_rings(geom):
    """Yield each polygon as a list of rings, for Polygon and MultiPolygon alike."""
    coords = geom.get("coordinates", [])
    gtype = geom.get("type", "")
    if gtype == "Polygon":
        yield coords
    elif gtype == "MultiPolygon":
        for poly in coords:
            yield poly


def octopus_rock_density(props: dict):
    """Rock density (g cm^-3) OCTOPUS used, read off its OWN St rate pair.

    EBEL/EBELNS are published in mm/kyr only.  Rather than hard-code a density to
    convert them, recover it from the basin's own EBE_GCMYR / EBE_MMKYR, which
    OCTOPUS publishes as the same rate in both units:

        EBE_GCMYR [g cm^-2 yr^-1] = EBE_MMKYR [mm/kyr] * 1e-4 [cm/yr per mm/kyr]
                                    * rho [g cm^-3]

    Over the 50 basins in the three tracked samples CSVs this gives 2.6474-2.6521
    (mean 2.6499), i.e. rho = 2.650 to the precision of the published 4-significant
    -figure EBE_GCMYR.  Returns None when either field is missing/zero.
    """
    try:
        g, m = float(props["EBE_GCMYR"]), float(props["EBE_MMKYR"])
    except (KeyError, TypeError, ValueError):
        return None
    return g / (m * 1.0e-4) if m else None


def harvest(args) -> int:
    scheme = SCALING_FIELDS[args.scaling]
    print(f"Harvesting OCTOPUS study {args.study} -> dataset '{args.dataset}' "
          f"(scaling {args.scaling})")
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
        # Selected rate pair (--scaling).  For LSDn/LSDnNS OCTOPUS publishes only
        # mm/kyr, so the g cm^-2 yr^-1 column is derived with OCTOPUS's own density.
        erate_mmkyr = p.get(scheme["mmkyr"], "")
        erate_err = p.get(scheme["err"], "")
        if scheme["gcmyr"]:
            erate_gcmyr = p.get(scheme["gcmyr"], "")
        else:
            rho = octopus_rock_density(p)
            erate_gcmyr = (f"{float(erate_mmkyr) * 1.0e-4 * rho:.6g}"
                           if rho and erate_mmkyr not in (None, "") else "")
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
            "published_erate_gcm2yr": erate_gcmyr,
            "published_erate_mm_kyr": erate_mmkyr,
            "published_erate_mm_kyr_err": erate_err,
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
    # The calculator differs BY SCHEME, and saying so matters: EBE_* is OCTOPUS's
    # own CAIRN run, but EBEL/EBELNS are Balco v3 point evaluations at
    # (CENTR_LAT, ATM_PRESS) reached through RIVERSAND -- a different code, a
    # different reduction of the basin, not a CAIRN scheme switch.
    calc = (f"OCTOPUS_{db_version}_CAIRN_St_Braucher" if args.scaling == "St"
            else f"OCTOPUS_{db_version}_BalcoV3_{args.scaling}_point_RIVERSAND")
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
        f"# SCALING SCHEME: {args.scaling}. Run lsdtt-cosmo-tool in the matching "
        f"scheme for a",
        "# like-for-like compare. St = EBE_* (OCTOPUS's own CAIRN, Stone/Lal -- the "
        "only",
        "# scheme CAIRN implements); LSDn = EBEL and LSDnNS = EBELNS, which are "
        "Balco v3",
        "# POINT evaluations at (CENTR_LAT, ATM_PRESS) via RIVERSAND, not CAIRN "
        "runs -- so",
        "# an LSDn dataset here mixes the kernel residual with a per-pixel/lumped "
        "difference",
        "# (see the DEM-free kernel oracle, octopus_lsdn_point_reference.csv).",
        "# Columns: latitude/longitude = OCTOPUS basin OUTLET point "
        "(crn_int_outlets);",
        f"# published_erate_mm_kyr = {scheme['mmkyr']} "
        f"(+/- published_erate_mm_kyr_err = {scheme['err']});",
        (f"# published_erate_gcm2yr = {scheme['gcmyr']} (the compared field);"
         if scheme["gcmyr"] else
         f"# published_erate_gcm2yr = {scheme['mmkyr']} x 1e-4 x rho (the compared "
         f"field), rho"),
        ("" if scheme["gcmyr"] else
         "# recovered per basin from OCTOPUS's own EBE_GCMYR/EBE_MMKYR pair "
         "(~2.650 g/cm3);"),
        "# octopus_topo_shield = BETOPO; snow_shielding = BESNOW; "
        "octopus_area_km2 = AREA;",
        "# atm_pressure_hpa = ATM_PRESS.",
    ]
    header = [h for h in header if h != ""]
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
                if args.scaling == "St" else
                f"OCTOPUS {db_version} {scheme['mmkyr']}: Balco et al. (2008) v3 "
                "POINT calculator via RIVERSAND, evaluated at the CAIRN-reduced "
                "(CENTR_LAT, ATM_PRESS) scalars -- NOT a CAIRN scheme (Saktura et "
                "al. 2025, GChron 7:113)."
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
            "scaling_scheme": args.scaling,
            "scaling_scheme_note": (
                "St (Stone/Lal) -- EBE_MMKYR / EBE_GCMYR, OCTOPUS's own CAIRN "
                "recalculation. OCTOPUS also stores EBEL (LSDn, topographically "
                "shielded) and EBELNS (LSDn, no topographic shielding), which are "
                "Balco v3 point results."
                if args.scaling == "St" else
                f"{args.scaling} -- {scheme['mmkyr']} (+/- {scheme['err']}), a Balco "
                "v3 POINT evaluation at (CENTR_LAT, ATM_PRESS). g cm^-2 yr^-1 "
                "derived with OCTOPUS's own EBE_GCMYR/EBE_MMKYR density (~2.650). "
                "This dataset therefore mixes the LSDn kernel residual with the "
                "per-pixel/lumped difference and the DEM difference: a monitoring "
                "band, not an oracle. The oracle is the DEM-free "
                "octopus_lsdn_point_reference.csv (test L12)."),
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
            "sample_snap_max_search_radius_pixels": args.snap_radius,
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


# --------------------------------------------------------------------------- #
# SECOND MODE -- --point-reference: the DEM-free L12 kernel oracle
# --------------------------------------------------------------------------- #

# Relative tolerance for the frozen-artefact regression check. Every harvested
# quantity is a stored OCTOPUS database field, so a fresh query is bit-identical
# unless the database itself moves; this absorbs only float formatting.
POINT_REF_CHECK_REL_TOL = 1e-9

# Compared as strings, not numbers: the AMS standard is a label, and a change in
# it changes what the concentration MEANS (§9b). It must be reported, not averaged.
POINT_REF_EXACT_COLS = ["sample", "dataset", "study_id", "BESTND"]


def read_frozen_point_reference(path: Path) -> dict:
    """Read the tracked artefact, skipping its '#' header. Keyed by dataset|sample."""
    if not path.is_file():
        return {}
    with open(path, newline="") as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    return {f"{r['dataset']}|{r['sample']}": r
            for r in csv.DictReader(lines) if r.get("sample")}


def check_against_frozen_point_reference(rows: list, out_path: Path) -> int:
    """Compare a fresh harvest against the tracked artefact WITHOUT modifying it.

    Same discipline as run_cosmo_production_isolation.py: the '# HARVESTED ...'
    provenance line moves on every run, so rewriting the file by default would
    dirty a tracked artefact without ever regression-checking its numbers. Only
    the numbers are compared here; --freeze regenerates deliberately, for review.
    """
    frozen = read_frozen_point_reference(out_path)
    if not frozen:
        print(f"CHECK: no frozen artefact at {out_path}; run once with --freeze to "
              "create a reference first.")
        return 1

    ok = True
    for r in rows:
        key = f"{r['dataset']}|{r['sample']}"
        ref = frozen.get(key)
        if ref is None:
            print(f"[CHECK][FAIL] {key}: not in the frozen artefact")
            ok = False
            continue
        for col in POINT_REFERENCE_FIELDS:
            new_s, ref_s = str(r.get(col, "")), str(ref.get(col, ""))
            if col in POINT_REF_EXACT_COLS:
                if new_s != ref_s:
                    print(f"[CHECK][FAIL] {key}.{col}: new={new_s!r} "
                          f"frozen={ref_s!r} (must match exactly)")
                    ok = False
                continue
            try:
                new_v, ref_v = float(new_s), float(ref_s)
            except ValueError:
                if new_s != ref_s:
                    print(f"[CHECK][FAIL] {key}.{col}: non-numeric and differing "
                          f"(new={new_s!r} frozen={ref_s!r})")
                    ok = False
                continue
            denom = abs(ref_v)
            rel = abs(new_v - ref_v) / denom if denom else abs(new_v - ref_v)
            if rel > POINT_REF_CHECK_REL_TOL:
                print(f"[CHECK][FAIL] {key}.{col}: new={new_v:.9g} "
                      f"frozen={ref_v:.9g} rel={rel:.3g}")
                ok = False

    missing = set(frozen) - {f"{r['dataset']}|{r['sample']}" for r in rows}
    if missing:
        print(f"[CHECK][FAIL] frozen basins absent from this harvest: "
              f"{', '.join(sorted(missing))}")
        ok = False

    print(f"CHECK: {'PASS' if ok else 'FAIL'} (rel tol {POINT_REF_CHECK_REL_TOL:g}, "
          f"{len(rows)} basins, artefact NOT modified)")
    if not ok:
        print("  If the change is intended, re-run with --freeze and review the diff.")
    return 0 if ok else 1


def read_selected_basins(published_dir: Path, dataset: str):
    """(study_id, [OBSID, ...]) for the basins already selected by a tracked dataset.

    The point reference deliberately reuses the SAME basin selection as the St
    datasets, so the kernel oracle (L12) and the per-pixel product tier are
    talking about the same landscapes and the two residuals are comparable.
    """
    manifest_path = published_dir / f"{dataset}_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    study = manifest["octopus_source"]["study_id"]
    samples_path = published_dir / manifest["sample_table"]
    with open(samples_path, newline="") as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    obsids = [r["octopus_obsid"] for r in csv.DictReader(lines)
              if r.get("octopus_obsid")]
    return study, obsids


def harvest_point_reference(args) -> int:
    """Harvest CENTR_LAT/ATM_PRESS/EBEL*/EBELNS* for the already-selected basins."""
    datasets = [d.strip() for d in args.datasets.split(",") if d.strip()]
    pdir = Path(args.published_dir)
    out_path = Path(args.point_reference_out)

    rows, studies = [], []
    for dataset in datasets:
        study, obsids = read_selected_basins(pdir, dataset)
        studies.append(f"{dataset}={study}")
        print(f"{dataset}: STUDYID {study}, {len(obsids)} selected basins")
        # The basin GEOMETRY is kept, not just the properties: CENTR_LON and
        # CENTR_LAT_DERIVED are computed from the polygon itself.
        by_oid = {f["properties"]["OBSID1"]: f
                  for f in wfs_features(BASINS_TYPE, study)}
        out_by_oid = {f["properties"]["OBSID1"]: f["geometry"]["coordinates"]
                      for f in wfs_features(OUTLETS_TYPE, study)}
        missing = [o for o in obsids if o not in by_oid or o not in out_by_oid]
        if missing:
            print(f"ERROR: OBSIDs absent from {study}: {missing}", file=sys.stderr)
            return 1
        for oid in obsids:
            p = by_oid[oid]["properties"]
            c_lon, c_lat = polygon_centroid(by_oid[oid]["geometry"])
            if c_lon is None:
                print(f"ERROR: degenerate basin polygon for {oid}", file=sys.stderr)
                return 1
            row = {"sample": oid, "dataset": dataset, "study_id": study,
                   "CENTR_LON": f"{c_lon:.5f}",
                   "CENTR_LAT_DERIVED": f"{c_lat:.5f}",
                   "OUTLET_LON": f"{out_by_oid[oid][0]:.5f}"}
            for col in POINT_REFERENCE_FIELDS:
                if col in row:
                    continue
                v = p.get(col)
                row[col] = "" if v is None else v
            rows.append(row)

    if not rows:
        print("ERROR: no basins harvested", file=sys.stderr)
        return 1

    if not args.freeze:
        return check_against_frozen_point_reference(rows, out_path)

    harvested = datetime.date.today().isoformat()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        f.write("# OCTOPUS LSDn POINT reference -- the DEM-free KERNEL oracle for "
                "test L12.\n")
        f.write("#\n")
        f.write("# This file validates the LSDn SCALING KERNEL, not lsdtt3's "
                "product. OCTOPUS does\n"
                "# not compute LSDn from a DEM: per Saktura et al. (2025, GChron "
                "7:113) CAIRN\n"
                "# implements only time-independent Stone-Lal, and for LSDn it "
                "reduces each basin\n"
                "# to two scalars -- CENTR_LAT (centroid latitude) and ATM_PRESS "
                "(the effective\n"
                "# basin-averaged pressure, itself iterated to match the spatially "
                "averaged\n"
                "# STONE-LAL production) -- and feeds them to the Balco et al. "
                "(2008) v3 point\n"
                "# calculator via RIVERSAND. EBEL is that point result. So the "
                "matched comparison\n"
                "# is DEM-FREE: evaluate our LSDn kernel at exactly (CENTR_LAT, "
                "ATM_PRESS), invert,\n"
                "# compare to EBEL. No DEM, no delineation, no COP30-vs-original "
                "confound; the\n"
                "# residual isolates the scaling kernel alone.\n")
        f.write("#\n")
        f.write("# CAVEAT, to quantify and NOT to fix here: ATM_PRESS is calibrated "
                "on Stone-Lal,\n"
                "# so it is strictly the WRONG effective pressure for LSDn, which "
                "has a different\n"
                "# pressure dependence. That is a systematic to measure (test L13, "
                "per-pixel vs\n"
                "# lumped), not truth.\n")
        f.write("#\n")
        f.write(f"# HARVESTED {harvested} from the OCTOPUS WFS by "
                f"tests/scripts/harvest_octopus.py --point-reference\n")
        f.write(f"#   WFS: {WFS_ENDPOINT}\n")
        f.write(f"#        typeName={BASINS_TYPE}, outputFormat=application/json, "
                f"srsName=EPSG:4326\n")
        f.write(f"#        CQL_FILTER STUDYID='...' per dataset: "
                f"{', '.join(studies)}\n")
        f.write("#   Basins are exactly those already selected by the tracked "
                "<dataset>_samples.csv\n#   files, so L12 and the per-pixel tier "
                "describe the same landscapes.\n")
        f.write("#\n")
        f.write("# FIELD NAMES confirmed against the live WFS DescribeFeatureType on "
                "crn_int_basins\n"
                "# (2026-07-26); the published field dictionary is incomplete, the "
                "WFS schema is\n"
                "# authoritative. EBEL/EBELNS are in mm/kyr (same units as "
                "EBE_MMKYR); _INT/_EXT\n"
                "# are Balco v3's internal/external uncertainty pair; NS = NO "
                "TOPOGRAPHIC SHIELDING\n"
                "# (EBELNS >= EBEL always; EBELNS/EBEL tracks 1/BETOPO to mean 4e-4, "
                "max 4.5e-3).\n"
                "# CORRECTED 2026-07-26: an earlier revision said 'no SNOW shielding'. "
                "BESNOW is\n"
                "# identically 1.0000 on all 50 basins here, so it explains nothing; "
                "the rows where\n"
                "# EBELNS == EBEL are those with BETOPO >= 0.996, below EBEL's "
                "3-sig-fig precision.\n")
        f.write("# EBE_GCMYR/EBE_MMKYR are carried so the mm/kyr -> g cm^-2 yr^-1 "
                "conversion uses\n"
                "# OCTOPUS's OWN density (their ratio gives 2.6474-2.6521, i.e. "
                "2.650 g/cm3) rather\n"
                "# than an assumed one. ELEV_AVE and AREA are context only and never "
                "enter the\n"
                "# inversion. BESTND is recorded because BE10NC is already "
                "renormalised (plan §9b).\n")
        f.write("#\n")
        f.write("# CENTR_LON IS DERIVED BY US, NOT A PUBLISHED OCTOPUS FIELD. "
                "OCTOPUS publishes\n"
                "# CENTR_LAT with no centroid longitude, and lsdn_scaling() needs a "
                "longitude because\n"
                "# Holocene Rc is the great-circle distance to a MOVING paleopole. "
                "CENTR_LON is the\n"
                "# area-weighted shoelace centroid longitude of OCTOPUS's OWN basin "
                "polygon from\n"
                "# crn_int_basins (exterior rings positive, interior rings "
                "subtracted, MultiPolygon\n"
                "# parts area-weighted). CENTR_LAT_DERIVED is the latitude from the "
                "same computation.\n")
        f.write("#\n")
        f.write("# RECORDED DISCREPANCY -- CENTR_LAT_DERIVED does NOT reproduce the "
                "published\n"
                "# CENTR_LAT, and nothing in this repository asserts that it does. "
                "Over these\n"
                "# basins CENTR_LAT_DERIVED - CENTR_LAT reaches 0.00863 deg with an "
                "RMS of\n"
                "# 0.00200 deg, and on S055WTS044 the published CENTR_LAT (34.3838) "
                "lies OUTSIDE\n"
                "# the polygon's own latitude range [34.37989, 34.38226] -- which a "
                "geometric\n"
                "# centroid cannot do. The geometry is nevertheless demonstrably "
                "correct: the same\n"
                "# rings reproduce the published AREA to 9.2e-5 relative on "
                "S014WTS036 (22.3779 km2\n"
                "# computed vs 22.38 published), and to better than 8e-4 on every "
                "basin large enough\n"
                "# that AREA's own 0.01 km2 quantisation does not dominate. The "
                "LEADING HYPOTHESIS,\n"
                "# UNCONFIRMED, is that CENTR_LAT is an EFFECTIVE latitude iterated "
                "against the\n"
                "# spatially averaged Stone-Lal production in the same way ATM_PRESS "
                "is, rather than\n"
                "# a geometric centroid; a separate research task is running to "
                "settle it. Until it\n"
                "# does, CENTR_LAT_DERIVED is a frozen RECORD of the discrepancy, "
                "not a check that\n"
                "# passed. L12 evaluates at the PUBLISHED CENTR_LAT and the derived "
                "CENTR_LON.\n")
        f.write("#\n")
        f.write("# OUTLET_LON is the basin OUTLET longitude from crn_int_outlets, "
                "kept as a DIAGNOSTIC:\n"
                "# it was the longitude L12 used before CENTR_LON existed, and "
                "CENTR_LON - OUTLET_LON\n"
                "# (max 0.0329 deg, RMS 0.0129 deg here) is the size of that former "
                "approximation.\n")
        w = csv.DictWriter(f, fieldnames=POINT_REFERENCE_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nFROZE artefact (--freeze): {out_path} ({len(rows)} basins)")
    print("  This file is tracked. Review the diff before committing: a numeric "
          "change here is a change in the ORACLE, not bookkeeping.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", help="OCTOPUS STUDYID, e.g. S014")
    ap.add_argument("--dataset",
                    help="dataset name for the tracked files, e.g. clapp2002")
    ap.add_argument("--out", default="tests/expected/cosmo/published",
                    help="output directory for the tracked artefacts")
    ap.add_argument("--scaling", choices=sorted(SCALING_FIELDS), default="St",
                    help="which OCTOPUS rate fields populate published_erate_*: "
                         "St = EBE_* (CAIRN Stone/Lal, the default so the tracked "
                         "St artefacts never move), LSDn = EBEL, LSDnNS = EBELNS "
                         "(Balco v3 point results, with / without topographic "
                         "shielding)")
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

    # --- second mode: the DEM-free LSDn point reference (L12 kernel oracle) ---
    ap.add_argument("--point-reference", action="store_true",
                    help="harvest CENTR_LAT/ATM_PRESS/EBEL*/EBELNS* for the basins "
                         "already selected by the tracked <dataset>_samples.csv "
                         "files, into the frozen L12 kernel oracle")
    ap.add_argument("--datasets", default="dethier2014,clapp2002,dibiase2010",
                    help="comma-separated tracked datasets to harvest "
                         "(--point-reference only)")
    ap.add_argument("--published-dir", default="tests/expected/cosmo/published",
                    help="directory holding <dataset>_samples.csv + _manifest.json")
    ap.add_argument("--point-reference-out",
                    default="tests/expected/cosmo/octopus_lsdn_point_reference.csv",
                    help="the frozen point-reference artefact")
    ap.add_argument("--freeze", action="store_true",
                    help="REGENERATE the tracked frozen artefact instead of "
                         "checking against it. Deliberate, reviewed act: the "
                         "artefact is tracked and a numeric change in it is a "
                         "change in the oracle. Default is check-only.")
    args = ap.parse_args()

    if args.point_reference:
        try:
            return harvest_point_reference(args)
        except Exception as exc:  # noqa: BLE001
            print(f"POINT-REFERENCE HARVEST FAILED: {exc}", file=sys.stderr)
            return 1
    if not args.study or not args.dataset:
        ap.error("--study and --dataset are required (unless --point-reference)")
    try:
        return harvest(args)
    except Exception as exc:  # noqa: BLE001
        print(f"HARVEST FAILED: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
