#!/usr/bin/env python3
"""Fetch DEM fixtures for lsdtt3 tests.

A general fetcher shared by every test suite that needs real topography. No DEM
is ever committed to the repository; fixtures are cut on demand from public
keyless sources and cached under tests/data/dem_cache/.

Two ways to use it:

  1. Named fixtures from the REGISTRY below -- deterministic, reproducible, and
     SHA256-pinnable. Use these for anything whose result is asserted.

  2. Ad-hoc bounding boxes (--bbox) -- for exploratory work and for the
     download-run-delete validation tier, where the DEM is transient and only a
     small summary CSV is kept (see docs/cosmo_cli_testing_plan.md Sec5b).
     Combine with --dest and --no-cache to write outside the cache.

Source: Copernicus GLO-30 on the AWS Open Data registry. Public, keyless, stable
URL scheme, stored as Cloud-Optimised GeoTIFFs, so a sub-window arrives via HTTP
range requests rather than a full 1-degree tile download. Bounding boxes that
span several tiles are mosaicked through an in-memory VRT.

Fixtures are reprojected to a metric CRS here, at fetch time, rather than by the
CLI under test. That matters whenever two implementations are compared: both
must read a byte-identical projected raster, or the differences are
uninterpretable. UTM zone is derived from the bounding box centre unless an EPSG
code is given.

Usage:
    pixi run fetch-test-dem --list
    pixi run fetch-test-dem --suite flow
    pixi run fetch-test-dem --fixture alpine --size medium
    pixi run fetch-test-dem --bbox 10.3 46.4 10.5 46.6 --dest /tmp/work/dem.tif --no-cache
"""

import argparse
import hashlib
import json
import math
import pathlib
import sys

try:
    from osgeo import gdal
except ImportError:  # pragma: no cover
    sys.exit("FATAL: GDAL Python bindings not found. Run under pixi, e.g. "
             "`pixi run fetch-test-dem --list`")

gdal.UseExceptions()

# /vsicurl/ tuning. GDAL_DISABLE_READDIR_ON_OPEN suppresses a bucket listing on
# every open (a large latency win against S3); the timeouts turn an offline
# machine into a prompt, bounded failure instead of a hung test.
gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
gdal.SetConfigOption("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
gdal.SetConfigOption("GDAL_HTTP_TIMEOUT", "60")
gdal.SetConfigOption("GDAL_HTTP_CONNECTTIMEOUT", "15")
gdal.SetConfigOption("GDAL_HTTP_MAX_RETRY", "3")
gdal.SetConfigOption("GDAL_HTTP_RETRY_DELAY", "2")

# Exit code 77 is CTest's "skip". A network fetch that cannot reach the source
# is an unavailable dependency, not a test failure.
SKIP_EXIT = 77


class FetchUnavailable(RuntimeError):
    """The source could not be reached -- callers should skip, not fail."""

COP30_URL = (
    "/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_{tile}_DEM/Copernicus_DSM_COG_10_{tile}_DEM.tif"
)

DEFAULT_RESOLUTION_M = 30.0

# NoData sentinel written into every fixture. Matches lsdtt_raster_meta's
# default so the C++ side needs no special handling. An explicit NoData value is
# part of the fixture contract, not a detail: comparisons are meaningless unless
# both sides agree which cells are valid, and some readers (RichDEM's LoadGDAL)
# refuse a raster that has none.
NODATA_VALUE = -9999.0

SIZES = {"small": 256, "medium": 720, "large": 2400}


# -----------------------------------------------------------------------------
# Fixture registry
#
# Add an entry here rather than writing a new fetch script. `suite` groups
# fixtures so a test suite can fetch only what it needs.
#
# A fixture is defined either by centre + pixel size (square, resolution-driven)
# or by an explicit WGS84 bbox. `epsg = None` means "derive the UTM zone".
# -----------------------------------------------------------------------------
REGISTRY = {
    "alpine": {
        "suite": "flow",
        "centre": (10.40, 46.50),
        "epsg": 32632,          # UTM 32N
        "description": "high relief, Italian Alps (Val Venosta / Ortler)",
    },
    "lakeland": {
        "suite": "flow",
        "centre": (25.40, 61.50),
        "epsg": 32635,          # UTM 35N
        "description": "low relief, Finnish Lakeland (many closed basins)",
    },
}

# SHA256 of each cached fixture, keyed "<fixture>_<size>". A mismatch means the
# upstream tile or the fetch parameters changed, which invalidates every stored
# expectation downstream -- so it is a hard failure, never a silent
# re-download. Entries are absent until pinned; run --print-hashes to fill in.
#
# SCOPE OF THE PIN: fetching is deterministic for a given GDAL build (verified:
# repeated --force runs reproduce the hash exactly), but warp output is
# GDAL-version sensitive -- a pixi re-solve that moved GDAL changed the alpine
# fixture by up to 2.6 m per cell while leaving lakeland byte-identical. So these
# hashes are only meaningful together with the committed pixi.lock, which pins
# the GDAL build. If you regenerate pixi.lock, expect to re-pin, and record the
# GDAL version alongside the hashes when you do.
PINNED_SHA256 = {}


def cache_dir() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parent.parent / "data" / "dem_cache"


def sha256_of(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def utm_epsg_for(lon: float, lat: float) -> int:
    """EPSG code of the UTM zone containing (lon, lat)."""
    zone = int((lon + 180.0) / 6.0) + 1
    return (32600 if lat >= 0 else 32700) + zone


def tiles_for_bbox(bbox):
    """Copernicus GLO-30 tile names covering a WGS84 bbox (west, south, east, north).

    Tiles are 1 degree squares named by their south-west corner, so the range is
    taken over floor(south)..floor(north) and likewise in longitude.
    """
    west, south, east, north = bbox
    names = []
    for lat in range(math.floor(south), math.floor(north) + 1):
        for lon in range(math.floor(west), math.floor(east) + 1):
            ns = "N" if lat >= 0 else "S"
            ew = "E" if lon >= 0 else "W"
            names.append(f"{ns}{abs(lat):02d}_00_{ew}{abs(lon):03d}_00")
    return names


def source_dataset(bbox):
    """Open the source as a single dataset, mosaicking tiles when the bbox spans several."""
    tiles = tiles_for_bbox(bbox)
    urls = [COP30_URL.format(tile=t) for t in tiles]

    try:
        if len(urls) == 1:
            return gdal.Open(urls[0]), tiles

        # A bbox crossing a tile boundary is normal for real study areas (the
        # cosmo validation basins in particular), so mosaic rather than failing.
        return gdal.BuildVRT("/vsimem/mosaic.vrt", urls), tiles
    except RuntimeError as exc:
        # Distinguish "no network" from "wrong tile name": the first is a skip,
        # the second is a genuine error worth surfacing.
        msg = str(exc)
        if "CURL" in msg or "resolve host" in msg or "timed out" in msg:
            raise FetchUnavailable(
                f"could not reach the DEM source ({msg.strip()}). "
                "Check network access; tests depending on this should skip."
            ) from exc
        raise


def fetch_bbox(bbox, out_path: pathlib.Path, epsg=None,
               resolution=DEFAULT_RESOLUTION_M, size_px=None):
    """Cut, reproject and write one DEM. Returns a dict of what was produced.

    bbox       (west, south, east, north) in WGS84 degrees
    epsg       target EPSG; None derives the UTM zone from the bbox centre
    size_px    if set, crop the exact square of this many pixels from the centre
    """
    west, south, east, north = bbox
    if epsg is None:
        epsg = utm_epsg_for(0.5 * (west + east), 0.5 * (south + north))

    src, tiles = source_dataset(bbox)
    if src is None:
        raise RuntimeError(f"could not open source tiles {tiles}")

    # Step 1: windowed read over /vsicurl (HTTP range requests).
    cut = gdal.Translate("/vsimem/cut.tif", src, projWin=[west, north, east, south])
    src = None

    # Step 2: reproject to a metric grid. Cubic matches the CLI default
    # resampling_method. dstNodata matters: without it, cells outside the source
    # footprint default to 0, which reads as sea level and becomes a giant
    # spurious depression.
    warped = gdal.Warp(
        "/vsimem/warp.tif", cut,
        dstSRS=f"EPSG:{epsg}",
        xRes=resolution, yRes=resolution,
        resampleAlg="cubic",
        targetAlignedPixels=True,
        dstNodata=NODATA_VALUE,
    )
    cut = None

    out_path.parent.mkdir(parents=True, exist_ok=True)
    translate_kwargs = dict(noData=NODATA_VALUE,
                            creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])

    if size_px is not None:
        # Crop an exact square from the centre so fixture dimensions are
        # deterministic rather than a function of GDAL's warp extent arithmetic.
        w, h = warped.RasterXSize, warped.RasterYSize
        if w < size_px or h < size_px:
            raise RuntimeError(
                f"warped raster {w}x{h} is smaller than the requested "
                f"{size_px}x{size_px}; widen the bbox")
        translate_kwargs["srcWin"] = [(w - size_px) // 2, (h - size_px) // 2,
                                      size_px, size_px]

    gdal.Translate(str(out_path), warped, **translate_kwargs)
    warped = None
    for tmp in ("/vsimem/cut.tif", "/vsimem/warp.tif", "/vsimem/mosaic.vrt"):
        try:
            gdal.Unlink(tmp)
        except RuntimeError:
            pass

    return {"path": str(out_path), "epsg": epsg, "tiles": tiles}


def bbox_for_fixture(spec, size_px, resolution):
    """WGS84 bbox for a registry entry, wide enough to survive reprojection."""
    if "bbox" in spec:
        return spec["bbox"]

    lon, lat = spec["centre"]
    # Generous margin so the post-warp centre crop never runs off the edge onto
    # nodata.
    half_m = 0.5 * size_px * resolution * 1.25
    dlat = half_m / 111_320.0
    dlon = dlat / max(math.cos(math.radians(lat)), 1e-6)
    return (lon - dlon, lat - dlat, lon + dlon, lat + dlat)


def ensure_fixture(name, size_name, force=False, dest=None, use_cache=True):
    spec = REGISTRY[name]
    size_px = SIZES[size_name]
    key = f"{name}_{size_name}"

    out_path = pathlib.Path(dest) if dest else cache_dir() / f"{key}.tif"

    if out_path.exists() and use_cache and not force:
        print(f"  cached: {out_path}")
    else:
        bbox = bbox_for_fixture(spec, size_px, DEFAULT_RESOLUTION_M)
        print(f"  cutting {key} from bbox "
              f"({bbox[0]:.3f},{bbox[1]:.3f}) - ({bbox[2]:.3f},{bbox[3]:.3f})")
        info = fetch_bbox(bbox, out_path, epsg=spec.get("epsg"),
                          size_px=size_px)
        print(f"  wrote:  {out_path}  [EPSG:{info['epsg']}, "
              f"tiles: {', '.join(info['tiles'])}]")

    digest = sha256_of(out_path)
    pinned = PINNED_SHA256.get(key)
    if pinned is None:
        print(f"  sha256: {digest}  (UNPINNED -- paste into PINNED_SHA256)")
    elif pinned != digest:
        raise SystemExit(
            f"FATAL: {key} sha256 mismatch.\n"
            f"  expected {pinned}\n  got      {digest}\n"
            "The upstream tile changed, or the fetch parameters were edited.\n"
            "Every downstream expectation is void until this is resolved.\n"
            "Do not simply re-pin without understanding why it moved.")
    else:
        print(f"  sha256: {digest}  (pinned, OK)")

    return out_path, digest


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixture", help="registry fixture name, or 'all'")
    ap.add_argument("--suite", help="fetch every fixture in this suite (e.g. flow)")
    ap.add_argument("--size", choices=list(SIZES) + ["all"], default="medium")
    ap.add_argument("--bbox", nargs=4, type=float,
                    metavar=("WEST", "SOUTH", "EAST", "NORTH"),
                    help="ad-hoc WGS84 bounding box instead of a registry fixture")
    ap.add_argument("--dest", help="explicit output path (implies a single fixture/bbox)")
    ap.add_argument("--epsg", type=int, help="target EPSG (default: derived UTM zone)")
    ap.add_argument("--resolution", type=float, default=DEFAULT_RESOLUTION_M,
                    help="target grid spacing in metres (default 30)")
    ap.add_argument("--no-cache", action="store_true",
                    help="always re-fetch; for transient download-run-delete use")
    ap.add_argument("--force", action="store_true", help="re-fetch even if cached")
    ap.add_argument("--list", action="store_true", help="list registry fixtures and exit")
    ap.add_argument("--print-hashes", action="store_true",
                    help="print a PINNED_SHA256 block for what was fetched")
    args = ap.parse_args()

    if args.list:
        print(f"{'fixture':<14}{'suite':<10}description")
        for name, spec in sorted(REGISTRY.items()):
            print(f"{name:<14}{spec['suite']:<10}{spec['description']}")
        print(f"\nsizes: {', '.join(f'{k} ({v}px)' for k, v in SIZES.items())}")
        return 0

    # Ad-hoc bbox: one raster, no registry, no pinning.
    if args.bbox:
        if not args.dest:
            return ap.error("--bbox requires --dest")
        dest = pathlib.Path(args.dest)
        print(f"fetching ad-hoc bbox {tuple(args.bbox)}")
        info = fetch_bbox(tuple(args.bbox), dest, epsg=args.epsg,
                          resolution=args.resolution)
        print(f"  wrote: {dest}  [EPSG:{info['epsg']}, "
              f"tiles: {', '.join(info['tiles'])}]")
        return 0

    if args.suite:
        names = [n for n, s in REGISTRY.items() if s["suite"] == args.suite]
        if not names:
            return ap.error(f"no fixtures in suite {args.suite!r}; try --list")
    elif args.fixture in (None, "all"):
        names = sorted(REGISTRY)
    else:
        if args.fixture not in REGISTRY:
            return ap.error(f"unknown fixture {args.fixture!r}; try --list")
        names = [args.fixture]

    sizes = list(SIZES) if args.size == "all" else [args.size]
    if args.dest and (len(names) > 1 or len(sizes) > 1):
        return ap.error("--dest requires exactly one fixture and one size")

    results = {}
    for name in sorted(names):
        print(f"{name}: {REGISTRY[name]['description']}")
        for size_name in sizes:
            _, digest = ensure_fixture(name, size_name, force=args.force,
                                       dest=args.dest,
                                       use_cache=not args.no_cache)
            results[f"{name}_{size_name}"] = digest
        print()

    if args.print_hashes:
        print("PINNED_SHA256 = {")
        for key, digest in sorted(results.items()):
            print(f'    "{key}": "{digest}",')
        print("}")

    if not args.dest:
        manifest = cache_dir() / "manifest.json"
        existing = {}
        if manifest.exists():
            existing = json.loads(manifest.read_text())
        existing.update(results)
        manifest.write_text(json.dumps(existing, indent=2, sort_keys=True) + "\n")
        print(f"wrote {manifest}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except FetchUnavailable as exc:
        print(f"SKIP: {exc}", file=sys.stderr)
        sys.exit(SKIP_EXIT)
