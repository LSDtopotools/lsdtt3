#!/usr/bin/env python3
"""
download_dem_tiles.py  --  Download DEM data from three sources:

  1. Microsoft Planetary Computer (default)
       No API key required. Downloads full tiles intersecting a point or bbox.
       Useful for global datasets such as Copernicus GLO-30.

  2. OpenTopography REST API  (--ot-dataset)
       Requires a free API key from https://opentopography.org
       The key is resolved in this order (first found wins):
         a) Environment variable:  export OT_API_KEY="your_key"
         b) Key file:              --ot-api-key-file /path/to/keyfile.txt
                                   (file must contain the key as the first non-blank line;
                                    set permissions with: chmod 600 keyfile.txt
                                    and never commit the file to a repository)
         c) Interactive prompt:    the script will ask if neither of the above is set

  3. WCS endpoint  (--wcs-url)
        For services that require an API key the same resolution order applies,
        using the environment variable WCS_API_KEY or --api-key-file.
        Many public WCS services (e.g. national mapping agencies) need no key at all.

  4. Copernicus Data Space Ecosystem EU-DEM  (--cdse-eu-dem)
        Uses the CDSE STAC endpoint (collection default: cop-dem-eea-10-laea-tif).
        Downloads either the authenticated zipped product asset (default) or
        another chosen asset from matching STAC items.

        Credentials are resolved in this order (first found wins):
          a) Command line:          --cdse-username / --cdse-password
          b) Environment variables: CDSE_USERNAME / CDSE_PASSWORD
          c) Interactive prompt (only if needed)

Credential examples (bash):
  export OT_API_KEY="your_ot_key"
  export WCS_API_KEY="your_wcs_key"
  export CDSE_USERNAME="your_cdse_username"
  export CDSE_PASSWORD="your_cdse_password"
"""

from pathlib import Path
from urllib.parse import parse_qsl, urlencode, urljoin, urlparse, urlunparse
import argparse
import getpass
import os
import re
import requests

try:
    import planetary_computer
except ImportError:
    planetary_computer = None

try:
    from pystac_client import Client
except ImportError:
    Client = None

CATALOG_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
CDSE_STAC_URL = "https://stac.dataspace.copernicus.eu/v1"
CDSE_TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)


def parse_latlon(text):
    """
    Parse a string like:
        "55.95244789111162, -3.2193971081922665"
    Returns (lat, lon)
    """
    parts = text.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"Invalid coordinate '{text}'. Expected format: 'lat, lon'"
        )

    try:
        lat = float(parts[0].strip())
        lon = float(parts[1].strip())
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid coordinate '{text}'. Could not parse numbers."
        )

    if not (-90 <= lat <= 90):
        raise argparse.ArgumentTypeError(f"Latitude out of range in '{text}'")
    if not (-180 <= lon <= 180):
        raise argparse.ArgumentTypeError(f"Longitude out of range in '{text}'")

    return lat, lon


def to_bbox(c1, c2):
    """
    Convert two opposite corners from (lat, lon) to STAC bbox:
    [west, south, east, north]
    """
    lat1, lon1 = c1
    lat2, lon2 = c2

    west = min(lon1, lon2)
    south = min(lat1, lat2)
    east = max(lon1, lon2)
    north = max(lat1, lat2)

    return [west, south, east, north]


def point_geojson(latlon):
    """
    Convert (lat, lon) to GeoJSON Point with coordinates [lon, lat]
    """
    lat, lon = latlon
    return {
        "type": "Point",
        "coordinates": [lon, lat],
    }


def infer_suffix_from_url(url):
    path = urlparse(url).path.lower()
    if path.endswith(".tif"):
        return ".tif"
    if path.endswith(".tiff"):
        return ".tiff"
    return ".tif"


def _follow_redirects_with_auth(session, url, headers, timeout=300, max_redirects=10):
    """GET url, preserving Authorization headers across cross-host redirects."""
    current = url
    for _ in range(max_redirects):
        r = session.get(
            current,
            stream=True,
            timeout=timeout,
            headers=headers,
            allow_redirects=False,
        )
        if r.is_redirect or r.is_permanent_redirect:
            loc = r.headers.get("location")
            if not loc:
                return r
            current = urljoin(current, loc)
            continue
        return r
    raise RuntimeError(f"Too many redirects while downloading: {url}")


def _candidate_download_urls(url):
    """Return fallback URL candidates for CDSE product downloads."""
    parsed = urlparse(url)
    host = parsed.netloc.lower()
    out = [url]
    if host == "download.dataspace.copernicus.eu":
        out.append(urlunparse((parsed.scheme, "zipper.dataspace.copernicus.eu",
                               parsed.path, parsed.params, parsed.query, parsed.fragment)))
        out.append(urlunparse((parsed.scheme, "catalogue.dataspace.copernicus.eu",
                               parsed.path, parsed.params, parsed.query, parsed.fragment)))
    return out


def download_file(url, output_path, headers=None):
    output_path = Path(output_path)
    headers = headers or {}

    errors = []
    with requests.Session() as session:
        for candidate in _candidate_download_urls(url):
            r = _follow_redirects_with_auth(
                session, candidate, headers=headers, timeout=300
            )
            with r:
                if r.status_code in (401, 403):
                    snippet = r.text[:400].replace("\n", " ") if r.text else ""
                    errors.append(f"{candidate} -> {r.status_code} {snippet}")
                    continue

                r.raise_for_status()
                total = int(r.headers.get("content-length", 0))
                downloaded = 0

                with open(output_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if not chunk:
                            continue
                        f.write(chunk)
                        downloaded += len(chunk)

                        if total:
                            percent = downloaded * 100 / total
                            print(f"\r  {output_path.name}: {percent:5.1f}% ", end="", flush=True)

                if total:
                    print()
                print(f"  saved: {output_path.resolve()}")
                return

    detail = " | ".join(errors) if errors else f"URL: {url}"
    if any("DAT-ZIP-608" in e for e in errors):
        detail += (
            "\nWARNING: CDSE returned DAT-ZIP-608 (Access forbidden). "
            "Your login is valid but the account likely lacks COP-DEM/CCM "
            "download entitlement (license/terms acceptance). "
            "Enable COP-DEM access in CDSE account licensing, or use the "
            "EU-DEM WCS mode as a fallback."
        )

    raise RuntimeError("Download failed on all candidate URLs. " + detail)


def build_output_name(prefix, index, total, suffix):
    """
    If only one file is being downloaded, save as:
        prefix.tif
    If multiple files are being downloaded, save as:
        prefix_001.tif, prefix_002.tif, ...
    """
    if total == 1:
        return f"{prefix}{suffix}"
    return f"{prefix}_{index:03d}{suffix}"


def read_api_key(key_file):
    """Return the first non-empty, non-comment line from key_file."""
    with open(key_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                return line
    raise ValueError(f"No API key found in '{key_file}'")


def get_cdse_access_token(username, password):
    payload = {
        "client_id": "cdse-public",
        "grant_type": "password",
        "username": username,
        "password": password,
    }
    r = requests.post(CDSE_TOKEN_URL, data=payload, timeout=60)
    r.raise_for_status()
    body = r.json()
    token = body.get("access_token")
    if not token:
        raise RuntimeError("CDSE token response did not include access_token")
    return token


def bbox_intersects(a, b):
    """Return True if bbox a intersects bbox b, where bboxes are [w,s,e,n]."""
    aw, as_, ae, an = a
    bw, bs, be, bn = b
    return (aw <= be) and (ae >= bw) and (as_ <= bn) and (an >= bs)


def grid_code_to_bbox(grid_code):
    """Parse CDEM grid code (e.g. CDEM-N58W005) to [w,s,e,n] bbox."""
    m = re.match(r"^CDEM-([NS])(\d{2})([EW])(\d{2,3})$", str(grid_code))
    if not m:
        return None
    ns, lat_s, ew, lon_s = m.groups()
    lat = int(lat_s)
    lon = int(lon_s)
    if ns == "S":
        lat = -lat
    if ew == "W":
        lon = -lon
    return [float(lon), float(lat), float(lon + 1), float(lat + 1)]


def filter_items_by_bbox(items, bbox):
    """Drop STAC items whose bbox does not intersect the requested bbox."""
    kept = []
    dropped = []
    for item in items:
        props = item.properties or {}
        gb = grid_code_to_bbox(props.get("grid:code"))
        if gb is not None and not bbox_intersects(gb, bbox):
            dropped.append(item)
            continue
        ib = getattr(item, "bbox", None)
        if ib and len(ib) >= 4:
            ib4 = [ib[0], ib[1], ib[2], ib[3]]
            if not bbox_intersects(ib4, bbox):
                dropped.append(item)
                continue
        kept.append(item)
    return kept, dropped


def latest_by_grid_code(items):
    """Deduplicate CDSE DEM items by grid code, keep most recently updated."""
    picked = {}
    for item in items:
        props = item.properties or {}
        key = props.get("grid:code") or item.id
        stamp = props.get("updated") or props.get("published") or ""
        cur = picked.get(key)
        if cur is None:
            picked[key] = item
            continue
        cur_stamp = (
            (cur.properties or {}).get("updated")
            or (cur.properties or {}).get("published")
            or ""
        )
        if stamp > cur_stamp:
            picked[key] = item
    out = list(picked.values())
    out.sort(key=lambda it: ((it.properties or {}).get("grid:code", ""), it.id))
    return out


def append_key_to_url(url, api_key, param_name="key"):
    """Append api_key as a query parameter to url."""
    parsed = urlparse(url)
    query = dict(parse_qsl(parsed.query, keep_blank_values=True))
    query[param_name] = api_key
    return urlunparse(
        (parsed.scheme, parsed.netloc, parsed.path, parsed.params, urlencode(query), parsed.fragment)
    )


def append_query_params(url, params):
    """Append/overwrite query parameters on a URL."""
    parsed = urlparse(url)
    query = dict(parse_qsl(parsed.query, keep_blank_values=True))
    query.update(params)
    return urlunparse(
        (parsed.scheme, parsed.netloc, parsed.path, parsed.params, urlencode(query), parsed.fragment)
    )


def wcs_preflight_getcapabilities(url, version):
    """Preflight WCS endpoint before OWSLib calls.

    Raises RuntimeError with a clear message when the endpoint is not serving
    WCS capabilities at the requested URL/version.
    """
    caps_url = append_query_params(
        url,
        {
            "SERVICE": "WCS",
            "REQUEST": "GetCapabilities",
            "VERSION": version,
        },
    )

    r = requests.get(caps_url, timeout=60)
    if r.status_code != 200:
        snippet = (r.text or "").strip().replace("\n", " ")[:500]
        raise RuntimeError(
            f"WCS preflight failed (HTTP {r.status_code}) for GetCapabilities. "
            f"URL: {caps_url}. Response snippet: {snippet}"
        )

    body = (r.text or "")
    if (
        "WCS_Capabilities" not in body
        and "wcs:Capabilities" not in body
        and "http://www.opengis.net/wcs" not in body
    ):
        snippet = body.strip().replace("\n", " ")[:500]
        raise RuntimeError(
            "WCS preflight failed: GetCapabilities did not return a WCS "
            f"capabilities document. URL: {caps_url}. Response snippet: {snippet}"
        )


def download_wcs_coverage(url, coverage, bbox, out_file, crs="EPSG:4326",
                          version="1.0.0", fmt="image/tiff", res=None):
    """Download a single coverage from a WCS endpoint.

    bbox must be [west, south, east, north] in the given crs.
    """
    from owslib.wcs import WebCoverageService

    west, south, east, north = bbox

    # Preflight first so endpoint/version problems fail with a clear message.
    wcs_preflight_getcapabilities(url, version)

    print(f"Connecting to WCS: {url}")
    wcs = WebCoverageService(url, version=version)

    kwargs = {
        "identifier": coverage,
        "bbox": (west, south, east, north),
        "crs": crs,
        "format": fmt,
    }
    if res is not None:
        kwargs["resx"] = res
        kwargs["resy"] = res

    print(f"Coverage: {coverage}")
    print(f"BBox:     {bbox}")
    print(f"CRS:      {crs}")
    if res is not None:
        print(f"Res:      {res}")
    print()

    response = wcs.getCoverage(**kwargs)
    data = response.read()

    if not data:
        raise RuntimeError("Server returned an empty response.")

    probe = data[:200].lstrip().lower()
    if probe.startswith(b"<"):
        snippet = data[:1000].decode(errors="replace")
        raise RuntimeError(
            f"Server returned XML/text instead of raster data.\nSnippet:\n{snippet}"
        )

    out_path = Path(out_file)
    with open(out_path, "wb") as f:
        f.write(data)

    print(f"Saved: {out_path.resolve()}  ({len(data)} bytes)")


OT_USGS_DATASETS = {"USGS30m", "USGS10m", "USGS1m"}
OT_GLOBAL_URL  = "https://portal.opentopography.org/API/globaldem"
OT_USGS_URL    = "https://portal.opentopography.org/API/usgsdem"


def download_opentopography(dataset, bbox, api_key, out_file):
    """Download a DEM from the OpenTopography REST API.

    dataset  -- OT demtype/datasetName string, e.g. 'SRTMGL1', 'COP30', 'USGS10m'
    bbox     -- [west, south, east, north] in WGS84 decimal degrees
    api_key  -- OT API key string
    out_file -- output file path
    """
    west, south, east, north = bbox

    if dataset in OT_USGS_DATASETS:
        base = OT_USGS_URL
        params = {
            "datasetName": dataset,
            "south": south, "north": north,
            "west": west,   "east": east,
            "outputFormat": "GTiff",
            "API_Key": api_key,
        }
    else:
        base = OT_GLOBAL_URL
        params = {
            "demtype": dataset,
            "south": south, "north": north,
            "west": west,   "east": east,
            "outputFormat": "GTiff",
            "API_Key": api_key,
        }

    url = f"{base}?{urlencode(params)}"
    safe_url = url.replace(api_key, "***")

    print(f"Dataset:  {dataset}")
    print(f"BBox:     {bbox}")
    print(f"URL (key hidden): {safe_url}")
    print()

    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()

        # OT returns an error as XML/JSON on a 200 if the bbox is bad, check it.
        content_type = r.headers.get("content-type", "")
        if "xml" in content_type or "json" in content_type or "text" in content_type:
            body = r.text[:1000]
            raise RuntimeError(
                f"OT returned text instead of raster data ({content_type}):\n{body}"
            )

        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        out_path = Path(out_file)

        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    print(f"\r  {out_path.name}: {downloaded * 100 / total:5.1f}% ", end="", flush=True)

    if total:
        print()
    print(f"Saved: {out_path.resolve()}  ({downloaded} bytes)")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download DEM tiles from Microsoft Planetary Computer.\n\n"
            "Pass either:\n"
            "  - one coordinate: download the full tile(s) containing that point\n"
            "  - two coordinates: download all full tiles intersecting the bbox "
            "defined by those two corners\n\n"
            "Coordinates must be quoted and given as 'lat, lon'."
        ),
        epilog=(
            "Help and credential setup\n"
            "--------------------------\n"
            "This script supports four modes:\n"
            "  1) Planetary Computer (default, no credentials)\n"
            "  2) OpenTopography (--ot-dataset, needs API key)\n"
            "  3) WCS (--wcs-url, may need API key)\n"
            "  4) CDSE EU-DEM (--cdse-eu-dem, product downloads need login)\n\n"
            "How to add credentials with export (bash):\n"
            "  export OT_API_KEY=\"your_ot_key\"\n"
            "  export WCS_API_KEY=\"your_wcs_key\"\n"
            "  export CDSE_USERNAME=\"your_cdse_username\"\n"
            "  export CDSE_PASSWORD=\"your_cdse_password\"\n\n"
            "Run examples:\n"
            "  # EU-DEM 10 m (CDSE) by bbox\n"
            "  python download_dem_tiles.py \"57.0, -4.5\" \"58.0, -3.0\" --cdse-eu-dem -o eudem\n"
            "\n"
            "  # OpenTopography COP30 by bbox\n"
            "  python download_dem_tiles.py \"57.0, -4.5\" \"58.0, -3.0\" --ot-dataset COP30 -o cop30\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "coords",
        nargs="+",
        type=parse_latlon,
        help=(
            "One or two coordinates as 'lat, lon'\n"
            "Examples:\n"
            "  \"55.95, -3.21\"\n"
            "  \"55.95, -3.21\" \"55.96, -3.20\""
        ),
    )
    parser.add_argument(
        "-o", "--output-prefix",
        default="my_DEM",
        help="Prefix for output files (default: my_DEM)",
    )
    parser.add_argument(
        "-d", "--output-dir",
        default=".",
        help="Directory to save files into (default: current directory)",
    )
    parser.add_argument(
        "-c", "--collection",
        default="cop-dem-glo-30",
        help="Planetary Computer collection (default: cop-dem-glo-30)",
    )
    parser.add_argument(
        "-a", "--asset",
        default="data",
        help="Asset key to download from each STAC item (default: data)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=200,
        help="Maximum number of matching items to fetch (default: 200)",
    )

    wcs_group = parser.add_argument_group(
        "WCS options",
        "Use these to download from a WCS endpoint instead of Planetary Computer.",
    )
    wcs_group.add_argument(
        "--wcs-url",
        help="WCS endpoint URL.  If provided, WCS mode is used.",
    )
    wcs_group.add_argument(
        "--wcs-coverage",
        help="Coverage ID to request from the WCS service.",
    )
    wcs_group.add_argument(
        "--wcs-version",
        default="1.0.0",
        help="WCS protocol version (default: 1.0.0)",
    )
    wcs_group.add_argument(
        "--wcs-crs",
        default="EPSG:4326",
        help="CRS to pass in the WCS GetCoverage request (default: EPSG:4326)",
    )
    wcs_group.add_argument(
        "--wcs-res",
        type=float,
        default=None,
        help="Cell size in CRS units (e.g. 0.0001 for ~10 m in EPSG:4326)",
    )
    wcs_group.add_argument(
        "--wcs-format",
        default="image/tiff",
        help="Output format requested from WCS (default: image/tiff)",
    )
    wcs_group.add_argument(
        "--api-key-file",
        help="Path to a file whose first non-blank line is the API key.  "
             "The key is appended to the WCS URL as a query parameter.",
    )
    wcs_group.add_argument(
        "--api-key-param",
        default="key",
        help="Query-parameter name used to pass the API key (default: key)",
    )

    ot_group = parser.add_argument_group(
        "OpenTopography options",
        "Use these to download from the OpenTopography REST API instead of Planetary Computer.",
    )
    ot_group.add_argument(
        "--ot-dataset",
        help=(
            "OpenTopography dataset name.  If provided, OT mode is used.\n"
            "Global datasets: SRTMGL1 (30m), SRTMGL3 (90m), AW3D30, COP30, COP90,\n"
            "                 SRTM15Plus, GEBCOIceTopo, GEBCOSubIceTopo\n"
            "USGS datasets:   USGS30m, USGS10m, USGS1m"
        ),
    )
    ot_group.add_argument(
        "--ot-api-key-file",
        help="Path to a file whose first non-blank line is your OpenTopography API key.",
    )

    cdse_group = parser.add_argument_group(
        "CDSE EU-DEM options",
        "Use these to download from Copernicus Data Space Ecosystem STAC.",
    )
    cdse_group.add_argument(
        "--cdse-eu-dem",
        action="store_true",
        help="Enable CDSE mode (defaults to EU-DEM 10 m collection).",
    )
    cdse_group.add_argument(
        "--cdse-collection",
        default="cop-dem-eea-10-laea-tif",
        help=(
            "CDSE STAC collection (default: cop-dem-eea-10-laea-tif). "
            "Examples: cop-dem-glo-30-dged-cog, cop-dem-glo-90-dged-cog"
        ),
    )
    cdse_group.add_argument(
        "--cdse-asset",
        default="product",
        help=(
            "Asset key to download in CDSE mode (default: product). "
            "Use 'product' for zipped download endpoint, or 'data' for tile raster asset."
        ),
    )
    cdse_group.add_argument(
        "--cdse-username",
        default=None,
        help="CDSE username (or set CDSE_USERNAME).",
    )
    cdse_group.add_argument(
        "--cdse-password",
        default=None,
        help="CDSE password (or set CDSE_PASSWORD).",
    )

    args = parser.parse_args()

    if len(args.coords) not in (1, 2):
        parser.error("You must provide either 1 coordinate or 2 coordinates.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ WCS --
    if args.wcs_url:
        if not args.wcs_coverage:
            parser.error("--wcs-url requires --wcs-coverage")

        if len(args.coords) != 2:
            parser.error("WCS mode requires exactly 2 coordinates to define a bounding box.")

        bbox = to_bbox(args.coords[0], args.coords[1])

        wcs_url = args.wcs_url
        wcs_api_key = (
            os.environ.get("WCS_API_KEY")
            or (read_api_key(args.api_key_file) if args.api_key_file else None)
        )
        if wcs_api_key:
            wcs_url = append_key_to_url(
                wcs_url, wcs_api_key, param_name=args.api_key_param
            )
        else:
            print(
                "No WCS API key provided (WCS_API_KEY / --api-key-file). "
                "Proceeding without key."
            )

        suffix = ".tif"
        out_file = out_dir / f"{args.output_prefix}{suffix}"

        try:
            download_wcs_coverage(
                url=wcs_url,
                coverage=args.wcs_coverage,
                bbox=bbox,
                out_file=out_file,
                crs=args.wcs_crs,
                version=args.wcs_version,
                fmt=args.wcs_format,
                res=args.wcs_res,
            )
        except Exception as e:
            raise SystemExit(f"WCS download failed: {e}")

        print("\nDone.")
        return
    # --------------------------------------------------------- end WCS block --

    # ---------------------------------------------------------- OT block ----
    if args.ot_dataset:
        if len(args.coords) != 2:
            parser.error("OpenTopography mode requires exactly 2 coordinates to define a bounding box.")

        bbox = to_bbox(args.coords[0], args.coords[1])
        ot_api_key = (
            os.environ.get("OT_API_KEY")
            or (read_api_key(args.ot_api_key_file) if args.ot_api_key_file else None)
            or getpass.getpass("OpenTopography API key: ")
        )
        out_file = out_dir / f"{args.output_prefix}_{args.ot_dataset}.tif"

        try:
            download_opentopography(
                dataset=args.ot_dataset,
                bbox=bbox,
                api_key=ot_api_key,
                out_file=out_file,
            )
        except Exception as e:
            # Requests can include the full query URL in connection errors;
            # never allow the API key to reach a pipeline log.
            detail = str(e).replace(ot_api_key, "***")
            raise SystemExit(f"OpenTopography download failed: {detail}")

        print("\nDone.")
        return
    # ------------------------------------------------------- end OT block --

    # ------------------------------------------------------- CDSE block ----
    if args.cdse_eu_dem:
        if Client is None:
            raise SystemExit(
                "CDSE mode requires pystac-client; install the Python helper dependencies."
            )
        catalog = Client.open(CDSE_STAC_URL)

        print(f"CDSE collection: {args.cdse_collection}")
        print(f"Asset key:        {args.cdse_asset}")

        req_bbox = None
        if len(args.coords) == 1:
            lat, lon = args.coords[0]
            print("Mode:             point search")
            print(f"Point:            lat={lat}, lon={lon}")
            print()

            search = catalog.search(
                collections=[args.cdse_collection],
                intersects=point_geojson(args.coords[0]),
                limit=args.limit,
            )
        else:
            bbox = to_bbox(args.coords[0], args.coords[1])
            req_bbox = bbox
            print("Mode:             bbox search")
            print(f"BBOX:             {bbox}")
            print()

            search = catalog.search(
                collections=[args.cdse_collection],
                bbox=bbox,
                limit=args.limit,
            )

        raw_items = list(search.items())
        if req_bbox is not None:
            raw_items, dropped = filter_items_by_bbox(raw_items, req_bbox)
            if dropped:
                print(
                    f"Dropped {len(dropped)} item(s) with non-intersecting bboxes "
                    "(catalog anomaly)."
                )
        if not raw_items:
            raise SystemExit("No matching CDSE DEM tile(s) found.")

        items = latest_by_grid_code(raw_items)
        if not items:
            raise SystemExit("No usable CDSE DEM tile(s) after deduplication.")

        print(
            f"Found {len(raw_items)} item(s), using {len(items)} latest unique grid tile(s):"
        )
        for i, item in enumerate(items, start=1):
            grid = (item.properties or {}).get("grid:code", "")
            stamp = (item.properties or {}).get("updated", "")
            label = f"{item.id} [{grid}]" if grid else item.id
            if stamp:
                label += f" updated={stamp}"
            print(f"  {i:03d}: {label}")
        print()

        need_auth = args.cdse_asset == "product"
        cdse_headers = {}
        if need_auth:
            username = args.cdse_username or os.environ.get("CDSE_USERNAME")
            password = args.cdse_password or os.environ.get("CDSE_PASSWORD")
            if not username:
                username = input("CDSE username: ").strip()
            if not password:
                password = getpass.getpass("CDSE password: ")
            token = get_cdse_access_token(username, password)
            cdse_headers = {"Authorization": f"Bearer {token}"}

        total = len(items)
        for i, item in enumerate(items, start=1):
            if args.cdse_asset not in item.assets:
                print(f"Skipping {item.id}: missing asset '{args.cdse_asset}'")
                print(f"Available assets: {list(item.assets.keys())}")
                continue

            url = item.assets[args.cdse_asset].href
            if args.cdse_asset == "product":
                suffix = ".zip"
            else:
                suffix = infer_suffix_from_url(url)

            filename = build_output_name(args.output_prefix, i, total, suffix)
            output_path = out_dir / filename

            print(f"Downloading {item.id} -> {filename}")
            download_file(url, output_path, headers=cdse_headers)

        print()
        print("Done.")
        return
    # --------------------------------------------------- end CDSE block ----

    if Client is None or planetary_computer is None:
        raise SystemExit(
            "Planetary Computer mode requires pystac-client and planetary-computer; "
            "install the Python helper dependencies."
        )
    catalog = Client.open(CATALOG_URL)

    print(f"Collection: {args.collection}")
    print(f"Asset key:  {args.asset}")

    if len(args.coords) == 1:
        lat, lon = args.coords[0]
        print(f"Mode:       point search")
        print(f"Point:      lat={lat}, lon={lon}")
        print()

        search = catalog.search(
            collections=[args.collection],
            intersects=point_geojson(args.coords[0]),
            limit=args.limit,
        )

    else:
        bbox = to_bbox(args.coords[0], args.coords[1])
        print(f"Mode:       bbox search")
        print(f"BBOX:       {bbox}")
        print()

        search = catalog.search(
            collections=[args.collection],
            bbox=bbox,
            limit=args.limit,
        )

    items = list(search.items())
    if not items:
        raise SystemExit("No matching tile(s) found.")

    items.sort(key=lambda item: item.id)

    print(f"Found {len(items)} tile(s):")
    for i, item in enumerate(items, start=1):
        print(f"  {i:03d}: {item.id}")
    print()

    total = len(items)

    for i, raw_item in enumerate(items, start=1):
        item = planetary_computer.sign(raw_item)

        if args.asset not in item.assets:
            print(f"Skipping {item.id}: missing asset '{args.asset}'")
            print(f"Available assets: {list(item.assets.keys())}")
            continue

        url = item.assets[args.asset].href
        suffix = infer_suffix_from_url(url)
        filename = build_output_name(args.output_prefix, i, total, suffix)
        output_path = out_dir / filename

        print(f"Downloading {item.id} -> {filename}")
        download_file(url, output_path)

    print()
    print("Done.")


if __name__ == "__main__":
    main()
