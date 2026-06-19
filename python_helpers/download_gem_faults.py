#!/usr/bin/env python3

import argparse
import os
import sys
import urllib.error
import urllib.request

# -----------------------------------------------------------------------------
# Download the GEM Global Active Faults Database (GEM GAF-DB).
#
# Source repository (GEM Foundation, Styron & Pagani 2020):
#   https://github.com/GEMScienceTools/gem-global-active-faults
#
# Licence: Creative Commons CC BY-SA 4.0 (attribute GEM / Styron & Pagani; any
# redistributed derivative must also be share-alike). See the repo LICENSE.
#
# This script simply DOWNLOADS the published file(s) and KEEPS GEM's own
# nomenclature unchanged. No fields are renamed, dropped, or reclassified — the
# kinematic field `slip_type` and its controlled vocabulary (Normal, Reverse,
# Dextral, Sinistral, oblique compounds such as Normal-Dextral, etc.) are left
# exactly as GEM publishes them. Map them to thrust / normal / transform later
# in your own workflow if you wish.
#
# Geometry: fault-trace LineStrings, EPSG:4326 (WGS84). Both GeoJSON and
# GeoPackage drop straight into the lsdtt3 GDAL-based vector tools.
#
# Examples:
#   # global harmonized GeoPackage into the current directory:
#   python3 download_gem_faults.py
#
#   # GeoJSON instead, to a chosen directory:
#   python3 download_gem_faults.py --format geojson -o gem_faults/
#
#   # if GEM has renamed/moved the file, point at the exact repo-relative path:
#   python3 download_gem_faults.py --path gpkg/gem_active_faults_harmonized.gpkg
#
# To clip to a region afterwards while preserving all attributes, use GDAL, e.g.
#   ogr2ogr -clipsrc <xmin> <ymin> <xmax> <ymax> faults_clip.fgb gem_active_faults_harmonized.gpkg
#
# NOTE: the exact repo paths/filenames below reflect the GEM repo layout at time
# of writing and were not live-verified in this environment. If a download
# 404s, check the repository and override with --path / --branch.
# -----------------------------------------------------------------------------

DEFAULT_REPO   = "GEMScienceTools/gem-global-active-faults"
DEFAULT_BRANCH = "master"

# Known repo-relative paths for the harmonized global dataset, keyed by format.
DATASET_PATHS = {
    "gpkg":    "gpkg/gem_active_faults_harmonized.gpkg",
    "geojson": "geojson/gem_active_faults_harmonized.geojson",
}


def raw_url(repo, branch, repo_path):
    return f"https://raw.githubusercontent.com/{repo}/{branch}/{repo_path}"


def download_file(url, path, overwrite=False, timeout=120):
    if os.path.exists(path) and not overwrite:
        print(f"SKIP  {path}  (already exists; use --overwrite to replace)")
        return "skip"

    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    req = urllib.request.Request(
        url,
        headers={"User-Agent": "gem-faults-downloader/1.0"}
    )

    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            total = resp.length or 0
            read = 0
            chunk = 1 << 20  # 1 MiB
            with open(path, "wb") as f:
                while True:
                    buf = resp.read(chunk)
                    if not buf:
                        break
                    f.write(buf)
                    read += len(buf)
                    if total:
                        pct = 100.0 * read / total
                        print(f"\r  {read/1e6:7.1f} MB / {total/1e6:7.1f} MB "
                              f"({pct:5.1f}%)", end="", file=sys.stderr)
                    else:
                        print(f"\r  {read/1e6:7.1f} MB", end="", file=sys.stderr)
            print("", file=sys.stderr)

        print(f"SAVED {path}")
        return "saved"

    except urllib.error.HTTPError as e:
        print(f"FAIL  {url}  HTTP {e.code}", file=sys.stderr)
        if e.code == 404:
            print("      The file path or branch may have changed. Check the "
                  "repository and rerun with --path and/or --branch:\n"
                  f"      https://github.com/{DEFAULT_REPO}", file=sys.stderr)
        return "failed"
    except urllib.error.URLError as e:
        print(f"FAIL  {url}  {e}", file=sys.stderr)
        return "failed"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download the GEM Global Active Faults Database, keeping GEM's own\n"
            "nomenclature (slip_type etc.) unchanged.\n\n"
            "Source: https://github.com/GEMScienceTools/gem-global-active-faults\n"
            "Licence: CC BY-SA 4.0 (attribution + share-alike).\n\n"
            "The kinematic field is `slip_type` with values such as Normal,\n"
            "Reverse, Dextral, Sinistral and oblique compounds. These are left\n"
            "exactly as published; classify to thrust/normal/transform yourself."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-o", "--output",
        default=".",
        help="Output file path, or a directory (the published filename is "
             "used). Default: current directory."
    )
    parser.add_argument(
        "--format",
        choices=sorted(DATASET_PATHS.keys()),
        default="gpkg",
        help="File format to fetch: gpkg (default) or geojson."
    )
    parser.add_argument(
        "--path",
        default=None,
        help="Override the repo-relative path of the file to download "
             "(e.g. gpkg/gem_active_faults_harmonized.gpkg). Takes precedence "
             "over --format."
    )
    parser.add_argument(
        "--repo",
        default=DEFAULT_REPO,
        help=f"GitHub owner/repo (default: {DEFAULT_REPO})."
    )
    parser.add_argument(
        "--branch",
        default=DEFAULT_BRANCH,
        help=f"Branch or tag (default: {DEFAULT_BRANCH}; try 'main' if master 404s)."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing output file."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the URL and output path without downloading."
    )

    args = parser.parse_args()

    repo_path = args.path if args.path else DATASET_PATHS[args.format]
    url = raw_url(args.repo, args.branch, repo_path)

    # Resolve the output path: if --output is a directory (existing or with a
    # trailing slash), append the published filename.
    out = args.output
    if out.endswith(os.sep) or os.path.isdir(out):
        out = os.path.join(out, os.path.basename(repo_path))

    print(f"Repository : {args.repo} @ {args.branch}")
    print(f"File       : {repo_path}")
    print(f"URL        : {url}")
    print(f"Output     : {out}")

    if args.dry_run:
        return

    result = download_file(url, out, overwrite=args.overwrite)

    if result == "saved":
        print("\nDone. GEM nomenclature preserved (see the `slip_type` field).")
        print("Attribution required: GEM Global Active Faults Database "
              "(Styron & Pagani, 2020), CC BY-SA 4.0.")
    elif result == "skip":
        print("\nNothing downloaded (file already present).")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
