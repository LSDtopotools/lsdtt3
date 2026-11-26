#!/usr/bin/env python3
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Install PROJ datum-shift and geoid grids where lsdtt3 looks for them.

lsdtt3 carries its own PROJ data (``share/proj``) and makes its tools use it,
so grids installed elsewhere (by ``projsync``, conda or Homebrew) are not seen.
This script downloads grids from the PROJ content delivery network
(https://cdn.proj.org) into that directory. It needs only Python 3; download
it on its own if you do not have the lsdtt3 repository:

    curl -O https://raw.githubusercontent.com/LSDtopotools/lsdtt3/main/tools/fetch_proj_grids.py

    python3 fetch_proj_grids.py where                       # the directory it uses
    python3 fetch_proj_grids.py list                        # grids installed there
    python3 fetch_proj_grids.py fetch us_nga_egm08_25.tif   # add one (or several)
    python3 fetch_proj_grids.py fetch --defaults            # the default set

The directory is found in this order: ``--dir``; the release package that
holds the ``lsdtt-raster-preprocessing`` on your PATH (its ``share/proj``); the
installed Python package (``lsdtt3/share/proj``, for ``pip install lsdtt3``).
Run it in the shell or environment where you use lsdtt3. When an lsdtt3 tool
warns that a transformation needs a grid that is not installed, the warning
names it. See docs/what_is_lsdtt3/note_on_grids.md.

Standard library only: the release packaging scripts run this same file.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

CDN = "https://cdn.proj.org"
USER_AGENT = "lsdtt3-fetch-proj-grids"

# The grids every release package and wheel carries, with their SHA-256. The
# same names are in src/core/lsdtt_proj_grid_check.cpp (default_proj_grids);
# tests/scripts/check_proj_grid_lists.py checks the two lists agree.
DEFAULT_GRIDS = {
    # OSGB36 (British National Grid) <-> ETRS89, Ordnance Survey OSTN15. 2.9 MB.
    "uk_os_OSTN15_NTv2_OSGBtoETRS.tif":
        "5d6ed64d2119952c4c559fa1fccbc594b6520fc3ec3ef2fc10be13202c4384fa",
    # EGM96 geoid, 15': lsdtt3's vertical_datum egm96. 2.6 MB.
    "us_nga_egm96_15.tif":
        "db493027562c9b004d7220fa881f5603adada4e1c5029b933fa7de4547b0e78d",
    # NADCON NAD27 <-> NAD83, contiguous US. 0.2 MB.
    "us_noaa_conus.tif":
        "44611d823c48e5347500ee6afe40ff33d2b88cf817bf59f705ed4a4c3bd687d7",
    # NADCON5 NAD27 <-> NAD83(1986), contiguous US. 0.1 MB.
    "us_noaa_nadcon5_nad27_nad83_1986_conus.tif":
        "c7d587e0d0b39b9f46c7de850b9a6a468c17b7139f82a313aa43aa4a64d94fa8",
}

GRID_SUFFIXES = (".tif", ".tiff", ".gsb", ".gtx", ".byn")


def _version() -> str:
    return "1"


def _is_pip_launcher(path: Path) -> bool:
    try:
        with path.open("rb") as f:
            head = f.read(2)
        if head == b"#!":
            return True
        # pip's Windows launchers sit in a Python environment's Scripts folder.
        return path.parent.name.lower() == "scripts" and (path.parent.parent / "python.exe").is_file()
    except OSError:
        return False


def find_grid_dir(explicit: str | None = None) -> Path:
    """The share/proj directory lsdtt3 uses here (see the module docstring)."""
    if explicit:
        return Path(explicit).expanduser().resolve()
    # A release package: <root>/bin/lsdtt-raster-preprocessing and <root>/share/proj.
    # A pip install puts a Python launcher (a "#!" script, or a .exe launcher
    # on Windows) on PATH under the same name; that is not a release package.
    tool = shutil.which("lsdtt-raster-preprocessing")
    if tool and not _is_pip_launcher(Path(tool)):
        bundle = Path(tool).resolve().parent.parent / "share" / "proj"
        if (bundle / "proj.db").is_file():
            return bundle
    # The Python package (pip install lsdtt3): lsdtt3/share/proj. Found
    # without importing lsdtt3.
    spec = importlib.util.find_spec("lsdtt3")
    if spec and spec.submodule_search_locations:
        for location in spec.submodule_search_locations:
            wheel = Path(location) / "share" / "proj"
            if (wheel / "proj.db").is_file():
                return wheel
    raise SystemExit(
        "fetch_proj_grids.py: cannot find lsdtt3's share/proj. Run it where lsdtt3 is set up "
        "(after source setenv.sh, or in the Python environment with lsdtt3), or give the "
        "directory with --dir <lsdtt3 folder>/share/proj."
    )


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch_grid(name: str, directory: Path, force: bool = False, base_url: str = CDN) -> Path:
    """Download one grid into directory; verify its checksum when known."""
    if "/" in name or "\\" in name or not name.endswith(GRID_SUFFIXES):
        raise ValueError(f"not a PROJ grid file name: {name!r} (e.g. us_nga_egm08_25.tif)")
    target = directory / name
    expected = DEFAULT_GRIDS.get(name)
    if target.is_file() and not force:
        if expected is None or _sha256(target) == expected:
            print(f"  {name}: already installed")
            return target
    directory.mkdir(parents=True, exist_ok=True)
    url = f"{base_url.rstrip('/')}/{name}"
    request = urllib.request.Request(url, headers={"User-Agent": f"{USER_AGENT}/{_version()}"})
    fd, tmp_name = tempfile.mkstemp(prefix=".download-", dir=directory)
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as out, urllib.request.urlopen(request, timeout=120) as response:
            while True:
                block = response.read(1 << 20)
                if not block:
                    break
                out.write(block)
        with tmp.open("rb") as f:
            magic = f.read(4)
        if name.endswith((".tif", ".tiff")) and magic not in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+"):
            raise RuntimeError(f"{url} did not return a GeoTIFF")
        if expected is not None and _sha256(tmp) != expected:
            raise RuntimeError(f"{name}: checksum does not match the expected value")
        tmp.replace(target)
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"{url}: HTTP {e.code} (is the grid name right? see https://cdn.proj.org)") from e
    finally:
        if tmp.exists():
            tmp.unlink()
    print(f"  {name}: installed ({target.stat().st_size / 1048576:.1f} MB)")
    return target


def installed_grids(directory: Path) -> list[Path]:
    return sorted(p for p in directory.iterdir() if p.is_file() and p.name.endswith(GRID_SUFFIXES))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="fetch_proj_grids.py", description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", help="the share/proj directory to use (default: lsdtt3's own)")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("where", help="print the directory grids are installed in")
    sub.add_parser("list", help="list the grids installed there")
    f = sub.add_parser("fetch", help="download grids from https://cdn.proj.org")
    f.add_argument("names", nargs="*", help="grid file names, e.g. us_nga_egm08_25.tif")
    f.add_argument("--defaults", action="store_true", help="the default set carried by every release")
    f.add_argument("--force", action="store_true", help="download again even if installed")
    f.add_argument("--base-url", default=CDN, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    directory = find_grid_dir(args.dir)
    if args.command == "where":
        print(directory)
        return 0
    if args.command == "list":
        grids = installed_grids(directory)
        print(f"{directory}: {len(grids)} grid(s)")
        for p in grids:
            tag = "default" if p.name in DEFAULT_GRIDS else "added"
            print(f"  {p.name}  {p.stat().st_size / 1048576:.1f} MB  ({tag})")
        missing = [g for g in DEFAULT_GRIDS if not (directory / g).is_file()]
        if missing:
            print("missing defaults: " + ", ".join(missing))
        return 0
    names = list(args.names) + (list(DEFAULT_GRIDS) if args.defaults else [])
    if not names:
        parser.error("fetch: name at least one grid, or --defaults")
    print(f"Installing into {directory}")
    failed = 0
    for name in dict.fromkeys(names):
        try:
            fetch_grid(name, directory, force=args.force, base_url=args.base_url)
        except Exception as e:  # report every failure, then exit non-zero
            print(f"  {name}: FAILED: {e}", file=sys.stderr)
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
