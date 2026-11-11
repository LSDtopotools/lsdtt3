# SPDX-License-Identifier: AGPL-3.0-or-later
"""Check a directory of release wheels before they are uploaded to a package index.

Run by the PyPI upload jobs (.github/workflows/release.yml and
python-wheels.yml) after the wheels are downloaded, before anything is sent.
An index accepts a version only once, so a wrong set of files is refused here,
where a re-run costs nothing:

  - exactly one wheel for each release platform, and nothing else;
  - every wheel is cp312-abi3 (pyproject.toml wheel.py-api) and repaired
    (a manylinux tag on Linux, not linux_x86_64);
  - every wheel carries the version of the tag being released;
  - every file is under PyPI's default 100 MB per-file limit.

    python check_dist_wheels.py --dir dist --tag v0.9.4

Exits 0 if every check passes, 1 otherwise. Uses only the standard library.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# One wheel per python-wheels.yml leg. A pattern on the platform tag.
PLATFORMS = {
    "Linux x86_64": re.compile(r"^manylinux_\d+_\d+_x86_64$"),
    "macOS ARM64": re.compile(r"^macosx_\d+_\d+_arm64$"),
    "macOS Intel": re.compile(r"^macosx_\d+_\d+_x86_64$"),
    "Windows x86_64": re.compile(r"^win_amd64$"),
}
PYPI_FILE_LIMIT = 100 * 1024 * 1024
# name-version-python-abi-platform.whl (no build tag in these wheels)
WHEEL_RE = re.compile(
    r"^(?P<name>[^-]+)-(?P<version>[^-]+)-(?P<py>[^-]+)-(?P<abi>[^-]+)-(?P<plat>[^-]+)\.whl$"
)


def check(directory: Path, tag: str) -> list[str]:
    errors: list[str] = []
    version = tag[1:] if tag.startswith("v") else tag
    files = sorted(p for p in directory.iterdir() if p.is_file())
    if not files:
        return [f"no files in {directory}"]

    found: dict[str, list[str]] = {name: [] for name in PLATFORMS}
    for path in files:
        print(f"  {path.stat().st_size / 1048576:7.1f} MB  {path.name}")
        m = WHEEL_RE.match(path.name)
        if not m:
            errors.append(f"{path.name}: not a wheel file name")
            continue
        if m["name"] != "lsdtt3":
            errors.append(f"{path.name}: distribution is {m['name']}, not lsdtt3")
        if m["version"] != version:
            errors.append(f"{path.name}: version {m['version']}, but the tag is {tag}")
        if (m["py"], m["abi"]) != ("cp312", "abi3"):
            errors.append(f"{path.name}: {m['py']}-{m['abi']}, expected cp312-abi3")
        if path.stat().st_size >= PYPI_FILE_LIMIT:
            errors.append(f"{path.name}: {path.stat().st_size} bytes, over PyPI's 100 MB limit")
        # A wheel can carry several platform tags joined by '.'; any match counts.
        plats = m["plat"].split(".")
        legs = [leg for leg, pat in PLATFORMS.items() if any(pat.match(t) for t in plats)]
        if not legs:
            errors.append(f"{path.name}: platform {m['plat']} is not a release platform "
                          "(a linux_x86_64 tag means auditwheel did not run)")
        for leg in legs:
            found[leg].append(path.name)

    for leg, names in found.items():
        if len(names) != 1:
            errors.append(f"{leg}: {len(names)} wheels, expected 1 ({', '.join(names) or 'none'})")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dir", required=True, type=Path, help="directory holding the wheels")
    parser.add_argument("--tag", required=True, help="the release tag, vX.Y.Z or vX.Y.Z.postN")
    args = parser.parse_args()
    print(f"Wheels in {args.dir} for {args.tag}:")
    errors = check(args.dir, args.tag)
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)
    if errors:
        return 1
    print(f"PASS: one cp312-abi3 wheel for each of the {len(PLATFORMS)} platforms, version {args.tag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
