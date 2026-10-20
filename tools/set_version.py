#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Set the lsdtt3 version number in every file that carries it.

lsdtt3 has one version number. project(lsdtt3 VERSION X.Y.Z) in CMakeLists.txt
is the truth; this script writes the same X.Y.Z into the checked copies:

  CMakeLists.txt                          project(lsdtt3 VERSION X.Y.Z ...)
  python_bindings/pyproject.toml          version = "X.Y.Z"
  qgis_plugin/lsdtt_qgis/metadata.txt     version=X.Y.Z
  conda-recipe/meta.yaml                  {% set version = "X.Y.Z" %}
  pixi.toml                               version = "X.Y.Z"
  debian/changelog                        a new top entry

and then runs the same check as the unit_version_consistency CTest test.

A new X.Y.Z also resets the conda build number to 0, clears any .postN from
pyproject.toml, and adds a "Release X.Y.Z." changelog entry.

Rebuilds of an unchanged source keep X.Y.Z and add only an ecosystem suffix:

  --post N              pyproject.toml becomes X.Y.Z.postN (a repackaged
                        wheel; tag vX.Y.Z.postN). Nothing else changes.
  --debian-revision N   a new changelog entry X.Y.Z+pkgN. Nothing else changes.

Both need X.Y.Z equal to the current CMake version. A conda-only rebuild is
an edit of build: number: in conda-recipe/meta.yaml.

The changelog entry is written in the format `dch` produces. The distribution
and maintainer are copied from the previous entry, unless DEBFULLNAME and
DEBEMAIL are set.

Usage:
    tools/set_version.py 0.7.0
    tools/set_version.py 0.7.0 --post 1
    tools/set_version.py 0.7.0 --debian-revision 2
    tools/set_version.py 0.7.0 --dry-run
"""

from __future__ import annotations

import argparse
import email.utils
import os
from pathlib import Path
import re
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tests" / "scripts"))
import check_version_consistency  # noqa: E402


class Edit:
    """One file rewritten in memory, written only if every edit succeeds."""

    def __init__(self, root: Path, relpath: str):
        self.path = root / relpath
        self.relpath = relpath
        self.original = self.path.read_text(encoding="utf-8")
        self.text = self.original

    def sub(self, pattern: str, replacement, what: str, flags: int = re.MULTILINE) -> None:
        new, count = re.subn(pattern, replacement, self.text, count=1, flags=flags)
        if count != 1:
            raise SystemExit(f"ERROR: {self.relpath}: could not find {what}")
        self.text = new

    @property
    def changed(self) -> bool:
        return self.text != self.original


def _in_table(edit: Edit, table: str, key_pattern: str, value: str, what: str) -> None:
    """Replace the quoted value of a key inside one TOML table only."""
    match = re.search(rf"^\[{re.escape(table)}\]\s*$(.*?)(?=^\[|\Z)", edit.text,
                      re.MULTILINE | re.DOTALL)
    if match is None:
        raise SystemExit(f"ERROR: {edit.relpath}: no [{table}] table")
    body = match.group(1)
    new_body, count = re.subn(rf'^({key_pattern}\s*=\s*")[^"]*(")', rf"\g<1>{value}\g<2>",
                              body, count=1, flags=re.MULTILINE)
    if count != 1:
        raise SystemExit(f"ERROR: {edit.relpath}: could not find {what} in [{table}]")
    edit.text = edit.text[:match.start(1)] + new_body + edit.text[match.end(1):]


def changelog_entry(previous: str, version: str, message: str) -> str:
    first = previous.splitlines()[0]
    match = re.match(r"^lsdtt3 \([^)]*\) (\S+); (urgency=\S+)", first)
    if match is None:
        raise SystemExit(f"ERROR: debian/changelog: unexpected top line {first!r}")
    distribution, urgency = match.groups()
    name, mail = os.environ.get("DEBFULLNAME"), os.environ.get("DEBEMAIL")
    if name and mail:
        maintainer = f"{name} <{mail}>"
    else:
        trailer = re.search(r"^ -- (.+?>)  ", previous, re.MULTILINE)
        if trailer is None:
            raise SystemExit("ERROR: debian/changelog: no maintainer line in the top "
                             "entry; set DEBFULLNAME and DEBEMAIL")
        maintainer = trailer.group(1)
    date = email.utils.formatdate(localtime=True)
    return (f"lsdtt3 ({version}) {distribution}; {urgency}\n\n"
            f"  * {message}\n\n"
            f" -- {maintainer}  {date}\n\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(__doc__.splitlines()[2:]))
    parser.add_argument("version", help="X.Y.Z, without a leading v")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--post", type=int, metavar="N",
                       help="wheel-only rebuild: pyproject.toml becomes X.Y.Z.postN")
    group.add_argument("--debian-revision", type=int, metavar="N",
                       help="Debian-only rebuild: changelog entry X.Y.Z+pkgN")
    parser.add_argument("--dry-run", action="store_true",
                        help="show what would change, write nothing")
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    root = args.repo_root.resolve()

    version = args.version
    if version.startswith("v"):
        raise SystemExit(f"ERROR: '{version}': give the bare number; the 'v' belongs "
                         "to the git tag only")
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise SystemExit(f"ERROR: '{version}' is not X.Y.Z (digits only; suffixes "
                         "come from --post / --debian-revision)")
    for value, flag in ((args.post, "--post"), (args.debian_revision, "--debian-revision")):
        if value is not None and value < 1:
            raise SystemExit(f"ERROR: {flag} must be 1 or more")

    current = check_version_consistency.read_versions(root)["CMakeLists.txt"]
    rebuild = args.post is not None or args.debian_revision is not None
    if rebuild and version != current:
        raise SystemExit(f"ERROR: --post / --debian-revision repackage the current "
                         f"version {current}; got {version}. Set a new version first.")

    edits: list[Edit] = []
    if args.post is not None:
        pyproject = Edit(root, "python_bindings/pyproject.toml")
        _in_table(pyproject, "project", "version", f"{version}.post{args.post}", "version")
        edits.append(pyproject)
    elif args.debian_revision is not None:
        changelog = Edit(root, "debian/changelog")
        full = f"{version}+pkg{args.debian_revision}"
        if changelog.text.startswith(f"lsdtt3 ({full})"):
            raise SystemExit(f"ERROR: debian/changelog already starts with {full}")
        changelog.text = changelog_entry(changelog.text, full,
                                         f"Repackaging of {version}.") + changelog.text
        edits.append(changelog)
    else:
        cmake = Edit(root, "CMakeLists.txt")
        cmake.sub(r"^(project\(lsdtt3\s+VERSION\s+)\S+?(?=[\s)])", rf"\g<1>{version}",
                  "project(lsdtt3 VERSION ...)")
        pyproject = Edit(root, "python_bindings/pyproject.toml")
        _in_table(pyproject, "project", "version", version, "version")
        metadata = Edit(root, "qgis_plugin/lsdtt_qgis/metadata.txt")
        metadata.sub(r"^(version\s*=\s*)\S+", rf"\g<1>{version}", "version=")
        conda = Edit(root, "conda-recipe/meta.yaml")
        conda.sub(r'^(\{%\s*set\s+version\s*=\s*")[^"]*(")', rf"\g<1>{version}\g<2>",
                  "{% set version = ... %}")
        pixi = Edit(root, "pixi.toml")
        _in_table(pixi, "workspace", "version", version, "version")
        edits += [cmake, pyproject, metadata, conda, pixi]
        if version != current:
            conda.sub(r"^(build:\s*\n\s+number:\s*)\d+", r"\g<1>0", "build: number:")
            changelog = Edit(root, "debian/changelog")
            if not changelog.text.startswith(f"lsdtt3 ({version})"):
                changelog.text = changelog_entry(changelog.text, version,
                                                 f"Release {version}.") + changelog.text
            edits.append(changelog)

    for edit in edits:
        state = "changed" if edit.changed else "unchanged"
        print(f"  {edit.relpath:40s} {state}")
    if args.dry_run:
        print("Dry run: nothing written.")
        return 0
    for edit in edits:
        if edit.changed:
            edit.path.write_text(edit.text, encoding="utf-8")

    return check_version_consistency.main(["--repo-root", str(root)])


if __name__ == "__main__":
    sys.exit(main())
