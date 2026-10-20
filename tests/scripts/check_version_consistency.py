#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Check that every copy of the lsdtt3 version number agrees.

The source of truth is project(lsdtt3 VERSION X.Y.Z) in the root
CMakeLists.txt. These files carry checked copies:

  python_bindings/pyproject.toml          version = "X.Y.Z" (may add .postN)
  qgis_plugin/lsdtt_qgis/metadata.txt     version=X.Y.Z
  conda-recipe/meta.yaml                  {% set version = "X.Y.Z" %}
  pixi.toml                               [workspace] version = "X.Y.Z"
  debian/changelog                        top entry, X.Y.Z (may add +pkgN)

Only the X.Y.Z stems are compared. A wheel-only rebuild is X.Y.Z.postN in
pyproject.toml and a Debian-only rebuild is X.Y.Z+pkgN; conda rebuilds use the
recipe's build number, which this check ignores. A leading "v" is rejected
everywhere: it belongs to the git tag, not to the files.

It also checks that the PyPI page states why some versions have no wheel:
pyproject.toml must declare readme = "README.md", and python_bindings/README.md
must have a "## Versioning" section.

tools/set_version.py rewrites all of these together.

Standard library only (tier_python_stdlib).

Usage:
    check_version_consistency.py --repo-root <dir>
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

STEM = r"(\d+\.\d+\.\d+)"


def _search(path: Path, pattern: str, what: str, flags: int = re.MULTILINE) -> re.Match:
    text = path.read_text(encoding="utf-8")
    match = re.search(pattern, text, flags)
    if match is None:
        raise ValueError(f"{path}: could not find {what}")
    return match


def read_versions(root: Path) -> dict[str, str]:
    """Return {file: full version string} for every copy, CMake first."""
    found: dict[str, str] = {}
    found["CMakeLists.txt"] = _search(
        root / "CMakeLists.txt",
        r"^project\(lsdtt3\s+VERSION\s+(\S+?)[\s)]",
        "project(lsdtt3 VERSION ...)").group(1)

    pyproject = root / "python_bindings" / "pyproject.toml"
    project_table = _search(pyproject, r"^\[project\]\s*$(.*?)(?=^\[|\Z)",
                            "[project] table", re.MULTILINE | re.DOTALL).group(1)
    match = re.search(r'^version\s*=\s*"([^"]*)"', project_table, re.MULTILINE)
    if match is None:
        raise ValueError(f"{pyproject}: no static version in [project]")
    found["python_bindings/pyproject.toml"] = match.group(1)

    found["qgis_plugin/lsdtt_qgis/metadata.txt"] = _search(
        root / "qgis_plugin" / "lsdtt_qgis" / "metadata.txt",
        r"^version\s*=\s*(\S+)\s*$", "version=").group(1)

    found["conda-recipe/meta.yaml"] = _search(
        root / "conda-recipe" / "meta.yaml",
        r'^\{%\s*set\s+version\s*=\s*"([^"]*)"\s*%\}', "{% set version = ... %}").group(1)

    pixi = root / "pixi.toml"
    workspace = _search(pixi, r"^\[workspace\]\s*$(.*?)(?=^\[|\Z)",
                        "[workspace] table", re.MULTILINE | re.DOTALL).group(1)
    match = re.search(r'^version\s*=\s*"([^"]*)"', workspace, re.MULTILINE)
    if match is None:
        raise ValueError(f"{pixi}: no version in [workspace]")
    found["pixi.toml"] = match.group(1)

    first_line = (root / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0]
    match = re.match(r"^lsdtt3 \(([^)]*)\)", first_line)
    if match is None:
        raise ValueError(f"debian/changelog: top line is not 'lsdtt3 (<version>) ...': {first_line!r}")
    found["debian/changelog"] = match.group(1)
    return found


# Allowed shape of each copy: the stem plus the one suffix its ecosystem uses.
ALLOWED = {
    "CMakeLists.txt": rf"^{STEM}$",
    "python_bindings/pyproject.toml": rf"^{STEM}(\.post\d+)?$",
    "qgis_plugin/lsdtt_qgis/metadata.txt": rf"^{STEM}$",
    "conda-recipe/meta.yaml": rf"^{STEM}$",
    "pixi.toml": rf"^{STEM}$",
    "debian/changelog": rf"^{STEM}(\+pkg\d+)?$",
}


def check_versions(found: dict[str, str]) -> list[str]:
    errors = []
    stems = {}
    for name, value in found.items():
        match = re.match(ALLOWED[name], value)
        if match is None:
            errors.append(f"{name}: '{value}' is not of the form {ALLOWED[name]}")
        else:
            stems[name] = match.group(1)
    truth = stems.get("CMakeLists.txt")
    if truth is not None:
        for name, stem in stems.items():
            if stem != truth:
                errors.append(f"{name}: version '{found[name]}' differs from "
                              f"CMakeLists.txt '{truth}'")
    return errors


def check_pypi_readme(root: Path) -> list[str]:
    errors = []
    pyproject = (root / "python_bindings" / "pyproject.toml").read_text(encoding="utf-8")
    if not re.search(r'^readme\s*=\s*"README\.md"\s*$', pyproject, re.MULTILINE):
        errors.append('python_bindings/pyproject.toml: missing readme = "README.md" '
                      "(PyPI would show no description)")
    readme = (root / "python_bindings" / "README.md").read_text(encoding="utf-8")
    if not re.search(r"^## Versioning\s*$", readme, re.MULTILINE):
        errors.append("python_bindings/README.md: missing the '## Versioning' section "
                      "that explains why some lsdtt3 versions have no wheel")
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repo-root", type=Path, required=True)
    args = parser.parse_args(argv)
    root = args.repo_root.resolve()

    try:
        found = read_versions(root)
    except (OSError, ValueError, IndexError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    errors = check_versions(found) + check_pypi_readme(root)
    for name, value in found.items():
        print(f"  {name:40s} {value}")
    if errors:
        for error in errors:
            print(f"FAIL: {error}", file=sys.stderr)
        print("Fix with: tools/set_version.py <X.Y.Z>", file=sys.stderr)
        return 1
    print("PASS: all version copies agree")
    return 0


if __name__ == "__main__":
    sys.exit(main())
