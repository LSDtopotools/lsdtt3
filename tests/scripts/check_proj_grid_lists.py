#!/usr/bin/env python3
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Check tools/fetch_proj_grids.py against the C++ default grid list, offline.

1. DEFAULT_GRIDS in tools/fetch_proj_grids.py names exactly the grids in
   default_proj_grids() (src/core/lsdtt_proj_grid_check.cpp), each with a
   SHA-256 checksum.
2. fetch, list and where work against a local file:// "CDN": a grid is
   installed, a second fetch keeps it, a default grid whose content does not
   match its checksum is refused and not left behind, a non-GeoTIFF and a bad
   name are refused, and a missing grid exits non-zero.

Standard library only (tier_python_stdlib).

Usage:
    check_proj_grid_lists.py --repo-root <dir>
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import re
import sys
import tempfile
from pathlib import Path


def load_tool(root: Path):
    spec = importlib.util.spec_from_file_location("fetch_proj_grids", root / "tools" / "fetch_proj_grids.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cpp_defaults(root: Path) -> list[str]:
    text = (root / "src" / "core" / "lsdtt_proj_grid_check.cpp").read_text(encoding="utf-8")
    m = re.search(r"default_proj_grids\(\)\s*\{.*?grids\s*=\s*\{(.*?)\};", text, re.DOTALL)
    if m is None:
        raise SystemExit("ERROR: could not find default_proj_grids() in lsdtt_proj_grid_check.cpp")
    return re.findall(r'"([^"]+)"', m.group(1))


def run(tool, *argv) -> tuple[int, str]:
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = tool.main(list(argv))
        except SystemExit as e:
            code = e.code if isinstance(e.code, int) else 1
    return code, out.getvalue() + err.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    root = parser.parse_args().repo_root.resolve()
    tool = load_tool(root)
    failures = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        print(("  PASS  " if ok else "  FAIL  ") + name + (f"  ({detail})" if detail else ""))
        if not ok:
            failures.append(name)

    cpp = cpp_defaults(root)
    py = list(tool.DEFAULT_GRIDS)
    check("C++ and Python default grid lists agree", cpp == py, f"{cpp} vs {py}")
    check("every default grid has a SHA-256",
          all(isinstance(v, str) and re.fullmatch(r"[0-9a-f]{64}", v) for v in tool.DEFAULT_GRIDS.values()))

    with tempfile.TemporaryDirectory() as tmp:
        cdn = Path(tmp) / "cdn"
        target = Path(tmp) / "share" / "proj"
        cdn.mkdir()
        target.mkdir(parents=True)
        (target / "proj.db").write_bytes(b"")
        (cdn / "xx_test_grid.tif").write_bytes(b"II*\x00" + b"\x00" * 64)
        (cdn / "not_tiff.tif").write_bytes(b"<html>not found</html>")
        default = py[0]
        (cdn / default).write_bytes(b"II*\x00 wrong content")
        base = cdn.as_uri()

        code, text = run(tool, "--dir", str(target), "fetch", "xx_test_grid.tif", "--base-url", base)
        check("fetch installs a grid", code == 0 and (target / "xx_test_grid.tif").is_file(), text.strip())
        code, text = run(tool, "--dir", str(target), "fetch", "xx_test_grid.tif", "--base-url", base)
        check("a second fetch keeps it", code == 0 and "already installed" in text)
        code, text = run(tool, "--dir", str(target), "fetch", default, "--base-url", base)
        check("a default grid failing its checksum is refused and not left behind",
              code != 0 and "checksum" in text and not (target / default).exists())
        code, text = run(tool, "--dir", str(target), "fetch", "not_tiff.tif", "--base-url", base)
        check("a non-GeoTIFF is refused", code != 0 and not (target / "not_tiff.tif").exists())
        code, text = run(tool, "--dir", str(target), "fetch", "../evil.tif", "--base-url", base)
        check("a name with a path is refused", code != 0)
        code, text = run(tool, "--dir", str(target), "fetch", "missing_grid.tif", "--base-url", base)
        check("a missing grid exits non-zero", code != 0)
        code, text = run(tool, "--dir", str(target), "list")
        check("list shows the grid and the missing defaults",
              code == 0 and "xx_test_grid.tif" in text and "missing defaults" in text)
        code, text = run(tool, "--dir", str(target), "where")
        check("where prints the directory", code == 0 and text.strip() == str(target.resolve()))
        check("no partial downloads left", not list(target.glob(".download-*")))

    print(f"RESULT: {'PASS' if not failures else 'FAIL'}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
