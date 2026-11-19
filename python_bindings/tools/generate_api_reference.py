#!/usr/bin/env python3

"""Generate complete domain-grouped Markdown from the runtime API manifest."""

import argparse
import importlib.util
import json
from pathlib import Path
import re
import sys


def _load_topics():
    """Load lsdtt3/_topics.py by path: CI runs this without the extension built."""
    path = Path(__file__).resolve().parents[1] / "src" / "lsdtt3" / "_topics.py"
    spec = importlib.util.spec_from_file_location("lsdtt3_topics", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_topics = _load_topics()
DOMAINS = _topics.DOMAINS
GEOSPATIAL_RUNTIME = _topics.GEOSPATIAL_RUNTIME
domain = _topics.domain


DOMAIN_GUIDANCE = {
    "graphflood": "Canonical convergence names are used throughout. Legacy tol_* keywords warn and old/new conflicts fail, including accelerator backend options.",
    "landscape": "LandscapeModel thresholds use specific catchment area A/dx in metres. Flat relief windows use a full side length in pixels.",
    "fluvial": "Valley identification is upstream of centreline tracing. extract_valley_centrelines requires an explicit, co-registered valley_mask.",
}


def materialize(signature):
    definition = signature["definition"]
    for index in range(len(signature["defaults"]) - 1, -1, -1):
        definition = definition.replace(f"\\{index}", repr(signature["defaults"][index]))
    return definition


def render_export(export, deprecations):
    lines = [f"## `{export['name']}`", "", export["doc"] or "No summary is available."]
    if deprecations:
        lines.extend(["", "**Deprecated compatibility names:**"])
        for item in deprecations:
            lines.append(f"- `{item['old']}` -> `{item['new']}`; earliest removal: {item['earliest_removal']}.")
    if export["signatures"]:
        lines.extend(["", "```python"])
        lines.extend(materialize(signature) for signature in export["signatures"])
        lines.append("```")
    if export["members"]:
        lines.extend(["", "### Members", ""])
        for member in export["members"]:
            lines.extend([f"#### `{member['name']}`", "", member["doc"] or "No summary is available."])
            if member["signatures"]:
                lines.extend(["", "```python"])
                lines.extend(materialize(signature) for signature in member["signatures"])
                lines.append("```")
            lines.append("")
    return "\n".join(lines).rstrip()


def render_page(key, title, exports, deprecated):
    lines = [
        f"# {title} API",
        "",
        "This page is generated from the built Python runtime manifest. Do not edit signatures here.",
        "",
    ]
    if key in DOMAIN_GUIDANCE:
        lines.extend([DOMAIN_GUIDANCE[key], ""])
    lines.append("\n\n---\n\n".join(
        render_export(export, deprecated.get(export["name"], [])) for export in exports
    ))
    rendered = "\n".join(lines).rstrip() + "\n"
    return "\n".join(line.rstrip() for line in rendered.splitlines()) + "\n"


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--update", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    manifest = json.loads((root / "python_bindings/api/public_api.json").read_text(encoding="utf-8"))
    decisions = json.loads((root / "python_bindings/api/decisions.json").read_text(encoding="utf-8"))
    deprecated = {}
    for item in decisions["deprecated_names"]:
        deprecated.setdefault(item["symbol"], []).append(item)
    grouped = {key: [] for key in DOMAINS}
    for export in manifest["exports"]:
        grouped[domain(export["name"])].append(export)
    assigned = [export["name"] for exports in grouped.values() for export in exports]
    if sorted(assigned) != sorted(export["name"] for export in manifest["exports"]):
        raise RuntimeError("API reference assignment is not exact")
    docs = root / "docs/python_api"
    changed = []
    for key, title in DOMAINS.items():
        path = docs / f"{key}.md"
        rendered = render_page(key, title, grouped[key], deprecated)
        if not path.exists() or path.read_text(encoding="utf-8") != rendered:
            changed.append(path)
            if args.update:
                path.write_text(rendered, encoding="utf-8")
    if args.check and changed:
        raise RuntimeError(f"generated Python API pages are stale: {[path.name for path in changed]}")
    print(f"{'updated' if args.update else 'checked'} {len(DOMAINS)} pages for {len(assigned)} exports")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
