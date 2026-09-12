#!/usr/bin/env python3

"""Generate deterministic CLI parameter reference sections from `-h json`."""

import argparse
import json
import pathlib
import re
import subprocess
import sys


BEGIN_MARKER = "<!-- BEGIN GENERATED CLI PARAMETERS -->"
END_MARKER = "<!-- END GENERATED CLI PARAMETERS -->"
PARAMETER_FIELDS = (
    "name",
    "type",
    "default",
    "description",
    "guidance",
    "group",
    "role",
    "parent",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--update", action="store_true", help="Rewrite generated sections")
    mode.add_argument("--check", action="store_true", help="Fail if generated sections differ")
    parser.add_argument("--bin-dir", required=True, type=pathlib.Path)
    parser.add_argument(
        "--repo-root",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parents[1],
    )
    return parser.parse_args()


def fail(message: str) -> None:
    raise RuntimeError(message)


def expected_page(source: pathlib.Path, docs_dir: pathlib.Path) -> pathlib.Path:
    suffix = source.stem.replace("-", "_")
    return docs_dir / f"cli_reference_{suffix}.md"


def validate_inventory(cli_dir: pathlib.Path, docs_dir: pathlib.Path) -> list[pathlib.Path]:
    sources = sorted(cli_dir.glob("lsdtt-*.cpp"), key=lambda path: path.name)
    pages = set(docs_dir.glob("cli_reference_lsdtt_*.md"))
    expected = {expected_page(source, docs_dir) for source in sources}
    if pages != expected:
        missing = sorted(str(path) for path in expected - pages)
        extra = sorted(str(path) for path in pages - expected)
        fail(f"CLI source/reference inventory mismatch\nmissing: {missing}\nextra: {extra}")
    return sources


def load_help(executable: pathlib.Path, source: pathlib.Path) -> dict:
    if not executable.is_file():
        fail(f"missing CLI binary: {executable}")
    if executable.stat().st_mtime < source.stat().st_mtime:
        fail(f"stale CLI binary: {executable} is older than {source}")

    result = subprocess.run(
        [str(executable), "-h", "json"], capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        fail(
            f"JSON help failed for {executable.name} ({result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    start = result.stdout.find("{")
    if start < 0:
        fail(f"JSON help for {executable.name} contains no object")
    try:
        payload = json.loads(result.stdout[start:])
    except json.JSONDecodeError as error:
        fail(f"invalid JSON help for {executable.name}: {error}")

    for field in ("tool", "version", "citation", "features", "parameters"):
        if field not in payload:
            fail(f"JSON help for {executable.name} lacks top-level field {field!r}")
    if not all(isinstance(payload[field], str) for field in ("tool", "version", "citation")):
        fail(f"JSON help for {executable.name} has invalid string metadata")
    if not isinstance(payload["features"], list) or not all(
        isinstance(value, str) for value in payload["features"]
    ):
        fail(f"JSON help for {executable.name} has invalid features")
    if not isinstance(payload["parameters"], list):
        fail(f"JSON help for {executable.name} has invalid parameters")

    names = []
    for index, parameter in enumerate(payload["parameters"]):
        if not isinstance(parameter, dict):
            fail(f"parameter {index} for {executable.name} is not an object")
        missing = [field for field in PARAMETER_FIELDS if field not in parameter]
        if missing:
            fail(f"parameter {index} for {executable.name} lacks fields {missing}")
        if not all(isinstance(parameter[field], str) for field in PARAMETER_FIELDS):
            fail(f"parameter {index} for {executable.name} has non-string fields")
        names.append(parameter["name"])
    if len(names) != len(set(names)):
        fail(f"JSON help for {executable.name} contains duplicate parameter names")
    return payload


def markdown_text(value: str) -> str:
    text = " ".join(value.split())
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return text.replace("|", "\\|") or "-"


def code_text(value: str) -> str:
    return f"`{value.replace('`', '')}`" if value else "-"


def render(payload: dict) -> str:
    lines = [
        BEGIN_MARKER,
        "## Parameters",
        "",
        "This section is generated from the executable's `-h json` registry. Do not edit it manually.",
    ]
    groups = []
    for parameter in payload["parameters"]:
        group = parameter["group"] or "General"
        if group not in groups:
            groups.append(group)
    for group in groups:
        lines.extend(
            [
                "",
                f"### {markdown_text(group)}",
                "",
                "| Parameter | Type | Default | Role | Parent | Description | Guidance |",
                "|---|---|---|---|---|---|---|",
            ]
        )
        for parameter in payload["parameters"]:
            if (parameter["group"] or "General") != group:
                continue
            lines.append(
                "| "
                + " | ".join(
                    (
                        code_text(parameter["name"]),
                        code_text(parameter["type"]),
                        code_text(parameter["default"]),
                        code_text(parameter["role"]),
                        code_text(parameter["parent"]),
                        markdown_text(parameter["description"]),
                        markdown_text(parameter["guidance"]),
                    )
                )
                + " |"
            )
    lines.extend(["", END_MARKER])
    return "\n".join(lines)


def replace_generated_section(original: str, generated: str, page: pathlib.Path) -> str:
    marker_pattern = re.compile(
        re.escape(BEGIN_MARKER) + r".*?" + re.escape(END_MARKER), re.DOTALL
    )
    matches = list(marker_pattern.finditer(original))
    if len(matches) > 1:
        fail(f"multiple generated parameter sections in {page}")
    if matches:
        return marker_pattern.sub(generated, original)
    return original.rstrip() + "\n\n" + generated + "\n"


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    bin_dir = args.bin_dir.resolve()
    cli_dir = root / "cli"
    docs_dir = root / "docs" / "cli_reference"
    sources = validate_inventory(cli_dir, docs_dir)

    changed = []
    total_parameters = 0
    for source in sources:
        executable = bin_dir / source.stem
        payload = load_help(executable, source)
        total_parameters += len(payload["parameters"])
        page = expected_page(source, docs_dir)
        original = page.read_text(encoding="utf-8")
        updated = replace_generated_section(original, render(payload), page)
        if updated != original:
            changed.append(page)
            if args.update:
                page.write_text(updated, encoding="utf-8")

    if args.check and changed:
        relative = "\n".join(str(path.relative_to(root)) for path in changed)
        fail(f"CLI reference sections are stale or missing:\n{relative}")

    action = "updated" if args.update else "checked"
    print(f"{action} {len(sources)} CLI reference pages ({total_parameters} parameters)")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
