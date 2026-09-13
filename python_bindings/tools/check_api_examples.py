#!/usr/bin/env python3

"""Compile Python examples and reject calls to unknown lsdtt3 exports."""

import ast
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "api" / "public_api.json").read_text(encoding="utf-8"))
    decisions = json.loads((root / "api" / "decisions.json").read_text(encoding="utf-8"))
    export_entries = {entry["name"]: entry for entry in manifest["exports"]}
    exports = set(export_entries)
    contracts = {}
    for entry in manifest["exports"]:
        definitions = [signature["definition"] for signature in entry["signatures"]]
        contracts[entry["name"]] = {
            "parameters": set().union(*(
                set(re.findall(r"\b([A-Za-z_]\w*)(?=:)", definition))
                for definition in definitions
            )),
            "variadic_keywords": any("**" in definition for definition in definitions),
        }
    banned = {
        decision["old"] for decision in decisions["completed_hard_renames"]
    }
    banned.update(decision["old"] for decision in decisions["deprecated_names"])

    unknown = []
    stale = []
    invalid_keywords = []

    def check_tree(tree, label):
        inferred_types = {}
        for node in ast.walk(tree):
            if (isinstance(node, ast.Assign) and len(node.targets) == 1
                    and isinstance(node.targets[0], ast.Name)
                    and isinstance(node.value, ast.Call)
                    and isinstance(node.value.func, ast.Attribute)
                    and isinstance(node.value.func.value, ast.Name)
                    and node.value.func.value.id == "lsdtt3"
                    and node.value.func.attr in exports
                    and export_entries[node.value.func.attr]["kind"] == "class"):
                inferred_types[node.targets[0].id] = node.value.func.attr
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)):
                continue
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "lsdtt3":
                symbol = node.func.attr
                callable_contract = contracts.get(symbol)
            elif isinstance(node.func.value, ast.Name) and node.func.value.id in inferred_types:
                class_name = inferred_types[node.func.value.id]
                class_entry = export_entries[class_name]
                member = next(
                    (item for item in class_entry["members"] if item["name"] == node.func.attr),
                    None,
                )
                if member is None:
                    unknown.append(f"{label}:{node.lineno}:{class_name}.{node.func.attr}")
                    continue
                definitions = [signature["definition"] for signature in member["signatures"]]
                callable_contract = {
                    "parameters": set().union(*(
                        set(re.findall(r"\b([A-Za-z_]\w*)(?=:)", definition))
                        for definition in definitions
                    )),
                    "variadic_keywords": any("**" in definition for definition in definitions),
                }
                symbol = f"{class_name}.{node.func.attr}"
            else:
                continue
            if symbol not in exports:
                if "." not in symbol:
                    unknown.append(f"{label}:{node.lineno}:{symbol}")
                    continue
            if callable_contract is None or callable_contract["variadic_keywords"]:
                continue
            for keyword in node.keywords:
                if keyword.arg is not None and keyword.arg not in callable_contract["parameters"]:
                    invalid_keywords.append(
                        f"{label}:{node.lineno}:{symbol}.{keyword.arg}"
                    )

    for path in sorted((root / "examples").glob("*.py")):
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        compile(tree, str(path), "exec")
        check_tree(tree, path.name)
        for name in banned:
            if re.search(rf"\b{re.escape(name)}\b", source):
                stale.append(f"{path.name}:{name}")

    docs_root = root.parents[0] / "docs" / "python"
    for path in sorted(docs_root.glob("*.md")):
        source = path.read_text(encoding="utf-8")
        for index, block in enumerate(re.findall(r"```python\n(.*?)```", source, re.DOTALL), 1):
            try:
                tree = ast.parse(block, filename=str(path))
            except SyntaxError as error:
                raise RuntimeError(f"invalid Python block in {path.name}: {error}") from error
            check_tree(tree, f"{path.name}:python-block-{index}")
    if unknown or stale or invalid_keywords:
        raise RuntimeError(
            "invalid example API usage\n"
            f"unknown: {unknown}\nstale: {stale}\ninvalid keywords: {invalid_keywords}"
        )
    if args.run:
        with tempfile.TemporaryDirectory(prefix="lsdtt3-api-example-") as work:
            subprocess.run(
                [sys.executable, str(root / "examples" / "api_smoke.py")],
                cwd=work, check=True,
            )
    print("Python examples compile and use known public exports")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
