#!/usr/bin/env python3

"""Generate/check source-facade signatures against the runtime API manifest."""

import argparse
import ast
import json
from pathlib import Path
import re
import sys


def signature_node(node):
    if isinstance(node, ast.ClassDef):
        node = next(
            (child for child in node.body
             if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
             and child.name == "__init__"),
            None,
        )
        if node is None:
            return None
    return node


def parameter_records(node):
    node = signature_node(node)
    if node is None:
        return []
    positional = list(node.args.posonlyargs) + list(node.args.args)
    default_offset = len(positional) - len(node.args.defaults)
    records = []
    for index, argument in enumerate(positional):
        if argument.arg in {"self", "cls"}:
            continue
        default = None
        if index >= default_offset:
            default = ast.unparse(node.args.defaults[index - default_offset])
        kind = "positional_only" if index < len(node.args.posonlyargs) else "positional_or_keyword"
        records.append({
            "name": argument.arg,
            "kind": kind,
            "annotation": ast.unparse(argument.annotation) if argument.annotation else "",
            "default": default,
        })
    for argument, default_node in zip(node.args.kwonlyargs, node.args.kw_defaults):
        records.append({
            "name": argument.arg,
            "kind": "keyword_only",
            "annotation": ast.unparse(argument.annotation) if argument.annotation else "",
            "default": ast.unparse(default_node) if default_node is not None else None,
        })
    if node.args.vararg is not None:
        records.append({"name": node.args.vararg.arg, "kind": "var_positional", "annotation": "", "default": None})
    if node.args.kwarg is not None:
        records.append({"name": node.args.kwarg.arg, "kind": "var_keyword", "annotation": "", "default": None})
    return records


def runtime_signature_names(export):
    signatures = []
    for signature in export["signatures"]:
        names = re.findall(r"\b([A-Za-z_]\w*)(?=:)", signature["definition"])
        signatures.append([name for name in names if name not in {"self", "cls"}])
    return signatures


def runtime_parameter_records(signature, symbol):
    definition = signature["definition"]
    for index in range(len(signature["defaults"]) - 1, -1, -1):
        definition = definition.replace(f"\\{index}", repr(signature["defaults"][index]))
    if definition.startswith("("):
        definition = f"def {symbol}{definition}"
    try:
        node = ast.parse(definition + ":\n    pass").body[0]
    except SyntaxError:
        return []
    return parameter_records(node)


def build_contract(root: Path):
    manifest = json.loads((root / "api" / "public_api.json").read_text(encoding="utf-8"))
    facade_path = root / "src" / "lsdtt3" / "_api_docs.py"
    tree = ast.parse(facade_path.read_text(encoding="utf-8"), filename=str(facade_path))
    facade = {
        node.name: node for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
        and not node.name.startswith("_")
    }
    exports = {entry["name"]: entry for entry in manifest["exports"]}
    records = []
    for name in sorted(set(facade).intersection(exports)):
        facade_node = signature_node(facade[name])
        facade_records = parameter_records(facade[name])
        facade_names = [record["name"] for record in facade_records]
        facade_definition = ast.unparse(facade_node.args) if facade_node is not None else ""
        runtime_names = runtime_signature_names(exports[name])
        runtime_records = [
            runtime_parameter_records(signature, name)
            for signature in exports[name]["signatures"]
        ]
        records.append(
            {
                "symbol": name,
                "facade_parameters": facade_names,
                "facade_definition": facade_definition,
                "facade_parameter_records": facade_records,
                "runtime_overloads": runtime_names,
                "runtime_signatures": exports[name]["signatures"],
                "matches_runtime_overload": facade_names in runtime_names,
                "matches_runtime_signature": facade_records in runtime_records,
            }
        )
    return {"schema_version": 1, "symbols": records}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--update", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = root / "api" / "documentation_signatures.json"
    exceptions_output = root / "api" / "documentation_signature_exceptions.json"
    contract = build_contract(root)
    exceptions = {
        "schema_version": 1,
        "reason": "Explicit Phase 3 facade-signature backlog",
        "symbols": [
            record["symbol"] for record in contract["symbols"]
            if not record["matches_runtime_signature"]
        ],
    }
    rendered = json.dumps(contract, indent=2) + "\n"
    exceptions_rendered = json.dumps(exceptions, indent=2) + "\n"
    if args.update:
        output.write_text(rendered, encoding="utf-8")
        exceptions_output.write_text(exceptions_rendered, encoding="utf-8")
        print(f"updated {output} and {exceptions_output}")
        return 0
    if output.read_text(encoding="utf-8") != rendered:
        raise RuntimeError("Python documentation signature contract is stale")
    if exceptions_output.read_text(encoding="utf-8") != exceptions_rendered:
        raise RuntimeError("Python documentation signature exceptions are stale")
    print("Python documentation signature contract is current")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
