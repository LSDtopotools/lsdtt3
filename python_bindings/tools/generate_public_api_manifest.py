#!/usr/bin/env python3

"""Generate or check the deterministic lsdtt3 Python public API manifest."""

from __future__ import annotations

import argparse
import importlib
import inspect
import json
import math
from pathlib import Path
import re
import sys
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--update", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "api" / "public_api.json",
    )
    return parser.parse_args()


def normalize_text(value: str) -> str:
    return value.replace("lsdtt3._core.", "lsdtt3.")


def normalize_default(value: Any) -> Any:
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            return repr(value)
        return value
    if isinstance(value, (tuple, list)):
        return [normalize_default(item) for item in value]
    rendered = repr(value)
    if " at 0x" in rendered:
        raise RuntimeError(f"unstable default representation: {rendered}")
    return normalize_text(rendered)


def first_doc_line(value: Any) -> str:
    doc = inspect.getdoc(value) or ""
    return normalize_text(doc.splitlines()[0]) if doc else ""


def full_doc(value: Any) -> str:
    return normalize_text(inspect.getdoc(value) or "")


def callable_signatures(value: Any) -> list[dict[str, Any]]:
    native = inspect.unwrap(value)
    nb_signatures = getattr(native, "__nb_signature__", None)
    if nb_signatures is not None:
        return [
            {
                "definition": normalize_text(definition),
                "defaults": [normalize_default(item) for item in (defaults or ())],
            }
            for definition, _doc, defaults in nb_signatures
        ]

    try:
        signature = str(inspect.signature(value))
    except (TypeError, ValueError):
        return []
    return [{"definition": normalize_text(signature), "defaults": []}]


def public_class_members(cls: type[Any]) -> list[dict[str, Any]]:
    names = {name for name in dir(cls) if not name.startswith("_")}
    names.add("__init__")
    members = []
    for name in sorted(names):
        try:
            value = inspect.getattr_static(cls, name)
        except AttributeError:
            continue
        if isinstance(value, property):
            members.append(
                {
                    "name": name,
                    "kind": "property",
                    "doc": full_doc(value),
                    "signatures": [],
                }
            )
            continue

        bound_value = getattr(cls, name, value)
        native_signatures = getattr(inspect.unwrap(bound_value), "__nb_signature__", None)
        if name == "__init__" and native_signatures is None and type(value).__name__ == "wrapper_descriptor":
            continue
        signatures = callable_signatures(bound_value)
        if signatures or name == "__init__":
            kind = "static_method" if type(value).__name__ == "nb_func" else "method"
            members.append(
                {
                    "name": name,
                    "kind": kind,
                    "doc": full_doc(bound_value),
                    "signatures": signatures,
                }
            )
    return members


def build_manifest() -> dict[str, Any]:
    package = importlib.import_module("lsdtt3")
    package_root = Path(__file__).resolve().parents[1]
    project = tomllib.loads(
        (package_root / "pyproject.toml").read_text(encoding="utf-8")
    )

    exports = []
    export_names = list(package.__all__)
    if len(export_names) != len(set(export_names)):
        raise RuntimeError("lsdtt3.__all__ contains duplicate names")
    for name in sorted(export_names):
        if not hasattr(package, name):
            raise RuntimeError(f"lsdtt3.__all__ names missing attribute {name!r}")
        value = getattr(package, name)
        if inspect.isclass(value):
            members = public_class_members(value)
            if value.__module__ == "lsdtt3._core":
                class_signatures = next(
                    (member["signatures"] for member in members if member["name"] == "__init__"),
                    [],
                )
            else:
                class_signatures = callable_signatures(value)
            exports.append(
                {
                    "name": name,
                    "kind": "class",
                    "doc": full_doc(value),
                    "members": members,
                    "signatures": class_signatures,
                }
            )
        elif callable(value):
            exports.append(
                {
                    "name": name,
                    "kind": "function",
                    "doc": full_doc(value),
                    "members": [],
                    "signatures": callable_signatures(value),
                }
            )
        else:
            exports.append(
                {
                    "name": name,
                    "kind": "value",
                    "doc": full_doc(value),
                    "members": [],
                    "signatures": [],
                }
            )

    return {
        "schema_version": 1,
        "package": "lsdtt3",
        "signature_source": "nanobind __nb_signature__ v3",
        "export_count": len(exports),
        "exports": exports,
        "console_scripts": sorted(project["project"]["scripts"]),
        "declared_drivers": sorted(importlib.import_module("lsdtt3._cli")._DRIVERS),
    }


def validate_decisions(manifest: dict[str, Any], decisions: dict[str, Any]) -> None:
    required = {
        "schema_version", "api_release_status", "release_evidence",
        "approved_name_differences", "completed_hard_renames",
        "approved_default_differences", "deprecated_names",
        "known_runtime_test_failures", "optional_features",
    }
    if set(decisions) != required:
        raise RuntimeError(
            f"decision manifest fields differ: expected {sorted(required)}, "
            f"got {sorted(decisions)}"
        )
    if decisions["api_release_status"] not in {"unresolved", "unshipped", "shipped"}:
        raise RuntimeError("invalid api_release_status")
    exports = {entry["name"]: entry for entry in manifest["exports"]}
    for decision in decisions["approved_name_differences"]:
        symbol, parameter = decision["python"].split(".", 1)
        if symbol not in exports:
            raise RuntimeError(f"name decision references missing symbol {symbol!r}")
        definitions = "\n".join(
            signature["definition"] for signature in exports[symbol]["signatures"]
        )
        if parameter not in definitions:
            raise RuntimeError(f"name decision references missing parameter {decision['python']!r}")
    for decision in decisions["completed_hard_renames"]:
        symbol = decision["symbol"]
        if symbol not in exports:
            raise RuntimeError(f"hard rename references missing symbol {symbol!r}")
        definitions = "\n".join(
            signature["definition"] for signature in exports[symbol]["signatures"]
        )
        if decision["new"] not in definitions or decision["old"] in definitions:
            raise RuntimeError(f"stale hard-rename decision for {symbol}")
    for decision in decisions["approved_default_differences"]:
        if set(decision) != {"python_contracts", "cli", "status"} or "unresolved" in decision["status"]:
            raise RuntimeError("default decisions must have a characterized status")
        for contract in decision["python_contracts"]:
            entry = exports.get(contract["symbol"])
            if entry is None:
                raise RuntimeError("default decision references a missing symbol")
            signatures = entry["signatures"]
            if "member" in contract:
                signatures = next(
                    (member["signatures"] for member in entry["members"]
                     if member["name"] == contract["member"]),
                    [],
                )
            found = False
            for signature in signatures:
                match = re.search(
                    rf"\b{re.escape(contract['parameter'])}\b[^,)]* = \\(\d+)",
                    signature["definition"],
                )
                if match is None:
                    continue
                actual = signature["defaults"][int(match.group(1))]
                expected_default = contract["default"]
                if isinstance(actual, float) or isinstance(expected_default, float):
                    found = math.isclose(float(actual), float(expected_default), rel_tol=1e-6)
                else:
                    found = actual == expected_default
                if found:
                    break
            if not found:
                raise RuntimeError(
                    f"stale default decision for {contract['symbol']}."
                    f"{contract['parameter']}"
                )
    for decision in decisions["deprecated_names"]:
        expected = {"symbol", "old", "new", "introduced", "earliest_removal"}
        if set(decision) != expected or decision["symbol"] not in exports:
            raise RuntimeError("deprecated-name decision has an invalid or stale schema")
    feature_names = set()
    for feature in decisions["optional_features"]:
        if not {"name", "symbols"}.issubset(feature) or feature["name"] in feature_names:
            raise RuntimeError("optional feature has an invalid schema")
        feature_names.add(feature["name"])
        if not set(feature["symbols"]).issubset(exports):
            raise RuntimeError("optional feature references a missing symbol")
        driver = feature.get("conditional_driver")
        if driver is not None and driver not in manifest["declared_drivers"]:
            raise RuntimeError(f"optional feature references missing driver {driver!r}")
    for failure in decisions["known_runtime_test_failures"]:
        if not failure.get("tests") or not failure.get("reason"):
            raise RuntimeError("known runtime failure lacks tests or reason")
    if decisions["api_release_status"] == "shipped":
        policy = decisions["release_evidence"].get("compatibility_policy", "")
        if "minor release" not in policy:
            raise RuntimeError("shipped API lacks a release-based compatibility policy")


def main() -> int:
    args = parse_args()
    output = args.output.resolve()
    manifest = build_manifest()
    decisions_path = Path(__file__).resolve().parents[1] / "api" / "decisions.json"
    decisions = json.loads(decisions_path.read_text(encoding="utf-8"))
    validate_decisions(manifest, decisions)
    rendered = json.dumps(manifest, indent=2, sort_keys=False) + "\n"
    if args.update:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
        print(f"updated {output}")
        return 0
    if not output.is_file():
        raise RuntimeError(f"manifest does not exist: {output}")
    if output.read_text(encoding="utf-8") != rendered:
        raise RuntimeError(f"public API manifest is stale: {output}")
    print(f"public API manifest is current: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
