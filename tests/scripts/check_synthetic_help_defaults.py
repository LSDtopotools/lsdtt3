#!/usr/bin/env python3

import argparse
import json
import math
import pathlib
import subprocess
import tempfile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check synthetic CLI defaults against its JSON help registry"
    )
    parser.add_argument("--exe", required=True, type=pathlib.Path)
    return parser.parse_args()


def run_checked(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def json_payload(stdout: str) -> dict:
    start = stdout.find("{")
    if start < 0:
        raise RuntimeError("JSON help output contains no object")
    return json.loads(stdout[start:])


def ingested_values(path: pathlib.Path) -> dict[str, str]:
    values = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line and not line.startswith("#") and ":" in line:
            key, value = line.split(":", 1)
            values[key.strip()] = value.strip()
    return values


def defaults_equal(parameter_type: str, expected: str, actual: str) -> bool:
    if parameter_type == "bool":
        return expected.lower() == actual.lower()
    if parameter_type == "int":
        return int(expected) == int(actual)
    if parameter_type == "float":
        return math.isclose(float(expected), float(actual), rel_tol=1e-6, abs_tol=1e-9)
    return expected == actual


def main() -> int:
    args = parse_args()
    exe = args.exe.resolve()

    help_data = json_payload(run_checked([str(exe), "-h", "json"]).stdout)
    help_names = [entry["name"] for entry in help_data["parameters"]]
    if len(help_names) != len(set(help_names)):
        raise RuntimeError("JSON help contains duplicate parameter names")

    with tempfile.TemporaryDirectory(prefix="lsdtt-synthetic-help-") as temp_name:
        temp_dir = pathlib.Path(temp_name)
        param_path = temp_dir / "synthetic_help.param"
        param_path.write_text(
            f"write_path: {temp_dir.as_posix()}/\n"
            "write_fname: synthetic_help\n",
            encoding="utf-8",
        )
        run_checked([str(exe), str(temp_dir), param_path.name])
        defaults = ingested_values(
            temp_dir / "synthetic_help_ingested_params.param"
        )

    help_set = set(help_names)
    defaults_set = set(defaults)
    if help_set != defaults_set:
        missing = sorted(defaults_set - help_set)
        extra = sorted(help_set - defaults_set)
        raise RuntimeError(
            "synthetic defaults/help registry mismatch\n"
            f"missing from help: {missing}\nextra in help: {extra}"
        )

    parser_io = {"read_path", "read_fname", "write_path", "write_fname"}
    mismatches = []
    for parameter in help_data["parameters"]:
        if parameter["name"] in parser_io:
            continue
        if not defaults_equal(
            parameter["type"], parameter["default"], defaults[parameter["name"]]
        ):
            mismatches.append(
                f"{parameter['name']}: help={parameter['default']!r}, "
                f"ingested={defaults[parameter['name']]!r}"
            )
    if mismatches:
        raise RuntimeError("synthetic default-value mismatch\n" + "\n".join(mismatches))

    print(f"synthetic defaults/help registry exact match: {len(help_names)} keys")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
