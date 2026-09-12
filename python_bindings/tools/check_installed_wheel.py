#!/usr/bin/env python3

"""Validate an installed lsdtt3 wheel's API and bundled console drivers."""

import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess

import lsdtt3


def load_manifest_generator(root: Path):
    path = root / "tools" / "generate_public_api_manifest.py"
    spec = importlib.util.spec_from_file_location("public_api_manifest", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=root / "api" / "public_api.json")
    parser.add_argument("--smoke-console-scripts", action="store_true")
    args = parser.parse_args()

    expected = json.loads(args.manifest.read_text(encoding="utf-8"))
    actual = load_manifest_generator(root).build_manifest()
    if actual != expected:
        raise RuntimeError("installed wheel public API differs from the manifest")

    package_bin = Path(lsdtt3.__file__).resolve().parent / "bin"
    decisions = json.loads((root / "api" / "decisions.json").read_text(encoding="utf-8"))
    conditional_drivers = {
        feature["conditional_driver"]: feature["name"]
        for feature in decisions["optional_features"]
        if "conditional_driver" in feature
    }
    missing = []
    for name in expected["declared_drivers"]:
        executable_name = name + (".exe" if __import__("os").name == "nt" else "")
        executable = package_bin / executable_name
        if not executable.is_file():
            if name in conditional_drivers:
                feature = conditional_drivers[name]
                if feature == "fftw" and not lsdtt3.spectral_available():
                    continue
            missing.append(name)
            continue
        if args.smoke_console_scripts:
            script = shutil.which(name)
            if script is None:
                raise RuntimeError(f"console script is not installed: {name}")
            result = subprocess.run(
                [script, "-v"], capture_output=True, text=True, check=False
            )
            if result.returncode != 0:
                raise RuntimeError(f"console script failed: {name}\n{result.stderr}")
    if missing:
        raise RuntimeError(f"wheel lacks declared native drivers: {missing}")
    print(f"installed wheel API and {len(expected['declared_drivers'])} drivers are current")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
