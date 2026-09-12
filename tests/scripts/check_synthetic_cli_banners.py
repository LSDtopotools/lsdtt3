#!/usr/bin/env python3

import argparse
import pathlib
import subprocess


EXPECTED = {
    "surfaces": (
        "one-shot surfaces and ground-truth fixtures",
        "lsdtt-synthetic-networks",
        "lsdtt-lem",
    ),
    "networks": (
        "optimised drainage networks and chi-based topography",
        "lsdtt-synthetic-surfaces",
        "lsdtt-lem",
    ),
    "lem": (
        "Evolves loaded or synthetic topography through geomorphic processes",
        "lsdtt-synthetic-surfaces",
        "lsdtt-synthetic-networks",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check the three synthetic CLI startup guidance banners"
    )
    parser.add_argument("--surfaces", required=True, type=pathlib.Path)
    parser.add_argument("--networks", required=True, type=pathlib.Path)
    parser.add_argument("--lem", required=True, type=pathlib.Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for role in EXPECTED:
        executable = getattr(args, role)
        result = subprocess.run(
            [str(executable), "-h"], capture_output=True, text=True, check=False
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"{executable} -h failed ({result.returncode})\n{result.stderr}"
            )
        missing = [text for text in EXPECTED[role] if text not in result.stdout]
        if missing:
            raise RuntimeError(f"{role} banner lacks required text: {missing}")
    print("synthetic CLI startup guidance is complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
