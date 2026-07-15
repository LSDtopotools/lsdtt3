#!/usr/bin/env python3
import argparse
import configparser
from pathlib import Path
import zipfile


EXCLUDED_DIRS = {
    "__pycache__",
    ".pytest_cache",
}

EXCLUDED_SUFFIXES = {
    ".pyc",
    ".pyo",
}


def read_plugin_metadata(metadata_path: Path):
    parser = configparser.ConfigParser()
    with metadata_path.open("r", encoding="utf-8") as src:
        parser.read_file(src)

    if not parser.has_section("general"):
        raise RuntimeError("metadata.txt is missing [general] section")

    name = parser.get("general", "name", fallback="lsdtt_qgis").strip()
    version = parser.get("general", "version", fallback="0.0.0").strip()
    return name, version


def should_include(path: Path, plugin_dir: Path, include_tests: bool):
    relative = path.relative_to(plugin_dir)

    for part in relative.parts:
        if part in EXCLUDED_DIRS:
            return False

    if path.suffix in EXCLUDED_SUFFIXES:
        return False

    if not include_tests and relative.parts and relative.parts[0] == "tests":
        return False

    return path.is_file()


def collect_files(plugin_dir: Path, include_tests: bool):
    files = []
    for path in plugin_dir.rglob("*"):
        if should_include(path, plugin_dir, include_tests):
            files.append(path)
    return sorted(files)


def main():
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Package the lsdtt3 QGIS plugin into an installable zip.",
    )
    parser.add_argument(
        "--plugin-dir",
        default=str(script_dir / "lsdtt_qgis"),
        help="Path to plugin folder (default: qgis_plugin/lsdtt_qgis)",
    )
    parser.add_argument(
        "--out-dir",
        default=str(script_dir / "dist"),
        help="Output directory for zip artifact (default: qgis_plugin/dist)",
    )
    parser.add_argument(
        "--include-tests",
        action="store_true",
        help="Include tests/ files in the packaged plugin",
    )
    args = parser.parse_args()

    plugin_dir = Path(args.plugin_dir).resolve()
    out_dir = Path(args.out_dir).resolve()

    metadata_path = plugin_dir / "metadata.txt"
    if not plugin_dir.is_dir():
        raise RuntimeError("Plugin directory not found: {}".format(plugin_dir))
    if not metadata_path.exists():
        raise RuntimeError("metadata.txt not found in plugin directory: {}".format(metadata_path))

    plugin_name, version = read_plugin_metadata(metadata_path)
    files = collect_files(plugin_dir, include_tests=args.include_tests)
    if not files:
        raise RuntimeError("No files found to package under {}".format(plugin_dir))

    out_dir.mkdir(parents=True, exist_ok=True)
    zip_name = "{}-{}.zip".format(plugin_name, version)
    zip_path = out_dir / zip_name

    if zip_path.exists():
        zip_path.unlink()

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for file_path in files:
            archive_name = plugin_dir.name + "/" + str(file_path.relative_to(plugin_dir)).replace("\\", "/")
            zf.write(file_path, archive_name)

    print("Created plugin zip: {}".format(zip_path))
    print("Packaged files: {}".format(len(files)))


if __name__ == "__main__":
    main()
