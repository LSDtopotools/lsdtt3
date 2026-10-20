import configparser
import io
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import zipfile


PLUGIN_DIR = Path(__file__).resolve().parents[1]
QGIS_PLUGIN_DIR = PLUGIN_DIR.parent

PHASE_4_ALGORITHMS = (
    "basin_extraction.py",
    "channel_extraction.py",
    "channel_network_analysis.py",
    "chi_analysis.py",
    "side_tributary_basins.py",
)

EXCLUDED_DIRS = {"tests", "__pycache__", ".pytest_cache"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo"}


def _metadata(contents):
    parser = configparser.ConfigParser()
    parser.read_file(io.StringIO(contents))
    return parser["general"]


def _expected_archive_names():
    names = []
    for path in PLUGIN_DIR.rglob("*"):
        relative = path.relative_to(PLUGIN_DIR)
        if (
            path.is_file()
            and not EXCLUDED_DIRS.intersection(relative.parts)
            and path.suffix not in EXCLUDED_SUFFIXES
        ):
            names.append("{}/{}".format(PLUGIN_DIR.name, relative.as_posix()))
    return sorted(names)


def test_distributable_matches_plugin_source_contract(tmp_path):
    source_metadata = _metadata(
        (PLUGIN_DIR / "metadata.txt").read_text(encoding="utf-8")
    )
    assert source_metadata["name"] == "lsdtt3"
    # One version number for the repository: the plugin carries the CMake one.
    cmake_version = re.search(
        r"^project\(lsdtt3\s+VERSION\s+(\d+\.\d+\.\d+)",
        (PLUGIN_DIR.parents[1] / "CMakeLists.txt").read_text(encoding="utf-8"),
        re.MULTILINE,
    ).group(1)
    assert source_metadata["version"] == cmake_version
    assert source_metadata["hasProcessingProvider"] == "yes"

    subprocess.run(
        [
            sys.executable,
            str(QGIS_PLUGIN_DIR / "package_plugin.py"),
            "--out-dir",
            str(tmp_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    expected_zip = tmp_path / "lsdtt3-qgis-plugin-{}.zip".format(
        source_metadata["version"]
    )
    assert sorted(tmp_path.glob("*.zip")) == [expected_zip]

    with zipfile.ZipFile(expected_zip) as archive:
        names = archive.namelist()
        assert names == _expected_archive_names()
        assert {PurePosixPath(name).parts[0] for name in names} == {PLUGIN_DIR.name}

        packaged_metadata = _metadata(
            archive.read(PLUGIN_DIR.name + "/metadata.txt").decode("utf-8")
        )
        assert packaged_metadata["name"] == source_metadata["name"]
        assert packaged_metadata["version"] == source_metadata["version"]

        for name in names:
            path = PurePosixPath(name)
            assert not EXCLUDED_DIRS.intersection(path.parts)
            assert path.suffix not in EXCLUDED_SUFFIXES

        for filename in PHASE_4_ALGORITHMS:
            contents = archive.read(
                "{}/algorithms/{}".format(PLUGIN_DIR.name, filename)
            ).decode("utf-8")
            assert 'THRESHOLD_CONTRIBUTING_PIXELS = "THRESHOLD_CONTRIBUTING_PIXELS"' in contents
            assert '"threshold_contributing_pixels"' in contents
            assert "SOURCE_THRESHOLD" not in contents
            assert "source_pixels_threshold" not in contents

        surface_metrics = archive.read(
            PLUGIN_DIR.name + "/algorithms/surface_metrics.py"
        ).decode("utf-8")
        assert 'self.SLOPE_THRESHOLD_LAYERS: "slope_thresholds"' in surface_metrics
        assert '"slope_thresholds": self.parameterAsString(' in surface_metrics

        manifest = archive.read(
            PLUGIN_DIR.name + "/backend/manifest.py"
        ).decode("utf-8")
        assert "def _series(version):" in manifest
        assert "Upgrade the plugin to" in manifest
