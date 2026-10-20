import configparser
from pathlib import Path

import pytest

from lsdtt_qgis.backend.manifest import (
    backend_compatibility_error,
    current_plugin_version,
    manifest_supports_plugin,
)


def test_manifest_requires_range():
    with pytest.raises(ValueError, match="no plugin_version_range"):
        manifest_supports_plugin({}, "0.1.0")


def test_manifest_respects_version_range_tokens():
    manifest = {"plugin_version_range": ">=0.1.0 <0.3.0"}
    assert manifest_supports_plugin(manifest, "0.1.0")
    assert manifest_supports_plugin(manifest, "0.2.5")
    assert not manifest_supports_plugin(manifest, "1.0.0")
    assert not manifest_supports_plugin(manifest, "0.0.9")


@pytest.mark.parametrize(
    "version_range",
    ["0.1.0", ">=", ">=development", ">=0.1.0,", ">=0.2.0 nope", 12],
)
def test_manifest_rejects_malformed_version_ranges(version_range):
    with pytest.raises(ValueError):
        manifest_supports_plugin({"plugin_version_range": version_range}, "0.1.0")


def test_current_plugin_version_comes_from_metadata():
    parser = configparser.ConfigParser()
    parser.read(Path(__file__).parents[1] / "metadata.txt", encoding="utf-8")
    assert current_plugin_version() == parser.get("general", "version")


def _bundle(version):
    """A manifest as the packaging scripts write it for backend `version`."""
    major, minor, _ = (int(part) for part in version.split("."))
    return {
        "backend_version": version,
        "plugin_version_range": ">={}.{}.0 <{}.{}.0".format(major, minor, major, minor + 1),
    }


def test_missing_manifest_or_version_is_refused():
    assert "no manifest.json" in backend_compatibility_error(None, "0.7.0")
    assert "Install the lsdtt3 0.7 backend bundle" in backend_compatibility_error(
        None, "0.7.0"
    )
    assert "not a JSON object" in backend_compatibility_error([], "0.7.0")
    assert "no backend_version" in backend_compatibility_error({}, "0.7.0")


def test_same_major_minor_is_accepted_at_any_patch_level():
    assert backend_compatibility_error(_bundle("0.7.0"), "0.7.0") is None
    assert backend_compatibility_error(_bundle("0.7.3"), "0.7.0") is None
    assert backend_compatibility_error(_bundle("0.7.0"), "0.7.5") is None
    assert backend_compatibility_error(_bundle("1.2.0"), "1.2.9") is None


def test_older_backend_asks_for_the_matching_backend():
    error = backend_compatibility_error(_bundle("0.6.0"), "0.7.0")
    assert "backend is version 0.6.0" in error
    assert "Install the lsdtt3 0.7 backend bundle" in error


def test_newer_backend_asks_for_a_plugin_upgrade():
    error = backend_compatibility_error(_bundle("0.8.1"), "0.7.0")
    assert "newer than this plugin" in error
    assert "Upgrade the plugin to 0.8.x" in error
    assert "Upgrade the plugin to 1.0.x" in backend_compatibility_error(
        _bundle("1.0.0"), "0.7.0"
    )


def test_installed_plugin_accepts_its_own_series():
    version = current_plugin_version()
    assert backend_compatibility_error(_bundle(version)) is None


def test_current_wrappers_reject_invalid_versions():
    assert "invalid backend_version" in backend_compatibility_error(
        {"backend_version": "development"}, "0.7.0"
    )
    assert "invalid version" in backend_compatibility_error(
        _bundle("0.7.0"), "0.7.0-rc1"
    )


def test_current_wrappers_reject_missing_or_malformed_plugin_range():
    assert "invalid plugin_version_range" in backend_compatibility_error(
        {"backend_version": "0.7.0"}, "0.7.0"
    )
    assert "invalid plugin_version_range" in backend_compatibility_error(
        {"backend_version": "0.7.0", "plugin_version_range": ">=banana"},
        "0.7.0",
    )


def test_current_wrappers_reject_plugin_outside_bundle_range():
    # Same series, but a bundle packaged with a narrower --plugin-version-range.
    error = backend_compatibility_error(
        {"backend_version": "0.7.2", "plugin_version_range": ">=0.7.2 <0.8.0"},
        "0.7.0",
    )
    assert "supports plugin versions" in error
    assert "current plugin version is 0.7.0" in error


# --- every backend CLI the plugin uses -------------------------------------

import re
from pathlib import Path

from lsdtt_qgis.backend.health import (
    REQUIRED_CLIS,
    backend_cli_error,
    missing_clis,
)

PLUGIN_DIR = Path(__file__).resolve().parents[1]
CLI_PROFILES = PLUGIN_DIR.parents[1] / "tools" / "cli_profiles.txt"


class _FakeManager:
    """Stand-in for BackendManager: a backend holding only `present` CLIs."""

    def __init__(self, present, path="/backend"):
        self._present = set(present)
        self._path = path

    def backend_path(self):
        return self._path

    def executable_path(self, cli_name):
        if cli_name not in self._present:
            raise RuntimeError("Executable not found for '{}'".format(cli_name))
        return "{}/bin/{}".format(self._path, cli_name)


def test_complete_backend_has_no_missing_clis():
    manager = _FakeManager(REQUIRED_CLIS)
    assert missing_clis(manager) == []
    assert backend_cli_error(manager) is None


def test_backend_cli_error_lists_every_missing_cli():
    present = [cli for cli in REQUIRED_CLIS
               if cli not in ("lsdtt-swath-profile", "lsdtt-basin-stats")]
    error = backend_cli_error(_FakeManager(present))
    assert "missing 2 of the {} CLIs".format(len(REQUIRED_CLIS)) in error
    assert "lsdtt-swath-profile" in error
    assert "lsdtt-basin-stats" in error
    assert "/backend" in error


def test_backend_cli_error_without_backend_path():
    assert "not configured" in backend_cli_error(_FakeManager([], path=""))


def test_required_clis_match_the_clis_the_algorithms_run():
    # Any quoted "lsdtt-..." name in an algorithm module is a CLI it runs or
    # probes for help. A new wrapper must add its CLI to REQUIRED_CLIS.
    used = set()
    for path in (PLUGIN_DIR / "algorithms").glob("*.py"):
        used.update(re.findall(r'"(lsdtt-[a-z-]+)"', path.read_text(encoding="utf-8")))
    assert used == set(REQUIRED_CLIS)
    assert len(REQUIRED_CLIS) == len(set(REQUIRED_CLIS))


def test_required_clis_match_the_qgis_bundle_profile():
    # The qgis backend bundle must ship exactly what the plugin needs.
    profile = set()
    for line in CLI_PROFILES.read_text(encoding="utf-8").splitlines():
        fields = line.split("#", 1)[0].split()
        if len(fields) == 2 and fields[0] == "qgis":
            profile.add(fields[1])
    assert profile == set(REQUIRED_CLIS)
