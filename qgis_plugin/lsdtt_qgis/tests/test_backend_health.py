import pytest

from lsdtt_qgis.backend.manifest import (
    MINIMUM_BACKEND_VERSION,
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
    assert current_plugin_version() == "0.4.0"


def test_current_wrappers_require_versioned_phase_4_backend():
    assert MINIMUM_BACKEND_VERSION == "0.6.0"
    assert "no manifest.json" in backend_compatibility_error(None)
    assert "not a JSON object" in backend_compatibility_error([])
    assert "no backend_version" in backend_compatibility_error({})
    assert "version 0.6.0 or newer" in backend_compatibility_error(
        {"backend_version": "0.5.9"}
    )
    compatible = {
        "backend_version": "0.6.0",
        "plugin_version_range": ">=0.1.0 <0.3.0",
    }
    assert backend_compatibility_error(compatible, "0.1.0") is None
    compatible["backend_version"] = "1.0.0"
    assert backend_compatibility_error(compatible, "0.1.0") is None


def test_current_wrappers_reject_invalid_backend_version():
    assert "invalid backend_version" in backend_compatibility_error(
        {"backend_version": "development"}
    )


def test_current_wrappers_reject_missing_or_malformed_plugin_range():
    assert "invalid plugin_version_range" in backend_compatibility_error(
        {"backend_version": "0.6.0"}, "0.1.0"
    )
    assert "invalid plugin_version_range" in backend_compatibility_error(
        {"backend_version": "0.6.0", "plugin_version_range": ">=banana"},
        "0.1.0",
    )


def test_current_wrappers_reject_incompatible_plugin_version():
    error = backend_compatibility_error(
        {"backend_version": "0.6.0", "plugin_version_range": ">=0.2.0 <0.3.0"},
        "0.1.0",
    )
    assert "supports plugin versions" in error
    assert "current plugin version is 0.1.0" in error


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
