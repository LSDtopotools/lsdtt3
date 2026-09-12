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
    assert current_plugin_version() == "0.2.0"


def test_current_wrappers_require_versioned_phase_4_backend():
    assert MINIMUM_BACKEND_VERSION == "0.5.0"
    assert "no manifest.json" in backend_compatibility_error(None)
    assert "not a JSON object" in backend_compatibility_error([])
    assert "no backend_version" in backend_compatibility_error({})
    assert "version 0.5.0 or newer" in backend_compatibility_error(
        {"backend_version": "0.4.9"}
    )
    compatible = {
        "backend_version": "0.5.0",
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
        {"backend_version": "0.5.0"}, "0.1.0"
    )
    assert "invalid plugin_version_range" in backend_compatibility_error(
        {"backend_version": "0.5.0", "plugin_version_range": ">=banana"},
        "0.1.0",
    )


def test_current_wrappers_reject_incompatible_plugin_version():
    error = backend_compatibility_error(
        {"backend_version": "0.5.0", "plugin_version_range": ">=0.2.0 <0.3.0"},
        "0.1.0",
    )
    assert "supports plugin versions" in error
    assert "current plugin version is 0.1.0" in error
