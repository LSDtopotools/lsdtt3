from lsdtt_qgis.backend.manifest import manifest_supports_plugin


def test_manifest_accepts_missing_range():
    assert manifest_supports_plugin({}, "0.1.0")


def test_manifest_respects_version_range_tokens():
    manifest = {"plugin_version_range": ">=0.1.0 <0.3.0"}
    assert manifest_supports_plugin(manifest, "0.1.0")
    assert manifest_supports_plugin(manifest, "0.2.5")
    assert not manifest_supports_plugin(manifest, "1.0.0")
    assert not manifest_supports_plugin(manifest, "0.0.9")
