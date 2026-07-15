import json
import os


def load_manifest(backend_dir):
    manifest_path = os.path.join(backend_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        return None
    with open(manifest_path, "r", encoding="utf-8") as src:
        return json.load(src)


def manifest_supports_plugin(manifest_data, plugin_version):
    if manifest_data is None:
        return True
    version_range = manifest_data.get("plugin_version_range")
    if not version_range:
        return True

    version = _parse_version(plugin_version)
    for token in version_range.split():
        if token.startswith(">="):
            if version < _parse_version(token[2:]):
                return False
        elif token.startswith(">"):
            if version <= _parse_version(token[1:]):
                return False
        elif token.startswith("<="):
            if version > _parse_version(token[2:]):
                return False
        elif token.startswith("<"):
            if version >= _parse_version(token[1:]):
                return False
        elif token.startswith("=="):
            if version != _parse_version(token[2:]):
                return False
    return True


def _parse_version(value):
    parts = value.strip().split(".")
    out = []
    for piece in parts[:3]:
        out.append(int(piece))
    while len(out) < 3:
        out.append(0)
    return tuple(out)
