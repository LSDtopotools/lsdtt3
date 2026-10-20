import configparser
import json
import os
import re


# lsdtt3 has one version number for the backend, the plugin and the Python
# package. A plugin X.Y.* accepts a backend X.Y.*: same major and minor, any
# patch level. A minor bump is how a release marks a change one component
# depends on (a parameter rename, a removed output, a new CLI behaviour the
# wrappers rely on), so a mismatch there means the wrappers may emit
# parameters the backend does not understand.
_VERSION_PATTERN = re.compile(r"^[0-9]+(?:\.[0-9]+){0,2}$")
_RANGE_OPERATORS = (">=", "<=", "==", ">", "<")


def load_manifest(backend_dir):
    manifest_path = os.path.join(backend_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        return None
    with open(manifest_path, "r", encoding="utf-8") as src:
        return json.load(src)


def current_plugin_version():
    """Read the installed plugin version from its QGIS metadata."""
    metadata_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)), "metadata.txt"
    )
    parser = configparser.ConfigParser()
    with open(metadata_path, "r", encoding="utf-8") as src:
        parser.read_file(src)
    if not parser.has_option("general", "version"):
        raise ValueError("metadata.txt has no [general] version")
    return parser.get("general", "version").strip()


def manifest_supports_plugin(manifest_data, plugin_version):
    """Return whether plugin_version satisfies a required manifest range.

    Missing and malformed ranges are invalid rather than unconstrained: every
    supported bundle must state which plugin releases it supports.
    """
    if not isinstance(manifest_data, dict):
        raise ValueError("manifest data is not an object")
    version_range = manifest_data.get("plugin_version_range")
    if not version_range:
        raise ValueError("manifest has no plugin_version_range")
    if not isinstance(version_range, str):
        raise ValueError("plugin_version_range must be a string")

    version = _parse_version(plugin_version)
    constraints = []
    for token in version_range.split():
        operator = next(
            (item for item in _RANGE_OPERATORS if token.startswith(item)), None
        )
        if operator is None or token == operator:
            raise ValueError("invalid plugin_version_range token: {!r}".format(token))
        boundary = _parse_version(token[len(operator):])
        constraints.append((operator, boundary))

    for operator, boundary in constraints:
        if operator == ">=":
            if version < boundary:
                return False
        elif operator == ">":
            if version <= boundary:
                return False
        elif operator == "<=":
            if version > boundary:
                return False
        elif operator == "<":
            if version >= boundary:
                return False
        elif operator == "==":
            if version != boundary:
                return False
    return True


def backend_compatibility_error(manifest_data, plugin_version=None):
    """Return why a backend cannot run current wrappers, or None.

    The backend must have the plugin's major.minor version, and the bundle's
    own plugin_version_range must include the plugin.
    """
    if plugin_version is None:
        try:
            plugin_version = current_plugin_version()
        except (OSError, configparser.Error, ValueError) as exc:
            return "Could not determine the current plugin version: {}.".format(exc)
    try:
        plugin_series = _series(_parse_version(plugin_version))
    except (TypeError, ValueError):
        return "The plugin has an invalid version: {!r}.".format(plugin_version)

    if manifest_data is None:
        return (
            "The configured backend has no manifest.json, so its version cannot "
            "be verified. Install the lsdtt3 {} backend bundle."
        ).format(plugin_series)
    if not isinstance(manifest_data, dict):
        return "The configured backend manifest is not a JSON object."

    backend_version = manifest_data.get("backend_version")
    if not backend_version:
        return (
            "The configured backend manifest has no backend_version. Install the "
            "lsdtt3 {} backend bundle."
        ).format(plugin_series)
    try:
        backend = _parse_version(backend_version)
    except (TypeError, ValueError):
        return "The configured backend has an invalid backend_version: {!r}.".format(
            backend_version
        )

    plugin = _parse_version(plugin_version)
    if backend[:2] < plugin[:2]:
        return (
            "The configured lsdtt3 backend is version {}, but this plugin is version "
            "{} and needs a {} backend. Install the lsdtt3 {} backend bundle."
        ).format(backend_version, plugin_version, plugin_series, plugin_series)
    if backend[:2] > plugin[:2]:
        return (
            "The configured lsdtt3 backend is version {}, which is newer than this "
            "plugin (version {}). Upgrade the plugin to {}.x."
        ).format(backend_version, plugin_version, _series(backend))

    try:
        supported = manifest_supports_plugin(manifest_data, plugin_version)
    except (TypeError, ValueError) as exc:
        return (
            "The configured backend has an invalid plugin_version_range: {}."
        ).format(exc)
    if not supported:
        return (
            "The configured backend supports plugin versions {!r}, but the current "
            "plugin version is {}. Install a compatible backend bundle."
        ).format(manifest_data.get("plugin_version_range"), plugin_version)
    return None


def _series(version):
    """The major.minor series of a parsed version, as text: (0, 7, 2) -> '0.7'."""
    return "{}.{}".format(version[0], version[1])


def _parse_version(value):
    if not isinstance(value, str) or not _VERSION_PATTERN.fullmatch(value.strip()):
        raise ValueError("invalid version {!r}".format(value))
    parts = value.strip().split(".")
    out = []
    for piece in parts[:3]:
        out.append(int(piece))
    while len(out) < 3:
        out.append(0)
    return tuple(out)
