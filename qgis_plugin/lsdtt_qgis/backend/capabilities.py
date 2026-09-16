"""Backend capability probing via the lsdtt --help-json interface.

Backends from 2026-07-15 onwards report their compiled optional
dependencies (fftw, eigen, netcdf, gpu_kernels) and full parameter
metadata through `<tool> --help-json`. Wrappers use this to fail fast with
a clear message when a requested option needs a feature the backend lacks.

Older backends do not support --help-json; probing then returns None and
callers must treat capabilities as unknown (i.e. not block the run).
"""

from ..naming import extract_json_payload
from .manager import BackendManager
from .runner import run_cli

_probe_cache = {}


def probe_tool(tool_name):
    """Return the parsed --help-json dict for a backend tool, or None.

    None means the probe failed (missing executable, old backend without
    --help-json, or unparseable output); capabilities are then unknown.
    Results are cached per session.
    """
    if tool_name in _probe_cache:
        return _probe_cache[tool_name]

    payload = None
    try:
        backend = BackendManager()
        executable = backend.executable_path(tool_name)
        result = run_cli(
            [executable, "--help-json"],
            timeout=30,
            env=backend.runtime_environment(),
        )
        if result.returncode == 0:
            payload = extract_json_payload(result.stdout or "")
    except Exception:
        payload = None

    _probe_cache[tool_name] = payload
    return payload


def backend_features(tool_name):
    """Return the backend tool's compiled feature set, or None if unknown."""
    payload = probe_tool(tool_name)
    if payload is None:
        return None
    features = payload.get("features")
    if not isinstance(features, list):
        return None
    return set(features)


def require_feature(tool_name, feature, option_description, feedback=None):
    """Check that the backend tool was built with a feature.

    Returns an error message string when the backend definitely lacks the
    feature, or None when it has it or capabilities are unknown (old
    backend). Callers raise QgsProcessingException on a non-None result.
    """
    features = backend_features(tool_name)
    if features is None:
        if feedback is not None:
            feedback.pushInfo(
                "Backend capabilities unknown ({} does not support "
                "--help-json); proceeding without checking for {}.".format(
                    tool_name, feature
                )
            )
        return None
    if feature in features:
        return None
    return (
        "{} requires a backend built with {} support, but the configured "
        "backend reports features {}. Install a backend bundle that "
        "includes {}.".format(
            option_description, feature,
            sorted(features) if features else "none", feature,
        )
    )


def clear_cache():
    """Reset cached probes (used after the backend path changes)."""
    _probe_cache.clear()
