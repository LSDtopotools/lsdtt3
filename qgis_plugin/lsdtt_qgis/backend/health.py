from .runner import run_cli


# Every backend CLI an algorithm in this plugin runs. Keep in step with the
# `qgis` profile in tools/cli_profiles.txt; tests/test_backend_health.py checks
# both that list and the CLIs named in algorithms/*.py.
REQUIRED_CLIS = [
    "lsdtt-raster-preprocessing",
    "lsdtt-channel-extraction",
    "lsdtt-chi-analysis",
    "lsdtt-basin-extraction",
    "lsdtt-basin-stats",
    "lsdtt-surface-metrics",
    "lsdtt-channel-network-analysis",
    "lsdtt-swath-profile",
]


def missing_clis(backend_manager):
    """Return the REQUIRED_CLIS that the configured backend does not contain."""
    missing = []
    for cli in REQUIRED_CLIS:
        try:
            backend_manager.executable_path(cli)
        except RuntimeError:
            missing.append(cli)
    return missing


def backend_cli_error(backend_manager):
    """Return why the configured backend cannot run every algorithm, or None.

    All missing CLIs are reported at once, so an old or partial bundle is
    flagged when the backend is saved in the settings dialog. Running an
    algorithm checks only that algorithm's CLI (support.run_lsdtt_tool).
    """
    if not backend_manager.backend_path():
        return "Backend path is not configured. Set it in plugin settings."
    missing = missing_clis(backend_manager)
    if not missing:
        return None
    return (
        "The configured backend ({}) is missing {} of the {} CLIs this plugin "
        "uses: {}. Install the qgis backend bundle that matches this plugin."
    ).format(
        backend_manager.backend_path(), len(missing), len(REQUIRED_CLIS),
        ", ".join(missing),
    )


def fast_health_check(backend_manager=None, feedback=None):
    if backend_manager is None:
        from .manager import BackendManager
        backend_manager = BackendManager()
    error = backend_cli_error(backend_manager)
    if error:
        raise RuntimeError(error)

    version_cli = backend_manager.executable_path("lsdtt-raster-preprocessing")
    result = run_cli(
        [version_cli, "--version"],
        feedback=feedback,
        timeout=20,
        env=backend_manager.runtime_environment(),
    )
    return result.returncode == 0
