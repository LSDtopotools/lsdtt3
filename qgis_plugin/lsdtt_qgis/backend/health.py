from .manager import BackendManager
from .runner import run_cli


REQUIRED_CLIS = [
    "lsdtt-raster-preprocessing",
    "lsdtt-channel-extraction",
    "lsdtt-chi-analysis",
    "lsdtt-basin-extraction",
]


def fast_health_check(backend_manager=None, feedback=None):
    manager = backend_manager or BackendManager()
    for cli in REQUIRED_CLIS:
        manager.executable_path(cli)

    version_cli = manager.executable_path("lsdtt-raster-preprocessing")
    result = run_cli(
        [version_cli, "--version"],
        feedback=feedback,
        timeout=20,
        env=manager.runtime_environment(),
    )
    return result.returncode == 0
