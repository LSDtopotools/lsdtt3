"""Launch bundled vanilla LSDTopoTools command-line drivers.

The functions in this module are installed as normal ``pip`` console scripts.
They forward every command-line argument directly to the matching native
driver, so its parameter-file interface and return code remain unchanged.

Authors
-------
B. Gailleton
"""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
from typing import Mapping, Optional


_DRIVERS = frozenset(
    {
        "lsdtt-raster-preprocessing", "lsdtt-synthetic-surfaces",
        "lsdtt-synthetic-networks", "lsdtt-surface-metrics", "lsdtt-slope-area",
        "lsdtt-channel-extraction", "lsdtt-internal-drainage",
        "lsdtt-hillslope-channel-coupling", "lsdtt-drainage-analysis",
        "lsdtt-distributed-hydrology", "lsdtt-feature-finder", "lsdtt-mountain-front", "lsdtt-outlet-matching",
        "lsdtt-basin-stats",
        "lsdtt-basin-extraction", "lsdtt-flat-surface-finder",
        "lsdtt-fluvial-surfaces",
        "lsdtt-valley-centreline", "lsdtt-graphflood",
        "lsdtt-channel-network-analysis", "lsdtt-cross-sections",
        "lsdtt-stream-burning", "lsdtt-water-line-repair",
        "lsdtt-water-centreline", "lsdtt-cosmo-tool",
        "lsdtt-cosmo-prediction", "lsdtt-cosmo-landslides",
        "lsdtt-cosmo-column-tool", "lsdtt-marsh-topography",
        "lsdtt-linear-stats",
        "lsdtt-swath-profile", "lsdtt-lem", "lsdtt-lithology-erodibility",
        "lsdtt-advection-diffusion",
        "lsdtt-sediment-routing",
        "lsdtt-ridgeline-analysis", "lsdtt-basin-surface-metrics",
        "lsdtt-divide-migration", "lsdtt-raster-calculator",
        "lsdtt-fetch-raster", "lsdtt-fetch-lithology", "lsdtt-fetch-erosion-rates",
        "lsdtt-fetch-river-lines",
        "lsdtt-bankfull-width",
        "lsdtt-spectral-analysis",
        "lsdtt-chi-analysis",
    }
)


def _driver_path(name: str) -> Path:
    """Return the packaged executable path for one driver.

    Parameters
    ----------
    name : str
        Native LSDTopoTools driver name, without a platform suffix.

    Returns
    -------
    pathlib.Path
        Absolute path to the packaged executable.

    Raises
    ------
    RuntimeError
        If this wheel was built without the requested driver.

    Authors
    -------
    B. Gailleton
    """
    suffix = ".exe" if os.name == "nt" else ""
    path = Path(__file__).resolve().parent / "bin" / f"{name}{suffix}"
    if not path.is_file():
        raise RuntimeError(
            f"This lsdtt3 wheel does not include '{name}'. "
            "Install a wheel built with LSDTT_PYTHON_CLI=ON."
        )
    return path


def _child_environment(base: Optional[Mapping[str, str]] = None) -> Optional[dict]:
    """Return the environment a bundled driver must be started with.

    On Windows a driver is a separate ``.exe``. delvewheel stages the wheel's
    DLLs beside ``lsdtt3`` or in a sibling ``lsdtt3.libs`` directory, and a child
    executable does not inherit Python's extension loader hooks, so those
    directories are put at the front of ``PATH``; without them the driver dies
    at start-up with STATUS_DLL_NOT_FOUND (0xC0000135). Elsewhere the drivers
    find their libraries through their rpath and need nothing added.

    Parameters
    ----------
    base : mapping, optional
        Environment to extend; ``os.environ`` when omitted.

    Returns
    -------
    dict or None
        The environment to pass to ``subprocess``; ``None`` (inherit the
        current one) off Windows when no ``base`` is given.

    Authors
    -------
    B. Gailleton
    """
    if os.name != "nt":
        return None if base is None else dict(base)
    package_dir = Path(__file__).resolve().parent
    library_dirs = (package_dir, package_dir.with_name(f"{package_dir.name}.libs"),
                    package_dir / ".libs")
    environment = dict(os.environ if base is None else base)
    environment["PATH"] = os.pathsep.join(
        [*(str(path) for path in library_dirs if path.is_dir()),
         environment.get("PATH", "")]
    )
    return environment


def _run(name: str) -> int:
    """Run a bundled native driver and return its exit status.

    Parameters
    ----------
    name : str
        Native LSDTopoTools driver name.

    Returns
    -------
    int
        Native driver's exit status.

    Authors
    -------
    B. Gailleton
    """
    executable = _driver_path(name)
    return subprocess.run(
        [os.fspath(executable), *sys.argv[1:]], check=False, env=_child_environment()
    ).returncode


def main() -> int:
    """Dispatch a pip console script to its bundled vanilla driver.

    Parameters
    ----------
    None

    Returns
    -------
    int
        Native driver's exit status.

    Raises
    ------
    RuntimeError
        If invoked under a console-script name that is not an LSDTopoTools
        driver distributed by this package.

    Authors
    -------
    B. Gailleton
    """
    name = Path(sys.argv[0]).stem
    if name not in _DRIVERS:
        raise RuntimeError(f"Unknown LSDTopoTools console script: '{name}'.")
    return _run(name)


def raster_preprocessing() -> int:
    """Run the bundled ``lsdtt-raster-preprocessing`` driver.

    Parameters
    ----------
    None

    Returns
    -------
    int
        Native driver's exit status.

    Authors
    -------
    B. Gailleton
    """
    return _run("lsdtt-raster-preprocessing")


def channel_extraction() -> int:
    """Run the bundled ``lsdtt-channel-extraction`` driver.

    Parameters
    ----------
    None

    Returns
    -------
    int
        Native driver's exit status.

    Authors
    -------
    B. Gailleton
    """
    return _run("lsdtt-channel-extraction")


def chi_analysis() -> int:
    """Run the bundled ``lsdtt-chi-analysis`` driver.

    Parameters
    ----------
    None

    Returns
    -------
    int
        Native driver's exit status.

    Authors
    -------
    B. Gailleton
    """
    return _run("lsdtt-chi-analysis")


def graphflood() -> int:
    """Run the bundled ``lsdtt-graphflood`` driver.

    Parameters
    ----------
    None

    Returns
    -------
    int
        Native driver's exit status.

    Authors
    -------
    B. Gailleton
    """
    return _run("lsdtt-graphflood")


def fetch_raster() -> int:
    """Run the bundled ``lsdtt-fetch-raster`` driver.

    Parameters
    ----------
    None

    Returns
    -------
    int
        Native driver's exit status.

    Authors
    -------
    B. Gailleton
    """
    return _run("lsdtt-fetch-raster")
