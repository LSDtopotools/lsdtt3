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
        [os.fspath(executable), *sys.argv[1:]], check=False
    ).returncode


def raster_preprocessing() -> int:
    """Run the bundled ``lsdtt-raster-preprocessing`` driver.

    Authors
    -------
    B. Gailleton
    """
    _run("lsdtt-raster-preprocessing")


def channel_extraction() -> int:
    """Run the bundled ``lsdtt-channel-extraction`` driver.

    Authors
    -------
    B. Gailleton
    """
    _run("lsdtt-channel-extraction")


def chi_analysis() -> int:
    """Run the bundled ``lsdtt-chi-analysis`` driver.

    Authors
    -------
    B. Gailleton
    """
    _run("lsdtt-chi-analysis")


def graphflood() -> int:
    """Run the bundled ``lsdtt-graphflood`` driver.

    Authors
    -------
    B. Gailleton
    """
    _run("lsdtt-graphflood")


def fetch_raster() -> int:
    """Run the bundled ``lsdtt-fetch-raster`` driver.

    Authors
    -------
    B. Gailleton
    """
    _run("lsdtt-fetch-raster")
