"""Shared helpers for lsdtt Processing algorithm wrappers.

The lsdtt CLIs report per-step failures as warnings and still exit 0, so
output-existence checks in the wrappers are the primary error detection.

Pure path/format logic lives in lsdtt_qgis.naming (qgis-free, unit-testable);
this module adds the qgis-dependent wrappers around it.
"""

import datetime
import os
import tempfile

from qgis.core import QgsProcessingException

from ..backend.manager import BackendManager
from ..backend.runner import run_cli, write_log
from ..naming import (
    build_backend_help_index,
    expected_output_path,
    extension_for_format,
    lookup_param_help,
    next_qgis_driver_name,
    split_source,
    strip_layer_suffix,
)
from ..templates.lsdtt_params import write_parameter_file

__all__ = [
    "attach_backend_help",
    "extension_for_format",
    "layer_file_path",
    "require_output",
    "run_lsdtt_tool",
    "split_layer_source",
]


def attach_backend_help(alg, tool_name, overrides=None):
    """Attach backend-derived tooltips to an algorithm's parameters.

    Looks up each QGIS parameter in the cached `<tool> --help-json` payload
    (see backend/capabilities.py) and sets its help text from the backend's
    own parameter declarations (description + guidance). Call at the end of
    initAlgorithm, after all addParameter calls.

    Silently no-ops when the backend is missing, old (no --help-json), or a
    parameter has no match — help presence never affects whether an
    algorithm runs. `overrides` maps QGIS parameter names to lsdtt keys for
    cases the naming conventions cannot resolve.
    """
    from ..backend.capabilities import probe_tool

    index = build_backend_help_index(probe_tool(tool_name))
    if not index:
        return
    overrides = overrides or {}
    for pdef in alg.parameterDefinitions():
        text = lookup_param_help(index, pdef.name(), overrides.get(pdef.name()))
        if text and hasattr(pdef, "setHelp"):
            pdef.setHelp(text)


def split_layer_source(source, kind="layer"):
    """Return (path, fname) for a QGIS layer source string.

    Vector sources may carry a "|layername=..." suffix; only the file part
    is passed to the lsdtt CLIs (they read the first layer).
    """
    parts = split_source(source)
    if parts is None:
        raise QgsProcessingException(
            "Could not resolve {} path/name from source: {}".format(kind, source)
        )
    return parts


def layer_file_path(source):
    """Return the plain file path of a layer source (pipe suffix stripped)."""
    return strip_layer_suffix(source)


def run_lsdtt_tool(tool_name, param_values, output_folder, output_prefix, feedback):
    """Write a parameter file, run the named lsdtt CLI, and log the result.

    Reproducibility contract: the parameter file is written next to the
    input DEM as <dem_stem>_qgis_analysis_NNN.driver and kept after the
    run, so any QGIS analysis can be re-run identically from the command
    line. A temporary directory is used only when there is no DEM input
    (or its directory is not writable), in which case the driver is
    discarded with the directory.

    Raises QgsProcessingException on a non-zero exit code.
    """
    log_stem = tool_name.replace("lsdtt-", "").replace("-", "_")

    driver_dir = param_values.get("read_path")
    dem_fname = param_values.get("read_fname")
    temp_ctx = None
    if (
        driver_dir
        and dem_fname
        and os.path.isdir(driver_dir)
        and os.access(driver_dir, os.W_OK)
    ):
        dem_stem = os.path.splitext(dem_fname)[0]
        params_name = next_qgis_driver_name(driver_dir, dem_stem)
        persistent = True
    else:
        temp_ctx = tempfile.TemporaryDirectory(prefix="lsdtt_qgis_")
        driver_dir = temp_ctx.name
        params_name = "{}.driver".format(log_stem)
        persistent = False

    try:
        params_file = os.path.join(driver_dir, params_name)
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        header = [
            "This parameter file was written by the lsdtt3 QGIS plugin",
            "for a {} analysis run on {}.".format(tool_name, timestamp),
            "Reproduce the same analysis from the command line with:",
            "  {} {} {}".format(tool_name, driver_dir, params_name),
        ]
        write_parameter_file(params_file, param_values, header_lines=header)
        if persistent:
            feedback.pushInfo(
                "Saved reproducible parameter file: {}".format(params_file)
            )

        backend = BackendManager()
        executable = backend.executable_path(tool_name)
        command = [executable, driver_dir, params_name]
        log_path = os.path.join(
            output_folder, "{}_{}.log".format(output_prefix, log_stem)
        )
        feedback.pushInfo("Running {}...".format(tool_name))
        result = run_cli(command, feedback=feedback, env=backend.runtime_environment())
        write_log(log_path, command, result)
        if result.returncode != 0:
            raise QgsProcessingException(
                "{} failed with exit code {}. See log: {}".format(
                    tool_name, result.returncode, log_path
                )
            )
    finally:
        if temp_ctx is not None:
            temp_ctx.cleanup()
    return log_path


def require_output(output_folder, output_prefix, suffix, description):
    """Return the expected output path, raising if the file was not created."""
    path = expected_output_path(output_folder, output_prefix, suffix)
    if not os.path.exists(path):
        raise QgsProcessingException(
            "Expected {} output not found: {}. "
            "Check the run log in the output folder for warnings.".format(
                description, path
            )
        )
    return path
