"""Shared helpers for lsdtt Processing algorithm wrappers.

The lsdtt CLIs report per-step failures as warnings and still exit 0, so
output-existence checks in the wrappers are the primary error detection.

Pure path/format logic lives in lsdtt_qgis.naming (qgis-free, unit-testable);
this module adds the qgis-dependent wrappers around it.
"""

import os
import tempfile

from qgis.core import QgsProcessingException

from ..backend.manager import BackendManager
from ..backend.runner import run_cli, write_log
from ..naming import expected_output_path, extension_for_format, split_source, strip_layer_suffix
from ..templates.lsdtt_params import write_parameter_file

__all__ = [
    "extension_for_format",
    "layer_file_path",
    "require_output",
    "run_lsdtt_tool",
    "split_layer_source",
]


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

    Raises QgsProcessingException on a non-zero exit code.
    """
    log_stem = tool_name.replace("lsdtt-", "").replace("-", "_")
    with tempfile.TemporaryDirectory(prefix="lsdtt_qgis_") as temp_dir:
        params_name = "{}.driver".format(log_stem)
        params_file = os.path.join(temp_dir, params_name)
        write_parameter_file(params_file, param_values)

        backend = BackendManager()
        executable = backend.executable_path(tool_name)
        command = [executable, temp_dir, params_name]
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
