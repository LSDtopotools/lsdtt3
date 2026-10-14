"""Shared helpers for lsdtt Processing algorithm wrappers.

The lsdtt CLIs report per-step failures as warnings and still exit 0, so
output-existence checks in the wrappers are the primary error detection.

Pure path/format logic lives in lsdtt_qgis.naming (qgis-free, unit-testable);
this module adds the qgis-dependent wrappers around it.
"""

import csv
import datetime
import os
import tempfile

from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsGeometry,
    QgsProcessingException,
)

from ..backend.manager import BackendManager
from ..backend.manifest import backend_compatibility_error, load_manifest
from ..backend.runner import run_cli, write_log
from ..naming import (
    build_backend_help_index,
    expected_output_path,
    extension_for_format,
    lookup_param_help,
    next_qgis_driver_name,
    split_source,
    strip_layer_suffix,
    write_line_geojson,
)
from ..templates.lsdtt_params import write_parameter_file

__all__ = [
    "attach_backend_help",
    "export_baseline_line",
    "export_points_as_latlon_csv",
    "extension_for_format",
    "layer_file_path",
    "require_output",
    "run_lsdtt_tool",
    "split_layer_source",
]

# Attribute columns passed through to an exported outlet CSV when the source
# layer carries them: lsdtt-basin-stats picks these up by name.
PASSTHROUGH_POINT_FIELDS = ("gauge_id", "gauge_name")


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


def export_points_as_latlon_csv(layer, csv_path, context, feedback=None):
    """Write a QGIS point layer to a latitude/longitude CSV in EPSG:4326.

    Some lsdtt CLIs (lsdtt-basin-stats) read outlet points only as a CSV
    with `latitude`/`longitude` columns, so a layer in any CRS is
    reprojected here. Any PASSTHROUGH_POINT_FIELDS present on the layer are
    copied through as extra columns. Returns the CSV path.
    """
    wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
    transform = QgsCoordinateTransform(layer.crs(), wgs84, context.transformContext())

    field_names = [f.name() for f in layer.fields()]
    passthrough = [name for name in PASSTHROUGH_POINT_FIELDS if name in field_names]

    written = 0
    with open(csv_path, "w", encoding="utf-8", newline="") as dst:
        writer = csv.writer(dst)
        writer.writerow(["latitude", "longitude"] + passthrough)
        for feature in layer.getFeatures():
            geometry = feature.geometry()
            if geometry is None or geometry.isEmpty():
                continue
            point = geometry.centroid().asPoint()
            try:
                point = transform.transform(point)
            except Exception as exc:
                raise QgsProcessingException(
                    "Could not reproject outlet points to EPSG:4326: {}".format(exc)
                )
            writer.writerow(
                ["{:.10f}".format(point.y()), "{:.10f}".format(point.x())]
                + [feature[name] for name in passthrough]
            )
            written += 1

    if written == 0:
        raise QgsProcessingException(
            "Outlet point layer contains no usable point geometries."
        )
    if feedback is not None:
        feedback.pushInfo(
            "Wrote {} outlet point(s) to {}".format(written, csv_path)
        )
    return csv_path


def export_baseline_line(layer, target_crs, path, context, feedback=None):
    """Write one line feature of a QGIS layer, reprojected, as GeoJSON.

    lsdtt-swath-profile uses the first LineString of its baseline file and
    reads its coordinates as being in the DEM's CRS, without reprojecting.
    So the baseline is taken from the first selected feature (or the first
    feature, when nothing is selected), transformed to `target_crs`, and
    written as a single LineString. A multipart feature with more than one
    part is refused rather than silently cut to its first part.
    Returns `path`.
    """
    transform = QgsCoordinateTransform(layer.crs(), target_crs, context.transformContext())
    use_selection = layer.selectedFeatureCount() > 0
    features = layer.selectedFeatures() if use_selection else layer.getFeatures()

    for feature in features:
        geometry = feature.geometry()
        if geometry is None or geometry.isEmpty():
            continue
        geometry = QgsGeometry(geometry)
        try:
            geometry.transform(transform)
        except Exception as exc:
            raise QgsProcessingException(
                "Could not reproject the baseline to the DEM CRS: {}".format(exc)
            )
        if geometry.isMultipart():
            parts = geometry.asMultiPolyline()
            if len(parts) != 1:
                raise QgsProcessingException(
                    "The baseline feature has {} parts; the swath needs one "
                    "continuous line. Merge or split it first (for example "
                    "with Multipart to Singleparts).".format(len(parts))
                )
            line = parts[0]
        else:
            line = geometry.asPolyline()
        if len(line) < 2:
            continue
        write_line_geojson(
            path, [(point.x(), point.y()) for point in line], target_crs.authid()
        )
        if feedback is not None:
            feedback.pushInfo(
                "Wrote baseline ({} vertices, {} feature id {}) to {}".format(
                    len(line),
                    "selected" if use_selection else "first",
                    feature.id(),
                    path,
                )
            )
        return path

    raise QgsProcessingException(
        "The baseline layer has no line feature with at least two vertices."
    )


def run_lsdtt_tool(tool_name, param_values, output_folder, output_prefix, feedback):
    """Write a parameter file, run the named lsdtt CLI, and log the result.

    Reproducibility contract: the parameter file is written next to the
    input DEM as <dem_stem>_qgis_analysis_NNN.driver and kept after the
    run, so any QGIS analysis can be re-run identically from the command
    line. A temporary directory is used only when there is no DEM input
    (or its directory is not writable), in which case the driver is
    discarded with the directory.

    Raises QgsProcessingException, before anything runs, if the backend is
    incompatible or lacks this algorithm's CLI (only this one: an older
    bundle still runs the algorithms it has); and on a non-zero exit code.
    """
    backend = BackendManager()
    compatibility_error = backend_compatibility_error(
        load_manifest(backend.backend_path())
    )
    if compatibility_error:
        raise QgsProcessingException(compatibility_error)
    try:
        executable = backend.executable_path(tool_name)
    except RuntimeError:
        raise QgsProcessingException(
            "The configured backend ({}) has no {}. Install the qgis backend "
            "bundle that matches this plugin.".format(backend.backend_path(), tool_name)
        )

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
