"""Pure path/format helpers shared by the algorithm wrappers.

This module must stay free of qgis imports so the unit tests can run
without a QGIS installation.
"""

import os


def strip_layer_suffix(source):
    """Return the file path of a QGIS layer source, dropping any
    "|layername=..." style suffix."""
    return source.split("|")[0]


def split_source(source):
    """Return (directory, filename) for a layer source string, or None if
    either part cannot be resolved."""
    path = strip_layer_suffix(source)
    directory, fname = os.path.split(path)
    if not directory or not fname:
        return None
    return directory, fname


def extension_for_format(format_name):
    """Map an lsdtt vector format name to its file extension."""
    if format_name == "flatgeobuf":
        return "fgb"
    return format_name


def expected_output_path(output_folder, output_prefix, suffix):
    """Build the path of an lsdtt output file: <folder>/<prefix><suffix>."""
    return os.path.join(output_folder, "{}{}".format(output_prefix, suffix))
