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


def next_qgis_driver_name(directory, dem_stem):
    """Return the next free parameter-file name for a QGIS-driven run.

    Files follow <dem_stem>_qgis_analysis_NNN.driver with NNN starting at
    001 and incrementing per existing file in `directory`, so every QGIS
    run leaves a reproducible driver next to the input DEM.
    """
    prefix = "{}_qgis_analysis_".format(dem_stem)
    highest = 0
    try:
        entries = os.listdir(directory)
    except OSError:
        entries = []
    for entry in entries:
        if not entry.startswith(prefix) or not entry.endswith(".driver"):
            continue
        counter = entry[len(prefix):-len(".driver")]
        if counter.isdigit():
            highest = max(highest, int(counter))
    return "{}{:03d}.driver".format(prefix, highest + 1)


# Aliases from common QGIS wrapper parameter names to lsdtt parameter keys,
# for cases where the lowercased QGIS name is not itself an lsdtt key.
_HELP_KEY_ALIASES = {
    "input_dem": "read_fname",
    "output_folder": "write_path",
    "output_prefix": "write_fname",
    "output_format": "write_output_format",
    "threads": "n_threads",
}


def build_backend_help_index(payload):
    """Return {lsdtt_key: parameter_record} from a --help-json payload.

    Returns {} for None or malformed payloads so callers can no-op.
    """
    if not isinstance(payload, dict):
        return {}
    params = payload.get("parameters")
    if not isinstance(params, list):
        return {}
    index = {}
    for record in params:
        if isinstance(record, dict) and record.get("name"):
            # Lowercased keys: lsdtt names are matched case-insensitively
            # (e.g. QGIS "CHI_A0" vs backend "chi_A0").
            index[record["name"].lower()] = record
    return index


def lookup_param_help(index, qgis_name, override=None):
    """Return backend help text for a QGIS parameter name, or None.

    Matching order: explicit override key, lowercased QGIS name, alias
    table, then the lsdtt naming-convention suffixes `_raster` (bool
    output switches), `_choice` (algorithm enums), and the
    `*_format` -> `*_output_format` vector-format convention. All
    comparisons are case-insensitive.
    """
    key = override.lower() if override else None
    if key is None:
        base = qgis_name.lower()
        candidates = [
            base,
            _HELP_KEY_ALIASES.get(base),
            base + "_raster",
            base + "_choice",
        ]
        if base.endswith("_format"):
            candidates.append(base[: -len("_format")] + "_output_format")
        for candidate in candidates:
            if candidate and candidate in index:
                key = candidate
                break
    record = index.get(key) if key else None
    if not isinstance(record, dict):
        return None
    description = (record.get("description") or "").strip()
    guidance = (record.get("guidance") or "").strip()
    parts = [p for p in (description, guidance) if p]
    if not parts:
        return None
    return "\n\n".join(parts)


def extract_json_payload(text):
    """Return the JSON object embedded in CLI output, or None.

    lsdtt CLIs print a banner before the --help-json payload, so parse from
    the first opening brace (the banner contains none).
    """
    import json

    start = text.find("{")
    if start < 0:
        return None
    try:
        return json.loads(text[start:])
    except ValueError:
        return None


def write_line_geojson(path, points, crs_authid=None):
    """Write one ordered polyline as a single-feature GeoJSON LineString.

    `points` is a sequence of (x, y) pairs in the target CRS. The lsdtt CLIs
    do not reproject vectors, so the caller must already have transformed
    them; `crs_authid` (e.g. "EPSG:32633") is recorded only so the file
    displays in the right place when opened in QGIS. Returns `path`.
    """
    import json

    coordinates = [[float(x), float(y)] for x, y in points]
    if len(coordinates) < 2:
        raise ValueError("a line needs at least two points")
    collection = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"id": 1},
                "geometry": {"type": "LineString", "coordinates": coordinates},
            }
        ],
    }
    if crs_authid:
        collection["crs"] = {"type": "name", "properties": {"name": crs_authid}}
    with open(path, "w", encoding="utf-8") as dst:
        json.dump(collection, dst)
        dst.write("\n")
    return path


def source_layer_name(source):
    """Return the layer a QGIS vector source refers to.

    The "|layername=..." option when present (GeoPackage and other
    multi-layer files); otherwise the file stem, which is the layer name
    GDAL gives single-layer formats such as FlatGeobuf and Shapefile.
    """
    for option in source.split("|")[1:]:
        key, _, value = option.partition("=")
        if key.strip().lower() == "layername" and value.strip():
            return value.strip()
    return os.path.splitext(os.path.basename(strip_layer_suffix(source)))[0]
