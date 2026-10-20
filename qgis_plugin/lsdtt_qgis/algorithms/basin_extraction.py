import glob
import os
import re

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputFile,
    QgsProcessingOutputMultipleLayers,
    QgsProcessingOutputRasterLayer,
    QgsProcessingOutputVectorLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterDefinition,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFile,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
    QgsProcessingParameterVectorLayer,
)

from .support import (
    OutputOverwriteMixin,
    attach_backend_help,
    extension_for_format,
    layer_file_path,
    require_output,
    run_lsdtt_tool,
    split_layer_source,
)


def _grouped_outline_outputs(output_folder, output_prefix, suffix, output_format):
    extension = extension_for_format(output_format)
    if output_format == "gpkg":
        return [require_output(
            output_folder, output_prefix, "{}.{}".format(suffix, extension),
            "grouped basin outlines",
        )]

    pattern = os.path.join(
        output_folder, "{}{}_order_*.{}".format(output_prefix, suffix, extension)
    )
    matches = glob.glob(pattern)
    if not matches:
        raise QgsProcessingException(
            "Expected grouped basin outlines output not found: {}. "
            "Check the run log in the output folder for warnings.".format(pattern)
        )

    def order_number(path):
        match = re.search(r"_order_(-?\d+)\.{}$".format(re.escape(extension)), path)
        return int(match.group(1)) if match else 0

    return sorted(matches, key=lambda path: (order_number(path), path))


class BasinExtractionAlgorithm(OutputOverwriteMixin, QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    RASTER_IS_FILLED = "RASTER_IS_FILLED"
    FILL_MIN_SLOPE = "FILL_MIN_SLOPE"
    THRESHOLD_CONTRIBUTING_PIXELS = "THRESHOLD_CONTRIBUTING_PIXELS"
    SELECTION_MODE = "SELECTION_MODE"
    BASIN_OUTLET_FILE = "BASIN_OUTLET_FILE"
    CROSSING_LINE = "CROSSING_LINE"
    SELECTED_STREAM_ORDER = "SELECTED_STREAM_ORDER"
    WRITE_BASIN_OUTLINES = "WRITE_BASIN_OUTLINES"
    BASIN_OUTLINES_FORMAT = "BASIN_OUTLINES_FORMAT"
    WRITE_BASIN_OUTLETS = "WRITE_BASIN_OUTLETS"
    WRITE_BASINS_RASTER = "WRITE_BASINS_RASTER"
    CALCULATE_SPATIAL_INDEX = "CALCULATE_SPATIAL_INDEX"

    # Advanced: routing
    THREADS = "THREADS"
    BOUNDARY_CONDITIONS = "BOUNDARY_CONDITIONS"

    # Advanced: outlet snapping
    SNAPPING_METHOD = "SNAPPING_METHOD"
    SNAPPING_THRESHOLD = "SNAPPING_THRESHOLD"
    SNAPPING_MIN_RADIUS = "SNAPPING_MIN_RADIUS"
    SNAPPING_MAX_RADIUS = "SNAPPING_MAX_RADIUS"

    # Advanced: line crossings
    CROSSING_MIN_PIXELS = "CROSSING_MIN_PIXELS"

    # Advanced: basin filters
    FILTER_BY_STREAM_ORDER = "FILTER_BY_STREAM_ORDER"
    MIN_STREAM_ORDER = "MIN_STREAM_ORDER"
    MAX_STREAM_ORDER = "MAX_STREAM_ORDER"
    FILTER_BY_AREA = "FILTER_BY_AREA"
    MIN_AREA_PIXELS = "MIN_AREA_PIXELS"
    MAX_AREA_PIXELS = "MAX_AREA_PIXELS"
    FILTER_INCOMPLETE = "FILTER_INCOMPLETE"
    FILTER_ALL_BUT_LARGEST = "FILTER_ALL_BUT_LARGEST"
    FILTER_NESTED = "FILTER_NESTED"
    FILTER_NESTED_BY_ORDER = "FILTER_NESTED_BY_ORDER"

    # Advanced: extra outputs
    OUTLINES_AS_POLYGON = "OUTLINES_AS_POLYGON"
    OUTLINES_PER_ORDER = "OUTLINES_PER_ORDER"
    WRITE_BASINS_DATAMAP = "WRITE_BASINS_DATAMAP"
    WRITE_OUTLINES_CSV = "WRITE_OUTLINES_CSV"
    WRITE_SELECTED_OUTLETS_CSV = "WRITE_SELECTED_OUTLETS_CSV"
    WRITE_CHANNEL_LINES = "WRITE_CHANNEL_LINES"
    WRITE_CHANNEL_POINTS = "WRITE_CHANNEL_POINTS"
    WRITE_JUNCTIONS = "WRITE_JUNCTIONS"
    QA_VECTOR_FORMAT = "QA_VECTOR_FORMAT"

    OUTPUT_BASIN_OUTLINES = "OUTPUT_BASIN_OUTLINES"
    OUTPUT_BASIN_OUTLETS = "OUTPUT_BASIN_OUTLETS"
    OUTPUT_BASINS_RASTER = "OUTPUT_BASINS_RASTER"
    OUTPUT_BASINS_DATAMAP = "OUTPUT_BASINS_DATAMAP"
    OUTPUT_OUTLINES_CSV = "OUTPUT_OUTLINES_CSV"
    OUTPUT_SELECTED_OUTLETS_CSV = "OUTPUT_SELECTED_OUTLETS_CSV"
    OUTPUT_LINE_CROSSING_CSV = "OUTPUT_LINE_CROSSING_CSV"
    OUTPUT_CHANNEL_LINES = "OUTPUT_CHANNEL_LINES"
    OUTPUT_CHANNEL_POINTS = "OUTPUT_CHANNEL_POINTS"
    OUTPUT_JUNCTIONS = "OUTPUT_JUNCTIONS"

    SELECTION_OPTIONS = [
        "All basins",
        "From outlets",
        "From outlets (nested)",
        "From line crossings",
        "All basins of stream order",
    ]
    SNAPPING_OPTIONS = ["downstream", "expanding_square", "none"]
    VECTOR_FORMAT_OPTIONS = ["flatgeobuf", "gpkg"]

    def _add_advanced(self, param):
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_DEM,
                "Input DEM",
                [QgsProcessing.TypeRaster],
            )
        )
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, "Output folder"))
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_basin"))
        self.add_overwrite_parameter()
        self.addParameter(QgsProcessingParameterBoolean(self.RASTER_IS_FILLED, "Raster already filled", False))
        self.addParameter(
            QgsProcessingParameterNumber(
                self.FILL_MIN_SLOPE,
                "Fill min slope",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0001,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.THRESHOLD_CONTRIBUTING_PIXELS,
                "Threshold contributing pixels",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1000,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.SELECTION_MODE,
                "Selection mode",
                options=self.SELECTION_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterFile(
                self.BASIN_OUTLET_FILE,
                "Basin outlet file (required for outlet modes; CSV needs latitude/longitude columns)",
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.CROSSING_LINE,
                "Crossing line layer (required for line-crossing mode; same CRS as DEM)",
                [QgsProcessing.TypeVectorLine],
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.SELECTED_STREAM_ORDER,
                "Selection mode: stream order (used by 'All basins of stream order')",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=3,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_BASIN_OUTLINES, "Write basin outlines", True))
        self.addParameter(
            QgsProcessingParameterEnum(
                self.BASIN_OUTLINES_FORMAT,
                "Basin outlines format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_BASIN_OUTLETS,
                "Write requested outlet points (outlet modes only)",
                True,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_BASINS_RASTER, "Write basins raster", False))
        self.addParameter(QgsProcessingParameterBoolean(self.CALCULATE_SPATIAL_INDEX, "Calculate spatial index", True))

        # --- Advanced: routing ---
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.THREADS,
                "Threads (0 = auto)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterString(
                self.BOUNDARY_CONDITIONS,
                "Boundary conditions (four characters NESW; o=open, n=no-flux, p=periodic)",
                "oooo",
            )
        )

        # --- Advanced: outlet snapping ---
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.SNAPPING_METHOD,
                "Outlet snapping method",
                options=self.SNAPPING_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SNAPPING_THRESHOLD,
                "Outlet snapping contributing pixels threshold",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1000,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SNAPPING_MIN_RADIUS,
                "Outlet snapping minimum search radius (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SNAPPING_MAX_RADIUS,
                "Outlet snapping maximum search radius (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=25,
            )
        )

        # --- Advanced: line crossings ---
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CROSSING_MIN_PIXELS,
                "Line crossings: minimum contributing pixels at a crossing",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=500,
            )
        )

        # --- Advanced: basin filters ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FILTER_BY_STREAM_ORDER,
                "Filter: keep only basins within a stream-order range",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MIN_STREAM_ORDER,
                "Filter: minimum stream order (-1 = ignore)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=-1,
                defaultValue=-1,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MAX_STREAM_ORDER,
                "Filter: maximum stream order (-1 = ignore)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=-1,
                defaultValue=-1,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FILTER_BY_AREA,
                "Filter: keep only basins within an area range (pixels)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MIN_AREA_PIXELS,
                "Filter: minimum basin area in pixels (-1 = ignore)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=-1,
                defaultValue=-1,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MAX_AREA_PIXELS,
                "Filter: maximum basin area in pixels (-1 = ignore)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=-1,
                defaultValue=-1,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FILTER_INCOMPLETE,
                "Filter: remove incomplete basins (touching nodata or edge)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FILTER_ALL_BUT_LARGEST,
                "Filter: keep only the largest basin",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FILTER_NESTED,
                "Filter: remove nested basins",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FILTER_NESTED_BY_ORDER,
                "Filter: remove nested basins of the same stream order",
                False,
            )
        )

        # --- Advanced: extra outputs ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.OUTLINES_AS_POLYGON,
                "Write basin outlines as polygons (recommended)",
                True,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.OUTLINES_PER_ORDER,
                "Write basin outlines grouped per stream order",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_BASINS_DATAMAP,
                "Write basins datamap CSV (junction number, outlet index)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_OUTLINES_CSV,
                "Write basin outlines CSV",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_SELECTED_OUTLETS_CSV,
                "Write selected basin outlets CSV (lat/lon)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_CHANNEL_LINES,
                "QA: write channel network lines",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_CHANNEL_POINTS,
                "QA: write channel network points",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_JUNCTIONS,
                "QA: write channel network junctions",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.QA_VECTOR_FORMAT,
                "QA vector format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )

        self.addOutput(
            QgsProcessingOutputMultipleLayers(self.OUTPUT_BASIN_OUTLINES, "Basin outlines")
        )
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_BASIN_OUTLETS, "Requested outlet points"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_BASINS_RASTER, "Basins raster"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_BASINS_DATAMAP, "Basins datamap CSV"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_OUTLINES_CSV, "Basin outlines CSV"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_SELECTED_OUTLETS_CSV, "Selected basin outlets CSV"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_LINE_CROSSING_CSV, "Line-crossing basins CSV"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_CHANNEL_LINES, "Channel network lines"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_CHANNEL_POINTS, "Channel network points"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_JUNCTIONS, "Channel network junctions"))

        attach_backend_help(
            self,
            "lsdtt-basin-extraction",
            overrides={
                self.FILTER_BY_AREA: "filter_basins_by_area",
                self.FILTER_INCOMPLETE: "filter_incomplete_basins",
                self.FILTER_NESTED: "filter_nested_basins",
                self.MAX_AREA_PIXELS: "basins_max_area_pixels",
                self.MAX_STREAM_ORDER: "basins_max_stream_order",
                self.MIN_AREA_PIXELS: "basins_min_area_pixels",
                self.MIN_STREAM_ORDER: "basins_min_stream_order",
                self.OUTLINES_AS_POLYGON: "write_basin_outlines_as_polygon",
                self.THRESHOLD_CONTRIBUTING_PIXELS: "threshold_contributing_pixels",
                self.WRITE_BASIN_OUTLETS: "write_basin_outlet_points",
                self.WRITE_CHANNEL_LINES: "write_channel_network_lines",
                self.WRITE_CHANNEL_POINTS: "write_channel_network_points",
                self.WRITE_JUNCTIONS: "write_channel_network_junctions",
                self.WRITE_OUTLINES_CSV: "write_basin_outlines_to_csv",
                self.WRITE_SELECTED_OUTLETS_CSV: "write_selected_basin_outlets_csv",
                self.CROSSING_LINE: "crossing_line_fname",
                self.CROSSING_MIN_PIXELS: "crossing_minimum_contributing_pixels",
                self.SELECTED_STREAM_ORDER: "selected_stream_order",
                self.FILTER_ALL_BUT_LARGEST: "filter_all_but_largest_basin",
                self.FILTER_BY_STREAM_ORDER: "filter_basins_by_stream_order",
                self.FILTER_NESTED_BY_ORDER: "filter_nested_by_stream_order",
                self.OUTLINES_PER_ORDER: "write_basin_outlines_per_order",
                self.SNAPPING_MAX_RADIUS: "basin_outlet_snap_max_search_radius_pixels",
                self.SNAPPING_METHOD: "basin_outlet_snapping_method",
                self.SNAPPING_MIN_RADIUS: "basin_outlet_snap_min_search_radius_pixels",
                self.SNAPPING_THRESHOLD: "basin_outlet_contributing_pixels_threshold",
            },
        )

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        if dem_layer is None:
            raise QgsProcessingException("Input DEM is required.")
        read_path, read_fname = split_layer_source(dem_layer.source(), "DEM")

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        if not output_folder or not output_prefix:
            raise QgsProcessingException("Output folder and prefix are required.")
        os.makedirs(output_folder, exist_ok=True)

        selection = self.parameterAsEnum(parameters, self.SELECTION_MODE, context)
        outlet_file = self.parameterAsFile(parameters, self.BASIN_OUTLET_FILE, context)
        crossing_layer = self.parameterAsVectorLayer(parameters, self.CROSSING_LINE, context)
        selected_stream_order = self.parameterAsInt(parameters, self.SELECTED_STREAM_ORDER, context)

        outlet_mode = selection in (1, 2)
        crossing_mode = selection == 3
        order_mode = selection == 4
        if outlet_mode and not outlet_file:
            raise QgsProcessingException("Basin outlet file is required for outlet-based selection modes.")
        if crossing_mode and crossing_layer is None:
            raise QgsProcessingException("A crossing line layer is required for line-crossing mode.")
        if order_mode and selected_stream_order <= 0:
            raise QgsProcessingException("Selected stream order must be > 0 for this selection mode.")

        boundary_conditions = self.parameterAsString(parameters, self.BOUNDARY_CONDITIONS, context)
        if len(boundary_conditions) != 4 or any(c not in "onpb" for c in boundary_conditions.lower()):
            raise QgsProcessingException(
                "Boundary conditions must be four characters from o, n, p, b (e.g. oooo)."
            )

        write_outlines = self.parameterAsBool(parameters, self.WRITE_BASIN_OUTLINES, context)
        write_outlets = self.parameterAsBool(parameters, self.WRITE_BASIN_OUTLETS, context)
        write_basins_raster = self.parameterAsBool(parameters, self.WRITE_BASINS_RASTER, context)
        write_datamap = self.parameterAsBool(parameters, self.WRITE_BASINS_DATAMAP, context)
        write_outlines_csv = self.parameterAsBool(parameters, self.WRITE_OUTLINES_CSV, context)
        write_selected_csv = self.parameterAsBool(parameters, self.WRITE_SELECTED_OUTLETS_CSV, context)
        write_channel_lines = self.parameterAsBool(parameters, self.WRITE_CHANNEL_LINES, context)
        write_channel_points = self.parameterAsBool(parameters, self.WRITE_CHANNEL_POINTS, context)
        write_junctions = self.parameterAsBool(parameters, self.WRITE_JUNCTIONS, context)
        outlines_as_polygon = write_outlines and self.parameterAsBool(
            parameters, self.OUTLINES_AS_POLYGON, context
        )
        outlines_per_order = write_outlines and self.parameterAsBool(
            parameters, self.OUTLINES_PER_ORDER, context
        )

        if not any([
            write_outlines, write_basins_raster, write_datamap, write_outlines_csv,
            write_selected_csv, write_channel_lines, write_channel_points, write_junctions,
            outlet_mode and write_outlets, crossing_mode,
        ]):
            raise QgsProcessingException("Enable at least one basin output.")

        outlines_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.BASIN_OUTLINES_FORMAT, context)
        ]
        qa_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.QA_VECTOR_FORMAT, context)
        ]

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "n_threads": self.parameterAsInt(parameters, self.THREADS, context),
            "boundary_conditions": boundary_conditions,
            "raster_is_filled": self.parameterAsBool(parameters, self.RASTER_IS_FILLED, context),
            "fill_min_slope": self.parameterAsDouble(parameters, self.FILL_MIN_SLOPE, context),
            "threshold_contributing_pixels": self.parameterAsInt(
                parameters, self.THRESHOLD_CONTRIBUTING_PIXELS, context
            ),
            # selection
            "select_all_basins": selection == 0,
            "select_basins_from_outlets": selection == 1,
            "select_basins_from_outlets_nested": selection == 2,
            "select_basins_from_line_crossings": crossing_mode,
            "select_all_basins_of_stream_order": order_mode,
            "selected_stream_order": selected_stream_order,
            # filters
            "filter_basins_by_stream_order": self.parameterAsBool(
                parameters, self.FILTER_BY_STREAM_ORDER, context
            ),
            "basins_min_stream_order": self.parameterAsInt(parameters, self.MIN_STREAM_ORDER, context),
            "basins_max_stream_order": self.parameterAsInt(parameters, self.MAX_STREAM_ORDER, context),
            "filter_basins_by_area": self.parameterAsBool(parameters, self.FILTER_BY_AREA, context),
            "basins_min_area_pixels": self.parameterAsInt(parameters, self.MIN_AREA_PIXELS, context),
            "basins_max_area_pixels": self.parameterAsInt(parameters, self.MAX_AREA_PIXELS, context),
            "filter_incomplete_basins": self.parameterAsBool(parameters, self.FILTER_INCOMPLETE, context),
            "filter_all_but_largest_basin": self.parameterAsBool(
                parameters, self.FILTER_ALL_BUT_LARGEST, context
            ),
            "filter_nested_basins": self.parameterAsBool(parameters, self.FILTER_NESTED, context),
            "filter_nested_by_stream_order": self.parameterAsBool(
                parameters, self.FILTER_NESTED_BY_ORDER, context
            ),
            # outputs
            "write_basin_outlines": write_outlines,
            "write_basin_outlines_as_polygon": outlines_as_polygon,
            "write_basin_outlines_per_order": outlines_per_order,
            "basin_outlines_write_vector_format": outlines_format,
            "write_basins_raster": write_basins_raster,
            "write_basins_datamap": write_datamap,
            "write_basin_outlines_to_csv": write_outlines_csv,
            "write_selected_basin_outlets_csv": write_selected_csv,
            "write_channel_network_lines": write_channel_lines,
            "write_vector_format": qa_format,
            "write_channel_network_points": write_channel_points,
            "write_channel_network_junctions": write_junctions,
            "calculate_spatial_index": self.parameterAsBool(parameters, self.CALCULATE_SPATIAL_INDEX, context),
        }

        if outlet_mode:
            param_values.update({
                "basin_outlet_fname": outlet_file,
                "basin_outlet_snapping_method": self.SNAPPING_OPTIONS[
                    self.parameterAsEnum(parameters, self.SNAPPING_METHOD, context)
                ],
                "basin_outlet_contributing_pixels_threshold": self.parameterAsInt(
                    parameters, self.SNAPPING_THRESHOLD, context
                ),
                "basin_outlet_snap_min_search_radius_pixels": self.parameterAsInt(
                    parameters, self.SNAPPING_MIN_RADIUS, context
                ),
                "basin_outlet_snap_max_search_radius_pixels": self.parameterAsInt(
                    parameters, self.SNAPPING_MAX_RADIUS, context
                ),
                "write_basin_outlet_points": write_outlets,
            })

        if crossing_mode:
            param_values.update({
                "crossing_line_fname": layer_file_path(crossing_layer.source()),
                "crossing_minimum_contributing_pixels": self.parameterAsInt(
                    parameters, self.CROSSING_MIN_PIXELS, context
                ),
            })

        run_lsdtt_tool(
            "lsdtt-basin-extraction", param_values, output_folder, output_prefix, feedback
        )

        # Basin outlines are written as _junction_basin_polygons when polygons
        # are requested and _junction_basin_outlines when written as lines.
        outlines_suffix = (
            "_junction_basin_polygons" if outlines_as_polygon else "_junction_basin_outlines"
        )

        outputs = {}
        if write_outlines:
            if outlines_per_order:
                outline_outputs = _grouped_outline_outputs(
                    output_folder, output_prefix, outlines_suffix, outlines_format
                )
                feedback.pushInfo(
                    "Grouped basin outlines: {} output file(s): {}".format(
                        len(outline_outputs), ", ".join(outline_outputs)
                    )
                )
            else:
                outline_outputs = [require_output(
                    output_folder, output_prefix,
                    "{}.{}".format(outlines_suffix, extension_for_format(outlines_format)),
                    "basin outlines",
                )]
            outputs[self.OUTPUT_BASIN_OUTLINES] = outline_outputs
        if outlet_mode and write_outlets:
            outputs[self.OUTPUT_BASIN_OUTLETS] = require_output(
                output_folder, output_prefix, "_basin_outlet_points.fgb",
                "requested outlet points",
            )
        if write_basins_raster:
            outputs[self.OUTPUT_BASINS_RASTER] = require_output(
                output_folder, output_prefix, "_junction_basins.tif", "basins raster"
            )
        if write_datamap:
            outputs[self.OUTPUT_BASINS_DATAMAP] = require_output(
                output_folder, output_prefix, "_junction_basins_datamap.csv",
                "basins datamap",
            )
        if write_outlines_csv:
            outputs[self.OUTPUT_OUTLINES_CSV] = require_output(
                output_folder, output_prefix, "_basin_outlines.csv", "basin outlines CSV"
            )
        if write_selected_csv:
            outputs[self.OUTPUT_SELECTED_OUTLETS_CSV] = require_output(
                output_folder, output_prefix, "_selected_basin_outlets.csv",
                "selected basin outlets CSV",
            )
        if crossing_mode:
            outputs[self.OUTPUT_LINE_CROSSING_CSV] = require_output(
                output_folder, output_prefix, "_line_crossing_basins.csv",
                "line-crossing basins CSV",
            )
        if write_channel_lines:
            outputs[self.OUTPUT_CHANNEL_LINES] = require_output(
                output_folder, output_prefix,
                "_channel_network_lines.{}".format(extension_for_format(qa_format)),
                "channel network lines",
            )
        if write_channel_points:
            outputs[self.OUTPUT_CHANNEL_POINTS] = require_output(
                output_folder, output_prefix,
                "_channel_network_points.{}".format(extension_for_format(qa_format)),
                "channel network points",
            )
        if write_junctions:
            outputs[self.OUTPUT_JUNCTIONS] = require_output(
                output_folder, output_prefix,
                "_channel_network_junctions.{}".format(extension_for_format(qa_format)),
                "channel network junctions",
            )
        return outputs

    def name(self):
        return "basin_extraction"

    def displayName(self):
        return "Basin Extraction"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Runs lsdtt-basin-extraction: builds D8 flow routing and a channel "
            "network, selects basins (all junctions, supplied outlets, nested "
            "outlets, or channels crossing a line), optionally filters them "
            "(stream order, area, completeness, nesting), and writes basin "
            "outlines, rasters, and CSV summaries. Backend path is managed in "
            "plugin settings."
        )

    def createInstance(self):
        return BasinExtractionAlgorithm()
