import os

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputVectorLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterDefinition,
    QgsProcessingParameterEnum,
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


class SideTributaryBasinsAlgorithm(OutputOverwriteMixin, QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    INPUT_POINTS = "INPUT_POINTS"
    POINTS_LAYER_NAME = "POINTS_LAYER_NAME"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    MODE = "MODE"
    RASTER_IS_FILLED = "RASTER_IS_FILLED"
    FILL_MIN_SLOPE = "FILL_MIN_SLOPE"
    THRESHOLD_CONTRIBUTING_PIXELS = "THRESHOLD_CONTRIBUTING_PIXELS"
    MIN_CONTRIBUTING_PIXELS = "MIN_CONTRIBUTING_PIXELS"
    REJECT_INCOMPLETE = "REJECT_INCOMPLETE"
    WRITE_MAIN_CHANNEL = "WRITE_MAIN_CHANNEL"
    WRITE_OUTLETS = "WRITE_OUTLETS"
    WRITE_OUTLINES = "WRITE_OUTLINES"
    VECTOR_FORMAT = "VECTOR_FORMAT"

    # Advanced
    USE_STOPPING_POINT = "USE_STOPPING_POINT"
    STOPPING_WINDOW = "STOPPING_WINDOW"
    CALCULATE_SPATIAL_INDEX = "CALCULATE_SPATIAL_INDEX"

    OUTPUT_MAIN_CHANNEL = "OUTPUT_MAIN_CHANNEL"
    OUTPUT_OUTLETS = "OUTPUT_OUTLETS"
    OUTPUT_OUTLINES = "OUTPUT_OUTLINES"

    MODE_OPTIONS = [
        "Side tributaries from seed points",
        "Basins along a single channel (1-2 points)",
    ]
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
        self.addParameter(
            QgsProcessingParameterEnum(
                self.MODE,
                "Extraction mode",
                options=self.MODE_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT_POINTS,
                "Seed points (any CRS; single-channel mode uses 1-2 points)",
                [QgsProcessing.TypeVectorPoint],
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.POINTS_LAYER_NAME,
                "Point layer name (empty = first layer; seed-points mode only)",
                "",
                optional=True,
            )
        )
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, "Output folder"))
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_side_trib"))
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
            QgsProcessingParameterNumber(
                self.MIN_CONTRIBUTING_PIXELS,
                "Minimum contributing pixels for a side tributary",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1000,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.REJECT_INCOMPLETE,
                "Reject side basins touching nodata or the raster edge",
                True,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_MAIN_CHANNEL, "Write main channel", True))
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_OUTLETS, "Write tributary outlet points", True))
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_OUTLINES, "Write tributary basin outlines", True))
        self.addParameter(
            QgsProcessingParameterEnum(
                self.VECTOR_FORMAT,
                "Vector output format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )

        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.USE_STOPPING_POINT,
                "Treat the last input point as a stopping point (seed-points mode)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.STOPPING_WINDOW,
                "Stopping-point window half-width (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=5,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(self.CALCULATE_SPATIAL_INDEX, "Calculate spatial index", True)
        )

        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_MAIN_CHANNEL, "Main channel"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_OUTLETS, "Tributary outlet points"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_OUTLINES, "Tributary basin outlines"))

        attach_backend_help(
            self,
            "lsdtt-basin-extraction",
            overrides={
                self.THRESHOLD_CONTRIBUTING_PIXELS: "threshold_contributing_pixels",
                self.WRITE_MAIN_CHANNEL: "write_side_tributary_main_channel",
                self.WRITE_OUTLETS: "write_side_tributary_outlets",
                self.WRITE_OUTLINES: "write_side_tributary_basin_outlines",
                self.INPUT_POINTS: "side_tributary_points_fname",
                self.MIN_CONTRIBUTING_PIXELS: "side_tributary_minimum_contributing_pixels",
                self.POINTS_LAYER_NAME: "side_tributary_points_layer",
                self.REJECT_INCOMPLETE: "reject_incomplete_side_tributary_basins",
                self.STOPPING_WINDOW: "side_tributary_stopping_point_window_pixels",
                self.USE_STOPPING_POINT: "include_side_tributary_stopping_point",
                self.VECTOR_FORMAT: "side_tributary_write_vector_format",
            },
        )

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        points_layer = self.parameterAsVectorLayer(parameters, self.INPUT_POINTS, context)
        if dem_layer is None:
            raise QgsProcessingException("Input DEM is required.")
        if points_layer is None:
            raise QgsProcessingException("Seed point layer is required.")
        read_path, read_fname = split_layer_source(dem_layer.source(), "DEM")

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        if not output_folder or not output_prefix:
            raise QgsProcessingException("Output folder and prefix are required.")
        os.makedirs(output_folder, exist_ok=True)

        mode = self.parameterAsEnum(parameters, self.MODE, context)
        seed_points_mode = mode == 0
        points_path = layer_file_path(points_layer.source())

        write_main_channel = self.parameterAsBool(parameters, self.WRITE_MAIN_CHANNEL, context)
        write_outlets = self.parameterAsBool(parameters, self.WRITE_OUTLETS, context)
        write_outlines = self.parameterAsBool(parameters, self.WRITE_OUTLINES, context)
        if not any([write_main_channel, write_outlets, write_outlines]):
            raise QgsProcessingException(
                "Enable at least one output (main channel, outlets, or outlines)."
            )

        vector_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.VECTOR_FORMAT, context)
        ]

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "raster_is_filled": self.parameterAsBool(parameters, self.RASTER_IS_FILLED, context),
            "fill_min_slope": self.parameterAsDouble(parameters, self.FILL_MIN_SLOPE, context),
            "threshold_contributing_pixels": self.parameterAsInt(
                parameters, self.THRESHOLD_CONTRIBUTING_PIXELS, context
            ),
            "side_tributary_minimum_contributing_pixels": self.parameterAsInt(
                parameters, self.MIN_CONTRIBUTING_PIXELS, context
            ),
            "reject_incomplete_side_tributary_basins": self.parameterAsBool(
                parameters, self.REJECT_INCOMPLETE, context
            ),
            "write_side_tributary_main_channel": write_main_channel,
            "write_side_tributary_outlets": write_outlets,
            "write_side_tributary_basin_outlines": write_outlines,
            "side_tributary_write_vector_format": vector_format,
            "calculate_spatial_index": self.parameterAsBool(
                parameters, self.CALCULATE_SPATIAL_INDEX, context
            ),
        }

        if seed_points_mode:
            param_values["extract_side_tributary_basins_from_points"] = True
            param_values["side_tributary_points_fname"] = points_path
            layer_name = self.parameterAsString(parameters, self.POINTS_LAYER_NAME, context)
            if layer_name:
                param_values["side_tributary_points_layer"] = layer_name
            param_values["include_side_tributary_stopping_point"] = self.parameterAsBool(
                parameters, self.USE_STOPPING_POINT, context
            )
            param_values["side_tributary_stopping_point_window_pixels"] = self.parameterAsInt(
                parameters, self.STOPPING_WINDOW, context
            )
        else:
            param_values["extract_basins_from_single_channel"] = True
            param_values["single_channel_points_fname"] = points_path

        run_lsdtt_tool(
            "lsdtt-basin-extraction", param_values, output_folder, output_prefix, feedback
        )

        ext = extension_for_format(vector_format)
        outputs = {}
        if write_main_channel:
            outputs[self.OUTPUT_MAIN_CHANNEL] = require_output(
                output_folder, output_prefix,
                "_side_tributary_main_channel.{}".format(ext), "main channel",
            )
        if write_outlets:
            outputs[self.OUTPUT_OUTLETS] = require_output(
                output_folder, output_prefix,
                "_side_tributary_outlets.{}".format(ext), "tributary outlets",
            )
        if write_outlines:
            outputs[self.OUTPUT_OUTLINES] = require_output(
                output_folder, output_prefix,
                "_side_tributary_basin_outlines.{}".format(ext), "tributary basin outlines",
            )
        return outputs

    def name(self):
        return "side_tributary_basins"

    def displayName(self):
        return "Side Tributary Basins"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Traces D8 channels downstream from seed points and extracts the "
            "side-tributary basins draining into them "
            "(lsdtt-basin-extraction side-tributary workflow). Seed-points "
            "mode traces each point to the DEM edge (optionally stopping at "
            "the last point); single-channel mode builds one main channel "
            "from one or two points and extracts the basins draining into it."
        )

    def createInstance(self):
        return SideTributaryBasinsAlgorithm()
