import os

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputFile,
    QgsProcessingOutputRasterLayer,
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
    attach_backend_help,
    extension_for_format,
    layer_file_path,
    require_output,
    run_lsdtt_tool,
    split_layer_source,
)


class ChannelExtractionAlgorithm(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    RASTER_IS_FILLED = "RASTER_IS_FILLED"
    FILL_MIN_SLOPE = "FILL_MIN_SLOPE"
    SOURCE_THRESHOLD = "SOURCE_THRESHOLD"
    SOURCE_ALGORITHM = "SOURCE_ALGORITHM"
    WRITE_CHANNEL_LINES = "WRITE_CHANNEL_LINES"
    LINES_FORMAT = "LINES_FORMAT"
    WRITE_CHANNEL_POINTS = "WRITE_CHANNEL_POINTS"
    POINTS_FORMAT = "POINTS_FORMAT"
    WRITE_JUNCTIONS = "WRITE_JUNCTIONS"
    JUNCTIONS_FORMAT = "JUNCTIONS_FORMAT"
    CALCULATE_SPATIAL_INDEX = "CALCULATE_SPATIAL_INDEX"

    # Advanced: routing
    THREADS = "THREADS"
    BOUNDARY_CONDITIONS = "BOUNDARY_CONDITIONS"

    # Advanced: wiener source extraction
    WIENER_PRUNING_AREA = "WIENER_PRUNING_AREA"
    WIENER_FITTING_RADIUS = "WIENER_FITTING_RADIUS"
    WIENER_CC_THRESHOLD = "WIENER_CC_THRESHOLD"
    WIENER_WRITE_FILTERED_DEM = "WIENER_WRITE_FILTERED_DEM"
    WIENER_WRITE_CURVATURE = "WIENER_WRITE_CURVATURE"
    WIENER_WRITE_QQ = "WIENER_WRITE_QQ"

    # Advanced: flow-routing outputs
    WRITE_FLOW_LINES = "WRITE_FLOW_LINES"
    FLOW_LINES_FORMAT = "FLOW_LINES_FORMAT"
    WRITE_D8_FLOWDIR = "WRITE_D8_FLOWDIR"

    # Advanced: burned stream profiles
    STREAM_VECTOR = "STREAM_VECTOR"
    WRITE_BURNED_PROFILES = "WRITE_BURNED_PROFILES"
    BURNED_PROFILES_FORMAT = "BURNED_PROFILES_FORMAT"

    OUTPUT_LINES = "OUTPUT_LINES"
    OUTPUT_POINTS = "OUTPUT_POINTS"
    OUTPUT_JUNCTIONS = "OUTPUT_JUNCTIONS"
    OUTPUT_FLOW_LINES = "OUTPUT_FLOW_LINES"
    OUTPUT_D8_FLOWDIR = "OUTPUT_D8_FLOWDIR"
    OUTPUT_BURNED_PROFILES = "OUTPUT_BURNED_PROFILES"
    OUTPUT_WIENER_FILTERED_DEM = "OUTPUT_WIENER_FILTERED_DEM"
    OUTPUT_WIENER_CURVATURE = "OUTPUT_WIENER_CURVATURE"
    OUTPUT_WIENER_QQ = "OUTPUT_WIENER_QQ"

    SOURCE_ALGORITHM_OPTIONS = ["threshold", "wiener"]
    VECTOR_FORMAT_OPTIONS = ["flatgeobuf", "gpkg"]
    POINT_FORMAT_OPTIONS = ["flatgeobuf", "gpkg", "csv"]
    BURNED_PROFILES_FORMAT_OPTIONS = ["csv", "flatgeobuf", "gpkg"]

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
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_channel"))
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
                self.SOURCE_THRESHOLD,
                "Source threshold (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1000,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.SOURCE_ALGORITHM,
                "Source extraction algorithm (wiener requires an FFTW-enabled backend)",
                options=self.SOURCE_ALGORITHM_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_CHANNEL_LINES, "Write channel lines", True))
        self.addParameter(
            QgsProcessingParameterEnum(
                self.LINES_FORMAT,
                "Lines format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_CHANNEL_POINTS, "Write channel points", False))
        self.addParameter(
            QgsProcessingParameterEnum(
                self.POINTS_FORMAT,
                "Points format",
                options=self.POINT_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_JUNCTIONS, "Write junction points", False))
        self.addParameter(
            QgsProcessingParameterEnum(
                self.JUNCTIONS_FORMAT,
                "Junction format",
                options=self.POINT_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
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

        # --- Advanced: wiener source extraction ---
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.WIENER_PRUNING_AREA,
                "Wiener: pruning drainage area (m2)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1000.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.WIENER_FITTING_RADIUS,
                "Wiener: surface fitting radius (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=6.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.WIENER_CC_THRESHOLD,
                "Wiener: connected components threshold (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=100,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WIENER_WRITE_FILTERED_DEM,
                "Wiener: write filtered DEM diagnostic",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WIENER_WRITE_CURVATURE,
                "Wiener: write tangential curvature diagnostic",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WIENER_WRITE_QQ,
                "Wiener: write QQ table CSV",
                False,
            )
        )

        # --- Advanced: flow-routing outputs ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_FLOW_LINES,
                "Write flow lines (large output)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.FLOW_LINES_FORMAT,
                "Flow lines format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_D8_FLOWDIR,
                "Write ArcMap-style D8 flow-direction raster",
                False,
            )
        )

        # --- Advanced: burned stream profiles ---
        self._add_advanced(
            QgsProcessingParameterVectorLayer(
                self.STREAM_VECTOR,
                "Stream polyline layer for burned profiles (same CRS as DEM)",
                [QgsProcessing.TypeVectorLine],
                optional=True,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_BURNED_PROFILES,
                "Write burned stream profiles",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.BURNED_PROFILES_FORMAT,
                "Burned profiles format",
                options=self.BURNED_PROFILES_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )

        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_LINES, "Channel lines"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_POINTS, "Channel points output"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_JUNCTIONS, "Junctions output"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_FLOW_LINES, "Flow lines"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_D8_FLOWDIR, "D8 flow-direction raster"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_BURNED_PROFILES, "Burned stream profiles"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_WIENER_FILTERED_DEM, "Wiener filtered DEM"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_WIENER_CURVATURE, "Wiener tangential curvature"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_WIENER_QQ, "Wiener QQ table"))

        attach_backend_help(
            self,
            "lsdtt-channel-extraction",
            overrides={
                self.JUNCTIONS_FORMAT: "channel_network_junctions_output_format",
                self.LINES_FORMAT: "channel_network_lines_output_format",
                self.POINTS_FORMAT: "channel_network_points_output_format",
                self.SOURCE_THRESHOLD: "source_pixels_threshold",
                self.WIENER_WRITE_CURVATURE: "wiener_write_tangential_curvature",
                self.WIENER_WRITE_QQ: "wiener_write_qq_table",
                self.WRITE_CHANNEL_LINES: "write_channel_network_lines",
                self.WRITE_CHANNEL_POINTS: "write_channel_network_points",
                self.WRITE_D8_FLOWDIR: "write_arcmap_d8_flowdirection_raster",
                self.WRITE_JUNCTIONS: "write_channel_network_junctions",
                self.BURNED_PROFILES_FORMAT: "burned_stream_profiles_output_format",
                self.SOURCE_ALGORITHM: "source_extraction_algorithm_choice",
                self.STREAM_VECTOR: "stream_vector_path",
                self.WIENER_CC_THRESHOLD: "wiener_connected_components_threshold_pixels",
                self.WIENER_FITTING_RADIUS: "wiener_surface_fitting_radius_m",
                self.WIENER_PRUNING_AREA: "wiener_pruning_drainage_area_m2",
                self.WRITE_BURNED_PROFILES: "write_burned_stream_profiles",
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

        source_algorithm = self.SOURCE_ALGORITHM_OPTIONS[
            self.parameterAsEnum(parameters, self.SOURCE_ALGORITHM, context)
        ]
        if source_algorithm == "wiener":
            from ..backend.capabilities import require_feature
            error = require_feature(
                "lsdtt-channel-extraction", "fftw",
                "Wiener source extraction", feedback,
            )
            if error:
                raise QgsProcessingException(error)

        boundary_conditions = self.parameterAsString(parameters, self.BOUNDARY_CONDITIONS, context)
        if len(boundary_conditions) != 4 or any(c not in "onpb" for c in boundary_conditions.lower()):
            raise QgsProcessingException(
                "Boundary conditions must be four characters from o, n, p, b (e.g. oooo)."
            )

        lines_format = self.VECTOR_FORMAT_OPTIONS[self.parameterAsEnum(parameters, self.LINES_FORMAT, context)]
        points_format = self.POINT_FORMAT_OPTIONS[self.parameterAsEnum(parameters, self.POINTS_FORMAT, context)]
        junctions_format = self.POINT_FORMAT_OPTIONS[self.parameterAsEnum(parameters, self.JUNCTIONS_FORMAT, context)]
        flow_lines_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.FLOW_LINES_FORMAT, context)
        ]
        burned_profiles_format = self.BURNED_PROFILES_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.BURNED_PROFILES_FORMAT, context)
        ]

        write_lines = self.parameterAsBool(parameters, self.WRITE_CHANNEL_LINES, context)
        write_points = self.parameterAsBool(parameters, self.WRITE_CHANNEL_POINTS, context)
        write_junctions = self.parameterAsBool(parameters, self.WRITE_JUNCTIONS, context)
        write_flow_lines = self.parameterAsBool(parameters, self.WRITE_FLOW_LINES, context)
        write_d8_flowdir = self.parameterAsBool(parameters, self.WRITE_D8_FLOWDIR, context)
        write_burned_profiles = self.parameterAsBool(parameters, self.WRITE_BURNED_PROFILES, context)
        wiener_write_filtered = self.parameterAsBool(parameters, self.WIENER_WRITE_FILTERED_DEM, context)
        wiener_write_curvature = self.parameterAsBool(parameters, self.WIENER_WRITE_CURVATURE, context)
        wiener_write_qq = self.parameterAsBool(parameters, self.WIENER_WRITE_QQ, context)

        if not any([write_lines, write_points, write_junctions, write_flow_lines, write_d8_flowdir]):
            raise QgsProcessingException(
                "Enable at least one channel output (lines, points, junctions, "
                "flow lines, or flow-direction raster)."
            )

        stream_layer = self.parameterAsVectorLayer(parameters, self.STREAM_VECTOR, context)
        if write_burned_profiles and stream_layer is None:
            raise QgsProcessingException(
                "A stream polyline layer is required to write burned stream profiles."
            )

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "n_threads": self.parameterAsInt(parameters, self.THREADS, context),
            "boundary_conditions": boundary_conditions,
            "raster_is_filled": self.parameterAsBool(parameters, self.RASTER_IS_FILLED, context),
            "fill_min_slope": self.parameterAsDouble(parameters, self.FILL_MIN_SLOPE, context),
            "source_pixels_threshold": self.parameterAsInt(parameters, self.SOURCE_THRESHOLD, context),
            "source_extraction_algorithm_choice": source_algorithm,
            "write_channel_network_lines": write_lines,
            "channel_network_lines_output_format": lines_format,
            "write_channel_network_points": write_points,
            "channel_network_points_output_format": points_format,
            "write_channel_network_junctions": write_junctions,
            "channel_network_junctions_output_format": junctions_format,
            "write_flow_lines": write_flow_lines,
            "flow_lines_output_format": flow_lines_format,
            "write_arcmap_d8_flowdirection_raster": write_d8_flowdir,
            "calculate_spatial_index": self.parameterAsBool(parameters, self.CALCULATE_SPATIAL_INDEX, context),
        }

        if source_algorithm == "wiener":
            param_values.update({
                "wiener_pruning_drainage_area_m2": self.parameterAsDouble(
                    parameters, self.WIENER_PRUNING_AREA, context
                ),
                "wiener_surface_fitting_radius_m": self.parameterAsDouble(
                    parameters, self.WIENER_FITTING_RADIUS, context
                ),
                "wiener_connected_components_threshold_pixels": self.parameterAsInt(
                    parameters, self.WIENER_CC_THRESHOLD, context
                ),
                "wiener_write_filtered_dem": wiener_write_filtered,
                "wiener_write_tangential_curvature": wiener_write_curvature,
                "wiener_write_qq_table": wiener_write_qq,
                # cog guarantees .tif diagnostics regardless of the input format
                "wiener_diagnostics_raster_format": "cog",
            })

        if write_burned_profiles:
            param_values.update({
                "write_burned_stream_profiles": True,
                "stream_vector_path": layer_file_path(stream_layer.source()),
                "stream_burn_mode": "burn_in_place",
                "burned_stream_profiles_output_format": burned_profiles_format,
            })

        run_lsdtt_tool(
            "lsdtt-channel-extraction", param_values, output_folder, output_prefix, feedback
        )

        outputs = {}
        if write_lines:
            outputs[self.OUTPUT_LINES] = require_output(
                output_folder, output_prefix,
                "_channel_network_lines.{}".format(extension_for_format(lines_format)),
                "channel lines",
            )
        if write_points:
            outputs[self.OUTPUT_POINTS] = require_output(
                output_folder, output_prefix,
                "_channel_network_points.{}".format(extension_for_format(points_format)),
                "channel points",
            )
        if write_junctions:
            outputs[self.OUTPUT_JUNCTIONS] = require_output(
                output_folder, output_prefix,
                "_channel_network_junctions.{}".format(extension_for_format(junctions_format)),
                "channel junctions",
            )
        if write_flow_lines:
            outputs[self.OUTPUT_FLOW_LINES] = require_output(
                output_folder, output_prefix,
                "_flow_lines.{}".format(extension_for_format(flow_lines_format)),
                "flow lines",
            )
        if write_d8_flowdir:
            outputs[self.OUTPUT_D8_FLOWDIR] = require_output(
                output_folder, output_prefix, "_arcmap_d8_flowdirection.tif",
                "D8 flow direction",
            )
        if write_burned_profiles:
            outputs[self.OUTPUT_BURNED_PROFILES] = require_output(
                output_folder, output_prefix,
                "_burned_stream_profiles.{}".format(extension_for_format(burned_profiles_format)),
                "burned stream profiles",
            )
        if source_algorithm == "wiener":
            if wiener_write_filtered:
                outputs[self.OUTPUT_WIENER_FILTERED_DEM] = require_output(
                    output_folder, output_prefix, "_wiener_filtered_dem.tif",
                    "Wiener filtered DEM",
                )
            if wiener_write_curvature:
                outputs[self.OUTPUT_WIENER_CURVATURE] = require_output(
                    output_folder, output_prefix, "_wiener_tangential_curvature.tif",
                    "Wiener tangential curvature",
                )
            if wiener_write_qq:
                outputs[self.OUTPUT_WIENER_QQ] = require_output(
                    output_folder, output_prefix, "_wiener_qq.csv", "Wiener QQ table"
                )
        return outputs

    def name(self):
        return "channel_extraction"

    def displayName(self):
        return "Channel Extraction"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Runs lsdtt-channel-extraction: builds D8 flow routing, extracts "
            "channel sources (contributing-pixel threshold or Wiener-filtered "
            "curvature), and writes the channel network as lines, points, and "
            "junctions. Advanced options add flow lines, an ArcMap-style D8 "
            "flow-direction raster, and burned stream profiles from a supplied "
            "polyline network. Backend path is managed in plugin settings."
        )

    def createInstance(self):
        return ChannelExtractionAlgorithm()
