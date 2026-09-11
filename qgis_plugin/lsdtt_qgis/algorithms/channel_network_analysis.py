import os

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputFile,
    QgsProcessingOutputVectorLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterDefinition,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFile,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
)

from .support import (
    attach_backend_help,
    extension_for_format,
    require_output,
    run_lsdtt_tool,
    split_layer_source,
)


class ChannelNetworkAnalysisAlgorithm(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    RASTER_IS_FILLED = "RASTER_IS_FILLED"
    FILL_MIN_SLOPE = "FILL_MIN_SLOPE"
    THRESHOLD_CONTRIBUTING_PIXELS = "THRESHOLD_CONTRIBUTING_PIXELS"
    SELECTION_MODE = "SELECTION_MODE"
    BASIN_OUTLET_FILE = "BASIN_OUTLET_FILE"
    VECTOR_FORMAT = "VECTOR_FORMAT"
    WRITE_LINK_BEARINGS = "WRITE_LINK_BEARINGS"
    WRITE_JUNCTION_ANGLES = "WRITE_JUNCTION_ANGLES"
    WRITE_STRAHLER_CHANNELS = "WRITE_STRAHLER_CHANNELS"
    WRITE_TOKUNAGA_CHANNELS = "WRITE_TOKUNAGA_CHANNELS"
    WRITE_TOKUNAGA_CHANNEL_DATA = "WRITE_TOKUNAGA_CHANNEL_DATA"
    WRITE_TOKUNAGA_SUMMARY = "WRITE_TOKUNAGA_SUMMARY"
    WRITE_AREA_RATIO_DATA = "WRITE_AREA_RATIO_DATA"
    WRITE_AREA_RATIO_SUMMARY = "WRITE_AREA_RATIO_SUMMARY"
    WRITE_COM_BEARINGS = "WRITE_COM_BEARINGS"

    # Advanced
    THREADS = "THREADS"
    BOUNDARY_CONDITIONS = "BOUNDARY_CONDITIONS"
    CALCULATE_SPATIAL_INDEX = "CALCULATE_SPATIAL_INDEX"
    SNAPPING_METHOD = "SNAPPING_METHOD"
    SNAPPING_THRESHOLD = "SNAPPING_THRESHOLD"
    SNAPPING_MIN_RADIUS = "SNAPPING_MIN_RADIUS"
    SNAPPING_MAX_RADIUS = "SNAPPING_MAX_RADIUS"
    LINK_SAMPLING_MODE = "LINK_SAMPLING_MODE"
    LINK_SAMPLING_FIXED_LENGTH = "LINK_SAMPLING_FIXED_LENGTH"
    LINK_SAMPLING_FIXED_N_PIXELS = "LINK_SAMPLING_FIXED_N_PIXELS"
    LINK_SAMPLING_MIN_NODES = "LINK_SAMPLING_MIN_NODES"
    LINK_BEARING_ANCHOR = "LINK_BEARING_ANCHOR"
    WRITE_LINK_BEARING_VECTORS = "WRITE_LINK_BEARING_VECTORS"
    WRITE_JUNCTION_ANGLE_VECTORS = "WRITE_JUNCTION_ANGLE_VECTORS"
    AREA_RATIO_BOOTSTRAP = "AREA_RATIO_BOOTSTRAP"
    AREA_RATIO_SEED = "AREA_RATIO_SEED"
    ONLY_COMPLETE_BASINS = "ONLY_COMPLETE_BASINS"
    MINIMUM_BASIN_PIXELS = "MINIMUM_BASIN_PIXELS"

    OUTPUT_LINK_BEARINGS = "OUTPUT_LINK_BEARINGS"
    OUTPUT_JUNCTION_ANGLES = "OUTPUT_JUNCTION_ANGLES"
    OUTPUT_STRAHLER_CHANNELS = "OUTPUT_STRAHLER_CHANNELS"
    OUTPUT_TOKUNAGA_CHANNELS = "OUTPUT_TOKUNAGA_CHANNELS"
    OUTPUT_TOKUNAGA_CHANNEL_DATA = "OUTPUT_TOKUNAGA_CHANNEL_DATA"
    OUTPUT_TOKUNAGA_SUMMARY = "OUTPUT_TOKUNAGA_SUMMARY"
    OUTPUT_TOKUNAGA_PAIR_COUNTS = "OUTPUT_TOKUNAGA_PAIR_COUNTS"
    OUTPUT_TOKUNAGA_OFFSET_STATS = "OUTPUT_TOKUNAGA_OFFSET_STATS"
    OUTPUT_AREA_RATIOS = "OUTPUT_AREA_RATIOS"
    OUTPUT_AREA_RATIO_SUMMARY = "OUTPUT_AREA_RATIO_SUMMARY"
    OUTPUT_AREA_RATIO_DECADES = "OUTPUT_AREA_RATIO_DECADES"
    OUTPUT_COM_BEARINGS = "OUTPUT_COM_BEARINGS"
    OUTPUT_LINK_BEARING_VECTORS = "OUTPUT_LINK_BEARING_VECTORS"
    OUTPUT_JUNCTION_ANGLE_VECTORS = "OUTPUT_JUNCTION_ANGLE_VECTORS"

    VECTOR_FORMAT_OPTIONS = ["flatgeobuf", "gpkg", "csv"]
    SELECTION_MODE_OPTIONS = [
        "All junctions",
        "Upstream of all complete-basin outlets",
        "Upstream of outlet points from a file",
    ]
    SELECTION_MODE_KEYS = [
        "analyse_all_junctions",
        "analyse_upstream_of_all_outlets",
        "analyse_upstream_of_outlet_csv",
    ]
    SNAPPING_OPTIONS = ["downstream", "expanding_square", "none"]
    LINK_SAMPLING_OPTIONS = ["full_link", "fixed_length", "fixed_n_pixels"]
    ANCHOR_OPTIONS = ["upstream_end", "downstream_end"]

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
        self.addParameter(
            QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_network")
        )
        self.addParameter(
            QgsProcessingParameterBoolean(self.RASTER_IS_FILLED, "Raster already filled", False)
        )
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
                "Junction selection mode",
                options=self.SELECTION_MODE_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterFile(
                self.BASIN_OUTLET_FILE,
                "Basin outlet file (CSV/FlatGeobuf/GeoPackage; outlet-file mode only)",
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.VECTOR_FORMAT,
                "Output format for the bearing / angle / channel layers",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_LINK_BEARINGS, "Write channel-link bearings", True
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(self.WRITE_JUNCTION_ANGLES, "Write junction angles", False)
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_STRAHLER_CHANNELS, "Write Strahler channels", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_TOKUNAGA_CHANNELS, "Write Tokunaga channels", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_TOKUNAGA_CHANNEL_DATA, "Write Tokunaga channel attributes (CSV)", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_TOKUNAGA_SUMMARY, "Write Tokunaga summary (CSV)", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_AREA_RATIO_DATA, "Write junction area ratios (CSV)", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_AREA_RATIO_SUMMARY, "Write junction area-ratio summary (CSV)", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_COM_BEARINGS, "Write junction centre-of-mass bearings (CSV)", False
            )
        )

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
            QgsProcessingParameterString(self.BOUNDARY_CONDITIONS, "Boundary conditions", "oooo")
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(self.CALCULATE_SPATIAL_INDEX, "Calculate spatial index", True)
        )
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
                "Outlet snapping contributing-pixel threshold",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1000,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SNAPPING_MIN_RADIUS,
                "Minimum snapping search radius (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SNAPPING_MAX_RADIUS,
                "Maximum snapping search radius (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=25,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.LINK_SAMPLING_MODE,
                "Link sampling mode",
                options=self.LINK_SAMPLING_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.LINK_SAMPLING_FIXED_LENGTH,
                "Fixed sampling length (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1000.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.LINK_SAMPLING_FIXED_N_PIXELS,
                "Fixed sampling length (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=10,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.LINK_SAMPLING_MIN_NODES,
                "Minimum sampled nodes per link",
                type=QgsProcessingParameterNumber.Integer,
                minValue=2,
                defaultValue=2,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.LINK_BEARING_ANCHOR,
                "Link-bearing sampling anchor",
                options=self.ANCHOR_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_LINK_BEARING_VECTORS, "Write channel-link bearing QA vectors", False
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_JUNCTION_ANGLE_VECTORS, "Write junction-angle bearing QA vectors", False
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.AREA_RATIO_BOOTSTRAP,
                "Junction area-ratio bootstrap replicates",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=1000,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.AREA_RATIO_SEED,
                "Junction area-ratio bootstrap seed",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=1,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.ONLY_COMPLETE_BASINS, "Restrict analyses to complete basins", True
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MINIMUM_BASIN_PIXELS,
                "Minimum upslope basin pixels per junction",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=0,
            )
        )

        self.addOutput(
            QgsProcessingOutputVectorLayer(self.OUTPUT_LINK_BEARINGS, "Channel link bearings")
        )
        self.addOutput(
            QgsProcessingOutputVectorLayer(self.OUTPUT_JUNCTION_ANGLES, "Junction angles")
        )
        self.addOutput(
            QgsProcessingOutputVectorLayer(self.OUTPUT_STRAHLER_CHANNELS, "Strahler channels")
        )
        self.addOutput(
            QgsProcessingOutputVectorLayer(self.OUTPUT_TOKUNAGA_CHANNELS, "Tokunaga channels")
        )
        self.addOutput(
            QgsProcessingOutputFile(self.OUTPUT_TOKUNAGA_CHANNEL_DATA, "Tokunaga channel data")
        )
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_TOKUNAGA_SUMMARY, "Tokunaga fit summary"))
        self.addOutput(
            QgsProcessingOutputFile(self.OUTPUT_TOKUNAGA_PAIR_COUNTS, "Tokunaga pair counts")
        )
        self.addOutput(
            QgsProcessingOutputFile(self.OUTPUT_TOKUNAGA_OFFSET_STATS, "Tokunaga offset statistics")
        )
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_AREA_RATIOS, "Junction area ratios"))
        self.addOutput(
            QgsProcessingOutputFile(self.OUTPUT_AREA_RATIO_SUMMARY, "Junction area-ratio summary")
        )
        self.addOutput(
            QgsProcessingOutputFile(self.OUTPUT_AREA_RATIO_DECADES, "Junction area-ratio decades")
        )
        self.addOutput(
            QgsProcessingOutputFile(self.OUTPUT_COM_BEARINGS, "Junction centre-of-mass bearings")
        )
        self.addOutput(
            QgsProcessingOutputVectorLayer(
                self.OUTPUT_LINK_BEARING_VECTORS, "Channel-link bearing QA vectors"
            )
        )
        self.addOutput(
            QgsProcessingOutputVectorLayer(
                self.OUTPUT_JUNCTION_ANGLE_VECTORS, "Junction-angle bearing QA vectors"
            )
        )

        attach_backend_help(
            self,
            "lsdtt-channel-network-analysis",
            overrides={
                self.WRITE_LINK_BEARINGS: "write_channel_link_bearings",
                self.WRITE_JUNCTION_ANGLES: "write_junction_angles",
                self.WRITE_STRAHLER_CHANNELS: "write_strahler_channels",
                self.WRITE_TOKUNAGA_CHANNELS: "write_tokunaga_channels",
                self.WRITE_TOKUNAGA_CHANNEL_DATA: "write_tokunaga_channel_data",
                self.WRITE_TOKUNAGA_SUMMARY: "write_tokunaga_summary_data",
                self.WRITE_AREA_RATIO_DATA: "write_junction_area_ratio_data",
                self.WRITE_AREA_RATIO_SUMMARY: "write_junction_area_ratio_summary",
                self.WRITE_COM_BEARINGS: "write_junction_centre_of_mass_bearings",
                self.WRITE_LINK_BEARING_VECTORS: "write_channel_link_bearing_vectors",
                self.WRITE_JUNCTION_ANGLE_VECTORS: "write_junction_angle_bearing_vectors",
                self.VECTOR_FORMAT: "channel_link_bearings_output_format",
                self.SNAPPING_METHOD: "basin_outlet_snapping_method",
                self.SNAPPING_THRESHOLD: "basin_outlet_contributing_pixels_threshold",
                self.SNAPPING_MIN_RADIUS: "basin_outlet_min_search_radius_pixels",
                self.SNAPPING_MAX_RADIUS: "basin_outlet_max_search_radius_pixels",
                self.LINK_SAMPLING_MODE: "link_sampling_mode",
                self.LINK_SAMPLING_FIXED_LENGTH: "link_sampling_fixed_length_metres",
                self.LINK_SAMPLING_FIXED_N_PIXELS: "link_sampling_fixed_n_pixels",
                self.LINK_SAMPLING_MIN_NODES: "link_sampling_min_nodes",
                self.LINK_BEARING_ANCHOR: "channel_link_bearing_sample_anchor",
                self.AREA_RATIO_BOOTSTRAP: "junction_area_ratio_bootstrap",
                self.AREA_RATIO_SEED: "junction_area_ratio_seed",
                self.THRESHOLD_CONTRIBUTING_PIXELS: "threshold_contributing_pixels",
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

        selection_mode = self.parameterAsEnum(parameters, self.SELECTION_MODE, context)
        outlet_file = self.parameterAsFile(parameters, self.BASIN_OUTLET_FILE, context)
        if selection_mode == 2 and not outlet_file:
            raise QgsProcessingException(
                "A basin outlet file is required for the outlet-file selection mode."
            )

        write_link_bearings = self.parameterAsBool(parameters, self.WRITE_LINK_BEARINGS, context)
        write_junction_angles = self.parameterAsBool(parameters, self.WRITE_JUNCTION_ANGLES, context)
        write_strahler = self.parameterAsBool(parameters, self.WRITE_STRAHLER_CHANNELS, context)
        write_tokunaga = self.parameterAsBool(parameters, self.WRITE_TOKUNAGA_CHANNELS, context)
        write_tokunaga_data = self.parameterAsBool(
            parameters, self.WRITE_TOKUNAGA_CHANNEL_DATA, context
        )
        write_tokunaga_summary = self.parameterAsBool(
            parameters, self.WRITE_TOKUNAGA_SUMMARY, context
        )
        write_area_ratio_data = self.parameterAsBool(parameters, self.WRITE_AREA_RATIO_DATA, context)
        write_area_ratio_summary = self.parameterAsBool(
            parameters, self.WRITE_AREA_RATIO_SUMMARY, context
        )
        write_com_bearings = self.parameterAsBool(parameters, self.WRITE_COM_BEARINGS, context)
        write_link_bearing_vectors = self.parameterAsBool(
            parameters, self.WRITE_LINK_BEARING_VECTORS, context
        )
        write_junction_angle_vectors = self.parameterAsBool(
            parameters, self.WRITE_JUNCTION_ANGLE_VECTORS, context
        )
        if not any([
            write_link_bearings, write_junction_angles, write_strahler, write_tokunaga,
            write_tokunaga_data, write_tokunaga_summary, write_area_ratio_data,
            write_area_ratio_summary, write_com_bearings, write_link_bearing_vectors,
            write_junction_angle_vectors,
        ]):
            raise QgsProcessingException("Enable at least one analysis output.")

        vector_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.VECTOR_FORMAT, context)
        ]
        # The QA bearing-vector layers are geometry-only: CSV is not a valid
        # format for them, so they fall back to FlatGeobuf.
        qa_format = "flatgeobuf" if vector_format == "csv" else vector_format

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "n_threads": self.parameterAsInt(parameters, self.THREADS, context),
            "boundary_conditions": self.parameterAsString(
                parameters, self.BOUNDARY_CONDITIONS, context
            ),
            "calculate_spatial_index": self.parameterAsBool(
                parameters, self.CALCULATE_SPATIAL_INDEX, context
            ),
            "raster_is_filled": self.parameterAsBool(parameters, self.RASTER_IS_FILLED, context),
            "fill_min_slope": self.parameterAsDouble(parameters, self.FILL_MIN_SLOPE, context),
            "threshold_contributing_pixels": self.parameterAsInt(
                parameters, self.THRESHOLD_CONTRIBUTING_PIXELS, context
            ),
            self.SELECTION_MODE_KEYS[selection_mode]: True,
            "only_complete_basins": self.parameterAsBool(
                parameters, self.ONLY_COMPLETE_BASINS, context
            ),
            "minimum_basin_pixels": self.parameterAsInt(
                parameters, self.MINIMUM_BASIN_PIXELS, context
            ),
            # link sampling for bearings / junction angles
            "link_sampling_mode": self.LINK_SAMPLING_OPTIONS[
                self.parameterAsEnum(parameters, self.LINK_SAMPLING_MODE, context)
            ],
            "link_sampling_fixed_length_metres": self.parameterAsDouble(
                parameters, self.LINK_SAMPLING_FIXED_LENGTH, context
            ),
            "link_sampling_fixed_n_pixels": self.parameterAsInt(
                parameters, self.LINK_SAMPLING_FIXED_N_PIXELS, context
            ),
            "link_sampling_min_nodes": self.parameterAsInt(
                parameters, self.LINK_SAMPLING_MIN_NODES, context
            ),
            "channel_link_bearing_sample_anchor": self.ANCHOR_OPTIONS[
                self.parameterAsEnum(parameters, self.LINK_BEARING_ANCHOR, context)
            ],
            # analysis outputs
            "write_channel_link_bearings": write_link_bearings,
            "channel_link_bearings_output_format": vector_format,
            "write_junction_angles": write_junction_angles,
            "junction_angles_output_format": vector_format,
            "write_strahler_channels": write_strahler,
            "strahler_channels_output_format": vector_format,
            "write_tokunaga_channels": write_tokunaga,
            "tokunaga_channels_output_format": vector_format,
            "write_tokunaga_channel_data": write_tokunaga_data,
            "write_tokunaga_summary_data": write_tokunaga_summary,
            "write_junction_area_ratio_data": write_area_ratio_data,
            "write_junction_area_ratio_summary": write_area_ratio_summary,
            "junction_area_ratio_bootstrap": self.parameterAsInt(
                parameters, self.AREA_RATIO_BOOTSTRAP, context
            ),
            "junction_area_ratio_seed": self.parameterAsInt(
                parameters, self.AREA_RATIO_SEED, context
            ),
            "write_junction_centre_of_mass_bearings": write_com_bearings,
            # QA bearing vectors
            "write_channel_link_bearing_vectors": write_link_bearing_vectors,
            "channel_link_bearing_vectors_output_format": qa_format,
            "write_junction_angle_bearing_vectors": write_junction_angle_vectors,
            "junction_angle_bearing_vectors_output_format": qa_format,
        }

        if outlet_file:
            param_values.update({
                "basin_outlet_file": outlet_file,
                "basin_outlet_snapping_method": self.SNAPPING_OPTIONS[
                    self.parameterAsEnum(parameters, self.SNAPPING_METHOD, context)
                ],
                "basin_outlet_contributing_pixels_threshold": self.parameterAsInt(
                    parameters, self.SNAPPING_THRESHOLD, context
                ),
                "basin_outlet_min_search_radius_pixels": self.parameterAsInt(
                    parameters, self.SNAPPING_MIN_RADIUS, context
                ),
                "basin_outlet_max_search_radius_pixels": self.parameterAsInt(
                    parameters, self.SNAPPING_MAX_RADIUS, context
                ),
            })

        run_lsdtt_tool(
            "lsdtt-channel-network-analysis", param_values, output_folder, output_prefix, feedback
        )

        ext = extension_for_format(vector_format)
        outputs = {}
        if write_link_bearings:
            outputs[self.OUTPUT_LINK_BEARINGS] = require_output(
                output_folder, output_prefix,
                "_channel_link_bearings.{}".format(ext), "channel link bearings",
            )
        if write_junction_angles:
            outputs[self.OUTPUT_JUNCTION_ANGLES] = require_output(
                output_folder, output_prefix,
                "_junction_angles.{}".format(ext), "junction angles",
            )
        if write_strahler:
            outputs[self.OUTPUT_STRAHLER_CHANNELS] = require_output(
                output_folder, output_prefix,
                "_strahler_channels.{}".format(ext), "Strahler channels",
            )
        if write_tokunaga:
            outputs[self.OUTPUT_TOKUNAGA_CHANNELS] = require_output(
                output_folder, output_prefix,
                "_tokunaga_channels.{}".format(ext), "Tokunaga channels",
            )
        if write_tokunaga_data:
            outputs[self.OUTPUT_TOKUNAGA_CHANNEL_DATA] = require_output(
                output_folder, output_prefix,
                "_tokunaga_channel_data.csv", "Tokunaga channel data",
            )
        if write_tokunaga_summary:
            outputs[self.OUTPUT_TOKUNAGA_PAIR_COUNTS] = require_output(
                output_folder, output_prefix,
                "_tokunaga_summary_pair_counts.csv", "Tokunaga pair counts",
            )
            outputs[self.OUTPUT_TOKUNAGA_OFFSET_STATS] = require_output(
                output_folder, output_prefix,
                "_tokunaga_summary_offset_stats.csv", "Tokunaga offset statistics",
            )
            outputs[self.OUTPUT_TOKUNAGA_SUMMARY] = require_output(
                output_folder, output_prefix,
                "_tokunaga_summary_fit.csv", "Tokunaga fit summary",
            )
        if write_area_ratio_data:
            outputs[self.OUTPUT_AREA_RATIOS] = require_output(
                output_folder, output_prefix,
                "_junction_area_ratios.csv", "junction area ratios",
            )
        if write_area_ratio_summary:
            outputs[self.OUTPUT_AREA_RATIO_SUMMARY] = require_output(
                output_folder, output_prefix,
                "_junction_area_ratio_summary.csv", "junction area-ratio summary",
            )
            outputs[self.OUTPUT_AREA_RATIO_DECADES] = require_output(
                output_folder, output_prefix,
                "_junction_area_ratio_decades.csv", "junction area-ratio decades",
            )
        if write_com_bearings:
            outputs[self.OUTPUT_COM_BEARINGS] = require_output(
                output_folder, output_prefix,
                "_junction_com_bearings.csv", "junction centre-of-mass bearings",
            )
        if write_link_bearing_vectors:
            outputs[self.OUTPUT_LINK_BEARING_VECTORS] = require_output(
                output_folder, output_prefix,
                "_channel_link_bearing_vectors.{}".format(extension_for_format(qa_format)),
                "channel-link bearing QA vectors",
            )
        if write_junction_angle_vectors:
            outputs[self.OUTPUT_JUNCTION_ANGLE_VECTORS] = require_output(
                output_folder, output_prefix,
                "_junction_angle_bearing_vectors.{}".format(extension_for_format(qa_format)),
                "junction-angle bearing QA vectors",
            )
        return outputs

    def name(self):
        return "channel_network_analysis"

    def displayName(self):
        return "Channel Network Analysis"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Runs lsdtt-channel-network-analysis: builds a threshold-defined "
            "channel network and analyses its topology and geometry. Outputs "
            "channel-link bearings, junction angles, Strahler and Tokunaga "
            "channel layers with their summary tables, junction "
            "drainage-area ratios, and junction centre-of-mass bearings. "
            "Junctions can be selected network-wide, upstream of all "
            "complete-basin outlets, or upstream of supplied outlet points. "
            "Backend path is managed in plugin settings."
        )

    def createInstance(self):
        return ChannelNetworkAnalysisAlgorithm()
