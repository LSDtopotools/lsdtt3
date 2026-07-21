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
    export_points_as_latlon_csv,
    extension_for_format,
    require_output,
    run_lsdtt_tool,
    split_layer_source,
)


class BasinStatsAlgorithm(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    INPUT_OUTLETS = "INPUT_OUTLETS"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    CALCULATE_CHANNEL_STATS = "CALCULATE_CHANNEL_STATS"
    THRESHOLD_CONTRIBUTING_PIXELS = "THRESHOLD_CONTRIBUTING_PIXELS"
    CALCULATE_SHAPE_METRICS = "CALCULATE_SHAPE_METRICS"
    CALCULATE_COM_BEARING = "CALCULATE_COM_BEARING"
    CALCULATE_JUNCTION_AREA_RATIOS = "CALCULATE_JUNCTION_AREA_RATIOS"
    WRITE_BASIN_POLYGONS = "WRITE_BASIN_POLYGONS"
    BASIN_POLYGON_FORMAT = "BASIN_POLYGON_FORMAT"
    REJECT_NODATA_EDGE_BASINS = "REJECT_NODATA_EDGE_BASINS"
    SNAPPING_METHOD = "SNAPPING_METHOD"
    SNAPPING_THRESHOLD = "SNAPPING_THRESHOLD"

    # Advanced
    THREADS = "THREADS"
    OUTPUT_FORMAT = "OUTPUT_FORMAT"
    SNAPPING_MIN_RADIUS = "SNAPPING_MIN_RADIUS"
    SNAPPING_MAX_RADIUS = "SNAPPING_MAX_RADIUS"
    MINIMUM_ELEVATION = "MINIMUM_ELEVATION"
    FILL_BEFORE_ROUTING = "FILL_BEFORE_ROUTING"
    FILL_MIN_SLOPE = "FILL_MIN_SLOPE"
    BASIN_LENGTH_DEFINITION = "BASIN_LENGTH_DEFINITION"
    SHAPE_PERIMETER_SIMPLIFY_CELLS = "SHAPE_PERIMETER_SIMPLIFY_CELLS"
    JUNCTION_AREA_RATIO_BOOTSTRAP = "JUNCTION_AREA_RATIO_BOOTSTRAP"
    WRITE_BASIN_PDF_CDF = "WRITE_BASIN_PDF_CDF"
    N_HISTOGRAM_BINS = "N_HISTOGRAM_BINS"
    WRITE_HILLSHADE = "WRITE_HILLSHADE"
    HILLSHADE_AZIMUTH = "HILLSHADE_AZIMUTH"
    HILLSHADE_ALTITUDE = "HILLSHADE_ALTITUDE"
    HILLSHADE_Z_FACTOR = "HILLSHADE_Z_FACTOR"

    OUTPUT_STATS = "OUTPUT_STATS"
    OUTPUT_POLYGONS = "OUTPUT_POLYGONS"
    OUTPUT_AREA_RATIOS = "OUTPUT_AREA_RATIOS"
    OUTPUT_HILLSHADE = "OUTPUT_HILLSHADE"

    VECTOR_FORMAT_OPTIONS = ["flatgeobuf", "gpkg"]
    OUTPUT_FORMAT_OPTIONS = ["cog", "gtiff"]
    SNAPPING_OPTIONS = ["none", "downstream", "expanding_square"]
    BASIN_LENGTH_OPTIONS = ["flow_path", "geometric"]

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
            QgsProcessingParameterVectorLayer(
                self.INPUT_OUTLETS,
                "Basin outlet points (any CRS; gauge_id / gauge_name fields "
                "are passed through)",
                [QgsProcessing.TypeVectorPoint],
            )
        )
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, "Output folder"))
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_basin_stats"))
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.CALCULATE_CHANNEL_STATS,
                "Calculate channel statistics (channel length by stream order)",
                True,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.THRESHOLD_CONTRIBUTING_PIXELS,
                "Source threshold (contributing pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=500,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.CALCULATE_SHAPE_METRICS, "Calculate basin planform shape metrics", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.CALCULATE_COM_BEARING, "Calculate outlet-to-centre-of-mass bearing", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.CALCULATE_JUNCTION_AREA_RATIOS, "Calculate junction drainage-area ratios", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_BASIN_POLYGONS,
                "Write basin polygons with the metrics as attributes",
                True,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.BASIN_POLYGON_FORMAT,
                "Basin polygon format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.REJECT_NODATA_EDGE_BASINS,
                "Reject basins touching nodata or the raster edge",
                True,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.SNAPPING_METHOD,
                "Outlet snapping method",
                options=self.SNAPPING_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
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
                self.THREADS,
                "Threads (0 = auto)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.OUTPUT_FORMAT,
                "Output raster format",
                options=self.OUTPUT_FORMAT_OPTIONS,
                defaultValue=0,
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
            QgsProcessingParameterNumber(
                self.MINIMUM_ELEVATION,
                "Minimum elevation (pixels at or below are set to nodata)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FILL_BEFORE_ROUTING, "Fill pits before routing", True
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FILL_MIN_SLOPE,
                "Fill min slope",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0001,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.BASIN_LENGTH_DEFINITION,
                "Basin length definition",
                options=self.BASIN_LENGTH_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SHAPE_PERIMETER_SIMPLIFY_CELLS,
                "Perimeter simplification tolerance (cells)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.JUNCTION_AREA_RATIO_BOOTSTRAP,
                "Junction area-ratio bootstrap replicates",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=1000,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_BASIN_PDF_CDF, "Write per-basin elevation PDF/CDF CSVs", False
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.N_HISTOGRAM_BINS,
                "Histogram bins for PDF/CDF output",
                type=QgsProcessingParameterNumber.Integer,
                minValue=2,
                defaultValue=100,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(self.WRITE_HILLSHADE, "Write hillshade raster", False)
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HILLSHADE_AZIMUTH,
                "Hillshade azimuth (degrees)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=315.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HILLSHADE_ALTITUDE,
                "Hillshade altitude (degrees)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=45.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HILLSHADE_Z_FACTOR,
                "Hillshade z factor",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
            )
        )

        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_STATS, "Basin statistics table"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_POLYGONS, "Basin polygons"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_AREA_RATIOS, "Junction area ratios"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_HILLSHADE, "Hillshade"))

        attach_backend_help(
            self,
            "lsdtt-basin-stats",
            overrides={
                self.INPUT_OUTLETS: "basin_outlet_csv",
                self.THRESHOLD_CONTRIBUTING_PIXELS: "threshold_contributing_pixels",
                self.CALCULATE_COM_BEARING: "calculate_centre_of_mass_bearing",
                self.BASIN_POLYGON_FORMAT: "basin_polygon_format",
                self.SNAPPING_METHOD: "basin_outlet_snapping_method",
                self.SNAPPING_THRESHOLD: "basin_outlet_contributing_pixels_threshold",
                self.SNAPPING_MIN_RADIUS: "basin_outlet_min_search_radius_pixels",
                self.SNAPPING_MAX_RADIUS: "basin_outlet_max_search_radius_pixels",
                self.WRITE_BASIN_PDF_CDF: "write_basin_pdf_cdf_csv",
                self.JUNCTION_AREA_RATIO_BOOTSTRAP: "junction_area_ratio_bootstrap",
            },
        )

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        outlet_layer = self.parameterAsVectorLayer(parameters, self.INPUT_OUTLETS, context)
        if dem_layer is None:
            raise QgsProcessingException("Input DEM is required.")
        if outlet_layer is None:
            raise QgsProcessingException("Basin outlet point layer is required.")
        read_path, read_fname = split_layer_source(dem_layer.source(), "DEM")

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        if not output_folder or not output_prefix:
            raise QgsProcessingException("Output folder and prefix are required.")
        os.makedirs(output_folder, exist_ok=True)

        # lsdtt-basin-stats reads outlets only as a latitude/longitude CSV.
        outlet_csv = os.path.join(
            output_folder, "{}_qgis_outlets.csv".format(output_prefix)
        )
        export_points_as_latlon_csv(outlet_layer, outlet_csv, context, feedback)

        calculate_channel_stats = self.parameterAsBool(
            parameters, self.CALCULATE_CHANNEL_STATS, context
        )
        calculate_area_ratios = self.parameterAsBool(
            parameters, self.CALCULATE_JUNCTION_AREA_RATIOS, context
        )
        write_polygons = self.parameterAsBool(parameters, self.WRITE_BASIN_POLYGONS, context)
        write_hillshade = self.parameterAsBool(parameters, self.WRITE_HILLSHADE, context)
        if calculate_area_ratios and not calculate_channel_stats:
            raise QgsProcessingException(
                "Junction drainage-area ratios require channel statistics to be enabled."
            )
        polygon_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.BASIN_POLYGON_FORMAT, context)
        ]

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "n_threads": self.parameterAsInt(parameters, self.THREADS, context),
            "write_output_format": self.OUTPUT_FORMAT_OPTIONS[
                self.parameterAsEnum(parameters, self.OUTPUT_FORMAT, context)
            ],
            "basin_outlet_csv": outlet_csv,
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
            "minimum_elevation": self.parameterAsDouble(
                parameters, self.MINIMUM_ELEVATION, context
            ),
            "fill_before_routing": self.parameterAsBool(
                parameters, self.FILL_BEFORE_ROUTING, context
            ),
            "fill_min_slope": self.parameterAsDouble(parameters, self.FILL_MIN_SLOPE, context),
            "reject_nodata_edge_basins": self.parameterAsBool(
                parameters, self.REJECT_NODATA_EDGE_BASINS, context
            ),
            "calculate_channel_stats": calculate_channel_stats,
            "threshold_contributing_pixels": self.parameterAsInt(
                parameters, self.THRESHOLD_CONTRIBUTING_PIXELS, context
            ),
            "calculate_junction_area_ratios": calculate_area_ratios,
            "junction_area_ratio_bootstrap": self.parameterAsInt(
                parameters, self.JUNCTION_AREA_RATIO_BOOTSTRAP, context
            ),
            "calculate_shape_metrics": self.parameterAsBool(
                parameters, self.CALCULATE_SHAPE_METRICS, context
            ),
            "basin_length_definition": self.BASIN_LENGTH_OPTIONS[
                self.parameterAsEnum(parameters, self.BASIN_LENGTH_DEFINITION, context)
            ],
            "shape_perimeter_simplify_cells": self.parameterAsDouble(
                parameters, self.SHAPE_PERIMETER_SIMPLIFY_CELLS, context
            ),
            "calculate_centre_of_mass_bearing": self.parameterAsBool(
                parameters, self.CALCULATE_COM_BEARING, context
            ),
            "write_basin_polygons": write_polygons,
            "basin_polygon_format": polygon_format,
            "write_basin_pdf_cdf_csv": self.parameterAsBool(
                parameters, self.WRITE_BASIN_PDF_CDF, context
            ),
            "n_histogram_bins": self.parameterAsInt(parameters, self.N_HISTOGRAM_BINS, context),
            "write_hillshade_raster": write_hillshade,
            "hillshade_azimuth": self.parameterAsDouble(
                parameters, self.HILLSHADE_AZIMUTH, context
            ),
            "hillshade_altitude": self.parameterAsDouble(
                parameters, self.HILLSHADE_ALTITUDE, context
            ),
            "hillshade_z_factor": self.parameterAsDouble(
                parameters, self.HILLSHADE_Z_FACTOR, context
            ),
        }

        run_lsdtt_tool(
            "lsdtt-basin-stats", param_values, output_folder, output_prefix, feedback
        )

        outputs = {
            self.OUTPUT_STATS: require_output(
                output_folder, output_prefix, "_basin_stats.csv", "basin statistics",
            )
        }
        if write_polygons:
            outputs[self.OUTPUT_POLYGONS] = require_output(
                output_folder, output_prefix,
                "_basin_polygons.{}".format(extension_for_format(polygon_format)),
                "basin polygons",
            )
        if calculate_area_ratios:
            outputs[self.OUTPUT_AREA_RATIOS] = require_output(
                output_folder, output_prefix, "_junction_area_ratios.csv",
                "junction area ratios",
            )
        if write_hillshade:
            outputs[self.OUTPUT_HILLSHADE] = require_output(
                output_folder, output_prefix, "_hillshade.tif", "hillshade",
            )
        return outputs

    def name(self):
        return "basin_stats"

    def displayName(self):
        return "Basin Statistics"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Runs lsdtt-basin-stats: per-basin hydromorphic statistics for "
            "basins draining to supplied outlet points (relief, hypsometry, "
            "channel length by stream order, optional planform shape metrics, "
            "centre-of-mass bearing, and junction drainage-area ratios). The "
            "outlet layer is exported to the latitude/longitude CSV the CLI "
            "expects, and the metrics are written both to "
            "<prefix>_basin_stats.csv and, when basin polygons are enabled, "
            "as attributes on the basin polygon layer. Backend path is "
            "managed in plugin settings."
        )

    def createInstance(self):
        return BasinStatsAlgorithm()
