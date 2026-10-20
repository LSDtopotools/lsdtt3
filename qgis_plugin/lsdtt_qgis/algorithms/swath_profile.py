import os

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputFile,
    QgsProcessingOutputRasterLayer,
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
    export_baseline_line,
    layer_file_path,
    require_output,
    run_lsdtt_tool,
    split_layer_source,
)


class SwathProfileAlgorithm(OutputOverwriteMixin, QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    INPUT_BASELINE = "INPUT_BASELINE"
    VALUE_RASTER = "VALUE_RASTER"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    HALF_WIDTH_METRES = "HALF_WIDTH_METRES"
    BIN_WIDTH_METRES = "BIN_WIDTH_METRES"
    WRITE_DISTANCE_RASTERS = "WRITE_DISTANCE_RASTERS"

    # Advanced
    THREADS = "THREADS"
    OUTPUT_FORMAT = "OUTPUT_FORMAT"

    OUTPUT_PROFILE = "OUTPUT_PROFILE"
    OUTPUT_ALONG_AXIS = "OUTPUT_ALONG_AXIS"
    OUTPUT_PERP_DISTANCE = "OUTPUT_PERP_DISTANCE"
    OUTPUT_SIGNED_PERP = "OUTPUT_SIGNED_PERP"

    OUTPUT_FORMAT_OPTIONS = ["cog", "gtiff"]

    # (output id, CLI suffix, description) for the write_distance_rasters set.
    DISTANCE_RASTERS = (
        (OUTPUT_ALONG_AXIS, "_swath_along_axis.tif", "along-axis distance"),
        (OUTPUT_PERP_DISTANCE, "_swath_perp_distance.tif", "perpendicular distance"),
        (OUTPUT_SIGNED_PERP, "_swath_signed_perp.tif", "signed perpendicular distance"),
    )

    def _add_advanced(self, param):
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_DEM,
                "Input DEM (projected, metres)",
                [QgsProcessing.TypeRaster],
            )
        )
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT_BASELINE,
                "Swath baseline line (any CRS; first selected feature, or the "
                "first feature if none is selected)",
                [QgsProcessing.TypeVectorLine],
            )
        )
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.VALUE_RASTER,
                "Raster to profile (optional; must share the DEM's grid; "
                "default is the DEM's elevation)",
                optional=True,
            )
        )
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, "Output folder"))
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_swath"))
        self.add_overwrite_parameter()
        self.addParameter(
            QgsProcessingParameterNumber(
                self.HALF_WIDTH_METRES,
                "Swath half-width (m; 0 = every valid pixel)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.BIN_WIDTH_METRES,
                "Along-axis bin width (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.001,
                defaultValue=1000.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_DISTANCE_RASTERS,
                "Write along-axis and perpendicular distance rasters",
                True,
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

        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_PROFILE, "Swath profile table"))
        for output_id, _, description in self.DISTANCE_RASTERS:
            self.addOutput(
                QgsProcessingOutputRasterLayer(output_id, description.capitalize())
            )

        attach_backend_help(
            self,
            "lsdtt-swath-profile",
            overrides={
                self.INPUT_BASELINE: "baseline_vector_fname",
                self.VALUE_RASTER: "value_raster_fname",
            },
        )

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        baseline_layer = self.parameterAsVectorLayer(parameters, self.INPUT_BASELINE, context)
        value_layer = self.parameterAsRasterLayer(parameters, self.VALUE_RASTER, context)
        if dem_layer is None:
            raise QgsProcessingException("Input DEM is required.")
        if baseline_layer is None:
            raise QgsProcessingException("Swath baseline layer is required.")
        read_path, read_fname = split_layer_source(dem_layer.source(), "DEM")

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        if not output_folder or not output_prefix:
            raise QgsProcessingException("Output folder and prefix are required.")
        os.makedirs(output_folder, exist_ok=True)

        # lsdtt-swath-profile reads baseline coordinates as DEM-CRS metres and
        # does not reproject, so export the chosen line in the DEM's CRS.
        baseline_path = os.path.join(
            output_folder, "{}_qgis_baseline.geojson".format(output_prefix)
        )
        export_baseline_line(
            baseline_layer, dem_layer.crs(), baseline_path, context, feedback
        )

        write_distance_rasters = self.parameterAsBool(
            parameters, self.WRITE_DISTANCE_RASTERS, context
        )
        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "n_threads": self.parameterAsInt(parameters, self.THREADS, context),
            "write_output_format": self.OUTPUT_FORMAT_OPTIONS[
                self.parameterAsEnum(parameters, self.OUTPUT_FORMAT, context)
            ],
            "baseline_vector_fname": baseline_path,
            "value_raster_fname": (
                layer_file_path(value_layer.source()) if value_layer is not None else "NULL"
            ),
            "half_width_metres": self.parameterAsDouble(
                parameters, self.HALF_WIDTH_METRES, context
            ),
            "bin_width_metres": self.parameterAsDouble(
                parameters, self.BIN_WIDTH_METRES, context
            ),
            "write_distance_rasters": write_distance_rasters,
        }

        run_lsdtt_tool(
            "lsdtt-swath-profile", param_values, output_folder, output_prefix, feedback
        )

        outputs = {
            self.OUTPUT_PROFILE: require_output(
                output_folder, output_prefix, "_swath_profile.csv", "swath profile table",
            )
        }
        if write_distance_rasters:
            for output_id, suffix, description in self.DISTANCE_RASTERS:
                outputs[output_id] = require_output(
                    output_folder, output_prefix, suffix, description + " raster",
                )
        return outputs

    def name(self):
        return "swath_profile"

    def displayName(self):
        return "Swath Profile"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Runs lsdtt-swath-profile: every valid pixel within the half-width "
            "of a baseline line is projected onto it and binned by distance "
            "along it. Each bin reports the pixel count, the mean and 21 "
            "percentiles (p0, p5, ... p100) of the DEM elevation, or of another "
            "raster on the same grid, in <prefix>_swath_profile.csv. The "
            "baseline can be in any CRS: it is reprojected to the DEM's CRS "
            "and saved as <prefix>_qgis_baseline.geojson. Optional rasters "
            "give each pixel's distance along and across the baseline. "
            "Backend path is managed in plugin settings."
        )

    def createInstance(self):
        return SwathProfileAlgorithm()
