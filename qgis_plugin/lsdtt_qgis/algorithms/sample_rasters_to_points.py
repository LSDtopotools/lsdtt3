import os

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputVectorLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterDefinition,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFile,
    QgsProcessingParameterFolderDestination,
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


class SampleRastersToPointsAlgorithm(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    INPUT_POINTS = "INPUT_POINTS"
    RASTER_ATTRIBUTE_CSV = "RASTER_ATTRIBUTE_CSV"
    LAYER_NAME = "LAYER_NAME"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    POINTS_FORMAT = "POINTS_FORMAT"
    ENFORCE_PROJECTION = "ENFORCE_PROJECTION"

    OUTPUT_POINTS = "OUTPUT_POINTS"

    POINTS_FORMAT_OPTIONS = ["flatgeobuf", "gpkg"]

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_DEM,
                "Reference raster (defines run CRS; sampled rasters listed in CSV)",
                [QgsProcessing.TypeRaster],
            )
        )
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT_POINTS,
                "Input point layer",
                [QgsProcessing.TypeVectorPoint],
            )
        )
        self.addParameter(
            QgsProcessingParameterFile(
                self.RASTER_ATTRIBUTE_CSV,
                "CSV listing rasters to sample (columns: raster,attribute)",
                extension="csv",
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.LAYER_NAME,
                "Point layer name (empty = first layer)",
                "",
                optional=True,
            )
        )
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, "Output folder"))
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_sample"))
        self.addParameter(
            QgsProcessingParameterEnum(
                self.POINTS_FORMAT,
                "Output point format",
                options=self.POINTS_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        proj = QgsProcessingParameterBoolean(
            self.ENFORCE_PROJECTION,
            "Enforce projection (reproject geographic reference raster)",
            True,
        )
        proj.setFlags(proj.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(proj)

        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_POINTS, "Sampled points"))

        attach_backend_help(
            self,
            "lsdtt-raster-preprocessing",
            overrides={
                self.INPUT_POINTS: "point_data_fname",
                self.RASTER_ATTRIBUTE_CSV: "raster_attribute_fname",
                self.LAYER_NAME: "point_data_layer_name",
            },
        )

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        points_layer = self.parameterAsVectorLayer(parameters, self.INPUT_POINTS, context)
        csv_path = self.parameterAsString(parameters, self.RASTER_ATTRIBUTE_CSV, context)
        if dem_layer is None:
            raise QgsProcessingException("Reference raster is required.")
        if points_layer is None:
            raise QgsProcessingException("Input point layer is required.")
        if not csv_path or not os.path.exists(csv_path):
            raise QgsProcessingException("Raster/attribute CSV not found: {}".format(csv_path))
        read_path, read_fname = split_layer_source(dem_layer.source(), "raster")

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        if not output_folder or not output_prefix:
            raise QgsProcessingException("Output folder and prefix are required.")
        os.makedirs(output_folder, exist_ok=True)

        points_format = self.POINTS_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.POINTS_FORMAT, context)
        ]

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "enforce_projection": self.parameterAsBool(parameters, self.ENFORCE_PROJECTION, context),
            # Disable the default sea-level masking: this algorithm only samples.
            "min_elevation": -99999.0,
            "sample_rasters_to_points": True,
            "point_data_fname": layer_file_path(points_layer.source()),
            "raster_attribute_fname": csv_path,
            "write_vector_format": points_format,
        }

        layer_name = self.parameterAsString(parameters, self.LAYER_NAME, context)
        if layer_name:
            param_values["point_data_layer_name"] = layer_name

        run_lsdtt_tool(
            "lsdtt-raster-preprocessing", param_values, output_folder, output_prefix, feedback
        )

        return {
            self.OUTPUT_POINTS: require_output(
                output_folder, output_prefix,
                "_sampled_points.{}".format(extension_for_format(points_format)),
                "sampled points",
            )
        }

    def name(self):
        return "sample_rasters_to_points"

    def displayName(self):
        return "Sample Rasters at Points"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Samples a list of rasters (from a raster,attribute CSV) at each "
            "point of a point layer and writes a new point layer with the "
            "sampled values appended as attributes "
            "(lsdtt-raster-preprocessing sample_rasters_to_points workflow)."
        )

    def createInstance(self):
        return SampleRastersToPointsAlgorithm()
