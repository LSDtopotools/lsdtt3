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


class ClipVectorToRasterAlgorithm(OutputOverwriteMixin, QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    INPUT_VECTOR = "INPUT_VECTOR"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    VECTOR_FORMAT = "VECTOR_FORMAT"
    ENFORCE_PROJECTION = "ENFORCE_PROJECTION"

    OUTPUT_CLIPPED = "OUTPUT_CLIPPED"

    VECTOR_FORMAT_OPTIONS = ["flatgeobuf", "gpkg"]

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_DEM,
                "Raster defining the clip extent",
                [QgsProcessing.TypeRaster],
            )
        )
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT_VECTOR,
                "Vector layer to clip",
                [QgsProcessing.TypeVectorAnyGeometry],
            )
        )
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, "Output folder"))
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_vclip"))
        self.add_overwrite_parameter()
        self.addParameter(
            QgsProcessingParameterEnum(
                self.VECTOR_FORMAT,
                "Output vector format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        proj = QgsProcessingParameterBoolean(
            self.ENFORCE_PROJECTION,
            "Enforce projection (reproject geographic rasters before clipping)",
            True,
        )
        proj.setFlags(proj.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(proj)

        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_CLIPPED, "Clipped vector"))

        attach_backend_help(
            self,
            "lsdtt-raster-preprocessing",
            overrides={self.INPUT_VECTOR: "vector_data_fname"},
        )

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        vector_layer = self.parameterAsVectorLayer(parameters, self.INPUT_VECTOR, context)
        if dem_layer is None:
            raise QgsProcessingException("Extent raster is required.")
        if vector_layer is None:
            raise QgsProcessingException("Input vector layer is required.")
        read_path, read_fname = split_layer_source(dem_layer.source(), "raster")

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        if not output_folder or not output_prefix:
            raise QgsProcessingException("Output folder and prefix are required.")
        os.makedirs(output_folder, exist_ok=True)

        vector_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.VECTOR_FORMAT, context)
        ]

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "enforce_projection": self.parameterAsBool(parameters, self.ENFORCE_PROJECTION, context),
            # Disable the default sea-level masking: this algorithm only clips.
            "min_elevation": -99999.0,
            "clip_vector_data_to_raster": True,
            "vector_data_fname": layer_file_path(vector_layer.source()),
            "write_vector_format": vector_format,
        }

        run_lsdtt_tool(
            "lsdtt-raster-preprocessing", param_values, output_folder, output_prefix, feedback
        )

        return {
            self.OUTPUT_CLIPPED: require_output(
                output_folder, output_prefix,
                "_clipped_vector.{}".format(extension_for_format(vector_format)),
                "clipped vector",
            )
        }

    def name(self):
        return "clip_vector_to_raster"

    def displayName(self):
        return "Clip Vector to Raster Extent"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Clips a vector layer to the bounding box of a raster. The vector "
            "CRS need not match the raster CRS (lsdtt-raster-preprocessing "
            "clip_vector_data_to_raster workflow)."
        )

    def createInstance(self):
        return ClipVectorToRasterAlgorithm()
