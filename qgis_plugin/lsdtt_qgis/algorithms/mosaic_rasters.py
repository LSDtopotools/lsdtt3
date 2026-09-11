import os

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputRasterLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterDefinition,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterString,
)

from .support import attach_backend_help, layer_file_path, require_output, run_lsdtt_tool


class MosaicRastersAlgorithm(QgsProcessingAlgorithm):
    INPUT_RASTERS = "INPUT_RASTERS"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    OUTPUT_FORMAT = "OUTPUT_FORMAT"
    ENFORCE_PROJECTION = "ENFORCE_PROJECTION"

    OUTPUT_MOSAIC = "OUTPUT_MOSAIC"

    OUTPUT_FORMAT_OPTIONS = ["cog", "gtiff"]

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterMultipleLayers(
                self.INPUT_RASTERS,
                "Input rasters (later layers win in overlap regions)",
                QgsProcessing.TypeRaster,
            )
        )
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, "Output folder"))
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_mosaic"))
        self.addParameter(
            QgsProcessingParameterEnum(
                self.OUTPUT_FORMAT,
                "Output raster format",
                options=self.OUTPUT_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        proj = QgsProcessingParameterBoolean(
            self.ENFORCE_PROJECTION,
            "Enforce projection (reproject a geographic mosaic afterwards)",
            True,
        )
        proj.setFlags(proj.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(proj)

        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_MOSAIC, "Mosaic raster"))

        attach_backend_help(
            self,
            "lsdtt-raster-preprocessing",
            overrides={self.INPUT_RASTERS: "mosaic_file_list"},
        )

    def processAlgorithm(self, parameters, context, feedback):
        raster_layers = self.parameterAsLayerList(parameters, self.INPUT_RASTERS, context)
        if not raster_layers or len(raster_layers) < 2:
            raise QgsProcessingException("At least two input rasters are required.")

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        if not output_folder or not output_prefix:
            raise QgsProcessingException("Output folder and prefix are required.")
        os.makedirs(output_folder, exist_ok=True)

        # The CLI reads the mosaic inputs from a text file (one path per line).
        # Write it to the output folder so failed runs are easy to debug.
        list_path = os.path.join(output_folder, "{}_mosaic_file_list.txt".format(output_prefix))
        with open(list_path, "w", encoding="utf-8") as dst:
            for layer in raster_layers:
                dst.write(layer_file_path(layer.source()) + "\n")

        param_values = {
            "write_path": output_folder,
            "write_fname": output_prefix,
            "write_output_format": self.OUTPUT_FORMAT_OPTIONS[
                self.parameterAsEnum(parameters, self.OUTPUT_FORMAT, context)
            ],
            "enforce_projection": self.parameterAsBool(parameters, self.ENFORCE_PROJECTION, context),
            # Disable the default sea-level masking: this algorithm only mosaics.
            "minimum_elevation": -99999.0,
            "mosaic_file_list": list_path,
        }

        run_lsdtt_tool(
            "lsdtt-raster-preprocessing", param_values, output_folder, output_prefix, feedback
        )

        return {
            self.OUTPUT_MOSAIC: require_output(
                output_folder, output_prefix, "_mosaic.tif", "mosaic"
            )
        }

    def name(self):
        return "mosaic_rasters"

    def displayName(self):
        return "Mosaic Rasters"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Merges rasters with matching CRS, cell size, and pixel-aligned "
            "origins into a single raster. In overlap regions the last listed "
            "raster wins (lsdtt-raster-preprocessing mosaic workflow). Rasters "
            "with a mismatched CRS or grid are skipped with a warning in the log."
        )

    def createInstance(self):
        return MosaicRastersAlgorithm()
