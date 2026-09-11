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
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
    QgsProcessingParameterVectorLayer,
)

from .support import attach_backend_help, layer_file_path, require_output, run_lsdtt_tool, split_layer_source


class ClipRasterToPolygonAlgorithm(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    INPUT_POLYGON = "INPUT_POLYGON"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    FORCE_RECTANGULAR = "FORCE_RECTANGULAR"
    CLIP_RULE = "CLIP_RULE"
    OUTPUT_FORMAT = "OUTPUT_FORMAT"
    ENFORCE_PROJECTION = "ENFORCE_PROJECTION"

    OUTPUT_CLIPPED = "OUTPUT_CLIPPED"

    CLIP_RULE_OPTIONS = ["majority_coverage", "any_coverage", "full_coverage"]
    OUTPUT_FORMAT_OPTIONS = ["cog", "gtiff"]

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_DEM,
                "Input raster",
                [QgsProcessing.TypeRaster],
            )
        )
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT_POLYGON,
                "Clip polygon (first polygon of first layer is used)",
                [QgsProcessing.TypeVectorPolygon],
            )
        )
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, "Output folder"))
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_clip"))
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.FORCE_RECTANGULAR,
                "Keep all pixels in bounding box (skip polygon-shape masking)",
                False,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.CLIP_RULE,
                "Pixel inclusion rule",
                options=self.CLIP_RULE_OPTIONS,
                defaultValue=0,
            )
        )
        fmt = QgsProcessingParameterEnum(
            self.OUTPUT_FORMAT,
            "Output raster format",
            options=self.OUTPUT_FORMAT_OPTIONS,
            defaultValue=0,
        )
        fmt.setFlags(fmt.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(fmt)
        proj = QgsProcessingParameterBoolean(
            self.ENFORCE_PROJECTION,
            "Enforce projection (reproject geographic rasters before clipping)",
            True,
        )
        proj.setFlags(proj.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(proj)

        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_CLIPPED, "Clipped raster"))

        attach_backend_help(
            self,
            "lsdtt-raster-preprocessing",
            overrides={self.INPUT_POLYGON: "polygon_clip_fname"},
        )

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        polygon_layer = self.parameterAsVectorLayer(parameters, self.INPUT_POLYGON, context)
        if dem_layer is None:
            raise QgsProcessingException("Input raster is required.")
        if polygon_layer is None:
            raise QgsProcessingException("Clip polygon is required.")
        read_path, read_fname = split_layer_source(dem_layer.source(), "raster")

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        if not output_folder or not output_prefix:
            raise QgsProcessingException("Output folder and prefix are required.")
        os.makedirs(output_folder, exist_ok=True)

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "write_output_format": self.OUTPUT_FORMAT_OPTIONS[
                self.parameterAsEnum(parameters, self.OUTPUT_FORMAT, context)
            ],
            "enforce_projection": self.parameterAsBool(parameters, self.ENFORCE_PROJECTION, context),
            # Disable the default sea-level masking: this algorithm only clips.
            "minimum_elevation": -99999.0,
            "clip_raster_to_polygon": True,
            "polygon_clip_fname": layer_file_path(polygon_layer.source()),
            "force_rectangular_raster": self.parameterAsBool(parameters, self.FORCE_RECTANGULAR, context),
            "polygon_clip_rule": self.CLIP_RULE_OPTIONS[
                self.parameterAsEnum(parameters, self.CLIP_RULE, context)
            ],
        }

        run_lsdtt_tool(
            "lsdtt-raster-preprocessing", param_values, output_folder, output_prefix, feedback
        )

        return {
            self.OUTPUT_CLIPPED: require_output(
                output_folder, output_prefix, "_polygon_clipped.tif", "clipped raster"
            )
        }

    def name(self):
        return "clip_raster_to_polygon"

    def displayName(self):
        return "Clip Raster to Polygon"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Crops a raster to the bounding box of a polygon and optionally "
            "masks pixels outside the polygon shape (lsdtt-raster-preprocessing "
            "clip_raster_to_polygon workflow)."
        )

    def createInstance(self):
        return ClipRasterToPolygonAlgorithm()
