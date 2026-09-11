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
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
)

from .support import (
    attach_backend_help,
    require_output,
    run_lsdtt_tool,
    split_layer_source,
)

# Derived polyfit metrics, in the order the CLI lists them. The metric name
# is also its output suffix stem: <prefix>_direct_polyfit_<metric>.tif
POLYFIT_METRIC_OPTIONS = [
    "fitted_surface",
    "dzdx",
    "dzdy",
    "d2zdx2",
    "d2zdxdy",
    "d2zdy2",
    "slope",
    "slope_degrees",
    "aspect_radians",
    "aspect_degrees",
    "laplacian",
    "hessian_determinant",
    "mean_curvature",
    "gaussian_curvature",
    "maximum_curvature",
    "minimum_curvature",
    "curvedness",
    "shape_index",
    "profile_curvature",
    "plan_curvature",
    "tangential_curvature",
    "difference_curvature",
    "horizontal_excess_curvature",
    "vertical_excess_curvature",
    "accumulation_curvature",
    "ring_curvature",
    "rotor",
]

# The wrapper pins polyfit_method to the CLI default: the direct fit needs no
# FFTW and is the preferred method, and pinning it keeps output suffixes
# predictable.
POLYFIT_METHOD = "direct"


class SurfaceMetricsAlgorithm(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    WRITE_HILLSHADE = "WRITE_HILLSHADE"
    WRITE_LOCAL_RELIEF = "WRITE_LOCAL_RELIEF"
    LOCAL_RELIEF_RADIUS = "LOCAL_RELIEF_RADIUS"
    POLYFIT_METRICS = "POLYFIT_METRICS"
    POLYFIT_LAMBDA = "POLYFIT_LAMBDA"
    WRITE_SHAPE_CLASSIFICATION = "WRITE_SHAPE_CLASSIFICATION"
    SHAPE_CURVATURE_THRESHOLD = "SHAPE_CURVATURE_THRESHOLD"
    WRITE_MADK2 = "WRITE_MADK2"
    MADK2_WINDOW_RADIUS = "MADK2_WINDOW_RADIUS"
    WRITE_CLOTH_ENVELOPE = "WRITE_CLOTH_ENVELOPE"
    WRITE_CLOTH_INCISION = "WRITE_CLOTH_INCISION"
    WRITE_CLOTH_CREST_MASK = "WRITE_CLOTH_CREST_MASK"

    # Advanced
    THREADS = "THREADS"
    OUTPUT_FORMAT = "OUTPUT_FORMAT"
    HILLSHADE_AZIMUTH = "HILLSHADE_AZIMUTH"
    HILLSHADE_ALTITUDE = "HILLSHADE_ALTITUDE"
    HILLSHADE_Z_FACTOR = "HILLSHADE_Z_FACTOR"
    LOCAL_RELIEF_WINDOW_SHAPE = "LOCAL_RELIEF_WINDOW_SHAPE"
    WRITE_POLYFIT_RASTERS = "WRITE_POLYFIT_RASTERS"
    POLYFIT_Z_FACTOR = "POLYFIT_Z_FACTOR"
    POLYFIT_FLAT_SLOPE_TOLERANCE = "POLYFIT_FLAT_SLOPE_TOLERANCE"
    WRITE_SLOPE_THRESHOLD_LAYERS = "WRITE_SLOPE_THRESHOLD_LAYERS"
    SLOPE_THRESHOLD_LAYERS = "SLOPE_THRESHOLD_LAYERS"
    CLAMP_SLOPE_TO_THRESHOLD = "CLAMP_SLOPE_TO_THRESHOLD"
    SHAPE_CURVATURE_METHOD = "SHAPE_CURVATURE_METHOD"
    WRITE_SHAPE_CURVATURE_RASTERS = "WRITE_SHAPE_CURVATURE_RASTERS"
    MADK2_LAG_PIXELS = "MADK2_LAG_PIXELS"
    CLOTH_RIGIDITY = "CLOTH_RIGIDITY"
    CLOTH_GRAVITY = "CLOTH_GRAVITY"
    CLOTH_CLASS_THRESHOLD = "CLOTH_CLASS_THRESHOLD"
    CLOTH_MAX_ITERATIONS = "CLOTH_MAX_ITERATIONS"
    CLOTH_CONVERGENCE_TOL = "CLOTH_CONVERGENCE_TOL"

    OUTPUT_HILLSHADE = "OUTPUT_HILLSHADE"
    OUTPUT_LOCAL_RELIEF = "OUTPUT_LOCAL_RELIEF"
    OUTPUT_FIRST_POLYFIT_METRIC = "OUTPUT_FIRST_POLYFIT_METRIC"
    OUTPUT_SHAPE_CLASS = "OUTPUT_SHAPE_CLASS"
    OUTPUT_MADK2_RA = "OUTPUT_MADK2_RA"
    OUTPUT_CLOTH_INCISION = "OUTPUT_CLOTH_INCISION"

    OUTPUT_FORMAT_OPTIONS = ["cog", "gtiff"]
    WINDOW_SHAPE_OPTIONS = ["circular", "square"]
    SHAPE_CURVATURE_METHOD_OPTIONS = ["polyfit", "fundamental_form"]

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
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_surface"))
        self.addParameter(
            QgsProcessingParameterBoolean(self.WRITE_HILLSHADE, "Write hillshade raster", True)
        )
        self.addParameter(
            QgsProcessingParameterBoolean(self.WRITE_LOCAL_RELIEF, "Write local relief raster", False)
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.LOCAL_RELIEF_RADIUS,
                "Local relief window radius (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=10.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.POLYFIT_METRICS,
                "Polynomial-fit metrics to write",
                options=POLYFIT_METRIC_OPTIONS,
                allowMultiple=True,
                defaultValue=[],
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.POLYFIT_LAMBDA,
                "Polyfit lengthscale lambda (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=70.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_SHAPE_CLASSIFICATION,
                "Write surface shape classification (SMAP)",
                False,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.SHAPE_CURVATURE_THRESHOLD,
                "Shape classification curvature threshold (kt)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(self.WRITE_MADK2, "Write MADk2 roughness rasters", False)
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MADK2_WINDOW_RADIUS,
                "MADk2 window radius (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=10.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_CLOTH_ENVELOPE, "Write cloth-drape envelope raster", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_CLOTH_INCISION, "Write cloth-drape incision raster", False
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_CLOTH_CREST_MASK, "Write cloth-drape crest mask raster", False
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
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.LOCAL_RELIEF_WINDOW_SHAPE,
                "Local relief window shape",
                options=self.WINDOW_SHAPE_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_POLYFIT_RASTERS,
                "Write the six polyfit coefficient rasters",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.POLYFIT_Z_FACTOR,
                "Polyfit z factor",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.POLYFIT_FLAT_SLOPE_TOLERANCE,
                "Polyfit flat-slope tolerance",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1e-12,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_SLOPE_THRESHOLD_LAYERS, "Write slope threshold layers", False
            )
        )
        self._add_advanced(
            QgsProcessingParameterString(
                self.SLOPE_THRESHOLD_LAYERS,
                "Slope thresholds (degrees, comma separated)",
                "3,7,12",
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.CLAMP_SLOPE_TO_THRESHOLD, "Clamp slope to the threshold value", True
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.SHAPE_CURVATURE_METHOD,
                "Shape classification curvature method",
                options=self.SHAPE_CURVATURE_METHOD_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_SHAPE_CURVATURE_RASTERS,
                "Write the shape Gaussian/mean curvature rasters",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MADK2_LAG_PIXELS,
                "MADk2 lag (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CLOTH_RIGIDITY,
                "Cloth rigidity",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CLOTH_GRAVITY,
                "Cloth gravity (m/step^2)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.3,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CLOTH_CLASS_THRESHOLD,
                "Cloth contact tolerance (m)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.5,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CLOTH_MAX_ITERATIONS,
                "Cloth maximum iterations",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=500,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CLOTH_CONVERGENCE_TOL,
                "Cloth convergence tolerance (m/step)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.005,
            )
        )

        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_HILLSHADE, "Hillshade"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_LOCAL_RELIEF, "Local relief"))
        self.addOutput(
            QgsProcessingOutputRasterLayer(
                self.OUTPUT_FIRST_POLYFIT_METRIC, "First requested polyfit metric"
            )
        )
        self.addOutput(
            QgsProcessingOutputRasterLayer(self.OUTPUT_SHAPE_CLASS, "Surface shape classification")
        )
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_MADK2_RA, "MADk2 roughness (RA)"))
        self.addOutput(
            QgsProcessingOutputRasterLayer(self.OUTPUT_CLOTH_INCISION, "Cloth-drape incision")
        )

        attach_backend_help(
            self,
            "lsdtt-surface-metrics",
            overrides={
                self.LOCAL_RELIEF_RADIUS: "local_relief_radius_metres",
                self.LOCAL_RELIEF_WINDOW_SHAPE: "local_relief_window_shape",
                self.POLYFIT_METRICS: "polyfit_metrics",
                self.POLYFIT_LAMBDA: "polyfit_lambda_metres",
                self.WRITE_SHAPE_CLASSIFICATION: "write_shape_classification",
                self.SHAPE_CURVATURE_METHOD: "shape_curvature_method",
                self.WRITE_MADK2: "write_madk2_rasters",
                self.MADK2_WINDOW_RADIUS: "madk2_window_radius_metres",
                self.WRITE_POLYFIT_RASTERS: "write_polyfit_rasters",
                self.SLOPE_THRESHOLD_LAYERS: "slope_thresholds",
                self.WRITE_SHAPE_CURVATURE_RASTERS: "write_shape_curvature_rasters",
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

        metrics = [
            POLYFIT_METRIC_OPTIONS[index]
            for index in self.parameterAsEnums(parameters, self.POLYFIT_METRICS, context)
        ]
        write_hillshade = self.parameterAsBool(parameters, self.WRITE_HILLSHADE, context)
        write_local_relief = self.parameterAsBool(parameters, self.WRITE_LOCAL_RELIEF, context)
        write_shape_class = self.parameterAsBool(
            parameters, self.WRITE_SHAPE_CLASSIFICATION, context
        )
        write_madk2 = self.parameterAsBool(parameters, self.WRITE_MADK2, context)
        write_polyfit_rasters = self.parameterAsBool(
            parameters, self.WRITE_POLYFIT_RASTERS, context
        )
        write_cloth_envelope = self.parameterAsBool(parameters, self.WRITE_CLOTH_ENVELOPE, context)
        write_cloth_incision = self.parameterAsBool(parameters, self.WRITE_CLOTH_INCISION, context)
        write_cloth_crest = self.parameterAsBool(parameters, self.WRITE_CLOTH_CREST_MASK, context)
        write_slope_thresholds = self.parameterAsBool(
            parameters, self.WRITE_SLOPE_THRESHOLD_LAYERS, context
        )
        if not any([
            write_hillshade, write_local_relief, write_shape_class, write_madk2,
            write_polyfit_rasters, write_cloth_envelope, write_cloth_incision,
            write_cloth_crest, write_slope_thresholds, metrics,
        ]):
            raise QgsProcessingException(
                "Select at least one output: hillshade, local relief, polyfit "
                "metrics, shape classification, MADk2, slope thresholds, or a "
                "cloth-drape raster."
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
            "write_local_relief_raster": write_local_relief,
            "local_relief_radius_metres": self.parameterAsDouble(
                parameters, self.LOCAL_RELIEF_RADIUS, context
            ),
            "local_relief_window_shape": self.WINDOW_SHAPE_OPTIONS[
                self.parameterAsEnum(parameters, self.LOCAL_RELIEF_WINDOW_SHAPE, context)
            ],
            # Polynomial surface fit. The method is pinned to the CLI default
            # so output suffixes stay predictable (see POLYFIT_METHOD).
            "polyfit_method": POLYFIT_METHOD,
            "polyfit_metrics": ",".join(metrics),
            "write_polyfit_rasters": write_polyfit_rasters,
            "polyfit_lambda_metres": self.parameterAsDouble(
                parameters, self.POLYFIT_LAMBDA, context
            ),
            "polyfit_z_factor": self.parameterAsDouble(parameters, self.POLYFIT_Z_FACTOR, context),
            "polyfit_flat_slope_tolerance": self.parameterAsDouble(
                parameters, self.POLYFIT_FLAT_SLOPE_TOLERANCE, context
            ),
            "write_slope_threshold_layers": write_slope_thresholds,
            "slope_thresholds": self.parameterAsString(
                parameters, self.SLOPE_THRESHOLD_LAYERS, context
            ),
            "clamp_slope_to_threshold": self.parameterAsBool(
                parameters, self.CLAMP_SLOPE_TO_THRESHOLD, context
            ),
            # Surface shape classification (SMAP)
            "write_shape_classification": write_shape_class,
            "shape_curvature_method": self.SHAPE_CURVATURE_METHOD_OPTIONS[
                self.parameterAsEnum(parameters, self.SHAPE_CURVATURE_METHOD, context)
            ],
            "shape_curvature_threshold": self.parameterAsDouble(
                parameters, self.SHAPE_CURVATURE_THRESHOLD, context
            ),
            "write_shape_curvature_rasters": self.parameterAsBool(
                parameters, self.WRITE_SHAPE_CURVATURE_RASTERS, context
            ),
            # MADk2 roughness
            "write_madk2_rasters": write_madk2,
            "madk2_window_radius_metres": self.parameterAsDouble(
                parameters, self.MADK2_WINDOW_RADIUS, context
            ),
            "madk2_lag_pixels": self.parameterAsInt(parameters, self.MADK2_LAG_PIXELS, context),
            # Cloth drape
            "write_cloth_envelope_raster": write_cloth_envelope,
            "write_cloth_incision_raster": write_cloth_incision,
            "write_cloth_crest_mask_raster": write_cloth_crest,
            "cloth_rigidity": self.parameterAsDouble(parameters, self.CLOTH_RIGIDITY, context),
            "cloth_gravity": self.parameterAsDouble(parameters, self.CLOTH_GRAVITY, context),
            "cloth_class_threshold": self.parameterAsDouble(
                parameters, self.CLOTH_CLASS_THRESHOLD, context
            ),
            "cloth_max_iterations": self.parameterAsInt(
                parameters, self.CLOTH_MAX_ITERATIONS, context
            ),
            "cloth_convergence_tol": self.parameterAsDouble(
                parameters, self.CLOTH_CONVERGENCE_TOL, context
            ),
        }

        run_lsdtt_tool(
            "lsdtt-surface-metrics", param_values, output_folder, output_prefix, feedback
        )

        outputs = {}
        if write_hillshade:
            outputs[self.OUTPUT_HILLSHADE] = require_output(
                output_folder, output_prefix, "_hillshade.tif", "hillshade",
            )
        if write_local_relief:
            outputs[self.OUTPUT_LOCAL_RELIEF] = require_output(
                output_folder, output_prefix, "_local_relief.tif", "local relief",
            )
        if metrics:
            outputs[self.OUTPUT_FIRST_POLYFIT_METRIC] = require_output(
                output_folder, output_prefix,
                "_{}_polyfit_{}.tif".format(POLYFIT_METHOD, metrics[0]),
                "polyfit {}".format(metrics[0]),
            )
        if write_shape_class:
            outputs[self.OUTPUT_SHAPE_CLASS] = require_output(
                output_folder, output_prefix,
                "_{}_shape_class.tif".format(POLYFIT_METHOD),
                "shape classification",
            )
        if write_madk2:
            outputs[self.OUTPUT_MADK2_RA] = require_output(
                output_folder, output_prefix, "_madk2_ra.tif", "MADk2 roughness",
            )
        if write_cloth_incision:
            outputs[self.OUTPUT_CLOTH_INCISION] = require_output(
                output_folder, output_prefix, "_cloth_incision.tif", "cloth incision",
            )
        return outputs

    def name(self):
        return "surface_metrics"

    def displayName(self):
        return "Surface Metrics"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Runs lsdtt-surface-metrics: DEM-only surface analysis with no "
            "flow routing. Writes hillshade and local relief, any of the "
            "derived polynomial-fit metrics (slope, aspect, the curvature "
            "family) from a direct quadratic fit at a chosen lengthscale, "
            "the SMAP surface-shape classification, MADk2 roughness, and the "
            "cloth-drape envelope/incision/crest rasters. The OpenCV-gated "
            "image-analysis outputs (Frangi, Otsu, Hough, adaptive "
            "thresholding, skeletons) are not exposed here because the "
            "Windows backend bundle is built without OpenCV. Backend path is "
            "managed in plugin settings."
        )

    def createInstance(self):
        return SurfaceMetricsAlgorithm()
