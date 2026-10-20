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
    OutputOverwriteMixin,
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


class SurfaceMetricsAlgorithm(OutputOverwriteMixin, QgsProcessingAlgorithm):
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

    # Advanced: compute backend and image-analysis outputs
    COMPUTE_BACKEND = "COMPUTE_BACKEND"
    WRITE_FRANGI = "WRITE_FRANGI"
    WRITE_FRANGI_SCALE = "WRITE_FRANGI_SCALE"
    FRANGI_SIGMA_MIN = "FRANGI_SIGMA_MIN"
    FRANGI_SIGMA_MAX = "FRANGI_SIGMA_MAX"
    FRANGI_N_SCALES = "FRANGI_N_SCALES"
    FRANGI_BETA = "FRANGI_BETA"
    FRANGI_C = "FRANGI_C"
    FRANGI_TARGET = "FRANGI_TARGET"
    WRITE_OTSU = "WRITE_OTSU"
    WRITE_HOUGH_LINES = "WRITE_HOUGH_LINES"
    HOUGH_RHO = "HOUGH_RHO"
    HOUGH_THETA_DEGREES = "HOUGH_THETA_DEGREES"
    HOUGH_THRESHOLD = "HOUGH_THRESHOLD"
    HOUGH_MIN_LINE_LENGTH = "HOUGH_MIN_LINE_LENGTH"
    HOUGH_MAX_LINE_GAP = "HOUGH_MAX_LINE_GAP"
    CANNY_THRESHOLD = "CANNY_THRESHOLD"
    CANNY_RATIO = "CANNY_RATIO"
    CANNY_GAUSSIAN_SIGMA = "CANNY_GAUSSIAN_SIGMA"
    WRITE_CURVATURE_SKELETONS = "WRITE_CURVATURE_SKELETONS"
    CURVATURE_SKELETON_LAMBDA = "CURVATURE_SKELETON_LAMBDA"
    VEC_SKELETON_THRESHOLD = "VEC_SKELETON_THRESHOLD"
    VEC_SKELETON_USE_CDF_THRESHOLD = "VEC_SKELETON_USE_CDF_THRESHOLD"
    VEC_SKELETON_WRITE_CDF = "VEC_SKELETON_WRITE_CDF"
    PROFILE_RELIEF_WINDOW = "PROFILE_RELIEF_WINDOW"
    PROFILE_RELIEF_THRESHOLD = "PROFILE_RELIEF_THRESHOLD"
    CURVATURE_SKELETON_SPUR_MAX = "CURVATURE_SKELETON_SPUR_MAX"
    CURVATURE_SKELETON_MIN_COMPONENT = "CURVATURE_SKELETON_MIN_COMPONENT"
    WRITE_VEC_ADAPTIVE = "WRITE_VEC_ADAPTIVE"
    VEC_ADAPTIVE_BLOCK_SIZE = "VEC_ADAPTIVE_BLOCK_SIZE"
    VEC_ADAPTIVE_C = "VEC_ADAPTIVE_C"
    VEC_ADAPTIVE_GAUSSIAN_SIGMA = "VEC_ADAPTIVE_GAUSSIAN_SIGMA"
    WRITE_VEC_ADAPTIVE_SKELETONS = "WRITE_VEC_ADAPTIVE_SKELETONS"

    OUTPUT_HILLSHADE = "OUTPUT_HILLSHADE"
    OUTPUT_LOCAL_RELIEF = "OUTPUT_LOCAL_RELIEF"
    OUTPUT_FIRST_POLYFIT_METRIC = "OUTPUT_FIRST_POLYFIT_METRIC"
    OUTPUT_SHAPE_CLASS = "OUTPUT_SHAPE_CLASS"
    OUTPUT_MADK2_RA = "OUTPUT_MADK2_RA"
    OUTPUT_CLOTH_INCISION = "OUTPUT_CLOTH_INCISION"
    OUTPUT_FRANGI = "OUTPUT_FRANGI"
    OUTPUT_FRANGI_SCALE = "OUTPUT_FRANGI_SCALE"
    OUTPUT_OTSU = "OUTPUT_OTSU"
    OUTPUT_HOUGH_LINES = "OUTPUT_HOUGH_LINES"
    OUTPUT_VEC_SKELETON = "OUTPUT_VEC_SKELETON"
    OUTPUT_PROFILE_SKELETON = "OUTPUT_PROFILE_SKELETON"
    OUTPUT_VEC_ADAPTIVE_MEAN = "OUTPUT_VEC_ADAPTIVE_MEAN"
    OUTPUT_VEC_ADAPTIVE_GAUSSIAN = "OUTPUT_VEC_ADAPTIVE_GAUSSIAN"
    OUTPUT_VEC_ADAPTIVE_MEAN_SKELETON = "OUTPUT_VEC_ADAPTIVE_MEAN_SKELETON"
    OUTPUT_VEC_ADAPTIVE_GAUSSIAN_SKELETON = "OUTPUT_VEC_ADAPTIVE_GAUSSIAN_SKELETON"

    OUTPUT_FORMAT_OPTIONS = ["cog", "gtiff"]
    WINDOW_SHAPE_OPTIONS = ["circular", "square"]
    SHAPE_CURVATURE_METHOD_OPTIONS = ["polyfit", "fundamental_form"]
    COMPUTE_BACKEND_OPTIONS = ["auto", "cpu", "gpu"]
    FRANGI_TARGET_OPTIONS = ["ridges", "valleys"]

    # Image-analysis outputs registered as QGIS layers: (switch, output id,
    # file suffix, description). The curvature-skeleton and adaptive
    # thresholding runs also write intermediate rasters (and optionally a CDF
    # table) to the output folder; only the skeletons and masks are loaded.
    IMAGE_OUTPUTS = [
        ("write_frangi_raster", OUTPUT_FRANGI, "_frangi.tif", "Frangi response"),
        ("write_frangi_scale_raster", OUTPUT_FRANGI_SCALE, "_frangi_scale.tif",
         "Frangi scale at maximum response"),
        ("write_otsu_raster", OUTPUT_OTSU, "_otsu.tif", "Otsu segmentation"),
        ("write_hough_lines_raster", OUTPUT_HOUGH_LINES, "_hough_lines.tif", "Hough lines"),
        ("write_curvature_skeleton_rasters", OUTPUT_VEC_SKELETON, "_vec_skeleton.tif",
         "VEC skeleton"),
        ("write_curvature_skeleton_rasters", OUTPUT_PROFILE_SKELETON,
         "_profile_skeleton.tif", "profile-curvature relief skeleton"),
        ("write_vec_adaptive_threshold_rasters", OUTPUT_VEC_ADAPTIVE_MEAN,
         "_vec_adaptive_mean.tif", "VEC adaptive mean threshold"),
        ("write_vec_adaptive_threshold_rasters", OUTPUT_VEC_ADAPTIVE_GAUSSIAN,
         "_vec_adaptive_gaussian.tif", "VEC adaptive Gaussian threshold"),
        ("write_vec_adaptive_skeleton_rasters", OUTPUT_VEC_ADAPTIVE_MEAN_SKELETON,
         "_vec_adaptive_mean_skeleton.tif", "VEC adaptive mean skeleton"),
        ("write_vec_adaptive_skeleton_rasters", OUTPUT_VEC_ADAPTIVE_GAUSSIAN_SKELETON,
         "_vec_adaptive_gaussian_skeleton.tif", "VEC adaptive Gaussian skeleton"),
    ]

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
        self.add_overwrite_parameter()
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
                "Threads (0 = all hardware; set a physical-core limit on shared hosts)",
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

        # --- Advanced: compute backend ---
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.COMPUTE_BACKEND,
                "Compute backend for kernels with a GPU port (only the direct polyfit has one)",
                options=self.COMPUTE_BACKEND_OPTIONS,
                defaultValue=0,
            )
        )

        # --- Advanced: Frangi filter ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_FRANGI,
                "Write Frangi filter response raster",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_FRANGI_SCALE,
                "Write Frangi scale-at-maximum-response raster",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FRANGI_SIGMA_MIN,
                "Frangi: minimum scale (metres)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FRANGI_SIGMA_MAX,
                "Frangi: maximum scale (metres)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=5.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FRANGI_N_SCALES,
                "Frangi: number of scales",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=4,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FRANGI_BETA,
                "Frangi: beta",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=0.5,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FRANGI_C,
                "Frangi: c",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.FRANGI_TARGET,
                "Frangi: feature polarity",
                options=self.FRANGI_TARGET_OPTIONS,
                defaultValue=0,
            )
        )

        # --- Advanced: Otsu, Hough and Canny ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_OTSU,
                "Write Otsu threshold segmentation raster",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_HOUGH_LINES,
                "Write Hough straight-line raster (from Canny edges)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HOUGH_RHO,
                "Hough: distance resolution (pixels)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HOUGH_THETA_DEGREES,
                "Hough: angle resolution (degrees)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HOUGH_THRESHOLD,
                "Hough: minimum accumulator votes",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=50,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HOUGH_MIN_LINE_LENGTH,
                "Hough: minimum line length (pixels)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=30.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HOUGH_MAX_LINE_GAP,
                "Hough: maximum gap bridged (pixels)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=10.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CANNY_THRESHOLD,
                "Canny: lower hysteresis threshold (0-255)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                maxValue=255.0,
                defaultValue=50.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CANNY_RATIO,
                "Canny: upper / lower threshold ratio",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=2.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CANNY_GAUSSIAN_SIGMA,
                "Canny: Gaussian blur sigma",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )

        # --- Advanced: curvature skeletons ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_CURVATURE_SKELETONS,
                "Write slope-break skeletons from VEC and profile-curvature relief",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CURVATURE_SKELETON_LAMBDA,
                "Curvature skeletons: polyfit half-wavelength (metres)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=100.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.VEC_SKELETON_THRESHOLD,
                "Curvature skeletons: VEC threshold",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.VEC_SKELETON_USE_CDF_THRESHOLD,
                "Curvature skeletons: choose the VEC threshold by CDF divergence",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.VEC_SKELETON_WRITE_CDF,
                "Curvature skeletons: write the VEC CDF table",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.PROFILE_RELIEF_WINDOW,
                "Curvature skeletons: profile-curvature relief window radius (metres)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=50.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.PROFILE_RELIEF_THRESHOLD,
                "Curvature skeletons: profile-curvature relief threshold",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.1,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CURVATURE_SKELETON_SPUR_MAX,
                "Curvature skeletons: maximum spur pruned (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=10,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.CURVATURE_SKELETON_MIN_COMPONENT,
                "Curvature skeletons: minimum component kept (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1,
            )
        )

        # --- Advanced: adaptive VEC thresholding ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_VEC_ADAPTIVE,
                "Write adaptive (mean and Gaussian) VEC threshold masks",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.VEC_ADAPTIVE_BLOCK_SIZE,
                "Adaptive VEC: neighbourhood block size (pixels, odd)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=3,
                defaultValue=51,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.VEC_ADAPTIVE_C,
                "Adaptive VEC: constant subtracted from the local mean",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.VEC_ADAPTIVE_GAUSSIAN_SIGMA,
                "Adaptive VEC: Gaussian pre-blur sigma (0 = none)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=0.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_VEC_ADAPTIVE_SKELETONS,
                "Adaptive VEC: also thin the masks to skeletons",
                False,
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

        for _, output_id, _, description in self.IMAGE_OUTPUTS:
            self.addOutput(QgsProcessingOutputRasterLayer(output_id, description[0].upper() + description[1:]))

        attach_backend_help(
            self,
            "lsdtt-surface-metrics",
            overrides={
                self.LOCAL_RELIEF_RADIUS: "local_relief_radius_metres",
                self.LOCAL_RELIEF_WINDOW_SHAPE: "local_relief_window_shape",
                self.POLYFIT_METRICS: "polyfit_metrics",
                self.POLYFIT_LAMBDA: "polyfit_lambda_metres",
                self.WRITE_SHAPE_CLASSIFICATION: "write_shape_classification_raster",
                self.SHAPE_CURVATURE_METHOD: "shape_curvature_method",
                self.WRITE_MADK2: "write_madk2_rasters",
                self.MADK2_WINDOW_RADIUS: "madk2_window_radius_metres",
                self.WRITE_POLYFIT_RASTERS: "write_polyfit_rasters",
                self.SLOPE_THRESHOLD_LAYERS: "slope_thresholds",
                self.WRITE_SHAPE_CURVATURE_RASTERS: "write_shape_curvature_rasters",
                self.WRITE_FRANGI: "write_frangi_raster",
                self.WRITE_FRANGI_SCALE: "write_frangi_scale_raster",
                self.FRANGI_SIGMA_MIN: "frangi_sigma_min_metres",
                self.FRANGI_SIGMA_MAX: "frangi_sigma_max_metres",
                self.WRITE_OTSU: "write_otsu_raster",
                self.WRITE_HOUGH_LINES: "write_hough_lines_raster",
                self.WRITE_CURVATURE_SKELETONS: "write_curvature_skeleton_rasters",
                self.CURVATURE_SKELETON_LAMBDA: "curvature_skeleton_lambda_metres",
                self.PROFILE_RELIEF_WINDOW: "profile_relief_window_metres",
                self.CURVATURE_SKELETON_SPUR_MAX: "curvature_skeleton_spur_max_pixels",
                self.CURVATURE_SKELETON_MIN_COMPONENT: "curvature_skeleton_min_component_pixels",
                self.WRITE_VEC_ADAPTIVE: "write_vec_adaptive_threshold_rasters",
                self.VEC_ADAPTIVE_BLOCK_SIZE: "vec_adaptive_block_size_pixels",
                self.WRITE_VEC_ADAPTIVE_SKELETONS: "write_vec_adaptive_skeleton_rasters",
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
        image_switches = {
            "write_frangi_raster": self.WRITE_FRANGI,
            "write_frangi_scale_raster": self.WRITE_FRANGI_SCALE,
            "write_otsu_raster": self.WRITE_OTSU,
            "write_hough_lines_raster": self.WRITE_HOUGH_LINES,
            "write_curvature_skeleton_rasters": self.WRITE_CURVATURE_SKELETONS,
            "write_vec_adaptive_threshold_rasters": self.WRITE_VEC_ADAPTIVE,
            "write_vec_adaptive_skeleton_rasters": self.WRITE_VEC_ADAPTIVE_SKELETONS,
        }
        image_values = {
            key: self.parameterAsBool(parameters, qgis_id, context)
            for key, qgis_id in image_switches.items()
        }
        if (image_values["write_vec_adaptive_skeleton_rasters"]
                and not image_values["write_vec_adaptive_threshold_rasters"]):
            raise QgsProcessingException(
                "Adaptive VEC skeletons thin the adaptive threshold masks; "
                "turn on the adaptive VEC threshold masks as well."
            )
        block_size = self.parameterAsInt(parameters, self.VEC_ADAPTIVE_BLOCK_SIZE, context)
        if image_values["write_vec_adaptive_threshold_rasters"] and block_size % 2 == 0:
            raise QgsProcessingException(
                "The adaptive VEC block size must be an odd number of pixels."
            )
        if not any([
            write_hillshade, write_local_relief, write_shape_class, write_madk2,
            write_polyfit_rasters, write_cloth_envelope, write_cloth_incision,
            write_cloth_crest, write_slope_thresholds, metrics,
        ] + list(image_values.values())):
            raise QgsProcessingException(
                "Select at least one output: hillshade, local relief, polyfit "
                "metrics, shape classification, MADk2, slope thresholds, a "
                "cloth-drape raster, or an image-analysis raster."
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
            "write_shape_classification_raster": write_shape_class,
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
            "cloth_convergence_tolerance_metres": self.parameterAsDouble(
                parameters, self.CLOTH_CONVERGENCE_TOL, context
            ),
            "compute_backend": self.COMPUTE_BACKEND_OPTIONS[
                self.parameterAsEnum(parameters, self.COMPUTE_BACKEND, context)
            ],
            # Frangi filter
            "frangi_sigma_min_metres": self.parameterAsDouble(parameters, self.FRANGI_SIGMA_MIN, context),
            "frangi_sigma_max_metres": self.parameterAsDouble(parameters, self.FRANGI_SIGMA_MAX, context),
            "frangi_n_scales": self.parameterAsInt(parameters, self.FRANGI_N_SCALES, context),
            "frangi_beta": self.parameterAsDouble(parameters, self.FRANGI_BETA, context),
            "frangi_c": self.parameterAsDouble(parameters, self.FRANGI_C, context),
            "frangi_target": self.FRANGI_TARGET_OPTIONS[
                self.parameterAsEnum(parameters, self.FRANGI_TARGET, context)
            ],
            # Hough lines and the Canny edges they are found in
            "hough_rho": self.parameterAsDouble(parameters, self.HOUGH_RHO, context),
            "hough_theta_degrees": self.parameterAsDouble(parameters, self.HOUGH_THETA_DEGREES, context),
            "hough_threshold": self.parameterAsInt(parameters, self.HOUGH_THRESHOLD, context),
            "hough_min_line_length": self.parameterAsDouble(
                parameters, self.HOUGH_MIN_LINE_LENGTH, context
            ),
            "hough_max_line_gap": self.parameterAsDouble(parameters, self.HOUGH_MAX_LINE_GAP, context),
            "canny_threshold": self.parameterAsDouble(parameters, self.CANNY_THRESHOLD, context),
            "canny_ratio": self.parameterAsDouble(parameters, self.CANNY_RATIO, context),
            "canny_gaussian_sigma": self.parameterAsDouble(
                parameters, self.CANNY_GAUSSIAN_SIGMA, context
            ),
            # Curvature skeletons
            "curvature_skeleton_lambda_metres": self.parameterAsDouble(
                parameters, self.CURVATURE_SKELETON_LAMBDA, context
            ),
            "vec_skeleton_threshold": self.parameterAsDouble(
                parameters, self.VEC_SKELETON_THRESHOLD, context
            ),
            "vec_skeleton_use_cdf_threshold": self.parameterAsBool(
                parameters, self.VEC_SKELETON_USE_CDF_THRESHOLD, context
            ),
            "vec_skeleton_write_cdf": self.parameterAsBool(
                parameters, self.VEC_SKELETON_WRITE_CDF, context
            ),
            "profile_relief_window_metres": self.parameterAsDouble(
                parameters, self.PROFILE_RELIEF_WINDOW, context
            ),
            "profile_relief_threshold": self.parameterAsDouble(
                parameters, self.PROFILE_RELIEF_THRESHOLD, context
            ),
            "curvature_skeleton_spur_max_pixels": self.parameterAsInt(
                parameters, self.CURVATURE_SKELETON_SPUR_MAX, context
            ),
            "curvature_skeleton_min_component_pixels": self.parameterAsInt(
                parameters, self.CURVATURE_SKELETON_MIN_COMPONENT, context
            ),
            # Adaptive VEC thresholding
            "vec_adaptive_block_size_pixels": block_size,
            "vec_adaptive_C": self.parameterAsDouble(parameters, self.VEC_ADAPTIVE_C, context),
            "vec_adaptive_gaussian_sigma": self.parameterAsDouble(
                parameters, self.VEC_ADAPTIVE_GAUSSIAN_SIGMA, context
            ),
        }
        param_values.update(image_values)

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
        for switch, output_id, suffix, description in self.IMAGE_OUTPUTS:
            if image_values[switch]:
                outputs[output_id] = require_output(
                    output_folder, output_prefix, suffix, description,
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
            "cloth-drape envelope/incision/crest rasters, and (advanced) the "
            "image-analysis outputs: Frangi filter, Otsu segmentation, Hough "
            "lines, curvature skeletons and adaptive VEC thresholding. Backend "
            "path is managed in plugin settings."
        )

    def createInstance(self):
        return SurfaceMetricsAlgorithm()
