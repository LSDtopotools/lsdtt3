import os

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
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
    OutputOverwriteMixin,
    attach_backend_help,
    extension_for_format,
    layer_file_path,
    require_output,
    run_lsdtt_tool,
    split_layer_source,
)


class RasterPreprocessingAlgorithm(OutputOverwriteMixin, QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    OUTPUT_FORMAT = "OUTPUT_FORMAT"
    ENFORCE_PROJECTION = "ENFORCE_PROJECTION"
    MINIMUM_ELEVATION = "MINIMUM_ELEVATION"
    MAXIMUM_ELEVATION = "MAXIMUM_ELEVATION"
    WRITE_FILL = "WRITE_FILL"
    WRITE_BREACHED = "WRITE_BREACHED"
    WRITE_HILLSHADE = "WRITE_HILLSHADE"

    # Advanced: threading and reprojection
    THREADS = "THREADS"
    RESAMPLING_METHOD = "RESAMPLING_METHOD"
    GRID_SPACING = "GRID_SPACING"
    FORCE_REPROJECT = "FORCE_REPROJECT"
    FORCE_RESAMPLING = "FORCE_RESAMPLING"
    FORCED_GRID_SPACING = "FORCED_GRID_SPACING"
    ALLOW_LAMBERT = "ALLOW_LAMBERT"
    WARP_MEM_MB = "WARP_MEM_MB"

    # Advanced: masking and nodata repair
    BINARY_MASK = "BINARY_MASK"
    WRITE_NODATA_FILLED = "WRITE_NODATA_FILLED"
    NODATA_FILL_METHOD = "NODATA_FILL_METHOD"
    NODATA_FILL_MAX_ITER = "NODATA_FILL_MAX_ITER"
    WRITE_CLIPPED_RECT = "WRITE_CLIPPED_RECT"

    # Advanced: breaching
    BREACH_BEFORE_FILL = "BREACH_BEFORE_FILL"
    FILL_REMAINDER_AFTER_BREACH = "FILL_REMAINDER_AFTER_BREACH"
    BREACH_ALGORITHM = "BREACH_ALGORITHM"
    BREACH_LC_MAX_DIST = "BREACH_LC_MAX_DIST"
    BREACH_LC_MAX_COST = "BREACH_LC_MAX_COST"
    BREACH_LC_FLAT_INCREMENT = "BREACH_LC_FLAT_INCREMENT"
    BREACH_LC_MINIMIZE_DIST = "BREACH_LC_MINIMIZE_DIST"
    BREACH_LC_EDGE_POLICY = "BREACH_LC_EDGE_POLICY"
    WRITE_CARVE_DIFF = "WRITE_CARVE_DIFF"

    # Advanced: filling
    FILL_ALGORITHM = "FILL_ALGORITHM"
    FILL_MIN_SLOPE = "FILL_MIN_SLOPE"
    WRITE_FILL_DIFF = "WRITE_FILL_DIFF"

    # Advanced: flat resolution and fill-spill-merge
    PRIORITY_FLOOD_QUEUE = "PRIORITY_FLOOD_QUEUE"
    FLAT_RESOLUTION = "FLAT_RESOLUTION"
    FSM_RAINFALL_RATE = "FSM_RAINFALL_RATE"
    FSM_BOUNDARY_CONDITIONS = "FSM_BOUNDARY_CONDITIONS"

    # Advanced: vertical units
    VERTICAL_UNITS = "VERTICAL_UNITS"
    CHECK_VERTICAL_UNITS_USGS = "CHECK_VERTICAL_UNITS_USGS"

    # Advanced: edge-preserving smoothing (applied before breaching and filling)
    SMOOTH_ALGORITHM = "SMOOTH_ALGORITHM"
    PM_ITERATIONS = "PM_ITERATIONS"
    PM_LAMBDA = "PM_LAMBDA"
    PM_CONDUCTANCE = "PM_CONDUCTANCE"
    PM_K_QUANTILE = "PM_K_QUANTILE"
    PM_K_METRES_PER_CELL = "PM_K_METRES_PER_CELL"
    WRITE_STENCIL_DIFF = "WRITE_STENCIL_DIFF"
    COHERENCE_ITERATIONS = "COHERENCE_ITERATIONS"
    COHERENCE_TIME_STEP = "COHERENCE_TIME_STEP"
    COHERENCE_DERIVATIVE_RADIUS = "COHERENCE_DERIVATIVE_RADIUS"
    COHERENCE_INTEGRATION_RADIUS = "COHERENCE_INTEGRATION_RADIUS"
    COHERENCE_CONTRAST = "COHERENCE_CONTRAST"
    COHERENCE_MIN_DIFFUSIVITY = "COHERENCE_MIN_DIFFUSIVITY"

    # Advanced: pixel replacement
    REPLACE_PIXELS = "REPLACE_PIXELS"
    REPLACE_PIXELS_VALUE_FIELD = "REPLACE_PIXELS_VALUE_FIELD"
    WRITE_REPLACED_PIXELS = "WRITE_REPLACED_PIXELS"

    # Advanced: sink candidates
    WRITE_SINK_CANDIDATES = "WRITE_SINK_CANDIDATES"
    SINK_CANDIDATE_MIN_PIXELS = "SINK_CANDIDATE_MIN_PIXELS"

    # Advanced: depression polygons
    WRITE_DEPRESSION_POLYGONS = "WRITE_DEPRESSION_POLYGONS"
    DEPRESSION_MIN_PIXELS = "DEPRESSION_MIN_PIXELS"
    DEPRESSION_MIN_DEPTH = "DEPRESSION_MIN_DEPTH"
    DEPRESSION_FORMAT = "DEPRESSION_FORMAT"

    # Advanced: hillshade
    HILLSHADE_AZIMUTH = "HILLSHADE_AZIMUTH"
    HILLSHADE_ALTITUDE = "HILLSHADE_ALTITUDE"
    HILLSHADE_Z_FACTOR = "HILLSHADE_Z_FACTOR"
    WRITE_MD_HILLSHADE = "WRITE_MD_HILLSHADE"
    MD_N_DIRECTIONS = "MD_N_DIRECTIONS"
    MD_START_AZIMUTH = "MD_START_AZIMUTH"

    OUTPUT_FILL = "OUTPUT_FILL"
    OUTPUT_BREACHED = "OUTPUT_BREACHED"
    OUTPUT_HILLSHADE = "OUTPUT_HILLSHADE"
    OUTPUT_MD_HILLSHADE = "OUTPUT_MD_HILLSHADE"
    OUTPUT_NODATA_FILLED = "OUTPUT_NODATA_FILLED"
    OUTPUT_CLIPPED_RECT = "OUTPUT_CLIPPED_RECT"
    OUTPUT_FILL_DIFF = "OUTPUT_FILL_DIFF"
    OUTPUT_CARVE_DIFF = "OUTPUT_CARVE_DIFF"
    OUTPUT_BINARY_MASKED = "OUTPUT_BINARY_MASKED"
    OUTPUT_DEPRESSIONS = "OUTPUT_DEPRESSIONS"
    OUTPUT_SMOOTHED = "OUTPUT_SMOOTHED"
    OUTPUT_STENCIL_DIFF = "OUTPUT_STENCIL_DIFF"
    OUTPUT_REPLACED_PIXELS = "OUTPUT_REPLACED_PIXELS"
    OUTPUT_SINK_CANDIDATES = "OUTPUT_SINK_CANDIDATES"

    OUTPUT_FORMAT_OPTIONS = ["cog", "gtiff"]
    RESAMPLING_OPTIONS = [
        "cubic", "near", "bilinear", "cubicspline", "lanczos",
        "average", "mode", "max", "min", "med", "q1", "q3", "sum",
    ]
    NODATA_FILL_OPTIONS = ["laplacian", "iterative_weighted_average"]
    BREACH_ALGORITHM_OPTIONS = ["lindsay2016", "least_cost"]
    BREACH_LC_EDGE_POLICY_OPTIONS = ["outlet", "whitebox"]
    FILL_ALGORITHM_OPTIONS = ["barnes2014", "barnesetal2021"]
    VECTOR_FORMAT_OPTIONS = ["flatgeobuf", "gpkg"]
    PRIORITY_FLOOD_QUEUE_OPTIONS = ["binary", "radix"]
    FLAT_RESOLUTION_OPTIONS = ["inline_gradient", "barnes2014b_gradient", "barnes2014b_float_step"]
    VERTICAL_UNITS_OPTIONS = ["metres", "feet", "us_survey_feet"]
    SMOOTH_ALGORITHM_OPTIONS = ["none", "perona_malik", "coherence_diffusion"]
    PM_CONDUCTANCE_OPTIONS = ["exponential", "quadratic"]

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
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_preprocess"))
        self.add_overwrite_parameter()
        self.addParameter(
            QgsProcessingParameterEnum(
                self.OUTPUT_FORMAT,
                "Output raster format",
                options=self.OUTPUT_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.ENFORCE_PROJECTION, "Enforce projection", True))
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MINIMUM_ELEVATION,
                "Minimum elevation (pixels at or below become nodata; negative disables)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAXIMUM_ELEVATION,
                "Maximum elevation (pixels above become nodata)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=30000.0,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_FILL, "Write fill raster", True))
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_BREACHED, "Write breached raster", False))
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_HILLSHADE, "Write hillshade raster", False))

        # --- Advanced: threading and reprojection ---
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
                self.RESAMPLING_METHOD,
                "Reprojection resampling method",
                options=self.RESAMPLING_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.GRID_SPACING,
                "Reprojection grid spacing in metres (negative = auto)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=-1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FORCE_REPROJECT,
                "Force reprojection to the recommended analysis CRS",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FORCE_RESAMPLING,
                "Force a user-specified grid spacing",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FORCED_GRID_SPACING,
                "Forced grid spacing in metres (used when forcing resampling)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=-1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.ALLOW_LAMBERT,
                "Allow Lambert Azimuthal Equal Area for very wide rasters",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.WARP_MEM_MB,
                "GDALWarp working memory (MB)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=64,
                defaultValue=512,
            )
        )

        # --- Advanced: masking and nodata repair ---
        self._add_advanced(
            QgsProcessingParameterRasterLayer(
                self.BINARY_MASK,
                "Binary mask raster (pixels > 0.5 kept; others set to nodata)",
                optional=True,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_NODATA_FILLED,
                "Write raster with interior nodata filled",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.NODATA_FILL_METHOD,
                "Interior nodata fill method",
                options=self.NODATA_FILL_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.NODATA_FILL_MAX_ITER,
                "Maximum iterations for nodata filling",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=10000,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_CLIPPED_RECT,
                "Write raster clipped to largest all-valid rectangle",
                False,
            )
        )

        # --- Advanced: breaching ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.BREACH_BEFORE_FILL,
                "Breach depressions before filling pits",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.FILL_REMAINDER_AFTER_BREACH,
                "Fill pits that could not be fully breached",
                True,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.BREACH_ALGORITHM,
                "Breach algorithm",
                options=self.BREACH_ALGORITHM_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.BREACH_LC_MAX_DIST,
                "Least-cost breaching: maximum breach path length (cells; "
                "0 = edge-policy default: 1000 outlet, 100 whitebox)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                maxValue=32767,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.BREACH_LC_MAX_COST,
                "Least-cost breaching: maximum total lowering (negative = unlimited)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=-1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.BREACH_LC_FLAT_INCREMENT,
                "Least-cost breaching: elevation drop per cell along a breach (0 = auto)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.BREACH_LC_MINIMIZE_DIST,
                "Least-cost breaching: weight path cost by length",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.BREACH_LC_EDGE_POLICY,
                "Least-cost breaching: raster edge policy",
                options=self.BREACH_LC_EDGE_POLICY_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_CARVE_DIFF,
                "Write carve difference raster (breached minus original)",
                False,
            )
        )

        # --- Advanced: filling ---
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.FILL_ALGORITHM,
                "Fill algorithm",
                options=self.FILL_ALGORITHM_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FILL_MIN_SLOPE,
                "Fill minimum slope",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0001,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_FILL_DIFF,
                "Write fill difference raster (filled minus original)",
                False,
            )
        )

        self._add_advanced(
            QgsProcessingParameterEnum(
                self.PRIORITY_FLOOD_QUEUE,
                "Barnes 2014 fill: priority-flood queue",
                options=self.PRIORITY_FLOOD_QUEUE_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.FLAT_RESOLUTION,
                "Flat resolution (how filled flats are made drainable)",
                options=self.FLAT_RESOLUTION_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FSM_RAINFALL_RATE,
                "Fill-spill-merge: uniform rainfall depth per cell (elevation units)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=0.1,
            )
        )
        self._add_advanced(
            QgsProcessingParameterString(
                self.FSM_BOUNDARY_CONDITIONS,
                "Fill-spill-merge: boundary conditions N E S W (o open, b baselevel, n no-flux, p periodic)",
                "oooo",
            )
        )

        # --- Advanced: vertical units ---
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.VERTICAL_UNITS,
                "Vertical units of the input elevations",
                options=self.VERTICAL_UNITS_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.CHECK_VERTICAL_UNITS_USGS,
                "If the DEM looks like feet, check against USGS 3DEP (sends 5 coordinates)",
                True,
            )
        )

        # --- Advanced: edge-preserving smoothing ---
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.SMOOTH_ALGORITHM,
                "Edge-preserving smoothing before breaching and filling",
                options=self.SMOOTH_ALGORITHM_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.PM_ITERATIONS,
                "Perona-Malik: iterations",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=10,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.PM_LAMBDA,
                "Perona-Malik: timestep lambda (0 < lambda <= 0.25)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                maxValue=0.25,
                defaultValue=0.2,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.PM_CONDUCTANCE,
                "Perona-Malik: conductance function",
                options=self.PM_CONDUCTANCE_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.PM_K_QUANTILE,
                "Perona-Malik: edge threshold K as a quantile of neighbour differences (0 = state K below)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                maxValue=1.0,
                defaultValue=0.9,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.PM_K_METRES_PER_CELL,
                "Perona-Malik: edge threshold K in metres per cell (used when the quantile is 0)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=0.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_STENCIL_DIFF,
                "Perona-Malik: write the neighbour-difference raster K was read from",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.COHERENCE_ITERATIONS,
                "Coherence diffusion: iterations",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=10,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.COHERENCE_TIME_STEP,
                "Coherence diffusion: timestep tau (0 < tau <= 0.25)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                maxValue=0.25,
                defaultValue=0.15,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.COHERENCE_DERIVATIVE_RADIUS,
                "Coherence diffusion: gradient pre-smoothing radius in metres (0 = one cell)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=0.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.COHERENCE_INTEGRATION_RADIUS,
                "Coherence diffusion: structure-tensor integration radius in metres (0 = three cells)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=0.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.COHERENCE_CONTRAST,
                "Coherence diffusion: contrast C",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=0.2,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.COHERENCE_MIN_DIFFUSIVITY,
                "Coherence diffusion: across-feature diffusivity floor (0 to 1)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                maxValue=1.0,
                defaultValue=0.001,
            )
        )

        # --- Advanced: pixel replacement ---
        self._add_advanced(
            QgsProcessingParameterVectorLayer(
                self.REPLACE_PIXELS,
                "Pixel replacement: point layer of pixels to overwrite",
                [QgsProcessing.TypeVectorPoint],
                optional=True,
            )
        )
        self._add_advanced(
            QgsProcessingParameterString(
                self.REPLACE_PIXELS_VALUE_FIELD,
                "Pixel replacement: field holding the new value",
                "value",
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_REPLACED_PIXELS,
                "Pixel replacement: write the DEM after replacement",
                False,
            )
        )

        # --- Advanced: sink candidates ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_SINK_CANDIDATES,
                "Write candidate internally drained sinks (point layer)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SINK_CANDIDATE_MIN_PIXELS,
                "Sink candidates: minimum sink size (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1,
            )
        )

        # --- Advanced: depression polygons ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_DEPRESSION_POLYGONS,
                "Write depression polygons (vector layer with area and volume)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.DEPRESSION_MIN_PIXELS,
                "Minimum depression size (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.DEPRESSION_MIN_DEPTH,
                "Minimum depression depth (metres)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=0.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.DEPRESSION_FORMAT,
                "Depression polygons format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )

        # --- Advanced: hillshade ---
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HILLSHADE_AZIMUTH,
                "Hillshade azimuth (degrees clockwise from north)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                maxValue=360.0,
                defaultValue=315.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HILLSHADE_ALTITUDE,
                "Hillshade sun altitude (degrees)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                maxValue=90.0,
                defaultValue=45.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.HILLSHADE_Z_FACTOR,
                "Hillshade vertical exaggeration",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_MD_HILLSHADE,
                "Write multidirectional hillshade raster",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MD_N_DIRECTIONS,
                "Multidirectional hillshade: number of directions",
                type=QgsProcessingParameterNumber.Integer,
                minValue=2,
                defaultValue=8,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MD_START_AZIMUTH,
                "Multidirectional hillshade: starting azimuth (degrees)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                maxValue=360.0,
                defaultValue=0.0,
            )
        )

        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_FILL, "Fill raster"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_BREACHED, "Breached raster"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_HILLSHADE, "Hillshade raster"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_MD_HILLSHADE, "Multidirectional hillshade raster"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_NODATA_FILLED, "Nodata-filled raster"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_CLIPPED_RECT, "Clipped rectangle raster"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_FILL_DIFF, "Fill difference raster"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_CARVE_DIFF, "Carve difference raster"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_BINARY_MASKED, "Binary-masked raster"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_DEPRESSIONS, "Depression polygons"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_SMOOTHED, "Smoothed DEM"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_STENCIL_DIFF, "Neighbour-difference raster"))
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT_REPLACED_PIXELS, "DEM after pixel replacement"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_SINK_CANDIDATES, "Sink candidates"))

        attach_backend_help(
            self,
            "lsdtt-raster-preprocessing",
            overrides={
                self.ALLOW_LAMBERT: "allow_lambert_projection",
                self.WRITE_CLIPPED_RECT: "write_clipped_rectangle_raster",
                self.WRITE_MD_HILLSHADE: "write_multidirectional_hillshade_raster",
                self.WRITE_NODATA_FILLED: "write_filled_nodata_raster",
                self.BINARY_MASK: "binary_mask_fname",
                self.DEPRESSION_FORMAT: "write_vector_format",
                self.MD_N_DIRECTIONS: "multidirectional_hillshade_n_directions",
                self.MD_START_AZIMUTH: "multidirectional_hillshade_start_azimuth",
                self.FLAT_RESOLUTION: "fill_flat_resolution_choice",
                self.FSM_RAINFALL_RATE: "fill_spill_merge_rainfall_rate",
                self.CHECK_VERTICAL_UNITS_USGS: "check_vertical_units_against_usgs",
                self.SMOOTH_ALGORITHM: "smooth_algorithm_choice",
                self.PM_ITERATIONS: "perona_malik_iterations",
                self.PM_LAMBDA: "perona_malik_lambda",
                self.PM_CONDUCTANCE: "perona_malik_conductance",
                self.PM_K_QUANTILE: "perona_malik_k_quantile",
                self.PM_K_METRES_PER_CELL: "perona_malik_k_metres_per_cell",
                self.WRITE_STENCIL_DIFF: "write_stencil_difference_raster",
                self.COHERENCE_DERIVATIVE_RADIUS: "coherence_derivative_radius_metres",
                self.COHERENCE_INTEGRATION_RADIUS: "coherence_integration_radius_metres",
                self.REPLACE_PIXELS: "replace_pixels_fname",
                self.WRITE_REPLACED_PIXELS: "write_replaced_pixels_raster",
            },
        )

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        if dem_layer is None:
            raise QgsProcessingException("Input DEM is required.")
        read_path, read_fname = split_layer_source(dem_layer.source(), "DEM")

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        if not output_folder:
            raise QgsProcessingException("Output folder is required.")
        if not output_prefix:
            raise QgsProcessingException("Output prefix is required.")
        os.makedirs(output_folder, exist_ok=True)

        output_format = self.OUTPUT_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.OUTPUT_FORMAT, context)
        ]

        write_fill = self.parameterAsBool(parameters, self.WRITE_FILL, context)
        write_breached = self.parameterAsBool(parameters, self.WRITE_BREACHED, context)
        write_hillshade = self.parameterAsBool(parameters, self.WRITE_HILLSHADE, context)
        write_md_hillshade = self.parameterAsBool(parameters, self.WRITE_MD_HILLSHADE, context)
        write_nodata_filled = self.parameterAsBool(parameters, self.WRITE_NODATA_FILLED, context)
        write_clipped_rect = self.parameterAsBool(parameters, self.WRITE_CLIPPED_RECT, context)
        write_fill_diff = self.parameterAsBool(parameters, self.WRITE_FILL_DIFF, context)
        write_carve_diff = self.parameterAsBool(parameters, self.WRITE_CARVE_DIFF, context)
        write_depressions = self.parameterAsBool(parameters, self.WRITE_DEPRESSION_POLYGONS, context)
        breach_before_fill = self.parameterAsBool(parameters, self.BREACH_BEFORE_FILL, context)

        force_resampling = self.parameterAsBool(parameters, self.FORCE_RESAMPLING, context)
        forced_grid_spacing = self.parameterAsDouble(parameters, self.FORCED_GRID_SPACING, context)
        if force_resampling and forced_grid_spacing <= 0.0:
            raise QgsProcessingException(
                "Forced grid spacing must be positive when forcing resampling."
            )

        mask_layer = self.parameterAsRasterLayer(parameters, self.BINARY_MASK, context)

        smooth_algorithm = self.SMOOTH_ALGORITHM_OPTIONS[
            self.parameterAsEnum(parameters, self.SMOOTH_ALGORITHM, context)
        ]
        write_stencil_diff = self.parameterAsBool(parameters, self.WRITE_STENCIL_DIFF, context)
        if write_stencil_diff and smooth_algorithm != "perona_malik":
            raise QgsProcessingException(
                "The neighbour-difference raster is written by Perona-Malik smoothing; "
                "choose perona_malik or turn the option off."
            )
        fsm_boundary_conditions = self.parameterAsString(
            parameters, self.FSM_BOUNDARY_CONDITIONS, context
        ).strip().lower()
        if len(fsm_boundary_conditions) != 4 or any(c not in "obnp" for c in fsm_boundary_conditions):
            raise QgsProcessingException(
                "Fill-spill-merge boundary conditions must be four letters from o, b, n, p "
                "(north, east, south, west), e.g. oooo."
            )
        replace_layer = self.parameterAsVectorLayer(parameters, self.REPLACE_PIXELS, context)
        write_replaced = self.parameterAsBool(parameters, self.WRITE_REPLACED_PIXELS, context)
        if write_replaced and replace_layer is None:
            raise QgsProcessingException(
                "Writing the replaced-pixel DEM needs a pixel-replacement point layer."
            )
        write_sink_candidates = self.parameterAsBool(parameters, self.WRITE_SINK_CANDIDATES, context)

        depression_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.DEPRESSION_FORMAT, context)
        ]

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "write_output_format": output_format,
            "write_vector_format": depression_format,
            "n_threads": self.parameterAsInt(parameters, self.THREADS, context),
            # projection
            "enforce_projection": self.parameterAsBool(parameters, self.ENFORCE_PROJECTION, context),
            "resampling_method": self.RESAMPLING_OPTIONS[
                self.parameterAsEnum(parameters, self.RESAMPLING_METHOD, context)
            ],
            "grid_spacing": self.parameterAsDouble(parameters, self.GRID_SPACING, context),
            "force_reproject": self.parameterAsBool(parameters, self.FORCE_REPROJECT, context),
            "force_resampling": force_resampling,
            "forced_grid_spacing": forced_grid_spacing,
            "allow_lambert_projection": self.parameterAsBool(parameters, self.ALLOW_LAMBERT, context),
            "warp_mem_mb": self.parameterAsInt(parameters, self.WARP_MEM_MB, context),
            # masking
            "min_elevation": self.parameterAsDouble(parameters, self.MINIMUM_ELEVATION, context),
            "maximum_elevation": self.parameterAsDouble(parameters, self.MAXIMUM_ELEVATION, context),
            # nodata repair
            "write_filled_nodata_raster": write_nodata_filled,
            "nodata_fill_method_choice": self.NODATA_FILL_OPTIONS[
                self.parameterAsEnum(parameters, self.NODATA_FILL_METHOD, context)
            ],
            "nodata_fill_max_iter": self.parameterAsInt(parameters, self.NODATA_FILL_MAX_ITER, context),
            "write_clipped_rectangle_raster": write_clipped_rect,
            # breaching
            "write_breached_raster": write_breached,
            "write_carve_diff_raster": write_carve_diff,
            "breach_before_fill": breach_before_fill,
            "fill_remainder_after_breach": self.parameterAsBool(
                parameters, self.FILL_REMAINDER_AFTER_BREACH, context
            ),
            "breach_algorithm_choice": self.BREACH_ALGORITHM_OPTIONS[
                self.parameterAsEnum(parameters, self.BREACH_ALGORITHM, context)
            ],
            "breach_lc_flat_increment": self.parameterAsDouble(
                parameters, self.BREACH_LC_FLAT_INCREMENT, context
            ),
            "breach_lc_minimize_dist": self.parameterAsBool(
                parameters, self.BREACH_LC_MINIMIZE_DIST, context
            ),
            "breach_lc_edge_policy": self.BREACH_LC_EDGE_POLICY_OPTIONS[
                self.parameterAsEnum(parameters, self.BREACH_LC_EDGE_POLICY, context)
            ],
            # filling
            "write_fill_raster": write_fill,
            "write_fill_diff_raster": write_fill_diff,
            "fill_algorithm_choice": self.FILL_ALGORITHM_OPTIONS[
                self.parameterAsEnum(parameters, self.FILL_ALGORITHM, context)
            ],
            "fill_min_slope": self.parameterAsDouble(parameters, self.FILL_MIN_SLOPE, context),
            # depression polygons
            "write_depression_polygons": write_depressions,
            "depression_min_pixels": self.parameterAsInt(parameters, self.DEPRESSION_MIN_PIXELS, context),
            "depression_min_depth": self.parameterAsDouble(parameters, self.DEPRESSION_MIN_DEPTH, context),
            # hillshade
            "write_hillshade_raster": write_hillshade,
            "hillshade_azimuth": self.parameterAsDouble(parameters, self.HILLSHADE_AZIMUTH, context),
            "hillshade_altitude": self.parameterAsDouble(parameters, self.HILLSHADE_ALTITUDE, context),
            "hillshade_z_factor": self.parameterAsDouble(parameters, self.HILLSHADE_Z_FACTOR, context),
            "write_multidirectional_hillshade_raster": write_md_hillshade,
            "multidirectional_hillshade_n_directions": self.parameterAsInt(
                parameters, self.MD_N_DIRECTIONS, context
            ),
            "multidirectional_hillshade_start_azimuth": self.parameterAsDouble(
                parameters, self.MD_START_AZIMUTH, context
            ),
            # flat resolution and fill-spill-merge
            "priority_flood_queue": self.PRIORITY_FLOOD_QUEUE_OPTIONS[
                self.parameterAsEnum(parameters, self.PRIORITY_FLOOD_QUEUE, context)
            ],
            "fill_flat_resolution_choice": self.FLAT_RESOLUTION_OPTIONS[
                self.parameterAsEnum(parameters, self.FLAT_RESOLUTION, context)
            ],
            "fill_spill_merge_rainfall_rate": self.parameterAsDouble(
                parameters, self.FSM_RAINFALL_RATE, context
            ),
            "fsm_boundary_conditions": fsm_boundary_conditions,
            # vertical units
            "vertical_units": self.VERTICAL_UNITS_OPTIONS[
                self.parameterAsEnum(parameters, self.VERTICAL_UNITS, context)
            ],
            "check_vertical_units_against_usgs": self.parameterAsBool(
                parameters, self.CHECK_VERTICAL_UNITS_USGS, context
            ),
            # smoothing
            "smooth_algorithm_choice": smooth_algorithm,
            "perona_malik_iterations": self.parameterAsInt(parameters, self.PM_ITERATIONS, context),
            "perona_malik_lambda": self.parameterAsDouble(parameters, self.PM_LAMBDA, context),
            "perona_malik_conductance": self.PM_CONDUCTANCE_OPTIONS[
                self.parameterAsEnum(parameters, self.PM_CONDUCTANCE, context)
            ],
            "perona_malik_k_quantile": self.parameterAsDouble(parameters, self.PM_K_QUANTILE, context),
            "perona_malik_k_metres_per_cell": self.parameterAsDouble(
                parameters, self.PM_K_METRES_PER_CELL, context
            ),
            "write_stencil_difference_raster": write_stencil_diff,
            "coherence_iterations": self.parameterAsInt(parameters, self.COHERENCE_ITERATIONS, context),
            "coherence_time_step": self.parameterAsDouble(parameters, self.COHERENCE_TIME_STEP, context),
            "coherence_derivative_radius_metres": self.parameterAsDouble(
                parameters, self.COHERENCE_DERIVATIVE_RADIUS, context
            ),
            "coherence_integration_radius_metres": self.parameterAsDouble(
                parameters, self.COHERENCE_INTEGRATION_RADIUS, context
            ),
            "coherence_contrast": self.parameterAsDouble(parameters, self.COHERENCE_CONTRAST, context),
            "coherence_min_diffusivity": self.parameterAsDouble(
                parameters, self.COHERENCE_MIN_DIFFUSIVITY, context
            ),
            # sink candidates (a point layer in write_vector_format)
            "write_sink_candidates": write_sink_candidates,
            "sink_candidate_min_pixels": self.parameterAsInt(
                parameters, self.SINK_CANDIDATE_MIN_PIXELS, context
            ),
        }

        if replace_layer is not None:
            param_values["replace_pixels_fname"] = layer_file_path(replace_layer.source())
            param_values["replace_pixels_value_field"] = self.parameterAsString(
                parameters, self.REPLACE_PIXELS_VALUE_FIELD, context
            )
            param_values["write_replaced_pixels_raster"] = write_replaced

        # The CLI default for breach_lc_max_cost is unlimited (float max);
        # only pass the key when the user sets a positive limit.
        breach_lc_max_cost = self.parameterAsDouble(parameters, self.BREACH_LC_MAX_COST, context)
        if breach_lc_max_cost > 0.0:
            param_values["breach_lc_max_cost"] = breach_lc_max_cost

        # 0 leaves breach_lc_max_dist out, so the CLI applies the edge
        # policy's default (1000 outlet, 100 whitebox).
        breach_lc_max_dist = self.parameterAsInt(parameters, self.BREACH_LC_MAX_DIST, context)
        if breach_lc_max_dist > 0:
            param_values["breach_lc_max_dist"] = breach_lc_max_dist

        if mask_layer is not None:
            param_values["run_binary_mask"] = True
            param_values["binary_mask_fname"] = mask_layer.source().split("|")[0]

        run_lsdtt_tool(
            "lsdtt-raster-preprocessing", param_values, output_folder, output_prefix, feedback
        )

        outputs = {}
        if write_fill:
            outputs[self.OUTPUT_FILL] = require_output(
                output_folder, output_prefix, "_fill.tif", "fill"
            )
        if write_breached or breach_before_fill:
            outputs[self.OUTPUT_BREACHED] = require_output(
                output_folder, output_prefix, "_breached.tif", "breached"
            )
        if write_hillshade:
            outputs[self.OUTPUT_HILLSHADE] = require_output(
                output_folder, output_prefix, "_hillshade.tif", "hillshade"
            )
        if write_md_hillshade:
            outputs[self.OUTPUT_MD_HILLSHADE] = require_output(
                output_folder, output_prefix, "_hillshade_multidirectional.tif",
                "multidirectional hillshade",
            )
        if write_nodata_filled:
            outputs[self.OUTPUT_NODATA_FILLED] = require_output(
                output_folder, output_prefix, "_nodata_filled.tif", "nodata-filled"
            )
        if write_clipped_rect:
            outputs[self.OUTPUT_CLIPPED_RECT] = require_output(
                output_folder, output_prefix, "_clipped_rect.tif", "clipped rectangle"
            )
        if write_fill_diff or write_depressions:
            outputs[self.OUTPUT_FILL_DIFF] = require_output(
                output_folder, output_prefix, "_fill_diff.tif", "fill difference"
            )
        if write_carve_diff:
            outputs[self.OUTPUT_CARVE_DIFF] = require_output(
                output_folder, output_prefix, "_carve_diff.tif", "carve difference"
            )
        if mask_layer is not None:
            outputs[self.OUTPUT_BINARY_MASKED] = require_output(
                output_folder, output_prefix, "_binary_masked.tif", "binary-masked"
            )
        if smooth_algorithm != "none":
            outputs[self.OUTPUT_SMOOTHED] = require_output(
                output_folder, output_prefix, "_smoothed.tif", "smoothed DEM"
            )
        if write_stencil_diff:
            outputs[self.OUTPUT_STENCIL_DIFF] = require_output(
                output_folder, output_prefix, "_stencil_difference.tif",
                "neighbour-difference",
            )
        if write_replaced:
            outputs[self.OUTPUT_REPLACED_PIXELS] = require_output(
                output_folder, output_prefix, "_pixels_replaced.tif", "replaced-pixel DEM"
            )
        if write_sink_candidates:
            outputs[self.OUTPUT_SINK_CANDIDATES] = require_output(
                output_folder, output_prefix,
                "_sink_candidates.{}".format(extension_for_format(depression_format)),
                "sink candidates",
            )
        if write_depressions:
            outputs[self.OUTPUT_DEPRESSIONS] = require_output(
                output_folder, output_prefix,
                "_depressions.{}".format(extension_for_format(depression_format)),
                "depression polygons",
            )
        return outputs

    def name(self):
        return "raster_preprocessing"

    def displayName(self):
        return "Raster Preprocessing"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Runs lsdtt-raster-preprocessing using a generated parameter file. "
            "Covers reprojection, elevation and binary masking, interior nodata "
            "filling, edge-preserving smoothing (Perona-Malik or coherence "
            "diffusion), pixel replacement, depression breaching and pit filling "
            "(with flat resolution, difference rasters, depression polygons and "
            "sink candidates), vertical-unit checks, and single or "
            "multidirectional hillshade. Backend path is managed in plugin settings."
        )

    def createInstance(self):
        return RasterPreprocessingAlgorithm()
