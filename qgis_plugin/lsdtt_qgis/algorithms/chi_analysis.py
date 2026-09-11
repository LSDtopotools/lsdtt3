import glob
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


class ChiAnalysisAlgorithm(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    RASTER_IS_FILLED = "RASTER_IS_FILLED"
    FILL_MIN_SLOPE = "FILL_MIN_SLOPE"
    THRESHOLD_CONTRIBUTING_PIXELS = "THRESHOLD_CONTRIBUTING_PIXELS"
    MINIMUM_CHANNEL_PIXELS = "MINIMUM_CHANNEL_PIXELS"
    N_LARGEST_BASINS = "N_LARGEST_BASINS"
    BASIN_OUTLET_FILE = "BASIN_OUTLET_FILE"
    CHI_A0 = "CHI_A0"
    MN_MIN = "MN_MIN"
    MN_MAX = "MN_MAX"
    MN_STEP = "MN_STEP"
    FIT_NETWORKS = "FIT_NETWORKS"
    WRITE_FIT_TABLES = "WRITE_FIT_TABLES"
    WRITE_KSN_SUMMARY = "WRITE_KSN_SUMMARY"
    WRITE_SEGMENT_LINES = "WRITE_SEGMENT_LINES"
    SEGMENT_LINES_FORMAT = "SEGMENT_LINES_FORMAT"
    PRINT_DATA_MAPS = "PRINT_DATA_MAPS"
    DATA_MAPS_FORMAT = "DATA_MAPS_FORMAT"
    WRITE_DATA_MAP_LINES = "WRITE_DATA_MAP_LINES"
    DATA_MAP_LINES_FORMAT = "DATA_MAP_LINES_FORMAT"

    # Advanced: run control
    THREADS = "THREADS"
    BOUNDARY_CONDITIONS = "BOUNDARY_CONDITIONS"
    OUTPUT_FORMAT = "OUTPUT_FORMAT"
    CALCULATE_SPATIAL_INDEX = "CALCULATE_SPATIAL_INDEX"

    # Advanced: outlet snapping
    SNAPPING_METHOD = "SNAPPING_METHOD"
    SNAPPING_THRESHOLD = "SNAPPING_THRESHOLD"
    SNAPPING_MIN_RADIUS = "SNAPPING_MIN_RADIUS"
    SNAPPING_MAX_RADIUS = "SNAPPING_MAX_RADIUS"
    MOVE_OUTLET_UPSTREAM = "MOVE_OUTLET_UPSTREAM"

    # Advanced: fitting detail
    M_OVER_N = "M_OVER_N"
    FIT_CRITICAL_DIVERGENCE = "FIT_CRITICAL_DIVERGENCE"

    # Advanced: disorder uncertainty
    RUN_SUBSET_SENSITIVITY = "RUN_SUBSET_SENSITIVITY"
    SUBSET_N_SIDE_BRANCHES = "SUBSET_N_SIDE_BRANCHES"
    RUN_JACKKNIFE = "RUN_JACKKNIFE"
    RUN_BOOTSTRAP = "RUN_BOOTSTRAP"
    BOOTSTRAP_REPLICATES = "BOOTSTRAP_REPLICATES"
    BOOTSTRAP_FULL_TRIBUTARY_COUNT = "BOOTSTRAP_FULL_TRIBUTARY_COUNT"
    BOOTSTRAP_N_SIDE_BRANCHES = "BOOTSTRAP_N_SIDE_BRANCHES"
    BOOTSTRAP_WITH_REPLACEMENT = "BOOTSTRAP_WITH_REPLACEMENT"
    BOOTSTRAP_SEED = "BOOTSTRAP_SEED"
    WRITE_DISORDER_SAMPLE_TABLES = "WRITE_DISORDER_SAMPLE_TABLES"

    # Advanced: inversion
    RUN_INVERSION = "RUN_INVERSION"
    INVERSION_COARSE_CELL_SIZE = "INVERSION_COARSE_CELL_SIZE"
    INVERSION_ROUGHNESS_WEIGHT = "INVERSION_ROUGHNESS_WEIGHT"
    INVERSION_DIAGONAL_JITTER = "INVERSION_DIAGONAL_JITTER"
    INVERSION_USE_MAE = "INVERSION_USE_MAE"
    WRITE_INVERSION_RASTERS = "WRITE_INVERSION_RASTERS"

    OUTPUT_DISORDER_SWEEP = "OUTPUT_DISORDER_SWEEP"
    OUTPUT_SEGMENT_LINES = "OUTPUT_SEGMENT_LINES"
    OUTPUT_DATA_MAPS = "OUTPUT_DATA_MAPS"
    OUTPUT_DATA_MAP_LINES = "OUTPUT_DATA_MAP_LINES"
    OUTPUT_INVERSION_BEST = "OUTPUT_INVERSION_BEST"
    OUTPUT_INVERSION_TRIALS = "OUTPUT_INVERSION_TRIALS"

    VECTOR_FORMAT_OPTIONS = ["flatgeobuf", "gpkg"]
    POINT_FORMAT_OPTIONS = ["csv", "flatgeobuf", "gpkg"]
    OUTPUT_FORMAT_OPTIONS = ["cog", "gtiff"]
    SNAPPING_OPTIONS = ["downstream", "expanding_square", "none"]

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
        self.addParameter(QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", "lsdtt_chi"))
        self.addParameter(QgsProcessingParameterBoolean(self.RASTER_IS_FILLED, "Raster already filled", False))
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
            QgsProcessingParameterNumber(
                self.MINIMUM_CHANNEL_PIXELS,
                "Minimum channel pixels",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=20,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.N_LARGEST_BASINS,
                "Keep N largest basins (0 = all complete basins)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=5,
            )
        )
        self.addParameter(
            QgsProcessingParameterFile(
                self.BASIN_OUTLET_FILE,
                "Basin outlet file (optional; CSV needs latitude/longitude columns)",
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.CHI_A0,
                "Chi A0 (reference drainage area)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MN_MIN,
                "m/n min",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.3,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MN_MAX,
                "m/n max",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.6,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MN_STEP,
                "m/n step",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.02,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.FIT_NETWORKS, "Fit chi networks", True))
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_FIT_TABLES, "Write fit tables", True))
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_KSN_SUMMARY, "Write ksn summary", True))
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_SEGMENT_LINES, "Write segment lines", False))
        self.addParameter(
            QgsProcessingParameterEnum(
                self.SEGMENT_LINES_FORMAT,
                "Segment lines format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.PRINT_DATA_MAPS,
                "Write chi data maps (per-node points, outlets, sources, basin polygons)",
                False,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.DATA_MAPS_FORMAT,
                "Data maps format",
                options=self.POINT_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_DATA_MAP_LINES,
                "Write chi channels as polylines (one per source)",
                False,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.DATA_MAP_LINES_FORMAT,
                "Data map lines format",
                options=self.VECTOR_FORMAT_OPTIONS,
                defaultValue=0,
            )
        )

        # --- Advanced: run control ---
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
            QgsProcessingParameterString(
                self.BOUNDARY_CONDITIONS,
                "Boundary conditions (four characters NESW; o=open, n=no-flux, p=periodic)",
                "oooo",
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
            QgsProcessingParameterBoolean(self.CALCULATE_SPATIAL_INDEX, "Calculate spatial index", False)
        )

        # --- Advanced: outlet snapping ---
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
                "Outlet snapping contributing pixels threshold",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1000,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SNAPPING_MIN_RADIUS,
                "Outlet snapping minimum search radius (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SNAPPING_MAX_RADIUS,
                "Outlet snapping maximum search radius (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=25,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MOVE_OUTLET_UPSTREAM,
                "Move auto-detected outlets upstream (pixels)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=10,
            )
        )

        # --- Advanced: fitting detail ---
        self._add_advanced(
            QgsProcessingParameterString(
                self.M_OVER_N,
                "Explicit m/n value(s), comma-separated (overrides the sweep)",
                "",
                optional=True,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FIT_CRITICAL_DIVERGENCE,
                "Critical elevation divergence for branch splitting (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=5.0,
            )
        )

        # --- Advanced: disorder uncertainty ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.RUN_SUBSET_SENSITIVITY,
                "Run subset sensitivity analysis on the disorder metric",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.SUBSET_N_SIDE_BRANCHES,
                "Subset sensitivity: side branches per sample",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=3,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.RUN_JACKKNIFE,
                "Run jackknife uncertainty on the disorder metric",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.RUN_BOOTSTRAP,
                "Run bootstrap uncertainty on the disorder metric",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.BOOTSTRAP_REPLICATES,
                "Bootstrap: number of replicates",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1000,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.BOOTSTRAP_FULL_TRIBUTARY_COUNT,
                "Bootstrap: sample the full tributary count",
                True,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.BOOTSTRAP_N_SIDE_BRANCHES,
                "Bootstrap: side branches per replicate",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=3,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.BOOTSTRAP_WITH_REPLACEMENT,
                "Bootstrap: sample with replacement",
                True,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.BOOTSTRAP_SEED,
                "Bootstrap: random seed",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=5489,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_DISORDER_SAMPLE_TABLES,
                "Write per-sample disorder tables",
                False,
            )
        )

        # --- Advanced: inversion ---
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.RUN_INVERSION,
                "Run k_sn/u* inversion (requires an Eigen-enabled backend)",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.INVERSION_COARSE_CELL_SIZE,
                "Inversion: coarse cell size (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=2000.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.INVERSION_ROUGHNESS_WEIGHT,
                "Inversion: roughness weight",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0e4,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.INVERSION_DIAGONAL_JITTER,
                "Inversion: diagonal jitter",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0e-8,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.INVERSION_USE_MAE,
                "Inversion: use MAE for model selection",
                True,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.WRITE_INVERSION_RASTERS,
                "Inversion: write per-basin u* and support-count rasters",
                False,
            )
        )

        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_DISORDER_SWEEP, "Chi disorder sweep CSV"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_SEGMENT_LINES, "Chi fit segment lines"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_DATA_MAPS, "Chi data maps (per-node points)"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_DATA_MAP_LINES, "Chi data map lines"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_INVERSION_BEST, "Inversion best-model CSV"))
        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_INVERSION_TRIALS, "Inversion trials CSV"))

        attach_backend_help(
            self,
            "lsdtt-chi-analysis",
            overrides={
                self.BOOTSTRAP_SEED: "bootstrap_random_seed",
                self.FIT_NETWORKS: "fit_chi_networks",
                self.PRINT_DATA_MAPS: "print_chi_data_maps",
                self.THRESHOLD_CONTRIBUTING_PIXELS: "threshold_contributing_pixels",
                self.WRITE_DATA_MAP_LINES: "write_chi_data_map_lines",
                self.WRITE_KSN_SUMMARY: "write_ksn_basin_summary",
                self.WRITE_SEGMENT_LINES: "write_chi_fit_segment_lines",
                self.BOOTSTRAP_FULL_TRIBUTARY_COUNT: "bootstrap_sample_full_tributary_count",
                self.BOOTSTRAP_N_SIDE_BRANCHES: "bootstrap_n_side_branches_per_replicate",
                self.DATA_MAPS_FORMAT: "chi_data_maps_output_format",
                self.DATA_MAP_LINES_FORMAT: "chi_data_map_lines_output_format",
                self.FIT_CRITICAL_DIVERGENCE: "fit_critical_divergence_m",
                self.INVERSION_USE_MAE: "inversion_use_mae_for_model_selection",
                self.MOVE_OUTLET_UPSTREAM: "n_pixels_move_outlet_upstream",
                self.SEGMENT_LINES_FORMAT: "chi_fit_segment_lines_output_format",
                self.SNAPPING_MAX_RADIUS: "basin_outlet_max_search_radius_pixels",
                self.SNAPPING_METHOD: "basin_outlet_snapping_method",
                self.SNAPPING_MIN_RADIUS: "basin_outlet_min_search_radius_pixels",
                self.SNAPPING_THRESHOLD: "basin_outlet_contributing_pixels_threshold",
                self.SUBSET_N_SIDE_BRANCHES: "subset_n_side_branches_per_sample",
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

        boundary_conditions = self.parameterAsString(parameters, self.BOUNDARY_CONDITIONS, context)
        if len(boundary_conditions) != 4 or any(c not in "onpb" for c in boundary_conditions.lower()):
            raise QgsProcessingException(
                "Boundary conditions must be four characters from o, n, p, b (e.g. oooo)."
            )

        fit_networks = self.parameterAsBool(parameters, self.FIT_NETWORKS, context)
        write_segment_lines = self.parameterAsBool(parameters, self.WRITE_SEGMENT_LINES, context)
        if write_segment_lines and not fit_networks:
            fit_networks = True
            feedback.pushInfo(
                "Chi fit segment lines require network fitting; enabling fit_chi_networks."
            )
        print_data_maps = self.parameterAsBool(parameters, self.PRINT_DATA_MAPS, context)
        write_data_map_lines = self.parameterAsBool(parameters, self.WRITE_DATA_MAP_LINES, context)
        run_inversion = self.parameterAsBool(parameters, self.RUN_INVERSION, context)

        if not any([fit_networks, print_data_maps, run_inversion]):
            raise QgsProcessingException(
                "Enable at least one analysis (fit chi networks, chi data maps, or inversion)."
            )

        if run_inversion:
            from ..backend.capabilities import require_feature
            error = require_feature(
                "lsdtt-chi-analysis", "eigen",
                "The chi k_sn/u* inversion", feedback,
            )
            if error:
                raise QgsProcessingException(error)

        segment_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.SEGMENT_LINES_FORMAT, context)
        ]
        data_maps_format = self.POINT_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.DATA_MAPS_FORMAT, context)
        ]
        data_map_lines_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.DATA_MAP_LINES_FORMAT, context)
        ]

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "n_threads": self.parameterAsInt(parameters, self.THREADS, context),
            "boundary_conditions": boundary_conditions,
            "write_output_format": self.OUTPUT_FORMAT_OPTIONS[
                self.parameterAsEnum(parameters, self.OUTPUT_FORMAT, context)
            ],
            "calculate_spatial_index": self.parameterAsBool(
                parameters, self.CALCULATE_SPATIAL_INDEX, context
            ),
            # DEM preparation and network extraction
            "raster_is_filled": self.parameterAsBool(parameters, self.RASTER_IS_FILLED, context),
            "fill_min_slope": self.parameterAsDouble(parameters, self.FILL_MIN_SLOPE, context),
            "threshold_contributing_pixels": self.parameterAsInt(
                parameters, self.THRESHOLD_CONTRIBUTING_PIXELS, context
            ),
            "minimum_channel_pixels": self.parameterAsInt(parameters, self.MINIMUM_CHANNEL_PIXELS, context),
            "n_largest_basins": self.parameterAsInt(parameters, self.N_LARGEST_BASINS, context),
            "n_pixels_move_outlet_upstream": self.parameterAsInt(
                parameters, self.MOVE_OUTLET_UPSTREAM, context
            ),
            # chi extraction / fitting
            "chi_A0": self.parameterAsDouble(parameters, self.CHI_A0, context),
            "mn_min": self.parameterAsDouble(parameters, self.MN_MIN, context),
            "mn_max": self.parameterAsDouble(parameters, self.MN_MAX, context),
            "mn_step": self.parameterAsDouble(parameters, self.MN_STEP, context),
            "fit_chi_networks": fit_networks,
            "fit_critical_divergence_m": self.parameterAsDouble(
                parameters, self.FIT_CRITICAL_DIVERGENCE, context
            ),
            "write_fit_tables": self.parameterAsBool(parameters, self.WRITE_FIT_TABLES, context),
            "write_ksn_basin_summary": self.parameterAsBool(parameters, self.WRITE_KSN_SUMMARY, context),
            "write_chi_fit_segment_lines": write_segment_lines,
            "chi_fit_segment_lines_output_format": segment_format,
            # diagnostic data maps
            "print_chi_data_maps": print_data_maps,
            "chi_data_maps_output_format": data_maps_format,
            "write_chi_data_map_lines": write_data_map_lines,
            "chi_data_map_lines_output_format": data_map_lines_format,
            # disorder uncertainty
            "run_subset_sensitivity": self.parameterAsBool(
                parameters, self.RUN_SUBSET_SENSITIVITY, context
            ),
            "subset_n_side_branches_per_sample": self.parameterAsInt(
                parameters, self.SUBSET_N_SIDE_BRANCHES, context
            ),
            "run_jackknife": self.parameterAsBool(parameters, self.RUN_JACKKNIFE, context),
            "run_bootstrap": self.parameterAsBool(parameters, self.RUN_BOOTSTRAP, context),
            "bootstrap_replicates": self.parameterAsInt(parameters, self.BOOTSTRAP_REPLICATES, context),
            "bootstrap_sample_full_tributary_count": self.parameterAsBool(
                parameters, self.BOOTSTRAP_FULL_TRIBUTARY_COUNT, context
            ),
            "bootstrap_n_side_branches_per_replicate": self.parameterAsInt(
                parameters, self.BOOTSTRAP_N_SIDE_BRANCHES, context
            ),
            "bootstrap_with_replacement": self.parameterAsBool(
                parameters, self.BOOTSTRAP_WITH_REPLACEMENT, context
            ),
            "bootstrap_random_seed": self.parameterAsInt(parameters, self.BOOTSTRAP_SEED, context),
            "write_disorder_sample_tables": self.parameterAsBool(
                parameters, self.WRITE_DISORDER_SAMPLE_TABLES, context
            ),
            # inversion
            "run_inversion": run_inversion,
            "inversion_coarse_cell_size": self.parameterAsDouble(
                parameters, self.INVERSION_COARSE_CELL_SIZE, context
            ),
            "inversion_roughness_weight": self.parameterAsDouble(
                parameters, self.INVERSION_ROUGHNESS_WEIGHT, context
            ),
            "inversion_diagonal_jitter": self.parameterAsDouble(
                parameters, self.INVERSION_DIAGONAL_JITTER, context
            ),
            "inversion_use_mae_for_model_selection": self.parameterAsBool(
                parameters, self.INVERSION_USE_MAE, context
            ),
            "write_inversion_rasters": self.parameterAsBool(
                parameters, self.WRITE_INVERSION_RASTERS, context
            ),
        }

        outlet_file = self.parameterAsFile(parameters, self.BASIN_OUTLET_FILE, context)
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

        m_over_n = self.parameterAsString(parameters, self.M_OVER_N, context)
        if m_over_n:
            param_values["m_over_n"] = m_over_n

        run_lsdtt_tool(
            "lsdtt-chi-analysis", param_values, output_folder, output_prefix, feedback
        )

        outputs = {}
        if fit_networks:
            outputs[self.OUTPUT_DISORDER_SWEEP] = require_output(
                output_folder, output_prefix, "_chi_disorder_sweep.csv",
                "chi disorder sweep",
            )
        if write_segment_lines:
            segment_pattern = os.path.join(
                output_folder,
                "{}_chi_fit_segment_lines_mn_*.{}".format(
                    output_prefix, extension_for_format(segment_format)
                ),
            )
            segment_matches = sorted(glob.glob(segment_pattern))
            if not segment_matches:
                raise QgsProcessingException(
                    "Expected chi fit segment lines output not found: {}. "
                    "Check the run log in the output folder for warnings.".format(segment_pattern)
                )
            outputs[self.OUTPUT_SEGMENT_LINES] = segment_matches[0]
        if print_data_maps:
            outputs[self.OUTPUT_DATA_MAPS] = require_output(
                output_folder, output_prefix,
                "_chi_data_maps.{}".format(extension_for_format(data_maps_format)),
                "chi data maps",
            )
        if write_data_map_lines and print_data_maps:
            outputs[self.OUTPUT_DATA_MAP_LINES] = require_output(
                output_folder, output_prefix,
                "_chi_data_map_lines.{}".format(extension_for_format(data_map_lines_format)),
                "chi data map lines",
            )
        if run_inversion:
            outputs[self.OUTPUT_INVERSION_BEST] = require_output(
                output_folder, output_prefix, "_chi_ksn_inversion_best.csv",
                "inversion best-model summary",
            )
            outputs[self.OUTPUT_INVERSION_TRIALS] = require_output(
                output_folder, output_prefix, "_chi_ksn_inversion_trials.csv",
                "inversion trials",
            )
        return outputs

    def name(self):
        return "chi_analysis"

    def displayName(self):
        return "Chi Analysis"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def shortHelpString(self):
        return (
            "Runs lsdtt-chi-analysis: chi extraction on threshold-defined "
            "channel networks (largest basins or supplied outlets with "
            "snapping), piecewise-linear chi-network fitting with an m/n "
            "disorder sweep (optional subset/jackknife/bootstrap "
            "uncertainty), diagnostic chi data maps, and optional spatial "
            "k_sn/u* inversion (Eigen-enabled backend required). Backend "
            "path is managed in plugin settings."
        )

    def createInstance(self):
        return ChiAnalysisAlgorithm()
