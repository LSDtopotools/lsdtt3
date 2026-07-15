import glob
import os
import tempfile

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputFile,
    QgsProcessingOutputVectorLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
)

from ..backend.manager import BackendManager
from ..backend.runner import run_cli, write_log
from ..templates.lsdtt_params import write_parameter_file


class ChiAnalysisAlgorithm(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    RASTER_IS_FILLED = "RASTER_IS_FILLED"
    FILL_MIN_SLOPE = "FILL_MIN_SLOPE"
    SOURCE_THRESHOLD = "SOURCE_THRESHOLD"
    MINIMUM_CHANNEL_PIXELS = "MINIMUM_CHANNEL_PIXELS"
    N_LARGEST_BASINS = "N_LARGEST_BASINS"
    CHI_A0 = "CHI_A0"
    MN_MIN = "MN_MIN"
    MN_MAX = "MN_MAX"
    MN_STEP = "MN_STEP"
    FIT_NETWORKS = "FIT_NETWORKS"
    WRITE_FIT_TABLES = "WRITE_FIT_TABLES"
    WRITE_KSN_SUMMARY = "WRITE_KSN_SUMMARY"
    WRITE_SEGMENT_LINES = "WRITE_SEGMENT_LINES"
    SEGMENT_LINES_FORMAT = "SEGMENT_LINES_FORMAT"

    OUTPUT_DISORDER_SWEEP = "OUTPUT_DISORDER_SWEEP"
    OUTPUT_SEGMENT_LINES = "OUTPUT_SEGMENT_LINES"

    VECTOR_FORMAT_OPTIONS = ["flatgeobuf", "gpkg"]

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
                self.SOURCE_THRESHOLD,
                "Source threshold (pixels)",
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
                "Keep N largest basins (0 = all)",
                type=QgsProcessingParameterNumber.Integer,
                minValue=0,
                defaultValue=5,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.CHI_A0,
                "Chi A0",
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

        self.addOutput(QgsProcessingOutputFile(self.OUTPUT_DISORDER_SWEEP, "Chi disorder sweep CSV"))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT_SEGMENT_LINES, "Chi fit segment lines"))

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        if dem_layer is None:
            raise QgsProcessingException("Input DEM is required.")

        dem_path = dem_layer.source()
        read_path, read_fname = os.path.split(dem_path)

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        output_prefix = self.parameterAsString(parameters, self.OUTPUT_PREFIX, context)
        os.makedirs(output_folder, exist_ok=True)
        write_segment_lines = self.parameterAsBool(parameters, self.WRITE_SEGMENT_LINES, context)
        segment_format = self.VECTOR_FORMAT_OPTIONS[
            self.parameterAsEnum(parameters, self.SEGMENT_LINES_FORMAT, context)
        ]

        param_values = {
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "raster_is_filled": self.parameterAsBool(parameters, self.RASTER_IS_FILLED, context),
            "fill_min_slope": self.parameterAsDouble(parameters, self.FILL_MIN_SLOPE, context),
            "source_pixels_threshold": self.parameterAsInt(parameters, self.SOURCE_THRESHOLD, context),
            "minimum_channel_pixels": self.parameterAsInt(parameters, self.MINIMUM_CHANNEL_PIXELS, context),
            "n_largest_basins": self.parameterAsInt(parameters, self.N_LARGEST_BASINS, context),
            "chi_A0": self.parameterAsDouble(parameters, self.CHI_A0, context),
            "mn_min": self.parameterAsDouble(parameters, self.MN_MIN, context),
            "mn_max": self.parameterAsDouble(parameters, self.MN_MAX, context),
            "mn_step": self.parameterAsDouble(parameters, self.MN_STEP, context),
            "fit_chi_networks": self.parameterAsBool(parameters, self.FIT_NETWORKS, context),
            "write_fit_tables": self.parameterAsBool(parameters, self.WRITE_FIT_TABLES, context),
            "write_ksn_basin_summary": self.parameterAsBool(parameters, self.WRITE_KSN_SUMMARY, context),
            "write_chi_fit_segment_lines": write_segment_lines,
            "chi_fit_segment_lines_output_format": segment_format,
            "run_inversion": False,
        }

        with tempfile.TemporaryDirectory(prefix="lsdtt_qgis_") as temp_dir:
            params_name = "chi_analysis.driver"
            params_file = os.path.join(temp_dir, params_name)
            write_parameter_file(params_file, param_values)

            backend = BackendManager()
            executable = backend.executable_path("lsdtt-chi-analysis")
            command = [executable, temp_dir, params_name]
            log_path = os.path.join(output_folder, "{}_chi_analysis.log".format(output_prefix))
            feedback.pushInfo("Running lsdtt-chi-analysis...")
            result = run_cli(command, feedback=feedback, env=backend.runtime_environment())
            write_log(log_path, command, result)
            if result.returncode != 0:
                raise QgsProcessingException(
                    "lsdtt-chi-analysis failed with exit code {}. See log: {}".format(
                        result.returncode,
                        log_path,
                    )
                )

        disorder_path = os.path.join(output_folder, "{}_chi_disorder_sweep.csv".format(output_prefix))
        if not os.path.exists(disorder_path):
            raise QgsProcessingException("Expected chi disorder sweep output not found: {}".format(disorder_path))

        outputs = {self.OUTPUT_DISORDER_SWEEP: disorder_path}
        if write_segment_lines:
            segment_pattern = os.path.join(
                output_folder,
                "{}_chi_fit_segment_lines_mn_*.{}".format(output_prefix, _extension_for_format(segment_format)),
            )
            segment_matches = sorted(glob.glob(segment_pattern))
            if segment_matches:
                outputs[self.OUTPUT_SEGMENT_LINES] = segment_matches[0]
        return outputs

    def name(self):
        return "chi_analysis"

    def displayName(self):
        return "Chi Analysis"

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def createInstance(self):
        return ChiAnalysisAlgorithm()


def _extension_for_format(format_name):
    if format_name == "flatgeobuf":
        return "fgb"
    return format_name
