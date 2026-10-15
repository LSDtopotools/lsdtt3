"""Processing algorithms for lsdtt-drainage-analysis.

The CLI computes several independent products from one routed DEM. Each is
its own algorithm here, so the Processing toolbox shows a short, focused form
per product:

- Drainage Area (D8, D-infinity and/or MFD)
- Topographic Wetness Index
- Sediment Connectivity Index (Cavalli et al. 2013)
- Height Above Channel
- Soil Topographic Index (Scottish soil and land-use data)

They share the DEM, routing and output options through
DrainageAnalysisBase. Every write_* switch is written to the parameter file,
False unless the algorithm turns it on, so a run only produces what its form
asked for.
"""

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
    QgsProcessingParameterFile,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
    QgsProcessingParameterVectorLayer,
)

from .support import (
    attach_backend_help,
    layer_file_path,
    require_output,
    run_lsdtt_tool,
    source_layer_name,
    split_layer_source,
)


TOOL = "lsdtt-drainage-analysis"

# Every output switch the CLI has. Each run writes all of them, so the CLI's
# own defaults never add an output the form did not ask for.
WRITE_SWITCHES = (
    "write_d8_drainage_area_raster",
    "write_dinf_drainage_area_raster",
    "write_mfd_drainage_area_raster",
    "write_twi_raster",
    "write_connectivity_raster",
    "write_height_above_channel_raster",
    "write_sti_raster",
    "write_sti_proxy_raster",
    "write_sti_intermediate_rasters",
    "write_sti_lookup_csv",
    "write_sti_95th_percentile_rasters",
    "write_sti_p95_jenks_rasters",
)


class DrainageAnalysisBase(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    OUTPUT_PREFIX = "OUTPUT_PREFIX"
    RASTER_IS_FILLED = "RASTER_IS_FILLED"

    # Advanced
    FILL_MIN_SLOPE = "FILL_MIN_SLOPE"
    FILL_FLAT_RESOLUTION = "FILL_FLAT_RESOLUTION"
    BOUNDARY_CONDITIONS = "BOUNDARY_CONDITIONS"
    THREADS = "THREADS"
    OUTPUT_FORMAT = "OUTPUT_FORMAT"

    OUTPUT_FORMAT_OPTIONS = ["cog", "gtiff"]
    FILL_FLAT_RESOLUTION_OPTIONS = [
        "inline_gradient", "barnes2014b_gradient", "barnes2014b_float_step",
    ]

    # Set by each subclass.
    DEFAULT_PREFIX = "lsdtt_drainage"
    # (output id, CLI suffix, description); see expected_outputs().
    OUTPUTS = ()
    HELP_OVERRIDES = {}

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
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, "Output folder"))
        self.addParameter(
            QgsProcessingParameterString(self.OUTPUT_PREFIX, "Output prefix", self.DEFAULT_PREFIX)
        )
        self.add_product_parameters()
        self.addParameter(
            QgsProcessingParameterBoolean(self.RASTER_IS_FILLED, "Raster already filled", False)
        )

        self.add_advanced_product_parameters()
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.FILL_MIN_SLOPE,
                "Fill min slope",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0001,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.FILL_FLAT_RESOLUTION,
                "How flats are made drainable",
                options=self.FILL_FLAT_RESOLUTION_OPTIONS,
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

        for output_id, suffix, description in self.OUTPUTS:
            if suffix.endswith(".csv"):
                self.addOutput(QgsProcessingOutputFile(output_id, description))
            else:
                self.addOutput(QgsProcessingOutputRasterLayer(output_id, description))

        overrides = {
            self.FILL_FLAT_RESOLUTION: "fill_flat_resolution_choice",
        }
        overrides.update(self.HELP_OVERRIDES)
        attach_backend_help(self, TOOL, overrides=overrides)

    # --- hooks for subclasses -------------------------------------------------

    def add_product_parameters(self):
        """Add the main parameters of this product."""

    def add_advanced_product_parameters(self):
        """Add this product's advanced parameters."""

    def product_values(self, parameters, context):
        """Return this product's CLI keys, including its write_* switches."""
        return {}

    def expected_outputs(self, values):
        """Return the (output id, suffix, description) the run must produce."""
        return list(self.OUTPUTS)

    # --- shared run -----------------------------------------------------------

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

        param_values = {switch: False for switch in WRITE_SWITCHES}
        param_values.update({
            "read_path": read_path,
            "read_fname": read_fname,
            "write_path": output_folder,
            "write_fname": output_prefix,
            "n_threads": self.parameterAsInt(parameters, self.THREADS, context),
            "write_output_format": self.OUTPUT_FORMAT_OPTIONS[
                self.parameterAsEnum(parameters, self.OUTPUT_FORMAT, context)
            ],
            "boundary_conditions": boundary_conditions,
            "raster_is_filled": self.parameterAsBool(parameters, self.RASTER_IS_FILLED, context),
            "fill_min_slope": self.parameterAsDouble(parameters, self.FILL_MIN_SLOPE, context),
            "fill_flat_resolution_choice": self.FILL_FLAT_RESOLUTION_OPTIONS[
                self.parameterAsEnum(parameters, self.FILL_FLAT_RESOLUTION, context)
            ],
        })
        param_values.update(self.product_values(parameters, context))

        expected = self.expected_outputs(param_values)
        if not expected:
            raise QgsProcessingException("Select at least one output.")

        run_lsdtt_tool(TOOL, param_values, output_folder, output_prefix, feedback)

        return {
            output_id: require_output(output_folder, output_prefix, suffix, description)
            for output_id, suffix, description in expected
        }

    def group(self):
        return "lsdtt3"

    def groupId(self):
        return "lsdtt"

    def createInstance(self):
        return type(self)()


class DrainageAreaAlgorithm(DrainageAnalysisBase):
    WRITE_D8 = "WRITE_D8"
    WRITE_DINF = "WRITE_DINF"
    WRITE_MFD = "WRITE_MFD"
    MFD_EXPONENT = "MFD_EXPONENT"

    OUTPUT_D8 = "OUTPUT_D8"
    OUTPUT_DINF = "OUTPUT_DINF"
    OUTPUT_MFD = "OUTPUT_MFD"

    DEFAULT_PREFIX = "lsdtt_drainage_area"
    OUTPUTS = (
        (OUTPUT_D8, "_d8_drainage_area.tif", "D8 drainage area"),
        (OUTPUT_DINF, "_dinf_drainage_area.tif", "D-infinity drainage area"),
        (OUTPUT_MFD, "_mfd_drainage_area.tif", "MFD drainage area"),
    )
    _SWITCHES = (
        (WRITE_D8, "write_d8_drainage_area_raster", OUTPUT_D8),
        (WRITE_DINF, "write_dinf_drainage_area_raster", OUTPUT_DINF),
        (WRITE_MFD, "write_mfd_drainage_area_raster", OUTPUT_MFD),
    )
    HELP_OVERRIDES = {
        WRITE_D8: "write_d8_drainage_area_raster",
        WRITE_DINF: "write_dinf_drainage_area_raster",
        WRITE_MFD: "write_mfd_drainage_area_raster",
    }

    def add_product_parameters(self):
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_D8, "D8 drainage area", True))
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_DINF, "D-infinity drainage area", False))
        self.addParameter(
            QgsProcessingParameterBoolean(self.WRITE_MFD, "Multiple-flow-direction (MFD) drainage area", False)
        )

    def add_advanced_product_parameters(self):
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MFD_EXPONENT,
                "MFD slope exponent",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )

    def product_values(self, parameters, context):
        values = {
            key: self.parameterAsBool(parameters, param, context)
            for param, key, _ in self._SWITCHES
        }
        values["mfd_exponent"] = self.parameterAsDouble(parameters, self.MFD_EXPONENT, context)
        return values

    def expected_outputs(self, values):
        enabled = {output for _, key, output in self._SWITCHES if values[key]}
        return [entry for entry in self.OUTPUTS if entry[0] in enabled]

    def name(self):
        return "drainage_area"

    def displayName(self):
        return "Drainage Area"

    def shortHelpString(self):
        return (
            "Runs lsdtt-drainage-analysis: drainage area in square metres from "
            "D8 (all flow to the steepest neighbour), D-infinity (flow split "
            "between two neighbours) and/or multiple flow direction (flow "
            "split between all lower neighbours, weighted by slope). The DEM "
            "is pit-filled first unless it is marked as already filled. "
            "Backend path is managed in plugin settings."
        )


class TopographicWetnessIndexAlgorithm(DrainageAnalysisBase):
    TWI_DRAINAGE_AREA_METHOD = "TWI_DRAINAGE_AREA_METHOD"
    TWI_WINDOW_RADIUS_METRES = "TWI_WINDOW_RADIUS_METRES"
    MFD_EXPONENT = "MFD_EXPONENT"

    OUTPUT_TWI = "OUTPUT_TWI"

    METHOD_OPTIONS = ["Dinf", "MFD"]
    DEFAULT_PREFIX = "lsdtt_twi"
    OUTPUTS = ((OUTPUT_TWI, "_twi.tif", "Topographic wetness index"),)

    def add_product_parameters(self):
        self.addParameter(
            QgsProcessingParameterEnum(
                self.TWI_DRAINAGE_AREA_METHOD,
                "Drainage area method",
                options=self.METHOD_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.TWI_WINDOW_RADIUS_METRES,
                "Slope window radius (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=100.0,
            )
        )

    def add_advanced_product_parameters(self):
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MFD_EXPONENT,
                "MFD slope exponent (MFD method only)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )

    def product_values(self, parameters, context):
        return {
            "write_twi_raster": True,
            "twi_drainage_area_method": self.METHOD_OPTIONS[
                self.parameterAsEnum(parameters, self.TWI_DRAINAGE_AREA_METHOD, context)
            ],
            "twi_window_radius_metres": self.parameterAsDouble(
                parameters, self.TWI_WINDOW_RADIUS_METRES, context
            ),
            "mfd_exponent": self.parameterAsDouble(parameters, self.MFD_EXPONENT, context),
        }

    def name(self):
        return "topographic_wetness_index"

    def displayName(self):
        return "Topographic Wetness Index"

    def shortHelpString(self):
        return (
            "Runs lsdtt-drainage-analysis: the topographic wetness index "
            "ln(a / tan b), from specific catchment area a (D-infinity or MFD "
            "drainage area per unit contour width) and the local slope b, "
            "fitted over a moving window. High values mark places where water "
            "collects. Backend path is managed in plugin settings."
        )


class ConnectivityIndexAlgorithm(DrainageAnalysisBase):
    CONNECTIVITY_DRAINAGE_AREA_METHOD = "CONNECTIVITY_DRAINAGE_AREA_METHOD"
    CONNECTIVITY_AREA_THRESHOLD_M2 = "CONNECTIVITY_AREA_THRESHOLD_M2"
    CONNECTIVITY_WINDOW_RADIUS_METRES = "CONNECTIVITY_WINDOW_RADIUS_METRES"
    CONNECTIVITY_SLOPE_METHOD = "CONNECTIVITY_SLOPE_METHOD"
    CONNECTIVITY_ROUGHNESS_METHOD = "CONNECTIVITY_ROUGHNESS_METHOD"
    MFD_EXPONENT = "MFD_EXPONENT"

    OUTPUT_CONNECTIVITY = "OUTPUT_CONNECTIVITY"
    OUTPUT_D_UP = "OUTPUT_D_UP"
    OUTPUT_D_DOWN = "OUTPUT_D_DOWN"

    METHOD_OPTIONS = ["d8", "dinf", "mfd"]
    SLOPE_OPTIONS = ["polyfit_direct", "wavelet"]
    ROUGHNESS_OPTIONS = ["madk2", "iqr_residual"]
    DEFAULT_PREFIX = "lsdtt_connectivity"
    OUTPUTS = (
        (OUTPUT_CONNECTIVITY, "_connectivity_index.tif", "Connectivity index"),
        (OUTPUT_D_UP, "_connectivity_d_up.tif", "Upslope component D_up"),
        (OUTPUT_D_DOWN, "_connectivity_d_down.tif", "Downslope component D_down"),
    )

    def add_product_parameters(self):
        self.addParameter(
            QgsProcessingParameterNumber(
                self.CONNECTIVITY_AREA_THRESHOLD_M2,
                "Channel drainage-area threshold (m2)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1e6,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.CONNECTIVITY_WINDOW_RADIUS_METRES,
                "Slope and roughness window radius (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=100.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.CONNECTIVITY_DRAINAGE_AREA_METHOD,
                "Drainage area method",
                options=self.METHOD_OPTIONS,
                defaultValue=0,
            )
        )

    def add_advanced_product_parameters(self):
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.CONNECTIVITY_SLOPE_METHOD,
                "Slope method",
                options=self.SLOPE_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.CONNECTIVITY_ROUGHNESS_METHOD,
                "Roughness method",
                options=self.ROUGHNESS_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MFD_EXPONENT,
                "MFD slope exponent (MFD method only)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )

    def product_values(self, parameters, context):
        return {
            "write_connectivity_raster": True,
            "connectivity_drainage_area_method": self.METHOD_OPTIONS[
                self.parameterAsEnum(parameters, self.CONNECTIVITY_DRAINAGE_AREA_METHOD, context)
            ],
            "connectivity_area_threshold_m2": self.parameterAsDouble(
                parameters, self.CONNECTIVITY_AREA_THRESHOLD_M2, context
            ),
            "connectivity_window_radius_metres": self.parameterAsDouble(
                parameters, self.CONNECTIVITY_WINDOW_RADIUS_METRES, context
            ),
            "connectivity_slope_method": self.SLOPE_OPTIONS[
                self.parameterAsEnum(parameters, self.CONNECTIVITY_SLOPE_METHOD, context)
            ],
            "connectivity_roughness_method": self.ROUGHNESS_OPTIONS[
                self.parameterAsEnum(parameters, self.CONNECTIVITY_ROUGHNESS_METHOD, context)
            ],
            "mfd_exponent": self.parameterAsDouble(parameters, self.MFD_EXPONENT, context),
        }

    def name(self):
        return "connectivity_index"

    def displayName(self):
        return "Sediment Connectivity Index"

    def shortHelpString(self):
        return (
            "Runs lsdtt-drainage-analysis: the sediment connectivity index of "
            "Cavalli et al. (2013), IC = log10(D_up / D_down). D_up grows with "
            "the area, slope and smoothness upslope of a pixel; D_down with "
            "the length and roughness of its D8 path down to the nearest "
            "channel, where channels are pixels whose drainage area reaches "
            "the threshold. Higher IC means better connected to the channel "
            "network. D_up and D_down are written too. "
            "Backend path is managed in plugin settings."
        )


class HeightAboveChannelAlgorithm(DrainageAnalysisBase):
    THRESHOLD_CONTRIBUTING_PIXELS = "THRESHOLD_CONTRIBUTING_PIXELS"

    OUTPUT_HEIGHT_ABOVE_CHANNEL = "OUTPUT_HEIGHT_ABOVE_CHANNEL"

    DEFAULT_PREFIX = "lsdtt_hac"
    OUTPUTS = (
        (OUTPUT_HEIGHT_ABOVE_CHANNEL, "_height_above_channel.tif", "Height above channel"),
    )

    def add_product_parameters(self):
        self.addParameter(
            QgsProcessingParameterNumber(
                self.THRESHOLD_CONTRIBUTING_PIXELS,
                "Threshold contributing pixels",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=1000,
            )
        )

    def product_values(self, parameters, context):
        return {
            "write_height_above_channel_raster": True,
            "threshold_contributing_pixels": self.parameterAsInt(
                parameters, self.THRESHOLD_CONTRIBUTING_PIXELS, context
            ),
        }

    def name(self):
        return "height_above_channel"

    def displayName(self):
        return "Height Above Channel"

    def shortHelpString(self):
        return (
            "Runs lsdtt-drainage-analysis: for each pixel, the height above "
            "the first channel pixel reached along its D8 flow path. Channel "
            "pixels are those draining at least the threshold number of "
            "contributing pixels. Backend path is managed in plugin settings."
        )


class SoilTopographicIndexAlgorithm(DrainageAnalysisBase):
    SOIL_LAYER = "SOIL_LAYER"
    LANDUSE_LAYER = "LANDUSE_LAYER"
    KSAT_CSV = "KSAT_CSV"
    SPR_CSV = "SPR_CSV"
    WRITE_STI = "WRITE_STI"
    WRITE_STI_PROXY = "WRITE_STI_PROXY"
    STI_TERM_MODE = "STI_TERM_MODE"
    STI_DRAINAGE_AREA_METHOD = "STI_DRAINAGE_AREA_METHOD"
    STI_WINDOW_RADIUS_METRES = "STI_WINDOW_RADIUS_METRES"
    WRITE_P95 = "WRITE_P95"
    WRITE_P95_JENKS = "WRITE_P95_JENKS"
    STI_P95_JENKS_N_CLASSES = "STI_P95_JENKS_N_CLASSES"
    WRITE_LOOKUP_CSV = "WRITE_LOOKUP_CSV"

    # Advanced
    STI_SPR_TRANSFORM = "STI_SPR_TRANSFORM"
    STI_DEPTH_CAP_CENTIMETRES = "STI_DEPTH_CAP_CENTIMETRES"
    STI_MIN_VALID_DEPTH_CENTIMETRES = "STI_MIN_VALID_DEPTH_CENTIMETRES"
    STI_KS_UNIT_SCALE = "STI_KS_UNIT_SCALE"
    STI_P95_JENKS_MAX_HISTOGRAM_BINS = "STI_P95_JENKS_MAX_HISTOGRAM_BINS"
    STI_RESAMPLE_TO_1M = "STI_RESAMPLE_TO_1M"
    MFD_EXPONENT = "MFD_EXPONENT"

    OUTPUT_STI = "OUTPUT_STI"
    OUTPUT_STI_P95 = "OUTPUT_STI_P95"
    OUTPUT_STI_P95_JENKS = "OUTPUT_STI_P95_JENKS"
    OUTPUT_STI_PROXY = "OUTPUT_STI_PROXY"
    OUTPUT_STI_PROXY_P95 = "OUTPUT_STI_PROXY_P95"
    OUTPUT_STI_PROXY_P95_JENKS = "OUTPUT_STI_PROXY_P95_JENKS"
    OUTPUT_LOOKUP = "OUTPUT_LOOKUP"

    TERM_MODE_OPTIONS = ["ksd", "spr"]
    METHOD_OPTIONS = ["Dinf", "MFD"]
    SPR_TRANSFORM_OPTIONS = ["ln_1_minus_spr", "ln_spr"]
    DEFAULT_PREFIX = "lsdtt_sti"
    # The Jenks suffixes carry the class count (_sti_p95_jenks5.tif); the
    # descriptions here are used for the declared outputs only.
    OUTPUTS = (
        (OUTPUT_STI, "_sti.tif", "Soil topographic index"),
        (OUTPUT_STI_P95, "_sti_p95.tif", "STI, top 5%"),
        (OUTPUT_STI_P95_JENKS, "_sti_p95_jenks.tif", "STI top 5%, Jenks classes"),
        (OUTPUT_STI_PROXY, "_sti_proxy.tif", "STI proxy (runoff-based)"),
        (OUTPUT_STI_PROXY_P95, "_sti_proxy_p95.tif", "STI proxy, top 5%"),
        (OUTPUT_STI_PROXY_P95_JENKS, "_sti_proxy_p95_jenks.tif", "STI proxy top 5%, Jenks classes"),
        (OUTPUT_LOOKUP, "_sti_lookup.csv", "STI lookup diagnostics"),
    )
    HELP_OVERRIDES = {
        SOIL_LAYER: "sti_soil_fname",
        LANDUSE_LAYER: "sti_landuse_fname",
        KSAT_CSV: "sti_ksat_fname",
        SPR_CSV: "sti_spr_fname",
        WRITE_STI: "write_sti_raster",
        WRITE_STI_PROXY: "write_sti_proxy_raster",
        WRITE_P95: "write_sti_95th_percentile_rasters",
        WRITE_P95_JENKS: "write_sti_p95_jenks_rasters",
        WRITE_LOOKUP_CSV: "write_sti_lookup_csv",
        STI_KS_UNIT_SCALE: "sti_ks_unit_scale_to_metres_per_day",
        STI_RESAMPLE_TO_1M: "sti_resample_to_1m_before_calculation",
    }

    def add_product_parameters(self):
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.SOIL_LAYER,
                "Soil series polygons",
                [QgsProcessing.TypeVectorPolygon],
            )
        )
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.LANDUSE_LAYER,
                "Land-use polygons",
                [QgsProcessing.TypeVectorPolygon],
            )
        )
        self.addParameter(
            QgsProcessingParameterFile(
                self.KSAT_CSV,
                "Saturated conductivity profiles CSV (needed for the STI raster and lookup)",
                extension="csv",
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterFile(
                self.SPR_CSV,
                "Standard percentage runoff CSV (needed for the proxy, lookup, or spr mode)",
                extension="csv",
                optional=True,
            )
        )
        self.addParameter(QgsProcessingParameterBoolean(self.WRITE_STI, "Write STI raster", True))
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_STI_PROXY, "Write runoff-based STI proxy raster", False
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.STI_TERM_MODE,
                "Soil term for the STI raster",
                options=self.TERM_MODE_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.STI_DRAINAGE_AREA_METHOD,
                "Drainage area method",
                options=self.METHOD_OPTIONS,
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.STI_WINDOW_RADIUS_METRES,
                "Slope window radius (m)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.001,
                defaultValue=100.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(self.WRITE_P95, "Also write the top 5% of values", False)
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_P95_JENKS, "Classify the top 5% into Jenks classes", False
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.STI_P95_JENKS_N_CLASSES,
                "Number of Jenks classes",
                type=QgsProcessingParameterNumber.Integer,
                minValue=2,
                defaultValue=5,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.WRITE_LOOKUP_CSV, "Write the per-series / land-use lookup CSV", False
            )
        )

    def add_advanced_product_parameters(self):
        self._add_advanced(
            QgsProcessingParameterEnum(
                self.STI_SPR_TRANSFORM,
                "Runoff (SPR) transform",
                options=self.SPR_TRANSFORM_OPTIONS,
                defaultValue=0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.STI_DEPTH_CAP_CENTIMETRES,
                "Profile depth cap (cm)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.001,
                defaultValue=100.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.STI_MIN_VALID_DEPTH_CENTIMETRES,
                "Shortest profile accepted (cm)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.001,
                defaultValue=30.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.STI_KS_UNIT_SCALE,
                "Conductivity unit scale to metres per day",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.STI_P95_JENKS_MAX_HISTOGRAM_BINS,
                "Jenks: maximum histogram bins",
                type=QgsProcessingParameterNumber.Integer,
                minValue=1,
                defaultValue=2048,
            )
        )
        self._add_advanced(
            QgsProcessingParameterBoolean(
                self.STI_RESAMPLE_TO_1M,
                "Resample the DEM to 1 m (cubic) before calculating",
                False,
            )
        )
        self._add_advanced(
            QgsProcessingParameterNumber(
                self.MFD_EXPONENT,
                "MFD slope exponent (MFD method only)",
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0,
                defaultValue=1.0,
            )
        )

    def product_values(self, parameters, context):
        soil_layer = self.parameterAsVectorLayer(parameters, self.SOIL_LAYER, context)
        landuse_layer = self.parameterAsVectorLayer(parameters, self.LANDUSE_LAYER, context)
        if soil_layer is None or landuse_layer is None:
            raise QgsProcessingException("Soil and land-use polygon layers are required.")

        write_sti = self.parameterAsBool(parameters, self.WRITE_STI, context)
        write_proxy = self.parameterAsBool(parameters, self.WRITE_STI_PROXY, context)
        term_mode = self.TERM_MODE_OPTIONS[
            self.parameterAsEnum(parameters, self.STI_TERM_MODE, context)
        ]
        ksat_csv = self.parameterAsFile(parameters, self.KSAT_CSV, context)
        spr_csv = self.parameterAsFile(parameters, self.SPR_CSV, context)
        write_lookup = self.parameterAsBool(parameters, self.WRITE_LOOKUP_CSV, context)
        # Fail here, with the form's wording, rather than in the CLI log. The
        # CLI builds the conductivity table for any STI raster or lookup, and
        # the runoff table for the proxy, spr mode or the lookup.
        if (write_sti or write_lookup) and not ksat_csv:
            raise QgsProcessingException(
                "The STI raster and the lookup CSV need the saturated "
                "conductivity CSV."
            )
        if (write_proxy or write_lookup or (write_sti and term_mode == "spr")) and not spr_csv:
            raise QgsProcessingException(
                "The STI proxy, the lookup CSV and the STI raster in spr mode "
                "need the standard percentage runoff CSV."
            )

        soil_source = soil_layer.source()
        landuse_source = landuse_layer.source()
        values = {
            "write_sti_raster": write_sti,
            "write_sti_proxy_raster": write_proxy,
            "write_sti_95th_percentile_rasters": self.parameterAsBool(
                parameters, self.WRITE_P95, context
            ),
            "write_sti_p95_jenks_rasters": self.parameterAsBool(
                parameters, self.WRITE_P95_JENKS, context
            ),
            "write_sti_lookup_csv": write_lookup,
            "sti_soil_fname": layer_file_path(soil_source),
            "sti_soil_layer_name": source_layer_name(soil_source),
            "sti_landuse_fname": layer_file_path(landuse_source),
            "sti_landuse_layer_name": source_layer_name(landuse_source),
            "sti_term_mode": term_mode,
            "sti_drainage_area_method": self.METHOD_OPTIONS[
                self.parameterAsEnum(parameters, self.STI_DRAINAGE_AREA_METHOD, context)
            ],
            "sti_window_radius_metres": self.parameterAsDouble(
                parameters, self.STI_WINDOW_RADIUS_METRES, context
            ),
            "sti_p95_jenks_n_classes": self.parameterAsInt(
                parameters, self.STI_P95_JENKS_N_CLASSES, context
            ),
            "sti_spr_transform": self.SPR_TRANSFORM_OPTIONS[
                self.parameterAsEnum(parameters, self.STI_SPR_TRANSFORM, context)
            ],
            "sti_depth_cap_centimetres": self.parameterAsDouble(
                parameters, self.STI_DEPTH_CAP_CENTIMETRES, context
            ),
            "sti_min_valid_depth_centimetres": self.parameterAsDouble(
                parameters, self.STI_MIN_VALID_DEPTH_CENTIMETRES, context
            ),
            "sti_ks_unit_scale_to_metres_per_day": self.parameterAsDouble(
                parameters, self.STI_KS_UNIT_SCALE, context
            ),
            "sti_p95_jenks_max_histogram_bins": self.parameterAsInt(
                parameters, self.STI_P95_JENKS_MAX_HISTOGRAM_BINS, context
            ),
            "sti_resample_to_1m_before_calculation": self.parameterAsBool(
                parameters, self.STI_RESAMPLE_TO_1M, context
            ),
            "mfd_exponent": self.parameterAsDouble(parameters, self.MFD_EXPONENT, context),
        }
        # Only pass the tables that were given. The checks above guarantee the
        # CLI never reads a table it was not given, so its defaults (data
        # paths that do not exist next to a user's DEM) are left untouched.
        if ksat_csv:
            values["sti_ksat_fname"] = ksat_csv
        if spr_csv:
            values["sti_spr_fname"] = spr_csv
        return values

    def expected_outputs(self, values):
        p95 = values["write_sti_95th_percentile_rasters"]
        jenks = p95 and values["write_sti_p95_jenks_rasters"]
        jenks_suffix = "_jenks{}.tif".format(values["sti_p95_jenks_n_classes"])
        outputs = []
        for enabled, base, ids in (
            (values["write_sti_raster"], "_sti",
             (self.OUTPUT_STI, self.OUTPUT_STI_P95, self.OUTPUT_STI_P95_JENKS)),
            (values["write_sti_proxy_raster"], "_sti_proxy",
             (self.OUTPUT_STI_PROXY, self.OUTPUT_STI_PROXY_P95,
              self.OUTPUT_STI_PROXY_P95_JENKS)),
        ):
            if not enabled:
                continue
            outputs.append((ids[0], base + ".tif", "STI raster"))
            if p95:
                outputs.append((ids[1], base + "_p95.tif", "STI top 5% raster"))
            if jenks:
                outputs.append((ids[2], base + "_p95" + jenks_suffix, "STI Jenks raster"))
        if values["write_sti_lookup_csv"] and outputs:
            outputs.append((self.OUTPUT_LOOKUP, "_sti_lookup.csv", "STI lookup CSV"))
        return outputs

    def name(self):
        return "soil_topographic_index"

    def displayName(self):
        return "Soil Topographic Index"

    def shortHelpString(self):
        return (
            "Runs lsdtt-drainage-analysis: the soil topographic index, "
            "ln(a / tan b) less a soil term, with the soil term looked up from "
            "soil-series and land-use polygons. In ksd mode the term is "
            "ln(Ks*D) from a table of saturated conductivity profiles; in spr "
            "mode, and for the proxy raster, it comes from the standard "
            "percentage runoff of each soil series. The tables follow the "
            "Scottish soil data (Hutton series codes, SSKIB conductivities); "
            "see the lsdtt-drainage-analysis reference for their columns. "
            "Backend path is managed in plugin settings."
        )
