import os

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from .algorithms.basin_extraction import BasinExtractionAlgorithm
from .algorithms.basin_stats import BasinStatsAlgorithm
from .algorithms.channel_extraction import ChannelExtractionAlgorithm
from .algorithms.channel_network_analysis import ChannelNetworkAnalysisAlgorithm
from .algorithms.chi_analysis import ChiAnalysisAlgorithm
from .algorithms.clip_raster_to_polygon import ClipRasterToPolygonAlgorithm
from .algorithms.clip_vector_to_raster import ClipVectorToRasterAlgorithm
from .algorithms.drainage_analysis import (
    ConnectivityIndexAlgorithm,
    DrainageAreaAlgorithm,
    HeightAboveChannelAlgorithm,
    SoilTopographicIndexAlgorithm,
    TopographicWetnessIndexAlgorithm,
)
from .algorithms.mosaic_rasters import MosaicRastersAlgorithm
from .algorithms.raster_preprocessing import RasterPreprocessingAlgorithm
from .algorithms.sample_rasters_to_points import SampleRastersToPointsAlgorithm
from .algorithms.side_tributary_basins import SideTributaryBasinsAlgorithm
from .algorithms.surface_metrics import SurfaceMetricsAlgorithm
from .algorithms.swath_profile import SwathProfileAlgorithm


class LSDTTProvider(QgsProcessingProvider):
    def loadAlgorithms(self):
        self.addAlgorithm(RasterPreprocessingAlgorithm())
        self.addAlgorithm(ChannelExtractionAlgorithm())
        self.addAlgorithm(ChiAnalysisAlgorithm())
        self.addAlgorithm(BasinExtractionAlgorithm())
        self.addAlgorithm(BasinStatsAlgorithm())
        self.addAlgorithm(SurfaceMetricsAlgorithm())
        self.addAlgorithm(ChannelNetworkAnalysisAlgorithm())
        self.addAlgorithm(ClipRasterToPolygonAlgorithm())
        self.addAlgorithm(ClipVectorToRasterAlgorithm())
        self.addAlgorithm(MosaicRastersAlgorithm())
        self.addAlgorithm(SampleRastersToPointsAlgorithm())
        self.addAlgorithm(SideTributaryBasinsAlgorithm())
        self.addAlgorithm(SwathProfileAlgorithm())
        self.addAlgorithm(DrainageAreaAlgorithm())
        self.addAlgorithm(TopographicWetnessIndexAlgorithm())
        self.addAlgorithm(ConnectivityIndexAlgorithm())
        self.addAlgorithm(HeightAboveChannelAlgorithm())
        self.addAlgorithm(SoilTopographicIndexAlgorithm())

    def id(self):
        return "lsdtt"

    def name(self):
        return "lsdtt3"

    def longName(self):
        return "lsdtt3 (lsdtt)"

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(__file__), "icon.png"))
