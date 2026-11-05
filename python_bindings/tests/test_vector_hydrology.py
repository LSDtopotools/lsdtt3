"""Focused tests for vector hydrology bindings.

Authors
-------
B. Gailleton
"""

import unittest

import geopandas as gpd
import numpy as np
import pytest
from pyproj import CRS
from shapely.geometry import LineString, Polygon

import lsdtt3

pytestmark = [pytest.mark.unit]


CRS_UTM = "EPSG:32633"


class TestVectorHydrology(unittest.TestCase):
    @pytest.mark.exact_solution
    def test_repair_water_lines(self):
        water = gpd.GeoDataFrame(
            geometry=[LineString([(0, 0), (10, 0)]), LineString([(15, 0), (25, 0)])], crs=CRS_UTM
        )
        roads = gpd.GeoDataFrame(geometry=[LineString([(12.5, -5), (12.5, 5)])], crs=CRS_UTM)
        result = lsdtt3.repair_water_lines(water, roads, max_gap_length_metres=10.0)
        self.assertEqual(result.summary["n_bridges_accepted"], 1)
        self.assertEqual(len(result.lines), 3)
        self.assertEqual(len(result.bridges), 1)

    @pytest.mark.exact_solution
    @pytest.mark.smoke
    def test_water_centrelines(self):
        areas = gpd.GeoDataFrame(
            geometry=[Polygon([(0, 0), (100, 0), (100, 20), (0, 20)])], crs=CRS_UTM
        )
        tributaries = gpd.GeoDataFrame(
            geometry=[LineString([(-20, 10), (0, 10)]), LineString([(100, 10), (120, 10)])], crs=CRS_UTM
        )
        result = lsdtt3.water_centrelines(
            areas, tributaries, cell_size_metres=2.0, min_component_length_metres=0.0
        )
        self.assertEqual(result.summary["n_polygons"], 1)
        self.assertFalse(result.lines.empty)
        self.assertIn("role", result.lines.columns)

    @pytest.mark.property
    @pytest.mark.smoke
    def test_burn_streams_returns_conditioned_copy(self):
        crs_wkt = CRS.from_epsg(32633).to_wkt()
        meta = lsdtt3.RasterMeta.from_shape(5, 5, cell_size=10.0, crs_wkt=crs_wkt, epsg_code=32633)
        values = np.array([
            [100, 100, 100, 100, 100],
            [100, 100, 100, 100, 100],
            [100, 90, 95, 70, 60],
            [100, 100, 100, 100, 100],
            [100, 100, 100, 100, 100],
        ], dtype=np.float32)
        dem = lsdtt3.Raster.from_numpy(values, np.ones_like(values, dtype=bool), meta)
        lines = gpd.GeoDataFrame(geometry=[LineString([(5, 25), (45, 25)])], crs=CRS_UTM)
        result = lsdtt3.burn_streams(dem, lines)
        self.assertEqual(result.raster.shape, dem.shape)
        self.assertFalse(result.profiles.empty)
        self.assertTrue(np.all(result.raster.to_numpy() <= values))


if __name__ == "__main__":
    unittest.main()
