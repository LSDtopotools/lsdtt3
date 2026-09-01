import unittest

import geopandas as gpd
import numpy as np
from pyproj import CRS
from shapely.geometry import LineString

import lsdtt3


class LinearStatsTests(unittest.TestCase):
    def test_line_sampling_preserves_geometry_and_adds_statistics(self):
        meta = lsdtt3.RasterMeta.from_shape(
            5, 5, x_min=0.0, y_min=0.0, cell_size=10.0,
            crs_wkt=CRS.from_epsg(32633).to_wkt(), epsg_code=32633,
        )
        raster = lsdtt3.Raster.from_numpy(
            np.arange(25, dtype=np.float32).reshape(5, 5),
            np.ones((5, 5), dtype=bool),
            meta,
        )
        lines = gpd.GeoDataFrame(
            {"name": ["centre"]},
            geometry=[LineString([(5, 25), (45, 25)])],
            crs="EPSG:32633",
        )
        result = lsdtt3.linear_stats(raster, lines, spacing_metres=10.0)
        self.assertEqual(result.loc[0, "name"], "centre")
        self.assertEqual(result.loc[0, "n_values"], 5)
        self.assertAlmostEqual(result.loc[0, "length_metres"], 40.0)
        self.assertAlmostEqual(result.loc[0, "mean"], 12.0)


class BearingComparisonTests(unittest.TestCase):
    def test_rotation_invariant_and_axial_comparisons(self):
        first = np.array([5.0, 25.0, 60.0, 150.0])
        shifted = first + 30.0
        ordinary = lsdtt3.compare_bearings(first, shifted)
        invariant = lsdtt3.compare_bearings(first, shifted, rotation_invariant=True)
        axial = lsdtt3.compare_bearings(np.array([10.0, 50.0]), np.array([190.0, 230.0]), axial=True)
        self.assertGreater(ordinary["circular_wasserstein_degrees"], 0.0)
        self.assertAlmostEqual(invariant["circular_wasserstein_degrees"], 0.0, places=10)
        self.assertAlmostEqual(axial["circular_wasserstein_degrees"], 0.0, places=10)


class BasinOutlineRidgelineTests(unittest.TestCase):
    def test_selected_basin_outlines_are_split_into_georeferenced_ridges(self):
        meta = lsdtt3.RasterMeta.from_shape(9, 9, cell_size=10.0)
        rows, cols = np.indices((9, 9))
        dem = lsdtt3.Raster.from_numpy(
            (100.0 - 5.0 * rows + np.abs(cols - 4)).astype(np.float32),
            np.ones((9, 9), dtype=bool), meta,
        )
        flow = lsdtt3.FlowInfo(
            dem, boundary_conditions="nnnn", forced_outlets=np.array([[8, 4]]),
            allow_pits_as_outlets=True,
        )
        ridges = flow.basins(np.array([[8, 4]])).outline_ridgelines(dem)
        self.assertEqual(len(ridges), 2)
        self.assertTrue(
            {"basin_id", "ridgeline_id", "length_metres", "geometry"}
            <= set(ridges.columns)
        )

    def test_ridgeline_extrema_include_network_topology_diagnostics(self):
        meta = lsdtt3.RasterMeta.from_shape(21, 21, cell_size=10.0)
        rows, cols = np.indices((21, 21))
        dem = lsdtt3.Raster.from_numpy(
            (200.0 - 4.0 * rows + 2.0 * np.abs(cols - 10)).astype(np.float32),
            np.ones((21, 21), dtype=bool), meta,
        )
        flow = lsdtt3.FlowInfo(
            dem, boundary_conditions="nnnn", forced_outlets=np.array([[20, 10]]),
            allow_pits_as_outlets=True,
        )
        network = flow.channel_network(flow.channel_scope(5))
        extrema = lsdtt3.basin_ridgeline_extrema(
            flow.basins(np.array([[20, 10]])), network, dem,
            detrend_smoothing_pixels=1,
        )
        self.assertGreater(len(extrema), 0)
        self.assertTrue(
            {"type", "distance_to_tdj_metres", "ridgeline_slope_fit", "geometry"}
            <= set(extrema.columns)
        )


if __name__ == "__main__":
    unittest.main()
