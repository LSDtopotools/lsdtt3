import unittest

import numpy as np
import pytest

import lsdtt3

pytestmark = [pytest.mark.unit, pytest.mark.smoke]


class TerrainFeatureTest(unittest.TestCase):
    def dem(self, size=64):
        metadata = lsdtt3.RasterMeta.from_shape(size, size, cell_size=1.0)
        rows, cols = np.indices((size, size), dtype=np.float32)
        values = 0.01 * rows + 0.005 * cols + 0.03 * np.random.default_rng(8).standard_normal((size, size))
        return lsdtt3.Raster.from_numpy(values, np.ones((size, size), dtype=bool), metadata)

    def test_marsh_and_road_feature_results(self):
        dem = self.dem()
        marsh = lsdtt3.find_marsh_platform(
            dem, max_scarp_growth_iterations=2, max_platform_growth_iterations=2,
            n_threads=1,
        )
        roads = lsdtt3.find_roads(dem, smoothing_window_metres=3.0, n_threads=1)
        smooth_roads = lsdtt3.extract_roads_smooth(
            dem, roughness_radius_metres=1.0, path_min_length_pixels=2,
            min_component_pixels=1, n_threads=1,
        )

        self.assertEqual(marsh.platform_mask.shape, dem.shape)
        self.assertEqual(marsh.scarp_mask.shape, dem.shape)
        self.assertEqual(roads.mask.shape, dem.shape)
        self.assertEqual(smooth_roads.road_mask.shape, dem.shape)

    @pytest.mark.property
    @pytest.mark.exact_solution
    def test_ridgeline_and_valley_centreline(self):
        metadata = lsdtt3.RasterMeta.from_shape(48, 48, cell_size=1.0)
        rows, cols = np.indices((48, 48), dtype=np.float32)
        valley_dem = lsdtt3.Raster.from_numpy(
            0.1 * rows + 0.04 * (cols - 24.0) ** 2,
            np.ones((48, 48), dtype=bool), metadata,
        )
        pixels = np.column_stack((np.full(20, 10), np.arange(10, 30)))
        ridge = lsdtt3.Ridgeline(valley_dem, pixels)
        ridge.fill_elevation(valley_dem)
        ridge.fit_profile()

        channel = np.column_stack((np.arange(4, 44), np.full(40, 24)))
        valley_values = (np.abs(cols - 24) <= 6).astype(np.float32)
        valley_mask = lsdtt3.Raster.from_numpy(
            valley_values, np.ones(valley_values.shape, dtype=bool), metadata,
        )
        valleys = lsdtt3.extract_valley_centrelines(
            valley_dem, [channel], corridor_width_metres=15.0,
            centreline_loops=1, n_threads=1,
            valley_mask=valley_mask,
        )

        self.assertTrue(ridge.is_four_connected())
        self.assertEqual(len(ridge.elevation), len(pixels))
        self.assertEqual(valleys.valley_mask.shape, valley_dem.shape)
        self.assertEqual(valleys.centreline_pixels(0).shape[1], 2)
        np.testing.assert_allclose(
            lsdtt3.fit_valley_baseline([0.0, 1.0, 2.0], [10.0, np.nan, 7.0]),
            [10.0, 8.5, 7.0],
        )


if __name__ == "__main__":
    unittest.main()
