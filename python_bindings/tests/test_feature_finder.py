import unittest

import numpy as np
import pytest

import lsdtt3

pytestmark = [pytest.mark.unit, pytest.mark.smoke]


class FeatureFinderTest(unittest.TestCase):
    def raster(self):
        metadata = lsdtt3.RasterMeta.from_shape(64, 64, cell_size=1.0)
        rows, cols = np.indices((64, 64), dtype=np.float32)
        values = (
            0.03 * cols + 0.01 * rows
            + 10.0 * np.exp(-((cols - 32.0) ** 2 + (rows - 32.0) ** 2) / 200.0)
            - 4.0 * np.exp(-(cols - 30.0) ** 2 / 10.0)
        )
        return lsdtt3.Raster.from_numpy(
            values, np.ones_like(values, dtype=bool), metadata
        )

    def test_feature_detection_results(self):
        """relief_ratio, frangi_filter, find_channel_candidates (3 scales), find_drains and
        find_valleys return their rasters at the DEM's shape.
        """
        raster = self.raster()
        ratio = lsdtt3.relief_ratio(raster, 2.0, 8.0)
        frangi = lsdtt3.frangi_filter(raster, 1.0, 3.0, 3, target="valleys", n_threads=1)
        channels = lsdtt3.find_channel_candidates(raster, 3.0, 10.0, 3, n_threads=1)
        drains = lsdtt3.find_drains(
            raster, 1.0, 3.0, 3, accumulation_curvature_lambda=4.0, n_threads=1
        )
        valleys = lsdtt3.find_valleys(raster, 10.0, n_threads=1)

        self.assertEqual(ratio.shape, raster.shape)
        self.assertEqual(frangi.response.shape, raster.shape)
        self.assertEqual(frangi.best_scale_metres.shape, raster.shape)
        self.assertEqual(channels.score.shape, raster.shape)
        self.assertEqual(channels.mask.shape, raster.shape)
        self.assertEqual(len(channels.scales_metres), 3)
        self.assertEqual(len(channels.curvature_raw_by_scale), 3)
        self.assertEqual(drains.drain_mask.shape, raster.shape)
        self.assertEqual(drains.frangi_response.shape, raster.shape)
        self.assertEqual(valleys.valley_mask.shape, raster.shape)
        self.assertEqual(valleys.shape_class.shape, raster.shape)

    def test_mountain_front_and_outlets(self):
        """identify_mountain_front with two relief scales returns its masks, fitted
        surface, residual and one relief raster per scale; mountain_front_outlets
        returns a raster at the DEM's shape.
        """
        raster = self.raster()
        result = lsdtt3.identify_mountain_front(
            raster,
            relief_scales=[[3.0, 3.0], [8.0, 8.0]],
            fit_method="plane_ols",
            return_diagnostics=True,
            n_threads=1,
        )

        self.assertEqual(result.low_ground_mask.shape, raster.shape)
        self.assertEqual(result.fitted_surface.shape, raster.shape)
        self.assertEqual(result.residual.shape, raster.shape)
        self.assertEqual(result.mountain_mask.shape, raster.shape)
        self.assertEqual(len(result.relief_rasters), 2)

        contributing = lsdtt3.Raster.from_numpy(
            np.arange(64 * 64, dtype=np.float32).reshape(64, 64),
            np.ones(raster.shape, dtype=bool), raster.metadata,
        )
        outlets = lsdtt3.mountain_front_outlets(
            result.mountain_mask, contributing,
            min_contributing_pixels=100
        )
        self.assertEqual(outlets.shape, raster.shape)

    @pytest.mark.negative
    def test_road_flow_canonical_keywords(self):
        """find_road_flow_discrimination runs with the current keyword names, and the old
        names flow_accumulation_threshold_pixels and dem_is_filled are refused
        (TypeError).
        """
        raster = self.raster()
        result = lsdtt3.find_road_flow_discrimination(
            raster,
            regional_smoothing_metres=5.0,
            threshold_contributing_pixels=5,
            intermediate_window_metres=3.0,
            roughness_radius_metres=1.0,
            raster_is_filled=True,
            n_threads=1,
        )
        self.assertEqual(result.contributing_pixels.shape, raster.shape)
        self.assertEqual(result.channel_mask.shape, raster.shape)

        with self.assertRaises(TypeError):
            lsdtt3.find_road_flow_discrimination(
                raster, flow_accumulation_threshold_pixels=5
            )
        with self.assertRaises(TypeError):
            lsdtt3.find_road_flow_discrimination(raster, dem_is_filled=True)


if __name__ == "__main__":
    unittest.main()
