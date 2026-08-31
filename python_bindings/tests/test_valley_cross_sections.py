import unittest

import numpy as np

import lsdtt3


class ValleyCrossSectionDiagnosticsTests(unittest.TestCase):
    def test_optional_diagnostic_tables_are_exposed(self):
        """Feature and hydraulic fields survive the native-to-Python boundary."""
        meta = lsdtt3.RasterMeta.from_shape(41, 41, cell_size=10.0)
        rows, cols = np.indices((41, 41))
        data = (np.abs(cols - 20) * 2.0 + rows * 0.05).astype(np.float32)
        dem = lsdtt3.Raster.from_numpy(data, np.ones(data.shape, dtype=bool), meta)
        channel = np.column_stack((np.arange(4, 37), np.full(33, 20, dtype=int)))

        centrelines = lsdtt3.extract_valley_centrelines(
            dem,
            [channel],
            corridor_width_metres=120.0,
            fallback_valley_hand_metres=4.0,
            valley_dilation_pixels=1,
            centreline_loops=1,
        )
        sections = lsdtt3.extract_valley_cross_sections(
            dem,
            centrelines,
            spacing_metres=30.0,
            min_half_length_metres=30.0,
            max_half_length_metres=80.0,
            channel_mask=lsdtt3.scale(dem, 0.0),
            compute_hydraulic_depth=True,
        )

        self.assertGreater(sections.n_sections, 0)
        self.assertTrue(
            {"curvature", "curvature_valid", "hillshade", "hillshade_valid", "is_channel"}
            <= set(sections.samples.columns)
        )
        self.assertEqual(
            list(sections.features.columns),
            ["source_id", "station_id", "feature", "offset_metres"],
        )
        self.assertTrue(
            {"break_id", "stage_metres", "depth_above_thalweg_metres"}
            <= set(sections.hydraulic_breaks.columns)
        )
        self.assertTrue(
            {"width_metres", "area_metres2", "hydraulic_depth_metres"}
            <= set(sections.hydraulic_curve.columns)
        )


if __name__ == "__main__":
    unittest.main()
