import unittest

import numpy as np
import pytest

import lsdtt3

pytestmark = [pytest.mark.unit]


class FlatSurfaceTest(unittest.TestCase):
    @pytest.mark.exact_solution
    def test_relief_window_radius_sets_the_window_side(self):
        """relief_window_radius_metres sets the relief window: with radius 1 m a cell two
        cells from a 10 m spike sees no relief, with radius 2 m it sees the spike's 10
        m.
        """
        metadata = lsdtt3.RasterMeta.from_shape(9, 9, cell_size=1.0)
        values = np.zeros((9, 9), dtype=np.float32)
        values[4, 4] = 10.0
        dem = lsdtt3.Raster.from_numpy(
            values, np.ones(values.shape, dtype=bool), metadata
        )
        common = dict(
            use_madk2=True, madk2_window_radius_metres=1.0,
            use_residual_iqr=False, use_nvd=False, min_component_pixels=1,
            max_relief_metres=20.0, n_threads=1,
        )
        relief3 = lsdtt3.find_flat_surfaces(
            dem, relief_window_radius_metres=1.0, **common
        ).relief.to_numpy()
        relief5 = lsdtt3.find_flat_surfaces(
            dem, relief_window_radius_metres=2.0, **common
        ).relief.to_numpy()
        self.assertEqual(relief3[4, 2], np.float32(0.0))
        self.assertEqual(relief5[4, 2], np.float32(10.0))

    @pytest.mark.smoke
    def test_find_flat_surfaces_returns_diagnostics_and_labels(self):
        """find_flat_surfaces returns its roughness, relief, combined mask and labelled
        mask at the DEM's shape.
        """
        metadata = lsdtt3.RasterMeta.from_shape(64, 64, cell_size=1.0)
        rows, cols = np.indices((64, 64), dtype=np.float32)
        dem = lsdtt3.Raster.from_numpy(
            0.02 * rows + 0.01 * cols + 0.02 * np.sin(rows / 4.0),
            np.ones((64, 64), dtype=bool), metadata,
        )

        result = lsdtt3.find_flat_surfaces(
            dem,
            use_madk2=True,
            madk2_window_radius_metres=3.0,
            use_residual_iqr=False,
            use_nvd=False,
            min_component_pixels=1,
            relief_window_radius_metres=2.0,
            max_relief_metres=1.0,
            n_threads=1,
        )

        self.assertEqual(result.madk2_roughness.shape, dem.shape)
        self.assertEqual(result.relief.shape, dem.shape)
        self.assertEqual(result.combined_flat_mask.shape, dem.shape)
        self.assertEqual(result.labelled_flat_mask.shape, dem.shape)

    @pytest.mark.exact_solution
    def test_classify_flat_surfaces_uses_hand_and_channel_ids(self):
        """classify_flat_surfaces classes a low-HAND flat patch containing a channel cell
        as floodplain: one floodplain label in the table.
        """
        metadata = lsdtt3.RasterMeta.from_shape(7, 7, cell_size=1.0)
        valid = np.ones((7, 7), dtype=bool)
        labels = np.zeros((7, 7), dtype=np.float32)
        labels[2:5, 1:6] = 1.0
        hand = np.full((7, 7), 4.0, dtype=np.float32)
        hand[2:5, 1:6] = 0.5
        source_ids = np.full((7, 7), -1.0, dtype=np.float32)
        source_ids[3, 3] = 0.0

        classified = lsdtt3.classify_flat_surfaces(
            lsdtt3.Raster.from_numpy(labels, valid, metadata),
            lsdtt3.Raster.from_numpy(hand, valid, metadata),
            lsdtt3.Raster.from_numpy(source_ids, valid, metadata),
            n_threads=1,
        )

        self.assertEqual(classified.surface_class.shape, (7, 7))
        self.assertEqual(classified.floodplain_mask.shape, (7, 7))
        self.assertEqual(classified.n_floodplain_labels, 1)
        self.assertEqual(classified.table.loc[0, "surface_class"], "floodplain")


if __name__ == "__main__":
    unittest.main()
