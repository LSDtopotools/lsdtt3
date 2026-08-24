import unittest

import numpy as np

import lsdtt3


class RasterTest(unittest.TestCase):
    def raster(self):
        metadata = lsdtt3.RasterMeta.from_shape(
            5, 5, x_min=100.0, y_min=200.0, cell_size=10.0
        )
        data = np.arange(25, dtype=np.float32).reshape(5, 5)
        return lsdtt3.Raster.from_numpy(data, np.ones((5, 5), dtype=bool), metadata)

    def test_from_numpy_copies_data_mask_and_metadata(self):
        raster = self.raster()

        np.testing.assert_array_equal(raster.to_numpy(), np.arange(25, dtype=np.float32).reshape(5, 5))
        np.testing.assert_array_equal(raster.mask_numpy(), np.ones((5, 5), dtype=bool))
        self.assertEqual(raster.shape, (5, 5))
        self.assertEqual(raster.metadata.x_min, 100.0)
        self.assertEqual(raster.metadata.y_min, 200.0)
        self.assertEqual(raster.metadata.cell_size, 10.0)

    def test_fill_interior_nodata_preserves_input_and_edge_nodata(self):
        raster = self.raster()
        with_holes = lsdtt3.set_nodata_rectangle(raster, 2, 3, 2, 3)
        with_holes = lsdtt3.set_nodata_rectangle(with_holes, 0, 1, 0, 1)

        filled = lsdtt3.fill_interior_nodata(with_holes, method="iterative")

        self.assertFalse(with_holes.mask_numpy()[2, 2])
        self.assertFalse(with_holes.mask_numpy()[0, 0])
        self.assertTrue(filled.mask_numpy()[2, 2])
        self.assertFalse(filled.mask_numpy()[0, 0])

    def test_clip_to_maximum_rectangle_updates_shape_and_extent(self):
        clipped = lsdtt3.clip_to_maximum_rectangle(
            lsdtt3.set_nodata_rectangle(self.raster(), 0, 5, 0, 1)
        )

        self.assertEqual(clipped.shape, (5, 4))
        self.assertEqual(clipped.metadata.x_min, 110.0)
        self.assertEqual(clipped.metadata.y_min, 200.0)
        self.assertEqual(clipped.valid_pixel_count, 20)


if __name__ == "__main__":
    unittest.main()
