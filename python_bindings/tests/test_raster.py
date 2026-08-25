import unittest

import geopandas as gpd
import numpy as np
from shapely.geometry import Point, Polygon

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

    def test_masking_and_algebra_return_new_rasters(self):
        raster = self.raster()
        mask = lsdtt3.threshold_above(raster, 12.0)
        kept = lsdtt3.mask_by_binary(raster, mask)
        low_values_removed = lsdtt3.mask_to_nodata(raster, 12.0, below=True)

        self.assertEqual(mask.to_numpy()[2, 2], 1.0)
        self.assertEqual(mask.to_numpy()[0, 0], 0.0)
        self.assertFalse(kept.mask_numpy()[0, 0])
        self.assertTrue(kept.mask_numpy()[2, 2])
        self.assertFalse(low_values_removed.mask_numpy()[0, 0])
        np.testing.assert_array_equal(
            lsdtt3.scale(raster, 2.0).to_numpy(), raster.to_numpy() * 2.0
        )
        np.testing.assert_array_equal(
            lsdtt3.add(raster, raster).to_numpy(), raster.to_numpy() * 2.0
        )
        normalised = lsdtt3.normalize_unit_interval(raster).to_numpy()
        self.assertEqual(normalised.min(), 0.0)
        self.assertEqual(normalised.max(), 1.0)
        np.testing.assert_array_equal(raster.to_numpy(), np.arange(25, dtype=np.float32).reshape(5, 5))

    def test_preprocessing_and_surface_metric_outputs(self):
        metadata = lsdtt3.RasterMeta.from_shape(5, 5)
        depression_values = np.full((5, 5), 10.0, dtype=np.float32)
        depression_values[2, 2] = 0.0
        depression = lsdtt3.Raster.from_numpy(
            depression_values, np.ones((5, 5), dtype=bool), metadata
        )

        filled = lsdtt3.fill_pits(depression, min_slope=0.0, n_threads=1)
        depth = lsdtt3.fill_depth(depression, min_slope=0.0, n_threads=1)
        breached = lsdtt3.breach_depressions(depression, n_threads=1)
        carve = lsdtt3.carve_depth(depression, n_threads=1)

        self.assertEqual(filled.to_numpy()[2, 2], 10.0)
        self.assertEqual(depth.to_numpy()[2, 2], 10.0)
        self.assertEqual(depression.to_numpy()[2, 2], 0.0)
        self.assertEqual(breached.shape, depression.shape)
        self.assertTrue(np.all(carve.to_numpy()[carve.mask_numpy()] >= 0.0))

        plane = lsdtt3.Raster.from_numpy(
            np.add.outer(np.arange(5, dtype=np.float32), np.arange(5, dtype=np.float32)),
            np.ones((5, 5), dtype=bool), metadata,
        )
        self.assertTrue(lsdtt3.hillshade(plane).mask_numpy()[2, 2])
        self.assertTrue(lsdtt3.multidirectional_hillshade(plane).mask_numpy()[2, 2])
        self.assertGreater(lsdtt3.slope(plane).to_numpy()[2, 2], 0.0)
        self.assertTrue(lsdtt3.aspect(plane, n_threads=1).mask_numpy()[2, 2])
        self.assertTrue(lsdtt3.plan_curvature(plane).mask_numpy()[2, 2])
        self.assertTrue(lsdtt3.total_curvature(plane, radius_metres=1.0).mask_numpy()[2, 2])

    def test_boundary_sinks_and_least_cost_breaching(self):
        metadata = lsdtt3.RasterMeta.from_shape(5, 5)
        values = np.full((5, 5), 10.0, dtype=np.float32)
        values[2, 2] = 0.0
        depression = lsdtt3.Raster.from_numpy(
            values, np.ones((5, 5), dtype=bool), metadata
        )

        self.assertEqual(
            lsdtt3.make_boundary_conditions("closed", "periodic", "closed", "periodic"),
            "npnp",
        )
        outlets = lsdtt3.make_outlet_mask(depression, south=True)
        np.testing.assert_array_equal(outlets.to_numpy()[-1], np.ones(5, dtype=np.float32))
        self.assertEqual(
            lsdtt3.fill_pits_with_outlets(
                depression, outlets, min_slope=0.0, n_threads=1
            ).to_numpy()[2, 2],
            10.0,
        )

        points = np.array([[2, 2]], dtype=np.int64)
        sink_mask = lsdtt3.mask_from_points(depression, points)
        self.assertEqual(
            lsdtt3.fill_pits_with_sinks(
                depression, points, min_slope=0.0, n_threads=1
            ).to_numpy()[2, 2],
            0.0,
        )
        self.assertEqual(
            lsdtt3.fill_pits_with_sinks(
                depression, sink_mask, min_slope=0.0, n_threads=1
            ).to_numpy()[2, 2],
            0.0,
        )
        self.assertEqual(
            lsdtt3.fill_pits_boundary_aware(
                depression, "onnn", min_slope=0.0, n_threads=1
            ).to_numpy()[2, 2],
            10.0,
        )
        self.assertEqual(
            lsdtt3.breach_least_cost(depression, n_threads=1).shape,
            depression.shape,
        )
        self.assertTrue(np.all(
            lsdtt3.least_cost_carve_depth(depression, n_threads=1).to_numpy() >= 0.0
        ))

    def test_geopandas_reprojection_mosaic_clipping_and_sampling(self):
        crs = gpd.GeoSeries([], crs=32612).crs
        metadata = lsdtt3.RasterMeta.from_shape(
            5, 5, x_min=500000.0, y_min=4400000.0, cell_size=10.0,
            crs_wkt=crs.to_wkt(),
        )
        raster = lsdtt3.Raster.from_numpy(
            np.arange(25, dtype=np.float32).reshape(5, 5),
            np.ones((5, 5), dtype=bool), metadata,
        )
        polygon = gpd.GeoDataFrame(
            geometry=[Polygon([
                (500000.0, 4400000.0), (500050.0, 4400000.0),
                (500050.0, 4400050.0), (500000.0, 4400050.0),
            ])],
            crs=crs,
        )

        polygon_mask = lsdtt3.rasterize_polygons(raster, polygon)
        self.assertEqual(polygon_mask.valid_pixel_count, 25)
        np.testing.assert_array_equal(polygon_mask.to_numpy(), np.ones((5, 5), dtype=np.float32))
        self.assertEqual(lsdtt3.clip_to_polygons(raster, polygon).shape, (5, 5))

        sampled = lsdtt3.sample_points(
            raster,
            gpd.GeoDataFrame(geometry=[Point(500005.0, 4400045.0)], crs=crs),
            column="elevation",
        )
        self.assertEqual(sampled.loc[0, "elevation"], 0.0)
        np.testing.assert_array_equal(
            lsdtt3.mosaic([raster, raster]).to_numpy(), raster.to_numpy()
        )
        regridded = lsdtt3.reproject_to_grid(raster, metadata)
        np.testing.assert_array_equal(regridded.to_numpy(), raster.to_numpy())


if __name__ == "__main__":
    unittest.main()
