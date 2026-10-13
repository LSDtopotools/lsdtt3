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

    def test_raster_hierarchy_restricts_and_prolongs(self):
        metadata = lsdtt3.RasterMeta.from_shape(8, 8, cell_size=10.0)
        raster = lsdtt3.Raster.from_numpy(
            np.arange(64, dtype=np.float32).reshape(8, 8),
            np.ones((8, 8), dtype=bool), metadata,
        )
        hierarchy = lsdtt3.RasterHierarchy(raster)

        self.assertEqual(hierarchy.n_levels, 4)
        self.assertEqual(hierarchy.level_factor(2), 4)
        self.assertEqual(hierarchy.level_shape(2), (2, 2))
        coarse = hierarchy.restrict_to(2)
        self.assertEqual(coarse.shape, (2, 2))
        self.assertEqual(hierarchy.active_level, 2)
        restored = hierarchy.prolong_to(0)
        self.assertEqual(restored.shape, (8, 8))
        factor = hierarchy.add_level_for_factor(8)
        self.assertEqual(hierarchy.level_factor(factor), 8)

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

    def test_polyfit_metrics_returns_requested_rasters_from_one_fit(self):
        metadata = lsdtt3.RasterMeta.from_shape(7, 7)
        rows, cols = np.indices((7, 7), dtype=np.float32)
        paraboloid = lsdtt3.Raster.from_numpy(
            (rows - 3.0) ** 2 + (cols - 3.0) ** 2,
            np.ones((7, 7), dtype=bool), metadata,
        )
        names = [
            "fitted_surface", "dzdx", "dzdy", "d2zdx2", "d2zdxdy", "d2zdy2",
            "slope", "slope_degrees", "aspect_radians", "aspect_degrees",
            "laplacian", "hessian_determinant", "mean_curvature",
            "gaussian_curvature", "maximum_curvature", "minimum_curvature",
            "curvedness", "shape_index", "profile_curvature", "plan_curvature",
            "tangential_curvature", "difference_curvature",
            "horizontal_excess_curvature", "vertical_excess_curvature",
            "accumulation_curvature", "ring_curvature", "rotor",
        ]

        metrics = lsdtt3.polyfit_metrics(paraboloid, names, lambda_metres=3.0, n_threads=1)

        self.assertEqual(set(metrics), set(names))
        self.assertTrue(all(metric.shape == paraboloid.shape for metric in metrics.values()))
        self.assertTrue(metrics["mean_curvature"].mask_numpy()[3, 3])
        self.assertGreater(metrics["mean_curvature"].to_numpy()[3, 3], 0.0)
        self.assertGreater(metrics["gaussian_curvature"].to_numpy()[3, 3], 0.0)
        self.assertEqual(
            metrics["aspect_degrees"].to_numpy()[2, 3],
            lsdtt3.aspect(paraboloid, lambda_metres=3.0, n_threads=1).to_numpy()[2, 3],
        )
        with self.assertRaisesRegex(ValueError, "unknown metric"):
            lsdtt3.polyfit_metrics(paraboloid, ["not_a_metric"])

    def test_surface_shape_and_fundamental_form_principal_curvatures(self):
        metadata = lsdtt3.RasterMeta.from_shape(7, 7)
        rows, cols = np.indices((7, 7), dtype=np.float32)
        paraboloid = lsdtt3.Raster.from_numpy(
            (rows - 3.0) ** 2 + (cols - 3.0) ** 2,
            np.ones((7, 7), dtype=bool), metadata,
        )

        computed = lsdtt3.classify_surface_shape(
            paraboloid, lambda_metres=3.0, n_threads=1
        )
        polyfit_principal = lsdtt3.polyfit_metrics(
            paraboloid, ["maximum_curvature", "minimum_curvature"],
            lambda_metres=3.0, n_threads=1,
        )
        reused = lsdtt3.classify_surface_shape(
            polyfit_principal["maximum_curvature"],
            polyfit_principal["minimum_curvature"],
        )
        self.assertEqual(set(computed), {"shape_class", "gaussian_curvature", "mean_curvature"})
        np.testing.assert_array_equal(
            computed["shape_class"].to_numpy(), reused["shape_class"].to_numpy()
        )
        self.assertEqual(computed["shape_class"].to_numpy()[3, 3], -3.0)

        fundamental = lsdtt3.fundamental_form_principal_curvatures(paraboloid, n_threads=1)
        self.assertEqual(set(fundamental), {"maximum_curvature", "minimum_curvature"})
        self.assertTrue(fundamental["maximum_curvature"].mask_numpy()[3, 3])
        self.assertGreater(fundamental["minimum_curvature"].to_numpy()[3, 3], 0.0)

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

    def test_breach_least_cost_edge_policy_and_errors(self):
        # The 7x7 bowl of tests/data/tiny_breach_lc_dem.asc: pit 1.0 at (3, 3),
        # a 4 saddle at (3, 5) and a 3 on the east border at (3, 6).
        values = np.array([
            [9, 9, 9, 9, 9, 9, 9],
            [9, 6, 6, 6, 6, 6, 9],
            [9, 6, 2, 2, 2, 6, 9],
            [9, 6, 2, 1, 2, 4, 3],
            [9, 6, 2, 2, 2, 6, 9],
            [9, 6, 6, 6, 6, 6, 9],
            [9, 9, 9, 9, 9, 9, 9],
        ], dtype=np.float32)
        metadata = lsdtt3.RasterMeta.from_shape(7, 7, x_min=0.0, y_min=0.0, cell_size=10.0)
        dem = lsdtt3.Raster.from_numpy(values, np.ones((7, 7), dtype=bool), metadata)

        # Default edge_policy "outlet": breached off the east edge.
        outlet = lsdtt3.breach_least_cost(dem, flat_increment=0.125).to_numpy()
        expected = values.copy()
        expected[3, 3] = 1.875
        expected[2, 4] = 1.75
        expected[3, 5] = 1.625
        expected[3, 6] = 1.5
        np.testing.assert_array_equal(outlet, expected)

        # "whitebox": the pit is unsolved; without the remainder fill only the
        # pit raise remains.
        wb = lsdtt3.breach_least_cost(
            dem, flat_increment=0.125, fill_remainder=False, edge_policy="whitebox"
        ).to_numpy()
        expected_wb = values.copy()
        expected_wb[3, 3] = 1.875
        np.testing.assert_array_equal(wb, expected_wb)

        carve = lsdtt3.least_cost_carve_depth(dem, flat_increment=0.125).to_numpy()
        self.assertEqual(float(carve[3, 5]), 2.375)

        with self.assertRaisesRegex(ValueError, "edge_policy must be one of outlet, whitebox"):
            lsdtt3.breach_least_cost(dem, edge_policy="Outlet")
        with self.assertRaisesRegex(ValueError, "max_dist must be between 1 and 32767"):
            lsdtt3.breach_least_cost(dem, max_distance=0)
        with self.assertRaisesRegex(ValueError, "fill_min_slope"):
            lsdtt3.breach_least_cost(dem, fill_min_slope=-1.0)

    def test_breach_least_cost_max_distance_default_follows_edge_policy(self):
        # A pit at (2, 2) in a 2 m corridor walled by 9 m; the only cell low
        # enough to drain to is a nodata hole 120 steps along the corridor.
        # In whitebox mode the edge is no outlet, so max_distance decides:
        # None must behave as 100 (unsolved) there, and as 1000 in outlet mode.
        rows, cols = 5, 130
        values = np.full((rows, cols), 9.0, dtype=np.float32)
        values[2, 1:125] = 2.0
        values[2, 2] = 1.0
        valid = np.ones((rows, cols), dtype=bool)
        valid[2, 122] = False
        metadata = lsdtt3.RasterMeta.from_shape(rows, cols, x_min=0.0, y_min=0.0, cell_size=10.0)
        dem = lsdtt3.Raster.from_numpy(values, valid, metadata)

        def run(**kw):
            return lsdtt3.breach_least_cost(
                dem, flat_increment=0.125, fill_remainder=False, **kw).to_numpy()

        wb_default = run(edge_policy="whitebox")
        np.testing.assert_array_equal(wb_default, run(edge_policy="whitebox", max_distance=100))
        self.assertFalse(np.array_equal(
            wb_default, run(edge_policy="whitebox", max_distance=1000)))
        np.testing.assert_array_equal(run(), run(max_distance=1000))

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

    def test_in_memory_resampling(self):
        metadata = lsdtt3.RasterMeta.from_shape(4, 4, cell_size=1.0)
        values = np.arange(1, 17, dtype=np.float32).reshape(4, 4)
        raster = lsdtt3.Raster.from_numpy(values, np.ones((4, 4), bool), metadata)

        expected = {
            "nearest": np.array([[6, 8], [14, 16]], dtype=np.float32),
            "linear": np.array([[3.5, 5.5], [11.5, 13.5]], dtype=np.float32),
            "mean": np.array([[3.5, 5.5], [11.5, 13.5]], dtype=np.float32),
            "sum": np.array([[14, 22], [46, 54]], dtype=np.float32),
            "min": np.array([[1, 3], [9, 11]], dtype=np.float32),
            "max": np.array([[6, 8], [14, 16]], dtype=np.float32),
            "mode": np.array([[1, 3], [9, 11]], dtype=np.float32),
            "conservative_mean": np.array([[3.5, 5.5], [11.5, 13.5]], dtype=np.float32),
        }
        for method, wanted in expected.items():
            with self.subTest(method=method):
                result = lsdtt3.halve_raster(raster, method)
                np.testing.assert_allclose(result.to_numpy(), wanted)

        cubic = lsdtt3.halve_raster(raster, "cubic")
        self.assertEqual(cubic.shape, (2, 2))
        self.assertTrue(cubic.mask_numpy().all())
        self.assertEqual(lsdtt3.double_raster(raster, "nearest").shape, (8, 8))

        arbitrary = lsdtt3.RasterMeta.from_shape(3, 3, cell_size=4.0 / 3.0)
        self.assertEqual(lsdtt3.resample_raster(raster, arbitrary, "linear").shape, (3, 3))

    def test_power_of_two_helpers(self):
        metadata = lsdtt3.RasterMeta.from_shape(3, 3, cell_size=2.0)
        values = np.arange(9, dtype=np.float32).reshape(3, 3)
        raster = lsdtt3.Raster.from_numpy(values, np.ones((3, 3), bool), metadata)

        padded = lsdtt3.pad_to_nearest_power_of_two(raster)
        self.assertEqual(padded.shape, (4, 4))
        self.assertFalse(padded.mask_numpy()[0].any())
        np.testing.assert_array_equal(padded.to_numpy()[1:, :3], values)

        self.assertEqual(
            lsdtt3.interpolate_to_nearest_power_of_two(
                raster, direction="lower", method="linear"
            ).shape,
            (2, 2),
        )
        self.assertEqual(
            lsdtt3.interpolate_to_nearest_power_of_two(
                raster, direction="higher", method="linear"
            ).shape,
            (4, 4),
        )


if __name__ == "__main__":
    unittest.main()
