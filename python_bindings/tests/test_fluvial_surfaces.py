import unittest

import numpy as np

import lsdtt3


def raster(values, cell_size=1.0, valid=None):
    values = np.asarray(values, dtype=np.float32)
    meta = lsdtt3.RasterMeta.from_shape(values.shape[0], values.shape[1], cell_size=cell_size)
    if valid is None:
        valid = np.ones(values.shape, dtype=bool)
    return lsdtt3.Raster.from_numpy(values, valid, meta)


def flat_floored_valley(n=60, cell_size=30.0, half_floor=5, floor_gradient=0.01,
                        wall_gradient=0.3):
    """A valley draining to the last row, floor 2 * half_floor + 1 cells wide."""
    rows, cols = np.indices((n, n), dtype=np.float64)
    centre = n // 2
    across = np.maximum(np.abs(cols - centre) - half_floor, 0.0) * cell_size
    down = (n - 1 - rows) * cell_size
    return (100.0 + floor_gradient * down + wall_gradient * across).astype(np.float32)


class ThresholdAndLevelTests(unittest.TestCase):
    def test_select_threshold_fixed_and_qq(self):
        fixed = lsdtt3.select_threshold([1.0, 2.0, 3.0], "fixed", fixed_value=2.5)
        self.assertEqual(fixed["threshold"], 2.5)
        self.assertTrue(fixed["found"])
        self.assertEqual(fixed["n_values"], 3)

        rng = np.random.default_rng(3)
        sample = np.concatenate([rng.normal(10.0, 1.0, 5000), rng.uniform(0.0, 3.0, 800)])
        picked = lsdtt3.select_threshold(sample, "qq_gaussian")
        self.assertIn("threshold", picked)
        if picked["found"]:
            self.assertGreater(picked["threshold"], 0.0)
            self.assertLess(picked["threshold"], 10.0)

        with self.assertRaises(ValueError):
            lsdtt3.select_threshold(sample, "no_such_method")

    def test_select_threshold_reads_a_raster_ignoring_nodata(self):
        values = np.array([[1.0, 2.0], [3.0, 99.0]])
        valid = np.array([[True, True], [True, False]])
        result = lsdtt3.select_threshold(raster(values, valid=valid), "fixed", fixed_value=0.0)
        self.assertEqual(result["n_values"], 3)

    def test_cluster_surface_levels_splits_at_the_gap(self):
        result = lsdtt3.cluster_surface_levels(
            [1.0, 1.1, 1.2, 5.0, 5.1, 9.0, 9.2], required_gap_metres=1.0, n_null=99)
        clusters = result["clusters"]
        self.assertEqual(len(clusters), 3)
        self.assertEqual(list(clusters["n_values"]), [3, 2, 2])
        self.assertAlmostEqual(clusters["centre_metres"].iloc[1], 5.05)
        self.assertGreater(result["empty_fraction"], 0.5)
        self.assertTrue(0.0 < result["p_value"] <= 1.0)
        self.assertEqual(result["n_null"], 99)


class FloodplainThresholdTests(unittest.TestCase):
    def test_floodplain_must_touch_the_channel(self):
        hand = np.full((9, 9), 10.0)
        hand[:, 3:6] = 0.5            # low ground along the channel
        hand[0:2, 7:9] = 0.5          # low, but not connected to it
        channel = np.zeros((9, 9))
        channel[:, 4] = 1.0
        result = lsdtt3.classify_floodplain_terraces_by_thresholds(
            raster(np.ones((9, 9))), raster(hand), raster(channel),
            max_floodplain_hand_metres=1.0)
        fp = result.floodplain_mask.to_numpy()
        self.assertTrue(np.all(fp[:, 3:6] == 1.0))
        self.assertTrue(np.all(fp[0:2, 7:9] == 0.0))
        self.assertEqual(result.n_floodplain_dropped_no_channel, 1)
        self.assertEqual(result.n_floodplain_pixels, 27)

    def test_terraces_start_above_the_floodplain_bound(self):
        hand = np.full((9, 9), 20.0)
        hand[:, 3:6] = 0.5
        hand[:, 0:2] = 4.0
        channel = np.zeros((9, 9))
        channel[:, 4] = 1.0
        result = lsdtt3.classify_floodplain_terraces_by_thresholds(
            raster(np.ones((9, 9))), raster(hand), raster(channel),
            terrace_field=raster(hand), max_floodplain_hand_metres=1.0,
            find_terraces=True, min_terrace_height_above_floodplain_metres=1.0,
            max_terrace_height_metres=10.0)
        band = result.band_index.to_numpy()
        self.assertTrue(np.all(band[:, 0:2] == 1.0))
        self.assertTrue(np.all(band[:, 3:6] == 0.0))
        self.assertTrue(np.all(band[:, 6:9] == -1.0))

    def test_gradient_bound_needs_a_gradient(self):
        ones = raster(np.ones((5, 5)))
        with self.assertRaisesRegex(ValueError, "gradient"):
            lsdtt3.classify_floodplain_terraces_by_thresholds(
                ones, ones, ones, max_floodplain_hand_metres=1.0,
                max_floodplain_gradient=0.1)


class ComparisonTests(unittest.TestCase):
    def test_identical_methods_agree_and_bits_mark_each_method(self):
        a = np.array([[0, 1, 1], [2, 2, 0], [1, 0, 2]], dtype=np.float32)
        b = a.copy()
        b[0, 0] = 1.0
        classes = {"ladder": raster(a), "clubb2017": raster(b)}
        table = lsdtt3.compare_surface_classes(classes)
        self.assertEqual(len(table), 1)
        self.assertEqual(table.loc[0, "method_a"], "ladder")
        self.assertAlmostEqual(table.loc[0, "terrace_iou"], 1.0)
        self.assertAlmostEqual(table.loc[0, "floodplain_iou"], 3.0 / 4.0)

        bits = lsdtt3.surface_class_agreement(classes, "floodplain").to_numpy()
        self.assertEqual(bits[0, 0], 2.0)   # only clubb2017
        self.assertEqual(bits[0, 1], 3.0)   # both
        self.assertEqual(bits[1, 0], 0.0)

    def test_class_values_are_checked(self):
        bad = raster(np.full((3, 3), 4.0))
        with self.assertRaises(ValueError):
            lsdtt3.compare_surface_classes({"a": bad, "b": bad})

    def test_baseline_relief_on_a_plane(self):
        dem = np.tile(np.arange(10, dtype=np.float32) * 2.0, (10, 1))
        baseline = np.zeros((10, 10))
        baseline[:, 0] = 1.0
        relief = lsdtt3.relief_above_nearest_baseline(
            raster(dem, cell_size=10.0), raster(baseline, cell_size=10.0)).to_numpy()
        np.testing.assert_allclose(relief[:, 4], 8.0)


class RiserTests(unittest.TestCase):
    def test_a_step_between_tilted_treads_is_a_riser(self):
        n = 80
        rows, cols = np.indices((n, n), dtype=np.float64)
        dem = 0.02 * cols * 2.0 + np.where(cols >= 40, 5.0, 0.0)
        result = lsdtt3.find_risers(
            raster(dem, cell_size=2.0), polyfit_lambda_metres=10.0,
            vec_adaptive_block_size_pixels=31, min_step_height_metres=1.0)
        segments = result.segments
        self.assertGreater(len(segments), 0)
        self.assertGreater(segments["step_height_metres"].max(), 2.5)
        self.assertEqual(result.centreline.shape, (n, n))
        self.assertTrue({"segment_id", "length_metres", "face_gradient", "geometry"}
                        <= set(segments.columns))

    def test_bad_gap_closing_choice_is_refused(self):
        with self.assertRaises(ValueError):
            lsdtt3.find_risers(raster(np.zeros((20, 20))), gap_closing_algorithm_choice="nope")


class ScopeBasedTests(unittest.TestCase):
    def setUp(self):
        self.cell = 30.0
        self.dem = raster(flat_floored_valley(cell_size=self.cell), cell_size=self.cell)
        self.filled = lsdtt3.fill_pits(self.dem)
        self.flow = lsdtt3.FlowInfo(self.filled, boundary_conditions="oooo")

    def test_stream_order_definition_keeps_fewer_pixels(self):
        every = self.flow.channel_scope(threshold_contributing_pixels=10)
        ordered = self.flow.channel_scope(
            threshold_contributing_pixels=10,
            channel_definition_algorithm_choice="stream_order", min_channel_stream_order=2)
        n_every = int(every.mask.to_numpy().sum())
        n_ordered = int(ordered.mask.to_numpy().sum())
        self.assertGreater(n_every, 0)
        self.assertLess(n_ordered, n_every)
        with self.assertRaises(ValueError):
            self.flow.channel_scope(threshold_contributing_pixels=10,
                                    channel_definition_algorithm_choice="nope")

    def test_perimeter_source_where_the_floor_enters(self):
        # A 25-cell (750 m) floor: wide enough that the 300 m minimum survives
        # the gradient fit's reach into the walls.
        dem = raster(flat_floored_valley(cell_size=self.cell, half_floor=12), cell_size=self.cell)
        filled = lsdtt3.fill_pits(dem)
        flow = lsdtt3.FlowInfo(filled, boundary_conditions="oooo")
        result = lsdtt3.find_perimeter_sources(dem, filled, flow)
        sources = result.sources
        self.assertEqual(len(sources), 1)
        self.assertEqual(int(sources["row"].iloc[0]), 0)
        self.assertEqual(sorted(result.crossings["verdict"]), ["outlet_crossing", "source"])
        self.assertEqual(result.pixels.shape, (len(sources), 2))
        self.assertTrue(result.outlets["selected"].any())
        self.assertGreater(result.report["profile_length_metres"], 0.0)
        self.assertIn("verdict", result.crossings.columns)

    def test_ladder_keywords_and_tables(self):
        scope = self.flow.channel_scope(threshold_contributing_pixels=10)
        labels = raster(np.where(np.abs(np.indices((60, 60))[1] - 30) <= 5, 1.0, 0.0),
                        cell_size=self.cell)
        result = lsdtt3.classify_flat_surfaces_with_scope(
            labels, self.filled, scope, max_floodplain_hand_metres=2.0,
            ladder_break_algorithm_choice="jenks", ladder_min_window_pixels=10,
            ladder_unresolved_window_algorithm_choice="nearest_resolved",
            n_threads=1)
        self.assertEqual(result.surface_class.shape, (60, 60))
        self.assertEqual(result.banded_labels.shape, (60, 60))
        self.assertTrue({"source_id", "status", "breaks_metres"} <= set(result.ladder_windows.columns))
        self.assertTrue({"label", "band", "median_field"} <= set(result.patches.columns))
        self.assertEqual(len(result.training_polygons), 0)
        with self.assertRaises(ValueError):
            lsdtt3.classify_flat_surfaces_with_scope(
                labels, self.filled, scope, ladder_gap_rank="nope")
        with self.assertRaises(ValueError):
            lsdtt3.classify_flat_surfaces_with_scope(
                labels, self.filled, scope, riser_mask=labels)

    def test_legacy_classify_keywords_warn(self):
        scope = self.flow.channel_scope(threshold_contributing_pixels=10)
        labels = raster(np.zeros((60, 60)), cell_size=self.cell)
        with self.assertWarns(DeprecationWarning):
            lsdtt3.classify_flat_surfaces_with_scope(
                labels, self.filled, scope, floodplain_max_height_metres=2.0)

    def test_width_network_relative_elevation_and_profile(self):
        scope = self.flow.channel_scope(threshold_contributing_pixels=10)
        floor = raster(np.where(np.abs(np.indices((60, 60))[1] - 30) <= 8, 1.0, 0.0),
                       cell_size=self.cell)
        network = lsdtt3.build_valley_width_network(self.filled, scope, floor)
        self.assertGreater(len(network.segments), 0)
        centrelines = network.centrelines
        rem, owner = centrelines.relative_elevation(self.filled)
        self.assertEqual(rem.shape, (60, 60))
        on_floor = floor.to_numpy() > 0.5
        self.assertGreater(rem.mask_numpy()[on_floor].mean(), 0.9)
        self.assertFalse(rem.mask_numpy()[~on_floor].any())
        patches = raster(np.where(np.abs(np.indices((60, 60))[1] - 30) <= 2, 1.0, 0.0),
                         cell_size=self.cell)
        bins, table = centrelines.terrace_profile(self.filled, patches,
                                                  terrace_profile_bin_metres=300.0)
        self.assertGreater(len(bins), 0)
        self.assertIn("relative_elevation_median", bins.columns)
        self.assertEqual(list(table["patch_label"]), [1])

        from_scope = lsdtt3.extract_valley_centrelines(
            self.filled, scope, valley_mask=floor,
            centreline_baseline_algorithm_choice="channel_crossings")
        rem2, _ = from_scope.relative_elevation(self.filled, "whole_dem")
        self.assertEqual(rem2.shape, (60, 60))

    def test_relative_elevation_needs_flow_ranked_centrelines(self):
        channel = np.column_stack((np.arange(2, 58), np.full(56, 30)))
        floor = raster(np.where(np.abs(np.indices((60, 60))[1] - 30) <= 8, 1.0, 0.0),
                       cell_size=self.cell)
        centrelines = lsdtt3.extract_valley_centrelines(
            self.filled, [channel], valley_mask=floor, centreline_loops=1)
        with self.assertRaisesRegex(ValueError, "drainage area"):
            centrelines.relative_elevation(self.filled)


class RidgelineMovingAverageTests(unittest.TestCase):
    def test_moving_average_fits_are_kept_per_window(self):
        values = np.tile(np.linspace(0.0, 20.0, 30, dtype=np.float32), (5, 1))
        values[2, ::4] += 1.0
        dem = raster(values, cell_size=10.0)
        pixels = np.column_stack((np.full(30, 2), np.arange(30)))
        ridge = lsdtt3.Ridgeline(dem, pixels)
        ridge.fill_elevation(dem)
        ridge.fit_profile(detrend_moving_average_window_pixels=[3, 5])
        fits = ridge.moving_average_fits
        self.assertEqual(sorted(fits), [3, 5])
        fitted, detrended = fits[3]
        self.assertEqual(len(fitted), 30)
        self.assertAlmostEqual(fitted[10] + detrended[10], ridge.elevation[10], places=4)
        with self.assertRaises(ValueError):
            ridge.fit_profile(detrend_moving_average_window_pixels=[4])


if __name__ == "__main__":
    unittest.main()
