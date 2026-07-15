# lsdtt3 smoke tests

This directory contains a lightweight CLI smoke-test suite wired into CTest.

## What it tests

- `smoke_raster_preprocessing`: runs `lsdtt-raster-preprocessing` on a tiny DEM and checks `tiny_fill.tif` exists and is non-empty.
- `smoke_raster_preprocessing_tif`: runs `lsdtt-raster-preprocessing` on `tiny30_dem.tif` and checks `tiny30_fill.tif` exists and is non-empty.
- `smoke_channel_extraction`: runs `lsdtt-channel-extraction` on the same DEM and checks `tiny_channel_network_lines.fgb` exists and is non-empty.
- `smoke_basin_extraction`: runs `lsdtt-basin-extraction` and checks `tiny_junction_basin_outlines.fgb` exists and is non-empty.
- `smoke_nodata_fill`: runs `lsdtt-raster-preprocessing` on a tiny DEM with an interior nodata pixel and checks `tiny_nodata_nodata_filled.tif` exists and reports `STATISTICS_VALID_PERCENT=100` via `gdalinfo -stats`.
- `smoke_nodata_edge_unfilled`: runs `lsdtt-raster-preprocessing` on a tiny DEM with an edge nodata pixel and checks `STATISTICS_VALID_PERCENT=96` (edge nodata remains unfilled).
- `smoke_drainage_analysis`: runs `lsdtt-drainage-analysis` and checks `tiny_d8_drainage_area.tif` exists.
- `smoke_chi_analysis`: runs `lsdtt-chi-analysis` and checks `tiny_chi_fit_nodes_mn_0p500.csv` exists.
- `smoke_raster_calculator_difference`: runs `lsdtt-raster-calculator` and checks `tiny_difference.tif` exists.

More rigorous checks currently included:

- `rigorous_channel_lines_feature_count`: requires channel lines output to contain at least 1 feature.
- `rigorous_basin_outlines_feature_count`: requires basin outlines output to contain at least 1 feature.
- `rigorous_tif_fill_raster_size`: requires the `tiny30` fill raster to be exactly `32x39` pixels.
- `rigorous_raster_calculator_difference_checksum`: checks the GDAL checksum of `tiny_difference.tif`.
- `rigorous_drainage_d8_checksum`: checks the GDAL checksum of `tiny_d8_drainage_area.tif`.
- `rigorous_chi_nodes_sha256`: checks SHA256 of `tiny_chi_fit_nodes_mn_0p500.csv`.
- `failure_channel_extraction_invalid_threshold`: expected-failure test for invalid `source_pixels_threshold`.
- `failure_raster_calculator_missing_second`: expected-failure test for missing `second_raster_fname`.

## Run locally

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Run by label:

```bash
ctest --test-dir build -L smoke --output-on-failure
ctest --test-dir build -L rigorous --output-on-failure
ctest --test-dir build -L failure --output-on-failure
```

Label conventions:

- `smoke_*` tests are labeled `smoke;cli`
- `rigorous_*` tests are labeled `rigorous;cli`
- `failure_*` tests are labeled `failure;cli`

## Add a new smoke test

1. Add a param template in `tests/params/`.
2. Add expected output filename(s).
3. Register it in `tests/CMakeLists.txt` with `add_lsdtt_cli_smoke_test(...)`.

The runner script is `tests/scripts/run_cli_smoke.py`.

## Fixtures

- `tests/data/tiny30_dem.tif` is a 30x30 GeoTIFF fixture clipped from `../data/scotland/scotland.tif` for future larger smoke/regression tests.
