# Image-operation fixtures

Frozen behaviour of OpenCV 4.13 for the operations in `src/cpu/lsdtt_raster_image_ops_cpu`, captured on
2026-09-14 (`capture_info.txt`) before OpenCV was removed as a dependency. The capture tool was deleted with
the dependency, so none of the OpenCV values here can be regenerated. Read by `unit_image_ops`
(`tests/cpp/test_image_ops.cpp`).

The fixtures hold hashes, not images. Until 2026-10-08 the raw OpenCV outputs were kept as well (`*_full.bin`,
`*_nodata.bin` and `hough_opencv_segments.csv`, 2.9 MB). Every check now runs on hashes. The raw files are
in git history if a failing check ever needs a pixel-level diff.

| File | Contents | Used by |
|---|---|---|
| `manifest.csv` | One row per (input, operation, parameters). `reference_hash` is FNV-1a 64 over OpenCV's output, encoded as one byte per pixel (`u8`) or as 0/1 bits packed eight to a byte, least significant bit first (`bits`). | O2: native output must hash the same on inputs without nodata |
| `nodata_native.csv` | Pinned FNV-1a hash of the NATIVE output for every record on an input with nodata, same encoding. Written by `lsdtt-test-image-ops <dir> --record-nodata`. | O3 |
| `hough_native.csv` | Pinned hash of the native Hough output for the 32 Hough records. Written with `--record-hough`. | H3 |
| `hough_opencv_segment_hashes.csv` | Per case, OpenCV's `HoughLinesP` segment count and an FNV-1a 64 hash of its segments in order, each as x0, y0, x1, y1 little-endian int32. | H5 |
| `primitives.csv` | Hashes of OpenCV's 8-bit / 16-bit primitive outputs on a 512x487 random image. | P |
| `planted_segments.csv` | The segments planted in the `lines` input. | H1 |
| `capture_info.txt` | OpenCV version, IPP state and CPU of the capture. | - |
| `e2e_hashes.csv` | Output hashes for `tests/scripts/capture_image_ops_e2e_hashes.py`. | that script |

Re-record `nodata_native.csv` or `hough_native.csv` only for a deliberate, reviewed change to the native
output. The OpenCV hashes are never re-recorded.

## Measurements at capture (OpenCV vs native)

Pixel-level comparisons that needed the raw OpenCV images. They were measured by `unit_image_ops` on
2026-10-08, the last run that had them, and cannot be re-measured. The test still reports how many records
differ, from the hashes.

Inputs with nodata. The native code normalises over valid pixels only, while OpenCV used every pixel, so they
differ on purpose. Pixel counts include nodata pixels for `u8`.

| Operation | Records that differ | Pixels that differ | Valid pixels in total |
|---|---|---|---|
| u8 | 7 of 8 | 64078 | 145058 |
| otsu | 6 of 8 | 36324 | 145058 |
| adaptive_mean | 118 of 144 | 109236 | 2611044 |
| adaptive_gaussian | 112 of 144 | 104468 | 2611044 |
| canny | 82 of 96 | 67792 | 1740696 |

Hough on inputs without nodata. The native Hough uses its own random-number generator, so its output is not
OpenCV's. H5 shows that the port is otherwise exact.

| Input | threshold | Line pixels, OpenCV | Line pixels, native | OpenCV px within 1 px of native | native px within 1 px of OpenCV |
|---|---|---|---|---|---|
| ramp_full | 50 | 0 | 0 | 1.000 | 1.000 |
| ramp_full | 80 | 0 | 0 | 1.000 | 1.000 |
| peaks_full | 50 | 0 | 0 | 1.000 | 1.000 |
| peaks_full | 80 | 0 | 0 | 1.000 | 1.000 |
| noise_full | 50 | 5419 | 5266 | 0.946 | 0.983 |
| noise_full | 80 | 1391 | 1365 | 0.922 | 0.934 |
| steps_full | 50 | 0 | 0 | 1.000 | 1.000 |
| steps_full | 80 | 1298 | 1299 | 1.000 | 1.000 |
| signed_full | 50 | 316 | 278 | 0.759 | 0.867 |
| signed_full | 80 | 911 | 806 | 0.780 | 0.893 |
| lines_full | 50 | 1821 | 2159 | 1.000 | 0.986 |
| lines_full | 80 | 1786 | 1763 | 0.994 | 1.000 |
| entropy_full | 50 | 0 | 0 | 1.000 | 1.000 |
| entropy_full | 80 | 71 | 0 | 0.000 | 1.000 |
