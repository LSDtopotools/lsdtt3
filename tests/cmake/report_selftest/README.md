# Test-report self-test fixture

Input for `unit_test_report_selftest`. `canned/` holds a miniature build tree,
stored FLAT with a `.in` suffix: its real names (`build/`, `Testing/`,
`CTestTestfile.cmake`) are all in `.gitignore`, and under them the fixture
was never committed. `make_test_report.cmake` lays the files out in a scratch
directory before building the report. The tree has
five tests (`sel_pass`, `sel_fail` -- a Python-harness test with its own work
directory and logs -- `sel_timeout`, `sel_skip`, and `sel_running`, which the
console log shows starting with no result, as when a run is interrupted), the console log that
`run_tests.cmake` would have written, `last_run.txt`, CTest's `LastTest.log`,
and a build-info file. `@HOME@` and `@BUILD@` are replaced with the real home
directory and the scratch copy's path before the report is built, so the
self-test can check that the home directory is masked.

Nothing here is executed; the paths only have to look real.
