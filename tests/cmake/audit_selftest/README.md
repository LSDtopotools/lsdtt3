# Audit self-test fixture

Input for `unit_test_inventory_audit_selftest`. `canned_CTestTestfile.cmake.in`
is a hand-written stand-in for a generated CTest file (not named
`CTestTestfile.cmake`, which `.gitignore` ignores everywhere), and `scripts/` for
`tests/scripts/`. Each check in `../audit_test_inventory.cmake` is violated on
purpose, and `expected.txt` lists exactly those violations. The self-test passes
only if the audit reports that set and nothing else, so the audit cannot pass
vacuously. Tests named `sel_ok_*` (and the benchmark) must NOT be reported.

None of this is run: the paths are never executed, only read.
