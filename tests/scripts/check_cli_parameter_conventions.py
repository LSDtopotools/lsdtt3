#!/usr/bin/env python3
"""CLI parameter-naming conformance check (Phase 0d).

Implements item **0d** of `docs/plan_docs/cli_consistency_remediation_plan.md`:
a CTest case that walks every built `lsdtt-*` CLI, reads its declared
parameters out of the help registry (`<cli> -h json`), and asserts the
mechanical naming rules of plan section 3 ("The standard to enforce"):

  R1 unit_suffix        no parameter name ends in `_m`, `_px` or `_deg`;
                        the spellings are `_metres`, `_pixels`, `_degrees`.
                        `<what>_area_threshold_m2` is endorsed by section 3
                        and is NOT matched (the name does not end in `_m`).
  R2 print_writes_file  no `print_*` parameter that writes a file; output
                        switches are `write_<output>`.
  R3 bool_raster_suffix a bool parameter that writes a raster is named
                        `write_<name>_raster` (`_rasters` for a set).
  R4 duplicate_format   at most one raster-format parameter
                        (`write_output_format`) and at most one vector-format
                        parameter (`write_vector_format`) per CLI.
  R5 minimum_prefix     bounds use `min_`/`max_`, not `minimum_`.  The family
                        `<subject>_minimum_contributing_pixels` is a
                        DELIBERATE exception recorded in sections 3 and 10b
                        (basin_, crossing_, side_tributary_) and is not
                        flagged.
  R6 tol_abbreviation   no `tol` abbreviation; metric tolerances are
                        `_tolerance_metres`, snap radii are
                        `<what>_snap_max_search_radius_pixels`.  Only a whole
                        underscore-token `tol` (or a `tol...` token that is not
                        `tolerance...`) is matched, so `tolerance` and words
                        that merely contain the letters are safe.

A CLI whose `-h json` cannot be parsed is itself reported, as rule
`help_json_unparseable`, through the same baseline mechanism.

REPORT-ONLY WITH A BASELINE ALLOW-LIST
--------------------------------------
Today's tree violates these rules in many places; the phases that fix them come
later.  So every violation present when the baseline was written is recorded in
the allow-list and does not fail the test.  The test fails only on a violation
that is NOT in the baseline -- i.e. a newly introduced one.

A baseline entry that no longer matches any violation ("stale") is REPORTED but
does NOT fail the test, so that a rename landing in a later phase is never
blocked by the allow-list; the operator simply deletes the stale lines.  The
allow-list is expected to shrink monotonically and never to grow.

Entries are keyed on (cli, parameter, rule), which is tight enough that
renaming a parameter always makes its entry stale rather than silently
transferring the exemption to the new name.

Usage:
    check_cli_parameter_conventions.py --bin-dir <dir> --baseline <json>
    check_cli_parameter_conventions.py --bin-dir <dir> --baseline <json> \
        --write-baseline        # regenerate the allow-list wholesale

Regenerating wholesale DEFEATS THE PURPOSE of the guard: it launders every
violation introduced since the last refresh into the allow-list.  Normally the
file should only ever be edited by DELETING lines.  `--write-baseline` exists
for the one legitimate case: the first creation of the file.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# Rule identifiers.  These strings are part of the baseline key: renaming one
# invalidates every baseline entry that uses it, which is deliberate.
# --------------------------------------------------------------------------
R_UNIT = "unit_suffix"
R_PRINT = "print_writes_file"
R_BOOL_RASTER = "bool_raster_suffix"
R_FORMAT = "duplicate_format"
R_MINIMUM = "minimum_prefix"
R_TOL = "tol_abbreviation"
R_UNPARSEABLE = "help_json_unparseable"

# Section 3 / 10b deliberate exception to the `min_` rule.
MINIMUM_EXCEPTION = re.compile(r"^.+_minimum_contributing_pixels$")

BAD_UNIT_SUFFIX = re.compile(r"_(m|px|deg)$")

CANONICAL_RASTER_FORMAT = "write_output_format"
CANONICAL_VECTOR_FORMAT = "write_vector_format"

# Reviewed Phase 5 exceptions: these selectors govern a genuinely different
# output domain or stage and inherit the canonical selector unless documented
# otherwise. Any new override must be reviewed and added explicitly.
REVIEWED_FORMAT_OVERRIDES = {
    "lsdtt-basin-extraction": {
        "basin_outlines_write_vector_format",
        "side_tributary_write_vector_format",
    },
    "lsdtt-channel-extraction": {
        "burned_stream_profiles_write_vector_format",
        "channel_network_junctions_write_vector_format",
        "channel_network_points_write_vector_format",
        "flow_lines_write_vector_format",
        "wiener_diagnostics_write_output_format",
    },
    "lsdtt-channel-network-analysis": {"qa_write_vector_format"},
    "lsdtt-chi-analysis": {"chi_data_maps_write_vector_format"},
    "lsdtt-cross-sections": {"cross_section_sample_points_write_vector_format"},
    "lsdtt-hillslope-channel-coupling": {"hilltop_points_write_vector_format"},
    "lsdtt-marsh-topography": {"marsh_points_write_vector_format"},
    "lsdtt-stream-burning": {
        "burned_profiles_write_vector_format",
        "junction_comparison_write_vector_format",
        "vector_sources_write_vector_format",
    },
}

# Legacy input selectors whose names predate the `_input_format` suffix. They
# are not output-selector duplicates and are outside the Phase 5 migration.
REVIEWED_INPUT_FORMATS = {"sink_candidate_input_format"}

# Heuristic vocabularies.  Deliberately narrow: a false positive here is a
# test that fails on correct code, which is worse than a missed violation that
# the later phases will catch by hand anyway.
WRITE_VERB = re.compile(
    r"\b(write|writes|writing|written|save|saves|saved|export|exports|exported)\b"
)
FILE_EXTENSION = re.compile(r"\.(csv|tif|tiff|txt|json|geojson|gpkg|fgb|png|asc|bil)\b")
RASTER_WORD = re.compile(r"\b(raster|rasters|geotiff|cog)\b|\.(tif|tiff|bil|asc)\b")
VECTOR_WORD = re.compile(
    r"\b(csv|gpkg|geopackage|flatgeobuf|fgb|shapefile|geojson|vector|point|points|line|lines|polygon|polygons|table)\b"
)


def discover_clis(bin_dir: Path, cli_src_dir: Path | None = None) -> list[Path]:
    """Executable regular files named lsdtt-* in bin_dir that are CURRENT CLIs.

    The build tree also holds strays such as
    `lsdtt-raster-preprocessing_help.csv`, so anything with a filename
    extension is skipped as well as anything non-executable.

    `lsdtt-test-*` binaries are skipped too: those are C++ unit-test drivers,
    not CLIs.  They have no parameter registry and they IGNORE their arguments,
    so `lsdtt-test-hillslope-implicit -h json` does not print help -- it runs
    the nine-minute implicit solver.

    STALE BINARIES ARE SKIPPED.  A build tree accumulates executables for CLIs
    that have since been deleted from CMakeLists.txt -- `lsdtt-landslides`,
    `lsdtt-hydraulic-routing` and `lsdtt-raster-basin-stats` were all present
    when this check was written, the oldest dating from April 2026.  Scanning
    them baselines violations for code that no longer exists, and the parameter
    registry of a binary built months ago says nothing about the current tree.
    So `cli/<name>.cpp` must exist for the binary to be scanned; the source
    directory is the authority on which CLIs there are, not the build tree.
    """
    found = []
    for entry in sorted(bin_dir.glob("lsdtt-*")):
        if not entry.is_file():
            continue
        if entry.suffix:  # .csv, .html, .exe is not expected on this platform
            continue
        if entry.name.startswith("lsdtt-test-"):
            continue
        if not os.access(entry, os.X_OK):
            continue
        if cli_src_dir is not None and not (cli_src_dir / (entry.name + ".cpp")).is_file():
            continue  # stale binary: no source in cli/
        found.append(entry)
    return found


def unmeasurable_clis(bin_dir: Path, cli_src_dir: Path) -> tuple[list[str], list[str]]:
    """CLI sources whose binary is missing, or older than the source.

    Returns (missing, stale).

    WHY THIS IS FATAL RATHER THAN A NOTE.  This check reads each CLI's parameter
    registry out of the BUILT BINARY via `-h json`.  A binary older than its
    source reports the parameters that CLI had when it was last compiled, not
    the ones it has now -- so a violation added since then is invisible and the
    check passes for the wrong reason.

    That is not hypothetical.  The Phase 0 baseline was generated against a tree
    in which 10 of 38 binaries were stale and one was absent, and it therefore
    recorded two violations that the source had already fixed while missing nine
    that were really there.  Both errors pointed the same way: the guard looked
    greener than the tree was.

    A guard that silently under-measures is worse than no guard, so an
    incomplete build fails here with the exact list of targets to build.
    """
    missing, stale = [], []
    for src in sorted(cli_src_dir.glob("lsdtt-*.cpp")):
        exe = bin_dir / src.stem
        if not exe.is_file():
            missing.append(src.stem)
        elif exe.stat().st_mtime < src.stat().st_mtime:
            stale.append(src.stem)
    return missing, stale


def read_parameters(exe: Path, timeout: float) -> list[dict]:
    """Run `<exe> -h json` and parse the JSON document that follows the banner.

    Raises RuntimeError with a short reason if the document cannot be read.
    """
    try:
        proc = subprocess.run(
            [str(exe), "-h", "json"],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError("timed out after %g s" % timeout)
    except OSError as exc:
        raise RuntimeError("could not execute: %s" % exc)

    text = proc.stdout
    brace = text.find("{")
    if brace < 0:
        raise RuntimeError(
            "no JSON document in stdout (exit code %d)" % proc.returncode
        )
    try:
        doc = json.loads(text[brace:])
    except json.JSONDecodeError as exc:
        raise RuntimeError("JSON parse error: %s" % exc)

    params = doc.get("parameters")
    if not isinstance(params, list):
        raise RuntimeError("document has no 'parameters' array")
    return params


def helptext(param: dict) -> str:
    return "{} {}".format(
        param.get("description", "") or "", param.get("guidance", "") or ""
    ).lower()


def looks_like_file_write(param: dict) -> bool:
    """R2 heuristic: does this parameter's help say it writes a file?

    Conservative: requires an explicit write/save/export verb, OR a literal
    file extension, in the description or guidance.  A `print_*` parameter
    whose help only says e.g. "report progress to the terminal" is not flagged.
    """
    text = helptext(param)
    return bool(WRITE_VERB.search(text) or FILE_EXTENSION.search(text))


def looks_like_raster_write(param: dict) -> bool:
    """R3 heuristic: is this a bool switch that writes a RASTER?

    Requires all of:
      * type is bool;
      * the name is an output switch (`write_`/`print_` prefix, or a `_write_`
        segment, which catches stage-prefixed forms like
        `wiener_write_filtered_dem`);
      * the help contains a raster word (raster/geotiff/cog or a raster file
        extension);
      * the help contains NO vector/table word (csv, gpkg, flatgeobuf, vector,
        points, lines, ...), so a switch that writes both, or writes a vector
        derived from a raster, is left alone rather than mis-flagged.
    """
    if param.get("type") != "bool":
        return False
    name = param.get("name", "")
    if not (name.startswith("write_") or name.startswith("print_") or "_write_" in name):
        return False
    text = helptext(param)
    if not RASTER_WORD.search(text):
        return False
    if VECTOR_WORD.search(text):
        return False
    return True


def classify_format_parameter(param: dict) -> str | None:
    """R4: is this a raster-format or a vector-format selector?

    Only names ending in `_format` are considered, and `_input_format` is
    excluded (section 3 makes input format a separate, per-input concept).
    A parameter is vector if its name mentions `vector`, or if its help
    mentions vector containers without mentioning raster containers; raster
    if its help mentions raster containers.  Anything undecidable returns
    None and is not counted.
    """
    name = param.get("name", "")
    if not name.endswith("_format"):
        return None
    if name.endswith("_input_format") or name in REVIEWED_INPUT_FORMATS:
        return None
    text = helptext(param)
    if "vector" in name:
        return "vector"
    if "raster" in name:
        return "raster"
    has_raster = bool(re.search(r"\b(cog|gtiff|geotiff|png|tif|tiff|bil|asc)\b", text))
    has_vector = bool(
        re.search(r"\b(gpkg|geopackage|flatgeobuf|fgb|shapefile|geojson)\b", text)
    )
    if has_vector and not has_raster:
        return "vector"
    if has_raster and not has_vector:
        return "raster"
    return None


def check_cli(cli_name: str, params: list[dict]) -> list[dict]:
    """Return every violation for one CLI, as {cli, parameter, rule, detail}."""
    violations = []

    def add(param_name, rule, detail):
        violations.append(
            {"cli": cli_name, "parameter": param_name, "rule": rule, "detail": detail}
        )

    formats = {"raster": [], "vector": []}

    for param in params:
        name = param.get("name", "")
        if not name:
            continue

        # R1
        if BAD_UNIT_SUFFIX.search(name):
            suffix = BAD_UNIT_SUFFIX.search(name).group(1)
            add(
                name,
                R_UNIT,
                "name ends in _%s; units are spelled out (_metres/_pixels/_degrees)"
                % suffix,
            )

        # R2
        if name.startswith("print_") and looks_like_file_write(param):
            add(
                name,
                R_PRINT,
                "print_* parameter whose help describes writing a file; "
                "output switches are write_<output>",
            )

        # R3
        if looks_like_raster_write(param) and not (
            name.endswith("_raster") or name.endswith("_rasters")
        ):
            add(
                name,
                R_BOOL_RASTER,
                "bool switch that writes a raster must be named "
                "write_<name>_raster (or _rasters for a set)",
            )

        # R5
        if "minimum_" in name and not MINIMUM_EXCEPTION.match(name):
            add(
                name,
                R_MINIMUM,
                "bounds use min_/max_, not minimum_ (the "
                "<subject>_minimum_contributing_pixels family is exempt)",
            )

        # R6
        for token in name.split("_"):
            if token.startswith("tol") and not token.startswith("tolerance"):
                add(
                    name,
                    R_TOL,
                    "'tol' abbreviation; metric tolerances are _tolerance_metres "
                    "and snap radii are <what>_snap_max_search_radius_pixels",
                )
                break

        # R4 (collected, judged below)
        kind = classify_format_parameter(param)
        if (kind is not None and
                name not in REVIEWED_FORMAT_OVERRIDES.get(cli_name, set())):
            formats[kind].append(name)

    # R4: at most one of each kind.  The canonical name keeps its place when
    # present; otherwise the first in sorted order is kept, so the choice of
    # which duplicates are flagged is deterministic across runs.
    for kind, canonical in (
        ("raster", CANONICAL_RASTER_FORMAT),
        ("vector", CANONICAL_VECTOR_FORMAT),
    ):
        names = sorted(set(formats[kind]))
        if len(names) <= 1:
            continue
        keep = canonical if canonical in names else names[0]
        for name in names:
            if name == keep:
                continue
            add(
                name,
                R_FORMAT,
                "%d %s-format parameters in this CLI (%s); there should be one, "
                "named %s" % (len(names), kind, ", ".join(names), canonical),
            )

    return violations


def key_of(violation: dict) -> tuple[str, str, str]:
    return (violation["cli"], violation["parameter"], violation["rule"])


BASELINE_README = [
    "Baseline allow-list for tests/scripts/check_cli_parameter_conventions.py",
    "(Phase 0d of docs/plan_docs/cli_consistency_remediation_plan.md).",
    "",
    "Every entry is a CLI parameter that violates a section-3 naming",
    "convention TODAY and is knowingly tolerated until the phase that renames",
    "it lands.  The test passes on these and fails only on violations that are",
    "NOT listed here.",
    "",
    "THIS LIST MUST ONLY EVER SHRINK.  Edit it by DELETING lines.  Do not",
    "regenerate it wholesale with --write-baseline: doing so launders every",
    "violation introduced since the last refresh into the allow-list and the",
    "guard stops guarding.  A stale entry (one whose violation no longer",
    "occurs) is reported by the test but does not fail it, so deleting stale",
    "lines is never urgent and never blocks a phase.",
]


def load_baseline(path: Path) -> set[tuple[str, str, str]]:
    with path.open() as handle:
        doc = json.load(handle)
    entries = doc.get("violations", [])
    return {(e["cli"], e["parameter"], e["rule"]) for e in entries}


def write_baseline(path: Path, violations: list[dict]) -> None:
    doc = {
        "README": BASELINE_README,
        "violations": [
            {"cli": v["cli"], "parameter": v["parameter"], "rule": v["rule"]}
            for v in sorted(violations, key=key_of)
        ],
    }
    with path.open("w") as handle:
        json.dump(doc, handle, indent=2)
        handle.write("\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--bin-dir",
        required=True,
        type=Path,
        help="directory holding the built lsdtt-* CLI executables",
    )
    ap.add_argument(
        "--baseline",
        required=True,
        type=Path,
        help="JSON allow-list of currently tolerated violations",
    )
    ap.add_argument(
        "--write-baseline",
        action="store_true",
        help="regenerate the allow-list from today's tree (see the module "
        "docstring: this normally defeats the purpose of the guard)",
    )
    ap.add_argument(
        "--cli-src-dir",
        type=Path,
        default=None,
        help="the repository's cli/ directory. When given, it is the authority "
        "on which CLIs exist: a built lsdtt-* executable with no matching "
        "cli/<name>.cpp is a STALE binary from a deleted CLI and is skipped, "
        "and any source with no binary is reported as unchecked.",
    )
    ap.add_argument("--timeout", type=float, default=60.0)
    args = ap.parse_args()

    if not args.bin_dir.is_dir():
        print("ERROR: --bin-dir does not exist: %s" % args.bin_dir)
        return 2
    if args.cli_src_dir is not None and not args.cli_src_dir.is_dir():
        print("ERROR: --cli-src-dir does not exist: %s" % args.cli_src_dir)
        return 2

    clis = discover_clis(args.bin_dir, args.cli_src_dir)
    if not clis:
        print("ERROR: no lsdtt-* executables found in %s" % args.bin_dir)
        return 2

    if args.cli_src_dir is not None:
        missing, stale = unmeasurable_clis(args.bin_dir, args.cli_src_dir)
        if missing or stale:
            print("ERROR: the build tree does not match cli/, so this check "
                  "cannot measure what is actually in the source.")
            if missing:
                print("  %d CLI(s) not built at all:" % len(missing))
                for name in missing:
                    print("      %s" % name)
            if stale:
                print("  %d binary/binaries older than their source:" % len(stale))
                for name in stale:
                    print("      %s" % name)
            print("  Build them and re-run, e.g.:")
            print("      cmake --build %s -j" % args.bin_dir)
            print("  See unmeasurable_clis() for why this is fatal rather than "
                  "a warning.")
            return 2

    violations: list[dict] = []
    n_parsed = 0
    n_params = 0
    for exe in clis:
        try:
            params = read_parameters(exe, args.timeout)
        except RuntimeError as exc:
            violations.append(
                {
                    "cli": exe.name,
                    "parameter": "-",
                    "rule": R_UNPARSEABLE,
                    "detail": "`-h json` could not be read: %s" % exc,
                }
            )
            continue
        n_parsed += 1
        n_params += len(params)
        violations.extend(check_cli(exe.name, params))

    print("CLI parameter-convention check (plan 0d)")
    print("  binaries found        : %d" % len(clis))
    print("  help registries parsed: %d" % n_parsed)
    print("  parameters inspected  : %d" % n_params)

    by_rule: dict[str, int] = {}
    for v in violations:
        by_rule[v["rule"]] = by_rule.get(v["rule"], 0) + 1
    print("  violations found      : %d" % len(violations))
    for rule in sorted(by_rule):
        print("      %-22s %d" % (rule, by_rule[rule]))

    if args.write_baseline:
        write_baseline(args.baseline, violations)
        print("\nWrote baseline with %d entries to %s" % (len(violations), args.baseline))
        return 0

    if not args.baseline.is_file():
        print("\nERROR: baseline file not found: %s" % args.baseline)
        print("Create it once with --write-baseline.")
        return 2

    baseline = load_baseline(args.baseline)
    seen = {key_of(v) for v in violations}

    new = sorted((v for v in violations if key_of(v) not in baseline), key=key_of)
    stale = sorted(entry for entry in baseline if entry not in seen)

    if stale:
        print(
            "\n%d stale baseline entr%s (no longer matching any violation) "
            "can be removed from %s:"
            % (len(stale), "y" if len(stale) == 1 else "ies", args.baseline)
        )
        for cli, param, rule in stale:
            print("  STALE  %s  %s  %s" % (cli, param, rule))

    if new:
        print(
            "\nFAIL: %d new parameter-convention violation%s "
            "(not in the baseline allow-list):"
            % (len(new), "" if len(new) == 1 else "s")
        )
        for v in new:
            print(
                "  %s  %s  [%s]  %s"
                % (v["cli"], v["parameter"], v["rule"], v["detail"])
            )
        print(
            "\nRename the parameter to the section-3 convention "
            "(docs/plan_docs/cli_consistency_remediation_plan.md section 3). "
            "Adding it to the baseline is not the fix: that list only shrinks."
        )
        return 1

    print("\nPASS: no violations outside the baseline allow-list (%d entries)."
          % len(baseline))
    return 0


if __name__ == "__main__":
    sys.exit(main())
