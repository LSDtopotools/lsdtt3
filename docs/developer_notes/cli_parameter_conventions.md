# CLI Parameter Naming Conventions

This document establishes the standard naming conventions for CLI parameters in `lsdtt3`.
If you are a developer and make a pull request failing to follow these instructions the lead developer will wag his finger at you and say "bad!". 

---

## Overview

Every lsdtt3 CLI declares its parameters through `lsdtt_param_parser`. The naming
standard below ensures that:

- A parameter's name tells the developer what it does and what units it uses.
- Renaming a parameter is non-breaking for existing `.param` files through the
  deprecated-alias mechanism.
- New contributors follow the existing patterns rather than inventing new ones.

The standard applies to all parameters. 
A CTest guard checks every help-registered parameter in every CLI, and fails on
any violation (see [Conformance testing](#conformance-testing)).

---

## The standard

| Concept | Canonical form | Notes |
|---|---|---|
| Channel initiation | `threshold_contributing_pixels` (int, pixels, 1000) | Matches lsdtopotools2 for cross-version continuity |
| Genuine areas | `<what>_area_threshold_m2` (float) | Never spelled the same as the pixel form |
| Input already filled, fill is mandatory | `raster_is_filled` (bool, false) | **Concept A** — a compute-saving assertion about the input |
| Fill is optional; false routes into sinks | `fill_before_routing` (bool, true) | **Concept B** — an analytical choice |
| Fill minimum slope | `fill_min_slope` (float, m/m, 1e-4) | Stage prefix only where a CLI fills more than one surface |
| Raster output format | `write_output_format` | One per CLI |
| Vector output format | `write_vector_format` | One per CLI; per-output overrides only where a container genuinely differs |
| Boolean raster output | `write_<name>_raster` | Plural `_rasters` where a set is written |
| Input file | `<what>_fname` | Never `_file`, `_filename`, `_csv`; format is a separate `<what>_input_format` |
| Fitting / moving window | `<stage>_<quantity>_radius_metres` | Radius, always; `polyfit_lambda_metres` reserved for half-wavelength |
| Snap radius | `<what>_snap_max_search_radius_pixels` (+ `_min_`) | No `tol`; metric tolerances `_tolerance_metres` |
| Units | `_metres`, `_pixels`, `_degrees` spelled out | Never `_m`, `_px`, `_deg` |
| Stage switch | `run_<stage>` | Not `do_`, `use_`, `apply_`, `compute_` |
| Output switch | `write_<output>` | Never `print_` for something that writes a file |
| Bounds | `min_` / `max_` prefix | Not `minimum_`, not `threshold_<x>` |
| Subject-qualified pixel-count filters | `<subject>_minimum_contributing_pixels` | **Deliberate exception to the row above.** The subject prefix (`basin_`, `crossing_`, `side_tributary_`) names the object being filtered and reads as plain English. Do not shorten to `min_` |

---

## Two concepts: filling a raster

The two fill-related parameters below ask different questions and must never be merged.
Each help text must state which question it is asking.

### Concept A: `raster_is_filled` — input assertion

**The question:** Has the input DEM already been filled (sinks removed)?

**The purpose:** A compute-saving assertion, not an analytical choice. The operation
requires a filled surface. If you set this to `true` but supply an unfilled DEM,
you get **bad routing** — divergences, loops, or sinks that should have been
breached. This is silent data corruption, not just a different answer.

**Default:** `false` (do the work, be safe)

**Implementation pattern:**
```cpp
if (!bf["raster_is_filled"])
{
    dem = lsdtt::fill_pits(dem);
}
```

**Legacy alias:** `flow_dem_is_filled` in `feature-finder` maps to
`raster_is_filled` with a deprecation warning.

### Concept B: `fill_before_routing` — analytical choice

**The question:** Should flow routing terminate at sinks, or should the surface
be filled first?

**The purpose:** A genuine analysis choice. Setting to `false` means flow is
routed on the unfilled surface and terminates at real sinks. This is intentional,
not an error.

**Default:** `true` (fill before routing)

**Implementation pattern:**
```cpp
dem = bf.at("fill_before_routing") ? lsdtt::fill_pits(dem_raw) : dem_raw;
```

**Legacy alias:** `fill_dem_for_flow_routing` in `fluvial-surfaces` maps to
`fill_before_routing` with a deprecation warning.

### Note: `fill_initial_surface`

The parameter `fill_initial_surface` in `lsdtt-lem` is neither concept. It applies
to the synthetic initial surface, not to routing. It stays as it is.

---

## Shared parameter declaration layer

The shared declaration layer in `src/core/lsdtt_common_parameters.hpp`
eliminates duplication of common parameter names and help text.

The layer exposes functions that populate a CLI's default maps and help registry
in one call. This ensures there is exactly one place each shared name is spelled:

```cpp
lsdtt::declare_routing_core_parameters(int_map, string_map, help_map);
    // Declares: boundary_conditions, n_threads

lsdtt::declare_raster_is_filled_parameter(bool_map, help_map);
lsdtt::declare_fill_before_routing_parameter(bool_map, help_map);
lsdtt::declare_fill_min_slope_parameter(float_map, help_map);
lsdtt::declare_channel_initiation_parameters(int_map, help_map);
    // Each declaration is independent; use only the concepts the CLI supports.
    // Every function accepts an optional prefix for multi-stage CLIs.

lsdtt::declare_raster_output_format_parameter(string_defaults, help);
lsdtt::declare_vector_output_format_parameter(string_defaults, help);
lsdtt::declare_spatial_index_parameter(bool_defaults, help);
    // Declares: write_output_format, write_vector_format,
    //           calculate_spatial_index

lsdtt::declare_basin_outlet_parameters(...);
    // Declares: basin_outlet_fname, basin_outlet_snapping_method,
    //           basin_outlet_snap_min_search_radius_pixels,
    //           basin_outlet_snap_max_search_radius_pixels,
    //           basin_outlet_contributing_pixels_threshold

lsdtt::declare_surface_fitting_parameters(...);
    // Declares: polyfit_lambda_metres, polyfit_metrics, compute_backend

lsdtt::declare_synthetic_domain_parameters(...);
    // Declares: n_rows, n_cols, cell_size_metres, random_seed,
    //           synthetic_crs_wkt

lsdtt::declare_synthetic_noise_parameters(...);
    // Declares: noise_type, noise_target_stddev_metres, noise_relief_metres,
    //           diamond_square_roughness, diamond_square_initial_scale

// Synthetic concepts are also independently selectable:
lsdtt::declare_synthetic_dimensions_parameters(int_map, help_map);
lsdtt::declare_synthetic_cell_size_parameter(float_map, help_map);
lsdtt::declare_synthetic_random_seed_parameter(int_map, help_map);
lsdtt::declare_synthetic_crs_wkt_parameter(string_map, help_map);
lsdtt::declare_synthetic_noise_type_parameter(string_map, help_map);
lsdtt::declare_synthetic_noise_target_stddev_parameter(float_map, help_map);
lsdtt::declare_synthetic_noise_relief_parameter(float_map, help_map);
lsdtt::declare_diamond_square_parameters(float_map, help_map);
```

CLIs use these declarations for the shared concepts they expose.

The flow and synthetic declarations are intentionally granular. `raster_is_filled` and
`fill_before_routing` are different questions, and many CLIs support only one.
Similarly, a CLI may consume a pre-filled surface without owning a fill-slope
parameter, or may define channels without exposing routing controls. Do not add
unsupported keys merely to use a larger declaration group. In particular, do
not add `synthetic_crs_wkt` to a CLI with an EPSG/centre-coordinate CRS model or
`noise_type` to a harness whose independent write switches can run together.

---

## Deprecated-alias adoption pattern

When a parameter is renamed or retired, the old name must be handled gracefully
so existing `.param` files do not silently break. The `lsdtt_param_parser` provides
five mechanisms, and the choice depends on what the parameter did.

**Most important rule:** The old key keeps working with a warning naming the new one.
Every rename is non-breaking for users. The cost of a rename is one line in an
alias map. There is no flag day.

However, the alias map only rescues a key that is *read by the CLI*. A key string
hard-coded in a script, notebook, or the QGIS plugin is not covered. A rename still
needs a tree sweep for those cases.

### Five maps: the decision rule

#### `set_deprecated_aliases()` — the parameter still exists under a new name

Use this when renaming a parameter. The value is carried across automatically,
and the old name is erased so the parser does not report it as unknown.

**When to use it:** The parameter exists; only the name changed. The user's
`.param` file continues to work, and the warning tells them the new name.

**Example:** `source_pixels_threshold` renamed to `threshold_contributing_pixels`.
Both names work; if both are supplied, the new one wins and the conflict is reported.

```cpp
pp.set_deprecated_aliases({
    {"source_pixels_threshold", "threshold_contributing_pixels"},
    {"old_fill_param", "fill_before_routing"},
});
```

**Warning message:**
```
WARNING: '<old_key>' is deprecated; use '<new_key>'. Its value has been applied
to '<new_key>' for this run.
```

If both old and new are present:
```
WARNING: '<old_key>' is deprecated and '<new_key>' was also supplied. Using
'<new_key>'; the deprecated key is ignored.
```

The warning contract names both keys; it does not name the executable. The
parser owns the parameter filename but is not given CLI identity, and adding
tool identity solely to decorate this warning would couple otherwise generic
key migration to presentation.

`set_deprecated_aliases()` replaces the parser's complete alias map. A CLI with
more than one migration must therefore pass all mappings in one call. Do not
call it repeatedly: a later call would discard the earlier aliases.

#### `set_deprecated_numeric_aliases()` — the name and numeric meaning changed

Use this narrow typed facility when a rename also requires an affine numeric
conversion. It is intended for unit or geometric convention migrations, such as
an old window diameter becoming a canonical radius. It is not a general callback
mechanism.

```cpp
pp.set_deprecated_numeric_aliases({
    {"slope_window_diameter_metres",
     {"slope_window_radius_metres", 0.5, 0.0, true}},
});
```

The policy fields are `new_key`, `scale`, `offset`, and `preserve_zero`; the
rewrite is `new = old * scale + offset`, except that an old zero remains exactly
zero when `preserve_zero` is true. The old value must be a complete finite number.
Malformed values, `nan`, `inf`, a non-finite result, and invalid policies stop
with a fatal diagnostic before typed parsing and analysis.

The canonical key wins when both are supplied, as for an ordinary alias. The old
value is still validated, then erased so it is not reported as unknown. Warnings
name both keys and print the transformation. Rewritten values use classic-locale
`double` `max_digits10` formatting, preserving sufficient text precision for
subsequent float or double parsing.

`set_deprecated_numeric_aliases()` replaces the complete transformed-alias map,
independently of `set_deprecated_aliases()`: multiple transformed entries and the
ordinary alias map coexist. Transformed aliases are applied before ordinary
aliases. Registering the same old key in both maps is a fatal configuration
error; do not rely on ordering to choose a policy.

#### `set_deprecated_compatibility_aliases()` — consolidated selectors

Use this only when several old per-output selectors are consolidated into one
canonical selector and old files may set those selectors differently. Each old
key maps to a distinct hidden storage key plus the public replacement named in
the warning. The old value remains attached to its original output for that run
without exposing fake public overrides.

Compatibility storage keys default to `same_as_write_vector_format` or
`same_as_write_output_format` and are resolved after parsing with
`resolve_output_format_override()`. They are implementation details and must not
be added to help or user documentation.

#### `set_removed_keys()` — the parameter is gone; fatal

Use this only for parameters whose whole purpose was to **turn something on**.
Silently ignoring such a key leaves the user believing an analysis ran when it did not.

**When to use it:** The parameter controlled an analysis mode or enabled a feature.
Removing it without the user noticing is a silent failure. Stop the run with a
clear migration message instead.

**Example:** If a parameter like `run_feature_analysis: true` is removed, a user
with an old `.param` file would expect that analysis to run. Silently skipping it
is wrong. Make it fatal so they fix the file.

```cpp
pp.set_removed_keys({
    {"run_obsolete_analysis",
     "The 'run_obsolete_analysis' parameter has been removed. "
     "Use 'run_new_analysis: true' instead."},
});
```

**Behavior:** The run stops (exit 1) and prints all removed-key messages together,
so a user migrating an old file gets the whole list in one attempt rather than
one key per run.

#### `set_ignored_keys()` — the parameter is gone; accepted and ignored

Use this only for parameters that merely **tuned or disabled** a mechanism that
no longer exists. Ignoring such a key cannot change what the user believes happened.

**When to use it:** The parameter optimised or disabled something, and the mechanism
no longer exists. The analysis still runs, and runs better. Ignoring the key is fine.

**Example:** If a parameter like `legacy_optimization_level` is removed, ignoring it
does not hide a missing analysis — it just means a tuning knob is gone. The analysis
still produces the same output.

```cpp
pp.set_ignored_keys({
    {"legacy_optimization_level",
     "'legacy_optimization_level' is no longer used; the code now auto-optimizes. "
     "This parameter is ignored and may be removed from your .param file."},
});
```

**Behavior:** The key is accepted, ignored, and the run continues. The message
is printed so the user knows the key was not lost in a typo.

#### `set_report_unknown_keys()` — default on

This is on by default and is separate from the three maps above. Any key in the
`.param` file that is not declared in any default map and is not reserved
(read_path, read_fname, write_path, write_fname, boundary_conditions) is reported
as unknown. This catches typos.

A warning is printed but the run continues — a CLI may legitimately read a key
by another route, or a key may be intended for a plugin or workflow wrapper.

### Implementation in a CLI

Every CLI gains a deprecated-aliases block in its `main()`, immediately after
the default maps are filled but before `parse_all_parameters()` is called:

```cpp
// --- Parameter defaults and help ---
int_default_map["threshold_contributing_pixels"] = 1000;
help_map["threshold_contributing_pixels"] = {"int","1000",
    "Channel-initiation threshold (upstream pixels).",
    "...", "Flow routing", "param", ""};

// ... other parameters ...

// --- Deprecated aliases, removed keys, ignored keys ---
pp.set_deprecated_aliases({
    {"source_pixels_threshold", "threshold_contributing_pixels"},
    {"old_fill_param", "fill_before_routing"},
});
pp.set_removed_keys({
    {"run_obsolete_analysis",
     "The 'run_obsolete_analysis' parameter has been removed. "
     "Use 'run_new_analysis: true' instead."},
});
pp.set_ignored_keys({
    {"legacy_optimization_level",
     "'legacy_optimization_level' is no longer used; the code now auto-optimizes."},
});

// --- Parse all parameters ---
pp.parse_all_parameters(float_default_map, int_default_map,
                        bool_default_map, string_default_map);
```

For reference implementations, see:
- `cli/lsdtt-flat-surface-finder.cpp`
- `cli/lsdtt-hillslope-channel-coupling.cpp`
- `cli/lsdtt-internal-drainage.cpp`
- `cli/lsdtt-valley-centreline.cpp`

Every future rename adds one line to the alias map. That is the cost and the
benefit.

---

## Conformance testing

The CTest case `unit_cli_parameter_conventions` mechanically enforces the
checkable subset of these rules. It runs
`tests/scripts/check_cli_parameter_conventions.py`, which runs every built
`lsdtt-*` binary with `-h json` and checks each declared parameter:

| Rule id | Checks |
|---|---|
| `unit_suffix` | no name ends in `_m`, `_px` or `_deg` (`_area_threshold_m2` is allowed) |
| `print_writes_file` | no `print_*` parameter whose help says it writes a file |
| `bool_raster_suffix` | a bool switch that writes a raster is `write_<name>_raster` / `_rasters` |
| `duplicate_format` | at most one raster-format and one vector-format parameter per CLI |
| `minimum_prefix` | bounds use `min_` / `max_`, except `<subject>_minimum_contributing_pixels` |
| `tol_abbreviation` | no `tol` token; use `_tolerance_metres` or `_snap_max_search_radius_pixels` |
| `help_json_unparseable` | the CLI's `-h json` output could be parsed |

The rules are heuristics over parameter names and help text. For example,
`bool_raster_suffix` only fires when the help text mentions a raster and no
vector or table output.

**The test fails on every violation.** The baseline allow-list
(`tests/data/cli_parameter_conventions_baseline.json`) was used while the
remediation phases were in progress, and is now empty. Do not add entries to
it. Fix the name instead. The only sanctioned exceptions are recorded in the
script itself:

- `REVIEWED_FORMAT_OVERRIDES`: per-CLI format selectors that govern a genuinely
  different output domain or stage. Examples are
  `wiener_diagnostics_write_output_format` in `lsdtt-channel-extraction` and
  `diagnostics_write_output_format` in `lsdtt-advection-diffusion`, for
  visualisation frames. Add a new one only after review, with a comment
  saying why.
- `REVIEWED_INPUT_FORMATS`: legacy input selectors that predate the
  `_input_format` suffix.

The test also refuses to run if the build tree does not match `cli/`: a CLI
source with no binary, or a binary older than its source, is a fatal error,
because checking stale binaries would measure parameters that no longer exist.
If it fails this way, rebuild everything (`cmake --build build -j`, re-running
`cmake` first if a new CLI was added) and run it again.

A companion test, `unit_cli_reference_is_current`, checks that the generated
parameter tables in `docs/cli_reference/` match the binaries. After adding or
renaming a parameter, regenerate them with
`python3 tools/generate_cli_reference.py --update --bin-dir build`.

---

## The exception: `<subject>_minimum_contributing_pixels`

Three parameters form a deliberate exception to the `min_` / `max_` prefix rule:

- `basin_minimum_contributing_pixels`
- `crossing_minimum_contributing_pixels`
- `side_tributary_minimum_contributing_pixels`

The subject prefix (`basin_`, `crossing_`, `side_tributary_`) says what is being
filtered and the full word `minimum` reads as plain English: "the minimum
contributing pixels for a [subject]". Do not shorten these to `min_` — the
readability of the English matters more than the prefix consistency in this case.

This is an explicit exception, not a conflict with the row above it in the standard
table.

---

## See also

- Architecture document: `docs/developer_notes/lsdtt3_architecture.md` — parameter
  declaration in CLIs
- Developer onboarding: `docs/developer_notes/developer_onboarding.md` — checklist
  for adding a parameter
- Parameter parser: `src/core/lsdtt_param_parser.hpp` — the engine that enforces
  these rules and applies aliases
