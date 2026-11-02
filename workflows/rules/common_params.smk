# Helpers shared by the per-item modules (valley_widths, chi_region,
# lithology_channels).
#
# Those modules pass CLI parameters through from the config verbatim: a block
# such as `valley_widths.fluvial_surfaces` holds `key: value` pairs written
# straight into the CLI's parameter file. Any parameter in the CLI's reference
# can therefore be set from the config without touching a rule, and the CLI's
# strict parser rejects a misspelt key with its own message. The workflow
# itself sets only the file I/O and the keys that chain one CLI to the next;
# setting one of those in the config is an error rather than a silent override.

import glob as _glob
import os as _os
from pathlib import Path

IO_KEYS = ("read_path", "read_fname", "write_path", "write_fname")


def workflow_path(path):
    """Absolute path; relative paths are relative to the workflows/ directory."""
    return _os.path.abspath(_os.path.expanduser(str(path)))


def lsdtt_exe(key):
    exes = config.get("executables", {})
    if key not in exes:
        raise ValueError(
            f"executables.{key} is not set in the config: give the path to the "
            f"CLI, or its bare name if it is on your PATH.")
    return exes[key]


def _param_value(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple)):
        return ",".join(_param_value(v) for v in value)
    return str(value)


def merged_cli_parameters(where, *blocks):
    """Merge CLI parameter blocks left to right; later blocks win.

    `where` names the config location for error messages. Each block is a dict
    or None. Keys with a value of None are dropped, so an item can unset a
    default.
    """
    merged = {}
    for block in blocks:
        if block is None:
            continue
        if not isinstance(block, dict):
            raise ValueError(f"{where}: expected a mapping of CLI parameters, got {block!r}")
        merged.update(block)
    return {k: v for k, v in merged.items() if v is not None}


def render_param_file(where, io, user, forced=None, comment=None):
    """Parameter-file text: io keys, then the forced chaining keys, then user keys.

    A user key that collides with an io or forced key is an error, because the
    workflow depends on those values to find the next stage's inputs.
    """
    forced = forced or {}
    clashes = sorted(k for k in user if k in io or k in forced)
    if clashes:
        raise ValueError(
            f"{where}: {', '.join(clashes)} set by the workflow and cannot be set "
            f"in the config (it links one stage to the next).")
    lines = []
    if comment:
        lines.extend(f"# {line}" for line in comment.splitlines())
    lines.extend(f"{k}: {_param_value(v)}" for k, v in io.items())
    if forced:
        lines.append("")
        lines.append("# set by the workflow")
        lines.extend(f"{k}: {_param_value(v)}" for k, v in forced.items())
    if user:
        lines.append("")
        lines.append("# from the config")
        lines.extend(f"{k}: {_param_value(v)}" for k, v in user.items())
    return "\n".join(lines) + "\n"


def io_block(read_file, write_dir, write_fname):
    read_file = workflow_path(read_file)
    return {
        "read_path": _os.path.dirname(read_file) + _os.sep,
        "read_fname": _os.path.basename(read_file),
        "write_path": workflow_path(write_dir) + _os.sep,
        "write_fname": write_fname,
    }


def items_from_config(section_name, item_key, dir_key, glob_key):
    """Items of a per-item module as {name: entry}.

    Entries listed under `item_key` win. With none listed, every file matching
    `glob_key` in `dir_key` becomes an item named after its file stem.
    """
    section = config[section_name]
    listed = section.get(item_key) or {}
    if listed:
        out = {}
        for name, entry in listed.items():
            entry = dict(entry or {})
            if "file" not in entry and "fetch" not in entry:
                entry["file"] = f"{name}.tif"
            out[str(name)] = entry
        return out
    directory = section.get(dir_key)
    pattern = section.get(glob_key, "*.tif")
    if not directory:
        raise ValueError(f"{section_name}: list {item_key}, or set {dir_key} to a folder of DEMs")
    files = sorted(_glob.glob(_os.path.join(workflow_path(directory), pattern)))
    if not files:
        raise ValueError(
            f"{section_name}: no file matches {pattern} in {workflow_path(directory)}")
    return {_os.path.splitext(_os.path.basename(f))[0]: {"file": f} for f in files}


def item_file(section_name, dir_key, entry):
    path = entry["file"]
    if _os.path.isabs(str(path)):
        return path
    directory = config[section_name].get(dir_key) or "."
    return _os.path.join(workflow_path(directory), str(path))


def fetch_param_text(where, fetch, write_dir, write_fname):
    """lsdtt-fetch-raster parameters for one item's `fetch:` block."""
    user = merged_cli_parameters(where, fetch)
    if "bounds_fname" in user:
        user["bounds_fname"] = workflow_path(user["bounds_fname"])
    io = {"write_path": workflow_path(write_dir) + _os.sep, "write_fname": write_fname}
    return render_param_file(where, io, user)


def item_name_pattern(items):
    """Wildcard regex matching exactly these item names (nothing when empty).

    Keeps one module's per-item rules from competing with another's for the
    same output path.
    """
    import re
    if not items:
        return "(?!x)x"
    return "|".join(re.escape(name) for name in sorted(items, key=len, reverse=True))


# Each CLI's log is kept; on failure its tail is printed so the reason shows
# in the Snakemake output, not only in the log.
CLI_SHELL = ("mkdir -p {params.run_dir} && "
             "{params.exe} {params.run_dir} {params.param_name} > {log} 2>&1 "
             "|| {{ tail -n 25 {log}; exit 1; }}")




# --- figures (lsdviztools3) -------------------------------------------------
#
# Figures come from lsdviztools3's lsdviz3-* commands, which are a separate
# install with a much larger environment than the workflows one. The config's
# `figures` block says whether to draw them:
#
#   figures:
#     enabled: auto          # auto: draw them if lsdviztools3 is found, else skip
#                            # true: lsdviztools3 is required; stop if it is missing
#                            # false: never draw them, and do not look
#     lsdviz3_command_prefix: ""   # put before every lsdviz3-* command, e.g.
#                                  # "pixi run --manifest-path ../../lsdviztools3/pixi.toml"
#
# lsdviztools3 is looked for by running `lsdviz3-drape --help` (behind the
# prefix) once per Snakemake invocation, in the environment the figure rules
# use: the GDAL/PROJ variables of the workflows environment removed and no
# lsdtt setenv script, so a Python stack from another environment is not
# pointed at the wrong data.

import shlex as _shlex
import subprocess as _subprocess
import sys as _sys

_GDAL_PROJ_VARS = ("GDAL_DATA", "GDAL_DRIVER_PATH", "PROJ_DATA", "PROJ_LIB", "PROJ_NETWORK")
_FIGURES_STATE = {}


def _figure_env():
    return {k: v for k, v in _os.environ.items() if k not in _GDAL_PROJ_VARS}


def lsdviz3_argv(command, *args):
    prefix = (config.get("figures") or {}).get("lsdviz3_command_prefix") or ""
    return _shlex.split(prefix) + [command] + [str(a) for a in args]


def _figures_mode():
    mode = (config.get("figures") or {}).get("enabled", "auto")
    if isinstance(mode, bool):
        return "true" if mode else "false"
    mode = str(mode).strip().lower()
    if mode not in ("auto", "true", "false"):
        raise ValueError("figures.enabled must be auto, true or false")
    return mode


def figures_enabled():
    """Whether this run draws figures; looks for lsdviztools3 once."""
    if "enabled" in _FIGURES_STATE:
        return _FIGURES_STATE["enabled"]
    mode = _figures_mode()
    if mode == "false":
        _FIGURES_STATE["enabled"] = False
        return False
    # Snakemake re-reads the workflow in a child process for some jobs; the
    # first look is passed down through the environment instead of repeated.
    cache_key = "LSDTT_WORKFLOW_LSDVIZ3_FOUND"
    cached = _os.environ.get(cache_key)
    if cached in ("0", "1") and _os.environ.get(cache_key + "_ARGV") == " ".join(
            lsdviz3_argv("lsdviz3-drape")):
        found = cached == "1"
        if mode == "true" and not found:
            raise ValueError("figures.enabled is true but lsdviztools3 was not found")
        _FIGURES_STATE["enabled"] = found
        return found
    argv = lsdviz3_argv("lsdviz3-drape", "--help")
    try:
        probe = _subprocess.run(argv, env=_figure_env(), capture_output=True,
                                text=True, timeout=600)
        found = probe.returncode == 0
        reason = "" if found else (probe.stderr or probe.stdout).strip().splitlines()[-1:]
    except (OSError, _subprocess.TimeoutExpired) as exc:
        found, reason = False, [str(exc)]
    if not found:
        detail = f"`{' '.join(argv)}` failed: {' '.join(reason) or 'no output'}"
        hint = ("Install lsdviztools3 so its lsdviz3-* commands are on your PATH, or set "
                "figures.lsdviz3_command_prefix (e.g. \"pixi run --manifest-path "
                "/path/to/lsdviztools3/pixi.toml\"), or set figures.enabled: false.")
        if mode == "true":
            raise ValueError(f"figures.enabled is true but lsdviztools3 was not found: {detail}. {hint}")
        print(f"Figures skipped: lsdviztools3 not found ({detail}). {hint}", file=_sys.stderr)
    else:
        print("Figures: lsdviztools3 found; figures will be drawn.", file=_sys.stderr)
    _os.environ[cache_key] = "1" if found else "0"
    _os.environ[cache_key + "_ARGV"] = " ".join(lsdviz3_argv("lsdviz3-drape"))
    _FIGURES_STATE["enabled"] = found
    return found


def run_lsdviz3(log, command, *args):
    """Run one lsdviz3-* command, its output to `log`; print the tail on failure."""
    argv = lsdviz3_argv(command, *args)
    Path(log).parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w", encoding="utf-8") as f:
        f.write("$ " + " ".join(_shlex.quote(a) for a in argv) + "\n")
        f.flush()
        result = _subprocess.run(argv, env=_figure_env(), stdout=f, stderr=_subprocess.STDOUT)
    if result.returncode != 0:
        tail = Path(log).read_text(encoding="utf-8").splitlines()[-25:]
        print("\n".join(tail), file=_sys.stderr)
        raise RuntimeError(f"{command} failed (exit {result.returncode}); log: {log}")


def lsdviz3_has(command):
    """Whether this lsdviztools3 provides `command` (checked once, cached like the probe).

    For commands newer than lsdviz3-drape: an older lsdviztools3 passes the main
    probe but lacks them. auto skips those figures with a message; true stops.
    """
    if not figures_enabled():
        return False
    key = f"has:{command}"
    if key in _FIGURES_STATE:
        return _FIGURES_STATE[key]
    env_key = "LSDTT_WORKFLOW_HAS_" + command.upper().replace("-", "_")
    cached = _os.environ.get(env_key)
    if cached in ("0", "1"):
        found = cached == "1"
    else:
        try:
            found = _subprocess.run(lsdviz3_argv(command, "--help"), env=_figure_env(),
                                    capture_output=True, timeout=600).returncode == 0
        except (OSError, _subprocess.TimeoutExpired):
            found = False
        _os.environ[env_key] = "1" if found else "0"
        if not found:
            message = (f"lsdviztools3 was found but has no {command}; update lsdviztools3 "
                       f"(git pull, then pixi run install) for these figures")
            if _figures_mode() == "true":
                raise ValueError(f"figures.enabled is true but {message}")
            print(f"Figures skipped: {message}.", file=_sys.stderr)
    _FIGURES_STATE[key] = found
    return found
