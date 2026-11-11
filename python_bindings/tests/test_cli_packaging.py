"""Tests for the pip-installed vanilla CLI dispatch layer."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import re

import pytest

pytestmark = [pytest.mark.unit]


_CLI_PATH = Path(__file__).parents[1] / "src" / "lsdtt3" / "_cli.py"
_SPEC = importlib.util.spec_from_file_location("lsdtt3_cli_test", _CLI_PATH)
assert _SPEC is not None and _SPEC.loader is not None
_cli = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_cli)


@pytest.mark.negative
def test_driver_path_reports_missing_driver(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """A wheel without the executable fails with a direct useful error."""
    monkeypatch.setattr(_cli, "__file__", str(tmp_path / "_cli.py"))
    with pytest.raises(RuntimeError, match="does not include 'lsdtt-graphflood'"):
        _cli._driver_path("lsdtt-graphflood")


@pytest.mark.exact_solution
def test_driver_path_uses_platform_layout(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Drivers are discovered underneath the installed package directory."""
    monkeypatch.setattr(_cli, "__file__", str(tmp_path / "_cli.py"))
    suffix = ".exe" if _cli.os.name == "nt" else ""
    executable = tmp_path / "bin" / f"lsdtt-chi-analysis{suffix}"
    executable.parent.mkdir()
    executable.touch()
    assert _cli._driver_path("lsdtt-chi-analysis") == executable


@pytest.mark.exact_solution
def test_run_forwards_arguments_and_exit_code(monkeypatch: pytest.MonkeyPatch) -> None:
    """The Python console script preserves native CLI arguments and status."""
    executable = Path("/tmp/lsdtt-chi-analysis")
    monkeypatch.setattr(_cli, "_driver_path", lambda _: executable)
    monkeypatch.setattr(_cli.sys, "argv", ["lsdtt-chi-analysis", "params.txt"])

    calls: list[tuple[list[str], bool, object]] = []

    class Completed:
        returncode = 7

    def run(command: list[str], *, check: bool, env: object) -> Completed:
        calls.append((command, check, env))
        return Completed()

    monkeypatch.setattr(_cli.subprocess, "run", run)
    assert _cli._run("lsdtt-chi-analysis") == 7
    assert len(calls) == 1
    command, check, env = calls[0]
    assert command == [str(executable), "params.txt"] and check is False
    if _cli.os.name == "nt":
        # the wheel's DLL directories are added to PATH for the child .exe
        assert isinstance(env, dict) and "PATH" in env
    else:
        assert env is None


@pytest.mark.exact_solution
def test_main_derives_driver_from_console_script_name(monkeypatch: pytest.MonkeyPatch) -> None:
    """One generic console entry point dispatches every declared driver."""
    monkeypatch.setattr(_cli.sys, "argv", ["/venv/bin/lsdtt-surface-metrics", "p.txt"])
    monkeypatch.setattr(_cli, "_run", lambda name: 3 if name == "lsdtt-surface-metrics" else 1)
    assert _cli.main() == 3


@pytest.mark.conformance
def test_driver_inventory_matches_the_vanilla_cmake_build() -> None:
    """Every root-CMake vanilla driver has a wheel dispatcher declaration."""
    root = Path(__file__).parents[2]
    vanilla_cmake = (root / "CMakeLists.txt").read_text(encoding="utf-8")
    expected = set(re.findall(r"lsdtt_cli\((lsdtt-[\w-]+)", vanilla_cmake))
    assert _cli._DRIVERS == expected


@pytest.mark.conformance
def test_driver_inventory_matches_wheel_declarations() -> None:
    """CMake targets and console scripts cannot drift from the dispatcher."""
    root = Path(__file__).parents[1]
    wheel_cmake = (root / "CMakeLists.txt").read_text(encoding="utf-8")
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    cmake_drivers = set(re.findall(r"lsdtt_python_cli\((lsdtt-[\w-]+)", wheel_cmake))
    console_scripts = set(
        re.findall(r"^(lsdtt-[\w-]+)\s*=\s*\"lsdtt3\._cli:main\"$", pyproject, re.MULTILINE)
    )
    assert cmake_drivers == _cli._DRIVERS
    assert console_scripts == _cli._DRIVERS
    assert "lsdtt-synthetic-surfaces" in _cli._DRIVERS
    assert "lsdtt-raster-synthetic-tests" not in _cli._DRIVERS
