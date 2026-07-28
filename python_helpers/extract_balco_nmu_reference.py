#!/usr/bin/env python3
"""Extract Balco CRONUS v3 muon inventory grids into a frozen CSV fixture.

Reads `Nmu.mat`, reports discovered structure, then writes a flat CSV table used
by lsdtt3 tests to compare muon-inventory interpolation against the published v3
reference data.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
from datetime import datetime, timezone
from typing import Any, Dict, Tuple

import numpy as np
import scipy.io as sio


NUCLIDE_FIELDS = {
    "Be10": "N10quartz",
    "Al26": "N26quartz",
    "C14": "N14quartz",
}


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _first_scalar_string(value: Any) -> str:
    arr = np.asarray(value)
    if arr.size == 0:
        return ""
    x = arr.flat[0]
    if isinstance(x, bytes):
        return x.decode("utf-8", errors="replace")
    return str(x)


def _to_float_array(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=np.float64)


def load_nmu(path: str) -> Tuple[Dict[str, np.ndarray], str, str]:
    """Load Nmu struct from a MATLAB file.

    Returns: (field dict, version string, loader_name)
    """
    try:
        d = sio.loadmat(path, squeeze_me=False, struct_as_record=False)
        keys = [k for k in d.keys() if not k.startswith("__")]
        if keys != ["Nmu"]:
            raise ValueError(f"{path}: expected only 'Nmu', found {keys}")
        nmu = d["Nmu"][0, 0]
        fields: Dict[str, np.ndarray] = {}
        for name in nmu._fieldnames:
            fields[name] = np.asarray(getattr(nmu, name))
        version = _first_scalar_string(fields.get("version", ""))
        return fields, version, "scipy.io.loadmat"
    except Exception as scipy_err:
        try:
            import h5py  # type: ignore
        except Exception:
            raise RuntimeError(
                "scipy.io.loadmat failed and h5py is not available for MATLAB v7.3 fallback"
            ) from scipy_err

        try:
            with h5py.File(path, "r") as h5:
                if "Nmu" not in h5:
                    raise ValueError("h5py fallback: missing top-level key 'Nmu'")
                nmu_obj = h5["Nmu"]
                if not isinstance(nmu_obj, h5py.Group):
                    raise ValueError("h5py fallback expects Nmu to be a group")

                fields = {}
                for key in nmu_obj.keys():
                    obj = nmu_obj[key]
                    if isinstance(obj, h5py.Dataset):
                        fields[key] = np.asarray(obj)
                    else:
                        fields[key] = np.asarray(obj)

                version = _first_scalar_string(fields.get("version", ""))
                return fields, version, "h5py"
        except Exception as h5_err:
            raise RuntimeError(
                "Failed to load Nmu.mat via scipy.io.loadmat and h5py fallback"
            ) from h5_err


def require_regular_axes(p_axis: np.ndarray, loge_axis: np.ndarray) -> None:
    p1 = np.asarray(p_axis, dtype=np.float64).reshape(-1)
    l1 = np.asarray(loge_axis, dtype=np.float64).reshape(-1)
    if p1.size < 2 or l1.size < 2:
        raise ValueError("axes must have at least 2 points each")
    pd = np.diff(p1)
    ld = np.diff(l1)
    if not np.allclose(pd, pd[0], rtol=0.0, atol=1e-12):
        raise ValueError("p axis is not regular")
    if not np.allclose(ld, ld[0], rtol=0.0, atol=1e-12):
        raise ValueError("log10_erate axis is not regular")


def print_structure(fields: Dict[str, np.ndarray], version: str, loader_name: str) -> None:
    print(f"loader: {loader_name}")
    print(f"Nmu.version: {version}")
    print("discovered fields:")
    for name in sorted(fields.keys()):
        arr = np.asarray(fields[name])
        line = f"  {name}: shape={arr.shape}, dtype={arr.dtype}"
        if arr.size and arr.dtype.kind in "US":
            line += f", first={arr.flat[0]}"
        print(line)


def format_f(x: float) -> str:
    return format(float(x), ".17g")


def write_csv(
    out_path: str,
    source_path: str,
    source_sha: str,
    version: str,
    p_axis: np.ndarray,
    log10_erate_axis: np.ndarray,
    grids: Dict[str, np.ndarray],
    loader_name: str,
) -> int:
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    p_vals = np.asarray(p_axis, dtype=np.float64).reshape(-1)
    l_vals = np.asarray(log10_erate_axis, dtype=np.float64).reshape(-1)

    rows = 0
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        fh.write("# Balco CRONUS v3 Nmu reference table for lsdtt3 tests\n")
        fh.write("# FROZEN reference - refreshing it is a reviewable event.\n")
        fh.write(f"# source_path: {source_path}\n")
        fh.write(f"# source_sha256: {source_sha}\n")
        fh.write(f"# Nmu.version: {version}\n")
        fh.write(f"# extracted_utc: {timestamp}\n")
        fh.write(f"# loader: {loader_name}\n")
        fh.write("# grid_layout: rows=log10_erate, cols=p_hPa\n")
        fh.write("# interpolation_space: bilinear in (p_hPa, log10_erate)\n")
        fh.write("# p_hPa_axis_csv: " + ",".join(format_f(v) for v in p_vals) + "\n")
        fh.write(
            "# log10_erate_axis_csv: " + ",".join(format_f(v) for v in l_vals) + "\n"
        )
        fh.write("nuclide,p_hPa,log10_erate,N_mu_atoms_g\n")

        for nuclide in ("Be10", "Al26", "C14"):
            grid = np.asarray(grids[nuclide], dtype=np.float64)
            if grid.shape != (l_vals.size, p_vals.size):
                raise ValueError(
                    f"{nuclide} grid shape {grid.shape} does not match "
                    f"(n_log10_erate={l_vals.size}, n_p={p_vals.size})"
                )
            for iy, loge in enumerate(l_vals):
                for ix, p_hpa in enumerate(p_vals):
                    fh.write(
                        f"{nuclide},{format_f(p_hpa)},{format_f(loge)},{format_f(grid[iy, ix])}\n"
                    )
                    rows += 1
    return rows


def main(argv: list[str]) -> int:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nmu-mat",
        default="/mnt/c/LSDTopoTools/crn_reference/v3/data/Nmu.mat",
        help="Path to Balco v3 Nmu.mat",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(repo_root, "tests", "expected", "cosmo", "balco_v3_nmu_reference.csv"),
        help="Output CSV path",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow overwriting existing output",
    )
    args = parser.parse_args(argv)

    source_path = os.path.abspath(args.nmu_mat)
    out_path = os.path.abspath(args.output)

    if not os.path.isfile(source_path):
        raise FileNotFoundError(source_path)
    out_dir = os.path.dirname(out_path)
    if not os.path.isdir(out_dir):
        raise FileNotFoundError(out_dir)
    if os.path.exists(out_path) and not args.force:
        raise FileExistsError(f"{out_path} exists; pass --force to overwrite")

    fields, version, loader_name = load_nmu(source_path)
    print_structure(fields, version, loader_name)

    p_axis = _to_float_array(fields["p"])
    log10_erate_axis = _to_float_array(fields["logeg"])
    p_axis = p_axis.reshape(-1)
    log10_erate_axis = log10_erate_axis.reshape(-1)

    require_regular_axes(p_axis, log10_erate_axis)

    grids: Dict[str, np.ndarray] = {}
    for nuclide, field_name in NUCLIDE_FIELDS.items():
        if field_name not in fields:
            raise KeyError(f"missing field '{field_name}' for {nuclide}")
        grids[nuclide] = np.asarray(fields[field_name], dtype=np.float64)

    source_sha = sha256_of(source_path)
    rows = write_csv(
        out_path=out_path,
        source_path=source_path,
        source_sha=source_sha,
        version=version,
        p_axis=p_axis,
        log10_erate_axis=log10_erate_axis,
        grids=grids,
        loader_name=loader_name,
    )

    print(f"wrote: {out_path}")
    print(f"rows: {rows}")
    print(
        "axes: "
        f"p_hPa[{p_axis.size}]={format_f(p_axis[0])}..{format_f(p_axis[-1])}, "
        f"log10_erate[{log10_erate_axis.size}]={format_f(log10_erate_axis[0])}..{format_f(log10_erate_axis[-1])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
