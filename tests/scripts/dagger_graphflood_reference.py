#!/usr/bin/env python3
"""DAGGER GraphFlood adapter for lsdtt3 hydraulic parity harness.

Implements the CLI contract expected by
tests/scripts/run_hydraulic_dagger_benchmark.py and writes:
- <out-prefix>_water_depth.tif
- <out-prefix>_water_surface.tif
- <out-prefix>_qout.tif

When DAGGER python bindings are unavailable, exits 77 so CTest records skip.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
from osgeo import gdal


SKIP_EXIT = 77


def die_skip(msg: str) -> int:
    print(f"SKIP: {msg}", file=sys.stderr)
    return SKIP_EXIT


def load_dem(path: pathlib.Path) -> tuple[
    np.ndarray,
    np.ndarray,
    tuple,
    str,
    float,
    float,
    float,
    float,
    float,
]:
    gdal.UseExceptions()
    ds = gdal.Open(str(path))
    if ds is None:
        raise RuntimeError(f"could not open DEM: {path}")

    band = ds.GetRasterBand(1)
    z = band.ReadAsArray().astype("float64")
    ndv = band.GetNoDataValue()

    valid = np.isfinite(z)
    if ndv is not None:
        valid &= z != ndv

    gt = ds.GetGeoTransform()
    proj = ds.GetProjectionRef() or ""
    dx = float(gt[1])
    dy = float(abs(gt[5]))
    x_min = float(gt[0])
    y_max = float(gt[3])
    return z, valid, tuple(gt), proj, ndv if ndv is not None else -9999.0, dx, dy, x_min, y_max


def write_raster(path: pathlib.Path, arr: np.ndarray, valid: np.ndarray,
                 gt: tuple, proj: str, ndv: float) -> None:
    driver = gdal.GetDriverByName("GTiff")
    if driver is None:
        raise RuntimeError("GDAL GTiff driver not available")

    ny, nx = arr.shape
    ds = driver.Create(str(path), nx, ny, 1, gdal.GDT_Float32)
    if ds is None:
        raise RuntimeError(f"failed to create output raster: {path}")

    ds.SetGeoTransform(gt)
    ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(float(ndv))

    out = np.array(arr, dtype="float32", copy=True)
    out[~valid] = np.float32(ndv)
    band.WriteArray(out)
    band.FlushCache()
    ds.FlushCache()
    ds = None


def import_dagger(dagger_root: pathlib.Path):
    primary = dagger_root / "wrappers" / "python"
    if primary.exists():
        sys.path.insert(0, str(primary))

    try:
        import dagger  # type: ignore
    except Exception as exc:
        raise RuntimeError(f"cannot import dagger module: {exc}") from exc

    if not (hasattr(dagger, "graph") and hasattr(dagger, "graphflood")):
        raise RuntimeError("imported dagger module does not expose graph/graphflood API")

    return dagger


def instantiate_graph(dagger, con):
    if hasattr(dagger, "graph"):
        try:
            return dagger.graph(con)
        except Exception:
            pass
        try:
            g = dagger.graph(con.nxy(), 8)
            g.init_graph(con)
            return g
        except Exception:
            pass
    raise RuntimeError("unsupported dagger.graph constructor/API")


def compute_graph(graph, flat_topo: np.ndarray) -> np.ndarray:
    for args in (
        (flat_topo, False, False),
        (flat_topo, False),
        (flat_topo,),
    ):
        try:
            out = graph.compute_graph(*args)
            if out is None:
                return flat_topo
            return np.asarray(out, dtype="float64")
        except Exception:
            continue
    return flat_topo


def as_bool(v: str) -> bool:
    return str(v).strip().lower() in {"1", "true", "yes", "on"}


def run_model(dagger, args: argparse.Namespace, z2d: np.ndarray, valid2d: np.ndarray,
              gt: tuple, proj: str, ndv: float, dx: float, dy: float,
              x_min: float, y_max: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ny, nx = z2d.shape
    topo = np.array(z2d, copy=True)
    topo[~valid2d] = ndv
    flat_topo = topo.ravel().astype("float64")

    con = dagger.D8N(nx, ny, dx, dy, x_min, y_max)
    if hasattr(con, "set_default_boundaries"):
        try:
            con.set_default_boundaries("4edges")
        except Exception:
            pass

    graph = instantiate_graph(dagger, con)
    processed_topo = compute_graph(graph, flat_topo)

    if not hasattr(dagger, "graphflood"):
        raise RuntimeError("dagger module has no graphflood class")
    gf = dagger.graphflood(graph, con)

    if args.flow_method == "mfd" and hasattr(gf, "enable_MFD"):
        gf.enable_MFD()
    elif hasattr(gf, "enable_SFD"):
        gf.enable_SFD()

    if hasattr(gf, "set_mannings"):
        gf.set_mannings(float(args.mannings_n))
    if hasattr(gf, "set_dt_hydro"):
        gf.set_dt_hydro(float(args.dt))
    if hasattr(gf, "set_topological_number"):
        gf.set_topological_number(float(args.mfd_exponent))

    gf.set_topo(processed_topo)
    gf.set_hw(np.zeros_like(processed_topo))

    if hasattr(gf, "set_water_input_by_constant_precipitation_rate"):
        gf.set_water_input_by_constant_precipitation_rate(float(args.rainfall_rate))

    if hasattr(gf, "enable_Qwout_recording"):
        gf.enable_Qwout_recording()

    dynamic = as_bool(args.dynamic_routing)
    for _ in range(int(args.n_steps)):
        if dynamic and hasattr(gf, "dynarun"):
            gf.dynarun()
        elif hasattr(gf, "run_hydro_only"):
            gf.run_hydro_only()
        elif hasattr(gf, "run"):
            gf.run()
        else:
            raise RuntimeError("graphflood has no runnable step method")

    hw = np.asarray(gf.get_hw(), dtype="float64").reshape((ny, nx))
    ws = np.asarray(gf.get_surface_topo(), dtype="float64").reshape((ny, nx))

    qout_flat = None
    if hasattr(gf, "get_Qwout_recording"):
        try:
            qout_flat = np.asarray(gf.get_Qwout_recording(), dtype="float64")
        except Exception:
            qout_flat = None
    if qout_flat is None and hasattr(gf, "get_Qwin"):
        qout_flat = np.asarray(gf.get_Qwin(), dtype="float64")
    if qout_flat is None:
        qout_flat = np.zeros(ny * nx, dtype="float64")
    qout = qout_flat.reshape((ny, nx))

    return hw, ws, qout


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dagger-root", required=True)
    ap.add_argument("--dem", required=True)
    ap.add_argument("--scenario", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out-prefix", default="dagger")
    ap.add_argument("--n-steps", required=True)
    ap.add_argument("--dt", required=True)
    ap.add_argument("--mannings-n", required=True)
    ap.add_argument("--rainfall-rate", required=True)
    ap.add_argument("--flow-method", required=True)
    ap.add_argument("--mfd-exponent", required=True)
    ap.add_argument("--dynamic-routing", required=True)
    ap.add_argument("--reroute-every", required=True)
    ap.add_argument("--run-mode", required=True)
    return ap.parse_args()


def main() -> int:
    args = parse_args()

    dagger_root = pathlib.Path(args.dagger_root).resolve()
    dem_path = pathlib.Path(args.dem).resolve()
    out_dir = pathlib.Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not dagger_root.exists():
        return die_skip(f"dagger root not found: {dagger_root}")
    if not dem_path.exists():
        return die_skip(f"dem not found: {dem_path}")

    try:
        dagger = import_dagger(dagger_root)
    except Exception as exc:
        return die_skip(str(exc))

    try:
        z2d, valid2d, gt, proj, ndv, dx, dy, x_min, y_max = load_dem(dem_path)
        hw, ws, qout = run_model(
            dagger, args, z2d, valid2d, gt, proj, ndv, dx, dy, x_min, y_max
        )

        write_raster(out_dir / f"{args.out_prefix}_water_depth.tif", hw, valid2d, gt, proj, ndv)
        write_raster(out_dir / f"{args.out_prefix}_water_surface.tif", ws, valid2d, gt, proj, ndv)
        write_raster(out_dir / f"{args.out_prefix}_qout.tif", qout, valid2d, gt, proj, ndv)
    except Exception as exc:
        msg = str(exc)
        dep_markers = (
            "_ARRAY_API",
            "numpy.core.multiarray failed to import",
            "cannot import dagger module",
            "No module named",
        )
        if any(marker in msg for marker in dep_markers):
            return die_skip(msg)
        print(f"FAIL: DAGGER adapter failed: {exc}", file=sys.stderr)
        return 1

    print("DAGGER adapter run completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
