#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask, shapes
from shapely.geometry import shape

matplotlib.use("Agg")


def outlines_from_excavation(exc_path: Path, out_fgb: Path) -> None:
    with rasterio.open(exc_path) as ds:
        arr = ds.read(1)
        scar_mask = np.isfinite(arr) & (arr > 0.0)
        geoms = []
        mean_depths = []
        for geom, val in shapes(scar_mask.astype("uint8"), mask=scar_mask, transform=ds.transform):
            if int(val) != 1:
                continue
            poly = shape(geom)
            if poly.is_empty:
                continue
            geoms.append(poly)
            pix_mask = geometry_mask([geom], out_shape=arr.shape, transform=ds.transform, invert=True)
            d = arr[pix_mask]
            d = d[np.isfinite(d) & (d > 0.0)]
            mean_depths.append(float(d.mean()) if d.size else np.nan)

        gdf = gpd.GeoDataFrame(
            {
                "polygon_id": np.arange(1, len(geoms) + 1),
                "mean_exc_m": mean_depths,
            },
            geometry=geoms,
            crs=ds.crs,
        )

    if len(gdf) == 0:
        raise RuntimeError("No excavated pixels found in excavation raster.")

    gdf["area_m2"] = gdf.geometry.area
    out_fgb.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(out_fgb, driver="FlatGeobuf")


def chart_larsen_comparison(inv_path: Path, out_png: Path, alpha: float, gamma: float) -> None:
    inv = pd.read_csv(inv_path)
    inv["area_m2"] = pd.to_numeric(inv["area_m2"], errors="coerce")
    inv["volume_m3"] = pd.to_numeric(inv["volume_m3"], errors="coerce")
    inv = inv[np.isfinite(inv["area_m2"]) & np.isfinite(inv["volume_m3"]) & (inv["area_m2"] > 0.0) & (inv["volume_m3"] > 0.0)]
    if inv.empty:
        raise RuntimeError("No valid landslides in inventory for Larsen comparison chart.")

    a_min = max(float(inv["area_m2"].min()), 1.0)
    a_max = float(inv["area_m2"].max())
    a_line = np.logspace(np.log10(a_min), np.log10(a_max), 200)
    v_line = alpha * np.power(a_line, gamma)

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.scatter(inv["area_m2"], inv["volume_m3"], s=18, alpha=0.7, label="Realised landslides")
    ax.plot(a_line, v_line, color="black", linewidth=2.0, label=f"Reference: V={alpha:g} A^{gamma:g}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Area (m^2)")
    ax.set_ylabel("Volume (m^3)")
    ax.set_title("Single-event landslide size vs Larsen power law")
    ax.grid(True, which="both", linestyle=":", linewidth=0.6, alpha=0.6)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def chart_apparent_vs_true(sample_effect_path: Path, out_png: Path, true_rate_mm_kyr: float) -> None:
    df = pd.read_csv(sample_effect_path)
    for col in ("erate_control_mm_kyr", "erate_median_mm_kyr"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[np.isfinite(df["erate_control_mm_kyr"]) & np.isfinite(df["erate_median_mm_kyr"])]
    if df.empty:
        raise RuntimeError("No valid rows in sample-effect CSV for apparent-vs-true chart.")

    x = df["erate_control_mm_kyr"].to_numpy(dtype=float)
    y = df["erate_median_mm_kyr"].to_numpy(dtype=float)
    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    pad = 0.08 * (hi - lo if hi > lo else 1.0)
    lims = (lo - pad, hi + pad)

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.scatter(x, y, s=36, alpha=0.85)
    ax.plot(lims, lims, "k--", linewidth=1.3, label="1:1")
    ax.axhline(true_rate_mm_kyr, color="tab:red", linewidth=1.5, label=f"True long-term rate ({true_rate_mm_kyr:g} mm/kyr)")
    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_xlabel("Apparent control erosion rate (mm/kyr)")
    ax.set_ylabel("Apparent landslide erosion rate (mm/kyr)")
    ax.set_title("Apparent vs true erosion rates (sample ensemble)")
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.6)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def chart_bias_ratio(sample_effect_path: Path, out_png: Path) -> None:
    df = pd.read_csv(sample_effect_path)
    for col in ("ratio_to_control_median", "ratio_to_control_p05", "ratio_to_control_p95"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "ratio_to_control_median" not in df.columns:
        raise RuntimeError("Sample-effect CSV is missing ratio_to_control_median.")

    med = df["ratio_to_control_median"].to_numpy(dtype=float)
    keep = np.isfinite(med)
    if not keep.any():
        raise RuntimeError("No finite ratio values in sample-effect CSV.")
    med = med[keep]
    names = df["sample_name"].astype(str).to_numpy()[keep] if "sample_name" in df.columns else np.array([str(i + 1) for i in range(med.size)])
    x = np.arange(med.size)

    has_ci = "ratio_to_control_p05" in df.columns and "ratio_to_control_p95" in df.columns
    fig, ax = plt.subplots(figsize=(max(7.0, 0.9 * med.size), 4.6))
    if has_ci:
        p05 = df["ratio_to_control_p05"].to_numpy(dtype=float)[keep]
        p95 = df["ratio_to_control_p95"].to_numpy(dtype=float)[keep]
        ok = np.isfinite(p05) & np.isfinite(p95)
        yerr_lo = np.where(ok, med - p05, 0.0)
        yerr_hi = np.where(ok, p95 - med, 0.0)
        ax.errorbar(x, med, yerr=[yerr_lo, yerr_hi], fmt="o", capsize=3)
    else:
        ax.plot(x, med, "o")

    ax.axhline(1.0, color="k", linestyle="--", linewidth=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Apparent rate ratio (landslide/control)")
    ax.set_title("Sample-wise landslide bias (Niemi-style ratio)")
    ax.grid(True, axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate landslide cosmo diagnostics.")
    parser.add_argument("--single-excavation", required=True)
    parser.add_argument("--single-inventory", required=True)
    parser.add_argument("--series-sample-effect", required=True)
    parser.add_argument("--out-outlines", required=True)
    parser.add_argument("--out-larsen-png", required=True)
    parser.add_argument("--out-apparent-vs-true-png", required=True)
    parser.add_argument("--out-ratio-png", required=True)
    parser.add_argument("--larsen-alpha", type=float, default=0.2344)
    parser.add_argument("--larsen-gamma", type=float, default=1.410)
    parser.add_argument("--true-rate-mm-kyr", type=float, default=50.0)
    args = parser.parse_args()

    outlines_from_excavation(Path(args.single_excavation), Path(args.out_outlines))
    chart_larsen_comparison(Path(args.single_inventory), Path(args.out_larsen_png), args.larsen_alpha, args.larsen_gamma)
    chart_apparent_vs_true(Path(args.series_sample_effect), Path(args.out_apparent_vs_true_png), args.true_rate_mm_kyr)
    chart_bias_ratio(Path(args.series_sample_effect), Path(args.out_ratio_png))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
