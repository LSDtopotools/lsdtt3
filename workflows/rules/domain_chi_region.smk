# chi_region: chi analysis of the basins in a region.
#
# Per region (rules in common_regions.smk for the first two):
#   1. (optional) lsdtt-fetch-raster fetches the DEM, for regions with a
#      `fetch:` block instead of a `file:`;
#   2. lsdtt-basin-extraction picks the basins (by default every complete,
#      non-nested basin above a minimum area) and writes their outlets;
#   3. lsdtt-chi-analysis runs on those basins: chi, k_sn, m/n disorder and
#      the piecewise chi fits;
#   4. when lsdviztools3 is found (see `figures`), a map of the chi
#      coordinate draped on the shaded DEM.
#
# Configuration: config/chi_region.yaml. Outputs:
#   results/<run_id>/<region>/chi/       lsdtt-chi-analysis outputs
#   results/<run_id>/<region>/figures/   PNGs

RUN_ID = config["run_id"]
ROOT = config["io"]["root"]


def mn_suffix(wc, block):
    """The _mn_<m/n> suffix of chi outputs, from the first m_over_n given."""
    mn = rg_block(wc, block).get("m_over_n")
    if mn is None:
        raise ValueError(
            f"{RG_SECTION} ({wc.region}): figures need m_over_n in the {block} block "
            f"(the chi outputs are named after it); set it, or set figures.enabled: false")
    first = mn[0] if isinstance(mn, (list, tuple)) else str(mn).split(",")[0]
    return f"{float(first):.3f}".replace(".", "p")


def chi_region_text(wc):
    forced = {
        "basin_outlet_fname": workflow_path(rg_outlets(wc.region)),
        "basin_outlet_snapping_method": "none",
    }
    user = rg_block(wc, "chi_analysis")
    if figures_enabled():
        user["write_chi_coordinate_raster"] = True   # the chi map draws it
    return render_param_file(
        f"chi_region ({wc.region}, chi_analysis)",
        io_block(rg_dem(wc), f"{ROOT}/{RUN_ID}/{wc.region}/chi", wc.region),
        user,
        forced,
        comment=f"chi_region workflow, run {RUN_ID}: lsdtt-chi-analysis",
    )


rule chi_region_write_params:
    input:
        dem=rg_dem,
        outlets=lambda wc: rg_outlets(wc.region)
    output:
        param=f"{RG_DIR}/chi/chi_analysis.param"
    params:
        text=chi_region_text
    run:
        Path(output.param).write_text(params.text, encoding="utf-8")


rule chi_region_chi:
    input:
        param=f"{RG_DIR}/chi/chi_analysis.param",
        dem=rg_dem,
        outlets=lambda wc: rg_outlets(wc.region)
    output:
        readme=f"{RG_DIR}/chi/{{region}}_README.csv"
    log:
        f"{RG_DIR}/logs/chi_analysis.log"
    params:
        exe=lambda wc: lsdtt_exe("chi_analysis"),
        run_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/chi",
        param_name="chi_analysis.param"
    message:
        "chi_region [{wildcards.region}]: lsdtt-chi-analysis"
    shell:
        CLI_SHELL


rule chi_region_figures:
    input:
        readme=f"{RG_DIR}/chi/{{region}}_README.csv",
        dem=rg_dem
    output:
        chi=f"{RG_DIR}/figures/{{region}}_chi_coordinate.png"
    log:
        f"{RG_DIR}/logs/figures_chi.log"
    params:
        chi=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/chi/{wc.region}_chi_coordinate_mn_{mn_suffix(wc, 'chi_analysis')}.tif"
    message:
        "chi_region [{wildcards.region}]: figures (lsdviztools3)"
    run:
        run_lsdviz3(log[0], "lsdviz3-drape", "--base", input.dem, "--shade-base",
                    "--drape", params.chi, "--out", output.chi, "--cmap", "viridis",
                    "--colorbar", "--colorbar-label", "chi (m)")


def chi_region_targets():
    if RG_SECTION != "chi_region":
        return []
    outs = [f"{ROOT}/{RUN_ID}/{r}/chi/{r}_README.csv" for r in RG_ITEMS]
    if figures_enabled():
        outs += [f"{ROOT}/{RUN_ID}/{r}/figures/{r}_chi_coordinate.png" for r in RG_ITEMS]
    return outs


def chi_region_manifest():
    return rg_common_manifest() + [
        f"{ROOT}/{RUN_ID}/{r}/chi/chi_analysis.param" for r in RG_ITEMS] + [
        t for t in chi_region_targets() if not t.endswith(".png")]
