# valley_widths: valley-floor (floodplain) widths over a channel network, for
# one DEM or a folder of them.
#
# Per DEM:
#   1. (optional) fetch the DEM with lsdtt-fetch-raster, for items with a
#      `fetch:` block instead of a `file:`;
#   2. lsdtt-fluvial-surfaces finds the floodplain (or valley) mask;
#   3. lsdtt-valley-centreline, run_valley_width_network, traces a valley
#      centreline down every channel segment and casts cross sections across
#      it against that mask, giving widths, centrelines and sections;
#   4. lsdviztools3 maps and profiles, when lsdviztools3 is found.
# Then every DEM's width table is gathered into one table and a summary.
#
# Configuration: config/valley_widths.yaml. Outputs:
#   results/<run_id>/<dem>/floodplain/   lsdtt-fluvial-surfaces outputs
#   results/<run_id>/<dem>/widths/       lsdtt-valley-centreline outputs
#   results/<run_id>/<dem>/figures/      PNGs, when lsdviztools3 is found (see `figures`)
#   results/<run_id>/<dem>/logs/         one log per CLI run
#   results/<run_id>/reports/            all_valley_widths.csv, valley_width_summary.csv

RUN_ID = config["run_id"]
ROOT = config["io"]["root"]

VW = config.get("valley_widths", {})
VW_ITEMS = (items_from_config("valley_widths", "dems", "dem_dir", "dem_glob")
            if config["analysis"]["module"] == "valley_widths" else {})
VW_DIR = f"{ROOT}/{RUN_ID}/{{dem}}"

wildcard_constraints:
    dem=item_name_pattern(VW_ITEMS)
VW_MASK_SUFFIX = {"floodplain": "floodplain_mask", "valley": "valley_mask"}


def vw_mask_kind():
    kind = VW.get("width_mask", "floodplain")
    if kind not in VW_MASK_SUFFIX:
        raise ValueError("valley_widths.width_mask must be floodplain or valley")
    return kind


def vw_dem(wc):
    entry = VW_ITEMS[wc.dem]
    if "fetch" in entry:
        return f"{ROOT}/{RUN_ID}/{wc.dem}/dem/{wc.dem}_DEM.tif"
    return item_file("valley_widths", "dem_dir", entry)


def vw_block(wc, block):
    """Shared parameters, then this CLI's block, then the DEM's own overrides."""
    entry = VW_ITEMS[wc.dem]
    return merged_cli_parameters(
        f"valley_widths ({wc.dem}, {block})",
        VW.get("shared_parameters"),
        VW.get(block),
        entry.get("shared_parameters"),
        entry.get(block),
    )


def vw_fetch_text(wc):
    return fetch_param_text(
        f"valley_widths.dems.{wc.dem}.fetch",
        VW_ITEMS[wc.dem]["fetch"],
        f"{ROOT}/{RUN_ID}/{wc.dem}/dem",
        wc.dem,
    )


def vw_floodplain_text(wc):
    forced = {f"write_{VW_MASK_SUFFIX[vw_mask_kind()]}_raster": True}
    return render_param_file(
        f"valley_widths ({wc.dem}, fluvial_surfaces)",
        io_block(vw_dem(wc), f"{ROOT}/{RUN_ID}/{wc.dem}/floodplain", wc.dem),
        vw_block(wc, "fluvial_surfaces"),
        forced,
        comment=f"valley_widths workflow, run {RUN_ID}: lsdtt-fluvial-surfaces",
    )


def vw_mask_path(dem):
    return f"{ROOT}/{RUN_ID}/{dem}/floodplain/{dem}_{VW_MASK_SUFFIX[vw_mask_kind()]}.tif"


def vw_widths_text(wc):
    forced = {
        "valley_mask_fname": workflow_path(vw_mask_path(wc.dem)),
        "run_valley_width_network": True,
    }
    if figures_enabled():
        forced["write_hillshade_raster"] = True
    return render_param_file(
        f"valley_widths ({wc.dem}, valley_centreline)",
        io_block(vw_dem(wc), f"{ROOT}/{RUN_ID}/{wc.dem}/widths", wc.dem),
        vw_block(wc, "valley_centreline"),
        forced,
        comment=f"valley_widths workflow, run {RUN_ID}: lsdtt-valley-centreline",
    )


rule vw_write_fetch_params:
    output:
        param=f"{VW_DIR}/dem/fetch.param"
    params:
        text=vw_fetch_text
    run:
        Path(output.param).write_text(params.text, encoding="utf-8")


rule vw_fetch_dem:
    input:
        param=f"{VW_DIR}/dem/fetch.param"
    output:
        dem=f"{VW_DIR}/dem/{{dem}}_DEM.tif"
    log:
        f"{VW_DIR}/logs/fetch.log"
    params:
        exe=lambda wc: lsdtt_exe("fetch_raster"),
        run_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.dem}/dem",
        param_name="fetch.param"
    message:
        "valley_widths [{wildcards.dem}]: fetching the DEM"
    shell:
        CLI_SHELL


rule vw_write_floodplain_params:
    input:
        dem=vw_dem
    output:
        param=f"{VW_DIR}/floodplain/fluvial_surfaces.param"
    params:
        text=vw_floodplain_text
    run:
        Path(output.param).write_text(params.text, encoding="utf-8")


rule vw_floodplain:
    input:
        param=f"{VW_DIR}/floodplain/fluvial_surfaces.param",
        dem=vw_dem
    output:
        mask=f"{VW_DIR}/floodplain/{{dem}}_{VW_MASK_SUFFIX.get(VW.get('width_mask', 'floodplain'), 'floodplain_mask')}.tif",
        readme=f"{VW_DIR}/floodplain/{{dem}}_README.csv"
    log:
        f"{VW_DIR}/logs/fluvial_surfaces.log"
    params:
        exe=lambda wc: lsdtt_exe("fluvial_surfaces"),
        run_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.dem}/floodplain",
        param_name="fluvial_surfaces.param"
    message:
        "valley_widths [{wildcards.dem}]: lsdtt-fluvial-surfaces (floodplain mask)"
    shell:
        CLI_SHELL


rule vw_write_widths_params:
    input:
        dem=vw_dem,
        mask=lambda wc: vw_mask_path(wc.dem)
    output:
        param=f"{VW_DIR}/widths/valley_centreline.param"
    params:
        text=vw_widths_text
    run:
        Path(output.param).write_text(params.text, encoding="utf-8")


rule vw_widths:
    input:
        param=f"{VW_DIR}/widths/valley_centreline.param",
        dem=vw_dem,
        mask=lambda wc: vw_mask_path(wc.dem)
    output:
        widths=f"{VW_DIR}/widths/{{dem}}_valley_widths.csv",
        segments=f"{VW_DIR}/widths/{{dem}}_valley_segments.csv",
        readme=f"{VW_DIR}/widths/{{dem}}_README.csv"
    log:
        f"{VW_DIR}/logs/valley_centreline.log"
    params:
        exe=lambda wc: lsdtt_exe("valley_centreline"),
        run_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.dem}/widths",
        param_name="valley_centreline.param"
    message:
        "valley_widths [{wildcards.dem}]: lsdtt-valley-centreline (valley width network)"
    shell:
        CLI_SHELL


rule vw_figures:
    input:
        widths=f"{VW_DIR}/widths/{{dem}}_valley_widths.csv",
        mask=lambda wc: vw_mask_path(wc.dem)
    output:
        map=f"{VW_DIR}/figures/{{dem}}_valley_widths_map.png",
        profile=f"{VW_DIR}/figures/{{dem}}_valley_widths_profile.png",
        channels=f"{VW_DIR}/figures/{{dem}}_valley_widths_channels.png",
        floodplain=f"{VW_DIR}/figures/{{dem}}_floodplain_map.png",
        qq=f"{VW_DIR}/figures/{{dem}}_floodplain_hand_qq.png"
    log:
        widths=f"{VW_DIR}/logs/figures_valley_widths.log",
        floodplain=f"{VW_DIR}/logs/figures_floodplain.log"
    params:
        widths_stem=lambda wc: workflow_path(f"{ROOT}/{RUN_ID}/{wc.dem}/widths/{wc.dem}"),
        floodplain_stem=lambda wc: workflow_path(f"{ROOT}/{RUN_ID}/{wc.dem}/floodplain/{wc.dem}"),
        hillshade=lambda wc: workflow_path(f"{ROOT}/{RUN_ID}/{wc.dem}/widths/{wc.dem}_hillshade.tif")
    message:
        "valley_widths [{wildcards.dem}]: figures (lsdviztools3)"
    run:
        run_lsdviz3(log.widths, "lsdviz3-valley-widths",
                    "--stem", params.widths_stem, "--base", params.hillshade,
                    "--valley-mask", input.mask, "--title", wildcards.dem,
                    "--out-map", output.map, "--out-profile", output.profile,
                    "--out-channels", output.channels)
        # The floodplain polygons over the hillshade, and the HAND QQ plots
        # its bounds were read from (clubb2017 / qq_gaussian runs).
        run_lsdviz3(log.floodplain, "lsdviz3-floodplain",
                    "--stem", params.floodplain_stem, "--base", params.hillshade,
                    "--out-map", output.floodplain, "--out-qq", output.qq)


rule vw_collect:
    input:
        widths=expand(f"{ROOT}/{RUN_ID}/{{dem}}/widths/{{dem}}_valley_widths.csv", dem=list(VW_ITEMS))
    output:
        table=f"{ROOT}/{RUN_ID}/reports/all_valley_widths.csv",
        summary=f"{ROOT}/{RUN_ID}/reports/valley_width_summary.csv"
    params:
        names=list(VW_ITEMS)
    message:
        "valley_widths: gathering the width tables"
    shell:
        "python scripts/summarise_valley_widths.py --names {params.names} "
        "--widths {input.widths} --table-out {output.table} --summary-out {output.summary}"


def valley_widths_targets():
    root = f"{ROOT}/{RUN_ID}"
    outs = [f"{root}/reports/all_valley_widths.csv",
            f"{root}/reports/valley_width_summary.csv"]
    for dem in VW_ITEMS:
        outs.append(f"{root}/{dem}/widths/{dem}_valley_segments.csv")
        if figures_enabled():
            outs.append(f"{root}/{dem}/figures/{dem}_valley_widths_map.png")
    return outs


def valley_widths_manifest():
    root = f"{ROOT}/{RUN_ID}"
    outs = []
    for dem in VW_ITEMS:
        outs += [f"{root}/{dem}/floodplain/fluvial_surfaces.param",
                 vw_mask_path(dem),
                 f"{root}/{dem}/widths/valley_centreline.param",
                 f"{root}/{dem}/widths/{dem}_valley_widths.csv",
                 f"{root}/{dem}/widths/{dem}_valley_segments.csv"]
    return outs + valley_widths_targets()[:2]
