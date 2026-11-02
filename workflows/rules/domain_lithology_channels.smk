# lithology_channels: channel steepness and channel-network statistics by
# lithologic unit, with the geology from Macrostrat.
#
# Per region (rules in common_regions.smk for the first three):
#   1. (optional) lsdtt-fetch-raster fetches the DEM;
#   2. lsdtt-fetch-lithology fetches the Macrostrat geologic map over the same
#      area, in the DEM's CRS;
#   3. lsdtt-basin-extraction picks the basins and writes their outlets;
#   4. lsdtt-chi-analysis: k_sn by lithologic unit in those basins
#      (per-unit summaries, histograms, bootstrap intervals);
#   5. lsdtt-channel-network-analysis: channel-link bearings, junction angles
#      and Strahler statistics grouped by lithologic unit;
#   6. when lsdviztools3 is found (see `figures`): the lithology map, k_sn by
#      unit, and bearing roses by unit.
#
# Configuration: config/lithology_channels.yaml. Outputs:
#   results/<run_id>/<region>/chi_lithology/  lsdtt-chi-analysis outputs
#   results/<run_id>/<region>/network/        lsdtt-channel-network-analysis outputs
#   results/<run_id>/<region>/figures/        PNGs

RUN_ID = config["run_id"]
ROOT = config["io"]["root"]


def lc_units(region):
    return f"{ROOT}/{RUN_ID}/{region}/lithology/{region}_macrostrat_units.fgb"


def lc_chi_text(wc):
    forced = {
        "basin_outlet_fname": workflow_path(rg_outlets(wc.region)),
        "basin_outlet_snapping_method": "none",
        "lithology_fname": workflow_path(lc_units(wc.region)),
    }
    user = rg_block(wc, "chi_analysis")
    if figures_enabled():
        user["write_lithology_code_raster"] = True   # the lithology map draws it
    return render_param_file(
        f"lithology_channels ({wc.region}, chi_analysis)",
        io_block(rg_dem(wc), f"{ROOT}/{RUN_ID}/{wc.region}/chi_lithology", wc.region),
        user,
        forced,
        comment=f"lithology_channels workflow, run {RUN_ID}: lsdtt-chi-analysis",
    )


def lc_network_text(wc):
    forced = {"lithology_fname": workflow_path(lc_units(wc.region))}
    return render_param_file(
        f"lithology_channels ({wc.region}, channel_network_analysis)",
        io_block(rg_dem(wc), f"{ROOT}/{RUN_ID}/{wc.region}/network", wc.region),
        rg_block(wc, "channel_network_analysis"),
        forced,
        comment=f"lithology_channels workflow, run {RUN_ID}: lsdtt-channel-network-analysis",
    )


rule lc_write_chi_params:
    input:
        dem=rg_dem,
        outlets=lambda wc: rg_outlets(wc.region),
        units=lambda wc: lc_units(wc.region)
    output:
        param=f"{RG_DIR}/chi_lithology/chi_analysis.param"
    params:
        text=lc_chi_text
    run:
        Path(output.param).write_text(params.text, encoding="utf-8")


rule lc_chi:
    input:
        param=f"{RG_DIR}/chi_lithology/chi_analysis.param",
        dem=rg_dem,
        outlets=lambda wc: rg_outlets(wc.region),
        units=lambda wc: lc_units(wc.region)
    output:
        readme=f"{RG_DIR}/chi_lithology/{{region}}_README.csv"
    log:
        f"{RG_DIR}/logs/chi_analysis.log"
    params:
        exe=lambda wc: lsdtt_exe("chi_analysis"),
        run_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/chi_lithology",
        param_name="chi_analysis.param"
    message:
        "lithology_channels [{wildcards.region}]: lsdtt-chi-analysis (k_sn by lithology)"
    shell:
        CLI_SHELL


rule lc_write_network_params:
    input:
        dem=rg_dem,
        units=lambda wc: lc_units(wc.region)
    output:
        param=f"{RG_DIR}/network/channel_network_analysis.param"
    params:
        text=lc_network_text
    run:
        Path(output.param).write_text(params.text, encoding="utf-8")


rule lc_network:
    input:
        param=f"{RG_DIR}/network/channel_network_analysis.param",
        dem=rg_dem,
        units=lambda wc: lc_units(wc.region)
    output:
        readme=f"{RG_DIR}/network/{{region}}_README.csv"
    log:
        f"{RG_DIR}/logs/channel_network_analysis.log"
    params:
        exe=lambda wc: lsdtt_exe("channel_network_analysis"),
        run_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/network",
        param_name="channel_network_analysis.param"
    message:
        "lithology_channels [{wildcards.region}]: lsdtt-channel-network-analysis (network statistics by lithology)"
    shell:
        CLI_SHELL


# lsdviz3-ksn-by-lithology and lsdviz3-channel-bearings write several PNGs
# whose names depend on the data, so each rule leaves a marker file.
rule lc_figures:
    input:
        chi=f"{RG_DIR}/chi_lithology/{{region}}_README.csv",
        network=f"{RG_DIR}/network/{{region}}_README.csv",
        dem=rg_dem
    output:
        lithology_map=f"{RG_DIR}/figures/{{region}}_lithology_map.png",
        done=f"{RG_DIR}/figures/{{region}}_figures.done"
    log:
        lithology_map=f"{RG_DIR}/logs/figures_lithology_map.log",
        ksn=f"{RG_DIR}/logs/figures_ksn_by_lithology.log",
        bearings=f"{RG_DIR}/logs/figures_channel_bearings.log"
    params:
        chi_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/chi_lithology",
        network_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/network",
        out_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/figures",
        mn=lambda wc: mn_suffix(wc, "chi_analysis"),
        # the roses are drawn for the stratification the network run used
        stratum=lambda wc: rg_block(wc, "channel_network_analysis").get(
            "channel_bearing_stratify_by", "none")
    message:
        "lithology_channels [{wildcards.region}]: figures (lsdviztools3)"
    run:
        import csv
        codes = f"{params.chi_dir}/{wildcards.region}_lithology_code.tif"
        # A map unit whose lithology field is empty comes out of lsdtt3 as a
        # unit with no name, which lsdviz3-lithology-map refuses to draw.
        names = f"{params.out_dir}/{wildcards.region}_lithology_names.csv"
        Path(params.out_dir).mkdir(parents=True, exist_ok=True)
        with open(f"{params.chi_dir}/{wildcards.region}_lithology_units.csv", newline="") as fin, \
                open(names, "w", newline="") as fout:
            writer = csv.writer(fout)
            writer.writerow(["code", "name"])
            for row in csv.DictReader(fin):
                writer.writerow([row["code"], row["name"] or "no lithology in map"])
        run_lsdviz3(log.lithology_map, "lsdviz3-lithology-map", "--codes", codes,
                    "--dem", input.dem, "--names", names, "--out", output.lithology_map)
        run_lsdviz3(log.ksn, "lsdviz3-ksn-by-lithology", "--data-dir", params.chi_dir,
                    "--prefix", wildcards.region, "--m-over-n", params.mn,
                    "--out-dir", params.out_dir, "--title", wildcards.region)
        run_lsdviz3(log.bearings, "lsdviz3-channel-bearings", "--data-dir", params.network_dir,
                    "--prefix", wildcards.region, "--out-dir", params.out_dir,
                    "--group-key", "lithology", "--stratum-key", params.stratum)
        Path(output.done).write_text("figures drawn\n", encoding="utf-8")


def lithology_channels_targets():
    if RG_SECTION != "lithology_channels":
        return []
    outs = []
    for r in RG_ITEMS:
        outs.append(f"{ROOT}/{RUN_ID}/{r}/chi_lithology/{r}_README.csv")
        outs.append(f"{ROOT}/{RUN_ID}/{r}/network/{r}_README.csv")
        if figures_enabled():
            outs.append(f"{ROOT}/{RUN_ID}/{r}/figures/{r}_figures.done")
    return outs


def lithology_channels_manifest():
    root = f"{ROOT}/{RUN_ID}"
    outs = rg_common_manifest()
    for r in RG_ITEMS:
        outs += [f"{root}/{r}/lithology/fetch_lithology.param", lc_units(r),
                 f"{root}/{r}/chi_lithology/chi_analysis.param",
                 f"{root}/{r}/network/channel_network_analysis.param"]
    return outs + [t for t in lithology_channels_targets() if not t.endswith(".done")]
