# Stages shared by the region modules (chi_region, lithology_channels): one
# DEM per region, fetched or read from a file, and the basins analysed in it.
#
#   results/<run_id>/<region>/dem/        lsdtt-fetch-raster output
#   results/<run_id>/<region>/basins/     lsdtt-basin-extraction output (the
#                                         region's own outlets, or basins picked
#                                         automatically)
#   results/<run_id>/<region>/lithology/  lsdtt-fetch-lithology output
#   results/<run_id>/<region>/logs/       one log per CLI run
#
# The active module's config section holds the regions, so a run reads
# config["chi_region"] or config["lithology_channels"].

RUN_ID = config["run_id"]
ROOT = config["io"]["root"]

REGION_MODULES = ("chi_region", "lithology_channels")
RG_SECTION = config["analysis"]["module"] if config["analysis"]["module"] in REGION_MODULES else None
RG = config.get(RG_SECTION, {}) if RG_SECTION else {}
RG_ITEMS = (items_from_config(RG_SECTION, "regions", "dem_dir", "dem_glob")
            if RG_SECTION else {})
RG_DIR = f"{ROOT}/{RUN_ID}/{{region}}"

wildcard_constraints:
    region=item_name_pattern(RG_ITEMS)


def rg_dem(wc):
    entry = RG_ITEMS[wc.region]
    if "fetch" in entry:
        return f"{ROOT}/{RUN_ID}/{wc.region}/dem/{wc.region}_DEM.tif"
    return item_file(RG_SECTION, "dem_dir", entry)


def rg_block(wc, block):
    """Shared parameters, then this CLI's block, then the region's overrides."""
    entry = RG_ITEMS[wc.region]
    return merged_cli_parameters(
        f"{RG_SECTION} ({wc.region}, {block})",
        RG.get("shared_parameters"),
        RG.get(block),
        entry.get("shared_parameters"),
        entry.get(block),
    )


def rg_outlets(region):
    return f"{ROOT}/{RUN_ID}/{region}/basins/{region}_selected_basin_outlets.csv"


def rg_fetch_text(wc):
    return fetch_param_text(
        f"{RG_SECTION}.regions.{wc.region}.fetch",
        RG_ITEMS[wc.region]["fetch"],
        f"{ROOT}/{RUN_ID}/{wc.region}/dem",
        wc.region,
    )


def rg_outlets_file(wc):
    """The region's own outlet points (CSV, FlatGeobuf or GeoPackage), or None.

    `basin_outlet_fname` under the region wins over one for the whole section.
    """
    entry = RG_ITEMS[wc.region]
    path = entry.get("basin_outlet_fname", RG.get("basin_outlet_fname"))
    return workflow_path(path) if path else None


def rg_basins_text(wc):
    """Basin selection: the supplied outlets if there are any, else the auto block.

    With outlets the `basin_extraction_from_outlets` block applies (snapping);
    without, `basin_extraction` (which basins to pick automatically).
    """
    outlets = rg_outlets_file(wc)
    forced = {"write_selected_basin_outlets_csv": True}
    if outlets:
        block = "basin_extraction_from_outlets"
        forced.update({"select_basins_from_outlets": True, "basin_outlet_fname": outlets})
    else:
        block = "basin_extraction"
    return render_param_file(
        f"{RG_SECTION} ({wc.region}, {block})",
        io_block(rg_dem(wc), f"{ROOT}/{RUN_ID}/{wc.region}/basins", wc.region),
        rg_block(wc, block),
        forced,
        comment=f"{RG_SECTION} workflow, run {RUN_ID}: lsdtt-basin-extraction",
    )


def rg_lithology_text(wc):
    """lsdtt-fetch-lithology over the same area, in the DEM's CRS.

    lsdtt reads a lithology polygon map onto the DEM grid without reprojecting
    it, so the map must come back in the DEM's EPSG. A fetched region passes
    its own box and target_epsg (0 = the UTM zone of the box centre in both
    CLIs); a region read from a file must state its EPSG as `epsg`.
    """
    entry = RG_ITEMS[wc.region]
    where = f"{RG_SECTION} ({wc.region}, fetch_lithology)"
    user = rg_block(wc, "fetch_lithology")
    if "fetch" in entry:
        fetch = entry["fetch"]
        if "bounds_fname" in fetch:
            area = {"bounds_fname": workflow_path(fetch["bounds_fname"])}
        else:
            area = {k: fetch[k] for k in ("west", "south", "east", "north")}
        epsg = fetch.get("target_epsg", 0)
    else:
        if "epsg" not in entry:
            raise ValueError(
                f"{where}: a region read from a file needs `epsg:` (the DEM's EPSG code) "
                f"so the lithology map is fetched in the same CRS")
        area = {"bounds_fname": workflow_path(rg_dem(wc))}
        epsg = entry["epsg"]
    io = {"write_path": workflow_path(f"{ROOT}/{RUN_ID}/{wc.region}/lithology") + "/",
          "write_fname": wc.region}
    return render_param_file(where, io, user, {**area, "output_epsg": epsg},
                             comment=f"{RG_SECTION} workflow, run {RUN_ID}: lsdtt-fetch-lithology")


rule rg_write_fetch_params:
    output:
        param=f"{RG_DIR}/dem/fetch.param"
    params:
        text=rg_fetch_text
    run:
        Path(output.param).write_text(params.text, encoding="utf-8")


rule rg_fetch_dem:
    input:
        param=f"{RG_DIR}/dem/fetch.param"
    output:
        dem=f"{RG_DIR}/dem/{{region}}_DEM.tif"
    log:
        f"{RG_DIR}/logs/fetch.log"
    params:
        exe=lambda wc: lsdtt_exe("fetch_raster"),
        run_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/dem",
        param_name="fetch.param"
    message:
        f"{RG_SECTION} [{{wildcards.region}}]: fetching the DEM"
    shell:
        CLI_SHELL


rule rg_write_basins_params:
    input:
        dem=rg_dem,
        outlets=lambda wc: [rg_outlets_file(wc)] if rg_outlets_file(wc) else []
    output:
        param=f"{RG_DIR}/basins/basin_extraction.param"
    params:
        text=rg_basins_text
    run:
        Path(output.param).write_text(params.text, encoding="utf-8")


rule rg_basins:
    input:
        param=f"{RG_DIR}/basins/basin_extraction.param",
        dem=rg_dem
    output:
        outlets=f"{RG_DIR}/basins/{{region}}_selected_basin_outlets.csv"
    log:
        f"{RG_DIR}/logs/basin_extraction.log"
    params:
        exe=lambda wc: lsdtt_exe("basin_extraction"),
        run_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/basins",
        param_name="basin_extraction.param"
    message:
        f"{RG_SECTION} [{{wildcards.region}}]: lsdtt-basin-extraction (basins to analyse)"
    shell:
        CLI_SHELL


rule rg_write_lithology_params:
    output:
        param=f"{RG_DIR}/lithology/fetch_lithology.param"
    params:
        text=rg_lithology_text
    run:
        Path(output.param).write_text(params.text, encoding="utf-8")


rule rg_fetch_lithology:
    input:
        param=f"{RG_DIR}/lithology/fetch_lithology.param",
        dem=rg_dem
    output:
        units=f"{RG_DIR}/lithology/{{region}}_macrostrat_units.fgb"
    log:
        f"{RG_DIR}/logs/fetch_lithology.log"
    params:
        exe=lambda wc: lsdtt_exe("fetch_lithology"),
        run_dir=lambda wc: f"{ROOT}/{RUN_ID}/{wc.region}/lithology",
        param_name="fetch_lithology.param"
    message:
        f"{RG_SECTION} [{{wildcards.region}}]: lsdtt-fetch-lithology (Macrostrat map)"
    shell:
        CLI_SHELL


def rg_common_manifest():
    root = f"{ROOT}/{RUN_ID}"
    outs = []
    for region, entry in RG_ITEMS.items():
        if "fetch" in entry:
            outs += [f"{root}/{region}/dem/fetch.param", f"{root}/{region}/dem/{region}_DEM.tif"]
        outs += [f"{root}/{region}/basins/basin_extraction.param", rg_outlets(region)]
    return outs
