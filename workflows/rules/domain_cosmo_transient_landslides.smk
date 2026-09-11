### Cosmogenic transient step-change + landslide ensemble.
###
### This is the Snakemake translation of the hand-run experiment documented in
### docs/plan_docs/plan_cosmo_transient_landslides.md and executed by
### p5_gen_ls.sh / p5_run_scenario.sh / p5_analyse.py in
### data/cosmo_transient_test/lsdtt3_version/ls_phase3/. The `.param` files in
### that directory are the ground truth: every key written below is copied from
### them verbatim, only the values that define a point in the ensemble
### (magnitude, seed, coupling, report time) are substituted.
###
### The ensemble is (magnitude) x (seed) x (coupling), expressed as wildcards:
###
###   mag   : a key of cosmo_transient_landslides.magnitudes, e.g. 5x, 10x
###   seed  : sNN for an ensemble member, or CTRL for the per-magnitude control
###   scen  : a key of cosmo_transient_landslides.scenarios (S1, S1f, S2, S3),
###           or the control scenario name (S0)
###
### Everything for one (mag, seed) lives in a single flat directory, exactly as
### the hand runs did. That is not cosmetic: `crn_nc_file` and the landslide
### `*_raster_prefix` parameters are resolved relative to `read_path`, never to
### the working directory, so the CRN netCDF and the landslide rasters must sit
### next to the DEM the run reads. `mask_fname`, `sample_points_filename` and
### `step_change_basin_outlet_csv` go through `add_read_path_if_needed` and may
### therefore be absolute; the others may not.

import os

RUN_ID = config["run_id"]
ROOT = config["io"]["root"]

# .get() throughout: the other workflow config files (e.g.
# config/cosmo_landslides_synthetic.yaml) replace config.yaml wholesale and do
# not carry this block, and every rule file is parsed for every module.
CTL = config.get("cosmo_transient_landslides", {})

CT_MAGNITUDES = list(CTL.get("magnitudes", {}).keys())
CT_SEEDS = [f"s{int(s):02d}" for s in CTL.get("seeds", [])]
CT_SCENARIOS = list(CTL.get("scenarios", {}).keys())
CT_CONTROL = CTL.get("control_scenario", {}).get("name", "S0")
CT_NUCLIDES = list(CTL.get("nuclides", ["Be10", "Al26", "C14"]))
CT_RATIOS = list(CTL.get("analysis", {}).get("ratios", []))

CT_SRC_DIR = os.path.abspath(os.path.expanduser(CTL.get("source_dir", "SOURCE_DIR")))
CT_SRC_DEM = CTL.get("source_dem_fname", "SOURCE_DEM")
CT_SRC_STEM = os.path.splitext(CT_SRC_DEM)[0]
CT_CRN_NC = CTL.get("crn_nc_file", "crn_parameters.nc")
CT_OUTLET_CSV = CTL.get("basin_outlet_csv", "basin_outlet.csv")

CT_RUNS = f"{ROOT}/{RUN_ID}/runs"
CT_PREP = f"{ROOT}/{RUN_ID}/prep"
CT_REPORTS = f"{ROOT}/{RUN_ID}/reports/cosmo_transient_landslides"


wildcard_constraints:
    mag=r"[0-9]+x",
    seed=r"(s[0-9]{2}|CTRL)",
    seed_num=r"[0-9]{2}",
    scen=r"S[0-9][a-z]?",
    nuc=r"(Be10|Al26|C14)",


# ---------------------------------------------------------------------------
# Helpers. These only assemble parameter text; no physics happens here.
# ---------------------------------------------------------------------------


def _ct_bool(value):
    return "true" if value else "false"


def _ct_coupling(scen):
    """(landslide_coupling, landslide_report_time_yr) for a scenario name."""
    if scen == CT_CONTROL:
        cfg = CTL["control_scenario"]
    else:
        cfg = CTL["scenarios"][scen]
    return cfg["coupling"], cfg["report_time_yr"]


def _ct_uplift_new(mag):
    return CTL["magnitudes"][mag]["uplift_new_m_yr"]


def _ct_event_volume_per_area_m(mag):
    """Geometric scar thickness for one event, in metres.

    plan_cosmo_transient_landslides.md section 5: the event is expressed as
    years of post-step erosion, converted from a ROCK-equivalent depth to a
    GEOMETRIC scar depth by rho_rock / rho_scar. The two densities are not
    interchangeable -- Larsen's V = alpha * A^gamma is fitted to measured scar
    volumes, so event_volume_per_area_m is a soil thickness while the erosion
    budget it balances is rock.

    Formatted to 5 decimal places because the hand runs used 0.17667 / 0.35333
    and the realised landslide field depends on the exact value.
    """
    ls = CTL["landslides"]
    override = ls.get("event_volume_per_area_m")
    if override is not None:
        return f"{float(override):.5f}"
    rho_rock = float(CTL["crn"]["rock_density_kg_m3"])
    rho_scar = float(ls["scar_density_kg_m3"])
    depth = _ct_uplift_new(mag) * float(ls["recurrence_yr"]) * (rho_rock / rho_scar)
    return f"{depth:.5f}"


def _ct_step_change_block(mag):
    sc = CTL["step_change"]
    return (
        "make_step_change_rasters: true\n"
        "write_step_change_rasters: true\n"
        f"step_change_production_from_generated_dem: {_ct_bool(sc['production_from_generated_dem'])}\n"
        # The generated step-change DEM carries a ~56 m elevation cliff along
        # drainage divides (plan section 7.2), so 5.4 % of tag-4 pixels drain
        # across it when D8 is rebuilt on that surface. Routing on the ORIGINAL
        # DEM removes the contamination entirely; the final runs used false.
        f"step_change_route_on_generated_dem: {_ct_bool(sc['route_on_generated_dem'])}\n"
        f"step_change_uplift_old_m_yr: {sc['uplift_old_m_yr']}\n"
        f"step_change_uplift_new_m_yr: {_ct_uplift_new(mag)}\n"
        f"step_change_age_yr: {sc['age_yr']}\n"
        f"step_change_target_relief_m: {sc['target_relief_m']}\n"
        f"step_change_m_over_n: {sc['m_over_n']}\n"
        f"step_change_n: {sc['n']}\n"
        f"step_change_A_0: {sc['A_0']}\n"
        f"step_change_channel_pixels: {sc['channel_pixels']}\n"
        f"step_change_Sc_old: {sc['Sc_old']}\n"
        f"step_change_Sc_new: {sc['Sc_new']}\n"
        f"step_change_basin_outlet_csv: {CT_OUTLET_CSV}\n"
    )


def _ct_crn_block():
    crn = CTL["crn"]
    return (
        # crn_nc_file is resolved as read_path + crn_nc_file, so it MUST be a
        # bare filename staged into read_path. Without it, atmosphere: ERA40
        # silently degrades to standard.
        f"crn_nc_file: {CT_CRN_NC}\n"
        f"scaling_scheme: {crn['scaling_scheme']}\n"
        f"muon_scheme: {crn['muon_scheme']}\n"
        f"atmosphere: {crn['atmosphere']}\n"
        f"rock_density_kg_m3: {crn['rock_density_kg_m3']}\n"
        f"sample_age_yr: {crn['sample_age_yr']}\n"
        f"topo_shield_value: {crn['topo_shield_value']}\n"
        f"snow_shield_eff_thickness: {crn['snow_shield_eff_thickness']}\n"
        f"self_shield_eff_thickness: {crn['self_shield_eff_thickness']}\n"
    )


def _ct_routing_block():
    rt = CTL["routing"]
    return (
        f"fill_before_routing: {_ct_bool(rt['fill_before_routing'])}\n"
        f"min_slope_for_fill: {rt['min_slope_for_fill']}\n"
        f"boundary_conditions: {rt['boundary_conditions']}\n"
        f"threshold_contributing_pixels: {rt['threshold_contributing_pixels']}\n"
        f"snap_search_radius_pixels: {rt['snap_search_radius_pixels']}\n"
    )


# ---------------------------------------------------------------------------
# Prep stage: generate the step-change landscape once per magnitude, then
# extract the channel-network points that become the sample set.
# ---------------------------------------------------------------------------


rule write_cosmo_transient_prep_params:
    output:
        param=f"{CT_PREP}/{{mag}}/sc_{{mag}}_gen.param"
    params:
        write_dir=f"{CT_PREP}/{{mag}}"
    run:
        from pathlib import Path

        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_dir = os.path.abspath(params.write_dir)

        # Copied from mega_sc<mag>_gen.param. A prediction mode must be active
        # or the CLI exits before writing anything, hence the Be10 concentration
        # raster; only the SC_DEM and SC_TAG products are consumed downstream.
        text = (
            f"read_path: {CT_SRC_DIR}/\n"
            f"write_path: {write_dir}/\n"
            f"read_fname: {CT_SRC_DEM}\n"
            f"write_fname: sc_{wildcards.mag}\n"
            "write_output_format: gtiff\n\n"
            + _ct_step_change_block(wildcards.mag)
            + "\n"
            "calculate_step_change_crn_concentration_raster: true\n"
            "write_crn_concentration_raster: true\n"
            "nuclide: Be10\n"
            + _ct_crn_block()
            + "\n"
            + _ct_routing_block()
        )
        out_path.write_text(text, encoding="utf-8")


rule run_cosmo_transient_prep:
    input:
        param=f"{CT_PREP}/{{mag}}/sc_{{mag}}_gen.param",
        dem=f"{CT_SRC_DIR}/{CT_SRC_DEM}",
        crn=f"{CT_SRC_DIR}/{CT_CRN_NC}",
        outlet=f"{CT_SRC_DIR}/{CT_OUTLET_CSV}"
    output:
        dem=f"{CT_PREP}/{{mag}}/sc_{{mag}}_SC_DEM.tif",
        tag=f"{CT_PREP}/{{mag}}/sc_{{mag}}_SC_TAG.tif",
        erate=f"{CT_PREP}/{{mag}}/sc_{{mag}}_SC_ERATE.tif",
        steptime=f"{CT_PREP}/{{mag}}/sc_{{mag}}_SC_STEPTIME.tif",
        chi=f"{CT_PREP}/{{mag}}/sc_{{mag}}_SC_CHI.tif"
    params:
        exe=lambda wc: config["executables"]["cosmo_prediction"],
        run_dir=f"{CT_PREP}/{{mag}}",
        param_name=lambda wc: f"sc_{wc.mag}_gen.param"
    log:
        f"{CT_PREP}/{{mag}}/sc_{{mag}}_gen.log"
    message:
        "Prep stage: generate the {wildcards.mag} step-change landscape"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name} > {log} 2>&1"


rule write_cosmo_transient_channel_params:
    input:
        dem=f"{CT_PREP}/{{mag}}/sc_{{mag}}_SC_DEM.tif"
    output:
        param=f"{CT_PREP}/{{mag}}/sc_{{mag}}_chan.param"
    params:
        work_dir=f"{CT_PREP}/{{mag}}"
    run:
        from pathlib import Path

        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        work_dir = os.path.abspath(params.work_dir)
        rt = CTL["routing"]

        text = (
            f"read_path: {work_dir}/\n"
            f"write_path: {work_dir}/\n"
            f"read_fname: sc_{wildcards.mag}_SC_DEM.tif\n"
            f"write_fname: sc_{wildcards.mag}_chan\n"
            "write_output_format: gtiff\n\n"
            f"fill_before_routing: {_ct_bool(rt['fill_before_routing'])}\n"
            f"min_slope_for_fill: {rt['min_slope_for_fill']}\n"
            f"boundary_conditions: {rt['boundary_conditions']}\n"
            f"threshold_contributing_pixels: {CTL['sampling']['threshold_contributing_pixels']}\n\n"
            "write_channel_network_points: true\n"
            "channel_network_points_output_format: csv\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule run_cosmo_transient_channel_points:
    input:
        param=f"{CT_PREP}/{{mag}}/sc_{{mag}}_chan.param",
        dem=f"{CT_PREP}/{{mag}}/sc_{{mag}}_SC_DEM.tif"
    output:
        points=f"{CT_PREP}/{{mag}}/sc_{{mag}}_chan_channel_network_points.csv"
    params:
        exe=lambda wc: config["executables"]["channel_extraction"],
        run_dir=f"{CT_PREP}/{{mag}}",
        param_name=lambda wc: f"sc_{wc.mag}_chan.param"
    log:
        f"{CT_PREP}/{{mag}}/sc_{{mag}}_chan.log"
    message:
        "Prep stage: channel-network sample points for {wildcards.mag}"
    shell:
        "{params.exe} {params.run_dir} {params.param_name} > {log} 2>&1"


# ---------------------------------------------------------------------------
# Staging. One flat directory per (magnitude, seed) holding everything the
# forward and inverse runs resolve against read_path.
# ---------------------------------------------------------------------------


rule stage_cosmo_transient_run_dir:
    input:
        dem=f"{CT_SRC_DIR}/{CT_SRC_DEM}",
        crn=f"{CT_SRC_DIR}/{CT_CRN_NC}",
        points=f"{CT_PREP}/{{mag}}/sc_{{mag}}_chan_channel_network_points.csv"
    output:
        dem=f"{CT_RUNS}/{{mag}}/{{seed}}/{CT_SRC_DEM}",
        crn=f"{CT_RUNS}/{{mag}}/{{seed}}/{CT_CRN_NC}"
    params:
        run_dir=f"{CT_RUNS}/{{mag}}/{{seed}}"
    message:
        "Stage: read_path inputs for {wildcards.mag}/{wildcards.seed}"
    run:
        import glob
        import shutil
        from pathlib import Path

        run_dir = Path(params.run_dir)
        run_dir.mkdir(parents=True, exist_ok=True)

        # The DEM may carry sidecars (.hdr for ENVI/bil, .prj, .aux.xml); the
        # reader needs all of them, so stage the whole stem.
        sources = sorted(glob.glob(os.path.join(CT_SRC_DIR, CT_SRC_STEM + ".*")))
        sources.append(os.path.join(CT_SRC_DIR, CT_CRN_NC))
        for src in sources:
            dst = run_dir / os.path.basename(src)
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            try:
                # Symlinks keep the tree small; NTFS mounts sometimes refuse
                # them, so fall back to a copy rather than failing the run.
                dst.symlink_to(os.path.abspath(src))
            except OSError:
                shutil.copyfile(src, dst)


# ---------------------------------------------------------------------------
# Landslide field. One realisation per (magnitude, seed), placed only on
# SC_TAG == 4 (post-knickpoint hillslopes).
# ---------------------------------------------------------------------------


rule write_cosmo_transient_landslide_params:
    input:
        dem=f"{CT_PREP}/{{mag}}/sc_{{mag}}_SC_DEM.tif"
    output:
        param=f"{CT_RUNS}/{{mag}}/s{{seed_num}}/ls.param"
    params:
        read_dir=f"{CT_PREP}/{{mag}}",
        write_dir=f"{CT_RUNS}/{{mag}}/s{{seed_num}}"
    run:
        from pathlib import Path

        ls = CTL["landslides"]
        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        read_dir = os.path.abspath(params.read_dir)
        write_dir = os.path.abspath(params.write_dir)

        # The mask is SC_TAG from a single magnitude for every run: the
        # knickpoint travels identically in 5x and 10x, so the tag rasters are
        # identical, and the hand runs used the 5x tags throughout.
        mask_mag = CTL["mask_magnitude"]
        mask_path = os.path.abspath(
            f"{CT_PREP}/{mask_mag}/sc_{mask_mag}_SC_TAG.tif")

        text = (
            f"read_path: {read_dir}/\n"
            f"write_path: {write_dir}/\n"
            f"read_fname: sc_{wildcards.mag}_SC_DEM.tif\n"
            "write_fname: ls\n"
            f"mask_fname: {mask_path}\n"
            # mask_fname alone treats ANY non-zero valid cell as active, which
            # would host landslides on tags 1-3 (unadjusted ground and
            # channels). mask_active_values narrows it to tag 4 only.
            f"mask_active_values: {ls['mask_active_values']}\n"
            f"landslide_mode: {ls['landslide_mode']}\n"
            f"event_time_yr: {ls['event_time_yr']}\n"
            f"event_volume_per_area_m: {_ct_event_volume_per_area_m(wildcards.mag)}\n"
            f"area_volume_preset: {ls['area_volume_preset']}\n"
            f"landslide_size_model: {ls['landslide_size_model']}\n"
            f"threshold_contributing_pixels: {ls['threshold_contributing_pixels']}\n"
            f"random_seed: {int(wildcards.seed_num)}\n"
            "write_excavation_depth_raster: true\n"
            "write_time_of_failure_raster: true\n"
            "write_landslide_csv: true\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule run_cosmo_transient_landslides:
    input:
        param=f"{CT_RUNS}/{{mag}}/s{{seed_num}}/ls.param",
        dem=f"{CT_PREP}/{{mag}}/sc_{{mag}}_SC_DEM.tif",
        mask=lambda wc: f"{CT_PREP}/{CTL['mask_magnitude']}/sc_{CTL['mask_magnitude']}_SC_TAG.tif"
    output:
        tfail=f"{CT_RUNS}/{{mag}}/s{{seed_num}}/ls_LS_TFAIL.tif",
        excavation=f"{CT_RUNS}/{{mag}}/s{{seed_num}}/ls_LS_EXCAVATION.tif",
        inventory=f"{CT_RUNS}/{{mag}}/s{{seed_num}}/ls_LS_INVENTORY.csv"
    params:
        exe=lambda wc: config["executables"]["cosmo_landslides"],
        run_dir=f"{CT_RUNS}/{{mag}}/s{{seed_num}}",
        param_name="ls.param"
    log:
        f"{CT_RUNS}/{{mag}}/s{{seed_num}}/ls.log"
    message:
        "Analysis stage: landslide field {wildcards.mag} seed s{wildcards.seed_num}"
    shell:
        "{params.exe} {params.run_dir} {params.param_name} > {log} 2>&1"


# ---------------------------------------------------------------------------
# Forward concentration. One run per (magnitude, seed, scenario); the control
# scenario runs once per magnitude under seed CTRL.
# ---------------------------------------------------------------------------


def _ct_conc_inputs(wc):
    base = f"{CT_RUNS}/{wc.mag}/{wc.seed}"
    ins = {
        "param": f"{base}/{wc.scen}_conc.param",
        "dem": f"{base}/{CT_SRC_DEM}",
        "crn": f"{base}/{CT_CRN_NC}",
        "points": f"{CT_PREP}/{wc.mag}/sc_{wc.mag}_chan_channel_network_points.csv",
        "outlet": f"{CT_SRC_DIR}/{CT_OUTLET_CSV}",
    }
    coupling, _ = _ct_coupling(wc.scen)
    if coupling != "none":
        ins["tfail"] = f"{base}/ls_LS_TFAIL.tif"
        ins["excavation"] = f"{base}/ls_LS_EXCAVATION.tif"
    return ins


rule write_cosmo_transient_conc_params:
    output:
        param=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc.param"
    params:
        run_dir=f"{CT_RUNS}/{{mag}}/{{seed}}"
    run:
        from pathlib import Path

        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        run_dir = os.path.abspath(params.run_dir)
        coupling, report_time = _ct_coupling(wildcards.scen)
        points = os.path.abspath(
            f"{CT_PREP}/{wildcards.mag}/sc_{wildcards.mag}_chan_channel_network_points.csv")

        # Written in the same key order as p5_run_scenario.sh so the generated
        # file can be diffed against the hand-run ground truth.
        text = (
            f"read_path: {run_dir}/\n"
            f"write_path: {run_dir}/\n"
            f"read_fname: {CT_SRC_DEM}\n"
            f"write_fname: {wildcards.scen}_conc\n"
            "write_output_format: gtiff\n"
            + _ct_step_change_block(wildcards.mag)
            + _ct_crn_block()
            + "calculate_multi_isotope_step_change: true\n"
            "calculate_accumulated_crn_concentration: true\n"
            "write_step_change_concentration_raster: true\n"
            "write_isotope_ratio_rasters: true\n"
            "write_accumulated_crn_raster: true\n"
            f"landslide_coupling: {coupling}\n"
        )
        if coupling != "none":
            # Both prefixes are resolved as read_path + prefix + ".tif", so they
            # must be bare stems in the run directory. Omitted for coupling
            # none, which never reads them.
            text += (
                "landslide_time_raster_prefix: ls_LS_TFAIL\n"
                "landslide_excavation_raster_prefix: ls_LS_EXCAVATION\n"
            )
        text += (
            f"landslide_report_time_yr: {report_time}\n"
            f"flush_time_threshold_yr: {CTL['coupling']['flush_time_threshold_yr']}\n"
            # scar_density_kg_m3 is the GEOMETRIC scar density (1500) and is not
            # interchangeable with rock_density_kg_m3 (2650) above.
            f"scar_density_kg_m3: {CTL['landslides']['scar_density_kg_m3']}\n"
            + _ct_routing_block()
            + f"sample_points_filename: {points}\n"
            f"sample_points_stride: {CTL['sampling']['sample_points_stride']}\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule run_cosmo_transient_conc:
    input:
        unpack(_ct_conc_inputs)
    output:
        sampled=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_sampled_crn_conc.csv",
        sampled_be10=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_sampled_crn_conc_Be10.csv",
        sampled_al26=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_sampled_crn_conc_Al26.csv",
        sampled_c14=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_sampled_crn_conc_C14.csv",
        sc_dem=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_SC_DEM.tif",
        sc_tag=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_SC_TAG.tif",
        sc_steptime=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_SC_STEPTIME.tif"
    params:
        exe=lambda wc: config["executables"]["cosmo_prediction"],
        run_dir=f"{CT_RUNS}/{{mag}}/{{seed}}",
        param_name=lambda wc: f"{wc.scen}_conc.param"
    log:
        f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc.log"
    message:
        "Analysis stage: forward CRN {wildcards.mag}/{wildcards.seed}/{wildcards.scen}"
    shell:
        "{params.exe} {params.run_dir} {params.param_name} > {log} 2>&1"


# ---------------------------------------------------------------------------
# Inversion. One per nuclide, reading the sampled concentrations back.
# ---------------------------------------------------------------------------


rule write_cosmo_transient_inversion_params:
    input:
        sampled=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_sampled_crn_conc_{{nuc}}.csv"
    output:
        param=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_inv_{{nuc}}.param"
    params:
        run_dir=f"{CT_RUNS}/{{mag}}/{{seed}}"
    run:
        from pathlib import Path

        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        run_dir = os.path.abspath(params.run_dir)

        text = (
            f"read_path: {run_dir}/\n"
            f"write_path: {run_dir}/\n"
            f"read_fname: {wildcards.scen}_conc_SC_DEM.tif\n"
            f"write_fname: {wildcards.scen}_inv_{wildcards.nuc}\n"
            f"nuclide: {wildcards.nuc}\n"
            f"points_filename: {wildcards.scen}_conc_sampled_crn_conc_{wildcards.nuc}.csv\n"
            f"concentration_column_name: {wildcards.nuc}_CONC\n"
            + _ct_crn_block()
            + f"self_consistent_age: {_ct_bool(CTL['inversion']['self_consistent_age'])}\n"
            + _ct_routing_block()
        )
        out_path.write_text(text, encoding="utf-8")


rule run_cosmo_transient_inversion:
    input:
        param=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_inv_{{nuc}}.param",
        sc_dem=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_SC_DEM.tif",
        sampled=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_sampled_crn_conc_{{nuc}}.csv",
        crn=f"{CT_RUNS}/{{mag}}/{{seed}}/{CT_CRN_NC}"
    output:
        erate=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_inv_{{nuc}}_CN_ERATE.csv"
    params:
        exe=lambda wc: config["executables"]["cosmo_tool"],
        run_dir=f"{CT_RUNS}/{{mag}}/{{seed}}",
        param_name=lambda wc: f"{wc.scen}_inv_{wc.nuc}.param"
    log:
        f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_inv_{{nuc}}.log"
    message:
        "Analysis stage: {wildcards.nuc} inversion {wildcards.mag}/{wildcards.seed}/{wildcards.scen}"
    shell:
        "{params.exe} {params.run_dir} {params.param_name} > {log} 2>&1"


# ---------------------------------------------------------------------------
# Join and ensemble statistics.
# ---------------------------------------------------------------------------


rule join_cosmo_transient_samples:
    input:
        param=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc.param",
        sampled=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_conc_sampled_crn_conc.csv",
        inversions=expand(
            f"{CT_RUNS}/{{{{mag}}}}/{{{{seed}}}}/{{{{scen}}}}_inv_{{nuc}}_CN_ERATE.csv",
            nuc=CT_NUCLIDES)
    output:
        joined=f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_joined_samples.csv"
    params:
        script="scripts/cosmo_transient_join_samples.py",
        run_dir=f"{CT_RUNS}/{{mag}}/{{seed}}",
        rel_tol=lambda wc: float(CTL["analysis"]["join_rel_tol"])
    log:
        f"{CT_RUNS}/{{mag}}/{{seed}}/{{scen}}_join.log"
    message:
        "Report stage: join samples {wildcards.mag}/{wildcards.seed}/{wildcards.scen}"
    shell:
        "python {params.script} "
        "--run-dir {params.run_dir} "
        "--stem {wildcards.scen} "
        "--rel-tol {params.rel_tol} "
        "--out {output.joined} > {log} 2>&1"


def _ct_ensemble_inputs(_wc):
    joined = []
    for mag in CT_MAGNITUDES:
        joined.append(f"{CT_RUNS}/{mag}/CTRL/{CT_CONTROL}_joined_samples.csv")
        for seed in CT_SEEDS:
            joined.append(f"{CT_RUNS}/{mag}/{seed}/ls_LS_INVENTORY.csv")
            for scen in CT_SCENARIOS:
                joined.append(f"{CT_RUNS}/{mag}/{seed}/{scen}_joined_samples.csv")
    return joined


rule cosmo_transient_ensemble:
    input:
        _ct_ensemble_inputs
    output:
        per_seed=f"{CT_REPORTS}/ensemble_per_seed.csv",
        summary=f"{CT_REPORTS}/ensemble_summary.csv",
        anchors=f"{CT_REPORTS}/ensemble_anchor_check.csv",
        report=f"{CT_REPORTS}/ensemble_report.txt"
    params:
        script="scripts/cosmo_transient_ensemble.py",
        runs_root=CT_RUNS,
        magnitudes=",".join(CT_MAGNITUDES),
        seeds=",".join(CT_SEEDS),
        scenarios=",".join(CT_SCENARIOS),
        control=CT_CONTROL,
        ratios=",".join(CT_RATIOS),
        anchor_eps=lambda wc: float(CTL["analysis"]["anchor_eps"]),
        anchor_alarm=lambda wc: float(CTL["analysis"]["anchor_alarm"]),
        expected_anchors=lambda wc: int(CTL["analysis"]["expected_anchor_count"]),
        equal_pairs=";".join(CTL["analysis"].get("assert_equal_scenarios", []))
    message:
        "Report stage: cosmogenic transient landslide ensemble statistics"
    shell:
        "python {params.script} "
        "--runs-root {params.runs_root} "
        "--magnitudes {params.magnitudes} "
        "--seeds {params.seeds} "
        "--scenarios {params.scenarios} "
        "--control {params.control} "
        "--ratios {params.ratios} "
        "--anchor-eps {params.anchor_eps} "
        "--anchor-alarm {params.anchor_alarm} "
        "--expected-anchors {params.expected_anchors} "
        "--assert-equal '{params.equal_pairs}' "
        "--out-per-seed {output.per_seed} "
        "--out-summary {output.summary} "
        "--out-anchors {output.anchors} "
        "--out-report {output.report}"
