RUN_ID = config["run_id"]
ROOT = config["io"]["root"]


rule write_synthetic_network_params:
    output:
        param=f"{ROOT}/{RUN_ID}/generate/synthetic_network.param"
    params:
        write_dir=f"{ROOT}/{RUN_ID}/generate"
    run:
        import os
        from pathlib import Path

        cfg = config["synthetic_network"]
        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_dir = os.path.abspath(params.write_dir)

        text = (
            f"write_path: {write_dir}\n"
            "write_fname: rn\n\n"
            f"n_rows: {cfg['n_rows']}\n"
            f"n_cols: {cfg['n_cols']}\n"
            f"cell_size_metres: {cfg['cell_size']}\n"
            f"domain_type: {cfg['domain_type']}\n"
            f"outlet_mode: {cfg['outlet_mode']}\n"
            f"island_fill_fraction: {cfg['island_fill_fraction']}\n"
            f"energy_exponent_gamma: {cfg['energy_exponent_gamma']}\n"
            f"random_seed: {cfg['random_seed']}\n\n"
            "write_dem_raster: true\n"
            "write_drainage_area_raster: true\n"
            "verify_routing: true\n"
            "write_output_format: gtiff\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule generate_random_network_dem:
    input:
        param=f"{ROOT}/{RUN_ID}/generate/synthetic_network.param"
    output:
        dem=f"{ROOT}/{RUN_ID}/generate/rn_ocn_dem.tif"
    params:
        exe=lambda wc: config["executables"]["synthetic_network"],
        run_dir=f"{ROOT}/{RUN_ID}/generate",
        param_name="synthetic_network.param"
    message:
        "Generate stage: lsdtt-synthetic-networks -> {output.dem}"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"


rule write_cosmo_step_change_params:
    input:
        dem=f"{ROOT}/{RUN_ID}/generate/rn_ocn_dem.tif"
    output:
        param=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/cosmo_step_change.param"
    params:
        read_dir=f"{ROOT}/{RUN_ID}/generate",
        write_dir=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change"
    run:
        import os
        from pathlib import Path

        cfg = config["cosmo_step_change"]
        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        read_dir = os.path.abspath(params.read_dir)
        write_dir = os.path.abspath(params.write_dir)

        text = (
            f"read_path: {read_dir}\n"
            "read_fname: rn_ocn_dem.tif\n"
            f"write_path: {write_dir}\n"
            "write_fname: sc\n\n"
            "make_step_change_rasters: true\n"
            "write_step_change_rasters: true\n"
            f"step_change_uplift_old_m_yr: {cfg['step_change_uplift_old_m_yr']}\n"
            f"step_change_uplift_new_m_yr: {cfg['step_change_uplift_new_m_yr']}\n"
            f"step_change_target_relief_metres: {cfg['step_change_target_relief_metres']}\n"
            f"step_change_age_yr: {cfg['step_change_age_yr']}\n"
            f"step_change_n: {cfg['step_change_n']}\n"
            f"step_change_channel_pixels: {cfg['step_change_channel_pixels']}\n\n"
            "calculate_multi_isotope_step_change: true\n"
            f"muon_scheme: {cfg['muon_scheme']}\n"
            f"scaling_scheme: {cfg['scaling_scheme']}\n"
            f"atmosphere: {cfg['atmosphere']}\n"
            f"topo_shield_value: {cfg['topo_shield_value']}\n"
            "write_output_format: gtiff\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule run_cosmo_step_change_random_network:
    input:
        param=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/cosmo_step_change.param"
    output:
        be10=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_Be10_step_crn_conc.tif",
        al26=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_Al26_step_crn_conc.tif",
        c14=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_C14_step_crn_conc.tif",
        r2610=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_R2610_step.tif",
        r1410=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_R1410_step.tif",
        steptime=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_SC_STEPTIME.tif"
    params:
        exe=lambda wc: config["executables"]["cosmo_prediction"],
        run_dir=f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change",
        param_name="cosmo_step_change.param"
    message:
        "Analysis stage: lsdtt-cosmo-prediction multi-isotope step-change on random-network DEM"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"
