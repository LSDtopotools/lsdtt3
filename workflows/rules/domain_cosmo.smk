RUN_ID = config["run_id"]
ROOT = config["io"]["root"]


rule write_cosmo_landslides_params:
    input:
        dem_fill=f"{ROOT}/{RUN_ID}/preprocess/dem_fill.tif"
    output:
        param=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/cosmo_landslides.param"
    params:
        read_dir=f"{ROOT}/{RUN_ID}/preprocess",
        write_dir=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides"
    run:
        import os
        from pathlib import Path

        cfg = config["cosmo_landslides"]
        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        read_dir = os.path.abspath(params.read_dir)
        write_dir = os.path.abspath(params.write_dir)

        text = (
            f"read_path: {read_dir}\n"
            "read_fname: dem_fill.tif\n"
            f"write_path: {write_dir}\n"
            "write_fname: ls\n\n"
            f"landslide_mode: {cfg['landslide_mode']}\n"
            f"threshold_contributing_pixels: {cfg['threshold_contributing_pixels']}\n"
            f"event_volume_m3: {cfg['event_volume_m3']}\n"
            f"event_time_yr: {cfg['event_time_yr']}\n"
            f"random_seed: {cfg['random_seed']}\n"
            "write_output_format: gtiff\n\n"
            f"calculate_crn_concentration: {'true' if cfg['calculate_crn_concentration'] else 'false'}\n"
            f"nuclide: {cfg['nuclide']}\n"
            f"muon_scheme: {cfg['muon_scheme']}\n"
            f"scaling_scheme: {cfg['scaling_scheme']}\n"
            f"background_eff_erate_g_cm2_yr: {cfg['background_eff_erate_g_cm2_yr']}\n"
            f"production_scaling: {cfg['production_scaling']}\n"
            f"total_shielding: {cfg['total_shielding']}\n"
            f"crn_report_time_yr: {cfg['crn_report_time_yr']}\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule run_cosmo_landslides:
    input:
        param=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/cosmo_landslides.param"
    output:
        inventory=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/ls_LS_INVENTORY.csv",
        tfail=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/ls_LS_TFAIL.tif",
        excavation=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/ls_LS_EXCAVATION.tif",
        conc=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/ls_LS_CONC.tif"
    params:
        exe=lambda wc: config["executables"]["cosmo_landslides"],
        run_dir=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides",
        param_name="cosmo_landslides.param"
    message:
        "Analysis stage: lsdtt-cosmo-landslides"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"
