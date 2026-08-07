RUN_ID = config["run_id"]
ROOT = config["io"]["root"]


rule write_synthetic_network_landslides_params:
    output:
        param=f"{ROOT}/{RUN_ID}/generate/ls_synthetic_network.param"
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
            "write_fname: ls_synth\n\n"
            f"n_rows: {cfg['n_rows']}\n"
            f"n_cols: {cfg['n_cols']}\n"
            f"cell_size: {cfg['cell_size']}\n"
            f"domain_type: {cfg['domain_type']}\n"
            f"outlet_mode: {cfg['outlet_mode']}\n"
            f"island_fill_fraction: {cfg['island_fill_fraction']}\n"
            f"energy_exponent_gamma: {cfg['energy_exponent_gamma']}\n"
            f"random_seed: {cfg['random_seed']}\n"
            "steady_state_mode: relief\n"
            f"target_relief: {config['cosmo_landslides_synthetic']['target_relief_m']}\n"
            f"m_exponent: {config['cosmo_landslides_synthetic']['m_exponent']}\n"
            f"n_exponent: {config['cosmo_landslides_synthetic']['n_exponent']}\n\n"
            "write_dem_raster: true\n"
            "write_drainage_area_raster: true\n"
            "verify_routing: true\n"
            "write_output_format: gtiff\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule generate_synthetic_landslide_dem:
    input:
        param=f"{ROOT}/{RUN_ID}/generate/ls_synthetic_network.param"
    output:
        dem=f"{ROOT}/{RUN_ID}/generate/ls_synth_ocn_dem.tif"
    params:
        exe=lambda wc: config["executables"]["synthetic_network"],
        run_dir=f"{ROOT}/{RUN_ID}/generate",
        param_name="ls_synthetic_network.param"
    message:
        "Generate stage: synthetic network DEM for landslide-cosmo workflow"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"


rule write_synthetic_landslide_preprocess_params:
    input:
        dem=f"{ROOT}/{RUN_ID}/generate/ls_synth_ocn_dem.tif"
    output:
        param=f"{ROOT}/{RUN_ID}/preprocess/ls_raster_preprocessing.param"
    params:
        read_dir=f"{ROOT}/{RUN_ID}/generate",
        write_dir=f"{ROOT}/{RUN_ID}/preprocess"
    run:
        import os
        from pathlib import Path

        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        read_dir = os.path.abspath(params.read_dir)
        write_dir = os.path.abspath(params.write_dir)
        text = (
            f"read_path: {read_dir}\n"
            "read_fname: ls_synth_ocn_dem.tif\n"
            f"write_path: {write_dir}\n"
            "write_fname: ls_synth\n\n"
            "enforce_projection: false\n"
            "write_output_format: gtiff\n"
            "write_fill_raster: true\n"
            "write_fill_diff_raster: false\n"
            "write_breached_raster: false\n"
            "write_hillshade_raster: true\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule preprocess_synthetic_landslide_dem:
    input:
        param=f"{ROOT}/{RUN_ID}/preprocess/ls_raster_preprocessing.param"
    output:
        dem_fill=f"{ROOT}/{RUN_ID}/preprocess/ls_synth_fill.tif",
        hillshade=f"{ROOT}/{RUN_ID}/preprocess/ls_synth_hillshade.tif"
    params:
        exe=lambda wc: config["executables"]["raster_preprocessing"],
        run_dir=f"{ROOT}/{RUN_ID}/preprocess",
        param_name="ls_raster_preprocessing.param"
    message:
        "Preprocess stage: fill + hillshade for synthetic landslide DEM"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"


rule write_landslide_single_event_params:
    input:
        dem_fill=f"{ROOT}/{RUN_ID}/preprocess/ls_synth_fill.tif"
    output:
        param=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/cosmo_landslides_single.param"
    params:
        read_dir=f"{ROOT}/{RUN_ID}/preprocess",
        write_dir=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single"
    run:
        import os
        from pathlib import Path

        cfg = config["cosmo_landslides_single"]
        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        read_dir = os.path.abspath(params.read_dir)
        write_dir = os.path.abspath(params.write_dir)
        text = (
            f"read_path: {read_dir}\n"
            "read_fname: ls_synth_fill.tif\n"
            f"write_path: {write_dir}\n"
            "write_fname: ls_single\n\n"
            "fill_before_routing: false\n"
            "write_output_format: gtiff\n\n"
            f"threshold_contributing_pixels: {cfg['threshold_contributing_pixels']}\n"
            "landslide_mode: single_event\n"
            f"event_time_yr: {cfg['event_time_yr']}\n"
            f"event_volume_m3: {cfg['event_volume_m3']}\n"
            f"random_seed: {cfg['random_seed']}\n"
            f"area_volume_preset: {cfg['area_volume_preset']}\n\n"
            f"calculate_crn_concentration: {'true' if cfg['calculate_crn_concentration'] else 'false'}\n"
            f"nuclide: {cfg['nuclide']}\n"
            f"muon_scheme: {cfg['muon_scheme']}\n"
            f"scaling_scheme: {cfg['scaling_scheme']}\n"
            f"background_eff_erate_g_cm2_yr: {cfg['background_eff_erate_g_cm2_yr']}\n"
            f"production_scaling: {cfg['production_scaling']}\n"
            f"total_shielding: {cfg['total_shielding']}\n"
            f"crn_report_time_yr: {cfg['crn_report_time_yr']}\n\n"
            "write_landslide_csv: true\n"
            "write_excavation_depth_raster: true\n"
            "write_time_of_failure_raster: true\n"
            "write_post_landslide_dem: true\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule run_landslide_single_event:
    input:
        param=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/cosmo_landslides_single.param"
    output:
        inventory=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_INVENTORY.csv",
        tfail=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_TFAIL.tif",
        excavation=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_EXCAVATION.tif",
        conc=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_CONC.tif",
        dem=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_DEM.tif"
    params:
        exe=lambda wc: config["executables"]["cosmo_landslides"],
        run_dir=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single",
        param_name="cosmo_landslides_single.param"
    message:
        "Analysis stage: single-event lsdtt-cosmo-landslides"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"


rule write_landslide_sample_points:
    input:
        dem_fill=f"{ROOT}/{RUN_ID}/preprocess/ls_synth_fill.tif"
    output:
        csv=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/sample_points.csv"
    params:
        n_points=lambda wc: int(config["cosmo_landslides_series"]["sample_points_count"])
    run:
        import numpy as np
        import rasterio
        from pathlib import Path
        from rasterio.warp import transform

        out_path = Path(output.csv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        n_points = max(1, int(params.n_points))

        with rasterio.open(input.dem_fill) as ds:
            arr = ds.read(1)
            mask = ds.read_masks(1) > 0
            valid = np.flatnonzero(mask & np.isfinite(arr))
            if valid.size < n_points:
                raise RuntimeError("Not enough valid pixels to generate sample points.")

            picks = []
            for i in range(n_points):
                frac = (i + 1) / (n_points + 1)
                j = int(round(frac * (valid.size - 1)))
                picks.append(int(valid[j]))

            rows = ["sample_name,latitude,longitude"]
            for i, flat in enumerate(picks, start=1):
                r, c = np.unravel_index(flat, arr.shape)
                x, y = ds.xy(int(r), int(c))
                lon, lat = transform(ds.crs, "EPSG:4326", [x], [y])
                rows.append(f"S{i},{lat[0]:.8f},{lon[0]:.8f}")

        out_path.write_text("\n".join(rows) + "\n", encoding="utf-8")


rule write_landslide_series_params:
    input:
        dem_fill=f"{ROOT}/{RUN_ID}/preprocess/ls_synth_fill.tif",
        samples=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/sample_points.csv"
    output:
        param=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/cosmo_landslides_series.param"
    params:
        read_dir=f"{ROOT}/{RUN_ID}/preprocess",
        write_dir=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series"
    run:
        import os
        from pathlib import Path

        cfg = config["cosmo_landslides_series"]
        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        read_dir = os.path.abspath(params.read_dir)
        write_dir = os.path.abspath(params.write_dir)
        text = (
            f"read_path: {read_dir}\n"
            "read_fname: ls_synth_fill.tif\n"
            f"write_path: {write_dir}\n"
            "write_fname: ls_series\n\n"
            "fill_before_routing: false\n"
            "write_output_format: gtiff\n\n"
            f"threshold_contributing_pixels: {cfg['threshold_contributing_pixels']}\n"
            "landslide_mode: series\n"
            f"series_start_yr: {cfg['series_start_yr']}\n"
            f"series_end_yr: {cfg['series_end_yr']}\n"
            f"mean_recurrence_yr: {cfg['mean_recurrence_yr']}\n"
            f"poisson_timing: {'true' if cfg['poisson_timing'] else 'false'}\n"
            f"long_term_erosion_rate_mm_kyr: {cfg['long_term_erosion_rate_mm_kyr']}\n"
            f"anneal_time_yr: {cfg['anneal_time_yr']}\n"
            f"random_seed: {cfg['random_seed']}\n"
            f"area_volume_preset: {cfg['area_volume_preset']}\n\n"
            f"calculate_crn_concentration: {'true' if cfg['calculate_crn_concentration'] else 'false'}\n"
            f"nuclide: {cfg['nuclide']}\n"
            f"muon_scheme: {cfg['muon_scheme']}\n"
            f"scaling_scheme: {cfg['scaling_scheme']}\n"
            f"background_eff_erate_g_cm2_yr: {cfg['background_eff_erate_g_cm2_yr']}\n"
            f"production_scaling: {cfg['production_scaling']}\n"
            f"total_shielding: {cfg['total_shielding']}\n"
            f"crn_report_time_yr: {cfg['crn_report_time_yr']}\n\n"
            "run_sample_point_ensemble: true\n"
            "sample_points_csv: sample_points.csv\n"
            "sample_name_column: sample_name\n"
            "sample_latitude_column: latitude\n"
            "sample_longitude_column: longitude\n"
            f"n_realisations: {cfg['n_realisations']}\n"
            f"ensemble_seed_stride: {cfg['ensemble_seed_stride']}\n"
            f"sample_snap_radius_pixels: {cfg['sample_snap_radius_pixels']}\n\n"
            "write_landslide_csv: true\n"
            "write_excavation_depth_raster: true\n"
            "write_time_of_failure_raster: true\n"
            "write_post_landslide_dem: true\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule run_landslide_series_event:
    input:
        param=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/cosmo_landslides_series.param"
    output:
        inventory=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/ls_series_LS_INVENTORY.csv",
        tfail=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/ls_series_LS_TFAIL.tif",
        excavation=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/ls_series_LS_EXCAVATION.tif",
        conc=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/ls_series_LS_CONC.tif",
        dem=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/ls_series_LS_DEM.tif",
        sample_effect=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/ls_series_LS_SAMPLE_EFFECT.csv"
    params:
        exe=lambda wc: config["executables"]["cosmo_landslides"],
        run_dir=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series",
        param_name="cosmo_landslides_series.param"
    message:
        "Analysis stage: series lsdtt-cosmo-landslides + sample ensemble"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"


rule postprocess_landslide_diagnostics:
    input:
        single_exc=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_EXCAVATION.tif",
        single_inv=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_INVENTORY.csv",
        sample_effect=f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/ls_series_LS_SAMPLE_EFFECT.csv"
    output:
        outlines=f"{ROOT}/{RUN_ID}/reports/cosmo_landslides/ls_single_landslide_outlines.fgb",
        larsen_png=f"{ROOT}/{RUN_ID}/reports/cosmo_landslides/ls_single_larsen_size_comparison.png",
        apparent_png=f"{ROOT}/{RUN_ID}/reports/cosmo_landslides/ls_series_apparent_vs_true_erosion.png",
        ratio_png=f"{ROOT}/{RUN_ID}/reports/cosmo_landslides/ls_series_sample_bias_ratio.png"
    params:
        script="scripts/landslide_cosmo_diagnostics.py",
        larsen_alpha=lambda wc: float(config["cosmo_landslides_diagnostics"]["larsen_alpha"]),
        larsen_gamma=lambda wc: float(config["cosmo_landslides_diagnostics"]["larsen_gamma"]),
        true_rate=lambda wc: float(config["cosmo_landslides_diagnostics"]["true_rate_mm_kyr"])
    message:
        "Report stage: outlines + landslide size + apparent-vs-true diagnostics"
    shell:
        "python {params.script} "
        "--single-excavation {input.single_exc} "
        "--single-inventory {input.single_inv} "
        "--series-sample-effect {input.sample_effect} "
        "--out-outlines {output.outlines} "
        "--out-larsen-png {output.larsen_png} "
        "--out-apparent-vs-true-png {output.apparent_png} "
        "--out-ratio-png {output.ratio_png} "
        "--larsen-alpha {params.larsen_alpha} "
        "--larsen-gamma {params.larsen_gamma} "
        "--true-rate-mm-kyr {params.true_rate}"
