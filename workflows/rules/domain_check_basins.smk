RUN_ID = config["run_id"]
ROOT = config["io"]["root"]


rule write_check_basins_params:
    input:
        dem_fill=f"{ROOT}/{RUN_ID}/preprocess/dem_fill.tif"
    output:
        param=f"{ROOT}/{RUN_ID}/extract/check_basins.param"
    params:
        read_dir=f"{ROOT}/{RUN_ID}/preprocess",
        write_dir=f"{ROOT}/{RUN_ID}/extract"
    run:
        import os
        from pathlib import Path

        cfg = config["check_basins"]
        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        read_dir = os.path.abspath(params.read_dir)
        write_dir = os.path.abspath(params.write_dir)
        outlets_file = os.path.abspath(cfg["outlets_file"])

        text = (
            f"read_path: {read_dir}\n"
            "read_fname: dem_fill.tif\n"
            f"write_path: {write_dir}\n"
            "write_fname: check\n\n"
            "raster_is_filled: true\n"
            f"threshold_contributing_pixels: {cfg['threshold_contributing_pixels']}\n\n"
            "select_basins_from_outlets: true\n"
            f"basin_outlet_fname: {outlets_file}\n"
            f"basin_outlet_snapping_method: {cfg['basin_outlet_snapping_method']}\n"
            f"basin_outlet_contributing_pixels_threshold: {cfg['basin_outlet_contributing_pixels_threshold']}\n"
            f"basin_outlet_snap_min_search_radius_pixels: {cfg['basin_outlet_snap_min_search_radius_pixels']}\n"
            f"basin_outlet_snap_max_search_radius_pixels: {cfg['basin_outlet_snap_max_search_radius_pixels']}\n"
            f"write_basin_outlet_points: {'true' if cfg['write_basin_outlet_points'] else 'false'}\n\n"
            "write_basin_outlines: true\n"
            "write_basin_outlines_as_polygon: false\n"
            "basin_outlines_write_vector_format: fgb\n"
            "write_channel_network_lines: true\n"
            "write_vector_format: fgb\n"
            "write_channel_network_points: false\n"
            "write_channel_network_junctions: false\n"
            "write_basins_raster: false\n"
            "calculate_spatial_index: false\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule run_check_basins_extract:
    input:
        param=f"{ROOT}/{RUN_ID}/extract/check_basins.param"
    output:
        basins=f"{ROOT}/{RUN_ID}/extract/check_junction_basin_outlines.fgb",
        channel_lines=f"{ROOT}/{RUN_ID}/extract/check_channel_network_lines.fgb"
    params:
        exe=lambda wc: config["executables"]["basin_extraction"],
        run_dir=f"{ROOT}/{RUN_ID}/extract",
        param_name="check_basins.param"
    message:
        "Extract stage: outlet-check basin + channel FGB outputs"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"
