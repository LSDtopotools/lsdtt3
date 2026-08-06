RUN_ID = config["run_id"]
ROOT = config["io"]["root"]


rule select_stream_order_for_minimum_basins:
    input:
        dem_fill=f"{ROOT}/{RUN_ID}/preprocess/dem_fill.tif"
    output:
        selected_order=f"{ROOT}/{RUN_ID}/extract/selected_stream_order.txt",
        report=f"{ROOT}/{RUN_ID}/extract/stream_order_scan.csv"
    params:
        exe=lambda wc: config["executables"]["basin_extraction"],
        min_basins=lambda wc: int(config["extract"]["minimum_basin_count"]),
        max_order=lambda wc: int(config["extract"]["max_search_stream_order"]),
        source_threshold=lambda wc: int(config["extract"]["source_pixels_threshold"]),
        work_dir=f"{ROOT}/{RUN_ID}/extract/order_probe",
        script="scripts/select_stream_order.py"
    message:
        "Extract stage: selecting largest stream order with at least {params.min_basins} basins"
    shell:
        "mkdir -p {params.work_dir} $(dirname {output.selected_order}) && "
        "python {params.script} "
        "--exe {params.exe} "
        "--dem {input.dem_fill} "
        "--work-dir {params.work_dir} "
        "--source-threshold {params.source_threshold} "
        "--min-basins {params.min_basins} "
        "--max-order {params.max_order} "
        "--selected-order-out {output.selected_order} "
        "--scan-report-out {output.report}"


rule write_basin_extraction_params:
    input:
        dem_fill=f"{ROOT}/{RUN_ID}/preprocess/dem_fill.tif",
        selected_order=f"{ROOT}/{RUN_ID}/extract/selected_stream_order.txt"
    output:
        param=f"{ROOT}/{RUN_ID}/extract/basin_extraction.param"
    params:
        read_dir=f"{ROOT}/{RUN_ID}/preprocess",
        write_dir=f"{ROOT}/{RUN_ID}/extract"
    run:
        import os
        from pathlib import Path

        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        read_dir = os.path.abspath(params.read_dir)
        write_dir = os.path.abspath(params.write_dir)
        selected_order = Path(input.selected_order).read_text(encoding="utf-8").strip()
        text = (
            f"read_path: {read_dir}\n"
            "read_fname: dem_fill.tif\n"
            f"write_path: {write_dir}\n"
            "write_fname: dem\n\n"
            "raster_is_filled: true\n"
            f"source_pixels_threshold: {config['extract']['source_pixels_threshold']}\n\n"
            "select_all_basins_of_stream_order: true\n"
            f"selected_stream_order: {selected_order}\n\n"
            "write_basin_outlines: true\n"
            "write_basin_outlines_as_polygon: false\n"
            "basin_outlines_output_format: fgb\n\n"
            "write_basin_outlet_points: false\n"
            "write_channel_network_lines: true\n"
            "channel_network_lines_output_format: fgb\n"
            "write_channel_network_points: false\n"
            "write_channel_network_junctions: true\n"
            "channel_network_junctions_output_format: fgb\n"
            "write_basins_raster: false\n"
            "calculate_spatial_index: false\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule extract_basins_and_channels:
    input:
        param=f"{ROOT}/{RUN_ID}/extract/basin_extraction.param"
    output:
        basins=f"{ROOT}/{RUN_ID}/extract/dem_junction_basin_outlines.fgb",
        channel_lines=f"{ROOT}/{RUN_ID}/extract/dem_channel_network_lines.fgb",
        channel_junctions=f"{ROOT}/{RUN_ID}/extract/dem_channel_network_junctions.fgb"
    params:
        exe=lambda wc: config["executables"]["basin_extraction"],
        run_dir=f"{ROOT}/{RUN_ID}/extract",
        param_name="basin_extraction.param"
    message:
        "Extract stage: basin outlines + channel network FGB outputs"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"
