RUN_ID = config["run_id"]
ROOT = config["io"]["root"]


rule write_preprocess_params:
    input:
        dem=f"{ROOT}/{RUN_ID}/fetch/dem_DEM.tif"
    output:
        param=f"{ROOT}/{RUN_ID}/preprocess/raster_preprocessing.param"
    params:
        read_dir=f"{ROOT}/{RUN_ID}/fetch",
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
            "read_fname: dem_DEM.tif\n"
            f"write_path: {write_dir}\n"
            "write_fname: dem\n\n"
            "enforce_projection: false\n"
            "write_output_format: gtiff\n\n"
            "write_fill_raster: true\n"
            "write_fill_diff_raster: false\n"
            "write_breached_raster: false\n"
            "write_hillshade_raster: false\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule preprocess_dem:
    input:
        param=f"{ROOT}/{RUN_ID}/preprocess/raster_preprocessing.param"
    output:
        dem_fill=f"{ROOT}/{RUN_ID}/preprocess/dem_fill.tif"
    params:
        exe=lambda wc: config["executables"]["raster_preprocessing"],
        run_dir=f"{ROOT}/{RUN_ID}/preprocess",
        param_name="raster_preprocessing.param"
    message:
        "Preprocess stage: lsdtt-raster-preprocessing -> {output.dem_fill}"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"
