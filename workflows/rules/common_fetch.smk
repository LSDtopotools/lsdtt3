RUN_ID = config["run_id"]
ROOT = config["io"]["root"]


rule write_fetch_params:
    output:
        param=f"{ROOT}/{RUN_ID}/fetch/fetch_bbox.param"
    params:
        write_dir=f"{ROOT}/{RUN_ID}/fetch",
        dem_source=lambda wc: config["fetch"]["dem_source"],
        west=lambda wc: config["fetch"]["west"],
        south=lambda wc: config["fetch"]["south"],
        east=lambda wc: config["fetch"]["east"],
        north=lambda wc: config["fetch"]["north"],
        target_epsg=lambda wc: config["fetch"]["target_epsg"],
        grid_spacing=lambda wc: config["fetch"]["grid_spacing"],
        size_pixels=lambda wc: config["fetch"]["size_pixels"],
        http_connect_timeout_seconds=lambda wc: config["fetch"]["http_connect_timeout_seconds"],
        http_timeout_seconds=lambda wc: config["fetch"]["http_timeout_seconds"]
    run:
        import os
        from pathlib import Path

        out_path = Path(output.param)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_dir = os.path.abspath(params.write_dir)
        text = (
            f"write_path: {write_dir}\n"
            "write_fname: dem\n\n"
            f"west: {params.west}\n"
            f"south: {params.south}\n"
            f"east: {params.east}\n"
            f"north: {params.north}\n\n"
            f"dem_source: {params.dem_source}\n"
            f"target_epsg: {params.target_epsg}\n"
            f"grid_spacing: {params.grid_spacing}\n"
            f"size_pixels: {params.size_pixels}\n"
            f"http_connect_timeout_seconds: {params.http_connect_timeout_seconds}\n"
            f"http_timeout_seconds: {params.http_timeout_seconds}\n"
        )
        out_path.write_text(text, encoding="utf-8")


rule fetch_dem_bbox:
    input:
        param=f"{ROOT}/{RUN_ID}/fetch/fetch_bbox.param"
    output:
        dem=f"{ROOT}/{RUN_ID}/fetch/dem_DEM.tif"
    params:
        exe=lambda wc: config["executables"]["fetch_raster"],
        run_dir=f"{ROOT}/{RUN_ID}/fetch",
        param_name="fetch_bbox.param"
    message:
        "Fetch stage (bbox): lsdtt-fetch-raster -> {output.dem}"
    shell:
        "mkdir -p {params.run_dir} && {params.exe} {params.run_dir} {params.param_name}"
