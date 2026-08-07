RUN_ID = config["run_id"]
ROOT = config["io"]["root"]
MODULE = config["analysis"]["module"]


def manifest_inputs(_wc):
    if MODULE == "cosmo_landslides":
        return [
            f"{ROOT}/{RUN_ID}/fetch/dem_DEM.tif",
            f"{ROOT}/{RUN_ID}/fetch/fetch_bbox.param",
            f"{ROOT}/{RUN_ID}/preprocess/raster_preprocessing.param",
            f"{ROOT}/{RUN_ID}/preprocess/dem_fill.tif",
            f"{ROOT}/{RUN_ID}/extract/selected_stream_order.txt",
            f"{ROOT}/{RUN_ID}/extract/stream_order_scan.csv",
            f"{ROOT}/{RUN_ID}/extract/basin_extraction.param",
            f"{ROOT}/{RUN_ID}/extract/dem_junction_basin_outlines.fgb",
            f"{ROOT}/{RUN_ID}/extract/dem_channel_network_lines.fgb",
            f"{ROOT}/{RUN_ID}/extract/dem_channel_network_junctions.fgb",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/cosmo_landslides.param",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/ls_LS_INVENTORY.csv",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/ls_LS_TFAIL.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/ls_LS_EXCAVATION.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides/ls_LS_CONC.tif",
        ]
    if MODULE == "cosmo_step_change_random_network":
        return [
            f"{ROOT}/{RUN_ID}/generate/synthetic_network.param",
            f"{ROOT}/{RUN_ID}/generate/rn_ocn_dem.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/cosmo_step_change.param",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_Be10_step_crn_conc.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_Al26_step_crn_conc.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_C14_step_crn_conc.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_R2610_step.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_R1410_step.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_step_change/sc_SC_STEPTIME.tif",
        ]
    if MODULE == "cosmo_landslides_synthetic":
        return [
            f"{ROOT}/{RUN_ID}/generate/ls_synthetic_network.param",
            f"{ROOT}/{RUN_ID}/generate/ls_synth_ocn_dem.tif",
            f"{ROOT}/{RUN_ID}/preprocess/ls_raster_preprocessing.param",
            f"{ROOT}/{RUN_ID}/preprocess/ls_synth_fill.tif",
            f"{ROOT}/{RUN_ID}/preprocess/ls_synth_hillshade.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/cosmo_landslides_single.param",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_INVENTORY.csv",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_TFAIL.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_EXCAVATION.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_single/ls_single_LS_CONC.tif",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/sample_points.csv",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/cosmo_landslides_series.param",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/ls_series_LS_INVENTORY.csv",
            f"{ROOT}/{RUN_ID}/analysis/cosmo_landslides_series/ls_series_LS_SAMPLE_EFFECT.csv",
            f"{ROOT}/{RUN_ID}/reports/cosmo_landslides/ls_single_landslide_outlines.fgb",
            f"{ROOT}/{RUN_ID}/reports/cosmo_landslides/ls_single_larsen_size_comparison.png",
            f"{ROOT}/{RUN_ID}/reports/cosmo_landslides/ls_series_apparent_vs_true_erosion.png",
            f"{ROOT}/{RUN_ID}/reports/cosmo_landslides/ls_series_sample_bias_ratio.png",
        ]
    if MODULE == "check_basins":
        return [
            f"{ROOT}/{RUN_ID}/fetch/fetch_bbox.param",
            f"{ROOT}/{RUN_ID}/fetch/dem_DEM.tif",
            f"{ROOT}/{RUN_ID}/preprocess/raster_preprocessing.param",
            f"{ROOT}/{RUN_ID}/preprocess/dem_fill.tif",
            f"{ROOT}/{RUN_ID}/extract/check_basins.param",
            f"{ROOT}/{RUN_ID}/extract/check_junction_basin_outlines.fgb",
            f"{ROOT}/{RUN_ID}/extract/check_channel_network_lines.fgb",
        ]
    return []


rule write_run_manifest:
    input:
        manifest_inputs
    output:
        manifest=f"{ROOT}/{RUN_ID}/reports/{RUN_ID}_RUN_MANIFEST.json"
    params:
        workflow_name=config["workflow_name"],
        run_id=RUN_ID,
        module=MODULE
    run:
        import hashlib
        import json
        from pathlib import Path

        def sha256(path: str) -> str:
            h = hashlib.sha256()
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(1024 * 1024), b""):
                    h.update(chunk)
            return h.hexdigest()

        out_path = Path(output.manifest)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        files = [{"path": p, "sha256": sha256(p)} for p in list(input)]
        payload = {
            "workflow_name": params.workflow_name,
            "run_id": params.run_id,
            "module": params.module,
            "files": files,
        }
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
