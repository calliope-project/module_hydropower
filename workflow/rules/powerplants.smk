"""Rules for processing powerstations."""

rule powerplants_adjust_location:
    message:
        "Adjusting hydro powerplant location to the nearest basin within buffer radius {params.buffer_radius}."
    params:
        geographic_crs = config["crs"]["geographic"],
        projected_crs = config["crs"]["projected"],
        buffer_radius = config["powerplants"]["basin_adjustment"]["buffer_radius"],
        max_dropped = config["powerplants"]["basin_adjustment"]["max_dropped"]
    input:
        powerplants = "resources/user/powerplants.parquet",
        basins = "results/hydrobasin_global.parquet"
    output:
        adjusted_powerplants = "results/adjusted_powerplants.parquet"
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/powerplants_adjust_location.py"
