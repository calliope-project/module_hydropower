"""Rules for processing powerstations."""


rule powerplants_adjust_location:
    message:
        "Adjusting hydro powerplant location to the nearest shape and basin."
    params:
        crs=config["crs"],
        basin_adjustment=config["powerplants"]["basin_adjustment"],
    input:
        basins=ancient(f"resources/automatic/hydrobasin_global_{config["pfafstetter_level"]}.parquet"),
        powerplants="resources/user/powerplants.parquet",
        shapes="resources/user/shapes.parquet"
    output:
        adjusted_powerplants="results/adjusted_powerplants.parquet",
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/powerplants_adjust_location.py"


rule powerplants_get_inflow_m3:
    message:
        "Calculating hydro powerplant inflow in m3."
    input:
        basins=ancient(f"resources/automatic/hydrobasin_global_{config["pfafstetter_level"]}.parquet"),
        shapes="resources/user/shapes.parquet",
        adjusted_powerplants="results/adjusted_powerplants.parquet",
        cutout=ancient("resources/automatic/cutout.nc"),
    output:
        inflow="results/by_powerplant_id/inflow_m3.parquet",
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/powerplants_get_inflow_m3.py"


rule powerplants_get_inflow_mwh:
    message:
        "Calculating powerplant generation in MWh and applying corrections using historical data."
    input:
        inflow_m3="results/by_powerplant_id/inflow_m3.parquet",
        adjusted_powerplants="results/adjusted_powerplants.parquet",
        generation="resources/user/national_generation.parquet"
    output:
        inflow_mwh="results/by_powerplant_id/inflow_mwh.parquet"
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/powerplants_get_inflow_mwh.py"


rule powerplants_get_cf_per_shape:
    message:
        "Calculating capacity factor timeseries per shape."
    input:
        adjusted_powerplants="results/adjusted_powerplants.parquet",
        inflow_mwh="results/by_powerplant_id/inflow_mwh.parquet"
    output:
        hydro_run_of_river="results/by_shape_id/hydro_run_of_river_cf.parquet",
        hydro_dam="results/by_shape_id/hydro_dam_cf.parquet"
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/powerplants_get_cf_per_shape.py"
