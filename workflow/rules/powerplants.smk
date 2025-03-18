"""Rules for processing powerstations."""


rule powerplants_adjust_location:
    message:
        "Adjusting hydro powerplant location to the nearest shape and basin."
    params:
        crs=config["crs"],
        basin_adjustment=config["powerplants"]["basin_adjustment"],
    input:
        basins="results/hydrobasin_global.parquet",
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
        shapes="resources/user/shapes.parquet",
        powerplants="results/adjusted_powerplants.parquet",
        basins="results/hydrobasin_global.parquet",
        cutout="resources/automatic/cutout.nc",
    output:
        inflow="results/inflow_m3.nc",
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/powerplants_get_inflow_m3.py"


# rule powerplants_get_inflow_mwh:
#     message:
#         "Adjusting powerplant generation using historical data."
#     params:
#         max_capacity_factor=config["powerplants"]["max_capacity_factor"],
#         start_year=config["years"]["start"],
#         end_year=config["years"]["end"],
#     input:
#         inflow_m3="results/inflow_m3.nc",
#         powerplants="results/adjusted_powerplants.parquet",
#         shapes="resources/user/shapes.parquet",
#         generation="resources/user/generation.parquet"
#     output:
#         inflow_mwh="results/inflow_mwh.nc"
#     conda:
#         "../envs/default.yaml"
#     script:
#         "../scripts/powerplants_get_inflow_mwh.py"
