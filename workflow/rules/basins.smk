"""Rules for processing the HydroBASINS dataset."""

rule basins_extract_pfafstetter_level:
    message: "Unzipping HydroBASINS file for '{wildcards.continent}' at pfafstetter level '{params.level}'."
    params:
        level = config["pfafstetter_level"],
        continent = lambda wc: wc.continent
    input:
        zip_file = "resources/automatic/hydrobasin_{continent}.zip"
    output:
        parquet_file = "resources/automatic/hydrobasin_{continent}.parquet"
    wildcard_constraints:
        continent = "|".join(internal["continent_codes"])
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/basins_extract_pfafstetter_level.py"


rule basins_combine_continents:
    message: "Combine all HydroBASINS into a single dataset."
    input:
        continent_files = expand("resources/automatic/hydrobasin_{continent}.parquet", continent=internal["continent_codes"])
    output:
        global_file = "results/hydrobasin_global.parquet"
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/basins_combine_continents.py"
