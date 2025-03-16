"""Rules to used to download automatic resource files."""


rule download_basin:
    message: "Downloading HydroBASINS file for '{wildcards.continent}'."
    params:
        url = lambda wc: internal["resources"]["automatic"]["HydroBASINS"].format(continent=wc.continent),
    output: temp("resources/automatic/hydrobasin_{continent}.zip")
    wildcard_constraints:
        continent = "|".join(internal["continent_codes"])
    conda: "../envs/shell.yaml"
    shell:
        "curl -sSLo {output} '{params.url}' "


rule download_cutout:
    message: "Downloading runoff cutout from {params.start_date} to {params.end_date}."
    params:
        geographic_crs = config["crs"]["geographic"],
        start_date = config["dates"]["start"],
        end_date = config["dates"]["end"]
    input:
        shapes = "resources/user/shapes.parquet"
    output:
        cutout = "resources/automatic/cutout.nc"
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/runoff_cutout.py"
