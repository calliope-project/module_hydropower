"""Rules to used to download automatic resource files."""


rule download_basin:
    message:
        "Downloading HydroBASINS file for '{wildcards.continent}'."
    params:
        url=lambda wc: internal["resources"]["automatic"]["HydroBASINS"].format(
            continent=wc.continent
        ),
    output:
        temp("resources/automatic/hydrobasin_{continent}.zip"),
    wildcard_constraints:
        continent="|".join(internal["continent_codes"]),
    conda:
        "../envs/shell.yaml"
    log:
        "logs/download_basin_{continent}.log",
    shell:
        "curl -sSLo {output} '{params.url}' "


rule download_cutout:
    message:
        "Downloading runoff cutout from {params.start_year}-01-01 to {params.end_year}-12-31."
    params:
        era5_crs=internal["era5_crs"],
        start_year=config["years"]["start"],
        end_year=config["years"]["end"],
    input:
        shapes="resources/user/shapes.parquet",
    output:
        cutout="resources/automatic/cutout.nc",
        plot=report(
            "resources/automatic/cutout.png",
            caption="../report/cutout.rst",
            category="Hydropower module",
        ),
    conda:
        "../envs/default.yaml"
    log:
        "logs/download_cutout.log",
    script:
        "../scripts/download_cutout.py"
