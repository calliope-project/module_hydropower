"""Rules to used to download automatic resource files."""


rule download_basin:
    message: "Downloading HydroBASINS file for '{wildcards.continent}'."
    conda: "../envs/shell.yaml"
    params:
        url = lambda wc: internal["resources"]["automatic"]["HydroBASINS"].format(continent=wc.continent),
    output: temp("resources/automatic/hydrobasin_{continent}.zip")
    wildcard_constraints:
        continent = "|".join(internal["continent_codes"])
    shell:
        "curl -sSLo {output} '{params.url}' "
