"""Combine all requested HydroBASIN datasets into one."""

import sys
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

if TYPE_CHECKING:
    snakemake: Any
sys.stderr = open(snakemake.log[0], "w")


def _plot_combined_basins(global_file, path):
    combined = gpd.read_parquet(global_file)
    ax = combined.plot(figsize=(20, 12))
    ax.set_title("Global hydro basins")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    plt.savefig(path, bbox_inches="tight")


def basins_combine_continents(continent_files, global_file):
    """Extract a specific pfafstetter level from a HydroBASINS zip file."""
    continents = [gpd.read_parquet(i) for i in continent_files]
    continent_crs = [i.crs for i in continents]
    assert len(set(continent_crs)) == 1

    combined = pd.concat(continents)
    combined = combined.reset_index(drop=True)
    combined.crs = continent_crs[0]

    combined.to_parquet(global_file)


if __name__ == "__main__":
    basins_combine_continents(
        continent_files=snakemake.input.continent_files,
        global_file=snakemake.output.global_file,
    )
    _plot_combined_basins(
        global_file=snakemake.output.global_file, path=snakemake.output.plot
    )
