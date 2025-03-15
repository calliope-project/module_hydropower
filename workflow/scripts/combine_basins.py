"""Combine all requested HydroBASIN datasets into one."""

from typing import TYPE_CHECKING, Any

import geopandas as gpd
import pandas as pd

if TYPE_CHECKING:
    snakemake: Any


def combine_basins(continent_files, global_file):
    """Extract a specific pfafstetter lavel from a HydroBASINS zip file."""
    continents = [gpd.read_parquet(i) for i in continent_files]
    continent_crs = [i.crs for i in continents]
    assert len(set(continent_crs)) == 1

    combined = pd.concat(continents)
    combined = combined.reset_index(drop=True)
    combined.crs = continent_crs[0]

    combined.to_parquet(global_file)


if __name__ == "__main__":
    combine_basins(
        continent_files=snakemake.input.continent_files,
        global_file=snakemake.output.global_file,
    )
