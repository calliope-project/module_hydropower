"""Transform from HydroBASINS zipped files to parquet."""

from typing import TYPE_CHECKING, Any

import geopandas as gpd

if TYPE_CHECKING:
    snakemake: Any

SHP_ZIP_PATH = "hybas_{continent}_lev{level}_v1c.shp"


def extract_basin_pfafstetter_level(zip_file, continent, level, parquet_file):
    """Extract a specific pfafstetter level from a HydroBASINS zip file."""
    inner_path = SHP_ZIP_PATH.format(continent=continent, level=level)
    gdf = gpd.read_file(f"{zip_file}!{inner_path}")
    gdf.to_parquet(parquet_file)


if __name__ == "__main__":
    extract_basin_pfafstetter_level(
        zip_file=snakemake.input.zip_file,
        continent=snakemake.params.continent,
        level=snakemake.params.level,
        parquet_file=snakemake.output.parquet_file,
    )
