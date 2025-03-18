"""Get a cutout with runoff data."""

from typing import TYPE_CHECKING, Any

import atlite
import geopandas as gpd
from pyproj import CRS

if TYPE_CHECKING:
    snakemake: Any


def runoff_cutout(input_shapes, era5_crs, start_year, end_year, output_netcdf):
    """Download runoff data from ERA5."""
    # Run quick checks before executing the download
    assert CRS(era5_crs).is_geographic

    shapes = gpd.read_parquet(input_shapes)
    shapes = shapes.to_crs(era5_crs)
    bounds = shapes.total_bounds

    cutout = atlite.Cutout(
        path=output_netcdf,
        module=["era5"],
        x=slice(bounds[0], bounds[2]),
        y=slice(bounds[1], bounds[3]),
        time=slice(f"{start_year}-01-01", f"{end_year}-12-31"),
    )
    cutout.prepare(features=["runoff"])
    assert cutout.crs == era5_crs


if __name__ == "__main__":
    runoff_cutout(
        input_shapes=snakemake.input.shapes,
        era5_crs=snakemake.params.era5_crs,
        start_year=snakemake.params.start_year,
        end_year=snakemake.params.end_year,
        output_netcdf=snakemake.output.cutout,
    )
