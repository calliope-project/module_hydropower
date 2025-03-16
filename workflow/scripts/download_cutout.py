"""Get a cutout with runoff data."""

from datetime import date
from typing import TYPE_CHECKING, Any

import atlite
import geopandas as gpd
from pyproj import CRS

if TYPE_CHECKING:
    snakemake: Any


def runoff_cutout(input_shapes, era5_crs, start_date, end_date, output_netcdf):
    """Download runoff data from ERA5."""
    # Run quick checks before executing the download
    assert CRS(era5_crs).is_geographic
    assert date.fromisoformat(start_date)
    assert date.fromisoformat(end_date)
    assert start_date != end_date, "Single-day requests are not supported."

    shapes = gpd.read_parquet(input_shapes)
    shapes = shapes.to_crs(era5_crs)
    bounds = shapes.total_bounds

    cutout = atlite.Cutout(
        path=output_netcdf,
        module=["era5"],
        x=slice(bounds[0], bounds[2]),
        y=slice(bounds[1], bounds[3]),
        time=slice(start_date, end_date),
    )
    cutout.prepare(features=["runoff"])
    assert cutout.crs == era5_crs


if __name__ == "__main__":
    runoff_cutout(
        input_shapes=snakemake.input.shapes,
        era5_crs=snakemake.params.era5_crs,
        start_date=snakemake.params.start_date,
        end_date=snakemake.params.end_date,
        output_netcdf=snakemake.output.cutout,
    )
