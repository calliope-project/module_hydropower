"""Calculate the water inflow of relevant powerplants."""

from pathlib import Path
from typing import TYPE_CHECKING, Any

import atlite
import geopandas as gpd

if TYPE_CHECKING:
    snakemake: Any

INFLOW_TYPES = ["hydro_run_of_river", "hydro_dam"]


def powerplants_get_inflow(
    shapes_file: Path,
    powerplants_file: Path,
    basins_file: Path,
    cutout_file: Path,
    inflow_file: Path,
):
    """Get powerplant inflow using 'atlite'.

    Args:
        shapes_file (Path): user shapes.
        powerplants_file (Path): adjusted powerplants file.
        basins_file (Path): global basins file.
        cutout_file (Path): atlite cutout fitting the user shape.
        inflow_file (Path): resulting water inflow per powerplant.
    """
    shapes = gpd.read_parquet(shapes_file)
    powerplants = gpd.read_parquet(powerplants_file)
    basins = gpd.read_parquet(basins_file)
    cutout = atlite.Cutout(path=cutout_file)

    # atlite does not support CRS conversion, adapt data to it instead
    era5_crs = cutout.crs
    assert era5_crs.is_geographic

    shapes = shapes.to_crs(era5_crs)
    powerplants = powerplants.to_crs(era5_crs)
    basins = basins.to_crs(era5_crs)

    # Process only relevant powerplants
    powerplants = powerplants[powerplants["powerplant_type"].isin(INFLOW_TYPES)]
    within = powerplants.apply(
        lambda x: shapes.contains(x.geometry).any(), axis="columns"
    )
    powerplants_within = powerplants[within].copy()

    # Re-format for atlite
    powerplants_within["lon"] = powerplants_within.geometry.apply(lambda point: point.x)
    powerplants_within["lat"] = powerplants_within.geometry.apply(lambda point: point.y)
    powerplants_within = powerplants_within.drop("geometry", axis="columns").set_index(
        "powerplant_id"
    )
    # Calculate inflow
    inflow = cutout.hydro(plants=powerplants_within, hydrobasins=basins)
    inflow.attrs = {
        "units": "cubic_meter",
        "long_name": "Water inflow"
    }
    # Reformat and add useful metadata
    inflow = inflow.rename(plant="powerplant_id")
    inflow.name = "inflow"

    inflow.to_netcdf(inflow_file)


if __name__ == "__main__":
    powerplants_get_inflow(
        shapes_file=snakemake.input.shapes,
        powerplants_file=snakemake.input.powerplants,
        basins_file=snakemake.input.basins,
        cutout_file=snakemake.input.cutout,
        inflow_file=snakemake.output.inflow
    )
