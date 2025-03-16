"""Adjustment of powerplant position to their nearest basin."""

from pathlib import Path
from typing import TYPE_CHECKING, Any

import _schema as schema
import geopandas as gpd
from pyproj import CRS

if TYPE_CHECKING:
    snakemake: Any


def powerplants_adjust_location(
    powerplants_path: Path,
    basins_path: Path,
    geographic_crs: str,
    projected_crs: str,
    buffer_radius: float,
    max_dropped: int,
    adjusted_powerplants_path: Path,
):
    """Adjust powerplants to their closest basin.

    Mismatched stations outside of the buffer radius will be dropped.

    Args:
        powerplants_path (Path): user inputted powerplant data.
        basins_path (Path): merged global basins.
        geographic_crs (str): geographic crs (e.g., EPSG:4326).
        projected_crs (str): projected crs to use for distance-based operations (e.g., EPSG:3857).
        buffer_radius (float): maximum distance allowed for the adjustment. Matches the unit of the projected CRS.
        max_dropped (int): maximum allowed number of stations to drop.
        adjusted_powerplants_path (Path): location to place the ajusted powerplant dataset.

    Raises:
        ValueError: dropped powerplants outside of bounds exceeded the configured maximum.
    """
    # Check the correctness of the requested crs
    assert CRS(projected_crs).is_projected
    assert CRS(geographic_crs).is_geographic

    # Read and validate input files
    powerplants = schema.powerplants_schema(gpd.read_parquet(powerplants_path))
    basins = gpd.read_parquet(basins_path)

    powerplants = powerplants.to_crs(geographic_crs)
    basins = basins.to_crs(geographic_crs)
    outside = powerplants.apply(
        lambda x: not basins.contains(x.geometry).any(), axis="columns"
    )

    # Identify and fix powerplants outside of bounds
    # Distance-based operations should use a projected CRS
    powerplants = powerplants.to_crs(projected_crs)
    basins = basins.to_crs(projected_crs)
    for index, data in powerplants[outside].iterrows():
        point = data["geometry"]
        distances = basins.distance(point)
        min_distance = distances.min()
        min_id = distances.idxmin()
        if min_distance <= buffer_radius:
            intersection = point.buffer(buffer_radius).intersection(
                basins.loc[min_id, "geometry"]
            )
            powerplants.loc[index, "geometry"] = intersection.representative_point()
            outside = outside.drop(index, axis="index")

    # Drop powerplants outside of basins
    if len(powerplants[outside]) > max_dropped:
        raise ValueError(
            f"Dropped powerplants ({len(outside)}) exceeds the configured maximum ({max_dropped})."
        )
    powerplants = powerplants[~outside]

    powerplants = powerplants.to_crs(geographic_crs)
    powerplants.to_parquet(adjusted_powerplants_path)


if __name__ == "__main__":
    powerplants_adjust_location(
        powerplants_path=snakemake.input.powerplants,
        basins_path=snakemake.input.basins,
        geographic_crs=snakemake.params.geographic_crs,
        projected_crs=snakemake.params.projected_crs,
        buffer_radius=snakemake.params.buffer_radius,
        max_dropped=snakemake.params.max_dropped,
        adjusted_powerplants_path=snakemake.output.adjusted_powerplants,
    )
