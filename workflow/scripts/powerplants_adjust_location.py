"""Adjustment of powerplant position to their nearest basin."""

from pathlib import Path
from typing import TYPE_CHECKING, Any

import geopandas as gpd
from _schema import PowerplantSchema, ShapeSchema
from pyproj import CRS

if TYPE_CHECKING:
    snakemake: Any


def powerplants_adjust_location(
    basins_path: Path,
    powerplants_path: Path,
    shapes_path: Path,
    crs: dict[str, str],
    basin_adjustment: dict[str, int],
    adjusted_powerplants_path: Path,
):
    """Adjust powerplants to their closest basin.

    Mismatched stations outside of the buffer radius will be dropped.

    Args:
        basins_path (Path): merged global basins.
        powerplants_path (Path): user inputted powerplant data.
        shapes_path (Path): user inputted shapes.
        crs (dict[str, str]): geographic and projected CRS codes to use for operations.
        basin_adjustment (dict[str, int]): settings to use during adjustment operation.
        adjusted_powerplants_path (Path): location to place the ajusted powerplant dataset.

    Raises:
        ValueError: dropped powerplants outside of bounds exceeded the configured maximum.
    """
    # Check the correctness of the requested crs
    assert CRS(crs["projected"]).is_projected
    assert CRS(crs["geographic"]).is_geographic

    # Read and validate input files
    basins = gpd.read_parquet(basins_path)
    powerplants = gpd.read_parquet(powerplants_path)
    PowerplantSchema.validate(powerplants)
    shapes = gpd.read_parquet(shapes_path)
    ShapeSchema.validate(shapes)

    # Coordinate-based operations must use a geographic CRS
    basins = basins.to_crs(crs["geographic"])
    powerplants = powerplants.to_crs(crs["geographic"])
    shapes = shapes.to_crs(crs["geographic"])

    # Identify powerplants outside of basins / shapes
    outside_basins = powerplants[
        powerplants.geometry.apply(lambda x: not basins.contains(x).any())
    ]
    outside_shapes = powerplants[
        powerplants.geometry.apply(lambda x: not shapes.contains(x).any())
    ]
    total_outside = len(set(outside_basins.index) | set(outside_shapes.index))
    max_outside = basin_adjustment["max_outside"]
    if total_outside > max_outside:
        raise ValueError(
            f"Powerplants outside basins and/or shapes ({total_outside}) exceeds the maximum ({max_outside})."
        )

    # Distance-based operations should use a projected CRS
    basins = basins.to_crs(crs["projected"])
    powerplants = powerplants.to_crs(crs["projected"])
    shapes = shapes.to_crs(crs["projected"])
    buffer_radius = basin_adjustment["buffer_radius"]

    # Shift powerplant location to the nearest basin if it is within the buffer
    to_drop = set(outside_basins.index)
    for index, data in outside_basins.iterrows():
        point = data["geometry"]
        distances = basins.distance(point)
        min_distance = distances.min()
        min_id = distances.idxmin()
        if min_distance <= buffer_radius:
            intersection = point.buffer(buffer_radius).intersection(
                basins.loc[min_id, "geometry"]
            )
            powerplants.loc[index, "geometry"] = intersection.representative_point()
            to_drop.remove(index)

    # Assign country / shape IDs to each adjusted plant if it is within the buffer
    powerplants = powerplants.sjoin_nearest(
        shapes[["country_id", "shape_id", "geometry"]],
        how="left",
        max_distance=buffer_radius,
    )

    # Powerplants outside the buffer will be dropped
    to_drop |= set(powerplants[powerplants["shape_id"].isna()].index)
    total_dropped = len(to_drop)
    max_dropped = basin_adjustment["max_dropped"]
    if total_dropped > max_dropped:
        raise ValueError(
            f"Dropped powerplants ({total_dropped}) exceeds the configured maximum of ({max_dropped})."
        )
    powerplants = powerplants.drop(to_drop, axis="index")

    # Re-validate and save
    powerplants = powerplants.to_crs(crs["geographic"])
    PowerplantSchema.validate(powerplants)
    powerplants.to_parquet(adjusted_powerplants_path)


if __name__ == "__main__":
    powerplants_adjust_location(
        basins_path=snakemake.input.basins,
        powerplants_path=snakemake.input.powerplants,
        shapes_path=snakemake.input.shapes,
        crs=snakemake.params.crs,
        basin_adjustment=snakemake.params.basin_adjustment,
        adjusted_powerplants_path=snakemake.output.adjusted_powerplants,
    )
