"""Capacity factor calculation for hydropower basin and run-of-river plants."""

from typing import TYPE_CHECKING, Any

import geopandas as gpd
import numpy as np
import pandas as pd
from _schema import POWERPLANT_TYPES, PowerplantSchema

if TYPE_CHECKING:
    snakemake: Any


def _get_capacity_factors_timeseries(
    tech: str, powerplants: pd.DataFrame, inflow_mwh: pd.DataFrame
) -> pd.DataFrame:
    """Calculate capacity factor timeseries within a shape for a given technology.

    Args:
        tech (str): name of the powerplant technology.
        powerplants (pd.DataFrame): powerplant data file.
        inflow_mwh (pd.DataFrame): timeseries of energy inflow per powerplant.

    Returns:
        pd.DataFrame: capacity factor timeseries (row: timestep, column: shape_id).
    """
    tech_powerplants = powerplants[powerplants["powerplant_type"] == tech]
    group = tech_powerplants.groupby(["shape_id"])
    shape_net_cap = group["net_generation_capacity_mw"].sum()
    shape_powerplants = group["powerplant_id"].apply(list)

    shape_ids = sorted(tech_powerplants.shape_id.unique())
    cf_timeseries = pd.DataFrame(np.nan, index=inflow_mwh.index, columns=shape_ids)
    for shape_id in shape_ids:
        cf_timeseries[shape_id] = (
            inflow_mwh[shape_powerplants[shape_id]].sum(axis="columns")
            / shape_net_cap[shape_id]
        )

    if cf_timeseries.isna().any().any():
        ValueError(
            f"Calculated capacity factor timeseries must not contain null values {tech}."
        )

    cf_timeseries.attrs = {
        "long_name": "Capacity factors",
        "units": None,
        "powerplant_type": tech,
    }
    return cf_timeseries


def powerplants_get_cf_per_shape(
    powerplants_file: gpd.GeoDataFrame,
    inflow_mwh_file: pd.DataFrame,
    output_paths: dict[str, str],
):
    """Construct a capacity factor file for each type of hydro plant."""
    powerplants = gpd.read_parquet(powerplants_file)
    inflow_mwh = pd.read_parquet(inflow_mwh_file)

    PowerplantSchema.validate(powerplants)

    for plant_type in POWERPLANT_TYPES:
        cap_factors = _get_capacity_factors_timeseries(
            plant_type, powerplants, inflow_mwh
        )
        cap_factors.to_parquet(output_paths.get(plant_type))


if __name__ == "__main__":
    powerplants_get_cf_per_shape(
        powerplants_file=snakemake.input.adjusted_powerplants,
        inflow_mwh_file=snakemake.input.inflow_mwh,
        output_paths=snakemake.output,
    )
