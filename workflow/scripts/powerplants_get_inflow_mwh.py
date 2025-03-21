"""Obtain magnitude-adjusted energy inflows for each powerplant in MWh."""

from calendar import isleap
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import numpy as np
import pandas as pd
from _schema import NationalGenerationSchema, PowerplantSchema
from scipy.optimize import minimize

if TYPE_CHECKING:
    snakemake: Any


# ---
# Taken from Euro-Calliope (MIT licensed)s
# https://github.com/calliope-project/euro-calliope/blob/e9cb908a5b4c6274148a16b59e5dd0b412aaf560/scripts/hydro/inflow_mwh.py
def _water_inflow_m3_to_mwh(inflow_m3: pd.Series, annual_generation: float, cap: float):
    def generation(scaling_factor):
        generation = inflow_m3 * scaling_factor
        generation[generation > cap] = cap
        return generation

    def residual(scaling_factor):
        return abs(generation(scaling_factor).sum() - annual_generation)

    x0 = annual_generation / inflow_m3.sum()
    res = minimize(residual, x0, method="nelder-mead", options={"xatol": 1e-10})
    assert res.success, print(res)
    assert res.fun < 1, print(res)  # error smaller than 1 MWh

    return generation(res.x)


# ---


def estimate_annual_powerplant_generation(
    powerplants: pd.DataFrame,
    inflow_m3: pd.DataFrame,
    national_generation: pd.DataFrame,
    year: int,
) -> pd.DataFrame:
    """Obtain magnitude-corrected hydropower timeseries dataset for a given year."""
    inflow_m3_yr = inflow_m3[inflow_m3.index.year == year]
    generation_yr = national_generation[national_generation.year == year]
    plants_by_id = powerplants.set_index("powerplant_id")

    # Re-scale capacity share to account for erroneous zero timesteps in the inflow data
    # TODO: identify the reasoning behind this
    scaling_factor = inflow_m3_yr.where(inflow_m3_yr > 1).count(
        axis="index"
    ) / inflow_m3_yr.count(axis="index")

    national_cap_share_per_powerplant = (
        plants_by_id["net_generation_capacity_mw"]
        .multiply(scaling_factor, axis="index")
        .groupby(plants_by_id["country_id"])
        .transform(lambda x: x / x.sum())
    )

    annual_national_generation = generation_yr.set_index("country_id")["generation_mwh"]
    assert annual_national_generation.index.is_unique

    annual_powerplant_mwh = plants_by_id.apply(
        lambda x: annual_national_generation[x.country_id]
        * national_cap_share_per_powerplant[x.name],
        axis="columns",
    )
    days_in_year = 366 if isleap(year) else 365
    inflow_mwh_yr = pd.DataFrame(
        np.nan, index=inflow_m3_yr.index, columns=inflow_m3_yr.columns
    )
    for powerplant_id, data in plants_by_id.iterrows():
        m3_timeseries = inflow_m3_yr[powerplant_id]
        annual_generation = annual_powerplant_mwh[powerplant_id]
        net_capacity = data.net_generation_capacity_mw
        max_generation = net_capacity * days_in_year * 24
        if annual_generation > max_generation:
            raise ValueError(
                f"Annual generation ({annual_generation}) exceeds the maximum ({max_generation})."
            )

        inflow_mwh_yr[powerplant_id] = _water_inflow_m3_to_mwh(
            m3_timeseries, annual_generation, net_capacity
        )
    return inflow_mwh_yr


def powerplants_get_inflow_mwh(
    inflow_m3_file: str,
    powerplants_file: str,
    national_generation_file: str,
    inflow_mwh_file: str,
):
    """Generate hydropower timeseries using unscaled water time series.

    The inflow timeseries will be scaled utilising national generation data,
    with the aim of correcting its magnitude while retaining its dynamics.

    General assumptions:
        - Annual generation never exceeds installed capacity.
        - Hydro dams have no spillage.

    Args:
        inflow_m3_file (str): Dataset with water inflow per-powerplant in m3.
        powerplants_file (str): Powerplants dataset (adjusted).
        national_generation_file (str): Annual national hydropower generation per country.
        inflow_mwh_file (str): Resulting file with energy inflow per powerplant in MWh.
    """
    inflow_m3 = pd.read_parquet(inflow_m3_file)
    powerplants = gpd.read_parquet(powerplants_file)
    generation = pd.read_parquet(national_generation_file)
    PowerplantSchema.validate(powerplants)
    NationalGenerationSchema.validate(generation)

    year_results = []
    for year in sorted(inflow_m3.index.year.unique()):
        inflow_mwh_yr = estimate_annual_powerplant_generation(
            powerplants, inflow_m3, generation, year
        )
        year_results.append(inflow_mwh_yr)

    inflow_mwh = pd.concat(year_results)
    inflow_mwh.attrs = {"long_name": "Energy inflow", "units": "Megawatt hour"}
    inflow_mwh.to_parquet(inflow_mwh_file)


if __name__ == "__main__":
    powerplants_get_inflow_mwh(
        inflow_m3_file=snakemake.input.inflow_m3,
        powerplants_file=snakemake.input.adjusted_powerplants,
        national_generation_file=snakemake.input.generation,
        inflow_mwh_file=snakemake.output.inflow_mwh,
    )
