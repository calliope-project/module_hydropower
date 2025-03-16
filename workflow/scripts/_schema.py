"""Schemas for user given resources and module results."""

import pandera as pan
from pandera import DataFrameModel, Field
from pandera.typing import Series
from pandera.typing.geopandas import GeoSeries

TYPES = ["hydro_run_of_river", "hydro_dam", "hydro_pumped_storage"]


class Powerplants(DataFrameModel):
    powerplant_id: Series[str] = Field(unique=True)
    net_generation_capacity_mw: Series[float] = Field(gt=0)
    net_pumping_capacity_mw: Series[float] = Field(nullable=True)
    storage_capacity_mwh: Series[float] = Field(nullable=True)
    powerplant_type: Series[str] = Field(isin=TYPES)
    geometry: GeoSeries


shapes_schema = pan.DataFrameSchema(
    {
        "shape_id": pan.Column(str, unique=True),
        "country_id": pan.Column(str),
        "class": pan.Column(str, pan.Check.isin(["land"])),
        "geometry": pan.Column("geometry"),
    },
    coerce=True,
    strict=False,
)
