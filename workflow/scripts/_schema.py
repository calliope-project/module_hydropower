"""Schemas for user given resources and module results."""

from pandera import DataFrameModel, Field
from pandera.typing import Index, Series
from pandera.typing.geopandas import GeoSeries

POWERPLANT_TYPES = ["hydro_run_of_river", "hydro_dam"]


class PowerplantSchema(DataFrameModel):
    class Config:
        coerce = True
        strict = False

    powerplant_id: Series[str] = Field(unique=True)
    "Unique powerplant ID."
    net_generation_capacity_mw: Series[float] = Field(ge=0)
    "Net Generation Capacity in Megawatts."
    storage_capacity_mwh: Series[float] = Field(nullable=True)
    "Storage capacity in Megawatt hour."
    powerplant_type: Series[str] = Field(isin=POWERPLANT_TYPES)
    "Type of hydropower plant."
    geometry: GeoSeries
    "Plant location (Point)."
    index: Index[int] = Field(unique=True)


class ShapeSchema(DataFrameModel):
    class Config:
        coerce = True
        strict = False

    shape_id: Series[str] = Field(unique=True)
    "Unique ID for this shape."
    country_id: Series[str]
    "ISO alpha-3 code."
    shape_class: Series[str] = Field(isin=["land"])
    "Shape class. Only 'land' is accepted."
    geometry: GeoSeries
    "Shape polygon."
    index: Index[int] = Field(unique=True)


class NationalGenerationSchema(DataFrameModel):
    class Config:
        coerce = True
        strict = False

    country_id: Series[str]
    "ISO alpha-3 code."
    year: Series[int]
    "Year of the data sample."
    generation_mwh: Series[float] = Field(ge=0)
    "Combined hydropower generation at the given year."
    index: Index[int] = Field(unique=True)
