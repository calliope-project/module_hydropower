"""Schemas for user given resources and module results."""

import pandera as pan

TYPES = ["hydro_run_of_river", "hydro_dam", "hydro_pumped_storage"]

powerplants_schema = pan.DataFrameSchema(
    {
        "powerplant_id": pan.Column(str, unique=True),
        "net_generation_capacity_mw": pan.Column(float, pan.Check.gt(0)),
        "net_pumping_capacity_mw": pan.Column(float, nullable=True),
        "storage_capacity_mwh": pan.Column(float, nullable=True),
        "powerplant_type": pan.Column(str, pan.Check.isin(TYPES)),
        "geometry": pan.Column("geometry"),
    },
    coerce=True,
    strict=True
)

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
