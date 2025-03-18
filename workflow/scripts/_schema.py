"""Schemas for user given resources and module results."""

import pandera as pan

TYPES = ["hydro_run_of_river", "hydro_dam", "hydro_pumped_storage"]

powerplants_schema = pan.DataFrameSchema(
    {
        "powerplant_id": pan.Column(
            str, unique=True, description="Unique powerplant ID."
        ),
        "net_generation_capacity_mw": pan.Column(
            float, pan.Check.gt(0), description="Net Generation Capacity in Megawatts."
        ),
        "net_pumping_capacity_mw": pan.Column(
            float, nullable=True, description="Net pumping capacity in Megawatts."
        ),
        "storage_capacity_mwh": pan.Column(
            float, nullable=True, description="Storage capacity in Megawatt hour."
        ),
        "powerplant_type": pan.Column(
            str,
            pan.Check.isin(TYPES),
            description=f"Type of hydro powerplant. One of: {TYPES}.",
        ),
        "geometry": pan.Column("geometry"),
    },
    coerce=True,
    strict=False,
)

shapes_schema = pan.DataFrameSchema(
    {
        "shape_id": pan.Column(
            str, unique=True, description="Unique ID of this shape."
        ),
        "country_id": pan.Column(str, description="ISO alpha-3 code."),
        "class": pan.Column(
            str,
            pan.Check.isin(["land"]),
            description="Shape class. Only land is accepted.",
        ),
        "geometry": pan.Column("geometry"),
    },
    coerce=True,
    strict=False,
)

generation_schema = pan.DataFrameSchema(
    {
        "country_id": pan.Column(str, description="ISO alpha-3 code."),
        "year": pan.Column(
            int, pan.Check.ge(0), description="Year of the data sample."
        ),
        "generation_mwh": pan.Column(
            float,
            pan.Check.ge(0),
            description="Combined hydropower generation at the given year.",
        ),
    },
    coerce=True,
    strict=False,
)
