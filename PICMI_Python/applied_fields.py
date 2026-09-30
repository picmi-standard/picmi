"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
from typing import ClassVar, Self

from pydantic import Field, model_validator

from .base import Expression, PICMI_AppliedField, _PICMIModel, PICMI_ExpressionParameters

# ---------------
# Applied fields
# ---------------


class PICMI_ConstantAppliedField(PICMI_AppliedField):
    """
    Describes a constant applied field
    """
    Ex: float | None = Field(
        default=None,
        description="Constant Ex field [V/m]"
    )
    Ey: float | None = Field(
        default=None,
        description="Constant Ey field [V/m]"
    )
    Ez: float | None = Field(
        default=None,
        description="Constant Ez field [V/m]"
    )
    Bx: float | None = Field(
        default=None,
        description="Constant Bx field [T]"
    )
    By: float | None = Field(
        default=None,
        description="Constant By field [T]"
    )
    Bz: float | None = Field(
        default=None,
        description="Constant Bz field [T]"
    )
    lower_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Lower bound of the region where the field is applied [m]."
    )
    upper_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Upper bound of the region where the field is applied [m]"
    )


class PICMI_AnalyticAppliedField(PICMI_AppliedField, PICMI_ExpressionParameters):
    """
    Describes an analytic applied field

    The expressions should be in terms of the position and time, written as 'x', 'y', 'z', 't'.
    Parameters can be used in the expression with the values given as additional keyword arguments.
    Expressions should be relative to the lab frame.
    """
    _expression_fields: ClassVar[tuple[str, ...]] = (
        "Ex_expression", "Ey_expression", "Ez_expression",
        "Bx_expression", "By_expression", "Bz_expression",
    )

    Ex_expression: Expression | None = Field(
        default=None,
        description="Analytic expression describing Ex field [V/m]"
    )
    Ey_expression: Expression | None = Field(
        default=None,
        description="Analytic expression describing Ey field [V/m]"
    )
    Ez_expression: Expression | None = Field(
        default=None,
        description="Analytic expression describing Ez field [V/m]"
    )
    Bx_expression: Expression | None = Field(
        default=None,
        description="Analytic expression describing Bx field [T]"
    )
    By_expression: Expression | None = Field(
        default=None,
        description="Analytic expression describing By field [T]"
    )
    Bz_expression: Expression | None = Field(
        default=None,
        description="Analytic expression describing Bz field [T]"
    )
    lower_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Lower bound of the region where the field is applied [m]."
    )
    upper_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Upper bound of the region where the field is applied [m]"
    )
    user_defined_kw: dict = Field(
        default_factory=dict,
        description="Constants referenced in the expressions, collected from otherwise-unrecognized keyword arguments."
    )


class PICMI_Mirror(PICMI_AppliedField):
    """
    Describes a perfectly reflecting mirror, where the E and B fields are zeroed
    out in a plane of finite thickness.

    Only one of the [x,y,z]_front_location should be specified. The mirror will be set
    perpendicular to the respective direction and infinite in the others.
    The depth of the mirror will be the maximum of the specified depth and number_of_cells,
    or the code's default value if neither are specified.
    """
    x_front_location: float | None = Field(
        default=None,
        description="Location in x of the front of the mirror [m]"
    )
    y_front_location: float | None = Field(
        default=None,
        description="Location in y of the front of the mirror [m]"
    )
    z_front_location: float | None = Field(
        default=None,
        description="Location in z of the front of the mirror [m]"
    )
    depth: float | None = Field(
        default=None,
        description="Depth of the mirror [m]"
    )
    number_of_cells: int | None = Field(
        default=None,
        description="Minimum number of cells zeroed out"
    )

    @model_validator(mode="after")
    def _one_front_location(self) -> Self:
        assert [self.x_front_location, self.y_front_location, self.z_front_location].count(None) == 2, (
            "At least one and only one of [x,y,z]_front_location should be specified."
        )
        return self


class PICMI_LoadAppliedField(PICMI_AppliedField):
    """
    The E and B fields read from file are applied to the particles directly. (They are not affected by the field solver.)
    The expected format is the file is OpenPMD with axes (x,y,z) in Cartesian, or (r,z) in Cylindrical geometry.
    """
    read_fields_from_path: str = Field(
        description="Path to file with field data"
    )
    load_B: bool = Field(
        default=True,
        description="If False, do not load magnetic field"
    )
    load_E: bool = Field(
        default=True,
        description="If False, do not load electric field"
    )


class PICMI_LoadGriddedField(PICMI_AppliedField):
    """
    The data read in is used to initialize the E and B fields on the grid at the start of the simulation.
    The expected format is the file is OpenPMD with axes (x,y,z) in Cartesian, or (r,z) in Cylindrical geometry.
    """
    read_fields_from_path: str = Field(
        description="Path to file with field data"
    )
    load_B: bool = Field(
        default=True,
        description="If False, do not load magnetic field"
    )
    load_E: bool = Field(
        default=True,
        description="If False, do not load electric field"
    )


PICMI_AnyAppliedField = (
    PICMI_ConstantAppliedField
    | PICMI_AnalyticAppliedField
    | PICMI_Mirror
    | PICMI_LoadAppliedField
    | PICMI_LoadGriddedField
    | PICMI_AppliedField
)
