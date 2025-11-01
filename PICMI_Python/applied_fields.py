"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
from __future__ import annotations
import sys
from typing import Any

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

from collections.abc import Sequence

import re

from pydantic import Field, model_validator

from .base import _ClassWithInit

# ---------------
# Applied fields
# ---------------


class PICMI_ConstantAppliedField(_ClassWithInit):
    """
    Describes a constant applied field.
    """
    Ex: float | None = Field(default=None, description="Constant Ex field [V/m]")
    Ey: float | None = Field(default=None, description="Constant Ey field [V/m]")
    Ez: float | None = Field(default=None, description="Constant Ez field [V/m]")
    Bx: float | None = Field(default=None, description="Constant Bx field [T]")
    By: float | None = Field(default=None, description="Constant By field [T]")
    Bz: float | None = Field(default=None, description="Constant Bz field [T]")
    lower_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Lower bound of the region where the field is applied [m]"
    )
    upper_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Upper bound of the region where the field is applied [m]"
    )


class PICMI_AnalyticAppliedField(_ClassWithInit):
    """
    Describes an analytic applied field.

    The expressions should be in terms of the position and time, written as 'x', 'y', 'z', 't'.
    Parameters can be used in the expression with the values given as additional keyword arguments.
    Expressions should be relative to the lab frame.
    """
    Ex_expression: str | None = Field(default=None, description="Analytic expression describing Ex field [V/m]")
    Ey_expression: str | None = Field(default=None, description="Analytic expression describing Ey field [V/m]")
    Ez_expression: str | None = Field(default=None, description="Analytic expression describing Ez field [V/m]")
    Bx_expression: str | None = Field(default=None, description="Analytic expression describing Bx field [T]")
    By_expression: str | None = Field(default=None, description="Analytic expression describing By field [T]")
    Bz_expression: str | None = Field(default=None, description="Analytic expression describing Bz field [T]")
    lower_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Lower bound of the region where the field is applied [m]"
    )
    upper_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Upper bound of the region where the field is applied [m]"
    )
    user_defined_kw: dict[str, Any] = Field(default_factory=dict, exclude=True, repr=False)
    
    @model_validator(mode='after')
    def _extract_user_defined_kw(self) -> Self:
        """Extract user-defined keywords from expressions"""
        # Extract user-defined keywords from extra fields that appear in expressions
        if self.model_extra:
            for k in list(self.model_extra.keys()):
                if ((self.Ex_expression is not None and re.search(rf'\b{k}\b', self.Ex_expression)) or
                    (self.Ey_expression is not None and re.search(rf'\b{k}\b', self.Ey_expression)) or
                    (self.Ez_expression is not None and re.search(rf'\b{k}\b', self.Ez_expression)) or
                    (self.Bx_expression is not None and re.search(rf'\b{k}\b', self.Bx_expression)) or
                    (self.By_expression is not None and re.search(rf'\b{k}\b', self.By_expression)) or
                    (self.Bz_expression is not None and re.search(rf'\b{k}\b', self.Bz_expression))):
                    self.user_defined_kw[k] = self.model_extra.pop(k)
        return self


class PICMI_Mirror(_ClassWithInit):
    """
    Describes a perfectly reflecting mirror, where the E and B fields are zeroed
    out in a plane of finite thickness.

    Only one of the [x,y,z]_front_location should be specified. The mirror will be set
    perpendicular to the respective direction and infinite in the others.
    The depth of the mirror will be the maximum of the specified depth and number_of_cells,
    or the code's default value if neither are specified.
    """
    x_front_location: float | None = Field(default=None, description="Location in x of the front of the mirror [m]")
    y_front_location: float | None = Field(default=None, description="Location in y of the front of the mirror [m]")
    z_front_location: float | None = Field(default=None, description="Location in z of the front of the mirror [m]")
    depth: float | None = Field(default=None, description="Depth of the mirror [m]")
    number_of_cells: int | None = Field(default=None, description="Minimum number of cells zeroed out")
    
    @model_validator(mode='after')
    def _validate_mirror_location(self) -> Self:
        """Validate that exactly one location is specified"""
        locations = [self.x_front_location, self.y_front_location, self.z_front_location]
        if locations.count(None) != 2:
            raise ValueError('Exactly one of [x,y,z]_front_location should be specified.')
        return self

class PICMI_LoadAppliedField(_ClassWithInit):
    """
    The E and B fields read from file are applied to the particles directly.
    (They are not affected by the field solver.)
    The expected format is the file is OpenPMD with axes (x,y,z) in Cartesian,
    or (r,z) in Cylindrical geometry.
    """
    read_fields_from_path: str = Field(description="Path to file with field data")
    load_B: bool = Field(default=True, description="If False, do not load magnetic field")
    load_E: bool = Field(default=True, description="If False, do not load electric field")


class PICMI_LoadGriddedField(_ClassWithInit):
    """
    The data read in is used to initialize the E and B fields on the grid at the start of the simulation.
    The expected format is the file is OpenPMD with axes (x,y,z) in Cartesian,
    or (r,z) in Cylindrical geometry.
    """
    read_fields_from_path: str = Field(description="Path to file with field data")
    load_B: bool = Field(default=True, description="If False, do not load magnetic field")
    load_E: bool = Field(default=True, description="If False, do not load electric field")
