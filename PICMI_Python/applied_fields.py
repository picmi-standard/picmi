"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import re
import typing

from autoclass import autoargs
from typeguard import typechecked

from .base import _ClassWithInit
from . import picmi_types

# ---------------
# Applied fields
# ---------------


@typechecked
class PICMI_ConstantAppliedField(_ClassWithInit):
    """
    Describes a constant applied field
      - Ex: Constant Ex field (float) [V/m]
      - Ey: Constant Ey field (float) [V/m]
      - Ez: Constant Ez field (float) [V/m]
      - Bx: Constant Bx field (float) [T]
      - By: Constant By field (float) [T]
      - Bz: Constant Bz field (float) [T]
      - lower_bound=[None,None,None]: Lower bound of the region where the field is applied (vector) [m]
      - upper_bound=[None,None,None]: Upper bound of the region where the field is applied (vector) [m]
    """
    @autoargs(exclude=['kw'])
    def __init__(self, Ex : float = None,
                       Ey : float = None,
                       Ez : float = None,
                       Bx : float = None,
                       By : float = None,
                       Bz : float = None,
                       lower_bound : picmi_types.VectorFloat3 = [None,None,None],
                       upper_bound : picmi_types.VectorFloat3 = [None,None,None],
                       **kw):

        self.handle_init(kw)


@typechecked
class PICMI_AnalyticAppliedField(_ClassWithInit):
    """
    Describes an analytic applied field
      - Ex_expression: Analytic expression describing Ex field (string) [V/m]
      - Ey_expression: Analytic expression describing Ey field (string) [V/m]
      - Ez_expression: Analytic expression describing Ez field (string) [V/m]
      - Bx_expression: Analytic expression describing Bx field (string) [T]
      - By_expression: Analytic expression describing By field (string) [T]
      - Bz_expression: Analytic expression describing Bz field (string) [T]
               Expressions should be in terms of the position and time, written as 'x', 'y', 'z', 't'.
               Parameters can be used in the expression with the values given as additional keyword arguments.
               Expressions should be relative to the lab frame.
      - lower_bound=[None,None,None]: Lower bound of the region where the field is applied (vector) [m]
      - upper_bound=[None,None,None]: Upper bound of the region where the field is applied (vector) [m]
    """
    @autoargs(exclude=['kw'])
    def __init__(self, Ex_expression : picmi_types.Expression = None,
                       Ey_expression : picmi_types.Expression = None,
                       Ez_expression : picmi_types.Expression = None,
                       Bx_expression : picmi_types.Expression = None,
                       By_expression : picmi_types.Expression = None,
                       Bz_expression : picmi_types.Expression = None,
                       lower_bound : picmi_types.VectorFloat3 = [None,None,None],
                       upper_bound : picmi_types.VectorFloat3 = [None,None,None],
                       **kw):

        self.handle_init(kw)


class PICMI_Mirror(_ClassWithInit):
    """
    Describes a perfectly reflecting mirror, where the E and B fields are zeroed
    out in a plane of finite thickness.
      - x_front_location: Location in x of the front of the nirror (float) [m]
      - y_front_location: Location in y of the front of the nirror (float) [m]
      - z_front_location: Location in z of the front of the nirror (float) [m]
      - depth: Depth of the mirror (float) [m]
      - number_of_cells: Minimum numer of cells zeroed out (integer)

    Only one of the [x,y,z]_front_location should be specified. The mirror will be set
    perpendicular to the respective direction and infinite in the others.
    The depth of the mirror will be the maximum of the specified depth and number_of_cells,
    or the code's default value if neither are specified.
    """

    def __init__(self, x_front_location : float = None,
                       y_front_location : float = None,
                       z_front_location : float = None,
                       depth : float = None,
                       number_of_cells : int = None,
                       **kw):

        assert [x_front_location,y_front_location,z_front_location].count(None) == 2,\
               Exception('At least one and only one of [x,y,z]_front_location should be specified.')

        self.handle_init(kw)

