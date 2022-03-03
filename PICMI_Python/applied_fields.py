"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import re

from .base import _ClassWithInit

# ---------------
# Applied fields
# ---------------


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
    def __init__(self, Ex=None, Ey=None, Ez=None, Bx=None, By=None, Bz=None,
                 lower_bound=[None,None,None], upper_bound=[None,None,None],
                 **kw):

        self.Ex = Ex
        self.Ey = Ey
        self.Ez = Ez
        self.Bx = Bx
        self.By = By
        self.Bz = Bz

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

        self.handle_init(kw)


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
    def __init__(self, Ex_expression=None, Ey_expression=None, Ez_expression=None,
                       Bx_expression=None, By_expression=None, Bz_expression=None,
                 lower_bound=[None,None,None], upper_bound=[None,None,None],
                 **kw):

        self.Ex_expression = Ex_expression
        self.Ey_expression = Ey_expression
        self.Ez_expression = Ez_expression
        self.Bx_expression = Bx_expression
        self.By_expression = By_expression
        self.Bz_expression = Bz_expression

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

        # --- Find any user defined keywords in the kw dictionary.
        # --- Save them and delete them from kw.
        # --- It's up to the code to make sure that all parameters
        # --- used in the expression are defined.
        self.user_defined_kw = {}
        for k in list(kw.keys()):
            if ((self.Ex_expression is not None and re.search(r'\b%s\b'%k, self.Ex_expression)) or
                (self.Ey_expression is not None and re.search(r'\b%s\b'%k, self.Ey_expression)) or
                (self.Ez_expression is not None and re.search(r'\b%s\b'%k, self.Ez_expression)) or
                (self.Bx_expression is not None and re.search(r'\b%s\b'%k, self.Bx_expression)) or
                (self.By_expression is not None and re.search(r'\b%s\b'%k, self.By_expression)) or
                (self.Bz_expression is not None and re.search(r'\b%s\b'%k, self.Bz_expression))):
                self.user_defined_kw[k] = kw[k]
                del kw[k]

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

    def __init__(self, x_front_location=None, y_front_location=None, z_front_location=None,
                 depth=None, number_of_cells=None, **kw):

        assert [x_front_location,y_front_location,z_front_location].count(None) == 2,\
               Exception('At least one and only one of [x,y,z]_front_location should be specified.')

        self.x_front_location = x_front_location
        self.y_front_location = y_front_location
        self.z_front_location = z_front_location
        self.depth = depth
        self.number_of_cells = number_of_cells

        self.handle_init(kw)

