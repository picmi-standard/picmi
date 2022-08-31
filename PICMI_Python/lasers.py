"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys
import re
import typing

from autoclass import autoargs
from typeguard import typechecked

from .base import _ClassWithInit, _get_constants
from . import picmi_types

# ---------------
# Physics objects
# ---------------

@typechecked
class PICMI_GaussianLaser(_ClassWithInit):
    """
    Specifies a Gaussian laser distribution.

    More precisely, the electric field **near the focal plane** is given by:

    .. math::

        E(\\boldsymbol{x},t) = a_0\\times E_0\,
        \exp\left( -\\frac{r^2}{w_0^2} - \\frac{(z-z_0-ct)^2}{c^2\\tau^2} \\right)
        \cos[ k_0( z - z_0 - ct ) - \phi_{cep} ]

    where :math:`k_0 = 2\pi/\\lambda_0` is the wavevector and where
    :math:`E_0 = m_e c^2 k_0 / q_e` is the field amplitude for :math:`a_0=1`.

    .. note::

        The additional terms that arise **far from the focal plane**
        (Gouy phase, wavefront curvature, ...) are not included in the above
        formula for simplicity, but are of course taken into account by
        the code, when initializing the laser pulse away from the focal plane.

    Parameters:
    -----------
      - wavelength: Laser wavelength [m], defined as :math:`\\lambda_0` in the above formula
      - waist: Waist of the Gaussian pulse at focus [m], defined as :math:`w_0` in the above formula
      - duration: Duration of the Gaussian pulse [s], defined as :math:`\\tau` in the above formula
      - focal_position=[0,0,0]: Position of the laser focus (vector) [m]
      - centroid_position=[0,0,0]: Position of the laser centroid at time 0 (vector) [m]
      - propagation_direction=[0,0,1]: Direction of propagation (unit vector) [1]
      - polarization_direction=[1,0,0]: Direction of polarization (unit vector) [1]
      - a0: Normalized vector potential at focus
            Specify either a0 or E0 (E0 takes precedence).
      - E0: Maximum amplitude of the laser field [V/m]
            Specify either a0 or E0 (E0 takes precedence).
      - phi0: Carrier envelope phase (CEP) [rad]
      - zeta: Spatial chirp at focus (in the lab frame) [m.s]
      - beta: Angular dispersion at focus (in the lab frame) [rad.s]
      - phi2: Temporal chirp at focus (in the lab frame) [s^2]
      - fill_in=True: Flags whether to fill in the empty spaced opened up when the grid moves
      - name=None: Optional name of the laser
    """
    @autoargs(exclude=['kw'])
    def __init__(self, wavelength : float,
                       waist : float,
                       duration : float,
                       focal_position : picmi_types.VectorFloat3 = [0., 0., 0.],
                       centroid_position : picmi_types.VectorFloat3 = [0., 0., 0.],
                       propagation_direction : picmi_types.VectorFloat3 = [0., 0., 1.],
                       polarization_direction : picmi_types.VectorFloat3 = [1., 0., 0.],
                       a0 : float = None,
                       E0 : float = None,
                       phi0 : float = None,
                       zeta : float = None,
                       beta : float = None,
                       phi2 : float = None,
                       fill_in : bool = True,
                       name : str = None,
                       **kw):

        assert E0 is not None or a0 is not None, 'One of E0 or a0 must be speficied'

        k0 = 2.*math.pi/wavelength
        if E0 is None:
            self.E0 = a0*_get_constants().m_e*_get_constants().c**2*k0/_get_constants().q_e
        if a0 is None:
            self.a0 = E0/(_get_constants().m_e*_get_constants().c**2*k0/_get_constants().q_e)

        self.handle_init(kw)


@typechecked
class PICMI_AnalyticLaser(_ClassWithInit):
    """
    Specifies a laser with an analytically described distribution
      - name=None: Optional name of the laser
      - field_expression: Analytic expression describing the electric field of the laser(string) [V/m]
                            Expression should be in terms of the position, 'X', 'Y', in the plane orthogonal
                            to the propagation direction, and 't' the time. The expression should describe
                            the full field, including the oscillitory component.
                            Parameters can be used in the expression with the values given as keyword arguments.
      - propagation_direction=[0,0,1]: Direction of propagation (unit vector) [1]
      - polarization_direction=[1,0,0]: Direction of polarization (unit vector) [1]

      Even though the parameters below should be built into the expression, some codes require
      specified values for numerical purposes.

      - wavelength: Laser wavelength
      - amax: Maximum normalized vector potential
      - Emax: Maximum amplitude of the laser field [V/m]
      Specify either amax or Emax (Emax takes precedence).
      - fill_in=True: Flags whether to fill in the empty spaced opened up when the grid moves
    """
    @autoargs(exclude=['kw'])
    def __init__(self, field_expression : picmi_types.Expression,
                       wavelength : float,
                       propagation_direction : picmi_types.VectorFloat3 = [0., 0., 1.],
                       polarization_direction : picmi_types.VectorFloat3 = [1., 0., 0.],
                       amax : float = None,
                       Emax : float = None,
                       name : str = None,
                       fill_in : bool = True,
                       **kw):

        assert Emax is not None or amax is not None, 'One of Emax or amax must be speficied'

        k0 = 2.*math.pi/wavelength
        if Emax is None:
            self.Emax = amax*_get_constants().m_e*_get_constants().c**2*k0/_get_constants().q_e
        if amax is None:
            self.amax = Emax/(_get_constants().m_e*_get_constants().c**2*k0/_get_constants().q_e)

        self.handle_init(kw)


# ------------------
# Numeric Objects
# ------------------


@typechecked
class PICMI_LaserAntenna(_ClassWithInit):
    """
    Specifies the laser antenna injection method
      - position: Position of antenna launching the laser (vector) [m]
      - normal_vector: Vector normal to antenna plane (vector) [1]
    """
    @autoargs(exclude=['kw'])
    def __init__(self, position : picmi_types.VectorFloat3,
                       normal_vector : picmi_types.VectorFloat3,
                       **kw):

        self.handle_init(kw)
