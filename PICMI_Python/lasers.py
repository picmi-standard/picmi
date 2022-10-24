"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys
import re

from .base import _ClassWithInit, _get_constants

# ---------------
# Physics objects
# ---------------

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

    Parameters
    ----------
    wavelength: float
        Laser wavelength [m], defined as :math:`\\lambda_0` in the above formula

    waist: float
        Waist of the Gaussian pulse at focus [m], defined as :math:`w_0` in the above formula

    duration: float
        Duration of the Gaussian pulse [s], defined as :math:`\\tau` in the above formula

    propagation_direction: unit vector of floats
        Direction of propagation [1]

    polarization_direction: unit vector of floats
        Direction of polarization [1]

    focal_position: vector of floats
        Position of the laser focus [m]

    centroid_position: vector of floats
        Position of the laser centroid at time 0 [m]

    a0: float
        Normalized vector potential at focus
        Specify either a0 or E0 (E0 takes precedence).

    E0: float
        Maximum amplitude of the laser field [V/m]
        Specify either a0 or E0 (E0 takes precedence).

    phi0: float
        Carrier envelope phase (CEP) [rad]

    zeta: float
        Spatial chirp at focus (in the lab frame) [m.s]

    beta: float
        Angular dispersion at focus (in the lab frame) [rad.s]

    phi2: float
        Temporal chirp at focus (in the lab frame) [s^2]

    fill_in: bool, default=True
        Flags whether to fill in the empty spaced opened up when the grid moves

    name: string, optional
        Optional name of the laser
    """
    def __init__(self, wavelength, waist, duration,
                 propagation_direction,
                 polarization_direction,
                 focal_position,
                 centroid_position,
                 a0 = None,
                 E0 = None,
                 phi0 = None,
                 zeta = None,
                 beta = None,
                 phi2 = None,
                 name = None,
                 fill_in = True,
                 **kw):

        assert E0 is not None or a0 is not None, 'One of E0 or a0 must be speficied'

        k0 = 2.*math.pi/wavelength
        if E0 is None:
            E0 = a0*_get_constants().m_e*_get_constants().c**2*k0/_get_constants().q_e
        if a0 is None:
            a0 = E0/(_get_constants().m_e*_get_constants().c**2*k0/_get_constants().q_e)

        self.wavelength = wavelength
        self.k0 = k0
        self.waist = waist
        self.duration = duration
        self.focal_position = focal_position
        self.centroid_position = centroid_position
        self.propagation_direction = propagation_direction
        self.polarization_direction = polarization_direction
        self.a0 = a0
        self.E0 = E0
        self.phi0 = phi0
        self.zeta = zeta
        self.beta = beta
        self.phi2 = phi2
        self.name = name
        self.fill_in = fill_in

        self.handle_init(kw)


class PICMI_AnalyticLaser(_ClassWithInit):
    """
    Specifies a laser with an analytically described distribution

    Parameters
    ----------
    name=None: string, optional
        Optional name of the laser

    field_expression: string
        Analytic expression describing the electric field of the laser [V/m]
        Expression should be in terms of the position, 'X', 'Y', in the plane orthogonal
        to the propagation direction, and 't' the time. The expression should describe
        the full field, including the oscillitory component.
        Parameters can be used in the expression with the values given as keyword arguments.

    wavelength: float
        Laser wavelength.
        This should be built into the expression, but some codes require a specified value for numerical purposes.

    propagation_direction: unit vector of floats
        Direction of propagation [1]

    polarization_direction: unit vector of floats
        Direction of polarization [1]

    amax: float, optional
        Maximum normalized vector potential.
        Specify either amax or Emax (Emax takes precedence).
        This should be built into the expression, but some codes require a specified value for numerical purposes.

    Emax: float, optional
        Maximum amplitude of the laser field [V/m].
        Specify either amax or Emax (Emax takes precedence).
        This should be built into the expression, but some codes require a specified value for numerical purposes.

    fill_in: bool, default=True
        Flags whether to fill in the empty spaced opened up when the grid moves
    """
    def __init__(self, field_expression,
                 wavelength,
                 propagation_direction,
                 polarization_direction,
                 amax = None,
                 Emax = None,
                 name = None,
                 fill_in = True,
                 **kw):

        assert Emax is not None or amax is not None, 'One of Emax or amax must be speficied'

        k0 = 2.*math.pi/wavelength
        if Emax is None:
            Emax = amax*_get_constants().m_e*_get_constants().c**2*k0/_get_constants().q_e
        if amax is None:
            amax = Emax/(_get_constants().m_e*_get_constants().c**2*k0/_get_constants().q_e)

        self.wavelength = wavelength
        self.field_expression = field_expression
        self.k0 = k0
        self.propagation_direction = propagation_direction
        self.polarization_direction = polarization_direction
        self.amax = amax
        self.Emax = Emax
        self.name = name
        self.fill_in = fill_in

        self.field_expression = '{}'.format(field_expression).replace('\n', '')

        # --- Find any user defined keywords in the kw dictionary.
        # --- Save them and delete them from kw.
        # --- It's up to the code to make sure that all parameters
        # --- used in the expression are defined.
        self.user_defined_kw = {}
        for k in list(kw.keys()):
            if re.search(r'\b%s\b'%k, self.field_expression):
                self.user_defined_kw[k] = kw[k]
                del kw[k]

        self.handle_init(kw)


# ------------------
# Numeric Objects
# ------------------


class PICMI_LaserAntenna(_ClassWithInit):
    """
    Specifies the laser antenna injection method

    Parameters
    ----------
    position: vector of strings
        Position of antenna launching the laser [m]

    normal_vector: vector of strings, optional
        Vector normal to antenna plane, defaults to the laser direction
        of propagation [1]
    """
    def __init__(self, position, normal_vector=None, **kw):

        self.position = position
        self.normal_vector = normal_vector

        self.handle_init(kw)
