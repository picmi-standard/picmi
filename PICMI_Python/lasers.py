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

import math
import re

from pydantic import Field, model_validator

from .base import _ClassWithInit, _get_constants

# ---------------
# Physics objects
# ---------------

class PICMI_GaussianLaser(_ClassWithInit):
    r"""
    Specifies a Gaussian laser distribution.

    More precisely, the electric field **near the focal plane** is given by:

    .. math::

        E(\boldsymbol{x},t) = a_0\times E_0\,
        \exp\left( -\frac{r^2}{w_0^2} - \frac{(z-z_0-ct)^2}{c^2\tau^2} \right)
        \cos[ k_0( z - z_0 - ct ) - \phi_{cep} ]

    where :math:`k_0 = 2\pi/\lambda_0` is the wavevector and where
    :math:`E_0 = m_e c^2 k_0 / q_e` is the field amplitude for :math:`a_0=1`.

    .. note::

        The additional terms that arise **far from the focal plane**
        (Gouy phase, wavefront curvature, ...) are not included in the above
        formula for simplicity, but are of course taken into account by
        the code, when initializing the laser pulse away from the focal plane.
    """
    wavelength: float = Field(
        description="Laser wavelength [m], defined as :math:`\\lambda_0` in the above formula"
    )
    waist: float = Field(
        description="Waist of the Gaussian pulse at focus [m], defined as :math:`w_0` in the above formula"
    )
    duration: float = Field(
        description="Duration of the Gaussian pulse [s], defined as :math:`\\tau` in the above formula"
    )
    propagation_direction: Sequence[float] = Field(
        description="Unit vector of length 3. Direction of propagation [1]"
    )
    polarization_direction: Sequence[float] = Field(
        description="Unit vector of length 3. Direction of polarization [1]"
    )
    focal_position: Sequence[float] = Field(
        description="Vector of length 3 of floats. Position of the laser focus [m]"
    )
    centroid_position: Sequence[float] = Field(
        description="Vector of length 3 of floats. Position of the laser centroid at time 0 [m]"
    )
    a0: float | None = Field(
        default=None,
        description="Normalized vector potential at focus. Specify either a0 or E0 (E0 takes precedence)."
    )
    E0: float | None = Field(
        default=None,
        description="Maximum amplitude of the laser field [V/m]. Specify either a0 or E0 (E0 takes precedence)."
    )
    phi0: float | None = Field(
        default=None,
        description="Carrier envelope phase (CEP) [rad]"
    )
    zeta: float | None = Field(
        default=None,
        description="Spatial chirp at focus (in the lab frame) [m.s]"
    )
    beta: float | None = Field(
        default=None,
        description="Angular dispersion at focus (in the lab frame) [rad.s]"
    )
    phi2: float | None = Field(
        default=None,
        description="Temporal chirp at focus (in the lab frame) [s^2]"
    )
    name: str | None = Field(
        default=None,
        description="Optional name of the laser"
    )
    fill_in: bool = Field(default=True, description="Flags whether to fill in the empty spaced opened up when the grid moves")
    k0: float = Field(default=0.0, exclude=True, init_var=False)
    
    @model_validator(mode='after')
    def _compute_a0_e0(self) -> Self:
        """Compute a0 and E0 from each other if needed"""
        if self.E0 is None and self.a0 is None:
            raise ValueError('One of E0 or a0 must be specified')
        
        k0 = 2.*math.pi/self.wavelength
        self.k0 = k0
        
        constants = _get_constants()
        if self.E0 is None:
            self.E0 = self.a0 * constants.m_e * constants.c**2 * k0 / constants.q_e
        if self.a0 is None:
            self.a0 = self.E0 / (constants.m_e * constants.c**2 * k0 / constants.q_e)
        
        return self


class PICMI_AnalyticLaser(_ClassWithInit):
    """
    Specifies a laser with an analytically described distribution.
    """
    field_expression: str = Field(
        description="Analytic expression describing the electric field of the laser [V/m]. Expression should be in terms of the position, 'X', 'Y', in the plane orthogonal to the propagation direction, and 't' the time. The expression should describe the full field, including the oscillitory component. Parameters can be used in the expression with the values given as keyword arguments."
    )
    wavelength: float = Field(
        description="Laser wavelength. This should be built into the expression, but some codes require a specified value for numerical purposes."
    )
    propagation_direction: Sequence[float] = Field(
        description="Unit vector of length 3 of floats. Direction of propagation [1]"
    )
    polarization_direction: Sequence[float] = Field(
        description="Unit vector of length 3 of floats. Direction of polarization [1]"
    )
    amax: float | None = Field(
        default=None,
        description="Maximum normalized vector potential. Specify either amax or Emax (Emax takes precedence). This should be built into the expression, but some codes require a specified value for numerical purposes."
    )
    Emax: float | None = Field(
        default=None,
        description="Maximum amplitude of the laser field [V/m]. Specify either amax or Emax (Emax takes precedence). This should be built into the expression, but some codes require a specified value for numerical purposes."
    )
    name: str | None = Field(
        default=None,
        description="Optional name of the laser"
    )
    fill_in: bool = Field(
        default=True,
        description="Flags whether to fill in the empty spaced opened up when the grid moves"
    )
    k0: float = Field(default=0.0, exclude=True, init_var=False)
    user_defined_kw: dict[str, Any] = Field(default_factory=dict, exclude=True, repr=False)
    
    @model_validator(mode='after')
    def _compute_amax_emax(self) -> Self:
        """Compute amax and Emax from each other if needed, and extract user-defined keywords"""
        if self.Emax is None and self.amax is None:
            raise ValueError('One of Emax or amax must be specified')
        
        # Normalize field expression
        self.field_expression = str(self.field_expression).replace('\n', '')
        
        k0 = 2.*math.pi/self.wavelength
        self.k0 = k0
        
        constants = _get_constants()
        if self.Emax is None:
            self.Emax = self.amax * constants.m_e * constants.c**2 * k0 / constants.q_e
        if self.amax is None:
            self.amax = self.Emax / (constants.m_e * constants.c**2 * k0 / constants.q_e)
        
        # Extract user-defined keywords from extra fields
        if self.model_extra:
            for k in list(self.model_extra.keys()):
                if re.search(rf'\b{k}\b', self.field_expression):
                    self.user_defined_kw[k] = self.model_extra.pop(k)
        
        return self


# ------------------
# Numeric Objects
# ------------------


class PICMI_LaserAntenna(_ClassWithInit):
    """
    Specifies the laser antenna injection method.
    """
    position: Sequence[float] = Field(
        description="Vector of floats. Position of antenna launching the laser [m]"
    )
    normal_vector: Sequence[float] | None = Field(
        default=None,
        description="Vector of floats, optional. Vector normal to antenna plane, defaults to the laser direction of propagation [1]"
    )
