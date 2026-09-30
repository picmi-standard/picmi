"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import math
from typing import ClassVar, Self

from pydantic import Field, PrivateAttr, model_validator

from .base import (
    Expression,
    PICMI_Laser,
    PICMI_LaserInjection,
    _PICMIModel,
    PICMI_ExpressionParameters,
    _get_constants,
    resolve_once,
)


# ---------------
# Physics objects
# ---------------
def _compute_E0_a0(wavelength, original_E0, original_a0, names=("E0", "a0")):
        """Return ``(E0, a0)``, computing the one that is ``None`` from the other one."""
        if original_E0 is None and original_a0 is None:
            raise ValueError(f'One of {names[0]} or {names[1]} must be specified')

        k0 = 2.*math.pi/wavelength
        constants = _get_constants()
        # field amplitude for a0 = 1
        factor = constants.m_e * constants.c**2 * k0 / constants.q_e

        # Note: compare to None explicitly, so that a zero amplitude is not treated as unset.
        E0 = original_a0 * factor if original_E0 is None else original_E0
        a0 = E0 / factor if original_a0 is None else original_a0

        # Both might have been given, so we check for consistency.
        # The relative tolerance is just an arbitrary cutoff for floating point errors.
        # Let's presume that the user did not purposefully choose a value
        # ever so slightly off to derail us here.
        if not math.isclose(E0, a0 * factor, rel_tol=1.e-6):
           raise ValueError(f"You provided inconsistent {original_a0=} and {original_E0=} which resulted in {a0=} and {E0=}.")
        return E0, a0




class _PICMILaser(PICMI_Laser):
    # Base of the lasers that are specified by either a normalized vector potential or a field
    # amplitude, which are computed from each other (without docstring, so that it is not
    # prepended to the ones of the derived classes).

    # (normalized vector potential, field amplitude) field names, defined by each laser
    _amplitude_fields: ClassVar[tuple[str, str]] = ("a0", "E0")

    # (wavelength, normalized vector potential, field amplitude) as of the last resolution
    # of the amplitudes, and the field name of the amplitude that the other one is derived from.
    _resolved_state: tuple | None = PrivateAttr(default=None)
    _amplitude_source: str | None = PrivateAttr(default=None)

    @property
    def k0(self) -> float:
        """Laser wavenumber [1/m], :math:`k_0 = 2\\pi/\\lambda_0`"""
        return 2.*math.pi/self.wavelength

    @model_validator(mode='after')
    @resolve_once
    def _compute_amplitudes(self) -> Self:
        """Compute the normalized vector potential and the field amplitude from each other if needed

        At construction, the amplitude that is not given is computed from the other one;
        if both are given, they must be consistent (the field amplitude then takes precedence
        later on).

        The validation re-runs on every later assignment and when the laser is passed to
        another PICMI object, at which point both amplitudes are set. To keep them
        consistent, the amplitude that was just assigned (or otherwise, the one the other
        was derived from) is kept and the other one is re-derived, e.g., when the wavelength
        changes.
        """
        a_name, E_name = self._amplitude_fields
        a0, E0 = getattr(self, a_name), getattr(self, E_name)
        source = self._amplitude_source
        if self._resolved_state is not None:
            if self._resolved_state == (self.wavelength, a0, E0):
                return self
            _, last_a0, last_E0 = self._resolved_state
            a0_changed, E0_changed = a0 != last_a0, E0 != last_E0
            if a0_changed and not E0_changed:
                # the normalized vector potential was assigned; if it was unset, derive it again
                source = a_name if a0 is not None else E_name
            elif E0_changed and not a0_changed:
                source = E_name if E0 is not None else a_name
            elif a0_changed and E0_changed:
                # both changed at once: treat like new input
                source = None
            if source == a_name:
                E0 = None
            elif source == E_name:
                a0 = None

        E0, a0 = _compute_E0_a0(self.wavelength, E0, a0, names=(E_name, a_name))
        if source is None:
            source = E_name if getattr(self, E_name) is not None else a_name

        setattr(self, E_name, E0)
        setattr(self, a_name, a0)
        self._amplitude_source = source
        self._resolved_state = (self.wavelength, a0, E0)
        return self


class PICMI_GaussianLaser(_PICMILaser):
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
        gt=0.,
        description="Laser wavelength [m], defined as :math:`\\lambda_0` in the above formula"
    )
    waist: float = Field(
        gt=0.,
        description="Waist of the Gaussian pulse at focus [m], defined as :math:`w_0` in the above formula"
    )
    duration: float = Field(
        gt=0.,
        description="Duration of the Gaussian pulse [s], defined as :math:`\\tau` in the above formula"
    )
    propagation_direction: list[float] = Field(
        description="Unit vector of length 3. Direction of propagation [1]"
    )
    polarization_direction: list[float] = Field(
        description="Unit vector of length 3. Direction of polarization [1]"
    )
    focal_position: list[float] = Field(
        description="Vector of length 3 of floats. Position of the laser focus [m]"
    )
    centroid_position: list[float] = Field(
        description="Vector of length 3 of floats. Position of the laser centroid at time 0 [m]"
    )
    a0: float | None = Field(
        default=None,
        description="Normalized vector potential at focus. Specify either a0 or E0 (if both are given, they must be consistent)."
    )
    E0: float | None = Field(
        default=None,
        description="Maximum amplitude of the laser field [V/m]. Specify either a0 or E0 (if both are given, they must be consistent)."
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

    _amplitude_fields: ClassVar[tuple[str, str]] = ("a0", "E0")


class PICMI_AnalyticLaser(_PICMILaser, PICMI_ExpressionParameters):
    """
    Specifies a laser with an analytically described distribution

    Parameters can be used in the field expression, with their values given as keyword arguments.
    """
    _amplitude_fields: ClassVar[tuple[str, str]] = ("amax", "Emax")
    _expression_fields: ClassVar[tuple[str, ...]] = ("field_expression",)

    field_expression: Expression = Field(
        description="Analytic expression describing the electric field of the laser [V/m]. Expression should be in terms of the position, 'X', 'Y', in the plane orthogonal to the propagation direction, and 't' the time. The expression should describe the full field, including the oscillitory component. Parameters can be used in the expression with the values given as keyword arguments."
    )
    wavelength: float = Field(
        gt=0.,
        description="Laser wavelength. This should be built into the expression, but some codes require a specified value for numerical purposes."
    )
    propagation_direction: list[float] = Field(
        description="Unit vector of length 3. Direction of propagation [1]"
    )
    polarization_direction: list[float] = Field(
        description="Unit vector of length 3. Direction of polarization [1]"
    )
    amax: float | None = Field(
        default=None,
        description="Maximum normalized vector potential. Specify either amax or Emax (if both are given, they must be consistent). This should be built into the expression, but some codes require a specified value for numerical purposes."
    )
    Emax: float | None = Field(
        default=None,
        description="Maximum amplitude of the laser field [V/m]. Specify either amax or Emax (if both are given, they must be consistent). This should be built into the expression, but some codes require a specified value for numerical purposes."
    )
    name: str | None = Field(
        default=None,
        description="Optional name of the laser"
    )
    fill_in: bool = Field(
        default=True,
        description="Flags whether to fill in the empty spaced opened up when the grid moves"
    )
    user_defined_kw: dict = Field(
        default_factory=dict,
        description="Constants referenced in the field expression, collected from otherwise-unrecognized keyword arguments."
    )


PICMI_AnyLaser = PICMI_GaussianLaser | PICMI_AnalyticLaser | PICMI_Laser


# ------------------
# Numeric Objects
# ------------------


class PICMI_LaserAntenna(PICMI_LaserInjection):
    """
    Specifies the laser antenna injection method
    """
    position: list[float] = Field(
        description="Vector of length 3. Position of antenna launching the laser [m]"
    )
    normal_vector: list[float] | None = Field(
        default=None,
        description="Vector of length 3. Vector normal to antenna plane, defaults to the laser direction of propagation [1]"
    )


PICMI_AnyLaserInjection = PICMI_LaserAntenna | PICMI_LaserInjection
