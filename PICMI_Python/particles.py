"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in the file are all particle related
"""

import numbers
from typing import ClassVar, Literal, Self

import numpy as np
from pydantic import Field, PrivateAttr, field_validator, model_validator

from .base import (
    Expression,
    PICMI_Distribution,
    PICMI_Layout,
    _PICMIModel,
    PICMI_ExpressionParameters,
    broadcast_validation,
    resolve_once,
    with_mutually_exclusive,
)
from .fields import PICMI_AnyGrid
from .interactions import PICMI_AnyInteraction, PICMI_FieldIonization

# ---------------
# Physics objects
# ---------------


class PICMI_GaussianBunchDistribution(PICMI_Distribution):
    """
    Describes a Gaussian distribution of particles
    """
    n_physical_particles: float = Field(ge=0.,
        description="Number of physical particles in the bunch"
    )
    rms_bunch_size: list[float] = Field(
        min_length=3,
        max_length=3,
        description="RMS bunch size at t=0 [m]"
    )
    rms_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        min_length=3,
        max_length=3,
        description="RMS velocity spread at t=0 [m/s]"
    )
    centroid_position: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        min_length=3,
        max_length=3,
        description="Position of the bunch centroid at t=0 [m]"
    )
    centroid_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        min_length=3,
        max_length=3,
        description="Velocity (gamma*V) of the bunch centroid at t=0 [m/s]"
    )
    velocity_divergence: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        min_length=3,
        max_length=3,
        description="Expansion rate of the bunch at t=0 [m/s/m]"
    )


class PICMI_UniformDistribution(PICMI_Distribution):
    """
    Describes a uniform density distribution of particles
    """
    density: float = Field(ge=0.,
        description="Physical number density [m^-3]"
    )
    lower_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Lower bound of the distribution [m]"
    )
    upper_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Upper bound of the distribution [m]"
    )
    rms_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Thermal velocity spread [m/s]"
    )
    directed_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Directed, average, proper velocity [m/s]"
    )
    fill_in: bool | None = Field(
        default=None,
        description="Flags whether to fill in the empty spaced opened up when the grid moves"
    )


class PICMI_FoilDistribution(PICMI_Distribution):
    """
    Describes a foil with optional exponential pre- and post-plasma ramps along the propagation direction.
    """
    density: float = Field(ge=0.,
        description="Physical number density [m^-3]"
    )
    front: float = Field(
        description="Position of front surface of foil [m]"
    )
    thickness: float = Field(ge=0.,
        description="Thickness of the foil [m]"
    )
    exponential_pre_plasma_length: float | None = Field(gt=0.,
        default=None,
        description="Length scale of exponential decay of pre-foil plasma density [m]"
    )
    exponential_pre_plasma_cutoff: float | None = Field(ge=0.,
        default=None,
        description="Cutoff length for exponential decay of pre-foil density [m]"
    )
    exponential_post_plasma_length: float | None = Field(gt=0.,
        default=None,
        description="Length scale of exponential decay of post-foil plasma density [m]"
    )
    exponential_post_plasma_cutoff: float | None = Field(ge=0.,
        default=None,
        description="Cutoff length for exponential decay of post-foil density [m]"
    )
    lower_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        min_length=3,
        max_length=3,
        description="Lower bound of the distribution [m]"
    )
    upper_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        min_length=3,
        max_length=3,
        description="Upper bound of the distribution [m]"
    )
    rms_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        min_length=3,
        max_length=3,
        description="Thermal velocity spread [m/s]"
    )
    directed_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        min_length=3,
        max_length=3,
        description="Directed, average, proper velocity [m/s]"
    )
    fill_in: bool | None = Field(
        default=None,
        description="Flags whether to fill in the empty spaced opened up when the grid moves"
    )

class PICMI_AnalyticFluxDistribution(PICMI_Distribution, PICMI_ExpressionParameters):
    """
    Describes a flux of particles emitted from a plane

    Parameters can be used in the flux expression, with their values given as keyword arguments.
    """
    _expression_fields: ClassVar[tuple[str, ...]] = ("flux",)

    flux: Expression = Field(
        description="Analytic expression describing flux of particles [m^-2.s^-1]. Expression should be in terms of the position and time, written as 'x', 'y', 'z', and 't'."
    )
    flux_normal_axis: str = Field(
        description="x, y, or z for 3D, x or z for 2D, or r, t, or z in RZ geometry"
    )
    surface_flux_position: float = Field(
        description="Location of the injection plane [m] along the direction specified by `flux_normal_axis`"
    )
    flux_direction: Literal[-1, 1] = Field(
        description="Direction of the flux relative to the plane: -1 or +1"
    )
    lower_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        min_length=3,
        max_length=3,
        description="Lower bound of the distribution [m]"
    )
    upper_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        min_length=3,
        max_length=3,
        description="Upper bound of the distribution [m]"
    )
    rms_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        min_length=3,
        max_length=3,
        description="Thermal velocity spread [m/s]"
    )
    directed_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        min_length=3,
        max_length=3,
        description="Directed, average, proper velocity [m/s]"
    )
    flux_tmin: float | None = Field(
        default=None,
        description="Time at which the flux injection will be turned on."
    )
    flux_tmax: float | None = Field(
        default=None,
        description="Time at which the flux injection will be turned off."
    )
    gaussian_flux_momentum_distribution: bool | None = Field(
        default=None,
        description="If True, the momentum distribution is v*Gaussian, in the direction normal to the plane. Otherwise, the momentum distribution is simply Gaussian."
    )
    user_defined_kw: dict = Field(
        default_factory=dict,
        description="Constants referenced in the flux expression, collected from otherwise-unrecognized keyword arguments."
    )

PICMI_UniformFluxDistribution = PICMI_AnalyticFluxDistribution

class PICMI_AnalyticDistribution(PICMI_Distribution, PICMI_ExpressionParameters):
    """
    Describes a plasma with density following a provided analytic expression

    Parameters can be used in the expressions, with their values given as keyword arguments.
    For example, this creates a distribution where the density is ``n0`` below ``rmax`` and
    zero elsewhere:

    .. code-block:: python

        dist = AnalyticDistribution(density_expression='((x**2+y**2)<rmax**2)*n0',
                                    rmax = 1.,
                                    n0 = 1.e20,
                                    ...)
    """
    _expression_fields: ClassVar[tuple[str, ...]] = (
        "density_expression", "momentum_expressions", "momentum_spread_expressions"
    )

    density_expression: Expression = Field(
        description="Analytic expression describing physical number density [m^-3]. Expression should be in terms of the position, written as 'x', 'y', and 'z'. Parameters can be used in the expression with the values given as keyword arguments."
    )
    momentum_expressions: list[Expression | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Analytic expressions describing the gamma*velocity for each axis [m/s]. Expressions should be in terms of the position, written as 'x', 'y', and 'z'. For any axis not supplied (set to None), directed_velocity will be used."
    )
    momentum_spread_expressions: list[Expression | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Analytic expressions describing the gamma*velocity Gaussian thermal spread sigma for each axis [m/s]. For any axis not supplied (set to None), zero will be used."
    )
    lower_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Lower bound of the distribution [m]"
    )
    upper_bound: list[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Upper bound of the distribution [m]"
    )
    rms_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Thermal velocity spread [m/s]"
    )
    directed_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Directed, average, proper velocity [m/s]"
    )
    fill_in: bool | None = Field(
        default=None,
        description="Flags whether to fill in the empty spaced opened up when the grid moves"
    )
    user_defined_kw: dict = Field(
        default_factory=dict,
        description="Constants referenced in the analytic expressions, collected from otherwise-unrecognized keyword arguments."
    )


class PICMI_ParticleListDistribution(PICMI_Distribution):
    """
    Load particles at the specified positions and velocities

    The positions and velocities can be given as lists, or as a single value that is used for
    all particles. All lists must have the same length.
    """
    x: list[float] = Field(
        default_factory=lambda: [0.],
        description="List of x positions of the particles [m]"
    )
    y: list[float] = Field(
        default_factory=lambda: [0.],
        description="List of y positions of the particles [m]"
    )
    z: list[float] = Field(
        default_factory=lambda: [0.],
        description="List of z positions of the particles [m]"
    )
    ux: list[float] = Field(
        default_factory=lambda: [0.],
        description="List of ux of the particles (ux = gamma*vx) [m/s]"
    )
    uy: list[float] = Field(
        default_factory=lambda: [0.],
        description="List of uy of the particles (uy = gamma*vy) [m/s]"
    )
    uz: list[float] = Field(
        default_factory=lambda: [0.],
        description="List of uz of the particles (uz = gamma*vz) [m/s]"
    )
    weight: float | list[float] = Field(
        default=0.,
        description="Particle weight or list of weights, number of real particles per simulation particle"
    )

    _per_particle_fields: ClassVar[tuple[str, ...]] = ("x", "y", "z", "ux", "uy", "uz")

    @field_validator(*_per_particle_fields, mode="before")
    @classmethod
    def _as_list(cls, value):
        # a single value (or an array) becomes a list
        return np.atleast_1d(value).tolist()

    @field_validator("weight", mode="before")
    @classmethod
    def _weight_as_float_or_list(cls, value):
        # note that the weight can be a scalar
        return value if np.ndim(value) == 0 else np.atleast_1d(value).tolist()

    @model_validator(mode="after")
    @resolve_once
    def _broadcast_to_the_number_of_particles(self) -> Self:
        lengths = {name: len(getattr(self, name)) for name in self._per_particle_fields}
        lengths["weight"] = np.size(self.weight)
        number_of_particles = max(lengths.values())
        for name, length in lengths.items():
            assert length in (number_of_particles, 1), f"Length of {name} doesn't match len of others"
        for name in self._per_particle_fields:
            if lengths[name] == 1 and number_of_particles > 1:
                setattr(self, name, getattr(self, name) * number_of_particles)
        return self


class PICMI_FromFileDistribution(PICMI_Distribution):
    """
    Load particles from an openPMD file.

    The openPMD file must contain the attributes `position`, `momentum`, `weighting`.
    """
    file_path: str = Field(
        description="Path to the openPMD file"
    )


PICMI_AnyDistribution = (
    PICMI_GaussianBunchDistribution
    | PICMI_UniformDistribution
    | PICMI_FoilDistribution
    | PICMI_AnalyticFluxDistribution
    | PICMI_AnalyticDistribution
    | PICMI_ParticleListDistribution
    | PICMI_FromFileDistribution
    | PICMI_Distribution
)


# ------------------
# Numeric Objects
# ------------------


class PICMI_ParticleDistributionPlanarInjector(_PICMIModel):
    """
    Describes the injection of particles from a plane
    """
    position: list[float] = Field(
        min_length=3,
        max_length=3,
        description="Position of the particle centroid [m]"
    )
    plane_normal: list[float] = Field(
        min_length=3,
        max_length=3,
        description="Vector normal to the plane of injection [1]"
    )
    plane_velocity: list[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        min_length=3,
        max_length=3,
        description="Velocity of the plane of injection [m/s]"
    )
    method: Literal["InPlace", "Plane"] = Field(
        default="InPlace",
        description="Method of injection"
    )


class PICMI_GriddedLayout(PICMI_Layout):
    """
    Specifies a gridded layout of particles
    """
    n_macroparticles_per_cell: list[int] = Field(
        min_length=1, max_length=3,
        description="Number of particles per cell along each axis (one entry per grid dimension, or a single number in 1D)"
    )
    grid: PICMI_AnyGrid | None = Field(
        default=None,
        description="Grid object specifying the grid to follow. If not specified, the underlying grid of the code is used."
    )

    def __init__(self, *args, **kwargs):
        if "n_macroparticle_per_cell" in kwargs and "n_macroparticles_per_cell" in kwargs:
            raise ValueError(
                f"You have given {kwargs['n_macroparticles_per_cell']=} and {kwargs['n_macroparticle_per_cell']=}. "
                    "Please only provide the former."
            )
        # Only translate the deprecated spelling if it was given, so that omitting both
        # reports the missing required argument instead of an invalid ``None``.
        if "n_macroparticle_per_cell" in kwargs:
            kwargs["n_macroparticles_per_cell"] = kwargs.pop("n_macroparticle_per_cell")
        return super().__init__(*args, **kwargs)

    @field_validator("n_macroparticles_per_cell", mode="before")
    @classmethod
    def _single_number_as_list(cls, value):
        # e.g., in 1D
        if isinstance(value, numbers.Number) and not isinstance(value, bool):
            return [value]
        return value

    # A validator of the field, instead of one in an Annotated type, which the documentation
    # would show as the type.
    @field_validator("n_macroparticles_per_cell", mode="after")
    @classmethod
    def _not_negative(cls, value):
        return broadcast_validation(
            value, condition=lambda v: v >= 0,
            message="All n_macroparticle_per_cell must be greater than or equal to 0.",
        )

    @property
    def n_macroparticle_per_cell(self):
        return self.n_macroparticles_per_cell

    @n_macroparticle_per_cell.setter
    def n_macroparticle_per_cell(self, value):
        self.n_macroparticles_per_cell = value


@with_mutually_exclusive("n_macroparticles_per_cell", "n_macroparticles", required=True)
class PICMI_PseudoRandomLayout(PICMI_Layout):
    """
    Specifies a pseudo-random layout of the particles
    """
    n_macroparticles: int | None = Field(ge=0,
        default=None,
        description="Total number of macroparticles to load. Either this argument or n_macroparticles_per_cell should be supplied (not both)."
    )
    n_macroparticles_per_cell: int | None = Field(ge=0,
        default=None,
        description="Number of macroparticles to load per cell. Either this argument or n_macroparticles should be supplied (not both)."
    )
    seed: int | None = Field(
        default=None,
        description="Pseudo-random number generator seed"
    )
    grid: PICMI_AnyGrid | None = Field(
        default=None,
        description="Grid object specifying the grid to follow for n_macroparticles_per_cell. If not specified, the underlying grid of the code is used."
    )


PICMI_AnyLayout = PICMI_GriddedLayout | PICMI_PseudoRandomLayout | PICMI_Layout


class PICMI_Species(_PICMIModel):
    """
    Sets up the species to be simulated.
    The species charge and mass can be specified by setting the particle type or by setting them directly.
    If the particle type is specified, the charge or mass can be set to override the value from the type.

    The particle advance method options:

    - 'Boris': Standard "leap-frog" Boris advance
    - 'Vay':
    - 'Higuera-Cary':
    - 'Li':
    - 'free-streaming': Advance with no fields
    - 'LLRK4': Landau-Lifschitz radiation reaction formula with RK-4)
    """

    methods_list: ClassVar[list[str]] = ["Boris", "Vay", "Higuera-Cary", "Li", "free-streaming", "LLRK4"]

    particle_type: str | None = Field(
        default=None,
        description="A string specifying an elementary particle, atom, or other, as defined in the openPMD 2 species type extension, openPMD-standard/EXT_SpeciesType.md",
    )
    name: str | None = Field(
        default=None, description="Name of the species. If not specified, it will be determined from the particle type."
    )
    method: str | None = Field(
        default=None,
        description="The particle advance method to use. Code-specific method can be specified using 'other:<method>'. The default is code dependent. Must be one of 'Boris', 'Vay', 'Higuera-Cary', 'Li', 'free-streaming', 'LLRK4', or start with 'other:'",
    )
    charge_state: float | None = Field(
        default=None, description="Charge state of the species (applies only to atoms) [1]"
    )
    charge: float | None = Field(
        default=None, description="Particle charge, if not specified, it will be determined from type [C]"
    )
    mass: float | None = Field(
        default=None, description="Particle mass, if not specified, it will be determined from type [kg]"
    )
    initial_distribution: PICMI_AnyDistribution | list[PICMI_AnyDistribution] | None = Field(
        default=None, description="The initial distribution loaded at t=0. A list of distributions may be given to superimpose several distributions on the same species."
    )
    density_scale: float | None = Field(
        default=None, description="A scale factor on the density given by the initial_distribution."
    )
    particle_shape: Literal["NGP", "linear", "quadratic", "cubic"] | int | None = Field(
        default=None,
        description="Particle shape used for deposition and gather. If not specified, the value from the Simulation object will be used. Other values maybe specified that are code dependent.",
    )
    interactions: list[PICMI_AnyInteraction] = Field(
        default_factory=list, description="List of interactions for this species"
    )

    @field_validator("method")
    @classmethod
    def _validate_method(cls, v):
        if v is not None and v not in PICMI_Species.methods_list and not v.startswith("other:"):
            raise ValueError(
                f'method must start with either "other:", or be one of the following: {", ".join(PICMI_Species.methods_list)}'
            )
        return v


# The species and the interactions reference each other
PICMI_FieldIonization.model_rebuild(_types_namespace={"PICMI_Species": PICMI_Species})
PICMI_Species.model_rebuild(force=True)


_ParticleShape = Literal["NGP", "linear", "quadratic", "cubic"] | int


class PICMI_MultiSpecies(_PICMIModel):
    """
    INCOMPLETE: proportions argument is not implemented
    Multiple species that are initialized with the same distribution.
    Each parameter can be list, giving a value for each species, or a single value which is given to all species.
    The species charge and mass can be specified by setting the particle type or by setting them directly.
    If the particle type is specified, the charge or mass can be set to override the value from the type.

    The species are created at construction, so their parameters cannot be changed afterwards.
    """

    # --- Note to developer: This class attribute needs to be set to the Species class
    # --- defined in the codes PICMI implementation.
    Species_class: ClassVar[type[PICMI_Species] | None] = None

    particle_types: str | list[str | None] | None = Field(
        default=None,
        frozen=True,
        description="A string specifying an elementary particle, atom, or other, as defined in the openPMD 2 species type extension, openPMD-standard/EXT_SpeciesType.md"
    )
    names: str | list[str | None] | None = Field(
        default=None,
        frozen=True,
        description="Names of the species"
    )
    charge_states: float | list[float | None] | None = Field(
        default=None,
        frozen=True,
        description="Charge states of the species (applies only to atoms)"
    )
    charges: float | list[float | None] | None = Field(
        default=None,
        frozen=True,
        description="Particle charges, required when type is not specified, otherwise determined from type [C]"
    )
    masses: float | list[float | None] | None = Field(
        default=None,
        frozen=True,
        description="Particle masses, required when type is not specified, otherwise determined from type [kg]"
    )
    proportions: float | list[float | None] | None = Field(
        default=None,
        frozen=True,
        description="Proportions of the initial distribution made up by each species"
    )
    initial_distribution: PICMI_AnyDistribution | list[PICMI_AnyDistribution] | None = Field(
        default=None,
        frozen=True,
        description="Initial particle distribution, applied to all species"
    )
    particle_shape: _ParticleShape | None = Field(
        default=None,
        frozen=True,
        description="Particle shape used for deposition and gather ('NGP', 'linear', 'quadratic', 'cubic'). If not specified, the value from the `Simulation` object will be used. Other values maybe specified that are code dependent."
    )

    _per_species_fields: ClassVar[tuple[str, ...]] = (
        "particle_types", "names", "charge_states", "charges", "masses", "proportions"
    )
    _species_instances_list: list[PICMI_Species] = PrivateAttr(default_factory=list)
    _species_instances_dict: dict[str, PICMI_Species] = PrivateAttr(default_factory=dict)

    @staticmethod
    def get_input_item(var, i):
        """The value for the i-th species of a parameter given for each species or for all"""
        if var is None or not isinstance(var, list):
            return var
        return var[i]

    @model_validator(mode="after")
    @resolve_once
    def _create_species(self) -> Self:
        if self._species_instances_list:
            # created at construction
            return self

        # The lists give a value per species, single values are given to all species.
        given = [getattr(self, name) for name in self._per_species_fields if getattr(self, name) is not None]
        if not given:
            raise ValueError(f"At least one of {', '.join(self._per_species_fields)} must be specified")
        lengths = {len(var) for var in given if isinstance(var, list)}
        if len(lengths) > 1:
            raise ValueError("All inputs must have the same length")
        nspecies = lengths.pop() if lengths else 1
        if PICMI_MultiSpecies.Species_class is None:
            raise TypeError("The implementing code must set PICMI_MultiSpecies.Species_class")

        for i in range(nspecies):
            name = self.get_input_item(self.names, i)
            species = PICMI_MultiSpecies.Species_class(
                particle_type=self.get_input_item(self.particle_types, i),
                name=name,
                charge=self.get_input_item(self.charges, i),
                charge_state=self.get_input_item(self.charge_states, i),
                mass=self.get_input_item(self.masses, i),
                initial_distribution=self.initial_distribution,
                density_scale=self.get_input_item(self.proportions, i),
                particle_shape=self.particle_shape,
            )
            self._species_instances_list.append(species)
            if name is not None:
                self._species_instances_dict[name] = species
        return self

    @property
    def nspecies(self):
        """Number of species"""
        return len(self._species_instances_list)

    @property
    def species_instances_list(self):
        """The species instances, in order"""
        return self._species_instances_list

    @property
    def species_instances_dict(self):
        """The species instances that have a name, by name"""
        return self._species_instances_dict

    def __len__(self):
        return self.nspecies

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._species_instances_dict[key]
        else:
            return self._species_instances_list[key]


PICMI_AnySpecies = PICMI_Species | PICMI_MultiSpecies
