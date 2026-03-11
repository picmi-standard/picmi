"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in the file are all particle related
"""
import math
import re
import sys
import numpy as np
from typing import Any, Literal

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

from collections.abc import Sequence

from pydantic import Field, field_validator, model_validator

from .base import _ClassWithInit

from .fields import (
    PICMI_Cartesian1DGrid,
    PICMI_Cartesian2DGrid,
    PICMI_Cartesian3DGrid,
    PICMI_CylindricalGrid,
)

# Type aliases using Python 3.10+ union syntax
PICMI_Grid = (
    PICMI_Cartesian1DGrid
    | PICMI_Cartesian2DGrid
    | PICMI_Cartesian3DGrid
    | PICMI_CylindricalGrid
)
# PICMI_Distribution will be defined at end of file after all distribution classes
# PICMI_Interaction uses forward reference to avoid circular import with interactions.py
PICMI_Interaction = "PICMI_FieldIonization"

# ---------------
# Physics objects
# ---------------


class PICMI_Species(_ClassWithInit):
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
    methods_list: list[str] = ['Boris', 'Vay', 'Higuera-Cary', 'Li', 'free-streaming', 'LLRK4']

    particle_type: str | None = Field(
        default=None,
        description="String, optional. A string specifying an elementary particle, atom, or other, as defined in the openPMD 2 species type extension, openPMD-standard/EXT_SpeciesType.md"
    )
    name: str | None = Field(
        default=None,
        description="String, optional. Name of the species"
    )
    method: str | None = Field(
        default=None,
        description="String, optional. The particle advance method to use. Code-specific method can be specified using 'other:<method>'. The default is code dependent. Must be one of 'Boris', 'Vay', 'Higuera-Cary', 'Li', 'free-streaming', 'LLRK4', or start with 'other:'"
    )
    charge_state: float | None = Field(
        default=None,
        description="Float, optional. Charge state of the species (applies only to atoms) [1]"
    )
    charge: float | None = Field(
        default=None,
        description="Float, optional. Particle charge, required when type is not specified, otherwise determined from type [C]"
    )
    mass: float | None = Field(
        default=None,
        description="Float, optional. Particle mass, required when type is not specified, otherwise determined from type [kg]"
    )
    initial_distribution: "PICMI_GaussianBunchDistribution | PICMI_UniformDistribution | PICMI_FoilDistribution | PICMI_AnalyticFluxDistribution | PICMI_AnalyticDistribution | PICMI_ParticleListDistribution | PICMI_FromFileDistribution | None" = Field(
        default=None,
        description="Distribution instance. The initial distribution loaded at t=0. Must be one of the standard distributions implemented."
    )
    density_scale: float | None = Field(
        default=None,
        description="Float, optional. A scale factor on the density given by the initial_distribution"
    )
    particle_shape: Literal['NGP', 'linear', 'quadratic', 'cubic'] | None = Field(
        default=None,
        description="String, optional. Particle shape used for deposition and gather. If not specified, the value from the Simulation object will be used. Other values maybe specified that are code dependent."
    )
    interactions: "list[PICMI_FieldIonization]" = Field(
        default_factory=list,
        exclude=True,
        description="List of interactions for this species"
    )

    @field_validator('method')
    @classmethod
    def _validate_method(cls, v):
        if v is not None and v not in PICMI_Species.methods_list and not v.startswith('other:'):
            raise ValueError(
                f'method must start with either "other:", or be one of the following: {", ".join(PICMI_Species.methods_list)}'
            )
        return v


class PICMI_MultiSpecies(_ClassWithInit):
    """
    INCOMPLETE: proportions argument is not implemented.
    Multiple species that are initialized with the same distribution.
    Each parameter can be list, giving a value for each species, or a single value which is given to all species.
    The species charge and mass can be specified by setting the particle type or by setting them directly.
    If the particle type is specified, the charge or mass can be set to override the value from the type.
    """
    # --- Note to developer: This class attribute needs to be set to the Species class
    # --- defined in the codes PICMI implementation.
    Species_class: type | None = None

    particle_types: str | Sequence[str] | None = Field(
        default=None,
        description="String or list of strings, optional. A string specifying an elementary particle, atom, or other, as defined in the openPMD 2 species type extension, openPMD-standard/EXT_SpeciesType.md"
    )
    names: str | Sequence[str] | None = Field(
        default=None,
        description="String or list of strings, optional. Names of the species"
    )
    charge_states: float | Sequence[float] | None = Field(
        default=None,
        description="Float or list of floats, optional. Charge states of the species (applies only to atoms)"
    )
    charges: float | Sequence[float] | None = Field(
        default=None,
        description="Float or list of floats, optional. Particle charges, required when type is not specified, otherwise determined from type [C]"
    )
    masses: float | Sequence[float] | None = Field(
        default=None,
        description="Float or list of floats, optional. Particle masses, required when type is not specified, otherwise determined from type [kg]"
    )
    proportions: float | Sequence[float] | None = Field(
        default=None,
        description="Float or list of floats, optional. Proportions of the initial distribution made up by each species"
    )
    initial_distribution: "PICMI_GaussianBunchDistribution | PICMI_UniformDistribution | PICMI_FoilDistribution | PICMI_AnalyticFluxDistribution | PICMI_AnalyticDistribution | PICMI_ParticleListDistribution | PICMI_FromFileDistribution | None" = Field(
        default=None,
        description="Distribution instance. Initial particle distribution, applied to all species"
    )
    particle_shape: Literal['NGP', 'linear', 'quadratic', 'cubic'] | None = Field(
        default=None,
        description="String, optional. Particle shape used for deposition and gather. If not specified, the value from the Simulation object will be used. Other values maybe specified that are code dependent."
    )
    nspecies: int | None = Field(
        default=None,
        exclude=True,
        description="Number of species (computed from input parameters)"
    )
    species_instances_list: list[Any] = Field(
        default_factory=list,
        exclude=True,
        description="List of species instances"
    )
    species_instances_dict: dict[str, Any] = Field(
        default_factory=dict,
        exclude=True,
        description="Dictionary of species instances keyed by name"
    )

    @model_validator(mode='after')
    def _create_species_instances(self) -> Self:
        """Check nspecies consistency and create species instances"""
        # Check nspecies consistency
        self.nspecies = None
        self.check_nspecies(self.particle_types)
        self.check_nspecies(self.names)
        self.check_nspecies(self.charges)
        self.check_nspecies(self.charge_states)
        self.check_nspecies(self.masses)
        self.check_nspecies(self.proportions)

        if self.nspecies is None:
            raise ValueError('At least one of particle_types, names, charges, charge_states, masses, or proportions must be specified')

        # Create the instances of each species
        self.species_instances_list = []
        self.species_instances_dict = {}
        for i in range(self.nspecies):
            particle_type = self.get_input_item(self.particle_types, i)
            name = self.get_input_item(self.names, i)
            charge = self.get_input_item(self.charges, i)
            charge_state = self.get_input_item(self.charge_states, i)
            mass = self.get_input_item(self.masses, i)
            proportion = self.get_input_item(self.proportions, i)
            specie = PICMI_MultiSpecies.Species_class(
                particle_type=particle_type,
                name=name,
                charge=charge,
                charge_state=charge_state,
                mass=mass,
                initial_distribution=self.initial_distribution,
                density_scale=proportion
            )
            self.species_instances_list.append(specie)
            if name is not None:
                self.species_instances_dict[name] = specie

        return self

    def check_nspecies(self, var):
        if var is not None:
            try:
                nvars = len(var)
            except TypeError:
                nvars = 1
            if self.nspecies is not None and self.nspecies != nvars:
                raise ValueError('All inputs must have the same length')
            self.nspecies = nvars

    def get_input_item(self, var, i):
        if var is None:
            return None
        else:
            try:
                len(var)
            except TypeError:
                return var
            else:
                return var[i]

    def __len__(self):
        return self.nspecies if self.nspecies is not None else 0

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.species_instances_dict[key]
        else:
            return self.species_instances_list[key]


class PICMI_GaussianBunchDistribution(_ClassWithInit):
    """
    Describes a Gaussian distribution of particles
    """
    n_physical_particles: int = Field(
        description="Integer. Number of physical particles in the bunch"
    )
    rms_bunch_size: Sequence[float] = Field(
        description="Vector of length 3 of floats. RMS bunch size at t=0 [m]"
    )
    rms_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. RMS velocity spread at t=0 [m/s]"
    )
    centroid_position: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Position of the bunch centroid at t=0 [m]"
    )
    centroid_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Velocity (gamma*V) of the bunch centroid at t=0 [m/s]"
    )
    velocity_divergence: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Expansion rate of the bunch at t=0 [m/s/m]"
    )


class PICMI_UniformDistribution(_ClassWithInit):
    """
    Describes a uniform density distribution of particles
    """
    density: float = Field(
        description="Float. Physical number density [m^-3]"
    )
    lower_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Vector of length 3 of floats, optional. Lower bound of the distribution [m]"
    )
    upper_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Vector of length 3 of floats, optional. Upper bound of the distribution [m]"
    )
    rms_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Thermal velocity spread [m/s]"
    )
    directed_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Directed, average, proper velocity [m/s]"
    )
    fill_in: bool | None = Field(
        default=None,
        description="Bool, optional. Flags whether to fill in the empty spaced opened up when the grid moves"
    )


class PICMI_FoilDistribution(_ClassWithInit):
    """
    Describes a foil with optional exponential pre- and post-plasma ramps along the propagation direction.
    """
    density: float = Field(
        description="Float. Physical number density [m^-3]"
    )
    front: float = Field(
        description="Float. Position of front surface of foil [m]"
    )
    thickness: float = Field(
        description="Float >= 0. Thickness of the foil [m]"
    )
    exponential_pre_plasma_length: float | None = Field(
        default=None,
        description="Float > 0, optional. Length scale of exponential decay of pre-foil plasma density [m]"
    )
    exponential_pre_plasma_cutoff: float | None = Field(
        default=None,
        description="Float >= 0, optional. Cutoff length for exponential decay of pre-foil density [m]"
    )
    exponential_post_plasma_length: float | None = Field(
        default=None,
        description="Float > 0, optional. Length scale of exponential decay of post-foil plasma density [m]"
    )
    exponential_post_plasma_cutoff: float | None = Field(
        default=None,
        description="Float >= 0, optional. Cutoff length for exponential decay of post-foil density [m]"
    )
    lower_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Vector of length 3 of floats, optional. Lower bound of the distribution [m]"
    )
    upper_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Vector of length 3 of floats, optional. Upper bound of the distribution [m]"
    )
    rms_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Thermal velocity spread [m/s]"
    )
    directed_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Directed, average, proper velocity [m/s]"
    )
    fill_in: bool | None = Field(
        default=None,
        description="Bool, optional. Flags whether to fill in the empty spaced opened up when the grid moves"
    )


class PICMI_AnalyticFluxDistribution(_ClassWithInit):
    """
    Describes a flux of particles emitted from a plane.

    The flux expression should be in terms of the position and time, written as 'x', 'y', 'z', and 't'.
    Parameters can be used in the expression with the values given as keyword arguments.
    """
    flux: str = Field(
        description="String. Analytic expression describing flux of particles [m^-2.s^-1]. Expression should be in terms of the position and time, written as 'x', 'y', 'z', and 't'."
    )
    flux_normal_axis: str = Field(
        description="String. x, y, or z for 3D, x or z for 2D, or r, t, or z in RZ geometry"
    )
    surface_flux_position: float = Field(
        description="Float. Location of the injection plane [m] along the direction specified by flux_normal_axis"
    )
    flux_direction: int = Field(
        description="Integer. Direction of the flux relative to the plane: -1 or +1"
    )
    lower_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Vector of floats, optional. Lower bound of the distribution [m]"
    )
    upper_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Vector of floats, optional. Upper bound of the distribution [m]"
    )
    rms_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of floats, default=[0.,0.,0.]. Thermal velocity spread [m/s]"
    )
    directed_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of floats, default=[0.,0.,0.]. Directed, average, proper velocity [m/s]"
    )
    flux_tmin: float | None = Field(
        default=None,
        description="Float, optional. Time at which the flux injection will be turned on"
    )
    flux_tmax: float | None = Field(
        default=None,
        description="Float, optional. Time at which the flux injection will be turned off"
    )
    gaussian_flux_momentum_distribution: bool | None = Field(
        default=None,
        description="Bool, optional. If True, the momentum distribution is v*Gaussian, in the direction normal to the plane. Otherwise, the momentum distribution is simply Gaussian."
    )
    user_defined_kw: dict[str, Any] = Field(
        default_factory=dict,
        exclude=True,
        description="Dictionary of user-defined keyword arguments extracted from extra fields"
    )

    @model_validator(mode='before')
    @classmethod
    def _normalize_flux(cls, data: Any) -> Any:
        """Normalize flux string and extract user-defined kwargs"""
        if not isinstance(data, dict):
            return data
        
        data = dict(data)
        
        # Normalize flux string
        if 'flux' in data:
            data['flux'] = f'{data["flux"]}'.replace('\n', '')
        
        return data
    
    @model_validator(mode='after')
    def _extract_user_defined_kw(self) -> Self:
        """Extract user-defined keywords from extra fields"""
        if self.model_extra:
            flux_str = self.flux
            self.user_defined_kw = {}
            keys_to_remove = []
            for k, v in self.model_extra.items():
                if re.search(rf'\b{k}\b', flux_str):
                    self.user_defined_kw[k] = v
                    keys_to_remove.append(k)
            # Remove extracted keys from model_extra
            for k in keys_to_remove:
                del self.model_extra[k]
        
        return self

PICMI_UniformFluxDistribution = PICMI_AnalyticFluxDistribution

class PICMI_AnalyticDistribution(_ClassWithInit):
    """
    Describes a plasma with density following a provided analytic expression.

    Expressions should be in terms of the position, written as 'x', 'y', and 'z'.
    Parameters can be used in the expression with the values given as keyword arguments.

    This example will create a distribution where the density is n0 below rmax and zero elsewhere::

    .. code-block:: python

      dist = AnalyticDistribution(density_expression='((x**2+y**2)<rmax**2)*n0',
                                  rmax = 1.,
                                  n0 = 1.e20,
                                  ...)
    """
    density_expression: str = Field(
        description="String. Analytic expression describing physical number density [m^-3]. Expression should be in terms of the position, written as 'x', 'y', and 'z'. Parameters can be used in the expression with the values given as keyword arguments."
    )
    momentum_expressions: Sequence[str | None] = Field(
        default_factory=lambda: [None, None, None],
        description="List of strings. Analytic expressions describing the gamma*velocity for each axis [m/s]. Expressions should be in terms of the position, written as 'x', 'y', and 'z'. Parameters can be used in the expression with the values given as keyword arguments. For any axis not supplied (set to None), directed_velocity will be used."
    )
    momentum_spread_expressions: Sequence[str | None] = Field(
        default_factory=lambda: [None, None, None],
        description="List of strings. Analytic expressions describing the gamma*velocity Gaussian thermal spread sigma for each axis [m/s]. Expressions should be in terms of the position, written as 'x', 'y', and 'z'. Parameters can be used in the expression with the values given as keyword arguments. For any axis not supplied (set to None), zero will be used."
    )
    lower_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Vector of length 3 of floats, optional. Lower bound of the distribution [m]"
    )
    upper_bound: Sequence[float | None] = Field(
        default_factory=lambda: [None, None, None],
        description="Vector of length 3 of floats, optional. Upper bound of the distribution [m]"
    )
    rms_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Thermal velocity spread [m/s]"
    )
    directed_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Directed, average, proper velocity [m/s]"
    )
    fill_in: bool | None = Field(
        default=None,
        description="Bool, optional. Flags whether to fill in the empty spaced opened up when the grid moves"
    )
    user_defined_kw: dict[str, Any] = Field(
        default_factory=dict,
        exclude=True,
        description="Dictionary of user-defined keyword arguments extracted from extra fields"
    )

    @model_validator(mode='before')
    @classmethod
    def _normalize_expressions(cls, data: Any) -> Any:
        """Normalize expression strings"""
        if not isinstance(data, dict):
            return data
        
        data = dict(data)
        
        # Normalize density expression
        if 'density_expression' in data:
            data['density_expression'] = f'{data["density_expression"]}'.replace('\n', '')
        
        # Normalize momentum expressions
        if 'momentum_expressions' in data and data['momentum_expressions']:
            momentum_exprs = data['momentum_expressions']
            for i in range(3):
                if i < len(momentum_exprs) and momentum_exprs[i] is not None:
                    momentum_exprs[i] = f'{momentum_exprs[i]}'.replace('\n', '')
        
        # Normalize momentum spread expressions
        if 'momentum_spread_expressions' in data and data['momentum_spread_expressions']:
            momentum_spread_exprs = data['momentum_spread_expressions']
            for i in range(3):
                if i < len(momentum_spread_exprs) and momentum_spread_exprs[i] is not None:
                    momentum_spread_exprs[i] = f'{momentum_spread_exprs[i]}'.replace('\n', '')
        
        return data
    
    @model_validator(mode='after')
    def _convert_expressions_and_extract_kw(self) -> Self:
        """Convert momentum expressions to strings and extract user-defined keywords"""
        # Convert momentum expressions to string if needed
        for idir in range(3):
            if self.momentum_expressions[idir] is not None:
                self.momentum_expressions[idir] = f'{self.momentum_expressions[idir]}'.replace('\n', '')
            if self.momentum_spread_expressions[idir] is not None:
                self.momentum_spread_expressions[idir] = f'{self.momentum_spread_expressions[idir]}'.replace('\n', '')
        
        # Extract user-defined keywords from extra fields
        if self.model_extra:
            self.user_defined_kw = {}
            keys_to_remove = []
            
            # Check density expression
            for k, v in self.model_extra.items():
                if re.search(rf'\b{k}\b', self.density_expression):
                    self.user_defined_kw[k] = v
                    keys_to_remove.append(k)
                # Check momentum expressions
                elif self.momentum_expressions[0] is not None and re.search(rf'\b{k}\b', self.momentum_expressions[0]):
                    self.user_defined_kw[k] = v
                    keys_to_remove.append(k)
                elif self.momentum_expressions[1] is not None and re.search(rf'\b{k}\b', self.momentum_expressions[1]):
                    self.user_defined_kw[k] = v
                    keys_to_remove.append(k)
                elif self.momentum_expressions[2] is not None and re.search(rf'\b{k}\b', self.momentum_expressions[2]):
                    self.user_defined_kw[k] = v
                    keys_to_remove.append(k)
            
            # Remove extracted keys from model_extra
            for k in keys_to_remove:
                del self.model_extra[k]
        
        return self


class PICMI_ParticleListDistribution(_ClassWithInit):
    """
    Load particles at the specified positions and velocities
    """
    x: float | Sequence[float] = Field(
        default=0.,
        description="Float or list of floats, default=0. List of x positions of the particles [m]"
    )
    y: float | Sequence[float] = Field(
        default=0.,
        description="Float or list of floats, default=0. List of y positions of the particles [m]"
    )
    z: float | Sequence[float] = Field(
        default=0.,
        description="Float or list of floats, default=0. List of z positions of the particles [m]"
    )
    ux: float | Sequence[float] = Field(
        default=0.,
        description="Float or list of floats, default=0. List of ux positions of the particles (ux = gamma*vx) [m/s]"
    )
    uy: float | Sequence[float] = Field(
        default=0.,
        description="Float or list of floats, default=0. List of uy positions of the particles (uy = gamma*vy) [m/s]"
    )
    uz: float | Sequence[float] = Field(
        default=0.,
        description="Float or list of floats, default=0. List of uz positions of the particles (uz = gamma*vz) [m/s]"
    )
    weight: float | Sequence[float] = Field(
        default=0.,
        description="Float or list of floats, default=0. Particle weight or list of weights, number of real particles per simulation particle"
    )

    @model_validator(mode='after')
    def _normalize_arrays(self) -> Self:
        """Normalize arrays: convert scalars to arrays and ensure all arrays have compatible lengths"""
        # Get length of arrays, set to one for scalars
        lenx = np.size(self.x)
        leny = np.size(self.y)
        lenz = np.size(self.z)
        lenux = np.size(self.ux)
        lenuy = np.size(self.uy)
        lenuz = np.size(self.uz)
        lenw = np.size(self.weight)

        maxlen = max(lenx, leny, lenz, lenux, lenuy, lenuz, lenw)
        
        if maxlen == 1:
            # All scalars, nothing to do
            return self
        
        # Validate lengths
        if not (lenx == maxlen or lenx == 1):
            raise ValueError("Length of x doesn't match len of others")
        if not (leny == maxlen or leny == 1):
            raise ValueError("Length of y doesn't match len of others")
        if not (lenz == maxlen or lenz == 1):
            raise ValueError("Length of z doesn't match len of others")
        if not (lenux == maxlen or lenux == 1):
            raise ValueError("Length of ux doesn't match len of others")
        if not (lenuy == maxlen or lenuy == 1):
            raise ValueError("Length of uy doesn't match len of others")
        if not (lenuz == maxlen or lenuz == 1):
            raise ValueError("Length of uz doesn't match len of others")
        if not (lenw == maxlen or lenw == 1):
            raise ValueError("Length of weight doesn't match len of others")

        # Convert scalars to arrays
        if lenx == 1:
            self.x = np.array(self.x) * np.ones(maxlen)
        if leny == 1:
            self.y = np.array(self.y) * np.ones(maxlen)
        if lenz == 1:
            self.z = np.array(self.z) * np.ones(maxlen)
        if lenux == 1:
            self.ux = np.array(self.ux) * np.ones(maxlen)
        if lenuy == 1:
            self.uy = np.array(self.uy) * np.ones(maxlen)
        if lenuz == 1:
            self.uz = np.array(self.uz) * np.ones(maxlen, 'd')
        # Note that weight can be a scalar (don't convert if lenw == 1)

        return self


class PICMI_FromFileDistribution(_ClassWithInit):
    """
    Load particles from an openPMD file.

    The openPMD file must contain the attributes `position`, `momentum`, `weighting`.
    """
    file_path: str = Field(
        description="String. Path to the openPMD file"
    )

# ------------------
# Numeric Objects
# ------------------


class PICMI_ParticleDistributionPlanarInjector(_ClassWithInit):
    """
    Describes the injection of particles from a plane
    """
    position: Sequence[float] = Field(
        description="Vector of length 3 of floats. Position of the particle centroid [m]"
    )
    plane_normal: Sequence[float] = Field(
        description="Vector of length 3 of floats. Vector normal to the plane of injection [1]"
    )
    plane_velocity: Sequence[float] = Field(
        default_factory=lambda: [0., 0., 0.],
        description="Vector of length 3 of floats, default=[0.,0.,0.]. Velocity of the plane of injection [m/s]"
    )
    method: Literal['InPlace', 'Plane'] = Field(
        default='InPlace',
        description="String, default='InPlace'. Method: 'InPlace' or 'Plane'"
    )


class PICMI_GriddedLayout(_ClassWithInit):
    """
    Specifies a gridded layout of particles
    """
    n_macroparticle_per_cell: Sequence[int] = Field(
        description="Vector of integers. Number of particles per cell along each axis"
    )
    grid: "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid | PICMI_CylindricalGrid | None" = Field(
        default=None,
        description="Grid instance, optional. Grid object specifying the grid to follow. If not specified, the underlying grid of the code is used."
    )


class PICMI_PseudoRandomLayout(_ClassWithInit):
    """
    Specifies a pseudo-random layout of the particles
    """
    n_macroparticles: int | None = Field(
        default=None,
        description="Integer, optional. Total number of macroparticles to load. Either this argument or n_macroparticles_per_cell should be supplied (not both)."
    )
    n_macroparticles_per_cell: int | None = Field(
        default=None,
        description="Integer, optional. Number of macroparticles to load per cell. Either this argument or n_macroparticles should be supplied (not both)."
    )
    seed: int | None = Field(
        default=None,
        description="Integer, optional. Pseudo-random number generator seed"
    )
    grid: "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid | PICMI_CylindricalGrid | None" = Field(
        default=None,
        description="Grid instance, optional. Grid object specifying the grid to follow for n_macroparticles_per_cell. If not specified, the underlying grid of the code is used."
    )

    @model_validator(mode='after')
    def _validate_either_or(self) -> Self:
        if (self.n_macroparticles is None) == (self.n_macroparticles_per_cell is None):
            raise ValueError('Exactly one of n_macroparticles or n_macroparticles_per_cell must be specified')
        return self

# Define PICMI_Distribution type alias after all distribution classes are defined
PICMI_Distribution = (
    PICMI_GaussianBunchDistribution
    | PICMI_UniformDistribution
    | PICMI_FoilDistribution
    | PICMI_AnalyticFluxDistribution
    | PICMI_AnalyticDistribution
    | PICMI_ParticleListDistribution
    | PICMI_FromFileDistribution
)
