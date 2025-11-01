"""Simulation class following the PICMI standard
This should be the base classes for Python implementation of the PICMI standard
"""
from __future__ import annotations
from typing import Any, Literal
from collections.abc import Sequence

from pydantic import Field

from .base import _ClassWithInit

# Import classes needed for type aliases (safe at runtime since these modules don't import simulation)
from .fields import (
    PICMI_ElectromagneticSolver,
    PICMI_ElectrostaticSolver,
    PICMI_MagnetostaticSolver,
)
from .particles import (
    PICMI_Species,
    PICMI_MultiSpecies,
    PICMI_GriddedLayout,
    PICMI_PseudoRandomLayout,
    PICMI_ParticleDistributionPlanarInjector,
)
from .lasers import (
    PICMI_GaussianLaser,
    PICMI_AnalyticLaser,
    PICMI_LaserAntenna,
)
from .applied_fields import (
    PICMI_ConstantAppliedField,
    PICMI_AnalyticAppliedField,
    PICMI_Mirror,
    PICMI_LoadAppliedField,
    PICMI_LoadGriddedField,
)
from .diagnostics import (
    PICMI_FieldDiagnostic,
    PICMI_ElectrostaticFieldDiagnostic,
    PICMI_ParticleDiagnostic,
    PICMI_ParticleBoundaryScrapingDiagnostic,
    PICMI_LabFrameFieldDiagnostic,
    PICMI_LabFrameParticleDiagnostic,
)
from .interactions import PICMI_FieldIonization

# Type unions for PICMI objects using Python 3.10+ union syntax
# (Note: Not used in type hints anymore, but kept for backwards compatibility)
PICMI_Solver = (
    PICMI_ElectromagneticSolver
    | PICMI_ElectrostaticSolver
    | PICMI_MagnetostaticSolver
)
PICMI_SpeciesType = PICMI_Species | PICMI_MultiSpecies
PICMI_Layout = (
    PICMI_GriddedLayout
    | PICMI_PseudoRandomLayout
    | PICMI_ParticleDistributionPlanarInjector
)
PICMI_Laser = PICMI_GaussianLaser | PICMI_AnalyticLaser
PICMI_AppliedField = (
    PICMI_ConstantAppliedField
    | PICMI_AnalyticAppliedField
    | PICMI_Mirror
    | PICMI_LoadAppliedField
    | PICMI_LoadGriddedField
)
PICMI_Diagnostic = (
    PICMI_FieldDiagnostic
    | PICMI_ElectrostaticFieldDiagnostic
    | PICMI_ParticleDiagnostic
    | PICMI_ParticleBoundaryScrapingDiagnostic
    | PICMI_LabFrameFieldDiagnostic
    | PICMI_LabFrameParticleDiagnostic
)
PICMI_Interaction = PICMI_FieldIonization

# ---------------------
# Main simulation object
# ---------------------

class PICMI_Simulation(_ClassWithInit):
    """
    Creates a Simulation object.
    """
    solver: "PICMI_ElectromagneticSolver | PICMI_ElectrostaticSolver | PICMI_MagnetostaticSolver | None" = Field(
        default=None,
        description="Field solver instance. This is the field solver to be used in the simulation. It should be an instance of field solver classes (PICMI_ElectromagneticSolver, PICMI_ElectrostaticSolver, or PICMI_MagnetostaticSolver)."
    )
    time_step_size: float | None = Field(
        default=None,
        description="Absolute time step size of the simulation [s]. Needed if the CFL is not specified elsewhere."
    )
    verbose: int | None = Field(
        default=None,
        description="Verbosity flag. A larger integer results in more verbose output"
    )
    max_steps: int | None = Field(
        default=None,
        description="Maximum number of time steps. Specify either this, or max_time, or use the step function directly."
    )
    max_time: float | None = Field(
        default=None,
        description="Maximum physical time to run the simulation [s]. Specify either this, or max_steps, or use the step function directly."
    )
    particle_shape: Literal['NGP', 'linear', 'quadratic', 'cubic'] = Field(
        default='linear',
        description="Default particle shape for species added to this simulation"
    )
    gamma_boost: float | None = Field(
        default=None,
        description="Lorentz factor of the boosted simulation frame. Note that all input values should be in the lab frame."
    )
    load_balancing: Any = Field(
        default=None, description="Load balancing configuration"
    )  # Keep as Any - no specific PICMI class for this yet

    # Runtime lists that are populated via methods
    species: "list[PICMI_Species | PICMI_MultiSpecies]" = Field(default_factory=list, exclude=True)
    layouts: "list[PICMI_GriddedLayout | PICMI_PseudoRandomLayout | PICMI_ParticleDistributionPlanarInjector]" = Field(default_factory=list, exclude=True)
    initialize_self_fields: list[bool | None] = Field(default_factory=list, exclude=True)
    injection_plane_positions: list[Sequence[float] | None] = Field(default_factory=list, exclude=True)
    injection_plane_normal_vectors: list[Sequence[float] | None] = Field(default_factory=list, exclude=True)
    lasers: "list[PICMI_GaussianLaser | PICMI_AnalyticLaser]" = Field(default_factory=list, exclude=True)
    laser_injection_methods: list[PICMI_LaserAntenna] = Field(default_factory=list, exclude=True)
    applied_fields: "list[PICMI_ConstantAppliedField | PICMI_AnalyticAppliedField | PICMI_Mirror | PICMI_LoadAppliedField | PICMI_LoadGriddedField]" = Field(default_factory=list, exclude=True)
    diagnostics: "list[PICMI_FieldDiagnostic | PICMI_ElectrostaticFieldDiagnostic | PICMI_ParticleDiagnostic | PICMI_ParticleBoundaryScrapingDiagnostic | PICMI_LabFrameFieldDiagnostic | PICMI_LabFrameParticleDiagnostic]" = Field(default_factory=list, exclude=True)
    interactions: "list[PICMI_FieldIonization]" = Field(
        default_factory=list,
        exclude=True,
        description="List of interaction objects"
    )

    def add_species(
        self,
        species: "PICMI_Species | PICMI_MultiSpecies",
        layout: "PICMI_GriddedLayout | PICMI_PseudoRandomLayout | PICMI_ParticleDistributionPlanarInjector",
        initialize_self_field: bool | None = None,
    ):
        """
        Add species to be used in the simulation

        Parameters
        ----------
        species: species instance
            An instance of one of the PICMI species objects.
            Defines species to be added from the *physical* point of view
            (e.g. charge, mass, initial distribution of particles).

        layout: layout instance
            An instance of one of the PICMI particle layout objects.
            Defines how particles are added into the simulation, from the *numerical* point of view.

        initialize_self_field: bool, optional
            Whether the initial space-charge fields of this species
            is calculated and added to the simulation
        """
        self.species.append(species)
        self.layouts.append(layout)
        self.initialize_self_fields.append(initialize_self_field)
        self.injection_plane_positions.append(None)
        self.injection_plane_normal_vectors.append(None)


    def add_species_through_plane(
        self,
        species: "PICMI_Species | PICMI_MultiSpecies",
        layout: "PICMI_GriddedLayout | PICMI_PseudoRandomLayout | PICMI_ParticleDistributionPlanarInjector",
        injection_plane_position: Sequence[float],
        injection_plane_normal_vector: Sequence[float],
        initialize_self_field: bool | None = None,
    ):
        """
        Add species to be used in the simulation that are injected through a plane
        during the simulation.

        Parameters
        ----------
        species: species instance
            An instance of one of the PICMI species objects.
            Defines species to be added from the *physical* point of view
            (e.g. charge, mass, initial distribution of particles).

        layout: layout instance
            An instance of one of the PICMI layout objects.
            Defines how particles are added into the simulation, from the *numerical* point of view.

        initialize_self_field: bool, optional
            Whether the initial space-charge fields of this species
            is calculated and added to the simulation

        injection_plane_position: vector of floats
            Position of one point of the injection plane

        injection_plane_normal_vector: vector of floats
            Vector normal to injection plane
        """
        self.species.append(species)
        self.layouts.append(layout)
        self.initialize_self_fields.append(initialize_self_field)
        self.injection_plane_positions.append(injection_plane_position)
        self.injection_plane_normal_vectors.append(injection_plane_normal_vector)


    def add_laser(
        self,
        laser: "PICMI_GaussianLaser | PICMI_AnalyticLaser",
        injection_method: PICMI_LaserAntenna,
    ):
        """
        Add a laser pulses that to be injected in the simulation

        Parameters
        ----------
        laser_profile: laser instance
            One of laser profile instances.
            Specifies the **physical** properties of the laser pulse
            (e.g. spatial and temporal profile, wavelength, amplitude, etc.).

        injection_method: laser injection instance, optional
            Specifies how the laser is injected (numerically) into the simulation
            (e.g. through a laser antenna, or directly added to the mesh).
            This argument describes an **algorithm**, not a physical object.
            It is up to each code to define the default method
            of injection, if the user does not provide injection_method.
        """
        self.lasers.append(laser)
        self.laser_injection_methods.append(injection_method)

    def add_applied_field(
        self,
        applied_field: "PICMI_ConstantAppliedField | PICMI_AnalyticAppliedField | PICMI_Mirror | PICMI_LoadAppliedField | PICMI_LoadGriddedField",
    ):
        """
        Add an applied field

        Parameters
        ----------
        applied_field: applied field instance
            One of the applied field instance.
            Specifies the properties of the applied field.
        """
        self.applied_fields.append(applied_field)

    def add_diagnostic(
        self,
        diagnostic: "PICMI_FieldDiagnostic | PICMI_ElectrostaticFieldDiagnostic | PICMI_ParticleDiagnostic | PICMI_ParticleBoundaryScrapingDiagnostic | PICMI_LabFrameFieldDiagnostic | PICMI_LabFrameParticleDiagnostic",
    ):
        """
        Add a diagnostic

        Parameters
        ----------
        diagnostic: diagnostic instance
            One of the diagnostic instances.
        """
        self.diagnostics.append(diagnostic)

    def add_interaction(self, interaction: "PICMI_FieldIonization"):
        """
        Add an interaction

        Parameters
        ----------
        interaction: interaction instance
            One of the interaction objects.
        """
        self.interactions.append(interaction)

    def set_max_step(self, max_steps: int):
        """
        Set the default number of steps for the simulation (i.e. the number
        of steps that gets written when calling `write_input_file`).

        Note: this is equivalent to passing `max_steps` as an argument,
        when initializing the `Simulation` object

        Parameter
        ---------
        max_steps: integer
            Maximum number of time steps
        """
        self.max_steps = max_steps

    def write_input_file(self, file_name: str):
        """
        Write the parameters of the simulation, as defined in the PICMI input,
        into a code-specific input file.

        This can be used for codes that are not Python-driven (e.g. compiled,
        pure C++ or Fortran codes) and expect a text input in a given format.

        Parameters
        ----------
        file_name: string
            The path to the file that will be created
        """
        raise NotImplementedError

    def step(self, nsteps=1):
        """
        Run the simulation for `nsteps` timesteps

        Parameters
        ----------
        nsteps: integer, default=1
            The number of timesteps
        """
        raise NotImplementedError

    def extension(self):
        """
        Reserved for code-specific extensions, for example returns a class instance
        that has further methods for manipulating a PIC simulation.
        """
        raise NotImplementedError
