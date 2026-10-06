"""Simulation class following the PICMI standard
This should be the base classes for Python implementation of the PICMI standard
"""
from typing import Any, Literal

from pydantic import Field

from .applied_fields import PICMI_AnyAppliedField
from .base import _PICMIModel
from .diagnostics import PICMI_AnyDiagnostic
from .fields import PICMI_AnySolver
from .interactions import PICMI_AnyInteraction
from .lasers import PICMI_AnyLaser, PICMI_AnyLaserInjection
from .particles import PICMI_AnyLayout, PICMI_AnySpecies

# ---------------------
# Main simulation object
# ---------------------

class PICMI_Simulation(_PICMIModel):
    """
    Creates a Simulation object
    """

    solver: PICMI_AnySolver | None = Field(
        default=None,
        description="This is the field solver to be used in the simulation. It should be an instance of field solver classes."
    )
    time_step_size: float | None = Field(
        default=None,
        description="Absolute time step size of the simulation [s]. Needed if the CFL is not specified elsewhere."
    )
    max_steps: int | None = Field(
        default=None,
        description="Maximum number of time steps. Specify either this, or max_time, or use the step function directly."
    )
    max_time: float | None = Field(
        default=None,
        description="Maximum physical time to run the simulation [s]. Specify either this, or max_steps, or use the step function directly."
    )
    verbose: int | None = Field(
        default=None,
        description="Verbosity flag. A larger integer results in more verbose output"
    )
    particle_shape: Literal["NGP", "linear", "quadratic", "cubic"] | int | None = Field(
        default="linear",
        description="Default particle shape for species added to this simulation. One of 'NGP', 'linear', 'quadratic', 'cubic', or the equivalent integer interpolation order."
    )
    gamma_boost: float | None = Field(
        default=None,
        description="Lorentz factor of the boosted simulation frame. Note that all input values should be in the lab frame."
    )
    # The standard leaves the meaning, and thus the type, of this parameter to the codes.
    load_balancing: Any | None = Field(
        default=None,
        description="Controls load balancing (code dependent)."
    )

    # The following lists are populated through the add_* methods rather than at construction.
    # The entries with the same index belong together, e.g., species[i] and layouts[i].
    species: list[PICMI_AnySpecies] = Field(
        default_factory=list,
        description="Species added with add_species or add_species_through_plane"
    )
    layouts: list[PICMI_AnyLayout | list[PICMI_AnyLayout] | None] = Field(
        default_factory=list,
        description="Layout (or list of layouts, one per initial distribution) of each species"
    )
    initialize_self_fields: list[bool | None] = Field(
        default_factory=list,
        description="Whether the initial space-charge fields of each species are calculated"
    )
    injection_plane_positions: list[float | list[float] | None] = Field(
        default_factory=list,
        description="Position of one point of the injection plane of each species"
    )
    injection_plane_normal_vectors: list[list[float] | None] = Field(
        default_factory=list,
        description="Vector normal to the injection plane of each species"
    )
    lasers: list[PICMI_AnyLaser] = Field(
        default_factory=list,
        description="Lasers added with add_laser"
    )
    laser_injection_methods: list[PICMI_AnyLaserInjection | None] = Field(
        default_factory=list,
        description="Injection method of each laser"
    )
    applied_fields: list[PICMI_AnyAppliedField] = Field(
        default_factory=list,
        description="Applied fields added with add_applied_field"
    )
    diagnostics: list[PICMI_AnyDiagnostic] = Field(
        default_factory=list,
        description="Diagnostics added with add_diagnostic"
    )
    interactions: list[PICMI_AnyInteraction] = Field(
        default_factory=list,
        description="Interactions added with add_interaction"
    )

    def _append(self, **entries):
        """Append to list fields. Assigning extended lists (instead of appending in place)
        validates the new entries; either all lists are extended or none."""
        with self._atomic_update():
            for name, entry in entries.items():
                setattr(self, name, [*getattr(self, name), entry])

    def add_species(self, species, layout, initialize_self_field=None):
        """
        Add species to be used in the simulation

        Parameters
        ----------
        species : PICMI_AnySpecies
            An instance of one of the PICMI species objects.
            Defines species to be added from the *physical* point of view
            (e.g. charge, mass, initial distribution of particles).

        layout : PICMI_AnyLayout, list of PICMI_AnyLayout, or None
            An instance of one of the PICMI particle layout objects (or a list of them, one per
            initial distribution). Defines how particles are added into the simulation, from the
            *numerical* point of view.

        initialize_self_field : bool, optional
            Whether the initial space-charge fields of this species
            is calculated and added to the simulation
        """
        self._append(
            species=species,
            layouts=layout,
            initialize_self_fields=initialize_self_field,
            injection_plane_positions=None,
            injection_plane_normal_vectors=None,
        )


    def add_species_through_plane(self, species, layout,
                                  injection_plane_position, injection_plane_normal_vector,
                                  initialize_self_field=None):
        """
        Add species to be used in the simulation that are injected through a plane
        during the simulation.

        Parameters
        ----------
        species : PICMI_AnySpecies
            An instance of one of the PICMI species objects.
            Defines species to be added from the *physical* point of view
            (e.g. charge, mass, initial distribution of particles).

        layout : PICMI_AnyLayout, list of PICMI_AnyLayout, or None
            An instance of one of the PICMI layout objects (or a list of them, one per initial
            distribution). Defines how particles are added into the simulation, from the
            *numerical* point of view.

        initialize_self_field : bool, optional
            Whether the initial space-charge fields of this species
            is calculated and added to the simulation

        injection_plane_position : float or list of float
            Position of one point of the injection plane

        injection_plane_normal_vector : list of float
            Vector normal to injection plane
        """
        self._append(
            species=species,
            layouts=layout,
            initialize_self_fields=initialize_self_field,
            injection_plane_positions=injection_plane_position,
            injection_plane_normal_vectors=injection_plane_normal_vector,
        )


    def add_laser(self, laser, injection_method):
        """
        Add a laser pulse that is injected into the simulation

        Parameters
        ----------
        laser : PICMI_AnyLaser
            One of the laser profile instances.
            Specifies the **physical** properties of the laser pulse
            (e.g. spatial and temporal profile, wavelength, amplitude, etc.).

        injection_method : PICMI_AnyLaserInjection or None
            Specifies how the laser is injected (numerically) into the simulation
            (e.g. through a laser antenna, or directly added to the mesh).
            This argument describes an **algorithm**, not a physical object.
            It is up to each code to define the default method
            of injection, if the user does not provide injection_method.
        """
        self._append(lasers=laser, laser_injection_methods=injection_method)

    def add_applied_field(self, applied_field):
        """
        Add an applied field

        Parameters
        ----------
        applied_field : PICMI_AnyAppliedField
            One of the applied field instances.
            Specifies the properties of the applied field.
        """
        self._append(applied_fields=applied_field)

    def add_diagnostic(self, diagnostic):
        """
        Add a diagnostic

        Parameters
        ----------
        diagnostic : PICMI_AnyDiagnostic
            One of the diagnostic instances.
        """
        self._append(diagnostics=diagnostic)

    def add_interaction(self, interaction):
        """
        Add an interaction

        Parameters
        ----------
        interaction : PICMI_AnyInteraction
            One of the interaction objects.
        """
        self._append(interactions=interaction)

    def set_max_step(self, max_steps):
        """
        Set the default number of steps for the simulation (i.e. the number
        of steps that gets written when calling `write_input_file`).

        Note: this is equivalent to passing `max_steps` as an argument,
        when initializing the `Simulation` object

        Parameters
        ----------
        max_steps : int
            Maximum number of time steps
        """
        self.max_steps = max_steps

    def write_input_file(self, file_name):
        """
        Write the parameters of the simulation, as defined in the PICMI input,
        into a code-specific input file.

        This can be used for codes that are not Python-driven (e.g. compiled,
        pure C++ or Fortran codes) and expect a text input in a given format.

        Parameters
        ----------
        file_name : str
            The path to the file that will be created
        """
        raise NotImplementedError

    def step(self, nsteps=1):
        """
        Run the simulation for `nsteps` timesteps

        Parameters
        ----------
        nsteps : int, default 1
            The number of timesteps
        """
        raise NotImplementedError

    def run(self):
        """
        Run the full simulation (up to max_time or max_step are reached)
        """
        raise NotImplementedError

    def extension(self):
        """
        Reserved for code-specific extensions, for example returns a class instance
        that has further methods for manipulating a PIC simulation.
        """
        raise NotImplementedError
