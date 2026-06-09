"""Simulation class following the PICMI standard
This should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys
from typing import Any

from pydantic import Field

from .base import _ClassWithInit, _PICMIModel

# ---------------------
# Main simulation object
# ---------------------

class PICMI_Simulation(_PICMIModel):
    """
    Creates a Simulation object
    """

    solver: Any | None = Field(
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
    particle_shape: str | int | None = Field(
        default="linear",
        description="Default particle shape for species added to this simulation. One of 'NGP', 'linear', 'quadratic', 'cubic', or the equivalent integer interpolation order."
    )
    gamma_boost: float | None = Field(
        default=None,
        description="Lorentz factor of the boosted simulation frame. Note that all input values should be in the lab frame."
    )
    load_balancing: Any | None = Field(
        default=None,
        description="Controls load balancing (code dependent)."
    )

    # The following lists are populated through the add_* methods rather than at construction.
    species: list = Field(default_factory=list)
    layouts: list = Field(default_factory=list)
    initialize_self_fields: list = Field(default_factory=list)
    injection_plane_positions: list = Field(default_factory=list)
    injection_plane_normal_vectors: list = Field(default_factory=list)
    lasers: list = Field(default_factory=list)
    laser_injection_methods: list = Field(default_factory=list)
    applied_fields: list = Field(default_factory=list)
    diagnostics: list = Field(default_factory=list)
    interactions: list = Field(default_factory=list)

    def add_species(self, species, layout, initialize_self_field=None):
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


    def add_species_through_plane(self, species, layout,
                                  injection_plane_position, injection_plane_normal_vector,
                                  initialize_self_field=None):
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


    def add_laser(self, laser, injection_method):
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

    def add_applied_field(self, applied_field):
        """
        Add an applied field

        Parameters
        ----------
        applied_field: applied field instance
            One of the applied field instance.
            Specifies the properties of the applied field.
        """
        self.applied_fields.append(applied_field)

    def add_diagnostic(self, diagnostic):
        """
        Add a diagnostic

        Parameters
        ----------
        diagnostic: diagnostic instance
            One of the diagnostic instances.
        """
        self.diagnostics.append(diagnostic)

    def add_interaction(self, interaction):
        """
        Add an interaction

        Parameters
        ----------
        interaction: interaction instance
            One of the interaction objects.
        """
        self.interactions.append(interaction)

    def set_max_step(self, max_steps):
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

    def write_input_file(self, file_name):
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
