"""Simulation class following the PICMI standard
This should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys
import typing

from autoclass import autoargs
from typeguard import typechecked

from .base import _ClassWithInit
from . import fields
from . import particles
from . import diagnostics
from . import lasers
from . import diagnostics
from . import picmi_types
from . import applied_fields

# ---------------------
# Main simulation object
# ---------------------

@typechecked
class PICMI_Simulation(_ClassWithInit):
    """
    Creates a Simulation object

    Parameters
    ----------
    solver: object
        An instance of one of the PICMI field solvers ; see :doc:`field`
        This is the field solver to be used in the simulation

    time_step_size: float
        Absolute time step size of the simulation [s]
        (needed if the CFL is not specified elsewhere)

    max_steps: int
        Maximum number of time steps
        (Specify either this, or `max_time`, or use the `step` function directly)

    max_time: float
        Maximum physical time to run the simulation [s]
        (Specify either this, or `max_steps`, or use the `step` function directly)

    verbose: int
        Verbosity flag (A larger integer results in more verbose output.)

    particle_shape: str
        Default particle shape for species added to this simulation.
        Possible values are 'NGP', 'linear', 'quadratic', 'cubic'.

    gamma_boost: float
        Lorentz factor of the boosted simulation frame.
        (Note that all input values should be in the lab frame)

    kw: additional code arguments
        Code specific arguments ; should be prefixed with the `codename`
    """

    @autoargs(exclude=['kw'])
    def __init__(self, solver : picmi_types.SolverType = None,
                       time_step_size : float = None,
                       max_steps : int = None,
                       max_time : float = None,
                       verbose : int = None,
                       particle_shape : str = 'linear',
                       gamma_boost : float = None,
                       cpu_split : bool = None,
                       load_balancing : bool = None,
                       **kw):

        self.species = []
        self.layouts = []
        self.initialize_self_fields = []
        self.injection_plane_positions = []
        self.injection_plane_normal_vectors = []

        self.lasers = []
        self.laser_injection_methods = []

        self.applied_fields = []

        self.diagnostics = []

        self.handle_init(kw)

    def add_species(self, species : picmi_types.SpeciesType,
                          layout : picmi_types.LayoutType = None,
                          initialize_self_field : bool = None):
        """
        Add species to be used in the simulation

        Parameters
        ----------
        - species : object
              an instance of one of the PICMI species object ; see :doc:`species`
              Defines added species from the *physical* point of view
              (e.g. charge, mass, initial distribution of particles)

        - layout : object
              an instance of one of the PICMI layout object ; see :doc:`layout`
              Defines how particles are added into the simulation, from the *numerical* point of view

        - initialize_self_field : bool
              Whether the initial space-charge fields of this species
              is calculated and added to the simulation
        """
        self.species.append(species)
        self.layouts.append(layout)
        self.initialize_self_fields.append(initialize_self_field)
        self.injection_plane_positions.append(None)
        self.injection_plane_normal_vectors.append(None)


    def add_species_through_plane(self, species : picmi_types.SpeciesType,
                                        layout : picmi_types.LayoutType,
                                        injection_plane_position : picmi_types.VectorFloat3,
                                        injection_plane_normal_vector : picmi_types.VectorFloat3,
                                        initialize_self_field : bool = None ):
        """
        Add species to be used in the simulation

        Parameters
        ----------
        Same as `add_species`, with the following two possible arguments:

        - injection_plane_position: Position of one point of the injection plane (vector of floats)
        - injection_plane_normal_vector: Vector normal to injection plane (vector of floats)
        """
        self.species.append(species)
        self.layouts.append(layout)
        self.initialize_self_fields.append(initialize_self_field)
        self.injection_plane_positions.append(injection_plane_position)
        self.injection_plane_normal_vectors.append(injection_plane_normal_vector)


    def add_laser(self, laser : picmi_types.LaserType,
                        injection_method : picmi_types.LaserInjectionType):
        """
        Add a laser pulses that to be injected in the simulation

        Parameters
        ----------
          - laser_profile: object
                one of laser profile objects
                Specifies the **physical** properties of the laser pulse.
                (e.g. spatial and temporal profile, wavelength, amplitude, etc.)

          - injection_method: object (optional)
                a laser injector object (optional)
                Specifies how the laser is injected (numerically) into the simulation
                (e.g. through a laser antenna, or directly added to the mesh).
                This argument describes an **algorithm**, not a physical object.
                It is optional. (It is up to each code to define the default method
                of injection, if the user does not provide injection_method)
        """
        self.lasers.append(laser)
        self.laser_injection_methods.append(injection_method)

    def add_applied_field(self, applied_field : picmi_types.AppliedFieldType):
        """
        Add an applied field

        Parameters
        ----------
          - applied_field: object
                one of the applied field objects
                Specifies the properties of the applied field.
        """
        self.applied_fields.append(applied_field)

    def add_diagnostic(self, diagnostic : picmi_types.DiagnosticType):
        """
        Add a diagnostic
          - diagnostic: object
                One of the diagnostic objects.
        """
        self.diagnostics.append(diagnostic)

    def set_max_step(self, max_steps : int):
        """
        Set the default number of steps for the simulation (i.e. the number
        of steps that gets written when calling `write_input_file`)

        Note: this is equivalent to passing `max_steps` as an argument,
        when initializing the `Simulation` object

        Parameter
        ---------
        max_steps: int
            Maximum number of time steps
        """
        self.max_steps = max_steps

    def write_input_file(self, file_name : str):
        """
        Write the parameters of the simulation, as defined in the PICMI input,
        into another, more code-specific input file.

        This can be used for codes that are not Python-driven (e.g. compiled,
        pure C++ or Fortran codes) and expect a text input in a given format.

        Parameters
        ----------
        file_name: str
            The path to the file that will be created.
        """
        raise NotImplementedError

    def step(self, nsteps : int = 1):
        """
        Run the simulation for `nsteps` timesteps

        Parameters
        ----------
        nsteps: int
            The number of timesteps
        """
        raise NotImplementedError

    def extension(self):
        """
        Reserved for code-specific extensions, for example a class instance
        that has further methods for manipulating a PIC simulation.
        """
        raise NotImplementedError
