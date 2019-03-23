"""Simulation class following the PICMI standard
This should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys

from .base import _ClassWithInit

# ---------------------
# Main simulation object
# ---------------------

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

    def __init__(self, solver=None, time_step_size=None, max_steps=None, max_time=None, verbose=None,
                particle_shape='linear', gamma_boost=None, cpu_split=None, load_balancing=None, **kw):

        self.solver = solver
        self.time_step_size = time_step_size
        self.verbose = verbose
        self.max_steps = max_steps
        self.max_time = max_time
        self.particle_shape = particle_shape
        self.gamma_boost = gamma_boost

        self.species = []
        self.layouts = []
        self.initialize_self_fields = []

        self.lasers = []
        self.laser_injection_methods = []

        self.diagnostics = []
        
        self.cpu_split = cpu_split
        self.load_balancing = load_balancing
        

        self.handle_init(kw)

    def add_species(self, species, layout, initialize_self_field=False):
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

    def add_laser(self, laser, injection_method):
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

    def add_diagnostic(self, diagnostic):
        """
        Add a diagnostic
          - diagnostic: object 
                One of the diagnostic objects. 
        """
        self.diagnostics.append(diagnostic)

    def write_input_file(self, file_name):
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

    def step(self, nsteps=1):
        """
        Run the simulation for `nsteps` timesteps

        Parameters
        ----------
        nsteps: int
            The number of timesteps
        """
        raise NotImplementedError
