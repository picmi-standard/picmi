"""Simulation class following the PICMI standard
This should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys

from .base import _ClassWithInit

# ---------------------
# Main simuation object
# ---------------------

class PICMI_Simulation(_ClassWithInit):
    """
    Simulation
      - solver: Field solver to be used in the simulation (an instance of one of the implemented solvers)
      - time_step_size: Absolute time step size of the simulation [s]
                        (needed if the CFL is not specified elsewhere)
      - max_steps: Maximum number of time steps
      - max_time: Maximum time to run the simulation [s]
      - verbose: Verbosity flag
      - particle_shape: Default particle shape for species added to this simulation. Possible values are 'NGP', 'linear', 'quadratic', 'cubic'
    """

    def __init__(self, solver=None, time_step_size=None, max_steps=None, max_time=None, verbose=None,
                 particle_shape='linear', **kw):

        self.solver = solver
        self.time_step_size = time_step_size
        self.verbose = verbose
        self.max_steps = max_steps
        self.max_time = max_time
        self.particle_shape = particle_shape

        self.species = []
        self.layouts = []
        self.initialize_self_fields = []

        self.lasers = []
        self.laser_injection_methods = []

        self.handle_init(kw)

    def add_species(self, species, layout, initialize_self_field=False):
        """
        Add species to be used in the simulation
        - species: species object
        - layout: particle layout for initial distribution
        - initialize_self_field=False: whether the initial space-charge fields
        of this species is calculated and added to the simulation.
        """
        self.species.append(species)
        self.layouts.append(layout)
        self.initialize_self_fields.append(initialize_self_field)

    def add_laser(self, laser, injection_method):
        """
        Add a laser pulses that to be injected in the simulation
          - laser_profile: one of laser profile objects
                           Specifies the **physical** properties of the laser pulse.
                           (e.g. spatial and temporal profile, wavelength, amplitude, etc.)
          - injection_method: a laser injector object (optional)
                              Specifies how the laser is injected (numerically) into the simulation
                              (e.g. through a laser antenna, or directly added to the mesh).
                              This argument describes an **algorithm**, not a physical object.
                              It is optional. (It is up to each code to define the default method
                              of injection, if the user does not provide injection_method)
        """
        self.lasers.append(laser)
        self.laser_injection_methods.append(injection_method)

    def write_input_file(self, file_name):
        raise NotImplementedError

    def step(self, nsteps=1):
        raise NotImplementedError
