"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in the file are all diagnostics related
"""
import numpy as np

from .base import _ClassWithInit

# ----------------------------
# Simulation frame diagnostics
# ----------------------------


class PICMI_FieldDiagnostic(_ClassWithInit):
    """
    Defines the electromagnetic field diagnostics in the simulation frame
      - period=1: Period of time steps that the diagnostic is performed
      - field_types=["rho", "E", "B", "J"]: List of field types to write out
      - write_dir='.': Directory where data is to be written
      - step_min=0: Minimum step at which diagnostics could be written
      - step_max=inf: Maximum step at which diagnostics could be written
    """
    def __init__(self, period = 1,
                 field_types = ["rho", "E", "B", "J"],
                 write_dir = None,
                 step_min = 0,
                 step_max = np.inf,
                 **kw):

        self.period = period
        self.field_types = field_types
        self.write_dir = write_dir
        self.step_min = step_min
        self.step_max = step_max

        self.handle_init(kw)


class PICMI_ElectrostaticFieldDiagnostic(_ClassWithInit):
    """
    Defines the electrostatic field diagnostics in the simulation frame
      - period=1: Period of time steps that the diagnostic is performed
      - write_dir='.': Directory where data is to be written
      - step_min=0: Minimum step at which diagnostics could be written
      - step_max=inf: Maximum step at which diagnostics could be written
    """
    def __init__(self, period = 1,
                 write_dir = None,
                 step_min = 0,
                 step_max = np.inf,
                 **kw):

        self.period = period
        self.write_dir = write_dir
        self.step_min = step_min
        self.step_max = step_max

        self.handle_init(kw)


class PICMI_ParticleDiagnostic(_ClassWithInit) :
    """
    Defines the particle diagnostics in the simulation frame
      - period=1: Period of time steps that the diagnostic is performed
      - species: Species or list of species to write out
                 Note that the name attribute must be defined for the species.
      - particle_data=["position", "momentum", "weighting"]: The data to be written out
      - write_dir='.': Directory where data is to be written
      - step_min=0: Minimum step at which diagnostics could be written
      - step_max=inf: Maximum step at which diagnostics could be written
    """

    def __init__(self, period = 1,
                 species = None,
                 particle_data = ["position", "momentum", "weighting"],
                 write_dir = None,
                 step_min = 0,
                 step_max = np.inf,
                 **kw):

        self.period = period
        self.species = species
        self.particle_data = particle_data
        self.write_dir = write_dir
        self.step_min = step_min
        self.step_max = step_max

        self.handle_init(kw)


# ----------------------------
# Lab frame diagnostics
# ----------------------------


class PICMI_LabFrameFieldDiagnostic(_ClassWithInit):
    """
    Defines the electromagnetic field diagnostics in the lab frame
      - num_snapshots: Number of lab frame snapshots to make
      - dt_snapshots: Time between each snapshot
      - field_types=["rho", "E", "B", "J"]: List of field types to write out
      - write_dir='.': Directory where data is to be written
    """
    def __init__(self, num_snapshots, dt_snapshots,
                 field_types = ["rho", "E", "B", "J"],
                 write_dir = None,
                 **kw):

        self.num_snapshots = num_snapshots
        self.dt_snapshots = dt_snapshots
        self.field_types = field_types
        self.write_dir = write_dir

        self.handle_init(kw)


class PICMI_LabFrameParticleDiagnostic(_ClassWithInit):
    """
    Defines the particle diagnostics in the lab frame
      - num_snapshots: Number of lab frame snapshots to make
      - dt_snapshots: Time between each snapshot
      - species: Species or list of species to write out
                 Note that the name attribute must be defined for the species.
      - particle_data=["position", "momentum", "weighting"]: The data to be written out
      - write_dir='.': Directory where data is to be written
    """
    def __init__(self, num_snapshots, dt_snapshots,
                 species = None,
                 particle_data = ["position", "momentum", "weighting"],
                 write_dir = None,
                 **kw):

        self.num_snapshots = num_snapshots
        self.dt_snapshots = dt_snapshots
        self.species = species
        self.particle_data = particle_data
        self.write_dir = write_dir

        self.handle_init(kw)
