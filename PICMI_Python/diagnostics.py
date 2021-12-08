"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in the file are all diagnostics related
"""

from .base import _ClassWithInit

# ----------------------------
# Simulation frame diagnostics
# ----------------------------


class PICMI_FieldDiagnostic(_ClassWithInit):
    """
    Defines the electromagnetic field diagnostics in the simulation frame
      - grid: Grid object for the diagnostic
      - period: Period of time steps that the diagnostic is performed
      - data_list: List of quantities to write out. Possible values 'rho', 'E', 'B', 'J', 'Ex' etc.
      - write_dir='.': Directory where data is to be written
      - step_min=None: Minimum step at which diagnostics could be written (optional)
                       Defaults to step 0.
      - step_max=None: Maximum step at which diagnostics could be written (optional)
                       Defaults to no limit.
      - number_of_cells=None: Number of cells in each dimension (optional)
                              If not given, will be obtained from grid.
      - lower_bound=None: Lower corner of diagnostics box in each direction (optional)
                          If not given, will be obtained from grid.
      - upper_bound=None: Higher corner of diagnostics box in each direction (optional)
                          If not given, will be obtained from grid.
      - parallelio=None: If set to True, field diagnostics are dumped in parallel (optional)
      - name: Sets the base name for the diagnostic output files (optional)

    """
    def __init__(self, grid, period, data_list,
                 write_dir = None,
                 step_min = None,
                 step_max = None,
                 number_of_cells = None,
                 lower_bound = None,
                 upper_bound = None,
                 parallelio = None,
                 name = None,
                 **kw):

        assert isinstance(data_list, list), 'FieldDiagnostic: data_list must be a list'

        self.grid = grid
        self.period = period
        self.data_list = data_list
        self.write_dir = write_dir
        self.step_min = step_min
        self.step_max = step_max
        self.number_of_cells = number_of_cells
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.parallelio = parallelio
        self.name = name

        self.handle_init(kw)


class PICMI_ElectrostaticFieldDiagnostic(_ClassWithInit):
    """
    Defines the electrostatic field diagnostics in the simulation frame
      - grid: Grid object for the diagnostic
      - period: Period of time steps that the diagnostic is performed
      - data_list: List of quantities to write out. Possible values 'rho', 'E', 'B', 'Ex' etc.
      - write_dir='.': Directory where data is to be written
      - step_min=None: Minimum step at which diagnostics could be written (optional)
                       Defaults to step 0.
      - step_max=None: Maximum step at which diagnostics could be written (optional)
                       Defaults to no limit.
      - number_of_cells=None: Number of cells in each dimension (optional)
                              If not given, will be obtained from grid.
      - lower_bound=None: Lower corner of diagnostics box in each direction (optional)
                          If not given, will be obtained from grid.
      - upper_bound=None: Higher corner of diagnostics box in each direction (optional)
                          If not given, will be obtained from grid.
      - parallelio=None: If set to True, field diagnostics are dumped in parallel (optional)
      - name: Sets the base name for the diagnostic output files (optional)
    """
    def __init__(self, grid, period, data_list,
                 write_dir = None,
                 step_min = None,
                 step_max = None,
                 number_of_cells = None,
                 lower_bound = None,
                 upper_bound = None,
                 parallelio = None,
                 name = None,
                 **kw):

        assert isinstance(data_list, list), 'ElectrostaticFieldDiagnostic: data_list must be a list'

        self.grid = grid
        self.period = period
        self.data_list = data_list
        self.write_dir = write_dir
        self.step_min = step_min
        self.step_max = step_max
        self.number_of_cells = number_of_cells
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.parallelio = parallelio
        self.name = name

        self.handle_init(kw)


class PICMI_ParticleDiagnostic(_ClassWithInit) :
    """
    Defines the particle diagnostics in the simulation frame
      - period: Period of time steps that the diagnostic is performed
      - species: Species or list of species to write out
                 Note that the name attribute must be defined for the species.
      - data_list: The data to be written out. Possible values 'position', 'momentum', 'weighting'.
      - write_dir='.': Directory where data is to be written
      - step_min=None: Minimum step at which diagnostics could be written (optional)
                       Defaults to step 0.
      - step_max=None: Maximum step at which diagnostics could be written (optional)
                       Defaults to no limit.
      - parallelio=None: If set to True, particle diagnostics are dumped in parallel (optional)
      - name: Sets the base name for the diagnostic output files (optional)
    """

    def __init__(self, period, species, data_list,
                 write_dir = None,
                 step_min = None,
                 step_max = None,
                 parallelio = None,
                 name = None,
                 **kw):

        assert isinstance(data_list, list), 'ParticleDiagnostic: data_list must be a list'

        self.period = period
        self.species = species
        self.data_list = data_list
        self.write_dir = write_dir
        self.step_min = step_min
        self.step_max = step_max
        self.parallelio = parallelio
        self.name = name

        self.handle_init(kw)


# ----------------------------
# Lab frame diagnostics
# ----------------------------


class PICMI_LabFrameFieldDiagnostic(_ClassWithInit):
    """
    Defines the electromagnetic field diagnostics in the lab frame
      - grid: Grid object for the diagnostic
      - num_snapshots: Number of lab frame snapshots to make
      - dt_snapshots: Time between each snapshot in lab frame
      - data_list: List of quantities to write out. Possible values 'rho', 'E', 'B', 'J', 'Ex' etc.
      - z_subsampling=1: A factor which is applied on the resolution of the lab frame reconstruction. (integer)
      - time_start=0.: Time for the first snapshot in lab frame
      - write_dir='.': Directory where data is to be written
      - parallelio=None: If set to True, field diagnostics are dumped in parallel (optional)
      - name: Sets the base name for the diagnostic output files (optional)
    """
    def __init__(self, grid, num_snapshots, dt_snapshots, data_list,
                 z_subsampling = 1, time_start = 0.,
                 write_dir = None,
                 parallelio = None,
                 name = None,
                 **kw):

        assert isinstance(data_list, list), 'LabFrameFieldDiagnostic: data_list must be a list'

        self.grid = grid
        self.num_snapshots = num_snapshots
        self.dt_snapshots = dt_snapshots
        self.z_subsampling = z_subsampling
        self.time_start = time_start
        self.data_list = data_list
        self.write_dir = write_dir
        self.parallelio = parallelio
        self.name = name

        self.handle_init(kw)


class PICMI_LabFrameParticleDiagnostic(_ClassWithInit):
    """
    Defines the particle diagnostics in the lab frame
      - grid: Grid object for the diagnostic
      - num_snapshots: Number of lab frame snapshots to make
      - dt_snapshots: Time between each snapshot in lab frame
      - data_list: The data to be written out. Possible values 'position', 'momentum', 'weighting'.
      - time_start=0.: Time for the first snapshot in lab frame
      - species: Species or list of species to write out
                 Note that the name attribute must be defined for the species.
      - write_dir='.': Directory where data is to be written
      - parallelio=None: If set to True, particle diagnostics are dumped in parallel (optional)
      - name: Sets the base name for the diagnostic output files (optional)
    """
    def __init__(self, grid, num_snapshots, dt_snapshots, data_list,
                 time_start = 0.,
                 species = None,
                 write_dir = None,
                 parallelio = None,
                 name = None,
                 **kw):

        assert isinstance(data_list, list), 'LabFrameParticleDiagnostic: data_list must be a list'

        self.grid = grid
        self.num_snapshots = num_snapshots
        self.dt_snapshots = dt_snapshots
        self.time_start = time_start
        self.species = species
        self.data_list = data_list
        self.write_dir = write_dir
        self.parallelio = parallelio
        self.name = name

        self.handle_init(kw)
