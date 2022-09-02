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

    Parameters
    ----------
    grid: grid instance
        Grid object for the diagnostic

    period: integer
        Period of time steps that the diagnostic is performed

    data_list: list of strings, optional
        List of quantities to write out. Possible values 'rho', 'E', 'B', 'J', 'Ex' etc.
        Defaults to the output list of the implementing code.

    write_dir: string, optional
        Directory where data is to be written

    step_min: integer, default=0
        Minimum step at which diagnostics could be written

    step_max: integer, default=unbounded
        Maximum step at which diagnostics could be written

    number_of_cells: vector of integers, optional
        Number of cells in each dimension.
        If not given, will be obtained from grid.

    lower_bound: vector of floats, optional
        Lower corner of diagnostics box in each direction.
        If not given, will be obtained from grid.

    upper_bound: vector of floats, optional
        Higher corner of diagnostics box in each direction.
        If not given, will be obtained from grid.

    parallelio: bool, optional
        If set to True, field diagnostics are dumped in parallel

    name: string, optional
        Sets the base name for the diagnostic output files
    """
    def __init__(self, grid, period, data_list=None,
                 write_dir = None,
                 step_min = None,
                 step_max = None,
                 number_of_cells = None,
                 lower_bound = None,
                 upper_bound = None,
                 parallelio = None,
                 name = None,
                 **kw):

        if data_list is not None:
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

    Parameters
    ----------
    grid: grid instance
        Grid object for the diagnostic

    period: integer
        Period of time steps that the diagnostic is performed

    data_list: list of strings, optional
        List of quantities to write out. Possible values 'rho', 'E', 'B', 'Ex' etc.
        Defaults to the output list of the implementing code.

    write_dir: string, optional
        Directory where data is to be written

    step_min: integer, default=0
        Minimum step at which diagnostics could be written

    step_max: integer, default=unbounded
        Maximum step at which diagnostics could be written

    number_of_cells: vector of integers, optional
        Number of cells in each dimension.
        If not given, will be obtained from grid.

    lower_bound: vector of floats, optional
        Lower corner of diagnostics box in each direction.
        If not given, will be obtained from grid.

    upper_bound: vector of floats, optional
        Higher corner of diagnostics box in each direction.
        If not given, will be obtained from grid.

    parallelio: bool, optional
        If set to True, field diagnostics are dumped in parallel

    name: string, optional
        Sets the base name for the diagnostic output files
    """
    def __init__(self, grid, period, data_list=None,
                 write_dir = None,
                 step_min = None,
                 step_max = None,
                 number_of_cells = None,
                 lower_bound = None,
                 upper_bound = None,
                 parallelio = None,
                 name = None,
                 **kw):

        if data_list is not None:
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

    Parameters
    ----------
    period: integer
        Period of time steps that the diagnostic is performed

    species: species instance or list of species instances
        Species to write out.
        Note that the name attribute must be defined for the species.

    data_list: list of strings, optional
        The data to be written out. Possible values 'position', 'momentum', 'weighting'.
        Defaults to the output list of the implementing code.

    write_dir: string, optional
        Directory where data is to be written

    step_min: integer, default=0
        Minimum step at which diagnostics could be written

    step_max: integer, default=unbounded
        Maximum step at which diagnostics could be written

    parallelio: bool, optional
        If set to True, particle diagnostics are dumped in parallel

    name: string, optional
        Sets the base name for the diagnostic output files
    """

    def __init__(self, period, species, data_list=None,
                 write_dir = None,
                 step_min = None,
                 step_max = None,
                 parallelio = None,
                 name = None,
                 **kw):

        if data_list is not None:
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

    Parameters
    ----------
    grid: grid instance
        Grid object for the diagnostic

    num_snapshots: integer
        Number of lab frame snapshots to make

    dt_snapshots: float
        Time between each snapshot in lab frame

    data_list: list of strings, optional
        List of quantities to write out. Possible values 'rho', 'E', 'B', 'J', 'Ex' etc.
        Defaults to the output list of the implementing code.

    z_subsampling: integer, default=1
        A factor which is applied on the resolution of the lab frame reconstruction

    time_start: float, default=0
        Time for the first snapshot in lab frame

    write_dir: string, optional
        Directory where data is to be written

    parallelio: bool, optional
        If set to True, field diagnostics are dumped in parallel

    name: string, optional
        Sets the base name for the diagnostic output files
    """
    def __init__(self, grid, num_snapshots, dt_snapshots, data_list=None,
                 z_subsampling = 1, time_start = 0.,
                 write_dir = None,
                 parallelio = None,
                 name = None,
                 **kw):

        if data_list is not None:
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

    Parameters
    ----------
    grid: grid instance
        Grid object for the diagnostic

    num_snapshots: integer
        Number of lab frame snapshots to make

    dt_snapshots: float
        Time between each snapshot in lab frame

    species: species instance or list of species instances
        Species to write out.
        Note that the name attribute must be defined for the species.

    data_list: list of strings, optional
        The data to be written out. Possible values 'position', 'momentum', 'weighting'.
        Defaults to the output list of the implementing code.

    time_start: float, default=0
        Time for the first snapshot in lab frame

    write_dir: string, optional
        Directory where data is to be written

    parallelio: bool, optional
        If set to True, particle diagnostics are dumped in parallel

    name: string, optional
        Sets the base name for the diagnostic output files
    """
    def __init__(self, grid, num_snapshots, dt_snapshots, data_list=None,
                 time_start = 0.,
                 species = None,
                 write_dir = None,
                 parallelio = None,
                 name = None,
                 **kw):

        if data_list is not None:
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
