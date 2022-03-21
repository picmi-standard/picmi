"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in the file are all diagnostics related
"""
import typing
from collections.abc import Sequence

from autoclass import autoargs
from typeguard import typechecked

from .base import _ClassWithInit
from . import fields
from . import particles
from . import picmi_types


# ----------------------------
# Simulation frame diagnostics
# ----------------------------


@typechecked
class PICMI_FieldDiagnostic(_ClassWithInit):
    """
    Defines the electromagnetic field diagnostics in the simulation frame
      - grid: Grid object for the diagnostic
      - period: Period of time steps that the diagnostic is performed
      - data_list=None: List of quantities to write out. Possible values 'rho', 'E', 'B', 'J', 'Ex' etc.
                        Defaults to the output list of the implementing code.
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
    @autoargs(exclude=['kw'])
    def __init__(self, grid : picmi_types.GridType,
                       period : int,
                       data_list : Sequence[str] = None,
                       write_dir : str = None,
                       step_min : int = None,
                       step_max : int = None,
                       number_of_cells : Sequence[int] = None,
                       lower_bound : Sequence[float] = None,
                       upper_bound : Sequence[float] = None,
                       parallelio : bool = None,
                       name : str = None,
                       **kw):

        if number_of_cells is not None:
            assert len(number_of_cells) == grid.number_of_dimensions, \
                'FieldDiagnostic: length of number_of_cells must be the same as the dimensionality of the grid'
        if lower_bound is not None:
            assert len(lower_bound) == grid.number_of_dimensions, \
                'FieldDiagnostic: length of lower_bound must be the same as the dimensionality of the grid'
        if upper_bound is not None:
            assert len(upper_bound) == grid.number_of_dimensions, \
                'FieldDiagnostic: length of upper_bound must be the same as the dimensionality of the grid'

        self.handle_init(kw)


@typechecked
class PICMI_ElectrostaticFieldDiagnostic(_ClassWithInit):
    """
    Defines the electrostatic field diagnostics in the simulation frame
      - grid: Grid object for the diagnostic
      - period: Period of time steps that the diagnostic is performed
      - data_list=None: List of quantities to write out. Possible values 'rho', 'E', 'B', 'Ex' etc.
                        Defaults to the output list of the implementing code.
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
    @autoargs(exclude=['kw'])
    def __init__(self, grid : picmi_types.GridType,
                       period : int,
                       data_list : Sequence[str] = None,
                       write_dir : str = None,
                       step_min : int = None,
                       step_max : int = None,
                       number_of_cells : Sequence[int] = None,
                       lower_bound : Sequence[float] = None,
                       upper_bound : Sequence[float] = None,
                       parallelio : bool = None,
                       name : str = None,
                       **kw):

        if number_of_cells is not None:
            assert len(number_of_cells) == grid.number_of_dimensions, \
                'ElectrostaticFieldDiagnostic: length of number_of_cells must be the same as the dimensionality of the grid'
        if lower_bound is not None:
            assert len(lower_bound) == grid.number_of_dimensions, \
                'ElectrostaticFieldDiagnostic: length of lower_bound must be the same as the dimensionality of the grid'
        if upper_bound is not None:
            assert len(upper_bound) == grid.number_of_dimensions, \
                'ElectrostaticFieldDiagnostic: length of upper_bound must be the same as the dimensionality of the grid'

        self.handle_init(kw)


@typechecked
class PICMI_ParticleDiagnostic(_ClassWithInit) :
    """
    Defines the particle diagnostics in the simulation frame
      - period: Period of time steps that the diagnostic is performed
      - species: Species or list of species to write out
                 Note that the name attribute must be defined for the species.
      - data_list=None: The data to be written out. Possible values 'position', 'momentum', 'weighting'.
                        Defaults to the output list of the implementing code.
      - write_dir='.': Directory where data is to be written
      - step_min=None: Minimum step at which diagnostics could be written (optional)
                       Defaults to step 0.
      - step_max=None: Maximum step at which diagnostics could be written (optional)
                       Defaults to no limit.
      - parallelio=None: If set to True, particle diagnostics are dumped in parallel (optional)
      - name: Sets the base name for the diagnostic output files (optional)
    """

    @autoargs(exclude=['kw'])
    def __init__(self, period : int,
                       species : picmi_types.SpeciesArgument,
                       data_list : Sequence[str] = None,
                       write_dir : str = None,
                       step_min : int = None,
                       step_max : int = None,
                       parallelio : bool = None,
                       name : str = None,
                       **kw):

        self.handle_init(kw)


# ----------------------------
# Lab frame diagnostics
# ----------------------------


@typechecked
class PICMI_LabFrameFieldDiagnostic(_ClassWithInit):
    """
    Defines the electromagnetic field diagnostics in the lab frame
      - grid: Grid object for the diagnostic
      - num_snapshots: Number of lab frame snapshots to make
      - dt_snapshots: Time between each snapshot in lab frame
      - data_list=None: List of quantities to write out. Possible values 'rho', 'E', 'B', 'J', 'Ex' etc.
                        Defaults to the output list of the implementing code.
      - z_subsampling=1: A factor which is applied on the resolution of the lab frame reconstruction. (integer)
      - time_start=0.: Time for the first snapshot in lab frame
      - write_dir='.': Directory where data is to be written
      - parallelio=None: If set to True, field diagnostics are dumped in parallel (optional)
      - name: Sets the base name for the diagnostic output files (optional)
    """
    @autoargs(exclude=['kw'])
    def __init__(self, grid : picmi_types.GridType,
                       num_snapshots : int,
                       dt_snapshots : float,
                       data_list : Sequence[str] = None,
                       z_subsampling : int = 1,
                       time_start : float = 0.,
                       write_dir : str = None,
                       parallelio : bool = None,
                       name : str = None,
                       **kw):

        self.handle_init(kw)


@typechecked
class PICMI_LabFrameParticleDiagnostic(_ClassWithInit):
    """
    Defines the particle diagnostics in the lab frame
      - grid: Grid object for the diagnostic
      - num_snapshots: Number of lab frame snapshots to make
      - dt_snapshots: Time between each snapshot in lab frame
      - data_list=None: The data to be written out. Possible values 'position', 'momentum', 'weighting'.
                        Defaults to the output list of the implementing code.
      - time_start=0.: Time for the first snapshot in lab frame
      - species: Species or list of species to write out
                 Note that the name attribute must be defined for the species.
      - write_dir='.': Directory where data is to be written
      - parallelio=None: If set to True, particle diagnostics are dumped in parallel (optional)
      - name: Sets the base name for the diagnostic output files (optional)
    """
    @autoargs(exclude=['kw'])
    def __init__(self, grid : picmi_types.GridType,
                       num_snapshots : int,
                       dt_snapshots : float,
                       data_list : Sequence[str] = None,
                       time_start : float = 0.,
                       species : picmi_types.SpeciesArgument = None,
                       write_dir : str = None,
                       parallelio : bool = None,
                       name : str = None,
                       **kw):

        self.handle_init(kw)
