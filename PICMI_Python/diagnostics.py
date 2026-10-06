"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in the file are all diagnostics related
"""

from pydantic import Field

from .base import PICMI_Diagnostic, _PICMIModel
from .fields import PICMI_AnyGrid
from .particles import PICMI_AnySpecies

# ----------------------------
# Simulation frame diagnostics
# ----------------------------


class PICMI_FieldDiagnostic(PICMI_Diagnostic):
    """
    Defines the electromagnetic field diagnostics in the simulation frame
    """
    grid: PICMI_AnyGrid = Field(
        description="Grid object for the diagnostic"
    )
    period: int = Field(
        description="Period of time steps that the diagnostic is performed"
    )
    data_list: list[str] | None = Field(
        default=None,
        description="List of quantities to write out. Possible values 'rho', 'E', 'B', 'J', 'Ex' etc. Defaults to the output list of the implementing code."
    )
    write_dir: str | None = Field(
        default=None,
        description="Directory where data is to be written"
    )
    step_min: int | None = Field(
        default=None,
        description="Minimum step at which diagnostics could be written (default 0)"
    )
    step_max: int | None = Field(
        default=None,
        description="Maximum step at which diagnostics could be written (default unbounded)"
    )
    number_of_cells: list[int] | None = Field(
        default=None,
        description="Number of cells in each dimension. If not given, will be obtained from grid."
    )
    lower_bound: list[float] | None = Field(
        default=None,
        description="Lower corner of diagnostics box in each direction. If not given, will be obtained from grid."
    )
    upper_bound: list[float] | None = Field(
        default=None,
        description="Higher corner of diagnostics box in each direction. If not given, will be obtained from grid."
    )
    parallelio: bool | None = Field(
        default=None,
        description="If set to True, field diagnostics are dumped in parallel"
    )
    name: str | None = Field(
        default=None,
        description="Sets the base name for the diagnostic output files"
    )


class PICMI_ElectrostaticFieldDiagnostic(PICMI_Diagnostic):
    """
    Defines the electrostatic field diagnostics in the simulation frame
    """
    grid: PICMI_AnyGrid = Field(
        description="Grid object for the diagnostic"
    )
    period: int = Field(
        description="Period of time steps that the diagnostic is performed"
    )
    data_list: list[str] | None = Field(
        default=None,
        description="List of quantities to write out. Possible values 'rho', 'E', 'B', 'Ex' etc. Defaults to the output list of the implementing code."
    )
    write_dir: str | None = Field(
        default=None,
        description="Directory where data is to be written"
    )
    step_min: int | None = Field(
        default=None,
        description="Minimum step at which diagnostics could be written (default 0)"
    )
    step_max: int | None = Field(
        default=None,
        description="Maximum step at which diagnostics could be written (default unbounded)"
    )
    number_of_cells: list[int] | None = Field(
        default=None,
        description="Number of cells in each dimension. If not given, will be obtained from grid."
    )
    lower_bound: list[float] | None = Field(
        default=None,
        description="Lower corner of diagnostics box in each direction. If not given, will be obtained from grid."
    )
    upper_bound: list[float] | None = Field(
        default=None,
        description="Higher corner of diagnostics box in each direction. If not given, will be obtained from grid."
    )
    parallelio: bool | None = Field(
        default=None,
        description="If set to True, field diagnostics are dumped in parallel"
    )
    name: str | None = Field(
        default=None,
        description="Sets the base name for the diagnostic output files"
    )


class PICMI_ParticleDiagnostic(PICMI_Diagnostic):
    """
    Defines the particle diagnostics in the simulation frame
    """
    period: int = Field(
        description="Period of time steps that the diagnostic is performed"
    )
    species: PICMI_AnySpecies | list[PICMI_AnySpecies] | None = Field(
        default=None,
        description="Species instance (or MultiSpecies instance) or list of them to write out. If not specified, all species are written. Note that the name attribute must be defined for the species."
    )
    data_list: list[str] | None = Field(
        default=None,
        description="The data to be written out. Possible values 'position', 'momentum', 'weighting'. Defaults to the output list of the implementing code."
    )
    write_dir: str | None = Field(
        default=None,
        description="Directory where data is to be written"
    )
    step_min: int | None = Field(
        default=None,
        description="Minimum step at which diagnostics could be written (default 0)"
    )
    step_max: int | None = Field(
        default=None,
        description="Maximum step at which diagnostics could be written (default unbounded)"
    )
    parallelio: bool | None = Field(
        default=None,
        description="If set to True, particle diagnostics are dumped in parallel"
    )
    name: str | None = Field(
        default=None,
        description="Sets the base name for the diagnostic output files"
    )


class PICMI_ParticleBoundaryScrapingDiagnostic(PICMI_Diagnostic):
    """
    Defines the particle diagnostics that are used to collect the particles that are absorbed at the boundaries, throughout the simulation.
    """
    period: int = Field(
        description="Period of time steps that the diagnostic is performed"
    )
    species: PICMI_AnySpecies | list[PICMI_AnySpecies] | None = Field(
        default=None,
        description="Species instance (or MultiSpecies instance) or list of them to write out. If not specified, all species are written. Note that the name attribute must be defined for the species."
    )
    data_list: list[str] | None = Field(
        default=None,
        description="The data to be written out. Possible values 'position', 'momentum', 'weighting'. Defaults to the output list of the implementing code."
    )
    write_dir: str | None = Field(
        default=None,
        description="Directory where data is to be written"
    )
    parallelio: bool | None = Field(
        default=None,
        description="If set to True, particle diagnostics are dumped in parallel"
    )
    name: str | None = Field(
        default=None,
        description="Sets the base name for the diagnostic output files"
    )


# ----------------------------
# Lab frame diagnostics
# ----------------------------


class PICMI_LabFrameFieldDiagnostic(PICMI_Diagnostic):
    """
    Defines the electromagnetic field diagnostics in the lab frame
    """
    grid: PICMI_AnyGrid = Field(
        description="Grid object for the diagnostic"
    )
    num_snapshots: int = Field(
        description="Number of lab frame snapshots to make"
    )
    dt_snapshots: float = Field(
        description="Time between each snapshot in lab frame"
    )
    data_list: list[str] | None = Field(
        default=None,
        description="List of quantities to write out. Possible values 'rho', 'E', 'B', 'J', 'Ex' etc. Defaults to the output list of the implementing code."
    )
    z_subsampling: int = Field(
        default=1,
        description="A factor which is applied on the resolution of the lab frame reconstruction"
    )
    time_start: float = Field(
        default=0.,
        description="Time for the first snapshot in lab frame"
    )
    write_dir: str | None = Field(
        default=None,
        description="Directory where data is to be written"
    )
    parallelio: bool | None = Field(
        default=None,
        description="If set to True, field diagnostics are dumped in parallel"
    )
    name: str | None = Field(
        default=None,
        description="Sets the base name for the diagnostic output files"
    )


class PICMI_LabFrameParticleDiagnostic(PICMI_Diagnostic):
    """
    Defines the particle diagnostics in the lab frame
    """
    grid: PICMI_AnyGrid = Field(
        description="Grid object for the diagnostic"
    )
    num_snapshots: int = Field(
        description="Number of lab frame snapshots to make"
    )
    dt_snapshots: float = Field(
        description="Time between each snapshot in lab frame"
    )
    data_list: list[str] | None = Field(
        default=None,
        description="The data to be written out. Possible values 'position', 'momentum', 'weighting'. Defaults to the output list of the implementing code."
    )
    time_start: float = Field(
        default=0.,
        description="Time for the first snapshot in lab frame"
    )
    species: PICMI_AnySpecies | list[PICMI_AnySpecies] | None = Field(
        default=None,
        description="Species instance (or MultiSpecies instance) or list of them to write out. If not specified, all species are written. Note that the name attribute must be defined for the species."
    )
    write_dir: str | None = Field(
        default=None,
        description="Directory where data is to be written"
    )
    parallelio: bool | None = Field(
        default=None,
        description="If set to True, particle diagnostics are dumped in parallel"
    )
    name: str | None = Field(
        default=None,
        description="Sets the base name for the diagnostic output files"
    )


PICMI_AnyDiagnostic = (
    PICMI_FieldDiagnostic
    | PICMI_ElectrostaticFieldDiagnostic
    | PICMI_ParticleDiagnostic
    | PICMI_ParticleBoundaryScrapingDiagnostic
    | PICMI_LabFrameFieldDiagnostic
    | PICMI_LabFrameParticleDiagnostic
    | PICMI_Diagnostic
)
