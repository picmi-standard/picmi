"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in the file are all diagnostics related
"""
from __future__ import annotations
from typing import Any
from collections.abc import Sequence

from pydantic import Field, field_validator

from .base import _ClassWithInit

from .fields import (
    PICMI_Cartesian1DGrid,
    PICMI_Cartesian2DGrid,
    PICMI_Cartesian3DGrid,
    PICMI_CylindricalGrid,
)
from .particles import PICMI_Species, PICMI_MultiSpecies

# Type aliases using Python 3.10+ union syntax
PICMI_Grid = (
    PICMI_Cartesian1DGrid
    | PICMI_Cartesian2DGrid
    | PICMI_Cartesian3DGrid
    | PICMI_CylindricalGrid
)
PICMI_SpeciesType = PICMI_Species | PICMI_MultiSpecies

# ----------------------------
# Simulation frame diagnostics
# ----------------------------


class PICMI_FieldDiagnostic(_ClassWithInit):
    """
    Defines the electromagnetic field diagnostics in the simulation frame.
    """
    grid: "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid | PICMI_CylindricalGrid" = Field(description="Grid object for the diagnostic")
    period: int = Field(description="Period of time steps that the diagnostic is performed")
    data_list: list[str] | None = Field(
        default=None,
        description="List of quantities to write out. Possible values 'rho', 'E', 'B', 'J', 'Ex' etc. Defaults to the output list of the implementing code."
    )
    write_dir: str | None = Field(default=None, description="Directory where data is to be written")
    step_min: int | None = Field(default=None, description="Minimum step at which diagnostics could be written")
    step_max: int | None = Field(default=None, description="Maximum step at which diagnostics could be written")
    number_of_cells: Sequence[int] | None = Field(
        default=None,
        description="Number of cells in each dimension. If not given, will be obtained from grid."
    )
    lower_bound: Sequence[float] | None = Field(
        default=None,
        description="Lower corner of diagnostics box in each direction. If not given, will be obtained from grid."
    )
    upper_bound: Sequence[float] | None = Field(
        default=None,
        description="Higher corner of diagnostics box in each direction. If not given, will be obtained from grid."
    )
    parallelio: bool | None = Field(default=None, description="If set to True, field diagnostics are dumped in parallel")
    name: str | None = Field(default=None, description="Sets the base name for the diagnostic output files")
    
    @field_validator('data_list')
    @classmethod
    def _validate_data_list(cls, v):
        if v is not None and not isinstance(v, list):
            raise ValueError('FieldDiagnostic: data_list must be a list')
        return v


class PICMI_ElectrostaticFieldDiagnostic(_ClassWithInit):
    """
    Defines the electrostatic field diagnostics in the simulation frame.
    """
    grid: "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid | PICMI_CylindricalGrid" = Field(description="Grid object for the diagnostic")
    period: int = Field(description="Period of time steps that the diagnostic is performed")
    data_list: list[str] | None = Field(
        default=None,
        description="List of quantities to write out. Possible values 'rho', 'E', 'B', 'Ex' etc. Defaults to the output list of the implementing code."
    )
    write_dir: str | None = Field(default=None, description="Directory where data is to be written")
    step_min: int | None = Field(default=None, description="Minimum step at which diagnostics could be written")
    step_max: int | None = Field(default=None, description="Maximum step at which diagnostics could be written")
    number_of_cells: Sequence[int] | None = Field(
        default=None,
        description="Number of cells in each dimension. If not given, will be obtained from grid."
    )
    lower_bound: Sequence[float] | None = Field(
        default=None,
        description="Lower corner of diagnostics box in each direction. If not given, will be obtained from grid."
    )
    upper_bound: Sequence[float] | None = Field(
        default=None,
        description="Higher corner of diagnostics box in each direction. If not given, will be obtained from grid."
    )
    parallelio: bool | None = Field(default=None, description="If set to True, field diagnostics are dumped in parallel")
    name: str | None = Field(default=None, description="Sets the base name for the diagnostic output files")
    
    @field_validator('data_list')
    @classmethod
    def _validate_data_list(cls, v):
        if v is not None and not isinstance(v, list):
            raise ValueError('ElectrostaticFieldDiagnostic: data_list must be a list')
        return v


class PICMI_ParticleDiagnostic(_ClassWithInit) :
    """
    Defines the particle diagnostics in the simulation frame.
    """
    period: int = Field(description="Period of time steps that the diagnostic is performed")
    species: "PICMI_Species | PICMI_MultiSpecies | list[PICMI_Species | PICMI_MultiSpecies] | None" = Field(
        default=None,
        description="Species instance or list of species instances. Species to write out. If not specified, all species are written. Note that the name attribute must be defined for the species."
    )
    data_list: list[str] | None = Field(
        default=None,
        description="The data to be written out. Possible values 'position', 'momentum', 'weighting'. Defaults to the output list of the implementing code."
    )
    write_dir: str | None = Field(default=None, description="Directory where data is to be written")
    step_min: int | None = Field(default=None, description="Minimum step at which diagnostics could be written")
    step_max: int | None = Field(default=None, description="Maximum step at which diagnostics could be written")
    parallelio: bool | None = Field(default=None, description="If set to True, particle diagnostics are dumped in parallel")
    name: str | None = Field(default=None, description="Sets the base name for the diagnostic output files")
    
    @field_validator('data_list')
    @classmethod
    def _validate_data_list(cls, v):
        if v is not None and not isinstance(v, list):
            raise ValueError('ParticleDiagnostic: data_list must be a list')
        return v


class PICMI_ParticleBoundaryScrapingDiagnostic(_ClassWithInit) :
    """
    Defines the particle diagnostics that are used to collect the particles that are absorbed at the boundaries, throughout the simulation.
    """
    period: int = Field(description="Period of time steps that the diagnostic is performed")
    species: "PICMI_Species | PICMI_MultiSpecies | list[PICMI_Species | PICMI_MultiSpecies] | None" = Field(
        default=None,
        description="Species instance or list of species instances. Species to write out. If not specified, all species are written. Note that the name attribute must be defined for the species."
    )
    data_list: list[str] | None = Field(
        default=None,
        description="The data to be written out. Possible values 'position', 'momentum', 'weighting'. Defaults to the output list of the implementing code."
    )
    write_dir: str | None = Field(default=None, description="Directory where data is to be written")
    parallelio: bool | None = Field(default=None, description="If set to True, particle diagnostics are dumped in parallel")
    name: str | None = Field(default=None, description="Sets the base name for the diagnostic output files")
    
    @field_validator('data_list')
    @classmethod
    def _validate_data_list(cls, v):
        if v is not None and not isinstance(v, list):
            raise ValueError('ParticleBoundaryScrapingDiagnostic: data_list must be a list')
        return v

# ----------------------------
# Lab frame diagnostics
# ----------------------------


class PICMI_LabFrameFieldDiagnostic(_ClassWithInit):
    """
    Defines the electromagnetic field diagnostics in the lab frame.
    """
    grid: "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid | PICMI_CylindricalGrid" = Field(description="Grid object for the diagnostic")
    num_snapshots: int = Field(description="Number of lab frame snapshots to make")
    dt_snapshots: float = Field(description="Time between each snapshot in lab frame")
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
    write_dir: str | None = Field(default=None, description="Directory where data is to be written")
    parallelio: bool | None = Field(default=None, description="If set to True, field diagnostics are dumped in parallel")
    name: str | None = Field(default=None, description="Sets the base name for the diagnostic output files")
    
    @field_validator('data_list')
    @classmethod
    def _validate_data_list(cls, v):
        if v is not None and not isinstance(v, list):
            raise ValueError('LabFrameFieldDiagnostic: data_list must be a list')
        return v


class PICMI_LabFrameParticleDiagnostic(_ClassWithInit):
    """
    Defines the particle diagnostics in the lab frame.
    """
    grid: "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid | PICMI_CylindricalGrid" = Field(description="Grid object for the diagnostic")
    num_snapshots: int = Field(description="Number of lab frame snapshots to make")
    dt_snapshots: float = Field(description="Time between each snapshot in lab frame")
    species: "PICMI_Species | PICMI_MultiSpecies | list[PICMI_Species | PICMI_MultiSpecies] | None" = Field(
        default=None,
        description="Species instance or list of species instances. Species to write out. If not specified, all species are written. Note that the name attribute must be defined for the species."
    )
    data_list: list[str] | None = Field(
        default=None,
        description="The data to be written out. Possible values 'position', 'momentum', 'weighting'. Defaults to the output list of the implementing code."
    )
    time_start: float = Field(default=0., description="Time for the first snapshot in lab frame")
    write_dir: str | None = Field(default=None, description="Directory where data is to be written")
    parallelio: bool | None = Field(default=None, description="If set to True, particle diagnostics are dumped in parallel")
    name: str | None = Field(default=None, description="Sets the base name for the diagnostic output files")
    
    @field_validator('data_list')
    @classmethod
    def _validate_data_list(cls, v):
        if v is not None and not isinstance(v, list):
            raise ValueError('LabFrameParticleDiagnostic: data_list must be a list')
        return v